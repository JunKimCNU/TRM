#!/usr/bin/env python3
"""Verify newly imported archives against a TSV manifest, without moving files.

This verifies compressed bytes and MD5 only. It does not validate mate pairing,
read quality, library identity, or biological suitability and does not update the
project's analysis/download ledgers. Run it on a separate incoming directory.
Schema 2 JSON records the local file identity verified by that MD5 pass. The
optional batch receipt check can reuse it only while path/stat identity remains
unchanged; copied files must be verified again on the destination filesystem.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import sys

REQUIRED = ('file_id', 'run_accession', 'archive_name', 'relative_path', 'url',
            'expected_bytes', 'expected_md5', 'scientific_name', 'task_id')
OPTIONAL = ('role', 'representation')
REPORTS = ('verification.tsv', 'verification.json', 'verified_files.tsv')
BLOCK_BYTES = 8 * 1024 * 1024


class ValidationError(ValueError):
    pass


def validate_relative_path(value, archive_name):
    path = PurePosixPath(value)
    parts = value.split('/')
    if (not value or path.is_absolute() or '\\' in value or ':' in value
            or any(part in ('', '.', '..') for part in parts)
            or str(path) != value or any(ord(ch) < 32 for ch in value)):
        raise ValidationError('Unsafe relative_path: ' + repr(value))
    if archive_name != parts[-1] or '/' in archive_name or '\\' in archive_name:
        raise ValidationError('archive_name must equal the final relative_path component')
    if any(re.search(r'(^|[._-])(part|partial|failed|quarantine)([._-]|$)', part, re.I)
           for part in parts):
        raise ValidationError('Partial or quarantined paths are not importable: ' + value)


def read_manifest(path):
    raw = Path(path).read_bytes()
    try:
        text = raw.decode('utf-8-sig')
    except UnicodeDecodeError as exc:
        raise ValidationError('Manifest must be UTF-8 TSV') from exc
    reader = csv.DictReader(io.StringIO(text, newline=''), delimiter='\t')
    fields = reader.fieldnames or []
    if len(fields) != len(set(fields)) or not set(REQUIRED).issubset(fields):
        raise ValidationError('Missing required or duplicate manifest columns')
    rows, by_id, by_path = [], {}, {}
    for line, raw_row in enumerate(reader, 2):
        if None in raw_row or any(value is None for value in raw_row.values()):
            raise ValidationError('Malformed TSV row at line %d' % line)
        row = {key: raw_row.get(key, '') for key in REQUIRED + OPTIONAL}
        if any(not row[key] for key in REQUIRED):
            raise ValidationError('Empty required value at line %d' % line)
        if any('\x00' in value for value in row.values()):
            raise ValidationError('NUL in manifest at line %d' % line)
        validate_relative_path(row['relative_path'], row['archive_name'])
        if not re.fullmatch(r'[0-9]+', row['expected_bytes']):
            raise ValidationError('Invalid expected_bytes at line %d' % line)
        row['expected_bytes'] = int(row['expected_bytes'])
        row['expected_md5'] = row['expected_md5'].lower()
        if not re.fullmatch(r'[0-9a-f]{32}', row['expected_md5']):
            raise ValidationError('Invalid expected_md5 at line %d' % line)
        for index, key in [(by_id, row['file_id']), (by_path, row['relative_path'])]:
            if key in index and index[key] != row:
                raise ValidationError('Conflicting duplicate identity at line %d' % line)
        if row['file_id'] in by_id:
            continue  # Repeated identical rows refer to one archive.
        by_id[row['file_id']] = row
        by_path[row['relative_path']] = row
        rows.append(row)
    if not rows:
        raise ValidationError('Manifest contains no archives')
    return rows, {'path': str(Path(path).resolve()), 'bytes': len(raw),
                  'sha256': hashlib.sha256(raw).hexdigest()}


def open_regular_file(root, relative):
    """Open below a directory FD, refusing symlinks at every path component.

    O_NONBLOCK ensures a malicious or accidental FIFO cannot hang a verifier.
    Directory-relative opens also prevent a renamed path escaping via symlinks.
    """
    directory = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        parts = relative.split('/')
        for part in parts[:-1]:
            child = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                            dir_fd=directory)
            os.close(directory)
            directory = child
        descriptor = os.open(parts[-1], os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK,
                             dir_fd=directory)
    finally:
        os.close(directory)
    if not stat.S_ISREG(os.fstat(descriptor).st_mode):
        os.close(descriptor)
        raise ValidationError('Input is not a regular file')
    return os.fdopen(descriptor, 'rb')


def fingerprint(info):
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns, info.st_mode


def verified_file_identity(info):
    return dict(device=info.st_dev, inode=info.st_ino, size_bytes=info.st_size,
                mtime_ns=info.st_mtime_ns, ctime_ns=info.st_ctime_ns, mode=info.st_mode)


def verify_one(root, row):
    result = dict(row, absolute_path=str(root / row['relative_path']),
                  observed_bytes=None, observed_md5=None, status='error', reason='')
    try:
        with open_regular_file(root, row['relative_path']) as source:
            before = os.fstat(source.fileno())
            result['observed_bytes'] = before.st_size
            if before.st_size != row['expected_bytes']:
                result.update(status='mismatch', reason='size_mismatch; MD5 not computed')
                return result
            digest = hashlib.md5(usedforsecurity=False)
            while True:
                block = source.read(BLOCK_BYTES)
                if not block:
                    break
                digest.update(block)
            after = os.fstat(source.fileno())
            if fingerprint(before) != fingerprint(after):
                result['reason'] = 'File changed while being verified; retry after copying finishes'
                return result
            # A copy tool can replace a pathname while the old open descriptor
            # still refers to unchanged bytes. Detect that case as well.
            with open_regular_file(root, row['relative_path']) as current:
                if fingerprint(before) != fingerprint(os.fstat(current.fileno())):
                    result['reason'] = 'File pathname changed while being verified; retry after copying finishes'
                    return result
            result['observed_md5'] = digest.hexdigest()
            if digest.hexdigest() == row['expected_md5']:
                result.update(status='valid', reason='size_and_md5_match',
                              verified_file_identity=verified_file_identity(before))
            else:
                result.update(status='mismatch', reason='md5_mismatch')
    except FileNotFoundError:
        result.update(status='missing', reason='Expected final archive is absent; .part files are not used')
    except (OSError, ValidationError) as exc:
        result['reason'] = '%s: %s' % (type(exc).__name__, exc)
    return result


def run(manifest, input_root, output, jobs=4, require_all=False):
    if not 1 <= jobs <= 144:
        raise ValidationError('--jobs must be between 1 and 144')
    rows, manifest_pin = read_manifest(manifest)
    root = Path(input_root).resolve(strict=True)
    if not root.is_dir():
        raise ValidationError('--input-root must be an existing directory')
    raw_output = Path(output)
    if raw_output.is_symlink():
        raise ValidationError('--output must not be a symlink')
    target = raw_output.resolve()
    if target == root or root in target.parents:
        raise ValidationError('--output must be outside --input-root to protect input archives')
    if target.exists() and (not target.is_dir() or any(target.iterdir())):
        raise ValidationError('--output must be a new or empty directory; existing reports are never overwritten')
    target.mkdir(parents=True, exist_ok=True)
    # Reserve all report names before hashing. Exclusive creation also rejects
    # another verifier racing to write into the same output directory.
    reservations = []
    try:
        for name in REPORTS:
            descriptor = os.open(target / name, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
            reservations.append((target / name, descriptor))
    except BaseException:
        for path, descriptor in reservations:
            os.close(descriptor)
            path.unlink()
        raise
    started = datetime.now(timezone.utc).isoformat()
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        results = list(executor.map(lambda row: verify_one(root, row), rows))
    counts = {status: sum(item['status'] == status for item in results)
              for status in ('valid', 'missing', 'mismatch', 'error')}
    failed = bool(counts['mismatch'] or counts['error'] or (require_all and counts['missing']))
    receipt = {'schema_version': 2, 'program': 'trm-external-download-verification',
               'status': 'complete', 'started_utc': started,
               'completed_utc': datetime.now(timezone.utc).isoformat(),
               'manifest': manifest_pin, 'input_root': str(root), 'jobs': jobs,
               'require_all': require_all, 'counts': counts, 'total_files': len(results),
               'all_manifest_files_present_and_valid': counts['valid'] == len(results),
               'exit_code': 1 if failed else 0, 'files_moved_or_promoted': False,
               'project_ledgers_modified': False,
               'input_identity_method': 'local_posix_stat_v1',
               'receipt_reuse_scope': 'User-trusted local verification report; path/device/inode/size/mtime/ctime/mode must still match. Copying to another path or filesystem requires new verification. Not a signed authenticity assertion.',
               'validation_scope': 'Compressed archive bytes and MD5 only; not mate pairing, read quality, library identity, or biological suitability',
               'files': results}
    fields = REQUIRED + OPTIONAL + ('absolute_path', 'observed_bytes', 'observed_md5', 'status', 'reason')
    # Write via the already reserved descriptors, never reopen report pathnames.
    for (path, descriptor), content in zip(reservations, [results, receipt,
            [item for item in results if item['status'] == 'valid']]):
        with os.fdopen(descriptor, 'w', encoding='utf-8', newline='') as handle:
            if path.name.endswith('.json'):
                json.dump(content, handle, indent=2, ensure_ascii=False)
                handle.write('\n')
            else:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
                writer.writeheader()
                writer.writerows(content)
    return receipt


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True, help='UTF-8 TSV from the external download package')
    parser.add_argument('--input-root', required=True, help='Directory containing only newly imported files')
    parser.add_argument('--output', required=True, help='New or empty report directory outside the input root')
    parser.add_argument('--jobs', type=int, default=4, help='Concurrent archive checks (default: 4; maximum: 144)')
    parser.add_argument('--require-all', action='store_true', help='Treat missing manifest files as failure')
    args = parser.parse_args(argv)
    try:
        receipt = run(args.manifest, args.input_root, args.output, args.jobs, args.require_all)
    except (OSError, ValidationError) as exc:
        print('Verification could not start: %s' % exc, file=sys.stderr)
        return 2
    print(json.dumps({'counts': receipt['counts'], 'exit_code': receipt['exit_code'],
                      'output': str(Path(args.output).resolve())}))
    return receipt['exit_code']


if __name__ == '__main__':
    sys.exit(main())
