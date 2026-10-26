"""Discovery batch writer locks and verified, isolated executable snapshots.

Locks are held by open file descriptions and inherited by native children.
Closing the coordinator descriptor therefore does not unlock an orphan scan.
These checks do not hash raw FASTQ input or change scanner observations.
"""
from __future__ import annotations

from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import tempfile

VERIFIED_IDENTITY_FIELDS = ('device', 'inode', 'size_bytes', 'mtime_ns', 'ctime_ns', 'mode')
MAX_VERIFICATION_REPORT_BYTES = 256 * 1024 * 1024


def file_identity(path: Path) -> list[int]:
    info = Path(path).stat()
    return [info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns,
            info.st_ctime_ns, info.st_mode]


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


@contextmanager
def output_lock(directory: Path, name: str = '.batch.lock'):
    """Refuse a second writer without changing its state, summary or lock data."""
    directory = Path(directory)
    if directory.is_symlink():
        raise ValueError(f'Batch output directory must not be a symlink: {directory}')
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    if path.is_symlink():
        raise ValueError(f'Batch lock must not be a symlink: {path}')
    descriptor = os.open(path, os.O_CREAT | os.O_RDWR | getattr(os, 'O_NOFOLLOW', 0), 0o600)
    try:
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            raise ValueError(f'Batch lock is not a regular file: {path}')
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError(f'Batch output is locked by another writer: {directory}') from error
        # Do not issue LOCK_UN: a surviving native child shares this file
        # description and must keep the lock after its coordinator exits.
        yield descriptor
    finally:
        os.close(descriptor)


def assert_snapshot(snapshot: dict, *, hash_contents: bool = False) -> None:
    path = Path(snapshot['path'])
    if path.is_symlink() or path.resolve() != path or file_identity(path) != snapshot['identity']:
        raise ValueError(f'Frozen scanner changed during batch execution: {path}')
    if hash_contents and digest(path) != snapshot['sha256']:
        raise ValueError(f'Frozen scanner SHA256 changed during batch execution: {path}')
    if file_identity(path) != snapshot['identity']:
        raise ValueError(f'Frozen scanner changed while verifying: {path}')


def freeze_binary(output: Path, executable: dict) -> dict:
    """Copy the selected scanner into a content-addressed, read-only snapshot.

    The caller owns the batch output lock. No hardlink to the installed binary
    is used: replacement or modification of the original must not affect jobs.
    Existing corrupt snapshots are refused instead of overwritten silently.
    """
    source = Path(executable['path']).resolve(strict=True)
    before = file_identity(source)
    expected = executable['sha256']
    if not stat.S_ISREG(before[-1]) or digest(source) != expected or file_identity(source) != before:
        raise ValueError(f'Selected scanner changed before snapshot creation: {source}')
    root = Path(output).resolve() / '.batch-binaries'
    directory = root / expected
    if root.is_symlink() or directory.is_symlink():
        raise ValueError('Frozen scanner directories must not be symlinks')
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / 'trm'
    if destination.is_symlink():
        raise ValueError(f'Frozen scanner must not be a symlink: {destination}')
    if destination.exists():
        if not destination.is_file() or not os.access(destination, os.X_OK) or digest(destination) != expected:
            raise ValueError(f'Existing frozen scanner is corrupt or not executable: {destination}')
        if destination.stat().st_mode & 0o222:
            raise ValueError(f'Existing frozen scanner is not read-only: {destination}')
    else:
        descriptor, temporary = tempfile.mkstemp(prefix='.snapshot-', dir=directory)
        try:
            with os.fdopen(descriptor, 'wb') as target, source.open('rb') as original:
                shutil.copyfileobj(original, target, length=1024 * 1024)
                target.flush()
                os.fsync(target.fileno())
            if (file_identity(source) != before or digest(Path(temporary)) != expected
                    or digest(source) != expected or file_identity(source) != before):
                raise ValueError(f'Selected scanner changed while copying: {source}')
            os.chmod(temporary, 0o555)
            os.replace(temporary, destination)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    snapshot = {'path': str(destination), 'source_path': str(source), 'sha256': expected,
                'identity': file_identity(destination)}
    assert_snapshot(snapshot, hash_contents=True)
    return snapshot


def verified_file_identity(info) -> dict:
    """Schema shared with the standalone external-archive verifier (v2)."""
    return dict(device=info.st_dev, inode=info.st_ino, size_bytes=info.st_size,
                mtime_ns=info.st_mtime_ns, ctime_ns=info.st_ctime_ns, mode=info.st_mode)


@contextmanager
def _open_regular_nofollow(path):
    """Open an absolute local path without following any symlink component."""
    path = Path(os.path.abspath(path))
    directory = os.open(path.anchor, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        parts = path.parts[1:]
        if not parts:
            raise ValueError('Expected a regular file, not a filesystem root')
        for component in parts[:-1]:
            child = os.open(component, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                            dir_fd=directory)
            os.close(directory)
            directory = child
        descriptor = os.open(parts[-1], os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK,
                             dir_fd=directory)
    finally:
        os.close(directory)
    with os.fdopen(descriptor, 'rb') as handle:
        if not stat.S_ISREG(os.fstat(handle.fileno()).st_mode):
            raise ValueError(f'Verified input is not a regular file: {path}')
        yield handle


def _current_verified_identity(path):
    with _open_regular_nofollow(path) as handle:
        return verified_file_identity(os.fstat(handle.fileno()))


def _unique_json_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'Duplicate field in input verification report: {key}')
        result[key] = value
    return result


def load_input_verification(report_path: Path, paths: list[Path]) -> dict:
    """Bind selected local files to a completed, user-trusted MD5 receipt.

    No archive content is opened for reading here: only the receipt is read.
    Local inode/ctime checks catch ordinary modification and timestamp-preserving
    replacement, but are not cryptographic authentication of a supplied receipt.
    """
    report_path = Path(os.path.abspath(report_path))
    with _open_regular_nofollow(report_path) as handle:
        before = verified_file_identity(os.fstat(handle.fileno()))
        if before['size_bytes'] > MAX_VERIFICATION_REPORT_BYTES:
            raise ValueError('Input verification report exceeds 256 MiB')
        payload = handle.read(MAX_VERIFICATION_REPORT_BYTES + 1)
        if (len(payload) > MAX_VERIFICATION_REPORT_BYTES
                or before != verified_file_identity(os.fstat(handle.fileno()))):
            raise ValueError('Input verification report changed while reading')
    if _current_verified_identity(report_path) != before:
        raise ValueError('Input verification report pathname changed while reading')
    document = json.loads(payload, object_pairs_hook=_unique_json_object)
    if (not isinstance(document, dict) or document.get('schema_version') != 2
            or document.get('program') != 'trm-external-download-verification'
            or document.get('status') != 'complete'
            or document.get('input_identity_method') != 'local_posix_stat_v1'
            or not isinstance(document.get('completed_utc'), str) or not document['completed_utc']
            or not isinstance(document.get('files'), list)):
        raise ValueError('Input verification requires a completed schema-2 local verifier JSON; rerun scripts/verify_external_download.py after copying')
    selected = {str(Path(os.path.abspath(path))) for path in paths}
    bound = {}
    for row in document['files']:
        if not isinstance(row, dict):
            raise ValueError('Invalid file record in input verification report')
        name = row.get('absolute_path')
        if not isinstance(name, str) or name not in selected:
            continue
        identity = row.get('verified_file_identity')
        expected = row.get('expected_bytes')
        expected_md5, observed_md5 = row.get('expected_md5'), row.get('observed_md5')
        if (row.get('status') != 'valid' or type(expected) is not int or expected < 0
                or type(row.get('observed_bytes')) is not int or row['observed_bytes'] != expected
                or not isinstance(expected_md5, str) or not re.fullmatch('[0-9a-f]{32}', expected_md5)
                or observed_md5 != expected_md5
                or not isinstance(identity, dict) or set(identity) != set(VERIFIED_IDENTITY_FIELDS)
                or any(type(identity[key]) is not int or identity[key] < 0 for key in VERIFIED_IDENTITY_FIELDS)
                or identity['size_bytes'] != expected or not stat.S_ISREG(identity['mode'])):
            raise ValueError(f'Selected input lacks valid size/MD5/local identity verification: {name}')
        value = {'path': name, 'size_bytes': expected, 'md5': observed_md5, 'identity': identity}
        if name in bound and bound[name] != value:
            raise ValueError(f'Conflicting input verification records: {name}')
        bound[name] = value
    missing = selected - bound.keys()
    if missing:
        raise ValueError('Input verification report does not cover selected files: ' + ', '.join(sorted(missing)))
    context = {'report': {'path': str(report_path), 'sha256': hashlib.sha256(payload).hexdigest(),
                          'identity': before}, 'files': bound}
    assert_input_verification(context)
    return context


def assert_input_verification(context: dict | None, paths=None, *, hash_report=False) -> None:
    if context is None:
        return
    report = context['report']
    if _current_verified_identity(report['path']) != report['identity']:
        raise ValueError('Input verification report changed during batch execution')
    if hash_report:
        with _open_regular_nofollow(report['path']) as handle:
            if verified_file_identity(os.fstat(handle.fileno())) != report['identity']:
                raise ValueError('Input verification report changed before checksum verification')
            digest_ = hashlib.sha256()
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                digest_.update(block)
            if (digest_.hexdigest() != report['sha256']
                    or verified_file_identity(os.fstat(handle.fileno())) != report['identity']):
                raise ValueError('Input verification report changed during checksum verification')
        if _current_verified_identity(report['path']) != report['identity']:
            raise ValueError('Input verification report pathname changed during checksum verification')
    selected = context['files'] if paths is None else [str(Path(os.path.abspath(path))) for path in paths]
    for path in selected:
        if path not in context['files']:
            raise ValueError(f'Input verification report does not cover selected file: {path}')
        if _current_verified_identity(path) != context['files'][path]['identity']:
            raise ValueError(f'Verified input identity is stale or replaced; reverify on this filesystem: {path}')


def input_verification_binding(context: dict, paths: list[Path]) -> dict:
    assert_input_verification(context, paths)
    return {'schema_version': 1, 'method': 'trusted_local_md5_receipt_and_posix_identity',
            'receipt': {key: context['report'][key] for key in ('path', 'sha256')},
            'files': [context['files'][str(Path(os.path.abspath(path)))] for path in paths]}
