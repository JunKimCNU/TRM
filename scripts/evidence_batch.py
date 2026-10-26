#!/usr/bin/env python3
"""Collect full-input junction evidence without repeating all-period discovery.

The saved discovery result remains authoritative for the original all-k scan.
Target selection is explicit, bounded and recorded; nonselected families are
never described as having been tested for junction support.
FASTQ paths are relative to the manifest directory unless absolute. path1 is
required; missing files are retained for per-sample failure reporting.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import time
import tempfile
import uuid


EXECUTION_VERSION = '1.2.0'


def sha256(path):
    """Hash a file without loading the whole artifact into coordinator memory."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def canonical(motif):
    motif = motif.upper()
    for size in range(1, len(motif) + 1):
        if len(motif) % size == 0 and motif[:size] * (len(motif) // size) == motif:
            motif = motif[:size]
            break
    rc = motif.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    return min(s[i:] + s[:i] for s in (motif, rc) for i in range(len(s)))


def read_tsv(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def index_samples(rows, label):
    indexed = {}
    for row in rows:
        sample = row.get('sample_id', '')
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}', sample) or sample in indexed:
            raise ValueError(f'{label}: invalid or duplicate sample_id {sample!r}')
        indexed[sample] = row
    if not indexed:
        raise ValueError(f'{label}: no samples')
    return indexed


def normalize_manifest_inputs(rows, manifest_path):
    """Validate scheduled file inputs without opening FASTQs or dropping rows."""
    parent = Path(manifest_path).resolve().parent
    normalized = []
    for row in rows:
        sample = row.get('sample_id', '')
        if not isinstance(row.get('path1'), str) or not row['path1'].strip():
            raise ValueError(f'manifest sample {sample}: path1 must not be empty')
        item, paths, identities = dict(row), set(), set()
        for key in ('path1', 'path2'):
            value = row.get(key)
            if value in (None, ''):
                item[key] = ''
                continue
            if not isinstance(value, str):
                raise ValueError(f'manifest sample {sample}: {key} must be a file path')
            path = Path(value).expanduser()
            path = (path if path.is_absolute() else parent / path).resolve()
            if path in paths:
                raise ValueError(f'manifest sample {sample}: path1 and path2 refer to the same file')
            paths.add(path)
            try:
                info = path.stat()
            except FileNotFoundError:
                pass  # A missing input must remain a recorded per-sample failure.
            else:
                if not stat.S_ISREG(info.st_mode):
                    raise ValueError(f'manifest sample {sample}: {key} must be a regular file, not a directory or stream: {path}')
                identity = (info.st_dev, info.st_ino)
                if identity in identities:
                    raise ValueError(f'manifest sample {sample}: path1 and path2 refer to the same file')
                identities.add(identity)
            item[key] = str(path)
        normalized.append(item)
    return normalized


def write_json(path, data, *, compact=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8') as stream:
            if compact:
                # One-shot encoding uses the standard library's C encoder.
                # No fields are dropped; callers opt in to whitespace removal.
                stream.write(json.dumps(data, ensure_ascii=False, allow_nan=False,
                                        separators=(',', ':')))
            else:
                json.dump(data, stream, indent=2, ensure_ascii=False, allow_nan=False)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def select_targets(report, known, max_novel=96):
    rows = {canonical(r['motif']): r for r in
            report.get('candidates', []) + report.get('rejected_candidates', [])}
    # Include unknown families only by generic observational criteria, never by
    # proximity to TTAGGG or a motif-specific blacklist.
    eligible = [r for family, r in rows.items() if family not in known
                and 3 <= len(family) <= 32 and r.get('reads', 0) >= 5
                and r.get('full_repeat_reads', 0) >= 2]
    eligible.sort(key=lambda r: (-r.get('repeat_bases', 0), -r.get('reads', 0), r['motif']))
    selected = eligible[:max_novel]
    motifs = sorted(known | {canonical(r['motif']) for r in selected})
    if len(motifs) > 256:
        raise ValueError('target family count exceeds evidence scanner limit 256')
    return {'motifs': motifs, 'selection': {
        'known_family_count': len(known), 'novel_screening_candidates': len(eligible),
        'novel_targets': len(selected), 'novel_target_limit': max_novel,
        'novel_candidates_not_targeted': len(eligible) - len(selected),
        'rule': 'All catalogue families <=32 bp plus unknown 3..32-bp families with >=5 reads and >=2 full-repeat records, top repeat-bases up to the stated cap.',
        'limitations': 'A target-limited second pass cannot rule out an unselected, low-count, long or irregular novel repeat.',
    }}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sample_shard(sample_id, shard_count):
    return int.from_bytes(hashlib.sha256(sample_id.encode('utf-8')).digest(), 'big') % shard_count


def available_cpus():
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)


def positive_int(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError('must be a positive integer')
    return number


def positive_seconds(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError('must be finite and greater than zero')
    return number


def native_options(args):
    """Only nondefault additions enter legacy-compatible sample identities."""
    from native_diagnostic_options import native_diagnostic_options
    result = native_diagnostic_options(args)
    if getattr(args, 'memory_limit_mib', None) is not None:
        result += ['--memory-limit-mib', str(args.memory_limit_mib)]
    if getattr(args, 'paired_layout', None):
        result += ['--paired-layout', args.paired_layout]
    if getattr(args, 'flank_period_annotation', None) == 'on':
        result += ['--flank-period-annotation', 'on']
    return result


def input_metadata(paths):
    result = []
    for path in paths:
        path = Path(path).resolve()
        info = path.stat()
        if not stat.S_ISREG(info.st_mode):
            raise ValueError(f'input is not a regular file: {path}')
        result.append({'path': str(path), 'bytes': info.st_size, 'mtime_ns': info.st_mtime_ns})
    return result


def file_fingerprint(path):
    info = Path(path).stat()
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns)


def binary_digest(args, binary):
    cached = getattr(args, '_frozen_binary', None)
    if cached is not None:
        if str(Path(binary).resolve()) != cached['path'] or file_fingerprint(binary) != cached['fingerprint']:
            raise ValueError('Frozen binary metadata changed during this execution')
        return cached['sha256']
    return sha256(binary)  # Backwards-compatible direct run_one callers.


def validate_summary(report, label):
    if not isinstance(report, dict) or not isinstance(report.get('summary'), dict):
        raise ValueError(f'{label}: result lacks summary object')
    for name in ('records', 'bases'):
        value = report['summary'].get(name)
        if type(value) is not int or value < 0:
            raise ValueError(f'{label}: summary requires nonnegative integer {name}')
    if 'sampled_prefix' in report['summary'] and type(report['summary']['sampled_prefix']) is not bool:
        raise ValueError(f'{label}: sampled_prefix must be boolean')
    if 'program' in report and report['program'] != 'trm':
        raise ValueError(f'{label}: unexpected program')
    if 'version' in report and (not isinstance(report['version'], str) or not report['version']):
        raise ValueError(f'{label}: invalid version metadata')


def validate_report_inputs(report, inputs, label):
    """Validate supplied native provenance; absent legacy fields are unmeasured."""
    if 'inputs' not in report:
        return 'legacy_input_metadata_absent'
    given = report['inputs']
    if not isinstance(given, list) or len(given) != len(inputs):
        raise ValueError(f'{label}: input count differs from scheduled files')
    for observed, expected in zip(given, inputs):
        if not isinstance(observed, dict):
            raise ValueError(f'{label}: invalid input metadata')
        if observed.get('path') is not None:
            if not isinstance(observed['path'], str) or not observed['path']:
                raise ValueError(f'{label}: invalid input path')
            if str(Path(observed['path']).expanduser().resolve()) != expected['path']:
                raise ValueError(f'{label}: input path differs from scheduled file')
        for field, expected_value in [('bytes', expected['bytes']), ('mtime_ns', expected['mtime_ns']),
                ('modified_unix_seconds', expected['mtime_ns'] // 1_000_000_000)]:
            if observed.get(field) is not None and observed[field] != expected_value:
                raise ValueError(f'{label}: input {field} differs from scheduled file')
    return 'supplied_input_metadata_validated'


def validate_evidence(evidence, discovery, targets, inputs, args):
    from native_diagnostic_options import validate_diagnostic_output
    validate_diagnostic_output(evidence, args)
    validate_summary(evidence, 'evidence')
    if evidence.get('analysis_kind', 'targeted_evidence') != 'targeted_evidence':
        raise ValueError('Evidence output is not targeted_evidence')
    rows = evidence.get('candidate_evidence')
    if not isinstance(rows, list) or any(not isinstance(row, dict) or not isinstance(row.get('motif'), str) for row in rows):
        raise ValueError('Evidence output lacks candidate_evidence motifs')
    motifs = [row['motif'] for row in rows]
    if len(motifs) != len(set(motifs)) or set(motifs) != set(targets['motifs']):
        raise ValueError('Target/output family set differs')
    if 'target_motifs' in evidence and (not isinstance(evidence['target_motifs'], list)
            or set(evidence['target_motifs']) != set(targets['motifs'])):
        raise ValueError('Target/output target_motifs differs')
    if not args.max_records:
        if evidence['summary'].get('sampled_prefix', False):
            raise ValueError('Full-input evidence unexpectedly reports a sampled prefix')
        for field in ('records', 'bases'):
            if evidence['summary'][field] != discovery['summary'][field]:
                raise ValueError(f'{field} differs from original full scan')
    elif evidence['summary']['records'] > args.max_records:
        raise ValueError('Evidence exceeds requested record prefix')
    validate_report_inputs(evidence, inputs, 'evidence')
    if getattr(args, 'paired_layout', None):
        paired = evidence.get('paired_end')
        if not isinstance(paired, dict) or paired.get('measured') is not True or paired.get('layout') != args.paired_layout:
            raise ValueError('Evidence paired-end provenance is missing or differs from requested options')
    if getattr(args, 'flank_period_annotation', None) == 'on':
        junction = evidence.get('junction_options')
        if not isinstance(junction, dict) or junction.get('flank_period_annotation') is not True:
            raise ValueError('Evidence flank annotation provenance is missing or differs from requested option')


def safe_sample_directory(output, sample):
    output = Path(output).resolve()
    directory = output / sample
    if directory.is_symlink() or directory.resolve() != directory:
        raise ValueError(f'sample output directory must not be a symlink: {directory}')
    if directory.exists() and not directory.is_dir():
        raise ValueError(f'sample output is not a directory: {directory}')
    for name in ('status.json', 'targets.json', 'evidence.json', 'evidence.tsv', 'evidence.rejected.tsv',
                 'stdout.log', 'stderr.log', 'resource.tsv', 'attempts'):
        if (directory / name).is_symlink():
            raise ValueError(f'managed output must not be a symlink: {directory / name}')
    return directory


def persist_status(output, sample, status):
    """A rejected sample symlink gets a failure receipt outside that directory."""
    try:
        directory = safe_sample_directory(output, sample)
    except ValueError:
        directory = Path(output) / '.failures'
        if directory.is_symlink():
            raise ValueError('failure receipt directory must not be a symlink')
        path = directory / f'{sample}.json'
    else:
        path = directory / 'status.json'
    write_json(path, status)
    return path


def resources(path):
    lines = Path(path).read_text().strip().splitlines()
    values = lines[-1].split('\t') if lines else []
    if len(values) not in (4, 5):
        raise ValueError(f'Invalid GNU time resource report: {path}')
    wall, user, system = map(float, values[:3])
    rss = int(values[3])
    if any(not math.isfinite(x) or x < 0 for x in (wall, user, system)) or rss < 0:
        raise ValueError(f'Invalid GNU time resource values: {path}')
    result = {'resource_wall_seconds': wall, 'user_cpu_seconds': user,
              'system_cpu_seconds': system, 'max_rss_kb': rss}
    if len(values) == 5:
        result['resource_exit_code'] = int(values[4])
    return result


def archive_previous(directory, attempt):
    previous = attempt / 'previous_published'
    for name in ('status.json', 'targets.json', 'evidence.json', 'evidence.tsv', 'evidence.rejected.tsv',
                 'stdout.log', 'stderr.log', 'resource.tsv'):
        old = directory / name
        if old.exists():
            if not old.is_file() or old.is_symlink():
                raise ValueError(f'Unexpected managed output: {old}')
            previous.mkdir(exist_ok=True)
            os.replace(old, previous / name)


class SourceFailure(ValueError):
    pass


def run_one(item, args, known, frozen_binary):
    """Execute one isolated attempt, or return an unchanged valid legacy status."""
    sample = item['sample_id']
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}', sample):
        raise ValueError(f'unsafe sample ID: {sample!r}')
    output_root, source_root = Path(args.output).resolve(), Path(args.source).resolve()
    if output_root == source_root or output_root in source_root.parents or source_root in output_root.parents:
        raise ValueError('Evidence output and source directories must be separate and not nested')
    started = time.monotonic()
    status = {'sample_id': sample, 'status': 'failed', 'execution_version': EXECUTION_VERSION,
              'execution_id': getattr(args, 'execution_id', None) or uuid.uuid4().hex,
              'started_utc': utc_now(), 'subprocess_started': False}
    attempt = None
    try:
        directory = safe_sample_directory(args.output, sample)
        directory.mkdir(parents=True, exist_ok=True)
        status_path = directory / 'status.json'
        source = Path(args.source) / sample / 'result.json'
        paths = [Path(item[key]).resolve() for key in ('path1', 'path2') if item.get(key)]
        status['inputs'] = list(map(str, paths))
        if not isinstance(item.get('path1'), str) or not item['path1'].strip():
            raise ValueError('path1 must not be empty; no FASTQ input was scheduled')
        if getattr(args, 'paired_layout', None) and (len(paths) != 2 or args.max_records % 2):
            raise ValueError('paired FR requires exactly two input files and an even record prefix')
        for path in [source, *paths]:
            if path.resolve() == Path(args.output).resolve() or Path(args.output).resolve() in path.resolve().parents:
                raise ValueError('Evidence output must be disjoint from source and input files')
        try:
            source_stamp = file_fingerprint(source)
            with source.open() as handle:
                report = json.load(handle)
            source_digest = sha256(source)
            if file_fingerprint(source) != source_stamp:
                raise ValueError('Source changed while reading discovery JSON')
            validate_summary(report, 'source discovery')
            if report.get('analysis_kind', 'de_novo_discovery') != 'de_novo_discovery':
                raise ValueError('Source is not a de novo discovery result')
            if report['summary'].get('sampled_prefix', False):
                raise ValueError('Source discovery is a sampled prefix, not a completed full-input discovery')
            for key in ('candidates', 'rejected_candidates'):
                if key in report and not isinstance(report[key], list):
                    raise ValueError(f'Source {key} is not an array')
        except (OSError, ValueError, TypeError) as error:
            raise SourceFailure(f'Original full discovery did not finish or is invalid: {error}; no biological nondetection inference.') from error
        targets = select_targets(report, known, args.max_novel_targets)
        inputs = input_metadata(paths)
        status['source_input_provenance'] = validate_report_inputs(report, inputs, 'source discovery')
        identity = {'source_sha256': source_digest, 'inputs': inputs,
                    'binary_sha256': binary_digest(args, frozen_binary), 'targets': targets,
                    'threads': args.threads, 'max_records': args.max_records}
        options = native_options(args)
        if options:
            identity['native_options'] = options
        status.update(identity=identity, source_result=str(source.resolve()))
        if args.resume and status_path.is_file():
            try:
                with status_path.open() as handle:
                    old = json.load(handle)
                output = directory / 'evidence.json'
                if (old.get('status') == 'complete' and old.get('identity') == identity
                        and output.is_file() and not output.is_symlink()
                        and sha256(output) == old.get('result_sha256')):
                    with output.open() as handle:
                        evidence = json.load(handle)
                    validate_evidence(evidence, report, targets, inputs, args)
                    if input_metadata(paths) != inputs or sha256(source) != identity['source_sha256']:
                        raise ValueError('Source/input changed during resume validation')
                    if binary_digest(args, frozen_binary) != identity['binary_sha256']:
                        raise ValueError('Binary changed during resume validation')
                    return old
            except (OSError, ValueError, TypeError, KeyError, AttributeError):
                pass  # Incomplete, corrupt or changed artifacts require a new attempt.
        attempts = directory / 'attempts'
        attempts.mkdir(exist_ok=True)
        attempt = Path(tempfile.mkdtemp(prefix='attempt-', dir=attempts))
        status['attempt_directory'] = str(attempt.resolve())
        archive_previous(directory, attempt)
        write_json(attempt / 'targets.json', targets)
        write_json(directory / 'targets.json', targets)
        command = [str(frozen_binary), 'evidence', *map(str, paths), '--targets', str(attempt / 'targets.json'),
                   '--output', str(attempt / 'evidence'), '--threads', str(args.threads)]
        if args.max_records:
            command += ['--max-records', str(args.max_records)]
        command += options
        status.update(command=command, status='running', subprocess_started=True)
        persist_status(args.output, sample, status)
        write_json(attempt / 'status.json', status)
        timed = ['/usr/bin/time', '-f', '%e\t%U\t%S\t%M\t%x', '-o', str(attempt / 'resource.tsv'), *command]
        scan_started = time.monotonic()
        with (attempt / 'stdout.log').open('w') as stdout, (attempt / 'stderr.log').open('w') as stderr:
            completed = subprocess.run(timed, stdout=stdout, stderr=stderr, check=False)
        status.update(wall_seconds=time.monotonic() - scan_started, exit_code=completed.returncode, status='failed')
        status.update(resources(attempt / 'resource.tsv'))
        if status.get('resource_exit_code', completed.returncode) != completed.returncode:
            raise ValueError('GNU time exit code differs from subprocess exit code')
        if completed.returncode:
            raise RuntimeError((attempt / 'stderr.log').read_text()[-4000:] or f'evidence exited {completed.returncode}')
        generated = attempt / 'evidence.json'
        with generated.open() as handle:
            evidence = json.load(handle)
        validate_evidence(evidence, report, targets, inputs, args)
        if input_metadata(paths) != inputs:
            raise ValueError('Input file metadata changed during scan')
        if sha256(source) != identity['source_sha256']:
            raise ValueError('Source discovery changed during scan')
        if binary_digest(args, frozen_binary) != identity['binary_sha256']:
            raise ValueError('Binary changed during scan')
        with (attempt / 'targets.json').open() as handle:
            if json.load(handle) != targets:
                raise ValueError('Target definition changed during scan')
        # Only newly generated, verified artifacts can become the published result.
        for name in ('evidence.json', 'evidence.tsv', 'evidence.rejected.tsv', 'stdout.log', 'stderr.log', 'resource.tsv'):
            artifact = attempt / name
            if artifact.exists():
                if not artifact.is_file() or artifact.is_symlink():
                    raise ValueError(f'Unexpected generated artifact: {artifact}')
                descriptor, temporary = tempfile.mkstemp(prefix='.publish-', dir=directory)
                os.close(descriptor)
                try:
                    shutil.copy2(artifact, temporary)
                    os.replace(temporary, directory / name)
                finally:
                    if os.path.exists(temporary):
                        os.unlink(temporary)
        output = directory / 'evidence.json'
        full = not bool(args.max_records) and not evidence['summary'].get('sampled_prefix', False)
        status.update(status='complete', result_json=str(output.resolve()), result_sha256=sha256(output),
                      records=evidence['summary']['records'], bases=evidence['summary']['bases'],
                      input_bytes=sum(row['bytes'] for row in inputs) if full else 0,
                      input_artifact_bytes_total=sum(row['bytes'] for row in inputs),
                      compressed_bytes_processed_measured=full,
                      input_byte_scope='complete_input_artifact_bytes' if full else 'prefix_processed_bytes_unknown')
    except Exception as error:
        status.update(status='source_failed' if isinstance(error, SourceFailure) else 'failed', reason=str(error))
        if attempt is None:
            # A failure before native launch must not leave an old published
            # result next to a new failed state for consumers that inspect paths.
            try:
                directory = safe_sample_directory(args.output, sample)
                if any((directory / name).exists() for name in ('evidence.json', 'evidence.tsv', 'evidence.rejected.tsv')):
                    attempts = directory / 'attempts'
                    attempts.mkdir(exist_ok=True)
                    attempt = Path(tempfile.mkdtemp(prefix='failed-before-launch-', dir=attempts))
                    status['attempt_directory'] = str(attempt.resolve())
                    archive_previous(directory, attempt)
            except (OSError, ValueError):
                pass  # The failure receipt may have to use the safe side directory.
    status.setdefault('wall_seconds', time.monotonic() - started)
    status['finished_utc'] = utc_now()
    if attempt is not None:
        write_json(attempt / 'status.json', status)
        # Failed attempts keep their logs in the immutable attempt directory.
    persist_status(args.output, sample, status)
    return status


def aggregate(rows):
    return {'samples': len(rows), 'records': sum(row.get('records') or 0 for row in rows),
            'bases': sum(row.get('bases') or 0 for row in rows),
            'input_bytes': sum(row.get('input_bytes') or 0 for row in rows),
            'sum_sample_wall_seconds': sum(row.get('wall_seconds') or 0 for row in rows),
            'user_cpu_seconds': sum(row.get('user_cpu_seconds') or 0 for row in rows),
            'system_cpu_seconds': sum(row.get('system_cpu_seconds') or 0 for row in rows),
            'max_sample_rss_kb': max((row.get('max_rss_kb') or 0 for row in rows), default=0),
            'resource_measurements': sum(row.get('user_cpu_seconds') is not None for row in rows)}


def metrics_report(statuses, args, started, total, finished, checkpoint_count, max_inflight):
    successful = [row for row in statuses if row['status'] == 'complete']
    scanned = [row for row in successful if row.get('execution_id') == args.execution_id]
    resumed = [row for row in successful if row.get('execution_id') != args.execution_id]
    attempted = [row for row in statuses if row.get('execution_id') == args.execution_id and row.get('subprocess_started')]
    historical, current, reused, attempts = map(aggregate, (successful, scanned, resumed, attempted))
    elapsed = time.monotonic() - started
    measured_bytes = bool(scanned) and all(row.get('compressed_bytes_processed_measured', not args.max_records) for row in scanned)
    failed = sum(row['status'] == 'failed' for row in statuses)
    source_failed = sum(row['status'] == 'source_failed' for row in statuses)
    report = {'execution_version': EXECUTION_VERSION, 'execution_id': args.execution_id,
              'wall_seconds': elapsed, 'sample_count': len(statuses), 'total_samples': total,
              'completed': len(successful), 'failed': failed, 'source_failed': source_failed,
              'records': historical['records'], 'bases': historical['bases'], 'input_bytes': historical['input_bytes'],
              'sum_sample_wall_seconds': historical['sum_sample_wall_seconds'],
              'sum_user_cpu_seconds': historical['user_cpu_seconds'], 'sum_system_cpu_seconds': historical['system_cpu_seconds'],
              'max_sample_rss_kb': historical['max_sample_rss_kb'],
              'legacy_counter_scope': 'completed sample artifacts, including resumed historical measurements; not current execution resource consumption',
              'compressed_MB_per_second': current['input_bytes'] / 1e6 / elapsed if measured_bytes else None,
              'compressed_throughput_scope': 'successful scans performed in this execution only; null for no scans or unknown prefix bytes; input artifact bytes, not physical OS I/O',
              'current_execution': {'scanned': current, 'resumed_artifacts': reused, 'attempted': attempts,
                  'resumed_native_cpu_seconds_this_execution': 0, 'resumed_input_bytes_scanned_this_execution': 0,
                  'attempt_resource_scope': 'native process resources of this execution, including failed attempts; coordinator resources excluded'},
              'scanned_samples': len(scanned), 'resumed_samples': len(resumed),
              'scanned_bases': current['bases'], 'scanned_input_bytes': current['input_bytes'],
              'all_succeeded': finished and len(successful) == total, 'finished': finished,
              'checkpoint_count': checkpoint_count, 'max_inflight': max_inflight,
              'schedule': args.schedule, 'jobs': args.jobs, 'threads': args.threads, 'cpu_budget': args.cpu_budget,
              'record_resources': True, 'checkpoint_seconds': args.checkpoint_seconds,
              'resource_plan': getattr(args, '_resource_plan', None),
              'paired_layout': args.paired_layout, 'flank_period_annotation': args.flank_period_annotation == 'on',
              'memory_limit_mib_per_scan': args.memory_limit_mib,
              'configured_jobs_times_memory_limit_mib': args.jobs * args.memory_limit_mib if args.memory_limit_mib else None,
              'memory_limit_note': 'Per-process virtual-address-space ceilings; their sum is not measured aggregate RSS and excludes coordinator/other programs.',
              'shard_index': args.shard_index, 'shard_count': args.shard_count,
              'shard_algorithm': 'sha256_utf8_sample_id_big_endian_modulo_count' if args.shard_count else None}
    return report


def run_scheduled(manifest, args, known, binary):
    """At most jobs live futures; sample states commit immediately, summaries coalesce."""
    def scheduled_bytes(row):
        total = 0
        for key in ('path1', 'path2'):
            if row.get(key):
                try:
                    total += Path(row[key]).stat().st_size
                except OSError:
                    pass
        return total
    scheduled = sorted(manifest, key=lambda row: -scheduled_bytes(row)) if args.schedule == 'largest-first' else manifest
    order = {row['sample_id']: i for i, row in enumerate(manifest)}
    started, statuses = time.monotonic(), []
    submitted = checkpoint_count = max_inflight = 0

    def checkpoint(pending, finished=False):
        nonlocal checkpoint_count
        checkpoint_count += 1
        ordered = sorted(statuses, key=lambda row: order[row['sample_id']])
        report = metrics_report(ordered, args, started, len(manifest), finished, checkpoint_count, max_inflight)
        write_json(args.output / 'metrics.json', report)
        samples = [{**{key: row.get(key) for key in ('sample_id', 'status', 'records', 'bases', 'reason', 'result_json', 'wall_seconds',
                            'input_bytes', 'user_cpu_seconds', 'system_cpu_seconds', 'max_rss_kb')},
                    'execution_status': 'resumed' if row['status'] == 'complete' and row.get('execution_id') != args.execution_id
                                        else 'scanned' if row['status'] == 'complete' else row['status']} for row in ordered]
        write_json(args.output / 'summary.json', {**report, 'samples': samples})
        running = sorted(sample for future, sample in pending.items() if future.running())
        write_json(args.output / 'progress.json', {
            'execution_version': EXECUTION_VERSION, 'execution_id': args.execution_id,
            'resource_plan': getattr(args, '_resource_plan', None),
            'updated_utc': utc_now(), 'finished': len(statuses), 'execution_finished': finished,
            'total': len(manifest), 'wall_seconds': report['wall_seconds'],
            'status_counts': {key: sum(row['status'] == key for row in statuses) for key in ('complete', 'source_failed', 'failed')},
            'submitted': submitted, 'inflight': len(pending), 'running': running, 'running_count': len(running),
            'pending': len(manifest) - len(statuses) - len(running), 'checkpoint_count': checkpoint_count,
            'max_inflight': max_inflight, 'resumed_samples': report['resumed_samples']})
        return report

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = {}
        def replenish():
            nonlocal submitted, max_inflight
            while submitted < len(scheduled) and len(pending) < args.jobs:
                row = scheduled[submitted]
                pending[pool.submit(run_one, row, args, known, binary)] = row['sample_id']
                submitted += 1
            max_inflight = max(max_inflight, len(pending))
        replenish()
        checkpoint(pending)
        next_checkpoint = time.monotonic() + args.checkpoint_seconds
        while pending:
            done, _ = wait(pending, timeout=min(10.0, max(0, next_checkpoint - time.monotonic())), return_when=FIRST_COMPLETED)
            for future in done:
                sample = pending.pop(future)
                try:
                    status = future.result()
                except Exception as error:
                    status = {'sample_id': sample, 'status': 'failed', 'reason': str(error),
                              'execution_id': args.execution_id, 'execution_version': EXECUTION_VERSION}
                    persist_status(args.output, sample, status)
                statuses.append(status)
                print(f"{len(statuses)}/{len(manifest)} {sample} {status['status']}", flush=True)
            replenish()
            if time.monotonic() >= next_checkpoint:
                checkpoint(pending)
                next_checkpoint = time.monotonic() + args.checkpoint_seconds
    return checkpoint({}, finished=True)


def freeze_binary(args):
    binary = args.output / 'trm-evidence'
    if binary.is_symlink():
        raise ValueError('Frozen output binary must not be a symlink')
    digest = sha256(args.binary)
    if binary.exists() and sha256(binary) != digest:
        archive = args.output / '.binaries' / sha256(binary)
        if (args.output / '.binaries').is_symlink() or archive.is_symlink():
            raise ValueError('Binary archive must not be a symlink')
        archive.mkdir(parents=True, exist_ok=True)
        archived = archive / 'trm-evidence'
        if archived.is_symlink():
            raise ValueError('Archived binary must not be a symlink')
        if not archived.exists():
            descriptor, temporary = tempfile.mkstemp(prefix='.binary-', dir=archive)
            os.close(descriptor)
            try:
                shutil.copy2(binary, temporary)
                os.replace(temporary, archived)
            finally:
                if os.path.exists(temporary):
                    os.unlink(temporary)
        elif sha256(archived) != sha256(binary):
            raise ValueError('Archived binary hash conflicts with its directory')
    if not binary.exists() or sha256(binary) != digest:
        descriptor, temporary = tempfile.mkstemp(prefix='.binary-', dir=args.output)
        os.close(descriptor)
        try:
            shutil.copy2(args.binary, temporary)
            if sha256(temporary) != digest or sha256(args.binary) != digest:
                raise ValueError('Binary changed while freezing a snapshot')
            os.replace(temporary, binary)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    args._frozen_binary = {'path': str(binary.resolve()), 'sha256': digest, 'fingerprint': file_fingerprint(binary)}
    return binary


def main():
    from runtime_binary import resolve_binary
    from runner_resources import add_arguments as add_resource_arguments, preflight
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=Path('manifests/existing_fastq.tsv'))
    parser.add_argument('--source', type=Path, default=Path('results/wgs_taxonomy_v030'))
    parser.add_argument('--catalog', type=Path, default=Path('data/telomere_motifs.tsv'))
    parser.add_argument("--binary", help="Explicit scanner override; default: TRM_BINARY or verified bin/trm")
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--jobs', type=positive_int, default=16)
    parser.add_argument('--threads', type=positive_int, default=2)
    parser.add_argument('--cpu-budget', type=positive_int, default=available_cpus(),
                        help='maximum concurrent scanner workers, excluding helpers (default: CPU affinity count)')
    add_resource_arguments(parser)
    parser.add_argument('--memory-limit-mib', type=positive_int)
    parser.add_argument('--paired-layout', choices=('fr',))
    parser.add_argument('--flank-period-annotation', choices=('on', 'off'), default='off')
    parser.add_argument('--checkpoint-seconds', type=positive_seconds, default=10.0)
    parser.add_argument('--schedule', choices=('manifest', 'largest-first'), default='largest-first')
    parser.add_argument('--shard-index', type=int)
    parser.add_argument('--shard-count', type=positive_int)
    parser.add_argument('--max-novel-targets', type=int, default=96)
    parser.add_argument('--max-records', type=int, default=0)
    parser.add_argument('--samples', help='Comma-separated sample IDs; omit for all')
    parser.add_argument('--resume', action='store_true')
    from native_diagnostic_options import add_arguments, native_diagnostic_options
    add_arguments(parser)
    args = parser.parse_args()
    if not 1 <= args.jobs <= 256 or not 1 <= args.threads <= 256 or not 0 <= args.max_novel_targets <= 220 or args.max_records < 0:
        parser.error('require jobs/threads 1..256, novel target cap 0..220, max-records >=0')
    if args.jobs * args.threads > args.cpu_budget:
        parser.error(f'jobs * threads ({args.jobs * args.threads}) exceeds cpu-budget ({args.cpu_budget})')
    if (args.shard_index is None) != (args.shard_count is None) or (args.shard_count is not None and not 0 <= args.shard_index < args.shard_count):
        parser.error('shard-index and shard-count must be supplied together, with 0 <= index < count')
    if args.paired_layout and args.max_records % 2:
        parser.error('paired --max-records must be even')
    try:
        native_diagnostic_options(args)
        args.manifest = args.manifest.resolve()
        args.source = args.source.resolve()
        root_output = args.output.resolve()
        if root_output == args.source or root_output in args.source.parents or args.source in root_output.parents:
            raise ValueError('Evidence output and archived source directories must be separate and not nested')
        manifest = list(index_samples(read_tsv(args.manifest), 'manifest').values())
        reserved = {'runner.lock', 'trm-evidence', 'run_metadata.json', 'metrics.json', 'summary.json', 'progress.json'}
        for row in manifest:
            if None in row or not isinstance(row.get('path1'), str) or not row['path1'].strip():
                raise ValueError(f"manifest sample {row['sample_id']}: path1 must not be empty and fields must be valid")
            if row['sample_id'] in reserved:
                raise ValueError(f"manifest sample ID is a reserved output name: {row['sample_id']}")
            # Protect every declared source path, even when this invocation
            # selects another sample/shard. Do not require unselected files to exist.
            for key in ('path1', 'path2'):
                value = row.get(key)
                if value in (None, ''):
                    continue
                if not isinstance(value, str):
                    raise ValueError(f"manifest sample {row['sample_id']}: invalid {key}")
                declared = Path(value).expanduser()
                declared = (declared if declared.is_absolute() else args.manifest.parent / declared).resolve()
                if declared == root_output or root_output in declared.parents:
                    raise ValueError('Evidence output must be disjoint from every declared input file, including unselected rows')
        if args.samples:
            requested = set(args.samples.split(','))
            if not requested <= {row['sample_id'] for row in manifest}:
                raise ValueError('Requested samples absent from manifest')
            manifest = [row for row in manifest if row['sample_id'] in requested]
        if args.shard_count is not None:
            manifest = [row for row in manifest if sample_shard(row['sample_id'], args.shard_count) == args.shard_index]
        manifest = normalize_manifest_inputs(manifest, args.manifest)
        if args.paired_layout:
            for row in manifest:
                if not row.get('path2'):
                    raise ValueError(f"sample {row['sample_id']}: --paired-layout requires exactly two input files")
        if preflight(parser, args, selected_samples=len(manifest), timed=True,
                     cpu_affinity=available_cpus()):
            return 0
        args.output = root_output if args.shard_count is None else root_output / f'shard-{args.shard_index:05d}-of-{args.shard_count:05d}'
        if args.output.is_symlink() or args.output.resolve() != args.output:
            raise ValueError('Shard output directory must not be a symlink')
        for row in manifest:
            for key in ('path1', 'path2'):
                if row.get(key) and (Path(row[key]) == root_output or root_output in Path(row[key]).parents):
                    raise ValueError('Evidence output must be disjoint from input files')
        args.binary = resolve_binary(args.binary)
        for path in (args.manifest, args.catalog.resolve(), args.binary.resolve()):
            if path == root_output or root_output in path.parents:
                raise ValueError('Evidence output must be disjoint from manifest, catalog and supplied binary')
        if not args.binary.is_file() or not os.access(args.binary, os.X_OK):
            raise ValueError(f'Binary is not an executable file: {args.binary}')
        args.binary = args.binary.resolve()
        if not Path('/usr/bin/time').is_file() or not os.access('/usr/bin/time', os.X_OK):
            raise ValueError('GNU /usr/bin/time is required for per-scan resource capture')
        catalog = read_tsv(args.catalog)
        if any(not isinstance(row.get('motif'), str) or not re.fullmatch('[ACGTacgt]+', row['motif']) for row in catalog):
            raise ValueError('Catalog motifs must be nonempty ACGT strings')
        known = {canonical(row['motif']) for row in catalog if len(canonical(row['motif'])) <= 32}
        if len(known) > 256:
            raise ValueError('Catalog target family count exceeds native limit 256')
        for name in ('runner.lock', 'trm-evidence', 'run_metadata.json', 'progress.json', 'metrics.json', 'summary.json'):
            if (args.output / name).is_symlink():
                raise ValueError(f'Managed output must not be a symlink: {args.output / name}')
    except (OSError, ValueError, TypeError) as error:
        parser.error(str(error))
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / 'runner.lock').open('a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error('Another evidence runner is using this output directory/shard')
        try:
            binary = freeze_binary(args)
        except (OSError, ValueError) as error:
            parser.error(str(error))
        args.execution_id = uuid.uuid4().hex
        metadata = {'execution_version': EXECUTION_VERSION, 'execution_id': args.execution_id,
                    'started_utc': utc_now(), 'operation': 'targeted_full_input_evidence',
                    'binary_sha256': args._frozen_binary['sha256'], 'catalog_sha256': sha256(args.catalog),
                    'script_sha256': sha256(__file__), 'jobs': args.jobs, 'threads': args.threads,
                    'cpu_budget': args.cpu_budget, 'record_resources': True, 'schedule': args.schedule,
                    'resource_plan': args._resource_plan,
                    'checkpoint_seconds': args.checkpoint_seconds, 'samples': len(manifest),
                    'max_records': args.max_records, 'source': str(args.source), 'manifest_sha256': sha256(args.manifest),
                    'native_options': native_options(args), 'shard_index': args.shard_index, 'shard_count': args.shard_count,
                    'binary_identity_policy': 'Frozen runner-owned snapshot hash verified once at start/end; full stat identity checked per job; direct run_one computes hash.',
                    'resume_identity_policy': 'Legacy identity preserved when no nondefault native options; scheduler settings do not invalidate sample evidence.'}
        write_json(args.output / 'run_metadata.json', metadata)
        report = run_scheduled(manifest, args, known, binary)
        if file_fingerprint(binary) != args._frozen_binary['fingerprint'] or sha256(binary) != args._frozen_binary['sha256']:
            report['all_succeeded'] = False
            report['execution_error'] = 'Frozen binary changed during execution'
            write_json(args.output / 'metrics.json', report)
            summary_path = args.output / 'summary.json'
            with summary_path.open() as handle:
                summary = json.load(handle)
            summary.update(all_succeeded=False, execution_error=report['execution_error'])
            write_json(summary_path, summary)
            with (args.output / 'progress.json').open() as handle:
                progress = json.load(handle)
            progress['execution_error'] = report['execution_error']
            write_json(args.output / 'progress.json', progress)
        metadata.update(finished_utc=utc_now(), all_succeeded=report['all_succeeded'])
        write_json(args.output / 'run_metadata.json', metadata)
        return int(not report['all_succeeded'])


if __name__ == '__main__':
    raise SystemExit(main())
