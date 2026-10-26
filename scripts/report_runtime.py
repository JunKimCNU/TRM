"""Bounded, restartable interpretation workers; never reads sequencing inputs.

Native/evidence measurements are immutable inputs. Only assessed JSON and small
artifact sidecars are read. Every complete cache receipt binds its functional
inputs and output bytes. Resource accounting belongs to this stage alone.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from contextlib import ExitStack, contextmanager
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
import uuid

CACHE_SCHEMA = 1
PATH_ARGS = ('source', 'evidence', 'expected_evidence', 'followup_diagnostics', 'technical_boundary_evidence', 'reference_localization_evidence', 'cross_period_evidence', 'saved_structure_evidence', 'assembly_context_evidence', 'quality_path_evidence', 'long_read_context_evidence', 'metadata',
             'library_context', 'library_evidence', 'previous_assessment', 'library_overrides',
             'published_context', 'assignments', 'catalog', 'taxon_references', 'assumptions',
             'duplicate_audit', 'output', 'comparison_output')
GLOBAL_ARGS = ('metadata', 'library_context', 'library_evidence', 'library_overrides',
               'published_context', 'assignments', 'catalog', 'taxon_references',
               'assumptions', 'duplicate_audit')
RESOURCE_FIELDS = ('sample_id', 'attempt_id', 'status', 'mode', 'wall_seconds',
                   'user_cpu_seconds', 'system_cpu_seconds', 'max_rss_kib',
                   'fingerprint_seconds', 'assessment_seconds', 'render_seconds',
                   'exit_code', 'log_file')


def utc():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def token(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                    separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def identity(path):
    path = Path(path)
    try:
        stat = path.stat()
    except FileNotFoundError:
        return None
    return [stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns,
            stat.st_mode]


def binding(path):
    path = Path(path).resolve()
    before = identity(path)
    result = {'path': str(path), 'kind': 'missing' if before is None else
              'file' if path.is_file() else 'directory' if path.is_dir() else 'other'}
    if result['kind'] == 'file':
        result.update(size_bytes=before[2], sha256=digest(path))
    if before != identity(path):
        raise ValueError('input changed while fingerprinting: ' + str(path))
    return result, before


def assert_unchanged(pins):
    for path, before in pins.items():
        if isinstance(before, dict):
            if str(Path(path).resolve()) != before['resolved_path']:
                raise ValueError('input path resolution changed during report: ' + path)
            before = before['identity']
        current = identity(path)
        # Directories are presence/type guards. Output publication or unrelated
        # sibling creation must not invalidate a diagnostic root directory.
        if before is not None and (before[5] & 0o170000) == 0o040000:
            if current is None or current[0:2] != before[0:2] or current[5] != before[5]:
                raise ValueError('input directory changed during report: ' + path)
        elif current != before:
            raise ValueError('input changed during report: ' + path)


def path_pin(path):
    lexical = str(Path(os.path.abspath(path)))
    resolved = str(Path(lexical).resolve())
    stamp = identity(lexical)
    if resolved != str(Path(lexical).resolve()):
        raise ValueError('input path resolution changed while pinning: ' + lexical)
    return lexical, {'resolved_path': resolved, 'identity': stamp}


def collect(paths):
    bindings, pins = {}, {}
    for path in paths:
        lexical, pin = path_pin(path)
        value, stamp = binding(lexical)
        if (value['path'] != pin['resolved_path'] or stamp != pin['identity']):
            raise ValueError('input path changed while fingerprinting: ' + lexical)
        assert_unchanged({lexical: pin})
        bindings[value['path']], pins[lexical] = value, pin
    return bindings, pins


def prepare_globals(args):
    from report_code_dependencies import local_dependency_paths
    paths = [getattr(args, name) for name in GLOBAL_ARGS if getattr(args, name)]
    values, pins = collect(paths)
    code, code_pins = collect(local_dependency_paths(Path(__file__).resolve().parent))
    return {'bindings': values, 'pins': pins, 'code': code, 'code_pins': code_pins}


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name('.' + path.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        with temporary.open('w') as handle:
            json.dump(value, handle, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def output_locks(*directories):
    with ExitStack() as stack:
        descriptors = []
        for directory in sorted({Path(d).resolve() for d in directories}):
            directory.mkdir(parents=True, exist_ok=True)
            lock = stack.enter_context((directory / '.report.lock').open('a+'))
            try:
                fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise ValueError('report output is locked by another writer: ' + str(directory)) from error
            descriptors.append(lock.fileno())
        yield descriptors


def validate_limits(args):
    if args.jobs < 1 or args.jobs > len(os.sched_getaffinity(0)):
        raise ValueError('--jobs must be positive and not exceed available CPU affinity')
    if args.memory_limit_mib is not None and args.memory_limit_mib < 64:
        raise ValueError('--memory-limit-mib must be at least 64 (per-worker address-space cap)')
    if args.memory_budget_mib is not None:
        if args.memory_limit_mib is None or args.memory_budget_mib <= 0:
            raise ValueError('--memory-budget-mib requires a positive per-worker --memory-limit-mib')
        if args.jobs * args.memory_limit_mib > args.memory_budget_mib:
            raise ValueError('jobs times per-worker memory limit exceeds --memory-budget-mib')


def sample_dependencies(args, sample, reviewed):
    paths = [args.source / sample, args.source / sample / 'result.json',
             args.source / sample / 'state.json',
             args.source / sample / 'original_scan_metrics.json', args.source / sample / 'metrics.json']
    for name in ('evidence', 'expected_evidence', 'followup_diagnostics'):
        root = getattr(args, name)
        if root:
            paths.extend([root, root / sample, root / sample / 'status.json', root / sample / 'evidence.json'])
    if getattr(args, 'technical_boundary_evidence', None):
        from technical_boundary_evidence import source_dependencies
        paths.extend(source_dependencies(args.technical_boundary_evidence))
    if getattr(args, 'reference_localization_evidence', None):
        from reference_localization_evidence import source_dependencies
        paths.extend(source_dependencies(args.reference_localization_evidence))
    if getattr(args, 'saved_structure_evidence', None):
        from saved_repeat_structure import source_dependencies
        paths.extend(source_dependencies(args.saved_structure_evidence))
    if getattr(args, 'assembly_context_evidence', None):
        from assembly_context_evidence import source_dependencies
        paths.extend(source_dependencies(args.assembly_context_evidence))
    if getattr(args, 'quality_path_evidence', None):
        from quality_path_evidence import source_dependencies
        paths.extend(source_dependencies(args.quality_path_evidence))
    if getattr(args, 'long_read_context_evidence', None):
        from long_read_context_evidence import source_dependencies
        paths.extend(source_dependencies(args.long_read_context_evidence))
    if getattr(args, 'cross_period_evidence', None):
        from cross_period_evidence import source_dependencies
        paths.extend(source_dependencies(args.cross_period_evidence))
    if reviewed:
        paths.append(Path(reviewed['detailed_evidence_path']))
        if reviewed.get('previous_assessment_path'):
            previous = Path(reviewed['previous_assessment_path'])
            if not previous.is_absolute():
                previous = args.library_evidence.parent / previous
            paths.append(previous)
    return collect(paths)


def presentation_dependencies(args, sample):
    paths = []
    if args.previous_assessment:
        paths = [args.previous_assessment / 'samples' / (sample + '.json'),
                 args.previous_assessment / (sample + '.json')]
    return collect(paths)


def verify_artifacts(receipt, output=None, sample=None, pins=None):
    artifacts = receipt.get('artifacts')
    verified_inodes = {}
    if not isinstance(artifacts, dict) or not artifacts:
        return False
    for path, expected in artifacts.items():
        try:
            stamp = identity(path)
            saved = verified_inodes.get(tuple(stamp)) if stamp is not None else None
            if saved is None:
                actual, stamp = binding(path)
                if stamp is not None:
                    verified_inodes[tuple(stamp)] = actual
            else:
                actual = {**saved, 'path': str(Path(path).resolve())}
                if stamp != identity(path):
                    return False
            if pins is not None and '.report-cache' in Path(path).parts:
                pins[str(path)] = stamp
        except (OSError, ValueError):
            return False
        if actual != expected:
            return False
    try:
        document, bundle_path = Path(receipt['document']), Path(receipt['bundle'])
        if document != bundle_path.parent / 'assessment.json':
            return False
        bundle = json.loads(bundle_path.read_text())
        allowed = {'input_applicability.tsv', 'flank_period_applicability.tsv',
                   'clade_variant_observations.tsv',
                   'candidate_quality_observations.tsv',
                   'technical_boundary_diagnostics.tsv', 'technical_boundary_events.tsv',
                   'reference_localization_diagnostics.tsv', 'reference_localization_families.tsv',
                   'reference_localization_references.tsv', 'reference_localization_events.tsv',
                   'cross_period_diagnostics.tsv', 'cross_period_families.tsv',
                   'cross_period_intervals.tsv', 'cross_period_events.tsv',
                   'cross_period_contexts.tsv', 'cross_period_positions.tsv',
                   'assembly_context_summary.tsv', 'assembly_context_scaffolds.tsv',
                   'assembly_context_intervals.tsv', 'assembly_context_events.tsv',
                   'quality_path_summary.tsv', 'quality_path_models.tsv', 'quality_path_read_comparisons.tsv',
                   'long_read_context_summary.tsv', 'long_read_context_reads.tsv', 'long_read_context_contexts.tsv',
                   'expected_family_observations.tsv', 'nearby_families.tsv',
                   'followup_diagnostics.tsv', 'runtime.tsv', 'adequacy.tsv',
                   'illustrative_power_scenarios.tsv', 'candidates.tsv'}
        chunks = bundle['chunks']
        localization_chunks = {'reference_localization_diagnostics.tsv', 'reference_localization_families.tsv',
                               'reference_localization_references.tsv', 'reference_localization_events.tsv'}
        cross_period_chunks = {'cross_period_diagnostics.tsv', 'cross_period_families.tsv',
                               'cross_period_intervals.tsv', 'cross_period_events.tsv',
                               'cross_period_contexts.tsv', 'cross_period_positions.tsv'}
        assembly_chunks = {'assembly_context_summary.tsv', 'assembly_context_scaffolds.tsv',
                           'assembly_context_intervals.tsv', 'assembly_context_events.tsv'}
        quality_path_chunks = {'quality_path_summary.tsv', 'quality_path_models.tsv', 'quality_path_read_comparisons.tsv'}
        long_read_chunks = {'long_read_context_summary.tsv', 'long_read_context_reads.tsv', 'long_read_context_contexts.tsv'}
        if (bundle['sample_id'] != sample or not isinstance(chunks, list)
                or len(chunks) != len(set(chunks)) or not set(chunks) <= allowed
                or not {'candidates.tsv', *localization_chunks, *cross_period_chunks, *assembly_chunks, *quality_path_chunks, *long_read_chunks} <= set(chunks)):
            return False
        required = {str(document), str(bundle_path), str(bundle_path.parent / 'assessment.html'),
                    str(bundle_path.parent / 'summary_projection.json')}
        required.update(str(bundle_path.parent / name) for name in chunks)
        if output is not None:
            required.update(str((Path(output) / 'samples' / (sample + suffix)).resolve()) for suffix in ('.json', '.html'))
            required.add(str((Path(output) / 'summaries' / (sample + '.json')).resolve()))
        return set(artifacts) == required
    except (OSError, ValueError, KeyError, TypeError):
        return False


def publish_file(source, destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name('.' + destination.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        try:
            os.link(source, temporary)
        except OSError:
            shutil.copyfile(source, temporary)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def worker(task_path, report):
    """Run in a fresh process so peak RSS and address-space limits are per sample."""
    task = json.loads(Path(task_path).read_text())
    args = argparse.Namespace(**task['args'])
    for name in PATH_ARGS:
        setattr(args, name, Path(getattr(args, name)) if getattr(args, name) else None)
    if args.memory_limit_mib:
        cap = args.memory_limit_mib * 1024 * 1024
        _, hard = resource.getrlimit(resource.RLIMIT_AS)
        if hard != resource.RLIM_INFINITY:
            cap = min(cap, hard)
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
    sample = task['sample_id']
    attempt = Path(task['attempt'])
    cache_path = Path(task['cache_path'])
    globals_ = task['globals']
    assert_unchanged(globals_['pins'])
    assert_unchanged(globals_['code_pins'])
    start = time.monotonic()
    reviewed = args.library_evidence_rows.get(sample)
    dependencies, pins = sample_dependencies(args, sample, reviewed)
    presentation, presentation_pins = presentation_dependencies(args, sample)
    if reviewed:
        expected = reviewed['evidence_sha256']
        actual = dependencies[str(Path(reviewed['detailed_evidence_path']).resolve())]
        if actual.get('sha256') != expected:
            raise ValueError('reviewed library evidence changed after validation: ' + sample)
        if reviewed.get('previous_assessment_path'):
            previous = Path(reviewed['previous_assessment_path'])
            if not previous.is_absolute():
                previous = args.library_evidence.parent / previous
            if dependencies[str(previous.resolve())].get('sha256') != reviewed.get('previous_assessment_sha256'):
                raise ValueError('reviewed library previous assessment changed after validation: ' + sample)
    reviewed_key = {k: v for k, v in (reviewed or {}).items() if k != 'library_evidence_table_sha256'}
    functional = token({'schema': CACHE_SCHEMA, 'code': globals_['code'],
        'python_version': sys.version, 'python_executable': str(Path(sys.executable).resolve()),
        'metadata': task['metadata'], 'assignment': task['assignment'],
        'references': task['references'], 'taxon_references': args.taxon_reference_rows,
        'assumptions': task['assumptions'], 'duplicate_groups': task['duplicate_groups'],
        'reviewed_library_evidence': reviewed_key, 'dependencies': dependencies})
    provenance_paths = {str(getattr(args, name).resolve()) for name in ('catalog', 'assignments', 'taxon_references', 'library_evidence') if getattr(args, name)}
    presentation_key = token({'globals': {p: b for p, b in globals_['bindings'].items() if p in provenance_paths}, 'previous_assessment': presentation,
        'output': str(args.output), 'comparison_output': str(args.comparison_output),
        'json_format': args.json_format})
    prior = None
    cache_read_pins = {}
    if args.resume and cache_path.is_file():
        try:
            possible = json.loads(cache_path.read_text())
            required = {'schema_version', 'sample_id', 'status', 'functional_key', 'presentation_key',
                        'dependencies', 'dynamic_dependencies', 'dynamic_aliases', 'document', 'bundle', 'artifacts'}
            if (required <= possible.keys() and possible.get('schema_version') == CACHE_SCHEMA and possible.get('sample_id') == sample
                    and possible.get('status') == 'complete' and possible.get('functional_key') == functional
                    and verify_artifacts(possible, args.output, sample, cache_read_pins)):
                alias_pins = {}
                for lexical, resolved in possible['dynamic_aliases'].items():
                    original, pin = path_pin(lexical)
                    if pin['resolved_path'] != resolved:
                        raise ValueError('cached diagnostic alias changed: ' + original)
                    alias_pins[original] = pin
                dynamic_paths = possible['dynamic_dependencies']
                if any(Path(p).is_file() and Path(p).stat().st_size > 64 * 1024 * 1024 for p in dynamic_paths):
                    raise ValueError('oversized cached diagnostic artifact')
                dynamic, dynamic_pins = collect(dynamic_paths)
                if dynamic == possible.get('dynamic_dependencies', {}):
                    prior = possible
                    pins.update(dynamic_pins)
                    pins.update(alias_pins)
        except (OSError, ValueError, TypeError, KeyError):
            prior = None
    if prior is None:
        cache_read_pins = {}
    assert_unchanged(cache_read_pins)
    fingerprint_seconds = time.monotonic() - start
    mode = 'recomputed'
    assessment_seconds = render_seconds = 0.0
    if prior and prior['presentation_key'] == presentation_key:
        mode = 'reused'
        result = dict(prior)
    else:
        # Reuse already computed SHA256 values inside sample_report; this avoids
        # a second full native JSON hash pass. File identities are checked again
        # after reading/assessment before anything is published.
        all_bindings = {**globals_['bindings'], **dependencies, **presentation}
        original_sha = report.sha256
        def pinned_sha(path):
            value = all_bindings.get(str(Path(path).resolve()))
            return value['sha256'] if value and value['kind'] == 'file' else original_sha(path)
        report.sha256 = pinned_sha
        dynamic_dependencies, dynamic_pins, dynamic_aliases = {}, {}, {}
        if prior:
            mode = 'refreshed'
            document = json.loads(Path(prior['document']).read_text())
            bundle = json.loads(Path(prior['bundle']).read_text())
            row, comparison = bundle['summary'], bundle['comparison']
            for argument, field in (('catalog', 'catalogue_sha256'), ('assignments', 'assignment_sha256'),
                                    ('taxon_references', 'taxon_reference_sha256')):
                value = getattr(args, argument)
                document['provenance'][field] = pinned_sha(value) if value else None
            document['reviewed_library_evidence'] = reviewed
            document['provenance']['library_evidence_table_sha256'] = reviewed.get('library_evidence_table_sha256') if reviewed else None
            dynamic_dependencies = prior['dynamic_dependencies']
            dynamic_aliases = prior['dynamic_aliases']
        else:
            import followup_diagnostics
            original_read = followup_diagnostics._read
            original_artifact = followup_diagnostics._artifact
            def tracked_read(path):
                before = identity(Path(path).resolve())
                content, actual = original_read(path)
                resolved = actual['path']
                if before != identity(resolved):
                    raise ValueError('diagnostic input changed while reading: ' + resolved)
                dynamic_dependencies[resolved] = {'kind': 'file', **actual}
                dynamic_pins[resolved] = before
                return content, actual
            def tracked_artifact(ref, table=False):
                # Capture only artifacts actually read. Historical raw-input
                # bindings normalized by _relation are deliberately not followed.
                lexical, pin = path_pin(ref['path'])
                value, actual = original_artifact(ref, table=table)
                assert_unchanged({lexical: pin})
                if actual['path'] != pin['resolved_path']:
                    raise ValueError('diagnostic path changed while reading: ' + lexical)
                dynamic_aliases[lexical] = pin['resolved_path']
                dynamic_pins[lexical] = pin
                return value, actual
            followup_diagnostics._read = tracked_read
            followup_diagnostics._artifact = tracked_artifact
            compute_start = time.monotonic()
            try:
                document, row, comparison = report.sample_report(sample, task['metadata'], task['assignment'],
                    args, task['references'], task['assumptions'], task['duplicate_groups'])
            finally:
                followup_diagnostics._read = original_read
                followup_diagnostics._artifact = original_artifact
            assessment_seconds = time.monotonic() - compute_start
        render_start = time.monotonic()
        bundle = report.sample_artifacts(document, row, comparison, args, attempt)
        render_seconds = time.monotonic() - render_start
        report.sha256 = original_sha
        pins.update(dynamic_pins)
        assert_unchanged(cache_read_pins)
        assert_unchanged(pins)
        assert_unchanged(presentation_pins)
        assert_unchanged(globals_['pins'])
        assert_unchanged(globals_['code_pins'])
        document_path = attempt / 'assessment.json'
        html_path = attempt / 'assessment.html'
        artifact_paths = [document_path, html_path, attempt / 'bundle.json', attempt / 'summary_projection.json']
        artifact_paths.extend(attempt / name for name in bundle['chunks'])
        artifacts, _ = collect(artifact_paths)
        result = {'schema_version': CACHE_SCHEMA, 'sample_id': sample, 'status': 'complete',
                  'functional_key': functional, 'presentation_key': presentation_key,
                  'dependencies': dependencies, 'dynamic_dependencies': dynamic_dependencies,
                  'dynamic_aliases': dynamic_aliases,
                  'document': str(document_path), 'bundle': str(attempt / 'bundle.json'),
                  'artifacts': artifacts, 'completed_utc': utc()}
    assert_unchanged(cache_read_pins)
    assert_unchanged(pins)
    assert_unchanged(presentation_pins)
    assert_unchanged(globals_['pins'])
    assert_unchanged(globals_['code_pins'])
    atomic_json(attempt / 'result.json', {'receipt': result, 'mode': mode,
        'fingerprint_seconds': fingerprint_seconds, 'assessment_seconds': assessment_seconds,
        'render_seconds': render_seconds, 'input_pins': pins, 'presentation_pins': presentation_pins,
        'cache_read_pins': cache_read_pins})


def run_samples(args, metadata, assignments, references, defaults, per_sample, duplicates, globals_):
    validate_limits(args)
    assert_unchanged(globals_['pins'])
    started = time.monotonic()
    execution_id = args._execution_id
    cache_root = args.output / '.report-cache'
    executions = args.output / '.report-executions' / execution_id
    executions.mkdir(parents=True, exist_ok=True)
    status_path = args.output / 'interpretation_status.json'
    status = {'schema_version': CACHE_SCHEMA, 'status': 'running', 'execution_id': execution_id,
              'started_utc': utc(), 'controller_pid': os.getpid(), 'jobs': args.jobs,
              'memory_limit_mib_per_worker': args.memory_limit_mib,
              'memory_budget_mib_workers': args.memory_budget_mib, 'samples': len(metadata),
              'sample_ids': sorted(metadata), 'completed': 0, 'failed': 0, 'resource_scope': 'interpretation workers only; coordinator memory additional'}
    atomic_json(status_path, status)
    task_args = {k: str(v.resolve()) if isinstance(v, Path) else v for k, v in vars(args).items()
                 if k not in {'runtime_globals'}}
    rows, results = [], {}
    script = Path(__file__).resolve().with_name('report_assessment.py')
    def launch(sample):
        attempt_id = execution_id + '-' + sample
        attempt = cache_root / sample / 'attempts' / attempt_id
        attempt.mkdir(parents=True)
        # Cohort-wide sample maps must not be copied into every task. Their
        # source files remain pinned globally, and this sample's exact values
        # participate in its functional cache identity.
        sample_args = {**task_args, 'library_evidence_rows': {
            sample: args.library_evidence_rows[sample]} if sample in args.library_evidence_rows else {}}
        task = {'sample_id': sample, 'attempt': str(attempt.resolve()),
                'cache_path': str((cache_root / sample / 'receipt.json').resolve()),
                'args': sample_args, 'globals': globals_, 'metadata': metadata[sample],
                'assignment': assignments[sample], 'references': references,
                'assumptions': {**defaults, **per_sample.get(sample, {})},
                'duplicate_groups': {sample: duplicates[sample]} if sample in duplicates else {}}
        atomic_json(attempt / 'task.json', task)
        log = attempt / 'worker.log'
        clock_file = attempt / 'gnu_time.tsv'
        command = ['/usr/bin/time', '-f', '%e\t%U\t%S\t%M\t%x', '-o', str(clock_file),
                   sys.executable, '-B', str(script), '--_worker-task', str(attempt / 'task.json')]
        begin = time.monotonic()
        with log.open('w') as handle:
            # Inherit the controller lock descriptions through GNU time. If
            # the controller is killed, surviving workers keep the namespace
            # locked until they exit, preventing overlapping restart budgets.
            process = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT,
                                     pass_fds=tuple(args._runtime_lock_fds))
        measured = {}
        try:
            values = clock_file.read_text().splitlines()[-1].split('\t')
            if len(values) == 5:
                measured = dict(zip(('wall_seconds', 'user_cpu_seconds', 'system_cpu_seconds', 'max_rss_kib', 'exit_code'),
                                    [float(v) for v in values]))
        except (OSError, ValueError, IndexError):
            pass
        row = {'sample_id': sample, 'attempt_id': attempt_id,
               'status': 'complete' if process.returncode == 0 else 'failed', 'mode': None,
               'wall_seconds': time.monotonic() - begin, 'exit_code': process.returncode,
               'log_file': str(log.resolve()), **measured}
        result = None
        if process.returncode == 0:
            try:
                result = json.loads((attempt / 'result.json').read_text())
                row.update({k: result[k] for k in ('mode', 'fingerprint_seconds', 'assessment_seconds', 'render_seconds')})
                assert_unchanged(result['cache_read_pins'])
                assert_unchanged(result['input_pins'])
                assert_unchanged(result['presentation_pins'])
                assert_unchanged(globals_['pins'])
                assert_unchanged(globals_['code_pins'])
                receipt = result['receipt']
                # Only the controller publishes while holding both output
                # locks. Orphan workers can finish private attempt files only.
                if row['mode'] != 'reused':
                    for name, folder, suffix in (('assessment.json', 'samples', '.json'),
                                                 ('assessment.html', 'samples', '.html'),
                                                 ('summary_projection.json', 'summaries', '.json')):
                        source = Path(receipt['bundle']).parent / name
                        destination = (args.output / folder / (sample + suffix)).resolve()
                        publish_file(source, destination)
                        receipt['artifacts'][str(destination)] = {
                            **receipt['artifacts'][str(source)], 'path': str(destination)}
                    atomic_json(cache_root / sample / 'receipt.json', receipt)
                receipt['_runtime_artifact_pins'] = {p: identity(p) for p in receipt['artifacts']}
                receipt['_runtime_input_pins'] = result['input_pins']
                receipt['_runtime_presentation_pins'] = result['presentation_pins']
            except (OSError, ValueError, KeyError) as error:
                with log.open('a') as handle:
                    handle.write('\nController refused worker publication: ' + str(error) + '\n')
                row.update(status='failed', exit_code=1)
        atomic_json(attempt / 'resources.json', row)
        return sample, row, result
    # At most jobs futures are live; no growing queue of sample documents/tasks.
    def input_size(sample):
        path = args.source / sample / 'result.json'
        return path.stat().st_size if path.is_file() else 0
    iterator = iter(sorted(metadata, key=lambda sample: (-input_size(sample), sample)))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = {}
        for _ in range(min(args.jobs, len(metadata))):
            sample = next(iterator)
            pending[pool.submit(launch, sample)] = sample
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                del pending[future]
                sample, row, result = future.result()
                rows.append(row)
                if row['status'] == 'complete':
                    results[sample] = result
                    status['completed'] += 1
                else:
                    status['failed'] += 1
                    print(Path(row['log_file']).read_text()[-12000:], file=sys.stderr)
                atomic_json(status_path, status)
                sample = next(iterator, None)
                if sample is not None:
                    pending[pool.submit(launch, sample)] = sample
    rows.sort(key=lambda row: row['sample_id'])
    resource_file = executions / 'interpretation_resources.tsv'
    with resource_file.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=RESOURCE_FIELDS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    publish_file(resource_file, args.output / 'interpretation_resources.tsv')
    status.update(status='failed' if status['failed'] else 'workers_complete',
                  worker_stage_wall_seconds=time.monotonic() - started,
                  resources_file=str(resource_file.resolve()), finished_workers_utc=utc(),
                  modes={mode: sum(row['mode'] == mode for row in rows) for mode in ('recomputed', 'refreshed', 'reused')})
    atomic_json(status_path, status)
    atomic_json(executions / 'status.json', status)
    if status['failed']:
        raise ValueError(f"{status['failed']} interpretation worker(s) failed; completed samples remain resumable; see {resource_file}")
    assert_unchanged(globals_['pins'])
    assert_unchanged(globals_['code_pins'])
    for result in results.values():
        assert_unchanged(result['input_pins'])
        assert_unchanged(result['presentation_pins'])
    receipts = [results[sample]['receipt'] for sample in sorted(results)]
    validate_receipts(receipts)
    return receipts, status, executions


def validate_receipts(receipts):
    for receipt in receipts:
        for key in ('_runtime_artifact_pins', '_runtime_input_pins', '_runtime_presentation_pins'):
            assert_unchanged(receipt[key])


def read_bundle(receipt):
    path = Path(receipt['bundle'])
    before = identity(path)
    content = path.read_bytes()
    if (hashlib.sha256(content).hexdigest() != receipt['artifacts'][str(path)]['sha256']
            or before != identity(path)):
        raise ValueError('cached bundle changed during aggregation: ' + str(path))
    return json.loads(content)


def merge_chunks(receipts, filename, output):
    """Stream complete TSV chunks in sample order; verify bytes while copying."""
    output = Path(output)
    temporary = output.with_name('.' + output.name + '.' + uuid.uuid4().hex + '.tmp')
    header = None
    written = False
    try:
        with temporary.open('wb') as target:
            for receipt in receipts:
                chunk = Path(receipt['bundle']).parent / filename
                expected = receipt['artifacts'].get(str(chunk))
                if expected is None:
                    continue
                hasher = hashlib.sha256()
                with chunk.open('rb') as source:
                    current = source.readline()
                    hasher.update(current)
                    if header is None:
                        header = current
                        target.write(header)
                    elif current != header:
                        raise ValueError('inconsistent chunk table header: ' + str(chunk))
                    for block in iter(lambda: source.read(1024 * 1024), b''):
                        target.write(block)
                        hasher.update(block)
                if hasher.hexdigest() != expected['sha256']:
                    raise ValueError('cached table chunk changed during merge: ' + str(chunk))
                written = True
            target.flush()
            os.fsync(target.fileno())
        if written:
            os.replace(temporary, output)
        else:
            output.unlink(missing_ok=True)
    finally:
        temporary.unlink(missing_ok=True)
