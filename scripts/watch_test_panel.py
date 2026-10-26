#!/usr/bin/env python3
"""Scan each fully downloaded, archive-verified pilot run and update a live report.

Never reads .part files. Downloading and scanning have separate timings. This
watcher makes no network requests and does not change the frozen detector.
"""
import argparse
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import html
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import re

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from trm_assessment import canonical


def sha(path, algorithm='sha256'):
    digest = hashlib.new(algorithm)
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def table(path):
    with path.open() as source:
        return list(csv.DictReader(source, delimiter='\t'))


def write_table(path, rows, fields=None):
    target = io.StringIO(newline='')
    writer = csv.DictWriter(target, fieldnames=fields or list(rows[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
    atomic(path, target.getvalue())


def atomic(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(value)
    os.replace(temporary, path)


def inspect_download(row, roots, cache, verified_metadata=None):
    paths, retained = [], 0
    urls, sizes, hashes = (row[k].split(';') for k in ('download_urls', 'fastq_bytes', 'fastq_md5'))
    if not len(urls) == len(sizes) == len(hashes) or not 1 <= len(urls) <= 2:
        raise ValueError('Unaligned archive file lists')
    for url, size, md5 in zip(urls, sizes, hashes):
        name, size, md5 = url.rsplit('/', 1)[1], int(size), md5.lower()
        if size <= 0 or not re.fullmatch(r'[a-f0-9]{32}', md5) or not re.fullmatch(r'[A-Za-z0-9_.-]+', name) or name in ('.', '..'):
            raise ValueError('Invalid archive size or MD5')
        valid = None
        observed = 0
        for root in roots:
            path = root / row['sample_id'] / name
            part = path.with_name(path.name + '.part')
            if part.exists():
                observed = max(observed, min(part.stat().st_size, size))
            if path.exists():
                stat = path.stat()
                observed = max(observed, min(stat.st_size, size))
                key = (str(path), stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, size, md5)
                # Recheck failed validation: filesystems may retain timestamps
                # during quick same-size writes, so a failed cache is not final.
                if not cache.get(key):
                    verified = stat.st_size == size and sha(path, 'md5') == md5
                    after = path.stat()
                    unchanged = (after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns) == (stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
                    cache[key] = verified and unchanged
                if cache[key]:
                    valid = path.resolve()
                    if verified_metadata is not None:
                        verified_metadata[str(valid)] = dict(path=str(valid), size=stat.st_size,
                            inode=stat.st_ino, mtime_ns=stat.st_mtime_ns, ctime_ns=stat.st_ctime_ns)
                    break
        retained += observed
        if valid:
            paths.append(str(valid))
    return paths if len(paths) == len(urls) else [], retained


def execute(command, log):
    started = time.monotonic()
    try:
        result = subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log.write_text(result.stdout)
        return dict(command=command, returncode=result.returncode, wall_seconds=time.monotonic()-started)
    except OSError as error:
        log.write_text(str(error) + '\n')
        return dict(command=command, returncode=127, error=str(error), wall_seconds=time.monotonic()-started)


def source_fingerprint(paths):
    result = []
    for name in paths:
        path = Path(name)
        stat = path.stat()
        result.append(dict(path=str(path.resolve()), size=stat.st_size, inode=stat.st_ino,
                           mtime_ns=stat.st_mtime_ns, ctime_ns=stat.st_ctime_ns))
    return result


def scan_artifacts(directory, paths, executable_sha):
    """Validate the batch-owned commit record, rather than a process exit alone."""
    raw_path, metrics_path, state_path = (directory / name for name in ('result.json', 'metrics.json', 'state.json'))
    raw = json.loads(raw_path.read_text())
    metrics = json.loads(metrics_path.read_text())
    state = json.loads(state_path.read_text())
    current = source_fingerprint(paths)
    if not current or state.get('status') != 'complete' or state.get('result_sha256') != sha(raw_path):
        raise ValueError('Missing/invalid batch completion fingerprint')
    expected_inputs = [{key: row[key] for key in ('path', 'size', 'mtime_ns')} for row in current]
    if state.get('fingerprint', {}).get('inputs') != expected_inputs:
        raise ValueError('Scan inputs differ from currently verified archive files')
    if state.get('fingerprint', {}).get('executable', {}).get('sha256') != executable_sha:
        raise ValueError('Scan executable differs from frozen binary')
    options = state.get('fingerprint', {}).get('options', [])
    if '--prior-policy' not in options or options[options.index('--prior-policy') + 1:][:1] != ['canonical']:
        raise ValueError('Scan does not establish canonical prior policy')
    if raw.get('summary', {}).get('sampled_prefix') is not False:
        raise ValueError('Result does not establish a complete, unsampled input scan')
    for field in ('records', 'bases'):
        if not isinstance(raw['summary'].get(field), int) or raw['summary'][field] < 0:
            raise ValueError('Invalid result summary')
    if metrics.get('exit_code') != 0 or not isinstance(metrics.get('wall_seconds'), (int, float)) or metrics['wall_seconds'] <= 0:
        raise ValueError('Invalid scanner timing/completion metrics')
    if metrics.get('input_compressed_bytes') != sum(row['size'] for row in current):
        raise ValueError('Scanner byte count differs from verified inputs')
    if state.get('metrics') != metrics:
        raise ValueError('Timing metrics differ from batch completion record')
    artifacts = {name: sha(directory / name) for name in ('result.json', 'metrics.json', 'state.json')}
    return raw, metrics, current, artifacts


def freeze(path, values, allow_additions=False):
    if path.exists():
        previous = json.loads(path.read_text())
        if previous != values and not (allow_additions and all(values.get(key) == value for key, value in previous.items())):
            raise ValueError('Frozen inputs/code changed; use a separate output directory')
    atomic(path, json.dumps(values, indent=2) + '\n')


def panel_metrics(successes):
    """Cumulative per-run metrics; never use overwritten one-run batch summary."""
    total_bytes = sum(metrics['input_compressed_bytes'] for _, metrics in successes.values())
    wall = sum(metrics['wall_seconds'] for _, metrics in successes.values())
    cpu_known = all(isinstance(metrics.get(field), (int, float)) for _, metrics in successes.values()
                    for field in ('user_cpu_seconds', 'system_cpu_seconds'))
    return dict(completed_samples=len(successes), sample_ids=sorted(successes),
        completed_input_compressed_bytes=total_bytes,
        total_records=sum(raw['summary']['records'] for raw, _ in successes.values()),
        total_bases=sum(raw['summary']['bases'] for raw, _ in successes.values()),
        sum_scan_wall_seconds=wall,
        scan_only_compressed_MBps=total_bytes / wall / 1e6 if wall else None,
        cpu_accounting_complete=cpu_known,
        sum_user_cpu_seconds=sum(metrics.get('user_cpu_seconds') or 0 for _, metrics in successes.values()) if cpu_known else None,
        sum_system_cpu_seconds=sum(metrics.get('system_cpu_seconds') or 0 for _, metrics in successes.values()) if cpu_known else None,
        max_sample_rss_kb=max((metrics.get('max_rss_kb') or 0 for _, metrics in successes.values()), default=0),
        interpretation='Sum of retained successful per-run scanner metrics; excludes downloads, archive checks, waiting and report generation. This is not end-to-end acquisition throughput. scans/summary.json describes only the latest batch invocation.')


def report_is_ready(report, signatures, output):
    if not (report and report.get('returncode') == 0 and report.get('scan_artifacts') == signatures):
        return False
    try:
        metrics = json.loads((output/'assessed/metrics.json').read_text())
        if metrics['samples'] != len(signatures) or not (output/'assessed/index.html').is_file():
            return False
        if {row['sample_id'] for row in table(output/'assessed/samples.tsv')} != set(signatures):
            return False
        paths = {'samples.tsv': output/'assessed/samples.tsv', 'candidates.tsv': output/'assessed/candidates.tsv',
                 'comparison.tsv': output/'raw_vs_assessed/comparison.tsv'}
        return all(metrics['files'][name] == sha(path) for name, path in paths.items())
    except (OSError, ValueError, KeyError, TypeError):
        return False


def changed_frozen_files(values):
    return [name for name, digest in values.items()
            if (sha(Path(name)) if Path(name).is_file() else None) != digest]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--download-roots', type=Path, nargs='+', required=True)
    p.add_argument('--metadata', type=Path, required=True)
    p.add_argument('--assignments', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--binary', type=Path, default=ROOT/'results/release_v040_20260921/trm')
    p.add_argument('--poll-seconds', type=int, default=30)
    p.add_argument('--max-wait-seconds', type=int, default=43200)
    p.add_argument('--stop-file', type=Path, help='Default: OUTPUT/STOP_REQUESTED. Finish the active command, then stop without launching another.')
    args = p.parse_args()
    if args.poll_seconds < 1 or args.max_wait_seconds < 1:
        p.error('Time limits must be positive')
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    stop_file = args.stop_file.resolve() if args.stop_file else args.output/'STOP_REQUESTED'
    with (args.output / 'watch.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return run_watch(args, stop_file)


def run_watch(args, stop_file):
    rows, meta, assignments = table(args.manifest), table(args.metadata), table(args.assignments)
    if not rows or any(len({r['sample_id'] for r in group}) != len(group) for group in (rows, meta, assignments)):
        raise ValueError('Panel tables must be nonempty and contain unique sample IDs')
    if any(not re.fullmatch(r'[A-Za-z0-9_-]+', r['sample_id']) for r in rows):
        raise ValueError('Unsafe sample identifier')
    if {r['sample_id'] for r in rows} != {r['sample_id'] for r in meta} or {r['sample_id'] for r in rows} != {r['sample_id'] for r in assignments}:
        raise ValueError('Panel, metadata and assignment sample sets must match')
    state_path = args.output/'completed_scans.json'
    completed = json.loads(state_path.read_text()) if state_path.exists() else {}
    if not isinstance(completed, dict) or set(completed) - {r['sample_id'] for r in rows}:
        raise ValueError('Invalid completed-scan state')
    input_paths = [args.manifest, args.metadata, args.assignments, args.binary,
        ROOT/'data/telomere_motifs.tsv', ROOT/'data/taxon_motif_references.tsv']
    duplicate_audit = ROOT/'results/motif_review_20260920/duplicate_audit.json'
    if duplicate_audit.exists():
        input_paths.append(duplicate_audit)
    config = {str(path.resolve()): sha(path) for path in input_paths}
    config.setdefault(str(duplicate_audit.resolve()), None)
    config_path = args.output/'frozen_inputs.json'
    freeze(config_path, config, allow_additions=not completed)
    code_paths = [Path(__file__), *(ROOT/'scripts'/name for name in
        ('batch.py', 'report_assessment.py', 'trm_assessment.py', 'assess_adequacy.py', 'evidence_batch.py'))]
    code = {str(path.resolve()): sha(path) for path in code_paths}
    if completed and not (args.output/'frozen_code.json').exists():
        raise ValueError('Legacy completed scans lack frozen code provenance; use a separate output directory')
    freeze(args.output/'frozen_code.json', code)
    metadata = {r['sample_id']: r for r in meta}
    for row in rows:
        # Preserve taxonomy while adding the reviewed protocol fields.
        for key in ('library_selection', 'library_construction_protocol', 'experiment_title', 'sample_description'):
            metadata[row['sample_id']][key] = row.get(key, '')
    metadata_fields = sorted({k for r in meta for k in r})
    cache, start = {}, time.monotonic()
    report_path = args.output/'report_execution.json'
    report = json.loads(report_path.read_text()) if report_path.exists() else None
    attempted_reports = set()
    while True:
        status, successes = [], {}
        stopped = stop_file.exists()
        timed_out = time.monotonic()-start >= args.max_wait_seconds
        changed_inputs = changed_frozen_files({**config, **code})
        fatal_error = 'Frozen files changed: ' + '; '.join(changed_inputs) if changed_inputs else ''
        for row in rows:
            sid = row['sample_id']
            stopped = stopped or stop_file.exists()
            timed_out = timed_out or time.monotonic()-start >= args.max_wait_seconds
            paths, retained, problem, verified_metadata = [], 0, '', {}
            if not stopped and not fatal_error:
                try:
                    paths, retained = inspect_download(row, args.download_roots, cache, verified_metadata)
                except (OSError, ValueError) as error:
                    problem = str(error)
            new_scan = False
            if paths and sid not in completed and not stopped and not timed_out and not fatal_error:
                verified_fingerprint = [verified_metadata[path] for path in paths]
                if source_fingerprint(paths) != verified_fingerprint:
                    problem, paths = 'Source changed after archive verification', []
            if paths and sid not in completed and not stopped and not timed_out and not fatal_error:
                manifest = args.output/(sid+'.scan_manifest.tsv')
                if len(paths) > 2:
                    raise ValueError('This paired pilot watcher expects at most two files')
                write_table(manifest, [dict(sample_id=sid, path1=paths[0], path2=paths[1] if len(paths)>1 else '')])
                command = [sys.executable, str(ROOT/'scripts/batch.py'), str(manifest), '--binary', str(args.binary.resolve()),
                    '--output-dir', str(args.output/'scans'), '--jobs', '1', '--threads', '4', '--cpu-budget', '4',
                    '--prior-policy', 'canonical', '--record-resources', '--resume']
                execution = execute(command, args.output/(sid+'.scan.log'))
                completed[sid] = execution
                new_scan = True
                batch_summary = args.output/'scans/summary.json'
                if batch_summary.exists():
                    (args.output/'scan_invocations').mkdir(exist_ok=True)
                    atomic(args.output/'scan_invocations'/(sid+'.batch_summary.json'), batch_summary.read_text())
            scan_state = ('not_checked_stopping' if stopped else 'configuration_changed' if fatal_error else
                          'scan_failed' if sid in completed and completed[sid].get('returncode') != 0 else
                          'waiting_for_verified_download')
            if sid in completed and completed[sid].get('returncode') == 0 and not stopped and not fatal_error:
                try:
                    raw, metrics, current, artifacts = scan_artifacts(args.output/'scans'/sid, paths,
                                                                    config[str(args.binary.resolve())])
                    if new_scan:
                        if current != verified_fingerprint:
                            raise ValueError('Archive source changed between verification and completed scan')
                        completed[sid].update(source_fingerprint=current, artifact_sha256=artifacts)
                    elif completed[sid].get('source_fingerprint') != current or completed[sid].get('artifact_sha256') != artifacts:
                        raise ValueError('Completed scan source/artifact fingerprint changed or is missing')
                    successes[sid] = (raw, metrics)
                    scan_state = 'scan_complete'
                except (OSError, ValueError, KeyError, TypeError, IndexError) as error:
                    scan_state, problem = 'scan_invalid', str(error)
                    if new_scan:
                        completed[sid]['validation_error'] = problem
            if new_scan:
                atomic(state_path, json.dumps(completed, indent=2)+'\n')
            entry = dict(sample_id=sid, scientific_name=row['scientific_name'],
                state=scan_state,
                expected_motifs=row.get('expected_motifs',''), expected_architecture=row['expected_architecture'],
                expected_bytes=int(row['total_fastq_bytes']), retained_bytes=retained,
                expected_family_reads='', records='', bases='', scan_wall_seconds='',
                interpretation=problem or ('Pending' if sid not in completed else ''))
            if sid in successes:
                raw, metrics = successes[sid]
                families = {canonical(m) for m in row.get('expected_motifs','').split(',') if m}
                observed = {canonical(c['motif']):c for c in raw.get('candidates',[])+raw.get('rejected_candidates',[])}
                entry.update(expected_family_reads=sum(observed.get(f,{}).get('reads',0) for f in families) if families else '',
                    records=raw['summary']['records'], bases=raw['summary']['bases'])
                entry['scan_wall_seconds'] = metrics['wall_seconds']
                entry['interpretation'] = 'Observed repeat-family support; inspect assessed report for terminal evidence. Absence not inferred.'
            status.append(entry)
        signatures = {sid: completed[sid]['artifact_sha256'] for sid in sorted(successes)}
        signature_key = json.dumps(signatures, sort_keys=True)
        changed_inputs = changed_frozen_files({**config, **code})
        if changed_inputs:
            fatal_error = 'Frozen files changed: ' + '; '.join(changed_inputs)
        stopped = stopped or stop_file.exists()
        ready = report_is_ready(report, signatures, args.output) if signatures and not fatal_error else False
        if signatures and not ready and signature_key not in attempted_reports and not stopped and not fatal_error:
                ids = set(successes)
                write_table(args.output/'completed_metadata.tsv', [metadata[s] for s in sorted(ids)], metadata_fields)
                write_table(args.output/'completed_assignments.tsv', [r for r in assignments if r['sample_id'] in ids])
                report = execute([sys.executable, str(ROOT/'scripts/report_assessment.py'), '--source', str(args.output/'scans'),
                    '--metadata', str(args.output/'completed_metadata.tsv'), '--assignments', str(args.output/'completed_assignments.tsv'),
                    '--duplicate-audit', str(duplicate_audit),
                    '--output', str(args.output/'assessed'), '--comparison-output', str(args.output/'raw_vs_assessed')], args.output/'report.log')
                report['scan_artifacts'] = signatures
                if report['returncode'] == 0 and not report_is_ready(report, signatures, args.output):
                    report.update(returncode=1, validation_error='Report command returned success without required output artifacts')
                atomic(report_path, json.dumps(report, indent=2)+'\n')
                attempted_reports.add(signature_key)
                ready = report_is_ready(report, signatures, args.output)
        aggregate = panel_metrics(successes)
        aggregate['updated_utc'] = datetime.now(timezone.utc).isoformat()
        # A stop can skip source verification; retain the last validated aggregate
        # rather than replace it with an artificial zero/partial-sample total.
        aggregate_path = args.output/'panel_metrics.json'
        if not ((stopped or fatal_error) and aggregate_path.exists()):
            atomic(aggregate_path, json.dumps(aggregate, indent=2)+'\n')
        now = datetime.now(timezone.utc).isoformat()
        all_done = all(s['state'] in {'scan_complete','scan_failed','scan_invalid'} for s in status)
        timed_out = time.monotonic()-start >= args.max_wait_seconds
        stopped = stopped or stop_file.exists()
        complete = all_done and all(s['state']=='scan_complete' for s in status) and ready and not stopped and not fatal_error
        report_state = 'ready' if ready else 'failed' if report and report.get('returncode') != 0 and report.get('scan_artifacts') == signatures else 'pending'
        state = dict(updated_utc=now, pid=os.getpid(), watcher_state='configuration_changed' if fatal_error else
            'stop_requested' if stopped else 'finished' if complete else 'finished_with_errors' if all_done else 'wait_limit_reached' if timed_out else 'watching',
            complete=complete, report_state=report_state, report_ready=ready, error=fatal_error,
            stop_request_path=str(stop_file),
            elapsed_seconds=time.monotonic()-start, samples=status,
            scope='Complete runs only. No raw read download, no AWS deployment, frozen canonical policy, no parameter tuning.')
        atomic(args.output/'status.json', json.dumps(state, indent=2)+'\n')
        write_table(args.output/'panel_status.tsv', status)
        fields = ['sample_id','scientific_name','state','expected_motifs','retained_bytes','expected_bytes','expected_family_reads','scan_wall_seconds']
        body = ''.join('<tr>'+''.join('<td>'+html.escape(str(r[k]))+'</td>' for k in fields)+'</tr>' for r in status)
        atomic(args.output/'index.html', '<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="60">'
            '<title>TRM 새 검증 panel 진행 상황</title><style>body{font:15px system-ui;margin:2em}td,th{padding:.6em;border:1px solid #ddd}table{border-collapse:collapse}</style>'
            '<h1>새 검증 panel</h1><p>Updated '+now+' · '+state['watcher_state']+'</p>'
            '<p>받는 중인 파일은 분석하지 않습니다. 두 FASTQ의 크기·MD5가 모두 일치한 run만 전체 분석합니다. '
            '기대 motif 관측은 염색체 말단 위치의 확정과 다릅니다. 미검출은 TRM 부재 판정이 아닙니다.</p>'+
            ('<p><a href="assessed/index.html">검증된 완료 샘플의 상세 해석</a> · ' if ready else
             '<p>상세 해석 보고서: '+html.escape(report_state)+' (현재 완료 목록에 대해 검증되지 않음). · ')+
            '<a href="panel_status.tsv">TSV</a> · <a href="panel_metrics.json">누적 분석 성능</a></p>'
            '<table><tr>'+''.join('<th>'+f+'</th>' for f in fields)+'</tr>'+body+'</table>')
        if all_done or timed_out or stopped or fatal_error:
            return 0 if complete or stopped else 1
        # Recheck stop requests at least once per second while idle.
        for _ in range(args.poll_seconds):
            if stop_file.exists():
                break
            time.sleep(1)


if __name__ == '__main__':
    raise SystemExit(main())
