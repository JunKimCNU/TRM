#!/usr/bin/env python3
"""Download complete, archive-verified pilot files with observable timings.

Downloads with bounded concurrency, archive integrity checks, and resumable
partial files. Optional single-first observations are not controlled network trials.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time
from urllib.parse import urlsplit
import uuid


VERIFIED = {'verified_existing', 'downloaded_verified', 'verified_existing_partial'}
BYTE_ACCOUNTING = ('transferred_bytes counts net retained file growth during this invocation, '
                   'not total network traffic; curl size_download may omit internal retries')


def save(path, data):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    os.replace(temporary, path)


def digests(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    size = 0
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            md5.update(block); sha.update(block); size += len(block)
    return size, md5.hexdigest(), sha.hexdigest()


def _fetch(item, out, phase, invocation_id):
    dest = out / item['sample_id']
    dest.mkdir(exist_ok=True)
    path = dest / item['name']
    started = time.monotonic()
    if path.exists():
        size, md5, sha = digests(path)
        if (size, md5) != (item['expected_bytes'], item['expected_md5']):
            return dict(item, path=str(path.resolve()), status='existing_archive_mismatch_preserved',
                        observed_bytes=size, observed_md5=md5, sha256=sha,
                        elapsed_seconds=time.monotonic() - started, phase=phase,
                        transferred_bytes=0)
        result = dict(item, path=str(path.resolve()), status='verified_existing', sha256=sha,
                      elapsed_seconds=time.monotonic() - started, phase=phase,
                      transferred_bytes=0)
        save(dest / (item['name'] + '.download.json'), result)
        return result
    part = dest / (item['name'] + '.part')
    resumed_from = part.stat().st_size if part.exists() else 0
    # A complete partial needs no HTTP request (servers may reject an EOF range).
    # Oversized/corrupt files are preserved for inspection, never silently reset.
    if resumed_from >= item['expected_bytes']:
        check_started = time.monotonic()
        size, md5, sha = digests(part)
        valid = (size, md5) == (item['expected_bytes'], item['expected_md5'])
        if valid:
            os.replace(part, path)
        return dict(item, path=str(path.resolve()), status=('verified_existing_partial' if valid
                    else 'archive_mismatch_partial_preserved'), observed_bytes=size,
                    observed_md5=md5, sha256=sha, resumed_from_bytes=resumed_from,
                    transferred_bytes=0, phase=phase,
                    verification_seconds=time.monotonic() - check_started,
                    elapsed_seconds=time.monotonic() - started)
    command = ['curl', '-L', '--fail', '--retry', '2', '--retry-delay', '2',
               '--connect-timeout', '30', '--max-time', '3600', '--speed-time', '120',
               '--speed-limit', '1024', '--continue-at', '-', '--output', str(part),
               '--write-out', '{"http_code":%{http_code},"size_download":%{size_download},"speed_download":%{speed_download},"time_total":%{time_total},"time_starttransfer":%{time_starttransfer}}',
               '--silent', '--show-error', item['url']]
    attempts = []
    attempt_path = dest / (item['name'] + '.attempts.json')
    if attempt_path.exists():
        attempts = json.loads(attempt_path.read_text())
        if not isinstance(attempts, list) or not all(isinstance(row, dict) for row in attempts):
            raise ValueError('Invalid attempt-history format; original history preserved')
    history_count = len(attempts)
    for attempt in range(4):
        before = part.stat().st_size if part.exists() else 0
        attempt_start = time.monotonic()
        started_utc = datetime.now(timezone.utc).isoformat()
        save(dest / (item['name'] + '.active.json'), dict(item, pid=os.getpid(),
            status='running', invocation_id=invocation_id,
            started_utc=started_utc, resumed_from_bytes=before, attempt=len(attempts) + 1))
        execution = subprocess.run(command, capture_output=True, text=True)
        with (dest / (item['name'] + '.curl.stderr.log')).open('a') as log:
            log.write(execution.stderr)
        try:
            transfer = json.loads(execution.stdout)
        except json.JSONDecodeError:
            transfer = {'raw_stdout': execution.stdout}
        after = part.stat().st_size if part.exists() else 0
        attempts.append(dict(started_utc=started_utc, returncode=execution.returncode,
            invocation_id=invocation_id, attempt_this_invocation=attempt + 1,
            before_bytes=before, after_bytes=after, retained_new_bytes=after-before,
            elapsed_seconds=time.monotonic()-attempt_start, curl=transfer))
        save(attempt_path, attempts)
        # Curl versions differ in which interruptions --retry considers transient.
        # Restart using the actual retained offset for premature EOF/reset/timeout.
        # A nominal curl success can still leave an incomplete archive (for
        # example an empty proxy response). Resume it within the same bounded
        # budget; only the expected archive size plus MD5 establishes success.
        if after >= item['expected_bytes'] or execution.returncode not in {0, 18, 28, 52, 56}:
            break
        if attempt < 3:
            time.sleep(2)
    result = dict(item, path=str(path.resolve()), command=command, returncode=execution.returncode,
                  phase=phase, resumed_from_bytes=resumed_from, curl=transfer, attempts=attempts,
                  prior_attempt_count=history_count, attempts_this_invocation=len(attempts)-history_count,
                  retained_new_bytes_this_invocation=after-resumed_from,
                  transferred_bytes=after-resumed_from,
                  transfer_wall_seconds=time.monotonic() - started)
    observed_size = part.stat().st_size if part.exists() else 0
    if not part.exists() or observed_size < item['expected_bytes']:
        # Incomplete transport is retryable regardless of curl's exit code.
        # Do not hash/rename an incomplete archive or discard its retained bytes.
        result.update(status='download_failed_partial_preserved',
                      failure_reason='archive_incomplete', observed_bytes=observed_size,
                      missing_bytes=item['expected_bytes'] - observed_size,
                      elapsed_seconds=time.monotonic() - started)
    else:
        check_started = time.monotonic()
        size, md5, sha = digests(part)
        result.update(observed_bytes=size, observed_md5=md5, sha256=sha,
                      verification_seconds=time.monotonic() - check_started)
        if (size, md5) != (item['expected_bytes'], item['expected_md5']):
            result['status'] = 'archive_mismatch_partial_preserved'
        else:
            os.replace(part, path)
            result.update(status='downloaded_verified', transferred_bytes=size - resumed_from)
        result['elapsed_seconds'] = time.monotonic() - started
    save(dest / (item['name'] + '.download.json'), result)
    print(f'{item["sample_id"]} {item["name"]}: {result["status"]} {result["elapsed_seconds"]:.1f}s', flush=True)
    return result


def fetch(item, out, phase):
    """One-file failures return evidence; other completed files remain reviewable."""
    dest = out / item['sample_id']
    dest.mkdir(parents=True, exist_ok=True)
    lock_path = dest / (item['name'] + '.lock')
    started = time.monotonic()
    invocation_id = uuid.uuid4().hex
    with lock_path.open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            # Do not overwrite the active owner process's per-file records.
            return dict(item, status='download_in_progress', phase=phase, transferred_bytes=0,
                        elapsed_seconds=time.monotonic() - started)
        try:
            identity = {key: item[key] for key in ('sample_id', 'name', 'url', 'expected_bytes', 'expected_md5')}
            identity_path = dest / (item['name'] + '.identity.json')
            for previous_path in (identity_path, dest / (item['name'] + '.active.json'),
                                  dest / (item['name'] + '.download.json')):
                if previous_path.exists():
                    previous = json.loads(previous_path.read_text())
                    if not isinstance(previous, dict):
                        raise ValueError('Invalid prior file provenance; file preserved')
                    if any(key in previous and previous[key] != value for key, value in identity.items()):
                        raise ValueError('Prior archive identity differs; existing bytes and identity preserved')
            if not identity_path.exists():
                save(identity_path, identity)
            result = _fetch(item, out, phase, invocation_id)
        except Exception as error:
            result = dict(item, path=str((dest / item['name']).resolve()),
                          status='download_exception_partial_preserved', phase=phase,
                          error_type=type(error).__name__, error=str(error),
                          elapsed_seconds=time.monotonic() - started,
                          transferred_bytes=None)
        result.update(invocation_id=invocation_id, byte_accounting=BYTE_ACCOUNTING)
        save(dest / (item['name'] + '.download.json'), result)
        save(dest / (item['name'] + '.active.json'), dict(item, pid=os.getpid(),
             status='finished', result_status=result['status'], invocation_id=invocation_id,
             finished_utc=datetime.now(timezone.utc).isoformat()))
        return result


def planned_items(rows):
    items = []
    for row in rows:
        if not re.fullmatch(r'[A-Za-z0-9_-]+', row['sample_id']):
            raise ValueError('Unsafe sample identifier')
        urls, sizes, md5s = (row[k].split(';') for k in ('download_urls', 'fastq_bytes', 'fastq_md5'))
        if not len(urls) == len(sizes) == len(md5s) or not urls:
            raise ValueError('Archive file lists do not align')
        if len(urls) > 2:
            raise ValueError('The scan manifest supports at most two FASTQ files per sample')
        for url, size, md5 in zip(urls, sizes, md5s):
            parsed = urlsplit(url)
            name = parsed.path.rsplit('/', 1)[-1]
            if (parsed.scheme != 'https' or parsed.netloc != 'ftp.sra.ebi.ac.uk'
                    or parsed.query or parsed.fragment
                    or not re.fullmatch(r'[A-Za-z0-9_.-]+\.fastq\.gz', name)
                    or not re.fullmatch(r'[a-fA-F0-9]{32}', md5) or int(size) <= 0):
                raise ValueError('Unexpected archive URL, size or MD5')
            items.append({'sample_id': row['sample_id'], 'name': name,
                          'url': url, 'expected_bytes': int(size), 'expected_md5': md5.lower()})
    if not items or len({(x['sample_id'], x['name']) for x in items}) != len(items):
        raise ValueError('Empty/duplicate file plan')
    return items


def safe_fetch(item, out, phase):
    try:
        return fetch(item, out, phase)
    except Exception as error:
        # Even an unwritable per-file destination must not lose other results.
        return dict(item, status='download_exception_partial_preserved', phase=phase,
                    error_type=type(error).__name__, error=str(error), transferred_bytes=None)


def supervise(args):
    """Bounded restart rounds after a prior invocation releases the output."""
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/'supervisor.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        started = time.monotonic()
        identity = uuid.uuid4().hex
        history = args.output/'rounds'/identity
        history.mkdir(parents=True)
        state = dict(invocation_id=identity, started_utc=datetime.now(timezone.utc).isoformat(),
                     rounds_planned=args.rounds, rounds_finished=0, status='waiting_for_previous_invocation')
        save(args.output/'supervisor.json', state)
        if args.wait_for_finished:
            while True:
                progress = args.output/'progress.json'
                if progress.exists() and json.loads(progress.read_text()).get('status') == 'finished':
                    break
                if time.monotonic()-started > 86400:
                    state['status'] = 'previous_invocation_wait_expired'
                    save(args.output/'supervisor.json', state)
                    return 1
                time.sleep(30)
        child_args = ['--manifest', str(args.manifest), '--output', str(args.output), '--jobs', str(args.jobs)]
        if args.samples:
            child_args += ['--samples', *args.samples]
        for number in range(1, args.rounds + 1):
            state.update(status='downloading', round=number)
            save(args.output/'supervisor.json', state)
            code = main(child_args)
            summary = json.loads((args.output/'download_summary.json').read_text())
            save(history/f'round_{number}.json', summary)
            state.update(rounds_finished=number, elapsed_seconds=time.monotonic()-started)
            if code == 0:
                state['status'] = 'complete'
                break
            # Retry transport failures only; integrity/identity errors need review.
            if any(r['status'] not in VERIFIED | {'download_failed_partial_preserved'} for r in summary['results']):
                state['status'] = 'stopped_for_integrity_or_local_error'
                break
            if number == args.rounds or time.monotonic()-started > 86400:
                state['status'] = 'retry_budget_exhausted'
                break
            state['status'] = 'retry_delay'
            save(args.output/'supervisor.json', state)
            time.sleep(30)
        state['finished_utc'] = datetime.now(timezone.utc).isoformat()
        save(args.output/'supervisor.json', state)
        return 0 if state['status']=='complete' else 1


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--samples', nargs='+')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--jobs', type=int, default=3)
    p.add_argument('--single-first', action='store_true', help='Finish one file before concurrent transfers')
    p.add_argument('--rounds', type=int, default=1, help='Bounded resumable transfer rounds (1..12)')
    p.add_argument('--wait-for-finished', action='store_true', help='Wait for an already running invocation of this output before resumable rounds')
    args = p.parse_args(argv)
    if args.jobs < 1 or not 1 <= args.rounds <= 12:
        p.error('jobs must be positive and rounds 1..12')
    with args.manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if args.samples:
        if set(args.samples) - {r['sample_id'] for r in rows}:
            p.error('Unknown samples')
        rows = [r for r in rows if r['sample_id'] in args.samples]
    items = planned_items(rows)
    if args.rounds > 1 or args.wait_for_finished:
        return supervise(args)
    args.output.mkdir(parents=True, exist_ok=True)
    source_manifest = args.manifest.read_bytes()
    (args.output / 'source_manifest.tsv').write_bytes(source_manifest)
    provenance = {'source_manifest_sha256': hashlib.sha256(source_manifest).hexdigest(),
                  'runner_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    start = time.monotonic()
    started_utc = datetime.now(timezone.utc).isoformat()
    save(args.output / 'progress.json', {'files_planned': len(items), 'files_finished': 0,
        'started_utc': started_utc, 'pid': os.getpid(), 'results': []})
    results = [safe_fetch(items[0], args.output, 'single_transfer')] if args.single_first else []
    parallel_start = time.monotonic()
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(safe_fetch, item, args.output, 'concurrent_batch'): item for item in (items[1:] if args.single_first else items)}
        for future in as_completed(futures):
            try:
                results.append(future.result())
            except Exception as error:
                results.append(dict(futures[future], status='download_exception_partial_preserved',
                                    error_type=type(error).__name__, error=str(error), transferred_bytes=None))
            save(args.output / 'progress.json', {'files_planned': len(items), 'files_finished': len(results),
                                               'started_utc': started_utc, 'pid': os.getpid(),
                                               'elapsed_seconds': time.monotonic() - start, 'results': results})
    verified = {(r['sample_id'], r['name']): r for r in results if r['status'] in VERIFIED}
    manifest_path = args.output / 'verified_fastq.tsv'
    with manifest_path.open('w') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=['sample_id', 'path1', 'path2'])
        writer.writeheader()
        for row in rows:
            names = [url.rsplit('/', 1)[1] for url in row['download_urls'].split(';')]
            if all((row['sample_id'], name) in verified for name in names) and len(names) <= 2:
                paths = [verified[(row['sample_id'], name)]['path'] for name in names]
                writer.writerow({'sample_id': row['sample_id'], 'path1': paths[0],
                                 'path2': paths[1] if len(paths) > 1 else ''})
    save(args.output / 'download_summary.json', {'complete': len(verified) == len(items),
        **provenance,
        'started_utc': started_utc, 'finished_utc': datetime.now(timezone.utc).isoformat(),
        'elapsed_seconds': time.monotonic() - start, 'concurrent_phase_seconds': time.monotonic() - parallel_start,
        'concurrent_jobs': args.jobs, 'planned_bytes': sum(i['expected_bytes'] for i in items),
        'verified_files': len(verified), 'planned_files': len(items), 'results': results,
        'retained_file_growth_bytes': sum(r.get('transferred_bytes') or 0 for r in results),
        'retained_file_growth_accounting_complete': all(r.get('transferred_bytes') is not None for r in results),
        'byte_accounting': BYTE_ACCOUNTING,
        'interpretation': 'Single-transfer and concurrent observations differ in files and time; do not infer a controlled parallel speedup. Existing-file checks transfer zero bytes.'})
    save(args.output / 'progress.json', {'files_planned': len(items), 'files_finished': len(results),
         'status': 'finished', 'complete': len(verified) == len(items), 'started_utc': started_utc,
         'pid': os.getpid(), 'elapsed_seconds': time.monotonic() - start, 'results': results})
    return 0 if len(verified) == len(items) else 1


if __name__ == '__main__':
    raise SystemExit(main())
