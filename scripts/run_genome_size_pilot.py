#!/usr/bin/env python3
"""Separate, reproducible KMC/GenomeScope pilot; never changes TRM adequacy.

Uses complete local FASTQs, not repeat-selected reads or prefixes. Each fit is
conditional on its explicitly supplied ploidy and high-frequency cutoff.
"""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def measured_run(command, stem, env):
    start = time.monotonic()
    stem = Path(stem)
    resource_file = stem.with_suffix('.resources.json')
    timed = ['/usr/bin/time', '-f', '{"wall_seconds":%e,"user_cpu_seconds":%U,"system_cpu_seconds":%S,"max_rss_kb":%M}', '-o', str(resource_file)] + command
    with stem.with_suffix('.stdout.log').open('w') as out, stem.with_suffix('.stderr.log').open('w') as err:
        result = subprocess.run(timed, stdout=out, stderr=err, env=env)
    metrics = {'command': command, 'returncode': result.returncode,
               'elapsed_seconds': time.monotonic() - start}
    # GNU time writes an extra diagnostic line for unsuccessful children.
    lines = resource_file.read_text().splitlines()
    metrics.update(json.loads(lines[-1]))
    stem.with_suffix('.execution.json').write_text(json.dumps(metrics, indent=2) + '\n')
    if result.returncode:
        raise RuntimeError(f'Failed command; see {stem}.stderr.log')
    return metrics


def summarize_histogram(path, cap, thresholds):
    distinct = occurrences = saturated = max_observed = 0
    tail = {str(t): 0 for t in thresholds}
    with path.open() as handle:
        for line in handle:
            multiplicity, count = map(int, line.split())
            distinct += count
            occurrences += multiplicity * count
            if count:
                max_observed = multiplicity
            if multiplicity == cap:
                saturated = count
            for threshold in thresholds:
                if multiplicity > threshold:
                    tail[str(threshold)] += multiplicity * count
    return {'distinct_kmers': distinct, 'counted_kmer_occurrences_lower_bound': occurrences,
            'max_observed_multiplicity': max_observed, 'counter_cap': cap,
            'kmers_at_counter_cap': saturated, 'saturated': bool(saturated),
            'occurrences_above_fit_cutoff': tail,
            'occurrence_fraction_above_fit_cutoff': {t: n / occurrences if occurrences else None for t, n in tail.items()}}


def model_converged(output_directory):
    """Require the model's success marker and positive finite size estimates.

    GenomeScope can write summary.txt even when fitting fails. A successful
    fit is still only a model result, not evidence that its assumptions hold.
    """
    output_directory = Path(output_directory)
    try:
        progress = (output_directory / 'progress.txt').read_text().split()
        model = (output_directory / 'model.txt').read_text()
        summary = (output_directory / 'summary.txt').read_text().splitlines()
    except OSError:
        return False
    if not progress or progress[-1] not in {'done', 'pass'} or not model.strip():
        return False
    size_lines = [line for line in summary if line.startswith('Genome Haploid Length')]
    if len(size_lines) != 1:
        return False
    numbers = re.findall(r'[-+]?(?:\d[\d,]*\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', size_lines[0])
    if len(numbers) != 2:
        return False
    estimates = [float(value.replace(',', '')) for value in numbers]
    return all(math.isfinite(value) and value > 0 for value in estimates)


def run_one(row, k, args, env):
    sid = row['sample_id']
    out = args.output / sid / f'k{k}'
    out.mkdir(parents=True, exist_ok=False)
    files = [Path(row[key]) for key in ('path1', 'path2') if row.get(key)]
    input_list = out / 'inputs.txt'
    input_list.write_text(''.join(str(p.resolve()) + '\n' for p in files))
    inputs = [{'path': str(p.resolve()), 'size_bytes': p.stat().st_size,
               'mtime_ns': p.stat().st_mtime_ns} for p in files]
    database = out / 'counts'
    histogram = out / 'counts.histo'
    start = time.monotonic()
    with tempfile.TemporaryDirectory(prefix=f'trm_genome_{sid}_k{k}_') as tmp:
        count = measured_run([str(args.tool_bin / 'kmc'), f'-k{k}', f'-t{args.threads}',
                              f'-m{args.memory_gb}', '-sm', '-ci1', f'-cs{args.counter_cap}',
                              '-cx1000000000', '-fq', '-hp', f'-j{out / "kmc_summary.json"}',
                              '@' + str(input_list), str(database), tmp], out / 'count', env)
    histogram_run = measured_run([str(args.tool_bin / 'kmc_tools'), f'-t{args.threads}',
                                  'transform', str(database), 'histogram', str(histogram),
                                  '-ci1', f'-cx{args.counter_cap}'], out / 'histogram', env)
    hist = summarize_histogram(histogram, args.counter_cap, args.cutoffs)
    # Retain the exact full histogram. Smaller dense files have an explicit
    # zero sentinel because this installed GenomeScope drops its final row.
    fits = []
    for cutoff in args.cutoffs:
        fit_hist = out / f'fit_input_max{cutoff}.histo'
        with histogram.open() as src, fit_hist.open('w') as dst:
            for line in src:
                if int(line.split()[0]) > cutoff:
                    break
                dst.write(line)
            dst.write(f'{cutoff + 1}\t0\n')
        for ploidy in args.ploidies:
            fit_out = out / f'p{ploidy}_max{cutoff}'
            fit_out.mkdir()
            fit = measured_run([str(args.tool_bin / 'Rscript'), str(args.tool_bin / 'genomescope2'),
                                '-i', str(fit_hist), '-o', str(fit_out), '-k', str(k),
                                '-p', str(ploidy), '-m', str(cutoff)], fit_out / 'fit', env)
            fits.append({'ploidy_assumption': ploidy, 'max_kmer_coverage': cutoff,
                         'output_directory': str(fit_out), 'execution': fit,
                         'model_converged': model_converged(fit_out)})
    result = {'sample_id': sid, 'k': k, 'inputs': inputs, 'canonical': True,
              'all_local_reads': True, 'quality_trimming': False,
              'adapter_trimming': False, 'input_hash_note': 'FASTQ identity uses path/size/mtime and prior archive provenance; no additional full-file checksum pass.',
              'count_execution': count, 'histogram_execution': histogram_run,
              'histogram': hist, 'histogram_sha256': sha256(histogram), 'fits': fits,
              'elapsed_seconds': time.monotonic() - start,
              'interpretation': 'Exploratory genome-size estimates only. Ploidy, tail contamination, read duplication, coverage bias and telomere representation are not established.'}
    (out / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(f'Completed {sid} k={k}: {result["elapsed_seconds"]:.1f}s', flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=Path('manifests/existing_fastq.tsv'))
    parser.add_argument('--samples', nargs='+', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--tool-bin', type=Path, required=True)
    parser.add_argument('--k', nargs='+', type=int, default=[21, 31])
    parser.add_argument('--ploidies', nargs='+', type=int, default=[1, 2])
    parser.add_argument('--cutoffs', nargs='+', type=int, default=[10000])
    parser.add_argument('--counter-cap', type=int, default=1000000)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--memory-gb', type=int, default=8)
    parser.add_argument('--jobs', type=int, default=3)
    args = parser.parse_args()
    if not (args.jobs > 0 and args.threads > 0 and args.memory_gb > 0 and
            all(1 <= k <= 256 for k in args.k) and
            all(1 <= p <= 6 for p in args.ploidies) and
            all(20 <= c < args.counter_cap <= 1000000000 for c in args.cutoffs)):
        parser.error('Invalid resource/model limits or cutoff >= counter cap')
    args.output = args.output.resolve()
    args.tool_bin = args.tool_bin.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    rows = {r['sample_id']: r for r in csv.DictReader(args.manifest.open(), delimiter='\t')}
    unknown = set(args.samples) - rows.keys()
    if unknown:
        parser.error(f'Unknown samples: {sorted(unknown)}')
    env = dict(os.environ)
    env['PATH'] = str(args.tool_bin) + os.pathsep + env.get('PATH', '')
    env['OPENBLAS_NUM_THREADS'] = '1'
    env['OMP_NUM_THREADS'] = '1'
    start = time.monotonic()
    tools = {name: {'path': str(args.tool_bin / name), 'sha256': sha256(args.tool_bin / name)}
             for name in ('kmc', 'kmc_tools', 'Rscript', 'genomescope2')}
    (args.output / 'run_configuration.json').write_text(json.dumps({
        'arguments': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        'tools': tools, 'script_sha256': sha256(__file__), 'started_unix': time.time()}, indent=2) + '\n')
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_one, rows[sid], k, args, env) for k in args.k for sid in args.samples]
        results = [future.result() for future in futures]
    (args.output / 'results.json').write_text(json.dumps({'elapsed_seconds': time.monotonic() - start,
                                                       'results': results}, indent=2) + '\n')


if __name__ == '__main__':
    main()
