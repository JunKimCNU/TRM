#!/usr/bin/env python3
"""Bounded inspection of existing FASTQs; sequence signatures do not diagnose RNA.

Extracts a prefix for exploratory mapping, explicitly not a random or complete
sample. Original inputs are never modified. RNA-seq cDNA and genomic DNA both
use A/C/G/T; library identity needs additional alignment/provenance evidence.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import json
from pathlib import Path
import time

PATTERNS = {
    'illumina_adapter_13bp': b'AGATCGGAAGAGC',
    'nextera_adapter_19bp': b'CTGTCTCTTATACACATCT',
    'ligation_GATCGATC': b'GATCGATC',
}


def terminal_run(seq, base, right=False):
    iterator = reversed(seq) if right else iter(seq)
    n = 0
    for char in iterator:
        if char != base:
            break
        n += 1
    return n


def record(handle):
    name = handle.readline()
    if not name:
        return None
    seq, plus, qual = (handle.readline().rstrip(b'\r\n') for _ in range(3))
    if not name.startswith(b'@') or not plus.startswith(b'+') or len(seq) != len(qual) or not seq:
        raise ValueError('Malformed/empty FASTQ record in inspected prefix')
    return name.rstrip(b'\r\n'), seq.upper(), qual


def inspect_file(path, output, limit):
    started = time.monotonic()
    stat = path.stat()
    counts = Counter()
    lengths = Counter()
    bases = Counter()
    signatures = Counter()
    quality_sum = 0
    q20_bases = q30_bases = 0
    first_name = last_name = None
    input_open = gzip.open if path.suffix == '.gz' else open
    extracted = output / (path.name.removesuffix('.gz') + '.prefix.fastq.gz')
    with input_open(path, 'rb') as source, gzip.open(extracted, 'wb', compresslevel=1) as target:
        for _ in range(limit):
            item = record(source)
            if item is None:
                break
            name, seq, qual = item
            if first_name is None:
                first_name = name.decode(errors='replace')
            last_name = name.decode(errors='replace')
            target.write(name + b'\n' + seq + b'\n+\n' + qual + b'\n')
            counts['records'] += 1
            counts['bases'] += len(seq)
            lengths[len(seq)] += 1
            bases.update(seq)
            signatures[hashlib.blake2b(seq, digest_size=16).digest()] += 1
            counts['records_with_N'] += b'N' in seq
            quality_sum += sum(qual) - 33 * len(qual)
            q20_bases += sum(c >= 53 for c in qual)
            q30_bases += sum(c >= 63 for c in qual)
            for base in b'ATGC':
                for side in ('left', 'right'):
                    run = terminal_run(seq, base, right=(side == 'right'))
                    for threshold in (12, 20):
                        counts[f'{side}_{chr(base)}_run_ge{threshold}'] += run >= threshold
            for label, pattern in PATTERNS.items():
                counts[label + '_records'] += pattern in seq
    n, b = counts['records'], counts['bases']
    result = {'input_path': str(path.resolve()), 'input_size_bytes': stat.st_size,
              'input_mtime_ns': stat.st_mtime_ns, 'requested_prefix_records': limit,
              'sampling': 'first_records_per_file_not_random', 'sampled_prefix': True,
              'first_record_name': first_name, 'last_record_name': last_name,
              'extracted_path': str(extracted.resolve()),
              'extracted_sha256': hashlib.sha256(extracted.read_bytes()).hexdigest(),
              'counts': dict(counts), 'length_histogram': dict(lengths),
              'gc_fraction': (bases[ord('G')] + bases[ord('C')]) / b if b else None,
              'mean_phred_assuming_phred33': quality_sum / b if b else None,
              'q20_fraction_assuming_phred33': q20_bases / b if b else None,
              'q30_fraction_assuming_phred33': q30_bases / b if b else None,
              'distinct_exact_sequences': len(signatures),
              'extra_identical_sequence_fraction': (n - len(signatures)) / n if n else None,
              'largest_identical_sequence_count': max(signatures.values(), default=0),
              'elapsed_seconds': time.monotonic() - started}
    if path.stat().st_size != stat.st_size or path.stat().st_mtime_ns != stat.st_mtime_ns:
        raise ValueError('Input changed during prefix inspection')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, default=Path('manifests/existing_fastq.tsv'))
    p.add_argument('--samples', nargs='+', required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--records-per-file', type=int, default=100000)
    p.add_argument('--jobs', type=int, default=4)
    args = p.parse_args()
    if args.records_per_file <= 0 or args.jobs <= 0:
        p.error('Record count and jobs must be positive')
    rows = {r['sample_id']: r for r in csv.DictReader(args.manifest.open(), delimiter='\t')}
    if set(args.samples) - rows.keys():
        p.error('Unknown sample IDs')
    args.output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    jobs = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for sid in args.samples:
            dest = args.output / sid
            dest.mkdir()
            for key in ('path1', 'path2'):
                if rows[sid].get(key):
                    jobs.append((sid, key, pool.submit(inspect_file, Path(rows[sid][key]), dest, args.records_per_file)))
        results = []
        for sid, mate, future in jobs:
            result = future.result()
            result.update(sample_id=sid, input_slot=mate)
            (args.output / sid / f'{mate}.json').write_text(json.dumps(result, indent=2) + '\n')
            results.append(result)
            print(f'{sid} {mate}: {result["counts"]["records"]} records', flush=True)
    (args.output / 'inspection.json').write_text(json.dumps({
        'schema_version': 1, 'elapsed_seconds': time.monotonic() - start,
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'results': results,
        'interpretation': 'Descriptive prefix QC only. Terminal polyA/T, duplicate sequences, adapters, GC or ligation sequence alone cannot distinguish RNA-seq from genomic WGS or establish WGA.'}, indent=2) + '\n')


if __name__ == '__main__':
    main()
