#!/usr/bin/env python3
"""Offline synthetic paired FASTQ -> scan -> resume -> assessment smoke check.

The output directory must not exist. No real specimen, taxonomic identification,
telomere sensitivity, chromosome-end confirmation or network access is implied.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
import random
import subprocess
import sys

from runtime_binary import ROOT, resolve_binary, sha256


def table(path, rows):
    with path.open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--binary', help='Explicit native binary override')
    args = parser.parse_args()
    binary = resolve_binary(args.binary)
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    rng = random.Random(260801)
    random_seq = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    sequences = ['CCCTAA' * 25 for _ in range(16)]
    sequences += ['CCCTAA' * 15 + random_seq(60) for _ in range(8)]
    sequences += [random_seq(150) for _ in range(8)]
    for mate in (1, 2):
        records = []
        for i, sequence in enumerate(sequences):
            if mate == 2:
                sequence = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
            records.append(f'@synthetic_{i}/{mate}\n{sequence}\n+\n' + 'I' * len(sequence) + '\n')
        (out / f'R{mate}.fastq.gz').write_bytes(gzip.compress(''.join(records).encode(), mtime=0))
    table(out / 'samples.tsv', [{'sample_id': 'synthetic', 'path1': 'R1.fastq.gz', 'path2': 'R2.fastq.gz'}])
    table(out / 'metadata.tsv', [{'sample_id': 'synthetic', 'status': 'unresolved',
                                 'scientific_name': 'Synthetic software fixture', 'tax_id': '', 'lineage': '',
                                 'library_strategy': '', 'library_source': '', 'library_layout': 'PAIRED'}])
    commands = []

    def run(label, command):
        with (out / f'{label}.log').open('w') as log:
            result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, timeout=120)
        commands.append({'label': label, 'argv': command, 'exit_code': result.returncode})
        (out / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
        if result.returncode:
            raise RuntimeError(f'{label} failed; see {out / (label + ".log")}')

    py = [sys.executable, '-B']
    batch = py + [str(ROOT / 'scripts/batch.py'), str(out / 'samples.tsv'), '--binary', str(binary),
                  '--output-dir', str(out / 'native'), '--jobs', '1', '--threads', '1', '--cpu-budget', '1',
                  '--thread-budget', '4', '--paired-layout', 'fr', '--memory-limit-mib', '512',
                  '--quality-diagnostics', 'phred33']
    run('scan', batch)
    result_path = out / 'native/synthetic/result.json'
    result = json.loads(result_path.read_text())
    require(result['summary']['records'] == 64, 'Unexpected input record count')
    require(result['summary']['bases'] == 9600, 'Unexpected input base count')
    observed = result['candidates'] + result.get('rejected_candidates', [])
    family = next(c for c in observed if c['motif'] == 'AACCCT')
    require(family['reads'] >= 48, 'Synthetic TTAGGG-family support missing')
    before = sha256(result_path)
    run('resume', batch + ['--resume'])
    summary = json.loads((out / 'native/summary.json').read_text())
    require(summary['resumed_samples'] == 1 and sha256(result_path) == before, 'Resume did not preserve the completed result')
    run('priors', py + [str(ROOT / 'scripts/prepare_taxon_priors.py'), '--metadata', str(out / 'metadata.tsv'),
                       '--references', str(ROOT / 'data/taxon_motif_references.tsv'), '--output-dir', str(out / 'priors')])
    run('assessment', py + [str(ROOT / 'scripts/report_assessment.py'), '--source', str(out / 'native'),
                           '--metadata', str(out / 'metadata.tsv'), '--assignments', str(out / 'priors/assignments.tsv'),
                           '--catalog', str(ROOT / 'data/telomere_motifs.tsv'),
                           '--taxon-references', str(ROOT / 'data/taxon_motif_references.tsv'),
                           '--output', str(out / 'assessment'), '--comparison-output', str(out / 'comparison'),
                           '--json-format', 'compact', '--jobs', '1'])
    assessment = json.loads((out / 'assessment/samples/synthetic.json').read_text())
    require(assessment['schema_version'] == 22, 'Unexpected assessment schema')
    require((out / 'assessment/index.html').is_file(), 'Assessment HTML missing')
    receipt = {'status': 'passed', 'synthetic_only': True, 'records': 64, 'bases': 9600,
               'expected_family': 'AACCCT', 'observed_family_reads': family['reads'],
               'result_sha256': before, 'resume_preserved_result': True,
               'assessment_schema': 22, 'binary_sha256': sha256(binary)}
    (out / 'smoke.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
