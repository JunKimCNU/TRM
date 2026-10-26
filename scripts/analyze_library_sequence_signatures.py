#!/usr/bin/env python3
"""Reference-free descriptive signatures from saved FASTQ prefixes only.

This audit does not assign RNA/DNA identity, read independence, or TRM adequacy.
The exact canonical nematode leader panel is intentionally incomplete.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import time

LEADERS = {'SL1': b'GGTTTAATTACCCAAGTTTGAG', 'SL2': b'GGTTTTAACCCAGTTACTCAAG'}
COMPLEMENT = bytes.maketrans(b'ACGT', b'TGCA')
LIGATION = b'GATCGATC'
VALID = re.compile(b'[ACGT]+')
RNA_IDS = {'ERR13319074', 'ERR13319101', 'ERR13623718', 'ERR13623736'}
CONTROL_IDS = {'ERR1233601', 'DRR674033'}
SOURCES = [
    {'purpose': 'Canonical nematode SL1/SL2 primer sequences in a primary experiment',
     'url': 'https://link.springer.com/article/10.1186/1471-2121-6-5'},
    {'purpose': 'SL2 variation and incomplete terminal signatures in nematode RNA data',
     'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC5833154/'},
    {'purpose': 'Noncanonical nematode spliced leaders: a negative canonical panel cannot rule out RNA',
     'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC2271357/'},
    {'purpose': 'Primary Arima Hi-C methods: GATC-GATC is one possible ligation junction',
     'url': 'https://pmc.ncbi.nlm.nih.gov/articles/7195841/'},
    {'purpose': 'Hi-C sequence signatures need read mapping and interaction QC',
     'url': 'https://academic.oup.com/gigascience/article/9/1/giz158/5695848'},
]


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def positions(seq, motif):
    start = 0
    while True:
        start = seq.find(motif, start)
        if start < 0:
            return
        yield start
        start += 1


def leader_signature(seq, motif, slack=3, flank=20):
    """Full motif at read 5' or its RC at 3', including overlapping occurrences."""
    forward, reverse = list(positions(seq, motif)), list(positions(seq, revcomp(motif)))
    left = [pos for pos in forward if pos <= slack]
    right = [pos for pos in reverse if len(seq) - pos - len(motif) <= slack]
    contexts = []
    for pos in left:
        contexts.append(seq[pos + len(motif):pos + len(motif) + flank])
    for pos in right:
        contexts.append(revcomp(seq[max(0, pos - flank):pos]))
    contexts = [c.decode() for c in contexts if len(c) == flank and VALID.fullmatch(c)]
    return {'anywhere_occurrences': len(forward) + len(reverse),
            'anywhere_records': int(bool(forward or reverse)),
            'anchored_5prime_records': int(bool(left)),
            'anchored_rc_3prime_records': int(bool(right)),
            'anchored_records': int(bool(left or right)),
            'contexts': contexts}


def valid_windows(seq, k):
    return sum(max(0, len(m.group()) - k + 1) for m in VALID.finditer(seq))


def expected_probability(bases, transition_counts, transition_denominators, markov=False):
    total = sum(bases.values())
    if not total:
        return 0.0
    p = bases.get('G', 0) / total
    for a, b in zip(LIGATION, LIGATION[1:]):
        if markov:
            key = bytes((a, b)).decode()
            denom = transition_denominators.get(chr(a), 0)
            p *= transition_counts.get(key, 0) / denom if denom else 0.0
        else:
            p *= bases.get(chr(b), 0) / total
    return p


def inspect_prefix(row):
    started = time.monotonic()
    path = Path(row['extracted_path'])
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    if sha != row['extracted_sha256']:
        raise ValueError(f'Prefix checksum mismatch: {path}')
    counts, bases, transitions, denominators = Counter(), Counter(), Counter(), Counter()
    leader_counts = {name: Counter() for name in LEADERS}
    contexts = {name: Counter() for name in LEADERS}
    needed_pairs = sorted(set(bytes(pair) for pair in zip(LIGATION, LIGATION[1:])))
    assert all(pair[0] != pair[1] for pair in needed_pairs), 'bytes.count assumes non-self-overlap'
    with gzip.open(path, 'rb') as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().rstrip(b'\r\n').upper()
            plus, quality = handle.readline(), handle.readline().rstrip(b'\r\n')
            if not header.startswith(b'@') or not plus.startswith(b'+') or not seq or len(seq) != len(quality):
                raise ValueError(f'Malformed saved prefix: {path}')
            counts['records'] += 1
            counts['bases'] += len(seq)
            counts['terminal_polyA_or_rc_polyT_ge20_records'] += seq.endswith(b'A' * 20) or seq.startswith(b'T' * 20)
            counts['terminal_polyA_or_rc_polyT_ge12_records'] += seq.endswith(b'A' * 12) or seq.startswith(b'T' * 12)
            all_leader = False
            for name, motif in LEADERS.items():
                hit = leader_signature(seq, motif)
                contexts[name].update(hit.pop('contexts'))
                leader_counts[name].update(hit)
                all_leader |= bool(hit['anchored_records'])
                suffix = motif[-15:]
                partial = any(seq.startswith(suffix, offset) for offset in range(4)) or any(
                    seq[max(0, len(seq)-offset-len(suffix)):len(seq)-offset] == revcomp(suffix) for offset in range(4))
                leader_counts[name]['suffix15_anchored_records_including_partial'] += partial
            counts['any_full_leader_anchored_records'] += all_leader
            found = list(positions(seq, LIGATION))
            internal = [pos for pos in found if pos >= 8 and pos + len(LIGATION) <= len(seq) - 8]
            counts['ligation_records'] += bool(found)
            counts['ligation_occurrences'] += len(found)
            counts['ligation_internal_records'] += bool(internal)
            counts['ligation_internal_occurrences'] += len(internal)
            counts['all_valid_8mer_windows'] += valid_windows(seq, 8)
            counts['internal_valid_8mer_windows'] += valid_windows(seq[8:-8], 8)
            # Segment boundaries exclude ambiguous bases and cross-read transitions.
            for match in VALID.finditer(seq):
                segment = match.group()
                for base in b'ACGT':
                    value = segment.count(base)
                    bases[chr(base)] += value
                    denominators[chr(base)] += value - (segment[-1] == base)
                for pair in needed_pairs:
                    transitions[pair.decode()] += segment.count(pair)
    if counts['records'] != row['counts']['records'] or counts['bases'] != row['counts']['bases']:
        raise ValueError(f'Prefix record/base totals differ from original inspection: {path}')
    expected = {}
    for model in ('iid', 'markov1'):
        probability = expected_probability(bases, transitions, denominators, model == 'markov1')
        for domain in ('all', 'internal'):
            windows = counts[f'{domain}_valid_8mer_windows']
            expected[f'{domain}_{model}_expected_occurrences'] = windows * probability
    leader_rows = {}
    for name in LEADERS:
        leader_rows[name] = dict(leader_counts[name]) | {
            'distinct_20bp_contexts': len(contexts[name]),
            'context_counts': dict(contexts[name]),
            'largest_context_count': max(contexts[name].values(), default=0)}
    return {'sample_id': row['sample_id'], 'input_slot': row['input_slot'],
            'prefix_path': str(path), 'prefix_sha256': sha, 'counts': dict(counts),
            'leaders': leader_rows, 'base_counts_acgt': dict(bases),
            'transition_counts_for_model': dict(transitions),
            'transition_denominators': dict(denominators), 'ligation_expectation': expected,
            'elapsed_seconds': time.monotonic() - started,
            'inspection_qc': {k: row[k] for k in ('gc_fraction', 'extra_identical_sequence_fraction',
                                                  'largest_identical_sequence_count')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspection', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    inspection = json.loads(args.inspection.read_text())
    with args.metadata.open() as stream:
        metadata = {row['sample_id']: row for row in csv.DictReader(stream, delimiter='\t')}
    files = []
    for row in inspection['results']:
        result = inspect_prefix(row)
        files.append(result)
        print(row['sample_id'], row['input_slot'], result['counts']['any_full_leader_anchored_records'], flush=True)
    groups = defaultdict(list)
    for row in files:
        groups[row['sample_id']].append(row)
    samples = []
    for sid, group in groups.items():
        counts, expected = Counter(), Counter()
        leader = {name: Counter() for name in LEADERS}
        for row in group:
            counts.update(row['counts'])
            expected.update(row['ligation_expectation'])
            for name in LEADERS:
                leader[name].update(row['leaders'][name]['context_counts'])
        meta = metadata[sid]
        sample = {'sample_id': sid, 'scientific_name': meta['scientific_name'],
                  'role': 'RNA_PolyA_metadata_conflict' if sid in RNA_IDS else (
                      'ordinary_genomic_comparison' if sid in CONTROL_IDS else 'HiC_protocol_annotation'),
                  'library_protocol': meta['library_construction_protocol'],
                  'prefix_files': len(group), 'prefix_records': counts['records'], 'prefix_bases': counts['bases'],
                  'sampling': 'first_100000_records_per_available_file_not_random',
                  'full_SL_anchored_records': counts['any_full_leader_anchored_records'],
                  'full_SL_anchored_per_million_reads': counts['any_full_leader_anchored_records'] / counts['records'] * 1e6,
                  'polyA_or_rc_polyT_ge20_records': counts['terminal_polyA_or_rc_polyT_ge20_records'],
                  'polyA_or_rc_polyT_ge20_per_million_reads': counts['terminal_polyA_or_rc_polyT_ge20_records'] / counts['records'] * 1e6,
                  'ligation_GATCGATC_records': counts['ligation_records'],
                  'ligation_read_fraction': counts['ligation_records'] / counts['records'],
                  'ligation_internal_occurrences': counts['ligation_internal_occurrences'],
                  'ligation_internal_valid_windows': counts['internal_valid_8mer_windows'],
                  'ligation_internal_iid_expected': expected['internal_iid_expected_occurrences'],
                  'ligation_internal_markov1_expected': expected['internal_markov1_expected_occurrences'],
                  'ligation_internal_enrichment_vs_iid': counts['ligation_internal_occurrences'] / expected['internal_iid_expected_occurrences'] if expected['internal_iid_expected_occurrences'] else None,
                  'ligation_internal_enrichment_vs_markov1': counts['ligation_internal_occurrences'] / expected['internal_markov1_expected_occurrences'] if expected['internal_markov1_expected_occurrences'] else None,
                  'RNA_DNA_identity_from_this_audit': 'unresolved',
                  'interpretation': 'Descriptive signature comparison; no RNA/DNA or telomere-adequacy diagnosis. Read pairs are not independent observations.'}
        for name in LEADERS:
            for key in ('anywhere_records', 'anywhere_occurrences', 'anchored_records', 'anchored_5prime_records',
                        'anchored_rc_3prime_records', 'suffix15_anchored_records_including_partial'):
                sample[f'{name}_{key}'] = sum(row['leaders'][name][key] for row in group)
            sample[f'{name}_distinct_anchored_20bp_contexts'] = len(leader[name])
            sample[f'{name}_largest_anchored_context_count'] = max(leader[name].values(), default=0)
        samples.append(sample)
    report = {'schema_version': 1, 'analysis': 'reference_free_saved_prefix_signature_audit',
              'source_inspection': str(args.inspection.resolve()),
              'source_inspection_sha256': hashlib.sha256(args.inspection.read_bytes()).hexdigest(),
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'elapsed_seconds': time.monotonic() - started, 'full_input_rescan': False,
              'leader_motifs': {k: v.decode() for k, v in LEADERS.items()}, 'sources': SOURCES,
              'methods': {'leader_anchor_slack_bp': 3, 'leader_context_bp': 20,
                          'partial_leader_suffix_bp': 15, 'ligation_internal_margin_each_end_bp': 8,
                          'ligation_panel': 'GATCGATC only; subset of possible Arima junctions, not a complete v2 panel',
                          'background': 'Per-file pooled ACGT IID and first-order Markov expectations × valid 8-mer windows; expectations summed across files. These ignore higher-order repeats, positional bias and dependency and are not p-values.'},
              'limitations': ['The two ordinary-genomic libraries are descriptive comparisons, not a sensitivity calibration with known RNA controls.',
                              'Leader genes themselves exist in genomic DNA. Terminal hits need context diversity and preferably splice-aware reference alignment.',
                              'Canonical-leader absence does not rule out RNA: fragmentation, species variants and protocol effects can remove the signature.',
                              'PolyA/T and duplicate counts alone do not identify molecular source.',
                              'GATCGATC also occurs naturally. Enrichment supports Hi-C protocol annotations but does not establish interaction validity.',
                              'Only saved nonrandom prefixes were inspected; mates and PCR copies may be dependent.'],
              'files': files, 'samples': samples}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / 'sequence_evidence.json').write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    with (args.output_dir / 'sequence_evidence.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(samples[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(samples)


if __name__ == '__main__':
    main()
