#!/usr/bin/env python3
"""Bounded library-identity evidence against a matched annotated reference.

This diagnostic never modifies TRM classifications. A splice-negative result is
inconclusive; synthetic controls are technical checks, not biological controls.
"""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import random
import re
import subprocess
import time


CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')
CANONICAL = {'GTAG': '+', 'GCAG': '+', 'ATAC': '+',
             'CTAC': '-', 'CTGC': '-', 'GTAT': '-'}
CATEGORIES = ('intergenic', 'gene_nonexonic', 'exonic')


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def read_fasta(path):
    records = {}
    name = None
    pieces = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    records[name] = ''.join(pieces).upper()
                name, pieces = line[1:].split()[0], []
            else:
                pieces.append(line.strip())
    if name is not None:
        records[name] = ''.join(pieces).upper()
    return records


def read_annotation(path, genome):
    """Use 0-based half-open coordinates; exon union overrides gene body."""
    masks = {name: bytearray(len(seq)) for name, seq in genome.items()}
    genes, exons, transcripts = [], [], defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            chrom, feature = fields[0], fields[2]
            if feature not in ('gene', 'exon'):
                continue
            start, end = int(fields[3]) - 1, int(fields[4])
            if chrom not in genome or not 0 <= start < end <= len(genome[chrom]):
                raise ValueError('Annotation coordinates outside the matched FASTA')
            if feature == 'gene':
                genes.append((chrom, start, end))
            else:
                exons.append((chrom, start, end))
                attrs = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
                for parent in attrs.get('Parent', '').split(','):
                    if parent:
                        transcripts[parent].append((chrom, start, end, fields[6]))
    for chrom, start, end in genes:
        masks[chrom][start:end] = b'\1' * (end - start)
    for chrom, start, end in exons:
        masks[chrom][start:end] = b'\2' * (end - start)
    introns = set()
    for parent, segments in transcripts.items():
        segments.sort(key=lambda s: (s[0], s[1], s[2]))
        for left, right in zip(segments, segments[1:]):
            if left[0] == right[0] and left[3] == right[3] and left[2] < right[1]:
                introns.add((left[0], left[2], right[1], left[3]))
    whole, acgt = Counter(), Counter()
    for chrom, seq in genome.items():
        whole.update(masks[chrom])
        acgt.update(mask for base, mask in zip(seq, masks[chrom]) if base in 'ACGT')
    return masks, transcripts, introns, {
        'all_reference_bases': {name: whole[i] for i, name in enumerate(CATEGORIES)},
        'reference_ACGT_bases': {name: acgt[i] for i, name in enumerate(CATEGORIES)},
        'definition': 'Exon union; annotated gene body minus exon union; remaining reference. Gene-nonexonic includes introns and any unannotated exonic/UTR bases. Exons take precedence across overlapping genes.',
        'annotated_intron_coordinates_stranded': len(introns),
    }


def cigar_events(cigar, start):
    """Return matched blocks and skipped intervals with SAM query positions."""
    ops = [(int(length), code) for length, code in CIGAR.findall(cigar)]
    if ''.join(str(length) + code for length, code in ops) != cigar:
        raise ValueError('Invalid CIGAR: ' + cigar)
    ref, query, blocks, junctions = start, 0, [], []
    for index, (length, code) in enumerate(ops):
        if code in 'M=X':
            blocks.append((ref, ref + length, query, query + length))
            ref += length
            query += length
        elif code == 'N':
            left = ops[index - 1][0] if index and ops[index - 1][1] in 'M=X' else 0
            right = ops[index + 1][0] if index + 1 < len(ops) and ops[index + 1][1] in 'M=X' else 0
            junctions.append((ref, ref + length, query, left, right))
            ref += length
        elif code == 'D':
            ref += length
        elif code in 'IS':
            query += length
        # H and P consume neither query nor reference.
    return blocks, junctions


def junction_details(chrom, seq, read, junction, introns, anchor=10):
    start, end, query, left, right = junction
    motif = seq[start:start + 2] + seq[end - 2:end]
    strand = CANONICAL.get(motif)
    long_enough = end - start >= 20 and left >= anchor and right >= anchor
    exact = (long_enough and query >= anchor and query + anchor <= len(read)
             and all(base in 'ACGT' for base in read[query - anchor:query + anchor])
             and read[query - anchor:query] == seq[start - anchor:start]
             and read[query:query + anchor] == seq[end:end + anchor])
    annotated = any((chrom, start, end, s) in introns for s in ('+', '-'))
    annotated_same_strand = bool(strand and (chrom, start, end, strand) in introns)
    return motif, strand, annotated, annotated_same_strand, bool(exact and strand)


def fraction(numerator, denominator):
    return numerator / denominator if denominator else None


def analyze_sam(path, genome, masks, introns, output, mapq=30, anchor=10):
    counts, base_all, base_acgt, motifs = Counter(), Counter(), Counter(), Counter()
    pairs_high, pairs_strong, starts_strong = set(), set(), set()
    junction_rows = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith('@'):
                continue
            fields = line.rstrip('\n').split('\t')
            flag = int(fields[1])
            counts['all_sam_alignment_records'] += 1
            if flag & (0x100 | 0x800):
                counts['secondary_or_supplementary_records'] += 1
                continue
            counts['primary_records'] += 1
            if flag & 4:
                counts['primary_unmapped_records'] += 1
                continue
            counts['primary_mapped_records'] += 1
            if flag & 2:
                counts['primary_proper_pair_records'] += 1
            blocks, skips = cigar_events(fields[5], int(fields[3]) - 1)
            if skips:
                counts['primary_spliced_records_any_mapq'] += 1
            quality = int(fields[4])
            if quality == 255 or quality < mapq:
                continue
            counts['high_mapq_primary_mapped_records'] += 1
            pairs_high.add(fields[0])
            if flag & 0x400:
                counts['high_mapq_marked_duplicate_records'] += 1
            chrom, read = fields[2], fields[9].upper()
            for start, end, _, _ in blocks:
                base_all.update(masks[chrom][start:end])
                base_acgt.update(mask for base, mask in zip(genome[chrom][start:end], masks[chrom][start:end]) if base in 'ACGT')
            if skips:
                counts['high_mapq_primary_spliced_records'] += 1
            strong_read = False
            annotated_strong_read = False
            for skip in skips:
                start, end, _, left, right = skip
                motif, strand, annotated, same_strand, strong = junction_details(
                    chrom, genome[chrom], read, skip, introns, anchor)
                counts['high_mapq_junction_observations'] += 1
                motifs[motif] += 1
                key = (chrom, start, end)
                row = junction_rows.setdefault(key, {
                    'chrom': chrom, 'intron_start_0based': start, 'intron_end_exclusive': end,
                    'genomic_dinucleotides': motif, 'inferred_intron_strand': strand,
                    'annotation_coordinate_match': annotated,
                    'annotation_coordinate_and_strand_match': same_strand,
                    'high_mapq_read_observations': 0, 'strong_read_observations': 0,
                    'pairs': set(), 'strong_pairs': set(), 'strong_placements': set(),
                })
                row['high_mapq_read_observations'] += 1
                row['pairs'].add(fields[0])
                if strong:
                    strong_read = True
                    annotated_strong_read |= same_strand
                    counts['strong_junction_observations'] += 1
                    counts['strong_annotation_matched_junction_observations'] += int(same_strand)
                    row['strong_read_observations'] += 1
                    row['strong_pairs'].add(fields[0])
                    placement = (int(fields[3]), bool(flag & 16), fields[5], bool(flag & 64))
                    row['strong_placements'].add(placement)
                    starts_strong.add((chrom, *placement))
            if strong_read:
                counts['strong_spliced_primary_records'] += 1
                pairs_strong.add(fields[0])
            if annotated_strong_read:
                counts['strong_annotation_matched_spliced_primary_records'] += 1
    rows = []
    for row in junction_rows.values():
        for field in ('pairs', 'strong_pairs', 'strong_placements'):
            row['distinct_' + field] = len(row.pop(field))
        rows.append(row)
    rows.sort(key=lambda row: (-row['strong_read_observations'], row['chrom'], row['intron_start_0based']))
    fields = list(rows[0]) if rows else ['chrom', 'intron_start_0based', 'intron_end_exclusive']
    with open(output / 'junctions.tsv', 'w') as handle:
        handle.write('\t'.join(fields) + '\n')
        for row in rows:
            handle.write('\t'.join(str(row[field]) for field in fields) + '\n')
    hq = counts['high_mapq_primary_mapped_records']
    return {
        'counts': dict(counts),
        'thresholds': {'mapq_min': mapq, 'mapq_255_excluded': True, 'min_exact_anchor_bases_each_side': anchor, 'min_intron_length': 20},
        'primary_mapped_fraction': fraction(counts['primary_mapped_records'], counts['primary_records']),
        'high_mapq_spliced_fraction': fraction(counts['high_mapq_primary_spliced_records'], hq),
        'high_mapq_strong_spliced_fraction': fraction(counts['strong_spliced_primary_records'], hq),
        'high_mapq_strong_annotation_matched_spliced_fraction': fraction(counts['strong_annotation_matched_spliced_primary_records'], hq),
        'distinct_high_mapq_pairs': len(pairs_high), 'distinct_strong_spliced_pairs': len(pairs_strong),
        'distinct_strong_read_placements': len(starts_strong),
        'distinct_high_mapq_junctions': len(rows),
        'distinct_strong_junctions': sum(row['strong_read_observations'] > 0 for row in rows),
        'distinct_strong_annotation_matched_junctions': sum(row['strong_read_observations'] > 0 and row['annotation_coordinate_and_strand_match'] for row in rows),
        'distinct_strong_junctions_with_multiple_placements': sum(row['distinct_strong_placements'] >= 2 for row in rows),
        'high_mapq_junction_motif_observations': dict(motifs),
        'high_mapq_aligned_reference_bases': {name: base_all[i] for i, name in enumerate(CATEGORIES)},
        'high_mapq_aligned_reference_ACGT_bases': {name: base_acgt[i] for i, name in enumerate(CATEGORIES)},
        'high_mapq_aligned_reference_ACGT_fractions': {name: fraction(base_acgt[i], sum(base_acgt.values())) for i, name in enumerate(CATEGORIES)},
        'counting_notes': 'Read-level primary SAM records; supplementary/secondary excluded. Paired overlapping bases counted twice. Junction requires a CIGAR N, canonical GT-AG/GC-AG/AT-AC or reverse-complement genomic motif, >=20 bp gap, and >=10 exact aligned bases immediately on each side. Distinct placements are not proof of independent molecules; no duplicate removal.',
    }


def reverse_complement(seq):
    return seq.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def simulate_controls(genome, transcripts, output, count=4000, seed=20260922, read_length=151, fragment_length=350):
    """Ideal error-free controls; RNA deliberately enriched for junction reads."""
    rng = random.Random(seed)
    from bisect import bisect_right
    chroms = [chrom for chrom, seq in genome.items() if len(seq) >= fragment_length]
    cumulative, total = [], 0
    for chrom in chroms:
        total += len(genome[chrom]) - fragment_length + 1
        cumulative.append(total)
    rna = []
    for parent, segments in transcripts.items():
        if len(segments) < 2 or len({s[0] for s in segments}) != 1 or len({s[3] for s in segments}) != 1:
            continue
        seq = ''.join(genome[c][start:end] for c, start, end, _ in segments)
        if len(seq) < fragment_length:
            continue
        positions, length = [], 0
        for _, start, end, _ in segments[:-1]:
            length += end - start
            positions.append(length)
        if segments[0][3] == '-':
            seq = reverse_complement(seq)
            positions = sorted(len(seq) - pos for pos in positions)
        rna.append((seq, positions))
    details = {'pairs_per_control': count, 'read_length': read_length, 'fragment_length': fragment_length,
               'random_seed': seed, 'sequencing_error_rate': 0,
               'limitations': 'Same-reference ideal synthetic reads, no polymorphism, sequencing bias or contamination. RNA fragments selected to have a >=10 bp anchored junction in at least one mate; not an expression model or estimate of biological sensitivity. Genomic fragments uniform over possible starts conditional on ACGT-only fragment.'}
    if not rna or not chroms:
        raise ValueError('Reference does not permit both synthetic controls')
    for label in ('synthetic_genomic_DNA', 'synthetic_junction_enriched_RNA'):
        directory = output / label
        directory.mkdir(exist_ok=True)
        attempts, accepted = 0, 0
        with gzip.open(directory / 'R1.fastq.gz', 'wt') as f1, gzip.open(directory / 'R2.fastq.gz', 'wt') as f2:
            while accepted < count:
                attempts += 1
                if attempts > count * 1000:
                    raise RuntimeError('Cannot obtain requested synthetic fragments')
                if label == 'synthetic_genomic_DNA':
                    draw = rng.randrange(total)
                    index = bisect_right(cumulative, draw)
                    start = draw - (cumulative[index - 1] if index else 0)
                    frag = genome[chroms[index]][start:start + fragment_length]
                else:
                    seq, boundaries = rng.choice(rna)
                    start = rng.randrange(len(seq) - fragment_length + 1)
                    relative = [pos - start for pos in boundaries]
                    if not any(10 <= pos <= read_length - 10 or fragment_length - read_length + 10 <= pos <= fragment_length - 10 for pos in relative):
                        continue
                    frag = seq[start:start + fragment_length]
                if any(base not in 'ACGT' for base in frag):
                    continue
                accepted += 1
                if rng.randrange(2):
                    frag = reverse_complement(frag)
                for handle, read in ((f1, frag[:read_length]), (f2, reverse_complement(frag[-read_length:]))):
                    handle.write(f'@{label}_{accepted}\n{read}\n+\n' + 'I' * read_length + '\n')
        details[label] = {'accepted_pairs': accepted, 'attempted_fragments': attempts}
    return details


def measured_run(command, output, label):
    start = time.perf_counter()
    with open(output / (label + '.stdout.log'), 'w') as stdout, open(output / (label + '.stderr.log'), 'w') as stderr:
        result = subprocess.run(['/usr/bin/time', '-v', '-o', str(output / (label + '.resources.txt')), *command], stdout=stdout, stderr=stderr)
    row = {'argv': command, 'returncode': result.returncode, 'elapsed_seconds': time.perf_counter() - start}
    (output / (label + '.command.json')).write_text(json.dumps(row, indent=2) + '\n')
    if result.returncode:
        raise RuntimeError('External command failed; see ' + str(output / (label + '.stderr.log')))
    return row


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genome', type=Path, required=True)
    parser.add_argument('--gff', type=Path, required=True)
    parser.add_argument('--read1', type=Path, required=True)
    parser.add_argument('--read2', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--index', type=Path, required=True)
    parser.add_argument('--tool-bin', type=Path, required=True,
                        help='Directory containing the external hisat2 executable')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--control-pairs', type=int, default=4000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    genome = read_fasta(args.genome)
    masks, transcripts, introns, space = read_annotation(args.gff, genome)
    simulation = simulate_controls(genome, transcripts, args.output, args.control_pairs)
    datasets = [('ERR13623736_prefix', args.read1, args.read2)]
    for label in ('synthetic_genomic_DNA', 'synthetic_junction_enriched_RNA'):
        datasets.append((label, args.output / label / 'R1.fastq.gz', args.output / label / 'R2.fastq.gz'))
    results = {}
    for label, read1, read2 in datasets:
        directory = args.output / label
        directory.mkdir(exist_ok=True)
        sam = directory / 'alignments.sam'
        command = [str(args.tool_bin / 'hisat2'), '-p', str(args.threads), '--seed', '20260922',
                   '-x', str(args.index), '-1', str(read1), '-2', str(read2), '-S', str(sam)]
        execution = measured_run(command, directory, 'hisat2')
        analysis_start = time.perf_counter()
        summary = analyze_sam(sam, genome, masks, introns, directory)
        summary['analysis_elapsed_seconds'] = time.perf_counter() - analysis_start
        summary['alignment_execution'] = execution
        summary['input_paths'] = [str(read1), str(read2)]
        summary['input_sha256'] = [sha256(read1), sha256(read2)]
        (directory / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        results[label] = summary
        print(label, json.dumps(summary['counts']), flush=True)
    result = {'created_utc': datetime.now(timezone.utc).isoformat(), 'reference_genome': str(args.genome),
              'reference_gff': str(args.gff), 'reference_sha256': [sha256(args.genome), sha256(args.gff)],
              'script_sha256': sha256(__file__), 'annotation_space': space, 'simulation': simulation,
              'datasets': results, 'total_elapsed_seconds_excluding_index': time.perf_counter() - start,
              'alignment_design': 'Genome-only HISAT2 index; splice-aware default; no annotation hints; no forced RNA strand; untrimmed paired prefix reads. No production classification edits.',
              'interpretation_guard': 'Strong supported splice junctions support an RNA component. Absence of junctions does not prove DNA; low mapping, old fragmented assembly, different BioSample, missing annotation and nonrandom prefix limit inference. Synthetic controls are technical only.'}
    (args.output / 'diagnostic.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
