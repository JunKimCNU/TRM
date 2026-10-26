#!/usr/bin/env python3
"""Join exported repeat coordinates to adjacent-anchor alignments, without TRM promotion.

PAF requires minimap2 base-level alignment (-c), recorded source hashes and a
complete successful execution receipt. A reference edge is not a physical end.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import html
import json
from pathlib import Path
import re
import time


DEFAULTS = dict(min_mapq=20, min_identity=0.85, min_query_coverage=0.80,
                min_aligned_bp=200, max_junction_gap=20, end_window=5000,
                read_end_window=50, competing_score_fraction=0.90, long_tract_bp=500)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fasta_lengths(path):
    lengths = {}
    opener = gzip.open if str(path).endswith('.gz') else open
    current = None
    with opener(path, 'rt', encoding='ascii') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                words = line[1:].split()
                if not words or words[0] in lengths:
                    raise ValueError('empty or duplicate FASTA identifier')
                current = words[0]
                lengths[current] = 0
            elif current is None:
                raise ValueError('sequence before FASTA identifier')
            else:
                lengths[current] += len(line)
    if any(n == 0 for n in lengths.values()):
        raise ValueError('empty FASTA sequence')
    return lengths


def parse_paf(line):
    """Validate minimap2 genomic base-alignment counts, without replaying bases.

    In minimap2, nn counts ambiguous columns (including inserted/deleted
    ambiguous bases). These consume CIGAR coordinates but are excluded from
    PAF matches/block length. An '=' can include N/N after --eqx conversion.
    Keep the existing PAF identity denominator and returned row contract.
    """
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) < 12:
        raise ValueError('PAF needs at least 12 columns')
    names = ['query', 'query_length', 'query_start', 'query_end', 'strand', 'target',
             'target_length', 'target_start', 'target_end', 'matches', 'block_length', 'mapq']
    row = dict(zip(names, fields[:12]))
    for key in names:
        if key not in {'query', 'strand', 'target'}:
            if not re.fullmatch(r'[0-9]+', row[key]):
                raise ValueError('invalid PAF integer: ' + key)
            row[key] = int(row[key])
    if (not row['query'] or not row['target'] or row['strand'] not in {'+', '-'} or
            not 0 <= row['query_start'] < row['query_end'] <= row['query_length'] or
            not 0 <= row['target_start'] < row['target_end'] <= row['target_length'] or
            not 0 <= row['matches'] <= row['block_length'] or row['block_length'] <= 0 or
            not 0 <= row['mapq'] <= 255):
        raise ValueError('invalid PAF coordinate/count')
    tags, tag_types = {}, {}
    for field in fields[12:]:
        pieces = field.split(':', 2)
        if (len(pieces) != 3 or not re.fullmatch(r'[A-Za-z][A-Za-z0-9]', pieces[0])
                or pieces[0] in tags):
            raise ValueError('invalid or duplicate PAF tag')
        tags[pieces[0]] = pieces[2]
        tag_types[pieces[0]] = pieces[1]
    if 'cg' not in tags or 'AS' not in tags or tags.get('tp') not in {'P', 'S', 'I', 'i'}:
        raise ValueError('PAF must include minimap2 -c cg/AS/tp alignment tags')
    for key, expected in [('cg', 'Z'), ('AS', 'i'), ('tp', 'A'), ('nn', 'i'), ('NM', 'i')]:
        if key in tags and tag_types[key] != expected:
            raise ValueError('invalid PAF tag type: ' + key)
    if not re.fullmatch(r'[+-]?[0-9]+', tags['AS']):
        raise ValueError('invalid PAF integer tag: AS')
    integer_tags = {}
    for key in ('nn', 'NM'):
        if key in tags:
            if not re.fullmatch(r'[+]?[0-9]+', tags[key]):
                raise ValueError('invalid nonnegative PAF integer tag: ' + key)
            integer_tags[key] = int(tags[key])
    ambiguous = integer_tags.get('nn', 0)
    cigar = re.findall(r'([1-9][0-9]*)([MID=X])', tags['cg'])
    if ''.join(n + op for n, op in cigar) != tags['cg'] or not cigar:
        raise ValueError('unsupported or malformed genomic CIGAR')
    query_span = sum(int(n) for n, op in cigar if op in 'MI=X')
    target_span = sum(int(n) for n, op in cigar if op in 'MD=X')
    columns = sum(int(n) for n, _ in cigar)
    equals = sum(int(n) for n, op in cigar if op == '=')
    minimum_matches = max(0, equals - ambiguous)
    maximum_matches = min(sum(int(n) for n, op in cigar if op in 'M='), row['block_length'])
    if (query_span != row['query_end'] - row['query_start'] or
            target_span != row['target_end'] - row['target_start'] or
            columns != row['block_length'] + ambiguous or
            not minimum_matches <= row['matches'] <= maximum_matches):
        raise ValueError('CIGAR spans disagree with PAF')
    if 'NM' in integer_tags and integer_tags['NM'] != columns - row['matches']:
        raise ValueError('PAF NM disagrees with CIGAR columns and matches')
    row['alignment_score'] = int(tags['AS'])
    row['alignment_type'] = tags['tp']
    row['identity'] = row['matches'] / row['block_length']
    row['aligned_query_bp'] = row['query_end'] - row['query_start']
    row['query_coverage'] = row['aligned_query_bp'] / row['query_length']
    return row


def select_alignment(rows, settings):
    primary = [r for r in rows if r['alignment_type'] == 'P']
    if not rows:
        return None, 'unmapped'
    if len(primary) != 1:
        return None, 'missing_or_multiple_primary_alignments'
    best = primary[0]
    if best['mapq'] == 255 or best['mapq'] < settings['min_mapq']:
        return None, 'low_or_unknown_mapping_quality'
    if best['identity'] < settings['min_identity']:
        return None, 'low_alignment_identity'
    if best['query_coverage'] < settings['min_query_coverage'] or best['aligned_query_bp'] < settings['min_aligned_bp']:
        return None, 'insufficient_anchor_alignment'
    for other in rows:
        if other is best:
            continue
        if other['alignment_score'] >= best['alignment_score'] * settings['competing_score_fraction']:
            return None, 'competitive_or_split_alignment'
    return best, 'qualified_within_reported_alignments'


def anchor_geometry(tract, side, alignment, settings):
    left = side == 'left'
    gap = alignment['query_length'] - alignment['query_end'] if left else alignment['query_start']
    outward_right = left == (alignment['strand'] == '+')
    end = 'right' if outward_right else 'left'
    boundary = alignment['target_end'] if outward_right else alignment['target_start']
    distance = alignment['target_length'] - boundary if outward_right else boundary
    repeat_strand = tract['strand']
    if repeat_strand in {'+', '-'}:
        projected = '+' if repeat_strand == alignment['strand'] else '-'
        polarity = 'expected' if projected == ('+' if outward_right else '-') else 'opposite'
    else:
        projected, polarity = 'ambiguous', 'ambiguous'
    # Use mapped query boundaries; never extrapolate through an unaligned gap.
    outer_read_gap = tract['right_end_distance'] if left else tract['left_end_distance']
    near_junction = gap <= settings['max_junction_gap']
    near_end = distance <= settings['end_window']
    consistent = near_junction and near_end and polarity == 'expected'
    return {'anchor_side': side, 'target': alignment['target'], 'mapping_strand': alignment['strand'],
            'assembly_repeat_strand': projected, 'outward_target_end': end,
            'mapped_junction_boundary0': boundary, 'anchor_boundary_to_assembly_end_bp': distance,
            'unaligned_anchor_junction_bp': gap, 'repeat_outer_read_gap_bp': outer_read_gap,
            'polarity': polarity, 'near_junction': near_junction, 'near_assembly_end': near_end,
            'assembly_end_anchor_concordant': consistent,
            'read_end_and_assembly_end_concordant': consistent and outer_read_gap <= settings['read_end_window'],
            'tract_ge_long_threshold': tract['tract_length'] >= settings['long_tract_bp'],
            **{key: alignment[key] for key in ['mapq', 'identity', 'query_coverage', 'aligned_query_bp',
                                              'target_start', 'target_end', 'query_start', 'query_end']}}


def evaluate_tract(tract, alignments, settings, compartments=None):
    result = dict(tract)
    anchors = []
    for side in ['left', 'right']:
        name = tract.get(side + '_anchor_id', '')
        if not name or name in {'NA', '.', 'none'}:
            continue
        chosen, status = select_alignment(alignments.get(name, []), settings)
        row = {'anchor_id': name, 'anchor_side': side, 'mapping_status': status,
               'reported_alignments': len(alignments.get(name, []))}
        if chosen:
            row.update(anchor_geometry(tract, side, chosen, settings))
            row['target_compartment'] = (compartments or {}).get(chosen['target'], 'unannotated')
            if row['target_compartment'] == 'organelle':
                row['assembly_end_anchor_concordant'] = False
                row['read_end_and_assembly_end_concordant'] = False
        anchors.append(row)
    qualified = [a for a in anchors if a['mapping_status'] == 'qualified_within_reported_alignments']
    concordant = [a for a in qualified if a['assembly_end_anchor_concordant']]
    joint = [a for a in qualified if a['read_end_and_assembly_end_concordant']]
    if len(qualified) == 2:
        status = 'two_sided_anchor_context'
    elif joint:
        status = 'read_end_and_assembly_end_concordant'
    elif concordant:
        status = 'assembly_end_anchor_only'
    elif qualified:
        status = 'qualified_anchor_without_terminal_concordance'
    elif anchors:
        status = 'anchor_mapping_unresolved'
    else:
        status = 'no_exportable_adjacent_anchor'
    result.update(localization_status=status, anchors=anchors, qualified_anchor_count=len(qualified),
                  two_sided_qualified_anchors=len(qualified) == 2,
                  # Two qualified sides remain a distinct geometry; do not count as terminal confirmation.
                  joint_terminal_context=len(qualified) == 1 and bool(joint),
                  joint_long_terminal_context=len(qualified) == 1 and bool(joint) and
                    tract['tract_length'] >= settings['long_tract_bp'],
                  joint_nuclear_terminal_context=len(qualified) == 1 and bool(joint) and
                    joint[0]['target_compartment'] == 'nuclear_chromosome',
                  physical_chromosome_end_confirmed=False)
    return result


def table(path, rows):
    with path.open('w', newline='') as f:
        if rows:
            writer = csv.DictWriter(f, list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)


def same_read_context(outcomes, settings):
    """Describe separate anchor and long-tract events; never merge their intervals.

    Added for interrupted arrays: this table is descriptive and does not change
    the strict same-tract localization flags, classification, or detection power.
    """
    grouped = defaultdict(list)
    for tract in outcomes:
        grouped[tract['read_ordinal']].append(tract)
    rows = []
    for tracts in grouped.values():
        for long_tract in tracts:
            if long_tract['tract_length'] < settings['long_tract_bp']:
                continue
            for anchor_tract in tracts:
                if anchor_tract['tract_id'] == long_tract['tract_id']:
                    continue
                for anchor in anchor_tract['anchors']:
                    if not anchor.get('assembly_end_anchor_concordant'):
                        continue
                    left_anchor = anchor['anchor_side'] == 'left'
                    order = (anchor_tract['end0'] <= long_tract['start0'] if left_anchor
                             else long_tract['end0'] <= anchor_tract['start0'])
                    gap = max(0, long_tract['start0'] - anchor_tract['end0'],
                              anchor_tract['start0'] - long_tract['end0'])
                    rows.append({'read_ordinal': long_tract['read_ordinal'], 'read_id': long_tract['read_id'],
                                 'long_tract_id': long_tract['tract_id'], 'long_start0': long_tract['start0'],
                                 'long_end0': long_tract['end0'], 'long_tract_bp': long_tract['tract_length'],
                                 'anchor_tract_id': anchor_tract['tract_id'], 'anchor_tract_start0': anchor_tract['start0'],
                                 'anchor_tract_end0': anchor_tract['end0'], 'anchor_tract_bp': anchor_tract['tract_length'],
                                 'anchor_id': anchor['anchor_id'], 'target': anchor['target'],
                                 'target_end': anchor['outward_target_end'],
                                 'target_compartment': anchor['target_compartment'],
                                 'intervening_read_gap_bp': gap, 'outward_read_order_consistent': order,
                                 'repeat_strands_agree': long_tract['strand'] == anchor_tract['strand'] and
                                   long_tract['strand'] in {'+', '-'},
                                 'strict_same_tract_joint_long': long_tract['joint_long_terminal_context'],
                                 'tracts_merged': False, 'physical_chromosome_end_confirmed': False})
    return rows


def load_exports(directory, anchor_lengths):
    with (directory / 'tracts.tsv').open(newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    seen, anchors = set(), {}
    for row in rows:
        if row['tract_id'] in seen:
            raise ValueError('duplicate tract identifier')
        seen.add(row['tract_id'])
        for key in ['read_ordinal', 'read_length', 'tract_index', 'start0', 'end0', 'tract_length',
                    'left_end_distance', 'right_end_distance']:
            row[key] = int(row[key])
        row['identity'] = float(row['identity'])
        if (not 0 <= row['start0'] < row['end0'] <= row['read_length'] or
                row['end0'] - row['start0'] != row['tract_length'] or
                row['left_end_distance'] != row['start0'] or
                row['right_end_distance'] != row['read_length'] - row['end0'] or
                row['strand'] not in {'+', '-', 'ambiguous'} or not 0 <= row['identity'] <= 1):
            raise ValueError('invalid tract geometry')
        for side in ['left', 'right']:
            name = row.get(side + '_anchor_id', '')
            if not name or name in {'NA', '.', 'none'}:
                continue
            start, end, length = [int(row[side + '_anchor_' + k]) for k in ['start0', 'end0', 'length']]
            if (name in anchors or not 0 <= start < end <= row['read_length'] or
                    end-start != length or anchor_lengths.get(name) != length or
                    (side == 'left' and end != row['start0']) or
                    (side == 'right' and start != row['end0'])):
                raise ValueError('invalid or non-adjacent anchor geometry')
            anchors[name] = row['tract_id']
    if set(anchors) != set(anchor_lengths):
        raise ValueError('anchor FASTA and tract identifiers differ')
    return rows


def validate_receipt(receipt, paths):
    if receipt.get('status') != 'complete' or receipt.get('exit_code') != 0:
        raise ValueError('alignment execution is incomplete or failed')
    command = receipt.get('command', [])
    if '-c' not in command or '--secondary=yes' not in command:
        raise ValueError('recorded mapping must request base alignments and secondary hits')
    for key, path in paths.items():
        saved = receipt.get(key, {})
        if Path(saved.get('path', '')).resolve() != path.resolve() or saved.get('sha256') != digest(path):
            raise ValueError('alignment source/hash mismatch: ' + key)
    # Requires mapping the exact reference FASTA and exact exported anchors, not an unrelated index.
    if command[-2:] != [str(paths['reference'].resolve()), str(paths['queries'].resolve())]:
        raise ValueError('alignment command input paths differ from receipt')


def validate_export(summary, directory):
    if summary.get('status') != 'complete' or summary.get('complete_input') is not True:
        raise ValueError('tract export must have completed full input')
    for name in ['tracts.tsv', 'anchors.fasta']:
        saved = summary.get('outputs', {}).get(name, {})
        if saved.get('sha256') != digest(directory / name):
            raise ValueError('export output hash mismatch: ' + name)


def assembly_compartments(path, lengths):
    if path is None:
        return {}
    result = {}
    with path.open() as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            row = line.rstrip('\r\n').split('\t')
            if len(row) != 10:
                raise ValueError('invalid NCBI assembly report row')
            location = row[3].lower()
            kind = ('organelle' if location in {'mitochondrion', 'chloroplast', 'plastid', 'apicoplast'}
                    else 'nuclear_chromosome' if location == 'chromosome' else 'other')
            for name in [row[0], row[4], row[6]]:
                if name in lengths:
                    if int(row[8]) != lengths[name] or name in result:
                        # Identical aliases in a single row are allowed.
                        if int(row[8]) != lengths[name] or result[name] != kind:
                            raise ValueError('assembly report sequence identity/length conflict')
                    result[name] = kind
    return result


def add_repeat_context(outcomes, index, annotate):
    """Add reference-coordinate diagnostics without changing any geometry flags."""
    statuses = Counter()
    joint_overlap, joint_no_overlap, joint_incomplete = set(), set(), set()
    for tract in outcomes:
        for anchor in tract['anchors']:
            if anchor['mapping_status'] != 'qualified_within_reported_alignments':
                anchor['repeat_context_status'] = 'not_assessed_unqualified_alignment'
                anchor['unique_nonrepeat_subtelomere_confirmed'] = False
                statuses[anchor['repeat_context_status']] += 1
                continue
            annotation = annotate(index, anchor['target'], anchor['target_start'], anchor['target_end'])
            anchor.update(annotation)
            statuses[annotation['repeat_context_status']] += 1
            if tract['joint_terminal_context'] and anchor['read_end_and_assembly_end_concordant']:
                if annotation['known_target_repeat_overlap_bp']:
                    joint_overlap.add(tract['read_ordinal'])
                elif annotation['audited_span_bp'] == annotation['reference_span_bp']:
                    joint_no_overlap.add(tract['read_ordinal'])
                else:
                    joint_incomplete.add(tract['read_ordinal'])
    return {'provenance': index['provenance'], 'anchor_status_counts': dict(statuses),
            'joint_geometry_records_with_known_repeat_overlap': len(joint_overlap),
            'joint_geometry_records_no_known_overlap_full_audit_scope': len(joint_no_overlap),
            'joint_geometry_records_no_known_overlap_incomplete_scope': len(joint_incomplete),
            'record_groups_can_overlap_across_distinct_tracts': True,
            'unique_nonrepeat_subtelomere_confirmed': False,
            'geometry_counts_changed': False,
            'interpretation': 'Overlap is measured on the reference span, including any reference-consuming alignment gaps, not on query bases. Only the selected family and audited windows are assessed. No observed overlap does not prove unique or nonrepetitive sequence.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['export-dir', 'assembly', 'paf', 'alignment-receipt', 'output-dir']:
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--assembly-report', type=Path)
    parser.add_argument('--repeat-audit', type=Path,
                        help='Optional exact assembly-end audit JSON; annotate selected-family reference overlap without changing geometry')
    for name, value in DEFAULTS.items():
        parser.add_argument('--' + name.replace('_', '-'), type=type(value), default=value)
    args = parser.parse_args()
    settings = {name: getattr(args, name) for name in DEFAULTS}
    if any(v < 0 for v in settings.values()) or any(not 0 < settings[k] <= 1 for k in
            ['min_identity', 'min_query_coverage', 'competing_score_fraction']):
        parser.error('invalid thresholds')
    if args.output_dir.exists():
        parser.error('output directory already exists; keep previous results')
    started = time.monotonic()
    anchors_path = args.export_dir / 'anchors.fasta'
    export_summary = json.loads((args.export_dir / 'summary.json').read_text())
    validate_export(export_summary, args.export_dir)
    receipt = json.loads(args.alignment_receipt.read_text())
    validate_receipt(receipt, {'reference': args.assembly, 'queries': anchors_path, 'paf': args.paf})
    anchor_lengths, assembly_lengths = fasta_lengths(anchors_path), fasta_lengths(args.assembly)
    if not assembly_lengths:
        raise ValueError('assembly FASTA is empty')
    compartments = assembly_compartments(args.assembly_report, assembly_lengths)
    tracts = load_exports(args.export_dir, anchor_lengths)
    exported_counts = {'tracts': len(tracts), 'anchors': len(anchor_lengths),
                       'target_reads': len({t['read_ordinal'] for t in tracts}),
                       'max_tract_bases': max((t['tract_length'] for t in tracts), default=0)}
    if any(export_summary.get(k) != value for k, value in exported_counts.items()):
        raise ValueError('export summary counts differ from parsed coordinates/anchors')
    alignments = defaultdict(list)
    with args.paf.open() as f:
        for line in f:
            if not line.strip():
                continue
            row = parse_paf(line)
            if anchor_lengths.get(row['query']) != row['query_length'] or assembly_lengths.get(row['target']) != row['target_length']:
                raise ValueError('PAF source identifier/length differs from FASTA')
            alignments[row['query']].append(row)
    outcomes = [evaluate_tract(t, alignments, settings, compartments) for t in tracts]
    repeat_context = None
    if args.repeat_audit:
        from anchor_repeat_context import load_repeat_audit, annotate_alignment
        repeat_index = load_repeat_audit(args.repeat_audit, digest(args.assembly),
                                        export_summary['family_key'], assembly_lengths)
        repeat_context = add_repeat_context(outcomes, repeat_index, annotate_alignment)
    read_contexts = same_read_context(outcomes, settings)
    counts = Counter(t['localization_status'] for t in outcomes)
    summary = {'schema_version': 2 if repeat_context is not None else 1, 'program': 'trm-terminal-anchor-validation',
               'created_utc': datetime.now(timezone.utc).isoformat(), 'settings': settings,
               'tracts': len(outcomes), 'exported_anchors': len(anchor_lengths),
               'repeat_supporting_records': len({t['read_ordinal'] for t in outcomes}),
               'status_counts': dict(counts),
               'joint_terminal_context_records': len({t['read_ordinal'] for t in outcomes if t['joint_terminal_context']}),
               'joint_nuclear_terminal_context_records': len({t['read_ordinal'] for t in outcomes if t['joint_nuclear_terminal_context']}),
               'joint_long_context_records': len({t['read_ordinal'] for t in outcomes if t['joint_long_terminal_context']}),
               'long_threshold_bp': settings['long_tract_bp'],
               'long_tracts_with_other_tract_anchor_context': len({r['long_tract_id'] for r in read_contexts}),
               'same_read_context_semantics': 'Descriptive links to other tracts on the same record; intervals are not merged and strict same-tract counts are unchanged.',
               'source_export_summary': export_summary,
               'source_sha256': {str(p.resolve()): digest(p) for p in
                   [args.export_dir / 'summary.json', args.export_dir / 'tracts.tsv', anchors_path,
                    args.assembly, args.paf, args.alignment_receipt] +
                    ([args.assembly_report] if args.assembly_report else []) +
                    ([args.repeat_audit] if args.repeat_audit else [])},
               'physical_chromosome_end_confirmed': False, 'changes_native_classification': False,
               'limitations': [
                   'Scaffold-end proximity plus an adjacent anchor is sequence context, not a physical chromosome-end experiment.',
                   'Same-sample or same-strain reads used in the assembly are not independent biological validation.',
                   'Anchor alignment allows indels; the exported repeat detector still tolerates substitutions only.',
                   'Missing anchors, low mapping quality, unassembled ends and fragmented tracts are unresolved, not biological negatives.',
                   'Reported alignment competitors are bounded by minimap2 settings; MAPQ is not proof of global uniqueness.',
                   'Two qualified flanks remain a separate context, including interrupted arrays; not automatically interstitial repeats.',
                   'Mapped boundary distance is not extrapolated through unaligned bases; terminal window is a declared heuristic.',
                   'No automatic catalogue, main-TRM, novelty, ALT, absence, karyotype or adequacy update.'
               ]}
    if repeat_context is not None:
        summary['reference_repeat_context'] = repeat_context
        summary['limitations'].append('Optional repeat overlap is exact-family annotation within recorded end windows; zero overlap is not a unique nonrepetitive subtelomere or absence claim.')
    args.output_dir.mkdir(parents=True)
    table(args.output_dir / 'tracts.tsv', [{k: v for k, v in row.items() if k != 'anchors'} for row in outcomes])
    table(args.output_dir / 'same_read_context.tsv', read_contexts)
    (args.output_dir / 'same_read_context.json').write_text(json.dumps(read_contexts, indent=2) + '\n')
    anchor_rows = [{'tract_id': row['tract_id'], **a} for row in outcomes for a in row['anchors']]
    if anchor_rows:
        keys = list(dict.fromkeys(k for row in anchor_rows for k in row))
        table(args.output_dir / 'anchors.tsv', [{k: row.get(k) for k in keys} for row in anchor_rows])
    summary['wall_seconds'] = time.monotonic() - started
    (args.output_dir / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')
    (args.output_dir / 'details.json').write_text(json.dumps(outcomes, ensure_ascii=False, indent=2) + '\n')
    content = ''.join('<tr>' + ''.join('<td>' + html.escape(str(row[k])) + '</td>' for k in
                     ['read_id', 'tract_length', 'strand', 'localization_status', 'joint_long_terminal_context']) + '</tr>'
                      for row in outcomes)
    repeat_html = ''
    if repeat_context is not None:
        context_rows = [row for row in anchor_rows if row['mapping_status'] == 'qualified_within_reported_alignments']
        context_table = ''.join('<tr>' + ''.join('<td>' + html.escape(str(row[key])) + '</td>' for key in
                               ['anchor_id', 'target', 'target_start', 'target_end', 'known_target_repeat_overlap_bp',
                                'reference_span_bp', 'audited_span_bp', 'repeat_context_status']) + '</tr>'
                                for row in context_rows)
        repeat_html = ('<h2>인접 anchor의 반복 중첩</h2><p>정렬의 reference 범위에 알려진 표적 반복이 겹치는지 별도로 확인했습니다. '
                       '검사 범위는 기록된 assembly 말단 window와 해당 family입니다. 중첩 0은 고유한 비반복 서열의 증명이 아닙니다. '
                       '아래 값으로 기존 위치·방향 계수를 바꾸지 않습니다.</p><p>기하 조건 지지 read 중 반복 중첩 관측: '
                       + str(repeat_context['joint_geometry_records_with_known_repeat_overlap']) + '; 검사 범위 전체에 속하고 중첩 미관측: '
                       + str(repeat_context['joint_geometry_records_no_known_overlap_full_audit_scope'])
                       + '. 같은 read의 서로 다른 구간이 두 그룹에 속할 수 있으므로 합산하지 않습니다.</p>'
                       '<table><tr><th>Anchor</th><th>Reference</th><th>Start0</th><th>End0</th><th>중첩 bp</th>'
                       '<th>Reference 범위 bp</th><th>검사 범위 bp</th><th>반복 맥락</th></tr>' + context_table + '</table>')
    (args.output_dir / 'index.html').write_text('<!doctype html><html lang="ko"><meta charset="utf-8">'
        '<title>TRM 인접 서열의 assembly 위치 검증</title><style>body{font:15px/1.6 system-ui;margin:30px}td,th{padding:7px;border-bottom:1px solid #ddd;text-align:left}</style>'
        '<h1>반복과 인접 서열의 assembly 위치</h1><p>Long-read 반복과 인접 anchor의 위치·방향을 연결했습니다. '
        'Assembly 끝은 실제 염색체 말단의 확정 증거가 아닙니다. 기존 TRM 판정과 검출력을 변경하지 않습니다.</p>'
        '<p><a href="summary.json">설정·요약·한계</a> · <a href="tracts.tsv">전체 구간 TSV</a> · '
        '<a href="details.json">정렬 상세 JSON</a> · <a href="same_read_context.json">같은 리드의 다른 구간 연결 근거</a></p>'
        + repeat_html + '<h2>전체 검출 구간</h2><table><tr><th>Read</th><th>Tract bp</th>'
        '<th>G-strand 방향</th><th>위치 근거</th><th>긴 구간·방향·말단 인접 동시 지지</th></tr>' + content + '</table></html>')
    print(json.dumps({k: summary[k] for k in ['tracts', 'status_counts', 'joint_terminal_context_records',
                                            'joint_long_context_records', 'wall_seconds']}))


if __name__ == '__main__':
    main()
