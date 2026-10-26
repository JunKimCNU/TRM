"""Explain observed sequence under a fixed vocabulary without joining gaps.

This evidence layer never changes the irregular detector's acceptance rules.
Intervals, coverage, local explanations and full paths are different objects.
"""
from bisect import bisect_right
from collections import defaultdict
import hashlib

import audit_irregular_arrays as audit

VERSION = '0.1.0'
COMPLEMENT = bytes.maketrans(b'ACGTRYSWKMBDHVN', b'TGCAYRSWMKVHDBN')


def natural(value, label, minimum=0):
    return audit.integer(value, label, minimum)


def window_coordinates(start, end, window_start, window_end, strand):
    for key, value in [('start', start), ('end', end), ('window_start', window_start), ('window_end', window_end)]:
        natural(value, key)
    if strand not in ('+', '-') or not 0 <= start < end <= window_end - window_start:
        raise ValueError('invalid window coordinates or strand')
    return (window_end - end, window_end - start) if strand == '-' else (window_start + start, window_start + end)


def solve_evidence(n, edges):
    """Earliest-start DAG, then minimum edits; count assignments only up to two.

    Edges: (start,end,min_cost,min_cost_variant_count,all_variant_count).
    A state retains counts for that SAME earliest-start span. Possible paths
    include higher-cost local explanations; optimal paths minimize total cost.
    All transitions must consume positive query length. There are no gap edges.
    """
    natural(n, 'window length')
    by_end, seen = defaultdict(list), set()
    for edge in edges:
        if len(edge) != 5:
            raise ValueError('edge requires five integers')
        a, b, cost, optimal, possible = edge
        for key, value in zip(('start', 'end', 'cost', 'optimal', 'possible'), edge):
            natural(value, key)
        if not 0 <= a < b <= n or not 1 <= optimal <= possible or (a, b) in seen:
            raise ValueError('invalid/duplicate edge or explanation counts')
        seen.add((a, b))
        by_end[b].append(edge)
    states = [None] * (n + 1)
    endpoints = []
    for end in range(1, n + 1):
        best = None
        for start, _, cost, optimal, possible in by_end.get(end, ()):
            previous = states[start]
            row = {'start0': previous['start0'] if previous else start, 'end0': end,
                   'minimum_edit_count': cost + (previous['minimum_edit_count'] if previous else 0),
                   'optimal_path_count_capped2': min(2, optimal * (previous['optimal_path_count_capped2'] if previous else 1)),
                   'possible_path_count_capped2': min(2, possible * (previous['possible_path_count_capped2'] if previous else 1))}
            if best is None or row['start0'] < best['start0']:
                best = row
            elif row['start0'] == best['start0']:
                all_ways = min(2, best['possible_path_count_capped2'] + row['possible_path_count_capped2'])
                if row['minimum_edit_count'] < best['minimum_edit_count']:
                    best = row
                elif row['minimum_edit_count'] == best['minimum_edit_count']:
                    best['optimal_path_count_capped2'] = min(2, best['optimal_path_count_capped2'] + row['optimal_path_count_capped2'])
                best['possible_path_count_capped2'] = all_ways
        if best is not None:
            states[end] = best
            endpoints.append(best)
    maximal, furthest = [], -1
    for row in sorted(endpoints, key=lambda x: (x['start0'], -x['end0'], x['minimum_edit_count'])):
        if row['end0'] > furthest:
            maximal.append(dict(row))
            furthest = row['end0']
    return maximal, endpoints


def coverage_segments(n, literal_intervals, token_intervals):
    """Partition bases by evidence coverage, never by biological repeat units."""
    natural(n, 'window length')
    arrays = []
    for intervals in (literal_intervals, token_intervals):
        delta = [0] * (n + 1)
        for a, b in intervals:
            natural(a, 'interval start'); natural(b, 'interval end')
            if not 0 <= a < b <= n:
                raise ValueError('coverage interval outside window')
            delta[a] += 1; delta[b] -= 1
        arrays.append(delta)
    exact = token = 0
    rows = []
    for pos in range(n):
        exact += arrays[0][pos]; token += arrays[1][pos]
        if exact and not token:
            raise ValueError('literal coverage must be contained in token coverage')
        kind = 'literal' if exact else 'edited_only' if token else 'uncovered'
        if rows and rows[-1]['class'] == kind:
            rows[-1]['end0'] = pos + 1
        else:
            rows.append({'start0': pos, 'end0': pos + 1, 'class': kind})
    return rows


def tokenize(sequence, units, max_edges=200000, max_explanations=1000000):
    """Enumerate explanations of observed slices, without synthesizing reads."""
    lengths = sorted({len(u) + d for u in units for d in ((-1, 0, 1) if len(u) >= 5 else (0,))})
    rows, cache, total = [], {}, 0
    for end in range(1, len(sequence) + 1):
        for size in lengths:
            start = end - size
            if start < 0:
                continue
            observed = sequence[start:end]
            if observed not in cache:
                cache[observed] = audit.explanations(observed, units) if audit.ACGT.fullmatch(observed) else []
            choices = cache[observed]
            if not choices:
                continue
            if len(rows) >= max_edges or total + len(choices) > max_explanations:
                raise ValueError('token/explanation computation limit exceeded; no complete result')
            total += len(choices)
            best = min(x[0] for x in choices)
            explanations = [{'edit_count': c, 'unit_index': u, 'edit': op, 'edit_offset': offset}
                            for c, u, op, offset in choices]
            rows.append({'start0': start, 'end0': end, 'minimum_edit_count': best,
                         'minimum_cost_variant_count': sum(x[0] == best for x in choices),
                         'all_variant_count': len(choices), 'explanations': explanations,
                         'unit_assignment_ambiguous': len({x[1] for x in choices}) > 1,
                         'unit_or_operation_ambiguous': len({(x[1], x[2]) for x in choices}) > 1,
                         'edit_position_ambiguous': any(sum(x[1:3] == key for x in choices) > 1
                                                        for key in {x[1:3] for x in choices})})
    return rows


def validated_units(units, allow_empty=False):
    if not isinstance(units, (list, tuple)) or (not units and not allow_empty):
        raise ValueError('explicit unit list required')
    if not units:
        return []
    strings = [u.decode('ascii') if isinstance(u, bytes) else u for u in units]
    if any(not isinstance(u, str) for u in strings):
        raise ValueError('unit must be an ASCII DNA string')
    audit.definitions_from_summary({'definitions': [{'id': 'validated_units', 'units': sorted(strings)}]})
    return [u.encode('ascii') for u in strings]


def explain_window(sequence, primary_units, explanatory_units, window, min_span=100,
                   min_exact_anchors=4, max_edges=200000, max_explanations=1000000):
    """Sequence is already in the declared window orientation (not synthesized)."""
    for key, value in [('min_span', min_span), ('min_exact_anchors', min_exact_anchors),
                       ('max_edges', max_edges), ('max_explanations', max_explanations)]:
        natural(value, key, 1)
    for key in ('record_length', 'start0', 'end0'):
        natural(window[key], key)
    if not isinstance(sequence, bytes) or not sequence:
        raise ValueError('nonempty observed bytes required')
    n = len(sequence)
    if (n != window['end0'] - window['start0'] or not audit.DNA.fullmatch(sequence) or
            not 0 <= window['start0'] < window['end0'] <= window['record_length'] or
            window['strand'] not in ('+', '-')):
        raise ValueError('window source length/alphabet mismatch')
    sequence = sequence.upper()
    primary = validated_units(primary_units)
    hypotheses = validated_units(explanatory_units, allow_empty=True)
    tokens = tokenize(sequence, primary, max_edges, max_explanations)
    literal = [(t['start0'], t['end0']) for t in tokens if t['minimum_edit_count'] == 0]
    intervals = [(t['start0'], t['end0']) for t in tokens]
    edges = [(t['start0'], t['end0'], t['minimum_edit_count'], t['minimum_cost_variant_count'], t['all_variant_count']) for t in tokens]
    maximal, endpoints = solve_evidence(n, edges)
    segment_rows = coverage_segments(n, literal, intervals)
    token_coverage = bytearray(n)
    for segment in segment_rows:
        if segment['class'] != 'uncovered':
            a, b = segment['start0'], segment['end0']
            token_coverage[a:b] = b'\x01' * (b-a)

    def locate(a, b):
        x, y = window_coordinates(a, b, window['start0'], window['end0'], window['strand'])
        return {'start0': a, 'end0': b, 'original_start0': x, 'original_end0': y, 'span_bp': b-a}

    starts = sorted({a for a, _ in intervals})
    left_crop = (window['start0'] != 0) if window['strand'] == '+' else (window['end0'] != window['record_length'])
    right_crop = (window['end0'] != window['record_length']) if window['strand'] == '+' else (window['start0'] != 0)
    maximum_token_length = max(map(len, primary)) + 1
    paths = []
    for row in maximal:
        a, b = row['start0'], row['end0']
        hit_count = audit.literal_hits(sequence[a:b], primary)
        next_index = bisect_right(starts, b)
        next_start = starts[next_index] if next_index < len(starts) else None
        x = {**row, **locate(a, b), 'literal_gate_hits': hit_count,
             'meets_existing_local_span_and_literal_gate': b-a >= min_span and hit_count >= min_exact_anchors,
             'successor_status': 'window_edge' if b == n else 'no_token_starts_at_exact_boundary',
             'next_token_start0': next_start,
             'next_base_has_token_coverage': bool(token_coverage[b]) if b < n else None,
             'left_crop_context_may_be_missing': left_crop and a < maximum_token_length,
             'right_crop_context_may_be_missing': right_crop and n-b < maximum_token_length,
             'intervals_may_overlap_do_not_sum_spans': True}
        if b < n and b in starts:
            raise AssertionError('maximal path cannot have an extendable successor')
        paths.append(x)
    for t in tokens:
        t.update(locate(t['start0'], t['end0']))
    for s in segment_rows:
        s.update(locate(s['start0'], s['end0']))
        s['basis'] = 'any_complete_primary_unit_token_coverage_not_contiguous_path'
        if s['class'] == 'uncovered':
            s['observed_source_slice'] = sequence[s['start0']:s['end0']].decode()
            s['contains_ambiguous_bases'] = not bool(audit.ACGT.fullmatch(sequence[s['start0']:s['end0']]))
    overlays = []
    for i, unit in enumerate(hypotheses):
        start = -1
        while True:
            start = sequence.find(unit, start + 1)
            if start < 0:
                break
            alternatives = audit.explanations(sequence[start:start+len(unit)], primary)
            overlays.append({**locate(start, start+len(unit)), 'hypothesis_unit_index': i,
                'hypothesis_literal': unit.decode(), 'hypothesis_role': 'separate_literal_not_in_primary_graph',
                'primary_same_interval_explanations': [
                    {'edit_count': c, 'unit_index': u, 'edit': op, 'edit_offset': offset}
                    for c, u, op, offset in alternatives],
                'local_unit_length_nonidentifiable': any(len(primary[u]) != len(unit) for _, u, _, _ in alternatives)})
    fragments = []
    for vocabulary, units in [('primary', primary), ('explanatory_literal', hypotheses)]:
        for i, unit in enumerate(units):
            for edge in ('left', 'right'):
                for size in range(1, min(n, len(unit)-1)+1):
                    a, b = (0, size) if edge == 'left' else (n-size, n)
                    target = unit[-size:] if edge == 'left' else unit[:size]
                    if sequence[a:b] != target:
                        continue
                    anchors = []
                    anchor_units = primary if vocabulary == 'primary' else [unit]
                    for anchor_unit in anchor_units:
                        c, d = (b, b+len(anchor_unit)) if edge == 'left' else (a-len(anchor_unit), a)
                        if 0 <= c < d <= n and sequence[c:d] == anchor_unit:
                            anchors.append({**locate(c, d), 'literal': anchor_unit.decode()})
                    original_edge = ('left' if window['strand'] == '+' else 'right') if edge == 'left' else ('right' if window['strand'] == '+' else 'left')
                    original_boundary = (window['start0'] if edge == 'left' else window['end0']) if window['strand'] == '+' else (window['end0'] if edge == 'left' else window['start0'])
                    genuine_record_edge = original_boundary in (0, window['record_length'])
                    complete_alternatives = audit.explanations(sequence[a:b], primary)
                    fragments.append({**locate(a, b), 'oriented_edge': edge, 'original_edge_side': original_edge,
                        'boundary_kind': 'original_record_edge' if genuine_record_edge else 'artificial_crop',
                        'vocabulary': vocabulary, 'unit_index': i, 'unit': unit.decode(),
                        'observed_source_slice': sequence[a:b].decode(), 'adjacent_literal_anchors': anchors,
                        'anchored_original_record_fragment': genuine_record_edge and bool(anchors),
                        'complete_primary_token_explanations_also_possible': [
                            {'edit_count': c, 'unit_index': u, 'edit': op, 'edit_offset': offset}
                            for c, u, op, offset in complete_alternatives],
                        'extends_detector_array': False, 'physical_telomere_confirmed': False})
    counts = {kind: sum(s['span_bp'] for s in segment_rows if s['class'] == kind)
              for kind in ('literal', 'edited_only', 'uncovered')}
    return {'window': window, 'oriented_slice_sha256': hashlib.sha256(sequence).hexdigest(),
        'primary_units': [u.decode() for u in primary], 'explanatory_units': [u.decode() for u in hypotheses],
        'tokens': tokens, 'coverage_segments': segment_rows, 'maximal_complete_paths': paths,
        'literal_hypotheses': overlays, 'boundary_fragments': fragments,
        'summary': {'window_bp': n, 'primary_token_intervals': len(tokens),
            'primary_local_explanations': sum(t['all_variant_count'] for t in tokens),
            'local_unit_ambiguous_intervals': sum(t['unit_assignment_ambiguous'] for t in tokens),
            'local_edit_position_ambiguous_intervals': sum(t['edit_position_ambiguous'] for t in tokens),
            'literal_union_bp': counts['literal'], 'edited_only_union_bp': counts['edited_only'],
            'allowed_token_union_bp': counts['literal']+counts['edited_only'], 'uncovered_bp': counts['uncovered'],
            'max_complete_path_bp': max((p['span_bp'] for p in paths), default=0),
            'maximal_complete_paths': len(paths),
            'optimal_parse_ambiguous_paths': sum(p['optimal_path_count_capped2'] == 2 for p in paths),
            'any_parse_ambiguous_paths': sum(p['possible_path_count_capped2'] == 2 for p in paths),
            'boundary_stops_inside_token_coverage': sum(p['next_base_has_token_coverage'] is True for p in paths),
            'explanatory_literal_occurrences': len(overlays),
            'local_length_nonidentifiable_occurrences': sum(x['local_unit_length_nonidentifiable'] for x in overlays),
            'anchored_original_record_fragments': sum(f['anchored_original_record_fragment'] for f in fragments),
            'anchored_artificial_crop_fragments': sum(f['boundary_kind'] == 'artificial_crop' and bool(f['adjacent_literal_anchors']) for f in fragments),
            'endpoint_states': len(endpoints)},
        'changes_detector_arrays': False, 'changes_biological_classification': False,
        'orientation_is_requested_not_joint_strand_disambiguation': True,
        'left_window_boundary_is_artificial_crop': left_crop,
        'right_window_boundary_is_artificial_crop': right_crop,
        'complete_tokens_crossing_crop_boundaries_not_scored': left_crop or right_crop,
        'path_counts_include_unit_operation_and_edit_position_assignments': True,
        'model_span_not_exact_repeat_length': True, 'full_path_sequences_not_enumerated': True}
