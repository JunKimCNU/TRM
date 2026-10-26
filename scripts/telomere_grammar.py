"""Exact finite-block language evidence, preserving published block phase.

This module explains supplied sequence windows. It does not classify telomeres,
bridge gaps, normalize each block's phase, or permit sequencing edits.
"""
from collections import defaultdict
import hashlib
import re

VERSION = '0.1.0'
COMPLEMENT = bytes.maketrans(b'ACGTRYSWKMBDHVN', b'TGCAYRSWMKVHDBN')


def natural(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(name + ' must be an integer >= ' + str(minimum))
    return value


def validate_units(units):
    # Contextual grammar blocks are not primitive, rotation-normalized motifs.
    if not isinstance(units, list) or not 1 <= len(units) <= 64:
        raise ValueError('explicit list of 1..64 phase-preserved blocks required')
    if any(not isinstance(u, str) or not re.fullmatch('[ACGT]{1,64}', u) for u in units):
        raise ValueError('each grammar block must be uppercase ACGT, length 1..64')
    if len(set(units)) != len(units):
        raise ValueError('duplicate literal grammar blocks')
    return [u.encode('ascii') for u in units]


def solve_paths(n, edges):
    """Return maximal and all earliest-start endpoint states of a labeled DAG.

    Each edge is (start, end, unit_index). All consume positive sequence length.
    At a fixed earliest start/end, min/max block counts and all-parse ambiguity
    are retained. The representative minimizes block count; it is not a claim
    about biological segmentation. Distinct labels are distinct interpretations.
    """
    natural(n, 'length')
    by_end, seen = defaultdict(list), set()
    for edge in edges:
        if len(edge) != 3:
            raise ValueError('three edge values required')
        a, b, unit = edge
        for v in edge:
            natural(v, 'edge value')
        if not 0 <= a < b <= n or (a, b, unit) in seen:
            raise ValueError('invalid or duplicate edge')
        seen.add((a, b, unit))
        by_end[b].append((a, b, unit))
    states = [None] * (n + 1)
    previous_edge = [None] * (n + 1)
    for end in range(1, n + 1):
        best, chosen = None, None
        for a, b, unit in sorted(by_end.get(end, [])):
            parent = states[a]
            candidate = {
                'start0': parent['start0'] if parent else a, 'end0': b,
                'minimum_block_count': 1 + (parent['minimum_block_count'] if parent else 0),
                'maximum_block_count': 1 + (parent['maximum_block_count'] if parent else 0),
                'path_count_capped2': parent['path_count_capped2'] if parent else 1,
            }
            edge = (a, b, unit)
            if best is None or candidate['start0'] < best['start0']:
                best, chosen = candidate, edge
            elif candidate['start0'] == best['start0']:
                count = min(2, best['path_count_capped2'] + candidate['path_count_capped2'])
                maximum = max(best['maximum_block_count'], candidate['maximum_block_count'])
                if (candidate['minimum_block_count'], edge) < (best['minimum_block_count'], chosen):
                    best, chosen = candidate, edge
                best['maximum_block_count'] = maximum
                best['path_count_capped2'] = count
        if best is not None:
            best['span_bp'] = end - best['start0']
            states[end], previous_edge[end] = best, chosen
    endpoints = [dict(s) for s in states if s is not None]
    maximal, furthest = [], -1
    for state in sorted(endpoints, key=lambda r: (r['start0'], -r['end0'])):
        if state['end0'] <= furthest:
            continue
        row = dict(state)
        trace, end = [], row['end0']
        while end > row['start0']:
            edge = previous_edge[end]
            if edge is None:
                raise AssertionError('missing exact path predecessor')
            trace.append(list(edge)); end = edge[0]
        if end != row['start0'] or len(trace) != row['minimum_block_count']:
            raise AssertionError('representative path inconsistent')
        row['representative_tokens'] = list(reversed(trace))
        maximal.append(row); furthest = row['end0']
    return maximal, endpoints


def analyze_window(sequence, units, window, min_span=100, min_blocks=4,
                   end_distance=50, max_edges=200000):
    for name, val in [('min_span', min_span), ('min_blocks', min_blocks), ('max_edges', max_edges)]:
        natural(val, name, 1)
    natural(end_distance, 'end_distance')
    if not isinstance(sequence, bytes) or not sequence or not re.fullmatch(b'[ACGTRYSWKMBDHVN]+', sequence):
        raise ValueError('nonempty uppercase observed DNA bytes required')
    blocks = validate_units(units)
    n = len(sequence)
    a, b, length, strand = (window[k] for k in ['start0', 'end0', 'record_length', 'strand'])
    for val in (a, b, length):
        natural(val, 'window coordinate')
    if not 0 <= a < b <= length or b-a != n or strand not in ('+', '-'):
        raise ValueError('window coordinates/strand disagree with supplied sequence')
    role = window.get('role')
    if role == 'deposited_terminal_window':
        if (a == 0) == (b == length) or strand != ('-' if a == 0 else '+'):
            raise ValueError('terminal window must face outward at one original end')
    elif role == 'coordinate_internal_unlabelled':
        if a == 0 or b == length:
            raise ValueError('internal window cannot touch an original end')
    else:
        raise ValueError('explicit supported window role required')

    def location(start, end):
        x, y = (b-end, b-start) if strand == '-' else (a+start, a+end)
        return {'start0': start, 'end0': end, 'original_start0': x, 'original_end0': y}

    edges = []
    for i, unit in enumerate(blocks):
        pos = -1
        while True:
            pos = sequence.find(unit, pos + 1)
            if pos < 0:
                break
            if len(edges) >= max_edges:
                raise ValueError('exact token limit exceeded; no complete result')
            edges.append((pos, pos+len(unit), i))
    edges.sort()
    paths, _ = solve_paths(n, edges)
    tokens = [{**location(x, y), 'unit_index': i} for x, y, i in edges]
    covered = bytearray(n)
    by_start, by_end = defaultdict(list), defaultdict(list)
    for x, y, i in edges:
        covered[x:y] = b'\x01' * (y-x)
        by_start[x].append([x, y, i]); by_end[y].append([x, y, i])
    segments = []
    for pos, value in enumerate(covered):
        kind = 'literal' if value else 'uncovered'
        if segments and segments[-1]['class'] == kind:
            segments[-1]['end0'] = pos + 1
        else:
            segments.append({'start0': pos, 'end0': pos+1, 'class': kind})
    for segment in segments:
        x, y = segment['start0'], segment['end0']
        segment.update(location(x, y)); segment['span_bp'] = y-x
        if segment['class'] == 'uncovered':
            segment['observed_source_slice'] = sequence[x:y].decode()
    is_terminal = window['role'] == 'deposited_terminal_window'
    for path in paths:
        x, y = path['start0'], path['end0']
        path.update(location(x, y))
        # Distance always refers to the ORIGINAL deposited end, never the crop.
        distance = (path['original_start0'] if a == 0 else length-path['original_end0']) if is_terminal else None
        path['original_end_distance_bp'] = distance
        path['meets_local_gate'] = path['span_bp'] >= min_span and path['minimum_block_count'] >= min_blocks
        path['near_original_end'] = distance <= end_distance if distance is not None else None
        path['span_at_least_500'] = path['span_bp'] >= 500
        path['next_base_covered'] = bool(covered[y]) if y < n else None
        path['termination'] = 'window_edge' if y == n else 'no_block_starts_at_exact_boundary'
        if y < n and by_start[y]:
            raise AssertionError('maximal path has an exact continuation')
    fragments = []
    for i, unit in enumerate(blocks):
        for side in ('left', 'right'):
            for size in range(1, min(n, len(unit)-1)+1):
                x, y = (0, size) if side == 'left' else (n-size, n)
                target = unit[-size:] if side == 'left' else unit[:size]
                if sequence[x:y] != target:
                    continue
                original_side = side if strand == '+' else ('right' if side == 'left' else 'left')
                genuine = (a == 0) if original_side == 'left' else (b == length)
                adjacent = by_start[y] if side == 'left' else by_end[x]
                fragments.append({**location(x, y), 'unit_index': i, 'oriented_side': side,
                    'original_side': original_side, 'boundary_kind': 'original_record_edge' if genuine else 'artificial_crop',
                    'adjacent_complete_tokens': adjacent, 'observed_source_slice': sequence[x:y].decode(),
                    'same_interval_complete_tokens': [edge for edge in by_start[x] if edge[1] == y],
                    'extends_path': False})
    near = [p for p in paths if p['near_original_end'] is True]
    summary = {'token_count': len(tokens), 'literal_union_bp': sum(covered),
               'uncovered_bp': n-sum(covered), 'maximal_paths': len(paths),
               'max_path_bp': max((p['span_bp'] for p in paths), default=0),
               'near_end_max_path_bp': max((p['span_bp'] for p in near), default=0) if is_terminal else None,
               'has_local_gate_near_end': any(p['meets_local_gate'] for p in near) if is_terminal else None,
               'has_500bp_path': any(p['span_at_least_500'] for p in paths),
               'has_500bp_path_near_end': any(p['span_at_least_500'] for p in near) if is_terminal else None,
               'ambiguous_maximal_paths': sum(p['path_count_capped2'] == 2 for p in paths),
               'boundary_fragments': len(fragments)}
    return {'window': window, 'units': units, 'oriented_slice_sha256': hashlib.sha256(sequence).hexdigest(),
            'tokens': tokens, 'coverage_segments': segments, 'paths': paths,
            'boundary_fragments': fragments, 'summary': summary}
