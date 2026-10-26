"""All optimal unit-cost paths against a cyclic repeat, as an optional diagnostic.

Both template strands and every initial phase are considered. Every query base
is consumed; subsequent phase changes require a paid edit. No traceback is
selected and no base quality participates in optimization. This module does
not classify variants, errors, telomeres, homology, or biological origin.
"""
from __future__ import annotations

SCHEMA_VERSION = 1
MAX_REGION_BP = 1000
MAX_MOTIF_BP = 64
MAX_TOTAL_STATES = 200_000
METHOD = 'global_unit_cost_cyclic_edit_all_optimal_paths'
OPERATIONS = ('match', 'substitution', 'insertion', 'deletion')
QUERY_OPERATIONS = OPERATIONS[:3]
_DNA = frozenset('ACGT')
_COMPLEMENT = str.maketrans('ACGT', 'TGCA')
LIMITATIONS = [
    'Distances and operation constraints describe a unit-cost sequence model, not sequencing-error probabilities or evolutionary homology.',
    'A selected native tract is an ascertained interval; preference for its selected family is not an independent test of variant accuracy.',
    'Both strands, every free initial phase, and all globally optimal paths are retained; phase ties are not independent evidence.',
    'Operation-count minima and maxima can come from different optimal paths and must not be combined into a single inferred alignment.',
    'Possible edit positions can be alternatives for one edit; their count is not the number of required edits.',
    'Only observed query bases have quality characters. Template-deletion gaps have no observed base and quality is null.',
    'Phred+33 values are conditional on unverified encoding and calibration; high values do not establish biological accuracy.',
    'No native count, candidate threshold, classification, adequacy, main TRM, or novel TRM is changed.',
]


def _primitive(motif):
    for period in range(1, len(motif) + 1):
        if len(motif) % period == 0 and motif == motif[:period] * (len(motif) // period):
            return motif[:period]
    raise AssertionError('nonempty bounded motif required')


def _forward(sequence, template):
    """Shortest distances from any initial phase; two sweeps close each cycle."""
    k = len(template)
    rows = [[0] * k]
    for base in sequence:
        previous = rows[-1]
        current = [min(previous[p] + 1,
                       previous[p - 1] + (base != template[p - 1])) for p in range(k)]
        for edge in range(2 * k):
            p = edge % k
            q = (p + 1) % k
            current[q] = min(current[q], current[p] + 1)
        rows.append(current)
    return rows


def _backward(sequence, template):
    """Shortest distances to any ending phase, with reverse deletion closure."""
    n, k = len(sequence), len(template)
    rows = [None] * (n + 1)
    rows[n] = [0] * k
    for i in range(n - 1, -1, -1):
        following = rows[i + 1]
        current = [min(following[p] + 1,
                       following[(p + 1) % k] + (sequence[i] != template[p])) for p in range(k)]
        for edge in range(2 * k - 1, -1, -1):
            p = edge % k
            current[p] = min(current[p], current[(p + 1) % k] + 1)
        rows[i] = current
    return rows


def _merge_vector(target, source, operation_index, choose):
    candidate = [v + (i == operation_index) for i, v in enumerate(source)]
    if target is None:
        return candidate
    return [choose(a, b) for a, b in zip(target, candidate)]


def _operation_ranges(n, k, query_edges, deletions, starts, ends):
    """Four independent min/max objectives on the complete optimal-edge DAG."""
    low = [[0] * 4 if p in starts else None for p in range(k)]
    high = [[0] * 4 if p in starts else None for p in range(k)]
    for i in range(n + 1):
        # A full positive-cost deletion cycle cannot lie on a shortest path.
        assert len(deletions[i]) < k
        for edge in range(2 * k):
            p = edge % k
            if p not in deletions[i] or low[p] is None:
                continue
            q = (p + 1) % k
            low[q] = _merge_vector(low[q], low[p], 3, min)
            high[q] = _merge_vector(high[q], high[p], 3, max)
        if i == n:
            break
        next_low, next_high = [None] * k, [None] * k
        for p, q, operation in query_edges[i]:
            if low[p] is None:
                raise AssertionError('optimal query edge lacks reachable prefix')
            op = OPERATIONS.index(operation)
            next_low[q] = _merge_vector(next_low[q], low[p], op, min)
            next_high[q] = _merge_vector(next_high[q], high[p], op, max)
        low, high = next_low, next_high
    if any(low[p] is None for p in ends):
        raise AssertionError('optimal ending phase not reachable')
    return {op: {'min': min(low[p][j] for p in ends),
                 'max': max(high[p][j] for p in ends)} for j, op in enumerate(OPERATIONS)}


def _boundary_deletion_ranges(n, k, query_edges, deletions, starts, ends):
    """Min/max deletions during each row visit, without enumerating tracebacks."""
    ranges = []
    for i in range(n + 1):
        entries = starts if i == 0 else {q for _, q, _ in query_edges[i - 1]}
        exits = ends if i == n else {p for p, _, _ in query_edges[i]}
        low = [0 if p in entries else None for p in range(k)]
        high = [0 if p in entries else None for p in range(k)]
        for edge in range(2 * k):
            p = edge % k
            if p not in deletions[i] or low[p] is None:
                continue
            q = (p + 1) % k
            low[q] = low[p] + 1 if low[q] is None else min(low[q], low[p] + 1)
            high[q] = high[p] + 1 if high[q] is None else max(high[q], high[p] + 1)
        assert exits and all(low[p] is not None for p in exits)
        ranges.append({'min': min(low[p] for p in exits), 'max': max(high[p] for p in exits)})
    return ranges


def _optimal_graph(sequence, template, forward, backward, cost, strand):
    n, k = len(sequence), len(template)
    starts = {p for p in range(k) if backward[0][p] == cost}
    ends = {p for p in range(k) if forward[n][p] == cost}
    query_edges = [[] for _ in range(n)]
    deletions = [set() for _ in range(n + 1)]
    edges = []
    query_ops = [set() for _ in range(n)]
    query_bases = [set() for _ in range(n)]
    deletion_bases = [set() for _ in range(n + 1)]

    def add(i, p, j, q, operation, edge_cost, template_base, query_base):
        if forward[i][p] + edge_cost + backward[j][q] != cost:
            return
        assert forward[i][p] + edge_cost == forward[j][q]
        edges.append(dict(from_query_boundary0=i, from_phase0=p,
            to_query_boundary0=j, to_phase0=q, operation=operation,
            template_base=template_base, query_base=query_base, cost=edge_cost))
        if operation == 'deletion':
            deletions[i].add(p)
            deletion_bases[i].add(template_base)
        else:
            query_edges[i].append((p, q, operation))
            query_ops[i].add(operation)
            if template_base is not None:
                query_bases[i].add(template_base)

    for i in range(n + 1):
        for p in range(k):
            q = (p + 1) % k
            add(i, p, i, q, 'deletion', 1, template[p], None)
            if i < n:
                operation = 'match' if sequence[i] == template[p] else 'substitution'
                add(i, p, i + 1, q, operation, int(operation != 'match'), template[p], sequence[i])
                add(i, p, i + 1, p, 'insertion', 1, None, sequence[i])
    assert all(query_ops)
    ranges = _operation_ranges(n, k, query_edges, deletions, starts, ends)
    boundary_ranges = _boundary_deletion_ranges(n, k, query_edges, deletions, starts, ends)
    public = dict(strand=strand, template=template, minimum_edit_cost=cost,
        globally_optimal=True, start_phases=sorted(starts), end_phases=sorted(ends),
        optimal_state_count=sum(forward[i][p] + backward[i][p] == cost
                                for i in range(n + 1) for p in range(k)),
        optimal_edge_count=len(edges), optimal_edges=edges, operation_count_ranges=ranges)
    return public, query_ops, query_bases, deletion_bases, boundary_ranges


def analyze_periodic_paths(sequence, motif, quality=None, start0=0):
    """Return all globally optimal cyclic-template operation constraints.

    Raises ValueError for malformed arguments/quality. Ambiguous sequence or
    motif, and fixed resource limits, return an explicit not_evaluated record.
    All coordinates are zero-based; query positions and deletion boundaries
    include relative and physical-read offsets. Empty query has score zero.
    """
    if not isinstance(sequence, str) or not isinstance(motif, str) or not motif:
        raise ValueError('sequence must be a string and motif a nonempty string')
    if type(start0) is not int or start0 < 0:
        raise ValueError('start0 must be a nonnegative integer')
    if quality is not None and (not isinstance(quality, str) or len(quality) != len(sequence)):
        raise ValueError('quality must be None or a string with the exact query length')
    n = len(sequence)
    result = dict(schema_version=SCHEMA_VERSION, method=METHOD,
        status='not_evaluated', evaluated=False, reason_codes=[],
        query_bp=n, query_start0=start0, query_end0=start0 + n,
        motif_input=motif, primitive_motif=None, primitive_period_bp=None,
        total_state_count=None,
        limits=dict(max_region_bp=MAX_REGION_BP, max_motif_bp=MAX_MOTIF_BP,
                    max_total_states=MAX_TOTAL_STATES),
        minimum_edit_cost=None, edit_cost_per_query_base=None,
        all_optimal_strands=[], operation_count_ranges=None, query_positions=[],
        deletion_boundaries=[], orientation_results=[],
        quality_available=quality is not None, quality_encoding_verified=False,
        quality_authenticity_verified=False,
        quality_interpretation='raw_ASCII_and_conditional_Phred33_no_calibration_claim',
        limitations=list(LIMITATIONS), native_counts_changed=False,
        classification_changed=False, biological_confirmation=None)
    if n > MAX_REGION_BP:
        result['reason_codes'].append('query_exceeds_region_bound')
    if len(motif) > MAX_MOTIF_BP:
        result['reason_codes'].append('motif_exceeds_length_bound')
    if result['reason_codes']:
        return result
    if quality is not None and any(not 33 <= ord(c) <= 126 for c in quality):
        raise ValueError('quality must contain printable ASCII FASTQ characters (33..126)')
    if set(sequence) - _DNA:
        result['reason_codes'].append('ambiguous_or_nonuppercase_query_DNA')
    if set(motif) - _DNA:
        result['reason_codes'].append('ambiguous_or_nonuppercase_motif_DNA')
    if result['reason_codes']:
        return result
    unit = _primitive(motif)
    k = len(unit)
    state_count = 2 * (n + 1) * k
    result.update(primitive_motif=unit, primitive_period_bp=k, total_state_count=state_count)
    if state_count > MAX_TOTAL_STATES:
        result['reason_codes'].append('query_motif_graph_exceeds_state_bound')
        return result
    templates = [('+', unit), ('-', unit.translate(_COMPLEMENT)[::-1])]
    forward = [_forward(sequence, template) for _, template in templates]
    costs = [min(rows[-1]) for rows in forward]
    best = min(costs)
    graphs = []
    orientations = []
    for (strand, template), rows, cost in zip(templates, forward, costs):
        if cost != best:
            orientations.append(dict(strand=strand, template=template,
                minimum_edit_cost=cost, globally_optimal=False,
                start_phases=[], end_phases=[], optimal_state_count=0,
                optimal_edge_count=0, optimal_edges=[], operation_count_ranges=None))
            continue
        graph = _optimal_graph(sequence, template, rows, _backward(sequence, template), best, strand)
        graphs.append(graph)
        orientations.append(graph[0])
    positions = []
    for i, base in enumerate(sequence):
        ops = set().union(*(g[1][i] for g in graphs))
        predicted = set().union(*(g[2][i] for g in graphs))
        character = quality[i] if quality is not None else None
        positions.append(dict(relative_position0=i, absolute_position0=start0 + i,
            query_base=base, operation_set=[op for op in QUERY_OPERATIONS if op in ops],
            template_base_set=sorted(predicted), mandatory_edit='match' not in ops,
            possible_edit=bool(ops - {'match'}), quality_character=character,
            quality_ascii=ord(character) if character is not None else None,
            conditional_phred33=ord(character) - 33 if character is not None else None))
    boundaries = []
    for i in range(n + 1):
        bases = set().union(*(g[3][i] for g in graphs))
        boundaries.append(dict(relative_boundary0=i, absolute_boundary0=start0 + i,
            template_base_set=sorted(bases),
            min_deletions=min(g[4][i]['min'] for g in graphs),
            max_deletions=max(g[4][i]['max'] for g in graphs), quality=None))
    ranges = {op: {'min': min(g[0]['operation_count_ranges'][op]['min'] for g in graphs),
                   'max': max(g[0]['operation_count_ranges'][op]['max'] for g in graphs)} for op in OPERATIONS}
    result.update(status='evaluated', evaluated=True, minimum_edit_cost=best,
        edit_cost_per_query_base=best / n if n else None,
        all_optimal_strands=[g[0]['strand'] for g in graphs],
        operation_count_ranges=ranges, query_positions=positions,
        deletion_boundaries=boundaries, orientation_results=orientations)
    return result
