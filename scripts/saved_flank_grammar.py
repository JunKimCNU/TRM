"""Exact block-language factors in bounded saved representative sequences.

Unlike complete-token assembly paths, query edges can cut a block at any phase.
Literal block phases are preserved between units. This is descriptive evidence,
not a telomere, nonrepeat, read continuity or genome-wide uniqueness classifier.
"""
from functools import lru_cache
from telomere_grammar import validate_units, natural

VERSION = '0.1.0'
MAX_SEQUENCE_BP = 80
MAX_TOTAL_UNIT_BP = 512
RC = str.maketrans('ACGT', 'TGCA')


def checked_units(units):
    validate_units(units)
    if sum(map(len, units)) > MAX_TOTAL_UNIT_BP:
        raise ValueError('sum of literal block lengths must be <=512')
    return tuple(units)


@lru_cache(maxsize=16)
def _machine(units):
    labels, predecessors, ends = [], [], []
    for unit in units:
        start = len(labels)
        labels.extend(unit)
        predecessors.extend([None] + list(range(start, start + len(unit) - 1)))
        ends.append(len(labels) - 1)
    return labels, predecessors, ends


def _factors(sequence, units):
    """Earliest start for each automaton state: O(n * total unit bases)."""
    labels, predecessors, ends = _machine(units)
    inactive = len(sequence) + 1
    previous = [inactive] * len(labels)
    endpoints = []
    for pos, base in enumerate(sequence):
        boundary = min(previous[i] for i in ends)
        current = [inactive] * len(labels)
        for state, label in enumerate(labels):
            if label == base:
                parent = boundary if predecessors[state] is None else previous[predecessors[state]]
                current[state] = min(pos, parent)
        earliest = min(current)
        if earliest <= pos:
            endpoints.append((earliest, pos + 1))
        previous = current
    maximal, farthest = [], -1
    for start, end in sorted(set(endpoints), key=lambda p: (p[0], -p[1])):
        if end > farthest:
            maximal.append((start, end))
            farthest = end
    return maximal


def analyze(sequence, units, min_span=20):
    """Annotate one 1..80bp uppercase ACGT query, with independent cutoff union.

    Full-query compatibility has no min_span gate. Outside the selected union
    means outside this exact model only. No result is an event-weighted estimate.
    """
    units = checked_units(units)
    natural(min_span, 'minimum union span', 1)
    if min_span > MAX_SEQUENCE_BP:
        raise ValueError('minimum union span must be <=80')
    if (not isinstance(sequence, str) or not 1 <= len(sequence) <= MAX_SEQUENCE_BP
            or set(sequence) - set('ACGT')):
        raise ValueError('saved sequence must contain 1..80 uppercase ACGT bases')
    n = len(sequence)
    views = [('+', sequence), ('-', sequence.translate(RC)[::-1])]
    all_intervals, full_strands = [], []
    for strand, view in views:
        factors = _factors(view, units)
        if (0, n) in factors:
            full_strands.append(strand)
        for a, b in factors:
            x, y = (a, b) if strand == '+' else (n-b, n-a)
            all_intervals.append((x, y, strand))
    pure = [unit for unit in units if any(
        all(base == unit[(pos+phase) % len(unit)] for pos, base in enumerate(view))
        for _, view in views for phase in range(len(unit)))]
    covered = [False] * n
    intervals = sorted((a, b, s) for a, b, s in all_intervals if b-a >= min_span)
    for a, b, _ in intervals:
        covered[a:b] = [True] * (b-a)
    prefix = next((i for i, v in enumerate(covered) if v), n)
    suffix = next((i for i, v in enumerate(reversed(covered)) if v), n)
    longest = current = 0
    for value in covered:
        current = 0 if value else current + 1
        longest = max(longest, current)
    return dict(sequence_bases=n, full_grammar_compatible=bool(full_strands),
        full_single_unit_compatible=bool(pure), full_compatible_units=pure,
        full_grammar_strands=full_strands,
        compatibility_class='single_unit_compatible' if pure else 'mixed_only' if full_strands else 'not_full_compatible',
        longest_grammar_span_bp=max((b-a for a, b, _ in all_intervals), default=0),
        minimum_union_span_bp=min_span, grammar_union_bp=sum(covered),
        prefix_outside_bp=prefix, suffix_outside_bp=suffix, max_contiguous_outside_bp=longest,
        intervals=[dict(start0=a, end0=b, strand=s) for a, b, s in intervals])
