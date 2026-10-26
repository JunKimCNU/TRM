#!/usr/bin/env python3
"""Slow independent exact-grammar oracle; deliberately separate from native code.

Intervals are 0-based half-open. A partial unit contributes span but never copies,
and is legal only at an actual record edge. All segmentations are enumerated with
counts capped at the qualifying count, so ambiguous parses cannot erase support.
"""
from __future__ import annotations


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _forward(sequence, units, min_span, min_copies):
    n = len(sequence)
    zero = (0,) * len(units)
    intervals = {}
    for start in range(n):
        initial = {(start, zero)}
        if start == 0:
            for unit in units:
                for partial in range(1, len(unit)):
                    if sequence.startswith(unit[-partial:]):
                        initial.add((partial, zero))
        pending = list(initial)
        seen = set(initial)
        while pending:
            end, counts = pending.pop()
            if end - start >= min_span:
                intervals.setdefault((start, end), set()).add(counts)
            remaining = n - end
            if remaining and n - start >= min_span:
                if any(remaining < len(unit) and sequence[end:] == unit[:remaining]
                       for unit in units):
                    intervals.setdefault((start, n), set()).add(counts)
            for index, unit in enumerate(units):
                if sequence.startswith(unit, end):
                    next_counts = list(counts)
                    next_counts[index] = min(min_copies, next_counts[index] + 1)
                    state = (end + len(unit), tuple(next_counts))
                    if state not in seen:
                        seen.add(state)
                        pending.append(state)
    # Retain all non-contained alternatives, independent of parsing tie order.
    maximal = []
    rightmost = -1
    for start, end in sorted(intervals, key=lambda item: (item[0], -item[1])):
        if end > rightmost:
            parses = intervals[start, end]
            minima = [min(counts[index] for counts in parses) for index in range(len(units))]
            maxima = [max(counts[index] for counts in parses) for index in range(len(units))]
            maximal.append({
                "start": start, "end": end,
                "guaranteed_mixed": sum(count >= min_copies for count in minima) >= 2,
                "possible_mixed": any(sum(count >= min_copies for count in counts) >= 2
                                      for counts in parses),
                "component_min_capped": minima, "component_max_capped": maxima,
                "count_cap": min_copies, "qualification_ambiguous": minima != maxima,
            })
            rightmost = end
    return maximal


def find_mixed(sequence, units, min_span=24, min_copies=2):
    sequence = sequence.upper()
    n = len(sequence)
    out = []
    for strand, oriented in (("+", sequence), ("-", reverse_complement(sequence))):
        for row in _forward(oriented, tuple(units), min_span, min_copies):
            start, end = row["start"], row["end"]
            if strand == "-":
                start, end = n - end, n - start
            out.append({**row, "start": start, "end": end, "strand": strand, "length": end - start})
    return sorted(out, key=lambda row: (row["start"], row["end"], row["strand"]))
