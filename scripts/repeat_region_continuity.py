"""Optional periodic-sequence continuity diagnostic; never a native TRM filter.

The score is a unit-cost edit distance to any finite segment of an infinite
periodic template. The query is consumed completely, the starting template
phase is free once, and subsequent phase changes require paid edits. The
distance divided by query length is neither native identity nor a biological
error-rate estimate. No alignment traceback or unique path is reported.
"""

from fractions import Fraction


VERSION = "0.1.0"
METHOD = "cyclic_unit_cost_edit_distance_complete_query"
DEFAULT_CONFIG = {
    "min_region_bp": 500,
    "max_region_bp": 100_000,
    "global_max_edit_ratio": 0.15,
    "window_bp": 100,
    "window_stride_bp": 25,
    "window_max_edit_ratio": 0.25,
    "gap_min_bp": 24,
    "gap_max_edit_ratio": 0.35,
    "gap_context_bp": 12,
    "gap_context_max_edit_ratio": 0.25,
}
_ACGT = frozenset("ACGT")


def _validate_motif(motif):
    if not isinstance(motif, str) or not motif or set(motif) - _ACGT:
        raise ValueError("motif must be a nonempty uppercase ACGT string")


def minimum_periodic_edit_score(sequence, motif):
    """Return the exact minimum cyclic-template edit score in O(n*k)/O(k).

    ``end_phases`` are all minimum-cost next-template-position states, numbered
    0 through len(motif)-1. They are not traceback paths or proof of alignment
    uniqueness. Empty query has cost zero, ratio None, and all end phases.

    A state p means motif[p] would be the next template base consumed. Initial
    costs are zero for every phase, allowing one arbitrary starting phase.
    Each row consumes one query base through insertion or substitution/match,
    then relaxes the forward template-deletion cycle at unit cost. Two sweeps
    propagate every possible minimum across the cycle; a cheapest deletion
    path never needs a complete positive-cost cycle.
    """
    _validate_motif(motif)
    if not isinstance(sequence, str) or set(sequence) - _ACGT:
        raise ValueError("sequence must contain only uppercase ACGT")
    period = len(motif)
    previous = [0] * period
    for base in sequence:
        current = [0] * period
        for phase in range(period):
            before = phase - 1
            insertion = previous[phase] + 1
            diagonal = previous[before] + (base != motif[before])
            current[phase] = min(insertion, diagonal)
        for edge in range(2 * period):
            phase = edge % period
            following = (phase + 1) % period
            deletion = current[phase] + 1
            if deletion < current[following]:
                current[following] = deletion
        previous = current
    cost = min(previous)
    return {
        "edit_cost": cost,
        "query_bp": len(sequence),
        "edit_cost_per_query_base": cost / len(sequence) if sequence else None,
        "end_phases": [phase for phase, value in enumerate(previous) if value == cost],
    }


def _positive_integer(value, label):
    if type(value) is not int or value <= 0:
        raise ValueError(label + " must be a positive integer")
    return value


def window_intervals(length, window_bp=100, stride_bp=25):
    """Return full sliding windows, appending the final window if necessary.

    A positive region shorter than one window is scored as one whole window.
    These are zero-based, half-open relative query coordinates. No sampled
    window is silently omitted. Empty query has no windows.
    """
    if type(length) is not int or length < 0:
        raise ValueError("length must be a nonnegative integer")
    _positive_integer(window_bp, "window_bp")
    _positive_integer(stride_bp, "stride_bp")
    if stride_bp > window_bp:
        raise ValueError("stride_bp must not leave uncovered sequence")
    if not length:
        return []
    if length <= window_bp:
        return [(0, length)]
    last_start = length - window_bp
    starts = list(range(0, last_start + 1, stride_bp))
    if starts[-1] != last_start:
        starts.append(last_start)
    return [(start, start + window_bp) for start in starts]


def _config(config):
    """Require versioned defaults, preventing accidental criterion changes."""
    result = dict(DEFAULT_CONFIG)
    if config is None:
        return result
    if not isinstance(config, dict):
        raise ValueError("config must be a dictionary or None")
    for key, value in config.items():
        if key not in result:
            raise ValueError("unknown config key: " + str(key))
        if type(value) is not type(result[key]) or value != result[key]:
            raise ValueError("version 0.1.0 uses fixed operational config: " + key)
    return result


def _gap_intervals(gaps, length):
    result = []
    previous_end = 0
    for interval in gaps:
        if not isinstance(interval, (tuple, list)) or len(interval) != 2:
            raise ValueError("each gap must be a pair of relative half-open coordinates")
        start, end = interval
        if type(start) is not int or type(end) is not int:
            raise ValueError("gap coordinates must be integers")
        if start < 0 or start > end or end > length:
            raise ValueError("gap coordinates are outside the region or reversed")
        if start == end:
            continue
        if start < previous_end:
            raise ValueError("gaps must be sorted and nonoverlapping")
        result.append((start, end))
        previous_end = end
    return result


def _passes(score, threshold):
    # Compare exact integer costs against the decimal operational threshold.
    limit = Fraction(str(threshold))
    return score["edit_cost"] * limit.denominator <= score["query_bp"] * limit.numerator


def _interval_score(sequence, motif, start, end):
    return {
        "start0": start,
        "end0": end,
        **minimum_periodic_edit_score(sequence[start:end], motif),
    }


def evaluate_region(sequence, motif, gaps=(), config=None):
    """Score a caller-selected complete region and every declared uncovered gap.

    ``gaps`` are sorted, nonoverlapping zero-based half-open coordinates relative
    to this sequence. Zero-length intervals are ignored. All nonempty gaps are
    scored; the standalone and context thresholds apply only to gaps >=24 bp.
    Context adds up to 12 bases on each side, clipped at the region boundary.

    Bounds and ACGT validation are explicit non-evaluation states, never scored
    negative evidence. The caller is responsible for tag/anchor coordinates,
    strand, primitive motif choice, provenance, and independent controls. This
    function never joins native tracts or alters the production detector.
    """
    _validate_motif(motif)
    if not isinstance(sequence, str):
        raise ValueError("sequence must be a string")
    settings = _config(config)
    intervals = _gap_intervals(gaps, len(sequence))
    result = {
        "schema_version": 1,
        "method": METHOD,
        "version": VERSION,
        "motif": motif,
        "query_bp": len(sequence),
        "coordinate_system": "zero-based half-open relative to supplied region",
        "config": settings,
        "evaluated": False,
        "status": None,
        "passes": None,
        "failed_reason_codes": [],
        "not_evaluated_reason_codes": [],
        "global": None,
        "windows": [],
        "gaps": [],
        "declared_nonempty_gap_count": len(intervals),
        "max_window_edit_ratio": None,
        "max_gap_edit_ratio": None,
        "max_gap_context_edit_ratio": None,
        "max_guarded_gap_edit_ratio": None,
        "max_guarded_gap_context_edit_ratio": None,
        "invalid_sequence_symbols": sorted(set(sequence) - _ACGT),
    }
    if result["invalid_sequence_symbols"]:
        result["not_evaluated_reason_codes"].append("ambiguous_or_non_acgt_sequence")
    if len(sequence) < settings["min_region_bp"]:
        result["not_evaluated_reason_codes"].append("region_below_minimum_scope")
    if len(sequence) > settings["max_region_bp"]:
        result["not_evaluated_reason_codes"].append("region_exceeds_length_limit")
    if result["not_evaluated_reason_codes"]:
        if result["invalid_sequence_symbols"]:
            result["status"] = "not_evaluated_ambiguous"
        elif len(sequence) < settings["min_region_bp"]:
            result["status"] = "outside_length_scope_too_short"
        else:
            result["status"] = "not_evaluated_length_limit"
        return result

    result["evaluated"] = True
    whole = _interval_score(sequence, motif, 0, len(sequence))
    whole["passes"] = _passes(whole, settings["global_max_edit_ratio"])
    result["global"] = whole
    if not whole["passes"]:
        result["failed_reason_codes"].append("global_edit_ratio_exceeded")

    for start, end in window_intervals(len(sequence), settings["window_bp"], settings["window_stride_bp"]):
        row = _interval_score(sequence, motif, start, end)
        row["passes"] = _passes(row, settings["window_max_edit_ratio"])
        result["windows"].append(row)
    result["max_window_edit_ratio"] = max(row["edit_cost_per_query_base"] for row in result["windows"])
    if any(not row["passes"] for row in result["windows"]):
        result["failed_reason_codes"].append("window_edit_ratio_exceeded")

    for start, end in intervals:
        row = _interval_score(sequence, motif, start, end)
        row["guard_applied"] = end - start >= settings["gap_min_bp"]
        row["passes"] = _passes(row, settings["gap_max_edit_ratio"]) if row["guard_applied"] else None
        context = settings["gap_context_bp"]
        context_row = _interval_score(sequence, motif, max(0, start - context), min(len(sequence), end + context))
        context_row["passes"] = _passes(context_row, settings["gap_context_max_edit_ratio"]) if row["guard_applied"] else None
        row["context"] = context_row
        result["gaps"].append(row)
    for key, rows, context in (
        ("max_gap_edit_ratio", result["gaps"], False),
        ("max_gap_context_edit_ratio", result["gaps"], True),
        ("max_guarded_gap_edit_ratio", [row for row in result["gaps"] if row["guard_applied"]], False),
        ("max_guarded_gap_context_edit_ratio", [row for row in result["gaps"] if row["guard_applied"]], True),
    ):
        if rows:
            result[key] = max((row["context"] if context else row)["edit_cost_per_query_base"] for row in rows)
    if any(row["passes"] is False for row in result["gaps"]):
        result["failed_reason_codes"].append("gap_edit_ratio_exceeded")
    if any(row["context"]["passes"] is False for row in result["gaps"]):
        result["failed_reason_codes"].append("gap_context_edit_ratio_exceeded")
    result["passes"] = not result["failed_reason_codes"]
    result["status"] = "passes_operational_continuity" if result["passes"] else "fails_operational_continuity"
    return result
