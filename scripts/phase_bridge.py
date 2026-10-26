"""Optional, bounded sequence bridges between registered repeat tracts.

This module evaluates a local hypothesis. It neither discovers nor merges
tracts, changes a native/policy classification, nor establishes a chromosome
end. Exact, adjacent flank anchors fix both template phases. Every observed gap
base must be consumed by a global unit-edit alignment; there are no free ends
or free phase changes. An insertion consumes query only, a deletion template
only, and a mismatch both. The edit budget is a declared heuristic, not an
error probability or evidence of the biological cause of a phase difference.

Only this small optional postprocessor imports tract_context; the original
fixed-phase model and production scan are unchanged.
"""

import math

try:
    from . import tract_context as context
except ImportError:  # Direct scripts-directory import / execution.
    import tract_context as context


VERSION = "0.1.0"
SCHEMA_VERSION = 1
DEFAULT_CONFIG = {
    "max_gap_bp": 1024,
    "max_span_bp": 10_000,
    "max_edits": 16,
    "max_edit_fraction": 0.05,
    "min_edits_floor": 1,
    "max_indel_run": 3,
    "min_anchor_copies": 2,
    "max_dp_states": 250_000,
}
MAX_EDIT_BUDGET = 256
MAX_DP_STATES = 2_000_000
SUCCESS = "compatible_under_bounded_phase_bridge"
REJECTED = "rejected_by_bounded_phase_bridge"
COMPUTATION_LIMIT = "not_evaluated_computation_limit"


def validate_config(config=None):
    """Return complete validated settings without mutating caller input.

The effective edit budget is min(max_edits, max(min_edits_floor,
floor(gap_length * max_edit_fraction))). The default floor of one deliberately
allows a one-base edit at a short or empty junction. Set it to zero to remove
that allowance. It is not a claim that all short gaps have <=5% edit rate.
"""
    if config is not None and not isinstance(config, dict):
        raise ValueError("Phase-bridge config must be a dictionary")
    config = {} if config is None else config
    unknown = set(config) - set(DEFAULT_CONFIG)
    if unknown:
        raise ValueError(f"Unknown phase-bridge config keys: {sorted(map(str, unknown))}")
    values = dict(DEFAULT_CONFIG, **config)
    context._integer(values["max_gap_bp"], "max_gap_bp", maximum=context.MAX_ANALYZED_INTERVAL_BP)
    context._integer(values["max_span_bp"], "max_span_bp", minimum=1, maximum=context.MAX_ANALYZED_INTERVAL_BP)
    context._integer(values["max_edits"], "max_edits", maximum=MAX_EDIT_BUDGET)
    context._integer(values["min_edits_floor"], "min_edits_floor", maximum=MAX_EDIT_BUDGET)
    context._integer(values["max_indel_run"], "max_indel_run", maximum=MAX_EDIT_BUDGET)
    context._integer(values["min_anchor_copies"], "min_anchor_copies", minimum=1, maximum=1000)
    context._integer(values["max_dp_states"], "max_dp_states", minimum=1, maximum=MAX_DP_STATES)
    fraction = values["max_edit_fraction"]
    if isinstance(fraction, bool) or not isinstance(fraction, (int, float)) or not math.isfinite(fraction) or not 0 <= fraction <= 1:
        raise ValueError("max_edit_fraction must be finite and between 0 and 1")
    return values


def _empty_alignment(sequence, motif, start_phase, end_phase, max_edits,
                     max_indel_run, max_dp_states):
    return dict(
        schema_version=SCHEMA_VERSION, module_version=VERSION,
        algorithm="global_unit_edit_cyclic_template_fixed_boundary_phases",
        motif=motif, motif_period_bp=len(motif), query_length_bp=len(sequence),
        start_phase=start_phase, end_phase=end_phase,
        phase_definition="Index of the next expected template base at each boundary",
        required_net_template_delta_mod_period=(end_phase - start_phase - len(sequence)) % len(motif),
        edit_budget=max_edits, max_indel_run=max_indel_run,
        max_dp_states=max_dp_states, dp_states=0,
        min_edits=None, optimal_path_count_capped2=0,
        optimal_path_ambiguous=None, path_ops=[], cigar=None,
        query_consumed_bp=None, template_consumed_bp=None,
        operation_counts=None, matched_query_intervals=[],
        mismatched_query_intervals=[], insertion_query_intervals=[],
        template_deletions=[], edit_runs=[],
        query_coordinate_system="zero-based half-open within supplied gap string",
        template_coordinate_system="unrolled repeat beginning at start_phase",
        all_query_bases_consumed=False, free_end_clipping_allowed=False,
        free_phase_changes_allowed=False, biological_cause_inferred=False,
        biological_tracts_merged=False, physical_chromosome_end_confirmed=False,
        changes_native_classification=False, changes_policy_classification=False,
    )


def align_bridge(sequence, motif, start_phase, end_phase, max_edits,
                 max_indel_run, max_dp_states=250_000):
    """Globally align a gap to a cyclic template with fixed boundary phases.

All phases index the *next* template base. Return the minimum edit cost among
paths meeting both the budget and the maximum consecutive I or D run. The run
limit is part of the DP state, so a forbidden unconstrained optimum does not
hide a valid constrained path. Consecutive I and D runs are bounded separately;
switching between I and D starts a new run and still pays for every base.

States are (query consumed, template consumed, last indel type, run length).
Template consumption is bounded by query length + edit budget. Lexicographic
query/template order is topological, including zero-cost matches. Each state
retains its minimum cost, deterministic first predecessor, and a count of all
optimal paths capped at two. All minimum-cost terminal states contribute to
ambiguity. No partially explored answer is returned if the state cap is hit.

Non-ACGT letters are reported as a rejected ambiguous gap, not aligned as free
mismatches. Invalid coordinates/types/settings raise ValueError. The public
kernel is usable without registered tracts for independent algorithm tests.
"""
    motif = context._motif(motif)
    sequence = context._chunk(sequence, 0, len(sequence) if isinstance(sequence, str) else 0)
    k, n = len(motif), len(sequence)
    context._integer(start_phase, "start_phase", maximum=k - 1)
    context._integer(end_phase, "end_phase", maximum=k - 1)
    context._integer(max_edits, "max_edits", maximum=MAX_EDIT_BUDGET)
    context._integer(max_indel_run, "max_indel_run", maximum=MAX_EDIT_BUDGET)
    context._integer(max_dp_states, "max_dp_states", minimum=1, maximum=MAX_DP_STATES)
    result = _empty_alignment(sequence, motif, start_phase, end_phase,
                              max_edits, max_indel_run, max_dp_states)
    if any(base not in "ACGT" for base in sequence):
        result.update(status=REJECTED, reason_codes=["ambiguous_gap_bases"], evaluated=False)
        return result

    # Each node is [cost, capped_path_count, predecessor_state, operation].
    # The run state is ('',0) after M/X and (I/D,length) after an indel.
    initial = (0, 0, "", 0)
    nodes = {initial: [0, 1, None, None]}
    current = {0: {("", 0): initial}}
    work_limit_hit = False

    def relax(state, cost, predecessor, operation, row):
        nonlocal work_limit_hit
        if cost > max_edits:
            return
        existing = nodes.get(state)
        if existing is None:
            if len(nodes) >= max_dp_states:
                work_limit_hit = True
                return
            nodes[state] = [cost, nodes[predecessor][1], predecessor, operation]
            row.setdefault(state[1], {})[(state[2], state[3])] = state
        elif cost < existing[0]:
            nodes[state] = [cost, nodes[predecessor][1], predecessor, operation]
        elif cost == existing[0]:
            existing[1] = min(2, existing[1] + nodes[predecessor][1])

    for i in range(n + 1):
        following = {}
        for j in range(max(0, i - max_edits), min(n + max_edits, i + max_edits) + 1):
            # Same-cell states never depend on one another. Sorting makes
            # tie-breaking independent of dictionary insertion details.
            for _, state in sorted(current.get(j, {}).items()):
                cost, _, _, _ = nodes[state]
                kind, run = state[2:]
                if i < n:
                    op = "M" if sequence[i] == motif[(start_phase + j) % k] else "X"
                    relax((i + 1, j + 1, "", 0), cost + (op == "X"), state, op, following)
                    next_run = run + 1 if kind == "I" else 1
                    if next_run <= max_indel_run:
                        relax((i + 1, j, "I", next_run), cost + 1, state, "I", following)
                next_run = run + 1 if kind == "D" else 1
                if j < n + max_edits and next_run <= max_indel_run:
                    relax((i, j + 1, "D", next_run), cost + 1, state, "D", current)
                if work_limit_hit:
                    result.update(status=COMPUTATION_LIMIT,
                                  reason_codes=["dp_state_cap_reached"],
                                  evaluated=False, dp_states=len(nodes))
                    return result
        if i < n:
            current = following

    result["dp_states"] = len(nodes)
    terminals = [state for j, states in current.items()
                 if (start_phase + j) % k == end_phase
                 for state in states.values()]
    if not terminals:
        result.update(status=REJECTED, reason_codes=["no_path_within_edit_budget_and_indel_run"], evaluated=True)
        return result
    minimum = min(nodes[state][0] for state in terminals)
    optimal = sorted(state for state in terminals if nodes[state][0] == minimum)
    path_count = min(2, sum(nodes[state][1] for state in optimal))
    selected = optimal[0]
    steps = []
    state = selected
    while nodes[state][2] is not None:
        predecessor, operation = nodes[state][2:]
        steps.append((operation, predecessor[0], state[0], predecessor[1], state[1]))
        state = predecessor
    steps.reverse()
    runs = []
    for operation, q0, q1, t0, t1 in steps:
        if runs and runs[-1]["operation"] == operation:
            runs[-1]["query_end0"] = q1
            runs[-1]["template_end0"] = t1
            runs[-1]["template_end_phase"] = (start_phase + t1) % k
            runs[-1]["length_bp"] += 1
        else:
            runs.append(dict(operation=operation, length_bp=1,
                             query_start0=q0, query_end0=q1,
                             template_start0=t0, template_end0=t1,
                             template_start_phase=(start_phase + t0) % k,
                             template_end_phase=(start_phase + t1) % k))
    counts = {op: sum(run["length_bp"] for run in runs if run["operation"] == op) for op in "MXID"}
    for run in runs:
        if run["operation"] != "M":
            run["query_bases"] = sequence[run["query_start0"]:run["query_end0"]]
            run["template_bases"] = "".join(motif[(start_phase + j) % k]
                                            for j in range(run["template_start0"], run["template_end0"]))
    result.update(
        status=SUCCESS, reason_codes=[], evaluated=True, min_edits=minimum,
        optimal_path_count_capped2=path_count, optimal_path_ambiguous=path_count > 1,
        deterministic_trace_policy="First optimal predecessor in lexicographic query/template/run traversal; smallest terminal state",
        path_ops=runs, cigar="".join(f"{run['length_bp']}{run['operation']}" for run in runs),
        query_consumed_bp=selected[0], template_consumed_bp=selected[1],
        operation_counts=counts, all_query_bases_consumed=selected[0] == n,
        matched_query_intervals=[[r["query_start0"], r["query_end0"]] for r in runs if r["operation"] == "M"],
        mismatched_query_intervals=[[r["query_start0"], r["query_end0"]] for r in runs if r["operation"] == "X"],
        insertion_query_intervals=[[r["query_start0"], r["query_end0"]] for r in runs if r["operation"] == "I"],
        template_deletions=[dict(r) for r in runs if r["operation"] == "D"],
        edit_runs=[dict(r) for r in runs if r["operation"] != "M"],
        max_observed_indel_run_bp=max((r["length_bp"] for r in runs if r["operation"] in "ID"), default=0),
        edit_rate_denominator_bp=max(n, 1), edit_cost_per_observed_base_floor1=minimum / max(n, 1),
        edit_rate_is_not_sequence_identity=True,
    )
    return result


def _exact_boundary_anchor(sequence, tract, motif, copies, side):
    length = copies * len(motif)
    if side == "left":
        start0, end0 = max(tract["start0"], tract["end0"] - length), tract["end0"]
    else:
        start0, end0 = tract["start0"], min(tract["end0"], tract["start0"] + length)
    chunk = context._chunk(sequence, start0, end0)
    template = motif if tract["strand"] == "+" else context.reverse_complement(motif)
    exact_phases = [phase for phase in range(len(motif))
                    if all(base == template[(i + phase) % len(motif)] for i, base in enumerate(chunk))]
    enough = len(chunk) >= length
    selected = exact_phases[0] if enough and len(exact_phases) == 1 else None
    return dict(
        tract_id=tract["tract_id"], side=side, declared_strand=tract["strand"],
        start0=start0, end0=end0, observed_bp=len(chunk), required_bp=length,
        required_exact_copies=copies, sufficient_length=enough,
        exact_start_phases=exact_phases, unique_exact_phase=selected is not None,
        selected_start_phase=selected,
        selected_absolute_phase=(selected - start0) % len(motif) if selected is not None else None,
        next_phase_at_gap_boundary=(selected + len(chunk)) % len(motif) if selected is not None and side == "left" else selected,
        exact_interval_sequence=chunk,
    )


def _offset_alignment_query_coordinates(alignment, offset):
    """Translate only query coordinates; cyclic template coordinates stay local."""
    for key in ("path_ops", "template_deletions", "edit_runs"):
        for row in alignment[key]:
            row["query_start0"] += offset
            row["query_end0"] += offset
    for key in ("matched_query_intervals", "mismatched_query_intervals", "insertion_query_intervals"):
        alignment[key] = [[a + offset, b + offset] for a, b in alignment[key]]
    alignment["query_coordinate_system"] = "zero-based half-open within supplied full record"
    alignment["gap_start0"] = offset
    alignment["gap_end0"] = offset + alignment["query_length_bp"]


def evaluate_pair(sequence, left, right, motif, config=None, alignments=None, expected_strand=None):
    """Evaluate two supplied tracts using exact flank anchors and a gap bridge.

Registered left/right tracts must be ordered, nonoverlapping and same-strand.
Only the suffix of the left tract and prefix of the right tract provide exact
anchors; their default lengths are two primitive copies each. Larger tract
interiors are not assumed exact and are not reanalyzed. The bridge is entirely
between their boundaries. Reverse-complement rotation-equivalent motifs and
non-ACGT gap bases are rejected. Unknown expected strand is disclosed, while
conflicting supplied expected strand is rejected. Alignment observations stay
descriptive and independent of sequence fit.
"""
    motif = context._motif(motif)
    values = validate_config(config)
    context._strand(expected_strand)
    left, right = context._tract(sequence, left, "left"), context._tract(sequence, right, "right")
    start0, end0 = min(left["start0"], right["start0"]), max(left["end0"], right["end0"])
    gap_bp = right["start0"] - left["end0"]
    direction_ambiguous = context.reverse_complement(motif) in motif + motif
    reasons = []
    if left["start0"] > right["start0"]:
        reasons.append("anchors_out_of_order")
    if max(left["start0"], right["start0"]) < min(left["end0"], right["end0"]):
        reasons.append("anchors_overlap")
    if left["strand"] != right["strand"]:
        reasons.append("anchor_orientation_conflict")
    if direction_ambiguous:
        reasons.append("anchor_template_strand_ambiguous")
    if expected_strand is not None and any(t["strand"] != expected_strand for t in (left, right)):
        reasons.append("expected_strand_conflict")
    if gap_bp > values["max_gap_bp"]:
        reasons.append("gap_exceeds_max_bp")
    if end0 - start0 > values["max_span_bp"]:
        reasons.append("span_exceeds_max_bp")
    terminal_status = ("unassessed_expected_strand_unknown" if expected_strand is None
                       else "direction_ambiguous" if direction_ambiguous
                       else "expected_strand_conflict" if "expected_strand_conflict" in reasons
                       else "expected_strand_consistent")
    result = dict(
        schema_version=SCHEMA_VERSION, module_version=VERSION,
        model="bounded_exact_flank_anchored_cyclic_edit_bridge",
        motif=motif, motif_period_bp=len(motif), config=values,
        left=dict(left), right=dict(right), start0=start0, end0=end0,
        length_bp=end0 - start0, observed_gap_bp=gap_bp,
        record_context=context._record_context(sequence, start0, end0),
        reverse_complement_rotation_equivalent=direction_ambiguous,
        terminal_context=dict(status=terminal_status, expected_strand=expected_strand,
                              record_edge_is_not_terminal_evidence=True, physical_chromosome_end_confirmed=False),
        alignment_evidence=context._alignment_evidence(sequence, alignments, start0, end0),
        anchor_phase_evidence=None, phase_compatibility=None, bridge=None,
        effective_edit_budget=None, diagnostics_skipped=bool(reasons),
        biological_tracts_merged=False, changes_native_classification=False,
        changes_policy_classification=False, physical_chromosome_end_confirmed=False,
        continuity_inference_permitted=False, biological_cause_inferred=False,
        limitations=[
            "A compatible bridge is an optional local sequence hypothesis, not a merged biological tract or physical chromosome end.",
            "Exact flank anchors reject errors touching the anchor windows; no boundary search or anchor rescue is performed.",
            "A one-edit budget floor permits short/empty junction edits and is not a 5% short-gap error-rate claim.",
            "The unit-edit model does not distinguish biological variation, sequencing error, or assembly error.",
            "Multiple optimum edit placements are counted up to two; the selected trace is illustrative and deterministic.",
            "Supplied mapping spans are descriptive, not a CIGAR continuity, chimera, or structural-variant assessment.",
            "No calibrated sensitivity, specificity, physical-end inference, or absence call is implied.",
        ],
    )
    if reasons:
        result.update(status=REJECTED, reason_codes=reasons)
        return result
    evidence = {side: _exact_boundary_anchor(sequence, tract, motif, values["min_anchor_copies"], side)
                for side, tract in (("left", left), ("right", right))}
    result["anchor_phase_evidence"] = evidence
    for side, item in evidence.items():
        if not item["sufficient_length"]:
            reasons.append(f"{side}_anchor_too_short")
        elif not item["exact_start_phases"]:
            reasons.append(f"{side}_anchor_not_exact")
        elif not item["unique_exact_phase"]:
            reasons.append(f"{side}_anchor_phase_ambiguous")
    gap = context._chunk(sequence, left["end0"], right["start0"])
    if any(base not in "ACGT" for base in gap):
        reasons.append("ambiguous_gap_bases")
    if reasons:
        result.update(status=REJECTED, reason_codes=reasons, diagnostics_skipped=True)
        return result
    budget = min(values["max_edits"], max(values["min_edits_floor"], math.floor(gap_bp * values["max_edit_fraction"])))
    result["effective_edit_budget"] = budget
    result["budget_definition"] = "min(max_edits,max(min_edits_floor,floor(observed_gap_bp*max_edit_fraction)))"
    start_phase, end_phase = [evidence[side]["next_phase_at_gap_boundary"] for side in ("left", "right")]
    result["phase_compatibility"] = dict(
        strand=left["strand"], start_phase=start_phase, end_phase=end_phase,
        left_absolute_phase=evidence["left"]["selected_absolute_phase"],
        right_absolute_phase=evidence["right"]["selected_absolute_phase"],
        phase_delta_mod_period=(evidence["right"]["selected_absolute_phase"] - evidence["left"]["selected_absolute_phase"]) % len(motif),
        required_net_template_delta_mod_period=(end_phase - start_phase - gap_bp) % len(motif),
        zero_edit_phase_consistent=(start_phase + gap_bp) % len(motif) == end_phase,
        empty_gap=gap_bp == 0,
    )
    template = motif if left["strand"] == "+" else context.reverse_complement(motif)
    bridge = align_bridge(gap, template, start_phase, end_phase, budget,
                          values["max_indel_run"], values["max_dp_states"])
    _offset_alignment_query_coordinates(bridge, left["end0"])
    result.update(bridge=bridge, status=bridge["status"], reason_codes=bridge["reason_codes"],
                  diagnostics_skipped=bridge["status"] == COMPUTATION_LIMIT)
    return result
