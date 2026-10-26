"""Additive repeat-unit specificity diagnostic, never a TRM decision/filter.

The supplied complete interval is compared with one supplied primitive unit.
Fixed-phase and globally consuming cyclic-edit evidence are kept separate.
The operational labels below are neither calibrated error probabilities nor
biological repeat/absence/novelty/terminal-location classifications.
"""
from collections import Counter
import hashlib

from repeat_region_continuity import minimum_periodic_edit_score
import tract_context

VERSION = "0.1.0"
SCHEMA_VERSION = 1
METHOD = "supplied_unit_vs_composition_and_homopolymer_complete_interval"
DEFAULT_CONFIG = {
    "max_region_bp": 10_000,
    "max_dp_cells": 5_000_000,
    "minimum_complete_unit_copies": 3,
    "minimum_discriminating_matches": 3,
    "minimum_fixed_identity": 0.90,
    "minimum_normalized_fixed_gain": 0.50,
    "minimum_normalized_edit_gain": 0.50,
    "minimum_discriminating_precision": 0.80,
    "minimum_discriminating_recall": 0.80,
}
# Public registration snapshot for callers; mutating it does not change policy.
DEFAULTS = dict(DEFAULT_CONFIG)
MAX_MOTIF_BP = 256
HARD_MAX_REGION_BP = 100_000
HARD_MAX_DP_CELLS = 50_000_000
BASES = "ACGT"


def validate_config(config=None):
    """Only resource bounds are adjustable in this version, before analysis."""
    result = dict(DEFAULT_CONFIG)
    if config is None:
        return result
    if not isinstance(config, dict):
        raise ValueError("config must be a dictionary or None")
    for key, value in config.items():
        if key not in result:
            raise ValueError("unknown config key: " + str(key))
        if key in ("max_region_bp", "max_dp_cells"):
            ceiling = HARD_MAX_REGION_BP if key == "max_region_bp" else HARD_MAX_DP_CELLS
            if type(value) is not int or not 1 <= value <= ceiling:
                raise ValueError(key + " must be a positive integer <= " + str(ceiling))
            result[key] = value
        elif type(value) is not type(result[key]) or value != result[key]:
            raise ValueError("version0.1.0 operational thresholds are fixed: " + key)
    return result


def primitive_unit(motif):
    for size in range(1, len(motif) + 1):
        if len(motif) % size == 0 and motif == motif[:size] * (len(motif) // size):
            return motif[:size]
    raise AssertionError("nonempty motif required")


def _residue_fixed_phase_scores(sequence, start0, motif):
    """Equivalent to tract_context fixed phases, extended to periods65..256.

    Residue/base counts make this O(n+k^2), without a query-sized template for
    every phase. The shared public helper is used for periods<=64; this version
    is independently tested against that helper on its full supported range.
    """
    k = len(motif)
    residue = [Counter() for _ in range(k)]
    for position, base in enumerate(sequence):
        residue[position % k][base] += 1
    scores = []
    for strand, template in [("+", motif), ("-", tract_context.reverse_complement(motif))]:
        for rotation in range(k):
            matches = sum(residue[r][template[(r + rotation) % k]] for r in range(k))
            scores.append({"strand": strand, "rotation_at_interval_start": rotation,
                           "absolute_template_phase": (rotation - start0) % k,
                           "motif_at_interval_start": template[rotation:] + template[:rotation],
                           "matches_bp": matches, "mismatches_bp": len(sequence) - matches,
                           "denominator_bp": len(sequence), "identity": matches / len(sequence) if sequence else None})
    best = max(x["matches_bp"] for x in scores)
    return {"model": "One fixed phase across the entire interval; no indels or phase switching",
            "scores": scores, "all_best_ties": [x for x in scores if x["matches_bp"] == best],
            "best_by_strand": {s: [x for x in scores if x["strand"] == s and x["matches_bp"] == max(y["matches_bp"] for y in scores if y["strand"] == s)] for s in ("+", "-")},
            "empty_interval": not sequence, "continuity_inference_permitted": False}


def _fixed(sequence, start0, motif):
    if len(motif) <= 64:
        return tract_context.fixed_phase_scores(sequence, start0, motif)
    return _residue_fixed_phase_scores(sequence, start0, motif)


def _finite_template_counts(unit, length):
    cycles, remainder = divmod(length, len(unit))
    counts = Counter({base: unit.count(base) * cycles for base in BASES})
    counts.update(unit[:remainder])
    return counts


def _discriminating(sequence, phase_unit, baseline):
    template_counts = _finite_template_counts(phase_unit, len(sequence))
    expected_positive = len(sequence) - template_counts[baseline]
    observed_positive = len(sequence) - sequence.count(baseline)
    tp = sum(base == phase_unit[i % len(phase_unit)] and base != baseline for i, base in enumerate(sequence))
    return {"majority_baseline_base": baseline, "template_nonmajority_positions": expected_positive,
            "query_nonmajority_positions": observed_positive, "correct_nonmajority_base_matches": tp,
            "precision": tp / expected_positive if expected_positive else None,
            "recall": tp / observed_positive if observed_positive else None,
            "precision_definition": "exact nonbaseline base matches / template nonbaseline positions",
            "recall_definition": "exact nonbaseline base matches / query nonbaseline positions",
            "zero_denominator_is_not_positive_evidence": True}


def evaluate_unit_specificity(sequence, motif, *, start0=0, config=None):
    """Return full metrics and an operational descriptive status.

    All strand/phase ties and all tied majority bases are retained. Precision
    and recall apply ONLY to the fixed-phase model, not an inferred edit path.
    The cyclic model has one free initial phase; all query bases are consumed
    and subsequent phase changes require paid edits. No traceback is claimed.

    ``start0`` records the original full-record interval start. It changes only
    coordinate metadata, not the fit. Lowercase is explicitly uppercased; N or
    other ambiguous symbols make the interval not_evaluated. No masking,
    prefix selection, caller-sequence mutation, or native reclassification.
    """
    if not isinstance(sequence, str) or not isinstance(motif, str):
        raise ValueError("sequence and motif must be strings")
    if type(start0) is not int or start0 < 0:
        raise ValueError("start0 must be a nonnegative integer")
    settings = validate_config(config)
    seq, supplied = sequence.upper(), motif.upper()
    valid_motif = bool(supplied) and set(supplied) <= set(BASES) and len(supplied) <= MAX_MOTIF_BP
    unit = primitive_unit(supplied) if valid_motif else None
    n = len(seq)
    result = {
        "schema_version": SCHEMA_VERSION, "version": VERSION, "method": METHOD,
        "status": "not_evaluated", "evaluated": False, "reason_codes": [],
        "supplied_motif": motif, "normalized_supplied_motif": supplied,
        "primitive_motif": unit, "primitive_period_bp": len(unit) if unit else None,
        "supplied_motif_reduced": unit != supplied if unit else None,
        "query_bp": n, "start0": start0, "end0": start0 + n,
        "coordinate_system": "original forward record, zero-based half-open",
        "sequence_normalization": "uppercase only; ambiguous bases are not masked or discarded",
        "normalized_sequence_sha256": hashlib.sha256(seq.encode()).hexdigest(),
        "config": settings, "operational_thresholds_calibrated": False,
        "fixed_phase": None, "majority_baseline": None, "cyclic_edit": None,
        "model_comparison": None, "changes_native_classification": False,
        "biological_repeat_unit_confirmed": False, "physical_terminal_confirmed": False,
        "biological_absence_or_novelty_inference_permitted": False,
        "scope": "Specificity to the supplied unit only; a mixed or degenerate biological array may lack specificity to this one literal unit.",
    }
    if not valid_motif:
        result["reason_codes"].append("motif_outside_nonempty_acgt_length1_to256_scope")
    if not seq:
        result["reason_codes"].append("empty_interval")
    if set(seq) - set(BASES):
        result["reason_codes"].append("ambiguous_or_non_acgt_sequence")
    if n > settings["max_region_bp"]:
        result["reason_codes"].append("interval_exceeds_registered_length_bound")
    dp_cells = 2 * n * len(unit) if unit else None
    result["estimated_cyclic_dp_cells_both_strands"] = dp_cells
    result["dp_work_definition"] = "2 * query_bp * primitive_period; each cell entails constant relaxation work"
    if dp_cells is not None and dp_cells > settings["max_dp_cells"]:
        result["reason_codes"].append("interval_exceeds_registered_dp_work_bound")
    if result["reason_codes"]:
        return result

    count = Counter(seq)
    maximum = max(count.values())
    majority = [base for base in BASES if count[base] == maximum]
    homopolymer_cost = n - maximum
    fixed = _fixed(seq, start0, unit)
    for row in fixed["scores"]:
        template = _finite_template_counts(row["motif_at_interval_start"], n)
        numerator = sum(count[b] * template[b] for b in BASES)
        row.update(template_base_counts=dict(template),
                   composition_null_expected_matches=numerator / n,
                   composition_null_expected_identity=numerator / (n * n),
                   observed_minus_composition_null_identity=row["identity"] - numerator / (n * n))
    fixed["all_best_tie_count"] = len(fixed["all_best_ties"])
    fixed["all_best_distinct_phase_templates"] = len({x["motif_at_interval_start"] for x in fixed["all_best_ties"]})
    fixed["phase_ties_are_not_independent_evidence_or_alignment_paths"] = True
    fixed["composition_null_scope"] = "Expected fixed-template match under exchangeable positions conditional on observed base counts. Not corrected for selecting the best strand/phase; not a p-value or calibrated null test."
    discriminating = []
    for row in fixed["all_best_ties"]:
        row["discriminating_bases"] = [_discriminating(seq, row["motif_at_interval_start"], base) for base in majority]
        discriminating.extend(row["discriminating_bases"])
    best_matches = fixed["all_best_ties"][0]["matches_bp"]
    fixed_gain = best_matches - maximum
    normalized_fixed_gain = fixed_gain / homopolymer_cost if homopolymer_cost else None
    edit = []
    for strand, template in [("+", unit), ("-", tract_context.reverse_complement(unit))]:
        edit.append({"strand": strand, "motif": template, **minimum_periodic_edit_score(seq, template)})
    edit_cost = min(x["edit_cost"] for x in edit)
    edit_gain = homopolymer_cost - edit_cost
    normalized_edit_gain = edit_gain / homopolymer_cost if homopolymer_cost else None
    adequate_copies = n >= len(unit) * settings["minimum_complete_unit_copies"]
    adequate_discriminating = all(x["correct_nonmajority_base_matches"] >= settings["minimum_discriminating_matches"] for x in discriminating)
    precision_recall = all(x["precision"] is not None and x["recall"] is not None
                           and x["precision"] >= settings["minimum_discriminating_precision"]
                           and x["recall"] >= settings["minimum_discriminating_recall"] for x in discriminating)
    identity_ok = best_matches / n >= settings["minimum_fixed_identity"]
    fixed_gain_ok = normalized_fixed_gain is not None and normalized_fixed_gain >= settings["minimum_normalized_fixed_gain"]
    edit_gain_ok = normalized_edit_gain is not None and normalized_edit_gain >= settings["minimum_normalized_edit_gain"]
    fixed_support = adequate_discriminating and precision_recall and identity_ok and fixed_gain_ok
    unit_support = adequate_copies and fixed_support and edit_gain_ok and len(unit) > 1
    status = "unit_specific_evidence" if unit_support else "composition_explains_fit" if fixed_gain <= 0 and edit_gain <= 0 else "ambiguous"
    reasons = []
    if unit_support:
        reasons.append("both_models_and_all_best_fixed_phase_discriminating_checks_support_unit")
    if fixed_gain <= 0:
        reasons.append("fixed_phase_does_not_improve_majority_baseline")
    if edit_gain <= 0:
        reasons.append("cyclic_edit_does_not_improve_best_homopolymer")
    if not adequate_copies:
        reasons.append("insufficient_complete_unit_copies_for_operational_support")
    if not adequate_discriminating:
        reasons.append("insufficient_discriminating_base_matches")
    if not precision_recall:
        reasons.append("best_phase_discriminating_precision_or_recall_below_operational_threshold")
    if not identity_ok:
        reasons.append("fixed_identity_below_operational_threshold")
    if not fixed_gain_ok:
        reasons.append("fixed_normalized_gain_below_threshold_or_undefined")
    if not edit_gain_ok:
        reasons.append("edit_normalized_gain_below_threshold_or_undefined")
    if fixed_support != edit_gain_ok:
        reasons.append("fixed_phase_and_indel_tolerant_operational_support_disagree")
    if len(unit) == 1:
        reasons.append("homopolymer_unit_not_distinguishable_from_composition_baseline")
    result.update(
        status=status, evaluated=True, reason_codes=reasons, fixed_phase=fixed,
        majority_baseline={"base_counts": {b: count[b] for b in BASES}, "all_tied_majority_bases": majority,
                           "matches_bp": maximum, "identity": maximum / n,
                           "best_homopolymer_edit_cost": homopolymer_cost,
                           "all_base_homopolymer_costs": {b: n - count[b] for b in BASES},
                           "homopolymer_edit_cost_identity": "Exact complete-query unit-cost edit distance to any constant template is query_bp minus its base count."},
        cyclic_edit={"method": "minimum_periodic_edit_score: globally consumes query; one free initial phase",
                     "by_strand": edit, "best_edit_cost": edit_cost,
                     "best_edit_cost_per_query_base": edit_cost / n,
                     "all_best_strands": [x["strand"] for x in edit if x["edit_cost"] == edit_cost],
                     "traceback_computed": False, "end_phases_are_not_fixed_start_phases": True},
        model_comparison={"fixed_match_gain_bp_over_majority": fixed_gain,
                          "fixed_match_gain_per_query_base": fixed_gain / n,
                          "fixed_normalized_gain": normalized_fixed_gain,
                          "edit_cost_gain_over_best_homopolymer": edit_gain,
                          "edit_cost_gain_per_query_base": edit_gain / n,
                          "edit_normalized_gain": normalized_edit_gain,
                          "normalization_denominator": homopolymer_cost,
                          "normalization_definition": "Improvement divided by best homopolymer mismatch/edit cost; zero denominator -> null, never positive evidence.",
                          "fixed_phase_beats_majority": fixed_gain > 0,
                          "cyclic_edit_beats_homopolymer": edit_gain > 0,
                          "fixed_operational_support": fixed_support,
                          "cyclic_operational_gain_support": edit_gain_ok,
                          "operational_model_disagreement": fixed_support != edit_gain_ok},
    )
    return result
