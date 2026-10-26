#!/usr/bin/env python3
"""Conditional short-read search power; never a diagnosis of telomere absence.

Genome size, chromosome number and library representation are not estimated
from read totals. Explicit, sourced assumptions are required for a supported
conditional-power classification. See docs/adequacy_model.md for equations.
"""
import argparse
import json
import math
from pathlib import Path
import re

REQUIRED_ASSUMPTIONS = (
    "genome_size_bp", "haploid_chromosome_count", "random_wgs_assumed",
    "telomere_retention_fraction", "effective_independent_read_fraction",
    "detector_sensitivity",
)
REFERENCE_URLS = [
    "https://www.sciencedirect.com/science/article/pii/0888754388900079",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC2940910/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC1560136/",
]
LIMITATIONS = [
    "Conditional power applies to raw exact-seed discovery at the scanner min_reads threshold, not passage through sequence-prior or background-enrichment filters.",
    "Read count, total bases, WGS metadata and missing WGA annotations do not prove random sampling, intact telomere representation or independent molecules.",
    "Effective independent-read fraction applies to seed-eligible reads; paired reads, PCR duplicates and repeated deposits cannot be assumed independent.",
    "The model assumes linear chromosomes, nonoverlapping terminal tracts, uniformly sampled read starts and the specified motif at the assumed chromosome ends.",
    "Nondetection cannot establish true motif absence, ALT or another maintenance mechanism; longer periods, short/interrupted tracts, few motif-bearing ends and library loss remain possible.",
    "The miss probability is conditional on a stated alternative model; it is not a posterior probability that the motif is absent, a p-value or an FDR.",
    "A missing ambiguous-record count is not zero: ambiguous-base totals can bound it conservatively; missing both counters prevents a supported adequacy classification.",
]


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _fraction(value, name, allow_zero=True):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite fraction")
    if not (0 <= value <= 1) or (not allow_zero and value == 0):
        raise ValueError(f"{name} must lie in {'[0,1]' if allow_zero else '(0,1]'}")
    return float(value)


def required_seed_bases(period, seed_bases, min_copies):
    """Matches src/detect.rs: whole motif copies, including ceiling to a period."""
    for value, name in ((period, "period"), (seed_bases, "seed_bases"), (min_copies, "min_copies")):
        _integer(value, name, 1)
    return period * max(min_copies, (seed_bases + period - 1) // period)


def eligible_read_bounds(records, bases, maximum, minimum_seed, histogram=None):
    """Rigorous length-only count bounds from totals and optional interval bins."""
    _integer(records, "records")
    _integer(bases, "bases")
    _integer(maximum, "max_record_bases")
    _integer(minimum_seed, "minimum_seed", 1)
    if bases > records * maximum:
        raise ValueError("bases exceeds records * max_record_bases")
    if not records or maximum < minimum_seed:
        return {"lower": 0, "upper": 0, "method": "read_lengths_cannot_reach_seed"}
    # At most x*M + (N-x)*(d-1) bases can be carried by x reads of length >= d.
    excess = max(0, bases - records * (minimum_seed - 1))
    denominator = maximum - minimum_seed + 1
    lower = (excess + denominator - 1) // denominator
    upper = min(records, bases // minimum_seed)
    method = "aggregate_length_bound"
    if histogram is not None:
        histogram_lower = sum(row["records"] for row in histogram if row["lower"] >= minimum_seed)
        histogram_upper = sum(row["records"] for row in histogram if row["upper"] >= minimum_seed)
        lower, upper = max(lower, histogram_lower), min(upper, histogram_upper)
        method = "histogram_and_aggregate_length_bounds"
    return {"lower": lower, "upper": upper, "method": method}


def binomial_miss_bound(trials, probability, min_hits):
    """P[X < min_hits]; exact finite sum for <= 10,000 terms, else upper bound."""
    _integer(trials, "trials")
    _integer(min_hits, "min_hits", 1)
    probability = _fraction(probability, "probability")
    if trials < min_hits or probability == 0:
        return {"miss_probability_upper_bound": 1.0, "log_miss_probability_upper_bound": 0.0, "method": "exact_binomial"}
    if probability == 1:
        return {"miss_probability_upper_bound": 0.0, "log_miss_probability_upper_bound": None, "method": "exact_binomial"}
    if min_hits > 10000:
        a = (min_hits - 1) / trials
        if a >= probability:
            log_value = 0.0
        else:
            divergence = (a * math.log(a / probability) if a else 0.0)
            divergence += (1 - a) * math.log((1 - a) / (1 - probability))
            log_value = min(0.0, -trials * divergence)
        method = "binomial_chernoff_upper_bound"
    else:
        log_term = trials * math.log1p(-probability)
        log_sum = log_term
        log_odds = math.log(probability) - math.log1p(-probability)
        for hits in range(1, min_hits):
            log_term += math.log(trials - hits + 1) - math.log(hits) + log_odds
            maximum = max(log_sum, log_term)
            log_sum = maximum + math.log1p(math.exp(-abs(log_sum - log_term)))
        log_value, method = min(0.0, log_sum), "exact_binomial"
    return {"miss_probability_upper_bound": math.exp(log_value),
            "log_miss_probability_upper_bound": log_value, "method": method,
            "floating_point_underflow": log_value < -745}


def _read_qc(report):
    summary, config = report.get("summary", {}), report.get("config", {})
    if not isinstance(summary, dict) or not isinstance(config, dict):
        raise ValueError("summary and config must be objects")
    records = _integer(summary.get("records"), "summary.records")
    bases = _integer(summary.get("bases"), "summary.bases")
    maximum = _integer(summary.get("max_record_bases"), "summary.max_record_bases")
    if bases > records * maximum or (not records and bases) or maximum > bases:
        raise ValueError("inconsistent read/base/maximum-length totals")
    qc = report.get("evidence_qc", {})
    if not isinstance(qc, dict):
        raise ValueError("evidence_qc must be an object")
    for key, expected in (("records", records), ("bases", bases), ("max_record_bases", maximum)):
        if key in qc and qc[key] != expected:
            raise ValueError(f"evidence_qc.{key} disagrees with summary")
    histogram = qc.get("length_histogram")
    if histogram is not None:
        if not isinstance(histogram, list):
            raise ValueError("length_histogram must be an array")
        last, count, low_bases, high_bases = 0, 0, 0, 0
        for row in histogram:
            lower = _integer(row.get("lower"), "histogram.lower", 1)
            upper = _integer(row.get("upper"), "histogram.upper", lower)
            size = _integer(row.get("records"), "histogram.records", 1)
            if lower <= last or lower > maximum:
                raise ValueError("histogram bins overlap, are unsorted or exceed observed maximum")
            last, count = upper, count + size
            low_bases += lower * size
            high_bases += min(upper, maximum) * size
        if count != records or not low_bases <= bases <= high_bases:
            raise ValueError("histogram contradicts total records/bases")
    ambiguous_records = qc.get("records_with_ambiguity")
    ambiguous_bases = qc.get("ambiguous_bases")
    if ambiguous_records is not None:
        _integer(ambiguous_records, "records_with_ambiguity")
        if ambiguous_records > records:
            raise ValueError("records_with_ambiguity exceeds records")
    if ambiguous_bases is not None:
        _integer(ambiguous_bases, "ambiguous_bases")
        if ambiguous_bases > bases:
            raise ValueError("ambiguous_bases exceeds bases")
    if ambiguous_records is not None and ambiguous_bases is not None:
        if not ambiguous_records <= ambiguous_bases <= ambiguous_records * maximum:
            raise ValueError("ambiguous record/base counts contradict each other")
    ambiguity = ambiguity_record_bounds(records, ambiguous_records, ambiguous_bases)
    from quality_observation import describe_quality
    quality = describe_quality(report, None)
    return {"records": records, "bases": bases, "max_record_bases": maximum,
            "mean_record_bases": bases / records if records else None,
            "length_histogram": histogram, "records_with_ambiguity": ambiguous_records,
            "ambiguous_bases": ambiguous_bases, "sampled_prefix": summary.get("sampled_prefix", None),
            "ambiguity_record_bounds": ambiguity,
            "length_distribution_available": histogram is not None,
            "fixed_length_proven_by_totals": bool(records and bases == records * maximum),
            "duplicate_fraction_measured": False, "terminal_library_retention_measured": False,
            "base_quality_distribution_available": quality['status'] in ('measured', 'partially_measured'),
            "base_quality_observation_status": quality['status'],
            "base_quality_used_for_power": False}


def ambiguity_record_bounds(records, ambiguous_records=None, ambiguous_bases=None):
    """A missing counter cannot certify ambiguity-free reads.

    Every ambiguous record contains at least one ambiguous base, so the number
    of such records is at most min(records, ambiguous_bases). The upper bound
    is appropriate when subtracting potentially unusable reads from a lower
    bound on seed-eligible reads. Unknown totals yield no proven clean reads.
    """
    if ambiguous_records is not None:
        return {"status": "measured_record_count", "upper": ambiguous_records,
                "record_count_measured": True, "bound_available": True}
    if ambiguous_bases is not None:
        return {"status": "bounded_from_ambiguous_bases", "upper": min(records, ambiguous_bases),
                "record_count_measured": False, "bound_available": True}
    return {"status": "unavailable", "upper": records,
            "record_count_measured": False, "bound_available": False}


def library_assessment(metadata):
    """Read explicit metadata; an absent amplification annotation proves nothing."""
    blockers, evidence = [], []
    strategy = str(metadata.get("library_strategy", "")).strip().upper()
    source = str(metadata.get("library_source", "")).strip().upper()
    selection = str(metadata.get("library_selection", "")).strip().upper()
    if strategy and strategy != "WGS":
        blockers.append("library_strategy_not_random_wgs")
        evidence.append({"field": "library_strategy", "value": metadata["library_strategy"]})
    if source and source != "GENOMIC":
        blockers.append("library_source_not_single_genome_dna")
        evidence.append({"field": "library_source", "value": metadata["library_source"]})
    if selection in {"PCR", "RANDOM PCR", "HYBRID SELECTION", "CDNA", "POLYA", "RT-PCR"}:
        blockers.append("library_selection_not_verified_uniform_genome_sampling")
        evidence.append({"field": "library_selection", "value": metadata["library_selection"]})
    if strategy in {"WGA", "MDA"} or selection in {"WGA", "MDA", "MALBAC", "DOP-PCR"}:
        blockers.append("explicit_whole_genome_amplification")
        if selection in {"WGA", "MDA", "MALBAC", "DOP-PCR"}:
            evidence.append({"field": "library_selection", "value": metadata["library_selection"]})
    if str(metadata.get("metagenome", "")).lower() in {"true", "yes", "1"}:
        blockers.append("mixed_metagenomic_target_not_modelled")
    for key in ("whole_genome_amplification", "wga", "mda"):
        if str(metadata.get(key, "")).lower() in {"true", "yes", "1"}:
            blockers.append("explicit_whole_genome_amplification")
            evidence.append({"field": key, "value": metadata[key]})
    method = str(metadata.get("amplification_method", "")).strip().upper()
    if method in {"MDA", "WGA", "MALBAC", "DOP-PCR", "REPLI-G", "GENOMIPHI"}:
        blockers.append("explicit_whole_genome_amplification")
        evidence.append({"field": "amplification_method", "value": metadata["amplification_method"]})
    for key in ("library_construction_protocol", "library_construction", "library_description"):
        text = str(metadata.get(key, ""))
        if re.search(r"\bWGA\b|\bMDA\b|whole[ -]genome amplification|multiple displacement amplification", text, re.I):
            blockers.append("amplification_protocol_requires_manual_review")
            evidence.append({"field": key, "value": text})
    return {"uniform_wgs_model_blockers": sorted(set(blockers)), "observed_metadata_evidence": evidence,
            "amplification_absence_verified": False,
            "note": "A protocol mention, including negated wording, triggers review rather than a positive WGA diagnosis; ordinary WGS labels do not establish absence of amplification."}


def _validate_assumptions(assumptions):
    if not isinstance(assumptions, dict):
        raise ValueError("assumptions must be an object")
    for key in ("genome_size_bp", "haploid_chromosome_count", "assumed_tract_bases", "motif_bearing_chromosome_ends"):
        if key in assumptions:
            _integer(assumptions[key], key, 1)
    for key in ("telomere_retention_fraction", "effective_independent_read_fraction", "detector_sensitivity"):
        if key in assumptions:
            _fraction(assumptions[key], key)
    if "random_wgs_assumed" in assumptions and not isinstance(assumptions["random_wgs_assumed"], bool):
        raise ValueError("random_wgs_assumed must be a boolean")
    threshold = assumptions.get("max_miss_probability", 0.05)
    _fraction(threshold, "max_miss_probability", allow_zero=False)
    if threshold >= 1:
        raise ValueError("max_miss_probability must be less than 1")
    provenance = assumptions.get("provenance", {})
    if not isinstance(provenance, dict) or any(not isinstance(v, str) or not v.strip() for v in provenance.values()):
        raise ValueError("assumptions.provenance must map parameters to nonempty source/assumption descriptions")


def assessing_power(qc, config, assumptions):
    """Calculate sensitivity under explicit assumptions, without assigning adequacy."""
    _validate_assumptions(assumptions)
    for key in REQUIRED_ASSUMPTIONS:
        if key not in assumptions:
            raise ValueError(f"missing power assumption {key}")
    genome = assumptions["genome_size_bp"]
    chromosomes = assumptions["haploid_chromosome_count"]
    tract = assumptions.get("assumed_tract_bases", 500)
    ends = assumptions.get("motif_bearing_chromosome_ends", 2 * chromosomes)
    if ends > 2 * chromosomes or genome < 2 * chromosomes * tract:
        raise ValueError("assumed chromosome ends/tracts are incompatible with genome size or chromosome number")
    minimum = _integer(config.get("min_k"), "config.min_k", 1)
    maximum = min(32, _integer(config.get("max_k"), "config.max_k", 1))
    if minimum > maximum:
        raise ValueError("no scanned periods in the supported 1..32 range")
    seed = _integer(config.get("seed_bases"), "config.seed_bases", 1)
    copies = _integer(config.get("min_copies"), "config.min_copies", 2)
    min_hits = _integer(config.get("min_reads"), "config.min_reads", 1)
    rows = []
    ambiguity = ambiguity_record_bounds(qc["records"], qc.get("records_with_ambiguity"), qc.get("ambiguous_bases"))
    for period in range(minimum, maximum + 1):
        required = required_seed_bases(period, seed, copies)
        eligible = eligible_read_bounds(qc["records"], qc["bases"], qc["max_record_bases"], required, qc.get("length_histogram"))
        clean_lower = max(0, eligible["lower"] - ambiguity["upper"])
        trials = min(clean_lower, math.floor(clean_lower * assumptions["effective_independent_read_fraction"]))
        exposure_per_end = max(0, tract - required + 1)
        geometry = ends * exposure_per_end / genome
        probability = geometry * assumptions["telomere_retention_fraction"] * assumptions["detector_sensitivity"]
        miss = binomial_miss_bound(trials, probability, min_hits)
        rows.append({"period": period, "required_exact_seed_bases": required,
                     "eligible_read_count_bounds": eligible, "ambiguity_free_eligible_reads_lower_bound": clean_lower,
                     "ambiguity_record_count_bound": dict(ambiguity),
                     "effective_independent_trials_assumed_lower_bound": trials,
                     "capture_start_positions_per_end": exposure_per_end,
                     "per_trial_capture_probability_assumed_lower_bound": probability,
                     "expected_independent_detected_reads_lower_bound": trials * probability,
                     "detection_probability_lower_bound": 1 - miss["miss_probability_upper_bound"], **miss})
    worst = max(rows, key=lambda row: row["miss_probability_upper_bound"])
    threshold = assumptions.get("max_miss_probability", 0.05)
    return {"model": "conditional_terminal_uniform_start_binomial_v1", "periods": rows,
            "assumed_terminal_tract_bases": tract, "assumed_motif_bearing_ends": ends,
            "min_independent_detected_reads_required": min_hits,
            "worst_case_period": worst["period"],
            "worst_case_miss_probability_upper_bound": worst["miss_probability_upper_bound"],
            "minimum_detection_probability_lower_bound": 1 - worst["miss_probability_upper_bound"],
            "max_miss_probability_target": threshold,
            "all_scanned_periods_up_to_32_meet_conditional_target": worst["miss_probability_upper_bound"] <= threshold,
            "mean_sequence_depth_given_genome_size": qc["bases"] / genome,
            "zero_probability_does_not_prove_absence": True}


def illustrative_grid(qc, config):
    rows = []
    for genome in (100_000_000, 1_000_000_000, 10_000_000_000):
        for chromosomes in (6, 20):
            for retention in (1.0, 0.1, 0.01):
                assumptions = {"genome_size_bp": genome, "haploid_chromosome_count": chromosomes,
                               "random_wgs_assumed": True, "telomere_retention_fraction": retention,
                               "effective_independent_read_fraction": 0.5, "detector_sensitivity": 0.5}
                power = assessing_power(qc, config, assumptions)
                rows.append({**assumptions, "assumed_tract_bases": 500,
                             "worst_case_miss_probability_upper_bound": power["worst_case_miss_probability_upper_bound"],
                             "minimum_detection_probability_lower_bound": power["minimum_detection_probability_lower_bound"]})
    return {"status": "illustrative_only_not_sample_parameters_or_adequacy_evidence", "scenarios": rows,
            "note": "Hypothetical values illustrate sensitivity. They are not estimated from this sample, not species-specific and never upgrade adequacy status."}


def assess_adequacy(report, metadata=None, assumptions=None, include_illustrative_grid=True):
    """Return qualified data adequacy independently of motif nomination/classification."""
    metadata, assumptions = metadata or {}, assumptions or {}
    _validate_assumptions(assumptions)
    result = {"schema_version": 1, "assessment": "qualified_short_read_nondetection",
              "status": "unassessable_library_or_genome", "reasons": [], "read_qc": None,
              "library": library_assessment(metadata), "assumptions": assumptions,
              "power": None, "limitations": LIMITATIONS, "references": REFERENCE_URLS,
              "true_absence_inferred": False, "alternative_maintenance_inferred": False,
              "motif_nomination_performed": False}
    failed = not isinstance(report, dict) or any(str(source.get(key, "")).lower() in {"failed", "error", "truncated"}
               for source in (report if isinstance(report, dict) else {}, metadata)
               for key in ("status", "state", "analysis_status"))
    if failed:
        result.update(status="observed_qc_failure", reasons=["scan_failure_or_missing_report"],
                      recommendation="Resolve the observed input/scan failure; no biological nondetection interpretation is justified.")
        return result
    try:
        qc = _read_qc(report)
    except (ValueError, AttributeError, TypeError) as error:
        result.update(status="observed_qc_failure", reasons=["invalid_or_inconsistent_report_qc"],
                      qc_error=str(error), recommendation="Validate the report and input integrity before interpreting nondetection.")
        return result
    result["read_qc"] = qc
    config = report.get("config", {})
    if not qc["records"] or not qc["bases"]:
        result.update(status="insufficient_data", reasons=["no_observed_reads_or_bases"],
                      recommendation="Acquire usable sequence data; zero observations do not imply motif absence.")
        return result
    if config.get("mode") != "reads":
        result.update(reasons=["assembly_not_short_read_sampling"],
                      recommendation="Use assembly-terminal evidence; this read-sampling model does not apply to contig counts.")
        return result
    try:
        min_k = _integer(config.get("min_k"), "config.min_k", 1)
        raw_max_k = _integer(config.get("max_k"), "config.max_k", 1)
        max_k = min(32, raw_max_k)
        seed = _integer(config.get("seed_bases"), "config.seed_bases", 1)
        copies = _integer(config.get("min_copies"), "config.min_copies", 2)
        min_reads = _integer(config.get("min_reads"), "config.min_reads", 1)
        result["detector_scope"] = {"min_k": min_k, "max_k": raw_max_k,
                                    "model_max_k": max_k, "seed_bases": seed,
                                    "min_copies": copies, "min_reads": min_reads,
                                    "model": "raw_single_period_exact_seed_discovery"}
        if min_k > max_k:
            result.update(reasons=["no_scanned_period_in_short_repeat_model_range"],
                          recommendation="This assessment covers scanned periods up to 32 bases; use a model appropriate to the actual longer-period search.")
            return result
    except ValueError as error:
        result.update(status="observed_qc_failure", reasons=["invalid_or_unsupported_detector_config"], qc_error=str(error))
        return result
    requirements = [required_seed_bases(k, seed, copies) for k in range(min_k, max_k + 1)]
    if qc["records"] < min_reads or qc["max_record_bases"] < min(requirements):
        result.update(status="insufficient_data", reasons=["read_count_or_length_cannot_reach_scanner_threshold"],
                      recommendation="Acquire more usable or longer reads before interpreting this short-repeat search.")
        return result
    missing = [key for key in REQUIRED_ASSUMPTIONS if key not in assumptions]
    missing_provenance = [key for key in REQUIRED_ASSUMPTIONS if key in assumptions and not assumptions.get("provenance", {}).get(key)]
    result["missing_assumptions"] = missing
    result["missing_assumption_provenance"] = missing_provenance
    blockers = list(result["library"]["uniform_wgs_model_blockers"])
    if not qc["ambiguity_record_bounds"]["bound_available"]:
        blockers.append("ambiguity_qc_unavailable")
    if qc["sampled_prefix"] is not False:
        blockers.append("input_completeness_or_prefix_sampling_not_verified")
    if assumptions.get("random_wgs_assumed") is False:
        blockers.append("random_wgs_model_explicitly_not_assumed")
    if missing:
        result["reasons"].append("missing_genome_or_library_assumptions")
    if missing_provenance:
        result["reasons"].append("assumptions_lack_explicit_provenance")
    result["reasons"].extend(blockers)
    if not missing:
        result["power"] = assessing_power(qc, config, assumptions)
    if not missing and not missing_provenance and not blockers:
        if result["power"]["all_scanned_periods_up_to_32_meet_conditional_target"]:
            result.update(status="conditional_power_supported", reasons=["stated_conditional_raw_discovery_power_target_met"])
        else:
            result.update(status="insufficient_data", reasons=["conditional_raw_discovery_power_below_target"])
    if include_illustrative_grid:
        result["illustrative_grid"] = illustrative_grid(qc, config)
    result["recommendation"] = (
        "If a separate evidence review finds no short-repeat candidate, report nondetection conditional on these assumptions; review longer motifs, terminal structure and library retention before any biological absence claim."
        if result["status"] == "conditional_power_supported" else
        "Document genome size/karyotype and library representation, independent-molecule and detector-sensitivity assumptions; inspect raw candidates and obtain terminal or complementary evidence where nondetection matters."
    )
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path)
    parser.add_argument("--metadata-json", type=Path)
    parser.add_argument("--assumptions-json", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    value = assess_adequacy(json.loads(args.result.read_text()),
                           json.loads(args.metadata_json.read_text()) if args.metadata_json else None,
                           json.loads(args.assumptions_json.read_text()) if args.assumptions_json else None)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
