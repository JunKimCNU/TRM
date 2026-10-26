"""Bounded, descriptive repeat-tract context and a conservative pair hypothesis.

Coordinates always refer to the supplied full record, zero-based and half-open.
This module does not discover, merge, relabel, or localize native repeat calls.
A compatible pair means only that its supplied sequence passes the explicitly
reported bounded model. It is not proof of a biological array or chromosome end.
There is no dependency on frozen, sample-specific result directories.
"""

from collections import Counter
import hashlib
import math
import statistics


VERSION = "0.1.0"
SCHEMA_VERSION = 1
MAX_ANALYZED_INTERVAL_BP = 100_000
MAX_AUXILIARY_RECORDS = 10_000
COMPLEMENT = str.maketrans("ACGT", "TGCA")
DEFAULT_PAIR_CONFIG = {
    "max_gap_bp": 64,
    "max_span_bp": 10_000,
    "max_uncovered_run_bp": 7,
    "min_gap_coverage": 0.8,
    "min_gap_fixed_phase_identity": 0.9,
    "min_anchor_fixed_phase_identity": 0.9,
    "min_exact_copies": 2,
    "require_same_phase": True,
}


def reverse_complement(motif):
    return motif.translate(COMPLEMENT)[::-1]


def _integer(value, name, minimum=0, maximum=None):
    if (isinstance(value, bool) or not isinstance(value, int)
            or value < minimum or (maximum is not None and value > maximum)):
        raise ValueError(f"{name} must be an integer between {minimum} and {maximum}")
    return value


def _motif(motif):
    if not isinstance(motif, str):
        raise ValueError("Motif must be an ACGT string")
    motif = motif.upper()
    if not 1 <= len(motif) <= 64 or any(x not in "ACGT" for x in motif):
        raise ValueError("Motif must contain 1 through 64 ACGT bases")
    for period in range(1, len(motif)):
        if len(motif) % period == 0 and motif == motif[:period] * (len(motif) // period):
            raise ValueError(f"Motif is nonprimitive; specify {motif[:period]!r} explicitly")
    return motif


def _strand(value, name="Expected strand", optional=True):
    if value not in ((None, "+", "-") if optional else ("+", "-")):
        raise ValueError(f"{name} must be + or -" + (" or unspecified" if optional else ""))


def _bounds(sequence, start0, end0):
    if not isinstance(sequence, str):
        raise ValueError("Sequence must be a full-record string")
    _integer(start0, "start0", maximum=len(sequence))
    _integer(end0, "end0", minimum=start0, maximum=len(sequence))


def _chunk(sequence, start0, end0):
    _bounds(sequence, start0, end0)
    if end0 - start0 > MAX_ANALYZED_INTERVAL_BP:
        raise ValueError(f"Requested interval exceeds {MAX_ANALYZED_INTERVAL_BP} bp work bound")
    chunk = sequence[start0:end0]
    if any(not ("A" <= base <= "Z" or "a" <= base <= "z") for base in chunk):
        raise ValueError("Analyzed sequence must contain ASCII letters only")
    return chunk.upper()


def _record_context(sequence, start0, end0):
    return dict(record_length_bp=len(sequence), left_record_distance_bp=start0,
                right_record_distance_bp=len(sequence) - end0,
                interval_left_cropped=start0 > 0, interval_right_cropped=end0 < len(sequence),
                record_start_reached=start0 == 0, record_end_reached=end0 == len(sequence),
                record_edge_is_not_physical_chromosome_end=True)


def union_intervals(intervals):
    """Numerical set union for coverage accounting; never a new tract call."""
    output = []
    for start, end in sorted(intervals):
        if start >= end:
            continue
        if output and start <= output[-1][1]:
            output[-1][1] = max(output[-1][1], end)
        else:
            output.append([start, end])
    return output


def complement_intervals(start, end, covered):
    islands, position = [], start
    for a, b in union_intervals(covered):
        a, b = max(a, start), min(b, end)
        if a >= b:
            continue
        if a > position:
            islands.append([position, a])
        position = max(position, b)
    if position < end:
        islands.append([position, end])
    return islands


def quality_summary(quality):
    if quality is None:
        return dict(available=False, encoding=None)
    if not isinstance(quality, str):
        raise ValueError("Phred+33 quality must be a string")
    values = [ord(char) - 33 for char in quality]
    if any(x < 0 or x > 93 for x in values):
        raise ValueError("Invalid Phred+33 character")
    return dict(available=True, encoding="Phred+33", bases=len(values),
                minimum=min(values) if values else None, maximum=max(values) if values else None,
                mean=statistics.mean(values) if values else None,
                median=statistics.median(values) if values else None,
                q_below_10=sum(x < 10 for x in values), q_below_20=sum(x < 20 for x in values),
                q_below_30=sum(x < 30 for x in values),
                interpretation="Base-call quality observations, not calibrated experiment-specific error probabilities")


def composition(sequence):
    counts = Counter(sequence)
    return dict(length_bp=len(sequence), base_counts=dict(sorted(counts.items())),
                acgt_bp=sum(counts[x] for x in "ACGT"),
                non_acgt_bp=sum(n for x, n in counts.items() if x not in "ACGT"))


def _coverage(sequence, quality, start, end, intervals):
    union = union_intervals((max(start, a), min(end, b)) for a, b in intervals)
    covered = sum(b - a for a, b in union)
    islands = complement_intervals(start, end, union)
    return dict(covered_bp=covered, uncovered_bp=end - start - covered,
                covered_fraction=covered / (end - start) if end > start else None,
                coverage_union_intervals=union,
                uncovered_islands=[dict(start0=a, end0=b, **composition(sequence[a:b].upper()),
                                        quality=quality_summary(quality[a:b] if quality is not None else None))
                                   for a, b in islands],
                union_is_base_accounting_only=True, biological_tracts_merged=False)


def _exact_segments(sequence, chunk, start0, end0, motif):
    k = len(motif)
    result = []
    for strand, template in [("+", motif), ("-", reverse_complement(motif))]:
        for rotation in range(k):
            begin = None
            for offset in range(len(chunk) + 1):
                match = offset < len(chunk) and chunk[offset] == template[(offset + rotation) % k]
                if match and begin is None:
                    begin = offset
                if not match and begin is not None:
                    if offset - begin >= k:
                        a, b = start0 + begin, start0 + offset
                        absolute_phase = (rotation - start0) % k
                        segment_rotation = (a + absolute_phase) % k
                        result.append(dict(start0=a, end0=b, length_bp=b - a, strand=strand,
                                           interval_start_rotation=rotation,
                                           absolute_template_phase=absolute_phase,
                                           rotation_at_segment_start=template[segment_rotation:] + template[:segment_rotation],
                                           complete_copies=(b - a) // k, partial_copy_bp=(b - a) % k,
                                           left_interval_censored=a == start0, right_interval_censored=b == end0,
                                           left_record_distance_bp=a, right_record_distance_bp=len(sequence) - b,
                                           exact_continuation_left_observed=a == start0 and a > 0 and sequence[a - 1].upper() == template[(a - 1 + absolute_phase) % k],
                                           exact_continuation_right_observed=b == end0 and b < len(sequence) and sequence[b].upper() == template[(b + absolute_phase) % k]))
                    begin = None
    result.sort(key=lambda r: (r["start0"], r["end0"], r["strand"], r["absolute_template_phase"]))
    return result


def exact_segments(sequence, start0, end0, motif):
    """Return maximal exact runs of at least one primitive period in the interval."""
    return _exact_segments(sequence, _chunk(sequence, start0, end0), start0, end0, _motif(motif))


def _fixed_phase_scores(chunk, start0, motif):
    scores = []
    for strand, template in [("+", motif), ("-", reverse_complement(motif))]:
        for rotation in range(len(motif)):
            matches = sum(base == template[(i + rotation) % len(motif)] for i, base in enumerate(chunk))
            scores.append(dict(strand=strand, rotation_at_interval_start=rotation,
                               absolute_template_phase=(rotation - start0) % len(motif),
                               motif_at_interval_start=template[rotation:] + template[:rotation],
                               matches_bp=matches, mismatches_bp=len(chunk) - matches,
                               denominator_bp=len(chunk), identity=matches / len(chunk) if chunk else None))
    best_count = max(s["matches_bp"] for s in scores)
    return dict(model="One fixed phase across the entire interval; substitutions/non-ACGT count as mismatches; no indels or phase switching",
                scores=scores, all_best_ties=[s for s in scores if s["matches_bp"] == best_count],
                best_by_strand={strand: [s for s in scores if s["strand"] == strand and s["matches_bp"] == max(x["matches_bp"] for x in scores if x["strand"] == strand)]
                                for strand in ["+", "-"]},
                empty_interval=not chunk, continuity_inference_permitted=False)


def fixed_phase_scores(chunk, start0, motif):
    """Score a supplied interval string with its original-record start coordinate."""
    _integer(start0, "start0")
    if not isinstance(chunk, str):
        raise ValueError("Interval sequence must be a string")
    return _fixed_phase_scores(_chunk(chunk, 0, len(chunk)), start0, _motif(motif))


def adjacent_phase_observations(segments, k):
    output = []
    for strand in ["+", "-"]:
        ordered = sorted((s for s in segments if s["strand"] == strand), key=lambda s: (s["start0"], s["end0"], s["absolute_template_phase"]))
        for left, right in zip(ordered, ordered[1:]):
            delta = (right["absolute_template_phase"] - left["absolute_template_phase"]) % k
            output.append(dict(strand=strand, left_start0=left["start0"], left_end0=left["end0"],
                               right_start0=right["start0"], right_end0=right["end0"],
                               left_absolute_phase=left["absolute_template_phase"], right_absolute_phase=right["absolute_template_phase"],
                               observed_phase_delta_mod_period=delta, observed_phase_changes=delta != 0,
                               gap_bp_negative_is_overlap=right["start0"] - left["end0"],
                               interpretation="Neighboring retained exact segments, not inferred indels or a continuity path"))
    return output


def _tract(sequence, tract, index):
    if not isinstance(tract, dict):
        raise ValueError("Each known tract must be a dictionary")
    if not {"start0", "end0", "strand"} <= tract.keys():
        raise ValueError("Known tract requires start0, end0, and strand")
    _bounds(sequence, tract["start0"], tract["end0"])
    if tract["start0"] == tract["end0"]:
        raise ValueError("Known tract must have positive length")
    _strand(tract["strand"], "Known tract strand", optional=False)
    if "tract_id" in tract and (not isinstance(tract["tract_id"], str) or not tract["tract_id"].strip()):
        raise ValueError("Supplied tract_id must be a nonempty string")
    if "identity" in tract and tract["identity"] is not None:
        identity = tract["identity"]
        if (isinstance(identity, bool) or not isinstance(identity, (int, float))
                or not math.isfinite(identity) or not 0 <= identity <= 1):
            raise ValueError("Known tract identity must be finite and between 0 and 1")
    return dict(tract, tract_id=tract.get("tract_id", f"tract_{index}"))


def analyze_interval(sequence, start0, end0, motif, quality=None, known_tracts=(), expected_strand=None, cutoffs=None):
    """Describe a bounded interval of a full record without changing tract calls.

Primitive motifs have 1--64 ACGT bases; repeated representations are rejected
instead of silently changing the requested period. Lowercase bases are accepted.
IUPAC/non-ACGT letters mismatch exact templates. Default cutoffs are one, two,
and four periods. Explicit cutoffs must be distinct integer base counts >= period.
Only analyzed sequence and quality characters are validated, to avoid rescanning
an entire chromosome for every small interval; full-record lengths are checked.
"""
    motif = _motif(motif)
    chunk = _chunk(sequence, start0, end0)
    _strand(expected_strand)
    if quality is not None and (not isinstance(quality, str) or len(quality) != len(sequence)):
        raise ValueError("Sequence and Phred+33 quality lengths differ")
    interval_quality = quality_summary(quality[start0:end0] if quality is not None else None)
    if cutoffs is None:
        cutoffs = (len(motif), 2 * len(motif), 4 * len(motif))
    if not isinstance(cutoffs, (tuple, list)) or not 1 <= len(cutoffs) <= 64:
        raise ValueError("cutoffs must contain 1 through 64 distinct integer thresholds")
    for cutoff in cutoffs:
        _integer(cutoff, "Exact-run cutoff", minimum=len(motif), maximum=MAX_ANALYZED_INTERVAL_BP)
    if len(set(cutoffs)) != len(cutoffs):
        raise ValueError("Exact-run cutoffs must be distinct")
    segments = _exact_segments(sequence, chunk, start0, end0, motif)
    views = {}
    for cutoff in cutoffs:
        selected = [s for s in segments if s["length_bp"] >= cutoff]
        view = dict(minimum_exact_run_bp=cutoff, segment_count=len(selected),
                    raw_segment_length_sum_bp=sum(s["length_bp"] for s in selected),
                    all_strands=_coverage(sequence, quality, start0, end0, [(s["start0"], s["end0"]) for s in selected]),
                    by_strand={strand: _coverage(sequence, quality, start0, end0, [(s["start0"], s["end0"]) for s in selected if s["strand"] == strand]) for strand in ["+", "-"]},
                    adjacent_phase_observations=adjacent_phase_observations(selected, len(motif)))
        view["double_counted_overlap_bp"] = view["raw_segment_length_sum_bp"] - view["all_strands"]["covered_bp"]
        views[str(cutoff)] = view
    intersections = []
    for index, supplied in enumerate(known_tracts):
        if index >= MAX_AUXILIARY_RECORDS:
            raise ValueError("Too many supplied known tracts")
        tract = _tract(sequence, supplied, index)
        a, b = tract["start0"], tract["end0"]
        lo, hi = max(a, start0), min(b, end0)
        if lo < hi:
            intersections.append(dict(tract_id=tract["tract_id"], original_start0=a, original_end0=b,
                                      original_length_bp=b - a, strand=tract["strand"], original_identity=tract.get("identity"),
                                      intersection_start0=lo, intersection_end0=hi, intersection_bp=hi - lo,
                                      clipped_by_interval=lo != a or hi != b,
                                      existing_joint_long_terminal_context=tract.get("joint_long_terminal_context"),
                                      existing_localization_status=tract.get("localization_status")))
    known_coverage = _coverage(sequence, quality, start0, end0, [(r["intersection_start0"], r["intersection_end0"]) for r in intersections])
    known_coverage["raw_intersection_length_sum_bp"] = sum(r["intersection_bp"] for r in intersections)
    known_coverage["double_counted_overlap_bp"] = known_coverage["raw_intersection_length_sum_bp"] - known_coverage["covered_bp"]
    direction_ambiguous = reverse_complement(motif) in motif + motif
    return dict(schema_version=SCHEMA_VERSION, module_version=VERSION,
                coordinate_system="original record, zero-based half-open",
                start0=start0, end0=end0, length_bp=end0 - start0, motif=motif, motif_period_bp=len(motif),
                expected_strand=expected_strand,
                expected_motif_in_original_record=motif if expected_strand == "+" else reverse_complement(motif) if expected_strand == "-" else None,
                reverse_complement_rotation_equivalent=direction_ambiguous,
                motif_strand_templates_distinguishable=not direction_ambiguous,
                record_context=_record_context(sequence, start0, end0),
                uppercase_interval_sequence_sha256=hashlib.sha256(chunk.encode("ascii")).hexdigest(),
                composition=composition(chunk), quality=interval_quality,
                exact_segments_ge_period=segments, exact_coverage_views=views,
                fixed_phase=_fixed_phase_scores(chunk, start0, motif),
                known_tract_intersections=intersections, known_tract_coverage=known_coverage,
                exact_cutoffs_bp=list(cutoffs), cutoff_interpretation="Descriptive within-interval exact-run lengths, not native detection thresholds",
                limitations=["Exact runs are censored by the requested interval; cutoffs apply to their within-interval lengths.",
                             "Coverage union counts bases; it does not establish one continuous repeat array.",
                             "Fixed-phase scoring allows mismatches, but no phase switching or indels.",
                             "Observed phase differences do not establish their biological or technical cause.",
                             "Expected direction is supplied metadata; both strand templates are evaluated.",
                             "A record edge or expected direction is not physical chromosome-end evidence.",
                             "No native, policy, adequacy, or physical-terminal classification changes."],
                biological_tracts_merged=False, changes_native_classification=False,
                changes_policy_classification=False, physical_chromosome_end_confirmed=False)


def validate_pair_config(config=None):
    """Return validated complete model settings without modifying caller input."""
    if config is not None and not isinstance(config, dict):
        raise ValueError("Pair config must be a dictionary")
    config = {} if config is None else config
    unknown = set(config) - set(DEFAULT_PAIR_CONFIG)
    if unknown:
        raise ValueError(f"Unknown pair config keys: {sorted(map(str, unknown))}")
    values = dict(DEFAULT_PAIR_CONFIG, **config)
    for name in ("max_gap_bp", "max_uncovered_run_bp"):
        _integer(values[name], name, maximum=MAX_ANALYZED_INTERVAL_BP)
    _integer(values["max_span_bp"], "max_span_bp", minimum=1, maximum=MAX_ANALYZED_INTERVAL_BP)
    _integer(values["min_exact_copies"], "min_exact_copies", minimum=1, maximum=1000)
    for name in ("min_gap_coverage", "min_gap_fixed_phase_identity", "min_anchor_fixed_phase_identity"):
        value = values[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"{name} must be finite and between 0 and 1")
    if not isinstance(values["require_same_phase"], bool):
        raise ValueError("require_same_phase must be boolean")
    return values


def _pair_config(config):
    """Compatibility alias for callers written before the public validator."""
    return validate_pair_config(config)


def _alignment_evidence(sequence, alignments, start0, end0):
    if alignments is None:
        return dict(available=False, status="not_provided", localization_confirmation=False,
                    physical_chromosome_end_confirmed=False)
    records = []
    for index, alignment in enumerate(alignments):
        if index >= MAX_AUXILIARY_RECORDS:
            raise ValueError("Too many supplied alignments")
        if not isinstance(alignment, dict) or not {"query_start0", "query_end0", "target", "strand", "mapq", "tp"} <= alignment.keys():
            raise ValueError("Alignment requires query_start0, query_end0, target, strand, mapq, and tp")
        a, b = alignment["query_start0"], alignment["query_end0"]
        _bounds(sequence, a, b)
        if a == b:
            raise ValueError("Alignment query interval must have positive length")
        _strand(alignment["strand"], "Alignment strand", optional=False)
        if not isinstance(alignment["target"], str) or not alignment["target"].strip():
            raise ValueError("Alignment target must be a nonempty string")
        _integer(alignment["mapq"], "Alignment MAPQ", maximum=255)
        if alignment["tp"] not in (None, "P", "S"):
            raise ValueError("Alignment tp must be P, S, or unspecified")
        row = dict(query_start0=a, query_end0=b, target=alignment["target"],
                   strand=alignment["strand"], mapq=alignment["mapq"], tp=alignment["tp"],
                   pair_overlap_bp=max(0, min(b, end0) - max(a, start0)),
                   covers_pair=a <= start0 and b >= end0,
                   mapq_available=alignment["mapq"] != 255)
        target_keys = {"target_start0", "target_end0", "target_length"}
        present_target_keys = target_keys.intersection(alignment)
        if present_target_keys and present_target_keys != target_keys:
            raise ValueError("Target coordinates require target_start0, target_end0, and target_length together")
        row["target_coordinates_available"] = bool(present_target_keys)
        if present_target_keys:
            target_length = _integer(alignment["target_length"], "target_length", minimum=1)
            target_start = _integer(alignment["target_start0"], "target_start0", maximum=target_length)
            target_end = _integer(alignment["target_end0"], "target_end0", minimum=target_start + 1, maximum=target_length)
            row.update(target_start0=target_start, target_end0=target_end, target_length=target_length)
        if "alignment_id" in alignment:
            if not isinstance(alignment["alignment_id"], str) or not alignment["alignment_id"].strip():
                raise ValueError("Supplied alignment_id must be a nonempty string")
            row["alignment_id"] = alignment["alignment_id"]
        records.append(row)
    unique_primary = {}
    for index, row in enumerate(records):
        if row["tp"] == "P":
            if row["target_coordinates_available"]:
                key = tuple(row[k] for k in ("query_start0", "query_end0", "target", "strand", "target_start0", "target_end0", "target_length"))
            else:
                # Identical query spans can map to different positions on the
                # same target. Missing positions cannot establish duplication.
                key = ("unpositioned_primary_row", index)
            unique_primary[key] = row
    primary = list(unique_primary.values())
    # Detect overlapping primary query spans in O(n log n), not pairwise O(n^2).
    competing, farthest_end = False, -1
    for row in sorted(primary, key=lambda r: (r["query_start0"], r["query_end0"])):
        competing |= row["query_start0"] < farthest_end
        farthest_end = max(farthest_end, row["query_end0"])
    split = len({(r["query_start0"], r["query_end0"]) for r in primary}) > 1
    if len(primary) > 1:
        status = "split_or_competing_primary_mapping"
    elif primary and not primary[0]["target_coordinates_available"]:
        status = "single_primary_query_span_only"
    elif primary and primary[0]["covers_pair"]:
        status = "single_primary_covers_pair"
    elif primary:
        status = "primary_mapping_does_not_cover_pair"
    else:
        status = "no_primary_mapping"
    return dict(available=True, status=status, alignments=records,
                distinct_primary_count=len(primary), duplicate_primary_rows=len([r for r in records if r["tp"] == "P"]) - len(primary),
                split_primary_mapping=split, competing_primary_mapping=competing,
                multiple_primary_targets=len({r["target"] for r in primary}) > 1,
                multiple_primary_strands=len({r["strand"] for r in primary}) > 1,
                primary_rows_without_target_coordinates=sum(not r["target_coordinates_available"] for r in primary),
                primary_target_positions_complete=all(r["target_coordinates_available"] for r in primary) if primary else None,
                primary_covers_pair_count=sum(r["covers_pair"] for r in primary),
                localization_confirmation=False, physical_chromosome_end_confirmed=False,
                interpretation="Supplied query and optional target spans only; missing target positions cannot establish duplicate mappings; no CIGAR continuity, reference endpoint, chimera, or structural-variant inference")


def _anchor_phase(sequence, anchor, motif):
    scores = fixed_phase_scores(sequence[anchor["start0"]:anchor["end0"]], anchor["start0"], motif)
    best = scores["best_by_strand"][anchor["strand"]]
    globally_best_strands = sorted({r["strand"] for r in scores["all_best_ties"]})
    return dict(tract_id=anchor["tract_id"], start0=anchor["start0"], end0=anchor["end0"],
                declared_strand=anchor["strand"], best_by_declared_strand=best,
                unique_best_phase=len(best) == 1,
                selected_absolute_phase=best[0]["absolute_template_phase"] if len(best) == 1 else None,
                best_identity=best[0]["identity"], globally_best_strands=globally_best_strands,
                template_strand_ambiguous=len(globally_best_strands) > 1,
                declared_strand_is_globally_best=anchor["strand"] in globally_best_strands,
                fixed_phase=scores)


def evaluate_pair(sequence, left, right, motif, config=None, alignments=None, expected_strand=None):
    """Evaluate two supplied tracts under a bounded sequence compatibility model.

The original anchors are never merged. Sequence compatibility and supplied
alignment observations are separate. Unknown expected direction leaves terminal
context unassessed, without itself rejecting sequence compatibility. An empty
gap has no denominator and bypasses gap score thresholds. Nonempty gaps use
exact runs wholly inside the gap: their cutoff does not borrow anchor bases.
With require_same_phase=False, differing anchor phases are permitted explicitly
and the gap's best same-strand phase is tested; no continuity inference follows.
"""
    motif = _motif(motif)
    values = _pair_config(config)
    _strand(expected_strand)
    left, right = _tract(sequence, left, "left"), _tract(sequence, right, "right")
    start0, end0 = min(left["start0"], right["start0"]), max(left["end0"], right["end0"])
    observed_gap = right["start0"] - left["end0"]
    reasons = []
    if left["start0"] > right["start0"]:
        reasons.append("anchors_out_of_order")
    if max(left["start0"], right["start0"]) < min(left["end0"], right["end0"]):
        reasons.append("anchors_overlap")
    same_strand = left["strand"] == right["strand"]
    if not same_strand:
        reasons.append("anchor_orientation_conflict")
    if expected_strand is not None and any(a["strand"] != expected_strand for a in (left, right)):
        reasons.append("expected_strand_conflict")
    direction_ambiguous = reverse_complement(motif) in motif + motif
    if direction_ambiguous:
        reasons.append("anchor_template_strand_ambiguous")
    if expected_strand is None:
        terminal_status = "unassessed_expected_strand_unknown"
    elif direction_ambiguous:
        terminal_status = "direction_ambiguous"
    elif "expected_strand_conflict" in reasons:
        terminal_status = "expected_strand_conflict"
    else:
        terminal_status = "expected_strand_consistent"
    if observed_gap > values["max_gap_bp"]:
        reasons.append("gap_exceeds_max_bp")
    if end0 - start0 > values["max_span_bp"]:
        reasons.append("span_exceeds_max_bp")
    result = dict(schema_version=SCHEMA_VERSION, module_version=VERSION,
                  model="bounded_local_repeat_sequence_compatibility", motif=motif, motif_period_bp=len(motif),
                  config=values, left=dict(left), right=dict(right), start0=start0, end0=end0,
                  length_bp=end0 - start0, observed_gap_bp=observed_gap,
                  record_context=_record_context(sequence, start0, end0),
                  terminal_context=dict(status=terminal_status, expected_strand=expected_strand,
                                        record_edge_is_not_terminal_evidence=True, physical_chromosome_end_confirmed=False),
                  reverse_complement_rotation_equivalent=direction_ambiguous,
                  alignment_evidence=_alignment_evidence(sequence, alignments, start0, end0),
                  anchor_phase_evidence=None, phase_compatibility=None, gap=None, bounding=None,
                  gap_same_strand_evidence=None, diagnostics_skipped=False,
                  biological_tracts_merged=False, changes_native_classification=False,
                  changes_policy_classification=False, physical_chromosome_end_confirmed=False,
                  continuity_inference_permitted=False,
                  limitations=["Compatibility is a local sequence hypothesis, not a merged biological array.",
                               "Fixed-phase mismatches do not identify their cause.",
                               "Exact gap coverage is censored at gap boundaries and can reject short perfect gaps.",
                               "Supplied alignment spans do not establish CIGAR continuity or physical terminal location.",
                               "No empirical sensitivity, false-positive rate, or absence claim is implied."])
    if "span_exceeds_max_bp" in reasons:
        result.update(status="rejected_by_bounded_model", reason_codes=list(dict.fromkeys(reasons)), diagnostics_skipped=True)
        return result
    anchor_evidence = {label: _anchor_phase(sequence, anchor, motif) for label, anchor in (("left", left), ("right", right))}
    result["anchor_phase_evidence"] = anchor_evidence
    for evidence in anchor_evidence.values():
        if not evidence["unique_best_phase"]:
            reasons.append("anchor_phase_ambiguous")
        if evidence["template_strand_ambiguous"]:
            reasons.append("anchor_template_strand_ambiguous")
        if not evidence["declared_strand_is_globally_best"]:
            reasons.append("anchor_declared_strand_not_best")
        if evidence["best_identity"] < values["min_anchor_fixed_phase_identity"]:
            reasons.append("anchor_identity_below_threshold")
    phases = [anchor_evidence[label]["selected_absolute_phase"] for label in ("left", "right")]
    same_phase = same_strand and None not in phases and phases[0] == phases[1]
    result["phase_compatibility"] = dict(same_strand=same_strand, left_absolute_phase=phases[0],
                                         right_absolute_phase=phases[1], same_absolute_phase=same_phase,
                                         phase_delta_mod_period=(phases[1] - phases[0]) % len(motif) if same_strand and None not in phases else None,
                                         require_same_phase=values["require_same_phase"])
    if values["require_same_phase"] and same_strand and None not in phases and not same_phase:
        reasons.append("anchor_absolute_phase_mismatch")
    cutoff = values["min_exact_copies"] * len(motif)
    cutoffs = tuple(sorted({len(motif), cutoff}))
    result["bounding"] = analyze_interval(sequence, start0, end0, motif, known_tracts=(left, right), expected_strand=expected_strand, cutoffs=cutoffs)
    if observed_gap >= 0:
        gap = analyze_interval(sequence, left["end0"], right["start0"], motif, expected_strand=expected_strand, cutoffs=cutoffs)
        result["gap"] = gap
        if observed_gap == 0:
            result["gap_same_strand_evidence"] = dict(empty_gap=True, strand=left["strand"] if same_strand else None,
                                                       minimum_exact_run_bp=cutoff, covered_fraction=None,
                                                       maximum_uncovered_run_bp=0, fixed_phase_identity=None,
                                                       score_model="empty_gap_thresholds_not_applicable")
        elif same_strand:
            coverage = gap["exact_coverage_views"][str(cutoff)]["by_strand"][left["strand"]]
            max_uncovered = max((x["length_bp"] for x in coverage["uncovered_islands"]), default=0)
            if same_phase:
                candidates = [s for s in gap["fixed_phase"]["scores"] if s["strand"] == left["strand"] and s["absolute_template_phase"] == phases[0]]
                score_model = "shared_anchor_absolute_phase"
            elif not values["require_same_phase"]:
                candidates = gap["fixed_phase"]["best_by_strand"][left["strand"]]
                score_model = "best_gap_phase_with_phase_switching_allowed_by_config"
            else:
                candidates = []
                score_model = "no_unique_shared_anchor_phase"
            identity = candidates[0]["identity"] if candidates else None
            result["gap_same_strand_evidence"] = dict(empty_gap=False, strand=left["strand"],
                                                       minimum_exact_run_bp=cutoff, covered_fraction=coverage["covered_fraction"],
                                                       covered_bp=coverage["covered_bp"], uncovered_bp=coverage["uncovered_bp"],
                                                       maximum_uncovered_run_bp=max_uncovered,
                                                       fixed_phase_identity=identity, fixed_phase_scores=candidates,
                                                       score_model=score_model, coverage=coverage)
            if coverage["covered_fraction"] < values["min_gap_coverage"]:
                reasons.append("gap_exact_coverage_below_threshold")
            if max_uncovered > values["max_uncovered_run_bp"]:
                reasons.append("gap_uncovered_run_exceeds_max_bp")
            if identity is not None and identity < values["min_gap_fixed_phase_identity"]:
                reasons.append("gap_fixed_phase_identity_below_threshold")
    result.update(status="rejected_by_bounded_model" if reasons else "compatible_under_bounded_model",
                  reason_codes=list(dict.fromkeys(reasons)))
    return result
