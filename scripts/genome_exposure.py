"""Pure numeric exposure under a fixed-fragment, uniform-start terminal model.

Inputs describe hypothetical geometry, not fitted sequencing-library properties.
This module reads no sequence files and makes no detector/classifier decisions.
Chromosome lengths are not supplied: the caller must explicitly assume that
each chromosome accommodates the fragment and that its two terminal supporting
start windows do not overlap. The aggregate genome-size check is necessary,
but cannot establish those chromosome-specific assumptions.
"""

import math


def _integer(value, name, minimum=1):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _overlaps(value):
    if isinstance(value, (tuple, list)):
        if len(value) != 2:
            raise ValueError("minimum_overlap_bp must be an integer or a pair of integers")
        return (_integer(value[0], "R1 minimum overlap"),
                _integer(value[1], "R2 minimum overlap"))
    overlap = _integer(value, "minimum_overlap_bp")
    return overlap, overlap


def full_repeat_minimum_overlaps(read1_bp, read2_bp, minimum_seed_overlap_bp):
    """Require both the seed threshold and ceil(90% of that mate's length).

    Pass the actual exact-seed length, including any whole-copy rounding. For
    a 7 bp motif, seed_bases=24 and min_copies=3 require 28 seed bases.
    This overlap threshold is a geometric approximation to the corresponding
    detector condition; it does not model errors, interruptions or seed phase.
    """
    r1 = _integer(read1_bp, "read1_bp")
    r2 = _integer(read2_bp, "read2_bp")
    d1, d2 = _overlaps(minimum_seed_overlap_bp)
    return max(d1, (9 * r1 + 9) // 10), max(d2, (9 * r2 + 9) // 10)


def paired_terminal_exposure(genome_bp, chromosome_count, tract_lengths_bp,
                             read1_bp, read2_bp, fragment_bp,
                             minimum_overlap_bp, number_pairs,
                             ends_nonoverlapping_assumed=True):
    """Expected positive records/pairs under declared terminal geometry.

    ``tract_lengths_bp`` contains exactly 2*C nonnegative integer lengths in
    chromosome left/right order. A zero length is an end with no modeled motif.
    ``minimum_overlap_bp`` is a positive integer or separate (R1, R2) thresholds;
    a threshold greater than its mate's length gives zero starts for that mate.
    ``number_pairs`` counts sequenced pairs, without molecule deduplication.

    Fragments are length F, wholly contained in one chromosome. R1 covers the
    first R1 bases and R2 the last R2 bases. Among D=G-C*(F-1) legal integer
    starts, all positions are equally likely. At a left end the near R1 width
    is max(L-d1+1,0), and the far R2 width max(L-(F-R2)-d2+1,0).
    The right end swaps near/far mates. Each width requires read_length>=d.
    Terminal windows are nested at a given end, so union=max, both=min.

    The output describes expected geometric exposure only. No genome coverage
    fit, telomere retention estimate or biological interpretation is performed.
    """
    genome = _integer(genome_bp, "genome_bp")
    chromosomes = _integer(chromosome_count, "chromosome_count")
    r1 = _integer(read1_bp, "read1_bp")
    r2 = _integer(read2_bp, "read2_bp")
    fragment = _integer(fragment_bp, "fragment_bp")
    pairs = _integer(number_pairs, "number_pairs", minimum=0)
    if fragment < max(r1, r2):
        raise ValueError("fragment_bp must be >= both mate lengths")
    if ends_nonoverlapping_assumed is not True:
        raise ValueError("Explicit nonoverlapping terminal-start windows assumption must be True")
    if not isinstance(tract_lengths_bp, (list, tuple)) or len(tract_lengths_bp) != 2 * chromosomes:
        raise ValueError("tract_lengths_bp must contain exactly 2*chromosome_count lengths in left/right order")
    lengths = [_integer(length, f"tract_lengths_bp[{index}]", minimum=0)
               for index, length in enumerate(tract_lengths_bp)]
    d1, d2 = _overlaps(minimum_overlap_bp)
    telomeric_bases = sum(lengths)
    excluded_starts = chromosomes * (fragment - 1)
    if genome <= telomeric_bases + excluded_starts:
        raise ValueError("Necessary nonoverlap condition requires genome_bp > sum(tract_lengths_bp) + C*(F-1)")
    denominator = genome - excluded_starts

    def width(length, read, required, offset):
        return max(length - offset - required + 1, 0) if read >= required else 0

    ends = []
    totals = {"positive_reads": 0, "supporting_pairs": 0,
              "one_positive_mate_pairs": 0, "both_positive_mate_pairs": 0}
    for index, length in enumerate(lengths):
        left = index % 2 == 0
        starts1 = width(length, r1, d1, 0 if left else fragment - r1)
        starts2 = width(length, r2, d2, fragment - r2 if left else 0)
        supporting = max(starts1, starts2)
        both = min(starts1, starts2)
        one = abs(starts1 - starts2)
        reads = starts1 + starts2
        ends.append({"chromosome_index": index // 2 + 1,
                     "side": "left" if left else "right", "tract_bp": length,
                     "zero_length_end": length == 0,
                     "r1_starts": starts1, "r2_starts": starts2,
                     "supporting_pair_starts": supporting, "both_mate_starts": both,
                     "one_mate_starts": one, "positive_read_starts": reads})
        totals["positive_reads"] += reads
        totals["supporting_pairs"] += supporting
        totals["one_positive_mate_pairs"] += one
        totals["both_positive_mate_pairs"] += both
    if totals["supporting_pairs"] > denominator:
        raise ValueError("Summed supporting-start probability exceeds 1 under supplied geometry")
    try:
        values = {"telomeric_base_fraction": telomeric_bases / genome,
                  "nominal_read_coverage": pairs * (r1 + r2) / genome,
                  **{f"expected_{key}": count * pairs / denominator
                     for key, count in totals.items()}}
    except OverflowError as exc:
        raise ValueError("Derived geometric exposure must be representable as finite numbers") from exc
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError("Derived geometric exposure must be finite")
    return {
        "schema_version": 1,
        "model": "fixed_fragment_uniform_start_terminal_exposure_v1",
        "inputs": {"genome_bp": genome, "chromosome_count": chromosomes,
                   "tract_lengths_bp": lengths, "read1_bp": r1, "read2_bp": r2,
                   "fragment_bp": fragment, "minimum_overlap_bp": [d1, d2],
                   "number_pairs": pairs},
        "assumptions": {
            "uniform_valid_fragment_starts": True, "fixed_fragment_length": True,
            "ends_nonoverlapping_assumed": True,
            "individual_chromosome_lengths_supplied": False,
            "chromosome_specific_fit_and_nonoverlap": "Assumed, not checked; only the necessary aggregate size condition is checked.",
            "zero_length_end": "An end with no modeled motif, not a statement about telomere maintenance.",
            "pair_scope": "Sequenced record pairs, not deduplicated independent molecules.",
            "scope": "Declared geometric exposure; library recovery and detector behavior are not fitted.",
        },
        "valid_fragment_starts": denominator, "start_counts": totals,
        "end_geometry": ends, **values,
    }
