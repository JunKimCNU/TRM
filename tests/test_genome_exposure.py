"""Numeric interval enumeration only: no sequencing reads or sequence fixtures."""

import importlib.util
import itertools
from pathlib import Path
import unittest


SPEC = importlib.util.spec_from_file_location(
    "genome_exposure", Path(__file__).resolve().parents[1] / "scripts/genome_exposure.py")
geometry = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(geometry)


def enumerate_coordinates(chromosome_lengths, terminal_lengths, r1, r2, fragment, d1, d2):
    """Independent overlap-of-intervals oracle for complete-fragment starts."""
    counts = {"positive_reads": 0, "supporting_pairs": 0,
              "one_positive_mate_pairs": 0, "both_positive_mate_pairs": 0}
    legal_starts = 0
    for index, chromosome in enumerate(chromosome_lengths):
        left, right = terminal_lengths[2 * index:2 * index + 2]
        ends = ((0, left), (chromosome - right, chromosome))
        for start in range(chromosome - fragment + 1):
            legal_starts += 1
            intervals = ((start, start + r1), (start + fragment - r2, start + fragment))
            positive = []
            for interval, required in zip(intervals, (d1, d2)):
                overlaps = [max(0, min(interval[1], end[1]) - max(interval[0], end[0])) for end in ends]
                positive.append(any(overlap >= required for overlap in overlaps))
            count = sum(positive)
            counts["positive_reads"] += count
            counts["supporting_pairs"] += count > 0
            counts["one_positive_mate_pairs"] += count == 1
            counts["both_positive_mate_pairs"] += count == 2
    return legal_starts, counts


class GeometricExposureTests(unittest.TestCase):
    def test_exhaustive_small_integer_coordinate_cases(self):
        cases = 0
        for r1, r2 in itertools.product(range(1, 5), repeat=2):
            for fragment in range(max(r1, r2), max(r1, r2) + 3):
                for left, right in itertools.product(range(0, 5), repeat=2):
                    chromosome = left + right + fragment + 3
                    for d1, d2 in ((1, 1), (r1, r2), (r1 + 1, 1), (1, r2 + 1), (2, 3)):
                        denominator, expected = enumerate_coordinates([chromosome], [left, right], r1, r2, fragment, d1, d2)
                        result = geometry.paired_terminal_exposure(chromosome, 1, [left, right], r1, r2, fragment, (d1, d2), 17)
                        self.assertEqual(result["valid_fragment_starts"], denominator)
                        self.assertEqual(result["start_counts"], expected,
                                         (chromosome, left, right, r1, r2, fragment, d1, d2))
                        for key, count in expected.items():
                            self.assertAlmostEqual(result[f"expected_{key}"], count * 17 / denominator)
                        cases += 1
        self.assertEqual(cases, 6000)

    def test_multiple_unequal_chromosomes_and_terminal_lengths(self):
        lengths = [0, 3, 6, 2, 5, 1]
        chromosome_lengths = [20, 25, 30]
        denominator, counts = enumerate_coordinates(chromosome_lengths, lengths, 4, 3, 8, 2, 3)
        result = geometry.paired_terminal_exposure(sum(chromosome_lengths), 3, lengths, 4, 3, 8, (2, 3), 100)
        self.assertEqual(result["valid_fragment_starts"], denominator)
        self.assertEqual(result["start_counts"], counts)
        self.assertAlmostEqual(result["telomeric_base_fraction"], sum(lengths) / sum(chromosome_lengths))
        self.assertAlmostEqual(result["nominal_read_coverage"], 100 * 7 / 75)
        self.assertTrue(result["end_geometry"][0]["zero_length_end"])

    def test_discrete_inclusive_seed_boundary(self):
        result = geometry.paired_terminal_exposure(100, 1, [5, 4], 5, 5, 10, 5, 1)
        self.assertEqual(result["start_counts"]["supporting_pairs"], 1)
        self.assertEqual(result["end_geometry"][0]["r1_starts"], 1)
        self.assertEqual(result["end_geometry"][1]["r2_starts"], 0)

    def test_seed_longer_than_one_or_both_mates(self):
        result = geometry.paired_terminal_exposure(100, 1, [20, 20], 4, 7, 12, (5, 6), 100)
        self.assertEqual(sum(end["r1_starts"] for end in result["end_geometry"]), 0)
        self.assertGreater(result["expected_positive_reads"], 0)
        result = geometry.paired_terminal_exposure(100, 1, [20, 20], 4, 7, 12, 8, 100)
        self.assertEqual(result["expected_positive_reads"], 0)

    def test_fragment_length_changes_near_and_far_exposure(self):
        outputs = [geometry.paired_terminal_exposure(1000, 1, [30, 30], 10, 10, fragment, 5, 100)
                   for fragment in (10, 20, 40)]
        self.assertEqual([x["end_geometry"][0]["r1_starts"] for x in outputs], [26, 26, 26])
        self.assertEqual([x["end_geometry"][0]["r2_starts"] for x in outputs], [26, 16, 0])
        self.assertEqual([x["end_geometry"][1]["r2_starts"] for x in outputs], [26, 26, 26])
        self.assertEqual([x["end_geometry"][1]["r1_starts"] for x in outputs], [26, 16, 0])
        for result in outputs:
            starts = result["start_counts"]
            self.assertEqual(starts["positive_reads"], starts["one_positive_mate_pairs"] + 2 * starts["both_positive_mate_pairs"])
            self.assertEqual(starts["supporting_pairs"], starts["one_positive_mate_pairs"] + starts["both_positive_mate_pairs"])

    def test_mate_and_end_swap_symmetry(self):
        first = geometry.paired_terminal_exposure(100, 1, [17, 11], 8, 5, 12, (3, 4), 50)
        second = geometry.paired_terminal_exposure(100, 1, [11, 17], 5, 8, 12, (4, 3), 50)
        self.assertEqual(first["start_counts"], second["start_counts"])

    def test_ploidy_consistent_genome_and_end_scaling(self):
        first = geometry.paired_terminal_exposure(200, 2, [3, 4, 5, 6], 7, 6, 12, 3, 100)
        copies = geometry.paired_terminal_exposure(400, 4, [3, 4, 5, 6] * 2, 7, 6, 12, 3, 100)
        self.assertEqual(first["telomeric_base_fraction"], copies["telomeric_base_fraction"])
        self.assertEqual(copies["valid_fragment_starts"], 2 * first["valid_fragment_starts"])
        for key in first["start_counts"]:
            self.assertEqual(copies["start_counts"][key], 2 * first["start_counts"][key])
            self.assertEqual(first[f"expected_{key}"], copies[f"expected_{key}"])
        self.assertEqual(first["nominal_read_coverage"], 2 * copies["nominal_read_coverage"])
        doubled_reads = geometry.paired_terminal_exposure(400, 4, [3, 4, 5, 6] * 2, 7, 6, 12, 3, 200)
        self.assertEqual(first["nominal_read_coverage"], doubled_reads["nominal_read_coverage"])
        self.assertEqual(doubled_reads["expected_positive_reads"], 2 * first["expected_positive_reads"])

    def test_zero_length_ends_and_zero_input_pairs(self):
        result = geometry.paired_terminal_exposure(100, 2, [0, 0, 0, 0], 5, 5, 10, 3, 100)
        self.assertEqual(result["telomeric_base_fraction"], 0)
        self.assertEqual(result["expected_positive_reads"], 0)
        self.assertTrue(all(end["zero_length_end"] for end in result["end_geometry"]))
        zero_pairs = geometry.paired_terminal_exposure(100, 1, [10, 10], 5, 5, 10, 3, 0)
        self.assertEqual(zero_pairs["nominal_read_coverage"], 0)
        self.assertEqual(zero_pairs["expected_positive_reads"], 0)

    def test_full_repeat_threshold_integer_rounding(self):
        self.assertEqual(geometry.full_repeat_minimum_overlaps(151, 151, 28), (136, 136))
        self.assertEqual(geometry.full_repeat_minimum_overlaps(10, 11, 1), (9, 10))
        self.assertEqual(geometry.full_repeat_minimum_overlaps(10, 11, (12, 4)), (12, 10))
        for r in range(1, 501):
            threshold = geometry.full_repeat_minimum_overlaps(r, r, 1)[0]
            self.assertGreaterEqual(10 * threshold, 9 * r)
            self.assertLess(10 * (threshold - 1), 9 * r)

    def test_outputs_have_no_biological_decision_fields(self):
        result = geometry.paired_terminal_exposure(100, 1, [10, 10], 5, 5, 10, 3, 10)
        for forbidden in ("sensitivity", "power", "absence", "adequacy", "classification"):
            self.assertFalse(any(forbidden in key for key in result))
        self.assertFalse(result["assumptions"]["individual_chromosome_lengths_supplied"])


class GeometryValidationTests(unittest.TestCase):
    def arguments(self, **overrides):
        values = dict(genome_bp=100, chromosome_count=1, tract_lengths_bp=[10, 10],
                      read1_bp=5, read2_bp=5, fragment_bp=10, minimum_overlap_bp=3,
                      number_pairs=10, ends_nonoverlapping_assumed=True)
        values.update(overrides)
        return values

    def test_invalid_integer_inputs(self):
        for key in ("genome_bp", "chromosome_count", "read1_bp", "read2_bp", "fragment_bp", "minimum_overlap_bp"):
            for invalid in (-1, 0, 1.0, True, float("nan"), float("inf"), "5"):
                with self.subTest(key=key, invalid=invalid), self.assertRaises(ValueError):
                    geometry.paired_terminal_exposure(**self.arguments(**{key: invalid}))
        for invalid in (-1, True, 1.0, float("nan"), float("inf")):
            with self.subTest(number_pairs=invalid), self.assertRaises(ValueError):
                geometry.paired_terminal_exposure(**self.arguments(number_pairs=invalid))

    def test_invalid_lengths_overlaps_and_fragment(self):
        cases = ({"tract_lengths_bp": [10]}, {"tract_lengths_bp": [10, 10, 10]},
                 {"tract_lengths_bp": [-1, 10]}, {"tract_lengths_bp": [1.0, 10]},
                 {"tract_lengths_bp": [True, 10]}, {"tract_lengths_bp": None},
                 {"fragment_bp": 4}, {"minimum_overlap_bp": (1,)},
                 {"minimum_overlap_bp": (1, 2, 3)}, {"minimum_overlap_bp": (0, 1)},
                 {"minimum_overlap_bp": (1, float("nan"))})
        for values in cases:
            with self.subTest(values=values), self.assertRaises(ValueError):
                geometry.paired_terminal_exposure(**self.arguments(**values))

    def test_nonoverlap_assumption_must_be_explicitly_true(self):
        for assumption in (False, None, "unknown", 1):
            with self.subTest(assumption=assumption), self.assertRaisesRegex(ValueError, "assumption"):
                geometry.paired_terminal_exposure(**self.arguments(ends_nonoverlapping_assumed=assumption))

    def test_necessary_aggregate_size_condition(self):
        # 20 modeled terminal bases plus 9 excluded fragment starts.
        for genome in (1, 9, 28, 29):
            with self.subTest(genome=genome), self.assertRaisesRegex(ValueError, "Necessary nonoverlap"):
                geometry.paired_terminal_exposure(**self.arguments(genome_bp=genome))
        result = geometry.paired_terminal_exposure(**self.arguments(genome_bp=30))
        self.assertLessEqual(result["start_counts"]["supporting_pairs"], result["valid_fragment_starts"])

    def test_overflowing_derived_output_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "finite"):
            geometry.paired_terminal_exposure(**self.arguments(number_pairs=10**1000))


if __name__ == "__main__":
    unittest.main()
