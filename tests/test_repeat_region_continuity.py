"""Algorithm and boundary tests for the optional periodic continuity scorer."""

import importlib.util
import itertools
from pathlib import Path
import random
import unittest


SPEC = importlib.util.spec_from_file_location(
    "repeat_region_continuity", Path(__file__).resolve().parents[1] / "scripts/repeat_region_continuity.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def oracle(sequence, motif):
    """Independent ordinary Levenshtein DP against each unrolled start phase.

    A query of n bases can always cost <=n by insertion only, so a better target
    cannot contain >2n bases. Standard global DP consumes the whole query; its
    last row covers every permitted finite target length, including zero.
    """
    size = len(sequence)
    best = size
    for phase in range(len(motif)):
        template = "".join(motif[(phase + i) % len(motif)] for i in range(2 * size))
        previous = list(range(len(template) + 1))
        for i, base in enumerate(sequence, 1):
            current = [i]
            for j, target_base in enumerate(template, 1):
                current.append(min(previous[j] + 1, current[j - 1] + 1,
                                   previous[j - 1] + (base != target_base)))
            previous = current
        best = min(best, min(previous))
    return best


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class PeriodicScoreTests(unittest.TestCase):
    def test_exhaustive_binary_queries_and_motifs_against_unrolled_global_oracle(self):
        queries = ["".join(x) for size in range(6) for x in itertools.product("AC", repeat=size)]
        motifs = ["".join(x) for size in range(1, 5) for x in itertools.product("AC", repeat=size)]
        for sequence in queries:
            for motif in motifs:
                with self.subTest(sequence=sequence, motif=motif):
                    self.assertEqual(M.minimum_periodic_edit_score(sequence, motif)["edit_cost"],
                                     oracle(sequence, motif))

    def test_random_four_base_queries_against_unrolled_global_oracle(self):
        rng = random.Random(938121)
        for _ in range(400):
            motif = "".join(rng.choices("ACGT", k=rng.randrange(1, 9)))
            sequence = "".join(rng.choices("ACGT", k=rng.randrange(0, 18)))
            self.assertEqual(M.minimum_periodic_edit_score(sequence, motif)["edit_cost"], oracle(sequence, motif))

    def test_rotations_and_reverse_complement_preserve_distance(self):
        sequence = "TTAGGGTTAAGGGTTTGGGTTAGCCCTTAGGG"
        motif = "TTAGGG"
        expected = M.minimum_periodic_edit_score(sequence, motif)["edit_cost"]
        for offset in range(len(motif)):
            rotated = motif[offset:] + motif[:offset]
            self.assertEqual(M.minimum_periodic_edit_score(sequence, rotated)["edit_cost"], expected)
            self.assertEqual(M.minimum_periodic_edit_score(reverse_complement(sequence), reverse_complement(rotated))["edit_cost"], expected)

    def test_empty_query_and_phase_definitions(self):
        self.assertEqual(M.minimum_periodic_edit_score("", "ACG"), {
            "edit_cost": 0, "query_bp": 0, "edit_cost_per_query_base": None, "end_phases": [0, 1, 2]})
        self.assertEqual(M.minimum_periodic_edit_score("CGAC", "ACG")["end_phases"], [2])

    def test_paid_indels_and_no_free_internal_phase_reset(self):
        self.assertEqual(M.minimum_periodic_edit_score("ACGTACGGTACGT", "ACGT")["edit_cost"], 1)
        self.assertEqual(M.minimum_periodic_edit_score("ACGTAGTACGT", "ACGT")["edit_cost"], 1)
        self.assertEqual(M.minimum_periodic_edit_score("ACACCCACAC", "AC")["edit_cost"], 1)
        # Every letter belongs to the motif, but free rephasing at each query
        # position would incorrectly score this string zero.
        self.assertGreater(M.minimum_periodic_edit_score("AAAACCCCGGGGTTTT", "ACGT")["edit_cost"], 0)

    def test_no_query_clipping(self):
        self.assertEqual(M.minimum_periodic_edit_score("CCCCAAAAACCCC", "A")["edit_cost"], 8)

    def test_strict_alphabet_and_nonempty_motif(self):
        for sequence, motif in [("N", "A"), ("acgt", "ACGT"), ("A", ""), ("A", "N"), ("A", "acgt"), (None, "A")]:
            with self.assertRaises(ValueError):
                M.minimum_periodic_edit_score(sequence, motif)


class RegionTests(unittest.TestCase):
    def test_clean_region_and_no_gap_denominators(self):
        row = M.evaluate_region("TTAGGG" * 100, "TTAGGG")
        self.assertTrue(row["passes"])
        self.assertTrue(row["evaluated"])
        self.assertEqual(row["global"]["edit_cost"], 0)
        self.assertEqual(row["failed_reason_codes"], [])
        self.assertEqual(len(row["windows"]), 21)
        self.assertIsNone(row["max_gap_edit_ratio"])

    def test_nonrepeat_island_cannot_be_diluted_by_long_flanks(self):
        sequence = "TTAGGG" * 300 + "C" * 100 + "TTAGGG" * 300
        row = M.evaluate_region(sequence, "TTAGGG", [(1800, 1900)])
        self.assertTrue(row["global"]["passes"])
        self.assertFalse(row["passes"])
        self.assertIn("window_edit_ratio_exceeded", row["failed_reason_codes"])
        self.assertIn("gap_edit_ratio_exceeded", row["failed_reason_codes"])
        self.assertIn("gap_context_edit_ratio_exceeded", row["failed_reason_codes"])

    def test_all_rows_scored_even_after_global_failure(self):
        row = M.evaluate_region("C" * 600, "A", [(0, 24), (100, 124), (550, 600)])
        self.assertFalse(row["global"]["passes"])
        self.assertEqual(len(row["windows"]), 21)
        self.assertEqual(len(row["gaps"]), 3)
        self.assertEqual(row["failed_reason_codes"], ["global_edit_ratio_exceeded", "window_edit_ratio_exceeded", "gap_edit_ratio_exceeded", "gap_context_edit_ratio_exceeded"])

    def test_local_gap_guard_detects_short_island(self):
        # 24 nonrepeat bases pass every 100 bp window at <=25%, but fail the
        # >=24 bp standalone and flanking-context gap criteria.
        row = M.evaluate_region("A" * 240 + "C" * 24 + "A" * 336, "A", [(240, 264)])
        self.assertTrue(row["global"]["passes"])
        self.assertTrue(all(window["passes"] for window in row["windows"]))
        self.assertFalse(row["gaps"][0]["passes"])
        self.assertFalse(row["gaps"][0]["context"]["passes"])

    def test_gap_context_is_independent_of_standalone_score(self):
        # Eight changes in 24 gap bases pass <=35%, while adjacent changes make
        # the 48 bp context exceed25%; all other sequence is still scored.
        sequence = "A" * 228 + "C" * 12 + "C" * 8 + "A" * 16 + "C" * 12 + "A" * 324
        row = M.evaluate_region(sequence, "A", [(240, 264)])
        self.assertTrue(row["gaps"][0]["passes"])
        self.assertFalse(row["gaps"][0]["context"]["passes"])

    def test_short_gaps_scored_but_not_guarded_empty_ignored(self):
        row = M.evaluate_region("A" * 240 + "C" * 23 + "A" * 337, "A", [(0, 0), (240, 263), (600, 600)])
        self.assertEqual(len(row["gaps"]), 1)
        self.assertFalse(row["gaps"][0]["guard_applied"])
        self.assertIsNone(row["gaps"][0]["passes"])
        self.assertIsNone(row["gaps"][0]["context"]["passes"])
        self.assertEqual(row["max_gap_edit_ratio"], 1.0)
        self.assertIsNone(row["max_guarded_gap_edit_ratio"])
        self.assertTrue(row["passes"])

    def test_boundary_context_clipped(self):
        row = M.evaluate_region("A" * 600, "A", [(0, 24), (576, 600)])
        contexts = [gap["context"] for gap in row["gaps"]]
        self.assertEqual([(r["start0"], r["end0"]) for r in contexts], [(0, 36), (564, 600)])

    def test_gaps_reject_invalid_order_overlap_and_coordinates(self):
        for gaps in [[(20, 50), (30, 60)], [(60, 70), (10, 20)], [(-1, 10)], [(0, 601)], [(10, 9)], [(False, 24)], [(0, 24.0)], [(0,)], ["xx"]]:
            with self.subTest(gaps=gaps), self.assertRaises(ValueError):
                M.evaluate_region("A" * 600, "A", gaps)

    def test_ambiguous_is_unscored_not_failure(self):
        row = M.evaluate_region("A" * 599 + "N", "A")
        self.assertEqual(row["status"], "not_evaluated_ambiguous")
        self.assertFalse(row["evaluated"])
        self.assertIsNone(row["passes"])
        self.assertIsNone(row["global"])
        self.assertEqual(row["failed_reason_codes"], [])
        self.assertEqual(row["invalid_sequence_symbols"], ["N"])

    def test_bounds_and_non_evaluation_reason_codes(self):
        for size, status in [(0, "outside_length_scope_too_short"), (499, "outside_length_scope_too_short"), (100_001, "not_evaluated_length_limit")]:
            row = M.evaluate_region("A" * size, "A")
            self.assertEqual(row["status"], status)
            self.assertFalse(row["evaluated"])
            self.assertIsNone(row["passes"])
        for size in (500, 100_000):
            self.assertTrue(M.evaluate_region("A" * size, "A")["passes"])
        row = M.evaluate_region("N", "A")
        self.assertEqual(row["not_evaluated_reason_codes"], ["ambiguous_or_non_acgt_sequence", "region_below_minimum_scope"])

    def test_window_tail_short_region_and_coverage(self):
        self.assertEqual(M.window_intervals(80), [(0, 80)])
        self.assertEqual(M.window_intervals(0), [])
        windows = M.window_intervals(511)
        self.assertEqual(windows[-2:], [(400, 500), (411, 511)])
        self.assertEqual(len(windows), len(set(windows)))
        self.assertEqual(set().union(*(set(range(a, b)) for a, b in windows)), set(range(511)))
        row = M.evaluate_region("A" * 485 + "C" * 26, "A")
        self.assertFalse(row["windows"][-1]["passes"])
        self.assertTrue(row["windows"][-2]["passes"])

    def test_exact_decimal_threshold_inclusive_and_config_fixed(self):
        self.assertTrue(M._passes({"edit_cost": 3, "query_bp": 20}, 0.15))
        self.assertFalse(M._passes({"edit_cost": 4, "query_bp": 20}, 0.15))
        M.evaluate_region("A" * 500, "A", config=dict(M.DEFAULT_CONFIG))
        for config in [{"window_bp": 101}, {"window_max_edit_ratio": 0.30}, {"unknown": 1}, {"window_bp": True}]:
            with self.assertRaises(ValueError):
                M.evaluate_region("A" * 500, "A", config=config)


if __name__ == "__main__":
    unittest.main()
