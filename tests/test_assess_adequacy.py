import importlib.util
import itertools
import json
import math
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location("assess_adequacy", Path(__file__).parents[1] / "scripts/assess_adequacy.py")
adequacy = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(adequacy)


def report(records=10_000_000, length=150):
    return {"config": {"mode": "reads", "min_k": 3, "max_k": 32, "min_copies": 3,
                       "seed_bases": 24, "min_reads": 2},
            "summary": {"records": records, "bases": records * length,
                        "max_record_bases": length if records else 0,
                        "sampled_prefix": False, "call": "no_call"},
            "evidence_qc": {"records_with_ambiguity": 0, "ambiguous_bases": 0},
            "candidates": [], "rejected_candidates": []}


def assumptions(**changes):
    value = {"genome_size_bp": 100_000_000, "haploid_chromosome_count": 6,
             "random_wgs_assumed": True, "telomere_retention_fraction": 0.1,
             "effective_independent_read_fraction": 0.5, "detector_sensitivity": 0.5,
             "provenance": {key: "explicit test hypothesis, not an observed species parameter"
                            for key in adequacy.REQUIRED_ASSUMPTIONS}}
    value.update(changes)
    return value


class AdequacyTests(unittest.TestCase):
    def test_actual_seed_rounding_not_unrounded_maximum(self):
        self.assertEqual(adequacy.required_seed_bases(5, 24, 3), 25)
        self.assertEqual(adequacy.required_seed_bases(7, 24, 3), 28)
        self.assertEqual(adequacy.required_seed_bases(32, 24, 3), 96)

    def test_length_bounds_hold_for_exhaustively_enumerated_small_data(self):
        for lengths in itertools.product(range(1, 6), repeat=4):
            for seed in range(1, 7):
                bound = adequacy.eligible_read_bounds(4, sum(lengths), max(lengths), seed)
                actual = sum(length >= seed for length in lengths)
                self.assertLessEqual(bound["lower"], actual)
                self.assertGreaterEqual(bound["upper"], actual)

    def test_histogram_improves_bounds_without_guessing_within_intervals(self):
        hist = [{"lower": 10, "upper": 10, "records": 5},
                {"lower": 90, "upper": 100, "records": 5}]
        bound = adequacy.eligible_read_bounds(10, 525, 100, 95, hist)
        self.assertLessEqual(bound["lower"], 5)
        self.assertEqual(bound["upper"], 5)
        exact = adequacy.eligible_read_bounds(10, 525, 100, 80, hist)
        self.assertEqual((exact["lower"], exact["upper"]), (5, 5))

    def test_exact_binomial_matches_hand_calculation(self):
        result = adequacy.binomial_miss_bound(10, 0.2, 2)
        self.assertAlmostEqual(result["miss_probability_upper_bound"], 0.8 ** 10 + 10 * 0.2 * 0.8 ** 9)
        self.assertEqual(adequacy.binomial_miss_bound(10, 0, 2)["miss_probability_upper_bound"], 1)
        self.assertEqual(adequacy.binomial_miss_bound(10, 1, 2)["miss_probability_upper_bound"], 0)
        self.assertEqual(adequacy.binomial_miss_bound(1, 0.5, 2)["miss_probability_upper_bound"], 1)

    def test_binomial_underflow_and_large_thresholds_are_labelled(self):
        underflow = adequacy.binomial_miss_bound(10**9, 0.01, 2)
        self.assertTrue(underflow["floating_point_underflow"])
        self.assertTrue(math.isfinite(underflow["log_miss_probability_upper_bound"]))
        bound = adequacy.binomial_miss_bound(10**7, 0.01, 10001)
        self.assertEqual(bound["method"], "binomial_chernoff_upper_bound")
        self.assertLessEqual(bound["miss_probability_upper_bound"], 1)

    def test_huge_read_count_cannot_supply_missing_genome_or_library_assumptions(self):
        value = adequacy.assess_adequacy(report(10**10), {"library_strategy": "WGS", "library_source": "GENOMIC"})
        self.assertEqual(value["status"], "unassessable_library_or_genome")
        self.assertIsNone(value["power"])
        self.assertFalse(value["library"]["amplification_absence_verified"])
        self.assertEqual(len(value["illustrative_grid"]["scenarios"]), 18)
        self.assertFalse(value["true_absence_inferred"])
        self.assertFalse(value["alternative_maintenance_inferred"])

    def test_power_supported_only_with_sourced_explicit_assumptions(self):
        value = adequacy.assess_adequacy(report(), assumptions=assumptions(), include_illustrative_grid=False)
        self.assertEqual(value["status"], "conditional_power_supported")
        self.assertLess(value["power"]["worst_case_miss_probability_upper_bound"], 0.05)
        self.assertEqual(value["power"]["mean_sequence_depth_given_genome_size"], 15)
        unsupported = adequacy.assess_adequacy(report(), assumptions=assumptions(provenance={}), include_illustrative_grid=False)
        self.assertEqual(unsupported["status"], "unassessable_library_or_genome")
        self.assertIsNotNone(unsupported["power"])

    def test_molecule_dependence_and_retention_reduce_power(self):
        high = adequacy.assess_adequacy(report(), assumptions=assumptions(), include_illustrative_grid=False)
        low = adequacy.assess_adequacy(report(), assumptions=assumptions(effective_independent_read_fraction=0.001), include_illustrative_grid=False)
        self.assertEqual(low["status"], "insufficient_data")
        self.assertGreater(low["power"]["worst_case_miss_probability_upper_bound"], high["power"]["worst_case_miss_probability_upper_bound"])
        zero = adequacy.assess_adequacy(report(), assumptions=assumptions(telomere_retention_fraction=0), include_illustrative_grid=False)
        self.assertEqual(zero["power"]["worst_case_miss_probability_upper_bound"], 1)

    def test_75bp_reads_do_not_support_96bp_exact_seed(self):
        value = adequacy.assess_adequacy(report(length=75), assumptions=assumptions(), include_illustrative_grid=False)
        self.assertEqual(value["status"], "insufficient_data")
        row = next(row for row in value["power"]["periods"] if row["period"] == 32)
        self.assertEqual(row["eligible_read_count_bounds"]["upper"], 0)
        self.assertEqual(row["miss_probability_upper_bound"], 1)

    def test_terminal_start_window_matches_small_enumerated_chromosome(self):
        chromosome, tract, read, seed = 100, 20, 15, 10
        starts = range(chromosome - read + 1)
        count = sum(min(start + read, tract) - start >= seed for start in starts)
        self.assertEqual(count, tract - seed + 1)
        config = {"min_k": 5, "max_k": 5, "seed_bases": 10, "min_copies": 2, "min_reads": 1}
        qc = adequacy._read_qc(report(100, read))
        power = adequacy.assessing_power(qc, config, assumptions(genome_size_bp=chromosome, haploid_chromosome_count=1, assumed_tract_bases=tract))
        self.assertEqual(power["periods"][0]["capture_start_positions_per_end"], count)

    def test_explicit_wga_rna_targeted_metadata_block_uniform_model(self):
        for metadata in ({"library_strategy": "WGA"}, {"amplification_method": "MDA"},
                         {"library_selection": "MDA"},
                         {"library_strategy": "RNA-Seq"}, {"library_strategy": "WXS"},
                         {"whole_genome_amplification": True}, {"metagenome": "true"}):
            with self.subTest(metadata=metadata):
                value = adequacy.assess_adequacy(report(), metadata, assumptions(), False)
                self.assertEqual(value["status"], "unassessable_library_or_genome")
                self.assertTrue(value["library"]["uniform_wgs_model_blockers"])

    def test_negated_wga_mention_is_review_not_positive_amplification_claim(self):
        value = adequacy.library_assessment({"library_construction_protocol": "Prepared without WGA"})
        self.assertIn("amplification_protocol_requires_manual_review", value["uniform_wgs_model_blockers"])
        self.assertNotIn("explicit_whole_genome_amplification", value["uniform_wgs_model_blockers"])

    def test_missing_prefix_completeness_and_actual_prefix_are_unassessable(self):
        for sampled in (None, True):
            data = report()
            data["summary"]["sampled_prefix"] = sampled
            value = adequacy.assess_adequacy(data, assumptions=assumptions(), include_illustrative_grid=False)
            self.assertEqual(value["status"], "unassessable_library_or_genome")

    def test_technical_failure_and_empty_input_are_distinguished(self):
        self.assertEqual(adequacy.assess_adequacy(None)["status"], "observed_qc_failure")
        self.assertEqual(adequacy.assess_adequacy({"state": "failed"})["status"], "observed_qc_failure")
        self.assertEqual(adequacy.assess_adequacy(report(0))["status"], "insufficient_data")
        self.assertEqual(adequacy.assess_adequacy(report(1))["status"], "insufficient_data")

    def test_histogram_mismatch_is_observed_qc_failure(self):
        data = report(10, 100)
        data["evidence_qc"] = {"length_histogram": [{"lower": 100, "upper": 100, "records": 9}]}
        value = adequacy.assess_adequacy(data)
        self.assertEqual(value["status"], "observed_qc_failure")
        self.assertIn("histogram", value["qc_error"])

    def test_ambiguous_read_records_subtracted_conservatively(self):
        data = report(100, 100)
        data["evidence_qc"] = {"records": 100, "bases": 10000, "max_record_bases": 100,
                               "records_with_ambiguity": 80, "ambiguous_bases": 80,
                               "length_histogram": [{"lower": 100, "upper": 100, "records": 100}]}
        value = adequacy.assess_adequacy(data, assumptions=assumptions(), include_illustrative_grid=False)
        row = value["power"]["periods"][0]
        self.assertEqual(row["ambiguity_free_eligible_reads_lower_bound"], 20)
        self.assertEqual(row["effective_independent_trials_assumed_lower_bound"], 10)

    def test_contradictory_ambiguity_qc_is_not_accepted(self):
        data = report(100, 100)
        data["evidence_qc"] = {"records_with_ambiguity": 80, "ambiguous_bases": 10}
        self.assertEqual(adequacy.assess_adequacy(data)["status"], "observed_qc_failure")

    def test_assembly_and_long_period_search_are_not_invalid_biological_data(self):
        data = report()
        data["config"]["mode"] = "assembly"
        self.assertEqual(adequacy.assess_adequacy(data)["status"], "unassessable_library_or_genome")
        data["config"].update(mode="reads", min_k=33, max_k=64)
        self.assertEqual(adequacy.assess_adequacy(data)["status"], "unassessable_library_or_genome")

    def test_invalid_assumption_values_raise_instead_of_silently_clamping(self):
        for changes in ({"genome_size_bp": 0}, {"detector_sensitivity": float("nan")},
                        {"telomere_retention_fraction": 2}, {"random_wgs_assumed": "true"},
                        {"haploid_chromosome_count": True}, {"max_miss_probability": 1}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                adequacy.assess_adequacy(report(), assumptions=assumptions(**changes))

    def test_json_output_has_no_nan_or_infinity_and_report_not_mutated(self):
        data = report()
        original = json.dumps(data, sort_keys=True)
        value = adequacy.assess_adequacy(data, assumptions=assumptions())
        json.dumps(value, allow_nan=False)
        self.assertEqual(json.dumps(data, sort_keys=True), original)
        self.assertFalse(value["motif_nomination_performed"])


if __name__ == "__main__":
    unittest.main()
