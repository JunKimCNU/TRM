"""Length applicability is an additive necessary-condition diagnostic only."""

import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from input_applicability import build_input_applicability


ZERO = "no_length_eligible_reads"
PRESENT = "length_eligible_reads_present_power_unknown"
UNRESOLVED = "length_eligibility_unresolved"


def report(length=150, records=10, **config_changes):
    config = {"mode": "reads", "min_k": 3, "max_k": 32,
              "seed_bases": 24, "min_copies": 3}
    config.update(config_changes)
    return {"config": config,
            "summary": {"records": records, "bases": records * length,
                        "max_record_bases": length if records else 0,
                        "sampled_prefix": False},
            "classification": "unchanged", "adequacy": {"status": "unknown"},
            "candidates": [{"motif": "TTAGGG", "reads": 5}],
            "rejected_candidates": [{"motif": "AAT", "reads": 100}]}


def period(result, value):
    return next(row for row in result["periods"] if row["period_bp"] == value)


class InputApplicabilityTests(unittest.TestCase):
    def test_seed_rounding_boundary_24_is_not_25(self):
        low = build_input_applicability(report(24, min_k=5, max_k=5))
        high = build_input_applicability(report(25, min_k=5, max_k=5))
        self.assertEqual(period(low, 5)["minimum_seed_bp"], 25)
        self.assertEqual(low["status"], ZERO)
        self.assertEqual(period(low, 5)["length_eligible_reads_upper"], 0)
        self.assertEqual(high["status"], PRESENT)
        self.assertEqual(period(high, 5)["length_eligible_reads_lower"], 10)
        self.assertEqual(period(high, 5)["length_eligible_reads_upper"], 10)

    def test_minimum_copies_boundary_and_max_k_64(self):
        for length, expected in ((95, ZERO), (96, PRESENT)):
            with self.subTest(length=length):
                result = build_input_applicability(report(length, min_k=32, max_k=32))
                self.assertEqual(result["status"], expected)
                self.assertEqual(period(result, 32)["minimum_seed_bp"], 96)
        result = build_input_applicability(report(192, min_k=64, max_k=64))
        self.assertEqual(len(result["periods"]), 1)
        self.assertEqual(period(result, 64)["minimum_seed_bp"], 192)
        self.assertEqual(result["status"], PRESENT)

    def test_whole_search_zero_requires_every_period_zero(self):
        short = build_input_applicability(report(23))
        self.assertEqual(short["status"], ZERO)
        self.assertTrue(all(row["status"] == ZERO for row in short["periods"]))
        mixed = build_input_applicability(report(24))
        self.assertEqual(mixed["status"], PRESENT)
        self.assertEqual(period(mixed, 3)["status"], PRESENT)
        self.assertEqual(period(mixed, 32)["status"], ZERO)

    def test_family_zero_does_not_make_whole_search_zero(self):
        result = build_input_applicability(report(40), ["TGTGTGGGTGTGGT"])
        self.assertEqual(result["status"], PRESENT)
        self.assertEqual(result["families"][0]["period_bp"], 14)
        self.assertEqual(result["families"][0]["minimum_seed_bp"], 42)
        self.assertEqual(result["families"][0]["status"], ZERO)

    def test_family_opportunity_does_not_overrule_impossible_discovery_range(self):
        result = build_input_applicability(report(30, min_k=32, max_k=32), ["TTAGGG"])
        self.assertEqual(result["status"], ZERO)
        row = result["families"][0]
        self.assertEqual(row["status"], PRESENT)
        self.assertTrue(row["outside_discovery_range"])
        self.assertEqual(row["minimum_seed_bp"], 24)

    def test_primitive_rotation_and_reverse_complement_deduplicate(self):
        motifs = ["TTAGGG", "CCCTAA", "TAGGGT", "TTAGGGTTAGGG", "ttaggg"]
        saved = copy.deepcopy(motifs)
        result = build_input_applicability(report(), motifs)
        self.assertEqual(motifs, saved)
        self.assertEqual(len(result["families"]), 1)
        family = result["families"][0]
        self.assertEqual(family["canonical_family"], "AACCCT")
        self.assertEqual(family["period_bp"], 6)
        self.assertEqual(family["expected_motifs"], sorted(set(m.upper() for m in motifs)))
        self.assertFalse(family["outside_discovery_range"])

    def test_exact_histogram_counts_and_overlapping_period_denominators(self):
        value = report()
        value["evidence_qc"] = {"records": 10, "bases": 1500, "max_record_bases": 150,
                                "length_histogram": [{"lower": 150, "upper": 150, "records": 10}]}
        result = build_input_applicability(value, ["TTAGGG"])
        self.assertEqual(result["records"], 10)
        self.assertTrue(result["length_distribution_available"])
        self.assertTrue(all(row["length_eligible_reads_lower"] == 10 for row in result["periods"]))
        self.assertEqual(result["families"][0]["bounds_method"], "histogram_and_aggregate_length_bounds")

    def coarse_report(self):
        value = report(min_k=19, max_k=19, min_copies=5)
        value["summary"].update(records=10, bases=525, max_record_bases=100)
        value["evidence_qc"] = {"length_histogram": [{"lower": 10, "upper": 10, "records": 5},
                                                     {"lower": 90, "upper": 100, "records": 5}]}
        return value

    def test_coarse_histogram_does_not_invent_exact_lengths(self):
        result = build_input_applicability(self.coarse_report())
        row = result["periods"][0]
        self.assertEqual(row["minimum_seed_bp"], 95)
        self.assertEqual((row["length_eligible_reads_lower"], row["length_eligible_reads_upper"]), (0, 5))
        self.assertEqual(result["status"], UNRESOLVED)

    def test_histogram_can_prove_eligible_count_despite_wide_aggregate_bounds(self):
        value = self.coarse_report()
        value["config"].update(min_k=16, max_k=16)
        result = build_input_applicability(value)
        row = result["periods"][0]
        self.assertEqual((row["length_eligible_reads_lower"], row["length_eligible_reads_upper"]), (5, 5))
        self.assertEqual(result["status"], PRESENT)

    def test_legacy_no_histogram_uses_aggregate_bounds(self):
        value = self.coarse_report()
        del value["evidence_qc"]
        result = build_input_applicability(value)
        self.assertFalse(result["length_distribution_available"])
        self.assertEqual(result["status"], UNRESOLVED)
        self.assertEqual(result["periods"][0]["bounds_method"], "aggregate_length_bound")
        self.assertEqual(result["periods"][0]["length_eligible_reads_upper"], 5)

    def test_prefix_zero_is_scoped_not_whole_input_zero(self):
        value = report(20)
        value["summary"]["sampled_prefix"] = True
        result = build_input_applicability(value, ["TTAGGG"])
        self.assertEqual(result["status"], "analysis_incomplete")
        self.assertEqual(result["scope"], "analyzed_records")
        for row in result["periods"] + result["families"]:
            self.assertEqual(row["status"], "analysis_incomplete")
            self.assertEqual(row["scope"], "analyzed_records")
            self.assertEqual(row["analyzed_records_length_status"], ZERO)
            self.assertEqual(row["length_eligible_reads_upper"], 0)

    def test_caller_and_saved_failure_gate_overrule_full_prefix_flag(self):
        for key in ("analysis_status", "status", "state"):
            value = report()
            value[key] = "failed"
            result = build_input_applicability(value)
            self.assertEqual(result["status"], "analysis_incomplete")
        result = build_input_applicability(report(), analysis_complete=False)
        self.assertEqual(result["status"], "analysis_incomplete")
        self.assertEqual(result["periods"][0]["analyzed_records_length_status"], PRESENT)

    def test_missing_or_non_boolean_completeness_is_unverified(self):
        for prefix in (None, "false", 0, 1, [], {}):
            with self.subTest(prefix=prefix):
                value = report()
                value["summary"]["sampled_prefix"] = prefix
                result = build_input_applicability(value)
                self.assertEqual(result["status"], "input_completeness_unverified")
                self.assertEqual(result["scope"], "unverified_input")
                self.assertEqual(result["periods"][0]["status"], "input_completeness_unverified")
        value = report()
        del value["summary"]["sampled_prefix"]
        self.assertEqual(build_input_applicability(value)["status"], "input_completeness_unverified")

    def test_missing_config_has_no_invented_defaults(self):
        for key in ("mode", "min_k", "max_k", "seed_bases", "min_copies"):
            with self.subTest(key=key):
                value = report()
                del value["config"][key]
                result = build_input_applicability(value, ["TTAGGG"])
                self.assertEqual(result["status"], "config_unavailable")
                self.assertEqual(result["periods"], [])
                self.assertIsNone(result["families"][0]["minimum_seed_bp"])
        value = report()
        del value["config"]
        self.assertEqual(build_input_applicability(value)["status"], "config_unavailable")

    def test_invalid_native_settings_are_rejected_and_work_is_bounded(self):
        changes = [{"min_k": 0}, {"max_k": 65}, {"max_k": 10**12}, {"min_k": 33},
                   {"seed_bases": 0}, {"seed_bases": 100_000_001}, {"min_copies": 1},
                   {"min_copies": 1_000_001}, {"min_copies": True}, {"max_k": 32.0},
                   {"seed_bases": "24"}, {"mode": None}]
        for change in changes:
            with self.subTest(change=change):
                result = build_input_applicability(report(**change))
                self.assertEqual(result["status"], "config_invalid")
                self.assertEqual(result["periods"], [])
        value = report()
        value["config"] = []
        self.assertEqual(build_input_applicability(value)["status"], "config_invalid")

    def test_large_valid_seed_needs_no_seed_string_allocation(self):
        result = build_input_applicability(report(seed_bases=100_000_000, min_copies=1_000_000))
        self.assertEqual(result["status"], ZERO)
        self.assertEqual(period(result, 3)["minimum_seed_bp"], 100_000_002)

    def test_missing_and_malformed_qc_remain_unavailable_not_zero(self):
        for key in ("records", "bases", "max_record_bases"):
            value = report()
            del value["summary"][key]
            result = build_input_applicability(value, ["TTAGGG"])
            self.assertEqual(result["status"], "qc_invalid")
            self.assertEqual(result["scope"], "unavailable")
            self.assertIsNone(result["families"][0]["length_eligible_reads_upper"])
        for malformed in (None, [], {"records": True}, {"bases": 1500.0},
                          {"length_histogram": [None]},
                          {"length_histogram": [{"lower": 150, "upper": 150, "records": 9}]}):
            value = report()
            value["evidence_qc"] = malformed
            self.assertEqual(build_input_applicability(value)["status"], "qc_invalid")

    def test_contradictory_totals_or_overlapping_histogram_are_invalid(self):
        for changes in ({"records": True}, {"bases": 1600}, {"max_record_bases": 2000}):
            value = report()
            value["summary"].update(changes)
            self.assertEqual(build_input_applicability(value)["status"], "qc_invalid")
        value = report()
        value["evidence_qc"] = {"length_histogram": [{"lower": 140, "upper": 150, "records": 5},
                                                     {"lower": 150, "upper": 160, "records": 5}]}
        self.assertEqual(build_input_applicability(value)["status"], "qc_invalid")

    def test_source_missing_and_non_read_mode_never_infer_zeros(self):
        for value in (None, [], "failed"):
            self.assertEqual(build_input_applicability(value)["status"], "source_unavailable")
        result = build_input_applicability(report(), ["TTAGGG"], source_available=False)
        self.assertEqual(result["families"][0]["status"], "source_unavailable")
        result = build_input_applicability(report(mode="assembly"), ["TTAGGG"])
        self.assertEqual(result["status"], "mode_not_supported")
        self.assertIsNone(result["families"][0]["length_eligible_reads_upper"])

    def test_empty_completed_input_is_length_zero_not_a_biological_decision(self):
        result = build_input_applicability(report(records=0), ["TTAGGG"])
        self.assertEqual(result["records"], 0)
        self.assertEqual(result["status"], ZERO)
        self.assertEqual(result["scope"], "complete_input")
        self.assertNotIn("classification", result)
        self.assertNotIn("adequacy", result)

    def test_invalid_expected_references_are_visible_and_do_not_erase_periods(self):
        result = build_input_applicability(report(), ["TTAGGG", "TG2-3(TG)1-6", "NNNN", "", 42, "A" * 65])
        self.assertEqual(result["status"], PRESENT)
        self.assertEqual(len(result["families"]), 1)
        self.assertEqual(len(result["invalid_expected_motifs"]), 5)
        self.assertEqual(len(result["periods"]), 30)
        self.assertEqual(build_input_applicability(report(), "TTAGGG")["families"][0]["period_bp"], 6)
        result = build_input_applicability(report(), None)
        self.assertEqual(result["status"], PRESENT)
        self.assertEqual(len(result["invalid_expected_motifs"]), 1)

    def test_input_and_all_existing_policy_observations_are_unchanged(self):
        value = self.coarse_report()
        before = copy.deepcopy(value)
        first = build_input_applicability(value, ["CCCTAA", "TTAGGG"])
        second = build_input_applicability(value, ["TTAGGG", "CCCTAA"])
        self.assertEqual(value, before)
        self.assertEqual(first, second)
        first["config"]["min_k"] = 1
        self.assertEqual(value, before)


if __name__ == "__main__":
    unittest.main()
