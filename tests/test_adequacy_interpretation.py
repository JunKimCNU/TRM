"""Scientific boundaries between observed repeats and conditional nondetection."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from assess_adequacy import REQUIRED_ASSUMPTIONS, assess_adequacy
from evidence_batch import canonical
from trm_assessment import assess_candidate, catalog_index, classify


def report(records=10_000_000, length=150):
    return {"config": {"mode": "reads", "min_k": 3, "max_k": 32,
                       "min_copies": 3, "seed_bases": 24, "min_reads": 2},
            "summary": {"records": records, "bases": records * length,
                        "max_record_bases": length, "sampled_prefix": False},
            # Power fixtures explicitly supply clean-read QC; absence is unknown.
            "evidence_qc": {"records_with_ambiguity": 0, "ambiguous_bases": 0}}


def assumptions():
    return {"genome_size_bp": 100_000_000, "haploid_chromosome_count": 6,
            "random_wgs_assumed": True, "telomere_retention_fraction": 0.1,
            "effective_independent_read_fraction": 0.5, "detector_sensitivity": 0.5,
            "provenance": {key: "Explicit test hypothesis; not a measurement"
                           for key in REQUIRED_ASSUMPTIONS}}


def observed(motif="TTAGGG", known=True, junctions=True):
    raw = {"motif": motif, "reads": 100, "full_repeat_reads": 60,
           "repeat_bases": 12000, "max_tract": 150}
    refs = catalog_index([{"motif": motif, "evidence": "experimental", "source": "fixture"}]) if known else {}
    target = {"motif": canonical(motif), "reads": 100, "repeat_bases": 12000,
              "junction_evidence": {
                  "left_forward_junction_reads": 12 if junctions else 0,
                  "left_reverse_junction_reads": 0,
                  "right_forward_junction_reads": 0,
                  "right_reverse_junction_reads": 12 if junctions else 0,
                  "junction_reads": 24 if junctions else 0,
                  "full_repeat_records": 60}}
    return assess_candidate(raw, target, refs, set(), 150)


class AdequacyInterpretationTests(unittest.TestCase):
    def test_missing_ambiguity_qc_cannot_support_conditional_nondetection(self):
        missing = report()
        del missing['evidence_qc']
        value = assess_adequacy(missing, assumptions=assumptions(), include_illustrative_grid=False)
        self.assertEqual(value['status'], 'unassessable_library_or_genome')
        self.assertIn('ambiguity_qc_unavailable', value['reasons'])
        self.assertEqual(classify([], value)[0], 'detection_inconclusive')
        self.assertFalse(value['true_absence_inferred'])

    def test_all_grid_rows_can_have_high_power_without_upgrading_actual_sample(self):
        value = assess_adequacy(report(10**10))
        self.assertEqual(value["status"], "unassessable_library_or_genome")
        self.assertIsNone(value["power"])
        self.assertTrue(all(row["worst_case_miss_probability_upper_bound"] < 0.05
                            for row in value["illustrative_grid"]["scenarios"]))
        self.assertEqual(classify([], value)[0], "detection_inconclusive")
        self.assertFalse(value["true_absence_inferred"])

    def test_each_required_assumption_matters_even_when_all_other_values_are_supplied(self):
        for omitted in REQUIRED_ASSUMPTIONS:
            with self.subTest(omitted=omitted):
                supplied = assumptions()
                del supplied[omitted]
                value = assess_adequacy(report(), assumptions=supplied, include_illustrative_grid=False)
                self.assertEqual(value["status"], "unassessable_library_or_genome")
                self.assertIn(omitted, value["missing_assumptions"])
                self.assertIsNone(value["power"])

    def test_all_end_power_does_not_generalize_to_motif_present_at_one_end(self):
        all_ends = assess_adequacy(report(), assumptions=assumptions(), include_illustrative_grid=False)
        one_end = assess_adequacy(report(), assumptions={**assumptions(), "motif_bearing_chromosome_ends": 1},
                                  include_illustrative_grid=False)
        self.assertEqual(all_ends["status"], "conditional_power_supported")
        self.assertEqual(one_end["status"], "insufficient_data")
        self.assertGreater(one_end["power"]["worst_case_miss_probability_upper_bound"],
                           all_ends["power"]["worst_case_miss_probability_upper_bound"])
        self.assertFalse(one_end["true_absence_inferred"])

    def test_a_numeric_power_calculation_is_not_an_adequacy_upgrade(self):
        cases = [({"library_strategy": "WGA"}, assumptions()),
                 ({}, {**assumptions(), "random_wgs_assumed": False}),
                 ({}, {**assumptions(), "provenance": {}})]
        for metadata, supplied in cases:
            with self.subTest(metadata=metadata, supplied=supplied):
                value = assess_adequacy(report(), metadata, supplied)
                self.assertIsNotNone(value["power"])
                self.assertTrue(value["power"]["all_scanned_periods_up_to_32_meet_conditional_target"])
                self.assertEqual(value["status"], "unassessable_library_or_genome")
                self.assertEqual(classify([], value)[0], "detection_inconclusive")

    def test_observed_known_and_novel_patterns_survive_missing_power_information(self):
        adequacy = assess_adequacy(report())
        self.assertEqual(adequacy["status"], "unassessable_library_or_genome")
        known = observed()
        full_only = observed(junctions=False)
        novel = observed("ACGATCG", known=False)
        self.assertEqual(classify([known], adequacy)[0], "known_trm_detected")
        self.assertEqual(classify([full_only], adequacy)[0], "known_repeat_detected_boundary_unresolved")
        self.assertEqual(classify([novel], adequacy)[0], "novel_geometry_nomination")
        self.assertTrue(novel['geometry_nomination'])
        self.assertFalse(novel['novel_candidate'])
        self.assertFalse(novel['confirmation_readiness']['validation_measured'])
        self.assertFalse(adequacy["motif_nomination_performed"])

    def test_wga_limits_nondetection_model_without_erasing_observed_known_pattern(self):
        adequacy = assess_adequacy(report(), {"amplification_method": "MDA"}, assumptions())
        self.assertEqual(adequacy["status"], "unassessable_library_or_genome")
        self.assertEqual(classify([observed()], adequacy)[0], "known_trm_detected")
        self.assertFalse(adequacy["alternative_maintenance_inferred"])

    def test_actual_qc_failure_and_non_eukaryotic_scope_override_positive_interpretation(self):
        inconsistent = report()
        inconsistent["summary"]["bases"] += 1
        failed = assess_adequacy(inconsistent)
        adequate = assess_adequacy(report(), assumptions=assumptions())
        for candidate in (observed(), observed("ACGATCG", known=False)):
            with self.subTest(candidate=candidate["candidate_status"]):
                self.assertEqual(classify([candidate], failed)[0], "detection_inconclusive_technical")
                self.assertEqual(classify([candidate], adequate, biological_scope=False)[0], "taxonomic_scope_unresolved")
                self.assertEqual(classify([candidate], adequate, analysis_complete=False)[0], "analysis_failed")

    def test_grid_is_optional_and_cannot_change_classification_or_supplied_power(self):
        for supplied, metadata in ((None, {}), (assumptions(), {}), (assumptions(), {"library_strategy": "WGA"})):
            with self.subTest(metadata=metadata, supplied=supplied):
                without = assess_adequacy(report(), metadata, supplied, include_illustrative_grid=False)
                with_grid = assess_adequacy(report(), metadata, supplied, include_illustrative_grid=True)
                copy_with = copy.deepcopy(with_grid)
                del copy_with["illustrative_grid"]
                self.assertEqual(without, copy_with)
                self.assertEqual(classify([], without), classify([], with_grid))


if __name__ == "__main__":
    unittest.main()
