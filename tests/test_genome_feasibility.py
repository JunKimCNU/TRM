"""Meaningful planning boundaries; no sequencing data are opened."""
from __future__ import annotations

import copy
import importlib.util
from pathlib import Path
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/assess_genome_feasibility.py"
SPEC = importlib.util.spec_from_file_location("assess_genome_feasibility", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FeasibilityTests(unittest.TestCase):
    def setUp(self):
        self.manifest = {"sample_id": "SRR1", "path1": "r1.fastq.gz", "path2": "r2.fastq.gz"}
        self.metadata = {"sample_id": "SRR1", "scientific_name": "Fixture", "tax_id": "1",
                         "lineage": "Eukaryota; Metazoa;", "instrument_platform": "ILLUMINA",
                         "library_strategy": "WGS", "library_source": "GENOMIC", "library_layout": "PAIRED"}
        self.library = {**self.metadata, "wga_mda_status": "not_reported", "library_selection": "RANDOM"}
        self.assessed = {"sample_id": "SRR1", "analysis_complete": True,
                         "input_summary": {"records": 10, "bases": 1500},
                         "data_adequacy": {"read_qc": {"records": 10, "bases": 1500,
                             "max_record_bases": 150, "sampled_prefix": False,
                             "length_histogram": [{"lower": 150, "upper": 150, "records": 10}]}}}

    def run_sample(self, flags=(), k=21):
        return MODULE.interpret_sample(self.manifest, self.metadata, self.library,
                                       self.assessed, list(flags), k, inspect_paths=False)

    def test_pilot_does_not_infer_genome_coverage_karyotype_or_no_wga(self):
        row = self.run_sample()
        self.assertEqual(row["genome_estimation_status"], "candidate_for_spectrum_pilot")
        self.assertIsNone(row["genome_size_bp"])
        self.assertIsNone(row["haploid_chromosome_count"])
        self.assertEqual(row["whole_genome_amplification_status"], "unknown")
        self.assertFalse(row["genome_estimation_performed"])
        self.assertFalse(row["local_pairing_verified"])
        self.assertEqual(row["eligible_records_lower_bound"], 10)

    def test_wgs_label_does_not_override_protocol_conflict(self):
        flags = [{"code": "rna_protocol_conflicts_with_wgs_label", "reason": "RNA PolyA"}]
        self.assertEqual(self.run_sample(flags)["genome_estimation_status"], "review_protocol_conflict")

    def test_failure_priority_preserves_protocol_reason(self):
        self.assessed["analysis_complete"] = False
        self.assessed["data_adequacy"] = {}
        flags = [{"code": "hic_protocol_conflicts_with_wgs_label", "reason": "Hi-C"}]
        row = self.run_sample(flags)
        self.assertEqual(row["genome_estimation_status"], "failed_input")
        self.assertIn("hic_protocol_conflicts_with_wgs_label", row["all_flags"])
        self.assertIsNone(row["records"])

    def test_no_pcr_does_not_resolve_wga(self):
        self.library["library_construction_protocol"] = "No PCR"
        row = self.run_sample()
        self.assertEqual(row["whole_genome_amplification_status"], "unknown")
        self.assertIn("not_evidence", row["library_pcr_status"])

    def test_one_file_is_not_assumed_missing_or_paired(self):
        self.manifest["path2"] = ""
        row = self.run_sample()
        self.assertEqual(row["genome_estimation_status"], "candidate_for_spectrum_pilot")
        self.assertIn("interleaved_unverified", row["local_layout_interpretation"])
        self.assertIn("archive_paired_local_one_file", row["all_flags"])

    def test_nonstandard_platform_not_hidden_by_read_volume(self):
        self.metadata["instrument_platform"] = "LS454"
        self.assertEqual(self.run_sample()["genome_estimation_status"], "nonstandard_platform")

    def test_qpcr_or_single_worm_is_not_automatic_wga(self):
        flags = [{"code": "sparse_qpcr_only_protocol", "reason": "qPCR only"},
                 {"code": "single_organism_or_cell_text", "reason": "single worm"}]
        row = self.run_sample(flags)
        self.assertEqual(row["genome_estimation_status"], "candidate_for_spectrum_pilot")
        self.assertEqual(row["whole_genome_amplification_status"], "unknown")
        self.assertEqual(len(row["all_flags"].split(";")), 2)

    def test_actual_wga_requires_amplification_review(self):
        self.library["wga_mda_status"] = "reported"
        row = self.run_sample()
        self.assertEqual(row["genome_estimation_status"], "review_amplification")
        self.assertEqual(row["whole_genome_amplification_status"], "reported")

    def test_metagenomic_and_negative_control_reasons_both_retained(self):
        self.metadata["library_source"] = self.library["library_source"] = "METAGENOMIC"
        flags = [{"code": "negative_control_description", "reason": "NegCon"}]
        row = self.run_sample(flags)
        self.assertEqual(row["genome_estimation_status"], "review_non_genomic_wgs")
        self.assertIn("negative_control_description", row["all_flags"])

    def test_bounded_histogram_does_not_pretend_exact_eligibility(self):
        self.assessed["data_adequacy"]["read_qc"]["length_histogram"] = [{"lower": 10, "upper": 200, "records": 10}]
        row = self.run_sample()
        self.assertEqual(row["eligible_records_lower_bound"], 0)
        self.assertEqual(row["eligible_records_upper_bound"], 10)

    def test_all_reads_shorter_than_k_requires_qc_review(self):
        self.assertEqual(self.run_sample(k=151)["genome_estimation_status"], "review_input_qc")

    def test_sample_identity_and_histogram_inconsistency_rejected(self):
        bad = copy.deepcopy(self.assessed)
        self.assessed["sample_id"] = "different"
        with self.assertRaisesRegex(ValueError, "identity"):
            self.run_sample()
        self.assessed = bad
        self.assessed["data_adequacy"]["read_qc"]["length_histogram"][0]["records"] = 9
        with self.assertRaisesRegex(ValueError, "histogram count"):
            self.run_sample()

    def test_metadata_disagreement_is_not_silently_overwritten(self):
        self.library["library_strategy"] = "RNA-Seq"
        with self.assertRaisesRegex(ValueError, "disagree"):
            self.run_sample()


if __name__ == "__main__":
    unittest.main()
