"""GenomeScope output can exist after failed or invalid model fitting."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    "report_genome_feasibility", Path(__file__).resolve().parents[1] / "scripts/report_genome_feasibility.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class GenomeReportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        (self.path / "summary.txt").write_text(
            "Genome Haploid Length         120,000,000 bp    130,000,000 bp\n"
            "Genome Repeat Length          10,000,000 bp     20,000,000 bp\n"
            "Model Fit                     12.1%             96.5%\n"
            "Read Error Rate               0.18%             0.20%\n")
        (self.path / "model.txt").write_text("Parameters:\n         Estimate Std. Error\nlength 1.2e+08 1.0\nkmercov 3.2e+01 0.1\n")
        (self.path / "progress.txt").write_text("starting\nround 0 converged\ndone\n")
        (self.path / "fit.execution.json").write_text(json.dumps({"returncode": 0}))

    def test_completed_model_is_exploratory_not_trusted_genome(self):
        row = MODULE.inspect_fit(self.path)
        self.assertEqual(row["fit_status"], "exploratory_fit")
        self.assertEqual(row["genome_haploid_length_min_bp"], 120000000)
        self.assertEqual(row["model_kmercov"], 32)
        self.assertEqual(row["all_kmer_fit_score_percent"], 12.1)
        self.assertEqual(row["full_model_fit_score_percent"], 96.5)

    def test_summary_after_failed_fit_is_rejected(self):
        (self.path / "progress.txt").write_text("starting\nfail\n")
        row = MODULE.inspect_fit(self.path)
        self.assertEqual(row["fit_status"], "rejected_model_output")
        self.assertIn("missing_success_marker", row["fit_guard_reasons"])

    def test_negative_model_parameter_rejected_despite_positive_summary(self):
        (self.path / "model.txt").write_text("length -120000000 1\nkmercov 30 1\n")
        row = MODULE.inspect_fit(self.path)
        self.assertEqual(row["fit_status"], "rejected_model_output")
        self.assertIn("nonpositive_or_missing_model_length", row["fit_guard_reasons"])

    def test_unrecorded_process_not_assumed_finished(self):
        (self.path / "fit.execution.json").unlink()
        self.assertEqual(MODULE.inspect_fit(self.path)["fit_status"], "pending_or_unrecorded")

    def test_nonzero_exit_rejects_otherwise_valid_output(self):
        (self.path / "fit.execution.json").write_text(json.dumps({"returncode": 1}))
        self.assertEqual(MODULE.inspect_fit(self.path)["fit_status"], "rejected_model_output")

    def test_nonfinite_negative_or_reversed_lengths_fail(self):
        for value in ("NaN bp 100 bp", "-1 bp 100 bp", "120 bp 100 bp", "Inf bp Inf bp"):
            with self.subTest(value=value):
                (self.path / "summary.txt").write_text("Genome Haploid Length " + value + "\n")
                self.assertEqual(MODULE.inspect_fit(self.path)["fit_status"], "rejected_model_output")

    def test_fit_scores_can_be_negative_and_are_not_ordered_ci(self):
        result = MODULE.summary_interval("Model Fit 98% -2%", "Model Fit", ordered=False)
        self.assertEqual(result, (98, -2))

    def test_duplicate_summary_or_parameter_is_not_silently_selected(self):
        self.assertIsNone(MODULE.summary_interval("Genome Haploid Length 1 bp 2 bp\nGenome Haploid Length 3 bp 4 bp", "Genome Haploid Length"))
        self.assertIsNone(MODULE.model_parameter("length 1 0\nlength 2 0", "length"))

    def test_saturated_counter_tail_recovered_from_exact_total(self):
        histogram = {"counted_kmer_occurrences_lower_bound": 1000,
                     "occurrences_above_fit_cutoff": {"10": 200}, "counter_cap": 100,
                     "distinct_kmers": 100}
        stats = {"#Total no. of k-mers": 1500, "#k-mers_below_min_threshold": 0,
                 "#k-mers_above_max_threshold": 0, "#Unique_counted_k-mers": 100}
        row = MODULE.exact_tail_audit(histogram, stats, 10)
        self.assertTrue(row["exact_tail_guard_passed"])
        self.assertEqual(row["exact_occurrences_above_fit_cutoff"], 700)
        self.assertAlmostEqual(row["exact_occurrence_fraction_above_fit_cutoff"], 700 / 1500)
        for key, value in (("#k-mers_below_min_threshold", 1), ("#Unique_counted_k-mers", 99), ("#Total no. of k-mers", 999)):
            changed = {**stats, key: value}
            self.assertFalse(MODULE.exact_tail_audit(histogram, changed, 10)["exact_tail_guard_passed"])
        self.assertFalse(MODULE.exact_tail_audit(histogram, stats, 100)["exact_tail_guard_passed"])

    def test_failed_summary_values_excluded_from_across_model_range(self):
        row = {"sample_id": "SRR1", "scientific_name": "species",
               **MODULE.inspect_fit(self.path)}
        bad = {**row, "fit_status": "rejected_model_output", "genome_haploid_length_max_bp": 9e20}
        ranges = MODULE.sample_model_ranges([row, bad])
        self.assertEqual(ranges[0]["maximum_over_models_bp"], 130000000)
        self.assertEqual(ranges[0]["exploratory_models"], 1)
        self.assertTrue(ranges[0]["review_required"])


if __name__ == "__main__":
    unittest.main()
