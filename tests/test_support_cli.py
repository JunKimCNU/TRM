"""Integration checks for explicit negative-control support and saved evidence."""

import copy
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

from test_reclassify_batch import candidate, canonical
from test_cli import fastq


from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]
CONFIG = {"min_k": 3, "max_k": 32, "seed_bases": 24, "min_copies": 3,
          "min_tract": 500, "min_identity": 0.85, "end_window": 50,
          "min_reads": 2, "max_candidates": 100000, "mode": "reads"}
DERIVED_ROW_KEYS = {"telomere_prior", "support", "candidate_eligible"}


class SupportCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binary = runtime_binary()
        if not cls.binary.is_file():
            raise RuntimeError(f"Build the release executable first: {cls.binary}")

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="trm-support-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.calls = 0

    def write_json(self, name, value):
        path = self.root / name
        path.write_text(json.dumps(value, indent=2) + "\n")
        return path

    def row(self, motif, reads, repeat_bases=None):
        result = candidate(motif, 100000.125)
        result.update(reads=reads, repeat_bases=reads * 96 if repeat_bases is None else repeat_bases,
                      max_tract=96, forward_reads=reads // 2, reverse_reads=reads - reads // 2,
                      left_end_reads=reads, right_end_reads=reads, full_repeat_reads=reads,
                      left_forward_reads=reads // 2, right_forward_reads=reads // 2,
                      left_reverse_reads=reads - reads // 2, right_reverse_reads=reads - reads // 2)
        result["examples"] = [{"identity": 0.9807692307692307, "record_id": "saved-example"}]
        return result

    def report(self, rows, config=None):
        return {"schema_version": 1, "program": "trm", "version": "0.2.0",
                "threads": 1, "inputs": [{"path": "/unavailable/original.fastq.gz", "bytes": 1234}],
                "config": copy.deepcopy(config or CONFIG),
                "summary": {"records": 100000, "bases": 10000000, "max_record_bases": 100,
                            "candidate_families_before_filter": len(rows),
                            "reported_candidates": len(rows), "sampled_prefix": False,
                            "elapsed_seconds": 1.2345678901234567,
                            "megabases_per_second": 8.1000000729, "call": "candidates_only"},
                "warnings": ["original scan warning"], "candidates": rows}

    def model(self, motifs=None, name="model.json", controls=None):
        config = copy.deepcopy(CONFIG)
        config["min_reads"] = 1
        control = {"control_id": "negative", "source_result": "/archived/control/result.json",
                   "source_result_sha256": "0" * 64, "config": config,
                   "records": 100000, "bases": 10000000, "max_record_bases": 100,
                   "sampled_prefix": False,
                   "motifs": {canonical(motif): {"reads": values[0], "repeat_bases": values[1]}
                              for motif, values in (motifs or {}).items()}}
        return self.write_json(name, {"schema_version": 1,
                                     "model_type": "empirical_repeat_background", "name": "test",
                                     "created_utc": "2026-09-20T00:00:00+00:00", "limitations": [],
                                     "controls": controls if controls is not None else [control]})

    def reclassify(self, report, model=None, *options, success=True):
        self.calls += 1
        source = report if isinstance(report, Path) else self.write_json(f"source{self.calls}.json", report)
        original_bytes = source.read_bytes()
        prefix = self.root / f"output{self.calls}"
        command = [str(self.binary), "reclassify", str(source), "--output", str(prefix),
                   "--prior-policy", "telomere", "--support-policy", "filter"]
        if model is not None:
            command.extend(["--background-model", str(model)])
        command.extend(options)
        completed = subprocess.run(command, capture_output=True, text=True, timeout=15, check=False)
        self.assertEqual(source.read_bytes(), original_bytes, "Reclassification changed its source")
        if not success:
            self.assertNotEqual(completed.returncode, 0, completed.stdout)
            self.assertTrue(completed.stderr.strip())
            self.assertFalse(prefix.with_suffix(".json").exists())
            return None, prefix, completed
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return json.loads(prefix.with_suffix(".json").read_text()), prefix, completed

    def all_rows(self, result):
        return {row["motif"]: row for row in result["candidates"] + result.get("rejected_candidates", [])}

    def test_abundant_repeat_matching_background_is_rejected(self):
        motif = "TTAGGG"
        result, _, _ = self.reclassify(self.report([self.row(motif, 10000)]),
                                       self.model({motif: (10000, 960000)}))
        row = self.all_rows(result)[canonical(motif)]
        self.assertTrue(row["telomere_prior"]["eligible"])
        self.assertFalse(row["support"]["eligible"])
        self.assertFalse(row["candidate_eligible"])
        self.assertTrue(row["support"]["reasons"])
        self.assertEqual(result["summary"]["sequence_prior_candidates"], 1)
        self.assertEqual(result["summary"]["joint_eligible_candidates"], 0)
        self.assertEqual(result["candidates"], [])
        self.assertEqual(result["summary"]["call"], "no_call")

    def test_strong_known_repeat_passes_with_unseen_background_motif(self):
        result, _, _ = self.reclassify(self.report([self.row("TTAGGG", 1000)]), self.model())
        row = result["candidates"][0]
        self.assertEqual(row["motif"], canonical("TTAGGG"))
        self.assertTrue(row["support"]["eligible"])
        self.assertTrue(row["candidate_eligible"])
        self.assertGreaterEqual(row["support"]["read_enrichment_lower"], 2)
        self.assertGreaterEqual(row["support"]["repeat_base_enrichment"], 2)
        self.assertEqual(result["summary"]["reported_candidates"], 1)
        self.assertEqual(result["summary"]["eligible_telomere_candidates"], 1)

    def test_few_known_motif_reads_do_not_automatically_pass(self):
        result, _, _ = self.reclassify(self.report([self.row("TTAGGG", 2)]),
                                       self.model({"TTAGGG": (20, 1920)}))
        self.assertTrue(result["rejected_candidates"][0]["telomere_prior"]["eligible"])
        self.assertFalse(result["rejected_candidates"][0]["support"]["eligible"])
        self.assertEqual(result["summary"]["eligible_telomere_candidates"], 0)

    def test_unseen_weak_motif_reports_limited_control_resolution(self):
        result, _, _ = self.reclassify(self.report([self.row("TTAGGG", 2)]), self.model())
        support = result["rejected_candidates"][0]["support"]
        self.assertFalse(support["eligible"])
        self.assertEqual(support["status"], "insufficient_background_resolution")
        self.assertIn("calibration_depth_limits_unseen_motif", support["reasons"])

    def test_two_reads_can_pass_when_the_sample_fraction_is_strong(self):
        report = self.report([self.row("TTAGGG", 2)])
        report["summary"].update(records=10, bases=1000)
        result, _, _ = self.reclassify(report, self.model())
        self.assertTrue(result["candidates"][0]["support"]["eligible"])
        self.assertEqual(result["candidates"][0]["reads"], 2)

    def test_read_enrichment_without_repeat_base_enrichment_fails(self):
        row = self.row("TTAGGG", 5000, repeat_bases=120000)
        row.update(max_tract=24, full_repeat_reads=0)
        result, _, _ = self.reclassify(self.report([row]), self.model({"TTAGGG": (1000, 96000)}))
        support = result["rejected_candidates"][0]["support"]
        self.assertGreaterEqual(support["read_enrichment_lower"], 2)
        self.assertLess(support["repeat_base_enrichment"], 2)
        self.assertFalse(support["eligible"])
        self.assertIn("repeat_base_enrichment_below_threshold", support["reasons"])
        self.assertNotIn("read_fraction_enrichment_below_threshold", support["reasons"])

    def test_maximum_control_background_is_used_instead_of_average(self):
        empty = json.loads(self.model().read_text())["controls"][0]
        high = copy.deepcopy(empty)
        high["control_id"] = "high_background"
        high["motifs"] = {canonical("TTAGGG"): {"reads": 10000, "repeat_bases": 960000}}
        model = self.model(name="two_controls.json", controls=[empty, high])
        result, _, _ = self.reclassify(self.report([self.row("TTAGGG", 15000)]), model)
        self.assertFalse(result["rejected_candidates"][0]["support"]["eligible"])

    def test_strong_unsupported_motif_does_not_bypass_prior(self):
        result, _, _ = self.reclassify(self.report([self.row("AAT", 1000)]), self.model())
        row = result["rejected_candidates"][0]
        self.assertTrue(row["support"]["eligible"])
        self.assertFalse(row["telomere_prior"]["eligible"])
        self.assertFalse(row["candidate_eligible"])
        self.assertEqual(result["candidates"], [])

    def test_related_motif_can_allow_aat_without_blacklist(self):
        result, _, _ = self.reclassify(self.report([self.row("AAT", 1000)]), self.model(),
                                       "--related-motif", "AAT")
        row = result["candidates"][0]
        self.assertEqual(row["motif"], canonical("AAT"))
        self.assertEqual(row["telomere_prior"]["class"], "related_exact")
        self.assertTrue(row["support"]["eligible"])
        self.assertTrue(row["candidate_eligible"])

    def test_annotate_keeps_failed_support_but_does_not_claim_eligibility(self):
        result, _, _ = self.reclassify(self.report([self.row("TTAGGG", 2)]),
                                       self.model({"TTAGGG": (20, 1920)}),
                                       "--support-policy", "annotate")
        self.assertEqual(len(result["candidates"]), 1)
        self.assertFalse(result["candidates"][0]["candidate_eligible"])
        self.assertFalse(result["candidates"][0]["support"]["eligible"])
        self.assertEqual(result["summary"]["eligible_telomere_candidates"], 0)
        self.assertEqual(result["summary"]["call"], "no_call")

    def test_original_evidence_precision_and_scan_metrics_are_preserved(self):
        original = self.report([self.row("TTAGGG", 1000), self.row("AAT", 3000)])
        result, _, _ = self.reclassify(original, self.model())
        self.assertFalse(result["reclassification"]["fastq_reread"])
        for key in ("records", "bases", "elapsed_seconds", "megabases_per_second", "max_record_bases"):
            self.assertEqual(result["summary"][key], original["summary"][key])
        self.assertEqual(result["inputs"], original["inputs"])
        self.assertEqual({motif: {k: v for k, v in row.items() if k not in DERIVED_ROW_KEYS}
                          for motif, row in self.all_rows(result).items()},
                         {row["motif"]: row for row in original["candidates"]})

    def test_reclassification_is_idempotent_and_no_model_resets_support(self):
        model = self.model({"TTAGGG": (20, 1920)})
        first, prefix, _ = self.reclassify(self.report([self.row("TTAGGG", 2)]), model)
        second, second_prefix, _ = self.reclassify(prefix.with_suffix(".json"), model)
        self.assertEqual(first["candidates"], second["candidates"])
        self.assertEqual(first["rejected_candidates"], second["rejected_candidates"])
        self.assertEqual(first["summary"], second["summary"])
        reset, _, _ = self.reclassify(second_prefix.with_suffix(".json"))
        self.assertEqual(len(reset["candidates"]), 1)
        self.assertNotIn("support", reset["candidates"][0])
        self.assertNotIn("candidate_eligible", reset["candidates"][0])
        self.assertEqual(reset["summary"]["eligible_telomere_candidates"], 1)

    def test_assembly_support_is_not_applicable_and_candidate_is_retained(self):
        config = copy.deepcopy(CONFIG)
        config["mode"] = "assembly"
        report = self.report([self.row("TTAGGG", 2)], config)
        report["summary"].update(records=10, bases=10000000, max_record_bases=1000000)
        result, _, _ = self.reclassify(report, self.model({"TTAGGG": (10000, 960000)}))
        row = result["candidates"][0]
        self.assertEqual(row["support"]["status"], "not_applicable")
        self.assertTrue(row["candidate_eligible"])

    def test_missing_model_and_detector_mismatch_fail_before_scan_input_is_opened(self):
        unavailable = self.root / "unavailable_input.fastq"
        missing_model = self.root / "missing_model.json"
        mismatch = self.model()
        for index, options in enumerate((
            ["--background-model", str(missing_model)],
            ["--background-model", str(mismatch), "--max-k", "16"],
        )):
            with self.subTest(options=options):
                prefix = self.root / f"unstarted{index}"
                completed = subprocess.run([str(self.binary), "scan", str(unavailable), "--output", str(prefix),
                                            *options], capture_output=True, text=True, timeout=10)
                self.assertNotEqual(completed.returncode, 0)
                self.assertNotIn(str(unavailable), completed.stderr)
                self.assertTrue(completed.stderr.strip())
                self.assertFalse(prefix.with_suffix(".json").exists())

    def test_invalid_options_are_rejected(self):
        source = self.report([self.row("TTAGGG", 1000)])
        model = self.model()
        for options in (("--support-policy", "invalid"), ("--support-min-enrichment", "0"),
                        ("--support-min-enrichment", "NaN"), ("--support-min-enrichment", "-1"),
                        ("--support-min-enrichment", "inf")):
            with self.subTest(options=options):
                self.reclassify(source, model, *options, success=False)

    def test_scan_can_apply_the_same_support_model(self):
        generator = random.Random(23)
        sequences = ["TTAGGG" * 16] * 20
        sequences += ["".join(generator.choices("ACGT", k=100)) for _ in range(1000)]
        path = self.root / "reads.fastq"
        path.write_bytes(fastq(sequences))
        prefix = self.root / "scanned"
        completed = subprocess.run([str(self.binary), "scan", str(path), "--output", str(prefix),
                                    "--threads", "1", "--background-model", str(self.model())],
                                   capture_output=True, text=True, timeout=20)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        result = json.loads(prefix.with_suffix(".json").read_text())
        row = self.all_rows(result)[canonical("TTAGGG")]
        self.assertTrue(row["support"]["eligible"])
        self.assertTrue(row["candidate_eligible"])

    def test_batch_resume_recomputes_after_model_contents_change_at_same_path(self):
        reads = self.root / "batch_reads.fastq"
        reads.write_bytes(fastq(["TTAGGG" * 16] * 4))
        manifest = self.root / "manifest.tsv"
        manifest.write_text(f"sample_id\tpath1\nsample1\t{reads}\n")
        output = self.root / "batch_output"
        model = self.model()
        command = [sys.executable, str(ROOT / "scripts/batch.py"), str(manifest),
                   "--binary", str(self.binary), "--output-dir", str(output),
                   "--jobs", "1", "--threads", "1", "--cpu-budget", "1",
                   "--background-model", str(model), "--support-policy", "filter"]

        def run(resume=False):
            completed = subprocess.run(command + (["--resume"] if resume else []),
                                       capture_output=True, text=True, timeout=20)
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            return (json.loads((output / "sample1/state.json").read_text()),
                    json.loads((output / "sample1/result.json").read_text()),
                    json.loads((output / "summary.json").read_text()))

        first_state, first_result, first_summary = run()
        self.assertEqual(first_summary["samples"][0]["status"], "complete")
        self.assertTrue(self.all_rows(first_result)[canonical("TTAGGG")]["candidate_eligible"])
        unchanged_state, _, unchanged_summary = run(resume=True)
        self.assertEqual(unchanged_summary["samples"][0]["status"], "resumed")
        self.assertEqual(first_state["fingerprint"], unchanged_state["fingerprint"])

        changed = json.loads(model.read_text())
        changed["controls"][0]["motifs"] = {
            canonical("TTAGGG"): {"reads": 100000, "repeat_bases": 9600000}}
        model.write_text(json.dumps(changed) + "\n")
        second_state, second_result, second_summary = run(resume=True)
        self.assertEqual(second_summary["samples"][0]["status"], "complete")
        self.assertNotEqual(first_state["fingerprint"]["background_model"]["sha256"],
                            second_state["fingerprint"]["background_model"]["sha256"])
        self.assertEqual(second_result["summary"]["eligible_telomere_candidates"], 0)
        self.assertFalse(self.all_rows(second_result)[canonical("TTAGGG")]["candidate_eligible"])


if __name__ == "__main__":
    unittest.main()
