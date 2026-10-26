#!/usr/bin/env python3
"""Background-profile provenance, completeness, and output-integrity checks."""

import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/calibrate_support.py"
SPEC = importlib.util.spec_from_file_location("calibrate_support", SCRIPT)
calibrate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(calibrate)


def fixture():
    return {
        "program": "trm", "schema_version": 2, "version": "0.2.0",
        "inputs": [{"path": "/example/negative.fastq.gz", "bytes": 1000}],
        "config": {"mode": "reads", "min_k": 3, "max_k": 32, "seed_bases": 24,
                   "min_copies": 3, "min_identity": 0.85, "min_tract": 500,
                   "end_window": 50, "min_reads": 1, "max_candidates": 100000},
        "summary": {"records": 1000, "bases": 100000, "max_record_bases": 100,
                    "sampled_prefix": False, "candidate_families_before_filter": 2,
                    "raw_reported_candidates": 2, "reported_candidates": 1},
        "candidates": [{"motif": "AACCCT", "period": 6, "reads": 10,
                        "repeat_bases": 700, "telomere_prior": {"eligible": True}}],
        "rejected_candidates": [{"motif": "AAAAAAAG", "period": 8, "reads": 5,
                                 "repeat_bases": 250, "telomere_prior": {"eligible": False}}],
    }


class CalibrationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="trm-calibration-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "negative.json"
        self.output = self.root / "profile.json"
        self.write(fixture())

    def write(self, value, path=None):
        path = path or self.source
        path.write_text(json.dumps(value), encoding="utf-8")
        return path

    def test_unknown_motif_included_and_source_counts_preserved(self):
        original = self.source.read_bytes()
        model = calibrate.write_model([self.source], self.output, "test background")
        self.assertEqual(model["model_type"], "empirical_repeat_background")
        self.assertEqual(model["name"], "test background")
        control = model["controls"][0]
        self.assertEqual(control["motifs"]["AAAAAAAG"], {"reads": 5, "repeat_bases": 250})
        self.assertEqual(control["config"], fixture()["config"])
        self.assertEqual(control["inputs"], fixture()["inputs"])
        self.assertEqual(control["source_result_sha256"], hashlib.sha256(original).hexdigest())
        self.assertEqual(control["source_result"], str(self.source.resolve()))
        self.assertEqual(model, json.loads(self.output.read_text()))
        self.assertEqual(original, self.source.read_bytes())

    def test_unseen_motifs_not_fabricated_and_censor_threshold_preserved(self):
        value = fixture()
        value["config"]["min_reads"] = 2
        value["summary"]["candidate_families_before_filter"] = 9
        self.write(value)
        control = calibrate.build_model([self.source])["controls"][0]
        self.assertEqual(control["config"]["min_reads"], 2)
        self.assertEqual(set(control["motifs"]), {"AAAAAAAG", "AACCCT"})
        self.assertNotIn("AAGCCT", control["motifs"])

    def test_duplicate_rows_ignore_only_annotations(self):
        value = fixture()
        duplicate = copy.deepcopy(value["candidates"][0])
        duplicate["telomere_prior"] = {"eligible": False}
        duplicate["support"] = {"eligible": False}
        duplicate["candidate_eligible"] = False
        value["rejected_candidates"].append(duplicate)
        self.write(value)
        self.assertEqual(len(calibrate.build_model([self.source])["controls"][0]["motifs"]), 2)
        duplicate["repeat_bases"] += 1
        self.write(value)
        with self.assertRaisesRegex(ValueError, "Conflicting duplicate"):
            calibrate.build_model([self.source])

    def test_invalid_counts_rejected(self):
        changes = [("records", 0), ("records", True), ("bases", -1),
                   ("bases", 100001), ("max_record_bases", 0)]
        for key, change in changes:
            with self.subTest(key=key, value=change):
                value = fixture()
                value["summary"][key] = change
                self.write(value)
                with self.assertRaises(ValueError):
                    calibrate.build_model([self.source])
        for key, change in [("reads", 1001), ("reads", True), ("reads", 0),
                            ("repeat_bases", 100001), ("repeat_bases", 1001),
                            ("repeat_bases", 0), ("period", 5)]:
            with self.subTest(key=key, value=change):
                value = fixture()
                value["candidates"][0][key] = change
                self.write(value)
                with self.assertRaises(ValueError):
                    calibrate.build_model([self.source])

    def test_missing_required_summary_or_counts_rejected(self):
        for container, key in [("summary", "records"), ("summary", "sampled_prefix"),
                               ("summary", "raw_reported_candidates"),
                               ("candidate", "reads"), ("candidate", "repeat_bases")]:
            with self.subTest(container=container, key=key):
                value = fixture()
                del (value["candidates"][0] if container == "candidate" else value[container])[key]
                self.write(value)
                with self.assertRaises(ValueError):
                    calibrate.build_model([self.source])

    def test_unknown_incompatible_and_invalid_configs_rejected(self):
        for key, change in [("mode", "assembly"), ("min_k", 0), ("max_k", 65),
                            ("min_identity", 1.1), ("min_identity", True),
                            ("min_reads", 0), ("seed_bases", 0), ("min_copies", 1),
                            ("extra_option", 1)]:
            with self.subTest(key=key, value=change):
                value = fixture()
                value["config"][key] = change
                self.write(value)
                with self.assertRaises(ValueError):
                    calibrate.build_model([self.source])
        self.write(fixture())
        second = fixture()
        second["config"]["seed_bases"] = 30
        other = self.write(second, self.root / "other.json")
        with self.assertRaisesRegex(ValueError, "incompatible detector"):
            calibrate.build_model([self.source, other])

    def test_partial_sampled_and_cap_limited_results_rejected(self):
        for key, change in [("sampled_prefix", True), ("candidate_families_before_filter", 3),
                            ("raw_reported_candidates", 3), ("reported_candidates", 0),
                            ("complete", False), ("status", "failed"),
                            ("candidate_families_before_filter", 100000)]:
            with self.subTest(key=key, value=change):
                value = fixture()
                value["summary"][key] = change
                self.write(value)
                with self.assertRaises(ValueError):
                    calibrate.build_model([self.source])
        for inputs in ([], None, [{"path": ""}], [{"path": "file", "bytes": -1}]):
            value = fixture()
            value["inputs"] = inputs
            self.write(value)
            with self.assertRaises(ValueError):
                calibrate.build_model([self.source])

    def test_invalid_noncanonical_and_nonprimitive_motifs_rejected(self):
        for motif in ("AAAN", "aaccc", "TTAGGG", "AATAAT", "", "A" * 65):
            with self.subTest(motif=motif):
                value = fixture()
                value["candidates"][0]["motif"] = motif
                value["candidates"][0]["period"] = len(motif)
                self.write(value)
                with self.assertRaises(ValueError):
                    calibrate.build_model([self.source])

    def test_duplicate_sources_hashes_and_control_ids_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate negative control"):
            calibrate.build_model([self.source, self.source])
        other = self.root / "other.json"
        other.write_bytes(self.source.read_bytes())
        with self.assertRaisesRegex(ValueError, "source_result_sha256"):
            calibrate.build_model([self.source, other])
        directory = self.root / "subdirectory"
        directory.mkdir()
        value = fixture()
        value["version"] = "0.3.0"
        duplicate_id = self.write(value, directory / self.source.name)
        with self.assertRaisesRegex(ValueError, "control_id"):
            calibrate.build_model([self.source, duplicate_id])

    def test_original_schema_and_different_censoring_thresholds_accepted(self):
        value = fixture()
        value["schema_version"] = 1
        value["candidates"] += value.pop("rejected_candidates")
        value["summary"]["reported_candidates"] = 2
        del value["summary"]["raw_reported_candidates"]
        self.write(value)
        other_value = fixture()
        other_value["config"]["min_reads"] = 2
        other_value["summary"]["candidate_families_before_filter"] = 9
        other = self.write(other_value, self.root / "other.json")
        self.assertEqual(len(calibrate.build_model([self.source, other])["controls"]), 2)

    def test_output_collision_and_failure_preserve_existing_files(self):
        original = self.source.read_bytes()
        with self.assertRaisesRegex(ValueError, "source result"):
            calibrate.write_model([self.source], self.source)
        self.output.write_text("existing model\n")
        with self.assertRaises(FileExistsError):
            calibrate.write_model([self.source], self.output)
        self.assertEqual(self.output.read_text(), "existing model\n")
        self.assertEqual(self.source.read_bytes(), original)
        self.output.unlink()
        self.source.write_text('{"partial":')
        with self.assertRaises(ValueError):
            calibrate.write_model([self.source], self.output)
        self.assertFalse(self.output.exists())
        self.assertFalse(list(self.root.glob(".background-*.part")))

    def test_nonfinite_and_duplicate_json_keys_rejected(self):
        self.source.write_text('{"program":"trm","program":"different"}')
        with self.assertRaisesRegex(ValueError, "Duplicate JSON key"):
            calibrate.build_model([self.source])
        self.source.write_text(json.dumps(fixture()).replace('0.85', 'NaN'))
        with self.assertRaisesRegex(ValueError, "Non-finite"):
            calibrate.build_model([self.source])

    def test_cli_requires_negative_inputs_and_builds_model(self):
        missing = subprocess.run([sys.executable, str(SCRIPT), "--output", str(self.output)],
                                 text=True, capture_output=True, timeout=10)
        self.assertNotEqual(missing.returncode, 0)
        completed = subprocess.run([sys.executable, str(SCRIPT), "--negative-result", str(self.source),
                                    "--output", str(self.output), "--name", "CLI check"],
                                   text=True, capture_output=True, timeout=10)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(json.loads(self.output.read_text())["name"], "CLI check")


if __name__ == "__main__":
    unittest.main()
