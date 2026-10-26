"""Failure publication and deterministic export checks; no synthetic DNA fixtures."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from archive_test_support import archive_only

ROOT = Path(__file__).resolve().parents[1]
DEFINITIONS = ROOT / "tests/fixtures/irregular/mixed_definitions.json"
ACTUAL_HIFI = ROOT / "data/matched_japonicus_20260926/SRR25243141/SRR25243141_subreads.fastq.gz"


class IrregularCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binary = Path(os.environ.get("TRM_IRREGULAR_BINARY", ROOT / "target/release/examples/scan_irregular_arrays"))
        if not cls.binary.exists() or not DEFINITIONS.exists():
            raise unittest.SkipTest("Build optional helper and retain published-unit definition first")

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="trm-irregular-cli-")
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)

    def run_cli(self, source, output, *options):
        return subprocess.run([str(self.binary), str(source), "--definitions", str(DEFINITIONS),
            "--output-dir", str(output), "--threads", "2", *options], capture_output=True, text=True, timeout=180)

    def test_empty_and_corrupt_inputs_publish_no_completed_summary(self):
        # These are empty/corrupt file containers, not constructed sequence data.
        for index, content in enumerate((b"", b"\x1f\x8b\x00invalid-container")):
            path = self.folder / f"input{index}"
            path.write_bytes(content)
            output = self.folder / f"out{index}"
            result = self.run_cli(path, output)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertFalse((output / "summary.json").exists())

    def test_existing_directory_is_not_overwritten(self):
        path = self.folder / "empty"
        path.write_bytes(b"")
        output = self.folder / "existing"
        output.mkdir()
        sentinel = output / "summary.json"
        sentinel.write_text('{"existing":true}\n')
        before = sentinel.read_bytes()
        self.assertNotEqual(self.run_cli(path, output).returncode, 0)
        self.assertEqual(sentinel.read_bytes(), before)

    def test_invalid_numeric_parameters_fail(self):
        path = self.folder / "empty"
        path.write_bytes(b"")
        for index, option in enumerate(("--threads", "--min-span", "--min-exact-anchors", "--max-nodes", "--max-record-bases")):
            output = self.folder / f"bad{index}"
            result = subprocess.run([str(self.binary), str(path), "--definitions", str(DEFINITIONS),
                "--output-dir", str(output), option, "0"], capture_output=True, text=True, timeout=30)
            self.assertNotEqual(result.returncode, 0, option)
            self.assertFalse((output / "summary.json").exists())

    @unittest.skipUnless(ACTUAL_HIFI.exists(), "Existing public HiFi archive needed")
    @archive_only
    def test_actual_record_limit_is_failure_not_zero_hit_completion(self):
        output = self.folder / "limited"
        result = self.run_cli(ACTUAL_HIFI, output, "--max-record-bases", "1")
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((output / "summary.json").exists())

    @unittest.skipUnless(ACTUAL_HIFI.exists(), "Existing public HiFi archive needed")
    @archive_only
    def test_actual_hifi_prefix_is_marked_and_parallel_export_is_deterministic(self):
        outputs = []
        for threads in (1, 2):
            output = self.folder / f"threads{threads}"
            result = subprocess.run([str(self.binary), str(ACTUAL_HIFI), "--definitions", str(DEFINITIONS),
                "--output-dir", str(output), "--threads", str(threads), "--max-records", "2000"],
                capture_output=True, text=True, timeout=180)
            self.assertEqual(result.returncode, 0, result.stderr)
            summary = json.loads((output / "summary.json").read_text())
            self.assertEqual(summary["records"], 2000)
            self.assertTrue(summary["sampled_prefix"])
            self.assertFalse(summary["complete_input"])
            for name in ("arrays.jsonl", "anchors.fasta"):
                self.assertEqual(hashlib.sha256((output / name).read_bytes()).hexdigest(), summary["outputs"][name]["sha256"])
            outputs.append(output)
        self.assertEqual((outputs[0] / "arrays.jsonl").read_bytes(), (outputs[1] / "arrays.jsonl").read_bytes())
        self.assertEqual((outputs[0] / "anchors.fasta").read_bytes(), (outputs[1] / "anchors.fasta").read_bytes())


if __name__ == "__main__":
    unittest.main()
