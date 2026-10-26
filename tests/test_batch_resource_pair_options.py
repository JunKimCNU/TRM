"""Resource and pairing options must survive batch orchestration and safe resume."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class BatchResourcePairOptionsTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.addCleanup(self.temporary.cleanup)
        self.binary = self.root / "scanner"
        # The stub isolates batch scheduling and option persistence. Native
        # memory enforcement, mate identities and biology have their own tests.
        self.binary.write_text(
            "#!/usr/bin/env python3\n"
            "import json, pathlib, sys\n"
            "args = sys.argv[1:]\n"
            "prefix = args[args.index('--output') + 1]\n"
            "result = {'summary': {'records': 2, 'bases': 48, 'elapsed_seconds': 0.001}, "
            "'candidates': [], 'received_options': args}\n"
            "pathlib.Path(prefix + '.json').write_text(json.dumps(result))\n"
        )
        self.binary.chmod(0o755)
        for mate in (1, 2):
            (self.root / f"r{mate}.fq").write_text(
                f"@fragment/{mate}\n{'TTAGGG' * 4}\n+\n{'I' * 24}\n")
        self.manifest = self.root / "samples.tsv"
        self.manifest.write_text("sample_id\tpath1\tpath2\nsample\tr1.fq\tr2.fq\n")
        self.output = self.root / "output"

    def run_batch(self, *options, success=True):
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts/batch.py"), str(self.manifest),
             "--binary", str(self.binary), "--output-dir", str(self.output),
             "--jobs", "2", "--threads", "1", "--cpu-budget", "2", *options],
            capture_output=True, text=True, timeout=15)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def state(self):
        return json.loads((self.output / "sample/state.json").read_text())

    def summary(self):
        return json.loads((self.output / "summary.json").read_text())

    def test_forwarding_defaults_and_per_process_limit_reporting(self):
        model = self.root / "background.json"
        model.write_text("{}\n")
        self.run_batch("--memory-limit-mib", "256", "--paired-layout", "fr",
                       "--background-model", str(model))
        options = self.state()["fingerprint"]["options"]
        for name, value in (("--memory-limit-mib", "256"), ("--paired-layout", "fr"),
                            ("--prior-policy", "canonical"), ("--support-policy", "annotate")):
            self.assertEqual(options[options.index(name) + 1], value)
        result = json.loads((self.output / "sample/result.json").read_text())
        self.assertEqual(result["received_options"][-len(options):], options)
        summary = self.summary()
        self.assertEqual(summary["memory_limit_mib_per_scan"], 256)
        self.assertEqual(summary["configured_jobs_times_memory_limit_mib"], 512)
        self.assertEqual(summary["paired_layout"], "fr")
        self.assertIn("not measured aggregate RSS", summary["memory_limit_note"])

    def test_unchanged_options_resume_but_limit_or_layout_changes_recompute(self):
        initial = ("--memory-limit-mib", "256", "--paired-layout", "fr")
        self.run_batch(*initial)
        self.run_batch(*initial, "--resume")
        self.assertEqual(self.summary()["samples"][0]["status"], "resumed")
        self.run_batch("--memory-limit-mib", "512", "--paired-layout", "fr", "--resume")
        self.assertEqual(self.summary()["samples"][0]["status"], "complete")
        self.run_batch("--memory-limit-mib", "512", "--resume")
        self.assertEqual(self.summary()["samples"][0]["status"], "complete")
        self.assertNotIn("--paired-layout", self.state()["fingerprint"]["options"])

    def test_default_independent_record_mode_accepts_single_file(self):
        self.manifest.write_text("sample_id\tpath1\nsample\tr1.fq\n")
        self.run_batch()
        options = self.state()["fingerprint"]["options"]
        self.assertNotIn("--paired-layout", options)
        self.assertNotIn("--memory-limit-mib", options)
        self.assertIsNone(self.summary()["memory_limit_mib_per_scan"])
        self.assertIsNone(self.summary()["configured_jobs_times_memory_limit_mib"])

    def test_any_single_file_row_rejects_paired_batch_before_launch(self):
        self.manifest.write_text(
            "sample_id\tpath1\tpath2\npaired\tr1.fq\tr2.fq\nsingle\tr1.fq\t\n")
        result = self.run_batch("--paired-layout", "fr", success=False)
        self.assertIn("sample single", result.stderr)
        self.assertIn("requires exactly two", result.stderr)
        self.assertFalse(self.output.exists())

    def test_assembly_rejects_paired_mode_before_launch(self):
        result = self.run_batch("--paired-layout", "fr", "--mode", "assembly", success=False)
        self.assertIn("requires --mode reads", result.stderr)
        self.assertFalse(self.output.exists())

    def test_nonpositive_memory_limits_rejected_before_launch(self):
        for limit in ("0", "-1"):
            with self.subTest(limit=limit):
                result = self.run_batch("--memory-limit-mib", limit, success=False)
                self.assertIn("positive integer", result.stderr)
                self.assertFalse(self.output.exists())

    def test_explicit_legacy_policies_remain_available(self):
        model = self.root / "background.json"
        model.write_text("{}\n")
        self.run_batch("--background-model", str(model), "--prior-policy", "telomere",
                       "--support-policy", "filter")
        options = self.state()["fingerprint"]["options"]
        self.assertEqual(options[options.index("--prior-policy") + 1], "telomere")
        self.assertEqual(options[options.index("--support-policy") + 1], "filter")


if __name__ == "__main__":
    unittest.main()
