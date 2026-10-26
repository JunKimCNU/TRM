"""Cohort scheduling must preserve independent results, resume, and shard scope."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
BATCH = ROOT / "scripts/batch.py"


def shard(sample_id, count):
    return int.from_bytes(hashlib.sha256(sample_id.encode()).digest(), "big") % count


class BatchScalingTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.manifest = self.root / "samples.tsv"
        self.output = self.root / "output"
        self.binary = self.root / "scanner"
        self.binary.write_text(
            "#!/usr/bin/env python3\n"
            "import json, pathlib, sys, time\n"
            "args = sys.argv[1:]\n"
            "source = pathlib.Path(args[1])\n"
            "with (source.parent / 'launches.txt').open('a') as stream:\n"
            "    stream.write(source.name + '\\n')\n"
            "time.sleep(0.12 if source.name.startswith('slow') else 0.005)\n"
            "if source.name.startswith('fail'):\n"
            "    sys.exit(7)\n"
            "prefix = args[args.index('--output') + 1]\n"
            "report = {'summary': {'records': 1, 'bases': source.stat().st_size},\n"
            "          'candidates': [{'motif': 'AACCCT', 'status': 'observed'}]}\n"
            "pathlib.Path(prefix + '.json').write_text(json.dumps(report))\n"
        )
        self.binary.chmod(0o755)

    def make_samples(self, names, missing=()):
        self.manifest.write_text("sample_id\tpath1\n" + "".join(
            f"{name}\t{name}.fq\n" for name in names))
        for index, name in enumerate(names):
            if name not in missing:
                (self.root / f"{name}.fq").write_text("A" * (index + 1))

    def run_batch(self, *options, expected=0, jobs=4):
        result = subprocess.run(
            [sys.executable, str(BATCH), str(self.manifest), "--binary", str(self.binary),
             "--output-dir", str(self.output), "--jobs", str(jobs), "--threads", "1",
             "--cpu-budget", str(jobs), *options],
            capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return result

    def summary(self, index=None, count=None):
        root = self.output if index is None else self.output / f"shard-{index:05d}-of-{count:05d}"
        return json.loads((root / "summary.json").read_text())

    def launches(self):
        path = self.root / "launches.txt"
        return path.read_text().splitlines() if path.exists() else []

    def test_bounded_jobs_coalesced_checkpoints_and_manifest_order(self):
        names = [f"sample{i:03}" for i in range(80)]
        self.make_samples(names)
        self.run_batch("--checkpoint-seconds", "300")
        report = self.summary()
        self.assertEqual(report["max_inflight"], 4)
        self.assertEqual(report["checkpoint_count"], 2)
        self.assertEqual(report["total_samples"], 80)
        self.assertEqual([item["sample_id"] for item in report["samples"]], names)
        self.assertTrue(report["all_succeeded"])
        self.assertEqual(len(self.launches()), 80)
        progress = json.loads((self.output / "progress.json").read_text())
        self.assertEqual((progress["pending"], progress["running_count"], progress["inflight"]), (0, 0, 0))
        self.assertEqual(progress["submitted"], 80)
        with (self.output / "summary.tsv").open() as source:
            rows = list(csv.DictReader(source, delimiter="\t"))
        self.assertEqual([row["sample_id"] for row in rows], names)

    def test_largest_first_changes_launch_order_but_not_results_or_resume(self):
        names = ["small", "medium", "large"]
        self.make_samples(names)
        self.run_batch("--schedule", "largest-first", "--checkpoint-seconds", "1e100", jobs=1)
        self.assertEqual(self.launches(), ["large.fq", "medium.fq", "small.fq"])
        report = self.summary()
        self.assertEqual([row["sample_id"] for row in report["samples"]], names)
        original = {name: (self.output / name / "result.json").read_bytes() for name in names}
        self.run_batch("--resume", "--schedule", "manifest", jobs=1)
        self.assertEqual(len(self.launches()), 3)
        self.assertEqual(self.summary()["resumed_samples"], 3)
        for name in names:
            self.assertEqual((self.output / name / "result.json").read_bytes(), original[name])

    def test_periodic_checkpoints_continue_while_long_sample_is_running(self):
        self.make_samples(["slow_sample"])
        self.run_batch("--checkpoint-seconds", "0.03", jobs=1)
        self.assertGreater(self.summary()["checkpoint_count"], 2)
        self.assertEqual(self.summary()["completed_samples"], 1)

    def test_failed_sample_and_corrupt_result_recompute_without_rescanning_good_sample(self):
        self.make_samples(["good", "corrupt", "fail_input"])
        self.run_batch(expected=1)
        (self.output / "corrupt/result.json").write_text("{}")
        self.run_batch("--resume", expected=1)
        self.assertEqual(self.launches().count("good.fq"), 1)
        self.assertEqual(self.launches().count("corrupt.fq"), 2)
        self.assertEqual(self.launches().count("fail_input.fq"), 2)
        report = self.summary()
        self.assertFalse(report["all_succeeded"])
        self.assertEqual((report["completed_samples"], report["failed_samples"], report["resumed_samples"]), (2, 1, 1))
        self.assertEqual(report["samples"][-1]["status"], "failed")

    def test_shards_cover_manifest_once_and_keep_results_separate(self):
        names = [f"sample{i:03}" for i in range(15)]
        self.make_samples(names)
        observed = []
        for index in range(4):
            self.run_batch("--shard-index", str(index), "--shard-count", "4")
            report = self.summary(index, 4)
            selected = [name for name in names if shard(name, 4) == index]
            self.assertEqual([row["sample_id"] for row in report["samples"]], selected)
            self.assertEqual(report["shard_index"], index)
            observed.extend(selected)
        self.assertEqual(sorted(observed), sorted(names))
        self.assertEqual(len(observed), len(set(observed)))
        self.assertEqual(len(self.launches()), 15)
        # Reordering the complete manifest does not change shard membership.
        self.make_samples(list(reversed(names)))
        self.run_batch("--shard-index", "1", "--shard-count", "4")
        self.assertEqual([row["sample_id"] for row in self.summary(1, 4)["samples"]],
                         [name for name in reversed(names) if shard(name, 4) == 1])

    def test_unselected_missing_inputs_are_not_statted_and_empty_shard_succeeds(self):
        names = ["only"]
        self.make_samples(names, missing=names)
        empty_index = 1 - shard("only", 2)
        self.run_batch("--shard-index", str(empty_index), "--shard-count", "2")
        report = self.summary(empty_index, 2)
        self.assertTrue(report["all_succeeded"])
        self.assertEqual(report["total_samples"], 0)
        self.assertEqual(report["max_inflight"], 0)
        self.assertEqual(self.launches(), [])
        table = self.output / f"shard-{empty_index:05d}-of-00002/summary.tsv"
        self.assertEqual(len(table.read_text().splitlines()), 1)
        self.run_batch("--shard-index", str(1 - empty_index), "--shard-count", "2", expected=2)

    def test_invalid_unselected_manifest_rows_still_rejected(self):
        self.make_samples(["only"], missing=["only"])
        empty_index = 1 - shard("only", 2)
        with self.manifest.open("a") as source:
            source.write("only\tabsent.fq\n")
        result = self.run_batch("--shard-index", str(empty_index), "--shard-count", "2", expected=2)
        self.assertIn("duplicate sample_id", result.stderr)
        self.assertFalse(self.output.exists())

    def test_invalid_checkpoint_and_partition_options_rejected_before_launch(self):
        self.make_samples(["sample"])
        for options in (("--checkpoint-seconds", "0"), ("--checkpoint-seconds", "nan"),
                        ("--checkpoint-seconds", "inf"), ("--shard-index", "0"),
                        ("--shard-count", "2"), ("--shard-index", "-1", "--shard-count", "2"),
                        ("--shard-index", "2", "--shard-count", "2")):
            with self.subTest(options=options):
                self.run_batch(*options, expected=2)
        self.assertEqual(self.launches(), [])
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
