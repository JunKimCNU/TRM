#!/usr/bin/env python3
"""Black-box regressions for seeds hidden by another family's extension.

Run after building the release executable:
  python3 tests/test_traversal_regression.py -v
  python3 tests/test_traversal_regression.py --binary /path/to/trm -v

TRM_BINARY is also supported. The small checked-in real context is always used;
the optional saved full record is exercised when available (or when selected by
--full-record / TRM_C24_FULL_RECORD). This is a regression fixture selected after
finding a discrepancy, not a new biological validation sample.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/c24_cross_family_context.fasta"
PROVENANCE = ROOT / "tests/fixtures/c24_cross_family_context.json"
FAMILY = "AAACCCT"


def reverse_complement(sequence: str) -> str:
    # The source assembly has soft masking. Normalize before complementing.
    return sequence.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


class TraversalRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.binary = runtime_binary()
        if not cls.binary.is_file():
            raise RuntimeError(f"build the release executable first: {cls.binary}")
        cls.provenance = json.loads(PROVENANCE.read_text())
        cls.sequence = "".join(
            line.strip() for line in FIXTURE.read_text().splitlines() if not line.startswith(">")
        ).upper()

    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="trm-traversal-")
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.target = self.directory / "targets.json"
        self.target.write_text(json.dumps({"motifs": ["TTTAGGG"]}))
        self.invocations = 0

    def fasta(self, sequences: list[str]) -> Path:
        path = self.directory / "input.fasta"
        path.write_text("".join(f">fixture_{index}\n{seq}\n" for index, seq in enumerate(sequences)))
        return path

    def run_cli(self, path: Path, action: str = "scan", *options: str) -> dict:
        self.invocations += 1
        prefix = self.directory / f"result_{self.invocations}"
        command = [str(self.binary), action, str(path), "--output", str(prefix),
                   "--mode", "assembly", "--min-reads", "1", "--threads", "1",
                   "--seed-bases", "24", "--min-copies", "3", "--min-identity", "0.85",
                   "--min-tract", "500"]
        if action == "evidence":
            command += ["--targets", str(self.target)]
        else:
            command += ["--min-k", "3", "--max-k", "32", "--prior-policy", "canonical"]
        completed = subprocess.run([*command, *options], capture_output=True, text=True, timeout=30)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return json.loads(prefix.with_suffix(".json").read_text())

    def families(self, report: dict) -> dict:
        # Prior membership must not hide observed, noncanonical repeat families.
        rows = report.get("candidate_evidence", report.get("candidates", []) + report.get("rejected_candidates", []))
        return {row["motif"]: row for row in rows}

    def candidate(self, report: dict) -> dict:
        candidates = self.families(report)
        self.assertIn(FAMILY, candidates)
        return candidates[FAMILY]

    def coordinates(self, candidate: dict) -> list[tuple[int, int, str]]:
        return [(tract["start"], tract["end"], tract["strand"]) for tract in candidate["examples"]]

    def test_fixture_has_recorded_provenance(self) -> None:
        expected = self.provenance["fixture"]
        self.assertEqual(hashlib.sha256(FIXTURE.read_bytes()).hexdigest(), expected["fasta_sha256"])
        self.assertEqual(hashlib.sha256(self.sequence.encode()).hexdigest(), expected["uppercase_sequence_sha256"])
        self.assertEqual(len(self.sequence), 1361)
        source = self.provenance["source"]
        self.assertEqual(source["source_end0"] - source["source_start0"], len(self.sequence))
        self.assertEqual(source["source_start0"] % 7, 0)

    def test_real_cross_family_tract_is_recovered_without_losing_existing_families(self) -> None:
        report = self.run_cli(FIXTURE)
        candidate = self.candidate(report)
        self.assertEqual(report["summary"]["records"], 1)
        self.assertEqual(report["summary"]["bases"], 1361)
        self.assertEqual(candidate["reads"], 1)
        self.assertEqual(candidate["max_tract"], 279)
        self.assertEqual(candidate["repeat_bases"], 307)
        self.assertEqual(self.coordinates(candidate), [(821, 849, "-"), (1005, 1284, "-")])
        self.assertEqual(candidate["long_tract_reads"], 0)
        self.assertEqual(candidate["right_end_reads"], 0)  # The longer tract is 77 bp inset.
        observed = self.families(report)
        for family, previous in self.provenance["existing_other_families"].items():
            with self.subTest(previous_family=family):
                self.assertIn(family, observed)
                self.assertEqual(observed[family]["reads"], 1)
                self.assertGreaterEqual(observed[family]["max_tract"], previous["max_tract"])

    def test_both_strands_and_all_seven_probe_phases_preserve_real_target_tracts(self) -> None:
        length = len(self.sequence)
        for is_reverse in (False, True):
            sequence = reverse_complement(self.sequence) if is_reverse else self.sequence
            expected = [(length - 1284, length - 1005, "+"), (length - 849, length - 821, "+")] if is_reverse else [(821, 849, "-"), (1005, 1284, "-")]
            for phase in range(7):
                with self.subTest(reverse_complement=is_reverse, phase=phase):
                    path = self.fasta(["N" * phase + sequence])
                    wanted = [(start + phase, end + phase, strand) for start, end, strand in expected]
                    # This exact observed context has independently recorded tract
                    # coordinates. The equality is not a promise for arbitrary arrays.
                    for action in ("scan", "evidence"):
                        with self.subTest(action=action):
                            report = self.run_cli(path, action)
                            candidate = self.candidate(report)
                            self.assertEqual(report["summary"]["bases"], length + phase)
                            self.assertEqual(candidate["max_tract"], 279)
                            self.assertEqual(candidate["repeat_bases"], 307)
                            self.assertEqual(self.coordinates(candidate), wanted)
                            self.assertEqual(candidate["long_tract_reads"], 0)

    def test_recovery_does_not_bridge_disjoint_tracts_or_separate_records_to_500bp(self) -> None:
        for sequences, records in (([self.sequence + "N" * 37 + self.sequence], 1),
                                   ([self.sequence, reverse_complement(self.sequence)], 2)):
            with self.subTest(records=records):
                report = self.run_cli(self.fasta(sequences))
                candidate = self.candidate(report)
                self.assertEqual(report["summary"]["records"], records)
                self.assertEqual(candidate["reads"], records)
                self.assertEqual(candidate["repeat_bases"], 614)
                self.assertGreater(candidate["repeat_bases"], 500)
                self.assertEqual(candidate["max_tract"], 279)
                self.assertEqual(candidate["long_tract_reads"], 0)

    def test_batch_and_worker_settings_preserve_recovered_observations(self) -> None:
        path = self.fasta([self.sequence, reverse_complement(self.sequence)] * 5)
        one = self.candidate(self.run_cli(path, "scan", "--batch-bases", "100"))
        many = self.candidate(self.run_cli(path, "scan", "--threads", "4", "--batch-bases", "4000000"))
        self.assertEqual({key: value for key, value in one.items() if key != "examples"},
                         {key: value for key, value in many.items() if key != "examples"})
        self.assertEqual(one["reads"], 10)
        self.assertEqual(one["max_tract"], 279)
        self.assertEqual(one["repeat_bases"], 3070)
        self.assertEqual((one["forward_reads"], one["reverse_reads"]), (5, 5))

    def test_optional_saved_full_record_keeps_family_support_and_consumes_every_base(self) -> None:
        source = self.provenance["source"]
        supplied = os.environ.get("TRM_C24_FULL_RECORD")
        path = Path(supplied) if supplied else ROOT / source["full_record_fixture"]
        if not path.is_file() and not supplied:
            self.skipTest("optional saved C24 full record is not installed")
        self.assertTrue(path.is_file(), path)
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), source["full_record_fixture_sha256"])
        report = self.run_cli(path)
        candidate = self.candidate(report)
        self.assertEqual(report["summary"]["records"], 1)
        self.assertEqual(report["summary"]["bases"], 29413772)
        self.assertFalse(report["summary"]["sampled_prefix"])
        self.assertEqual(candidate["reads"], 1)
        # Allow additional valid discovery tracts; only require preservation and
        # recovery of observations already established for this complete record.
        self.assertGreaterEqual(candidate["max_tract"], 279)
        self.assertGreaterEqual(candidate["repeat_bases"], 3995)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--full-record", type=Path)
    args, unittest_args = parser.parse_known_args()
    if args.binary:
        os.environ["TRM_BINARY"] = str(args.binary)
    if args.full_record:
        os.environ["TRM_C24_FULL_RECORD"] = str(args.full_record)
    unittest.main(argv=[__file__, *unittest_args])
