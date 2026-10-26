#!/usr/bin/env python3
"""Black-box checks of whole-read scope, positions and continuity accounting.

Run with the current release binary:
  python3 -m unittest discover -s tests -p 'test_full_read_scope.py' -v
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest


from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]
FAMILY = "AACCCT"  # Rotation / reverse-complement canonical form of TTAGGG.


def background(length: int, seed: int = 91037) -> str:
    return "".join(random.Random(seed).choices("ACGT", k=length))


class FullReadScopeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.binary = runtime_binary()
        if not cls.binary.is_file():
            raise RuntimeError(f"build the release executable first: {cls.binary}")

    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="trm-full-read-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.invocations = 0

    def fastq(self, name: str, sequences: list[str], quality: str = "I") -> Path:
        path = self.root / name
        with path.open("w", encoding="ascii") as stream:
            for index, sequence in enumerate(sequences):
                stream.write(f"@read_{index}\n{sequence}\n+\n{quality * len(sequence)}\n")
        return path

    def scan(self, paths: list[Path], *options: str) -> dict:
        self.invocations += 1
        prefix = self.root / f"result_{self.invocations}"
        completed = subprocess.run(
            [str(self.binary), "scan", *map(str, paths), "--output", str(prefix),
             "--threads", "1", "--min-reads", "1", *options],
            text=True, capture_output=True, check=False, timeout=15,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return json.loads(prefix.with_suffix(".json").read_text(encoding="utf-8"))

    def candidate(self, report: dict) -> dict:
        matches = [row for row in report["candidates"] if row["motif"] == FAMILY]
        self.assertEqual(len(matches), 1, report["candidates"])
        return matches[0]

    def test_signals_at_read_start_middle_and_end_have_observed_coordinates(self) -> None:
        repeat = "TTAGGG" * 10
        flank = background(257)
        fixtures = [
            (repeat + "N" + flank, 0),
            (flank + "N" + repeat + "N" + flank, 258),
            (flank + "N" + repeat, 258),
        ]
        for index, (sequence, start) in enumerate(fixtures):
            with self.subTest(position=index):
                report = self.scan([self.fastq(f"position_{index}.fastq", [sequence])])
                candidate = self.candidate(report)
                self.assertEqual(report["summary"]["bases"], len(sequence))
                self.assertEqual(candidate["repeat_bases"], 60)
                self.assertEqual(candidate["max_tract"], 60)
                self.assertEqual([(item["start"], item["end"]) for item in candidate["examples"]],
                                 [(start, start + 60)])

    def test_minimum_seed_is_seen_at_all_six_phases_beyond_base_150(self) -> None:
        for phase in range(6):
            with self.subTest(phase=phase):
                start = 180 + phase
                sequence = background(start - 1) + "N" + "TTAGGG" * 4 + "N" + background(53)
                report = self.scan([self.fastq(f"phase_{phase}.fastq", [sequence])])
                candidate = self.candidate(report)
                self.assertEqual(candidate["repeat_bases"], 24)
                self.assertEqual(candidate["max_tract"], 24)
                self.assertEqual(candidate["examples"][0]["start"], start)
                self.assertEqual(candidate["examples"][0]["end"], start + 24)

    def test_late_long_read_signal_is_not_cut_at_batch_or_read_length_boundary(self) -> None:
        start = 8193
        sequence = background(start - 1) + "N" + "TTAGGG" * 100 + "N" + background(511)
        path = self.fastq("long_read.fastq", [sequence])
        small_batch = self.scan([path], "--batch-bases", "100")
        large_batch = self.scan([path], "--batch-bases", "4000000")
        candidate = self.candidate(small_batch)
        self.assertEqual(candidate, self.candidate(large_batch))
        self.assertEqual(small_batch["summary"]["bases"], len(sequence))
        self.assertEqual(small_batch["summary"]["max_record_bases"], len(sequence))
        self.assertEqual(candidate["repeat_bases"], 600)
        self.assertEqual(candidate["max_tract"], 600)
        self.assertEqual(candidate["long_tract_reads"], 1)
        self.assertEqual(candidate["status"], "long_tract_supported")
        self.assertEqual(candidate["examples"][0]["start"], start)
        self.assertEqual(candidate["examples"][0]["end"], start + 600)

    def test_two_disjoint_tracts_in_one_full_read_do_not_create_500bp_continuity(self) -> None:
        sequence = "TTAGGG" * 50 + "N" * 37 + "TTAGGG" * 50
        report = self.scan([self.fastq("separate_tracts.fastq", [sequence])])
        candidate = self.candidate(report)
        self.assertEqual(report["summary"]["bases"], 637)
        self.assertEqual(candidate["reads"], 1)
        self.assertEqual(candidate["repeat_bases"], 600)
        self.assertEqual(candidate["max_tract"], 300)
        self.assertEqual(candidate["long_tract_reads"], 0)
        self.assertEqual(candidate["full_repeat_reads"], 0)
        self.assertEqual(candidate["status"], "repeat_candidate")
        self.assertEqual([(item["start"], item["end"]) for item in candidate["examples"]],
                         [(0, 300), (337, 637)])

    def test_isolated_motif_occurrences_are_not_a_tandem_repeat_census(self) -> None:
        tandem = self.scan([self.fastq("tandem.fastq", ["TTAGGG" * 100])])
        isolated = self.scan([self.fastq("isolated.fastq", [("TTAGGG" + "N") * 100])])
        self.assertEqual(self.candidate(tandem)["repeat_bases"], 600)
        self.assertFalse(any(row["motif"] == FAMILY for row in isolated["candidates"]))
        self.assertEqual(isolated["summary"]["bases"], 700)

    def test_all_records_from_all_inputs_are_used_unless_prefix_is_requested(self) -> None:
        first_sequences = [background(257, seed=index) for index in range(5)]
        second_sequences = [background(258, seed=20 + index) for index in range(5)]
        second_sequences.append("TTAGGG" * 100)
        paths = [self.fastq("sample_R1.fastq", first_sequences),
                 self.fastq("sample_R2.fastq", second_sequences)]
        complete = self.scan(paths)
        prefix = self.scan(paths, "--max-records", "10")
        self.assertEqual(complete["summary"]["records"], 11)
        self.assertEqual(complete["summary"]["bases"], sum(map(len, first_sequences + second_sequences)))
        self.assertFalse(complete["summary"]["sampled_prefix"])
        self.assertEqual(self.candidate(complete)["reads"], 1)
        self.assertEqual(self.candidate(complete)["examples"][0]["input_index"], 1)
        self.assertEqual(prefix["summary"]["records"], 10)
        self.assertTrue(prefix["summary"]["sampled_prefix"])
        self.assertFalse(any(row["motif"] == FAMILY for row in prefix["candidates"]))
        self.assertTrue(any("prefix" in warning for warning in prefix["warnings"]))

    def test_r1_r2_are_complete_separate_records_without_inferred_continuity(self) -> None:
        paths = [self.fastq("sample_R1.fastq", ["TTAGGG" * 50]),
                 self.fastq("sample_R2.fastq", ["CCCTAA" * 50])]
        report = self.scan(paths)
        candidate = self.candidate(report)
        self.assertEqual(report["summary"]["records"], 2)
        self.assertEqual(report["summary"]["bases"], 600)
        self.assertEqual(candidate["reads"], 2)
        self.assertEqual(candidate["repeat_bases"], 600)
        self.assertEqual(candidate["max_tract"], 300)
        self.assertEqual(candidate["long_tract_reads"], 0)
        self.assertEqual((candidate["forward_reads"], candidate["reverse_reads"]), (1, 1))

    def test_valid_quality_values_do_not_trim_or_weight_sequence_evidence(self) -> None:
        sequence = background(167) + "N" + "TTAGGG" * 10 + "N" + background(77)
        high = self.scan([self.fastq("high_quality.fastq", [sequence], quality="I")])
        low = self.scan([self.fastq("low_quality.fastq", [sequence], quality="!")])
        self.assertEqual(self.candidate(high), self.candidate(low))
        self.assertEqual(low["summary"]["bases"], len(sequence))
        self.assertEqual(self.candidate(low)["examples"][0]["start"], 168)


if __name__ == "__main__":
    unittest.main()
