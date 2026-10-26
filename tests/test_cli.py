#!/usr/bin/env python3
"""Subprocess integration checks for the compiled Rust CLI (Python stdlib only).

Install once with `python3 scripts/install_runtime.py`, then run:
  python3 -m unittest discover -s tests -p 'test_*.py' -v
Set TRM_BINARY to test a different executable. Each test uses a temporary
directory and subprocess timeouts, including the bounded-worker error path.
"""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest


from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def canonical(sequence: str) -> str:
    return min(value[index:] + value[:index]
               for value in (sequence, reverse_complement(sequence))
               for index in range(len(sequence)))


def fastq(sequences: list[str]) -> bytes:
    return "".join(f"@read_{index}\n{sequence}\n+\n{'I' * len(sequence)}\n"
                   for index, sequence in enumerate(sequences)).encode("ascii")


class CliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.binary = runtime_binary()
        if not cls.binary.is_file():
            raise RuntimeError(f"build the executable first with cargo build --release: {cls.binary}")

    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="trm-cli-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.invocations = 0

    def write(self, name: str, data: bytes) -> Path:
        path = self.root / name
        path.write_bytes(data)
        return path

    def run_scan(self, paths: list[Path | str], *options: str,
                 stdin: bytes | None = None, success: bool = True,
                 timeout: float = 10, prefix: Path | None = None) -> tuple[dict | None, subprocess.CompletedProcess]:
        self.invocations += 1
        if prefix is None:
            prefix = self.root / f"run_{self.invocations}"
        command = [str(self.binary), "scan", *map(str, paths), "--output", str(prefix),
                   "--threads", "1", "--min-reads", "1", *options]
        completed = subprocess.run(command, input=stdin, capture_output=True,
                                   check=False, timeout=timeout)
        diagnostic = completed.stderr.decode("utf-8", errors="replace")
        if success:
            self.assertEqual(completed.returncode, 0, diagnostic)
            self.assertTrue(Path(str(prefix) + ".tsv").is_file())
            result = json.loads(Path(str(prefix) + ".json").read_text())
            return result, completed
        self.assertNotEqual(completed.returncode, 0, "invalid invocation unexpectedly succeeded")
        self.assertTrue(diagnostic.strip(), "failure must have an actionable diagnostic")
        self.assertFalse(Path(str(prefix) + ".json").exists(), diagnostic)
        self.assertFalse(Path(str(prefix) + ".tsv").exists(), diagnostic)
        self.assertFalse(list(prefix.parent.glob(prefix.name + ".*.part")), diagnostic)
        return None, completed

    def candidate(self, result: dict, motif: str = "TTAGGG") -> dict:
        target = canonical(motif)
        matches = [candidate for candidate in result["candidates"] if candidate["motif"] == target]
        self.assertEqual(len(matches), 1, f"expected one canonical family {target}: {result['candidates']}")
        return matches[0]

    @staticmethod
    def biological_metrics(result: dict) -> dict:
        return {
            "summary": {key: value for key, value in result["summary"].items()
                        if key not in ("elapsed_seconds", "megabases_per_second")},
            "candidates": sorted(
                [{key: value for key, value in candidate.items() if key != "examples"}
                 for candidate in result["candidates"]], key=lambda item: item["motif"]),
            "config": result["config"], "warnings": result["warnings"],
        }

    def test_invalid_options_are_nonzero_without_output(self) -> None:
        path = self.write("valid.fastq", fastq(["TTAGGG" * 25]))
        for options in (("--threads", "0"), ("--threads", "257"),
                        ("--min-k", "0"), ("--min-k", "8", "--max-k", "7"),
                        ("--max-k", "65"), ("--min-identity", "NaN"),
                        ("--min-identity", "inf"), ("--min-identity", "1.1"),
                        ("--min-copies", "1"), ("--max-candidates", "0"),
                        ("--min-tract", "0"), ("--mode", "unknown"),
                        ("--not-an-option", "1")):
            with self.subTest(options=options):
                self.run_scan([path], *options, success=False)

    def test_empty_and_malformed_inputs_are_nonzero_without_output(self) -> None:
        fixtures = [b"", b"\n\n", b"ACGT\n", b">empty\n", b"@x\nACGT\n+\nIII\n",
                    b"@x\nACGT\n+\nIIIII\n", b"@x\nACGT\n+y\nIIII\n",
                    b"@x\nACZT\n+\nIIII\n", b"@x\nACGT\n+\nII I\n"]
        for index, fixture in enumerate(fixtures):
            with self.subTest(index=index):
                self.run_scan([self.write(f"bad_{index}.txt", fixture)], success=False)

    def test_missing_input_is_nonzero(self) -> None:
        self.run_scan([self.root / "missing.fastq"], success=False)

    def test_gzip_crc_truncation_and_garbage_are_rejected(self) -> None:
        compressed = gzip.compress(fastq(["TTAGGG" * 25] * 2), mtime=0)
        corrupt = bytearray(compressed)
        corrupt[-8] ^= 0xFF
        for label, data in (("crc", bytes(corrupt)), ("truncated", compressed[:-4]),
                            ("garbage", compressed + b"garbage")):
            with self.subTest(kind=label):
                self.run_scan([self.write(label + ".fastq.gz", data)], success=False)

    def test_gzip_magic_and_concatenated_members(self) -> None:
        first = gzip.compress(fastq(["TTAGGG" * 25]), mtime=0)
        second = gzip.compress(fastq(["CCCTAA" * 25] * 2), mtime=0)
        path = self.write("without_gzip_suffix.data", first + second)
        result, _ = self.run_scan([path])
        self.assertEqual(result["summary"]["records"], 3)
        self.assertEqual(result["summary"]["bases"], 450)
        self.assertEqual(self.candidate(result)["reads"], 3)

    def test_plain_and_gzip_stdin(self) -> None:
        content = fastq(["TTAGGG" * 25, "CCCTAA" * 25])
        plain, _ = self.run_scan(["-"], stdin=content)
        compressed, _ = self.run_scan(["-"], stdin=gzip.compress(content, mtime=0))
        self.assertEqual(self.biological_metrics(plain), self.biological_metrics(compressed))
        self.assertEqual(plain["summary"]["records"], 2)
        self.assertEqual(plain["inputs"][0]["path"], "-")
        self.assertIsNone(plain["inputs"][0]["bytes"])

    def test_many_batches_agree_between_one_and_four_threads(self) -> None:
        rng = random.Random(92471)
        motifs = ("TTAGGG", "TTAGG", "TTTAGGG", "TTAGGC", "CTCGGTTATGGG", "ACGACACC")
        sequences = []
        for index in range(960):
            if index % 3 == 0:
                sequence = "".join(rng.choices("ACGT", k=150))
            else:
                motif = motifs[index % len(motifs)]
                sequence = (motif * 151)[:150]
                if index % 2:
                    sequence = reverse_complement(sequence)
            sequences.append(sequence)
        path = self.write("many_batches.fastq", fastq(sequences))
        one, _ = self.run_scan([path], "--threads", "1", "--batch-bases", "1000")
        four, _ = self.run_scan([path], "--threads", "4", "--batch-bases", "1000")
        self.assertEqual(one["summary"]["records"], 960)
        self.assertEqual(one["summary"]["bases"], 144000)
        self.assertEqual(self.biological_metrics(one), self.biological_metrics(four))

    def test_record_limit_is_global_across_files_and_flagged(self) -> None:
        first_sequences = ["TTAGGG" * (25 + index) for index in range(11)]
        second_sequences = ["TTAGG" * (30 + index) for index in range(20)]
        paths = [self.write("first.fastq", fastq(first_sequences)),
                 self.write("second.fastq", fastq(second_sequences))]
        result, _ = self.run_scan(paths, "--max-records", "17", "--batch-bases", "1000")
        self.assertEqual(result["summary"]["records"], 17)
        self.assertEqual(result["summary"]["bases"],
                         sum(map(len, first_sequences + second_sequences[:6])))
        self.assertTrue(result["summary"]["sampled_prefix"])
        self.assertTrue(any("prefix" in warning for warning in result["warnings"]))
        full, _ = self.run_scan(paths)
        self.assertEqual(full["summary"]["records"], 31)
        self.assertFalse(full["summary"]["sampled_prefix"])

    def test_candidate_limit_failure_drains_workers_promptly(self) -> None:
        rng = random.Random(4102)
        motifs: set[str] = set()
        while len(motifs) < 100:
            motif = "".join(rng.choices("ACGT", k=11))
            if len(set(motif)) > 1:
                motifs.add(canonical(motif))
        sequences = [motif * 20 for motif in sorted(motifs)] * 4
        path = self.write("too_many_candidates.fastq", fastq(sequences))
        _, completed = self.run_scan([path], "--max-candidates", "1", "--threads", "4",
                                     "--batch-bases", "100", success=False, timeout=5)
        self.assertIn("candidate limit", completed.stderr.decode())

    def test_reverse_complement_swaps_orientation_counts(self) -> None:
        motif = canonical("TTAGGG")
        sequences = [motif * 25] * 3 + [reverse_complement(motif) * 25] * 2
        forward, _ = self.run_scan([self.write("original.fastq", fastq(sequences))])
        reverse, _ = self.run_scan([self.write("rc.fastq", fastq(list(map(reverse_complement, sequences))))])
        first = self.candidate(forward)
        second = self.candidate(reverse)
        self.assertEqual((first["forward_reads"], first["reverse_reads"]), (3, 2))
        self.assertEqual((second["forward_reads"], second["reverse_reads"]), (2, 3))
        self.assertEqual(first["ambiguous_strand_reads"], 0)
        for field in ("reads", "repeat_bases", "max_tract", "period", "long_tract_reads", "status"):
            self.assertEqual(first[field], second[field])

    def test_reverse_complement_symmetric_motif_has_ambiguous_strand(self) -> None:
        result, _ = self.run_scan([self.write("symmetric.fastq", fastq(["ATGC" * 40] * 2))],
                                  "--prior-policy", "discovery")
        candidate = self.candidate(result, "ATGC")
        self.assertEqual(candidate["ambiguous_strand_reads"], 2)
        self.assertEqual(candidate["forward_reads"], 0)
        self.assertEqual(candidate["reverse_reads"], 0)

    def test_short_reads_cannot_supply_a_contiguous_500_base_tract(self) -> None:
        result, _ = self.run_scan([self.write("short.fastq", fastq(["TTAGGG" * 25] * 20))])
        candidate = self.candidate(result)
        self.assertEqual(candidate["repeat_bases"], 3000)
        self.assertEqual(candidate["max_tract"], 150)
        self.assertEqual(candidate["long_tract_reads"], 0)
        self.assertEqual(candidate["status"], "repeat_candidate")

    def test_assembly_terminal_and_internal_evidence_are_distinct(self) -> None:
        rng = random.Random(740)
        middle = "".join(rng.choices("ACGT", k=2000))
        terminal = reverse_complement("TTAGGG" * 100) + middle + "TTAGGG" * 100
        internal = middle + "TTAGGG" * 100 + middle
        terminal_result, _ = self.run_scan([self.write("terminal.fa", f">contig\n{terminal}\n".encode())],
                                          "--mode", "assembly")
        internal_result, _ = self.run_scan([self.write("internal.fa", f">contig\n{internal}\n".encode())],
                                          "--mode", "assembly")
        ends = self.candidate(terminal_result)
        center = self.candidate(internal_result)
        self.assertEqual(ends["status"], "terminal_repeat_candidate")
        self.assertEqual((ends["long_left_end_reads"], ends["long_right_end_reads"]), (1, 1))
        self.assertEqual(ends["opposed_end_records"], 1)
        self.assertEqual(ends["left_forward_reads"] + ends["left_reverse_reads"], 1)
        self.assertEqual(ends["left_forward_reads"], ends["right_reverse_reads"])
        self.assertEqual(ends["left_reverse_reads"], ends["right_forward_reads"])
        self.assertEqual(center["status"], "long_tract_supported")
        self.assertEqual((center["left_end_reads"], center["right_end_reads"]), (0, 0))
        self.assertEqual(center["opposed_end_records"], 0)
        self.assertEqual(center["long_tract_reads"], 1)
        self.assertNotEqual(ends["status"], "confirmed_telomere")
        self.assertNotEqual(center["status"], "confirmed_telomere")

    def test_ambiguity_splits_tracts_without_summing_continuity(self) -> None:
        sequence = "TTAGGG" * 50 + "N" * 12 + "TTAGGG" * 50
        result, _ = self.run_scan([self.write("ambiguous.fastq", fastq([sequence]))])
        candidate = self.candidate(result)
        self.assertEqual(candidate["repeat_bases"], 600)
        self.assertEqual(candidate["max_tract"], 300)
        self.assertEqual(candidate["long_tract_reads"], 0)
        self.assertEqual(candidate["status"], "repeat_candidate")

    def test_exact_minimum_seed_after_an_off_grid_prefix(self) -> None:
        sequence = "N" + "TTAGGG" * 4 + "N"
        result, _ = self.run_scan([self.write("off_grid.fastq", fastq([sequence]))])
        candidate = self.candidate(result)
        self.assertEqual(candidate["reads"], 1)
        self.assertEqual(candidate["repeat_bases"], 24)
        self.assertEqual(candidate["max_tract"], 24)
        self.assertEqual(candidate["long_tract_reads"], 0)

    def test_similar_opposite_strands_do_not_form_a_long_tract(self) -> None:
        motif = "CACCGTATATACGGTT"
        sequence = motif * 20 + reverse_complement(motif) * 20
        self.assertEqual(len(sequence), 640)
        result, _ = self.run_scan([self.write("similar_strand_switch.fastq", fastq([sequence]))],
                                  "--prior-policy", "discovery")
        candidate = self.candidate(result, motif)
        self.assertEqual(candidate["period"], 16)
        self.assertEqual(candidate["reads"], 1)
        self.assertEqual(candidate["forward_reads"], 1)
        self.assertEqual(candidate["reverse_reads"], 1)
        self.assertLess(candidate["max_tract"], 500)
        self.assertEqual(candidate["long_tract_reads"], 0)
        self.assertEqual(candidate["status"], "repeat_candidate")

    def test_overlapping_minimum_seeds_recover_both_families(self) -> None:
        sequence = "NACGACGACGCCGCCGCN"
        result, _ = self.run_scan([self.write("overlapping_seeds.fastq", fastq([sequence]))],
                                  "--min-k", "3", "--max-k", "3", "--seed-bases", "1",
                                  "--min-identity", "1", "--prior-policy", "discovery")
        for motif in ("ACG", "GCC"):
            with self.subTest(motif=motif):
                candidate = self.candidate(result, motif)
                self.assertEqual(candidate["reads"], 1)
                self.assertEqual(candidate["max_tract"], 9)

    def test_output_cannot_overwrite_input(self) -> None:
        content = fastq(["TTAGGG" * 25])
        source = self.write("protected.json", content)
        command = [str(self.binary), "scan", str(source), "--output", str(self.root / "protected")]
        completed = subprocess.run(command, capture_output=True, timeout=10, check=False)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("overwrite", completed.stderr.decode())
        self.assertEqual(source.read_bytes(), content)
        self.assertFalse((self.root / "protected.tsv").exists())


if __name__ == "__main__":
    unittest.main()
