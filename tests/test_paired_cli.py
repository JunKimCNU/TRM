"""FR pair annotations, input integrity, and unchanged whole-read evidence.

Run against a compiled binary with TRM_BINARY, as for the other CLI tests.
Fixtures express molecule-level expectations without interpreting mates as
chromosome strands or treating two records as a 500 bp observed repeat.
"""

import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]
FLANK = "ACGTCCGATGCTAACGTCAGTGCATCAGACGTTCGACGTA"


def canonical(sequence):
    reverse = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return min(word[index:] + word[:index]
               for word in (sequence, reverse) for index in range(len(word)))


class PairedCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binary = runtime_binary()

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="trm-paired-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.counter = 0
        self.targets = self.root / "targets.json"
        self.targets.write_text(json.dumps({"motifs": ["TTAGGG", "TTAGGC"]}))

    def write_fastq(self, name, sequences, headers=None, compressed=False):
        path = self.root / name
        if headers is None:
            headers = [f"read{i}" for i in range(len(sequences))]
        self.assertEqual(len(headers), len(sequences))
        data = "".join(f"@{header}\n{seq}\n+\n{'I' * len(seq)}\n"
                       for header, seq in zip(headers, sequences)).encode()
        path.write_bytes(gzip.compress(data) if compressed else data)
        return path

    def write_pair(self, pairs, headers=None):
        left = [pair[0] for pair in pairs]
        right = [pair[1] for pair in pairs]
        left_headers, right_headers = headers if headers else (None, None)
        return [self.write_fastq("R1.fastq", left, left_headers),
                self.write_fastq("R2.fastq", right, right_headers)]

    def run_cli(self, paths, *options, mode="scan", paired=True, success=True):
        self.counter += 1
        output = self.root / f"result{self.counter}"
        command = [str(self.binary), mode, *map(str, paths), "--output", str(output),
                   "--threads", "1"]
        if mode != "mixed":
            command += ["--min-reads", "1", "--min-identity", "1"]
        if mode == "scan":
            command += ["--prior-policy", "discovery"]
        if mode == "evidence":
            command += ["--targets", str(self.targets)]
        if paired:
            command += ["--paired-layout", "fr"]
        completed = subprocess.run([*command, *options], text=True, capture_output=True, timeout=20)
        if not success:
            self.assertNotEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue(completed.stderr.strip())
            for suffix in (".json", ".tsv", ".rejected.tsv"):
                self.assertFalse(Path(str(output) + suffix).exists(), completed.stderr)
            self.assertFalse(list(self.root.glob(output.name + ".*.part")))
            return completed
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.output = output
        return json.loads(output.with_suffix(".json").read_text())

    @staticmethod
    def rows(report):
        return report.get("candidate_evidence", report.get("candidates", []))

    def row(self, report, motif="TTAGGG"):
        key = canonical(motif)
        return next(row for row in self.rows(report) if row["motif"] == key)

    def comparable(self, report, omit_pairs=False):
        excluded = {"examples"}
        if omit_pairs:
            excluded.add("paired_evidence")
        return {row["motif"]: {key: value for key, value in row.items() if key not in excluded}
                for row in self.rows(report)}

    def geometry_fixture(self):
        g, c = "TTAGGG" * 20, "CCCTAA" * 20
        return [(g, c), (g, g), ("TTAGGG" * 5 + "N" * 10 + "CCCTAA" * 5, c),
                (g, "N" * 120), ("N" * 120, c),
                ("GGTTAG" * 20, "CCCTAA" * 10 + FLANK),
                ("N" * 120, "N" * 120)]

    def test_opposite_same_ambiguous_and_one_sided_pairs_are_distinct(self):
        result = self.run_cli(self.write_pair(self.geometry_fixture()))
        stats = self.row(result)["paired_evidence"]
        self.assertEqual(stats, {
            "supporting_pairs": 6,
            "both_mates_repeat_pairs": 4,
            "r1_only_repeat_pairs": 1,
            "r2_only_repeat_pairs": 1,
            "opposite_strand_pairs": 2,
            "same_strand_pairs": 1,
            "ambiguous_strand_pairs": 1,
            "both_full_repeat_pairs": 2,
            "exactly_one_full_repeat_pairs": 4,
        })
        self.assertEqual(result["paired_end"]["validated_pairs"], 7)
        self.assertEqual(result["paired_end"]["scope"], "complete_inputs")
        self.assertEqual(result["summary"]["records"], 14)
        self.assertEqual(self.row(result)["reads"], 10)

    def test_paired_annotation_preserves_all_independent_read_metrics(self):
        paths = self.write_pair(self.geometry_fixture())
        for mode in ("scan", "evidence"):
            with self.subTest(mode=mode):
                independent = self.run_cli(paths, paired=False, mode=mode)
                paired = self.run_cli(paths, mode=mode)
                self.assertEqual(self.comparable(independent, omit_pairs=True),
                                 self.comparable(paired, omit_pairs=True))
                self.assertEqual(independent["evidence_qc"], paired["evidence_qc"])
                ignored = {"elapsed_seconds", "megabases_per_second"}
                self.assertEqual({k: v for k, v in independent["summary"].items() if k not in ignored},
                                 {k: v for k, v in paired["summary"].items() if k not in ignored})
                self.assertFalse(independent["paired_end"]["measured"])
                self.assertIsNone(independent["paired_end"]["validated_pairs"])
                for row in self.rows(independent):
                    self.assertIsNone(row.get("paired_evidence"))

    def test_phase_offset_and_repeat_only_mates_never_join_to_500bp(self):
        motif = "TTAGGG"
        pairs = []
        for shift in range(len(motif)):
            rotated = motif[shift:] + motif[:shift]
            reverse = rotated.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            pairs.append((rotated * 50, reverse * 50))
        row = self.row(self.run_cli(self.write_pair(pairs)))
        self.assertEqual(row["paired_evidence"]["opposite_strand_pairs"], 6)
        self.assertEqual(row["max_tract"], 300)
        self.assertEqual(row["long_tract_reads"], 0)
        self.assertEqual(row["junction_evidence"]["junction_reads"], 0)

    def test_symmetric_repeat_family_is_orientation_ambiguous(self):
        row = self.row(self.run_cli(self.write_pair([("ATGC" * 30, "ATGC" * 30)])), "ATGC")
        stats = row["paired_evidence"]
        self.assertEqual(stats["both_mates_repeat_pairs"], 1)
        self.assertEqual(stats["ambiguous_strand_pairs"], 1)
        self.assertEqual(stats["opposite_strand_pairs"], 0)
        self.assertEqual(stats["same_strand_pairs"], 0)

    def test_one_and_four_workers_preserve_pair_counts_qc_and_read_evidence(self):
        paths = self.write_pair(self.geometry_fixture() * 45)
        one = self.run_cli(paths, "--batch-bases", "100")
        four = self.run_cli(paths, "--batch-bases", "100", "--threads", "4")
        self.assertEqual(self.comparable(one), self.comparable(four))
        self.assertEqual(one["evidence_qc"], four["evidence_qc"])
        self.assertEqual(one["paired_end"], four["paired_end"])

    def test_evidence_retains_explicit_zero_pair_counts_and_tsv(self):
        report = self.run_cli(self.write_pair([("TTAGGG" * 20, "CCCTAA" * 20)]), mode="evidence")
        absent = self.row(report, "TTAGGC")
        self.assertEqual(absent["reads"], 0)
        self.assertTrue(absent["paired_evidence"])
        self.assertTrue(all(value == 0 for value in absent["paired_evidence"].values()))
        with self.output.with_suffix(".tsv").open() as handle:
            rows = {row["motif"]: row for row in csv.DictReader(handle, delimiter="\t")}
        for row in self.rows(report):
            for field, value in row["paired_evidence"].items():
                self.assertEqual(rows[row["motif"]]["pair_" + field], str(value))

    def test_valid_slash_illumina_and_partly_unlabelled_headers_with_gzip(self):
        headers = (["a/1", "b 1:N:0:ACGT", "c/1", "d", "e/1 1:Y:0:AG"],
                   ["a/2", "b 2:N:0:ACGT", "c", "d/2", "e/2 2:Y:0:AG"])
        paths = [self.write_fastq("one.data", ["TTAGGG" * 20] * 5, headers[0], compressed=True),
                 self.write_fastq("two.data", ["CCCTAA" * 20] * 5, headers[1], compressed=True)]
        report = self.run_cli(paths)
        self.assertEqual(report["paired_end"]["validated_pairs"], 5)
        self.assertEqual(self.row(report)["paired_evidence"]["opposite_strand_pairs"], 5)

    def test_explicit_wrong_mate_labels_and_conflicting_headers_fail(self):
        for left, right in [("a/1", "a/1"), ("a/2", "a/2"),
                            ("a 2:N:0:X", "a 1:N:0:X"),
                            ("a/1 2:N:0:X", "a/2"),
                            ("a/1", "a/2 1:N:0:X")]:
            with self.subTest(headers=(left, right)):
                paths = self.write_pair([("TTAGGG" * 20, "CCCTAA" * 20)], ([left], [right]))
                self.run_cli(paths, success=False)

    def test_archive_spot_comments_validate_mate_position_without_changing_identity(self):
        pairs = [("TTAGGG" * 20, "CCCTAA" * 20)] * 3
        left = ["SRR21149832.1 1/1", "ERR123.17 17/1", "DRR15.25 25/1"]
        right = ["SRR21149832.1 1/2", "ERR123.17 17/2", "DRR15.25 25/2"]
        report = self.run_cli(self.write_pair(pairs, (left, right)))
        self.assertEqual(report["paired_end"]["validated_pairs"], 3)
        error = self.run_cli(self.write_pair(pairs, (right, left)), success=False)
        self.assertIn("mate label", error.stderr)
        self.run_cli(self.write_pair(pairs, (left, left)), success=False)

    def test_mismatched_and_reordered_pair_ids_fail_without_rematching(self):
        pairs = [("TTAGGG" * 20, "CCCTAA" * 20)] * 2
        for left, right in [(["a/1", "b/1"], ["b/2", "a/2"]),
                            (["SRR1.1", "SRR1.2"], ["SRR1.2", "SRR1.3"])]:
            with self.subTest(headers=(left, right)):
                error = self.run_cli(self.write_pair(pairs, (left, right)), success=False)
                self.assertIn("out of order", error.stderr)

    def test_unequal_streams_and_truncated_second_mate_publish_no_result(self):
        paths = self.write_pair([("TTAGGG" * 20, "CCCTAA" * 20)] * 2)
        paths[1] = self.write_fastq("R2.fastq", ["CCCTAA" * 20])
        error = self.run_cli(paths, "--batch-bases", "1", "--threads", "4", success=False)
        self.assertIn("unequal record counts", error.stderr)
        with paths[1].open("ab") as handle:
            handle.write(b"@read1\nCCCTAACCCTAA\n+\nI\n")
        self.run_cli(paths, "--batch-bases", "1", success=False)

    def test_fasta_even_with_matching_ids_is_not_accepted_as_pairs(self):
        paths = [self.root / "one.fasta", self.root / "two.fasta"]
        for path in paths:
            path.write_text(">read0\n" + "TTAGGG" * 20 + "\n")
        error = self.run_cli(paths, success=False)
        self.assertIn("FASTQ", error.stderr)

    def test_prefix_is_pair_aligned_and_does_not_claim_unread_inputs_validated(self):
        paths = self.write_pair([("TTAGGG" * 20, "CCCTAA" * 20)] * 2,
                                (["a/1", "b/1"], ["a/2", "bad/2"]))
        error = self.run_cli(paths, "--max-records", "1", success=False)
        self.assertIn("even", error.stderr)
        report = self.run_cli(paths, "--max-records", "2")
        self.assertEqual(report["summary"]["records"], 2)
        self.assertTrue(report["summary"]["sampled_prefix"])
        self.assertEqual(report["paired_end"]["validated_pairs"], 1)
        self.assertEqual(report["paired_end"]["scope"], "synchronized_prefix_only")
        self.run_cli(paths, success=False)

    def test_pair_mode_rejects_unsupported_layout_input_count_stdin_and_same_file(self):
        paths = self.write_pair([("TTAGGG" * 20, "CCCTAA" * 20)])
        alias = self.root / "alias.fastq"
        alias.symlink_to(paths[0])
        hardlink = self.root / "hardlink.fastq"
        os.link(paths[0], hardlink)
        for inputs, options in [(paths, ["--paired-layout", "rf"]),
                                (paths[:1], []), (paths + paths[:1], []),
                                ([paths[0], "-"], []),
                                ([paths[0], paths[0]], []), ([paths[0], alias], []),
                                ([paths[0], hardlink], []),
                                (paths, ["--mode", "assembly"])]:
            with self.subTest(inputs=inputs, options=options):
                self.run_cli(inputs, *options, success=False)

    def test_pair_option_is_not_silently_ignored_in_mixed_mode(self):
        paths = self.write_pair([("TTAGGG" * 20, "CCCTAA" * 20)])
        definitions = self.root / "mixed.json"
        definitions.write_text(json.dumps([{"id": "test", "units": ["TTAGGG", "TTAGGC"]}]))
        error = self.run_cli(paths, "--mixed-definitions", str(definitions),
                             mode="mixed", success=False)
        self.assertIn("--paired-layout", error.stderr)


if __name__ == "__main__":
    unittest.main()
