"""Independent frozen exact-grammar fixtures and mixed CLI contract checks."""
import hashlib
import gzip
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests/fixtures/mixed"
VALIDATION = FIXTURES / "validation"
DATA = FIXTURES / "data"
spec = importlib.util.spec_from_file_location("mixed_reference_oracle", VALIDATION / "reference_oracle.py")
oracle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(oracle)


def frozen_fixture_path(original_path):
    """Relocate original pinned paths without changing receipts or assertions."""
    original = Path(original_path)
    locations = {
        'data/mixed_validation_20260923': DATA,
        'results/mixed_validation_20260923': VALIDATION,
    }
    parent = original.parent.as_posix()
    if parent not in locations:
        raise ValueError(f'Unregistered frozen fixture path: {original_path}')
    return locations[parent] / original.name


class MixedPreregistrationTests(unittest.TestCase):
    def test_frozen_fixture_and_oracle_hashes(self):
        for name, key in (("preregistered_settings.json", "fixture_hashes"),
                          ("api_validation_addendum.json", "hashes")):
            doc = json.loads((VALIDATION / name).read_text())
            for rel, digest in doc[key].items():
                with self.subTest(file=rel):
                    self.assertEqual(hashlib.sha256(frozen_fixture_path(rel).read_bytes()).hexdigest(), digest)

    def test_fixture_oracle_reproducibility(self):
        groups = json.loads((DATA / "definitions.json").read_text())
        groups += json.loads((DATA / "noncomposite_ambiguous_4_6_8.definitions.json").read_text())
        units = {row["id"]: row["units"] for row in groups}
        fixtures = json.loads((DATA / "fixtures.json").read_text())
        fixtures += json.loads((DATA / "api_addendum_fixtures.json").read_text())
        for row in fixtures:
            with self.subTest(record=row["id"]):
                self.assertEqual(oracle.find_mixed(row["sequence"], units[row["group"]]), row["oracle_intervals"])


class MixedCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binary = runtime_binary()

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="trm-mixed-test-")
        self.addCleanup(temporary.cleanup)
        self.temporary = Path(temporary.name)
        self.counter = 0

    def run_cli(self, fasta, definitions, *options, success=True):
        self.counter += 1
        output = self.temporary / f"result{self.counter}"
        command = [str(self.binary), "mixed", str(fasta), "--mixed-definitions", str(definitions),
                   "--output", str(output), "--threads", "1", "--mixed-max-examples", "10000", *options]
        completed = subprocess.run(command, text=True, capture_output=True, timeout=30)
        if not success:
            self.assertNotEqual(completed.returncode, 0, completed.stderr)
            self.assertFalse(output.with_suffix(".json").exists())
            return completed
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return json.loads(output.with_suffix(".json").read_text())

    def test_invalid_atomic_definitions_are_rejected(self):
        for units in (["ACG", "TCA", "ACGTCA"], ["ACG", "GAC"],
                      ["ACG", "CGT"], ["ACGACG", "TTCAGG"], ["ACN", "GTC"]):
            with self.subTest(units=units):
                definitions = self.temporary / "invalid.json"
                definitions.write_text(json.dumps([{"id": "invalid", "units": units}]))
                completed = self.run_cli(self.temporary / "missing.fasta", definitions, success=False)
                self.assertNotIn("failed to open", completed.stderr.lower())

    def test_mixed_mode_does_not_silently_accept_homogeneous_discovery_options(self):
        for options in (("--min-k", "3"), ("--seed-bases", "24"), ("--prior-policy", "canonical")):
            with self.subTest(options=options):
                self.run_cli(DATA / "sosm_9_11.fasta", DATA / "sosm_9_11.definitions.json",
                             *options, success=False)

    def test_134_frozen_valid_sequences_match_independent_oracle(self):
        fixtures = json.loads((DATA / "fixtures.json").read_text())
        fixtures += json.loads((DATA / "api_addendum_fixtures.json").read_text())
        valid = [row for row in fixtures if row["group"] != "ambiguous_3_3_6"]
        self.assertEqual(len(valid), 134)
        for fixture in valid:
            with self.subTest(record=fixture["id"]):
                fasta = self.temporary / "one.fasta"
                fasta.write_text(f'>{fixture["id"]}\n{fixture["sequence"]}\n')
                doc = self.run_cli(fasta, DATA / (fixture["group"] + ".definitions.json"))
                group = doc["mixed_arrays"]["groups"][0]
                expected = [hit for hit in fixture["oracle_intervals"] if hit["guaranteed_mixed"]]
                self.assertEqual(group["mixed_records"], int(bool(expected)))
                self.assertEqual(group["max_array_bp"], max((hit["length"] for hit in expected), default=0))
                self.assertEqual(group["max_any_array_bp"], max((hit["length"] for hit in fixture["oracle_intervals"]), default=0))
                expected_intervals = {(hit["start"], hit["end"]) for hit in expected}
                self.assertEqual({(hit["start"], hit["end"]) for hit in group["examples"]}, expected_intervals)

    def test_record_geometry_keeps_expected_and_opposite_polarity_separate(self):
        flank = "ACGTCCGATGCTAACGTCAGTGCATCAGACGTTCGACGTA"
        repeat = ("TTGGGTTAC" + "ACTTGGGTTAC") * 4
        reverse = oracle.reverse_complement(repeat)
        seqs = [repeat, reverse, flank + repeat, reverse + flank,
                repeat + flank, flank + reverse, flank + repeat + flank]
        fasta = self.temporary / "geometry.fasta"
        fasta.write_text("".join(f">record{i}\n{seq}\n" for i, seq in enumerate(seqs)))
        group = self.run_cli(fasta, DATA / "sosm_9_11.definitions.json")["mixed_arrays"]["groups"][0]
        self.assertEqual(group["mixed_records"], 7)
        self.assertEqual(group["full_mixed_records"], 2)
        self.assertEqual(group["entire_mixed_records"], 2)
        for field in ("g_right_boundary_records", "c_left_boundary_records",
                      "g_left_boundary_records", "c_right_boundary_records"):
            self.assertEqual(group[field], 1, (field, group))
        self.assertEqual(group["internal_array_records"], 1)
        self.assertEqual(group["clean_two_flank_records"], 1)

    def test_threads_preserve_groups_and_bounded_examples(self):
        fasta = self.temporary / "duplicate_ids.fasta"
        a, b = "TTGGGTTAC", "ACTTGGGTTAC"
        seqs = [a * 4 + b * 13, a * 15 + b * 4] * 120
        self.assertEqual(len(seqs[0]), len(seqs[1]))
        fasta.write_text("".join(f">same_id\n{seq}\n" for seq in seqs))
        definitions = DATA / "sosm_9_11.definitions.json"
        one = self.run_cli(fasta, definitions, "--batch-bases", "20", "--mixed-max-examples", "5")
        four = self.run_cli(fasta, definitions, "--threads", "4", "--batch-bases", "20", "--mixed-max-examples", "5")
        self.assertEqual(one["mixed_arrays"]["groups"], four["mixed_arrays"]["groups"])
        self.assertEqual(one["evidence_qc"], four["evidence_qc"])

    def test_malformed_gzip_and_truncated_fastq_publish_no_result(self):
        for name, content in (("corrupt.fastq.gz", b"\x1f\x8bdefinitely-not-a-gzip-stream"),
                              ("truncated.fastq.gz", gzip.compress(b"@read\nACGTACGT\n+\nII\n"))):
            with self.subTest(name=name):
                path = self.temporary / name
                path.write_bytes(content)
                self.run_cli(path, DATA / "sosm_9_11.definitions.json", success=False)
                self.assertFalse((self.temporary / f"result{self.counter}.tsv").exists())

    def test_definitions_cannot_be_overwritten_by_output_prefix(self):
        definitions = self.temporary / "protected.json"
        definitions.write_bytes((DATA / "sosm_9_11.definitions.json").read_bytes())
        before = definitions.read_bytes()
        completed = self.run_cli(DATA / "sosm_9_11.fasta", definitions,
                                 "--output", str(definitions.with_suffix("")), success=False)
        self.assertEqual(definitions.read_bytes(), before)
        self.assertIn("overwrite", completed.stderr.lower())

    def test_entire_read_and_input_files_used_and_prefix_marked(self):
        repeat = ("TTGGGTTAC" + "ACTTGGGTTAC") * 4
        first = self.temporary / "first.fasta"
        first.write_text(">noise\n" + "N" * 300 + "\n")
        second = self.temporary / "second.fasta"
        second.write_text(">late\n" + "N" * 200 + repeat + "\n")
        definitions = DATA / "sosm_9_11.definitions.json"
        all_reads = self.run_cli(first, definitions, str(second))
        prefix = self.run_cli(first, definitions, str(second), "--max-records", "1")
        self.assertEqual(all_reads["summary"]["records"], 2)
        self.assertEqual(all_reads["mixed_arrays"]["groups"][0]["mixed_records"], 1)
        self.assertEqual(all_reads["mixed_arrays"]["groups"][0]["examples"][0]["start"], 200)
        self.assertFalse(all_reads["summary"]["sampled_prefix"])
        self.assertEqual(prefix["summary"]["records"], 1)
        self.assertEqual(prefix["mixed_arrays"]["groups"][0]["mixed_records"], 0)
        self.assertTrue(prefix["summary"]["sampled_prefix"])


if __name__ == "__main__":
    unittest.main()
