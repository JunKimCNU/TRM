"""Independent literal-window checks for represented-reference seed opportunity."""

import contextlib
import csv
import gzip
import hashlib
import importlib.util
import io
import itertools
import json
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    "reference_seed_opportunity", Path(__file__).resolve().parents[1] / "scripts/reference_seed_opportunity.py")
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def literal_oracle(sequence, motif, read_length, seed_bases, copies, end_window):
    """Enumerate each read and every contained substring without interval helpers."""
    sequence = sequence.upper()
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    reverse = "".join(complement[base] for base in reversed(motif))
    required = len(motif) * max(copies, (seed_bases + len(motif) - 1) // len(motif))
    units = {part[offset:] + part[:offset] for part in (motif, reverse) for offset in range(len(motif))}
    words = {unit * (required // len(motif)) for unit in units}
    outcome = {scope: {field: 0 for field in module.COUNT_FIELDS}
               for scope in ("all_valid", "ambiguity_free")}
    for start in range(len(sequence) - read_length + 1):
        read = sequence[start:start + read_length]
        contexts = set()
        for offset in range(read_length - required + 1):
            if read[offset:offset + required] in words:
                a = start + offset
                context = "end" if a + required <= min(end_window, len(sequence)) or a >= max(0, len(sequence) - end_window) else "interior"
                contexts.add(context)
        for scope in ("all_valid", "ambiguity_free"):
            if scope == "ambiguity_free" and any(base not in "ACGT" for base in read):
                continue
            counts = outcome[scope]
            counts["valid_starts"] += 1
            counts["positive_starts"] += bool(contexts)
            counts["negative_starts"] += not contexts
            counts["end_proxy_only_starts"] += contexts == {"end"}
            counts["interior_only_starts"] += contexts == {"interior"}
            counts["both_context_starts"] += contexts == {"end", "interior"}
    return outcome


class OpportunityTests(unittest.TestCase):
    def test_exhaustive_short_sequences_and_literal_window_counts(self):
        # All 3^5 strings, including ambiguous boundaries, and every short read length.
        cases = 0
        for letters in itertools.product("ACN", repeat=5):
            sequence = "".join(letters)
            for motif in ("A", "AC"):
                target = module.make_targets([motif], seed_bases=1, min_copies=2)[0]
                for length in range(1, 7):
                    for end in (0, 2, 8):
                        expected = literal_oracle(sequence, motif, length, 1, 2, end)
                        actual = module.opportunities(sequence.encode(), target, length, end)
                        for scope in expected:
                            self.assertEqual(actual[scope], expected[scope], (sequence, motif, length, end, scope))
                        cases += 1
        self.assertEqual(cases, 8748)

    def test_larger_motifs_overlap_orientation_and_ambiguity(self):
        cases = 0
        for motif in ("ACG", "AT", "TTGGGTTAC", "ACTTGGGTTAC"):
            rc = "".join({"A": "T", "C": "G", "G": "C", "T": "A"}[b] for b in reversed(motif))
            for seq in (motif * 7, "N" + motif * 4 + "R" + rc * 4 + "N",
                        "CC" + motif * 4 + "GTC" + motif * 4 + "TT"):
                target = module.make_targets([motif])[0]
                for length in (1, 24, 33, 51, len(seq), len(seq) + 1):
                    for end in (0, 25, 50, len(seq)):
                        expected = literal_oracle(seq, motif, length, 24, 3, end)
                        result = module.opportunities(seq.encode(), target, length, end)
                        for scope in expected:
                            self.assertEqual(result[scope], expected[scope], (motif, length, end, scope))
                        cases += 1
        self.assertEqual(cases, 288)

    def test_terminal_vs_internal_opportunity_formula(self):
        motif, copies, length = "ACG", 5, 20
        target = module.make_targets([motif], seed_bases=1, min_copies=2)[0]
        terminal = (motif * copies + "N" * 80).encode()
        internal = ("N" * 40 + motif * copies + "N" * 40).encode()
        a = module.opportunities(terminal, target, length, 30)
        b = module.opportunities(internal, target, length, 30)
        self.assertEqual(a["all_valid"]["positive_starts"], 15 - 6 + 1)
        self.assertEqual(b["all_valid"]["positive_starts"], 15 + length - 2 * 6 + 1)
        self.assertEqual(a["ambiguity_free"]["positive_starts"], 0)

    def test_rotation_reverse_complement_and_harmonics_deduplicate(self):
        targets = module.make_targets(["TTAGGG", "TAGGGT", "CCCTAA", "TTAGGGTTAGGG"])
        self.assertEqual(len(targets), 1)
        self.assertEqual(targets[0]["period_bp"], 6)
        self.assertEqual(targets[0]["seed_length_bp"], 24)
        self.assertEqual(len(targets[0]["seed_words"]), 12)

    def test_reverse_complement_preserves_all_partition_counts(self):
        seq = b"NACGACGACGACGGTATACGACGACGNN"
        target = module.make_targets(["ACG"], seed_bases=1, min_copies=2)[0]
        rev = seq.translate(bytes.maketrans(b"ACGT", b"TGCA"))[::-1]
        for length in range(1, len(seq) + 2):
            for end in (0, 7, 100):
                a = module.opportunities(seq, target, length, end)
                b = module.opportunities(rev, target, length, end)
                for scope in ("all_valid", "ambiguity_free"):
                    self.assertEqual(a[scope], b[scope])

    def test_overlapping_end_regions_do_not_multiply_windows(self):
        target = module.make_targets(["AT"], seed_bases=1, min_copies=2)[0]
        actual = module.opportunities(b"ATATATATATAT", target, 8, 10)
        self.assertEqual(actual["all_valid"]["positive_starts"], 5)
        self.assertEqual(actual["all_valid"]["end_proxy_only_starts"], 5)
        self.assertEqual(actual["all_valid"]["both_context_starts"], 0)

    def test_ambiguity_between_seeds_has_matching_clean_numerator(self):
        target = module.make_targets(["AC"], seed_bases=1, min_copies=2)[0]
        seq = b"ACACACNACACAC"
        result = module.opportunities(seq, target, 8, 0)
        self.assertEqual(result["all_valid"]["positive_starts"], 6)
        self.assertEqual(result["ambiguity_free"]["valid_starts"], 0)
        self.assertEqual(result["ambiguity_free"]["positive_starts"], 0)

    def test_resource_and_input_guards(self):
        for motifs in ([], ["ACN"], [""], ["A" * 65], [4]):
            with self.assertRaises((ValueError, TypeError)):
                module.make_targets(motifs)
        with self.assertRaises(ValueError):
            module.clean_base_runs(b"ANANAN", max_intervals=2)
        with self.assertRaises(ValueError):
            list(module.fasta_records(io.BytesIO(b">a\nAAAA\n"), 3))
        for payload in (b"AAA\n", b">a\n", b">\nAA\n", b">a\nA C\n", b"@read\nAC\n"):
            with self.assertRaises(ValueError):
                list(module.fasta_records(io.BytesIO(payload), 100))

    def test_cli_gzip_full_eof_hash_case_records_and_no_cross_record_seed(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "source.fa.gz"
            payload = b">left chromosome\nacacAC\nN\n>right\nAC\n"
            source.write_bytes(gzip.compress(payload, mtime=0))
            output = directory / "output"
            with contextlib.redirect_stdout(io.StringIO()):
                status = module.main(["--fasta", str(source), "--motif", "AC", "--read-length", "4",
                                      "--seed-bases", "1", "--min-copies", "2", "--end-window-bp", "4", "--outdir", str(output)])
            self.assertEqual(status, 0)
            summary = json.loads((output / "summary.json").read_text())
            self.assertTrue(summary["full_fasta_eof_verified"])
            self.assertEqual(summary["source"]["records"], 2)
            self.assertEqual(summary["source"]["bases"], 9)
            self.assertEqual(summary["source"]["sha256"], hashlib.sha256(source.read_bytes()).hexdigest())
            with (output / "records.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["lowercase_bases"], "4")
            self.assertEqual(rows[0]["ambiguous_bases"], "1")
            with (output / "opportunities.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            right = [row for row in rows if row["record_id"] == "right"]
            self.assertTrue(all(row["valid_starts"] == "0" and row["positive_starts"] == "0" for row in right))
            self.assertFalse(list(output.glob("*.part")))
            with contextlib.redirect_stderr(io.StringIO()):
                repeated = module.main(["--fasta", str(source), "--motif", "AC", "--read-length", "4", "--outdir", str(output)])
            self.assertEqual(repeated, 1)
            self.assertFalse((output / "failure.json").exists())

    def test_duplicate_record_or_truncated_gzip_is_failure_not_biological_zero(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            cases = (b">same\nACACAC\n>same\nACACAC\n", gzip.compress(b">good\nACACAC\n", mtime=0)[:-5])
            for index, payload in enumerate(cases):
                source = directory / f"bad{index}"
                source.write_bytes(payload)
                output = directory / f"out{index}"
                with contextlib.redirect_stderr(io.StringIO()):
                    status = module.main(["--fasta", str(source), "--motif", "AC", "--read-length", "4", "--outdir", str(output)])
                self.assertEqual(status, 1)
                self.assertFalse((output / "summary.json").exists())
                failure = json.loads((output / "failure.json").read_text())
                self.assertFalse(failure["biological_result_available"])

    def test_targets_json_requires_array_and_primitive_seed_rounding(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "reference.fa"
            source.write_text(">a\n" + "TTGGGTTAC" * 4 + "\n")
            targets = directory / "targets.json"
            targets.write_text(json.dumps({"motifs": "TTGGGTTAC"}))
            with contextlib.redirect_stderr(io.StringIO()):
                status = module.main(["--fasta", str(source), "--targets-json", str(targets),
                                      "--read-length", "30", "--outdir", str(directory / "bad")])
            self.assertEqual(status, 1)
            targets.write_text(json.dumps({"motifs": ["TTGGGTTAC", "ACTTGGGTTAC"]}))
            with contextlib.redirect_stdout(io.StringIO()):
                status = module.main(["--fasta", str(source), "--targets-json", str(targets),
                                      "--read-length", "30", "--outdir", str(directory / "good")])
            self.assertEqual(status, 0)
            summary = json.loads((directory / "good" / "summary.json").read_text())
            self.assertEqual(sorted(item["seed_length_bp"] for item in summary["targets"]), [27, 33])


if __name__ == "__main__":
    unittest.main()
