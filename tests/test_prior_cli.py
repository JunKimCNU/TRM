"""Biological-prior regression tests using sequence inputs and saved evidence."""
import json
from pathlib import Path
import subprocess
import unittest

import test_cli
from test_cli import canonical, fastq, reverse_complement


class PriorCliTests(unittest.TestCase):
    setUpClass = classmethod(test_cli.CliTests.setUpClass.__func__)
    setUp = test_cli.CliTests.setUp
    write = test_cli.CliTests.write
    run_scan = test_cli.CliTests.run_scan
    candidate = test_cli.CliTests.candidate

    def test_user_negative_repeats_are_retained_only_as_rejected(self):
        motifs = ["AAAAAAAG", "AAT", "ACGCGGCG", "ACAT", "AAAT"]
        path = self.write("satellites.fastq", fastq([m * 210 for m in motifs] * 3))
        result, _ = self.run_scan([path])
        self.assertEqual(result["summary"]["call"], "no_call")
        self.assertEqual(result["candidates"], [])
        rejected = {c["motif"]: c for c in result["rejected_candidates"]}
        for motif in motifs:
            row = rejected[canonical(motif)]
            self.assertFalse(row["telomere_prior"]["eligible"])
            self.assertGreaterEqual(row["max_tract"], 500)
            self.assertEqual(row["telomere_prior"]["class"], "uncharacterized_repeat")

    def test_expected_variants_and_known_exception_are_eligible(self):
        motifs = ["TTAGGG", "TTTAGGG", "TTCAGGG", "TTAGGC", "TTAGAC", "TTAGG",
                  "CTCGGTTATGGG"]
        sequences = [m * 100 for m in motifs]
        sequences += [reverse_complement(m * 100) for m in motifs]
        result, _ = self.run_scan([self.write("positives.fastq", fastq(sequences))],
                                  "--prior-policy", "telomere")
        for motif in motifs:
            row = self.candidate(result, motif)
            self.assertTrue(row["telomere_prior"]["eligible"])
        self.assertEqual(self.candidate(result, "TTAGAC")["telomere_prior"]["confidence_tier"],
                         "known_reference")

    def test_large_satellite_does_not_outrank_short_known_motif(self):
        sequences = ["AAAT" * 1000] * 100 + ["TTAGGG" * 25] * 2
        result, _ = self.run_scan([self.write("mixture.fastq", fastq(sequences))])
        self.assertEqual(result["candidates"][0]["motif"], canonical("TTAGGG"))
        self.assertEqual(result["candidates"][0]["long_tract_reads"], 0)
        self.assertTrue(any(c["motif"] == canonical("AAAT") for c in result["rejected_candidates"]))

    def test_explicit_relative_context_can_support_a_divergent_family(self):
        # An invented reference checks mechanics; it is not a curated telomere.
        result, _ = self.run_scan([self.write("relative.fastq", fastq(
            ["ACGACACC" * 100, "TTAGGG" * 100]))],
            "--related-motif", "ACGACACC", "--lineage", "synthetic_relative",
            "--prior-policy", "telomere")
        row = self.candidate(result, "ACGACACC")
        self.assertEqual(row["telomere_prior"]["class"], "related_exact")
        self.assertEqual(row["telomere_prior"]["source_lineage"], "synthetic_relative")
        conflict = next(c for c in result["rejected_candidates"] if c["motif"] == canonical("TTAGGG"))
        self.assertEqual(conflict["telomere_prior"]["class"], "lineage_conflict")

    def test_reclassification_preserves_evidence_without_original_fastq(self):
        path = self.write("old.fastq", fastq(["AAAT" * 200, "TTAGGG" * 100]))
        source = self.root / "source"
        original, _ = self.run_scan([path], "--prior-policy", "discovery", prefix=source)
        original_bytes = source.with_suffix(".json").read_bytes()
        path.unlink()
        output = self.root / "updated"
        completed = subprocess.run([str(self.binary), "reclassify", str(source.with_suffix(".json")),
                                    "--output", str(output)], capture_output=True, timeout=10)
        self.assertEqual(completed.returncode, 0, completed.stderr.decode())
        updated = json.loads(output.with_suffix(".json").read_text())
        self.assertEqual(source.with_suffix(".json").read_bytes(), original_bytes)
        self.assertEqual(updated["summary"]["elapsed_seconds"], original["summary"]["elapsed_seconds"])
        self.assertFalse(updated["reclassification"]["fastq_reread"])
        combined = updated["candidates"] + updated["rejected_candidates"]
        old_rows = {c["motif"]: {k: v for k, v in c.items() if k != "telomere_prior"}
                    for c in original["candidates"]}
        new_rows = {c["motif"]: {k: v for k, v in c.items() if k != "telomere_prior"}
                    for c in combined}
        self.assertEqual(old_rows, new_rows)
        self.assertTrue(Path(str(output) + ".rejected.tsv").is_file())

    def test_prior_invalid_settings_and_source_overwrite_are_rejected(self):
        path = self.write("valid.fastq", fastq(["TTAGGG" * 25]))
        for options in (("--prior-policy", "unknown"), ("--prior-max-edits", "65"),
                        ("--prior-max-distance", "NaN"), ("--prior-max-distance", "1.1"),
                        ("--related-motif", "TTNGGG"), ("--lineage", "unmapped_species")):
            with self.subTest(options=options):
                self.run_scan([path], *options, success=False)
        source = self.root / "saved"
        self.run_scan([path], prefix=source)
        before = source.with_suffix(".json").read_bytes()
        completed = subprocess.run([str(self.binary), "reclassify", str(source.with_suffix(".json")),
                                    "--output", str(source)], capture_output=True, timeout=10)
        self.assertNotEqual(completed.returncode, 0)
        self.assertEqual(source.with_suffix(".json").read_bytes(), before)

    def test_reclassification_preserves_original_decimal_values_exactly(self):
        # Real archived identity value that default serde_json parsing rounded
        # upward by one ULP. Reclassification must not alter measured evidence.
        original = {
            "summary": {"records": 10, "bases": 1000, "elapsed_seconds": 1.2345678901234567},
            "candidates": [{"motif": "AAGCCT", "reads": 3, "rank_score": 223787.48683169714,
                            "examples": [{"identity": 0.9807692307692307}]}],
        }
        source = self.write("precision.json", json.dumps(original).encode())
        output = self.root / "precision_updated"
        completed = subprocess.run([str(self.binary), "reclassify", str(source), "--output", str(output)],
                                   capture_output=True, timeout=10)
        self.assertEqual(completed.returncode, 0, completed.stderr.decode())
        result = json.loads(output.with_suffix(".json").read_text())
        self.assertEqual(result["summary"]["elapsed_seconds"], original["summary"]["elapsed_seconds"])
        candidate = result["candidates"][0]
        self.assertEqual({k: v for k, v in candidate.items() if k != "telomere_prior"},
                         original["candidates"][0])


if __name__ == "__main__":
    unittest.main()
