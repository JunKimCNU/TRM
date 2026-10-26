#!/usr/bin/env python3
"""Saved-result reclassification integration tests; no FASTQ inputs are needed."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/reclassify_batch.py"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def canonical(motif: str) -> str:
    reverse = motif.translate(COMPLEMENT)[::-1]
    return min(sequence[index:] + sequence[:index] for sequence in (motif, reverse)
               for index in range(len(sequence)))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def candidate(motif: str, score: float) -> dict:
    motif = canonical(motif)
    return {"motif": motif, "reverse_complement": motif.translate(COMPLEMENT)[::-1],
            "period": len(motif), "distinct_bases": len(set(motif)), "reads": 10,
            "repeat_bases": 1500, "max_tract": 150, "long_tract_reads": 0,
            "forward_reads": 5, "reverse_reads": 5, "ambiguous_strand_reads": 0,
            "left_end_reads": 10, "right_end_reads": 10, "left_forward_reads": 5,
            "left_reverse_reads": 5, "right_forward_reads": 5, "right_reverse_reads": 5,
            "long_left_end_reads": 0, "long_right_end_reads": 0, "opposed_end_records": 0,
            "full_repeat_reads": 10, "status": "repeat_candidate", "rank_score": score,
            "examples": []}


class ReclassifyBatchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.binary = runtime_binary()
        if not cls.binary.is_file():
            raise RuntimeError(f"build the executable first: {cls.binary}")

    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="trm-reclassify-batch-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "source"
        self.source.mkdir()
        self.output = self.root / "reclassified"

    def make_batch(self, motifs: list[str], failed: bool = False) -> Path:
        destination = self.source / "sample1"
        destination.mkdir()
        result = {"schema_version": 1, "program": "trm", "version": "0.1.0", "threads": 1,
                  "command": ["trm", "scan", "/nonexistent/input.fastq.gz"],
                  "inputs": [{"path": "/nonexistent/input.fastq.gz", "bytes": 1000,
                              "modified_unix_seconds": 0}],
                  "config": {"min_k": 3, "max_k": 32, "seed_bases": 24, "min_copies": 3,
                             "min_tract": 500, "min_identity": 0.85, "end_window": 50,
                             "min_reads": 2, "max_candidates": 100000, "mode": "reads"},
                  "summary": {"records": 100, "bases": 15000, "max_record_bases": 150,
                              "candidate_families_before_filter": len(motifs), "reported_candidates": len(motifs),
                              "sampled_prefix": False, "elapsed_seconds": 42.0,
                              "megabases_per_second": 15000 / 42 / 1e6,
                              "call": "candidates_only" if motifs else "no_call"},
                  "warnings": ["original scan warning"],
                  "candidates": [candidate(motif, float(100000 - rank)) for rank, motif in enumerate(motifs)]}
        path = destination / "result.json"
        write_json(path, result)
        metrics = {"operation": "scan", "wall_seconds": 43.0, "exit_code": 0,
                   "input_compressed_bytes": 1000, "user_cpu_seconds": 20.0}
        write_json(destination / "metrics.json", metrics)
        write_json(destination / "state.json", {"status": "complete", "result_sha256": digest(path),
                                                 "metrics": metrics})
        samples = [{"sample_id": "sample1", "status": "complete"}]
        manifest = "sample_id\tpath1\n" + "sample1\t/nonexistent/input.fastq.gz\n"
        if failed:
            failed_dir = self.source / "failed1"
            failed_dir.mkdir()
            write_json(failed_dir / "state.json", {"status": "failed", "error": "unexpected end of file"})
            samples.append({"sample_id": "failed1", "status": "failed"})
            manifest += "failed1\t/nonexistent/truncated.fastq.gz\n"
        (self.source / "manifest.tsv").write_text(manifest, encoding="utf-8")
        write_json(self.source / "summary.json", {"wall_seconds": 45.0, "samples": samples})
        return path

    def run_batch(self, *options: str, success: bool = True, output: Path | None = None) -> subprocess.CompletedProcess:
        command = [sys.executable, str(SCRIPT), str(self.source), "--binary", str(self.binary),
                   "--output-dir", str(output or self.output), "--jobs", "2", *options]
        completed = subprocess.run(command, capture_output=True, text=True, timeout=30, check=False)
        if success:
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        else:
            self.assertNotEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        return completed

    def table(self, name: str) -> list[dict]:
        with (self.output / name).open(encoding="utf-8", newline="") as source:
            return list(csv.DictReader(source, delimiter="\t"))

    def priors(self, rows: list[str]) -> Path:
        path = self.root / "priors.tsv"
        path.write_text("sample_id\tlineage\trelated_motifs\tsource\n" + "\n".join(rows) + "\n",
                        encoding="utf-8")
        return path

    def snapshot_source(self) -> dict:
        return {str(path.relative_to(self.source)): (digest(path), path.stat().st_mtime_ns)
                for path in self.source.rglob("*") if path.is_file()}

    def test_saved_results_only_preserves_source_and_original_failures(self) -> None:
        original_path = self.make_batch(["AAAAAAAG", "TTAGGG"], failed=True)
        snapshot = self.snapshot_source()
        self.run_batch()
        self.assertEqual(self.snapshot_source(), snapshot)
        rows = {row["sample_id"]: row for row in self.table("sample_report.tsv")}
        row = rows["sample1"]
        self.assertEqual(row["status"], "complete")
        self.assertEqual(row["old_top_motif"], canonical("AAAAAAAG"))
        self.assertEqual(row["new_top_motif"], canonical("TTAGGG"))
        self.assertEqual(row["new_top_prior_class"], "exact_known")
        self.assertEqual(row["top_changed"], "True")
        self.assertEqual(row["accepted_candidate_count"], "1")
        self.assertEqual(row["rejected_candidate_count"], "1")
        self.assertEqual(row["original_result_sha256"], digest(original_path))
        self.assertEqual(rows["failed1"]["status"], "source_failed")
        self.assertIn("unexpected end of file", rows["failed1"]["error"])
        self.assertFalse((self.output / "failed1/result.json").exists())
        self.assertEqual([item["motif"] for item in self.table("qualified_candidates.tsv")], [canonical("TTAGGG")])
        self.assertEqual([item["motif"] for item in self.table("rejected_candidates.tsv")], [canonical("AAAAAAAG")])
        metrics = json.loads((self.output / "metrics.json").read_text())
        self.assertFalse(metrics["sequence_inputs_opened"])
        self.assertEqual(metrics["original_scan_batch_wall_seconds"], 45.0)
        self.assertEqual(metrics["original_scan_sum_sample_wall_seconds"], 43.0)
        self.assertEqual(metrics["source_failed_samples"], 1)
        self.assertEqual(metrics["input_artifact_bytes"], original_path.stat().st_size)
        self.assertGreater(metrics["wall_seconds"], 0)
        scan_metrics = json.loads((self.output / "sample1/original_scan_metrics.json").read_text())
        self.assertEqual(scan_metrics["wall_seconds"], 43.0)
        job = json.loads((self.output / "sample1/metrics.json").read_text())
        self.assertEqual(job["exit_code"], 0)
        self.assertIn("reclassify", job["command"])
        self.assertTrue((self.output / "sample1/stderr.log").is_file())

    def test_unsupported_repeats_produce_no_call(self) -> None:
        self.make_batch(["AAAAAAAG", "AAT", "ACGCGGCG", "ACAT", "AAAT"])
        self.run_batch()
        row = self.table("summary.tsv")[0]
        self.assertEqual(row["call"], "no_call")
        self.assertEqual(row["new_top_motif"], "")
        self.assertEqual(row["accepted_candidate_count"], "0")
        self.assertEqual(row["rejected_candidate_count"], "5")
        self.assertEqual(len(self.table("qualified_candidates.tsv")), 0)
        self.assertEqual(len(self.table("rejected_candidates.tsv")), 5)

    def test_related_motif_and_evidence_source_are_recorded(self) -> None:
        motif = "ACGTGCACTGATCGAT"
        self.make_batch([motif])
        priors = self.priors([f"sample1\tExample lineage\t{motif}\thttps://example.org/reference"])
        self.run_batch("--priors", str(priors))
        row = self.table("summary.tsv")[0]
        self.assertEqual(row["new_top_prior_class"], "related_exact")
        self.assertEqual(row["lineage"], "Example lineage")
        self.assertEqual(row["prior_source"], "https://example.org/reference")
        provenance = json.loads((self.output / "sample1/provenance.json").read_text())
        self.assertEqual(provenance["phylogenetic_prior"]["related_motifs"], [motif])
        job = json.loads((self.output / "sample1/metrics.json").read_text())
        self.assertIn("--lineage", job["command"])
        self.assertIn("--related-motif", job["command"])

    def test_duplicate_and_unknown_prior_samples_fail_before_output_creation(self) -> None:
        self.make_batch(["TTAGGG"])
        for rows, error in ((["absent\tExample\tTTAGGG\tref"], "unknown sample_id"),
                            (["sample1\tExample\tTTAGGG\tref"] * 2, "duplicate sample_id")):
            with self.subTest(error=error):
                completed = self.run_batch("--priors", str(self.priors(rows)), success=False)
                self.assertIn(error, completed.stderr)
                self.assertFalse(self.output.exists())

    def test_lineage_without_related_motifs_is_rejected(self) -> None:
        self.make_batch(["TTAGGG"])
        path = self.priors(["sample1\tExample lineage\t\tref"])
        completed = self.run_batch("--priors", str(path), success=False)
        self.assertIn("lineage requires related_motifs", completed.stderr)
        self.assertFalse(self.output.exists())

    def test_malformed_related_motif_is_rejected_by_engine(self) -> None:
        self.make_batch(["TTAGGG"])
        path = self.priors(["sample1\tExample lineage\tACGN\tref"])
        self.run_batch("--priors", str(path), success=False)
        row = self.table("summary.tsv")[0]
        self.assertEqual(row["status"], "reclassification_failed")
        job = json.loads((self.output / "sample1/metrics.json").read_text())
        self.assertNotEqual(job["exit_code"], 0)
        self.assertIn("ACGN", job["command"])
        self.assertTrue((self.output / "sample1/stderr.log").read_text().strip())

    def test_overlapping_output_roots_are_rejected_without_changes(self) -> None:
        self.make_batch(["TTAGGG"])
        snapshot = self.snapshot_source()
        for destination in (self.source, self.source / "new", self.root):
            with self.subTest(output=str(destination)):
                completed = self.run_batch(output=destination, success=False)
                self.assertIn("disjoint", completed.stderr)
        self.assertEqual(self.snapshot_source(), snapshot)

    def test_original_checksum_mismatch_is_retained_as_source_invalid(self) -> None:
        path = self.make_batch(["TTAGGG"])
        path.write_text(path.read_text() + "\n", encoding="utf-8")
        self.run_batch(success=False)
        row = self.table("summary.tsv")[0]
        self.assertEqual(row["status"], "source_invalid")
        self.assertIn("checksum", row["error"])
        self.assertFalse((self.output / "sample1/result.json").exists())

    def test_discovery_keeps_all_candidates_but_tables_show_prior_eligibility(self) -> None:
        self.make_batch(["AAAAAAAG", "TTAGGG"])
        self.run_batch("--prior-policy", "discovery")
        result = json.loads((self.output / "sample1/result.json").read_text())
        self.assertEqual(len(result["candidates"]), 2)
        self.assertEqual(len(self.table("qualified_candidates.tsv")), 1)
        self.assertEqual(len(self.table("rejected_candidates.tsv")), 1)


if __name__ == "__main__":
    unittest.main()
