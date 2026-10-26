"""Saved-artifact taxonomy joins must retain no-calls/failures and reject drift."""

from __future__ import annotations

import copy
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "report_taxonomy.py"
SPEC = importlib.util.spec_from_file_location("report_taxonomy", SCRIPT)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def fixtures() -> dict:
    metadata = [
        {"sample_id": "SRR1", "status": "resolved", "scientific_name": "Caenorhabditis elegans",
         "tax_id": "6239", "lineage": "Eukaryota; Nematoda; Caenorhabditis;",
         "library_strategy": "WGS", "instrument_platform": "ILLUMINA", "error": ""},
        {"sample_id": "SRR2", "status": "resolved", "scientific_name": "Example species",
         "tax_id": "99999", "lineage": "Eukaryota;", "library_strategy": "WGS",
         "instrument_platform": "ILLUMINA", "error": ""},
        {"sample_id": "SRR3", "status": "lookup_failed", "scientific_name": "", "tax_id": "",
         "lineage": "", "library_strategy": "", "instrument_platform": "", "error": "fixture lookup failure"},
    ]
    assignments = []
    for index, meta in enumerate(metadata):
        row = {field: "" for field in module.ASSIGNMENT_FIELDS}
        row.update({field: meta[field] for field in ("sample_id", "scientific_name", "tax_id", "lineage")})
        row.update(metadata_status=meta["status"], mapping_status="global_prior_no_curated_reference",
                   notes="Fixture assignment, not biological evidence.")
        if index == 0:
            row.update(mapping_status="exact_species_reference", related_motifs="TTAGGC",
                       reference_ids="fixture_ce", reference_evidence="experimental",
                       prior_source="https://example.org/reference")
        if index == 1:
            row.update(review_motifs="TTAGGG", review_reference_ids="fixture_computational")
        if index == 2:
            row["mapping_status"] = "global_prior_metadata_unresolved"
        assignments.append(row)
    current = []
    for index in range(3):
        row = {field: "" for field in module.BATCH_FIELDS}
        row.update(sample_id=f"SRR{index + 1}", status="complete", records="10", bases="1500",
                   raw_candidate_count="3", accepted_candidate_count="1", rejected_candidate_count="2",
                   sequence_prior_candidate_count="2", original_scan_wall_seconds="1.5",
                   reclassification_wall_seconds="0.01", call="candidates_only", new_top_motif="AAGCCT",
                   new_top_reads="3", new_top_max_tract="150", new_top_status="repeat_candidate",
                   new_top_prior_class="related_reference", result_json=f"SRR{index + 1}/result.json")
        if index == 0:
            row.update(related_motifs="TTAGGC", prior_source="https://example.org/reference",
                       lineage="Caenorhabditis elegans (TaxID 6239); exact_species_reference; fixture_ce")
        if index == 1:
            row.update(call="no_call", accepted_candidate_count="0", rejected_candidate_count="3",
                       new_top_motif="", new_top_reads="", new_top_max_tract="", new_top_status="",
                       new_top_prior_class="")
        if index == 2:
            row = {field: "" for field in module.BATCH_FIELDS}
            row.update(sample_id="SRR3", status="source_failed", error="truncated FASTQ in original scan",
                       original_scan_wall_seconds="1.8", reclassification_wall_seconds="0.001")
        current.append(row)
    previous = copy.deepcopy(current)
    for row in previous[:2]:
        row.update(new_top_motif="AACCCT", new_top_prior_class="exact_known", related_motifs="",
                   prior_source="", lineage="", reclassification_wall_seconds="0.02", call="candidates_only",
                   accepted_candidate_count="1", rejected_candidate_count="2", new_top_reads="3",
                   new_top_max_tract="150", new_top_status="repeat_candidate")

    def candidate(row: dict, matched: str) -> dict:
        return {"sample_id": row["sample_id"], "rank": "1", "motif": row["new_top_motif"],
                "reads": row["new_top_reads"], "max_tract": row["new_top_max_tract"],
                "status": row["new_top_status"], "prior_class": row["new_top_prior_class"],
                "prior_matched_reference": matched}

    return {"metadata": metadata, "assignments": assignments, "current": current, "previous": previous,
            "current_candidates": [candidate(current[0], "TTAGGC")],
            "previous_candidates": [candidate(row, "TTAGGG") for row in previous[:2]]}


def join(data: dict) -> list[dict]:
    return module.joined_rows(*(module.index_samples(data[key], key)
                                for key in ("metadata", "assignments", "current", "previous")),
                              data["current_candidates"], data["previous_candidates"])


def write_table(path: Path, rows: list[dict], fields: set[str]) -> None:
    columns = sorted(fields | {field for row in rows for field in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class TaxonomyReportTests(unittest.TestCase):
    def test_join_preserves_identity_review_only_no_call_and_failed_scan(self) -> None:
        rows = join(fixtures())
        self.assertEqual([row["sample_id"] for row in rows], ["SRR1", "SRR2", "SRR3"])
        first, no_call, failed = rows
        self.assertEqual(first["scientific_name"], "Caenorhabditis elegans")
        self.assertEqual(first["tax_id"], "6239")
        self.assertEqual(first["top_motif"], "AAGCCT")
        self.assertEqual(first["top_matched_reference"], "TTAGGC")
        self.assertEqual(first["previous_top_matched_reference"], "TTAGGG")
        self.assertEqual(first["top_changed_vs_previous"], "true")
        self.assertEqual(first["curated_expected_motifs"], "TTAGGC")
        self.assertEqual(first["instrument_platform"], "ILLUMINA")
        self.assertEqual(first["original_scan_wall_seconds"], "1.5")
        self.assertEqual(first["reclassification_wall_seconds"], "0.01")
        self.assertEqual(no_call["no_call"], "true")
        self.assertEqual(no_call["curated_expected_motifs"], "")
        self.assertEqual(no_call["review_only_motifs"], "TTAGGG")
        self.assertEqual(no_call["top_motif"], "")
        self.assertEqual(no_call["top_matched_reference"], "")
        self.assertEqual(failed["analysis_status"], "source_failed")
        self.assertEqual(failed["metadata_error"], "fixture lookup failure")
        self.assertIn("truncated", failed["error"])
        for field in ("call", "no_call", "top_changed_vs_previous", "records", "top_motif"):
            self.assertEqual(failed[field], "")

    def test_observed_motif_never_changes_unresolved_metadata(self) -> None:
        data = fixtures()
        for table in (data["metadata"], data["assignments"]):
            table[0].update(scientific_name="", tax_id="", lineage="")
        data["metadata"][0]["status"] = "lookup_failed"
        data["assignments"][0].update(metadata_status="lookup_failed", related_motifs="", prior_source="",
                                       mapping_status="global_prior_metadata_unresolved")
        data["current"][0].update(related_motifs="", prior_source="")
        row = join(data)[0]
        self.assertEqual(row["top_motif"], "AAGCCT")
        self.assertEqual(row["scientific_name"], "")
        self.assertEqual(row["tax_id"], "")

    def test_archive_sources_and_raw_taxonomy_stay_separate_from_resolved_identity(self) -> None:
        data = fixtures()
        data["metadata"][0].update(
            metadata_source="https://example.org/ena/run", taxonomy_source="https://example.org/taxonomy/6239",
            taxonomy_assignment_source="https://example.org/ncbi/sample", deposited_tax_id="1000",
            ena_tax_id="99999", ena_scientific_name="historical archive name", requested_tax_id="1000",
            tax_id_resolution="explicit_merge")
        row = join(data)[0]
        self.assertEqual(row["tax_id"], "6239")
        self.assertEqual(row["scientific_name"], "Caenorhabditis elegans")
        self.assertEqual(row["deposited_tax_id"], "1000")
        self.assertEqual(row["ena_tax_id"], "99999")
        self.assertEqual(row["ena_scientific_name"], "historical archive name")
        self.assertEqual(row["tax_id_resolution"], "explicit_merge")
        self.assertEqual(row["metadata_source"], "https://example.org/ena/run")
        self.assertEqual(row["taxonomy_source"], "https://example.org/taxonomy/6239")
        self.assertEqual(row["taxonomy_assignment_source_url"], "https://example.org/ncbi/sample")

    def test_missing_extra_or_duplicate_samples_in_any_sample_table_fail(self) -> None:
        for table in ("metadata", "assignments", "current", "previous"):
            for mutation in ("missing", "extra", "duplicate"):
                with self.subTest(table=table, mutation=mutation):
                    data = fixtures()
                    if mutation == "missing":
                        data[table].pop()
                    else:
                        extra = dict(data[table][0])
                        if mutation == "extra":
                            extra["sample_id"] = "SRR4"
                        data[table].append(extra)
                    with self.assertRaises(ValueError):
                        join(data)

    def test_assignment_metadata_identity_and_status_must_agree(self) -> None:
        for field in ("scientific_name", "tax_id", "lineage", "metadata_status"):
            with self.subTest(field=field):
                data = fixtures()
                data["assignments"][0][field] = "different"
                with self.assertRaisesRegex(ValueError, "disagrees with metadata"):
                    join(data)

    def test_actual_applied_prior_must_agree_with_assignment(self) -> None:
        for field in ("related_motifs", "prior_source"):
            data = fixtures()
            data["current"][0][field] = "wrong"
            with self.assertRaisesRegex(ValueError, "disagrees with assignment"):
                join(data)

    def test_original_measurements_must_be_preserved(self) -> None:
        for field in ("records", "bases", "original_scan_wall_seconds"):
            data = fixtures()
            data["current"][0][field] = "99"
            with self.assertRaisesRegex(ValueError, "changed between reclassifications"):
                join(data)

    def test_missing_duplicate_extra_and_failed_sample_candidates_fail(self) -> None:
        for table in ("current_candidates", "previous_candidates"):
            for mutation in ("missing", "duplicate", "extra", "failed"):
                with self.subTest(table=table, mutation=mutation):
                    data = fixtures()
                    if mutation == "missing":
                        data[table].pop(0)
                    else:
                        extra = dict(data[table][0])
                        if mutation == "extra":
                            extra["sample_id"] = "SRR4"
                        elif mutation == "failed":
                            extra["sample_id"] = "SRR3"
                        data[table].append(extra)
                    with self.assertRaises(ValueError):
                        join(data)

    def test_top_candidate_must_match_summary_and_call(self) -> None:
        for field in ("motif", "reads", "max_tract", "status", "prior_class"):
            data = fixtures()
            data["current_candidates"][0][field] = "wrong"
            with self.assertRaisesRegex(ValueError, "disagrees with sample report"):
                join(data)
        data = fixtures()
        data["current"][0]["call"] = "no_call"
        with self.assertRaisesRegex(ValueError, "call disagrees"):
            join(data)

    def test_unfinished_batch_inconsistent_counts_and_invalid_timing_fail(self) -> None:
        for changes in ({"status": "running"}, {"accepted_candidate_count": "2"},
                        {"sequence_prior_candidate_count": "0"}, {"records": "0"},
                        {"reclassification_wall_seconds": "NaN"}):
            data = fixtures()
            data["current"][0].update(changes)
            with self.assertRaises(ValueError):
                join(data)

    def make_cli(self, root: Path) -> tuple[list[str], Path]:
        data = fixtures()
        write_table(root / "metadata.tsv", data["metadata"], module.METADATA_FIELDS)
        write_table(root / "assignments.tsv", data["assignments"], module.ASSIGNMENT_FIELDS)
        for label in ("current", "previous"):
            directory = root / label
            directory.mkdir()
            write_table(directory / "sample_report.tsv", data[label], module.BATCH_FIELDS)
            write_table(directory / "qualified_candidates.tsv", data[label + "_candidates"], module.CANDIDATE_FIELDS)
        output = root / "taxon_sample_report.tsv"
        command = [sys.executable, str(SCRIPT), "--metadata", str(root / "metadata.tsv"),
                   "--assignments", str(root / "assignments.tsv"), "--batch", str(root / "current"),
                   "--previous-batch", str(root / "previous"), "--output", str(output)]
        return command, output

    def test_cli_hashes_metrics_and_no_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            command, output = self.make_cli(root)
            completed = subprocess.run(command, text=True, capture_output=True, timeout=10)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            sidecar = output.with_suffix(".summary.json")
            summary = json.loads(sidecar.read_text())
            self.assertEqual(summary["sample_count"], 3)
            self.assertEqual(summary["analysis_status_counts"], {"complete": 2, "source_failed": 1})
            self.assertEqual(summary["no_call_samples"], 1)
            self.assertEqual(summary["top_changed_samples"], 2)
            self.assertEqual(summary["qualified_candidate_count"], 1)
            self.assertEqual(summary["output_sha256"], module.sha256(output))
            self.assertEqual(summary["inputs"]["metadata"]["sha256"], module.sha256(root / "metadata.tsv"))
            self.assertFalse(summary["sequence_inputs_opened"])
            self.assertFalse(summary["motif_results_used_for_taxonomy"])
            before = (output.read_bytes(), sidecar.read_bytes())
            self.assertNotEqual(subprocess.run(command, capture_output=True, timeout=10).returncode, 0)
            self.assertEqual(before, (output.read_bytes(), sidecar.read_bytes()))
            self.assertFalse(list(root.glob(".*.tsv.*")))

    def test_existing_sidecar_prevents_creating_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            command, output = self.make_cli(root)
            sidecar = output.with_suffix(".summary.json")
            sidecar.write_text("keep me")
            self.assertNotEqual(subprocess.run(command, capture_output=True, timeout=10).returncode, 0)
            self.assertFalse(output.exists())
            self.assertEqual(sidecar.read_text(), "keep me")

    def test_invalid_tsv_headers_or_row_width_fail_before_output(self) -> None:
        for text in ("sample_id\tsample_id\nSRR1\tSRR1\n", "sample_id\nSRR1\n"):
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                command, output = self.make_cli(root)
                (root / "metadata.tsv").write_text(text)
                self.assertNotEqual(subprocess.run(command, capture_output=True, timeout=10).returncode, 0)
                self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
