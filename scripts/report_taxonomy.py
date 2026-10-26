#!/usr/bin/env python3
"""Join archive taxonomy and TRM reclassification tables without reading FASTQ.

All sample sets must match, including failed scans. Expected motifs come only
from curated assignments; observed motifs never determine taxonomic identity.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import tempfile
import time


METADATA_FIELDS = {"sample_id", "status", "scientific_name", "tax_id", "lineage"}
ASSIGNMENT_FIELDS = {"sample_id", "metadata_status", "scientific_name", "tax_id", "lineage",
                     "mapping_status", "related_motifs", "reference_ids", "reference_evidence",
                     "prior_source", "review_reference_ids", "review_motifs", "notes"}
BATCH_FIELDS = {"sample_id", "status", "new_top_motif", "new_top_reads", "new_top_max_tract",
                "new_top_status", "new_top_prior_class", "call", "raw_candidate_count",
                "accepted_candidate_count", "rejected_candidate_count", "sequence_prior_candidate_count",
                "records", "bases", "related_motifs", "prior_source", "lineage",
                "reclassification_wall_seconds", "original_scan_wall_seconds", "result_json", "error"}
CANDIDATE_FIELDS = {"sample_id", "rank", "motif", "reads", "max_tract", "status", "prior_class",
                    "prior_matched_reference"}
OUTPUT_FIELDS = [
    "sample_id", "scientific_name", "tax_id", "metadata_status", "lineage", "mapping_status",
    "curated_expected_motifs", "curated_reference_ids", "curated_reference_evidence",
    "curated_prior_source", "review_only_motifs", "review_only_reference_ids", "mapping_notes",
    "library_strategy", "library_source", "library_layout", "instrument_platform", "instrument_model",
    "run_accession", "sample_accession", "study_accession", "taxonomy_scientific_name", "taxonomy_rank",
    "metagenome", "metadata_source", "taxonomy_source", "metadata_error",
    "deposited_tax_id", "ena_scientific_name", "ena_tax_id", "requested_tax_id", "tax_id_resolution",
    "taxonomy_assignment_source_url",
    "analysis_status", "previous_analysis_status", "call", "no_call", "previous_call",
    "top_motif", "top_matched_reference", "top_prior_class", "top_status", "top_reads",
    "top_max_tract", "top_support_status", "top_reads_per_million", "top_read_enrichment_lower",
    "top_repeat_base_enrichment", "previous_top_motif", "previous_top_matched_reference",
    "top_changed_vs_previous", "records", "bases", "raw_candidate_count",
    "sequence_prior_candidate_count", "accepted_candidate_count", "rejected_candidate_count",
    "previous_accepted_candidate_count", "applied_related_motifs", "applied_prior_lineage",
    "applied_prior_source", "original_scan_wall_seconds", "reclassification_wall_seconds",
    "previous_reclassification_wall_seconds", "result_json", "previous_result_json", "error",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        if len(columns) != len(set(columns)) or not required.issubset(columns):
            raise ValueError(f"{path}: missing required or duplicate columns")
        rows = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{reader.line_num}: incorrect field count")
            rows.append({key: value.strip() for key, value in row.items()})
    return rows


def index_samples(rows: list[dict[str, str]], label: str) -> dict[str, dict[str, str]]:
    indexed = {}
    for row in rows:
        sample = row["sample_id"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", sample) or sample in indexed:
            raise ValueError(f"{label}: invalid or duplicate sample_id: {sample!r}")
        indexed[sample] = row
    if not indexed:
        raise ValueError(f"{label}: no samples")
    return indexed


def same_samples(expected: dict, observed: dict, label: str) -> None:
    if expected.keys() != observed.keys():
        missing = sorted(expected.keys() - observed.keys())
        extra = sorted(observed.keys() - expected.keys())
        raise ValueError(f"{label}: sample set mismatch; missing={missing}, extra={extra}")


def integer(value: str, label: str, minimum: int = 0) -> int:
    if not re.fullmatch(r"[0-9]+", value) or int(value) < minimum:
        raise ValueError(f"{label}: invalid integer {value!r}")
    return int(value)


def elapsed(value: str, label: str) -> None:
    if value and (not math.isfinite(float(value)) or float(value) < 0):
        raise ValueError(f"{label}: invalid elapsed time {value!r}")


def validate_batch(samples: dict, candidates: list[dict], label: str) -> dict[str, dict]:
    grouped: dict[str, dict[int, dict]] = {}
    for row in candidates:
        sample = row["sample_id"]
        if sample not in samples:
            raise ValueError(f"{label}: candidate belongs to unknown sample {sample!r}")
        if samples[sample]["status"] != "complete":
            raise ValueError(f"{label}: candidate belongs to incomplete sample {sample!r}")
        rank = integer(row["rank"], f"{label}/{sample}/rank", 1)
        ranks = grouped.setdefault(sample, {})
        if rank in ranks or any(other["motif"] == row["motif"] for other in ranks.values()):
            raise ValueError(f"{label}/{sample}: duplicate candidate rank or motif")
        ranks[rank] = row
    tops = {}
    for sample, row in samples.items():
        if row["status"] not in {"complete", "source_failed", "failed"}:
            raise ValueError(f"{label}/{sample}: unfinished or unrecognized status {row['status']!r}")
        elapsed(row["original_scan_wall_seconds"], f"{label}/{sample}/scan time")
        elapsed(row["reclassification_wall_seconds"], f"{label}/{sample}/reclassification time")
        if row["status"] != "complete":
            if row["call"] or row["new_top_motif"]:
                raise ValueError(f"{label}/{sample}: incomplete sample must not contain a call")
            continue
        raw = integer(row["raw_candidate_count"], f"{label}/{sample}/raw count")
        accepted = integer(row["accepted_candidate_count"], f"{label}/{sample}/accepted count")
        rejected = integer(row["rejected_candidate_count"], f"{label}/{sample}/rejected count")
        prior = integer(row["sequence_prior_candidate_count"], f"{label}/{sample}/prior count")
        integer(row["records"], f"{label}/{sample}/records", 1)
        integer(row["bases"], f"{label}/{sample}/bases", 1)
        if accepted + rejected != raw or not accepted <= prior <= raw:
            raise ValueError(f"{label}/{sample}: inconsistent candidate counts")
        ranks = grouped.get(sample, {})
        if set(ranks) != set(range(1, accepted + 1)):
            raise ValueError(f"{label}/{sample}: missing or inconsistent candidate ranks/counts")
        expected_call = "candidates_only" if accepted else "no_call"
        if row["call"] != expected_call:
            raise ValueError(f"{label}/{sample}: call disagrees with candidate count")
        if not accepted:
            if row["new_top_motif"] or row["new_top_reads"]:
                raise ValueError(f"{label}/{sample}: no_call must not contain a top candidate")
            continue
        top = ranks[1]
        for field in ("motif", "reads", "max_tract", "status", "prior_class"):
            if top[field] != row["new_top_" + field]:
                raise ValueError(f"{label}/{sample}: rank 1 {field} disagrees with sample report")
        tops[sample] = top
    return tops


def joined_rows(metadata: dict, assignments: dict, current: dict, previous: dict,
                current_candidates: list[dict], previous_candidates: list[dict]) -> list[dict]:
    for label, indexed in (("assignments", assignments), ("batch", current), ("previous batch", previous)):
        same_samples(metadata, indexed, label)
    top_current = validate_batch(current, current_candidates, "batch")
    top_previous = validate_batch(previous, previous_candidates, "previous batch")
    output = []
    for sample, meta in metadata.items():
        assigned, new, old = assignments[sample], current[sample], previous[sample]
        for field in ("scientific_name", "tax_id", "lineage"):
            if assigned[field] != meta[field]:
                raise ValueError(f"{sample}: assignment {field} disagrees with metadata")
        if assigned["metadata_status"] != meta["status"]:
            raise ValueError(f"{sample}: assignment metadata_status disagrees with metadata")
        if meta["status"] == "resolved":
            integer(meta["tax_id"], f"{sample}/TaxID", 1)
            if not meta["scientific_name"]:
                raise ValueError(f"{sample}: resolved metadata has no scientific name")
        complete = new["status"] == "complete"
        comparable = complete and old["status"] == "complete"
        if complete:
            for field in ("related_motifs", "prior_source"):
                if new[field] != assigned[field]:
                    raise ValueError(f"{sample}: applied {field} disagrees with assignment")
        if comparable:
            for field in ("records", "bases", "raw_candidate_count", "original_scan_wall_seconds"):
                if new[field] != old[field]:
                    raise ValueError(f"{sample}: original {field} changed between reclassifications")
        row = {field: "" for field in OUTPUT_FIELDS}
        row.update({key: meta.get(key, "") for key in OUTPUT_FIELDS if key in meta})
        row.update({
            "metadata_status": meta["status"], "metadata_error": meta.get("error", ""),
            "taxonomy_assignment_source_url": meta.get("taxonomy_assignment_source", ""),
            "mapping_status": assigned["mapping_status"],
            "curated_expected_motifs": assigned["related_motifs"],
            "curated_reference_ids": assigned["reference_ids"],
            "curated_reference_evidence": assigned["reference_evidence"],
            "curated_prior_source": assigned["prior_source"],
            "review_only_motifs": assigned["review_motifs"],
            "review_only_reference_ids": assigned["review_reference_ids"], "mapping_notes": assigned["notes"],
            "analysis_status": new["status"], "previous_analysis_status": old["status"],
            "call": new["call"], "no_call": str(new["call"] == "no_call").lower() if complete else "",
            "previous_call": old["call"], "top_motif": new["new_top_motif"],
            "top_matched_reference": top_current.get(sample, {}).get("prior_matched_reference", ""),
            "previous_top_motif": old["new_top_motif"],
            "previous_top_matched_reference": top_previous.get(sample, {}).get("prior_matched_reference", ""),
            "top_changed_vs_previous": str(new["new_top_motif"] != old["new_top_motif"]).lower() if comparable else "",
            "previous_accepted_candidate_count": old["accepted_candidate_count"],
            "applied_related_motifs": new["related_motifs"], "applied_prior_lineage": new["lineage"],
            "applied_prior_source": new["prior_source"],
            "previous_reclassification_wall_seconds": old["reclassification_wall_seconds"],
            "result_json": new["result_json"], "previous_result_json": old["result_json"], "error": new["error"],
        })
        for field in ("prior_class", "status", "reads", "max_tract", "support_status", "reads_per_million",
                      "read_enrichment_lower", "repeat_base_enrichment"):
            row["top_" + field] = new.get("new_top_" + field, "")
        for field in ("records", "bases", "raw_candidate_count", "sequence_prior_candidate_count",
                      "accepted_candidate_count", "rejected_candidate_count", "original_scan_wall_seconds",
                      "reclassification_wall_seconds"):
            row[field] = new[field]
        output.append(row)
    return output


def atomic_create(path: Path, contents: str) -> None:
    descriptor, temporary = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(contents)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, path)  # Atomic no-clobber publication, including races.
    finally:
        os.unlink(temporary)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--previous-batch", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    try:
        sidecar = args.output.with_suffix(".summary.json")
        if args.output == sidecar or args.output.exists() or sidecar.exists():
            raise ValueError("output or summary sidecar already exists; choose a new output path")
        sources = {
            "metadata": args.metadata, "assignments": args.assignments,
            "batch_samples": args.batch / "sample_report.tsv",
            "batch_candidates": args.batch / "qualified_candidates.tsv",
            "previous_samples": args.previous_batch / "sample_report.tsv",
            "previous_candidates": args.previous_batch / "qualified_candidates.tsv",
        }
        hashes = {name: sha256(path) for name, path in sources.items()}
        tables = {name: read_table(path, METADATA_FIELDS if name == "metadata" else
                                  ASSIGNMENT_FIELDS if name == "assignments" else
                                  CANDIDATE_FIELDS if name.endswith("candidates") else BATCH_FIELDS)
                  for name, path in sources.items()}
        indexed = {name: index_samples(rows, name) for name, rows in tables.items()
                   if not name.endswith("candidates")}
        rows = joined_rows(indexed["metadata"], indexed["assignments"], indexed["batch_samples"],
                           indexed["previous_samples"], tables["batch_candidates"], tables["previous_candidates"])
        if any(sha256(path) != hashes[name] for name, path in sources.items()):
            raise ValueError("an input artifact changed while generating the report")
        content = io.StringIO(newline="")
        writer = csv.DictWriter(content, OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        summary = {
            "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
            "script_sha256": sha256(Path(__file__)),
            "inputs": {name: {"path": str(path.resolve()), "sha256": hashes[name]}
                       for name, path in sources.items()},
            "output_path": str(args.output.resolve()),
            "output_sha256": hashlib.sha256(content.getvalue().encode("utf-8")).hexdigest(),
            "sample_count": len(rows), "analysis_status_counts": dict(Counter(row["analysis_status"] for row in rows)),
            "metadata_status_counts": dict(Counter(row["metadata_status"] for row in rows)),
            "mapping_status_counts": dict(Counter(row["mapping_status"] for row in rows)),
            "scientific_name_counts": dict(Counter(row["scientific_name"] for row in rows)),
            "no_call_samples": sum(row["no_call"] == "true" for row in rows),
            "top_changed_samples": sum(row["top_changed_vs_previous"] == "true" for row in rows),
            "qualified_candidate_count": sum(int(row["accepted_candidate_count"]) for row in rows
                                             if row["analysis_status"] == "complete"),
            "sequence_inputs_opened": False, "motif_results_used_for_taxonomy": False,
            "wall_seconds": time.perf_counter() - started,
            "limitations": [
                "Scientific names and TaxIDs are archive metadata, not sequence-based identity verification.",
                "Curated expected motifs are prior references, not confirmation of this sample's telomere location.",
                "Review-only motifs did not set the applied lineage prior.",
                "No call does not demonstrate motif absence or a telomere-maintenance mechanism.",
                "Per-sample scan and reclassification times are separate; sums do not equal parallel batch wall time.",
            ],
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        atomic_create(args.output, content.getvalue())
        atomic_create(sidecar, json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
        print(json.dumps({key: summary[key] for key in ("sample_count", "analysis_status_counts",
                                                       "no_call_samples", "top_changed_samples",
                                                       "qualified_candidate_count", "wall_seconds")}))
        return 0
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
