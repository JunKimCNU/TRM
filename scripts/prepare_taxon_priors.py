#!/usr/bin/env python3
"""Build auditable sample priors from archive taxonomy and curated references.

Species/isolate exceptions precede curated genus and nearest-clade defaults.
Unknown or complex references prevent broader inheritance without claiming
telomere absence. Taxonomy is never inferred from predicted repeat content.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import time

from reference_evidence import SUPPORTED_REFERENCE_EVIDENCE


DIRECT_EVIDENCE = SUPPORTED_REFERENCE_EVIDENCE
REFERENCE_STATUSES = {"supported", "candidate", "unknown", "out_of_scope"}
SCOPE_PRIORITY = {"clade": 1, "genus": 2, "species": 3, "isolate": 4}
REF_FIELDS = {"reference_id", "scientific_name", "reference_taxid", "scope",
              "related_motifs", "evidence", "source", "notes"}
META_FIELDS = {"sample_id", "status", "scientific_name", "tax_id", "lineage"}
ASSIGNMENT_FIELDS = [
    "sample_id", "metadata_status", "scientific_name", "tax_id", "lineage",
    "mapping_status", "related_motifs", "reference_ids", "reference_evidence",
    "prior_source", "review_reference_ids", "review_motifs", "notes",
    "reference_scope", "reference_status", "mapping_basis", "repeat_architecture",
    "suppressed_reference_ids",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path: Path, required: set[str]) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if len(set(fields)) != len(fields) or not required.issubset(fields):
            raise ValueError(f"{path}: missing required or duplicate columns")
        rows = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{reader.line_num}: incorrect field count")
            rows.append({key: value.strip() for key, value in row.items()})
    return rows


def motif_list(value: str) -> list[str]:
    motifs = value.upper().split(",")
    if any(not re.fullmatch(r"[ACGT]{1,64}", motif) for motif in motifs):
        raise ValueError(f"invalid related_motifs: {value!r}")
    return list(dict.fromkeys(motifs))


def load_references(path: Path) -> list[dict]:
    rows = table(path, REF_FIELDS)
    identifiers = set()
    for row in rows:
        rid = row["reference_id"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", rid) or rid in identifiers:
            raise ValueError(f"duplicate or invalid reference_id: {rid!r}")
        identifiers.add(rid)
        if row["scope"] not in SCOPE_PRIORITY:
            raise ValueError(f"{rid}: scope must be species, isolate, genus, or clade")
        if not re.fullmatch(r"[1-9][0-9]*", row["reference_taxid"]):
            raise ValueError(f"{rid}: require positive reference TaxID")
        if not row["scientific_name"] or not row["evidence"]:
            raise ValueError(f"{rid}: missing scientific name or evidence type")
        if not row["source"].startswith(("https://", "http://")):
            raise ValueError(f"{rid}: require a source URL")
        if row["scope"] in {"genus", "clade", "isolate"} and not row["notes"]:
            raise ValueError(f"{rid}: broader/isolate scope requires a curation note")
        status = row.get("reference_status") or (
            "supported" if row["evidence"] in DIRECT_EVIDENCE else "candidate")
        if status not in REFERENCE_STATUSES:
            raise ValueError(f"{rid}: invalid reference_status: {status!r}")
        if status == "supported" and row["evidence"] not in DIRECT_EVIDENCE:
            raise ValueError(f"{rid}: supported reference requires an eligible evidence type")
        row["reference_status"] = status
        row["repeat_architecture"] = row.get("repeat_architecture") or "exact_short_unit"
        aliases = [value.strip() for value in row.get("sample_ids", "").split(",") if value.strip()]
        if aliases and row["scope"] != "isolate":
            raise ValueError(f"{rid}: accession aliases require isolate scope")
        if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", value) for value in aliases):
            raise ValueError(f"{rid}: invalid sample_ids alias")
        row["sample_ids"] = aliases
        if not row["related_motifs"] and status in {"unknown", "out_of_scope"}:
            row["motifs"] = []
        else:
            row["motifs"] = motif_list(row["related_motifs"])
    return rows


def status_of(reference: dict) -> str:
    return reference.get("reference_status") or (
        "supported" if reference["evidence"] in DIRECT_EVIDENCE else "candidate")


def match_reference(sample: dict, reference: dict, names: list[str], ids: list[str]):
    """Return specificity and identity basis, using IDs ahead of name fallback."""
    scope, taxid = reference["scope"], reference["reference_taxid"]
    aliases = reference.get("sample_ids", [])
    if isinstance(aliases, str):
        aliases = [value.strip() for value in aliases.split(",") if value.strip()]
    if scope == "isolate" and (sample["sample_id"] in aliases
                                or sample.get("run_accession") in aliases):
        return (SCOPE_PRIORITY[scope], 0), "curated_accession_isolate_alias"
    if taxid == sample["tax_id"]:
        return (SCOPE_PRIORITY[scope], len(ids or names)), "exact_taxid"
    # An isolate name is not a species/genus-wide expectation.
    if scope == "isolate":
        return None
    if ids:
        if taxid in ids:
            return (SCOPE_PRIORITY[scope], ids.index(taxid)), "lineage_taxid"
        return None  # A conflicting/incomplete supplied ID list is not silently overridden.
    if scope in {"genus", "clade"} and reference["scientific_name"] in names:
        return (SCOPE_PRIORITY[scope], names.index(reference["scientific_name"])), "lineage_name_fallback"
    return None


def assign(metadata: list[dict], references: list[dict]) -> tuple[list[dict], list[dict]]:
    assignments, priors = [], []
    seen = set()
    for sample in metadata:
        sid = sample["sample_id"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", sid) or sid in seen:
            raise ValueError(f"invalid or duplicate metadata sample_id: {sid!r}")
        seen.add(sid)
        row = {key: "" for key in ASSIGNMENT_FIELDS}
        row.update({key: sample[key] for key in ("sample_id", "scientific_name", "tax_id", "lineage")})
        row["metadata_status"] = sample["status"]
        row["reference_status"] = "not_curated"
        row["mapping_status"] = "global_prior_metadata_unresolved"
        if sample["status"] != "resolved":
            row["notes"] = "Archive taxonomy unresolved; no species or motif inferred."
            assignments.append(row)
            continue
        if not re.fullmatch(r"[1-9][0-9]*", sample["tax_id"]) or not sample["scientific_name"]:
            raise ValueError(f"{sid}: resolved metadata requires valid tax_id and scientific_name")
        if (sample.get("metagenome", "false").lower() == "true"
                or sample.get("library_source", "").upper() == "METAGENOMIC"):
            row["mapping_status"] = "global_prior_mixed_sample"
            row["notes"] = "Archive marks metagenome; no single-organism prior applied."
            assignments.append(row)
            continue
        lineage = {name.strip() for name in sample["lineage"].split(";") if name.strip()}
        if lineage.intersection({"Bacteria", "Archaea"}):
            row["mapping_status"] = "global_prior_non_eukaryotic_metadata"
            row["notes"] = "Non-eukaryotic archive taxon: eukaryotic TRM interpretation needs review; global candidates retained."
            assignments.append(row)
            continue
        lineage_names = [item.strip() for item in sample["lineage"].split(";") if item.strip()]
        lineage_ids = [item.strip() for item in sample.get("lineage_tax_ids", "").split(";") if item.strip()]
        if any(not re.fullmatch(r"[1-9][0-9]*", item) for item in lineage_ids):
            raise ValueError(f"{sid}: malformed optional lineage_tax_ids")
        matches = []
        for reference in references:
            match = match_reference(sample, reference, lineage_names, lineage_ids)
            if match:
                matches.append((reference, *match))
        best_rank = max((rank for _, rank, _ in matches), default=None)
        specific = [(ref, basis) for ref, rank, basis in matches if rank == best_rank]
        # The nearest curated context determines the interpretation. Missing or
        # complex knowledge at that level must not inherit an ancestral motif.
        selected = [ref for ref, _ in specific if status_of(ref) == "supported"]
        blocking = [ref for ref, _ in specific if status_of(ref) in {"unknown", "out_of_scope"}]
        if blocking:
            selected = []
        selected_ids = {ref["reference_id"] for ref in selected}
        review = [ref for ref, _ in specific if ref["reference_id"] not in selected_ids]
        row["suppressed_reference_ids"] = ",".join(ref["reference_id"] for ref, rank, _ in matches
                                                   if rank != best_rank)
        row["review_reference_ids"] = ",".join(ref["reference_id"] for ref in review)
        row["review_motifs"] = ",".join(sorted({m for ref in review for m in ref["motifs"]}))
        row["mapping_status"] = "global_prior_no_curated_reference"
        row["notes"] = "No applicable curated reference; absence of curation is not telomere absence."
        if specific:
            row["reference_scope"] = ",".join(sorted({ref["scope"] for ref, _ in specific}))
            row["mapping_basis"] = ",".join(sorted({basis for _, basis in specific}))
            row["repeat_architecture"] = ",".join(sorted({ref.get("repeat_architecture", "exact_short_unit")
                                                          for ref, _ in specific}))
        if review and not selected:
            statuses = {status_of(ref) for ref in review}
            status = "out_of_scope" if "out_of_scope" in statuses else (
                "unknown" if "unknown" in statuses else "candidate")
            row["reference_status"] = status
            row["mapping_status"] = f"specific_reference_{status}"
            row["notes"] = ("More-specific curated context blocks broader default inheritance. "
                            "No short-motif prior emitted; raw discovery remains available. "
                            + " | ".join(ref["notes"] for ref in review))
        if selected:
            scope = selected[0]["scope"]
            row["mapping_status"] = {"species": "exact_species_reference", "isolate": "curated_isolate_reference",
                                     "genus": "curated_genus_reference", "clade": "clade_default_reference"}[scope]
            row["reference_status"] = "supported"
            row["related_motifs"] = ",".join(sorted({m for ref in selected for m in ref["motifs"]}))
            row["reference_ids"] = ",".join(ref["reference_id"] for ref in selected)
            row["reference_evidence"] = ",".join(sorted({ref["evidence"] for ref in selected}))
            row["prior_source"] = " | ".join(sorted({ref["source"] for ref in selected}))
            row["notes"] = ("All equally applicable supported motifs retained; reference use is not held-out validation. "
                            "Location remains unconfirmed in this sample. "
                            + ("Clade default is a weak expectation with known exceptions. " if scope == "clade" else "")
                            + " | ".join(ref["notes"] for ref in selected))
            if "lineage_name_fallback" in row["mapping_basis"]:
                row["notes"] += " Full lineage TaxIDs unavailable; matched exact archive lineage-name token."
            priors.append({
                "sample_id": sid,
                "lineage": f"{sample['scientific_name']} (TaxID {sample['tax_id']}); {row['mapping_status']}; {row['reference_ids']}",
                "related_motifs": row["related_motifs"], "source": row["prior_source"],
            })
        assignments.append(row)
    return assignments, priors


def write_table(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--references", type=Path, default=Path("data/taxon_motif_references.tsv"))
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    try:
        metadata = table(args.metadata, META_FIELDS)
        references = load_references(args.references)
        assignments, priors = assign(metadata, references)
        args.output_dir.mkdir(parents=True, exist_ok=False)
        write_table(args.output_dir / "assignments.tsv", assignments, ASSIGNMENT_FIELDS)
        write_table(args.output_dir / "priors.tsv", priors, ["sample_id", "lineage", "related_motifs", "source"])
        summary = {
            "schema_version": 2, "created_utc": datetime.now(timezone.utc).isoformat(),
            "metadata_path": str(args.metadata.resolve()), "metadata_sha256": sha256(args.metadata),
            "references_path": str(args.references.resolve()), "references_sha256": sha256(args.references),
            "script_sha256": sha256(Path(__file__)), "samples": len(metadata), "applied_priors": len(priors),
            "mapping_status_counts": dict(Counter(row["mapping_status"] for row in assignments)),
            "reference_count": len(references), "allowed_prior_evidence": sorted(DIRECT_EVIDENCE),
            "reference_status_counts": dict(Counter(row["reference_status"] for row in assignments)),
            "mapping_basis_counts": dict(Counter(row["mapping_basis"] or "none" for row in assignments)),
            "sequence_inputs_opened": False, "motif_results_used_for_taxonomy": False,
            "wall_seconds": time.perf_counter() - started,
            "output_sha256": {name: sha256(args.output_dir / name) for name in ("priors.tsv", "assignments.tsv")},
            "limitations": ["Archive taxonomy is submitter metadata, not a contamination check.",
                            "Taxonomic lineage is not a phylogeny with branch lengths.",
                            "Curated clade defaults are weak expectations; more-specific exceptions take precedence.",
                            "Exact lineage-name matching is used only when lineage TaxIDs are absent.",
                            "Unknown or complex references do not establish absence of telomeres or ALT.",
                            "Reference reuse and held-out validation independence are separate questions."]}
        with (args.output_dir / "summary.json").open("x", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        print(json.dumps({key: summary[key] for key in ("samples", "applied_priors", "mapping_status_counts", "wall_seconds")}))
        return 0
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
