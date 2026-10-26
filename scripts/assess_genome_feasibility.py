#!/usr/bin/env python3
"""Plan k-mer genome-spectrum pilots from saved read QC and public metadata.

This command neither opens FASTQ streams nor estimates a genome size.  It keeps
data completeness, protocol conflicts, library caveats, and biological scope as
independent flags.  Planning eligibility is not evidence of adequate coverage.
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


FIELDS = (
    "sample_id", "scientific_name", "tax_id", "genome_estimation_status",
    "input_status", "analysis_complete", "local_input_file_count",
    "local_existing_file_count", "local_input_compressed_bytes", "local_paths",
    "archive_library_layout", "local_layout_interpretation", "local_pairing_verified",
    "records", "bases", "min_read_length_lower_bound", "max_read_length",
    "mean_read_length", "records_with_ambiguity", "ambiguous_bases",
    "kmer_length", "eligible_records_lower_bound", "eligible_records_upper_bound",
    "sampled_prefix", "instrument_platform", "instrument_model",
    "library_strategy", "library_source", "library_selection",
    "library_construction_protocol", "library_protocol_flags", "all_flags",
    "review_reasons", "genome_size_status", "genome_size_bp", "karyotype_status",
    "haploid_chromosome_count", "whole_genome_amplification_status",
    "archive_wga_mda_status", "library_pcr_status", "coverage_status",
    "genome_estimation_performed", "read_quality_assessed", "duplicates_assessed",
    "metadata_provider", "library_source_url", "library_response_sha256",
    "assessment_json", "assessment_sha256",
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def indexed_tsv(path: Path) -> dict[str, dict]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "sample_id" not in reader.fieldnames:
            raise ValueError(f"sample_id column missing: {path}")
        rows = list(reader)
    result = {}
    for row in rows:
        sid = row["sample_id"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sid) or sid in result:
            raise ValueError(f"Unsafe or duplicate sample_id {sid!r}: {path}")
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"Malformed TSV row {sid}: {path}")
        result[sid] = row
    return result


def interpret_sample(manifest: dict, metadata: dict, library: dict, assessed: dict,
                     audit_flags: list[dict], k: int, *, inspect_paths=True) -> dict:
    sid = manifest["sample_id"]
    if assessed.get("sample_id") != sid:
        raise ValueError(f"Assessment sample identity mismatch: {sid}")
    for key in ("library_strategy", "library_source", "library_layout"):
        if metadata.get(key, "").strip() != library.get(key, "").strip():
            raise ValueError(f"Metadata tables disagree for {sid}: {key}")
    flags = {flag["code"]: dict(flag) for flag in audit_flags}

    def flag(code, reason):
        flags.setdefault(code, {"code": code, "reason": reason})

    paths = [Path(manifest[key]) for key in ("path1", "path2") if manifest.get(key)]
    if len(set(paths)) != len(paths):
        flag("duplicate_local_path", "The same local input path is listed more than once.")
    exists = [path.is_file() for path in paths] if inspect_paths else [True] * len(paths)
    complete = assessed.get("analysis_complete") is True
    qc = assessed.get("data_adequacy", {}).get("read_qc") or {}
    summary = assessed.get("input_summary") or {}
    if not complete:
        flag("prior_input_analysis_failed", "The saved complete-input TRM analysis failed; valid whole-input QC is unavailable.")
    if not paths or not all(exists):
        flag("missing_local_input", "One or more listed local sequence files are unavailable.")
    records, bases = qc.get("records"), qc.get("bases")
    qc_available = isinstance(records, int) and records >= 0 and isinstance(bases, int) and bases >= 0
    if complete and not qc_available:
        flag("missing_read_qc", "The completed assessment has no usable read/base totals.")
    if qc_available:
        for key in ("records", "bases"):
            if key in summary and summary[key] != qc[key]:
                raise ValueError(f"Saved summary/QC disagree for {sid}: {key}")
        if records == 0 or bases == 0:
            flag("empty_input", "The saved QC contains no sequence records or bases.")
    histogram = qc.get("length_histogram") or []
    lower = upper = None
    minimum = None
    if histogram:
        if not qc_available or sum(item["records"] for item in histogram) != records:
            raise ValueError(f"Length histogram count mismatch: {sid}")
        if any(not (0 <= item["lower"] <= item["upper"] and item["records"] >= 0)
               for item in histogram):
            raise ValueError(f"Invalid length histogram: {sid}")
        lower = sum(item["records"] for item in histogram if item["lower"] >= k)
        upper = sum(item["records"] for item in histogram if item["upper"] >= k)
        nonempty = [item for item in histogram if item["records"]]
        minimum = min((item["lower"] for item in nonempty), default=None)
        if upper == 0:
            flag("reads_shorter_than_k", "All recorded reads are shorter than the requested k-mer length.")
    sampled = qc.get("sampled_prefix", summary.get("sampled_prefix"))
    if sampled is True:
        flag("prefix_qc_only", "Saved QC covers a prefix; whole-input totals are not available.")
    platform = metadata.get("instrument_platform", "").upper()
    if platform not in {"ILLUMINA", "BGISEQ", "DNBSEQ"}:
        flag("nonstandard_platform", "The platform needs an explicit error-model and spectrum-fitting compatibility review.")
    source = library.get("library_source", "").upper()
    if source != "GENOMIC" or library.get("library_strategy", "").upper() != "WGS":
        flag("not_standard_genomic_wgs", "Deposited strategy/source does not describe genomic WGS.")
    if not metadata.get("lineage", "").strip().startswith("Eukaryota"):
        flag("non_eukaryote_taxon", "Outside this eukaryotic TRM panel; genome-size fitting requires an appropriate organism model.")
    wga = library.get("wga_mda_status", "not_reported")
    if wga == "reported":
        flag("reported_whole_genome_amplification", "Archive text reports WGA; nonuniform amplification requires model review.")
    elif wga in {"negated", "uncertain"}:
        flag("amplification_text_requires_review", "Negated or uncertain WGA language requires source review; absence is not established.")
    layout = library.get("library_layout", "")
    layout_state = "input_layout_not_verified"
    if layout == "PAIRED" and len(paths) == 2:
        layout_state = "two_files_consistent_with_archive_paired_label_not_pair_validated"
    elif layout == "PAIRED" and len(paths) == 1:
        layout_state = "one_file_may_be_single_mate_or_interleaved_unverified"
        flag("archive_paired_local_one_file", "Archive reports paired reads but one local file is listed; it may contain one mate or interleaved reads.")
    elif layout == "SINGLE" and len(paths) == 1:
        layout_state = "one_file_consistent_with_archive_single_label"
    elif layout == "SINGLE" and len(paths) > 1:
        flag("archive_single_multiple_local_files", "Multiple local files are listed for an archive single-read library.")

    codes = set(flags)
    status = "candidate_for_spectrum_pilot"
    if codes & {"prior_input_analysis_failed", "missing_local_input", "empty_input", "duplicate_local_path"}:
        status = "failed_input"
    elif codes & {"missing_read_qc", "prefix_qc_only", "reads_shorter_than_k"}:
        status = "review_input_qc"
    elif codes & {"hic_protocol_conflicts_with_wgs_label", "rna_protocol_conflicts_with_wgs_label"}:
        status = "review_protocol_conflict"
    elif "not_standard_genomic_wgs" in codes:
        status = "review_non_genomic_wgs"
    elif "non_eukaryote_taxon" in codes:
        status = "outside_eukaryotic_trm_scope"
    elif "nonstandard_platform" in codes:
        status = "nonstandard_platform"
    elif codes & {"reported_whole_genome_amplification", "amplification_text_requires_review"}:
        status = "review_amplification"
    elif codes & {"mate_pair_protocol", "linked_read_protocol"}:
        status = "review_library_structure"
    elif "pooled_population_title" in codes:
        status = "review_population_model"
    protocol_codes = {f["code"] for f in audit_flags}
    return {
        "sample_id": sid, "scientific_name": metadata.get("scientific_name", ""),
        "tax_id": metadata.get("tax_id", ""), "genome_estimation_status": status,
        "input_status": "failed" if status == "failed_input" else ("complete" if complete else "unknown"),
        "analysis_complete": complete, "local_input_file_count": len(paths),
        "local_existing_file_count": sum(exists),
        "local_input_compressed_bytes": sum(p.stat().st_size for p, ok in zip(paths, exists) if ok) if inspect_paths else None,
        "local_paths": json.dumps([str(p) for p in paths]), "archive_library_layout": layout,
        "local_layout_interpretation": layout_state, "local_pairing_verified": False,
        "records": records if qc_available else None, "bases": bases if qc_available else None,
        "min_read_length_lower_bound": minimum, "max_read_length": qc.get("max_record_bases"),
        "mean_read_length": bases / records if qc_available and records else None,
        "records_with_ambiguity": qc.get("records_with_ambiguity"), "ambiguous_bases": qc.get("ambiguous_bases"),
        "kmer_length": k, "eligible_records_lower_bound": lower, "eligible_records_upper_bound": upper,
        "sampled_prefix": sampled, "instrument_platform": platform,
        "instrument_model": metadata.get("instrument_model", ""),
        "library_strategy": library.get("library_strategy", ""), "library_source": source,
        "library_selection": library.get("library_selection", ""),
        "library_construction_protocol": library.get("library_construction_protocol", ""),
        "library_protocol_flags": ";".join(sorted(protocol_codes)), "all_flags": ";".join(sorted(codes)),
        "review_reasons": json.dumps(list(flags.values()), ensure_ascii=False, separators=(",", ":")),
        "genome_size_status": "unknown_not_estimated", "genome_size_bp": None,
        "karyotype_status": "unknown_not_inferred", "haploid_chromosome_count": None,
        "whole_genome_amplification_status": "reported" if wga == "reported" else "unknown",
        "archive_wga_mda_status": wga,
        "library_pcr_status": "archive_no_pcr_not_evidence_of_no_upstream_wga" if library.get("library_construction_protocol") == "No PCR" else "not_resolved",
        "coverage_status": "unknown_genome_size_and_spectrum_not_fitted", "genome_estimation_performed": False,
        "read_quality_assessed": qc.get("base_quality_distribution_available", False),
        "duplicates_assessed": qc.get("duplicate_fraction_measured", False),
        "metadata_provider": "ENA_INSDC_with_separately_sourced_NCBI_spot_checks",
        "library_source_url": library.get("source_url", ""),
        "library_response_sha256": library.get("response_sha256", ""),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("manifests/existing_fastq.tsv"))
    parser.add_argument("--assessed-dir", type=Path, default=Path("results/wgs_assessed_v040"))
    parser.add_argument("--metadata", type=Path, default=Path("data/sample_metadata_20260920/sample_metadata.tsv"))
    parser.add_argument("--library-context", type=Path, default=Path("data/library_context_20260921/library_context.tsv"))
    parser.add_argument("--library-audit", type=Path, default=Path("results/genome_feasibility_20260921/library_metadata_audit.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/genome_feasibility_20260921"))
    parser.add_argument("--k", type=int, default=21)
    args = parser.parse_args(argv)
    if not 3 <= args.k <= 255:
        parser.error("k must be between 3 and 255")
    try:
        manifest = indexed_tsv(args.manifest)
        metadata = indexed_tsv(args.metadata)
        library = indexed_tsv(args.library_context)
        if set(manifest) != set(metadata) or set(manifest) != set(library):
            raise ValueError("Manifest, taxonomy metadata and library context must identify the same samples")
        audit = json.loads(args.library_audit.read_text())
        if audit["record_count"] != len(manifest):
            raise ValueError("Library audit record count does not match manifest")
        expected_hashes = {item["sha256"] for item in audit["inputs"]}
        if digest(args.metadata) not in expected_hashes or digest(args.library_context) not in expected_hashes:
            raise ValueError("Library audit is not bound to these metadata files")
        flags = {}
        for entry in audit["flagged_samples"]:
            sid = entry["sample_id"]
            if sid in flags or sid not in manifest:
                raise ValueError(f"Duplicate or unknown audit sample: {sid}")
            flags[sid] = entry["flags"]
        output = args.output_dir.resolve()
        source = args.assessed_dir.resolve()
        if output == source or source in output.parents or output in source.parents:
            raise ValueError("Output must be separate from the existing assessed-results directory")
        rows, artifacts = [], []
        for sid, item in manifest.items():
            path = args.assessed_dir / "samples" / f"{sid}.json"
            assessed = json.loads(path.read_text())
            row = interpret_sample(item, metadata[sid], library[sid], assessed, flags.get(sid, []), args.k)
            row.update(assessment_json=str(path.resolve()), assessment_sha256=digest(path))
            rows.append(row)
            artifacts.append({"sample_id": sid, "path": str(path.resolve()), "sha256": row["assessment_sha256"]})
        args.output_dir.mkdir(parents=True, exist_ok=True)
        tsv = args.output_dir / "samples.tsv"
        with tsv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        summary = {
            "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
            "sample_count": len(rows), "genome_estimation_performed": False, "sequence_inputs_opened": False,
            "existing_assessment_results_modified": False,
            "status_counts": dict(Counter(r["genome_estimation_status"] for r in rows)),
            "input_status_counts": dict(Counter(r["input_status"] for r in rows)),
            "local_input_file_count_distribution": dict(Counter(r["local_input_file_count"] for r in rows)),
            "local_layout_counts": dict(Counter(r["local_layout_interpretation"] for r in rows)),
            "all_flag_counts": dict(Counter(code for r in rows for code in r["all_flags"].split(";") if code)),
            "total_records_with_complete_qc": sum(r["records"] or 0 for r in rows if r["analysis_complete"]),
            "total_bases_with_complete_qc": sum(r["bases"] or 0 for r in rows if r["analysis_complete"]),
            "unknown_genome_size_count": len(rows), "unknown_karyotype_count": len(rows),
            "unknown_wga_count": sum(r["whole_genome_amplification_status"] == "unknown" for r in rows),
            "kmer_length_for_read_length_eligibility_only": args.k,
            "limitations": [
                "candidate_for_spectrum_pilot is a planning status, not a successful fit or proof of adequate coverage.",
                "Eligible read counts only describe read length; they do not measure usable distinct k-mers, accuracy or genome coverage.",
                "Local file counts do not establish mate pairing, completeness of archive downloads, or read independence.",
                "No PCR and missing WGA annotations do not establish absence of whole-genome amplification.",
                "Unknown karyotype and telomere retention remain unresolved even if a later genome-size fit succeeds.",
                "Independent flags remain present when a higher-priority planning status is assigned.",
            ],
            "input_files": [{"path": str(p.resolve()), "sha256": digest(p)} for p in
                            (args.manifest, args.metadata, args.library_context, args.library_audit, Path(__file__))],
            "assessment_artifacts": artifacts,
            "samples_tsv": str(tsv.resolve()), "samples_tsv_sha256": digest(tsv),
        }
        (args.output_dir / "metadata_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
        print(json.dumps({key: summary[key] for key in ("sample_count", "status_counts", "input_status_counts")}, ensure_ascii=False))
        return 0
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    raise SystemExit(main())
