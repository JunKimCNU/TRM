#!/usr/bin/env python3
"""Archive ENA library protocols and flag explicitly reported WGA/MDA context.

No FASTQ is opened. Missing amplification text is not evidence of an unamplified
library. Raw responses, requested fields, hashes, and parsing decisions are kept.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import io
import json
from pathlib import Path
import re
import time
import urllib.parse

from fetch_sample_metadata import (
    ResponseCache, collect, json_bytes, load_accessions, sha256, utc_now, write_result,
)


FIELDS = (
    "run_accession", "experiment_accession", "sample_accession", "study_accession",
    "scientific_name", "library_strategy", "library_source", "library_layout",
    "library_selection", "library_name", "library_construction_protocol",
    "experimental_protocol", "extraction_protocol", "library_gen_protocol",
    "library_pcr_isolation_protocol", "pcr_isolation_protocol",
    "experiment_title", "sample_description",
)
TEXT_FIELDS = (
    "library_strategy", "library_selection", "library_name",
    "library_construction_protocol", "experimental_protocol", "extraction_protocol",
    "library_gen_protocol", "library_pcr_isolation_protocol", "pcr_isolation_protocol",
    "experiment_title", "sample_description",
)
OUTPUT_FIELDS = (
    "sample_id", *FIELDS, "status", "genomic_wgs_context", "wga_mda_status",
    "wga_mda_flag", "whole_genome_amplification", "wga_evidence",
    "library_context_notes", "source_url", "response_sha256", "error",
)
RETURN_FIELDS_URL = "https://www.ebi.ac.uk/ena/portal/api/returnFields?result=read_run"
WGA_PATTERN = re.compile(
    r"\b(?:WGA|MDA|MALBAC|PicoPLEX|GenomePlex|GenomiPhi|REPLI[\s-]*g)\b"
    r"|\bwhole[\s-]*genome[\s-]*amplif(?:ication|ied|ying)\b"
    r"|\bmultiple[\s-]*displacement[\s-]*amplif(?:ication|ied|ying)\b",
    re.IGNORECASE,
)


def run_url(accession: str) -> str:
    return "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode({
        "accession": accession, "result": "read_run", "fields": ",".join(FIELDS), "format": "tsv",
    })


def validate_fields(raw: bytes, _identifier: str) -> list[str]:
    rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t"))
    available = {row.get("columnId") for row in rows}
    missing = set(FIELDS) - available
    if missing:
        raise ValueError(f"ENA read_run returnFields lacks requested columns: {sorted(missing)}")
    return sorted(available)


def validate_run(raw: bytes, accession: str) -> dict:
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
    if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise ValueError("absent or duplicate ENA field names")
    if set(FIELDS) - set(reader.fieldnames):
        raise ValueError("ENA response lacks requested library context fields")
    rows = list(reader)
    if len(rows) != 1 or rows[0].get("run_accession") != accession:
        raise ValueError(f"ENA response must identify exactly requested run {accession}")
    if None in rows[0] or any(rows[0].get(field) is None for field in FIELDS):
        raise ValueError("malformed ENA library context TSV")
    return {field: rows[0][field] for field in FIELDS}


def mention_polarity(text: str, start: int, end: int) -> str:
    # This is a transparent text screen, not a biological determination. Keep
    # ambiguous/conflicting language visible for protocol review.
    prefix = re.split(r"[.;\n]", text[max(0, start - 100):start])[-1].lower()
    suffix = text[end:end + 75].lower()
    if (re.search(r"\b(?:no|without|not|never)\b(?:[\s-]+\w+){0,6}[\s:(-]*$", prefix)
            or re.match(r"\s*(?:\([^)]*\))?\s*(?:(?:was|were|is|are)\s+)?(?:not|never)\b", suffix)):
        return "negated"
    if re.search(r"\b(?:possible|possibly|may|might|whether|unknown)\b", prefix):
        return "uncertain"
    return "reported"


def interpret(row: dict) -> dict:
    evidence = []
    for field in TEXT_FIELDS:
        text = row.get(field, "")
        for match in WGA_PATTERN.finditer(text):
            evidence.append({
                "field": field, "term": match.group(),
                "polarity": mention_polarity(text, match.start(), match.end()),
                "context": text[max(0, match.start() - 100):match.end() + 130],
            })
    polarities = {item["polarity"] for item in evidence}
    status = "not_reported"
    if len(polarities) > 1 or "uncertain" in polarities:
        status = "uncertain"
    elif polarities:
        status = next(iter(polarities))
    strategy, source = row.get("library_strategy", "").upper(), row.get("library_source", "").upper()
    context = "other_or_missing"
    if source == "METAGENOMIC":
        context = "metagenomic"
    elif source == "GENOMIC":
        context = "genomic_wgs" if strategy == "WGS" else "genomic_non_wgs"
    notes = ["Archive text screen; absence of WGA/MDA terms is not proof of an unamplified library."]
    if status == "reported":
        notes.append("Reported whole-genome amplification can distort repeat abundance and coverage; retain protocol evidence.")
    elif status in {"negated", "uncertain"}:
        notes.append("Negated, uncertain, or conflicting amplification language requires manual protocol review.")
    if any(re.search(r"\bPCR[\s-]*free\b", row.get(field, ""), re.IGNORECASE) for field in TEXT_FIELDS):
        notes.append("PCR-free library preparation does not by itself rule out earlier whole-genome amplification.")
    return {
        "genomic_wgs_context": context, "wga_mda_status": status,
        "wga_mda_flag": "true" if status == "reported" else "",
        "whole_genome_amplification": "true" if status == "reported" else "",
        "wga_evidence": json.dumps(evidence, ensure_ascii=False, separators=(",", ":")),
        "library_context_notes": " ".join(notes),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=30)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--update", action="store_true")
    args = parser.parse_args(argv)
    if not 1 <= args.workers <= 16 or args.timeout <= 0 or not 0 <= args.retries <= 5:
        parser.error("require workers 1..16, timeout > 0, retries 0..5")
    try:
        started, monotonic = utc_now(), time.monotonic()
        accessions = load_accessions(args.manifest)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        cache = ResponseCache(args.output_dir / "raw", args.timeout, args.retries)
        _, field_provenance = cache.fetch("run", "READRUNFIELDS", RETURN_FIELDS_URL, validate_fields)
        runs, errors = collect(accessions, lambda key: cache.fetch("run", key, run_url(key), validate_run),
                               args.workers, "Library context")
        rows = []
        for accession in accessions:
            row = {field: "" for field in OUTPUT_FIELDS}
            row.update(sample_id=accession, run_accession=accession, source_url=run_url(accession))
            if accession in runs:
                data, provenance = runs[accession]
                row.update(data)
                row.update(interpret(data))
                row.update(status="resolved", response_sha256=provenance["sha256"])
            else:
                row.update(status="lookup_failed", error=errors[accession],
                           genomic_wgs_context="other_or_missing", wga_mda_status="metadata_missing",
                           wga_evidence="[]", library_context_notes="Metadata unavailable; amplification status unknown.")
            rows.append(row)
        stream = io.StringIO(newline="")
        writer = csv.DictWriter(stream, OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        output = stream.getvalue().encode()
        write_result(args.output_dir / "library_context.tsv", output, args.update)
        report = {
            "schema_version": 1, "started_utc": started, "finished_utc": utc_now(),
            "wall_seconds": time.monotonic() - monotonic, "samples": len(rows),
            "resolved": len(runs), "failed": len(errors), "requested_fields": FIELDS,
            "manifest": str(args.manifest.resolve()), "manifest_sha256": sha256(args.manifest.read_bytes()),
            "script_sha256": sha256(Path(__file__).read_bytes()),
            "cache_helper_sha256": sha256(Path(__file__).with_name("fetch_sample_metadata.py").read_bytes()),
            "output_sha256": sha256(output), "field_definition_provenance": field_provenance,
            "wga_mda_status_counts": dict(Counter(row["wga_mda_status"] for row in rows)),
            "genomic_wgs_context_counts": dict(Counter(row["genomic_wgs_context"] for row in rows)),
            "positive_wga_runs": [row["sample_id"] for row in rows if row["wga_mda_status"] == "reported"],
            "unresolved_runs": errors, "sequence_inputs_opened": False,
            "limitations": ["Archive protocol text may be incomplete or inconsistent.",
                            "A not_reported result does not prove absence of whole-genome amplification.",
                            "PCR-free library preparation does not rule out upstream WGA.",
                            "Negation and ambiguous text are heuristically flagged for review, not resolved automatically."],
        }
        write_result(args.output_dir / "provenance.json", json_bytes(report), args.update)
        print(json.dumps({key: report[key] for key in ("samples", "resolved", "failed", "wall_seconds",
                                                     "wga_mda_status_counts", "genomic_wgs_context_counts")}))
        return 0 if not errors else 1
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
