#!/usr/bin/env python3
"""Resolve manifest run accessions through ENA, retaining verifiable raw evidence.

Only API metadata are read. FASTQ paths in the manifest are never opened.
Cached responses must pass both SHA-256 and request/response identity checks.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import re
import sys
import time
import urllib.parse
import urllib.request

FIELDS = (
    "run_accession", "sample_accession", "study_accession", "scientific_name",
    "tax_id", "library_strategy", "library_source", "library_layout",
    "instrument_platform", "instrument_model",
)
OUTPUT_FIELDS = (
    "sample_id", *FIELDS, "taxonomy_scientific_name", "taxonomy_rank", "lineage",
    "metagenome", "status", "error", "run_response_sha256", "taxonomy_response_sha256",
    "metadata_source", "taxonomy_source",
    "ena_scientific_name", "ena_tax_id", "taxonomy_assignment_source",
    "fallback_response_sha256",
    "requested_tax_id", "tax_id_resolution",
)
ACCESSION = re.compile(r"(?:SRR|ERR|DRR)[0-9]+\Z")
MAX_RESPONSE_BYTES = 2 * 1024 * 1024


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def json_bytes(value):
    return (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode()


def run_url(accession):
    return "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode({
        "accession": accession, "result": "read_run", "fields": ",".join(FIELDS),
        "format": "tsv",
    })


def taxonomy_url(tax_id):
    return f"https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{tax_id}"


def ncbi_run_url(accession):
    return "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({
        "db": "sra", "id": accession, "rettype": "runinfo", "retmode": "text",
    })


def load_accessions(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or "sample_id" not in reader.fieldnames:
            raise ValueError("manifest requires a sample_id column")
        accessions = []
        seen = set()
        for line, row in enumerate(reader, 2):
            accession = row.get("sample_id", "")
            if not ACCESSION.fullmatch(accession):
                raise ValueError(f"manifest line {line}: invalid run accession {accession!r}")
            if accession in seen:
                raise ValueError(f"manifest line {line}: duplicate run accession {accession}")
            seen.add(accession)
            accessions.append(accession)
        if not accessions:
            raise ValueError("manifest contains no run accessions")
        return accessions


def validate_run(raw, accession):
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
    if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise ValueError("ENA run response has absent or duplicate column names")
    if set(FIELDS) - set(reader.fieldnames):
        raise ValueError("ENA run response lacks requested fields")
    rows = list(reader)
    if len(rows) != 1 or rows[0].get("run_accession") != accession:
        raise ValueError(f"ENA run response must contain exactly requested run {accession}")
    row = rows[0]
    if None in row or any(row.get(key) is None for key in FIELDS):
        raise ValueError("ENA run response has malformed TSV columns")
    if row["tax_id"] and (not row["tax_id"].isdigit() or int(row["tax_id"]) <= 0):
        raise ValueError("ENA run response tax_id must be a positive integer")
    return {key: row[key] for key in FIELDS}


def validate_ncbi_run(raw, accession):
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    required = {"Run", "TaxID", "ScientificName", "BioSample"}
    if not reader.fieldnames or required - set(reader.fieldnames):
        raise ValueError("NCBI runinfo lacks requested identity/taxonomy fields")
    rows = list(reader)
    if len(rows) != 1 or rows[0].get("Run") != accession:
        raise ValueError(f"NCBI runinfo must contain exactly requested run {accession}")
    row = rows[0]
    if (not row.get("TaxID", "").isdigit() or int(row["TaxID"]) <= 0
            or not row.get("ScientificName", "").strip() or not row.get("BioSample", "").strip()):
        raise ValueError("NCBI runinfo lacks valid TaxID, ScientificName or BioSample")
    return row


def validate_taxonomy(raw, tax_id):
    value = json.loads(raw)
    if not isinstance(value, dict):
        raise ValueError("ENA taxonomy response must be an object")
    returned_id = str(value.get("taxId", ""))
    if not returned_id.isdigit() or int(returned_id) <= 0:
        raise ValueError("ENA taxonomy response taxId must be a positive integer")
    merged = value.get("merged", [])
    if isinstance(merged, str):
        try:
            merged = json.loads(merged)
        except (ValueError, TypeError):
            raise ValueError("ENA taxonomy merged IDs must be an explicit integer list") from None
    if (not isinstance(merged, list)
            or any(not str(identifier).isdigit() or int(identifier) <= 0 for identifier in merged)):
        raise ValueError("ENA taxonomy merged IDs must be an explicit integer list")
    if returned_id != str(tax_id) and str(tax_id) not in {str(identifier) for identifier in merged}:
        raise ValueError(f"ENA taxonomy response must match requested tax ID {tax_id}")
    for key in ("scientificName", "rank", "lineage"):
        if not isinstance(value.get(key), str) or not value[key].strip():
            raise ValueError(f"ENA taxonomy response lacks {key}")
    if value.get("metagenome") not in ("true", "false", True, False):
        raise ValueError("ENA taxonomy response lacks boolean metagenome flag")
    return value


def download(url, timeout):
    request = urllib.request.Request(url, headers={"User-Agent": "trm-metadata/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        raw = response.read(MAX_RESPONSE_BYTES + 1)
    if len(raw) > MAX_RESPONSE_BYTES:
        raise ValueError("API response exceeds 2 MiB safety limit")
    return raw


def write_new(path, data):
    """Exclusive creation; never replace or silently repair an existing cache."""
    try:
        with path.open("xb") as stream:
            stream.write(data)
    except FileExistsError:
        if path.read_bytes() != data:
            raise ValueError(f"refusing to replace existing different file: {path}")


class ResponseCache:
    def __init__(self, root, timeout=30, retries=2, transport=None, sleeper=None):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.retries = retries
        self.transport = transport or download
        self.sleeper = sleeper or time.sleep

    def fetch(self, kind, identifier, url, validator):
        if kind not in ("run", "taxonomy", "ncbi") or not re.fullmatch(r"[A-Za-z0-9]+", identifier):
            raise ValueError("invalid cache request key")
        stem = f"{kind}_{identifier}"
        raw_path = self.root / (stem + ".response")
        meta_path = self.root / (stem + ".provenance.json")
        if raw_path.exists() or meta_path.exists():
            if not raw_path.is_file() or not meta_path.is_file():
                raise ValueError(f"incomplete cache entry {stem}; preserve and repair explicitly")
            metadata = json.loads(meta_path.read_text())
            raw = raw_path.read_bytes()
            if (metadata.get("schema_version") != 1 or metadata.get("url") != url
                    or metadata.get("identifier") != identifier or metadata.get("kind") != kind
                    or metadata.get("sha256") != sha256(raw)
                    or metadata.get("bytes") != len(raw)
                    or not metadata.get("retrieved_utc")):
                raise ValueError(f"cache provenance or SHA-256 mismatch: {stem}")
            parsed = validator(raw, identifier)
            return parsed, metadata
        for attempt in range(self.retries + 1):
            try:
                raw = self.transport(url, self.timeout)
                parsed = validator(raw, identifier)
                break
            except Exception:
                if attempt == self.retries:
                    raise
                self.sleeper(min(2 ** attempt, 8))
        metadata = {"schema_version": 1, "kind": kind, "identifier": identifier,
                    "url": url, "retrieved_utc": utc_now(), "sha256": sha256(raw),
                    "bytes": len(raw), "response_file": str(raw_path.resolve())}
        write_new(raw_path, raw)
        write_new(meta_path, json_bytes(metadata))
        return parsed, metadata


def collect(keys, worker, workers, label):
    resolved, errors = {}, {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(worker, key): key for key in keys}
        for index, future in enumerate(as_completed(futures), 1):
            key = futures[future]
            try:
                resolved[key] = future.result()
            except Exception as error:
                errors[key] = f"{type(error).__name__}: {error}"
            if index % 10 == 0 or index == len(keys):
                print(f"{label}: {index}/{len(keys)} resolved={len(resolved)} unresolved={len(errors)}", flush=True)
    return resolved, errors


def write_result(path, data, update):
    if not path.exists() or path.read_bytes() == data:
        write_new(path, data)
        return
    if not update:
        raise ValueError(f"different result already exists: {path}; use --update to preserve a backup and replace it")
    old = path.read_bytes()
    history = path.parent / "history"
    history.mkdir(exist_ok=True)
    backup = history / f"{path.name}.{sha256(old)}"
    write_new(backup, old)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    write_new(temporary, data)
    os.replace(temporary, path)
    print(f"Updated {path}; previous content preserved at {backup}", flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=30)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--update", action="store_true", help="replace changed summaries while preserving content-addressed backups")
    args = parser.parse_args(argv)
    if not 1 <= args.workers <= 16 or args.timeout <= 0 or not 0 <= args.retries <= 5:
        parser.error("require workers 1..16, timeout > 0 and retries 0..5")
    accessions = load_accessions(args.manifest)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cache = ResponseCache(args.output_dir / "raw", args.timeout, args.retries)
    started = utc_now()
    started_monotonic = time.monotonic()
    runs, run_errors = collect(accessions, lambda key: cache.fetch("run", key, run_url(key), validate_run),
                               args.workers, "Run metadata")
    original_runs = {key: dict(value[0]) for key, value in runs.items()}
    missing_ids = [key for key, value in runs.items() if not value[0]["tax_id"]]
    # A single worker keeps unauthenticated NCBI fallback traffic below the
    # three-requests-per-second limit with the explicit minimum interval.
    def fetch_ncbi(key):
        response = cache.fetch("ncbi", key, ncbi_run_url(key), validate_ncbi_run)
        time.sleep(0.35)
        return response
    fallback, fallback_errors = collect(missing_ids, fetch_ncbi, 1, "NCBI linked-run fallback")
    for accession, (response, source) in fallback.items():
        run = runs[accession][0]
        if run["sample_accession"] != response["BioSample"]:
            fallback_errors[accession] = "NCBI BioSample differs from the exact ENA run's sample accession; no assignment made"
            continue
        run["tax_id"] = response["TaxID"]
        if not run["scientific_name"]:
            run["scientific_name"] = response["ScientificName"]
    tax_ids = sorted({row[0]["tax_id"] for row in runs.values() if row[0]["tax_id"]}, key=int)
    distribution = Counter(row[0]["scientific_name"] or "[name pending taxonomy]" for row in runs.values())
    print("Run species distribution: " + json.dumps(dict(distribution.most_common()), ensure_ascii=False), flush=True)
    taxonomy, taxonomy_errors = collect(tax_ids, lambda key: cache.fetch("taxonomy", key, taxonomy_url(key), validate_taxonomy),
                                       args.workers, "Taxonomy")
    rows = []
    unresolved = []
    response_provenance = []
    warnings = []
    for accession in accessions:
        row = {field: "" for field in OUTPUT_FIELDS}
        row["sample_id"] = accession
        if accession in run_errors:
            row.update(status="unresolved_run", error=run_errors[accession])
        else:
            run, run_source = runs[accession]
            row.update(run)
            row["ena_scientific_name"] = original_runs[accession]["scientific_name"]
            row["ena_tax_id"] = original_runs[accession]["tax_id"]
            row["run_response_sha256"] = run_source["sha256"]
            row["metadata_source"] = run_source["url"]
            row["taxonomy_assignment_source"] = run_source["url"]
            response_provenance.append(run_source)
            tax_id = run["tax_id"]
            row["requested_tax_id"] = tax_id
            if accession in fallback and accession not in fallback_errors:
                row["fallback_response_sha256"] = fallback[accession][1]["sha256"]
                row["taxonomy_assignment_source"] = fallback[accession][1]["url"]
            if accession in fallback_errors:
                row.update(status="unresolved_taxonomy_assignment", error=fallback_errors[accession])
            elif tax_id in taxonomy_errors:
                row.update(status="unresolved_taxonomy", error=taxonomy_errors[tax_id])
            else:
                taxon, tax_source = taxonomy[tax_id]
                row.update(taxonomy_scientific_name=taxon["scientificName"], taxonomy_rank=taxon["rank"],
                           lineage=taxon["lineage"], metagenome=str(taxon["metagenome"]).lower(),
                           taxonomy_response_sha256=tax_source["sha256"], status="resolved")
                row["taxonomy_source"] = tax_source["url"]
                row["tax_id"] = str(taxon["taxId"])
                row["tax_id_resolution"] = "exact" if row["tax_id"] == tax_id else "explicit_merge"
                if row["tax_id_resolution"] == "explicit_merge":
                    row["taxonomy_assignment_source"] = tax_source["url"]
                if not row["scientific_name"]:
                    row["scientific_name"] = taxon["scientificName"]
                    row["taxonomy_assignment_source"] = tax_source["url"]
                if run["scientific_name"] and run["scientific_name"] != taxon["scientificName"]:
                    warnings.append({"sample_id": accession, "warning": "run and taxonomy scientific names differ; both retained",
                                     "run_name": run["scientific_name"], "taxonomy_name": taxon["scientificName"]})
        if row["status"] != "resolved":
            unresolved.append({"sample_id": accession, "status": row["status"], "error": row["error"]})
        rows.append(row)
    response_provenance.extend(value[1] for value in taxonomy.values())
    response_provenance.extend(value[1] for value in fallback.values())
    response_provenance.sort(key=lambda value: (value["kind"], value["identifier"]))
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    table_bytes = output.getvalue().encode("utf-8")
    provenance = {
        "schema_version": 1, "started_utc": started, "finished_utc": utc_now(),
        "elapsed_seconds": time.monotonic() - started_monotonic,
        "source": "European Nucleotide Archive run metadata and ENA taxonomy (NCBI taxonomy identifiers); NCBI SRA runinfo for explicitly linked missing TaxIDs",
        "manifest": str(args.manifest.resolve()), "manifest_sha256": sha256(args.manifest.read_bytes()),
        "script": str(Path(__file__).resolve()), "script_sha256": sha256(Path(__file__).read_bytes()),
        "sample_metadata_sha256": sha256(table_bytes), "requested_runs": len(accessions),
        "resolved_runs": len(rows) - len(unresolved), "unresolved": unresolved,
        "requested_taxonomy_ids": len(tax_ids),
        "unique_tax_ids": len({row["tax_id"] for row in rows if row["status"] == "resolved"}),
        "explicit_tax_id_merges": sum(row["tax_id_resolution"] == "explicit_merge" for row in rows),
        "species_distribution": dict(Counter(row["taxonomy_scientific_name"] for row in rows if row["status"] == "resolved").most_common()),
        "library_strategy_distribution": dict(Counter(row[0]["library_strategy"] for row in runs.values())),
        "warnings": warnings, "responses": response_provenance,
        "limitations": ["Deposited sample metadata are annotations, not an independently verified species identification.",
                        "Lineage strings are taxonomy classifications, not branch lengths or phylogenetic distance.",
                        "No species or lineage was inferred from read sequence or guessed for unresolved accessions.",
                        "TaxID changes are accepted only when the API explicitly lists the requested ID as merged into the returned current ID."],
        "fastq_files_opened": False,
    }
    write_result(args.output_dir / "sample_metadata.tsv", table_bytes, args.update)
    write_result(args.output_dir / "provenance.json", json_bytes(provenance), args.update)
    print(f"Saved {len(rows)} runs: {len(rows) - len(unresolved)} resolved, {len(unresolved)} unresolved; {len(tax_ids)} tax IDs", flush=True)
    return 1 if unresolved else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError) as error:
        print(f"metadata: {error}", file=sys.stderr)
        raise SystemExit(2)
