#!/usr/bin/env python3
"""Bounded, metadata-only ENA inventory for public eukaryotic WGS planning.

No sequence files, AWS resources, or paid queries are accessed.  The default is
a dry plan; --fetch collects archive-wide counts and at most 1000 run records.
The ENA search endpoint is partitioned by disjoint taxonomic queries; its first
N results are NOT a random sample and must not be extrapolated to the archive.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import time
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/ena/portal/api/"
BROAD_QUERY = 'tax_tree(2759) AND library_strategy="WGS"'
PRIMARY_QUERY = BROAD_QUERY + ' AND library_source="GENOMIC"'
PARTITIONS = (
    ("human", "tax_tree(9606)"),
    ("other_metazoa", "tax_tree(33208) AND NOT tax_tree(9606)"),
    ("viridiplantae", "tax_tree(33090)"),
    ("fungi", "tax_tree(4751)"),
    ("other_eukaryota", "NOT tax_tree(33208) AND NOT tax_tree(33090) AND NOT tax_tree(4751)"),
)
FIELDS = (
    "run_accession", "experiment_accession", "sample_accession", "study_accession",
    "tax_id", "scientific_name", "library_strategy", "library_source",
    "library_selection", "library_layout", "instrument_platform", "instrument_model",
    "read_count", "base_count", "fastq_bytes", "fastq_ftp", "fastq_md5",
    "sra_bytes", "sra_ftp", "submitted_bytes", "submitted_format", "submitted_ftp",
    "first_public", "last_updated",
)
DERIVED = ("inventory_partition", "fastq_total_bytes", "fastq_file_count", "sra_total_bytes",
           "submitted_total_bytes", "taxonomy_rank", "taxonomy_lookup_status")
SOURCES = {
    "ena_search": "https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/advanced-search.html",
    "ena_taxonomy_search": "https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/taxon-based-search.html",
    "ena_api": "https://www.ebi.ac.uk/ena/portal/api/swagger-ui/index.html",
    "aws_registry": "https://registry.opendata.aws/ncbi-sra/",
    "ncbi_aws_access": "https://www.ncbi.nlm.nih.gov/sra/docs/sra-aws-download/",
    "ncbi_metadata_units": "https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-metadata-table/",
}
RUN = re.compile(r"[SED]RR[0-9]+\Z")


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def api_url(endpoint, **parameters):
    return BASE + endpoint + "?" + urllib.parse.urlencode(parameters)


def query_for(partition):
    return PRIMARY_QUERY + " AND (" + dict(PARTITIONS)[partition] + ")"


def json_write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n")


def byte_sum(value):
    """Unknown is None, not zero; reject partial or malformed per-file lists."""
    if value == "":
        return None
    parts = value.split(";")
    if any(not part.isdigit() for part in parts):
        raise ValueError("invalid semicolon-delimited byte sizes")
    return sum(int(part) for part in parts)


def parse_count(raw):
    value = json.loads(raw)
    if not isinstance(value, dict) or not str(value.get("count", "")).isdigit():
        raise ValueError("invalid ENA count response")
    return int(value["count"])


def parse_group_counts(raw, field):
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
    if reader.fieldnames != [field, "count"]:
        raise ValueError("unexpected grouped count fields")
    result = {}
    for row in reader:
        key = row[field]
        if key in result or not row.get("count", "").isdigit() or None in row:
            raise ValueError("invalid or duplicate grouped count")
        result[key] = int(row["count"])
    return result


def parse_runs(raw, limit, partition):
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
    if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames) or set(FIELDS) - set(reader.fieldnames):
        raise ValueError("run response lacks unique requested fields")
    result, seen = [], set()
    for row in reader:
        if None in row or any(row.get(field) is None for field in FIELDS):
            raise ValueError("malformed TSV run row")
        accession = row["run_accession"]
        if not RUN.fullmatch(accession) or accession in seen:
            raise ValueError("invalid or duplicate run accession")
        if row["library_strategy"] != "WGS" or row["library_source"] != "GENOMIC":
            raise ValueError("run does not match requested WGS/GENOMIC query")
        if not row["tax_id"].isdigit() or int(row["tax_id"]) <= 0:
            raise ValueError("invalid tax_id")
        for field in ("read_count", "base_count"):
            if row[field] and not row[field].isdigit():
                raise ValueError("invalid count field")
        row = {field: row[field] for field in FIELDS}
        row.update(inventory_partition=partition, taxonomy_rank="", taxonomy_lookup_status="not_queried")
        for prefix in ("fastq", "sra", "submitted"):
            total = byte_sum(row[prefix + "_bytes"])
            row[prefix + "_total_bytes"] = "" if total is None else str(total)
        fastq_files = row["fastq_ftp"].split(";") if row["fastq_ftp"] else []
        if row["fastq_bytes"] and len(row["fastq_bytes"].split(";")) != len(fastq_files):
            raise ValueError("FASTQ byte/file list length mismatch")
        if row["fastq_md5"] and len(row["fastq_md5"].split(";")) != len(fastq_files):
            raise ValueError("FASTQ checksum/file list length mismatch")
        row["fastq_file_count"] = str(len(fastq_files))
        seen.add(accession)
        result.append(row)
        if len(result) > limit:
            raise ValueError("ENA returned more rows than bounded request")
    return result


class Fetcher:
    def __init__(self, root, timeout=60, max_response_bytes=2 * 1024 * 1024,
                 max_total_bytes=8 * 1024 * 1024):
        self.root = root
        self.timeout = timeout
        self.max_response_bytes = max_response_bytes
        self.max_total_bytes = max_total_bytes
        self.total = 0
        self.provenance = []

    def fetch(self, name, url):
        start = time.monotonic()
        entry = {"name": name, "url": url, "requested_at_utc": utc_now(), "status": "pending"}
        self.provenance.append(entry)
        cap = min(self.max_response_bytes, self.max_total_bytes - self.total)
        if cap <= 0:
            raise ValueError("metadata byte budget exhausted")
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "trm-public-inventory/1.0"})
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                raw = response.read(cap + 1)
                entry["response_url"] = response.url
                entry["content_type"] = response.headers.get("Content-Type", "")
            self.total += len(raw)
            if len(raw) > cap:
                raise ValueError("metadata response exceeds bounded byte budget")
            path = self.root / "raw" / (name + ".response")
            path.write_bytes(raw)
            entry.update(status="ok", bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest(),
                         local_path=str(path.relative_to(self.root)))
            return raw
        except Exception as error:
            entry.update(status="failed", error=f"{type(error).__name__}: {error}")
            raise
        finally:
            entry["elapsed_seconds"] = round(time.monotonic() - start, 6)
            json_write(self.root / "provenance.json", self.provenance)


def write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows, counts, partition_counts):
    """All size and taxon summaries describe only returned run records."""
    sampled = {"run_count": len(rows), "distinct_tax_ids": len({r["tax_id"] for r in rows}),
               "distinct_samples": len({r["sample_accession"] for r in rows if r["sample_accession"]}),
               "distinct_experiments": len({r["experiment_accession"] for r in rows}),
               "distinct_studies": len({r["study_accession"] for r in rows}),
               "rank_species_tax_ids": len({r["tax_id"] for r in rows if r["taxonomy_rank"] == "species"}),
               "by_platform": dict(Counter(r["instrument_platform"] for r in rows)),
               "by_library_source": dict(Counter(r["library_source"] for r in rows)),
               "by_partition": dict(Counter(r["inventory_partition"] for r in rows)),
               "total_bases_reported": sum(int(r["base_count"] or 0) for r in rows),
               "base_count_missing_runs": sum(not r["base_count"] for r in rows),
               "sampling_method": "first N archive-ordered results in each disjoint taxonomic query; not random",
               "archive_volume_extrapolation": None,
               "all_eukaryote_species_count": None}
    for prefix in ("fastq", "sra", "submitted"):
        field = prefix + "_total_bytes"
        sampled[prefix + "_bytes_sum_known_rows"] = sum(int(r[field]) for r in rows if r[field] != "")
        sampled[prefix + "_bytes_missing_runs"] = sum(r[field] == "" for r in rows)
    return {"archive_counts": counts, "partition_archive_counts": partition_counts,
            "pilot_only": sampled,
            "full_archive_species_count": None, "full_archive_data_bytes": None,
            "coverage": "Public ENA-indexed WGS records assigned within Eukaryota. GENOMIC subset is metadata-defined, not verified unbiased WGS.",
            "count_snapshot": "Independent live API requests, not an atomic database snapshot.",
            "unit_notes": {
                "run": "Archive run accession, not an independent individual, species, library, or molecule.",
                "tax_id": "Assigned archive taxon, which can be species, strain, subspecies, or an unresolved taxon.",
                "species": "Only tax_ids explicitly returned with taxonomy rank species; no complete species census or roll-up inferred.",
                "fastq_bytes": "Sum of per-file ENA FASTQ object sizes (bytes), normally compressed FASTQ; not bases or AWS SRA object bytes.",
                "sra_bytes": "ENA-reported per-file SRA byte sizes; not an independent HEAD measurement of AWS S3 objects.",
                "submitted_bytes": "Sizes of original submitted files; alternative representation, never added to FASTQ/SRA sizes.",
                "ncbi_mbytes": "NCBI cloud metadata defines megabytes of data in an SRA run; not FASTQ bytes, not a count of sequenced bases.",
                "read_count": "Archive-reported read_count retained without conversion to independent molecules or assumptions about paired spots.",
            }}


def run(args):
    root = Path(args.output)
    if root.exists() and any(root.iterdir()):
        raise ValueError("output directory is not empty; use a new path to preserve provenance")
    root.mkdir(parents=True, exist_ok=True)
    (root / "raw").mkdir()
    started = utc_now()
    plan = {"broad_query": BROAD_QUERY, "primary_query": PRIMARY_QUERY,
            "partitions": {name: query_for(name) for name, _ in PARTITIONS},
            "maximum_run_records": 5 * args.per_partition, "limit_per_partition": args.per_partition,
            "maximum_taxonomy_requests": args.taxonomy_limit,
            "network_sequence_downloads": False, "paid_cloud_queries": False,
            "sources": SOURCES, "started_at_utc": started}
    json_write(root / "plan.json", plan)
    if not args.fetch:
        print(json.dumps({"status": "dry_plan", "plan": str(root / "plan.json")}))
        return
    fetcher = Fetcher(root, timeout=args.timeout)
    counts, partition_counts, rows, grouped, errors = {}, {}, [], {}, []
    start = time.monotonic()

    def optional(name, operation):
        try:
            return operation()
        except Exception as error:
            errors.append({"operation": name, "error": f"{type(error).__name__}: {error}"})
            return None

    for name, query in (("wgs_all_sources", BROAD_QUERY), ("wgs_genomic", PRIMARY_QUERY)):
        counts[name] = optional(name, lambda name=name, query=query: parse_count(fetcher.fetch(
            "count_" + name, api_url("count", result="read_run", query=query, format="json"))))
    for name, query, field in (("all_wgs_by_source", BROAD_QUERY, "library_source"),
                               ("genomic_wgs_by_platform", PRIMARY_QUERY, "instrument_platform")):
        grouped[name] = optional(name, lambda name=name, query=query, field=field: parse_group_counts(fetcher.fetch(
            name, api_url("count", result="read_run", query=query, field=field, format="tsv")), field))
    seen = set()
    for name, _ in PARTITIONS:
        query = query_for(name)
        partition_counts[name] = optional("count_" + name, lambda name=name, query=query: parse_count(fetcher.fetch(
            "count_" + name, api_url("count", result="read_run", query=query, format="json"))))
        batch = optional("runs_" + name, lambda name=name, query=query: parse_runs(fetcher.fetch(
            "runs_" + name, api_url("search", result="read_run", query=query, fields=",".join(FIELDS),
                                    format="tsv", limit=args.per_partition)), args.per_partition, name))
        for row in batch or []:
            if row["run_accession"] in seen:
                raise ValueError("taxonomic partitions returned duplicate run; refusing double-counting")
            seen.add(row["run_accession"])
            rows.append(row)
    taxa = {}
    for tax_id in sorted({row["tax_id"] for row in rows}, key=int)[:args.taxonomy_limit]:
        def get_taxon(tax_id=tax_id):
            value = json.loads(fetcher.fetch("taxon_" + tax_id,
                "https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/" + tax_id))
            if not isinstance(value, dict) or str(value.get("taxId")) != tax_id or not isinstance(value.get("rank"), str):
                raise ValueError("taxonomy response identity/rank mismatch")
            return value
        taxa[tax_id] = optional("taxon_" + tax_id, get_taxon)
    for row in rows:
        taxon = taxa.get(row["tax_id"])
        if taxon:
            row["taxonomy_rank"] = taxon["rank"]
            row["taxonomy_lookup_status"] = "verified_current_id"
        elif row["tax_id"] in taxa:
            row["taxonomy_lookup_status"] = "lookup_failed"
    summary = summarize(rows, counts, partition_counts)
    summary.update(started_at_utc=started, finished_at_utc=utc_now(), elapsed_seconds=time.monotonic() - start,
                   metadata_response_bytes=fetcher.total, maximum_run_records=5 * args.per_partition,
                   full_query_group_counts=grouped, errors=errors, sources=SOURCES,
                   sequence_files_downloaded=0, status="complete_bounded_pilot" if not errors else "partial_pilot_with_errors")
    summary["count_consistency"] = {
        "partition_sum": sum(v for v in partition_counts.values() if v is not None),
        "partition_counts_all_available": all(v is not None for v in partition_counts.values()),
        "genomic_total": counts.get("wgs_genomic"),
        "platform_sum": sum((grouped.get("genomic_wgs_by_platform") or {}).values()),
        "source_sum": sum((grouped.get("all_wgs_by_source") or {}).values()),
    }
    write_tsv(root / "pilot_runs.tsv", (*FIELDS, *DERIVED), rows)
    tax_rows = defaultdict(list)
    for row in rows:
        tax_rows[(row["tax_id"], row["scientific_name"], row["taxonomy_rank"])].append(row)
    write_tsv(root / "pilot_taxa.tsv", ("tax_id", "scientific_name", "rank", "pilot_run_count", "pilot_fastq_bytes_known", "scope"),
              ({"tax_id": key[0], "scientific_name": key[1], "rank": key[2], "pilot_run_count": len(values),
                "pilot_fastq_bytes_known": sum(int(r["fastq_total_bytes"] or 0) for r in values),
                "scope": "pilot only; not archive-wide species frequency"} for key, values in sorted(tax_rows.items())))
    group_rows = []
    for label, values in grouped.items():
        for value, count in (values or {}).items():
            group_rows.append({"query_group": label, "value": value, "run_count": count, "scope": "full query count endpoint"})
    write_tsv(root / "archive_group_counts.tsv", ("query_group", "value", "run_count", "scope"), group_rows)
    json_write(root / "summary.json", summary)
    json_write(root / "taxonomy.json", taxa)
    readme = f"""# Public eukaryote WGS inventory pilot

Collected {summary['finished_at_utc']}. Status: {summary['status']}.

Archive-wide count endpoint: WGS all sources = {counts.get('wgs_all_sources')}; WGS + GENOMIC = {counts.get('wgs_genomic')} run records. These are public ENA-indexed records, not species, independent samples, or verified suitable libraries. The requests are live, not an atomic snapshot.

Detailed metadata: **{len(rows)} runs only**, capped at {5 * args.per_partition} across five disjoint taxonomy queries. The first archive-ordered records in each query are not a representative random sample. Species/taxon frequencies, file sizes and bases in pilot tables apply **only to this pilot**. No archive-wide species count, volume, cost, or runtime is extrapolated.

Files: `summary.json` contains full-query counts and explicitly nested pilot-only statistics; `archive_group_counts.tsv` contains full-query platform/source counts; `pilot_runs.tsv` contains per-run metadata and per-format sizes; `pilot_taxa.tsv` contains pilot coverage with explicit taxonomy rank; `taxonomy.json` stores taxon lookup results; `provenance.json` and `raw/` retain URLs, UTC times, byte counts and SHA-256 response evidence. `plan.json` records query limits. Errors, if any, are in `summary.json`.

## Size and identity rules

* ENA `fastq_bytes` lists individual FASTQ object sizes in bytes; semicolon entries are summed within a run. Missing sizes remain unknown. ENA `sra_bytes`, original `submitted_bytes`, sequenced `base_count`, and NCBI SRA metadata `mbytes` are different quantities. Never add alternative file representations or convert bases directly to download bytes.
* ENA size metadata do not verify AWS object size or availability. NCBI's `mbytes` describes megabytes of data in the SRA run, not FASTQ bytes. Cost estimates need the chosen actual representation and file locations.
* A tax_id can identify a strain, subspecies or unresolved taxon. Only verified `rank=species` IDs contribute to `rank_species_tax_ids`; this is a pilot lower-bound coverage measure, not a complete species census. The script does not roll strains up to species.
* Distinct run/BioSample/experiment IDs do not establish biological independence; repeated submissions and shared libraries require additional checking.
* GENOMIC/WGS metadata alone do not rule out Hi-C, WGA, capture, RNA mislabelling, contamination, or poor telomere retention. Protocol/QC review remains required.

## Future AWS planning

The official registry lists public SRA Normalized full-quality data in `sra-pub-run-odp` (us-east-1); these require SRA Toolkit conversion. Public original FASTQ/BAM/CRAM collections are subsets. ODP access is anonymous/free under the linked NCBI guidance; other original-file delivery paths have separate restrictions. No AWS deployment, paid Athena query, sequence download, or object availability probe was performed.

Next collect a complete, deduplicated metadata inventory, select suitable representative libraries per taxon, resolve taxonomic rank, verify target objects, and benchmark transfer + conversion + TRM + scratch storage on AWS. Keep novel/negative follow-up and optional genome estimation separate. Cost cannot be inferred from this biased 1000-run sample.

## Reproduce

`python3 scripts/inventory_public_wgs.py --fetch --output NEW_DIRECTORY --per-partition {args.per_partition} --taxonomy-limit {args.taxonomy_limit}`

Without `--fetch`, only a query plan is written. The per-partition cap cannot exceed 200 (1000 total); each response is capped at 2 MiB and all responses at 8 MiB. No offset pagination is assumed; disjoint queries are used. Failed operations are recorded, not converted into zero counts.

## Official sources

""" + "\n".join(f"* [{key}]({url})" for key, url in SOURCES.items()) + "\n"
    (root / "README.md").write_text(readme, encoding="utf-8")
    print(json.dumps({"status": summary["status"], "output": str(root), "archive_counts": counts,
                      "pilot_runs": len(rows), "response_bytes": fetcher.total, "errors": len(errors)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--fetch", action="store_true", help="make bounded public metadata HTTP requests")
    parser.add_argument("--per-partition", type=int, default=200)
    parser.add_argument("--taxonomy-limit", type=int, default=80)
    parser.add_argument("--timeout", type=int, default=60)
    args = parser.parse_args()
    if not 1 <= args.per_partition <= 200 or not 0 <= args.taxonomy_limit <= 100 or not 1 <= args.timeout <= 120:
        parser.error("require per-partition 1..200, taxonomy-limit 0..100, timeout 1..120")
    run(args)


if __name__ == "__main__":
    main()
