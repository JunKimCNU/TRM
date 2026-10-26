#!/usr/bin/env python3
"""Summarize completed TRM batch artifacts without reopening any FASTQ input.

Safe to run while a batch is in progress. Missing/failed/running samples remain
visible and their unfinished artifacts do not enter the candidate table.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import io
import json
import math
from pathlib import Path
import subprocess
import sys


try:
    from batch_common import atomic_text
except ModuleNotFoundError as error:
    if error.name != "batch_common":
        raise
    # Direct file loaders need not place scripts/ on sys.path.
    import importlib.util
    _spec = importlib.util.spec_from_file_location("_trm_batch_common", Path(__file__).with_name("batch_common.py"))
    _common = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_common)
    atomic_text = _common.atomic_text


CANDIDATE_FIELDS = [
    "sample_id", "rank", "motif", "reverse_complement", "period", "reads",
    "repeat_bases", "max_tract", "long_tract_reads", "status", "rank_score",
    "distinct_bases", "forward_reads", "reverse_reads", "ambiguous_strand_reads",
    "full_repeat_reads", "left_end_reads", "right_end_reads", "long_left_end_reads",
    "long_right_end_reads", "opposed_end_records", "result_json",
    "prior_class", "prior_confidence_tier", "prior_matched_reference", "prior_edit_distance",
    "prior_normalized_distance", "prior_eligible", "prior_source_lineage",
    "prior_reference_evidence",
    "candidate_eligible", "support_status", "support_reads_per_million",
    "support_repeat_bases_per_million", "support_read_enrichment_lower", "support_repeat_base_enrichment",
]
SAMPLE_FIELDS = [
    "sample_id", "run_status", "records", "bases", "input_compressed_bytes",
    "wall_seconds", "wall_source", "sequence_Mbps", "compressed_MBps",
    "user_cpu_seconds", "system_cpu_seconds", "max_rss_kb", "candidate_count",
    "top_motif", "top_reverse_complement", "top_status", "top_reads",
    "top_repeat_bases", "top_max_tract", "top_long_tract_reads", "top_rank_score",
    "max_record_bases", "sampled_prefix", "result_json", "metrics_json", "error",
    "eligible_telomere_candidates", "rejected_repeat_candidates", "top_prior_class",
    "top_prior_matched_reference", "top_prior_confidence_tier",
    "sequence_prior_candidates", "joint_eligible_candidates", "top_support_status",
    "top_support_read_enrichment_lower", "top_support_repeat_base_enrichment",
]


def read_json(path: Path) -> tuple[dict, str]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(value, dict):
            raise ValueError("expected a JSON object")
        return value, ""
    except FileNotFoundError:
        return {}, ""
    except (OSError, ValueError) as error:
        return {}, f"{path.name}: {error}"


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n",
                            extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    atomic_text(path, stream.getvalue())


def number(value: object) -> float | int | None:
    if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value):
        return value
    return None


def positive(value: object) -> float | int | None:
    value = number(value)
    return value if value is not None and value > 0 else None


def iso_time(value: object) -> datetime | None:
    if not isinstance(value, str):
        return None
    try:
        timestamp = datetime.fromisoformat(value.replace("Z", "+00:00"))
        return timestamp.replace(tzinfo=timezone.utc) if timestamp.tzinfo is None else timestamp
    except ValueError:
        return None


def display(value: object, digits: int = 0) -> str:
    if value is None or value == "":
        return "NA"
    numeric = number(value)
    if numeric is not None:
        return f"{numeric:,.{digits}f}"
    return str(value).replace("|", "\\|").replace("\n", " ")


def manifest_ids(path: Path) -> list[str]:
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if not reader.fieldnames or "sample_id" not in reader.fieldnames:
            raise ValueError(f"manifest lacks sample_id column: {path}")
        identifiers = [row["sample_id"] for row in reader]
    for sample_id in identifiers:
        if not sample_id or Path(sample_id).name != sample_id or sample_id in (".", ".."):
            raise ValueError(f"unsafe sample_id: {sample_id!r}")
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("manifest has duplicate sample_id values")
    return identifiers


def collect_sample(root: Path, sample_id: str, previous: dict) -> tuple[dict, list[dict]]:
    destination = root / sample_id
    state, state_error = read_json(destination / "state.json")
    metrics, metrics_error = read_json(destination / "metrics.json")
    if not metrics and isinstance(state.get("metrics"), dict):
        metrics = state["metrics"]
    errors = [error for error in (state_error, metrics_error) if error]
    status = state.get("status") or previous.get("status") or "pending"
    if status == "resumed":
        status = "complete"
    result = {}
    result_path = destination / "result.json"
    # A previous result may still exist during a fresh scan; never report it then.
    if status == "complete" or (status == "pending" and result_path.is_file() and not state_error):
        result, result_error = read_json(result_path)
        if result_error:
            errors.append(result_error)
        summary = result.get("summary")
        candidates = result.get("candidates")
        if (not isinstance(summary, dict) or not isinstance(candidates, list)
                or number(summary.get("records")) is None or number(summary.get("bases")) is None
                or any(not isinstance(candidate, dict) for candidate in candidates)):
            status = "invalid_result"
            errors.append("missing or invalid completed result.json")
            result = {}
        else:
            status = "complete"
    if state_error and status == "pending":
        status = "invalid_state"
    if status == "failed":
        errors.append(str(state.get("error") or previous.get("error") or "see stderr.log"))
    if status not in ("complete", "failed"):
        metrics = {}  # An old timing file can remain while a replacement scan runs.
    summary = result.get("summary", {})
    candidates = result.get("candidates", [])
    top = candidates[0] if candidates else {}
    input_bytes = number(metrics.get("input_compressed_bytes"))
    if input_bytes is None and result.get("inputs"):
        sizes = [number(item.get("bytes")) for item in result["inputs"]]
        if all(size is not None for size in sizes):
            input_bytes = sum(sizes)
    elapsed = positive(metrics.get("wall_seconds"))
    wall_source = "measured_subprocess" if elapsed is not None else ""
    if elapsed is None and status == "complete":
        elapsed = positive(summary.get("elapsed_seconds"))
        wall_source = "scanner_internal" if elapsed is not None else ""
    bases = number(summary.get("bases"))
    row = {
        "sample_id": sample_id, "run_status": status,
        "records": summary.get("records", ""), "bases": summary.get("bases", ""),
        "input_compressed_bytes": input_bytes, "wall_seconds": elapsed, "wall_source": wall_source,
        "sequence_Mbps": bases / elapsed / 1e6
            if status == "complete" and bases is not None and elapsed else None,
        "compressed_MBps": input_bytes / elapsed / 1e6
            if status == "complete" and input_bytes is not None and elapsed else None,
        "user_cpu_seconds": number(metrics.get("user_cpu_seconds")),
        "system_cpu_seconds": number(metrics.get("system_cpu_seconds")),
        "max_rss_kb": number(metrics.get("max_rss_kb")),
        "candidate_count": len(candidates) if result else "",
        "max_record_bases": summary.get("max_record_bases", ""),
        "sampled_prefix": summary.get("sampled_prefix", ""),
        "result_json": f"{sample_id}/result.json" if result else "",
        "metrics_json": f"{sample_id}/metrics.json" if (destination / "metrics.json").is_file() else "",
        "error": "; ".join(errors),
        "eligible_telomere_candidates": summary.get("eligible_telomere_candidates", ""),
        "rejected_repeat_candidates": summary.get("rejected_repeat_candidates", ""),
        "sequence_prior_candidates": summary.get("sequence_prior_candidates", ""),
        "joint_eligible_candidates": summary.get("joint_eligible_candidates", ""),
    }
    for field in ("motif", "reverse_complement", "status", "reads", "repeat_bases",
                  "max_tract", "long_tract_reads", "rank_score"):
        row["top_" + field] = top.get(field, "")
    all_candidates = [dict(candidate, sample_id=sample_id, rank=rank, result_json=row["result_json"])
                      for rank, candidate in enumerate(candidates, start=1)]
    for field in ("class", "matched_reference", "confidence_tier"):
        row["top_prior_" + field] = top.get("telomere_prior", {}).get(field, "")
    for field in ("status", "read_enrichment_lower", "repeat_base_enrichment"):
        row["top_support_" + field] = top.get("support", {}).get(field, "")
    for candidate in all_candidates:
        for field in ("class", "confidence_tier", "matched_reference", "edit_distance",
                      "normalized_distance", "eligible", "source_lineage", "reference_evidence"):
            candidate["prior_" + field] = candidate.get("telomere_prior", {}).get(field, "")
        for field in ("status", "reads_per_million", "repeat_bases_per_million",
                      "read_enrichment_lower", "repeat_base_enrichment"):
            candidate["support_" + field] = candidate.get("support", {}).get(field, "")
    return row, all_candidates


def generate(root: Path, manifest: Path | None = None) -> dict:
    generated = datetime.now(timezone.utc)
    summary, summary_error = read_json(root / "summary.json")
    metadata = {}
    metadata_errors = []
    for name in ("run.json", "batch_state.json", "run_metadata.json"):
        candidate, error = read_json(root / name)
        if candidate:
            metadata.update(candidate)
        if error:
            metadata_errors.append(error)
    identifiers = []
    manifest_value = manifest or metadata.get("manifest") or summary.get("manifest")
    if not manifest_value and (root / "manifest.tsv").is_file():
        manifest_value = root / "manifest.tsv"
    if manifest_value:
        manifest_path = Path(manifest_value)
        if not manifest_path.is_absolute():
            manifest_path = root / manifest_path
        try:
            identifiers = manifest_ids(manifest_path)
        except (OSError, ValueError) as error:
            if manifest is not None:
                raise
            metadata_errors.append(str(error))
    previous = {row["sample_id"]: row for row in summary.get("samples", [])
                if isinstance(row, dict) and isinstance(row.get("sample_id"), str)}
    if not identifiers:
        identifiers = list(previous)
    known = set(identifiers)
    for path in sorted(root.iterdir()):
        if (path.is_dir() and not path.name.startswith(".") and path.name not in known
                and any((path / name).is_file() for name in ("state.json", "result.json", "metrics.json"))):
            identifiers.append(path.name)
            known.add(path.name)
    rows, candidates = [], []
    for sample_id in identifiers:
        row, sample_candidates = collect_sample(root, sample_id, previous.get(sample_id, {}))
        rows.append(row)
        candidates.extend(sample_candidates)
    complete = [row for row in rows if row["run_status"] == "complete"]
    counts = Counter(row["run_status"] for row in rows)
    total_bases = sum(row["bases"] for row in complete)
    total_records = sum(row["records"] for row in complete)
    byte_rows = [row for row in complete if number(row["input_compressed_bytes"]) is not None]
    total_bytes = sum(row["input_compressed_bytes"] for row in byte_rows)
    timed = [row for row in complete if positive(row["wall_seconds"]) is not None]
    sample_seconds = sum(row["wall_seconds"] for row in timed)
    timing = summary if summary and not summary_error else metadata
    started = iso_time(timing.get("started_utc"))
    finished = iso_time(timing.get("finished_utc"))
    # summary.json describes a finished invocation; its wall includes resume checks.
    # Never divide previously scanned bases by an invocation that only resumed them.
    summary_has_resumed = any(row.get("status") == "resumed" for row in previous.values())
    batch_seconds = positive(timing.get("wall_seconds"))
    batch_time_source = "summary_invocation" if timing is summary else "run_metadata"
    if started and not finished and timing.get("finished") is False:
        batch_seconds = max(0, (generated - started).total_seconds())
        batch_time_source = "live_timestamp"
    elif batch_seconds is None and started:
        batch_seconds = max(0, ((finished or generated) - started).total_seconds())
        batch_time_source = "metadata_timestamps"
    if batch_seconds is None:
        batch_time_source = ""
    can_batch_rate = bool(batch_seconds and not summary_has_resumed)
    totals = {
        "generated_utc": generated.isoformat(), "output_dir": str(root),
        "manifest": str(manifest_value) if manifest_value else None,
        "sample_count": len(rows), "status_counts": dict(counts),
        "completed_samples": len(complete), "failed_samples": counts["failed"],
        "running_samples": counts["running"], "pending_samples": counts["pending"],
        "completed_records": total_records, "completed_bases": total_bases,
        "completed_input_compressed_bytes": total_bytes,
        "samples_with_input_bytes": len(byte_rows), "samples_with_elapsed_time": len(timed),
        "sum_sample_wall_seconds": sample_seconds, "batch_wall_seconds": batch_seconds,
        "batch_wall_source": batch_time_source,
        "batch_sequence_Mbps": total_bases / batch_seconds / 1e6 if can_batch_rate else None,
        "batch_compressed_MBps": total_bytes / batch_seconds / 1e6 if can_batch_rate else None,
        "sample_time_weighted_sequence_Mbps":
            sum(row["bases"] for row in timed) / sample_seconds / 1e6 if sample_seconds else None,
        "sample_time_weighted_compressed_MBps":
            sum(row["input_compressed_bytes"] or 0 for row in timed) / sample_seconds / 1e6
            if sample_seconds else None,
        "max_sample_rss_kb": max((row["max_rss_kb"] for row in complete
                                  if row["max_rss_kb"] is not None), default=None),
        "sum_sample_user_cpu_seconds": sum(row["user_cpu_seconds"] or 0 for row in complete),
        "sum_sample_system_cpu_seconds": sum(row["system_cpu_seconds"] or 0 for row in complete),
        "samples_with_rss": sum(row["max_rss_kb"] is not None for row in complete),
        "candidate_rows": len(candidates), "includes_resumed_samples": summary_has_resumed,
        "metadata_errors": ([summary_error] if summary_error else []) + metadata_errors,
    }
    lines = [
        "# TRM batch result and speed report", "",
        f"Generated: {generated.isoformat()}. Output: `{root}`.", "",
        f"Samples: **{len(complete)}/{len(rows)} complete**, {counts['failed']} failed, "
        f"{counts['running']} running, {counts['pending']} pending. "
        f"Other states: {sum(v for k, v in counts.items() if k not in {'complete', 'failed', 'running', 'pending'})}.", "",
        "This is a point-in-time snapshot. Only completed, readable result files enter the totals and "
        "candidate table; running samples may finish immediately after this report is generated.", "",
        "| Quantity | Value |", "| --- | ---: |",
        f"| Completed FASTQ records (mates counted separately) | {total_records:,} |",
        f"| Completed sequence bases | {total_bases:,} ({total_bases / 1e9:,.3f} Gb) |",
        f"| Completed input bytes | {total_bytes:,} ({total_bytes / 1e9:,.3f} GB) |",
        f"| Batch wall time (s) | {display(batch_seconds, 3)} |",
        f"| Sum of sample wall times (s) | {sample_seconds:,.3f} |",
        f"| Sum of sample user / system CPU times (s) | {totals['sum_sample_user_cpu_seconds']:,.3f} / {totals['sum_sample_system_cpu_seconds']:,.3f} |",
        f"| Batch sequence throughput (Mb/s) | {display(totals['batch_sequence_Mbps'], 3)} |",
        f"| Batch input throughput (MB/s) | {display(totals['batch_compressed_MBps'], 3)} |",
        f"| Sample-time weighted sequence throughput (Mb/s) | {display(totals['sample_time_weighted_sequence_Mbps'], 3)} |",
        f"| Highest individual sample RSS (MiB) | {display(totals['max_sample_rss_kb'] / 1024 if totals['max_sample_rss_kb'] is not None else None, 2)} |",
        f"| Reported candidate rows | {len(candidates):,} |", "",
        "Mb and MB use decimal units (10^6). Sequence throughput counts nucleotide bases; input throughput "
        "counts stored input bytes, compressed for .gz files. Batch throughput divides completed work by "
        "elapsed batch time and reflects parallelism. The sum of per-sample elapsed times can exceed batch "
        "wall time; its corresponding throughput measures a time-weighted individual scan. During execution, "
        "batch throughput excludes work within unfinished samples. RSS is the maximum of individual scanner "
        "peaks, not the simultaneous memory total. Missing measurements appear as NA.", "",
    ]
    if summary_has_resumed:
        lines += ["Batch throughput is NA because the finished invocation includes resumed samples; its "
                  "wall time does not measure their original scans. Per-sample metrics preserve original scans.", ""]
    if any(row["wall_source"] == "scanner_internal" for row in complete):
        lines += ["Some samples lack external timing; their recorded scanner-internal elapsed time is used "
                  "and identified in `sample_report.tsv` (`wall_source=scanner_internal`).", ""]
    lines += [
        "## Candidate interpretation", "",
        "Motifs are canonicalized across rotations and reverse complements. Their reported rank is retained "
        "from TRM; the top-ranked repeat is a candidate, not proof of telomeric location. Internal tandem "
        "repeats and satellites can dominate a WGS sample. Forward/reverse counts describe input-read "
        "orientation and do not infer chromosome-end direction.", "",
        "`long_tract_reads` counts records containing a single tract meeting the configured `min_tract` "
        "threshold (normally 500 bp). Reads shorter than this threshold cannot establish that contiguous "
        "length. A lack of candidates or long tracts does not establish telomeric-motif absence or ALT. "
        "Mates and PCR duplicates are counted as records and are not independent molecules. Motif period "
        "is repeat-unit length; distinct_bases separately describes its nucleotide composition.", "",
        "## Per-sample results", "",
        "| Sample | Run status | Top motif | Candidate status | Supporting records | Max tract (bp) | Long-tract records | Total records | Sequence (Mb) | Wall (s) | RSS (MiB) |",
        "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        label = f"[{row['sample_id']}]({row['result_json']})" if row["result_json"] else row["sample_id"]
        values = [label, row["run_status"], row["top_motif"], row["top_status"], row["top_reads"],
                  row["top_max_tract"], row["top_long_tract_reads"], row["records"]]
        lines.append("| " + " | ".join(display(value) for value in values) + " | " + " | ".join([
            display(row["bases"] / 1e6 if number(row["bases"]) is not None else None, 3),
            display(row["wall_seconds"], 3),
            display(row["max_rss_kb"] / 1024 if row["max_rss_kb"] is not None else None, 2),
        ]) + " |")
    errors = [f"- {row['sample_id']}: {row['error']}" for row in rows if row["error"]]
    errors += [f"- Report metadata: {error}" for error in totals["metadata_errors"]]
    if errors:
        lines += ["", "## Missing artifacts and errors", "", *errors]
    lines += [
        "", "## Saved artifacts", "",
        "- [All motif candidates](all_candidates.tsv): every reported candidate, original rank, orientation "
        "counts, tract evidence, and relative source path.",
        "- [Per-sample results and resource metrics](sample_report.tsv): includes pending and failed samples.",
        "- [Machine-readable totals](report_totals.json): counts, timestamp, speed definitions, missing metadata.",
        *(["- [Live input progress](live_progress.json): optional process-based byte-position estimates "
           "for unfinished samples; distinct from completed scan results."]
          if (root / "live_progress.json").is_file() else []),
        "- `SAMPLE/result.json` and `SAMPLE/result.tsv`: original scanner outputs; JSON includes scan "
        "configuration, source inputs, warnings, and example tracts.",
        "- `SAMPLE/metrics.json`, `SAMPLE/resources.tsv`, `SAMPLE/state.json`, `SAMPLE/stdout.log`, "
        "and `SAMPLE/stderr.log`: timing, resources, execution provenance, and logs when available.",
        "- `summary.json` and `summary.tsv`: batch launcher summary once its invocation finishes.", "",
        "These links are relative so this report can be moved together with its output directory. "
        "The original FASTQ data are not embedded in the report.", "",
    ]
    write_tsv(root / "all_candidates.tsv", CANDIDATE_FIELDS, candidates)
    write_tsv(root / "sample_report.tsv", SAMPLE_FIELDS, rows)
    atomic_text(root / "report_totals.json", json.dumps(totals, indent=2, sort_keys=True) + "\n")
    atomic_text(root / "run_report.md", "\n".join(lines))
    return totals


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True, help="existing batch output root")
    parser.add_argument("--manifest", type=Path, help="include all expected samples, even before they start")
    args = parser.parse_args()
    root = args.output_dir.resolve()
    if not root.is_dir():
        parser.error(f"batch output directory does not exist: {root}")
    monitor = Path(__file__).with_name("monitor_batch.py")
    if monitor.is_file():
        try:
            subprocess.run([sys.executable, str(monitor), "--output-dir", str(root)],
                           check=False, capture_output=True, timeout=10)
        except (OSError, subprocess.TimeoutExpired):
            pass  # Optional live telemetry must never prevent artifact reporting.
    try:
        totals = generate(root, args.manifest.resolve() if args.manifest else None)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Report: {root / 'run_report.md'}")
    print(f"Completed {totals['completed_samples']}/{totals['sample_count']} samples; "
          f"{totals['candidate_rows']} candidate rows; {totals['completed_bases']} bases")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
