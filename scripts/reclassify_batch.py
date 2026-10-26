#!/usr/bin/env python3
"""Apply telomere priors to saved batch results without opening sequence inputs.

The source batch stays intact. Output must be a new or empty directory disjoint
from the source batch. Optional priors TSV columns are sample_id, lineage,
related_motifs (comma separated), and source (citation or provenance).
Elapsed time and resource usage describe reclassification only; original scan
metrics are retained separately and are never used as postprocessing throughput.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time


try:
    from batch_common import positive_int, sha256
except ModuleNotFoundError as error:
    if error.name != "batch_common":
        raise
    # Direct file loaders need not place scripts/ on sys.path.
    import importlib.util
    _spec = importlib.util.spec_from_file_location("_trm_batch_common", Path(__file__).with_name("batch_common.py"))
    _common = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_common)
    positive_int = _common.positive_int
    sha256 = _common.sha256


SAFE_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}\Z")
SAMPLE_FIELDS = [
    "sample_id", "status", "old_top_motif", "new_top_motif", "new_top_prior_class",
    "new_top_reference_evidence",
    "new_top_support_status", "new_top_reads_per_million", "new_top_read_enrichment_lower",
    "new_top_repeat_base_enrichment", "sequence_prior_candidate_count",
    "new_top_reads", "new_top_max_tract", "new_top_status", "top_changed", "call",
    "raw_candidate_count", "accepted_candidate_count", "rejected_candidate_count",
    "records", "bases", "lineage", "related_motifs", "prior_source",
    "reclassification_wall_seconds", "user_cpu_seconds", "system_cpu_seconds", "max_rss_kb",
    "original_scan_wall_seconds", "input_artifact_bytes", "original_result_sha256",
    "source_result_json", "result_json", "error",
]
CANDIDATE_FIELDS = [
    "sample_id", "rank", "motif", "reverse_complement", "period", "reads", "repeat_bases",
    "max_tract", "long_tract_reads", "status", "rank_score", "prior_class", "prior_confidence_tier",
    "prior_matched_reference", "prior_edit_distance", "prior_normalized_distance",
    "prior_source_lineage", "prior_evidence_source", "prior_reference_evidence",
    "telomere_prior_json", "result_json",
    "candidate_eligible", "support_status", "support_reads_per_million",
    "support_repeat_bases_per_million", "support_read_enrichment_lower",
    "support_repeat_base_enrichment", "support_reasons",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, value: object) -> None:
    descriptor, temporary = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            json.dump(value, output, indent=2, sort_keys=True)
            output.write("\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def read_json(path: Path, optional: bool = False) -> dict:
    if optional and not path.exists():
        return {}
    with path.open(encoding="utf-8") as source:
        value = json.load(source)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def load_samples(root: Path) -> tuple[list[dict], dict]:
    summary = read_json(root / "summary.json", optional=True)
    previous = {row["sample_id"]: row for row in summary.get("samples", [])}
    if (root / "manifest.tsv").is_file():
        with (root / "manifest.tsv").open(encoding="utf-8-sig", newline="") as source:
            reader = csv.DictReader(source, delimiter="\t")
            if not reader.fieldnames or "sample_id" not in reader.fieldnames:
                raise ValueError("source manifest.tsv lacks sample_id column")
            identifiers = [row["sample_id"] for row in reader]
    else:
        identifiers = list(previous)
    if not identifiers:
        raise ValueError("source batch has no samples in manifest.tsv or summary.json")
    if any(not isinstance(item, str) or not SAFE_ID.fullmatch(item) for item in identifiers):
        raise ValueError("source batch contains an unsafe sample_id")
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("source batch contains duplicate sample_id values")
    samples = []
    for sample_id in identifiers:
        state = read_json(root / sample_id / "state.json", optional=True)
        status = state.get("status") or previous.get(sample_id, {}).get("status", "pending")
        samples.append({"sample_id": sample_id, "status": status, "state": state,
                        "previous": previous.get(sample_id, {})})
    return samples, summary


def load_priors(path: Path | None, identifiers: set[str]) -> dict:
    if path is None:
        return {}
    priors = {}
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if not reader.fieldnames or set(reader.fieldnames) != {
                "sample_id", "lineage", "related_motifs", "source"}:
            raise ValueError("priors TSV requires sample_id, lineage, related_motifs, source columns")
        if len(reader.fieldnames) != 4:
            raise ValueError("priors TSV contains duplicate columns")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"priors TSV line {reader.line_num}: incorrect number of fields")
            sample_id = row["sample_id"]
            if sample_id not in identifiers:
                raise ValueError(f"priors TSV has unknown sample_id: {sample_id!r}")
            if sample_id in priors:
                raise ValueError(f"priors TSV has duplicate sample_id: {sample_id!r}")
            motifs = [value.strip() for value in row["related_motifs"].split(",")] \
                if row["related_motifs"].strip() else []
            if any(not motif for motif in motifs):
                raise ValueError(f"priors TSV line {reader.line_num}: empty related motif")
            if row["lineage"].strip() and not motifs:
                raise ValueError(f"priors TSV line {reader.line_num}: lineage requires related_motifs; "
                                 "a taxon name alone does not supply phylogenetic evidence")
            # Sequence validity and supported lengths belong to the Rust engine.
            priors[sample_id] = {"lineage": row["lineage"].strip(), "related_motifs": motifs,
                                 "source": row["source"].strip()}
    return priors


def prior_class(candidate: dict) -> str:
    prior = candidate.get("telomere_prior") or {}
    return str(prior.get("class", prior.get("classification", prior.get("support", ""))))


def candidate_evidence(candidate: dict) -> dict:
    return {key: value for key, value in candidate.items()
            if key not in ("telomere_prior", "support", "candidate_eligible")}


def eligible(candidate: dict) -> bool:
    return bool(candidate.get("candidate_eligible", candidate.get("telomere_prior", {}).get("eligible", False)))


def all_candidates(result: dict):
    """Discovery output can repeat unsupported rows in its rejected subset."""
    seen = {}
    for key in ("candidates", "rejected_candidates"):
        for candidate in result.get(key, []):
            motif = candidate["motif"]
            if motif in seen:
                if candidate_evidence(candidate) != candidate_evidence(seen[motif]):
                    raise ValueError(f"conflicting duplicate evidence for motif {motif}")
                continue
            seen[motif] = candidate
            yield candidate


def read_resources(path: Path) -> dict:
    fields = path.read_text(encoding="utf-8").splitlines()[-1].split("\t")
    if len(fields) != 5:
        raise ValueError(f"invalid GNU time resource output: {path}")
    return dict(zip(("resource_wall_seconds", "user_cpu_seconds", "system_cpu_seconds",
                     "max_rss_kb", "resource_exit_code"),
                    (float(fields[0]), float(fields[1]), float(fields[2]), int(fields[3]), int(fields[4]))))


def reclassify_sample(sample: dict, args: argparse.Namespace, priors: dict) -> dict:
    sample_id = sample["sample_id"]
    source = args.source_batch / sample_id
    destination = args.output_dir / sample_id
    destination.mkdir()
    prior = priors.get(sample_id, {"lineage": "", "related_motifs": [], "source": ""})
    row = {"sample_id": sample_id, "lineage": prior["lineage"],
           "related_motifs": ",".join(prior["related_motifs"]), "prior_source": prior["source"],
           "source_result_json": str(source / "result.json")}
    metrics = {"operation": "reclassification", "started_utc": utc_now(),
               "sequence_inputs_opened": False, "command": [], "exit_code": None}
    provenance = {"source_batch": str(args.source_batch), "source_status": sample["status"],
                  "source_state": sample["state"], "phylogenetic_prior": prior}
    started = time.perf_counter()
    try:
        scan_metrics = read_json(source / "metrics.json", optional=True) \
            or sample["state"].get("metrics", {})
        atomic_json(destination / "original_scan_metrics.json", scan_metrics)
        row["original_scan_wall_seconds"] = scan_metrics.get("wall_seconds")
        if sample["status"] not in ("complete", "resumed"):
            row.update(status="source_failed" if sample["status"] == "failed" else "source_incomplete",
                       error=sample["state"].get("error") or sample["previous"].get("error")
                       or f"original sample status: {sample['status']}")
            return row
        result_path = source / "result.json"
        original_sha = sha256(result_path)
        row.update(original_result_sha256=original_sha, input_artifact_bytes=result_path.stat().st_size)
        provenance.update(original_result_sha256=original_sha,
                          source_result_json=str(result_path), input_artifact_bytes=row["input_artifact_bytes"])
        expected_sha = sample["state"].get("result_sha256")
        if expected_sha and expected_sha != original_sha:
            row["status"] = "source_invalid"
            raise ValueError("original result checksum differs from source state.json")
        original = read_json(result_path)
        original_candidates = original["candidates"]
        row["old_top_motif"] = original_candidates[0]["motif"] if original_candidates else ""
        original_evidence = {item["motif"]: candidate_evidence(item) for item in all_candidates(original)}
        row["raw_candidate_count"] = len(original_evidence)
        row["records"] = original["summary"]["records"]
        row["bases"] = original["summary"]["bases"]
        del original, original_candidates
        command = [args.binary, "reclassify", str(result_path), "--output", str(destination / "result"),
                   "--prior-policy", args.prior_policy, "--prior-max-edits", str(args.prior_max_edits),
                   "--prior-max-distance", str(args.prior_max_distance)]
        if prior["lineage"]:
            command.extend(("--lineage", prior["lineage"]))
        for motif in prior["related_motifs"]:
            command.extend(("--related-motif", motif))
        if args.background_model:
            command.extend(("--background-model", str(args.background_model),
                            "--support-policy", args.support_policy,
                            "--support-min-enrichment", str(args.support_min_enrichment)))
        metrics["command"] = command
        launched = command
        if args.time_binary:
            launched = [args.time_binary, "-f", "%e\t%U\t%S\t%M\t%x", "-o",
                        str(destination / "resources.tsv"), *command]
        command_started = time.perf_counter()
        with (destination / "stdout.log").open("w", encoding="utf-8") as stdout, \
                (destination / "stderr.log").open("w", encoding="utf-8") as stderr:
            completed = subprocess.run(launched, stdout=stdout, stderr=stderr, check=False)
        metrics.update(command_wall_seconds=time.perf_counter() - command_started,
                       exit_code=completed.returncode)
        if args.time_binary:
            metrics.update(read_resources(destination / "resources.tsv"))
        if completed.returncode:
            diagnostic = (destination / "stderr.log").read_text(encoding="utf-8", errors="replace")
            raise RuntimeError(f"reclassify exited {completed.returncode}: {diagnostic.strip()[-2000:]}")
        if sha256(result_path) != original_sha:
            raise ValueError("source result changed during reclassification")
        if args.background_model and sha256(args.background_model) != args.background_sha256:
            raise ValueError("background model changed during reclassification")
        result = read_json(destination / "result.json")
        accepted = result["candidates"]
        classified_count = 0
        eligible_count = 0
        for item in all_candidates(result):
            classified_count += 1
            eligible_count += eligible(item)
            if original_evidence.get(item["motif"]) != candidate_evidence(item):
                raise ValueError(f"original evidence changed for motif {item['motif']}")
        if classified_count != row["raw_candidate_count"]:
            raise ValueError("candidate counts changed: expected every raw candidate to be retained")
        if any(result["summary"].get(key) != row[key] for key in ("records", "bases")):
            raise ValueError("original read/base counts changed during reclassification")
        top = accepted[0] if accepted else {}
        row.update(status="complete", new_top_motif=top.get("motif", ""),
                   new_top_prior_class=prior_class(top), new_top_reads=top.get("reads", ""),
                   new_top_reference_evidence=json.dumps(top.get("telomere_prior", {}).get("reference_evidence", [])),
                   new_top_support_status=top.get("support", {}).get("status", "not_assessed"),
                   new_top_reads_per_million=top.get("support", {}).get("reads_per_million", ""),
                   new_top_read_enrichment_lower=top.get("support", {}).get("read_enrichment_lower", ""),
                   new_top_repeat_base_enrichment=top.get("support", {}).get("repeat_base_enrichment", ""),
                   sequence_prior_candidate_count=result["summary"].get("sequence_prior_candidates", eligible_count),
                   new_top_max_tract=top.get("max_tract", ""), new_top_status=top.get("status", ""),
                   top_changed=row["old_top_motif"] != top.get("motif", ""),
                   call=result["summary"].get("call", "candidates_only" if accepted else "no_call"),
                   accepted_candidate_count=eligible_count,
                   rejected_candidate_count=classified_count - eligible_count,
                   result_json=f"{sample_id}/result.json", error="")
        provenance["result_sha256"] = sha256(destination / "result.json")
    except (OSError, ValueError, KeyError, TypeError, RuntimeError) as error:
        row.setdefault("status", "reclassification_failed")
        row["error"] = str(error)
    finally:
        metrics.update(finished_utc=utc_now(), wall_seconds=time.perf_counter() - started)
        row["reclassification_wall_seconds"] = metrics["wall_seconds"]
        for key in ("user_cpu_seconds", "system_cpu_seconds", "max_rss_kb"):
            row[key] = metrics.get(key)
        atomic_json(destination / "metrics.json", metrics)
        atomic_json(destination / "provenance.json", provenance)
        atomic_json(destination / "state.json", {"status": row["status"], "error": row.get("error", ""),
                    "result_sha256": provenance.get("result_sha256")})
    return row


def write_reports(args: argparse.Namespace, rows: list[dict], metrics: dict) -> None:
    for name in ("summary.tsv", "sample_report.tsv"):
        with (args.output_dir / name).open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, SAMPLE_FIELDS, delimiter="\t", lineterminator="\n",
                                    extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    with (args.output_dir / "qualified_candidates.tsv").open("w", encoding="utf-8", newline="") as qualified, \
            (args.output_dir / "rejected_candidates.tsv").open("w", encoding="utf-8", newline="") as rejected:
        writers = [csv.DictWriter(handle, CANDIDATE_FIELDS, delimiter="\t", lineterminator="\n",
                                  extrasaction="ignore") for handle in (qualified, rejected)]
        for writer in writers:
            writer.writeheader()
        for row in rows:
            if row["status"] != "complete":
                continue
            result = read_json(args.output_dir / row["result_json"])
            ranks = [0, 0]
            for candidate in all_candidates(result):
                prior = candidate.get("telomere_prior", {})
                table = 0 if eligible(candidate) else 1
                ranks[table] += 1
                writers[table].writerow({**candidate, "sample_id": row["sample_id"], "rank": ranks[table],
                                        "prior_class": prior_class(candidate),
                                        **{"prior_" + key: prior.get(key) for key in
                                           ("confidence_tier", "matched_reference", "edit_distance",
                                            "normalized_distance", "source_lineage")},
                                        "prior_evidence_source": json.dumps(prior.get("evidence_source", [])),
                                        "prior_reference_evidence": json.dumps(prior.get("reference_evidence", [])),
                                        "telomere_prior_json": json.dumps(candidate.get("telomere_prior"),
                                                                           ensure_ascii=False, sort_keys=True),
                                        **{"support_" + key: candidate.get("support", {}).get(key) for key in
                                           ("status", "reads_per_million", "repeat_bases_per_million",
                                            "read_enrichment_lower", "repeat_base_enrichment")},
                                        "support_reasons": json.dumps(candidate.get("support", {}).get("reasons", [])),
                                        "result_json": row["result_json"]})
    counts = Counter(row["status"] for row in rows)
    complete = [row for row in rows if row["status"] == "complete"]
    metrics.update(status_counts=dict(counts), sample_count=len(rows), completed_samples=len(complete),
                   source_failed_samples=counts["source_failed"],
                   reclassification_failed_samples=counts["reclassification_failed"],
                   input_artifact_bytes=sum(row.get("input_artifact_bytes", 0) for row in complete),
                   raw_candidate_count=sum(row["raw_candidate_count"] for row in complete),
                   sequence_prior_candidate_count=sum(row["sequence_prior_candidate_count"] for row in complete),
                   accepted_candidate_count=sum(row["accepted_candidate_count"] for row in complete),
                   rejected_candidate_count=sum(row["rejected_candidate_count"] for row in complete),
                   top_changed_samples=sum(row["top_changed"] for row in complete),
                   no_call_samples=sum(not row["accepted_candidate_count"] for row in complete),
                   original_scan_sum_sample_wall_seconds=sum(row.get("original_scan_wall_seconds") or 0
                                                             for row in complete),
                   sum_user_cpu_seconds=sum(row.get("user_cpu_seconds") or 0 for row in complete),
                   sum_system_cpu_seconds=sum(row.get("system_cpu_seconds") or 0 for row in complete),
                   max_sample_rss_kb=max((row.get("max_rss_kb") or 0 for row in complete), default=0))
    metrics["wall_seconds"] = time.perf_counter() - args.started
    metrics["finished_utc"] = utc_now()
    metrics["input_artifact_MBps"] = metrics["input_artifact_bytes"] / metrics["wall_seconds"] / 1e6
    atomic_json(args.output_dir / "metrics.json", metrics)
    atomic_json(args.output_dir / "summary.json", {**metrics, "samples": rows})
    lines = ["# Telomeric motif candidate reclassification", "",
             f"Source batch: `{args.source_batch}`", "",
             f"Reclassified {len(complete)}/{len(rows)} samples in {metrics['wall_seconds']:.3f} seconds.",
             f"Original failed samples: {counts['source_failed']}; incomplete samples: {counts['source_incomplete']}; "
             f"invalid source artifacts: {counts['source_invalid']}; reclassification failures: {counts['reclassification_failed']}.",
             "", f"Raw candidates: {metrics['raw_candidate_count']:,}; accepted: {metrics['accepted_candidate_count']:,}; "
             f"rejected and retained: {metrics['rejected_candidate_count']:,}.",
             f"Top motif changed in {metrics['top_changed_samples']} samples; no accepted candidates in {metrics['no_call_samples']} samples.",
             "", "This operation reads saved JSON artifacts only. Original FASTQ inputs and batch outputs are unchanged.",
             "Postprocessing wall time includes artifact checks and report generation; it is not sequencing-data scan time.",
             f"Original batch scan wall time: {metrics.get('original_scan_batch_wall_seconds', 'unavailable')} seconds.",
             "Original per-sample scan metrics are saved as `original_scan_metrics.json` in each sample directory.",
             "", f"Sequence-prior candidates: {metrics['sequence_prior_candidate_count']:,}.",
             f"Background model: {metrics.get('background_model') or 'not supplied'}; support policy: {args.support_policy}; minimum enrichment: {args.support_min_enrichment}.",
             "Accepted candidates satisfy sequence-prior eligibility and, when an applicable background model was supplied, background support. These conditions do not prove chromosome-end location.",
             "In discovery or support-annotate mode, new_top_motif follows the primary result ranking and can be unsupported; candidate tables separate joint eligibility.",
             "Read/base fractions normalize sequence input amount, not genomic coverage. A no_call may reflect insufficient background-model resolution; use controls appropriate to the sample and read length.",
             "A no_call does not establish motif absence or a telomere maintenance mechanism.",
             "", "- [Sample comparison](sample_report.tsv)",
             "- [Qualified candidates](qualified_candidates.tsv)",
             "- [Rejected candidates](rejected_candidates.tsv)",
             "- [Reclassification metrics](metrics.json)", ""]
    (args.output_dir / "run_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    from runtime_binary import resolve_binary
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source_batch", type=Path, help="completed batch directory containing manifest.tsv and results")
    parser.add_argument("--binary", help="Explicit scanner override; default: TRM_BINARY or verified bin/trm")
    parser.add_argument("--output-dir", type=Path, required=True, help="new or empty directory disjoint from source batch")
    parser.add_argument("--priors", type=Path, help="optional TSV: sample_id, lineage, related_motifs, source")
    parser.add_argument("--jobs", type=positive_int, default=4)
    parser.add_argument("--prior-policy", choices=("canonical", "telomere", "discovery"), default="canonical")
    parser.add_argument("--prior-max-edits", type=int, default=2)
    parser.add_argument("--prior-max-distance", type=float, default=0.34)
    parser.add_argument("--background-model", type=Path)
    parser.add_argument("--support-policy", choices=("filter", "annotate"), default="annotate")
    parser.add_argument("--support-min-enrichment", type=float, default=2.0)
    args = parser.parse_args()
    args.started = time.perf_counter()
    try:
        if args.prior_max_edits < 0:
            raise ValueError("prior-max-edits must be nonnegative")
        if not math.isfinite(args.prior_max_distance) or not 0 <= args.prior_max_distance <= 1:
            raise ValueError("prior-max-distance must be a finite number between 0 and 1")
        if not math.isfinite(args.support_min_enrichment) or args.support_min_enrichment < 1:
            raise ValueError("support-min-enrichment must be finite and >=1")
        args.background_sha256 = None
        if args.background_model:
            args.background_model = args.background_model.resolve(strict=True)
            args.background_sha256 = sha256(args.background_model)
        args.source_batch = args.source_batch.resolve(strict=True)
        args.output_dir = args.output_dir.resolve()
        if (args.source_batch == args.output_dir or args.source_batch in args.output_dir.parents
                or args.output_dir in args.source_batch.parents):
            raise ValueError("output directory must be disjoint from the source batch (use a sibling directory)")
        if args.output_dir.exists() and (not args.output_dir.is_dir() or any(args.output_dir.iterdir())):
            raise ValueError("output directory must be new or empty; existing outputs are never overwritten")
        samples, source_summary = load_samples(args.source_batch)
        priors = load_priors(args.priors, {sample["sample_id"] for sample in samples})
        args.binary = str(resolve_binary(args.binary))
        args.time_binary = None
        if Path("/usr/bin/time").is_file():
            probe = subprocess.run(["/usr/bin/time", "--version"], capture_output=True, text=True, check=False)
            if probe.returncode == 0 and "GNU" in probe.stdout + probe.stderr:
                args.time_binary = "/usr/bin/time"
        args.output_dir.mkdir(parents=True, exist_ok=True)
        metrics = {"schema_version": 1, "operation": "reclassification", "started_utc": utc_now(),
                   "source_batch": str(args.source_batch), "output_dir": str(args.output_dir),
                   "sequence_inputs_opened": False, "original_outputs_modified": False,
                   "binary": args.binary, "binary_sha256": sha256(Path(args.binary)),
                   "script_sha256": sha256(Path(__file__).resolve()),
                   "jobs": args.jobs, "prior_policy": args.prior_policy, "prior_max_edits": args.prior_max_edits,
                   "prior_max_distance": args.prior_max_distance, "priors": priors,
                   "background_model": str(args.background_model) if args.background_model else None,
                   "background_model_sha256": args.background_sha256,
                   "support_policy": args.support_policy, "support_min_enrichment": args.support_min_enrichment,
                   "original_scan_batch_wall_seconds": source_summary.get("wall_seconds"),
                   "resource_measurement": "GNU time" if args.time_binary else "wall clock only"}
        atomic_json(args.output_dir / "run_metadata.json", metrics)
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            rows = list(executor.map(lambda sample: reclassify_sample(sample, args, priors), samples))
        write_reports(args, rows, metrics)
        print(f"Reclassified {metrics['completed_samples']}/{metrics['sample_count']} samples; "
              f"{metrics['accepted_candidate_count']} accepted; {metrics['rejected_candidate_count']} rejected; "
              f"{metrics['source_failed_samples']} original failures; {metrics['wall_seconds']:.3f} seconds")
        return 0 if all(row["status"] in ("complete", "source_failed") for row in rows) else 1
    except (OSError, ValueError, TypeError, KeyError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
