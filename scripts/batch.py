#!/usr/bin/env python3
"""Scan separate WGS samples with bounded resources and safe resume.

Manifest: tab-separated header sample_id, path1, and optional path2.
Input paths are relative to the manifest directory unless absolute. Each row is
one biological sample; rows are never pooled. Pair membership is supplied by the
user and cannot be inferred from FASTQ sequence content. --paired-layout fr
requires two synchronized mate files in every row and reads mode. The scanner
validates their pairing; filenames alone do not establish pair membership.

--memory-limit-mib is forwarded to each scanner process. With J concurrent jobs
and limit M, J*M is the sum of configured process limits, not measured total RSS.
The Python coordinator and other programs need additional system resources.

--flank-period-annotation on forwards the optional descriptive annotation to
the native scanner and changes resume identity. Default/explicit off forwards
no new native flag, preserving existing scanner commands and resume records.

Example:
  python3 scripts/batch.py samples.tsv --binary bin/trm \
      --output-dir results/batch --jobs 4 --threads 2 --resume

Resume requires identical input size/mtime/path, executable checksum, scan
options, and an intact previous result JSON. Input contents are not hashed.
Optionally, --input-verification-report binds every selected input to a local
schema-2 verify_external_download.py MD5 receipt and its device/inode/size/
mtime/ctime/mode. Stale or copied files require a new verification report;
batch checks receipt/stat identity without hashing the archive again.
Each output directory/shard permits one writer. Native children inherit the
lock, and execute a verified read-only binary snapshot rather than the mutable
installation path. Failed reruns retain old result bytes with an explicitly
failed current state; consumers must check state.json before using them.

Only --jobs samples are submitted at once. Full cohort summaries are refreshed
at --checkpoint-seconds intervals and at completion; each sample's durable
state is still committed immediately. --schedule largest-first uses total
input file bytes as a duration proxy without changing result order or content.
--shard-index I --shard-count N selects SHA256(sample_id) modulo N and writes
under output-dir/shard-IIIII-of-NNNNN. Only that shard's files must exist locally;
the complete manifest header and sample identifiers are always validated.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time
import uuid


try:
    from batch_common import atomic_text, positive_int, sha256
except ModuleNotFoundError as error:
    if error.name != "batch_common":
        raise
    # Direct file loaders need not place scripts/ on sys.path.
    import importlib.util
    _spec = importlib.util.spec_from_file_location("_trm_batch_common", Path(__file__).with_name("batch_common.py"))
    _common = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_common)
    atomic_text = _common.atomic_text
    positive_int = _common.positive_int
    sha256 = _common.sha256

try:
    from batch_integrity import (assert_input_verification, assert_snapshot, freeze_binary,
                                 input_verification_binding, load_input_verification, output_lock)
except ModuleNotFoundError as error:
    if error.name != "batch_integrity":
        raise
    import importlib.util
    _spec = importlib.util.spec_from_file_location("_trm_batch_integrity", Path(__file__).with_name("batch_integrity.py"))
    _integrity = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_integrity)
    assert_snapshot = _integrity.assert_snapshot
    freeze_binary = _integrity.freeze_binary
    output_lock = _integrity.output_lock
    assert_input_verification = _integrity.assert_input_verification
    input_verification_binding = _integrity.input_verification_binding
    load_input_verification = _integrity.load_input_verification


SAFE_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}\Z")
SUMMARY_FIELDS = ["sample_id", "status", "wall_seconds", "result_json", "records",
                  "bases", "candidate_count", "top_motif", "top_status", "error",
                  "input_compressed_bytes", "user_cpu_seconds", "system_cpu_seconds",
                  "max_rss_kb", "compressed_MBps", "sequence_Mbps", "metrics_json"]


def positive_seconds(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be finite and greater than zero")
    return number


def sample_shard(sample_id: str, shard_count: int) -> int:
    """Stable partition independent of manifest order and Python hash seeding."""
    return int.from_bytes(hashlib.sha256(sample_id.encode("utf-8")).digest(), "big") % shard_count


def available_cpus() -> int:
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        return max(1, os.cpu_count() or 1)


def atomic_json(path: Path, value: object) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def input_metadata(paths: list[Path]) -> list[dict]:
    result = []
    for path in paths:
        status = path.stat()
        if not path.is_file():
            raise ValueError(f"input is not a regular file: {path}")
        result.append({"path": str(path), "size": status.st_size, "mtime_ns": status.st_mtime_ns})
    return result


def load_manifest(path: Path, shard_index: int | None = None,
                  shard_count: int | None = None) -> list[dict]:
    samples = []
    identifiers: set[str] = set()
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames not in (["sample_id", "path1"], ["sample_id", "path1", "path2"]):
            raise ValueError("TSV header must be sample_id<TAB>path1[<TAB>path2]")
        for row in reader:
            line = reader.line_num
            if None in row or any(value is None for key, value in row.items() if key != "path2"):
                raise ValueError(f"manifest line {line}: incorrect number of fields")
            sample_id = row["sample_id"]
            if not SAFE_ID.fullmatch(sample_id):
                raise ValueError(f"manifest line {line}: unsafe sample_id {sample_id!r}")
            if sample_id in {"summary.json", "summary.tsv", "progress.json"}:
                raise ValueError(f"manifest line {line}: reserved output name {sample_id!r}")
            if sample_id in identifiers:
                raise ValueError(f"manifest line {line}: duplicate sample_id {sample_id!r}")
            identifiers.add(sample_id)
            if not row["path1"]:
                raise ValueError(f"manifest line {line}: path1 must not be empty")
            # Validate the complete manifest structure/IDs, but a cloud shard
            # only needs its own input files staged on the current host.
            if shard_count is not None and sample_shard(sample_id, shard_count) != shard_index:
                continue
            paths = []
            for value in (row["path1"], row.get("path2")):
                if value:
                    item = Path(value).expanduser()
                    paths.append((item if item.is_absolute() else path.parent / item).resolve(strict=True))
            if len(set(paths)) != len(paths):
                raise ValueError(f"manifest line {line}: path1 and path2 refer to the same file")
            input_metadata(paths)
            samples.append({"sample_id": sample_id, "paths": paths})
    if not identifiers:
        raise ValueError("manifest has no samples")
    return samples


def read_result(path: Path) -> dict:
    with path.open(encoding="utf-8") as source:
        result = json.load(source)
    if not isinstance(result, dict) or not isinstance(result.get("summary"), dict):
        raise ValueError("result lacks summary object")
    if not isinstance(result.get("candidates"), list):
        raise ValueError("result lacks candidates array")
    for key in ("records", "bases"):
        if not isinstance(result["summary"].get(key), int) or result["summary"][key] < 0:
            raise ValueError(f"result summary lacks nonnegative integer {key}")
    return result


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def resource_values(path: Path) -> dict:
    """GNU time can prepend a diagnostic on failure; the numeric row is last."""
    fields = path.read_text(encoding="utf-8").splitlines()[-1].split("\t")
    if len(fields) != 5:
        raise ValueError(f"invalid GNU time report: {path}")
    return {"resource_wall_seconds": float(fields[0]), "user_cpu_seconds": float(fields[1]),
            "system_cpu_seconds": float(fields[2]), "max_rss_kb": int(fields[3]),
            "resource_exit_code": int(fields[4])}


def summarize(sample_id: str, result: dict, path: Path, status: str, metrics: dict) -> dict:
    candidates = result["candidates"]
    top = candidates[0] if candidates else {}
    return {"sample_id": sample_id, "status": status, "wall_seconds": metrics.get("wall_seconds"),
            "result_json": str(path), "records": result["summary"]["records"],
            "bases": result["summary"]["bases"], "candidate_count": len(candidates),
            "top_motif": top.get("motif", ""), "top_status": top.get("status", ""), "error": "",
            **{key: metrics.get(key) for key in ("input_compressed_bytes", "user_cpu_seconds",
               "system_cpu_seconds", "max_rss_kb", "compressed_MBps", "sequence_Mbps")},
            "metrics_json": str(path.parent / "metrics.json") if (path.parent / "metrics.json").is_file() else ""}


def scan_sample(sample: dict, args: argparse.Namespace, executable: dict, options: list[str]) -> dict:
    """Protect direct callers as well as CLI cohorts against sample collisions."""
    sample_id = sample["sample_id"]
    started = time.perf_counter()
    try:
        if not SAFE_ID.fullmatch(sample_id):
            raise ValueError(f"unsafe sample_id {sample_id!r}")
        with output_lock(args.output_dir / sample_id, '.scan.lock') as sample_lock:
            return _scan_sample_locked(sample, args, executable, options, sample_lock)
    except (OSError, ValueError) as error:
        # A lock contender must not turn the active writer's state into failed.
        return {"sample_id": sample_id, "status": "failed", "error": str(error),
                "wall_seconds": time.perf_counter() - started, "result_json": "",
                "records": "", "bases": "", "candidate_count": "", "top_motif": "", "top_status": "",
                "input_compressed_bytes": None, "user_cpu_seconds": None,
                "system_cpu_seconds": None, "max_rss_kb": None, "compressed_MBps": None,
                "sequence_Mbps": None, "metrics_json": ""}


def _scan_sample_locked(sample: dict, args: argparse.Namespace, executable: dict,
                        options: list[str], sample_lock: int) -> dict:
    sample_id = sample["sample_id"]
    started = time.perf_counter()
    destination = args.output_dir / sample_id
    fingerprint = None
    metrics = {"schema_version": 1, "started_utc": utc_now(), "command": [],
               "user_cpu_seconds": None, "system_cpu_seconds": None, "max_rss_kb": None,
               "compressed_MBps": None, "sequence_Mbps": None, "exit_code": None}
    annotation_state = {"flank_period_annotation": True} if getattr(args, "flank_period_annotation", "off") == "on" else {}
    metrics.update(annotation_state)
    execution = {"schema_version": 2, "sample_id": sample_id,
                 "execution_id": getattr(args, "_execution_id", None) or uuid.uuid4().hex}
    try:
        if destination.is_symlink():
            raise ValueError(f"sample output directory must not be a symlink: {destination}")
        destination.mkdir(parents=True, exist_ok=True)
        paths = sample["paths"]
        fingerprint = {"format_version": 1, "inputs": input_metadata(paths),
                       "executable": executable, "options": options}
        verification = getattr(args, "_input_verification", None)
        if verification is not None:
            fingerprint["input_verification"] = input_verification_binding(verification, paths)
        if args.background_model:
            fingerprint["background_model"] = {"path": str(args.background_model),
                                               "sha256": sha256(args.background_model)}
        state_path = destination / "state.json"
        result_path = destination / "result.json"
        snapshot = getattr(args, "_execution_binary", None) or freeze_binary(args.output_dir, executable)
        assert_snapshot(snapshot)
        execution["execution_binary"] = snapshot
        if args.resume and state_path.is_file() and result_path.is_file():
            try:
                state = json.loads(state_path.read_text(encoding="utf-8"))
                if (state.get("status") == "complete" and state.get("fingerprint") == fingerprint
                        and state.get("result_sha256") == sha256(result_path)):
                    result = read_result(result_path)
                    from native_diagnostic_options import validate_diagnostic_output
                    validate_diagnostic_output(result, args)
                    previous_metrics = state.get("metrics") or {
                        "wall_seconds": result["summary"].get("elapsed_seconds"),
                        "input_compressed_bytes": sum(item["size"] for item in fingerprint["inputs"])}
                    assert_snapshot(snapshot)
                    assert_input_verification(verification, paths)
                    return summarize(sample_id, result, result_path, "resumed", previous_metrics)
            except (OSError, ValueError, TypeError, AttributeError):
                pass  # Incomplete/changed artifacts cause a fresh scan.
        atomic_json(state_path, {"status": "running", "fingerprint": fingerprint,
                    "retained_previous_result_is_not_current": result_path.exists(),
                    **execution, **annotation_state})
        with tempfile.TemporaryDirectory(prefix=".scan-", dir=destination) as work:
            prefix = Path(work) / "result"
            command = [snapshot["path"], "scan", *map(str, paths), "--output", str(prefix), *options]
            metrics.update({"command": command, "started_utc": utc_now(),
                            "input_compressed_bytes": sum(item["size"] for item in fingerprint["inputs"])})
            launched_command = command
            if args.record_resources:
                launched_command = [args.time_binary, "-f",
                    "wall_seconds\tuser_cpu_seconds\tsystem_cpu_seconds\tmax_rss_kb\texit_code\n%e\t%U\t%S\t%M\t%x",
                    "-o", str(destination / "resources.tsv"), *command]
            scan_started = time.perf_counter()
            with (destination / "stdout.log").open("w", encoding="utf-8") as stdout, \
                    (destination / "stderr.log").open("w", encoding="utf-8") as stderr:
                assert_snapshot(snapshot)
                assert_input_verification(verification, paths)
                completed = subprocess.run(launched_command, stdout=stdout, stderr=stderr, check=False,
                    pass_fds=tuple(getattr(args, "_runtime_lock_fds", ())) + (sample_lock,))
            metrics.update({"finished_utc": utc_now(), "wall_seconds": time.perf_counter() - scan_started,
                            "exit_code": completed.returncode})
            if args.record_resources:
                metrics.update(resource_values(destination / "resources.tsv"))
            metrics["compressed_MBps"] = metrics["input_compressed_bytes"] / metrics["wall_seconds"] / 1e6
            atomic_json(destination / "metrics.json", metrics)
            if completed.returncode:
                raise RuntimeError(f"scan exited {completed.returncode}; see {destination / 'stderr.log'}")
            generated = Path(str(prefix) + ".json")
            result = read_result(generated)
            from native_diagnostic_options import validate_diagnostic_output
            validate_diagnostic_output(result, args)
            metrics["sequence_Mbps"] = result["summary"]["bases"] / metrics["wall_seconds"] / 1e6
            metrics["scanner_wall_seconds"] = result["summary"].get("elapsed_seconds")
            atomic_json(destination / "metrics.json", metrics)
            if input_metadata(paths) != fingerprint["inputs"]:
                raise RuntimeError("input metadata changed during scan; result cannot be resumed safely")
            if args.background_model and sha256(args.background_model) != fingerprint["background_model"]["sha256"]:
                raise RuntimeError("background model changed during scan; result cannot be resumed safely")
            assert_snapshot(snapshot)
            assert_input_verification(verification, paths)
            for artifact in sorted(Path(work).iterdir()):
                if artifact.is_file() and artifact.name.startswith("result."):
                    os.replace(artifact, destination / artifact.name)
            atomic_json(state_path, {"status": "complete", "fingerprint": fingerprint,
                                    "result_sha256": sha256(result_path), "command": command,
                                    "metrics": metrics, "retained_previous_result_is_not_current": False,
                                    **execution, **annotation_state})
        return summarize(sample_id, result, result_path, "complete", metrics)
    except (OSError, ValueError, TypeError, KeyError, RuntimeError) as error:
        metrics.setdefault("wall_seconds", time.perf_counter() - started)
        metrics.setdefault("finished_utc", utc_now())
        if destination.is_dir() and not destination.is_symlink():
            try:
                atomic_json(destination / "metrics.json", metrics)
                atomic_json(destination / "state.json", {"status": "failed", "fingerprint": fingerprint,
                            "error": str(error), "metrics": metrics,
                            "retained_previous_result_is_not_current": (destination / "result.json").exists(),
                            **execution, **annotation_state})
            except OSError:
                pass
        return {"sample_id": sample_id, "status": "failed", "wall_seconds": metrics["wall_seconds"],
                "result_json": "", "records": "", "bases": "", "candidate_count": "",
                "top_motif": "", "top_status": "", "error": str(error),
                **{key: metrics.get(key) for key in ("input_compressed_bytes", "user_cpu_seconds",
                   "system_cpu_seconds", "max_rss_kb", "compressed_MBps", "sequence_Mbps")},
                "metrics_json": str(destination / "metrics.json") if (destination / "metrics.json").is_file() else ""}


def main() -> int:
    from runtime_binary import resolve_binary
    from runner_resources import add_arguments as add_resource_arguments, preflight
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--binary", help="Explicit scanner override; default: TRM_BINARY or verified bin/trm")
    parser.add_argument("--output-dir", type=Path, default=Path("results/batch"))
    parser.add_argument("--jobs", type=positive_int, default=1)
    parser.add_argument("--threads", type=positive_int, help="threads per scan; defaults to cpu-budget // jobs")
    parser.add_argument("--cpu-budget", type=positive_int, default=available_cpus(),
                        help="maximum concurrent scanner workers, excluding helpers (defaults to CPU affinity count)")
    add_resource_arguments(parser)
    parser.add_argument("--memory-limit-mib", type=positive_int,
                        help="memory limit per scanner process in MiB; not an aggregate RSS budget")
    parser.add_argument("--paired-layout", choices=("fr",),
                        help="validate and annotate inward-facing paired reads; requires two files per row")
    parser.add_argument("--flank-period-annotation", choices=("on", "off"), default="off",
                        help="optional native flank-period diagnostics; on requires a supporting scanner (default: off)")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--input-verification-report", type=Path,
                        help="optional completed schema-2 local archive verification JSON; binds selected inputs to prior MD5 and current stat identity without rehashing")
    parser.add_argument("--record-resources", action="store_true",
                        help="record GNU time CPU seconds and peak RSS for every scan")
    parser.add_argument("--checkpoint-seconds", type=positive_seconds, default=10.0,
                        help="interval between full summary checkpoints; final checkpoint is mandatory")
    parser.add_argument("--schedule", choices=("manifest", "largest-first"), default="manifest",
                        help="submission order; final summaries always retain manifest order")
    parser.add_argument("--shard-index", type=int,
                        help="zero-based sample-ID hash partition; requires --shard-count")
    parser.add_argument("--shard-count", type=positive_int,
                        help="number of disjoint sample-ID hash partitions; requires --shard-index")
    parser.add_argument("--mode", choices=("reads", "assembly"), default="reads")
    parser.add_argument("--min-k", type=positive_int, default=3)
    parser.add_argument("--max-k", type=positive_int, default=32)
    parser.add_argument("--min-tract", type=positive_int, default=500)
    parser.add_argument("--min-reads", type=positive_int, default=2)
    parser.add_argument("--prior-policy", choices=("canonical", "telomere", "discovery"), default="canonical")
    parser.add_argument("--prior-max-edits", type=int, default=2)
    parser.add_argument("--prior-max-distance", type=float, default=0.34)
    parser.add_argument("--background-model", type=Path)
    parser.add_argument("--support-policy", choices=("filter", "annotate"), default="annotate")
    parser.add_argument("--support-min-enrichment", type=float, default=2.0)
    from native_diagnostic_options import add_arguments, native_diagnostic_options
    add_arguments(parser)
    args = parser.parse_args()
    if (args.shard_index is None) != (args.shard_count is None):
        parser.error("--shard-index and --shard-count must be supplied together")
    if args.shard_index is not None and not 0 <= args.shard_index < args.shard_count:
        parser.error("--shard-index must be >= 0 and less than --shard-count")
    if args.threads is None:
        args.threads = max(1, args.cpu_budget // args.jobs)
    if args.jobs * args.threads > args.cpu_budget:
        parser.error(f"jobs * threads ({args.jobs * args.threads}) exceeds cpu-budget ({args.cpu_budget})")
    if args.min_k > args.max_k:
        parser.error("min-k must not exceed max-k")
    if args.paired_layout and args.mode != "reads":
        parser.error("--paired-layout requires --mode reads")
    try:
        native_diagnostic_options(args)
        args.time_binary = "/usr/bin/time"
        samples = load_manifest(args.manifest.resolve(strict=True), args.shard_index, args.shard_count)
        if args.paired_layout:
            for sample in samples:
                if len(sample["paths"]) != 2:
                    raise ValueError(
                        f"sample {sample['sample_id']}: --paired-layout requires exactly two input files")
        if preflight(parser, args, selected_samples=len(samples), timed=args.record_resources,
                     cpu_affinity=available_cpus()):
            return 0
        if args.record_resources:
            probe = subprocess.run([args.time_binary, "--version"], capture_output=True, text=True, check=False)
            if probe.returncode or "GNU" not in probe.stdout + probe.stderr:
                raise ValueError("--record-resources requires GNU /usr/bin/time")
        executable_path = resolve_binary(args.binary)
        executable = {"path": str(Path(executable_path).resolve()), "sha256": sha256(Path(executable_path))}
        if args.background_model:
            args.background_model = args.background_model.resolve(strict=True)
        args._input_verification = (load_input_verification(args.input_verification_report,
            [path for sample in samples for path in sample["paths"]])
            if args.input_verification_report else None)
        args.output_dir = args.output_dir.resolve()
        if args.shard_count is not None:
            args.output_dir /= f"shard-{args.shard_index:05d}-of-{args.shard_count:05d}"
        protected = [args.manifest.resolve(), *[path for sample in samples for path in sample["paths"]]]
        if args.input_verification_report:
            protected.append(args.input_verification_report.absolute())
        for source in protected:
            if source == args.output_dir or args.output_dir in source.parents:
                raise ValueError(f"Batch output must be disjoint from manifest, verification report and raw inputs: {source}")
    except (OSError, ValueError) as error:
        parser.error(str(error))
    try:
        with output_lock(args.output_dir) as lock:
            args._runtime_lock_fds = (lock,)
            args._execution_id = uuid.uuid4().hex
            assert_input_verification(args._input_verification)
            args._execution_binary = freeze_binary(args.output_dir, executable)
            return execute_batch(args, samples, executable)
    except (OSError, ValueError) as error:
        parser.error(str(error))


def execute_batch(args: argparse.Namespace, samples: list[dict], executable: dict) -> int:
    """Run while the caller holds the output namespace lock and binary pin."""
    options = ["--threads", str(args.threads), "--mode", args.mode,
               "--min-k", str(args.min_k), "--max-k", str(args.max_k),
               "--min-tract", str(args.min_tract), "--min-reads", str(args.min_reads),
               "--prior-policy", args.prior_policy,
               "--prior-max-edits", str(args.prior_max_edits),
               "--prior-max-distance", str(args.prior_max_distance)]
    from native_diagnostic_options import native_diagnostic_options
    options.extend(native_diagnostic_options(args))
    if args.memory_limit_mib is not None:
        options.extend(("--memory-limit-mib", str(args.memory_limit_mib)))
    if args.paired_layout:
        options.extend(("--paired-layout", args.paired_layout))
    if args.flank_period_annotation == "on":
        options.extend(("--flank-period-annotation", "on"))
    if args.background_model:
        options.extend(("--background-model", str(args.background_model),
                        "--support-policy", args.support_policy,
                        "--support-min-enrichment", str(args.support_min_enrichment)))
    started = time.perf_counter()
    started_utc = utc_now()
    results = []
    order = {sample["sample_id"]: index for index, sample in enumerate(samples)}
    input_sizes = {sample["sample_id"]: sum(item["size"] for item in input_metadata(sample["paths"]))
                   for sample in samples}
    total_input_bytes = sum(input_sizes.values())
    scheduled = sorted(samples, key=lambda sample: -input_sizes[sample["sample_id"]]) \
        if args.schedule == "largest-first" else samples
    submitted = 0
    checkpoint_count = 0
    max_inflight = 0
    execution_error = None

    def checkpoint(pending: dict, finished: bool = False) -> dict:
        nonlocal checkpoint_count
        checkpoint_count += 1
        ordered = sorted(results, key=lambda item: order[item["sample_id"]])
        succeeded = [item for item in ordered if item["status"] != "failed"]
        scanned = [item for item in succeeded if item["status"] == "complete"]
        elapsed = time.perf_counter() - started
        completed_bytes = sum(item.get("input_compressed_bytes") or 0 for item in succeeded)
        scanned_bytes = sum(item.get("input_compressed_bytes") or 0 for item in scanned)
        scanned_bases = sum(item["bases"] for item in scanned)
        report = {"manifest": str(args.manifest.resolve()), "jobs": args.jobs, "threads": args.threads,
                  "execution_id": args._execution_id, "execution_binary": args._execution_binary,
                  "execution_error": execution_error,
                  "binary_identity_policy": "Verified isolated snapshot; full SHA256 at start/end and stat identity before/after each scan",
                  "input_identity_policy": ("User-trusted local MD5 verification receipt plus current path/device/inode/size/mtime/ctime/mode; no raw rehash in batch"
                      if args._input_verification is not None else
                      "Path, size and mtime only; raw input content checksums are not verified by batch resume"),
                  "input_verification_report": args._input_verification['report'] if args._input_verification is not None else None,
                  "cpu_budget": args.cpu_budget, "record_resources": args.record_resources,
                  "resource_plan": getattr(args, "_resource_plan", None),
                  "schedule": args.schedule, "checkpoint_seconds": args.checkpoint_seconds,
                  "checkpoint_count": checkpoint_count, "max_inflight": max_inflight,
                  "shard_index": args.shard_index, "shard_count": args.shard_count,
                  "shard_algorithm": "sha256_utf8_sample_id_big_endian_modulo_count" if args.shard_count else None,
                  "paired_layout": args.paired_layout,
                  "flank_period_annotation": args.flank_period_annotation == "on",
                  "memory_limit_mib_per_scan": args.memory_limit_mib,
                  "configured_jobs_times_memory_limit_mib": (
                      args.jobs * args.memory_limit_mib if args.memory_limit_mib is not None else None),
                  "memory_limit_note": (
                      "Limits apply per scanner process; their sum is not measured aggregate RSS "
                      "and excludes the batch coordinator and other system processes."),
                  "started_utc": started_utc, "finished_utc": utc_now() if finished else None,
                  "wall_seconds": elapsed, "finished": finished,
                  "all_succeeded": finished and len(succeeded) == len(samples) and execution_error is None,
                  "total_samples": len(samples), "completed_samples": len(succeeded),
                  "failed_samples": len(ordered) - len(succeeded),
                  "resumed_samples": len(succeeded) - len(scanned),
                  "total_input_compressed_bytes": total_input_bytes,
                  "completed_input_compressed_bytes": completed_bytes,
                  "scanned_input_compressed_bytes": scanned_bytes,
                  "total_records": sum(item["records"] for item in succeeded),
                  "total_bases": sum(item["bases"] for item in succeeded),
                  "scanned_bases": scanned_bases,
                  "compressed_MBps": scanned_bytes / elapsed / 1e6,
                  "sequence_Mbps": scanned_bases / elapsed / 1e6,
                  "sum_scan_wall_seconds": sum(item.get("wall_seconds") or 0 for item in scanned),
                  "sum_user_cpu_seconds": sum(item.get("user_cpu_seconds") or 0 for item in scanned),
                  "sum_system_cpu_seconds": sum(item.get("system_cpu_seconds") or 0 for item in scanned),
                  "max_sample_rss_kb": max((item.get("max_rss_kb") or 0 for item in scanned), default=0),
                  "samples": ordered}
        atomic_json(args.output_dir / "summary.json", report)
        stream = io.StringIO(newline="")
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(ordered)
        atomic_text(args.output_dir / "summary.tsv", stream.getvalue())
        running = sorted(sample_id for future, sample_id in pending.items() if future.running())
        atomic_json(args.output_dir / "progress.json", {
            "updated_utc": utc_now(), "started_utc": started_utc,
            "execution_id": args._execution_id, "execution_error": execution_error,
            "resource_plan": getattr(args, "_resource_plan", None),
            "finished_utc": report["finished_utc"], "finished": finished,
            "total": len(samples), "completed": len(succeeded), "failed": report["failed_samples"],
            "resumed": report["resumed_samples"], "running": running,
            "flank_period_annotation": args.flank_period_annotation == "on",
            "running_count": len(running),
            "pending": len(samples) - len(results) - len(running),
            "submitted": submitted, "inflight": len(pending),
            "checkpoint_count": checkpoint_count,
            "elapsed_seconds": elapsed, "total_input_compressed_bytes": total_input_bytes,
            "completed_input_compressed_bytes": completed_bytes})
        return report

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = {}

        def replenish() -> None:
            nonlocal submitted, max_inflight
            while submitted < len(scheduled) and len(pending) < args.jobs:
                sample = scheduled[submitted]
                pending[pool.submit(scan_sample, sample, args, executable, options)] = sample["sample_id"]
                submitted += 1
            max_inflight = max(max_inflight, len(pending))

        replenish()
        checkpoint(pending)
        next_checkpoint = time.perf_counter() + args.checkpoint_seconds
        while pending:
            done, _ = wait(pending, timeout=min(10.0, max(0, next_checkpoint - time.perf_counter())),
                           return_when=FIRST_COMPLETED)
            for future in done:
                pending.pop(future)
                result = future.result()
                results.append(result)
                print(f"{result['status'].upper()} {result['sample_id']}" +
                      (f": {result['error']}" if result["error"] else ""), flush=True)
            replenish()
            if time.perf_counter() >= next_checkpoint:
                checkpoint(pending)
                next_checkpoint = time.perf_counter() + args.checkpoint_seconds
    try:
        assert_snapshot(args._execution_binary, hash_contents=True)
        assert_input_verification(args._input_verification, hash_report=True)
    except (OSError, ValueError) as error:
        execution_error = str(error)
    report = checkpoint({}, finished=True)
    print(f"Summary: {args.output_dir / 'summary.tsv'}")
    return 0 if report["all_succeeded"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
