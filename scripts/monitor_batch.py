#!/usr/bin/env python3
"""Estimate live compressed-input progress from /proc; never read FASTQ payloads.

Run once in the batch supervisor's PID namespace. File offsets include buffered
read-ahead, so rates and ETA describe estimated bytes read, not analyzed bases.
An inaccessible namespace produces unknown estimates rather than zero activity.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import tempfile


def read_json(path: Path, default):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return default


def observe_process(proc: Path, inputs: list[dict]) -> dict:
    offsets = {}
    for fd in (proc / "fd").iterdir():
        try:
            target = os.readlink(fd)
            for line in (proc / "fdinfo" / fd.name).read_text().splitlines():
                if line.startswith("pos:"):
                    offsets[target] = max(offsets.get(target, 0), int(line.split()[1]))
        except (OSError, ValueError):
            continue
    current = [i for i, item in enumerate(inputs) if item["path"] in offsets]
    read_bytes = None
    if current:
        index = max(current)  # The scanner consumes input paths sequentially.
        read_bytes = sum(item["bytes"] for item in inputs[:index])
        read_bytes += min(inputs[index]["bytes"], offsets[inputs[index]["path"]])
    rss = next((int(line.split()[1]) for line in (proc / "status").read_text().splitlines()
                if line.startswith("VmRSS:")), None)
    stat = (proc / "stat").read_text().rsplit(")", 1)[1].split()
    return {"pid": int(proc.name), "sample_id": inputs[0]["sample_id"],
            "estimated_input_bytes_read": read_bytes, "rss_kb": rss,
            "cpu_seconds": (int(stat[11]) + int(stat[12])) / os.sysconf("SC_CLK_TCK")}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    root = args.output_dir.resolve()
    inventory = {item["path"]: item for item in read_json(root / "input_inventory.json", [])}
    progress = read_json(root / "progress.json", {})
    previous = read_json(root / "live_progress.json", {})
    timestamp = datetime.now(timezone.utc)
    observed = []
    for proc in Path("/proc").iterdir():
        if not proc.name.isdecimal():
            continue
        try:
            command = (proc / "cmdline").read_bytes().split(b"\0")
            if not command or os.fsdecode(command[0]) != str(root / "code" / "trm"):
                continue
            inputs = [inventory[os.fsdecode(arg)] for arg in command if os.fsdecode(arg) in inventory]
            if inputs:
                observed.append(observe_process(proc, inputs))
        except (OSError, ValueError, IndexError):
            continue  # Processes may exit between proc reads.
    expected = set(progress.get("running", []))
    seen = {item["sample_id"] for item in observed}
    visible = bool(progress) and seen == expected and all(
        item["estimated_input_bytes_read"] is not None for item in observed)
    completed_bytes = progress.get("completed_input_compressed_bytes", 0)
    total_bytes = progress.get("total_input_compressed_bytes", 0)
    estimated = completed_bytes + sum(item["estimated_input_bytes_read"] or 0 for item in observed) if visible else None
    started = datetime.fromisoformat(progress.get("started_utc", timestamp.isoformat()))
    reference = previous.get("last_observed")
    rate = None
    if estimated is not None and reference:
        seconds = (timestamp - datetime.fromisoformat(reference["utc"])).total_seconds()
        delta = estimated - reference["estimated_input_bytes_read"]
        if seconds > 0 and delta > 0:
            rate = delta / seconds
    finished = bool(progress.get("finished"))
    if finished:
        rate = None  # Failed inputs do not imply that a finished batch has work queued.
    elapsed = progress.get("elapsed_seconds") if finished else (timestamp - started).total_seconds()
    record = {"utc": timestamp.isoformat(), "batch_elapsed_seconds": elapsed,
              "completed": progress.get("completed"), "running": progress.get("running_count"),
              "failed": progress.get("failed"), "finished": progress.get("finished"),
              "observed_worker_count": len(observed), "proc_visibility_complete": visible,
              "completed_input_compressed_bytes": completed_bytes, "total_input_compressed_bytes": total_bytes,
              "estimated_input_bytes_read": estimated,
              "estimated_fraction_read": estimated / total_bytes if estimated is not None and total_bytes else None,
              "estimated_rolling_bytes_per_second": rate,
              "estimated_eta_seconds": max(0, total_bytes - estimated) / rate if rate and not finished else None,
              "active_rss_kb": sum(item["rss_kb"] or 0 for item in observed) if visible else None,
              "active_cpu_seconds": sum(item["cpu_seconds"] for item in observed) if visible else None,
              "workers": observed,
              "note": "Estimates use compressed file offsets, include read-ahead, and do not measure analyzed bases. "
                      "Unknown values mean /proc visibility or sample transition prevented a complete observation."}
    if finished:
        record["note"] += " Batch has finished; ETA is not applicable. Failed inputs can leave validated input coverage below 100%."
    if estimated is not None:
        reference = {"utc": record["utc"], "estimated_input_bytes_read": estimated}
    record["last_observed"] = reference
    descriptor, temporary = tempfile.mkstemp(prefix=".live-progress-", dir=root)
    with os.fdopen(descriptor, "w") as output:
        json.dump(record, output, indent=2)
        output.write("\n")
    os.replace(temporary, root / "live_progress.json")
    with (root / "telemetry.jsonl").open("a") as output:
        output.write(json.dumps(record) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
