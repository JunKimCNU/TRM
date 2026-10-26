#!/usr/bin/env python3
"""Run TREW and optionally tidk; retain their outputs without calibrating scores.

Example:
  python scripts/compare_tools.py --r1 R1.fastq.gz --r2 R2.fastq.gz \
      --trew ./vendor/trew --threads 8 --output results/comparison
  python scripts/compare_tools.py --r1 R1.fastq.gz --r2 R2.fastq.gz \
      --tidk tidk --fasta reads.fasta --output results/comparison

tidk needs FASTA (preferably long reads or assemblies). Its run is a separate
input/algorithm comparison, not automatically a fair short-read speed benchmark.
"""

import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time


def atomic_json(path, document):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                         prefix=path.name + ".", suffix=".part", delete=False) as out:
            temporary = Path(out.name)
            json.dump(document, out, indent=2, ensure_ascii=False)
            out.write("\n")
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def resolve_executable(name):
    executable = shutil.which(name)
    if executable is None:
        raise ValueError("Executable not found or not executable: " + name)
    return str(Path(executable).resolve())


def input_metadata(path, hash_inputs):
    path = path.resolve(strict=True)
    if not path.is_file():
        raise ValueError("Input is not a regular file: " + str(path))
    stat = path.stat()
    result = {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if hash_inputs:
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        result["sha256"] = digest.hexdigest()
    return result


def run_capture(command, output, label, timeout):
    stdout_path = output / (label + ".stdout.txt")
    stderr_path = output / (label + ".stderr.txt")
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    begin = time.perf_counter()
    result = {"command": command, "command_display": shlex.join(command), "started_at_utc": started,
              "stdout_path": str(stdout_path), "stderr_path": str(stderr_path), "timeout_seconds": timeout,
              "exit_code": None, "timed_out": False}
    print(result["command_display"], file=sys.stderr)
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        try:
            process = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=timeout, check=False)
            result["exit_code"] = process.returncode
        except subprocess.TimeoutExpired:
            result["timed_out"] = True
        except OSError as error:
            result["error"] = str(error)
    result["elapsed_seconds"] = time.perf_counter() - begin
    return result


def captured_text(run):
    return "\n".join(Path(run[key]).read_text(encoding="utf-8", errors="replace")
                     for key in ("stdout_path", "stderr_path"))


def require_flags(help_run, flags, tool):
    if help_run["exit_code"] != 0:
        raise ValueError(tool + " --help failed; inspect captured help output")
    help_text = captured_text(help_run)
    missing = [flag for flag in flags if not re.search(re.escape(flag) + r"(?![\w-])", help_text)]
    if missing:
        raise ValueError(tool + " installed CLI lacks expected options: " + ", ".join(missing))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--r1", type=Path, required=True, help="Illumina paired-end R1 FASTQ[.gz]")
    parser.add_argument("--r2", type=Path, required=True, help="Illumina paired-end R2 FASTQ[.gz]")
    parser.add_argument("--trew", default="trew", help="TREW executable (default: trew)")
    parser.add_argument("--tidk", help="Optional tidk executable; requires --fasta")
    parser.add_argument("--fasta", type=Path, help="FASTA for optional tidk explore; no automatic read conversion")
    parser.add_argument("--k-min", type=int, default=3)
    parser.add_argument("--k-max", type=int, default=32)
    parser.add_argument("--threads", type=int, default=1, help="TREW threads; tidk explore does not advertise this option")
    parser.add_argument("--tidk-threshold", type=int, default=100, help="tidk sequential-repeat threshold, not a bp tract length")
    parser.add_argument("--tidk-distance", type=float, default=0.5, help="Fraction searched from each end; 0.5 scans entire FASTA records")
    parser.add_argument("--timeout", type=float, default=3600, help="Maximum seconds per analysis (default: 3600)")
    parser.add_argument("--output", type=Path, default=Path("results/comparison"), help="New or empty output directory")
    parser.add_argument("--hash-inputs", action="store_true", help="Record SHA256 by reading each input before timing")
    args = parser.parse_args()
    if not 3 <= args.k_min <= args.k_max <= 64:
        parser.error("TREW requires 3 <= --k-min <= --k-max <= 64")
    if args.threads < 1 or args.tidk_threshold < 1 or not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("Threads, tidk threshold and timeout must be positive")
    if not math.isfinite(args.tidk_distance) or not 0 <= args.tidk_distance <= 0.5:
        parser.error("--tidk-distance must be between 0 and 0.5")
    if bool(args.tidk) != bool(args.fasta):
        parser.error("--tidk and --fasta must be supplied together")
    report = {"schema_version": 1, "complete": False, "runs": [], "inputs": {},
              "interpretation": [
                  "TREW H/L outputs and Putative_TRM scores are tool-specific heuristics, not calibrated probabilities or confidence scores.",
                  "A short-read candidate cannot establish a contiguous >=500 bp tract or terminal chromosome localization.",
                  "No candidate is not proof of motif absence or ALT.",
                  "tidk consumes the separately supplied FASTA; differences in input or search thresholds preclude direct speed/accuracy claims.",
                  "Reported elapsed time is wall-clock subprocess time, excluding CLI probes and optional input hashing; cache state and system load are uncontrolled."
              ]}
    try:
        trew = resolve_executable(args.trew)
        tidk = resolve_executable(args.tidk) if args.tidk else None
        report["inputs"]["r1"] = input_metadata(args.r1, args.hash_inputs)
        report["inputs"]["r2"] = input_metadata(args.r2, args.hash_inputs)
        if args.r1.resolve() == args.r2.resolve():
            raise ValueError("R1 and R2 must be distinct files")
        if args.fasta:
            report["inputs"]["fasta"] = input_metadata(args.fasta, args.hash_inputs)
        output = args.output.resolve()
        if output.exists() and any(output.iterdir()):
            raise ValueError("Output directory must be empty to preserve earlier comparison records: " + str(output))
        output.mkdir(parents=True, exist_ok=True)
    except (OSError, ValueError) as error:
        print("Error: " + str(error), file=sys.stderr)
        return 1

    report_path = output / "comparison.json"
    atomic_json(report_path, report)

    def record(command, label, timeout=30):
        run = run_capture(command, output, label, timeout)
        run["label"] = label
        report["runs"].append(run)
        atomic_json(report_path, report)
        return run

    try:
        trew_help = record([trew, "short", "--help"], "trew.help")
        require_flags(trew_help, ["--paired_end", "--fq1", "--fq2", "-t"], "TREW short")
        version_run = record([trew, "--version"], "trew.version")
        report["trew_version"] = captured_text(version_run).strip() if version_run["exit_code"] == 0 else None
        trew_command = [trew, "short", str(args.k_min), str(args.k_max), "--paired_end",
                        "--fq1", report["inputs"]["r1"]["path"], "--fq2", report["inputs"]["r2"]["path"],
                        "-t", str(args.threads)]
        analysis = record(trew_command, "trew.analysis", args.timeout)
        success = analysis["exit_code"] == 0
        if tidk:
            tidk_help = record([tidk, "explore", "--help"], "tidk.help")
            require_flags(tidk_help, ["--minimum", "--maximum", "--distance", "--threshold"], "tidk explore")
            version_run = record([tidk, "--version"], "tidk.version")
            report["tidk_version"] = captured_text(version_run).strip() if version_run["exit_code"] == 0 else None
            tidk_command = [tidk, "explore", "--minimum", str(args.k_min), "--maximum", str(args.k_max),
                            "--distance", str(args.tidk_distance), "--threshold", str(args.tidk_threshold),
                            report["inputs"]["fasta"]["path"]]
            analysis = record(tidk_command, "tidk.analysis", args.timeout)
            success = success and analysis["exit_code"] == 0
        report["complete"] = success
        atomic_json(report_path, report)
    except (OSError, ValueError) as error:
        report["error"] = str(error)
        atomic_json(report_path, report)
        print("Error: " + str(error), file=sys.stderr)
        return 1
    print(str(report_path))
    return 0 if report["complete"] else 1


if __name__ == "__main__":
    sys.exit(main())
