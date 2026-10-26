#!/usr/bin/env python3
"""Evaluate a fixed background model on held-out reads and separate positives.

This is a smoke evaluation, not a cross-species specificity or FDR estimate.
The supplied held-out negative must not be an input used to build the model.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def raw_rows(report):
    return {row["motif"]: {k: v for k, v in row.items()
                          if k not in ("telomere_prior", "support", "candidate_eligible")}
            for key in ("candidates", "rejected_candidates") for row in report.get(key, [])}


def main():
    from runtime_binary import resolve_binary
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", help="Explicit scanner override; default: TRM_BINARY or verified bin/trm")
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--heldout-result", type=Path, required=True)
    parser.add_argument("--positive-result", type=Path, default=Path("results/real_v020_verified/arabidopsis.json"))
    parser.add_argument("--assembly-result", type=Path, default=Path("results/real_v020_verified/ce11_chrI.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/support_v030/validation"))
    args = parser.parse_args()
    try:
        args.binary = resolve_binary(args.binary)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error("output directory must be new or empty")
    model = json.loads(args.model.read_text())
    sources = [args.heldout_result, args.positive_result, args.assembly_result]
    calibration_hashes = {c["source_result_sha256"] for c in model["controls"]}
    if any(sha256(path) in calibration_hashes for path in sources):
        parser.error("evaluation result is also a calibration source")
    calibration_inputs = {str(Path(item["path"]).resolve()) for c in model["controls"]
                          for item in c.get("inputs", [])}
    for path in sources:
        source = json.loads(path.read_text())
        if any(str(Path(item["path"]).resolve()) in calibration_inputs for item in source.get("inputs", [])):
            parser.error("evaluation and calibration input paths overlap")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = {"binary_sha256": sha256(args.binary), "model_sha256": sha256(args.model),
              "model": str(args.model.resolve()), "min_enrichment": 2.0,
              "scope": "Drosophila paired-read hash holdout within the same library; separate Arabidopsis positive; assembly mode safeguard. Not independent negative taxa or population specificity.",
              "cases": []}
    for name, path in zip(("drosophila_heldout", "arabidopsis", "ce11_assembly"), sources):
        source_hash = sha256(path)
        source = json.loads(path.read_text())
        prefix = args.output_dir / name
        command = [str(args.binary.resolve()), "reclassify", str(path), "--background-model", str(args.model),
                   "--support-min-enrichment", "2", "--output", str(prefix)]
        started = time.perf_counter()
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        elapsed = time.perf_counter() - started
        Path(str(prefix) + ".stderr.log").write_text(completed.stderr)
        if completed.returncode:
            raise RuntimeError(completed.stderr)
        result = json.loads(Path(str(prefix) + ".json").read_text())
        candidates = result["candidates"]
        checks = {"original_file_preserved": source_hash == sha256(path),
                  "all_measured_evidence_preserved": raw_rows(source) == raw_rows(result),
                  "scan_timing_preserved": source["summary"]["elapsed_seconds"] == result["summary"]["elapsed_seconds"],
                  "all_primary_candidates_eligible": all(c["candidate_eligible"] for c in candidates)}
        if name == "drosophila_heldout":
            checks["no_eligible_candidates_in_this_holdout"] = len(candidates) == 0
        elif name == "arabidopsis":
            checks["expected_TTTAGGG_top_retained"] = bool(candidates) and candidates[0]["motif"] == "AAACCCT"
            checks["weak_three_record_TTAGGC_removed"] = all(c["motif"] != "AAGCCT" for c in candidates)
            checks["normalization_keeps_observed_read_count"] = bool(candidates) and candidates[0]["reads"] == 157
        else:
            checks["assembly_prior_candidates_retained"] = [c["motif"] for c in candidates] == [c["motif"] for c in source["candidates"]]
            checks["assembly_background_not_applied"] = all(c["support"]["status"] == "not_applicable" for c in candidates)
        entry = {"name": name, "source_result": str(path.resolve()), "source_sha256": source_hash,
                 "command": command, "wall_seconds": elapsed, "before_candidates": len(source["candidates"]),
                 "after_candidates": len(candidates), "top_motif": candidates[0]["motif"] if candidates else None,
                 "checks": checks, "passed": all(checks.values()), "summary": result["summary"]}
        report["cases"].append(entry)
        print(f"{'PASS' if entry['passed'] else 'FAIL'} {name}: {entry['before_candidates']} -> {len(candidates)} candidates, {elapsed:.3f}s")
    report["passed"] = all(c["passed"] for c in report["cases"])
    (args.output_dir / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
