#!/usr/bin/env python3
"""Build an empirical background profile from explicitly selected negative results.

This command copies measured counts; it does not fit thresholds or read positive
or held-out results. Runtime defaults z=1.96 and minimum enrichment=2 are scoring
heuristics, not calibrated confidence levels, significance tests, or FDR control.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile


CONFIG_KEYS = {"mode", "min_k", "max_k", "seed_bases", "min_copies", "min_identity",
               "min_tract", "end_window", "min_reads", "max_candidates"}
DETECTOR_KEYS = sorted(CONFIG_KEYS - {"min_reads", "max_candidates"})
ANNOTATIONS = {"telomere_prior", "support", "candidate_eligible"}
U64_MAX = (1 << 64) - 1
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def integer(value, label, minimum=0, maximum=U64_MAX):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(f"{label} must be an integer in {minimum}..{maximum}")
    return value


def unique_object(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"Duplicate JSON key: {key}")
        value[key] = item
    return value


def invalid_constant(value):
    raise ValueError(f"Non-finite JSON number: {value}")


def canonical(motif):
    reverse = motif.translate(COMPLEMENT)[::-1]
    return min(sequence[offset:] + sequence[:offset]
               for sequence in (motif, reverse) for offset in range(len(motif)))


def validate_config(config):
    if not isinstance(config, dict) or set(config) != CONFIG_KEYS:
        raise ValueError("Unknown or incomplete detector config; require exactly "
                         + ", ".join(sorted(CONFIG_KEYS)))
    if config["mode"] != "reads":
        raise ValueError("Negative controls must use reads mode")
    minimum = integer(config["min_k"], "min_k", 1, 64)
    integer(config["max_k"], "max_k", minimum, 64)
    integer(config["seed_bases"], "seed_bases", 1, 100_000_000)
    integer(config["min_copies"], "min_copies", 2, 1_000_000)
    for key in ("min_reads", "min_tract", "max_candidates"):
        integer(config[key], key, 1)
    integer(config["end_window"], "end_window")
    identity = config["min_identity"]
    if (type(identity) not in (float, int) or not math.isfinite(identity)
            or not 0.5 <= identity <= 1.0):
        raise ValueError("min_identity must be finite and between 0.5 and 1.0")


def load_control(path):
    path = Path(path).resolve(strict=True)
    contents = path.read_bytes()
    report = json.loads(contents, object_pairs_hook=unique_object,
                        parse_constant=invalid_constant)
    if not isinstance(report, dict) or report.get("program") != "trm":
        raise ValueError("Control must be a saved trm scan result object")
    if type(report.get("schema_version")) is not int or report["schema_version"] not in (1, 2, 3):
        raise ValueError("Unsupported scan-result schema_version")
    config = report.get("config")
    validate_config(config)
    summary = report.get("summary")
    if not isinstance(summary, dict):
        raise ValueError("Control lacks a complete scan summary")
    for obj in (report, summary):
        if obj.get("finished") is False or obj.get("complete") is False:
            raise ValueError("Partial scan results cannot be negative controls")
        if "status" in obj and obj["status"] not in ("complete", "completed", "success"):
            raise ValueError("Incomplete or failed scan status in negative control")
    if summary.get("sampled_prefix") is not False:
        raise ValueError("Negative control must have sampled_prefix=false")
    records = integer(summary.get("records"), "summary.records", 1)
    bases = integer(summary.get("bases"), "summary.bases", 1)
    longest = integer(summary.get("max_record_bases"), "summary.max_record_bases", 1, bases)
    if bases > records * longest:
        raise ValueError("summary.bases exceeds records * max_record_bases")
    before = integer(summary.get("candidate_families_before_filter"),
                     "summary.candidate_families_before_filter")
    if before >= config["max_candidates"]:
        raise ValueError("Control reaches candidate cap; completeness is uncertain")
    inputs = report.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        raise ValueError("Control requires nonempty input provenance")
    for item in inputs:
        if not isinstance(item, dict) or not isinstance(item.get("path"), str) or not item["path"]:
            raise ValueError("Each input requires a nonempty path")
        if "bytes" in item and item["bytes"] is not None:
            integer(item["bytes"], "input.bytes")
    candidates = report.get("candidates")
    rejected = report.get("rejected_candidates", [])
    if not isinstance(candidates, list) or not isinstance(rejected, list):
        raise ValueError("Control requires candidate and rejected-candidate arrays")
    reported = integer(summary.get("reported_candidates"), "summary.reported_candidates")
    if reported != len(candidates):
        raise ValueError("reported_candidates differs from candidates array length")
    evidence = {}
    motifs = {}
    for row in candidates + rejected:
        if not isinstance(row, dict):
            raise ValueError("Each candidate must be an object")
        motif = row.get("motif")
        if (not isinstance(motif, str) or not motif or len(motif) > 64
                or any(base not in "ACGT" for base in motif)):
            raise ValueError("Candidate motif must contain 1..64 uppercase DNA bases")
        if motif != canonical(motif):
            raise ValueError(f"Noncanonical candidate motif: {motif}")
        if any(len(motif) % size == 0 and motif == motif[:size] * (len(motif) // size)
               for size in range(1, len(motif))):
            raise ValueError(f"Candidate motif is not primitive: {motif}")
        period = integer(row.get("period"), "candidate.period", config["min_k"], config["max_k"])
        if period != len(motif):
            raise ValueError(f"Candidate period does not match motif length: {motif}")
        reads = integer(row.get("reads"), "candidate.reads", config["min_reads"], records)
        repeat_bases = integer(row.get("repeat_bases"), "candidate.repeat_bases", reads, bases)
        if repeat_bases > reads * longest:
            raise ValueError(f"Candidate bases exceed supported record lengths: {motif}")
        clean = {key: value for key, value in row.items() if key not in ANNOTATIONS}
        # JSON serialization retains numeric types, unlike Python's 1 == True.
        exact = json.dumps(clean, sort_keys=True, allow_nan=False, separators=(",", ":"))
        if motif in evidence and evidence[motif] != exact:
            raise ValueError(f"Conflicting duplicate candidate family: {motif}")
        evidence[motif] = exact
        motifs[motif] = {"reads": reads, "repeat_bases": repeat_bases}
    raw = len(motifs)
    if before < raw or (config["min_reads"] == 1 and before != raw):
        raise ValueError("Candidate-family summary indicates omitted or inconsistent raw rows")
    if "raw_reported_candidates" in summary:
        if integer(summary["raw_reported_candidates"], "summary.raw_reported_candidates") != raw:
            raise ValueError("raw_reported_candidates differs from unique candidate families")
    elif report["schema_version"] >= 2 or rejected:
        raise ValueError("Prior/support-classified reports require raw_reported_candidates")
    elif reported != raw:
        raise ValueError("Original scan report is missing raw candidate families")
    return {
        "control_id": path.stem,
        "source_result": str(path),
        "source_result_sha256": hashlib.sha256(contents).hexdigest(),
        "config": config,
        "records": records,
        "bases": bases,
        "max_record_bases": longest,
        "sampled_prefix": False,
        "inputs": inputs,
        "motifs": dict(sorted(motifs.items())),
    }


def build_model(negative_results, name=None):
    if not negative_results:
        raise ValueError("At least one explicitly selected negative result is required")
    if name is not None and (not isinstance(name, str) or not name.strip()):
        raise ValueError("Model name must not be blank")
    controls = []
    seen = {key: set() for key in ("control_id", "source_result", "source_result_sha256")}
    reference = None
    for source in negative_results:
        control = load_control(source)
        for key, values in seen.items():
            if control[key] in values:
                raise ValueError(f"Duplicate negative control {key}: {control[key]}")
            values.add(control[key])
        detector = {key: control["config"][key] for key in DETECTOR_KEYS}
        if reference is not None and detector != reference:
            raise ValueError("Negative controls have incompatible detector configurations")
        reference = detector
        controls.append(control)
    return {
        "schema_version": 1,
        "model_type": "empirical_repeat_background",
        "name": name.strip() if name else "explicit_negative_controls",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "limitations": [
            "Controls were explicitly designated negative by the caller; biological labels are not inferred or verified.",
            "Only specified negative results are read; thresholds are not fitted using positive or held-out results.",
            "Runtime z=1.96 and minimum enrichment=2 are configurable scoring heuristics, not calibrated confidence, p-values, or FDR.",
            "Reads are sequence records; paired reads and PCR duplicates are not independent biological molecules.",
            "Read and repeat-base fractions are not genomic coverage; read length and library preparation affect comparisons.",
            "Missing motifs are censored below each control's config.min_reads, not necessarily observed zero.",
            "Enrichment relative to these controls does not establish telomeric location or specificity across taxa.",
        ],
        "controls": controls,
    }


def write_model(negative_results, output, name=None):
    destination = Path(output).absolute()
    sources = [Path(path).resolve(strict=True) for path in negative_results]
    if destination.resolve() in sources:
        raise ValueError("Output must not overwrite a source result")
    # Also refuse empty existing files; this permits atomic, no-clobber publication.
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Refusing to overwrite existing output: {destination}")
    model = build_model(sources, name)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=destination.parent,
                                         prefix=".background-", suffix=".part", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(model, handle, indent=2, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return model


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--negative-result", type=Path, action="append", required=True,
                        help="Saved negative-control JSON; repeat to provide multiple controls")
    parser.add_argument("--output", type=Path, required=True, help="New background-model JSON")
    parser.add_argument("--name", help="Descriptive model name")
    args = parser.parse_args(argv)
    try:
        model = write_model(args.negative_result, args.output, args.name)
    except (ValueError, OSError) as error:
        print(f"calibrate_support: {error}", file=sys.stderr)
        return 1
    print(f"Wrote {len(model['controls'])} negative-control profiles to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
