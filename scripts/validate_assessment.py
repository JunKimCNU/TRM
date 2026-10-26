#!/usr/bin/env python3
"""Reproducible CLI-to-assessment integration validation, without threshold fitting.

Uses seeded artificial sequences and the existing Drosophila control. Writes
all fixtures, executed commands, measurements, code snapshots, and assertions
to a new output directory; refuses to overwrite a previous validation.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
SEED = 20260921
UNKNOWN = "ACGATCG"
CODE_FILES = (
    "scripts/validate_assessment.py", "scripts/evidence_batch.py", "scripts/runtime_binary.py",
    "scripts/trm_assessment.py", "scripts/assess_adequacy.py",
    "scripts/report_assessment.py", "data/telomere_motifs.tsv", "data/controls/controls_manifest.json",
    "src/main.rs", "src/detect.rs", "src/evidence.rs", "src/motif.rs", "src/input.rs",
    "src/prior.rs", "src/support.rs", "Cargo.lock", "Cargo.toml",
)


def digest(path):
    hasher = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def rc(motif):
    return motif.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def tandem(motif, length, phase=0):
    return (motif * ((length + phase) // len(motif) + 2))[phase:phase + length]


def random_bases(rng, length):
    return "".join(rng.choice("ACGT") for _ in range(length))


def repeat_fixture(rng, motif, mode, full=20, each=16):
    sequences = [tandem(motif if i % 2 else rc(motif), 150, i % len(motif)) for i in range(full)]
    if mode == "full_only":
        return sequences
    if mode == "internal":
        return sequences + [random_bases(rng, 40) + tandem(motif, 70, i % len(motif))
                            + random_bases(rng, 40) for i in range(each * 2)]
    for index in range(each):
        sequences.extend([
            tandem(rc(motif), 110, index % len(motif)) + random_bases(rng, 40),
            random_bases(rng, 40) + tandem(motif, 110, index % len(motif)),
        ])
    if mode in {"random_orientation", "partial_orientation"}:
        opposite = each if mode == "random_orientation" else each // 2
        for index in range(opposite):
            sequences.extend([
                tandem(motif, 110, index % len(motif)) + random_bases(rng, 40),
                random_bases(rng, 40) + tandem(rc(motif), 110, index % len(motif)),
            ])
    return sequences


def fastq(path, sequences):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for index, sequence in enumerate(sequences):
            handle.write(f"@seeded_{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def timed_command(command, log_prefix):
    started = time.monotonic()
    completed = subprocess.run(command, text=True, capture_output=True, check=False, timeout=120)
    wall = time.monotonic() - started
    log_prefix.parent.mkdir(parents=True, exist_ok=True)
    log_prefix.with_suffix(".stdout.log").write_text(completed.stdout)
    log_prefix.with_suffix(".stderr.log").write_text(completed.stderr)
    return {"command": command, "exit_code": completed.returncode, "wall_seconds": wall,
            "stderr_path": str(log_prefix.with_suffix(".stderr.log")),
            "stdout_path": str(log_prefix.with_suffix(".stdout.log"))}


def main():
    from runtime_binary import resolve_binary
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/assessment_validation_v050")
    parser.add_argument("--binary", help="Explicit scanner override; default: TRM_BINARY or verified bin/trm")
    parser.add_argument("--drosophila-discovery", type=Path,
                        default=ROOT / "results/evidence_engine_smoke_20260921/scan.json")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    try:
        args.binary = resolve_binary(args.binary)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    args.output = args.output.resolve()
    if args.output.exists():
        parser.error("Output already exists; choose a new directory to preserve the earlier audit")
    if not 1 <= args.threads <= 16:
        parser.error("threads must be 1..16")
    args.output.mkdir(parents=True)
    started = time.monotonic()
    result = {"schema_version": 1, "status": "running", "validation_kind": "cli_and_scientific_assessment_integration",
              "started_utc": datetime.now(timezone.utc).isoformat(), "random_seed": SEED,
              "fixtures": [], "checks": [], "scenarios": [], "source_snapshots": [],
              "limitations": [
                  "Synthetic fixtures test software behavior and explicit model assumptions, not biological sensitivity or FDR.",
                  "The same-control Drosophila check is not an independent validation cohort and cannot calibrate thresholds.",
                  "All timing is measured wall time under the current machine load, including concurrent WGS analysis.",
                  "No check infers true telomere absence, ALT, or physical chromosome-end localization.",
              ]}
    frozen = args.output / "snapshot"
    for relative in CODE_FILES:
        source = ROOT / relative
        if not source.is_file():
            continue
        destination = frozen / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        result["source_snapshots"].append({"original_path": str(source),
                                            "snapshot": str(destination), "sha256": digest(destination)})
    binary = frozen / "trm"
    shutil.copy2(args.binary, binary)
    result["binary"] = {"original_path": str(args.binary.resolve()), "snapshot": str(binary),
                        "sha256": digest(binary), "version": subprocess.check_output(
                            [str(binary), "--version"], text=True).strip()}
    # Import the frozen modules so ongoing edits cannot change this run midway.
    sys.path.insert(0, str(frozen / "scripts"))
    from evidence_batch import canonical, read_tsv, select_targets
    from trm_assessment import catalog_index, POLICY
    from assess_adequacy import REQUIRED_ASSUMPTIONS
    from report_assessment import sample_report
    catalog = frozen / "data/telomere_motifs.tsv"
    references = catalog_index(read_tsv(catalog))
    known = {family for family in references if len(family) <= 32}
    result["assessment_policy"] = POLICY
    result["known_target_family_count"] = len(known)
    result["unknown_fixture_family"] = canonical(UNKNOWN)
    if canonical(UNKNOWN) in references:
        raise ValueError("The predeclared novel fixture motif is in this catalogue; no silent replacement")
    rng = random.Random(SEED)
    generated = {
        "known_ttagac": repeat_fixture(rng, "TTAGAC", "positive"),
        "known_ttaggt": repeat_fixture(rng, "TTAGGT", "positive"),
        "novel_both_boundaries": repeat_fixture(rng, UNKNOWN, "positive"),
        "novel_internal_negative": repeat_fixture(rng, UNKNOWN, "internal"),
        "novel_random_orientation_negative": repeat_fixture(rng, UNKNOWN, "random_orientation"),
        "known_full_only": repeat_fixture(rng, "TTAGGT", "full_only"),
        "known_partial_with_secondary_novel": repeat_fixture(rng, "TTAGGT", "partial_orientation", full=80, each=48)
            + repeat_fixture(rng, UNKNOWN, "positive", full=10),
        "random_negative": [random_bases(rng, 150) for _ in range(1000)],
        "one_read_negative": [random_bases(rng, 150)],
    }
    metadata = {"status": "synthetic", "scientific_name": "Artificial sequence fixture", "tax_id": "",
                "lineage": "", "library_strategy": "WGS", "library_source": "GENOMIC",
                "library_layout": "SINGLE"}
    assignment_file = args.output / "assignments.tsv"
    assignment_file.write_text("sample_id\tscientific_name\ttax_id\tlineage\trelated_motifs\n")
    report_args = SimpleNamespace(source=args.output / "discovery", evidence=args.output / "evidence",
                                  catalog=catalog, assignments=assignment_file, output=args.output / "reports")

    def check(name, condition, actual=None, expected=None):
        result["checks"].append({"name": name, "passed": bool(condition), "actual": actual, "expected": expected})

    def analysis(sample, paths, saved_scan=None):
        directory = report_args.source / sample
        directory.mkdir(parents=True, exist_ok=True)
        scan_path = directory / "result.json"
        measurement = {"sample_id": sample, "inputs": [str(p) for p in paths]}
        if saved_scan:
            shutil.copy2(saved_scan, scan_path)
            measurement["discovery_reused"] = {"path": str(saved_scan), "sha256": digest(saved_scan)}
        else:
            command = [str(binary), "scan", *map(str, paths), "--output", str(directory / "result"),
                       "--threads", str(args.threads), "--prior-policy", "discovery"]
            measurement["scan"] = timed_command(command, directory / "scan")
            if measurement["scan"]["exit_code"]:
                result["fixtures"].append(measurement)
                raise RuntimeError(f"Unexpected scan failure for {sample}")
        raw = json.loads(scan_path.read_text())
        targets = select_targets(raw, known, max_novel=96)
        evidence_dir = report_args.evidence / sample
        evidence_dir.mkdir(parents=True, exist_ok=True)
        targets_path = evidence_dir / "targets.json"
        write_json(targets_path, targets)
        command = [str(binary), "evidence", *map(str, paths), "--targets", str(targets_path),
                   "--output", str(evidence_dir / "evidence"), "--threads", str(args.threads)]
        # Match the original detector including the control's min_reads=1.
        for field, flag in (("seed_bases", "--seed-bases"),
                            ("min_copies", "--min-copies"), ("min_identity", "--min-identity"), ("min_reads", "--min-reads")):
            command.extend([flag, str(raw["config"][field])])
        measurement["evidence"] = timed_command(command, evidence_dir / "run")
        measurement["detector_config"] = raw["config"]
        measurement["selection"] = targets["selection"]
        if measurement["evidence"]["exit_code"]:
            result["fixtures"].append(measurement)
            raise RuntimeError(f"Unexpected evidence failure for {sample}")
        evidence_path = evidence_dir / "evidence.json"
        targeted = json.loads(evidence_path.read_text())
        status = {"status": "complete", "identity": {"source_sha256": digest(scan_path), "targets": targets},
                  "result_sha256": digest(evidence_path), "command": command}
        write_json(evidence_dir / "status.json", status)
        measurement.update(records=raw["summary"]["records"], bases=raw["summary"]["bases"],
                           source_sha256=digest(scan_path), evidence_sha256=digest(evidence_path),
                           inputs_sha256=[{"path": str(p), "sha256": digest(p), "bytes": p.stat().st_size} for p in paths])
        result["fixtures"].append(measurement)
        for field in ("records", "bases"):
            check(f"{sample}: complete {field}", raw["summary"][field] == targeted["summary"][field],
                  targeted["summary"][field], raw["summary"][field])
        check(f"{sample}: complete input", targeted["summary"].get("sampled_prefix") is False)
        check(f"{sample}: exact target family set", set(targets["motifs"]) == {c["motif"] for c in targeted["candidate_evidence"]})
        return raw, targeted

    def assessment(sample, scenario=None, meta=None, assumptions=None, related="", native=False):
        name = scenario or sample
        meta = dict(meta or metadata)
        assignment = {key: meta.get(key, "") for key in ("scientific_name", "tax_id", "lineage")}
        assignment.update(sample_id=sample, related_motifs=related)
        local_args = SimpleNamespace(**vars(report_args))
        if native:
            local_args.evidence = None
        begin = time.monotonic()
        document, _, _ = sample_report(sample, meta, assignment, local_args, references, assumptions, {})
        path = report_args.output / f"{name}.json"
        write_json(path, document)
        result["scenarios"].append({"scenario": name, "fixture": sample, "classification": document["classification"],
            "data_adequacy": document["data_adequacy"]["status"], "detected_motif": document["detected_motif"],
            "novel_secondary_candidate_count": sum(c["novel_candidate"] for c in document["candidates"]),
            "wall_seconds": time.monotonic() - begin, "report": str(path), "sha256": digest(path),
            "assumptions": assumptions, "metadata": meta})
        check(f"{name}: no absence inferred", document["data_adequacy"]["true_absence_inferred"] is False)
        check(f"{name}: no ALT inferred", document["data_adequacy"]["alternative_maintenance_inferred"] is False)
        return document

    for sample, sequences in generated.items():
        path = args.output / "inputs" / f"{sample}.fastq"
        fastq(path, sequences)
        analysis(sample, [path])
    outputs = {sample: assessment(sample, related="TTAGAC" if sample == "known_ttagac" else
                                  "TTAGGT" if sample.startswith("known_") else "") for sample in generated}
    for sample, motif in (("known_ttagac", "TTAGAC"), ("known_ttaggt", "TTAGGT")):
        observed = outputs[sample]
        check(f"{sample}: established motif and boundary pattern", observed["classification"] == "known_trm_detected"
              and observed["top_candidate"]["canonical_family"] == canonical(motif), observed["classification"], "known_trm_detected")
    novel = outputs["novel_both_boundaries"]
    check("novel: uncatalogued family nominated", novel["classification"] == "novel_provisional_candidate"
          and novel["top_candidate"]["canonical_family"] == canonical(UNKNOWN), novel["classification"], "novel_provisional_candidate")
    for sample in ("novel_internal_negative", "novel_random_orientation_negative"):
        check(f"{sample}: no novel primary", outputs[sample]["classification"] not in {"novel_geometry_nomination", "novel_provisional_candidate"},
              outputs[sample]["classification"], "not novel_provisional_candidate")
        candidate = next(c for c in outputs[sample]["candidates"] if c["canonical_family"] == canonical(UNKNOWN))
        check(f"{sample}: same motif retained as unresolved", candidate["candidate_status"] == "uncatalogued_repeat",
              candidate["candidate_status"], "uncatalogued_repeat")
    full_only = outputs["known_full_only"]
    check("full-only known: boundary unresolved", full_only["classification"] == "known_repeat_detected_boundary_unresolved",
          full_only["classification"], "known_repeat_detected_boundary_unresolved")
    check("full-only known: zero junction records", full_only["top_candidate"]["unique_junction_reads"] == 0)
    mixed = outputs["known_partial_with_secondary_novel"]
    check("partial known dominates secondary novelty", mixed["classification"] == "known_repeat_detected_boundary_unresolved"
          and mixed["top_candidate"]["canonical_family"] == canonical("TTAGGT"), mixed["classification"],
          "known_repeat_detected_boundary_unresolved")
    check("partial known retains secondary novel evidence", any(c["novel_candidate"] for c in mixed["candidates"]))
    check("one-read negative: insufficient data", outputs["one_read_negative"]["data_adequacy"]["status"] == "insufficient_data")
    for sample in ("known_ttagac", "known_ttaggt", "novel_both_boundaries", "novel_internal_negative", "known_full_only"):
        native = assessment(sample, scenario=sample + "_native", native=True)
        check(f"{sample}: native and targeted report classifications agree", native["classification"] == outputs[sample]["classification"],
              native["classification"], outputs[sample]["classification"])
    assumptions = {"genome_size_bp": 100_000, "haploid_chromosome_count": 6, "random_wgs_assumed": True,
                   "telomere_retention_fraction": 1.0, "effective_independent_read_fraction": 1.0,
                   "detector_sensitivity": 1.0,
                   "provenance": {key: "Explicit synthetic scenario hypothesis; not a measured biological parameter"
                                  for key in REQUIRED_ASSUMPTIONS}}
    powerful = assessment("random_negative", "explicit_hypothetical_power", assumptions=assumptions)
    check("explicit assumptions: power supported without absence claim", powerful["classification"] == "canonical_not_detected_conditional_power",
          powerful["classification"], "canonical_not_detected_conditional_power")
    weak = assessment("random_negative", "large_genome_low_power", assumptions={**assumptions, "genome_size_bp": 10_000_000_000})
    check("large genome: insufficient conditional power", weak["data_adequacy"]["status"] == "insufficient_data",
          weak["data_adequacy"]["status"], "insufficient_data")
    wga = assessment("random_negative", "explicit_wga_metadata", assumptions=assumptions,
                     meta={**metadata, "amplification_method": "MDA", "library_strategy": "WGA"})
    check("explicit WGA: uniform WGS model unassessable", wga["data_adequacy"]["status"] == "unassessable_library_or_genome",
          wga["data_adequacy"]["status"], "unassessable_library_or_genome")
    observed_wga = assessment("known_ttaggt", "known_observed_despite_wga", assumptions=assumptions,
                              meta={**metadata, "amplification_method": "MDA", "library_strategy": "WGA"})
    check("WGA limitation does not erase known observed signal", observed_wga["classification"] == "known_trm_detected")
    damaged = args.output / "inputs/damaged.fastq"
    damaged.write_text("@broken\nTTAGGGTTAGGGTTAGGGTTAGGG\n+\nIIII\n")
    damaged_dir = report_args.source / "damaged"
    damaged_dir.mkdir(parents=True)
    for mode in ("scan", "evidence"):
        prefix = damaged_dir / ("result" if mode == "scan" else "evidence")
        command = [str(binary), mode, str(damaged), "--output", str(prefix), "--threads", str(args.threads)]
        if mode == "evidence":
            command.extend(["--targets", str(report_args.evidence / "known_ttaggt/targets.json")])
        executed = timed_command(command, damaged_dir / mode)
        result["fixtures"].append({"sample_id": "damaged_" + mode, "input_sha256": digest(damaged), "run": executed})
        check(f"damaged FASTQ: {mode} fails and has no final report", executed["exit_code"] != 0
              and not prefix.with_suffix(".json").exists(), executed["exit_code"], "nonzero")
    broken = assessment("damaged")
    check("damaged FASTQ: technical failure not biological nondetection", broken["classification"] == "analysis_failed",
          broken["classification"], "analysis_failed")
    write_json(args.output / "results.json", result)

    # A real negative architecture control, using full input and untuned selection.
    control_source = args.drosophila_discovery.resolve()
    control_raw = json.loads(control_source.read_text())
    controls = [Path(item["path"] if isinstance(item, dict) else item) for item in control_raw["inputs"]]
    controls = [path if path.is_absolute() else ROOT / path for path in controls]
    control_manifest = json.loads((frozen / "data/controls/controls_manifest.json").read_text())["controls"]["drosophila"]
    pinned_files = {entry["name"]: entry for entry in control_manifest["files"]}
    if {path.name for path in controls} != set(pinned_files):
        raise ValueError("Discovery inputs are not the two pinned Drosophila control files")
    if any(digest(path) != pinned_files[path.name]["sha256"] for path in controls):
        raise ValueError("Drosophila input hash differs from previously verified archive download")
    result["control_metadata_provenance"] = control_manifest
    check("real Drosophila: both FASTQ hashes match pinned archive downloads", True)
    analysis("drosophila_SRR1183635", controls, saved_scan=control_source)
    control_meta = {"status": "resolved", "scientific_name": control_manifest["scientific_name"], "tax_id": "7227",
                    "lineage": "Eukaryota; Metazoa; Arthropoda; Insecta; Diptera; Drosophilidae",
                    "library_strategy": control_manifest["library_strategy"], "library_source": control_manifest["library_source"],
                    "library_selection": control_manifest["library_selection"], "library_layout": "PAIRED",
                    "metadata_basis": "Existing downloaded SRR1183635 control metadata; no assumption of independent molecules or unamplified library"}
    control = assessment("drosophila_SRR1183635", meta=control_meta)
    check("real Drosophila: no novel primary", control["classification"] not in {"novel_geometry_nomination", "novel_provisional_candidate"},
          control["classification"], "not novel_provisional_candidate")
    check("real Drosophila: unknown assumptions prevent absence claim", control["data_adequacy"]["status"] == "unassessable_library_or_genome")
    result["finished_utc"] = datetime.now(timezone.utc).isoformat()
    result["status"] = "complete"
    result["total_wall_seconds"] = time.monotonic() - started
    result["passed"] = all(check["passed"] for check in result["checks"])
    result["passed_checks"] = sum(check["passed"] for check in result["checks"])
    result["total_checks"] = len(result["checks"])
    result["production_sources_changed_during_validation"] = [row["original_path"] for row in result["source_snapshots"]
        if digest(row["original_path"]) != row["sha256"]]
    write_json(args.output / "results.json", result)
    print(json.dumps({key: result[key] for key in ("passed", "passed_checks", "total_checks", "total_wall_seconds")}, indent=2))
    for item in result["checks"]:
        if not item["passed"]:
            print(json.dumps(item, ensure_ascii=False))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
