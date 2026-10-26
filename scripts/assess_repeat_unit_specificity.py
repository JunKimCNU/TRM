#!/usr/bin/env python3
"""Assess supplied motif hypotheses in registered actual FASTA/FASTQ intervals.

Schema1 plan: case_id, data_role, expected_source{sha256,size_bytes,records,bases},
regions[{id,record_id,record_length,start0,end0,motif,strand(optional)}],
config(optional resource bounds), references(optional hash-bound provenance).
Use --input for the source path, or source.path in the plan. Every source record
is streamed and the complete on-disk bytes, record count and base count checked.
This optional diagnostic never updates native results or classifies TRMs.
"""
import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
from pathlib import Path
import resource
import tempfile
import time

import audit_irregular_arrays as audit
from anchor_store import apply_memory_limit
from report_tract_context import binding, records_with_quality, verify_reference
import repeat_region_continuity
import repeat_unit_specificity as specificity
import tract_context

VERSION = specificity.VERSION
ROLES = {"development_previously_seen", "independent_validation", "unspecified"}
TSV_FIELDS = ["id", "record_id", "record_ordinal", "start0", "end0", "motif", "primitive_motif",
              "status", "reason_codes", "best_fixed_identity", "fixed_phase_ties",
              "majority_identity", "fixed_normalized_gain", "best_cyclic_edit_cost",
              "edit_normalized_gain", "operational_model_disagreement"]


def dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def text_hash(value):
    return hashlib.sha256(value.encode()).hexdigest()


def validate_plan(plan, args):
    if not isinstance(plan, dict) or type(plan.get("schema_version")) is not int or plan["schema_version"] != 1:
        raise ValueError("plan schema_version1 required")
    if not isinstance(plan.get("case_id"), str) or not plan["case_id"].strip():
        raise ValueError("nonempty case_id required")
    if plan.get("data_role") not in ROLES:
        raise ValueError("explicit data_role required")
    expected = plan.get("expected_source")
    if not isinstance(expected, dict) or not audit.SHA.fullmatch(str(expected.get("sha256", ""))):
        raise ValueError("expected source SHA256 required")
    for key in ("size_bytes", "records", "bases"):
        audit.integer(expected.get(key), "expected source " + key, 1)
    config = specificity.validate_config(plan.get("config"))
    regions = plan.get("regions")
    if not isinstance(regions, list) or not 1 <= len(regions) <= args.max_regions:
        raise ValueError("nonempty regions within max-regions required")
    ids, lengths = set(), {}
    total_bp = total_work = 0
    for row in regions:
        if not isinstance(row, dict) or not isinstance(row.get("id"), str) or not row["id"] or row["id"] in ids:
            raise ValueError("unique nonempty interval id required")
        ids.add(row["id"])
        rid = row.get("record_id")
        if not isinstance(rid, str) or not rid or any(c.isspace() for c in rid):
            raise ValueError("record_id must be a nonempty token")
        length = audit.integer(row.get("record_length"), "record length", 1)
        if length > args.max_record_bases:
            raise ValueError("registered source record exceeds max-record-bases")
        start, end = [audit.integer(row.get(k), k) for k in ("start0", "end0")]
        if not 0 <= start <= end <= length:
            raise ValueError("interval outside registered record")
        if rid in lengths and lengths[rid] != length:
            raise ValueError("conflicting registered lengths for one record")
        lengths[rid] = length
        if row.get("strand") not in (None, "+", "-"):
            raise ValueError("strand must be +, -, or omitted")
        if not isinstance(row.get("motif"), str):
            raise ValueError("supplied motif must be a string")
        # Invalid/unsupported motifs are preserved as not_evaluated diagnostics.
        m = row["motif"].upper()
        k = len(specificity.primitive_unit(m)) if m and set(m) <= set("ACGT") and len(m) <= specificity.MAX_MOTIF_BP else 0
        size = end - start
        total_bp += size
        if size <= config["max_region_bp"] and 2 * size * k <= config["max_dp_cells"]:
            total_work += 2 * size * k
    if total_bp > args.max_total_region_bases:
        raise ValueError("total requested interval bases exceed registered CLI bound")
    if total_work > args.max_total_dp_cells:
        raise ValueError("total evaluated DP cell estimate exceeds registered CLI bound")
    return lengths, config, total_work


def run(args):
    start = time.monotonic()
    initial = resource.getrusage(resource.RUSAGE_SELF)
    memory = apply_memory_limit(args.memory_limit_mib)
    for key in ("max_regions", "max_record_bases", "max_total_region_bases", "max_total_dp_cells"):
        audit.integer(getattr(args, key), key, 1)
    plan_path = Path(args.plan).resolve(strict=True)
    if plan_path.stat().st_size > 16 * 1024 * 1024:
        raise ValueError("plan exceeds16MiB")
    plan_pin = binding(plan_path)
    plan = audit.parse_json(plan_path.read_text())
    lengths, config, work = validate_plan(plan, args)
    source_name = args.input or plan.get("source", {}).get("path")
    if not source_name:
        raise ValueError("--input or plan source.path required")
    source = Path(source_name)
    if not source.is_absolute():
        source = (Path.cwd() if args.input else plan_path.parent) / source
    source = source.resolve(strict=True)
    output = Path(args.output).resolve()
    if output.exists():
        raise ValueError("output already exists; use a new directory")
    refs = [verify_reference(x, plan_path.parent) for x in plan.get("references", [])]
    modules = [Path(__file__), Path(specificity.__file__), Path(repeat_region_continuity.__file__),
               Path(tract_context.__file__), Path(audit.__file__),
               Path(__file__).with_name("report_tract_context.py"), Path(__file__).with_name("anchor_store.py")]
    code = [binding(x) for x in modules]
    by_record = defaultdict(list)
    for row in plan["regions"]:
        by_record[row["record_id"]].append(row)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="." + output.name + "-", dir=output.parent) as tmp:
        stage = Path(tmp) / "report"
        stage.mkdir()
        seen = set()
        count = bases = emitted = 0
        statuses = Counter()
        with ExitStack() as stack:
            verified = audit.VerifiedStream(stack, source, plan["expected_source"])
            details = stack.enter_context((stage / "details.jsonl").open("w"))
            table = stack.enter_context((stage / "regions.tsv").open("w", newline=""))
            writer = csv.DictWriter(table, fieldnames=TSV_FIELDS, delimiter="\t")
            writer.writeheader()
            for rid, sequence, quality in records_with_quality(verified.stream, args.max_record_bases):
                count += 1
                bases += len(sequence)
                if rid not in lengths:
                    continue
                if rid in seen:
                    raise ValueError("duplicate selected record_id: " + rid)
                seen.add(rid)
                if len(sequence) != lengths[rid]:
                    raise ValueError("record length differs from registration: " + rid)
                record_hash = text_hash(sequence)
                quality_hash = text_hash(quality) if quality is not None else None
                for spec in by_record[rid]:
                    lo, hi = spec["start0"], spec["end0"]
                    result = specificity.evaluate_unit_specificity(sequence[lo:hi], spec["motif"], start0=lo, config=config)
                    payload = {"id": spec["id"], "record_id": rid, "registered_region": spec,
                               "registered_strand": spec.get("strand"),
                               "registered_strand_role": "metadata only; both motif strands are evaluated on the original forward interval; this field does not orient or filter the sequence",
                               "source_record_ordinal": count, "record_length_bp": len(sequence),
                               "full_record_normalized_sequence_sha256": record_hash,
                               "full_record_quality_sha256": quality_hash,
                               "full_header_preserved": False,
                               "header_scope": "Existing parser preserves identifier token; full source bytes are SHA-bound. Header comments are not exported.",
                               "sequence_parser_normalization": "Existing audit parser uppercases DNA and strips FASTA sequence whitespace; coordinates refer to normalized biological sequence.",
                               "left_record_distance_bp": lo, "right_record_distance_bp": len(sequence) - hi,
                               "record_edge_is_not_physical_chromosome_end": True,
                               "quality": tract_context.quality_summary(quality[lo:hi] if quality is not None else None),
                               "analysis": result}
                    details.write(json.dumps(payload, ensure_ascii=False, allow_nan=False) + "\n")
                    fixed, majority, edit, comparison = [result[k] for k in ("fixed_phase", "majority_baseline", "cyclic_edit", "model_comparison")]
                    row = {"id": spec["id"], "record_id": rid, "record_ordinal": count, "start0": lo, "end0": hi,
                           "motif": spec["motif"], "primitive_motif": result["primitive_motif"], "status": result["status"],
                           "reason_codes": json.dumps(result["reason_codes"]),
                           "best_fixed_identity": fixed["all_best_ties"][0]["identity"] if fixed else None,
                           "fixed_phase_ties": fixed["all_best_tie_count"] if fixed else None,
                           "majority_identity": majority["identity"] if majority else None,
                           "fixed_normalized_gain": comparison["fixed_normalized_gain"] if comparison else None,
                           "best_cyclic_edit_cost": edit["best_edit_cost"] if edit else None,
                           "edit_normalized_gain": comparison["edit_normalized_gain"] if comparison else None,
                           "operational_model_disagreement": comparison["operational_model_disagreement"] if comparison else None}
                    writer.writerow(row)
                    emitted += 1
                    statuses[result["status"]] += 1
            input_pin = verified.finish()
        if seen != set(lengths):
            raise ValueError("registered source record missing")
        if count != plan["expected_source"]["records"] or bases != plan["expected_source"]["bases"]:
            raise ValueError("full source record/base count mismatch")
        if plan_pin != binding(plan_path):
            raise ValueError("plan changed during analysis")
        for ref in refs + code:
            verify_reference(ref, plan_path.parent)
        assert emitted == len(plan["regions"])
        usage = resource.getrusage(resource.RUSAGE_SELF)
        summary = {"status": "complete", "program": "trm-repeat-unit-specificity", "version": VERSION,
                   "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
                   "case_id": plan["case_id"], "data_role": plan["data_role"], "config": config,
                   "source": input_pin, "plan": plan_pin, "references": refs, "code": code,
                   "input_records": count, "input_bases": bases, "selected_records": len(seen),
                   "regions": emitted, "diagnostic_status_counts": dict(statuses),
                   "evaluated_DP_cell_estimate": work,
                   "limits": {k: getattr(args, k) for k in ("max_regions", "max_record_bases", "max_total_region_bases", "max_total_dp_cells")},
                   "runtime": {"wall_seconds": time.monotonic() - start,
                               "user_cpu_seconds": usage.ru_utime - initial.ru_utime,
                               "system_cpu_seconds": usage.ru_stime - initial.ru_stime,
                               "maximum_process_rss_kib": usage.ru_maxrss, **memory},
                   "changes_native_classification": False, "novel_trm_confirmed": False,
                   "registered_strand_is_metadata_only": True,
                   "physical_terminal_location_confirmed": False, "negative_detection_power_established": False,
                   "limitations": ["Fixed operational diagnostic settings, not calibrated sensitivity/specificity or p-values.",
                                   "Conditional on supplied motif and selected interval; not de novo unit discovery.",
                                   "Phase selection inflates best fixed matches; composition-null expectation is not selection-adjusted.",
                                   "Indel-tolerant edit cost has no traceback and its end states are not independent paths.",
                                   "Quality values are stored observations, not calibrated biological error rates.",
                                   "No automatic motif exclusion, TRM confirmation, biological absence or novelty call."]}
        summary["outputs"] = [{**binding(p), "path": str(output / p.name)} for p in sorted(stage.iterdir())]
        dump(stage / "summary.json", summary)
        os.rename(stage, output)
    return {"status": "complete", "output": str(output), "regions": emitted, "diagnostic_status_counts": dict(statuses)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="trm-repeat-unit-specificity " + VERSION)
    parser.add_argument("--input", help="Hash-bound FASTA/FASTQ, optionally gzip; or use plan source.path")
    parser.add_argument("--plan", required=True)
    parser.add_argument("--output", required=True, help="New directory; no overwrite")
    parser.add_argument("--memory-limit-mib", type=int, default=None)
    parser.add_argument("--max-record-bases", type=int, default=1_000_000)
    parser.add_argument("--max-regions", type=int, default=1000)
    parser.add_argument("--max-total-region-bases", type=int, default=2_000_000)
    parser.add_argument("--max-total-dp-cells", type=int, default=100_000_000)
    args = parser.parse_args(argv)
    try:
        print(json.dumps(run(args), ensure_ascii=False))
    except (ValueError, OSError, KeyError, TypeError, MemoryError, EOFError) as error:
        parser.exit(1, "repeat-unit-specificity: " + str(error) + "\n")


if __name__ == "__main__":
    main()
