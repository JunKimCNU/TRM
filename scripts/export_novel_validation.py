#!/usr/bin/env python3
"""Export traceable follow-up targets from a v0.5 assessment and native scan.

This reads JSON only. It does not map flanks, design primers, contact SRA, or
confirm chromosome-terminal localization. Missing measurements remain null.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import urlencode


LEVELS = {"geometry_nomination", "provisional_candidate", "externally_confirmed"}
QUADRANTS = ("left_forward", "left_reverse", "right_forward", "right_reverse", "ambiguous")
FIELDS = (
    "sample_id", "scientific_name", "canonical_family", "display_motif", "source_evidence_level",
    "input_scope", "validation_availability", "reads", "full_repeat_reads", "junction_reads",
    "left_forward_junction_reads", "left_reverse_junction_reads", "right_forward_junction_reads",
    "right_reverse_junction_reads", "left_ambiguous_junction_reads", "right_ambiguous_junction_reads",
    "clean_two_flank_repeat_records", "distinct_junction_sequences_lower_bound", "junction_sketch_saturated",
    "distinct_physical_flanks_lower_bound", "flank_sketch_saturated", "retained_flank_context_count",
    "exported_fasta_ids", "deduplicated_consistency_lower_bound", "deduplicated_consistency_upper_bound",
    "physical_flank_cluster_count", "source_missing_evidence", "required_proof_steps", "current_biosample",
    "sra_independent_search_url", "assessment_sha256", "raw_sha256",
)
LINK_FIELDS = ("sample_id", "canonical_family", "fasta_id", "sequence_kind", "anchor20", "signature",
               "events", "normalized_quadrant_counts", "source_flank_sketch_saturated")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def canonical(motif):
    if not isinstance(motif, str) or not motif or set(motif) - set("ACGT"):
        raise ValueError("motif must contain uppercase ACGT")
    for size in range(1, len(motif) + 1):
        if len(motif) % size == 0 and motif[:size] * (len(motif) // size) == motif:
            motif = motif[:size]
            break
    rc = motif.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return min(s[i:] + s[:i] for s in (motif, rc) for i in range(len(s)))


def nonnegative(value, label):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(label + " must be a nonnegative integer")
    return value


def validate_binding(assessment, raw, raw_sha, sample_id):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", sample_id):
        raise ValueError("invalid expected sample_id")
    if assessment.get("program") != "trm-assessment" or assessment.get("sample_id") != sample_id:
        raise ValueError("assessment sample/program identity mismatch")
    for source in (raw, assessment.get("metadata", {}), assessment.get("execution", {})):
        if source.get("sample_id") not in (None, "", sample_id):
            raise ValueError("inconsistent sample identity")
    provenance = assessment.get("provenance", {})
    if provenance.get("original_discovery_sha256") != raw_sha:
        raise ValueError("raw source SHA-256 does not match assessment provenance")
    original = provenance.get("original_discovery_result")
    if not isinstance(original, str) or not original:
        raise ValueError("assessment lacks native source path")
    if Path(original).parent.name != sample_id and raw.get("sample_id") != sample_id:
        raise ValueError("native source has no matching sample directory or explicit sample_id")
    if raw.get("program") != "trm" or raw.get("analysis_kind") != "de_novo_discovery":
        raise ValueError("expected native de novo scan, not targeted or unrelated output")
    # Separate targeted evidence would require its own hash and measurement file.
    evidence_sha = provenance.get("targeted_evidence_sha256")
    if evidence_sha is not None and evidence_sha != raw_sha:
        raise ValueError("assessment uses separate targeted evidence; native-only export is unsupported")
    raw_summary, assessed_summary = raw.get("summary", {}), assessment.get("input_summary", {})
    for field in ("records", "bases", "sampled_prefix"):
        if field in assessed_summary and assessed_summary[field] != raw_summary.get(field):
            raise ValueError("assessment/native input summary mismatch: " + field)
    prefix = raw_summary.get("sampled_prefix")
    return "sampled_prefix" if prefix is True else "full_input" if prefix is False else "not_available"


def validate_sketch(sketch, kind):
    observed_key, count_key = ("observed_events", "events") if kind == "flank_contexts" else ("observed_records", "records")
    required = ("capacity", "distinct_lower_bound", "saturated", "entries", observed_key)
    if not isinstance(sketch, dict) or any(key not in sketch for key in required):
        return False
    capacity = nonnegative(sketch["capacity"], kind + ".capacity")
    lower = nonnegative(sketch["distinct_lower_bound"], kind + ".distinct_lower_bound")
    observed = nonnegative(sketch[observed_key], kind + "." + observed_key)
    entries = sketch["entries"]
    if not capacity or not isinstance(entries, list) or not isinstance(sketch["saturated"], bool) or len(entries) > capacity:
        raise ValueError("invalid sketch shape: " + kind)
    total, seen = 0, set()
    for entry in entries:
        signature = entry.get("signature")
        if not isinstance(signature, str) or not signature or signature in seen:
            raise ValueError("invalid/duplicate retained signature: " + kind)
        seen.add(signature)
        count = nonnegative(entry.get(count_key), kind + "." + count_key)
        if not count:
            raise ValueError("retained signature has zero observations")
        total += count
        q = entry.get("normalized_quadrant_counts")
        if not isinstance(q, dict) or any(key not in q for key in QUADRANTS):
            return False
        for key in QUADRANTS:
            nonnegative(q[key], kind + ".quadrant." + key)
    if total > observed or lower < len(entries) or lower > observed:
        raise ValueError("contradictory sketch observation bounds")
    if not sketch["saturated"] and (lower != len(entries) or total != observed):
        raise ValueError("unsaturated sketch does not account for all observations")
    if sketch["saturated"] and (len(entries) != capacity or lower <= len(entries) or observed - total < lower - len(entries)):
        raise ValueError("saturated sketch lacks omitted-entry evidence")
    return True


def validation_status(native):
    junction = native.get("junction_evidence", {})
    validation = junction.get("validation")
    if not isinstance(validation, dict) or validation.get("measured") is not True:
        return "not_available", None
    if validation.get("schema_version") != 1 or "clean_two_flank_repeat_records" not in validation:
        return "incomplete", None
    nonnegative(validation["clean_two_flank_repeat_records"], "clean_two_flank_repeat_records")
    if not all(validate_sketch(validation.get(name), name) for name in ("junction_read_signatures", "flank_contexts")):
        return "incomplete", None
    if validation["flank_contexts"].get("flank_bases") != 20:
        return "incomplete", None
    if validation["junction_read_signatures"]["observed_records"] != junction.get("junction_reads"):
        raise ValueError("read sketch differs from native junction count")
    return "available", validation


def proof_steps(scope, availability):
    steps = [
        {"id": "independent_library_and_biosample", "required_before_confirmation": True,
         "status": "not_evaluated", "instruction": "Verify another biological specimen/BioSample and separately prepared genomic WGS or long-read library; compare ReadHash/RunHash and provenance. Another accession alone is not independent."},
        {"id": "local_assembly_or_long_read_terminal_localization", "required_before_confirmation": True,
         "status": "not_evaluated", "instruction": "Recover repeat-to-subtelomere continuity with local assembly or long reads, then establish physical chromosome-terminal position, outward polarity and array extent; inspect internal copies. A contig edge or short flank hit alone is insufficient."},
        {"id": "literature_and_taxonomy_verification", "required_before_confirmation": True,
         "status": "not_evaluated", "instruction": "Verify species/strain identity and original publications for this motif family and telomere architecture. Not in this catalogue does not establish literature-wide novelty."},
        {"id": "library_and_contamination_review", "required_before_confirmation": True,
         "status": "not_evaluated", "instruction": "Resolve WGS/library-prep/WGA/Hi-C/RNA annotations, contamination and PCR dependence using evidence; missing annotation does not demonstrate absence."},
    ]
    if scope != "full_input":
        steps.append({"id": "complete_full_input_scan", "required_before_confirmation": True,
                      "status": "missing", "instruction": "Run and verify complete input sequencing data; prefix observations cannot establish full-input nondetection power."})
    if availability != "available":
        steps.append({"id": "collect_v050_validation_evidence", "required_before_confirmation": True,
                      "status": "missing", "instruction": "Run v0.5 evidence collection to obtain distinct-read and physical-flank measurements; absent values are not zero."})
    return steps


def sra_search(assessment):
    meta = assessment.get("metadata", {})
    tax_id = str(assessment.get("tax_id") or meta.get("tax_id") or "")
    name = str(assessment.get("scientific_name") or meta.get("scientific_name") or "").replace('"', "")
    organism = "txid" + tax_id + "[Organism:exp]" if tax_id.isdigit() and int(tax_id) > 0 else '"' + name + '"[Organism]' if name else None
    query = organism + ' AND WGS[Strategy] AND genomic[Source]' if organism else None
    return {"query": query, "url": "https://www.ncbi.nlm.nih.gov/sra/?" + urlencode({"term": query}) if query else None,
            "current_biosample": assessment.get("independent_validation", {}).get("current_biosample") or meta.get("sample_accession"),
            "current_run": assessment["sample_id"], "network_request_performed": False,
            "alternative_dataset_verified": False,
            "instruction": "Search results require manual library/BioSample/provenance checks; the link does not select or verify an independent dataset."}


def build_export(assessment, raw, assessment_sha, raw_sha, sample_id):
    scope = validate_binding(assessment, raw, raw_sha, sample_id)
    search = sra_search(assessment)
    manifest = {"schema_version": 1, "program": "trm-novel-validation-export", "sample_id": sample_id,
                "assessment_sha256": assessment_sha, "raw_sha256": raw_sha,
                "assessment_version": assessment.get("version"), "native_version": raw.get("version"),
                "input_scope": scope, "sample_binding": "Explicit expected sample + assessment/metadata identity + native source hash and sample directory/ID; not an independent biological identity test.",
                "sra_search": search, "external_confirmation_verified": False,
                "reference_alignment_performed": False, "primer_design_performed": False,
                "limitations": ["Flanks and read signatures are sequence observations, not independent molecules, genomic loci or chromosome ends.",
                                "20-bp anchors and longer representatives may map nonuniquely; a short match alone cannot establish terminal localization.",
                                "Saturated sketches retain only part of the observed contexts; lower bounds are not exact distinct counts.",
                                "Outward-normalized flank orientation is preserved. FASTA deduplication uses exact strings, not reverse-complement or locus merging.",
                                "A representative flank is one retained sequence associated with a 20-bp anchor; it does not cover all alleles, errors or offsets.",
                                "Neither this export nor short-read geometry confirms a novel telomere or its maintenance mechanism."],
                "targets": [], "skipped_candidates": []}
    # v0.4 lacks the new evidence levels; do not reinterpret its legacy boolean.
    if not str(assessment.get("version", "")).startswith("0.5.") or assessment.get("schema_version") != 5:
        manifest.update(status="not_available", reason="v050_assessment_evidence_levels_required", target_count=None, exported_sequence_count=None)
        return manifest, [], [], {}
    native_families = {}
    for candidate in raw.get("candidates", []) + raw.get("rejected_candidates", []):
        family = canonical(candidate.get("motif"))
        if family in native_families:
            raise ValueError("duplicate native candidate family: " + family)
        native_families[family] = candidate
    seen, targets, links, sequences = set(), [], [], {}
    for candidate in assessment.get("candidates", []):
        family = candidate.get("canonical_family")
        if not isinstance(family, str) or canonical(family) != family or canonical(candidate.get("motif")) != family:
            raise ValueError("assessment candidate family/motif mismatch")
        if family in seen:
            raise ValueError("duplicate assessment family")
        seen.add(family)
        level = candidate.get("novel_evidence_level")
        if candidate.get("reference_status") != "not_in_catalogue" or level not in LEVELS:
            manifest["skipped_candidates"].append({"canonical_family": family, "reason": "not_unknown_geometry_or_stronger"})
            continue
        if family not in native_families:
            raise ValueError("selected assessment family absent from native scan: " + family)
        native = native_families[family]
        junction = native.get("junction_evidence", {})
        availability, validation = validation_status(native)
        direction = candidate.get("directionality", {})
        for key, value in direction.get("canonical_quadrants", {}).items():
            if value != junction.get(key):
                raise ValueError("assessment/native orientation count mismatch")
        readiness = candidate.get("confirmation_readiness", {})
        flanks = validation["flank_contexts"] if validation else None
        read_sketch = validation["junction_read_signatures"] if validation else None
        if validation:
            for key, value in (("distinct_junction_sequences_lower_bound", read_sketch["distinct_lower_bound"]),
                               ("distinct_physical_flanks_lower_bound", flanks["distinct_lower_bound"]),
                               ("flank_context_sketch_saturated", flanks["saturated"]),
                               ("clean_two_flank_repeat_records", validation["clean_two_flank_repeat_records"])):
                if readiness.get(key) is not None and readiness[key] != value:
                    raise ValueError("assessment/native validation count mismatch: " + key)
            if flanks["saturated"] and readiness.get("physical_flank_cluster_count") is not None:
                raise ValueError("saturated flank sketch cannot establish full-pool cluster count")
        fasta_ids, anchors = set(), set()
        for entry in flanks["entries"] if flanks else []:
            anchor = entry.get("sequence")
            if not isinstance(anchor, str) or len(anchor) != 20 or set(anchor) - set("ACGT") or anchor in anchors:
                raise ValueError("invalid or duplicate physical 20-bp anchor")
            anchors.add(anchor)
            packed = 0
            for base in anchor:
                packed = (packed << 2) | "ACGT".index(base)
            if entry["signature"] != format(packed, "010x"):
                raise ValueError("physical flank signature does not encode the reported anchor")
            representative = entry.get("representative_flank")
            if representative is None:
                sequence, sequence_kind = anchor, "anchor20_only"
            else:
                if not isinstance(representative, str) or not 20 <= len(representative) <= 80 or set(representative) - set("ACGT") or not representative.startswith(anchor):
                    raise ValueError("invalid representative physical flank")
                sequence, sequence_kind = representative, "representative_physical_flank"
            fasta_id = "flank_" + digest(sequence.encode())
            sequences[fasta_id] = sequence
            fasta_ids.add(fasta_id)
            links.append({"sample_id": sample_id, "canonical_family": family, "fasta_id": fasta_id,
                          "sequence_kind": sequence_kind, "anchor20": anchor, "signature": entry["signature"],
                          "events": entry["events"], "normalized_quadrant_counts": entry["normalized_quadrant_counts"],
                          "source_flank_sketch_saturated": flanks["saturated"]})
        steps = proof_steps(scope, availability)
        dedup = readiness.get("deduplicated_orientation", {}) if validation else {}
        target = {"sample_id": sample_id, "scientific_name": assessment.get("scientific_name"),
                  "canonical_family": family, "display_motif": candidate["motif"], "source_evidence_level": level,
                  "input_scope": scope, "validation_availability": availability,
                  "reads": native.get("reads"), "full_repeat_reads": native.get("full_repeat_reads"),
                  "junction_reads": junction.get("junction_reads"),
                  "clean_two_flank_repeat_records": validation["clean_two_flank_repeat_records"] if validation else None,
                  "distinct_junction_sequences_lower_bound": read_sketch["distinct_lower_bound"] if read_sketch else None,
                  "junction_sketch_saturated": read_sketch["saturated"] if read_sketch else None,
                  "distinct_physical_flanks_lower_bound": flanks["distinct_lower_bound"] if flanks else None,
                  "flank_sketch_saturated": flanks["saturated"] if flanks else None,
                  "retained_flank_context_count": len(flanks["entries"]) if flanks else None,
                  "exported_fasta_ids": sorted(fasta_ids),
                  "deduplicated_consistency_lower_bound": dedup.get("consistent_fraction_lower_bound"),
                  "deduplicated_consistency_upper_bound": dedup.get("consistent_fraction_upper_bound"),
                  "physical_flank_cluster_count": readiness.get("physical_flank_cluster_count") if validation else None,
                  "source_missing_evidence": readiness.get("missing_evidence", []),
                  "required_proof_steps": [step["id"] for step in steps],
                  "current_biosample": search["current_biosample"], "sra_independent_search_url": search["url"],
                  "assessment_sha256": assessment_sha, "raw_sha256": raw_sha}
        target.update({key: junction.get(key) for key in FIELDS if key.endswith("_junction_reads")})
        targets.append(target)
        manifest["targets"].append({**target, "missing_proof_checklist": steps,
            "source_directionality": direction, "source_confirmation_readiness": readiness,
            "native_validation": validation, "external_confirmation_verified": False})
    manifest.update(status="exported" if targets else "no_eligible_targets", target_count=len(targets), exported_sequence_count=len(sequences))
    return manifest, targets, links, sequences


def write_table(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "not_available" if value is None else json.dumps(value, separators=(",", ":"), sort_keys=True)
                             if isinstance(value, (dict, list, bool)) else value for key, value in row.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assessment", type=Path, required=True)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--sample-id", required=True, help="expected sample identifier, independently supplied by the caller")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    a_bytes, r_bytes = args.assessment.read_bytes(), args.raw.read_bytes()
    manifest, targets, links, sequences = build_export(json.loads(a_bytes), json.loads(r_bytes), digest(a_bytes), digest(r_bytes), args.sample_id)
    if args.output.exists() and any(args.output.iterdir()):
        parser.error("output directory must be empty; previous exports are never overwritten")
    args.output.mkdir(parents=True, exist_ok=True)
    write_table(args.output / "validation_targets.tsv", FIELDS, targets)
    write_table(args.output / "flank_associations.tsv", LINK_FIELDS, links)
    fasta = "".join(">" + name + "\n" + sequences[name] + "\n" for name in sorted(sequences))
    (args.output / "physical_flanks.fasta").write_text(fasta, encoding="ascii")
    manifest.update(created_utc=datetime.now(timezone.utc).isoformat(), input_files={"assessment": str(args.assessment.resolve()), "native_raw": str(args.raw.resolve())},
                    output_sha256={name: digest((args.output / name).read_bytes()) for name in ("validation_targets.tsv", "flank_associations.tsv", "physical_flanks.fasta")})
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": manifest["status"], "sample_id": args.sample_id, "targets": manifest["target_count"],
                      "exported_sequences": manifest["exported_sequence_count"], "output": str(args.output)}))


if __name__ == "__main__":
    main()
