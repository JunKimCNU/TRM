#!/usr/bin/env python3
"""Create and hash independent fixtures before observing mixed-engine output."""
from __future__ import annotations
import csv
import hashlib
import json
import random
from datetime import datetime, timezone
from pathlib import Path
from reference_oracle import find_mixed, reverse_complement

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/mixed_validation_20260923"
DATA = ROOT / "data/mixed_validation_20260923"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    settings = OUT / "preregistered_settings.json"
    if settings.exists():
        raise SystemExit("Already frozen; refusing to overwrite preregistration")
    groups = [
        {"id": "sosm_9_11", "units": ["TTGGGTTAC", "ACTTGGGTTAC"]},
        {"id": "plantlike_6_7", "units": ["TTCAGG", "TTTCAGG"]},
        {"id": "unrelated_5_5", "units": ["ACGTC", "GATCA"]},
        {"id": "ambiguous_3_3_6", "units": ["ACG", "TCA", "ACGTCA"]},
        {"id": "boundary_3_4", "units": ["AGT", "CCGA"]},
    ]
    rng = random.Random(20260923)
    fixtures = []

    def add(group, label, seq, category, expected=None):
        truth = find_mixed(seq, group["units"])
        detected = any(row["guaranteed_mixed"] for row in truth)
        if expected is not None:
            assert detected == expected, (group["id"], label, truth)
        fixtures.append({"id": group["id"] + "__" + label, "group": group["id"],
                         "category": category, "sequence": seq, "length": len(seq),
                         "expected_mixed": detected, "oracle_intervals": truth})

    for group in groups:
        a, b = group["units"][:2]
        body = a + b + a + b + b + a + a + b
        ambiguous = len(group["units"]) > 2
        positive = None if ambiguous else True
        add(group, "whole_plus", body, "whole_array", positive)
        add(group, "whole_minus", reverse_complement(body), "whole_array", positive)
        add(group, "right_g", "NNNNN" + body, "expected_polarity_boundary", positive)
        add(group, "left_c", reverse_complement(body) + "NNNNN", "expected_polarity_boundary", positive)
        add(group, "left_g", body + "NNNNN", "opposite_polarity_boundary", positive)
        add(group, "right_c", "NNNNN" + reverse_complement(body), "opposite_polarity_boundary", positive)
        add(group, "internal", "N" * 20 + body + "N" * 20, "internal_sequence_positive", positive)
        add(group, "pure_a", a * 12, "homogeneous_array", False)
        add(group, "pure_b", b * 12, "homogeneous_array", False)
        add(group, "one_b", a * 8 + b, "insufficient_component_copies", False)
        add(group, "separated", a * 6 + "N" + b * 6, "no_gap_bridging", False)
        add(group, "N_interrupt", (a + b) * 4 + "N" + (a + b) * 4,
            "separate_positive_arrays", positive)
        add(group, "wrong_phase_join", a * 5 + b[1:] + b[:1] + a * 5,
            "phase_reset_forbidden", None)
        add(group, "mixed_strands", a * 5 + reverse_complement(b) * 5,
            "orientation_switch_forbidden", False)
        shuffled = list(body)
        rng.shuffle(shuffled)
        add(group, "composition_shuffle", "".join(shuffled), "composition_matched_negative", False)
        for phase in range(1, len(a)):
            add(group, f"left_phase_{phase}", a[phase:] + body, "left_edge_partial", positive)
            add(group, f"right_phase_{phase}", body + a[:phase], "right_edge_partial", positive)
        for index in range(6):
            order = [rng.randrange(2) for _ in range(14)]
            order[:4] = [0, 0, 1, 1]
            rng.shuffle(order)
            sequence = "".join(group["units"][i] for i in order)
            add(group, f"random_order_{index}", sequence, "arbitrary_order", positive)
    small = groups[-1]
    a, b = small["units"]
    add(small, "span_23", a * 5 + b * 2, "minimum_span_below", False)
    add(small, "span_24", a * 4 + b * 3, "minimum_span_equal", True)
    add(small, "span_25", a * 3 + b * 4, "minimum_span_above", True)
    add(small, "partial_does_not_count", a * 8 + b + b[:3], "partial_copy_exclusion", False)
    add(groups[3], "ambiguous_all_composite", "ACGTCA" * 8,
        "possible_but_not_guaranteed_mixed", False)

    DATA.mkdir(exist_ok=True, parents=True)
    (DATA / "definitions.json").write_text(json.dumps(groups, indent=2) + "\n")
    (DATA / "fixtures.json").write_text(json.dumps(fixtures, indent=2) + "\n")
    with (DATA / "fixtures.fasta").open("w") as out:
        for fixture in fixtures:
            out.write(f'>{fixture["id"]}\n{fixture["sequence"]}\n')
    # One group per input prevents accidental cross-family expectations.
    for group in groups:
        (DATA / (group["id"] + ".definitions.json")).write_text(json.dumps([group], indent=2) + "\n")
        with (DATA / (group["id"] + ".fasta")).open("w") as out:
            for fixture in fixtures:
                if fixture["group"] == group["id"]:
                    out.write(f'>{fixture["id"]}\n{fixture["sequence"]}\n')

    with (ROOT / "data/next_holdout_20260923/S16.tsv").open(newline="") as inp:
        rows = list(csv.reader(inp, delimiter="\t"))
    truth_rows = [row for row in rows if row and row[0].startswith("m64033_")]
    assert len(truth_rows) == 87
    with (DATA / "published_87_segments.fasta").open("w") as out:
        for row in truth_rows:
            assert int(row[1]) == len(row[3])
            out.write(f'>{row[0]}\n{row[3]}\n')
    source_paths = [
        ROOT / "data/next_holdout_20260923/S16.tsv",
        ROOT / "results/holdout_osmophilus_20260923/acquisition_truth/published_read_mapping.tsv",
        ROOT / "results/holdout_osmophilus_20260923/acquisition_truth/sequence_unique_full_reads.fastq",
        ROOT / "results/holdout_osmophilus_20260923/truth_validation/per_read_summary.tsv",
        ROOT / "results/holdout_osmophilus_20260923/evaluation/dataset_summary.tsv",
    ]
    frozen = [*sorted(DATA.glob("*")), OUT / "reference_oracle.py", Path(__file__)]
    document = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "before_any_new_mixed_engine_output_inspection",
        "method": "reference_guided_exact_mixed_array_reconstruction_not_de_novo",
        "minimum_span_bp": 24, "minimum_complete_copies_per_component": 2,
        "required_distinct_components": 2,
        "grammar": "Exact supplied same-G-strand units in arbitrary order; both strands; no cyclic rotation reset at joins; no mismatch or gap bridging; record-edge partial units allowed and do not count full copies.",
        "ambiguity": "All exact parses for identical maximal span contribute per-unit count minima/maxima; mixed requires at least two unit minima >=2. Possible-only mixed is retained separately, never promoted by arbitrary tie resolution.",
        "primary_synthetic_gate": "Each record matches independent oracle guaranteed-mixed presence and maximum guaranteed-mixed span; intervals compared where native emits every interval. Pure arrays must not become mixed.",
        "fixtures": {"records": len(fixtures), "groups": groups,
                     "positive_records": sum(row["expected_mixed"] for row in fixtures)},
        "real_data": {"published_segments": 87, "sequence_unique_full_reads": 65,
                      "original_molecule_identity_confirmed": False,
                      "coordinate_truth": "Use saved unique-segment mapping; require native interval overlap with its published segment. Do not claim 65/87 read sensitivity or original ID recovery.",
                      "segment_diagnostic": "All87 published telomeric segments are diagnostic input, not independent molecules; report exact-grammar eligible subgroup separately from all87.",
                      "metrics": ["qualifying mixed interval presence", "max exact mixed span",
                                  "published array length", "frozen previous homogeneous tract max",
                                  "record-edge and orientation evidence", "runtime and peak RSS"]},
        "limitations": ["Synthetic sequence specificity is not biological telomere specificity.",
                        "Exact reconstruction cannot bridge substitutions or indels; a shorter recovered span is not the full biological telomere length.",
                        "No new downloads; matched sequencing modalities and assembly are not independent biological replicates.",
                        "Known-literature definitions are supplied explicitly; recovery does not measure de novo motif discovery."],
        "fixture_hashes": {str(path.relative_to(ROOT)): sha(path) for path in frozen},
        "source_hashes": {str(path.relative_to(ROOT)): sha(path) for path in source_paths if path.exists()},
    }
    settings.write_text(json.dumps(document, indent=2) + "\n")
    print(json.dumps({"settings": str(settings), "sha256": sha(settings),
                      "fixtures": len(fixtures), "positive": document["fixtures"]["positive_records"]}, indent=2))


if __name__ == "__main__":
    main()
