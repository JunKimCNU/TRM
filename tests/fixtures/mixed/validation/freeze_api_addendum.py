#!/usr/bin/env python3
"""API clarification frozen before any new native output is inspected."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from reference_oracle import find_mixed, reverse_complement

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = ROOT / "data/mixed_validation_20260923"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    destination = OUT / "api_validation_addendum.json"
    if destination.exists():
        raise SystemExit("Already frozen; refusing overwrite")
    group = {"id": "noncomposite_ambiguous_4_6_8", "units": ["ACGT", "ACGTAC", "GTACACGT"]}
    a, b, c = group["units"]
    ambiguous = b + c
    assert ambiguous == a + b + a
    sequences = [
        ("ambiguous_tandem", ambiguous * 4),
        ("ambiguous_tandem_rc", reverse_complement(ambiguous * 4)),
        ("ambiguous_internal", "NNNNN" + ambiguous * 4 + "NNNNN"),
        ("guaranteed_components", a * 5 + ambiguous * 4 + c * 5),
        ("guaranteed_components_rc", reverse_complement(a * 5 + ambiguous * 4 + c * 5)),
        ("left_clipped", b[1:] + ambiguous * 4),
        ("right_clipped", ambiguous * 4 + c[:-1]),
        ("both_clipped", b[1:] + ambiguous * 4 + c[:-1]),
    ]
    fixtures = []
    for label, sequence in sequences:
        truth = find_mixed(sequence, group["units"])
        fixtures.append({"id": group["id"] + "__" + label, "group": group["id"],
                         "category": "noncomposite_ambiguous_grammar", "sequence": sequence,
                         "length": len(sequence), "oracle_intervals": truth,
                         "expected_mixed": any(row["guaranteed_mixed"] for row in truth)})
    defpath = DATA / (group["id"] + ".definitions.json")
    fastapath = DATA / (group["id"] + ".fasta")
    truthpath = DATA / "api_addendum_fixtures.json"
    defpath.write_text(json.dumps([group], indent=2) + "\n")
    fastapath.write_text("".join(f'>{row["id"]}\n{row["sequence"]}\n' for row in fixtures))
    truthpath.write_text(json.dumps(fixtures, indent=2) + "\n")
    doc = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "before_any_new_mixed_engine_output_inspection",
        "reason": "Root confirmed an originally intended atomic-unit API constraint omitted in initial schema clarification: supplied units that are exact concatenations of other supplied units are rejected.",
        "original_preregistration_preserved": True,
        "original_preregistration_sha256": sha(OUT / "preregistered_settings.json"),
        "original_fixture_group_reclassified_as_invalid_definition": "ambiguous_3_3_6",
        "invalid_definition_records": 26,
        "remaining_original_valid_sequence_fixtures": 126,
        "additional_noncomposite_ambiguous_sequence_fixtures": len(fixtures),
        "new_group": group,
        "ambiguity_witness": {"sequence": ambiguous, "parse_1": [b, c], "parse_2": [a, b, a]},
        "definition_constraints": ["Every unit is primitive", "No rotational or reverse-complement family duplicate", "No unit is exactly representable as concatenation of the others"],
        "thresholds_changed": False,
        "hashes": {str(path.relative_to(ROOT)): sha(path) for path in [defpath, fastapath, truthpath, Path(__file__)]},
    }
    destination.write_text(json.dumps(doc, indent=2) + "\n")
    print(json.dumps({"addendum": str(destination), "sha256": sha(destination),
                      "valid_sequence_fixtures_total": 134}, indent=2))


if __name__ == "__main__":
    main()
