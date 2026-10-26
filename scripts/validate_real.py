#!/usr/bin/env python3
"""Check pinned real WGS controls and an optional C. elegans reference assembly.

Checksums are verified outside scan timing. These are biological smoke controls,
not a cross-taxon sensitivity/specificity benchmark or proof of telomere absence.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
from fetch_controls import CONTROLS


def digest(path, algorithm="sha256"):
    result = hashlib.new(algorithm)
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def main():
    from runtime_binary import resolve_binary
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", help="Explicit scanner override; default: TRM_BINARY or verified bin/trm")
    parser.add_argument("--data", type=Path, default=Path("data/controls"))
    parser.add_argument("--output", type=Path, default=Path("results/real"))
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    try:
        args.binary = resolve_binary(args.binary)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    args.output.mkdir(parents=True, exist_ok=True)
    cases = []
    for name in ("arabidopsis", "drosophila"):
        files = []
        for entry in CONTROLS[name]["files"]:
            path = args.data / entry["name"]
            if path.stat().st_size != entry["bytes"] or digest(path, "md5") != entry["md5"]:
                raise ValueError("size/MD5 mismatch: " + str(path))
            files.append(path)
        cases.append((name, files, []))
    ce11 = args.data / "ce11.chrI.fa.gz"
    if ce11.exists():
        if digest(ce11, "md5") != "40ce7b39a878f3d6a77bc99f7e2175be":
            raise ValueError("ce11 MD5 mismatch")
        cases.append(("ce11_chrI", [ce11], ["--mode", "assembly", "--min-reads", "1"]))
    report = {"binary_sha256": digest(args.binary), "threads": args.threads,
              "scope": "Two low-coverage WGS runs plus optional reference chromosome; not tree-of-life validation or a specificity estimate. The Drosophila checks cover four classical families, not absence of all eligible off-target repeats.",
              "cases": []}
    for name, files, extra in cases:
        prefix = args.output / name
        command = [str(args.binary.resolve()), "scan", *map(str, files),
                   "--threads", str(args.threads), "--output", str(prefix), *extra]
        started = time.perf_counter()
        result = subprocess.run(command, text=True, capture_output=True, check=False)
        elapsed = time.perf_counter() - started
        prefix.with_suffix(".stderr.log").write_text(result.stderr)
        if result.returncode:
            raise RuntimeError(result.stderr)
        data = json.loads(prefix.with_suffix(".json").read_text())
        candidates = {c["motif"]: c for c in data["candidates"]}
        checks = {}
        evidence = {}
        if name == "arabidopsis":
            checks["expected_motif_TTTAGGG_present"] = "AAACCCT" in candidates
            checks["expected_motif_TTTAGGG_ranked_first"] = bool(data["candidates"]) and data["candidates"][0]["motif"] == "AAACCCT"
            checks["expected_base_count"] = data["summary"]["bases"] == 88411707
            if "AAACCCT" in candidates:
                evidence = candidates["AAACCCT"]
                evidence = {**evidence, "rank": data["candidates"].index(candidates["AAACCCT"]) + 1}
        elif name == "drosophila":
            classical = sorted(set(candidates) & {"AACCCT", "AAACCCT", "AACCT", "AAGCCT"})
            checks["classical_telomeric_families_not_reported_in_this_run"] = not classical
            checks["expected_base_count"] = data["summary"]["bases"] == 52858552
            evidence = {"top_prior_candidate": data["candidates"][0] if data["candidates"] else None,
                        "classical_telomeric_families": classical,
                        "eligible_candidate_count": len(candidates),
                        "rejected_candidate_count": len(data.get("rejected_candidates", [])),
                        "interpretation": "Global sequence priors can admit off-target motifs in this control, including motifs known from other taxa. Passing these smoke checks does not establish specificity, motif absence, or ALT."}
        else:
            c = candidates.get("AAGCCT", {})
            checks["expected_motif_TTAGGC_present"] = bool(c)
            checks["expected_motif_TTAGGC_ranked_first"] = bool(data["candidates"]) and data["candidates"][0]["motif"] == "AAGCCT"
            checks["both_ends_opposite_strands"] = c.get("opposed_end_records", 0) > 0
            checks["short_terminal_tract_not_promoted_to_500bp"] = c.get("max_tract", 500) < 500 and c.get("long_tract_reads", 1) == 0
            evidence = c
        if name != "ce11_chrI":
            checks["short_reads_never_claim_500bp"] = all(c["long_tract_reads"] == 0 for c in candidates.values())
        entry = {"name": name, "command": command, "wall_seconds": elapsed,
                 "checks": checks, "passed": all(checks.values()), "summary": data["summary"], "evidence": evidence}
        report["cases"].append(entry)
        print(f"{'PASS' if entry['passed'] else 'FAIL'} {name}: {elapsed:.3f}s; {len(candidates)} candidate families")
    report["passed"] = all(c["passed"] for c in report["cases"])
    target = args.output / "validation.json"
    temporary = target.with_suffix(".json.part")
    temporary.write_text(json.dumps(report, indent=2) + "\n")
    temporary.replace(target)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
