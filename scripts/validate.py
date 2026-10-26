#!/usr/bin/env python3
"""Deterministic, independent synthetic controls for the trm command-line tool.

Run: python3 scripts/validate.py --binary bin/trm
An optional benchmark measures generated FASTQ scanning, not download time.
Synthetic controls test motif recovery and evidence semantics; they do not
establish that a repeat is biologically telomeric in an unknown organism.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import time


BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def canonical(motif: str) -> str:
    reverse = reverse_complement(motif)
    return min(sequence[i:] + sequence[:i]
               for sequence in (motif, reverse) for i in range(len(motif)))


def random_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choices(BASES, k=length))


def repeat(motif: str, length: int, phase: int = 0) -> str:
    return (motif * ((length + phase) // len(motif) + 1))[phase:phase + length]


def mutate(rng: random.Random, sequence: str, fraction: float) -> str:
    letters = list(sequence)
    for position in rng.sample(range(len(letters)), round(len(letters) * fraction)):
        letters[position] = rng.choice(BASES.replace(letters[position], ""))
    return "".join(letters)


def write_sequences(path: Path, sequences: list[str], fasta: bool = False) -> None:
    with path.open("w", encoding="ascii") as output:
        for index, sequence in enumerate(sequences):
            if fasta:
                output.write(f">control_{index}\n{sequence}\n")
            else:
                output.write(f"@control_{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as output:
        json.dump(value, output, indent=2, sort_keys=True)
        output.write("\n")
    os.replace(temporary, path)


def find_binary(value: str | None) -> str:
    from runtime_binary import resolve_binary
    return str(resolve_binary(value))


def run_scan(binary: str, inputs: list[Path], prefix: Path, threads: int,
             mode: str = "reads", extra: list[str] | None = None) -> tuple[dict, dict]:
    command = [binary, "scan", *map(str, inputs), "--threads", str(threads),
               "--output", str(prefix), "--min-k", "3", "--max-k", "32",
               "--min-tract", "500", "--min-reads", "1", "--mode", mode]
    if extra:
        command.extend(extra)
    started = time.perf_counter()
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    elapsed = time.perf_counter() - started
    prefix.with_suffix(".stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"scan exited {completed.returncode}: {completed.stderr[-2000:]}")
    with Path(str(prefix) + ".json").open(encoding="utf-8") as source:
        result = json.load(source)
    if not isinstance(result.get("summary"), dict) or not isinstance(result.get("candidates"), list):
        raise AssertionError("output must contain summary object and candidates array")
    metrics = {"wall_seconds": elapsed, "command": command,
               "input_bytes": sum(path.stat().st_size for path in inputs)}
    return result, metrics


def matching_candidate(result: dict, motif: str) -> dict:
    target = canonical(motif)
    matches = [candidate for candidate in result["candidates"]
               if canonical(candidate["motif"]) == target]
    if not matches:
        observed = [candidate.get("motif") for candidate in result["candidates"][:12]]
        raise AssertionError(f"missing {motif} (canonical {target}); observed {observed}")
    if len(matches) != 1:
        raise AssertionError(f"rotations/reverse complements were not collapsed for {motif}")
    candidate = matches[0]
    if candidate["motif"] != target:
        raise AssertionError(f"expected lexicographic canonical motif {target}, got {candidate['motif']}")
    if candidate.get("reverse_complement") != reverse_complement(target):
        raise AssertionError(f"wrong reverse complement for canonical motif {target}")
    if candidate.get("period") != len(motif):
        raise AssertionError(f"wrong primitive period for {motif}: {candidate.get('period')}")
    return candidate


def require_evidence_fields(candidate: dict) -> None:
    for name in ("reads", "repeat_bases", "max_tract", "long_tract_reads",
                 "forward_reads", "reverse_reads", "left_end_reads", "right_end_reads", "status"):
        if name not in candidate:
            raise AssertionError(f"candidate lacks evidence field: {name}")


def never_confirmed(result: dict) -> None:
    for candidate in result["candidates"]:
        if candidate.get("status") == "confirmed_telomere":
            raise AssertionError("sequence alone was labeled confirmed_telomere")


def build_cases(root: Path, seed: int) -> list[dict]:
    rng = random.Random(seed)
    cases: list[dict] = []

    def add(name: str, sequences: list[str], motif: str | None = None,
            mode: str = "reads", kind: str = "short") -> None:
        path = root / (name + (".fa" if mode == "assembly" else ".fastq"))
        write_sequences(path, sequences, fasta=mode == "assembly")
        cases.append({"name": name, "inputs": [path], "motif": motif,
                      "mode": mode, "kind": kind, "expected_records": len(sequences),
                      "expected_bases": sum(map(len, sequences))})

    # Distinct established telomeric repeat families plus one invented family.
    for motif in ("TTAGGG", "TTAGG", "TTTAGGG", "TTAGGC", "CTCGGTTATGGG", "ACGACACC"):
        sequences = []
        for index in range(48):
            sequence = repeat(motif, 150, index % len(motif))
            sequences.append(reverse_complement(sequence) if index % 2 else sequence)
        sequences.extend(random_dna(rng, 150) for _ in range(64))
        add("short_" + motif, sequences, motif)

    # Exact rotations and reverse complements must collapse into one family.
    motif = "TTAGGG"
    for index in range(len(motif)):
        rotated = motif[index:] + motif[:index]
        add(f"rotation_{index}", [repeat(rotated, 150)] * 4 +
            [reverse_complement(repeat(rotated, 150))] * 4, motif)

    for error_rate in (0.01, 0.03, 0.05):
        sequences = [mutate(rng, repeat(motif, 150, index % 6), error_rate)
                     for index in range(80)]
        sequences = [reverse_complement(value) if index % 2 else value
                     for index, value in enumerate(sequences)]
        sequences.extend(random_dna(rng, 150) for _ in range(64))
        add(f"substitution_{round(error_rate * 100)}pct", sequences, motif)

    long_tract = repeat(motif, 900)
    add("long_reads", [long_tract, reverse_complement(long_tract)] +
        [random_dna(rng, 1200) for _ in range(8)], motif, kind="long")
    add("terminal_assembly", [reverse_complement(long_tract) + random_dna(rng, 2000) + long_tract],
        motif, mode="assembly", kind="terminal")
    add("internal_str", [random_dna(rng, 2000) + long_tract + random_dna(rng, 2000)],
        motif, mode="assembly", kind="internal")
    add("negative_random", [random_dna(rng, 150) for _ in range(2000)], kind="random")
    add("negative_homopolymer", [base * 900 for base in BASES], kind="excluded_period")
    add("negative_dinucleotide", [repeat(motif, 900) for motif in ("AC", "AG", "AT", "CG")],
        kind="excluded_period")

    # Two files are one sample and the second file must contribute to totals.
    paired = [root / "paired_R1.fastq", root / "paired_R2.fastq"]
    write_sequences(paired[0], [repeat(motif, 150)] * 8)
    write_sequences(paired[1], [reverse_complement(repeat(motif, 150))] * 8)
    cases.append({"name": "paired_inputs", "inputs": paired, "motif": motif,
                  "mode": "reads", "kind": "short", "expected_records": 16,
                  "expected_bases": 2400})
    return cases


def evaluate(case: dict, result: dict) -> dict:
    summary = result["summary"]
    for actual, expected in (("records", "expected_records"), ("bases", "expected_bases")):
        if summary.get(actual) != case[expected]:
            raise AssertionError(f"summary.{actual}={summary.get(actual)}; expected {case[expected]}")
    never_confirmed(result)
    if case["kind"] == "excluded_period":
        if result["candidates"]:
            raise AssertionError("homopolymer/dinucleotide repeat passed the primitive period >= 3 filter")
        return {"candidate_count": 0}
    if case["motif"] is None:
        for candidate in result["candidates"]:
            require_evidence_fields(candidate)
            if candidate["long_tract_reads"] or candidate["max_tract"] >= 500:
                raise AssertionError("random 150-base reads were treated as a >=500-base tract")
        return {"candidate_count": len(result["candidates"])}
    candidate = matching_candidate(result, case["motif"])
    require_evidence_fields(candidate)
    kind = case["kind"]
    if kind == "short":
        if candidate["long_tract_reads"] != 0 or candidate["max_tract"] > 150:
            raise AssertionError("separate short reads were combined into long-tract evidence")
        if "long_tract_supported" in candidate["status"]:
            raise AssertionError("150-base reads were labeled long_tract_supported")
        if candidate["repeat_bases"] <= 500:
            raise AssertionError("positive short-read control should aggregate >500 repeat bases")
        if not candidate["forward_reads"] or not candidate["reverse_reads"]:
            raise AssertionError("both read orientations should support the canonical motif")
    if kind in ("long", "terminal", "internal"):
        if candidate["long_tract_reads"] < 1 or candidate["max_tract"] < 500:
            raise AssertionError("a contiguous 900-base tract was not recognized")
    if kind == "internal":
        if candidate["left_end_reads"] != 0 or candidate["right_end_reads"] != 0:
            raise AssertionError("internal STR was counted as assembly-terminal evidence")
    if kind == "terminal":
        if candidate["left_end_reads"] < 1 or candidate["right_end_reads"] < 1:
            raise AssertionError("assembly control did not recover both terminal tracts")
    return {key: candidate[key] for key in ("motif", "period", "reads", "repeat_bases",
            "max_tract", "long_tract_reads", "forward_reads", "reverse_reads", "status")}


def benchmark(binary: str, root: Path, size_mb: float, threads: int, seed: int) -> dict:
    path = root / "benchmark.fastq"
    rng = random.Random(seed + 1)
    target_bytes = math.ceil(size_mb * 1_000_000)
    records = 0
    # A deterministic, mostly nonrepetitive workload with 0.1% telomeric reads.
    with path.open("w", encoding="ascii") as output:
        while output.tell() < target_bytes:
            sequence = (repeat("TTAGGG", 150) if records % 1000 == 0 else random_dna(rng, 150))
            output.write(f"@benchmark_{records}\n{sequence}\n+\n{'I' * 150}\n")
            records += 1
    result, metrics = run_scan(binary, [path], root / "benchmark", threads)
    if result["summary"].get("records") != records:
        raise AssertionError("benchmark record count mismatch")
    matching_candidate(result, "TTAGGG")
    metrics.update({"records": records, "bases": records * 150, "threads": threads,
                    "input_MB_per_second": metrics["input_bytes"] / 1_000_000 / metrics["wall_seconds"],
                    "sequence_Mbases_per_second": records * 150 / 1_000_000 / metrics["wall_seconds"],
                    "notes": "Single generated uncompressed FASTQ run; wall time includes process startup and output. "
                             "Input may be cached. Excludes fixture generation; not a real WGS benchmark."})
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary", help="Explicit scanner override; default: TRM_BINARY or verified bin/trm")
    parser.add_argument("--output-dir", type=Path, default=Path("results/synthetic"))
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20260920)
    parser.add_argument("--benchmark-mb", type=float, default=0,
                        help="also generate and scan this many decimal MB of FASTQ")
    args = parser.parse_args()
    if args.threads < 1 or not math.isfinite(args.benchmark_mb) or args.benchmark_mb < 0:
        parser.error("threads must be positive and benchmark-mb must be finite and nonnegative")
    try:
        binary = find_binary(args.binary)
        root = args.output_dir.resolve()
        root.mkdir(parents=True, exist_ok=True)
        cases = build_cases(root, args.seed)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    report: dict = {"seed": args.seed, "binary": binary, "threads": args.threads,
                    "scope": "Synthetic recovery/evidence controls; not biological confirmation.", "tests": []}
    for case in cases:
        entry = {"name": case["name"]}
        try:
            result, metrics = run_scan(binary, case["inputs"], root / case["name"],
                                       args.threads, case["mode"],
                                       ["--prior-policy", "discovery"]
                                       if case["motif"] == "ACGACACC" else None)
            entry.update({"passed": True, "evidence": evaluate(case, result), **metrics})
            print(f"PASS {case['name']}", flush=True)
        except (OSError, ValueError, KeyError, TypeError, AssertionError, RuntimeError) as error:
            entry.update({"passed": False, "error": str(error)})
            print(f"FAIL {case['name']}: {error}", file=sys.stderr, flush=True)
        report["tests"].append(entry)
    if args.benchmark_mb:
        try:
            report["benchmark"] = benchmark(binary, root, args.benchmark_mb, args.threads, args.seed)
            print(f"BENCHMARK {report['benchmark']['input_MB_per_second']:.2f} MB/s input; "
                  f"{report['benchmark']['sequence_Mbases_per_second']:.2f} Mbases/s", flush=True)
        except (OSError, ValueError, KeyError, TypeError, AssertionError, RuntimeError) as error:
            report["tests"].append({"name": "benchmark", "passed": False, "error": str(error)})
            print(f"FAIL benchmark: {error}", file=sys.stderr)
    report["passed"] = all(item["passed"] for item in report["tests"])
    atomic_json(root / "validation.json", report)
    passed = sum(item["passed"] for item in report["tests"])
    print(f"{passed}/{len(report['tests'])} passed; report: {root / 'validation.json'}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
