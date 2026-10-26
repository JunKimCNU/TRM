#!/usr/bin/env python3
"""Exact-seed opportunity in linear represented FASTA records, without a power model."""

import argparse
import csv
import gzip
import hashlib
import heapq
import json
import os
from pathlib import Path
import re
import resource
import sys
import time

VERSION = "1.0.0"
DNA = frozenset(b"ACGT")
IUPAC = frozenset(b"ACGTRYSWKMBDHVNacgtryswkmbdhvn")
TRANSLATION = bytes.maketrans(b"ACGT", b"TGCA")
COUNT_FIELDS = ("valid_starts", "positive_starts", "end_proxy_only_starts",
                "interior_only_starts", "both_context_starts", "negative_starts")


def primitive(sequence):
    for size in range(1, len(sequence) + 1):
        if len(sequence) % size == 0 and sequence == sequence[:size] * (len(sequence) // size):
            return sequence[:size]
    raise AssertionError("unreachable")


def orbit(sequence):
    reverse = sequence.translate(TRANSLATION)[::-1]
    return sorted({part[offset:] + part[:offset] for part in (sequence, reverse)
                   for offset in range(len(part))})


def make_targets(motifs, seed_bases=24, min_copies=3):
    if seed_bases < 1 or min_copies < 1:
        raise ValueError("seed bases and minimum copies must be positive")
    grouped = {}
    for motif in motifs:
        if not isinstance(motif, str):
            raise ValueError("motifs must be strings")
        sequence = motif.upper().encode("ascii")
        if not 1 <= len(sequence) <= 64 or not set(sequence) <= DNA:
            raise ValueError("motifs require 1..64 A/C/G/T bases")
        unit = primitive(sequence)
        rotations = orbit(unit)
        family = rotations[0].decode()
        copies = max(min_copies, (seed_bases + len(unit) - 1) // len(unit))
        if len(unit) * copies > 100000:
            raise ValueError("exact seed length exceeds 100000 bp safety limit")
        target = grouped.setdefault(family, {
            "family": family, "period_bp": len(unit), "seed_length_bp": len(unit) * copies,
            "requested_motifs": [], "seed_words": [word * copies for word in rotations],
        })
        if motif.upper() not in target["requested_motifs"]:
            target["requested_motifs"].append(motif.upper())
    if not grouped or len(grouped) > 256:
        raise ValueError("require 1..256 unique primitive families")
    return [grouped[key] for key in sorted(grouped)]


def _append(intervals, start, stop, maximum):
    """Append a sorted half-open interval, coalescing overlap and adjacency."""
    if stop <= start:
        return
    if intervals and start <= intervals[-1][1]:
        intervals[-1] = (intervals[-1][0], max(stop, intervals[-1][1]))
    else:
        if len(intervals) >= maximum:
            raise ValueError(f"interval limit {maximum} exceeded")
        intervals.append((start, stop))


def seed_runs(sequence, target, end_window_bp, max_intervals=1000000):
    """Distinct literal seed starts; streaming merge avoids a per-hit position list."""
    def occurrences(word):
        position = sequence.find(word)
        while position >= 0:
            yield position
            position = sequence.find(word, position + 1)

    length, seed = len(sequence), target["seed_length_bp"]
    grouped = {"all": [], "end_proxy": [], "interior": []}
    previous = -1
    for start in heapq.merge(*(occurrences(word) for word in target["seed_words"])):
        if start == previous:
            continue
        previous = start
        context = ("end_proxy" if start + seed <= min(end_window_bp, length)
                   or start >= max(0, length - end_window_bp) else "interior")
        _append(grouped["all"], start, start + 1, max_intervals)
        _append(grouped[context], start, start + 1, max_intervals)
    return grouped


def read_start_intervals(runs, length, read_length, seed_length, max_intervals=1000000):
    result = []
    if read_length < seed_length or read_length > length:
        return result
    for start, stop in runs:
        _append(result, max(0, start + seed_length - read_length),
                min(stop, length - read_length + 1), max_intervals)
    return result


def intersection(first, second, max_intervals=1000000):
    result, i, j = [], 0, 0
    while i < len(first) and j < len(second):
        a, b = first[i], second[j]
        _append(result, max(a[0], b[0]), min(a[1], b[1]), max_intervals)
        if a[1] <= b[1]:
            i += 1
        else:
            j += 1
    return result


def cardinality(intervals):
    return sum(stop - start for start, stop in intervals)


def clean_base_runs(sequence, max_intervals=1000000):
    result = []
    for match in re.finditer(b"[ACGT]+", sequence):
        _append(result, match.start(), match.end(), max_intervals)
    return result


def opportunities(sequence, target, read_length, end_window_bp, max_intervals=1000000,
                  prepared_seed_runs=None, prepared_clean_runs=None):
    """Return interval sets and counts; caller supplies uppercase reference bytes."""
    if read_length < 1 or end_window_bp < 0 or max_intervals < 1:
        raise ValueError("read length/interval limit must be positive; end window nonnegative")
    seeds = prepared_seed_runs if prepared_seed_runs is not None else seed_runs(
        sequence, target, end_window_bp, max_intervals)
    clean = prepared_clean_runs if prepared_clean_runs is not None else clean_base_runs(
        sequence, max_intervals)
    length = len(sequence)
    clean_starts = [(a, b - read_length + 1) for a, b in clean if b - a >= read_length]
    sets = {key: read_start_intervals(value, length, read_length, target["seed_length_bp"],
                                     max_intervals) for key, value in seeds.items()}
    clean_sets = {key: intersection(value, clean_starts, max_intervals)
                  for key, value in sets.items()}

    def counts(selected, denominator):
        total, end, internal = (cardinality(selected[key]) for key in ("all", "end_proxy", "interior"))
        both = end + internal - total
        if not 0 <= both <= min(end, internal) or total > denominator:
            raise AssertionError("invalid opportunity partition")
        return dict(zip(COUNT_FIELDS, (denominator, total, total - internal,
                                     total - end, both, denominator - total)))

    return {"all_valid": counts(sets, max(0, length - read_length + 1)),
            "ambiguity_free": counts(clean_sets, cardinality(clean_starts)),
            "intervals": sets, "ambiguity_free_intervals": clean_sets, "seed_runs": seeds}


def fasta_records(handle, max_record_bases):
    name, description, sequence = None, None, bytearray()
    while True:
        line = handle.readline(max_record_bases + 2)
        if not line:
            break
        if len(line) > max_record_bases + 1:
            raise ValueError("FASTA line exceeds configured record safety limit")
        line = line.rstrip(b"\r\n")
        if not line:
            continue
        if line.startswith(b">"):
            if name is not None:
                if not sequence:
                    raise ValueError(f"empty FASTA record {name}")
                yield name, description, bytes(sequence)
            description = line[1:].decode("utf-8")
            fields = description.split()
            if not fields:
                raise ValueError("empty FASTA identifier")
            name, sequence = fields[0], bytearray()
        else:
            if name is None:
                raise ValueError("sequence before first FASTA header")
            if not set(line) <= IUPAC:
                raise ValueError(f"non-IUPAC sequence characters in {name}")
            if len(sequence) + len(line) > max_record_bases:
                raise ValueError(f"record {name} exceeds max-record-bases={max_record_bases}")
            sequence.extend(line)
    if name is None:
        raise ValueError("empty FASTA input")
    if not sequence:
        raise ValueError(f"empty FASTA record {name}")
    yield name, description, bytes(sequence)


def _ratio_fields(row):
    denominator = row["valid_starts"]
    row["q_reference_exact_seed"] = row["positive_starts"] / denominator if denominator else None
    return row


def _write_json(path, value):
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def run(args):
    started = time.monotonic()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=False)
    try:
        motifs = list(args.motif or [])
        if args.targets_json:
            supplied = json.loads(Path(args.targets_json).read_text())
            if isinstance(supplied, list):
                motifs.extend(supplied)
            elif isinstance(supplied, dict) and "motifs" in supplied:
                if not isinstance(supplied["motifs"], list):
                    raise ValueError("targets motifs must be an array")
                motifs.extend(supplied["motifs"])
            elif isinstance(supplied, dict) and "candidates" in supplied:
                if not isinstance(supplied["candidates"], list):
                    raise ValueError("targets candidates must be an array")
                motifs.extend(item["motif"] for item in supplied["candidates"])
            else:
                raise ValueError("targets JSON requires a motif list or motifs/candidates array")
        targets = make_targets(motifs, args.seed_bases, args.min_copies)
        lengths = sorted(set(args.read_length))
        if len(lengths) > 64 or any(r < 1 for r in lengths):
            raise ValueError("require 1..64 distinct positive read lengths")
        if args.end_window_bp < 0 or args.max_record_bases < 1 or args.max_intervals < 1:
            raise ValueError("invalid coordinate/resource limit")
        source = Path(args.fasta).resolve(strict=True)
        fields = ["record_id", "record_length_bp", "family", "period_bp", "seed_length_bp",
                  "read_length_bp", "denominator_scope", *COUNT_FIELDS, "q_reference_exact_seed"]
        interval_fields = ["record_id", "family", "read_length_bp", "denominator_scope", "context",
                           "start_inclusive", "stop_exclusive"]
        seed_fields = ["record_id", "family", "seed_length_bp", "context", "start_inclusive", "stop_exclusive"]
        record_fields = ["record_id", "description", "length_bp", "sequence_sha256_original_case",
                         "ambiguous_bases", "lowercase_bases", "reference_compartment"]
        aggregates, record_ids, total_bases, total_records = {}, set(), 0, 0
        handles = {}
        try:
            writers = {}
            for name, columns in (("opportunities.tsv", fields), ("read_start_intervals.tsv", interval_fields),
                                  ("seed_start_runs.tsv", seed_fields), ("records.tsv", record_fields)):
                handles[name] = (out / (name + ".part")).open("x", newline="")
                writers[name] = csv.DictWriter(handles[name], columns, delimiter="\t")
                writers[name].writeheader()
            with source.open("rb") as raw:
                before = os.fstat(raw.fileno())
                digest = hashlib.sha256()
                while chunk := raw.read(1 << 20):
                    digest.update(chunk)
                raw.seek(0)
                compressed = raw.read(2) == b"\x1f\x8b"
                raw.seek(0)
                stream = gzip.GzipFile(fileobj=raw) if compressed else raw
                for name, description, original in fasta_records(stream, args.max_record_bases):
                    if name in record_ids:
                        raise ValueError(f"duplicate FASTA identifier {name}")
                    record_ids.add(name)
                    sequence = original.upper()
                    total_records += 1
                    total_bases += len(sequence)
                    writers["records.tsv"].writerow({
                        "record_id": name, "description": description, "length_bp": len(sequence),
                        "sequence_sha256_original_case": hashlib.sha256(original).hexdigest(),
                        "ambiguous_bases": len(sequence) - sum(sequence.count(bytes([b])) for b in DNA),
                        "lowercase_bases": sum(97 <= b <= 122 for b in original),
                        "reference_compartment": "unassigned_preserved_record",
                    })
                    clean = clean_base_runs(sequence, args.max_intervals)
                    for target in targets:
                        seeds = seed_runs(sequence, target, args.end_window_bp, args.max_intervals)
                        for context, intervals in seeds.items():
                            for a, b in intervals:
                                writers["seed_start_runs.tsv"].writerow({
                                    "record_id": name, "family": target["family"],
                                    "seed_length_bp": target["seed_length_bp"], "context": context,
                                    "start_inclusive": a, "stop_exclusive": b})
                        for read_length in lengths:
                            result = opportunities(sequence, target, read_length, args.end_window_bp,
                                                   args.max_intervals, seeds, clean)
                            for scope, interval_key in (("all_valid", "intervals"),
                                                        ("ambiguity_free", "ambiguity_free_intervals")):
                                row = {"record_id": name, "record_length_bp": len(sequence),
                                       **{key: target[key] for key in ("family", "period_bp", "seed_length_bp")},
                                       "read_length_bp": read_length, "denominator_scope": scope,
                                       **_ratio_fields(result[scope].copy())}
                                writers["opportunities.tsv"].writerow(row)
                                key = target["family"], read_length, scope
                                accumulator = aggregates.setdefault(key, {name: 0 for name in COUNT_FIELDS})
                                for count in COUNT_FIELDS:
                                    accumulator[count] += row[count]
                                for context, intervals in result[interval_key].items():
                                    for a, b in intervals:
                                        writers["read_start_intervals.tsv"].writerow({
                                            "record_id": name, "family": target["family"],
                                            "read_length_bp": read_length, "denominator_scope": scope,
                                            "context": context, "start_inclusive": a, "stop_exclusive": b})
                after = os.fstat(raw.fileno())
                if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
                    raise ValueError("source changed during hashing/parsing")
                if compressed:
                    stream.close()
        finally:
            for handle in handles.values():
                handle.close()
        for name in handles:
            os.replace(out / (name + ".part"), out / name)
        summary = {
            "schema_version": 1, "version": VERSION, "status": "complete", "full_fasta_eof_verified": True,
            "source": {"path": str(source), "size_bytes": before.st_size, "sha256": digest.hexdigest(),
                       "compression": "gzip" if compressed else "plain", "records": total_records, "bases": total_bases},
            "configuration": {"read_lengths_bp": lengths, "seed_bases": args.seed_bases,
                              "min_copies": args.min_copies, "end_window_bp": args.end_window_bp,
                              "max_record_bases": args.max_record_bases, "max_intervals_per_list": args.max_intervals},
            "targets": [{key: value for key, value in item.items() if key != "seed_words"} for item in targets],
            "aggregate_reference_opportunities": [
                {"family": key[0], "read_length_bp": key[1], "denominator_scope": key[2],
                 **_ratio_fields(value)} for key, value in sorted(aggregates.items())],
            "scope": {
                "reference_geometry": "Linear represented records; no wraparound or cross-record windows.",
                "record_weights": "Every legal start across every supplied record counted once; no biological copy-number weights.",
                "compartments": "No inferred nuclear/organellar labels; all records preserved individually for external annotation.",
                "seed_context": "A seed is end_proxy only when wholly within either fixed end window; other seeds are interior.",
                "window_partition": "end_proxy_only/interior_only/both describe seed contexts contained by each distinct read-start window.",
                "orientation": "Rotation/reverse-complement family union, one coordinate per window; no mate or independent-molecule claim.",
                "q": "Conditional represented-reference exact-seed opportunity fraction, not detector sensitivity, library retention or nondetection power.",
                "zero": "Computed zero exact-seed opportunity under stated geometry, never biological absence.",
            },
            "resources": {"elapsed_seconds": time.monotonic() - started,
                          "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
            "outputs": {},
        }
        for name in handles:
            path = out / name
            checksum = hashlib.sha256()
            with path.open("rb") as handle:
                while chunk := handle.read(1 << 20):
                    checksum.update(chunk)
            summary["outputs"][name] = {"sha256": checksum.hexdigest(), "size_bytes": path.stat().st_size}
        summary["resources"] = {"elapsed_seconds": time.monotonic() - started,
                                "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        _write_json(out / "summary.json.part", summary)
        os.replace(out / "summary.json.part", out / "summary.json")
        return summary
    except Exception as exc:
        _write_json(out / "failure.json", {"status": "failed", "error_type": type(exc).__name__,
                                           "error": str(exc), "elapsed_seconds": time.monotonic() - started,
                                           "biological_result_available": False})
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--motif", action="append")
    parser.add_argument("--targets-json")
    parser.add_argument("--read-length", type=int, action="append", required=True)
    parser.add_argument("--end-window-bp", type=int, default=5000)
    parser.add_argument("--seed-bases", type=int, default=24)
    parser.add_argument("--min-copies", type=int, default=3)
    parser.add_argument("--max-record-bases", type=int, default=250000000)
    parser.add_argument("--max-intervals", type=int, default=1000000)
    parser.add_argument("--outdir", required=True, help="New directory; existing results are never overwritten")
    args = parser.parse_args(argv)
    try:
        result = run(args)
    except Exception as exc:
        print(f"reference_seed_opportunity: {exc}", file=sys.stderr)
        return 1
    print(json.dumps({"status": "complete", "records": result["source"]["records"],
                      "bases": result["source"]["bases"], "elapsed_seconds": result["resources"]["elapsed_seconds"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
