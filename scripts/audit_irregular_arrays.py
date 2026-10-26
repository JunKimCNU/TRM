#!/usr/bin/env python3
"""Replay exported irregular-array traces against their complete real sources.

This audit does not discover arrays or classify biological absence. It verifies
the exported evidence, using an independent FASTA/FASTQ reader and observed
source slices, never generated nucleotide fixtures. Memory is bounded by one
source record and its exported arrays/traces, plus the longest source line.
"""

import argparse
from collections import Counter
import contextlib
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import signal
import stat
import sys
import time


PROGRAM = "audit_irregular_arrays"
BUFFER = 1024 * 1024
DNA = re.compile(rb"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]*\Z")
ACGT = re.compile(rb"[ACGT]+\Z")
QUALITY = re.compile(rb"[!-~]+\Z")
SAFE_ID = re.compile(r"[A-Za-z0-9_.:-]+\Z")
SHA = re.compile(r"[a-f0-9]{64}\Z")
COMPLEMENT = bytes.maketrans(b"ACGT", b"TGCA")


def integer(value, label, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{label}: expected integer >= {minimum}")
    return value


def boolean(value, label):
    if type(value) is not bool:
        raise ValueError(f"{label}: expected boolean")
    return value


def same(actual, expected, label):
    if type(actual) is not type(expected) or actual != expected:
        raise ValueError(f"{label}: reported {actual!r}, independently obtained {expected!r}")


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON key: " + key)
        result[key] = value
    return result


def parse_json(value):
    return json.loads(value, object_pairs_hook=unique_object,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError("nonfinite JSON: " + x)))


def snapshot(path):
    value = Path(path).stat()
    if not stat.S_ISREG(value.st_mode):
        raise ValueError("regular file required: " + str(path))
    return {"device": value.st_dev, "inode": value.st_ino,
            "size_bytes": value.st_size, "mtime_ns": value.st_mtime_ns,
            "ctime_ns": value.st_ctime_ns}


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(BUFFER), b""):
            result.update(chunk)
    return result.hexdigest()


class HashingReader(io.RawIOBase):
    def __init__(self, raw):
        super().__init__()
        self.raw = raw
        self.hash = hashlib.sha256()
        self.bytes_read = 0

    def readable(self):
        return True

    def readinto(self, target):
        count = self.raw.readinto(target)
        if count:
            self.hash.update(memoryview(target)[:count])
            self.bytes_read += count
        return count


class VerifiedStream:
    """Hash on-disk bytes during one complete parsed traversal; check stability."""
    def __init__(self, stack, path, expected=None, compressed=True):
        self.path = Path(path).resolve(strict=True)
        self.before = snapshot(self.path)
        self.expected = expected
        if expected is not None:
            if not isinstance(expected, dict) or not SHA.fullmatch(str(expected.get("sha256", ""))):
                raise ValueError("valid expected SHA256 required: " + str(path))
            same(self.before["size_bytes"], integer(expected.get("size_bytes"), "expected file size"), str(path) + " size")
        self.raw = stack.enter_context(self.path.open("rb", buffering=0))
        opened = os.fstat(self.raw.fileno())
        if (opened.st_dev, opened.st_ino) != (self.before["device"], self.before["inode"]):
            raise ValueError("file changed before opening: " + str(path))
        self.hashing = stack.enter_context(HashingReader(self.raw))
        self.buffered = stack.enter_context(io.BufferedReader(self.hashing, BUFFER))
        self.compression = "gzip" if compressed and self.buffered.peek(2)[:2] == b"\x1f\x8b" else "plain"
        if self.compression == "gzip":
            decoder = stack.enter_context(gzip.GzipFile(fileobj=self.buffered, mode="rb"))
            self.stream = stack.enter_context(io.BufferedReader(decoder, BUFFER))
        else:
            self.stream = self.buffered

    def finish(self):
        if self.stream.read(1):
            raise ValueError("unconsumed content before EOF: " + str(self.path))
        same(snapshot(self.path), self.before, str(self.path) + " file stability")
        same(self.hashing.bytes_read, self.before["size_bytes"], str(self.path) + " complete raw byte count")
        value = self.hashing.hash.hexdigest()
        if self.expected is not None:
            same(value, self.expected["sha256"], str(self.path) + " SHA256")
        return {"path": str(self.path), "sha256": value,
                "size_bytes": self.hashing.bytes_read, "compression": self.compression,
                "validated_eof": True, "stat_before_after": self.before}


def without_eol(line):
    if line.endswith(b"\n"):
        line = line[:-1]
        if line.endswith(b"\r"):
            line = line[:-1]
    return line


def record_id(header):
    if any(x < 32 and x != 9 or x == 127 for x in header):
        raise ValueError("control character in record header")
    values = header.decode("utf-8").split()
    if not values:
        raise ValueError("record header has no identifier")
    return values[0]


def normalize_sequence(line, fasta):
    if fasta:
        line = line.replace(b" ", b"").replace(b"\t", b"")
    if not DNA.fullmatch(line):
        raise ValueError("invalid IUPAC DNA character")
    return line.upper()


def records(stream):
    """Independent wrapped FASTA/FASTQ reader, with native header/space rules."""
    pending = None
    file_format = None
    while True:
        line = pending if pending is not None else stream.readline()
        pending = None
        if not line:
            return
        line = without_eol(line)
        if file_format is None:
            if not line:
                continue
            file_format = line[:1]
            if file_format not in (b">", b"@"):
                raise ValueError("expected FASTA or FASTQ header")
        if not line.startswith(file_format):
            raise ValueError("invalid record header or blank line between FASTQ records")
        identifier = record_id(line[1:])
        parts = []
        if file_format == b">":
            while True:
                raw = stream.readline()
                if not raw:
                    break
                value = without_eol(raw)
                if value.startswith(b">"):
                    pending = raw
                    break
                parts.append(normalize_sequence(value, True))
            sequence = b"".join(parts)
        else:
            while True:
                raw = stream.readline()
                if not raw:
                    raise ValueError("truncated FASTQ sequence/separator")
                value = without_eol(raw)
                if value.startswith(b"+"):
                    if value[1:] and record_id(value[1:]) != identifier:
                        raise ValueError("FASTQ repeated separator ID differs")
                    break
                if not value:
                    raise ValueError("empty FASTQ sequence line")
                parts.append(normalize_sequence(value, False))
            sequence = b"".join(parts)
            quality_length = 0
            while quality_length < len(sequence):
                raw = stream.readline()
                if not raw:
                    raise ValueError("truncated FASTQ quality")
                value = without_eol(raw)
                if not QUALITY.fullmatch(value):
                    raise ValueError("invalid or empty FASTQ quality line")
                quality_length += len(value)
            if quality_length != len(sequence):
                raise ValueError("FASTQ quality/sequence lengths differ")
        if not sequence:
            raise ValueError("empty source sequence")
        yield identifier, sequence


def definitions_from_summary(summary):
    definitions = summary.get("definitions")
    if not isinstance(definitions, list) or not 1 <= len(definitions) <= 8:
        raise ValueError("require 1..8 definition groups")
    result = {}
    for group in definitions:
        name, units = group.get("id"), group.get("units")
        if not isinstance(name, str) or not name.strip() or name in result:
            raise ValueError("invalid or duplicate group ID")
        if not isinstance(units, list) or not 1 <= len(units) <= 8:
            raise ValueError("require 1..8 units per group")
        if any(not isinstance(x, str) or not re.fullmatch("[ACGT]{3,32}", x) for x in units):
            raise ValueError("invalid literal definition unit")
        same(units, sorted(units), "sorted definition units")
        forward = [x.encode("ascii") for x in units]
        families = set()
        for unit in forward:
            if any(len(unit) % k == 0 and unit == unit[:k] * (len(unit) // k)
                   for k in range(1, len(unit))):
                raise ValueError("nonprimitive definition unit")
            # Canonical comparison concerns only the supplied units, not fixtures.
            reverse = unit.translate(COMPLEMENT)[::-1]
            family = min(x[i:] + x[:i] for x in (unit, reverse) for i in range(len(x)))
            if family in families:
                raise ValueError("rotation/RC-equivalent duplicate definition")
            families.add(family)
        result[name] = {"+": forward,
                        "-": [x.translate(COMPLEMENT)[::-1] for x in forward]}
    return result


def literal_hits(sequence, units):
    """Independent earliest-finish interval selection, restricted to this span."""
    count, last_end = 0, 0
    while True:
        ends = [start + len(unit) for unit in units
                if (start := sequence.find(unit, last_end)) >= 0]
        if not ends:
            return count
        count += 1
        last_end = min(ends)


def explanations(observed, units):
    """All valid one-unit explanations of ACTUAL observed source bytes.

    Compare slices directly; do not construct edited nucleotide strings.
    Different permissible edit offsets are distinct ambiguous explanations.
    """
    result = []
    for index, unit in enumerate(units):
        k = len(unit)
        if observed == unit:
            result.append((0, index, "=", 0))
        if k < 5:
            continue
        if len(observed) == k:
            differences = [i for i, (a, b) in enumerate(zip(observed, unit)) if a != b]
            if len(differences) == 1:
                result.append((1, index, "S", differences[0]))
        elif len(observed) == k + 1:
            for offset in range(k + 1):
                if observed[:offset] == unit[:offset] and observed[offset + 1:] == unit[offset:]:
                    result.append((1, index, "I", offset))
        elif len(observed) == k - 1:
            for offset in range(k):
                if observed[:offset] == unit[:offset] and observed[offset:] == unit[offset + 1:]:
                    result.append((1, index, "D", offset))
    return result


def fixed_span_optimum(sequence, units):
    """Independent minimum-edit paths constrained to exactly this output span.

    Counts tied optimal explanations capped at two. This does not rediscover
    maximal intervals or test whether omitted source records contain arrays.
    """
    n = len(sequence)
    lengths = sorted({len(u) + d for u in units for d in ((-1, 0, 1) if len(u) >= 5 else (0,))})
    ring_length = max(lengths) + 1
    costs = [None] * ring_length
    ways = [0] * ring_length
    costs[0], ways[0] = 0, 1
    cache = {}
    for end in range(1, n + 1):
        slot = end % ring_length
        costs[slot], ways[slot] = None, 0
        for length in lengths:
            start = end - length
            previous = start % ring_length
            if start < 0 or costs[previous] is None:
                continue
            observed = sequence[start:end]
            if observed not in cache:
                if len(cache) >= 65_536:
                    cache.clear()
                variants = explanations(observed, units)
                if variants:
                    minimum = min(x[0] for x in variants)
                    cache[observed] = (minimum, min(2, sum(x[0] == minimum for x in variants)))
                else:
                    cache[observed] = None
            matched = cache[observed]
            if matched is None:
                continue
            cost, count = matched
            candidate = costs[previous] + cost
            count = min(2, ways[previous] * count)
            if costs[slot] is None or candidate < costs[slot]:
                costs[slot], ways[slot] = candidate, count
            elif candidate == costs[slot]:
                ways[slot] = min(2, ways[slot] + count)
    return costs[n % ring_length], ways[n % ring_length] > 1


def replay_trace(row, sequence, units, config):
    model = row.get("model")
    if not isinstance(model, dict) or not isinstance(model.get("trace"), list) or not model["trace"]:
        raise ValueError("nonempty model trace required")
    values = Counter(matches=0, substitutions=0, insertions=0, deletions=0,
                     model_bases=0, exact_unit_bases=0,
                     max_consecutive_edited_units=0, max_unmatched_query_run=0)
    unit_counts = [0] * len(units)
    position = row["start0"]
    edited_run = unmatched_run = 0
    for token in model["trace"]:
        start = integer(token.get("start0"), "token.start0")
        end = integer(token.get("end0"), "token.end0", 1)
        same(start, position, "trace contiguity (no free segmentation gaps)")
        if not start < end <= row["end0"]:
            raise ValueError("token outside array coordinates")
        index = integer(token.get("unit_index"), "token.unit_index")
        if index >= len(units):
            raise ValueError("token unit_index outside oriented definition")
        op = token.get("edit")
        offset = integer(token.get("edit_offset"), "token.edit_offset")
        if op not in {"=", "S", "I", "D"}:
            raise ValueError("invalid token edit operation")
        observed = sequence[start:end]
        if not ACGT.fullmatch(observed):
            raise ValueError("ambiguous nucleotide cannot be an edited token")
        possible = explanations(observed, units)
        selected = (0 if op == "=" else 1, index, op, offset)
        if selected not in possible:
            raise ValueError(f"token does not replay against source at [{start},{end})")
        minimum = min(x[0] for x in possible)
        same(selected[0], minimum, "minimum-cost explanation of this token")
        ambiguous = sum(x[0] == minimum for x in possible) > 1
        same(token.get("ambiguous"), ambiguous, "token explanation ambiguity")
        k = len(units[index])
        values["model_bases"] += k
        unit_counts[index] += 1
        values["matches"] += k - (op in {"S", "D"})
        for name, event in (("substitutions", "S"), ("insertions", "I"), ("deletions", "D")):
            values[name] += int(op == event)
        if op == "=":
            values["exact_unit_bases"] += k
            edited_run = 0
        else:
            edited_run += 1
            values["max_consecutive_edited_units"] = max(values["max_consecutive_edited_units"], edited_run)
        for q in range(end - start):
            unmatched_run = unmatched_run + 1 if op in {"S", "I"} and q == offset else 0
            values["max_unmatched_query_run"] = max(values["max_unmatched_query_run"], unmatched_run)
        position = end
    same(position, row["end0"], "trace covers complete array span")
    for key, value in values.items():
        same(model.get(key), value, "model." + key)
    for value in model.get("representative_unit_counts", []):
        integer(value, "representative unit count")
    same(model.get("representative_unit_counts"), unit_counts, "representative unit counts")
    span = sequence[row["start0"]:row["end0"]]
    gate = literal_hits(span, units)
    same(model.get("literal_gate_hits"), gate, "literal nonoverlapping gate count")
    if gate < config["min_exact_anchors"]:
        raise ValueError("array fails configured exact literal gate")
    columns = values["matches"] + values["substitutions"] + values["insertions"] + values["deletions"]
    identity = model.get("identity")
    if (type(identity) not in (int, float) or not math.isfinite(identity) or
            not math.isclose(identity, values["matches"] / columns, rel_tol=1e-12, abs_tol=1e-12)):
        raise ValueError("model identity differs from M/(M+S+I+D)")
    same(values["matches"] + values["substitutions"] + values["insertions"], row["span_bp"], "query-consuming operation count")
    same(values["matches"] + values["substitutions"] + values["deletions"], values["model_bases"], "model-consuming operation count")
    cost, ambiguous = fixed_span_optimum(span, units)
    same(values["substitutions"] + values["insertions"] + values["deletions"], cost, "fixed-span minimum edit cost")
    same(model.get("optimal_parse_ambiguous"), ambiguous, "fixed-span optimal parse ambiguity")
    return {"tokens": len(model["trace"]), "span_bp": len(span), **values}


def json_rows(stream):
    previous = (0, 0)
    for line_number, line in enumerate(stream, 1):
        if not line.strip():
            raise ValueError(f"blank arrays.jsonl line {line_number}")
        row = parse_json(line)
        if not isinstance(row, dict):
            raise ValueError("array row must be a JSON object")
        key = (integer(row.get("source_index"), "source_index", 1),
               integer(row.get("read_ordinal"), "read_ordinal", 1))
        if key < previous:
            raise ValueError("array rows not sorted by source index / record ordinal")
        previous = key
        yield row


def validate_row(row, sequence, identifier, definitions, config, array_index, array_ids):
    array_id = row.get("array_id")
    if not isinstance(array_id, str) or not SAFE_ID.fullmatch(array_id):
        raise ValueError("invalid safe array ID")
    if array_id in array_ids:
        raise ValueError("duplicate array ID within source record")
    array_ids.add(array_id)
    same(array_id, f"s{row['source_index']:03}_r{row['read_ordinal']:012}_a{array_index:06}",
         "globally unique ordinal-derived array ID")
    same(row.get("read_id"), identifier, "source record ID")
    same(row.get("read_length"), len(sequence), "source read length")
    group = row.get("group_id")
    if group not in definitions:
        raise ValueError("undefined array group")
    start = integer(row.get("start0"), "start0")
    end = integer(row.get("end0"), "end0", 1)
    if not start < end <= len(sequence):
        raise ValueError("invalid array interval")
    same(row.get("span_bp"), end - start, "span length")
    if end - start < config["min_span"]:
        raise ValueError("array shorter than configured minimum span")
    if row.get("alignment_strand") not in {"+", "-"}:
        raise ValueError("alignment_strand must select supplied or reverse-complement units")
    if row.get("strand") not in {"+", "-", "ambiguous"}:
        raise ValueError("invalid interpreted strand")
    same(row.get("left_end_distance"), start, "left end distance")
    same(row.get("right_end_distance"), len(sequence) - end, "right end distance")


def expected_anchor(row, side, sequence, record_rows, enabled):
    edge = row["start0"] if side == "left" else row["end0"]
    if not enabled:
        return {"id": "", "start0": edge, "end0": edge, "length": 0,
                "status": "disabled", "stop_reason": "export_not_requested"}
    if side == "left":
        touching = any(x["start0"] < edge <= x["end0"] for x in record_rows)
        hard = edge if touching else max([0] + [x["end0"] for x in record_rows if x["end0"] <= edge])
        start, end = max(edge - 1000, hard), edge
        stop = ("read_boundary" if hard == 0 else "reported_array") if start == hard else "maximum_length"
        for position in range(end - 1, start - 1, -1):
            if sequence[position] not in b"ACGT":
                start, stop = position + 1, "ambiguous_base"
                break
    else:
        touching = any(x["start0"] <= edge < x["end0"] for x in record_rows)
        hard = edge if touching else min([len(sequence)] + [x["start0"] for x in record_rows if x["start0"] >= edge])
        start, end = edge, min(edge + 1000, hard)
        stop = ("read_boundary" if hard == len(sequence) else "reported_array") if end == hard else "maximum_length"
        for position in range(start, end):
            if sequence[position] not in b"ACGT":
                end, stop = position, "ambiguous_base"
                break
    length = end - start
    return {"id": row["array_id"] + ("_L" if side == "left" else "_R") if length >= 200 else "",
            "start0": start, "end0": end, "length": length,
            "status": "exported" if length >= 200 else "below_minimum", "stop_reason": stop}


def verify_anchor(row, side, sequence, record_rows, enabled, anchors, exported_ids):
    anchor = row.get(side + "_anchor")
    if not isinstance(anchor, dict):
        raise ValueError("missing adjacent anchor descriptor")
    expected = expected_anchor(row, side, sequence, record_rows, enabled)
    for key, value in expected.items():
        same(anchor.get(key), value, side + " anchor " + key)
    start = integer(anchor.get("start0"), side + " anchor start")
    end = integer(anchor.get("end0"), side + " anchor end")
    if not 0 <= start <= end <= len(sequence):
        raise ValueError("anchor outside source")
    same(anchor.get("length"), end - start, "anchor length")
    same(end if side == "left" else start,
         row["start0"] if side == "left" else row["end0"], "anchor adjacency")
    if end - start > 1000:
        raise ValueError("anchor exceeds fixed 1000 bp maximum")
    name = anchor.get("id")
    if not isinstance(name, str):
        raise ValueError("anchor ID must be a string")
    if name:
        same(name, row["array_id"] + ("_L" if side == "left" else "_R"), "anchor ID")
        same(anchor.get("status"), "exported", "anchor status")
        if end - start < 200 or not ACGT.fullmatch(sequence[start:end]):
            raise ValueError("exported anchor shorter than 200 bp or contains ambiguity")
        found = next(anchors, None)
        if found is None:
            raise ValueError("missing FASTA anchor")
        same(found[0], name, "FASTA anchor identity/order")
        same(found[1], sequence[start:end], "FASTA anchor source sequence")
        if name in exported_ids:
            raise ValueError("duplicate anchor ID within source record")
        exported_ids.add(name)
        return 1
    if anchor.get("status") == "exported":
        raise ValueError("exported anchor has empty ID")
    return 0


def load_summary(directory):
    path = directory / "summary.json"
    before = snapshot(path)
    blob = path.read_bytes()
    same(snapshot(path), before, "summary stability")
    summary = parse_json(blob)
    for key, expected in (("schema_version", 1), ("program", "scan_irregular_arrays"),
                          ("model_version", "0.1.0"), ("status", "complete"),
                          ("complete_input", True), ("validated_eof", True),
                          ("sampled_prefix", False)):
        same(summary.get(key), expected, "summary." + key)
    if summary.get("export_truncated") is True or summary.get("resource_limit_hit") is True:
        raise ValueError("truncated or resource-capped export is not complete")
    sources = summary.get("inputs")
    if not isinstance(sources, list) or not sources:
        raise ValueError("complete input provenance required")
    for source in sources:
        if not isinstance(source.get("path"), str) or source["path"] == "-":
            raise ValueError("replay requires immutable file-backed source paths")
        integer(source.get("records"), "source records")
        integer(source.get("bases"), "source bases")
    config = summary.get("config")
    if not isinstance(config, dict):
        raise ValueError("model config required")
    integer(config.get("min_span"), "minimum span", 1)
    integer(config.get("min_exact_anchors"), "exact literal gate", 1)
    boolean(summary.get("export_anchors"), "export_anchors")
    outputs = summary.get("outputs")
    if not isinstance(outputs, dict):
        raise ValueError("output fingerprints required")
    for name in ("arrays.jsonl", "anchors.fasta"):
        if not isinstance(outputs.get(name), dict):
            raise ValueError("missing output fingerprint: " + name)
    return summary, {"path": str(path), "sha256": hashlib.sha256(blob).hexdigest(),
                     "size_bytes": len(blob), "stat_before_after": before}


def audit(directory, max_traces=None):
    directory = Path(directory).resolve(strict=True)
    summary, summary_fingerprint = load_summary(directory)
    definitions = definitions_from_summary(summary)
    totals = Counter(records=0, bases=0, arrays=0, anchors=0, traces_replayed=0,
                     arrays_not_replayed=0, records_with_arrays=0,
                     max_model_span_bp=0, arrays_model_span_ge_500=0,
                     arrays_covering_ge_90_percent_of_record=0,
                     arrays_model_span_ge_500_and_covering_ge_90_percent_of_record=0)
    group_counts = {name: Counter(records_with_arrays=0, arrays=0, modeled_span_bp_sum=0,
                   exact_unit_bases_sum=0, max_model_span_bp=0, arrays_model_span_ge_500=0,
                   arrays_covering_ge_90_percent_of_record=0,
                   arrays_model_span_ge_500_and_covering_ge_90_percent_of_record=0,
                   orientation_ambiguous_arrays=0, optimal_parse_ambiguous_arrays=0)
                   for name in definitions}
    evidence_totals = Counter()
    inputs = []
    with contextlib.ExitStack() as stack:
        array_file = VerifiedStream(stack, directory / "arrays.jsonl", summary["outputs"].get("arrays.jsonl"), False)
        anchor_file = VerifiedStream(stack, directory / "anchors.fasta", summary["outputs"].get("anchors.fasta"), False)
        rows = iter(json_rows(array_file.stream))
        next_row = next(rows, None)
        anchors = iter(records(anchor_file.stream))
        for source_index, source in enumerate(summary["inputs"], 1):
            with contextlib.ExitStack() as input_stack:
                started = time.monotonic()
                file = VerifiedStream(input_stack, source["path"], source)
                count = bases = 0
                for count, (identifier, sequence) in enumerate(records(file.stream), 1):
                    bases += len(sequence)
                    key = (source_index, count)
                    record_rows = []
                    while next_row is not None and (next_row["source_index"], next_row["read_ordinal"]) == key:
                        record_rows.append(next_row)
                        next_row = next(rows, None)
                    if next_row is not None and (next_row["source_index"], next_row["read_ordinal"]) < key:
                        raise ValueError("exported source ordinal was not found")
                    ids, anchor_ids = set(), set()
                    # Source ordinal remains the join key even when IDs repeat.
                    for array_index, row in enumerate(record_rows, 1):
                        validate_row(row, sequence, identifier, definitions, summary["config"], array_index, ids)
                    coordinates = [(x["group_id"], x["start0"], x["end0"], x["alignment_strand"]) for x in record_rows]
                    same(coordinates, sorted(coordinates), "within-record array order")
                    seen_groups = set()
                    for row in record_rows:
                        opposite_overlap = any(
                            row["group_id"] == other["group_id"] and
                            row["alignment_strand"] != other["alignment_strand"] and
                            row["start0"] < other["end0"] and other["start0"] < row["end0"]
                            for other in record_rows)
                        same(row["strand"], "ambiguous" if opposite_overlap else row["alignment_strand"], "opposite-orientation overlap ambiguity")
                        for side in ("left", "right"):
                            totals["anchors"] += verify_anchor(row, side, sequence, record_rows, summary["export_anchors"], anchors, anchor_ids)
                        if max_traces is None or totals["traces_replayed"] < max_traces:
                            values = replay_trace(row, sequence, definitions[row["group_id"]][row["alignment_strand"]], summary["config"])
                            for key, value in values.items():
                                if key.startswith("max_"):
                                    evidence_totals[key] = max(evidence_totals[key], value)
                                else:
                                    evidence_totals[key] += value
                            totals["traces_replayed"] += 1
                        else:
                            totals["arrays_not_replayed"] += 1
                        totals["arrays"] += 1
                        group = group_counts[row["group_id"]]
                        seen_groups.add(row["group_id"])
                        group["arrays"] += 1
                        group["modeled_span_bp_sum"] += row["span_bp"]
                        group["exact_unit_bases_sum"] += integer(row["model"].get("exact_unit_bases"), "exact unit bases")
                        group["orientation_ambiguous_arrays"] += int(row["strand"] == "ambiguous")
                        group["optimal_parse_ambiguous_arrays"] += int(boolean(row["model"].get("optimal_parse_ambiguous"), "optimal parse ambiguity"))
                        long_span = row["span_bp"] >= 500
                        full_record = row["span_bp"] * 10 >= len(sequence) * 9
                        for counts in (totals, group):
                            counts["max_model_span_bp"] = max(counts["max_model_span_bp"], row["span_bp"])
                            counts["arrays_model_span_ge_500"] += int(long_span)
                            counts["arrays_covering_ge_90_percent_of_record"] += int(full_record)
                            counts["arrays_model_span_ge_500_and_covering_ge_90_percent_of_record"] += int(long_span and full_record)
                    for name in seen_groups:
                        group_counts[name]["records_with_arrays"] += 1
                    totals["records_with_arrays"] += bool(record_rows)
                    if count % 500_000 == 0:
                        print(f"source {source_index}: {count:,} records, {bases:,} bp; {totals['traces_replayed']:,} traces replayed", file=sys.stderr, flush=True)
                if next_row is not None and next_row["source_index"] == source_index:
                    raise ValueError("array read ordinal exceeds complete input record count")
                same(count, source["records"], "complete source records")
                same(bases, source["bases"], "complete source bases")
                receipt = file.finish()
                receipt.update(source_index=source_index, records=count, bases=bases,
                               elapsed_seconds=time.monotonic() - started)
                inputs.append(receipt)
                totals["records"] += count
                totals["bases"] += bases
        if next_row is not None:
            raise ValueError("array source index exceeds declared source list")
        if next(anchors, None) is not None:
            raise ValueError("unreferenced or out-of-order FASTA anchor")
        outputs = {"arrays.jsonl": array_file.finish(), "anchors.fasta": anchor_file.finish()}
    for key in ("records", "bases", "arrays", "anchors"):
        same(totals[key], summary.get(key), "summary." + key)
    reported_counts = summary.get("counts")
    if not isinstance(reported_counts, dict):
        raise ValueError("export summary aggregate counts missing")
    for key in ("records", "bases", "arrays", "anchors", "records_with_arrays", "max_model_span_bp",
                "arrays_model_span_ge_500", "arrays_covering_ge_90_percent_of_record",
                "arrays_model_span_ge_500_and_covering_ge_90_percent_of_record"):
        same(totals[key], reported_counts.get(key), "summary.counts." + key)
    same(reported_counts.get("groups"), {name: dict(counts) for name, counts in group_counts.items()},
         "summary group array statistics")
    same(reported_counts.get("orientations_considered"), totals["records"] * 2 * len(definitions),
         "orientations considered")
    passed = integer(reported_counts.get("orientations_passing_literal_gate"), "orientations passing gate")
    unscored = integer(reported_counts.get("orientations_unscored_literal_gate"), "orientations unscored")
    same(passed + unscored, reported_counts["orientations_considered"], "orientation gate accounting")
    same(snapshot(summary_fingerprint["path"]), summary_fingerprint["stat_before_after"], "summary unchanged during audit")
    same(digest(summary_fingerprint["path"]), summary_fingerprint["sha256"], "summary final SHA256")
    return {"export_summary": summary_fingerprint, "inputs": inputs, "export_outputs": outputs,
            "counts": dict(totals), "group_counts": {name: dict(counts) for name, counts in group_counts.items()},
            "replayed_evidence_totals": dict(evidence_totals),
            "complete_input": True, "validated_eof": True,
            "all_exported_traces_replayed": totals["arrays_not_replayed"] == 0,
            "trace_selection": "all exported arrays" if max_traces is None else f"first {max_traces} exported arrays in source order",
            "checks": ["complete raw source SHA256 and stable file stats", "complete FASTA/FASTQ EOF and record/base counts",
                       "source ordinal/ID/length and array coordinates", "supplied orientation and reverse-complement literal units",
                       "actual source token replay with no free gaps", "M/S/I/D and query/model span identities",
                       "nonoverlapping literal-unit gate within each replayed span", "token and fixed-span optimal parse ambiguity",
                       "fixed-span minimum edit cost", "representative unit counts and edited/unmatched runs",
                       "opposite-orientation overlap ambiguity", "adjacent anchor source sequence, stopping boundaries and FASTA identity",
                       "aggregate/group array statistics and gate accounting arithmetic",
                       "export and summary SHA256 binding"],
            "limitations": ["Replay is conditional on exported intervals; no rescan of omitted records or test of detector sensitivity.",
                            "Whole-record gate pass/fail counts are checked for accounting consistency, not independently rediscovered on every negative record.",
                            "Fixed-span optimality is checked; full-record maximal-interval discovery is covered separately by numeric kernel review.",
                            "Model edits do not distinguish true unit variation from sequencing error.",
                            "Modeled span is not the count of exact-repeat bases; no biological classification or adequacy promotion.",
                            "Any bounded trace selection is an output-order prefix, not an independent or random validation sample."]}


def write_status(path, payload):
    temporary = path.with_suffix(".json.part")
    with temporary.open("w") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path, help="new audit directory")
    parser.add_argument("--max-traces", type=int, help="optional first-N exported trace replay; full source/hash/EOF and anchor checks still required")
    args = parser.parse_args(argv)
    if args.max_traces is not None and args.max_traces < 1:
        parser.error("--max-traces must be positive")
    if args.output_dir.exists():
        parser.error("output directory must be new")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    status_path = args.output_dir / "audit.json"
    started = time.monotonic()
    record = {"schema_version": 1, "program": PROGRAM, "status": "in_progress",
              "started_utc": datetime.now(timezone.utc).isoformat(),
              "command": [sys.executable, str(Path(__file__).resolve()), *(sys.argv[1:] if argv is None else argv)],
              "audit_source_sha256": digest(__file__), "max_traces": args.max_traces}
    write_status(status_path, record)

    def interrupted(signum, frame):
        raise KeyboardInterrupt(f"signal {signum}")

    old_term = signal.signal(signal.SIGTERM, interrupted)
    try:
        record.update(audit(args.export_dir, args.max_traces))
        record["status"] = "complete"
        record["completed_utc"] = datetime.now(timezone.utc).isoformat()
        record["elapsed_seconds"] = time.monotonic() - started
        write_status(status_path, record)
        print(json.dumps({"status": "complete", "audit": str(status_path), "counts": record["counts"]}, sort_keys=True))
        return 0
    except (Exception, KeyboardInterrupt) as exc:
        record.update(status="failed", error=f"{type(exc).__name__}: {exc}",
                      completed_utc=datetime.now(timezone.utc).isoformat(),
                      elapsed_seconds=time.monotonic() - started)
        write_status(status_path, record)
        print(record["error"], file=sys.stderr)
        return 1
    finally:
        signal.signal(signal.SIGTERM, old_term)


if __name__ == "__main__":
    raise SystemExit(main())
