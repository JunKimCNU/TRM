#!/usr/bin/env python3
"""Select nested paired FASTQ subsets by source pair ordinal in one input pass.

This samples records, not independent molecules. Different seeds are technical
resampling of one library, never independent biological replicates. Original
record bytes are retained; no motif or sequence-content selection is performed.
"""

import argparse
import contextlib
import csv
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import signal
import stat
import sys
import time
import zlib

try:  # Direct CLI execution and namespace-module/test imports.
    from split_control import normalized_id
except ModuleNotFoundError as exc:
    if exc.name != "split_control":
        raise
    from scripts.split_control import normalized_id


HASH_DOMAIN = b"trm-pair-ordinal-subsample-v1\x00"
UINT64_SPACE = 1 << 64
MAX_UINT64 = UINT64_SPACE - 1
PROGRESS_PAIRS = 500_000
BUFFER_SIZE = 1024 * 1024
DNA = re.compile(rb"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+")
QUALITY = re.compile(rb"[!-~]+")
SAMPLE_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")
SHA256 = re.compile(r"[0-9a-fA-F]{64}")
MANIFEST_FIELDS = (
    "sample_id", "source_sample_id", "path1", "path2", "seed", "fraction",
    "pairs", "records", "bases", "ordinals_path",
)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def uint64(value, name="seed", positive=False):
    """Validate decimal uint64 input without accepting booleans/floats."""
    if isinstance(value, bool) or not re.fullmatch(r"[0-9]+", str(value)):
        raise ValueError(f"{name} must be an unsigned decimal integer")
    number = int(value)
    if not (1 if positive else 0) <= number <= MAX_UINT64:
        raise ValueError(f"{name} must be in [{1 if positive else 0}, {MAX_UINT64}]")
    return number


def fraction_spec(value):
    """Return exact decimal label and floor(fraction * 2**64).

    Fraction(Decimal) avoids rounding by the active decimal context. Values too
    small to select any uint64 value are rejected rather than silently empty.
    """
    try:
        decimal = Decimal(str(value))
    except (InvalidOperation, ValueError) as exc:
        raise ValueError("fraction must be a finite decimal in (0, 1]") from exc
    if not decimal.is_finite() or not 0 < decimal <= 1:
        raise ValueError("fraction must be a finite decimal in (0, 1]")
    # Avoid giant fixed-point labels and unnecessarily huge integer conversion.
    if decimal.adjusted() < -20 or len(decimal.as_tuple().digits) > 100:
        raise ValueError("fraction is too small or has more than 100 decimal digits")
    exact = Fraction(decimal)
    threshold = exact.numerator * UINT64_SPACE // exact.denominator
    if threshold == 0:
        raise ValueError("fraction produces an empty uint64 selection interval")
    label = format(decimal, "f").rstrip("0").rstrip(".") if decimal < 1 else "1"
    return {"fraction": label, "label": label, "threshold": threshold}


def seed_prefix(seed, domain=HASH_DOMAIN):
    if not isinstance(domain, bytes) or not domain:
        raise ValueError("hash domain must be nonempty bytes")
    return domain + uint64(seed).to_bytes(8, "big")


def ordinal_value(seed, ordinal, domain=HASH_DOMAIN):
    """Uniform first64 digest bits, with fixed-width unsigned big-endian inputs."""
    ordinal = uint64(ordinal, "ordinal", positive=True)
    return int.from_bytes(hashlib.sha256(
        seed_prefix(seed, domain) + ordinal.to_bytes(8, "big")).digest()[:8], "big")


def selected(value, threshold):
    """Selection boundary helper, also used for numeric-only verification."""
    if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value < UINT64_SPACE:
        raise ValueError("value must be a uint64 integer")
    if isinstance(threshold, bool) or not isinstance(threshold, int) or not 0 <= threshold <= UINT64_SPACE:
        raise ValueError("threshold must be an integer in [0, 2**64]")
    return value < threshold


def validate_parameters(sample_id, fractions, seeds, expected_pairs, expected_bases,
                        r1_sha256, r2_sha256):
    if not SAMPLE_ID.fullmatch(sample_id) or sample_id in (".", ".."):
        raise ValueError("sample-id must start with an ASCII letter/digit and contain only letters, digits, _, -, .")
    specs = [fraction_spec(value) for value in fractions]
    if not specs or len({x["fraction"] for x in specs}) != len(specs):
        raise ValueError("fractions must be nonempty and distinct exact decimal values")
    seed_values = [uint64(value) for value in seeds]
    if not seed_values or len(set(seed_values)) != len(seed_values):
        raise ValueError("seeds must be nonempty and distinct uint64 values")
    pairs = uint64(expected_pairs, "expected-pairs", positive=True)
    bases = uint64(expected_bases, "expected-bases", positive=True)
    if bases < 2 * pairs:
        raise ValueError("expected-bases must allow at least one base per mate")
    hashes = [r1_sha256, r2_sha256]
    if any(not isinstance(x, str) or not SHA256.fullmatch(x) for x in hashes):
        raise ValueError("both source SHA256 values must contain exactly 64 hexadecimal characters")
    return sorted(specs, key=lambda x: Decimal(x["fraction"])), seed_values, pairs, bases, [x.lower() for x in hashes]


def stat_info(value):
    return {"device": value.st_dev, "inode": value.st_ino, "bytes": value.st_size,
            "mtime_ns": value.st_mtime_ns, "ctime_ns": value.st_ctime_ns}


def validate_source_paths(r1, r2):
    """Reject identical paths, symlink aliases, and hardlinked input mates."""
    paths = [Path(r1).resolve(strict=True), Path(r2).resolve(strict=True)]
    if os.path.samefile(paths[0], paths[1]):
        raise ValueError("R1 and R2 must be distinct files (including symlink/hardlink identity)")
    snapshots = []
    for path in paths:
        value = path.stat()
        if not stat.S_ISREG(value.st_mode):
            raise ValueError(f"Source must be a regular file: {path}")
        snapshots.append(stat_info(value))
    return paths, snapshots


class HashingReader(io.RawIOBase):
    """Hash original on-disk bytes while the parser consumes this one stream."""

    def __init__(self, raw):
        super().__init__()
        self.raw = raw
        self.digest = hashlib.sha256()
        self.bytes_read = 0

    def readable(self):
        return True

    def readinto(self, target):
        count = self.raw.readinto(target)
        if count:
            self.digest.update(memoryview(target)[:count])
            self.bytes_read += count
        return count


def open_source(stack, path, snapshot):
    raw = stack.enter_context(path.open("rb", buffering=0))
    if stat_info(os.fstat(raw.fileno())) != snapshot:
        raise ValueError(f"Input changed before opening: {path}")
    hashing = stack.enter_context(HashingReader(raw))
    buffered = stack.enter_context(io.BufferedReader(hashing, BUFFER_SIZE))
    compressed = buffered.peek(2)[:2] == b"\x1f\x8b"
    if compressed:
        decoder = stack.enter_context(gzip.GzipFile(fileobj=buffered, mode="rb"))
        reader = stack.enter_context(io.BufferedReader(decoder, BUFFER_SIZE))
    else:
        reader = buffered
    return reader, hashing, raw, "gzip" if compressed else "plain"


def records(reader, path, mate):
    """Strict four-line FASTQ, yielding unchanged record bytes and base counts."""
    ordinal = 0
    while True:
        header = reader.readline()
        if not header:
            return
        ordinal += 1
        sequence = reader.readline()
        plus = reader.readline()
        quality = reader.readline()
        try:
            if not sequence or not plus or not quality:
                raise ValueError("Truncated four-line FASTQ record")
            identifier = normalized_id(without_line_ending(header), mate)
            separator = without_line_ending(plus)
            if not separator.startswith(b"+"):
                raise ValueError("FASTQ separator must start with +")
            if separator[1:] and normalized_id(b"@" + separator[1:], mate) != identifier:
                raise ValueError("FASTQ separator ID does not match header ID")
            bases, scores = without_line_ending(sequence), without_line_ending(quality)
            if not bases or len(bases) != len(scores):
                raise ValueError("Empty sequence or sequence/quality length mismatch")
            if DNA.fullmatch(bases) is None:
                raise ValueError("Sequence contains non-IUPAC DNA characters")
            if QUALITY.fullmatch(scores) is None:
                raise ValueError("Quality contains non-printable FASTQ characters")
        except ValueError as exc:
            raise ValueError(f"{path}, record {ordinal}: {exc}") from exc
        yield identifier, header + sequence + plus + quality, len(bases)


def without_line_ending(line):
    """Remove one LF or CRLF terminator, leaving extra control bytes invalid."""
    if line.endswith(b"\n"):
        line = line[:-1]
        if line.endswith(b"\r"):
            line = line[:-1]
    return line


def atomic_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def hash_file(path):
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(BUFFER_SIZE), b""):
            digest.update(block)
            size += len(block)
    return {"bytes": size, "sha256": digest.hexdigest()}


def subsample(sample_id, r1, r2, output_dir, fractions, seeds, expected_pairs,
              expected_bases, r1_sha256, r2_sha256):
    started = time.monotonic()
    specs, seeds, expected_pairs, expected_bases, hashes = validate_parameters(
        sample_id, fractions, seeds, expected_pairs, expected_bases, r1_sha256, r2_sha256)
    paths, snapshots = validate_source_paths(r1, r2)
    output_dir = Path(output_dir).absolute()
    if output_dir.exists() or output_dir.is_symlink():
        raise FileExistsError(f"Output directory must be new: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=False)
    output_dir = output_dir.resolve()
    provenance_path = output_dir / "provenance.json"
    provenance = {
        "schema_version": 1, "status": "in_progress", "complete_input": False,
        "sample_id": sample_id, "started_at_utc": utc_now(),
        "expected_pairs": expected_pairs, "expected_records": 2 * expected_pairs,
        "expected_bases": expected_bases, "total_pairs": 0, "total_records": 0, "total_bases": 0,
        "inputs": [{"mate": i + 1, "path": str(path), "expected_sha256": hashes[i],
                    "stat_before": snapshots[i]} for i, path in enumerate(paths)],
        "selection": {
            "hash": "SHA256", "domain_hex": HASH_DOMAIN.hex(),
            "message": "domain bytes + seed uint64 big-endian + pair ordinal uint64 big-endian",
            "value": "first 8 digest bytes interpreted as unsigned big-endian integer",
            "decision": "value < floor(exact_decimal_fraction * 2**64)",
            "ordinal_start": 1, "duplicate_ids": "Separate pair ordinals; duplicate IDs are not grouped.",
            "fractions_nested_per_seed": True, "motif_independent": True,
            "unit": "paired FASTQ records, not independent molecules",
            "replicate_scope": "Seeds are technical resampling of the same library, not independent biological replicates.",
        },
        "compression": {"format": "gzip", "level": 1, "mtime": 0, "embedded_filename": "",
                        "zlib_runtime": zlib.ZLIB_RUNTIME_VERSION},
        "validation": {"full_eof": False, "mate_ids_and_labels": False,
                       "expected_counts": False, "source_sha256": False, "source_stat_unchanged": False},
        "partitions": [], "timings": {},
        "scope": "Full single-pass source traversal; unchanged selected record bytes; no biological interpretation or detector analysis.",
    }
    pending_files = []
    try:
        atomic_json(provenance_path, provenance)
        with contextlib.ExitStack() as stack:
            groups = []
            for seed in seeds:
                group = []
                for spec in specs:
                    destination = output_dir / f"seed_{seed}" / f"fraction_{spec['label']}"
                    destination.mkdir(parents=True, exist_ok=False)
                    partition = {"sample_id": f"{sample_id}__seed_{seed}__fraction_{spec['label']}",
                                 "source_sample_id": sample_id, "seed": seed, **spec,
                                 "pairs": 0, "records": 0, "bases": 0, "files": [],
                                 "ordinals": {"encoding": "uint64 big-endian, 1-based; exactly 8 bytes per ordinal",
                                              "count": 0}}
                    writers = []
                    for mate in (1, 2):
                        final = destination / f"R{mate}.fastq.gz"
                        partial = final.with_name(final.name + ".partial")
                        raw = stack.enter_context(partial.open("xb"))
                        writer = stack.enter_context(gzip.GzipFile(
                            filename="", mode="wb", fileobj=raw, compresslevel=1, mtime=0))
                        writers.append(writer)
                        info = {"mate": mate, "path": str(final), "partial_path": str(partial)}
                        partition["files"].append(info)
                        pending_files.append((partial, final, info))
                    final = destination / "selected_ordinals.u64be"
                    partial = final.with_name(final.name + ".partial")
                    ordinals = stack.enter_context(partial.open("xb"))
                    partition["ordinals"].update({"path": str(final), "partial_path": str(partial)})
                    pending_files.append((partial, final, partition["ordinals"]))
                    provenance["partitions"].append(partition)
                    group.append((spec["threshold"], writers[0], writers[1], ordinals, partition))
                groups.append((seed_prefix(seed), group))
            sources = [open_source(stack, path, snapshots[i]) for i, path in enumerate(paths)]
            for i, source in enumerate(sources):
                provenance["inputs"][i]["format"] = source[3]
            atomic_json(provenance_path, provenance)
            iterators = [records(source[0], paths[i], i + 1) for i, source in enumerate(sources)]
            traversal_started = time.monotonic()
            while True:
                first, second = next(iterators[0], None), next(iterators[1], None)
                if first is None and second is None:
                    break
                if first is None or second is None:
                    raise ValueError("R1/R2 record count mismatch")
                ordinal = provenance["total_pairs"] + 1
                if first[0] != second[0]:
                    raise ValueError(f"Mate ID mismatch at pair {ordinal}: {first[0]!r} versus {second[0]!r}")
                if ordinal > MAX_UINT64:
                    raise ValueError("Source pair ordinal exceeds uint64 capacity")
                ordinal_bytes = ordinal.to_bytes(8, "big")
                bases = first[2] + second[2]
                for prefix, group in groups:
                    value = int.from_bytes(hashlib.sha256(prefix + ordinal_bytes).digest()[:8], "big")
                    for threshold, writer1, writer2, ordinals, partition in group:
                        if value < threshold:
                            writer1.write(first[1])
                            writer2.write(second[1])
                            ordinals.write(ordinal_bytes)
                            partition["pairs"] += 1
                            partition["records"] += 2
                            partition["bases"] += bases
                            partition["ordinals"]["count"] += 1
                provenance["total_pairs"] = ordinal
                provenance["total_records"] += 2
                provenance["total_bases"] += bases
                if ordinal % PROGRESS_PAIRS == 0:
                    elapsed = time.monotonic() - traversal_started
                    print(f"{sample_id}: processed {ordinal:,}/{expected_pairs:,} pairs; "
                          f"{elapsed:.1f} s; {ordinal / max(elapsed, 1e-9):,.0f} pairs/s", file=sys.stderr, flush=True)
            provenance["timings"]["traversal_seconds"] = time.monotonic() - traversal_started
            provenance["validation"]["full_eof"] = True
            provenance["validation"]["mate_ids_and_labels"] = True
            # At decompressed EOF, gzip has checked the final CRC/size and consumed
            # all members. No second pass over original input bytes is performed.
            for i, (_, hashing, raw, _) in enumerate(sources):
                item = provenance["inputs"][i]
                item.update({"sha256": hashing.digest.hexdigest(), "bytes_read": hashing.bytes_read,
                             "stat_after": stat_info(paths[i].stat()),
                             "descriptor_stat_after": stat_info(os.fstat(raw.fileno()))})
                if item["stat_after"] != snapshots[i] or item["descriptor_stat_after"] != snapshots[i]:
                    raise ValueError(f"Source changed during traversal: {paths[i]}")
                if hashing.bytes_read != snapshots[i]["bytes"]:
                    raise ValueError(f"Raw source traversal did not consume full byte count: {paths[i]}")
                if item["sha256"] != hashes[i]:
                    raise ValueError(f"Source SHA256 mismatch: {paths[i]}")
            provenance["validation"]["source_sha256"] = True
            provenance["validation"]["source_stat_unchanged"] = True
            if provenance["total_pairs"] != expected_pairs or provenance["total_bases"] != expected_bases:
                raise ValueError(f"Source totals differ: pairs={provenance['total_pairs']} (expected {expected_pairs}), "
                                 f"bases={provenance['total_bases']} (expected {expected_bases})")
            provenance["validation"]["expected_counts"] = True
        # All gzip trailers and file buffers have closed successfully before hashing.
        hashing_started = time.monotonic()
        for partial, _, info in pending_files:
            info.update(hash_file(partial))
        for partition in provenance["partitions"]:
            if partition["ordinals"]["bytes"] != 8 * partition["pairs"]:
                raise ValueError("Internal ordinal-file length mismatch")
        provenance["timings"]["output_hash_seconds"] = time.monotonic() - hashing_started
        for partial, final, info in pending_files:
            partial.replace(final)
            info.pop("partial_path")
        manifest_path = output_dir / "manifest.tsv"
        temporary = manifest_path.with_name(manifest_path.name + ".tmp")
        with temporary.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for partition in provenance["partitions"]:
                row = {key: partition[key] for key in MANIFEST_FIELDS if key in partition}
                row.update(path1=partition["files"][0]["path"], path2=partition["files"][1]["path"],
                           ordinals_path=partition["ordinals"]["path"])
                writer.writerow(row)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(manifest_path)
        provenance["manifest"] = {"path": str(manifest_path), **hash_file(manifest_path)}
        provenance["status"] = "complete"
        provenance["complete_input"] = True
        provenance["finished_at_utc"] = utc_now()
        provenance["wall_seconds"] = time.monotonic() - started
        atomic_json(provenance_path, provenance)
    except BaseException as exc:
        provenance["status"] = "failed"
        provenance["complete_input"] = False
        provenance["error"] = {"type": type(exc).__name__, "message": str(exc)}
        provenance["finished_at_utc"] = utc_now()
        provenance["wall_seconds"] = time.monotonic() - started
        try:
            atomic_json(provenance_path, provenance)
        except OSError as save_error:
            print(f"Could not save failed provenance: {save_error}", file=sys.stderr)
        raise
    return provenance


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path, help="Must not already exist")
    parser.add_argument("--fractions", nargs="+", required=True, help="Exact decimal fractions in (0,1], nested within seed")
    parser.add_argument("--seeds", nargs="+", required=True, help="Distinct unsigned 64-bit integers")
    parser.add_argument("--expected-pairs", required=True)
    parser.add_argument("--expected-bases", required=True, help="Total sequence bases across both source mates")
    parser.add_argument("--r1-sha256", required=True, help="Expected SHA256 of original on-disk file bytes")
    parser.add_argument("--r2-sha256", required=True, help="Expected SHA256 of original on-disk file bytes")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

    def terminate(signum, _frame):
        raise InterruptedError(f"Interrupted by signal {signum}; partial outputs retained")

    previous_termination_handler = signal.signal(signal.SIGTERM, terminate)
    try:
        result = subsample(**vars(args))
    except (OSError, ValueError, EOFError, OverflowError, zlib.error) as exc:
        print(f"subsample_read_pairs: {exc}", file=sys.stderr)
        return 1
    finally:
        signal.signal(signal.SIGTERM, previous_termination_handler)
    print(json.dumps({"status": result["status"], "total_pairs": result["total_pairs"],
                      "total_bases": result["total_bases"], "subsets": len(result["partitions"]),
                      "wall_seconds": result["wall_seconds"],
                      "provenance": str(args.output_dir.resolve() / "provenance.json")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
