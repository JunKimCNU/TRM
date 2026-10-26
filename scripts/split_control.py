#!/usr/bin/env python3
"""Split paired FASTQ by read ID without separating mates or selecting motifs."""

import argparse
import contextlib
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import time
from datetime import datetime, timezone


HASH_DOMAIN = b"trm-control-split-v1\x00"


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def file_info(path):
    path = Path(path).resolve()
    stat = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    after = path.stat()
    if (stat.st_size, stat.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"Input changed during hashing: {path}")
    return {"path": str(path), "bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns, "sha256": digest.hexdigest()}


def normalized_id(header, expected_mate=None):
    """Use the first header token and remove a trailing /1 or /2."""
    if not header.startswith(b"@"): 
        raise ValueError("FASTQ header must start with @")
    tokens = header[1:].split()
    if not tokens:
        raise ValueError("FASTQ record ID is empty")
    identifier = tokens[0]
    mate = None
    if identifier.endswith((b"/1", b"/2")):
        mate = int(identifier[-1:])
        identifier = identifier[:-2]
    if not identifier:
        raise ValueError("FASTQ record ID is empty")
    if len(tokens) > 1:
        match = re.fullmatch(rb"\d+/([12])", tokens[1])
        casava = re.fullmatch(rb"([12]):[YN]:\d+:\S*", tokens[1])
        if match or casava:
            description_mate = int((match or casava).group(1))
            if mate is not None and mate != description_mate:
                raise ValueError("Conflicting mate indicators in FASTQ header")
            mate = description_mate
    if expected_mate is not None and mate is not None and mate != expected_mate:
        raise ValueError(f"Wrong mate indicator: expected R{expected_mate}, found R{mate}")
    return identifier


def records(path, expected_mate):
    """Read strict four-line FASTQ; gzip errors and partial records are fatal."""
    with gzip.open(path, "rb") as handle:
        index = 0
        while True:
            header = handle.readline()
            if not header:
                return
            index += 1
            lines = [header, handle.readline(), handle.readline(), handle.readline()]
            try:
                if not all(lines):
                    raise ValueError("Truncated FASTQ record")
                fields = [line.rstrip(b"\r\n") for line in lines]
                identifier = normalized_id(fields[0], expected_mate)
                if not fields[2].startswith(b"+"):
                    raise ValueError("FASTQ separator must start with +")
                if fields[2][1:] and normalized_id(b"@" + fields[2][1:]) != identifier:
                    raise ValueError("FASTQ separator ID does not match header ID")
                sequence, quality = fields[1], fields[3]
                if not sequence or len(sequence) != len(quality):
                    raise ValueError("Empty sequence or sequence/quality length mismatch")
                if any(base not in b"ACGTRYSWKMBDHVNacgtryswkmbdhvn" for base in sequence):
                    raise ValueError("Sequence contains non-IUPAC DNA characters")
                if any(score < 33 or score > 126 for score in quality):
                    raise ValueError("Quality contains non-printable FASTQ characters")
            except ValueError as exc:
                raise ValueError(f"{path}, record {index}: {exc}") from exc
            yield identifier, b"".join(lines), len(sequence)


def verify_manifest(inputs, manifest_path):
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    if not manifest.get("complete"):
        raise ValueError("Control manifest is not marked complete")
    entries = {str(Path(item["path"]).resolve()): item
               for control in manifest["controls"].values()
               for item in control["files"]}
    for item in inputs:
        expected = entries.get(item["path"])
        if expected is None:
            raise ValueError(f"Input missing from manifest: {item['path']}")
        if item["sha256"] != expected.get("sha256") or item["bytes"] != expected.get("bytes"):
            raise ValueError(f"Input checksum or size differs from manifest: {item['path']}")
    return file_info(manifest_path)


def split_control(r1, r2, output_dir, seed, calibration_fraction=0.5, manifest=None):
    start = time.monotonic()
    if not math.isfinite(calibration_fraction) or not 0 < calibration_fraction < 1:
        raise ValueError("Calibration fraction must be strictly between 0 and 1")
    if Path(r1).resolve() == Path(r2).resolve():
        raise ValueError("R1 and R2 must be distinct files")
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    inputs = [file_info(r1), file_info(r2)]
    manifest_info = verify_manifest(inputs, manifest) if manifest else None
    output_dir.mkdir(parents=True, exist_ok=False)
    threshold = int(calibration_fraction * (1 << 64))
    hash_prefix = HASH_DOMAIN + str(seed).encode("utf-8") + b"\x00"
    provenance = {
        "schema_version": 1, "status": "in_progress", "started_at_utc": utc_now(),
        "inputs": inputs, "verified_control_manifest": manifest_info,
        "seed": str(seed), "calibration_fraction": calibration_fraction,
        "assignment": {
            "hash": "SHA-256", "domain": "trm-control-split-v1",
            "message": "ASCII domain + NUL + UTF-8 seed + NUL + normalized FASTQ ID bytes",
            "decision": "first 8 digest bytes as unsigned big-endian integer < threshold => calibration",
            "threshold": threshold,
            "record_id_normalization": "first whitespace-delimited header token without @ or terminal /1 or /2",
            "pair_preserving": True, "motif_independent": True,
            "duplicate_ids": "Equal normalized IDs always receive the same partition.",
        },
        "partitions": {
            part: {"pairs": 0, "records": 0, "bases": 0, "files": []}
            for part in ("calibration", "heldout")
        },
        "note": "Both partitions are from the same biological library; held-out reads are not an independent negative cohort. No motif or sequence-content selection is applied.",
    }
    provenance_path = output_dir / "provenance.json"

    def save():
        temporary = provenance_path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(provenance, indent=2) + "\n")
        temporary.replace(provenance_path)

    save()
    try:
        with contextlib.ExitStack() as stack:
            writers = {}
            for part in provenance["partitions"]:
                for mate in (1, 2):
                    partial = output_dir / f"{part}_R{mate}.fastq.gz.partial"
                    raw = stack.enter_context(partial.open("xb"))
                    writers[part, mate] = stack.enter_context(
                        gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0))
            iterators = [records(r1, 1), records(r2, 2)]
            while True:
                pair = [next(iterator, None) for iterator in iterators]
                if pair == [None, None]:
                    break
                if None in pair:
                    raise ValueError("R1/R2 record count mismatch")
                first, second = pair
                if first[0] != second[0]:
                    raise ValueError(f"Mate ID mismatch: {first[0]!r} versus {second[0]!r}")
                value = int.from_bytes(hashlib.sha256(hash_prefix + first[0]).digest()[:8], "big")
                part = "calibration" if value < threshold else "heldout"
                stats = provenance["partitions"][part]
                stats["pairs"] += 1
                stats["records"] += 2
                stats["bases"] += first[2] + second[2]
                writers[part, 1].write(first[1])
                writers[part, 2].write(second[1])
        for item in inputs:
            stat = Path(item["path"]).stat()
            if (stat.st_size, stat.st_mtime_ns) != (item["bytes"], item["mtime_ns"]):
                raise ValueError(f"Input changed during splitting: {item['path']}")
        for part, stats in provenance["partitions"].items():
            for mate in (1, 2):
                final = output_dir / f"{part}_R{mate}.fastq.gz"
                final.with_suffix(".gz.partial").replace(final)
                stats["files"].append(file_info(final))
        provenance["total_pairs"] = sum(x["pairs"] for x in provenance["partitions"].values())
        provenance["total_records"] = 2 * provenance["total_pairs"]
        provenance["total_bases"] = sum(x["bases"] for x in provenance["partitions"].values())
        provenance["status"] = "complete"
    except Exception as exc:
        provenance["status"] = "failed"
        provenance["error"] = str(exc)
        raise
    finally:
        provenance["finished_at_utc"] = utc_now()
        provenance["wall_seconds"] = time.monotonic() - start
        save()
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("data/controls/support_split"))
    parser.add_argument("--seed", required=True, help="Explicit seed for deterministic ID hashing")
    parser.add_argument("--calibration-fraction", type=float, default=0.5)
    parser.add_argument("--manifest", type=Path, help="Require original SHA-256 and byte sizes to match this control manifest")
    args = parser.parse_args()
    try:
        result = split_control(args.r1, args.r2, args.output_dir, args.seed,
                               args.calibration_fraction, args.manifest)
    except (OSError, ValueError, EOFError, KeyError) as exc:
        print(f"split_control: {exc}", file=sys.stderr)
        return 1
    print(json.dumps({"status": result["status"], "total_pairs": result["total_pairs"],
                      "partitions": {name: stats["pairs"] for name, stats in result["partitions"].items()},
                      "wall_seconds": result["wall_seconds"],
                      "provenance": str(args.output_dir.resolve() / "provenance.json")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
