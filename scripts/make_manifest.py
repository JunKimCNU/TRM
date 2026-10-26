#!/usr/bin/env python3
"""Make one batch manifest row per accession; never pool different accessions.

Pair accession_1.fastq[.gz] with accession_2.fastq[.gz]. Unpaired files remain
single inputs and are listed on stderr. This groups filenames, not species.
"""
import argparse
import csv
from pathlib import Path
import re
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.directory.is_dir():
        parser.error("directory does not exist")
    files = sorted(p for p in args.directory.iterdir() if p.is_file()
                   and re.search(r"\.(?:fastq|fq)(?:\.gz)?$", p.name))
    if not files:
        parser.error("no FASTQ files found")
    rows, consumed, identifiers = [], set(), set()
    unpaired = []
    for path in files:
        if path in consumed:
            continue
        match = re.fullmatch(r"(.+?)(?:_([12]))?(\.(?:fastq|fq)(?:\.gz)?)", path.name)
        name, mate, suffix = match.groups()
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", name):
            parser.error(f"filename gives unsupported sample_id: {path.name}")
        if name in identifiers:
            parser.error(f"ambiguous repeated accession {name}; review files before grouping")
        identifiers.add(name)
        paths = [str(path.resolve()), ""]
        if mate == "1" and path.with_name(name + "_2" + suffix).is_file():
            second = path.with_name(name + "_2" + suffix)
            consumed.add(second)
            paths[1] = str(second.resolve())
        else:
            unpaired.append(path.name)
        rows.append([name, *paths])
        consumed.add(path)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + ".part")
    with temporary.open("w", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "path1", "path2"])
        writer.writerows(rows)
    temporary.replace(args.output)
    print(f"{len(rows)} samples: {len(rows)-len(unpaired)} paired, {len(unpaired)} single input")
    if unpaired:
        print("Single inputs (verify library/mate availability): " + ", ".join(unpaired), file=sys.stderr)


if __name__ == "__main__":
    main()
