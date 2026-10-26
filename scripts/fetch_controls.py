#!/usr/bin/env python3
"""Fetch public WGS controls and an optional assembly with Python's standard library.

Each complete gzip sequence file is checked against its archive byte count and MD5 before it
is atomically renamed. Existing valid downloads are reused; invalid existing
files are preserved and reported. No FASTQ is sampled or motif-enriched.
"""

import argparse
import copy
import datetime
import hashlib
import http.client
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import time
import urllib.error
import urllib.request


DEFAULT_CONTROLS = ("arabidopsis", "drosophila")

CONTROLS = {
    "ce11-chrI": {
        "assembly_accession": "GCF_000002985.6",
        "assembly_name": "WBcel235 / UCSC ce11",
        "scientific_name": "Caenorhabditis elegans",
        "role": "positive_assembly_chromosome",
        "expected_motifs": ["TTAGGC"],
        "data_type": "reference_assembly_FASTA_not_raw_WGS",
        "provenance_url": "https://hgdownload.soe.ucsc.edu/goldenPath/ce11/chromosomes/",
        "checksum_source_url": "https://hgdownload.soe.ucsc.edu/goldenPath/ce11/chromosomes/md5sum.txt",
        "description": "Chromosome I from the WBcel235 reference assembly served by UCSC as ce11.",
        "caveats": [
            "This is assembled chromosome sequence, not raw WGS reads; it validates motif and chromosome-end location only in this reference.",
            "Exact terminal repeat runs in this chromosome can be shorter than 500 bp; failure of the 500 bp gate must not erase weaker telomeric evidence.",
            "Included only when explicitly selected with --control ce11-chrI; default downloads remain the two paired WGS runs."
        ],
        "files": [
            {"name": "ce11.chrI.fa.gz", "bytes": 4772981,
             "md5": "40ce7b39a878f3d6a77bc99f7e2175be",
             "url": "https://hgdownload.soe.ucsc.edu/goldenPath/ce11/chromosomes/chrI.fa.gz"}
        ]
    },
    "arabidopsis": {
        "run_accession": "ERR12378605",
        "study_accession": "PRJEB71401",
        "experiment_accession": "ERX11755252",
        "scientific_name": "Arabidopsis thaliana",
        "role": "positive",
        "expected_motifs": ["TTTAGGG"],
        "library_strategy": "WGS",
        "library_source": "GENOMIC",
        "library_selection": "RANDOM",
        "instrument_model": "Illumina HiSeq 3000",
        "ena_read_count": 588296,
        "ena_base_count": 88411707,
        "provenance_url": "https://www.ebi.ac.uk/ena/browser/api/xml/PRJEB71401",
        "description": "Whole-genome resequencing of F2 Col x Ler chromosome-substitution populations, including wild type and recq4ab figl1 mutants.",
        "caveats": [
            "Low-coverage biological positive control; detection depends on the number and quality of sampled telomeric reads.",
            "Motifs are equivalent under cyclic rotation and reverse complement.",
            "Short-read support does not establish a contiguous 500 bp repeat tract or chromosome-end localization."
        ],
        "files": [
            {"name": "ERR12378605_1.fastq.gz", "bytes": 23794176,
             "md5": "0e10916b2863394f6bf1067d5a16b1fa",
             "url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/005/ERR12378605/ERR12378605_1.fastq.gz"},
            {"name": "ERR12378605_2.fastq.gz", "bytes": 28044007,
             "md5": "cd4e518006febd0a9c0e8704223cab2a",
             "url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/005/ERR12378605/ERR12378605_2.fastq.gz"}
        ]
    },
    "drosophila": {
        "run_accession": "SRR1183635",
        "study_accession": "PRJNA239534",
        "experiment_accession": "SRX481220",
        "scientific_name": "Drosophila melanogaster",
        "role": "negative_for_canonical_short_telomeric_repeat",
        "expected_motifs": [],
        "library_strategy": "WGS",
        "library_source": "GENOMIC",
        "library_selection": "RANDOM",
        "instrument_model": "Illumina HiSeq 2000",
        "ena_read_count": 261676,
        "ena_base_count": 52858552,
        "provenance_url": "https://www.ebi.ac.uk/ena/browser/api/xml/SRX481220",
        "description": "Pooled whole-genome sequencing from adults; mir-9a-one-copy-ancestral2 population.",
        "caveats": [
            "Drosophila maintains telomeres principally through HeT-A, TART and TAHRE retrotransposons, not a canonical short telomerase repeat.",
            "Other short tandem repeats can be abundant; detecting a candidate is not necessarily a false positive unless it is asserted to be telomeric.",
            "A no-call in this low-coverage run does not prove motif absence or an ALT mechanism.",
            "Repeat-positive contamination, interstitial repeats and sequencing artifacts remain possible."
        ],
        "files": [
            {"name": "SRR1183635_1.fastq.gz", "bytes": 18219353,
             "md5": "c824b06e7556bf4f495ee6c798d80c28",
             "url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/005/SRR1183635/SRR1183635_1.fastq.gz"},
            {"name": "SRR1183635_2.fastq.gz", "bytes": 18775150,
             "md5": "b6b1bf51bf602c8dc4edf9a2d96e0af5",
             "url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/005/SRR1183635/SRR1183635_2.fastq.gz"}
        ]
    }
}


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def atomic_json(path, document):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                         prefix=path.name + ".", suffix=".part", delete=False) as out:
            temporary = Path(out.name)
            json.dump(document, out, indent=2, ensure_ascii=False)
            out.write("\n")
            out.flush()
            os.fsync(out.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def file_digests(path):
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    size = 0
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            size += len(chunk)
            md5.update(chunk)
            sha256.update(chunk)
    return size, md5.hexdigest(), sha256.hexdigest()


def fetch_file(entry, destination, timeout, retries):
    path = destination / entry["name"]
    if path.exists():
        size, digest, sha256 = file_digests(path)
        if size != entry["bytes"] or digest != entry["md5"]:
            raise ValueError("Existing file fails archive size/MD5 validation; preserved: " + str(path))
        print("Verified existing " + path.name, file=sys.stderr)
        return {"path": str(path.resolve()), "status": "verified_existing", "sha256": sha256,
                "verified_at_utc": utc_now()}

    for attempt in range(retries + 1):
        temporary = None
        try:
            print("Downloading {} (attempt {}/{})".format(path.name, attempt + 1, retries + 1), file=sys.stderr)
            request = urllib.request.Request(entry["url"], headers={"User-Agent": "trm-control-fetch/1.0"})
            md5 = hashlib.md5()
            sha256 = hashlib.sha256()
            size = 0
            with urllib.request.urlopen(request, timeout=timeout) as response:
                with tempfile.NamedTemporaryFile(mode="wb", dir=destination,
                                                 prefix=path.name + ".", suffix=".part", delete=False) as out:
                    temporary = Path(out.name)
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        size += len(chunk)
                        if size > entry["bytes"]:
                            raise ValueError("Download exceeds expected archive size: " + path.name)
                        out.write(chunk)
                        md5.update(chunk)
                        sha256.update(chunk)
                    out.flush()
                    os.fsync(out.fileno())
            if size != entry["bytes"] or md5.hexdigest() != entry["md5"]:
                raise ValueError("Archive size/MD5 mismatch for {}: {} bytes, MD5 {}".format(
                    path.name, size, md5.hexdigest()))
            os.replace(temporary, path)
            print("Verified " + path.name, file=sys.stderr)
            return {"path": str(path.resolve()), "status": "downloaded_verified", "sha256": sha256.hexdigest(),
                    "verified_at_utc": utc_now()}
        except (OSError, urllib.error.URLError, http.client.HTTPException, ValueError) as error:
            if attempt == retries:
                raise RuntimeError("Failed after {} attempts: {}".format(retries + 1, error)) from error
            print("Retrying after error: " + str(error), file=sys.stderr)
            time.sleep(min(2 ** attempt, 10))
        finally:
            if temporary is not None and temporary.exists():
                temporary.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data/controls"), help="Download directory (default: data/controls)")
    parser.add_argument("--control", action="append", choices=sorted(CONTROLS), help="Select a control; repeat for multiple (default: arabidopsis and drosophila WGS)")
    parser.add_argument("--timeout", type=float, default=60, help="Socket timeout in seconds (default: 60)")
    parser.add_argument("--retries", type=int, default=2, help="Retries after the first failed attempt (default: 2)")
    parser.add_argument("--list", action="store_true", help="Print pinned metadata JSON without downloading or creating directories")
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0 or args.retries < 0:
        parser.error("--timeout must be positive and --retries must be nonnegative")
    names = list(dict.fromkeys(args.control or DEFAULT_CONTROLS))
    controls = {name: copy.deepcopy(CONTROLS[name]) for name in names}
    for control in controls.values():
        if "run_accession" in control:
            control["ena_metadata_url"] = (
                "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + control["run_accession"]
                + "&result=read_run&fields=run_accession,scientific_name,library_strategy,library_source,"
                "library_selection,read_count,base_count,fastq_ftp,fastq_bytes,fastq_md5&format=tsv")
    manifest = {"schema_version": 1, "metadata_verified_date": "2026-09-20",
                "downloaded_at_utc": None, "complete": False, "controls": controls,
                "total_expected_compressed_bytes": sum(f["bytes"] for c in controls.values() for f in c["files"]),
                "note": "Pinned whole WGS runs or an explicitly selected reference chromosome, not motif-selected read subsets. MD5 is used for archival content integrity, not authentication."}
    if args.list:
        print(json.dumps(manifest, indent=2, ensure_ascii=False))
        return 0
    args.output.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output / "controls_manifest.json"
    atomic_json(manifest_path, manifest)
    try:
        for control in controls.values():
            for entry in control["files"]:
                entry.update(fetch_file(entry, args.output, args.timeout, args.retries))
                atomic_json(manifest_path, manifest)
        manifest["complete"] = True
        manifest["downloaded_at_utc"] = utc_now()
        atomic_json(manifest_path, manifest)
    except (OSError, ValueError, RuntimeError) as error:
        manifest["error"] = str(error)
        atomic_json(manifest_path, manifest)
        print("Error: " + str(error), file=sys.stderr)
        return 1
    print(str(manifest_path.resolve()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
