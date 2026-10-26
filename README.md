# TRM — telomeric repeat screening from WGS

TRM searches complete FASTQ/FASTA records for tandem repeats and reports evidence for telomeric repeat motif candidates. A Rust scanner performs sequence analysis; Python tools organize samples and interpret saved observations using curated references and explicit library assumptions.

**A repeat call is a screening result, not confirmation of a physical telomere, a species' main motif, a novel motif, or telomere absence.** Read the [method and limits](docs/methodology.md) before interpreting nondetection.

| Component | Version |
|---|---|
| Native scanner | 0.8.5 |
| Interpretation policy | 0.7.2 |
| Assessment reporting / schema | 1.13.0 / 22 |

The source artifact is labeled **0.8.5-rc4**. This is a distribution candidate label; it does not change the native version. Source is available at [JunKimCNU/TRM](https://github.com/JunKimCNU/TRM); a source commit is distinct from a tagged GitHub Release. No standalone prebuilt binary is supplied. Builds embed the supplied catalog; record its hash and the executable hash because older local binaries with the same native version can contain earlier reference tables.

## Install and check

The initial supported workflow targets **Linux x86_64**. The release preparation environment uses **Rust/Cargo 1.75.0 and Python 3.11**. Other Python versions, ARM, macOS, Windows, and prebuilt-binary portability are not established by this release. Core Python commands use the standard library; resource recording requires GNU `/usr/bin/time`.

From the full source bundle's root:

```bash
python3 -B scripts/install_runtime.py --jobs 2
python3 -B scripts/runtime_binary.py
python3 -B scripts/release_smoke.py --output /tmp/trm-smoke
```

The installer builds from `Cargo.lock` into a temporary build directory and installs `bin/trm` with a version/hash receipt. Building may fetch Cargo dependencies; `--offline` requires an existing dependency cache. Use a fresh smoke output directory. The smoke command generates tiny synthetic paired reads and exercises the local pipeline without downloading biological data. It checks software behavior, not biological sensitivity.

The Cargo crate is a smaller **native-only** distribution. The full source bundle also contains the Python workflow, curated references, documentation and public tests. See [testing and release scope](docs/releasing.md).

## Analyze data

For synchronized conventional FR paired FASTQ from one sample:

```bash
bin/trm scan sample_R1.fastq.gz sample_R2.fastq.gz \
  --paired-layout fr --threads 2 --memory-limit-mib 2048 \
  --output sample_scan
```

This writes `sample_scan.json`, `sample_scan.tsv` and `sample_scan.rejected.tsv`. Read both primary and rejected candidate observations. Canonical family spelling groups rotations and reverse complements: the familiar `TTAGGG` family is written `AACCCT` in normalized output.

For multiple samples and taxonomy-aware reporting, follow the [complete workflow](docs/quickstart.md). Keep separate samples in separate manifest rows. Providing two paths without `--paired-layout fr` pools their records without requesting mate validation.

## What TRM provides

- Full-record seed-and-extend discovery, with default motif periods 3–32 bp and configurable native support up to 64 bp.
- Repeat-family counts, interval lengths, full-repeat reads, boundary orientation and optional FR mate evidence.
- Curated reference annotation, retained alternative candidates, and separate observed/unknown/failed states.
- Sample-level parallel execution, bounded task submission, hash-checked completed-result reuse and resource planning.
- Optional targeted and supplied mixed-unit searches; descriptive Phred+33 quality and producer timing diagnostics.

Exact seeds are required and the standard extension model handles substitutions, not general indel alignment. Short-read repeat totals do not establish one continuous long tract. Genome size, chromosome number and sequencing volume alone do not establish library telomere retention or rule out technical nondetection.

AWS orchestration, S3 staging, Spot recovery, GPU acceleration and FR diskless streaming are not provided by this release. Multiple independent samples can be distributed using existing jobs/shards; a single sample cannot be split across machines and its reports simply summed.

[Workflow](docs/quickstart.md) · [Method](docs/methodology.md) · [Resources and resume](docs/resources.md) · [Validation](docs/validation.md) · [CLI entry points](docs/cli.md) · [References](docs/references.md) · [Change log](CHANGELOG.md)

TRM is distributed under the [MIT license](LICENSE). Dependency information is in [THIRD_PARTY.md](THIRD_PARTY.md). No project publication DOI or author list is asserted here; record the release identifier, component versions, parameters and reference-table hashes when describing an analysis.
