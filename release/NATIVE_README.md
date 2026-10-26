# TRM native scanner

TRM 0.8.5 screens complete FASTQ/FASTA records for tandem repeats associated with telomeric repeat motif candidates. It is a research tool: a repeat observation is not proof of physical chromosome-end localization, a novel/main TRM, true absence or ALT.

This Cargo package contains the **native Rust program and its embedded motif catalog**. Python sample batching, taxonomy assignment and assessment reporting belong to the separate full source distribution. They are not installed by Cargo.

## Build and run

Rust 1.75 is the declared minimum and the initial tested build target is Linux x86_64. Other platform and prebuilt-binary compatibility is not asserted.

```bash
cargo build --release --locked
cargo test --locked
target/release/trm --help
target/release/trm scan R1.fastq.gz R2.fastq.gz \
  --paired-layout fr --threads 2 --memory-limit-mib 2048 \
  --output sample
```

All inputs in one invocation must represent the same sample. Gzip is detected from content. Single stdin input is supported; FR pairing requires two distinct synchronized FASTQ paths. Invalid FASTQ/gzip and resource failures are errors, not biological negative results.

Native subcommands are `scan` for discovery, `evidence` for supplied target families, `mixed` for explicitly supplied exact mixed-unit grammars, and `reclassify` for prior reinterpretation of saved JSON. Use `--help` for their actual argument contracts. No SRA container, BAM/CRAM or S3 URL reader is provided.

Default discovery searches primitive periods 3–32 bp with exact seeds of at least 24 bp and at least 3 copies, then permits substitution-tolerant extension at identity 0.85. Periods up to 64 bp are configurable. It does not implement general indel alignment or infer arbitrary mixed grammars. The whole record is used; batch size does not crop reads.

Output is `PREFIX.json`, `PREFIX.tsv`, and `PREFIX.rejected.tsv`. Inspect both primary and rejected observations. Rotation and reverse-complement families share a normalized spelling, while orientation is stored separately. Default canonical annotation uses established exact families. Discovery omission is not an explicitly measured zero, and rejected candidates are not biological negatives.

The 500 bp long-tract setting describes single-record evidence; summed short-read repeat bases cannot prove a contiguous 500 bp tract. FR pairs are checked but not merged, deduplicated or counted as independent biological replicates. Optional Phred+33 diagnostics summarize supplied qualities without altering discovery.

`--memory-limit-mib` is a Linux process address-space limit, not a total RSS budget. Workers share a producer/batch pipeline; extra threads or RAM do not guarantee proportional speedup. AWS orchestration and GPU execution are not implemented.

The default gzip backend uses zlib-rs. The alternative build is `cargo build --release --locked --no-default-features --features rust-gzip`. Dependency versions are pinned in Cargo.lock. TRM is MIT licensed; dependency license terms remain applicable.
