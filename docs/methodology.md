# Method and interpretation

TRM is a candidate-screening workflow. It separates the sequence that was observed, a software classification, a literature expectation, and direct evidence of chromosome-end localization.

```text
FASTA/FASTQ records from one sample
  -> parse and validate complete records; decode gzip when present
  -> bounded batches -> repeat-discovery workers
  -> exact tandem seeds across the full sequence
  -> substitution-tolerant extension; retain intervals and orientation
  -> primitive motif / rotation / reverse-complement family aggregation
  -> counts, QC, boundaries, optional paired and quality observations
  -> native JSON plus primary and rejected candidate tables
  -> taxonomy/reference assignments and separate Python assessment
  -> evidence axes, software interpretation, limitations and conditional power
```

## Discovery

The scanner examines each complete record rather than constructing a genome-wide k-mer frequency histogram. Here `k` is the repeat period: `TTAGGG` has period 6 even though it contains three distinct bases. No three-distinct-base filter is imposed.

For a tested period `k`, the minimum exact seed length is:

```text
k × max(min_copies, ceil(seed_bases / k))
```

Defaults are `seed_bases=24`, `min_copies=3`, and periods 3–32. Thus a 6 bp unit needs 24 exact tandem bases, a 5 bp unit needs 25, and a 32 bp unit needs 96. Seeds may start at any phase. The native maximum configurable period is 64; this does not automatically extend every downstream adequacy or targeted-selection model to that range.

The detector extends an exact seed in both directions against a periodic template, permitting substitutions under the configured score and final identity threshold, default 0.85. Ambiguous bases interrupt extension. It does not perform general insertion/deletion alignment. An error-rich repeat with no sufficient exact seed may be missed even when its average identity is high.

Repeats are reduced to primitive periods and equivalent rotations/reverse complements share a family. For example, `TTAGGG`, `TAGGGT` and `CCCTAA` have canonical spelling `AACCCT`. This lexical spelling is not the biological G-rich orientation. Coordinate intervals are zero-based and half-open, `[start,end)`; orientation remains a separate observation.

Candidate summaries include supporting records, repeat-base accounting and the longest observed tract. Separated tracts and different reads are not joined into an artificial long tract. The default `--min-tract 500` annotates single-record long-tract evidence; it is not a universal 500 bp exact requirement for every telomere or a rule discarding all short-read positives. Discovery defaults to reporting families with at least two supporting records.

## Reference policy and candidate visibility

The native default `canonical` policy promotes exact established or explicitly supplied related families. Other reportable repeats remain in `rejected_candidates` and the rejected table. They are not biological negatives. `--prior-policy discovery` retains all reportable families as primary while retaining reference annotation/tier ordering; it does not remove the reference system. The legacy `telomere` policy permits configured edit neighbours and should not be treated as a calibrated confidence model.

The global motif catalog and taxon-specific expectations answer different questions. An exact globally known family match does not establish that the analyzed organism uses it at chromosome ends. Species/isolate exceptions take precedence over broader defaults; unknown or complex references can block inappropriate inheritance. Taxonomic lineage matching is not a phylogenetic model with branch lengths.

Python reporting reads both primary and rejected observations. Native ordering, interpretation scores, the selected representative, known-family observation, and main/novel TRM confirmation remain separate. A representative outside the expected family does not imply the expected family was absent. Neither a particular string blacklist nor blanket threshold relaxation is used to force an expected answer.

## Boundaries, pairs and quality

Repeat-rich reads and repeat–exterior boundaries provide complementary evidence. Relative to a declared G-rich reference motif, a repeat toward a left endpoint is expected in the reverse-complement orientation and toward a right endpoint in the forward orientation. Counts at read boundaries are not direct measurements of physical chromosome termini. A scaffold boundary may be an assembly break.

With `--paired-layout fr`, exactly two synchronized FASTQ files are checked for mate identity, order and count. Pair orientation evidence is recorded alongside single-read observations. Mates are not concatenated, overlaps are not merged, PCR copies are not removed, and a pair is not two independent biological replicates. Without this option, input records are pooled without mate validation.

FASTQ structural and quality-length checks are always performed. Optional Phred+33 diagnostics summarize supplied qualities within the configured retention scope; they do not alter discovery. Mean quality does not establish accuracy at a motif-discriminating position, and constant Q30 does not prove calibrated error probabilities. Neither the scanner nor a checksum authenticates how qualities were produced.

Stored flank sketches are bounded representatives, not an exhaustive set of molecules or a uniqueness assay. A representative's longer sequence does not necessarily have support equal to every event associated with its shorter identity key. Selected reference alignments and mixed-unit annotations must retain their own selection, cap and specimen scope.

## Nondetection and independent validation

Keep these states distinct:

| State | Meaning |
|---|---|
| Not reported | Absent from saved discovery candidates under its detection/reporting rules; not an explicit zero measurement |
| Targeted zero | A completed specified target/model search measured zero support in the supplied input |
| Unmeasured | A required observation or supporting input was not collected |
| Technical failure | Invalid input, resource failure or incomplete execution; no biological negative conclusion |
| Conditional power | A nondetection calculation under supplied genome, chromosome, tract, retention, sensitivity and independence assumptions |

The adequacy model cannot infer chromosome-end retention from read count alone. Genome size, karyotype and repeat length describe potential exposure; library construction and independent molecules remain separate assumptions. The illustrative sensitivity grid is not a calibrated result for the sample. Synthetic insertion tests assess the algorithm, not the actual library's telomere recovery.

Missing a candidate cannot by itself distinguish true motif absence, a noncanonical maintenance mechanism, a longer/irregular repeat, incomplete data or library loss. Confirming a new or main TRM requires specimen-linked chromosome-end evidence and appropriate independent validation. Repeated accessions, mates, PCR copies and an assembly built from the same reads are not independent biological replications.

See [validation scope](validation.md) and [primary references](references.md). Implementation entry points are [input parsing](../src/input.rs), [repeat discovery](../src/detect.rs), [motif normalization](../src/motif.rs), [native aggregation](../src/main.rs), and [assessment](../scripts/report_assessment.py).
