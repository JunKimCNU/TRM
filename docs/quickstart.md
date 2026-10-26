# Complete local workflow

Run commands from the full source bundle's root. Install the managed scanner as described in the [README](../README.md). Keep the `scripts/` directory together: the commands share Python modules.

## A self-contained example

```bash
python3 -B scripts/release_smoke.py --output /tmp/trm-smoke
```

Use a new output directory. This creates synthetic paired input and metadata, runs scan/batch/report checks, and records the outcome. It does not download a control genome or demonstrate recovery of physical telomeres. Generated taxonomy labels are fixtures, not inferred organism identities.

Outputs include `native/`, `priors/`, `assessment/index.html`, `comparison/`, `commands.json` and `smoke.json`. Inspect the smoke receipt and command log rather than treating the presence of an HTML file alone as successful completion.

## Your own paired samples

Create a tab-separated `samples.tsv` with exactly these columns. Paths may be absolute or relative to the manifest directory. The example labels and paths below must be replaced with your data.

```tsv
sample_id	path1	path2
sample01	reads/sample01_R1.fastq.gz	reads/sample01_R2.fastq.gz
```

For single-input samples, use the two-column header `sample_id<TAB>path1`. Use separate batches for FR paired and single/unknown-layout libraries. IDs must be unique accession-safe labels using letters, digits, `_`, `.` and `-`.

```bash
python3 -B scripts/batch.py samples.tsv \
  --jobs 1 --threads 1 --cpu-budget 1 \
  --thread-budget 8 --thread-reserve 1 --memory-limit-mib 2048 \
  --paired-layout fr --record-resources \
  --output-dir work/native
```

These are modest starting settings, not universal minimum memory or optimal throughput. A job includes a producer and management threads in addition to scanner workers. See [resource planning](resources.md). Add `--resource-plan-only` to inspect the planned thread count without scanning; this is not a FASTQ integrity check.

The batch stores each sample under `work/native/SAMPLE/`, with `result.json` and execution state. To resume, repeat the same command with `--resume`. Only compatible completed results are reused. Changed inputs or options can require a new scan; a partial scan does not resume from a FASTQ byte offset.

## Metadata and reference assignments

Create a metadata TSV containing at least `sample_id`, `status`, `scientific_name`, `tax_id`, and `lineage`. The following is an illustration for a genuinely identified *C. elegans* sample, not a value to copy into an unidentified sample:

```tsv
sample_id	status	scientific_name	tax_id	lineage
sample01	resolved	Caenorhabditis elegans	6239	Eukaryota;Metazoa;Nematoda;Caenorhabditis
```

Use `status=unresolved` and leave unknown organism fields empty when identity is unavailable. Additional source metadata can be retained. Taxonomy must come from specimen/archive evidence; do not infer it from the motif. Lineage TaxIDs, when supplied, take precedence over name fallback. Clade defaults are weaker evidence than a matched species or isolate reference.

```bash
python3 -B scripts/prepare_taxon_priors.py \
  --metadata metadata.tsv \
  --references data/taxon_motif_references.tsv \
  --output-dir work/priors

python3 -B scripts/report_assessment.py \
  --source work/native --metadata metadata.tsv \
  --assignments work/priors/assignments.tsv \
  --catalog data/telomere_motifs.tsv \
  --taxon-references data/taxon_motif_references.tsv \
  --jobs 1 --json-format compact \
  --output work/assessment --comparison-output work/comparison
```

Use distinct input/output directories. The assignment and reporting steps do not reread FASTQ or alter the saved discovery observations. They connect those observations to explicit taxonomy and references. The comparison directory records interpretation differences; it is not automatically a comparison with a previous software release. A genuine previous assessment can be supplied with `--previous-assessment`.

Open `work/assessment/index.html` and inspect `samples.tsv`, `evidence_axes.tsv`, and the sample JSON files. Check execution completion and input scope before biological interpretation. Without sample-appropriate genome/library assumptions, nondetection power will often remain `unassessable_library_or_genome`; this does not invalidate a positive repeat observation.

## Optional diagnostics

If the FASTQ quality encoding is known to be Phred+33, `scan` and batch support `--quality-diagnostics phred33`. It describes quality distributions without trimming, filtering or reweighting reads. Do not use it to authenticate constant or transformed qualities. `--performance-profile on` records producer wall-clock components, not isolated CPU or physical disk costs. Both options default to off.

`--flank-period-annotation on` adds candidate-period flank descriptions without changing counts or calls. Its absence does not mean the flank is proven nonrepetitive. Advanced saved-evidence adapters require their own validated provenance contracts; they are not needed for the basic workflow.

For a separate full-input check of an explicitly selected motif:

```bash
printf '{"motifs":["TTAGGG"]}\n' > targets.json
bin/trm evidence sample_R1.fastq.gz sample_R2.fastq.gz \
  --targets targets.json --paired-layout fr --threads 2 \
  --output sample_targeted
```

Targeted analysis rereads the input and preserves measured zero-hit targets. Discovery omission and targeted zero are different observations. Do not overwrite discovery counts with targeted values. The complete interfaces are listed in [CLI entry points](cli.md).
