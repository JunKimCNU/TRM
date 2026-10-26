# Resources, throughput and resume

Start with explicit jobs and worker counts. More RAM or more workers does not guarantee faster scans: one producer decodes/parses input and dispatches bounded batches, while workers discover repeats. The default gzip decoder is not a parallel decoder. Independent samples often offer a more useful parallelization axis than adding workers to one file.

`--threads` means scanner workers per sample. `--jobs` means simultaneous samples. `--cpu-budget` checks `jobs × threads`; it does not account for every producer or management thread. The optional `--thread-budget` checks a broader plan and rejects an explicit overflow before native execution.

For the current batch topology:

```text
planned threads = 1 + min(jobs, selected samples) × (workers + 2 + H) + reserve
H = 1 when GNU time wraps each native job, otherwise 0
```

The plan includes native producer/main, Python executor/coordinator and GNU time. It is not an OS quota or a global scheduler. Other runners, decompression/transfer helpers and observers must fit the reserve or a separate total allocation. CPU affinity can differ from container CPU reservation; explicitly configure cloud/container budgets.

`--memory-limit-mib` is a per-native-process Linux virtual-address-space limit, not an RSS reservation or a total job-tree memory limit. The scanner can reduce its batch target to fit its plan, but individual long records, candidate tables and serialization still require space. A failed allocation is a technical failure. Batch size does not split long reads into shorter reads.

For Python reporting, `--memory-limit-mib` limits each interpretation worker and `--memory-budget-mib` checks jobs times that cap. Coordinator memory remains additional. Genome complexity, candidate counts and saved sketches affect reporting cost; native CPU/Gbp is not a substitute for measuring it.

## Scaling across samples

`batch.py` and `evidence_batch.py` bound pending futures and support periodic cohort checkpoints, `--schedule largest-first`, and deterministic `--shard-index/--shard-count`. Largest-first uses input bytes, which need not predict runtime exactly. Manifest/checkpoint state still scales with sample count. Shards partition samples; they do not implement within-sample distributed merging.

Keep paired mates in the same sample and shard. Separate layouts or incompatible library types rather than forcing them through FR validation. Determine practical settings with a fixed representative workload, completed-result equivalence, whole-batch wall time, actual CPU-hour/Gbp, compressed bytes, simultaneous process-tree memory and failed-attempt cost. Local measurements do not establish AWS, GPU or ARM performance.

## Completed-result reuse

`--resume` reuses compatible completed samples after state, input metadata, options and output hashes are checked. Default input identity uses path/size/mtime, not a newly computed whole-FASTQ checksum. A selected verified-download receipt can strengthen discovery identity; it must match the current local files. Moving data to another filesystem requires new local verification.

Batch output uses advisory filesystem locks and per-sample state. These are not distributed S3 locks. An unfinished sample restarts rather than continuing at its previous FASTQ offset. Failed/running state must not be interpreted as success merely because an older result file remains nearby.

Preserve native JSON, run states, parameters, input checksums, metadata, reference-table hashes and any evidence sidecars. Representative motif strings alone are insufficient for reproducible interpretation. Advanced sidecars may refer to validated external paths and need an explicit portability plan.

## Optional profiling

`--performance-profile on` reports source reading, decode/buffering, queue-send waits and a residual including parsing, mate checks, allocation and control work. These are overlapping-producer wall-clock observations, not isolated causal CPU costs or physical disk bytes. Shared-server load, cache and input order affect timing. Profile only when it answers a specified question and compare the same scientific options.

Use each command's [help](cli.md) for current options. AWS staging/orchestration and Spot recovery remain future work; this document describes the local implementation.
