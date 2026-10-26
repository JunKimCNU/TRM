# Validation scope

Validation has three separate purposes: software correctness, biological applicability and performance. Success on one does not prove the others.

## Public software checks

The source bundle provides native and Python tests plus a small generated end-to-end smoke workflow. Relevant contracts include full-record search, seed phases, family normalization, gzip/input errors, FR pairing, resource planning, failure/reuse states, reference assignment and separation of observed versus unmeasured evidence. Independent synthetic oracles are used where supplied in the test fixtures.

```bash
cargo test --locked
python3 -B scripts/run_tests.py
python3 -B scripts/release_smoke.py --output /tmp/trm-smoke
```

Install the managed binary before running Python checks. Run the smoke command in a fresh directory. See [release testing](releasing.md) for clean-source and alternate-backend checks. Use the actual test log and release artifact checksums to identify what passed; this document does not assert that every environment or optional analysis was tested.

Some historical tests require research archives that are not distributed. These are explicitly gated by `TRM_RUN_ARCHIVE_TESTS=1` and skipped in the public default suite. A skipped archive test is not a pass. The public suite excludes Windows downloader tests because that tool is outside this distribution. Synthetic fixture tests must not be reclassified as archive-only merely to hide a failure.

## Biological applicability

Curated references include different evidence types and scopes. A published telomere in one species is not automatic truth for its relatives or for every strain. For example, chromosome capping by the *C. elegans* TTAGGC repeat has direct experimental support, while transferring a clade default to an unstudied sample is only an expectation. [Wicky et al. (1996)](https://pubmed.ncbi.nlm.nih.gov/8799140/).

A valid biological evaluation should fix the reference answer and its terminal evidence, sample/strain linkage, library, searched motif range, depth and independent-replicate definition before calculating a detection rate. Include positive controls and plausible non-telomeric repeats. Distinguish recovery of a known family from correct localization or main-motif selection. A prior-curated development control is not automatically a held-out test.

Short reads can support abundant repeats, full-repeat reads and boundary orientation. Physical ends generally require appropriate long-read/assembly context or independent assays. A scaffold endpoint and a same-source assembly are useful context but not independent confirmation. An observed internal repeat may be real without being telomeric.

No tree-of-life sensitivity/specificity, universal minimum input size, absence/ALT classifier or novel-TRM confirmation rate is asserted. Nondetection remains conditional when terminal retention, genome/karyotype or detector-sensitivity assumptions are missing. See [method and limits](methodology.md).

## Performance evidence

Measure input bases, compressed bytes, wall time, actual user+system CPU, simultaneous memory, thread allocation and failed attempts separately. Report whether input is a prefix or complete supplied files, whether caches may be warm, and whether optional diagnostics are enabled. Sum of concurrent sample wall times is not cohort wall time; sum of per-process peak RSS is not simultaneous memory.

The implementation has been exercised on local WGS cohorts, but those observations are not a public benchmark promise for every taxon, hardware platform or cloud. This source release does not bundle raw cohort data or claim an independently reproduced cross-tool speed advantage. A publishable benchmark should provide accession/object identities, exact commands, hashes, workload/equivalence criteria and repeat measurements.
