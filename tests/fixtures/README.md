# Frozen release test fixtures

These small artifacts are byte-identical copies of existing independent controls,
oracles and diagnostic fixture code. `provenance.json` records their original
project-relative paths, exact sizes and SHA256 digests. The original research
archives remain unchanged. The release tests validate these hashes.

- `mixed/`: 2026-09-23 exact mixed-unit fixtures, independent oracle, original
  preregistration/addendum and generating scripts. The original receipts retain
  their original paths and hashes; the test resolves only those two registered
  path prefixes to the relocated files. `published_87_segments.fasta` is the
  original 19 KB published-segment artifact covered by the frozen receipt; the
  portable mixed CLI tests use the synthetic controls and do not claim
  biological sensitivity from them. Original real-data provenance remains in
  `mixed/validation/preregistered_settings.json`.
- `phase_bridge/` and `tract_context/`: frozen protocols and their receipts.
  Numeric/exhaustive oracles remain implemented independently in the tests.
- `volvariella/`: the original diagnostic runner, copied unchanged for its
  synthetic twelve-record contract tests; it is not a public workflow entrypoint.
- `irregular/`: the supplied published-unit definitions used by helper CLI tests.

Saved real WGS/assembly observations and historical executable comparisons are
archive-only integration tests. They are visible as explicit skips in the
portable suite. To run them in the original populated research workspace, set
`TRM_RUN_ARCHIVE_TESTS=1`; individual required-artifact checks still apply. No
flag downloads missing material. Optional native example helpers must also be
built or selected explicitly for their CLI tests.

Published-segment source: [G3 jkad028, Table S16](https://doi.org/10.1093/g3journal/jkad028).
Irregular unit source: [Yeast yea.3922](https://doi.org/10.1002/yea.3922).
