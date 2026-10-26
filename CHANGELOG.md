# Change log

## 0.8.5-rc2 — source distribution candidate, not yet published

- Prepare a standalone Linux research source distribution with English workflow, method, resource and validation documentation.
- Separate the native-only Cargo package from the full Python workflow source bundle.
- Preserve small public test fixtures and explicitly distinguish research-archive-only checks.
- Provide a generated synthetic paired workflow for installation and pipeline smoke checks.
- Synchronize fresh native builds with the supplied 23-entry motif catalog. An older locally managed binary embedded 22 entries: `TTATTGGG` now carries `assembly_terminal` rather than `comparative` evidence, and `AATGGGGGG` is included with `experimental_terminal` evidence. These are catalog updates already present in the supplied reference table, not new biological discoveries in this release.
- Correct a stale test expectation to reflect the current catalog evidence. Sequence detection counts and scanner thresholds are unchanged; reference eligibility and primary-versus-rejected interpretation can change when rebuilding against the current catalog. Record both the catalog and executable hashes when comparing results.

The packaged analysis components remain native **0.8.5**, interpretation policy **0.7.2**, reporting **1.13.0**, assessment schema **22**. Packaging and documentation do not constitute a new telomere detector or new biological validation. No public tag, release date, AWS capability or binary-platform support is implied by this preparation entry.
