# Change log

## 0.8.5-rc4 — portable invalid-input test

- Give the invalid-manifest CLI regression an explicit one-worker CPU budget. On the four-CPU GitHub runner, its former 32-worker defaults correctly stopped at the resource guard before reaching the intended manifest error. The test now isolates its intended contract on small machines while retaining all input/output assertions.
- Keep resource guards, runtime defaults, detector behavior and strict CI requirements unchanged. The related test module passes with actual one-CPU affinity; full hosted validation is recorded against the published commit.

## 0.8.5-rc3 — hosted CI compatibility follow-up

- Address current stable Clippy warnings in one Rust documentation comment and five test-only single-element slices. Detection logic, thresholds, reference tables and scientific output are unchanged.
- Preserve RC2 as the initial source publication; its hosted Rust tests passed, but stable Rust 1.99.1 rejected the new lint warnings. RC3 keeps the strict CI checks enabled and receives a separate source manifest and archive.

## 0.8.5-rc2 — initial source distribution candidate

- Prepare a standalone Linux research source distribution with English workflow, method, resource and validation documentation.
- Separate the native-only Cargo package from the full Python workflow source bundle.
- Preserve small public test fixtures and explicitly distinguish research-archive-only checks.
- Provide a generated synthetic paired workflow for installation and pipeline smoke checks.
- Synchronize fresh native builds with the supplied 23-entry motif catalog. An older locally managed binary embedded 22 entries: `TTATTGGG` now carries `assembly_terminal` rather than `comparative` evidence, and `AATGGGGGG` is included with `experimental_terminal` evidence. These are catalog updates already present in the supplied reference table, not new biological discoveries in this release.
- Correct a stale test expectation to reflect the current catalog evidence. Sequence detection counts and scanner thresholds are unchanged; reference eligibility and primary-versus-rejected interpretation can change when rebuilding against the current catalog. Record both the catalog and executable hashes when comparing results.

The packaged analysis components remain native **0.8.5**, interpretation policy **0.7.2**, reporting **1.13.0**, assessment schema **22**. Packaging and documentation do not constitute a new telomere detector or new biological validation. A source publication does not imply a public release tag, AWS capability or additional binary-platform support.
