# Testing and release process

The full source bundle and native Cargo crate have different scopes. The full bundle includes Python workflow modules, public docs, small curated data and test fixtures. The crate contains the native Rust program and its required embedded catalog. Neither distribution should contain raw WGS, downloaded assemblies, research result trees, credentials, installed runtime receipts or the separate Windows downloader workspace.

## Clean-source verification

Test the exact export in a fresh directory, not only the development workspace. The initial compatibility target is Linux x86_64 with Rust/Cargo1.75 and Python3.11. GNU `/usr/bin/time` is required for commands that record native resources. Network access may be needed once to populate the locked Cargo dependency cache.

```bash
cargo check --locked
cargo test --locked
cargo build --release --locked --examples
python3 -B scripts/install_runtime.py --jobs 2
python3 -B scripts/runtime_binary.py
python3 -B scripts/run_tests.py
python3 -B scripts/release_smoke.py --output /tmp/trm-smoke
```

Build the examples before Python testing so checks that use native exporters or irregular-analysis executables can run rather than silently skip for missing binaries. The public test runner skips explicitly identified archive-only tests unless `TRM_RUN_ARCHIVE_TESTS=1` is set with the required research inputs. Record skips separately from passes. Keep small independent oracles/fixtures needed by public tests in the export.

The default `fast-gzip` feature uses zlib-rs. A separate build can exercise the alternative Rust backend without replacing the managed binary:

```bash
CARGO_TARGET_DIR=/tmp/trm-rust-gzip-build cargo test --locked \
  --no-default-features --features rust-gzip
```

Backend tests do not alone prove end-to-end result equivalence. A release check should also run the same small scientific fixtures through freshly built binaries and compare scientific fields while allowing version/runtime provenance to differ as documented.

## Artifact checks

1. Inspect the Git export and `cargo package --list --allow-dirty` independently; their inclusion rules differ. Confirm the native catalog is present and no research data are included.
2. Extract the proposed source archive and repeat the clean-source checks. Verify public local links, mandatory fixtures, failed-input handling, FR pair validation and completed-result resume.
3. Record the final source manifest, component versions, Cargo.lock, reference hashes, build environment, test logs and SHA256 checksums. Never describe a source `cargo check` as testing a newly linked release executable.
4. Review the [dependency inventory](../THIRD_PARTY.md) and supply applicable upstream notices with any binary distribution. Declare architecture/libc requirements for each tested binary; do not infer portability from a successful local build.
5. Finalize repository identity, release tag, attribution/citation details and change log only when known. Publish the reviewed artifacts, not the working research directory. Keep AWS implementation as a later scoped addition.

This document specifies the release procedure; it is not a certificate that a tag has already been published. Consult the actual release checksums and CI logs for the artifacts being used. Small software fixtures establish execution contracts, not biological accuracy across the tree of life.
