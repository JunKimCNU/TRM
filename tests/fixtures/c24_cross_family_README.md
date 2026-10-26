# C24 cross-family traversal regression

`c24_cross_family_context.fasta` is an unedited, soft-masked 1,361 bp excerpt
from Arabidopsis thaliana C24 assembly **GCA_902705455.1 / araThaC24v1**,
record **CACSHJ010000087.1**, source interval **[29,412,411, 29,413,772)**
(zero-based, half-open). The JSON sidecar records archive and fixture checksums,
the exact parameters and the original diagnostic artifacts. The excerpt starts
at a multiple of seven, preserving the original period-seven scan grid.

This context was chosen *after* the complete matched-panel analysis exposed a
de novo/targeted discrepancy. It is a software regression fixture, not an
independent validation sample or a confirmed physical telomere endpoint.
Clipping creates an artificial left edge. The 279 bp tract ends 77 bp before the
original record's right edge, outside the registered 50 bp terminal window.

With frozen v0.8.0, an approximate extension of another seven-base family,
`ACCCTAT`, spans [957, 1361). Advancing the shared period-seven probe to that
extension's end skips a valid exact seed for `AAACCCT` (the rotation and reverse
complement family of `TTTAGGG`). The native scan retains only the earlier 28 bp
tract; targeted evidence independently reports two tracts [821, 849) and
[1005, 1284), totaling 307 bp with maximum 279 bp. Neither tract is 500 bp long.

The CLI regressions require recovery of both tracts, preservation of the other
observed families, all seven phase offsets in both orientations, worker/batch
consistency, and no artificial 500 bp tract from combined support. They use
the fixture's observed coordinates, not a universal assertion that de novo and
targeted algorithms must give identical outputs on every imperfect array.
The reverse-complement transformations uppercase the source before complementing
so that lowercase soft masking does not corrupt the control.

When the original saved 29,413,772 bp full-record fixture is available, an
additional test verifies complete record consumption and recovery of at least
the previously targeted 3,995 repeat bases and 279 bp maximum. The large record
is not copied into the test suite; its exact hash is checked before use.

Run with the current executable:

```sh
python3 tests/test_traversal_regression.py --binary target/release/trm -v
```

`TRM_BINARY` provides the same selection for `unittest discover`. The optional
full record can be selected using `--full-record` or `TRM_C24_FULL_RECORD`.
The original local-context recovery test is expected to fail against the frozen
v0.8.0 executable with `max_tract: 28` instead of `279`.
