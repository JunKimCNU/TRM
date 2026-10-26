# CLI entry points

The installed code is the authoritative option reference. Run these commands from the full source root; links identify the corresponding parser/implementation.

| Purpose | Help command | Source |
|---|---|---|
| Native scan/evidence/mixed/reclassify | `bin/trm --help` | [main.rs](../src/main.rs) |
| Build/install managed native binary | `python3 -B scripts/install_runtime.py --help` | [installer](../scripts/install_runtime.py) |
| Runtime identity | `python3 -B scripts/runtime_binary.py --help` | [runtime verifier](../scripts/runtime_binary.py) |
| Sample discovery batch | `python3 -B scripts/batch.py --help` | [batch](../scripts/batch.py) |
| Taxon reference assignment | `python3 -B scripts/prepare_taxon_priors.py --help` | [assignment](../scripts/prepare_taxon_priors.py) |
| Assessment/report | `python3 -B scripts/report_assessment.py --help` | [report](../scripts/report_assessment.py) |
| Selected-target batch | `python3 -B scripts/evidence_batch.py --help` | [targeted batch](../scripts/evidence_batch.py) |
| Completed-assessment summary | `python3 -B scripts/summarize_interpretation.py --help` | [summary](../scripts/summarize_interpretation.py) |
| Public Python tests | `python3 -B scripts/run_tests.py --help` | [test runner](../scripts/run_tests.py) |

These entry points do not define a stable Python import API. Advanced evidence options require the schema and provenance expected by their loaders. A directory containing arbitrary alignment or quality files is not interchangeable with a sealed evidence manifest.

The native binary reads FASTA/FASTQ, optionally gzip. It does not directly read SRA containers, BAM/CRAM or S3 URLs. Single stdin input is supported; FR paired mode requires two distinct FASTQ paths and rejects stdin.
