# Reference data and related work

TRM ships two small curated tables:

- [telomere_motifs.tsv](../data/telomere_motifs.tsv): motif families and supporting sources; embedded in the native build.
- [taxon_motif_references.tsv](../data/taxon_motif_references.tsv): species/isolate/genus/clade scope, evidence status, architecture and exceptions used by Python assignment/reporting.

These are selected references, not a complete catalog of eukaryotic telomeres. Preserve each row's source and scope when extending them. A published motif spelling describes an orientation; rotation/RC family equivalence does not remove the need to record biological orientation separately. Rebuild the native binary when changing its embedded catalog, and record both table hashes in an analysis.

Examples illustrate why simple “close to TTAGGG” or base-composition rules are insufficient:

- *Caenorhabditis elegans*: TTAGGC chromosome-end and capping evidence, [Wicky et al. (1996)](https://pubmed.ncbi.nlm.nih.gov/8799140/).
- *Allium*: the 12 bp repeat CTCGGTTATGGG and telomerase synthesis, [Fajkus et al. (2016)](https://pubmed.ncbi.nlm.nih.gov/26716914/). Onion should not be labeled motif-free or ALT solely because the common plant repeat was not found.

Source URLs for additional taxa and variants are retained directly in the tables. A reference's supporting experiment applies to its studied material, not automatically to a new WGS sample. Curated clade defaults remain weak expectations.

Related independently distributed software includes [TREW](https://github.com/Chemical118/TREW), [tidk](https://github.com/tolkit/telomeric-identifier), [TelFinder](https://github.com/bio-tombow/TelFinder), and [TeloSearchLR](https://github.com/gchchung/TeloSearchLR). They are not bundled or required for the basic TRM workflow. Cite the version and primary publication of any tool actually used; different inputs, thresholds and outputs prevent treating agreement as a controlled speed or sensitivity comparison.

SRA provenance also matters: [NCBI describes normalized, original and Lite representations](https://www.ncbi.nlm.nih.gov/sra/docs/sra-data-formats/). Lite-derived per-read quality values are not original per-base quality measurements. Record the representation and conversion commands when preparing FASTQ.
