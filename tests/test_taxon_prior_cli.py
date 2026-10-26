"""Archive identity can recover a curated distant motif without reading FASTQ again."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import unittest

import test_cli
from test_cli import canonical, fastq


class TaxonPriorIntegrationTests(unittest.TestCase):
    setUpClass = classmethod(test_cli.CliTests.setUpClass.__func__)
    setUp = test_cli.CliTests.setUp
    write = test_cli.CliTests.write
    run_scan = test_cli.CliTests.run_scan

    def test_metadata_reference_recovers_distant_family_and_keeps_raw_evidence(self):
        motif = "CCCAAGGTTTAAAAGG"
        path = self.write("sample.fastq", fastq([motif * 40] * 3 + ["TTAGGC" * 100] * 20))
        prefix = self.root / "before"
        before, _ = self.run_scan([path], prefix=prefix)
        family = canonical(motif)
        self.assertTrue(any(row["motif"] == family for row in before["rejected_candidates"]))
        self.assertFalse(any(row["motif"] == family for row in before["candidates"]))
        metadata = self.write("metadata.tsv", b"sample_id\tstatus\tscientific_name\ttax_id\tlineage\n"
                              b"SRR1\tresolved\tMeloidogyne hapla\t6305\tNematoda; Meloidogyne;\n")
        project = Path(__file__).resolve().parents[1]
        output = self.root / "priors"
        prepared = subprocess.run([sys.executable, str(project / "scripts/prepare_taxon_priors.py"),
                                   "--metadata", str(metadata), "--references",
                                   str(project / "data/taxon_motif_references.tsv"),
                                   "--output-dir", str(output)], capture_output=True, text=True, timeout=10)
        self.assertEqual(prepared.returncode, 0, prepared.stderr)
        with (output / "priors.tsv").open() as handle:
            prior = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(prior["related_motifs"], motif)
        path.unlink()
        after_path = self.root / "after"
        command = [str(self.binary), "reclassify", str(prefix) + ".json", "--related-motif", motif,
                   "--lineage", prior["lineage"], "--output", str(after_path)]
        result = subprocess.run(command, capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)
        after = json.loads(Path(str(after_path) + ".json").read_text())
        recovered = next(row for row in after["candidates"] if row["motif"] == family)
        self.assertEqual(recovered["telomere_prior"]["class"], "related_exact")
        self.assertEqual(recovered["reads"], 3)
        # Reference context must not overrule stronger observed repeat evidence.
        self.assertEqual(after["candidates"][0]["motif"], canonical("TTAGGC"))
        self.assertFalse(after["reclassification"]["fastq_reread"])
        def evidence(report):
            return {row["motif"]: {k: v for k, v in row.items()
                                   if k not in {"telomere_prior", "support", "candidate_eligible"}}
                    for row in report["candidates"] + report["rejected_candidates"]}
        self.assertEqual(evidence(before), evidence(after))


if __name__ == "__main__":
    unittest.main()
