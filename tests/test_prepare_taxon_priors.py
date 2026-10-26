"""Taxonomic identity must come from metadata, never predicted repeat content."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "prepare_taxon_priors.py"
sys.path.insert(0, str(SCRIPT.parent))
SPEC = importlib.util.spec_from_file_location("prepare_taxon_priors", SCRIPT)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def sample(taxid="6239", name="Caenorhabditis elegans", **changes):
    return {"sample_id": "SRR1", "status": "resolved", "tax_id": taxid,
            "scientific_name": name, "lineage": "Eukaryota; Metazoa; Nematoda; Caenorhabditis;",
            **changes}


def reference(taxid="6239", name="Caenorhabditis elegans", motif="TTAGGC", **changes):
    return {"reference_id": "ce_reference", "reference_taxid": taxid,
            "scientific_name": name, "scope": "species", "related_motifs": motif,
            "motifs": [motif], "evidence": "experimental", "source": "https://example.org/paper",
            "notes": "Test fixture, not biological evidence.", **changes}


class TaxonPriorTests(unittest.TestCase):
    def test_exact_taxid_uses_all_applicable_references(self):
        refs = [reference(), reference(motif="TTAGGT", reference_id="second_ref")]
        rows, priors = module.assign([sample()], refs)
        self.assertEqual(priors[0]["related_motifs"], "TTAGGC,TTAGGT")
        self.assertEqual(rows[0]["mapping_status"], "exact_species_reference")
        self.assertIn("TaxID 6239", priors[0]["lineage"])

    def test_related_genus_or_phylum_does_not_inherit_species_reference(self):
        rows, priors = module.assign([sample("2305860", "Caenorhabditis uteleia")], [reference()])
        self.assertEqual(priors, [])
        self.assertEqual(rows[0]["mapping_status"], "global_prior_no_curated_reference")

    def test_taxid_not_identical_name_is_identity_anchor(self):
        self.assertFalse(module.assign([sample("999999")], [reference()])[1])
        self.assertTrue(module.assign([sample(name="historical taxon synonym")], [reference()])[1])

    def test_computational_reference_remains_review_only(self):
        rows, priors = module.assign([sample()], [reference(evidence="computational")])
        self.assertEqual(priors, [])
        self.assertEqual(rows[0]["review_motifs"], "TTAGGC")

    def test_curated_genus_matches_exact_taxonomic_token_only(self):
        ref = reference("6237", "Caenorhabditis", scope="genus")
        rows, priors = module.assign([sample(lineage_tax_ids="1; 6237;")], [ref])
        self.assertEqual(rows[0]["mapping_status"], "curated_genus_reference")
        self.assertTrue(priors)
        fallback = module.assign([sample()], [ref])
        self.assertEqual(fallback[0][0]["mapping_basis"], "lineage_name_fallback")
        self.assertTrue(fallback[1])
        self.assertFalse(module.assign([sample(lineage="Nematoda; NotCaenorhabditis;")], [ref])[1])
        self.assertFalse(module.assign([sample(lineage_tax_ids="9999")], [ref])[1])

    def test_nearest_clade_default_and_exact_exception_precedence(self):
        refs = [reference("33208", "Metazoa", "TTAGGG", scope="clade", reference_id="animals",
                          evidence="clade_consensus"),
                reference("6231", "Nematoda", scope="clade", reference_id="nematodes",
                          evidence="clade_consensus")]
        unresolved_species = sample("999", "Uncurated nematode", lineage_tax_ids="2759;33208;6231")
        rows, priors = module.assign([unresolved_species], list(reversed(refs)))
        self.assertEqual(priors[0]["related_motifs"], "TTAGGC")
        self.assertEqual(rows[0]["mapping_status"], "clade_default_reference")
        self.assertEqual(rows[0]["mapping_basis"], "lineage_taxid")
        self.assertIn("animals", rows[0]["suppressed_reference_ids"])
        self.assertIn("weak expectation", rows[0]["notes"])
        refs.append(reference("999", "Uncurated nematode", "TTAGGT", reference_id="exception"))
        self.assertEqual(module.assign([unresolved_species], refs)[1][0]["related_motifs"], "TTAGGT")

    def test_species_reference_can_match_verified_descendant_taxid(self):
        rows, priors = module.assign([sample("100000", "C. elegans strain", lineage_tax_ids="2759;6239")],
                                    [reference()])
        self.assertEqual(priors[0]["related_motifs"], "TTAGGC")
        self.assertEqual(rows[0]["mapping_basis"], "lineage_taxid")

    def test_isolate_accession_alias_does_not_merge_taxids(self):
        isolate = reference("2950576", "Panagrellus sp. LJ2284", "TTAGAC", scope="isolate",
                            sample_ids=["DRR674033"], reference_id="lj2284")
        meta = sample("6233", "Panagrellus redivivus", sample_id="DRR674033")
        rows, priors = module.assign([meta], [reference("6233", "Panagrellus redivivus"), isolate])
        self.assertEqual(priors[0]["related_motifs"], "TTAGAC")
        self.assertEqual(rows[0]["tax_id"], "6233")
        self.assertEqual(rows[0]["mapping_basis"], "curated_accession_isolate_alias")
        meta["sample_id"] = "unlinked_run"
        self.assertFalse(module.assign([meta], [isolate])[1])

    def test_specific_unknown_and_complex_block_ancestor_default(self):
        broad = reference("6231", "Nematoda", scope="clade", reference_id="broad")
        for status in ["candidate", "unknown", "out_of_scope"]:
            specific = reference(evidence="comparative", reference_status=status,
                                 motifs=[] if status != "candidate" else ["TTAGGT"])
            rows, priors = module.assign([sample()], [broad, specific])
            self.assertFalse(priors)
            self.assertEqual(rows[0]["mapping_status"], "specific_reference_" + status)
            self.assertEqual(rows[0]["suppressed_reference_ids"], "broad")

    def test_bundled_clade_and_exception_catalogue(self):
        refs = module.load_references(SCRIPT.parents[1] / "data/taxon_motif_references.tsv")
        cases = [
            ("4679", "Allium cepa", "Eukaryota; Embryophyta; Asparagales; Allium;", "CTCGGTTATGGG"),
            ("103475", "Cestrum elegans", "Eukaryota; Embryophyta; Cestrum;", "TTTTTTAGGG"),
            ("193134", "Genlisea hispidula", "Eukaryota; Embryophyta; Genlisea;", "TTCAGG,TTTCAGG"),
            ("5476", "Candida albicans", "Eukaryota; Fungi; Saccharomycotina; Candida;", "ACGGATGTCTAACTTCTTGGTGT"),
            ("28985", "Kluyveromyces lactis", "Eukaryota; Fungi; Saccharomycotina; Kluyveromyces;", "ACGGATTTGATTAGGTATGTGGTGT"),
            ("5062", "Aspergillus oryzae", "Eukaryota; Fungi; Pezizomycotina; Aspergillus;", "TTAGGGTCAACA"),
            ("2305860", "Caenorhabditis uteleia", "Eukaryota; Metazoa; Nematoda; Caenorhabditis;", "TTAGGT"),
            ("6233", "Panagrellus redivivus", "Eukaryota; Metazoa; Nematoda; Panagrellus;", "TTAGAC"),
        ]
        for taxid, name, lineage, expected in cases:
            with self.subTest(name=name):
                rows, priors = module.assign([sample(taxid, name, lineage=lineage)], refs)
                self.assertEqual(priors[0]["related_motifs"], expected)
                self.assertEqual(rows[0]["reference_status"], "supported")
        for taxid, name, lineage in [
                ("4932", "Saccharomyces cerevisiae", "Eukaryota; Fungi; Saccharomycotina; Saccharomyces;"),
                ("4896", "Schizosaccharomyces pombe", "Eukaryota; Fungi; Schizosaccharomyces;"),
                ("6306", "Meloidogyne incognita", "Eukaryota; Metazoa; Nematoda; Meloidogyne;"),
                ("7227", "Drosophila melanogaster", "Eukaryota; Metazoa; Insecta; Diptera; Drosophila;")]:
            with self.subTest(name=name):
                rows, priors = module.assign([sample(taxid, name, lineage=lineage)], refs)
                self.assertFalse(priors)
                self.assertEqual(rows[0]["reference_status"], "out_of_scope")
        rows, priors = module.assign([sample("999", "Uncurated fungus", lineage="Eukaryota; Fungi;")], refs)
        self.assertFalse(priors)
        self.assertEqual(rows[0]["reference_status"], "unknown")

    def test_empty_motif_allowed_only_for_explicit_non_sequence_context(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "refs.tsv"
            ref = reference(related_motifs="", reference_status="out_of_scope")
            ref.pop("motifs")
            module.write_table(path, [ref], list(ref))
            self.assertEqual(module.load_references(path)[0]["motifs"], [])
            path.unlink()
            ref["reference_status"] = "supported"
            module.write_table(path, [ref], list(ref))
            with self.assertRaises(ValueError):
                module.load_references(path)

    def test_specific_exception_blocks_broader_rule_and_exact_wins(self):
        genus = reference("6237", "Caenorhabditis", scope="genus", reference_id="genus")
        exact = reference(motif="TTAGGT", evidence="computational")
        self.assertFalse(module.assign([sample(lineage_tax_ids="6237")], [genus, exact])[1])
        exact["evidence"] = "experimental_terminal"
        self.assertEqual(module.assign([sample()], [genus, exact])[1][0]["related_motifs"], "TTAGGT")

    def test_unresolved_or_metagenome_never_gets_species_prior(self):
        for changes in ({"status": "lookup_failed"}, {"metagenome": "true"},
                        {"library_source": "METAGENOMIC"}):
            with self.subTest(changes=changes):
                self.assertFalse(module.assign([sample(**changes)], [reference()])[1])

    def test_non_eukaryotic_metadata_flags_review_without_motif_blacklist(self):
        rows, priors = module.assign([sample(lineage="Bacteria; Bacteroidota;")], [reference()])
        self.assertFalse(priors)
        self.assertEqual(rows[0]["mapping_status"], "global_prior_non_eukaryotic_metadata")

    def test_duplicate_samples_and_invalid_resolved_taxid_fail(self):
        for samples in ([sample(), sample()], [sample(taxid="0")]):
            with self.assertRaises(ValueError):
                module.assign(samples, [reference()])

    def test_invalid_reference_scope_and_motif_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "refs.tsv"
            for changes in ({"scope": "phylum"}, {"related_motifs": "TTNGGG"},
                            {"reference_taxid": "0"}, {"source": ""}):
                ref = reference(**changes)
                ref.pop("motifs")
                module.write_table(path, [ref], sorted(module.REF_FIELDS))
                with self.assertRaises(ValueError):
                    module.load_references(path)
                path.unlink()

    def test_cli_provenance_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata, refs, output = root / "meta.tsv", root / "refs.tsv", root / "out"
            module.write_table(metadata, [sample()], sorted(module.META_FIELDS))
            ref = reference()
            ref.pop("motifs")
            module.write_table(refs, [ref], sorted(module.REF_FIELDS))
            command = [sys.executable, str(SCRIPT), "--metadata", str(metadata),
                       "--references", str(refs), "--output-dir", str(output)]
            result = subprocess.run(command, capture_output=True, text=True, timeout=10)
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads((output / "summary.json").read_text())
            self.assertFalse(report["motif_results_used_for_taxonomy"])
            self.assertEqual(report["metadata_sha256"], module.sha256(metadata))
            before = (output / "priors.tsv").read_bytes()
            self.assertNotEqual(subprocess.run(command, capture_output=True, timeout=10).returncode, 0)
            self.assertEqual((output / "priors.tsv").read_bytes(), before)
            with (output / "priors.tsv").open() as handle:
                priors = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(priors[0]["related_motifs"], "TTAGGC")


if __name__ == "__main__":
    unittest.main()
