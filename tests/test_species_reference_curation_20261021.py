"""Species reference evidence must not become genus or current-sample proof."""

import csv
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from prepare_taxon_priors import assign, load_references
from trm_assessment import catalog_index, reference_agreement


def sample(taxid, name, lineage_ids):
    return {"sample_id": "unrelated_archive", "status": "resolved", "tax_id": taxid,
            "scientific_name": name, "lineage": "Eukaryota",
            "lineage_tax_ids": lineage_ids}


class SpeciesReferenceCurationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.refs = load_references(ROOT / "data/taxon_motif_references.tsv")
        with (ROOT / "data/telomere_motifs.tsv").open() as handle:
            cls.catalog = list(csv.DictReader(handle, delimiter="\t"))

    def test_supported_species_replace_only_their_broader_expectation(self):
        cases = [
            ("45157", "Cyanidioschyzon merolae", "2759;2763", "AATGGGGGG"),
            ("6661", "Artemia franciscana", "2759;33208;6656;6660", "TTAGG"),
            ("36056", "Boletus edulis", "2759;4751;5204;36055", "TTAGG"),
            ("7425", "Nasonia vitripennis", "2759;33208;6656;50557;7424", "TTATTGGG"),
        ]
        for taxid, name, lineage, motif in cases:
            with self.subTest(name=name):
                rows, priors = assign([sample(taxid, name, lineage)], self.refs)
                self.assertEqual(rows[0]["mapping_status"], "exact_species_reference")
                self.assertEqual(rows[0]["reference_status"], "supported")
                self.assertEqual(rows[0]["related_motifs"], motif)
                self.assertEqual(priors[0]["related_motifs"], motif)
                # Taxonomic expectations cannot assert the current sample's end.
                self.assertIn("Location remains unconfirmed in this sample", rows[0]["notes"])

    def test_10d_descendant_matches_without_collapsing_strain_identity(self):
        rows, priors = assign([sample("280699", "Cyanidioschyzon merolae strain 10D",
                                     "2759;2763;45157")], self.refs)
        self.assertEqual(rows[0]["tax_id"], "280699")
        self.assertEqual(rows[0]["mapping_basis"], "lineage_taxid")
        self.assertEqual(priors[0]["related_motifs"], "AATGGGGGG")

    def test_species_rows_do_not_propagate_to_uncurated_relatives(self):
        cases = [
            ("999991", "Cyanidioschyzon sp.", "2759;2763", ""),
            ("999992", "Artemia sp.", "2759;33208;6656;6660", "TTAGGG"),
            ("999993", "Boletus sp.", "2759;4751;5204;36055", "TTAGGG"),
            ("999994", "Nasonia sp.", "2759;33208;6656;50557;7424", "TTAGG"),
        ]
        for taxid, name, lineage, expected in cases:
            with self.subTest(name=name):
                rows, _ = assign([sample(taxid, name, lineage)], self.refs)
                self.assertEqual(rows[0]["related_motifs"], expected)
                self.assertNotEqual(rows[0]["reference_scope"], "species")

    def test_nasonia_assembly_evidence_keeps_older_candidate_and_figure_limit(self):
        rows, _ = assign([sample("7425", "Nasonia vitripennis",
                                 "2759;33208;50557;7424")], self.refs)
        self.assertIn("nasonia_vitripennis_dallabenetta2020", rows[0]["reference_ids"])
        self.assertIn("nasonia_vitripennis_zhou2022", rows[0]["review_reference_ids"])
        self.assertEqual(rows[0]["reference_evidence"], "assembly_terminal")
        self.assertIn("not visually inspected", rows[0]["notes"])

    def test_catalog_recognizes_red_algal_unit_without_a_rhodophyta_default(self):
        entries = [r for r in self.catalog if r["motif"] == "AATGGGGGG"]
        self.assertEqual(len(entries), 1)
        self.assertEqual(entries[0]["evidence"], "experimental_terminal")
        self.assertTrue(catalog_index(entries))
        self.assertFalse(any(r["scope"] == "clade" and r["scientific_name"] == "Rhodophyta"
                             for r in self.refs))

    def test_boletus_assembly_evidence_is_not_experimental_terminal(self):
        ref = next(r for r in self.refs if r["reference_id"] == "boletus_edulis_brejon2025")
        self.assertEqual(ref["scope"], "species")
        self.assertEqual(ref["evidence"], "assembly_terminal")
        self.assertEqual(ref["sample_ids"], [])
        self.assertIn("not physical-terminal FISH", ref["notes"])

    def test_references_alone_do_not_confirm_current_sample_main_trm(self):
        rows, _ = assign([sample("45157", "Cyanidioschyzon merolae", "2759;2763")], self.refs)
        result = reference_agreement([], None, rows[0], [])
        self.assertEqual(result["status"], "expected_reference_unmeasured")
        self.assertIsNone(result["main_trm_motif"])
        self.assertEqual(result["main_trm_status"], "not_confirmed")

    def test_hericium_does_not_gain_species_confirmation_from_other_specimen(self):
        rows, _ = assign([sample("91752", "Hericium erinaceus", "2759;4751;5204")], self.refs)
        self.assertEqual(rows[0]["reference_scope"], "clade")
        self.assertEqual(rows[0]["related_motifs"], "TTAGGG")


if __name__ == "__main__":
    unittest.main()
