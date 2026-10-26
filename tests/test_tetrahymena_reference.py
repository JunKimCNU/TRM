"""Species expectation must not become identity verification or genus inheritance."""
import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location("tetrahymena_taxon_priors", ROOT / "scripts/prepare_taxon_priors.py")
priors_module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(priors_module)


class TetrahymenaReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.references = priors_module.load_references(ROOT / "data/taxon_motif_references.tsv")

    def sample(self, taxid="5911", name="Tetrahymena thermophila", lineage_ids="2759;5782;5908"):
        return {"sample_id": "metadata_identity_only", "status": "resolved", "tax_id": taxid,
                "scientific_name": name, "lineage": "Eukaryota; Ciliophora; Tetrahymena;",
                "lineage_tax_ids": lineage_ids}

    def test_species_gets_published_expectation_without_changing_identity(self):
        rows, priors = priors_module.assign([self.sample()], self.references)
        self.assertEqual(rows[0]["mapping_status"], "exact_species_reference")
        self.assertEqual(rows[0]["tax_id"], "5911")
        self.assertEqual(priors[0]["related_motifs"], "TTGGGG")
        self.assertIn("tetrahymena_thermophila_sheng1995", rows[0]["reference_ids"])

    def test_verified_species_descendant_can_use_expectation(self):
        rows, priors = priors_module.assign(
            [self.sample("312017", "Tetrahymena thermophila SB210", "2759;5908;5911")],
            self.references)
        self.assertEqual(priors[0]["related_motifs"], "TTGGGG")
        self.assertEqual(rows[0]["mapping_basis"], "lineage_taxid")

    def test_other_tetrahymena_species_does_not_inherit_species_row(self):
        rows, priors = priors_module.assign(
            [self.sample("5910", "Tetrahymena pyriformis")], self.references)
        self.assertEqual(priors, [])
        self.assertNotIn("tetrahymena_thermophila_sheng1995", rows[0]["reference_ids"])


if __name__ == "__main__":
    unittest.main()
