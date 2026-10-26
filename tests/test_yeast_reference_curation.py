"""Guard the biological scope and evidence boundaries of the yeast curation."""

import importlib.util
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location("yeast_curation_priors", ROOT / "scripts/prepare_taxon_priors.py")
PRIORS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PRIORS)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def rotations(sequence):
    return {sequence[i:] + sequence[:i] for i in range(len(sequence))}


def family(sequence):
    return min(rotations(sequence) | rotations(reverse_complement(sequence)))


def sample(taxid, name):
    return {"sample_id": "sample_" + taxid, "status": "resolved", "tax_id": taxid,
            "scientific_name": name,
            "lineage": "Eukaryota;Fungi;Saccharomycotina",
            "lineage_tax_ids": "2759;4751;147537;" + taxid}


class YeastReferenceCurationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.references = PRIORS.load_references(ROOT / "data/taxon_motif_references.tsv")
        cls.by_taxid = {row["reference_taxid"]: row for row in cls.references}

    def test_supported_species_emit_only_their_exact_reference(self):
        for taxid in ("4952", "4911"):
            with self.subTest(taxid=taxid):
                ref = self.by_taxid[taxid]
                assignments, priors = PRIORS.assign([sample(taxid, ref["scientific_name"])], self.references)
                self.assertEqual(ref["scope"], "species")
                self.assertEqual(ref["sample_ids"], [])
                self.assertEqual(assignments[0]["mapping_status"], "exact_species_reference")
                self.assertEqual(assignments[0]["reference_ids"], ref["reference_id"])
                self.assertEqual(assignments[0]["related_motifs"], ref["related_motifs"])
                self.assertEqual(len(priors), 1)

    def test_template_candidates_remain_review_only(self):
        for taxid in ("33165", "381046"):
            with self.subTest(taxid=taxid):
                ref = self.by_taxid[taxid]
                assignments, priors = PRIORS.assign([sample(taxid, ref["scientific_name"])], self.references)
                self.assertEqual(ref["scope"], "species")
                self.assertEqual(assignments[0]["reference_status"], "candidate")
                self.assertEqual(assignments[0]["review_motifs"], ref["related_motifs"])
                self.assertEqual(assignments[0]["related_motifs"], "")
                self.assertEqual(priors, [])

    def test_uncurated_relative_does_not_inherit_species_motif(self):
        relative = sample("51656", "Kluyveromyces dobzhanskii")
        assignments, priors = PRIORS.assign([relative], self.references)
        self.assertEqual(assignments[0]["reference_status"], "unknown")
        self.assertEqual(assignments[0]["related_motifs"], "")
        self.assertEqual(assignments[0]["review_motifs"], "")
        self.assertEqual(priors, [])

    def test_source_polarity_and_distinct_kluyveromyces_families(self):
        # 2021 S1 writes the C-rich strand for these species. Complementing
        # without reversing would quietly reverse the end-polarity convention.
        self.assertEqual(self.by_taxid["33165"]["related_motifs"],
                         reverse_complement("ATCCGTACACCACATACCTAATCAA"))
        self.assertEqual(self.by_taxid["381046"]["related_motifs"], reverse_complement("CCAACACC"))
        # 1994 Fig3 explicitly gives the strand pointing 5' -> chromosome end.
        self.assertIn(self.by_taxid["4911"]["related_motifs"],
                      rotations("ACGGATTTGATTAGTTATGTGGTGT"))
        self.assertNotEqual(family(self.by_taxid["4911"]["related_motifs"]),
                            family(self.by_taxid["28985"]["related_motifs"]))


if __name__ == "__main__":
    unittest.main()
