"""Metadata text must not turn missing protocols or ordinary PCR into WGA facts."""
import importlib.util
from pathlib import Path
import sys
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "fetch_library_context.py"
SPEC = importlib.util.spec_from_file_location("fetch_library_context", SCRIPT)
module = importlib.util.module_from_spec(SPEC)
sys.path.insert(0, str(SCRIPT.parent))
try:
    SPEC.loader.exec_module(module)
finally:
    sys.path.pop(0)


class LibraryContextTests(unittest.TestCase):
    def test_explicit_methods_and_kit_names(self):
        for description in ["DNA amplified by MDA.", "whole-genome amplification", "REPLI-g mini kit",
                            "multiple displacement amplification", "GenomiPhi", "MALBAC"]:
            with self.subTest(description=description):
                row = module.interpret({"library_construction_protocol": description})
                self.assertEqual(row["wga_mda_status"], "reported")
                self.assertEqual(row["whole_genome_amplification"], "true")

    def test_selection_wga_is_positive(self):
        self.assertEqual(module.interpret({"library_selection": "WGA"})["wga_mda_status"], "reported")

    def test_no_mention_does_not_mean_unamplified(self):
        for value in ["", "PCR amplification for eight cycles", "TruSeq PCR-free library"]:
            row = module.interpret({"library_construction_protocol": value})
            self.assertEqual(row["wga_mda_status"], "not_reported")
            self.assertEqual(row["whole_genome_amplification"], "")

    def test_explicit_negation_and_uncertainty(self):
        for description in ["No whole genome amplification was performed.", "without MDA", "WGA was not used"]:
            with self.subTest(description=description):
                self.assertEqual(module.interpret({"library_construction_protocol": description})["wga_mda_status"], "negated")
        self.assertEqual(module.interpret({"library_construction_protocol": "DNA may have undergone WGA"})["wga_mda_status"], "uncertain")

    def test_conflicting_fields_require_review(self):
        row = module.interpret({"library_selection": "WGA", "library_construction_protocol": "No WGA was performed."})
        self.assertEqual(row["wga_mda_status"], "uncertain")
        self.assertEqual(row["whole_genome_amplification"], "")

    def test_genomic_wgs_is_separate_from_amplification(self):
        row = module.interpret({"library_strategy": "WGS", "library_source": "GENOMIC", "library_selection": "WGA"})
        self.assertEqual(row["genomic_wgs_context"], "genomic_wgs")
        self.assertEqual(row["whole_genome_amplification"], "true")
        self.assertEqual(module.interpret({"library_source": "METAGENOMIC"})["genomic_wgs_context"], "metagenomic")

    def test_archive_identity_and_fields_required(self):
        header = "\t".join(module.FIELDS)
        row = "\t".join("SRR1" if field == "run_accession" else "" for field in module.FIELDS)
        raw = (header + "\n" + row + "\n").encode()
        self.assertEqual(module.validate_run(raw, "SRR1")["run_accession"], "SRR1")
        with self.assertRaises(ValueError):
            module.validate_run(raw, "SRR2")
        with self.assertRaises(ValueError):
            module.validate_fields(b"columnId\tdescription\nrun_accession\tx\n", "FIELDS")


if __name__ == "__main__":
    unittest.main()
