import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).parents[1] / "scripts/export_novel_validation.py"
SPEC = importlib.util.spec_from_file_location("exporter", SCRIPT)
exporter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(exporter)


def packed(sequence):
    value = 0
    for base in sequence:
        value = (value << 2) | "ACGT".index(base)
    return format(value, "010x")


def fixture():
    q = {key: int(key == "left_reverse") for key in exporter.QUADRANTS}
    anchors = ["ACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCA"]
    flanks = [{"signature": packed(a), "rank": "a" * 32, "sequence": a,
               "representative_flank": a + "ACGTTGCA", "events": 1,
               "quadrant_counts": q.copy(), "normalized_quadrant_counts": q.copy()} for a in anchors]
    validation = {"schema_version": 1, "measured": True, "clean_two_flank_repeat_records": 0,
                  "junction_read_signatures": {"capacity": 128, "observed_records": 2,
                      "distinct_lower_bound": 2, "saturated": False,
                      "entries": [{"signature": str(n) * 32, "records": 1, "quadrant_counts": q.copy(),
                                   "normalized_quadrant_counts": q.copy()} for n in (1, 2)]},
                  "flank_contexts": {"capacity": 128, "flank_bases": 20, "observed_events": 2,
                      "distinct_lower_bound": 2, "saturated": False, "entries": flanks}}
    native = {"motif": exporter.canonical("TTAGAC"), "reads": 12, "full_repeat_reads": 10,
              "junction_evidence": {"junction_reads": 2, "left_forward_junction_reads": 0,
                    "left_reverse_junction_reads": 2, "right_forward_junction_reads": 0,
                    "right_reverse_junction_reads": 0, "left_ambiguous_junction_reads": 0,
                    "right_ambiguous_junction_reads": 0, "validation": validation}}
    raw = {"program": "trm", "version": "0.5.0", "analysis_kind": "de_novo_discovery",
           "summary": {"records": 1000, "bases": 150000, "sampled_prefix": False},
           "candidates": [], "rejected_candidates": [native]}
    candidate = {"motif": "TTAGAC", "canonical_family": native["motif"], "reference_status": "not_in_catalogue",
                 "novel_evidence_level": "geometry_nomination", "confirmation_readiness": {"status": "measured",
                     "missing_evidence": ["minimum_distinct_junction_sequences"], "physical_flank_cluster_count": 2},
                 "directionality": {"canonical_quadrants": {"left_reverse_junction_reads": 2}}}
    assessment = {"program": "trm-assessment", "schema_version": 5, "version": "0.5.0", "sample_id": "SRR123",
                  "tax_id": "99999", "scientific_name": "Example organism", "metadata": {"sample_id": "SRR123", "sample_accession": "SAMN1"},
                  "input_summary": raw["summary"].copy(), "candidates": [candidate],
                  "provenance": {"original_discovery_result": "/source/SRR123/result.json"}}
    return assessment, raw


def export(assessment, raw, sample="SRR123"):
    raw_hash = exporter.digest(json.dumps(raw).encode())
    assessment = copy.deepcopy(assessment)
    assessment["provenance"]["original_discovery_sha256"] = raw_hash
    assessment["provenance"]["targeted_evidence_sha256"] = raw_hash
    return exporter.build_export(assessment, raw, "a" * 64, raw_hash, sample)


class ExportTests(unittest.TestCase):
    def test_valid_geometry_exports_flanks_and_required_proof(self):
        assessment, raw = fixture()
        manifest, rows, links, sequences = export(assessment, raw)
        self.assertEqual(manifest["status"], "exported")
        self.assertEqual(len(rows), 1)
        self.assertEqual(len(sequences), 2)
        self.assertEqual(rows[0]["distinct_physical_flanks_lower_bound"], 2)
        self.assertFalse(manifest["external_confirmation_verified"])
        self.assertIn("local_assembly_or_long_read_terminal_localization", rows[0]["required_proof_steps"])
        self.assertFalse(manifest["sra_search"]["network_request_performed"])
        self.assertIn("txid99999", manifest["sra_search"]["url"])
        self.assertEqual(links[0]["anchor20"], raw["rejected_candidates"][0]["junction_evidence"]["validation"]["flank_contexts"]["entries"][0]["sequence"])

    def test_only_unknown_geometry_or_stronger_is_selected(self):
        assessment, raw = fixture()
        assessment["candidates"][0]["novel_evidence_level"] = "uncatalogued_repeat"
        self.assertEqual(export(assessment, raw)[0]["target_count"], 0)
        assessment["candidates"][0]["novel_evidence_level"] = "provisional_candidate"
        assessment["candidates"][0]["reference_status"] = "established_reference"
        self.assertEqual(export(assessment, raw)[0]["target_count"], 0)

    def test_v040_is_unavailable_not_zero_evidence(self):
        assessment, raw = fixture()
        assessment.update(version="0.4.0", schema_version=4)
        manifest, rows, links, fasta = export(assessment, raw)
        self.assertEqual(manifest["status"], "not_available")
        self.assertIsNone(manifest["target_count"])
        self.assertEqual((rows, links, fasta), ([], [], {}))

    def test_v050_missing_native_validation_remains_null(self):
        assessment, raw = fixture()
        del raw["rejected_candidates"][0]["junction_evidence"]["validation"]
        manifest, rows, _, fasta = export(assessment, raw)
        self.assertEqual(rows[0]["validation_availability"], "not_available")
        self.assertIsNone(rows[0]["distinct_physical_flanks_lower_bound"])
        self.assertIsNone(rows[0]["retained_flank_context_count"])
        self.assertIsNone(rows[0]["physical_flank_cluster_count"])
        self.assertEqual(fasta, {})

    def test_source_hash_and_sample_identity_are_required(self):
        assessment, raw = fixture()
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            exporter.build_export(assessment, raw, "a", "b", "SRR123")
        with self.assertRaisesRegex(ValueError, "identity"):
            export(assessment, raw, "SRR456")
        assessment["metadata"]["sample_id"] = "SRR456"
        with self.assertRaisesRegex(ValueError, "identity"):
            export(assessment, raw)

    def test_source_path_must_bind_the_same_sample(self):
        assessment, raw = fixture()
        assessment["provenance"]["original_discovery_result"] = "/source/SRR456/result.json"
        with self.assertRaisesRegex(ValueError, "sample directory"):
            export(assessment, raw)

    def test_wrong_or_duplicate_family_is_rejected(self):
        assessment, raw = fixture()
        raw["rejected_candidates"][0]["motif"] = "TTAGGG"
        with self.assertRaisesRegex(ValueError, "absent"):
            export(assessment, raw)
        assessment, raw = fixture()
        raw["candidates"] = copy.deepcopy(raw["rejected_candidates"])
        with self.assertRaisesRegex(ValueError, "duplicate native"):
            export(assessment, raw)

    def test_assessed_motif_must_match_family(self):
        assessment, raw = fixture()
        assessment["candidates"][0]["motif"] = "TTAGGG"
        with self.assertRaisesRegex(ValueError, "family/motif"):
            export(assessment, raw)

    def test_orientation_contradiction_rejected(self):
        assessment, raw = fixture()
        assessment["candidates"][0]["directionality"]["canonical_quadrants"]["left_reverse_junction_reads"] = 9
        with self.assertRaisesRegex(ValueError, "orientation"):
            export(assessment, raw)

    def test_saturated_sketch_preserves_lower_bound(self):
        assessment, raw = fixture()
        j = raw["rejected_candidates"][0]["junction_evidence"]
        j["junction_reads"] = 4
        for name, observed in (("junction_read_signatures", "observed_records"), ("flank_contexts", "observed_events")):
            j["validation"][name].update(capacity=2, saturated=True, distinct_lower_bound=3, **{observed: 4})
        assessment["candidates"][0]["confirmation_readiness"]["physical_flank_cluster_count"] = None
        manifest, rows, _, fasta = export(assessment, raw)
        self.assertEqual(rows[0]["distinct_physical_flanks_lower_bound"], 3)
        self.assertTrue(rows[0]["flank_sketch_saturated"])
        self.assertEqual(rows[0]["retained_flank_context_count"], 2)
        self.assertIsNone(rows[0]["physical_flank_cluster_count"])
        self.assertEqual(len(fasta), 2)

    def test_bad_saturation_or_anchor_is_rejected(self):
        assessment, raw = fixture()
        sketch = raw["rejected_candidates"][0]["junction_evidence"]["validation"]["flank_contexts"]
        sketch["saturated"] = True
        with self.assertRaisesRegex(ValueError, "saturated"):
            export(assessment, raw)
        sketch["saturated"] = False
        sketch["entries"][0]["signature"] = "0" * 10
        with self.assertRaisesRegex(ValueError, "signature does not encode"):
            export(assessment, raw)

    def test_separate_targeted_measurement_cannot_be_exported_as_native(self):
        assessment, raw = fixture()
        raw_hash = exporter.digest(json.dumps(raw).encode())
        assessment["provenance"].update(original_discovery_sha256=raw_hash, targeted_evidence_sha256="different")
        with self.assertRaisesRegex(ValueError, "separate targeted"):
            exporter.build_export(assessment, raw, "a" * 64, raw_hash, "SRR123")

    def test_representative_must_begin_with_anchor(self):
        assessment, raw = fixture()
        entry = raw["rejected_candidates"][0]["junction_evidence"]["validation"]["flank_contexts"]["entries"][0]
        entry["representative_flank"] = "A" * 30
        with self.assertRaisesRegex(ValueError, "representative"):
            export(assessment, raw)

    def test_prefix_scope_requires_full_scan_followup(self):
        assessment, raw = fixture()
        raw["summary"]["sampled_prefix"] = True
        assessment["input_summary"]["sampled_prefix"] = True
        manifest, rows, _, _ = export(assessment, raw)
        self.assertEqual(manifest["input_scope"], "sampled_prefix")
        self.assertIn("complete_full_input_scan", rows[0]["required_proof_steps"])

    def test_global_fasta_dedup_preserves_associations(self):
        assessment, raw = fixture()
        other_raw = copy.deepcopy(raw["rejected_candidates"][0])
        other_raw["motif"] = exporter.canonical("TTAGGT")
        other_assessed = copy.deepcopy(assessment["candidates"][0])
        other_assessed.update(motif="TTAGGT", canonical_family=other_raw["motif"])
        raw["rejected_candidates"].append(other_raw)
        assessment["candidates"].append(other_assessed)
        manifest, rows, links, fasta = export(assessment, raw)
        self.assertEqual(len(rows), 2)
        self.assertEqual(len(links), 4)
        self.assertEqual(len(fasta), 2)

    def test_cli_writes_hash_manifest_and_protects_existing_output(self):
        assessment, raw = fixture()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw_path = root / "raw.json"
            raw_path.write_text(json.dumps(raw))
            assessment["provenance"]["original_discovery_sha256"] = exporter.digest(raw_path.read_bytes())
            apath = root / "assessment.json"
            apath.write_text(json.dumps(assessment))
            command = [sys.executable, str(SCRIPT), "--assessment", str(apath), "--raw", str(raw_path),
                       "--sample-id", "SRR123", "--output", str(root / "export")]
            process = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(process.returncode, 0, process.stderr)
            manifest = json.loads((root / "export/manifest.json").read_text())
            for name, digest in manifest["output_sha256"].items():
                self.assertEqual(exporter.digest((root / "export" / name).read_bytes()), digest)
            second = subprocess.run(command, text=True, capture_output=True)
            self.assertNotEqual(second.returncode, 0)
            self.assertIn("previous exports", second.stderr)


if __name__ == "__main__":
    unittest.main()
