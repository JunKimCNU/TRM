"""Mixed structural reports must not manufacture absence or decision upgrades."""
import copy
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    "report_mixed_arrays", Path(__file__).resolve().parents[1] / "scripts/report_mixed_arrays.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def native_fixture():
    group = {key: 0 for key in MODULE.GROUP_COUNTS}
    group.update(name="mixed", units=["TTGGGTTAC", "ACTTGGGTTAC"], records=8,
                 mixed_records=5, pure_only_records=3, total_mixed_arrays=6,
                 full_mixed_records=2, entire_mixed_records=1, repeat_bases=900,
                 max_array_bp=510, max_any_array_bp=600, long_array_records=1,
                 g_right_boundary_records=2, c_left_boundary_records=1,
                 g_left_boundary_records=1, clean_two_flank_records=1,
                 internal_array_records=2, right_g_end_records=2,
                 left_c_end_records=1,
                 examples=[{"record_id": "read<script>", "record_length": 150,
                            "input_index": 0, "start": 0, "end": 150, "strand": True,
                            "unit_counts": [8, 7], "unit_counts_max": [8, 7],
                            "core_start": 4, "core_end": 149, "partial_edge_bases": 5,
                            "segmentation_ambiguous": False}])
    return {"analysis_kind": "targeted_mixed_arrays", "version": "0.7.0",
            "summary": {"records": 10, "bases": 2000},
            "inputs": [{"path": "/data/test.fastq.gz", "bytes": 456, "modified_unix_seconds": 123}],
            "mixed_arrays": {"mode": "reads", "settings": {"full_fraction": 0.9,
                "long_array_bp": 500, "exact_only": True}, "groups": [group]}}


def previous_fixture():
    return {"version": "0.5.0", "summary": {"records": 10, "bases": 2000},
            "inputs": copy.deepcopy(native_fixture()["inputs"]), "candidates": [],
            "rejected_candidates": [{"motif": "TTGGGTTAC", "reads": 4,
                "full_repeat_reads": 0, "repeat_bases": 130, "max_tract": 40}]}


class MixedReportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)
        self.source = self.path / "mixed.json"
        self.previous = self.path / "previous.json"
        self.source.write_text(json.dumps(native_fixture()))
        self.previous.write_text(json.dumps(previous_fixture()))

    def test_positive_structure_is_separate_from_biological_call(self):
        report = MODULE.build_report(self.source, self.previous, "control")
        self.assertFalse(report["classification_changed"])
        self.assertFalse(report["data_adequacy_changed"])
        self.assertIsNone(report["biological_call"])
        self.assertEqual(report["groups"][0]["observation"], "exact_mixed_array_observed")
        self.assertEqual(report["groups"][0]["mixed_records"], 5)
        self.assertEqual(report["groups"][0]["total_mixed_arrays"], 6)
        self.assertFalse(report["previous"]["input_check"]["sequence_rehashed"])

    def test_missing_component_is_not_fabricated_zero(self):
        report = MODULE.build_report(self.source, self.previous)
        rows = {r["unit"]: r for r in report["previous_components"]}
        self.assertEqual(rows["TTGGGTTAC"]["previous_reads"], 4)
        self.assertEqual(rows["ACTTGGGTTAC"]["previous_observation"], "not_reported")
        self.assertIsNone(rows["ACTTGGGTTAC"]["previous_reads"])

    def test_explicit_targeted_zero_remains_distinct(self):
        previous = previous_fixture()
        previous["rejected_candidates"][0]["reads"] = 0
        self.previous.write_text(json.dumps(previous))
        row = MODULE.build_report(self.source, self.previous)["previous_components"][0]
        self.assertEqual(row["previous_observation"], "explicit_zero")
        self.assertEqual(row["previous_reads"], 0)

    def test_rotation_and_reverse_complement_previous_match(self):
        previous = previous_fixture()
        previous["rejected_candidates"][0]["motif"] = "CCCAAGTAA"
        self.assertEqual(MODULE.canonical("CCCAAGTAA"), MODULE.canonical("TTGGGTTAC"))
        self.previous.write_text(json.dumps(previous))
        self.assertEqual(MODULE.build_report(self.source, self.previous)["previous_components"][0]["previous_reads"], 4)

    def test_missing_measurement_rejected_not_defaulted_zero(self):
        for key in ("mixed_records", "full_mixed_records", "c_right_boundary_records"):
            with self.subTest(key=key):
                native = native_fixture()
                del native["mixed_arrays"]["groups"][0][key]
                with self.assertRaises(ValueError):
                    MODULE.validate_source(native)

    def test_record_subsets_and_mixed_union_bounds_enforced(self):
        for key, bad in (("entire_mixed_records", 3), ("mixed_records", 9),
                         ("total_mixed_arrays", 4), ("repeat_bases", 2001),
                         ("long_array_records", 6), ("g_right_boundary_records", 6)):
            with self.subTest(key=key):
                native = native_fixture()
                native["mixed_arrays"]["groups"][0][key] = bad
                with self.assertRaises(ValueError):
                    MODULE.validate_source(native)

    def test_negative_boolean_fractional_counts_rejected(self):
        for bad in (-1, True, 1.5, None):
            with self.subTest(bad=bad):
                native = native_fixture()
                native["mixed_arrays"]["groups"][0]["mixed_records"] = bad
                with self.assertRaises(ValueError):
                    MODULE.validate_source(native)

    def test_nonqualifying_records_are_not_reported_as_pure(self):
        native = native_fixture()
        group = native["mixed_arrays"]["groups"][0]
        group.update(pure_only_records=1, nonqualifying_records=2)
        self.source.write_text(json.dumps(native))
        row = MODULE.build_report(self.source)["groups"][0]
        self.assertEqual(row["pure_only_records"], 1)
        self.assertEqual(row["nonqualifying_records"], 2)
        self.assertEqual(row["records"], 8)

    def test_general_ambiguity_counts_can_include_nonmixed_records(self):
        native = native_fixture()
        group = native["mixed_arrays"]["groups"][0]
        group["segmentation_ambiguous_records"] = 7
        MODULE.validate_source(native)
        group["segmentation_ambiguous_records"] = 9
        with self.assertRaisesRegex(ValueError, "any-array"):
            MODULE.validate_source(native)

    def test_record_edges_report_censoring_not_complete_telomere(self):
        example = MODULE.build_report(self.source)["examples"][0]
        self.assertTrue(example["biological_length_censored"])
        self.assertTrue(example["touches_record_start"])
        self.assertTrue(example["touches_record_end"])
        self.assertEqual(example["span_bp"], 150)
        self.assertEqual(example["strand"], "G")
        self.assertEqual(example["native_example"]["partial_edge_bases"], 5)

    def test_invalid_example_coordinates_rejected(self):
        native = native_fixture()
        native["mixed_arrays"]["groups"][0]["examples"][0]["end"] = 151
        self.source.write_text(json.dumps(native))
        with self.assertRaisesRegex(ValueError, "coordinates"):
            MODULE.build_report(self.source)

    def test_previous_comparison_refuses_different_input_and_prefix(self):
        for target, key, value in (("inputs", "bytes", 457), ("inputs", "path", "/data/other.fq"),
                                   ("inputs", "modified_unix_seconds", 124),
                                   ("summary", "records", 9), ("summary", "bases", 1999)):
            with self.subTest(key=key):
                previous = previous_fixture()
                (previous[target][0] if target == "inputs" else previous[target])[key] = value
                with self.assertRaisesRegex(ValueError, "different"):
                    MODULE.check_previous_inputs(native_fixture(), previous)

    def test_duplicate_component_family_not_arbitrarily_selected(self):
        previous = previous_fixture()
        previous["candidates"] = copy.deepcopy(previous["rejected_candidates"])
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            MODULE.component_rows(native_fixture()["mixed_arrays"]["groups"], previous, "id", "hash")

    def test_zero_mixed_is_pattern_specific_and_pure_evidence_retained(self):
        native = native_fixture()
        group = native["mixed_arrays"]["groups"][0]
        group.update({key: 0 for key in MODULE.GROUP_COUNTS})
        group.update(records=3, pure_only_records=3, max_any_array_bp=600, examples=[])
        self.source.write_text(json.dumps(native))
        report = MODULE.build_report(self.source)
        row = report["groups"][0]
        self.assertEqual(row["observation"], "supplied_exact_mixed_pattern_not_observed")
        self.assertEqual(row["pure_only_records"], 3)
        self.assertIsNone(report["biological_call"])

    def test_outputs_have_blank_missing_counts_and_escaped_html(self):
        report = MODULE.build_report(self.source, self.previous, "<script>alert(1)</script>")
        output = self.path / "out"
        MODULE.write_report(report, output)
        rendered = (output / "index.html").read_text()
        self.assertNotIn("<script>", rendered)
        self.assertIn("&lt;script&gt;", rendered)
        self.assertIn("최대 관측 span", rendered)
        self.assertIn("기존 v0.6 TRM 판정", rendered)
        with (output / "previous_components.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[1]["previous_reads"], "")
        self.assertEqual(rows[0]["previous_full_repeat_reads"], "0")
        self.assertEqual(json.loads((output / "report.json").read_text())["biological_call"], None)

    def test_source_cannot_be_overwritten_by_report(self):
        source = self.path / "report.json"
        source.write_text(json.dumps(native_fixture()))
        with self.assertRaisesRegex(ValueError, "overwrite"):
            MODULE.write_report(MODULE.build_report(source), self.path)

    def test_duplicate_or_reverse_complement_units_rejected(self):
        native = native_fixture()
        native["mixed_arrays"]["groups"][0]["units"] = ["TTGGGTTAC", "CCCAAGTAA"]
        with self.assertRaisesRegex(ValueError, "distinct primitive"):
            MODULE.validate_source(native)


if __name__ == "__main__":
    unittest.main()
