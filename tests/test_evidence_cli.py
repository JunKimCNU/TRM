"""Whole-read targeted evidence, orientation geometry and bounded QC tests."""
import csv
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest

from archive_test_support import archive_only

from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]
FAMILY = "AACCCT"
FLANK = "ACGTCCGATGCTAACGTCAGTGCATCAGACGTTCGACGTA"


class EvidenceCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binary = runtime_binary()

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="trm-evidence-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.target = self.root / "targets.json"
        self.target.write_text(json.dumps({"motifs": ["TTAGGG", "TTAGGC"]}))
        self.counter = 0

    def fastq(self, sequences, name="reads.fastq"):
        path = self.root / name
        path.write_text("".join(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i, seq in enumerate(sequences)))
        return path

    def run_cli(self, paths, *options, mode="evidence", success=True):
        self.counter += 1
        output = self.root / f"out{self.counter}"
        command = [str(self.binary), mode, *map(str, paths), "--output", str(output),
                   "--threads", "1", "--min-reads", "1", "--min-identity", "1"]
        if mode == "evidence":
            command += ["--targets", str(self.target)]
        else:
            command += ["--prior-policy", "discovery"]
        completed = subprocess.run([*command, *options], text=True, capture_output=True, timeout=10)
        if not success:
            self.assertNotEqual(completed.returncode, 0)
            self.assertFalse(output.with_suffix(".json").exists())
            return completed
        self.assertEqual(completed.returncode, 0, completed.stderr)
        with output.with_suffix(".tsv").open() as handle:
            self.assertTrue(csv.DictReader(handle, delimiter="\t").fieldnames)
        return json.loads(output.with_suffix(".json").read_text())

    def candidate(self, report, family=FAMILY):
        return next(row for row in report.get("candidate_evidence", report.get("candidates", [])) if row["motif"] == family)

    def test_four_canonical_quadrants_and_full_repeat_separation(self):
        forward = FAMILY * 10
        reverse = "AGGGTT" * 10
        sequences = [forward + FLANK, reverse + FLANK, FLANK + forward, FLANK + reverse,
                     forward, reverse]
        result = self.run_cli([self.fastq(sequences)])
        junction = self.candidate(result)["junction_evidence"]
        self.assertEqual(junction["full_repeat_records"], 2)
        self.assertEqual(junction["junction_reads"], 4)
        for field in ("left_forward", "left_reverse", "right_forward", "right_reverse"):
            self.assertEqual(junction[field + "_junction_reads"], 1, junction)
        self.assertEqual(result["analysis_kind"], "targeted_evidence")
        self.assertEqual(result["summary"]["call"], "targeted_evidence_only")
        self.assertNotIn("candidates", result)

    def test_full_repeat_long_record_with_twenty_base_flank_is_not_junction(self):
        result = self.run_cli([self.fastq(["TTAGGG" * 100 + FLANK[:20]])])
        junction = self.candidate(result)["junction_evidence"]
        self.assertEqual(junction["full_repeat_records"], 1)
        self.assertEqual(junction["junction_reads"], 0)

    def test_internal_repeat_and_ambiguous_or_simple_flanks_are_not_junctions(self):
        sequences = [FLANK + "TTAGGG" * 10 + FLANK,
                     "TTAGGG" * 10 + "N" * 40,
                     "TTAGGG" * 10 + "A" * 40,
                     "TTAGGG" * 10 + "TTAGGC" * 10]
        report = self.run_cli([self.fastq(sequences)])
        junction = self.candidate(report)["junction_evidence"]
        self.assertEqual(junction["junction_reads"], 0)
        self.assertEqual(junction["internal_repeat_reads"], 1)
        self.assertEqual(junction["validation"]["clean_two_flank_repeat_records"], 1)
        self.assertEqual(report["evidence_qc"]["ambiguous_bases"], 40)
        self.assertEqual(report["evidence_qc"]["records_with_ambiguity"], 1)

    def test_disjoint_tracts_and_separate_records_never_make_500bp(self):
        sequences = ["TTAGGG" * 50 + "N" * 21 + "TTAGGG" * 50,
                     "TTAGGG" * 50, "CCCTAA" * 50]
        report = self.run_cli([self.fastq(sequences)])
        candidate = self.candidate(report)
        self.assertEqual(candidate["reads"], 3)
        self.assertEqual(candidate["repeat_bases"], 1200)
        self.assertEqual(candidate["max_tract"], 300)
        self.assertEqual(candidate["long_tract_reads"], 0)

    def test_late_minimum_seed_phases_match_full_discovery(self):
        sequences = ["N" * (160 + shift) + "TTAGGG" * 4 + "N" * 17 for shift in range(6)]
        path = self.fastq(sequences)
        targeted = self.run_cli([path], "--batch-bases", "10")
        exhaustive = self.run_cli([path], mode="scan")
        for field in ("reads", "repeat_bases", "max_tract", "forward_reads", "reverse_reads", "junction_evidence"):
            self.assertEqual(self.candidate(targeted)[field], self.candidate(exhaustive)[field], field)
        self.assertEqual(self.candidate(targeted)["reads"], 6)
        self.assertGreaterEqual(self.candidate(targeted)["examples"][0]["start"], 160)

    def test_entire_input_files_are_used_and_explicit_prefix_is_marked(self):
        paths = [self.fastq(["ACGT" * 25], "R1.fastq"), self.fastq(["TTAGGG" * 25], "R2.fastq")]
        report = self.run_cli(paths)
        prefix = self.run_cli(paths, "--max-records", "1")
        self.assertEqual(report["summary"]["records"], 2)
        self.assertFalse(report["summary"]["sampled_prefix"])
        self.assertEqual(self.candidate(report)["reads"], 1)
        self.assertTrue(prefix["summary"]["sampled_prefix"])
        self.assertEqual(self.candidate(prefix)["reads"], 0)

    def test_all_targets_including_zero_hits_and_saved_candidates_are_retained(self):
        self.target.write_text(json.dumps({"candidates": [{"motif": "TTAGGG"}, {"motif": "CCCTAA"},
                                                         {"motif": "TTAGGC"}]}))
        report = self.run_cli([self.fastq(["TTAGGG" * 25])])
        self.assertEqual(len(report["candidate_evidence"]), 2)
        self.assertEqual(self.candidate(report, "AAGCCT")["reads"], 0)
        self.assertEqual(self.candidate(report, "AAGCCT")["junction_evidence"]["junction_reads"], 0)

    def test_histogram_exact_short_and_bounded_long_bins_cover_all_records(self):
        sequences = ["N" * length for length in (20, 101, 1000, 1001, 1023, 1024, 2048)]
        report = self.run_cli([self.fastq(sequences)])
        qc = report["evidence_qc"]
        self.assertEqual(qc["records"], report["summary"]["records"])
        self.assertEqual(qc["bases"], sum(map(len, sequences)))
        self.assertEqual(qc["ambiguous_bases"], qc["bases"])
        self.assertEqual(qc["max_record_bases"], report["summary"]["max_record_bases"])
        self.assertEqual(sum(row["records"] for row in qc["length_histogram"]), len(sequences))
        self.assertIn({"lower": 1001, "upper": 1023, "records": 2}, qc["length_histogram"])
        self.assertIn({"lower": 1024, "upper": 2047, "records": 1}, qc["length_histogram"])

    def test_invalid_target_file_fails_before_opening_missing_input(self):
        for value in ({"motifs": []}, {"motifs": ["TTNGGG"]}, {"motifs": [""]}, {"motifs": [5]}, {}):
            with self.subTest(value=value):
                self.target.write_text(json.dumps(value))
                completed = self.run_cli([self.root / "missing.fastq"], success=False)
                self.assertNotIn("missing.fastq", completed.stderr)

    def test_target_limit_and_junction_options_validate(self):
        rng = random.Random(71)
        self.target.write_text(json.dumps({"motifs": ["".join(rng.choices("ACGT", k=24)) for _ in range(257)]}))
        completed = self.run_cli([self.root / "missing.fastq"], success=False)
        self.assertIn("target limit", completed.stderr)
        self.target.write_text(json.dumps({"motifs": ["TTAGGG"]}))
        for options in (("--junction-flank-bases", "19"), ("--junction-end-slack", "1001")):
            self.run_cli([self.root / "missing.fastq"], *options, success=False)

    def test_one_and_four_workers_preserve_evidence_and_qc(self):
        sequences = ["TTAGGG" * 10 + FLANK, "CCCTAA" * 10, FLANK + "TTAGGC" * 10] * 100
        path = self.fastq(sequences)
        one = self.run_cli([path], "--threads", "1", "--batch-bases", "100")
        many = self.run_cli([path], "--threads", "4", "--batch-bases", "100")
        self.assertEqual(one["evidence_qc"], many["evidence_qc"])
        left = {c["motif"]: {k: v for k, v in c.items() if k != "examples"} for c in one["candidate_evidence"]}
        right = {c["motif"]: {k: v for k, v in c.items() if k != "examples"} for c in many["candidate_evidence"]}
        self.assertEqual(left, right)

    def test_validation_collapses_rc_junctions_and_exports_same_physical_flank(self):
        sequence = FAMILY * 10 + FLANK
        reverse = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        result = self.run_cli([self.fastq([sequence, reverse, sequence])])
        evidence = self.candidate(result)['junction_evidence']
        validation = evidence['validation']
        self.assertTrue(validation['measured'])
        self.assertEqual(validation['schema_version'], 1)
        reads = validation['junction_read_signatures']
        self.assertEqual(reads['observed_records'], evidence['junction_reads'])
        self.assertEqual(reads['observed_records'], 3)
        self.assertEqual(reads['distinct_lower_bound'], 1)
        self.assertFalse(reads['saturated'])
        self.assertEqual(reads['entries'][0]['records'], 3)
        self.assertEqual(sum(reads['entries'][0]['normalized_quadrant_counts'].values()), 3)
        contexts = validation['flank_contexts']
        self.assertEqual(contexts['distinct_lower_bound'], 1)
        self.assertEqual(contexts['entries'][0]['events'], 3)
        entry = contexts['entries'][0]
        self.assertEqual(len(entry['sequence']), 20)
        self.assertTrue(entry['representative_flank'].startswith(entry['sequence']))
        self.assertLessEqual(len(entry['representative_flank']), 80)
        self.assertEqual(int(entry['signature'], 16), sum('ACGT'.index(base) << (2 * (19 - index))
                         for index, base in enumerate(entry['sequence'])))
        absent = self.candidate(result, 'AAGCCT')['junction_evidence']['validation']
        self.assertTrue(absent['measured'])
        self.assertEqual(absent['junction_read_signatures']['observed_records'], 0)

    def test_full_repeat_duplicates_do_not_enter_validation_sketches(self):
        result = self.run_cli([self.fastq([FAMILY * 20] * 100)])
        evidence = self.candidate(result)['junction_evidence']
        self.assertEqual(evidence['full_repeat_records'], 100)
        validation = evidence['validation']
        self.assertEqual(validation['junction_read_signatures']['entries'], [])
        self.assertEqual(validation['flank_contexts']['entries'], [])
        self.assertEqual(validation['clean_two_flank_repeat_records'], 0)

    def test_saturated_sketch_counts_match_across_worker_counts_and_batch_sizes(self):
        rng = random.Random(9457)
        originals = [FAMILY * 10 + ''.join(rng.choices('ACGT', k=60)) for _ in range(220)]
        sequences = []
        for sequence in originals:
            reverse = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
            sequences.extend([sequence] * 3 + [reverse] * 2)
        rng.shuffle(sequences)
        path = self.fastq(sequences)
        one = self.run_cli([path], '--threads', '1', '--batch-bases', '100000')
        many = self.run_cli([path], '--threads', '4', '--batch-bases', '100')
        left = self.candidate(one)['junction_evidence']
        self.assertEqual(left, self.candidate(many)['junction_evidence'])
        reads = left['validation']['junction_read_signatures']
        self.assertTrue(reads['saturated'])
        self.assertEqual(len(reads['entries']), 128)
        self.assertEqual(reads['distinct_lower_bound'], 129)
        self.assertTrue(all(entry['records'] == 5 for entry in reads['entries']))
        self.assertEqual(reads['observed_records'], left['junction_reads'])
        self.assertTrue(left['validation']['flank_contexts']['saturated'])

    @archive_only
    def test_new_validation_preserves_all_v040_candidate_fields(self):
        frozen = ROOT / 'results/release_v040_20260921/trm'
        if not frozen.is_file():
            self.skipTest('Optional preserved v0.4 binary unavailable')
        sequences = [FAMILY * 10 + FLANK, FLANK + 'TTAGGG' * 10,
                     FLANK + FAMILY * 10 + FLANK, FAMILY * 20,
                     'TTAGGC' * 12 + FLANK, 'ATC' * 25 + FLANK,
                     'GCTCAG' * 12 + 'N' * 50] * 3
        path = self.fastq(sequences)
        updated = self.run_cli([path], mode='scan')
        current = self.binary
        try:
            self.binary = frozen
            previous = self.run_cli([path], mode='scan')
        finally:
            self.binary = current
        for candidate in updated['candidates']:
            candidate['junction_evidence'].pop('validation')
        self.assertEqual(updated['candidates'], previous['candidates'])


if __name__ == "__main__":
    unittest.main()
