"""Opt-in candidate-period flank annotations, through the real native CLI.

The fixtures are deliberately smaller than a detector seed where a matching
flank must remain undetected. This tests the new annotation after real tract
detection, without injecting tract coordinates or changing biological calls.
Set TRM_BINARY to a separately built development executable.
"""

import copy
import csv
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest


from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]
MOTIF = "TTAGGGG"
FLANK = "ACGTCCGATGCTAACGTCAGTGCATCAGACGTTCGACGTA"
CLASS_NAMES = {"insufficient_window", "neither", "lag_only", "family_only", "both"}
QUADRANTS = {"left_forward", "left_reverse", "right_forward", "right_reverse", "ambiguous"}


def rc(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canonical(sequence):
    return min(word[offset:] + word[:offset]
               for word in (sequence, rc(sequence)) for offset in range(len(word)))


class FlankPeriodCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binary = runtime_binary()
        if not cls.binary.is_file():
            raise RuntimeError(f"Build the development executable or set TRM_BINARY: {cls.binary}")

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="trm-flank-period-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.counter = 0
        self.targets = self.root / "targets.json"
        self.targets.write_text(json.dumps({"motifs": [MOTIF, "TTAGGC"]}))

    def fastq(self, sequences, name="reads.fastq"):
        path = self.root / name
        path.write_text("".join(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n"
                                for i, seq in enumerate(sequences)))
        return path

    def run_cli(self, paths, *options, mode="evidence", annotation="on", success=True):
        self.counter += 1
        prefix = self.root / f"result{self.counter}"
        command = [str(self.binary), mode, *map(str, paths), "--output", str(prefix),
                   "--threads", "1", "--min-reads", "1", "--min-identity", "1"]
        if mode == "scan":
            command += ["--prior-policy", "discovery"]
        elif mode == "evidence":
            command += ["--targets", str(self.targets)]
        if annotation is not None:
            command += ["--flank-period-annotation", annotation]
        result = subprocess.run([*command, *options], capture_output=True, text=True, timeout=30)
        if not success:
            self.assertNotEqual(result.returncode, 0, result.stderr)
            self.assertTrue(result.stderr.strip())
            self.assertFalse(prefix.with_suffix(".json").exists())
            return result
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(prefix.with_suffix(".json").read_text())

    @staticmethod
    def rows(report):
        return report.get("candidate_evidence", report.get("candidates", []))

    def row(self, report, motif=MOTIF):
        return next(row for row in self.rows(report) if row["motif"] == canonical(motif))

    @staticmethod
    def legacy(report):
        """Remove only new fields and invocation timing, retaining row order/examples."""
        report = copy.deepcopy(report)
        report.pop("command", None)
        report.pop("flank_period_annotation_method", None)
        for name in ("elapsed_seconds", "megabases_per_second"):
            report["summary"].pop(name, None)
        def strip(value):
            if isinstance(value, dict):
                value.pop("flank_period_annotation", None)
                for child in value.values():
                    strip(child)
            elif isinstance(value, list):
                for child in value:
                    strip(child)
        strip(report)
        return report

    def annotation(self, report, motif=MOTIF):
        candidate = self.row(report, motif)
        value = candidate["junction_evidence"]["flank_period_annotation"]
        self.assertEqual(value["schema_version"], 1)
        self.assertEqual(value["candidate_motif"], canonical(motif))
        self.assertEqual(value["primitive_period_bases"], len(motif))
        self.assertEqual(value["required_window_bases"], max(20, 3 * len(motif)))
        self.assertEqual(value["threshold_percent"], 85)
        self.assertTrue(value["measured_all_events"])
        self.assertEqual(set(value["classes"]), CLASS_NAMES)
        self.assertEqual(sum(group["events"] for group in value["classes"].values()),
                         value["observed_events"])
        for group in value["classes"].values():
            self.assertEqual(set(group["quadrant_counts"]), QUADRANTS)
            self.assertEqual(sum(group["quadrant_counts"].values()), group["events"])
        return value

    def assert_classes(self, value, **expected):
        self.assertEqual({key: item["events"] for key, item in value["classes"].items()},
                         {key: expected.get(key, 0) for key in CLASS_NAMES})

    def fixtures(self):
        # Main seed requires 28 bp. These 21 bp flanks cannot seed this family,
        # and their first base deliberately breaks the main tract's phase.
        matching = "GGGGTTA" * 3
        mutated = matching[:8] + "C" + matching[9:10] + "C" + matching[11:]
        return [MOTIF * 8 + flank for flank in
                (matching, "ACCTTGA" * 3, mutated, FLANK[:21], FLANK[:20])]

    def test_default_off_and_explicit_off_preserve_all_legacy_fields(self):
        path = self.fastq(self.fixtures() + [MOTIF * 100, FLANK + MOTIF * 8 + FLANK])
        for mode in ("scan", "evidence"):
            with self.subTest(mode=mode):
                default = self.run_cli([path], mode=mode, annotation=None)
                off = self.run_cli([path], mode=mode, annotation="off")
                on = self.run_cli([path], mode=mode)
                for row in self.rows(default) + self.rows(off):
                    self.assertNotIn("flank_period_annotation", row["junction_evidence"])
                self.assertEqual(self.legacy(default), self.legacy(off))
                self.assertEqual(self.legacy(off), self.legacy(on))

    def test_five_classes_through_targeted_and_discovery_detection(self):
        path = self.fastq(self.fixtures())
        targeted = self.run_cli([path])
        discovery = self.run_cli([path], mode="scan")
        for report in (targeted, discovery):
            value = self.annotation(report)
            self.assertEqual(value["observed_events"], 5)
            self.assert_classes(value, insufficient_window=1, neither=1, lag_only=1,
                                family_only=1, both=1)
            self.assertEqual(self.row(report)["junction_evidence"]["junction_reads"], 5)
        self.assertEqual(self.annotation(targeted), self.annotation(discovery))

    def test_phase_and_reverse_complement_invariance_preserves_quadrants(self):
        sequences = []
        for phase in range(len(MOTIF)):
            unit = MOTIF[phase:] + MOTIF[:phase]
            flank_unit = next(MOTIF[offset:] + MOTIF[:offset]
                              for offset in range(len(MOTIF)) if MOTIF[offset] != unit[0])
            sequence = unit * 8 + flank_unit * 3
            sequences.extend([sequence, rc(sequence)])
        report = self.run_cli([self.fastq(sequences)])
        value = self.annotation(report)
        self.assert_classes(value, both=14)
        self.assertEqual(value["classes"]["both"]["quadrant_counts"],
                         {"left_forward": 0, "left_reverse": 7,
                          "right_forward": 7, "right_reverse": 0, "ambiguous": 0})

    def test_candidate_periods_32_and_64_require_full_bounded_window(self):
        rng = random.Random(421)
        for period in (32, 64):
            with self.subTest(period=period):
                unit = "T" + "".join(rng.choices("ACGT", k=period - 1))
                flank = "A" + "".join(rng.choices("ACGT", k=3 * period - 1))
                self.targets.write_text(json.dumps({"motifs": [unit]}))
                path = self.fastq([unit * 4 + flank, unit * 4 + flank[:-1]])
                for mode in ("evidence", "scan"):
                    options = ("--max-k", "64") if mode == "scan" else ()
                    value = self.annotation(self.run_cli([path], *options, mode=mode), unit)
                    self.assertEqual(value["required_window_bases"], 3 * period)
                    self.assert_classes(value, neither=1, insufficient_window=1)

    def test_ambiguity_and_same_family_tract_stop_context_without_bridging(self):
        # The final flank base is T, not the preceding phase's G, so the
        # downstream exact tract cannot extend backward into the 20 bp gap.
        flank = FLANK[:19] + "T"
        sequences = [MOTIF * 8 + flank + "N" + FLANK,
                     MOTIF * 8 + flank + MOTIF * 4]
        report = self.run_cli([self.fastq(sequences)])
        value = self.annotation(report)
        self.assert_classes(value, insufficient_window=3)
        self.assertEqual(self.row(report)["junction_evidence"]["junction_reads"], 2)
        self.assertEqual(value["observed_events"], 3)  # Events are not unique reads.

    def test_zero_target_and_ineligible_reads_are_measured_zero(self):
        sequences = [MOTIF * 100 + FLANK[:20], MOTIF * 20,
                     FLANK + MOTIF * 8 + FLANK, MOTIF * 8 + "N" * 40]
        report = self.run_cli([self.fastq(sequences)])
        for motif in (MOTIF, "TTAGGC"):
            value = self.annotation(report, motif)
            self.assertEqual(value["observed_events"], 0)
            self.assert_classes(value)
        self.assertEqual(self.row(report, "TTAGGC")["reads"], 0)

    def test_paired_metrics_and_worker_batch_merges_are_unchanged(self):
        left = self.fixtures() * 37
        paths = [self.fastq(left, "R1.fastq"), self.fastq([rc(s) for s in left], "R2.fastq")]
        for mode in ("scan", "evidence"):
            with self.subTest(mode=mode):
                one = self.run_cli(paths, "--paired-layout", "fr", "--batch-bases", "100000", mode=mode)
                many = self.run_cli(paths, "--paired-layout", "fr", "--threads", "4",
                                    "--batch-bases", "100", mode=mode)
                off = self.run_cli(paths, "--paired-layout", "fr", "--batch-bases", "100000",
                                   mode=mode, annotation="off")
                self.assertEqual(self.legacy(one), self.legacy(off))
                self.assertEqual(self.annotation(one), self.annotation(many))
                self.assertEqual(one["evidence_qc"], many["evidence_qc"])
                self.assertEqual(one["paired_end"], many["paired_end"])
                # Examples are bounded, and can be selected differently by batch.
                comparable = lambda report: [{k: v for k, v in row.items() if k != "examples"}
                                             for row in self.rows(report)]
                self.assertEqual(comparable(one), comparable(many))

    def test_invalid_value_and_unsupported_modes_fail_before_input(self):
        missing = self.root / "missing.fastq"
        result = self.run_cli([missing], annotation="yes", success=False)
        self.assertIn("flank-period-annotation", result.stderr)
        self.assertNotIn("failed to open", result.stderr.lower())
        for mode in ("mixed", "reclassify"):
            # Supply only otherwise generally accepted options to make the
            # annotation switch itself the rejection being checked.
            output = self.root / f"rejected-{mode}"
            result = subprocess.run([str(self.binary), mode, str(missing), "--output", str(output),
                                     "--flank-period-annotation", "on"],
                                    text=True, capture_output=True, timeout=10)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("flank-period-annotation", result.stderr)
            self.assertFalse(output.with_suffix(".json").exists())

    def test_readable_report_accepts_real_native_schema_and_preserves_missingness(self):
        path = self.fastq(self.fixtures())
        generated = []
        for annotation in ("on", "off"):
            self.run_cli([path], annotation=annotation)
            generated.append(self.root / f"result{self.counter}.json")
        output = self.root / "readable"
        result = subprocess.run([sys.executable, str(ROOT / "scripts/report_flank_period.py"),
                                 *map(str, generated), "--output", str(output)],
                                text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (output / "samples.tsv").open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual(len(rows), 4)
        statuses = [row["measurement_status"] for row in rows]
        self.assertEqual(statuses.count("unmeasured"), 2)
        self.assertEqual(statuses.count("measured_zero"), 1)
        measured = next(row for row in rows if row["observed_events"] == "5")
        self.assertEqual([measured[key] for key in CLASS_NAMES], ["1"] * 5)
        summary = json.loads((output / "summary.json").read_text())
        self.assertFalse(summary["changes_native_classification"])
        self.assertFalse(summary["reads_sequence_inputs"])


if __name__ == "__main__":
    unittest.main()
