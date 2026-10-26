"""Reporting contracts for optional flank annotations and their missingness."""
import copy
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import report_flank_period as report


def candidate(motif='AACCCCT', counts=None, measured=True):
    counts = counts or {}
    classes = {name: {'events': counts.get(name, 0),
                      'quadrant_counts': {q: counts.get(name, 0) if q == 'left_forward' else 0
                                          for q in report.QUADRANTS}}
               for name in report.CLASSES}
    observed = sum(counts.values())
    return {'motif': motif, 'period': len(motif), 'reads': observed,
            'status': 'repeat_candidate',
            'junction_evidence': {'junction_reads': observed,
                'flank_period_annotation': {'schema_version': 1, 'candidate_motif': motif,
                    'primitive_period_bases': len(motif), 'required_window_bases': max(20, 3 * len(motif)),
                    'threshold_percent': 85, 'observed_events': observed,
                    'measured_all_events': measured, 'classes': classes}}}


class ReporterTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.path = self.base / 'native.json'
        self.output = self.base / 'report'
        self.doc = {'sample_id': '<script>alert(1)</script>', 'analysis_kind': 'targeted_evidence',
                    'candidate_evidence': [candidate(counts={'both': 2, 'neither': 1})],
                    'rejected_candidates': []}

    def execute(self, document=None):
        self.path.write_text(json.dumps(document if document is not None else self.doc), encoding='utf-8')
        return subprocess.run([sys.executable, str(ROOT / 'scripts/report_flank_period.py'),
                               str(self.path), '--output', str(self.output)], capture_output=True, text=True)

    def test_success_source_hashes_html_escape_and_read_only_rerun(self):
        result = self.execute()
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads((self.output / 'summary.json').read_text())
        self.assertFalse(summary['changes_native_classification'])
        self.assertFalse(summary['reads_sequence_inputs'])
        self.assertEqual(summary['source_files'][0]['sha256'], hashlib.sha256(self.path.read_bytes()).hexdigest())
        for binding in summary['outputs']:
            data = Path(binding['path']).read_bytes()
            self.assertEqual(binding['size_bytes'], len(data))
            self.assertEqual(binding['sha256'], hashlib.sha256(data).hexdigest())
        page = (self.output / 'index.html').read_text()
        self.assertIn('&lt;script&gt;alert(1)&lt;/script&gt;', page)
        self.assertNotIn('<script>alert(1)</script>', page)
        self.assertIn('비반복 subtelomere', page)
        before = (self.output / 'summary.json').read_bytes()
        self.assertNotEqual(self.execute().returncode, 0)
        self.assertEqual((self.output / 'summary.json').read_bytes(), before)

    def test_missing_zero_partial_and_insufficient_are_distinct(self):
        missing = candidate('AACCCT')
        del missing['junction_evidence']['flank_period_annotation']
        zero = candidate('AACCT')
        partial = candidate('AACCCCT', counts={'insufficient_window': 1}, measured=False)
        partial['reads'] = partial['junction_evidence']['junction_reads'] = 5
        self.doc['candidate_evidence'] = [missing, zero, partial]
        self.assertEqual(self.execute().returncode, 0)
        with (self.output / 'samples.tsv').open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        self.assertEqual([r['measurement_status'] for r in rows], ['unmeasured', 'measured_zero', 'partial'])
        self.assertEqual([r['observed_events'] for r in rows], ['', '0', '1'])
        self.assertEqual([r['neither'] for r in rows], ['', '0', '0'])
        self.assertEqual(rows[2]['evaluated_events'], '0')
        self.assertEqual(rows[2]['native_junction_reads'], '5')

    def test_scan_includes_rejected_families_and_does_not_pool_files(self):
        primary = candidate('AACCCT', {'both': 1})
        rejected = candidate('AACCCCT', {'lag_only': 2})
        doc = {'candidates': [primary], 'rejected_candidates': [rejected]}
        self.path.write_text(json.dumps(doc))
        other = self.base / 'second.json'
        other.write_text(json.dumps(doc))
        summary = report.build_report([self.path, other], self.output)
        self.assertEqual(summary['candidate_result_rows'], 4)
        self.assertIsNone(summary['independent_biological_sample_count'])
        with (self.output / 'samples.tsv').open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        self.assertEqual([r['candidate_container'] for r in rows], ['candidates', 'rejected_candidates'] * 2)

    def test_corrupt_annotation_rejected_before_output_creation(self):
        mutations = [
            lambda a: a.update(schema_version=True),
            lambda a: a.update(schema_version=2),
            lambda a: a.update(candidate_motif='AACCCT'),
            lambda a: a.update(primitive_period_bases=6),
            lambda a: a.update(required_window_bases=20),
            lambda a: a.update(threshold_percent=84),
            lambda a: a.update(observed_events=True),
            lambda a: a.update(observed_events=9),
            lambda a: a.update(measured_all_events=1),
            lambda a: a['classes'].pop('both'),
            lambda a: a['classes'].update(unexpected={'events': 0}),
            lambda a: a['classes']['both'].update(events=False),
            lambda a: a['classes']['both']['quadrant_counts'].update(left_forward=1),
            lambda a: a['classes']['both']['quadrant_counts'].update(ambiguous=-1),
            lambda a: a['classes']['both']['quadrant_counts'].update(left_forward=float('nan')),
        ]
        for mutate in mutations:
            doc = copy.deepcopy(self.doc)
            mutate(doc['candidate_evidence'][0]['junction_evidence']['flank_period_annotation'])
            with self.subTest(mutation=mutate):
                result = self.execute(doc)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.output.exists())

    def test_complete_event_denominator_cannot_understate_junction_reads(self):
        row = self.doc['candidate_evidence'][0]
        row['reads'] = row['junction_evidence']['junction_reads'] = 4
        result = self.execute()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.output.exists())

    def test_duplicate_family_or_path_rejected_not_silently_double_counted(self):
        self.doc['rejected_candidates'] = copy.deepcopy(self.doc['candidate_evidence'])
        self.assertNotEqual(self.execute().returncode, 0)
        self.assertFalse(self.output.exists())
        self.doc['rejected_candidates'] = []
        self.path.write_text(json.dumps(self.doc))
        with self.assertRaisesRegex(ValueError, 'duplicate input path'):
            report.build_report([self.path, self.path], self.output)
        self.assertFalse(self.output.exists())

    def test_rotated_or_nonprimitive_motif_and_invalid_candidate_count_rejected(self):
        for field, value in [('motif', 'TTAGGGG'), ('motif', 'AACCCTAACCCT'), ('reads', True), ('period', True)]:
            with self.subTest(field=field, value=value):
                doc = copy.deepcopy(self.doc)
                doc['candidate_evidence'][0][field] = value
                self.assertNotEqual(self.execute(doc).returncode, 0)
                self.assertFalse(self.output.exists())

    def test_duplicate_json_keys_rejected(self):
        self.path.write_text('{"candidates":[],"candidates":[]}')
        with self.assertRaisesRegex(ValueError, 'duplicate JSON key'):
            report.build_report([self.path], self.output)
        self.assertFalse(self.output.exists())

    def test_no_candidate_results_is_successful_empty_report_not_one_negative_sample(self):
        summary_doc = {'candidates': [], 'rejected_candidates': []}
        self.assertEqual(self.execute(summary_doc).returncode, 0)
        summary = json.loads((self.output / 'summary.json').read_text())
        self.assertEqual(summary['candidate_result_rows'], 0)
        self.assertEqual(summary['measurement_status_counts'], {})
        self.assertIsNone(summary['independent_biological_sample_count'])


if __name__ == '__main__':
    unittest.main()
