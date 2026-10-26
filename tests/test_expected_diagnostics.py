"""Expected-only measurements cannot erase native candidates or create false zeros."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from expected_diagnostics import supplement_expected_observations


class ExpectedDiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.dest = self.folder / 'sample'
        self.dest.mkdir()
        self.raw = {'version': '0.8.0', 'inputs': [{'path': 'reads.fastq.gz', 'bytes': 123}],
                    'summary': {'records': 100, 'bases': 15000, 'sampled_prefix': False},
                    'config': {'mode': 'reads', 'seed_bases': 24, 'min_copies': 3,
                               'min_identity': .85, 'min_tract': 500, 'end_window': 50},
                    'junction_options': {'end_slack': 3, 'flank_bases': 20},
                    'candidates': [{'motif': 'AACCCT', 'reads': 30, 'full_repeat_reads': 12,
                                    'repeat_bases': 3000, 'max_tract': 150}], 'rejected_candidates': []}
        self.evidence = {**copy.deepcopy(self.raw), 'analysis_kind': 'targeted_evidence',
                         'candidate_evidence': [{'motif': 'AACCCC', 'reads': 0, 'full_repeat_reads': 0,
                                                 'repeat_bases': 0, 'max_tract': 0}]}
        self.status = {'sample_id': 'sample', 'status': 'complete', 'exit_code': 0,
                       'identity': {'source_sha256': 'source-hash', 'motifs': ['TTGGGG']}}
        self.observations = [{'motif': 'TTGGGG', 'canonical_family': 'AACCCC', 'recovery': 'not_reported',
                              'reads': None, 'full_repeat_reads': None, 'candidate_status': None}]

    def save(self):
        p = self.dest / 'evidence.json'
        p.write_text(json.dumps(self.evidence))
        self.status['result_sha256'] = hashlib.sha256(p.read_bytes()).hexdigest()
        (self.dest / 'status.json').write_text(json.dumps(self.status))

    def run_join(self, complete=True):
        return supplement_expected_observations('sample', self.observations, self.folder,
                                                self.raw, 'source-hash', complete)

    def test_zero_fills_only_expectation_table_and_preserves_all_native_observations(self):
        self.save()
        before = copy.deepcopy(self.raw)
        rows, provenance = self.run_join()
        self.assertEqual(rows[0]['recovery'], 'explicit_zero')
        self.assertEqual(rows[0]['reads'], 0)
        self.assertIsNone(self.observations[0]['reads'])
        self.assertEqual(self.raw, before)
        self.assertEqual(provenance['used_expected_families'], ['AACCCC'])

    def test_absent_sample_stays_unreported_not_zero(self):
        rows, provenance = supplement_expected_observations('absent', self.observations,
                                                            self.folder, self.raw, 'source-hash')
        self.assertEqual(rows[0]['recovery'], 'not_reported')
        self.assertIsNone(provenance)

    def test_one_read_trace_is_kept_without_inserting_or_promoting_a_native_candidate(self):
        self.evidence['candidate_evidence'][0].update(reads=1, repeat_bases=24, max_tract=24)
        self.save()
        rows, _ = self.run_join()
        self.assertEqual(rows[0]['reads'], 1)
        self.assertFalse(rows[0]['original_discovery_reported'])
        self.assertEqual(len(self.raw['candidates']), 1)

    def test_source_or_result_hash_change_is_rejected(self):
        self.save()
        self.status['identity']['source_sha256'] = 'different'
        self.save()
        with self.assertRaisesRegex(ValueError, 'source hash'):
            self.run_join()
        self.status['identity']['source_sha256'] = 'source-hash'
        self.save()
        with (self.dest / 'evidence.json').open('a') as f:
            f.write(' ')
        with self.assertRaisesRegex(ValueError, 'result hash'):
            self.run_join()

    def test_wrong_sample_failed_prefix_different_input_and_detector_are_rejected(self):
        changes = [('sample', lambda: self.status.update(sample_id='other')),
                   ('failed', lambda: self.status.update(status='failed')),
                   ('prefix', lambda: self.evidence['summary'].update(sampled_prefix=True)),
                   ('input', lambda: self.evidence.update(inputs=[{'path': 'other.fastq'}])),
                   ('version', lambda: self.evidence.update(version='0.9.0')),
                   ('seed', lambda: self.evidence['config'].update(seed_bases=18)),
                   ('totals', lambda: self.evidence['summary'].update(records=99))]
        original = copy.deepcopy((self.status, self.evidence))
        for label, mutate in changes:
            with self.subTest(label=label):
                self.status, self.evidence = copy.deepcopy(original)
                mutate()
                self.save()
                with self.assertRaises(ValueError):
                    self.run_join()

    def test_rotation_and_reverse_complement_target_aliases_are_same_family(self):
        self.status['identity']['motifs'] = ['GGGGTT', 'AACCCC']
        self.save()
        rows, _ = self.run_join()
        self.assertEqual(rows[0]['reads'], 0)

    def test_missing_target_and_duplicate_family_are_rejected(self):
        self.evidence['candidate_evidence'] = []
        self.save()
        with self.assertRaisesRegex(ValueError, 'target set'):
            self.run_join()
        self.evidence['candidate_evidence'] = [{'motif': 'AACCCC', 'reads': 0,
                                               'full_repeat_reads': 0, 'repeat_bases': 0, 'max_tract': 0}] * 2
        self.save()
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            self.run_join()

    def test_diagnostic_cannot_repair_failed_native_analysis(self):
        self.save()
        with self.assertRaisesRegex(ValueError, 'incomplete source'):
            self.run_join(complete=False)

    def test_positive_native_vs_measured_zero_is_explicit_disagreement(self):
        self.raw['candidates'].append({'motif': 'AACCCC', 'reads': 12, 'full_repeat_reads': 5,
                                       'repeat_bases': 1200, 'max_tract': 150})
        self.save()
        rows, _ = self.run_join()
        self.assertEqual(rows[0]['recovery'], 'measurement_disagreement')
        self.assertEqual(rows[0]['reads'], 0)
        self.assertEqual(rows[0]['original_discovery_reads'], 12)
        self.assertEqual(self.raw['candidates'][1]['reads'], 12)

    def test_unrequested_secondary_family_is_not_added_to_expectations(self):
        self.status['identity']['motifs'].append('TTAGGG')
        self.evidence['candidate_evidence'].append(copy.deepcopy(self.raw['candidates'][0]))
        self.save()
        rows, provenance = self.run_join()
        self.assertEqual(len(rows), 1)
        self.assertEqual(provenance['unused_registered_families'], ['AACCCT'])


if __name__ == '__main__':
    unittest.main()
