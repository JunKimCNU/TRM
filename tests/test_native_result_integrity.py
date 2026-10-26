"""Failed reruns and changed completion receipts must not become success caches."""
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'tests')]
import test_report_assessment as fixtures
from test_report_assessment import write_json, digest, read_tsv


class NativeResultIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.AssessmentReportTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.fixture.make_native_sources()
        self.sample = 'SRRknown'
        self.source = self.fixture.source / self.sample / 'result.json'
        self.state = self.source.with_name('state.json')

    def run_report(self, **kwargs):
        return self.fixture.run_report(native=True, **kwargs)

    def test_failed_or_running_rerun_does_not_interpret_retained_result_as_complete(self):
        for status in ('failed', 'running'):
            with self.subTest(status=status):
                write_json(self.state, {'status': status, 'result_sha256': digest(self.source)})
                self.run_report()
                document = self.fixture.document(self.sample)
                self.assertFalse(document['analysis_complete'])
                self.assertNotEqual(document['classification'], 'known_trm_detected')
                self.assertEqual(document['provenance']['targeted_scan_status'], 'source_not_complete')
                self.assertEqual(document['data_adequacy']['status'], 'observed_qc_failure')

    def test_completed_state_requires_matching_bytes_and_sample(self):
        for state in ({'status': 'complete', 'result_sha256': '0' * 64},
                      {'status': 'complete'},
                      {'status': 'complete', 'result_sha256': digest(self.source), 'sample_id': 'other'}):
            with self.subTest(state=state):
                write_json(self.state, state)
                self.run_report(success=False)
        write_json(self.state, {'status': 'complete', 'result_sha256': digest(self.source),
                                'sample_id': self.sample})
        self.run_report()
        document = self.fixture.document(self.sample)
        self.assertTrue(document['analysis_complete'])
        self.assertTrue(document['provenance']['native_result_validation']['execution_receipt_verified'])

    def test_standalone_native_remains_usable_without_fabricating_execution_receipt(self):
        self.run_report()
        document = self.fixture.document(self.sample)
        self.assertTrue(document['analysis_complete'])
        self.assertFalse(document['provenance']['native_result_validation']['execution_receipt_verified'])
        self.assertEqual(document['provenance']['native_result_validation']['status'], 'standalone_state_unavailable')

    def test_state_change_invalidates_warm_cache_even_if_result_bytes_stay_same(self):
        write_json(self.state, {'status': 'complete', 'result_sha256': digest(self.source)})
        self.run_report()
        self.run_report(extra=('--resume',))
        before = {r['sample_id']: r['mode'] for r in read_tsv(self.fixture.output / 'interpretation_resources.tsv')}
        self.assertEqual(before[self.sample], 'reused')
        write_json(self.state, {'status': 'failed', 'error': 'replacement input truncated'})
        self.run_report(extra=('--resume',))
        after = {r['sample_id']: r['mode'] for r in read_tsv(self.fixture.output / 'interpretation_resources.tsv')}
        self.assertEqual(after[self.sample], 'recomputed')
        self.assertFalse(self.fixture.document(self.sample)['analysis_complete'])
        self.assertEqual(after['SRRnovel'], 'reused')


if __name__ == '__main__':
    unittest.main()
