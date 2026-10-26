"""Fast producer summaries retain full-source integrity and independent audit."""
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'tests')]
import test_report_assessment as fixtures
import summarize_interpretation as summary


class SummaryProjectionTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.AssessmentReportTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.fixture.make_native_sources()
        self.fixture.run_report(native=True)

    def test_projection_and_full_json_have_identical_tables(self):
        fast = self.fixture.root / 'fast'
        full = self.fixture.root / 'full'
        sources = set((self.fixture.output / 'samples').glob('*.json'))
        original = Path.read_bytes
        def prohibit_full_json_materialization(path):
            if path in sources:
                raise AssertionError('Projection mode must stream source integrity, not materialize full JSON')
            return original(path)
        with patch.object(Path, 'read_bytes', prohibit_full_json_materialization):
            first = summary.summarize(self.fixture.output, fast)
        second = summary.summarize(self.fixture.output, full, use_projections=False)
        for name in ('samples.tsv', 'expected_families.tsv', 'review_queue.tsv'):
            self.assertEqual((fast / name).read_bytes(), (full / name).read_bytes())
        self.assertEqual(first['producer_projection_samples'], first['samples'])
        self.assertEqual(first['independently_recalculated_samples'], 0)
        self.assertEqual(second['producer_projection_samples'], 0)
        self.assertEqual(first['source_pins'], second['source_pins'])

    def test_tampered_full_source_is_rejected_even_with_valid_projection(self):
        source = self.fixture.output / 'samples/SRRknown.json'
        source.write_bytes(source.read_bytes().replace(b'TTAGGC', b'TTAGGG'))
        with self.assertRaisesRegex(ValueError, 'Assessment bytes changed'):
            summary.summarize(self.fixture.output, self.fixture.root / 'tampered-source')

    def test_projection_corruption_is_rejected_and_full_audit_remains_available(self):
        path = self.fixture.output / 'summaries/SRRknown.json'
        path.write_text('{}')
        with self.assertRaisesRegex(ValueError, 'projection bytes differ'):
            summary.summarize(self.fixture.output, self.fixture.root / 'tampered-projection')
        result = summary.summarize(self.fixture.output, self.fixture.root / 'audit', use_projections=False)
        self.assertEqual(result['independently_recalculated_samples'], result['samples'])

    def test_manifest_missing_or_wrong_cohort_cannot_silently_skip_samples(self):
        path = self.fixture.output / 'summary_projections.tsv'
        path.unlink()
        with self.assertRaisesRegex(ValueError, 'manifest checksum differs'):
            summary.summarize(self.fixture.output, self.fixture.root / 'missing-manifest')


if __name__ == '__main__':
    unittest.main()
