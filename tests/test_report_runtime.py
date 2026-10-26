"""End-to-end bounded reporting and semantic equivalence on offline observations."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'tests')]
import test_report_assessment as fixtures
import report_runtime
import report_assessment


class ReportRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.AssessmentReportTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.fixture.make_native_sources()
        self.command = list(self.fixture.command())
        index = self.command.index('--evidence')
        del self.command[index:index + 2]
        self.command += ['--jobs', '2', '--memory-limit-mib', '512', '--memory-budget-mib', '1024']

    def run_report(self, *extra):
        result = subprocess.run([*self.command, *extra], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return json.loads((self.fixture.output / 'interpretation_status.json').read_text())

    def resources(self):
        return {r['sample_id']: r for r in fixtures.read_tsv(self.fixture.output / 'interpretation_resources.tsv')}

    def test_warm_resume_reuses_all_and_has_separate_measured_resources(self):
        self.run_report()
        old_native = (self.fixture.output / 'runtime.tsv').read_bytes()
        documents = {p.name: p.read_bytes() for p in (self.fixture.output / 'samples').glob('*.json')}
        status = self.run_report('--resume')
        self.assertEqual(status['status'], 'complete')
        self.assertEqual(status['modes']['reused'], len(documents))
        self.assertEqual(status['modes']['recomputed'], 0)
        self.assertEqual((self.fixture.output / 'runtime.tsv').read_bytes(), old_native)
        self.assertEqual(documents, {p.name: p.read_bytes() for p in (self.fixture.output / 'samples').glob('*.json')})
        for row in self.resources().values():
            self.assertEqual(row['assessment_seconds'], '0.0')
            self.assertGreater(float(row['max_rss_kib']), 0)
            self.assertGreaterEqual(float(row['user_cpu_seconds']), 0)

    def test_selective_assignment_changes_recompute_one_and_refresh_others(self):
        self.run_report()
        rows = fixtures.read_tsv(self.fixture.assignments)
        rows[0]['prior_source'] = 'https://example.org/revised'
        changed = rows[0]['sample_id']
        fixtures.write_tsv(self.fixture.assignments, rows)
        self.run_report('--resume')
        for sample, row in self.resources().items():
            self.assertEqual(row['mode'], 'recomputed' if sample == changed else 'refreshed')

    def test_corrupted_public_sample_or_missing_receipt_field_recomputes(self):
        self.run_report()
        public = self.fixture.output / 'samples/SRRknown.json'
        public.write_text('{broken')
        self.run_report('--resume')
        self.assertEqual(self.resources()['SRRknown']['mode'], 'recomputed')
        receipt = self.fixture.output / '.report-cache/SRRknown/receipt.json'
        value = json.loads(receipt.read_text())
        value.pop('bundle')
        receipt.write_text(json.dumps(value))
        self.run_report('--resume')
        self.assertEqual(self.resources()['SRRknown']['mode'], 'recomputed')

    def test_source_symlink_retarget_during_compute_cannot_complete(self):
        source = self.fixture.source / 'SRRknown/result.json'
        target = source.with_name('original_target.json')
        source.rename(target)
        source.symlink_to(target)
        changed_target = source.with_name('other_target.json')
        raw = json.loads(target.read_text())
        raw['config']['fixture_changed_after_fingerprint'] = True
        fixtures.write_json(changed_target, raw)
        self.run_report()
        receipt = json.loads((self.fixture.output / '.report-cache/SRRknown/receipt.json').read_text())
        old_attempt = Path(receipt['bundle']).parent
        task = json.loads((old_attempt / 'task.json').read_text())
        attempt = self.fixture.root / 'symlink-race-attempt'
        attempt.mkdir()
        task['attempt'] = str(attempt)
        task['args']['resume'] = False
        task['args']['memory_limit_mib'] = None
        task_path = attempt / 'task.json'
        fixtures.write_json(task_path, task)
        original = report_assessment.sample_report
        def retarget_then_read(*args, **kwargs):
            source.unlink()
            source.symlink_to(changed_target)
            return original(*args, **kwargs)
        with patch.object(report_assessment, 'sample_report', side_effect=retarget_then_read):
            with self.assertRaisesRegex(ValueError, 'changed'):
                report_runtime.worker(task_path, report_assessment)
        self.assertFalse((attempt / 'result.json').exists())

    def test_job_count_does_not_change_documents_or_tables(self):
        self.run_report('--jobs', '1')
        names = ('samples.tsv', 'candidates.tsv', 'input_applicability.tsv', 'runtime.tsv',
                 'expected_family_observations.tsv', 'nearby_families.tsv', 'adequacy.tsv')
        before = {n: (self.fixture.output / n).read_bytes() for n in names}
        documents = {p.name: p.read_bytes() for p in (self.fixture.output / 'samples').glob('*.json')}
        self.run_report('--jobs', '2')
        self.assertEqual(before, {n: (self.fixture.output / n).read_bytes() for n in names})
        self.assertEqual(documents, {p.name: p.read_bytes() for p in (self.fixture.output / 'samples').glob('*.json')})


if __name__ == '__main__':
    unittest.main()
