"""Independent resume/cache boundary checks using small offline report fixtures."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'tests')]
import test_report_assessment as fixtures
import report_runtime
import report_assessment


class PostprocessCacheBoundaryTests(unittest.TestCase):
    def test_chunk_corruption_keeps_previous_published_table(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / 'published.tsv'
            output.write_text('previous complete table\n')
            attempt = root / 'attempt'
            attempt.mkdir()
            chunk = attempt / 'candidates.tsv'
            chunk.write_text('sample_id\nfixture\n')
            receipt = {'bundle': str(attempt / 'bundle.json'),
                       'artifacts': {str(chunk): {'sha256': '0' * 64}}}
            with self.assertRaisesRegex(ValueError, 'changed'):
                report_runtime.merge_chunks([receipt], 'candidates.tsv', output)
            self.assertEqual(output.read_text(), 'previous complete table\n')
            self.assertEqual(list(root.glob('.published.tsv.*.tmp')), [])

    def test_output_lock_refuses_second_writer(self):
        with tempfile.TemporaryDirectory() as temporary:
            output, comparison = Path(temporary) / 'current', Path(temporary) / 'comparison'
            with report_runtime.output_locks(output, comparison):
                with self.assertRaisesRegex(ValueError, 'locked'):
                    with report_runtime.output_locks(output, comparison):
                        self.fail('Second writer acquired the first writer\'s lock')

    def test_worker_memory_budget_cannot_be_oversubscribed(self):
        args = SimpleNamespace(jobs=2, memory_limit_mib=512, memory_budget_mib=512)
        with self.assertRaisesRegex(ValueError, 'memory'):
            report_runtime.validate_limits(args)

    def test_invalid_atomic_receipt_preserves_prior_complete_bytes(self):
        with tempfile.TemporaryDirectory() as temporary:
            receipt = Path(temporary) / 'receipt.json'
            report_runtime.atomic_json(receipt, {'status': 'complete'})
            before = receipt.read_bytes()
            with self.assertRaises(ValueError):
                report_runtime.atomic_json(receipt, {'status': 'complete', 'bad': float('nan')})
            self.assertEqual(receipt.read_bytes(), before)
            self.assertEqual(list(Path(temporary).glob('.*.tmp')), [])

    def test_worker_inherited_lock_survives_controller_descriptor_close(self):
        with tempfile.TemporaryDirectory() as temporary:
            output, comparison = Path(temporary) / 'current', Path(temporary) / 'comparison'
            process = None
            try:
                with report_runtime.output_locks(output, comparison) as descriptors:
                    process = subprocess.Popen([sys.executable, '-c',
                        'import sys; print("ready", flush=True); sys.stdin.read()'],
                        stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True,
                        pass_fds=tuple(descriptors))
                    self.assertEqual(process.stdout.readline().strip(), 'ready')
                with self.assertRaisesRegex(ValueError, 'locked'):
                    with report_runtime.output_locks(output, comparison):
                        self.fail('Surviving worker no longer retains inherited output locks')
            finally:
                if process is not None:
                    process.communicate(timeout=5)
            with report_runtime.output_locks(output, comparison):
                pass


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class PostprocessCacheReviewTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.AssessmentReportTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.fixture.make_native_sources()
        self.command = list(self.fixture.command())
        index = self.command.index('--evidence')
        del self.command[index:index + 2]
        self.command += ['--resume', '--jobs', '2', '--memory-limit-mib', '512',
                         '--memory-budget-mib', '1024']

    def run_report(self, *extra, success=True):
        completed = subprocess.run([*self.command, *extra], text=True, capture_output=True)
        self.assertEqual(completed.returncode == 0, success, completed.stdout + completed.stderr)
        return completed

    def document(self, sample='SRRknown'):
        return json.loads((self.fixture.output / 'samples' / (sample + '.json')).read_text())

    def test_new_original_metrics_supersede_cached_fallback(self):
        metrics = self.fixture.source / 'SRRknown/metrics.json'
        fixtures.write_json(metrics, {'wall_seconds': 7.0, 'max_rss_kb': 77})
        self.run_report()
        self.assertEqual(self.document()['execution']['original_scan_wall_seconds'], 7.0)
        original = metrics.with_name('original_scan_metrics.json')
        fixtures.write_json(original, {'wall_seconds': 23.0, 'max_rss_kb': 230})
        self.run_report()
        self.assertEqual(self.document()['execution']['original_scan_wall_seconds'], 23.0)

    def test_expected_directory_becoming_partial_cannot_reuse_cache(self):
        expected = self.fixture.root / 'expected'
        expected.mkdir()
        self.run_report('--expected-evidence', str(expected))
        (expected / 'SRRknown').mkdir()
        result = self.run_report('--expected-evidence', str(expected), success=False)
        self.assertIn('incomplete', result.stdout + result.stderr)

    def test_followup_root_removal_cannot_reuse_missing_sample_cache(self):
        followup = self.fixture.root / 'followup'
        followup.mkdir()
        self.run_report('--followup-diagnostics', str(followup))
        followup.rmdir()
        result = self.run_report('--followup-diagnostics', str(followup), success=False)
        self.assertIn('does not exist', result.stdout + result.stderr)

    def test_external_followup_failure_artifact_is_revalidated_on_hit(self):
        sample = 'SRRknown'
        source = self.fixture.source / sample / 'result.json'
        raw = json.loads(source.read_text())
        raw['version'] = '0.8.4'
        raw['inputs'] = [{'path': str(self.fixture.root / 'never_read.fastq.gz'), 'bytes': 999}]
        fixtures.write_json(source, raw)
        followup = self.fixture.root / 'followup'
        failure = self.fixture.root / 'external_failure.json'
        fixtures.write_json(failure, {'status': 'operational_failure', 'reason': 'fixture'})
        status = {'schema_version': 1, 'sample_id': sample, 'status': 'failed', 'exit_code': 1,
                  'failure_reason': 'Fixture diagnostic failure',
                  'identity': {'scientific_name': 'Fixture species', 'tax_id': '6239',
                               'native_result_sha256': digest(source), 'native_version': raw['version'],
                               'inputs': raw['inputs'], 'records': raw['summary']['records'],
                               'bases': raw['summary']['bases']},
                  'failure_artifact': {'path': str(failure), 'size_bytes': failure.stat().st_size,
                                       'sha256': digest(failure)}}
        fixtures.write_json(followup / sample / 'status.json', status)
        self.run_report('--followup-diagnostics', str(followup))
        self.assertEqual(self.document()['followup_diagnostics']['status'], 'failed')
        failure.write_text(failure.read_text() + '\n')
        result = self.run_report('--followup-diagnostics', str(followup), success=False)
        self.assertIn('hash or size mismatch', result.stdout + result.stderr)
        self.assertFalse((self.fixture.root / 'never_read.fastq.gz').exists())

    def test_assignment_edit_refreshes_whole_table_provenance_on_all_samples(self):
        self.run_report()
        rows = fixtures.read_tsv(self.fixture.assignments)
        next(row for row in rows if row['sample_id'] == 'SRRknown')['prior_source'] = 'https://example.org/updated-reference'
        fixtures.write_tsv(self.fixture.assignments, rows)
        current_hash = digest(self.fixture.assignments)
        self.run_report()
        for row in rows:
            with self.subTest(sample=row['sample_id']):
                self.assertEqual(self.document(row['sample_id'])['provenance']['assignment_sha256'], current_hash)
        modes = {row['sample_id']: row['mode'] for row in fixtures.read_tsv(
            self.fixture.output / 'interpretation_resources.tsv')}
        self.assertEqual(modes['SRRknown'], 'recomputed')
        self.assertEqual({mode for sample, mode in modes.items() if sample != 'SRRknown'}, {'refreshed'})

    def test_external_followup_symlink_retarget_cannot_reuse_warm_cache(self):
        sample = 'SRRknown'
        source = self.fixture.source / sample / 'result.json'
        raw = json.loads(source.read_text())
        raw['version'] = '0.8.4'
        raw['inputs'] = [{'path': str(self.fixture.root / 'never_read.fastq.gz'), 'bytes': 999}]
        fixtures.write_json(source, raw)
        followup = self.fixture.root / 'followup'
        original = self.fixture.root / 'failure_original.json'
        changed = self.fixture.root / 'failure_changed.json'
        fixtures.write_json(original, {'status': 'operational_failure', 'reason': 'first'})
        fixtures.write_json(changed, {'status': 'operational_failure', 'reason': 'other'})
        alias = self.fixture.root / 'failure_alias.json'
        alias.symlink_to(original)
        status = {'schema_version': 1, 'sample_id': sample, 'status': 'failed', 'exit_code': 1,
                  'failure_reason': 'Fixture diagnostic failure',
                  'identity': {'scientific_name': 'Fixture species', 'tax_id': '6239',
                               'native_result_sha256': digest(source), 'native_version': raw['version'],
                               'inputs': raw['inputs'], 'records': raw['summary']['records'],
                               'bases': raw['summary']['bases']},
                  'failure_artifact': {'path': str(alias), 'size_bytes': original.stat().st_size,
                                       'sha256': digest(original)}}
        fixtures.write_json(followup / sample / 'status.json', status)
        self.run_report('--followup-diagnostics', str(followup))
        self.assertEqual(self.document()['followup_diagnostics']['status'], 'failed')
        original_hash = digest(original)
        alias.unlink()
        alias.symlink_to(changed)
        result = self.run_report('--followup-diagnostics', str(followup), success=False)
        self.assertIn('hash or size mismatch', result.stdout + result.stderr)
        self.assertEqual(digest(original), original_hash)
        self.assertFalse((self.fixture.root / 'never_read.fastq.gz').exists())

    def test_reversed_metadata_order_keeps_deterministic_sample_and_candidate_order(self):
        self.run_report()
        before = {name: (self.fixture.output / name).read_bytes()
                  for name in ('samples.tsv', 'candidates.tsv', 'expected_family_observations.tsv')}
        rows = fixtures.read_tsv(self.fixture.metadata)
        fixtures.write_tsv(self.fixture.metadata, list(reversed(rows)))
        self.run_report()
        for name, content in before.items():
            with self.subTest(file=name):
                self.assertEqual((self.fixture.output / name).read_bytes(), content)

    def test_library_table_hash_refreshed_inside_reviewed_rows_and_provenance(self):
        table = self.fixture.root / 'reviewed.tsv'
        rows = []
        for sample in ('SRRknown', 'SRRnovel'):
            evidence = self.fixture.root / (sample + '_review.json')
            fixtures.write_json(evidence, {'sample_id': sample, 'findings': 'fixture'})
            rows.append({'sample_id': sample, 'scientific_name': 'Fixture species', 'tax_id': '6239',
                         'library_review_code': 'WGS_context_compatible', 'detailed_evidence_path': evidence.name,
                         'evidence_sha256': digest(evidence), 'review_notes': 'first'})
        fixtures.write_tsv(table, rows)
        self.run_report('--library-evidence', str(table))
        rows[0]['review_notes'] = 'updated'
        fixtures.write_tsv(table, rows)
        self.run_report('--library-evidence', str(table))
        for sample in ('SRRknown', 'SRRnovel'):
            with self.subTest(sample=sample):
                doc = self.document(sample)
                self.assertEqual(doc['provenance']['library_evidence_table_sha256'], digest(table))
                self.assertEqual(doc['reviewed_library_evidence']['library_evidence_table_sha256'], digest(table))

    def test_preferred_previous_assessment_appearing_refreshes_comparison(self):
        previous = self.fixture.root / 'previous'
        old = {'program': 'trm-assessment', 'version': '0.4.0', 'sample_id': 'SRRknown',
               'tax_id': '6239', 'classification': 'novel_trm_candidate', 'result': 'fixture old',
               'detected_motif': 'TTAGGC', 'top_candidate': {'canonical_family': 'AAGCCT'},
               'candidates': []}
        fixtures.write_json(previous / 'SRRknown.json', old)
        self.run_report('--previous-assessment', str(previous))
        nested = previous / 'samples/SRRknown.json'
        fixtures.write_json(nested, {**old, 'classification': 'known_trm_detected'})
        self.run_report('--previous-assessment', str(previous))
        rows = {row['sample_id']: row for row in fixtures.read_tsv(
            self.fixture.comparison / 'previous_assessment_comparison.tsv')}
        self.assertEqual(rows['SRRknown']['previous_classification'], 'known_trm_detected')
        self.assertEqual(rows['SRRknown']['previous_assessment_file'], str(nested.resolve()))

    def test_preflight_validation_failure_replaces_previous_completion(self):
        self.run_report()
        status_path = self.fixture.output / 'interpretation_status.json'
        previous = json.loads(status_path.read_text())
        rows = fixtures.read_tsv(self.fixture.assignments)
        fixtures.write_tsv(self.fixture.assignments, rows[:-1])
        self.run_report(success=False)
        current = json.loads(status_path.read_text())
        self.assertEqual(current['status'], 'failed')
        self.assertTrue(current.get('execution_id'))
        self.assertNotEqual(current['execution_id'], previous['execution_id'])
        failed_attempt = self.fixture.output / '.report-executions' / current['execution_id'] / 'status.json'
        self.assertEqual(json.loads(failed_attempt.read_text())['status'], 'failed')

    def test_lock_refusal_does_not_overwrite_active_writer_status(self):
        status_path = self.fixture.output / 'interpretation_status.json'
        fixtures.write_json(status_path, {'status': 'running', 'execution_id': 'active-writer'})
        before = status_path.read_bytes()
        with report_runtime.output_locks(self.fixture.output, self.fixture.comparison):
            result = self.run_report(success=False)
        self.assertIn('locked', result.stdout + result.stderr)
        self.assertEqual(status_path.read_bytes(), before)

    def test_keyboard_interrupt_marks_new_attempt_failed(self):
        status_path = self.fixture.output / 'interpretation_status.json'
        fixtures.write_json(status_path, {'status': 'complete', 'execution_id': 'prior-success'})
        with patch.object(sys, 'argv', self.command[1:]), patch.object(
                report_assessment, 'execute_report', side_effect=KeyboardInterrupt):
            with self.assertRaises(KeyboardInterrupt):
                report_assessment.main()
        current = json.loads(status_path.read_text())
        self.assertEqual(current['status'], 'failed')
        self.assertTrue(current.get('execution_id'))
        self.assertNotEqual(current['execution_id'], 'prior-success')

    def test_corrupted_cached_document_is_recomputed(self):
        self.run_report()
        receipt_path = self.fixture.output / '.report-cache/SRRknown/receipt.json'
        receipt = json.loads(receipt_path.read_text())
        cached_document = Path(receipt['document'])
        cached_document.write_text(cached_document.read_text() + '\n')
        self.run_report()
        modes = {row['sample_id']: row['mode'] for row in fixtures.read_tsv(
            self.fixture.output / 'interpretation_resources.tsv')}
        self.assertEqual(modes['SRRknown'], 'recomputed')
        self.assertEqual(self.document()['sample_id'], 'SRRknown')

    def test_cached_document_change_after_validation_cannot_seal_refresh(self):
        self.run_report()
        receipt_path = self.fixture.output / '.report-cache/SRRknown/receipt.json'
        receipt = json.loads(receipt_path.read_text())
        old_attempt = Path(receipt['bundle']).parent
        task = json.loads((old_attempt / 'task.json').read_text())
        rows = fixtures.read_tsv(self.fixture.assignments)
        next(row for row in rows if row['sample_id'] == 'SRRnovel')['prior_source'] = 'https://example.org/new-other-row'
        fixtures.write_tsv(self.fixture.assignments, rows)
        args = SimpleNamespace(**task['args'])
        for name in report_runtime.PATH_ARGS:
            setattr(args, name, Path(getattr(args, name)) if getattr(args, name) else None)
        task['globals'] = report_runtime.prepare_globals(args)
        # Run only the fixture worker in-process; never change the test runner's
        # resource limit or invoke a production interpretation process.
        task['args']['memory_limit_mib'] = None
        attempt = self.fixture.root / 'race-attempt'
        attempt.mkdir()
        task['attempt'] = str(attempt)
        task_path = attempt / 'task.json'
        fixtures.write_json(task_path, task)
        original = report_runtime.verify_artifacts
        def mutate_after_verification(value, *args, **kwargs):
            accepted = original(value, *args, **kwargs)
            if accepted:
                cached = Path(value['document'])
                document = json.loads(cached.read_text())
                document['metadata']['library_strategy'] = 'unexpected concurrent edit'
                cached.write_text(json.dumps(document))
            return accepted
        with patch.object(report_runtime, 'verify_artifacts', side_effect=mutate_after_verification):
            with self.assertRaisesRegex(ValueError, 'changed'):
                report_runtime.worker(task_path, report_assessment)
        self.assertFalse((attempt / 'result.json').exists())


if __name__ == '__main__':
    unittest.main()
