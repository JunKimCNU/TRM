"""A partial or mismatched archive must never enter the analysis queue."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

path = Path(__file__).resolve().parents[1]/'scripts/watch_test_panel.py'
spec = importlib.util.spec_from_file_location('watch_test_panel', path)
watch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(watch)


class WatcherTests(unittest.TestCase):
    def row(self):
        return dict(sample_id='SRR1', download_urls='https://archive/r1.gz;https://archive/r2.gz',
            fastq_bytes='4;4', fastq_md5=';'.join(hashlib.md5(b'ACGT').hexdigest() for _ in range(2)))

    def test_complete_sized_part_does_not_enter_queue(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'SRR1').mkdir()
            (root/'SRR1/r1.gz').write_bytes(b'ACGT')
            (root/'SRR1/r2.gz.part').write_bytes(b'ACGT')
            paths, retained = watch.inspect_download(self.row(), [root], {})
            self.assertEqual(paths, [])
            self.assertEqual(retained, 8)

    def test_size_alone_is_not_integrity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'SRR1').mkdir()
            (root/'SRR1/r1.gz').write_bytes(b'ACGT')
            (root/'SRR1/r2.gz').write_bytes(b'TTTT')
            cache = {}
            self.assertEqual(watch.inspect_download(self.row(), [root], cache)[0], [])
            (root/'SRR1/r2.gz').write_bytes(b'ACGT')
            self.assertEqual(len(watch.inspect_download(self.row(), [root], cache)[0]), 2)

    def test_duplicate_roots_do_not_inflate_retained_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'SRR1').mkdir()
            for name in ('r1.gz', 'r2.gz'):
                (root/'SRR1'/name).write_bytes(b'ACGT')
            paths, retained = watch.inspect_download(self.row(), [root, root], {})
            self.assertEqual(len(paths), 2)
            self.assertEqual(retained, 8)

    def test_bad_hash_is_rechecked_to_allow_same_size_replacement(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); (root/'SRR1').mkdir()
            (root/'SRR1/r1.gz').write_bytes(b'TTTT')
            cache = {}
            with patch.object(watch, 'sha', wraps=watch.sha) as digest:
                watch.inspect_download(self.row(), [root], cache)
                watch.inspect_download(self.row(), [root], cache)
            self.assertEqual(digest.call_count, 2)

    def test_file_changed_during_hash_is_not_verified(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); (root/'SRR1').mkdir()
            for name in ('r1.gz', 'r2.gz'):
                (root/'SRR1'/name).write_bytes(b'ACGT')
            original_sha = watch.sha
            def changing_sha(path, algorithm='sha256'):
                result = original_sha(path, algorithm)
                if path.name == 'r1.gz':
                    before = path.stat()
                    path.write_bytes(b'TTTT')
                    os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns + 1000000000))
                return result
            with patch.object(watch, 'sha', side_effect=changing_sha):
                paths, _ = watch.inspect_download(self.row(), [root], {})
            self.assertEqual(paths, [])

    def test_panel_metrics_sum_per_sample_not_latest_batch(self):
        values = {
            'a': ({'summary': {'records': 10, 'bases': 100}}, {'input_compressed_bytes': 1000000, 'wall_seconds': 2,
                'user_cpu_seconds': 3, 'system_cpu_seconds': 1, 'max_rss_kb': 100}),
            'b': ({'summary': {'records': 20, 'bases': 200}}, {'input_compressed_bytes': 2000000, 'wall_seconds': 4,
                'user_cpu_seconds': 6, 'system_cpu_seconds': 2, 'max_rss_kb': 200})}
        metrics = watch.panel_metrics(values)
        self.assertEqual(metrics['completed_samples'], 2)
        self.assertEqual(metrics['total_records'], 30)
        self.assertEqual(metrics['sum_scan_wall_seconds'], 6)
        self.assertEqual(metrics['scan_only_compressed_MBps'], .5)
        self.assertEqual(metrics['sum_user_cpu_seconds'], 9)
        self.assertEqual(metrics['max_sample_rss_kb'], 200)


class WatcherLifecycleTests(unittest.TestCase):
    def fixture(self, root, samples=('SRR1',)):
        (root/'data').mkdir()
        for filename in ('telomere_motifs.tsv', 'taxon_motif_references.tsv'):
            (root/'data'/filename).write_text('fixture\n')
        (root/'scripts').mkdir()
        for filename in ('batch.py', 'report_assessment.py', 'trm_assessment.py', 'assess_adequacy.py', 'evidence_batch.py'):
            (root/'scripts'/filename).write_text('# fixture\n')
        (root/'trm').write_bytes(b'frozen executable')
        rows = []
        for sid in samples:
            archive = root/'downloads'/sid; archive.mkdir(parents=True)
            for name in ('r1.gz', 'r2.gz'):
                (archive/name).write_bytes(b'ACGT')
            rows.append(dict(sample_id=sid, scientific_name='Test species',
                download_urls='https://archive/r1.gz;https://archive/r2.gz', fastq_bytes='4;4',
                fastq_md5=';'.join(hashlib.md5(b'ACGT').hexdigest() for _ in range(2)),
                total_fastq_bytes='8', expected_architecture='simple_repeat', expected_motifs='TTAGGG'))
        watch.write_table(root/'panel.tsv', rows)
        watch.write_table(root/'metadata.tsv', [dict(sample_id=sid, scientific_name='Test species') for sid in samples])
        watch.write_table(root/'assignments.tsv', [dict(sample_id=sid) for sid in samples])
        return ['watch', '--manifest', str(root/'panel.tsv'), '--metadata', str(root/'metadata.tsv'),
                '--assignments', str(root/'assignments.tsv'), '--download-roots', str(root/'downloads'),
                '--binary', str(root/'trm'), '--output', str(root/'out'), '--poll-seconds', '1', '--max-wait-seconds', '1']

    def scan_output(self, root, command):
        row = watch.table(Path(command[2]))[0]
        sid = row['sample_id']
        paths = [row[key] for key in ('path1', 'path2') if row.get(key)]
        directory = root/'out/scans'/sid; directory.mkdir(parents=True)
        raw = {'summary': {'records': 2, 'bases': 8, 'sampled_prefix': False}, 'candidates': []}
        metrics = {'exit_code': 0, 'input_compressed_bytes': 8, 'wall_seconds': 2.0,
                   'user_cpu_seconds': 3.0, 'system_cpu_seconds': .1, 'max_rss_kb': 10}
        (directory/'result.json').write_text(json.dumps(raw))
        (directory/'metrics.json').write_text(json.dumps(metrics))
        state = {'status': 'complete', 'result_sha256': watch.sha(directory/'result.json'), 'metrics': metrics,
                 'fingerprint': {'inputs': [{key: row[key] for key in ('path', 'size', 'mtime_ns')}
                                            for row in watch.source_fingerprint(paths)],
                                 'executable': {'sha256': watch.sha(root/'trm')},
                                 'options': ['--threads', '4', '--prior-policy', 'canonical']}}
        (directory/'state.json').write_text(json.dumps(state))
        (root/'out/scans/summary.json').write_text(json.dumps({'total_samples': 1, 'sample': sid}))

    def report_output(self, root, command):
        rows = watch.table(Path(command[command.index('--metadata') + 1]))
        assessed = root/'out/assessed'; assessed.mkdir(exist_ok=True)
        comparison = root/'out/raw_vs_assessed'; comparison.mkdir(exist_ok=True)
        (assessed/'index.html').write_text('fixture report')
        watch.write_table(assessed/'samples.tsv', [dict(sample_id=row['sample_id']) for row in rows])
        (assessed/'candidates.tsv').write_text('sample_id\n')
        (comparison/'comparison.tsv').write_text('sample_id\n')
        (assessed/'metrics.json').write_text(json.dumps({'samples': len(rows), 'files': {
            'samples.tsv': watch.sha(assessed/'samples.tsv'), 'candidates.tsv': watch.sha(assessed/'candidates.tsv'),
            'comparison.tsv': watch.sha(comparison/'comparison.tsv')}}))

    def runner(self, root, report_failure=False, stop_after_scan=False):
        def execute(command, log):
            if command[1].endswith('batch.py'):
                self.scan_output(root, command)
                if stop_after_scan:
                    (root/'out/STOP_REQUESTED').touch()
            elif report_failure:
                return {'command': command, 'returncode': 1, 'wall_seconds': .01}
            else:
                self.report_output(root, command)
            return {'command': command, 'returncode': 0, 'wall_seconds': .01}
        return execute

    def test_success_preserves_each_batch_and_cumulative_metrics(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); args = self.fixture(root, ('SRR1', 'SRR2'))
            with patch.object(watch, 'ROOT', root), patch.object(sys, 'argv', args), patch.object(watch, 'execute', side_effect=self.runner(root)):
                self.assertEqual(watch.main(), 0)
            self.assertTrue(json.loads((root/'out/status.json').read_text())['complete'])
            metrics = json.loads((root/'out/panel_metrics.json').read_text())
            self.assertEqual(metrics['completed_samples'], 2)
            self.assertEqual(metrics['completed_input_compressed_bytes'], 16)
            self.assertEqual(len(list((root/'out/scan_invocations').glob('*.json'))), 2)

    def test_report_failure_is_not_complete(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); args = self.fixture(root)
            with patch.object(watch, 'ROOT', root), patch.object(sys, 'argv', args), patch.object(watch, 'execute', side_effect=self.runner(root, report_failure=True)):
                self.assertEqual(watch.main(), 1)
            status = json.loads((root/'out/status.json').read_text())
            self.assertFalse(status['complete'])
            self.assertFalse(status['report_ready'])
            self.assertEqual(status['report_state'], 'failed')
            self.assertNotIn('href="assessed/index.html"', (root/'out/index.html').read_text())

    def test_restart_regenerates_failed_report_without_rescanning(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); args = self.fixture(root)
            with patch.object(watch, 'ROOT', root), patch.object(sys, 'argv', args):
                with patch.object(watch, 'execute', side_effect=self.runner(root, report_failure=True)):
                    self.assertEqual(watch.main(), 1)
                with patch.object(watch, 'execute', side_effect=self.runner(root)) as execute:
                    self.assertEqual(watch.main(), 0)
                self.assertEqual(execute.call_count, 1)
                self.assertTrue(execute.call_args.args[0][1].endswith('report_assessment.py'))

    def test_changed_source_invalidates_completed_scan(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); args = self.fixture(root)
            with patch.object(watch, 'ROOT', root), patch.object(sys, 'argv', args):
                with patch.object(watch, 'execute', side_effect=self.runner(root)):
                    self.assertEqual(watch.main(), 0)
                (root/'downloads/SRR1/r1.gz').write_bytes(b'TTTT')
                with patch.object(watch, 'execute') as execute:
                    self.assertEqual(watch.main(), 1)
                execute.assert_not_called()
            status = json.loads((root/'out/status.json').read_text())
            self.assertEqual(status['samples'][0]['state'], 'scan_invalid')
            self.assertFalse(status['report_ready'])

    def test_changed_result_invalidates_completed_scan(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); args = self.fixture(root)
            with patch.object(watch, 'ROOT', root), patch.object(sys, 'argv', args):
                with patch.object(watch, 'execute', side_effect=self.runner(root)):
                    self.assertEqual(watch.main(), 0)
                (root/'out/scans/SRR1/result.json').write_text('{"summary": {}}')
                with patch.object(watch, 'execute') as execute:
                    self.assertEqual(watch.main(), 1)
                execute.assert_not_called()
            self.assertEqual(json.loads((root/'out/status.json').read_text())['samples'][0]['state'], 'scan_invalid')

    def test_stop_after_active_scan_prevents_next_scan_and_report(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); args = self.fixture(root, ('SRR1', 'SRR2'))
            with patch.object(watch, 'ROOT', root), patch.object(sys, 'argv', args), patch.object(watch, 'execute', side_effect=self.runner(root, stop_after_scan=True)) as execute:
                self.assertEqual(watch.main(), 0)
            self.assertEqual(execute.call_count, 1)
            status = json.loads((root/'out/status.json').read_text())
            self.assertEqual(status['watcher_state'], 'stop_requested')
            self.assertFalse(status['complete'])
            self.assertTrue((root/'out/completed_scans.json').exists())

    def test_stop_before_work_launches_nothing(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); args = self.fixture(root)
            (root/'out').mkdir(); (root/'out/STOP_REQUESTED').touch()
            with patch.object(watch, 'ROOT', root), patch.object(sys, 'argv', args), patch.object(watch, 'execute') as execute:
                self.assertEqual(watch.main(), 0)
            execute.assert_not_called()
            self.assertFalse(json.loads((root/'out/status.json').read_text())['complete'])

    def test_zero_completed_legacy_freeze_can_add_code_provenance(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); args = self.fixture(root)
            (root/'out').mkdir(); (root/'out/STOP_REQUESTED').touch()
            legacy = {str(root/'panel.tsv'): watch.sha(root/'panel.tsv')}
            (root/'out/frozen_inputs.json').write_text(json.dumps(legacy))
            with patch.object(watch, 'ROOT', root), patch.object(sys, 'argv', args):
                self.assertEqual(watch.main(), 0)
            self.assertTrue((root/'out/frozen_code.json').exists())


if __name__ == '__main__':
    unittest.main()
