"""Input corruption and interrupted downloads must never become valid scans."""
import gzip
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch


def module(name):
    path = Path(__file__).resolve().parents[1] / 'scripts' / (name + '.py')
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


reads = module('inspect_library_reads')
download = module('fetch_test_panel')


class PrefixTests(unittest.TestCase):
    def test_truncated_quality_and_invalid_separator_rejected(self):
        for data in (b'@r\nACGT\n+\n!!', b'@r\nACGT\n-\n!!!!\n'):
            with self.subTest(data=data), self.assertRaises(ValueError):
                reads.record(io.BytesIO(data))

    def test_extracts_complete_records_without_changing_source(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            source = root / 'input.fastq.gz'
            data = b'@a\nACGT\n+\nIIII\n@b\nTTTT\n+\nIIII\n'
            source.write_bytes(gzip.compress(data))
            original = source.read_bytes()
            result = reads.inspect_file(source, root, 1)
            self.assertEqual(source.read_bytes(), original)
            self.assertEqual(result['counts']['records'], 1)
            self.assertEqual(gzip.open(result['extracted_path'], 'rb').read(), data[:15])
            self.assertEqual(result['sampling'], 'first_records_per_file_not_random')


class DownloadTests(unittest.TestCase):
    def setup_file(self, root, data=b'complete archive bytes'):
        item = dict(sample_id='SRR1', name='SRR1.fastq.gz', url='https://ftp.sra.ebi.ac.uk/file',
                    expected_bytes=len(data), expected_md5=hashlib.md5(data).hexdigest())
        return item, root / 'SRR1' / 'SRR1.fastq.gz'

    def test_supervisor_retries_transport_failure_and_stops_on_success(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            args = SimpleNamespace(output=root, rounds=3, wait_for_finished=False,
                manifest=root/'manifest.tsv', jobs=1, samples=['SRR1'])
            calls = []
            def attempt(argv):
                calls.append(argv)
                success = len(calls)==2
                download.save(root/'download_summary.json', {'results':[
                    {'status':'downloaded_verified' if success else 'download_failed_partial_preserved'}]})
                return 0 if success else 1
            with patch.object(download, 'main', attempt), patch.object(download.time, 'sleep'):
                self.assertEqual(download.supervise(args), 0)
            self.assertEqual(len(calls), 2)
            self.assertEqual(len(list((root/'rounds').glob('*/*.json'))), 2)

    def test_supervisor_does_not_retry_integrity_failure(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            args = SimpleNamespace(output=root, rounds=3, wait_for_finished=False,
                manifest=root/'manifest.tsv', jobs=1, samples=['SRR1'])
            def attempt(argv):
                download.save(root/'download_summary.json', {'results':[
                    {'status':'archive_mismatch_partial_preserved'}]})
                return 1
            with patch.object(download, 'main', side_effect=attempt) as child, patch.object(download.time, 'sleep'):
                self.assertEqual(download.supervise(args), 1)
                self.assertEqual(child.call_count, 1)

    def test_existing_archive_mismatch_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            item, target = self.setup_file(root)
            target.parent.mkdir()
            target.write_bytes(b'corrupt')
            result = download.fetch(item, root, 'test')
            self.assertEqual(result['status'], 'existing_archive_mismatch_preserved')
            self.assertEqual(target.read_bytes(), b'corrupt')

    def test_wrong_content_with_right_size_stays_partial(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            item, target = self.setup_file(root, b'correct')
            def transfer(cmd, **kwargs):
                Path(cmd[cmd.index('--output') + 1]).write_bytes(b'incorrect'[:7])
                return subprocess.CompletedProcess(cmd, 0, '{}', '')
            with patch.object(download.subprocess, 'run', transfer):
                result = download.fetch(item, root, 'test')
            self.assertEqual(result['status'], 'archive_mismatch_partial_preserved')
            self.assertFalse(target.exists())
            self.assertTrue(target.with_name(target.name + '.part').exists())

    def test_interrupted_transfer_resumes_and_verifies(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            data = b'complete archive bytes'
            item, target = self.setup_file(root, data)
            calls = []
            def transfer(cmd, **kwargs):
                part = Path(cmd[cmd.index('--output') + 1])
                calls.append(cmd)
                if len(calls) == 1:
                    part.write_bytes(data[:5])
                    return subprocess.CompletedProcess(cmd, 18, '{}', 'EOF')
                self.assertEqual(part.read_bytes(), data[:5])
                self.assertIn('--continue-at', cmd)
                with part.open('ab') as h:
                    h.write(data[5:])
                return subprocess.CompletedProcess(cmd, 0, '{}', '')
            with patch.object(download.subprocess, 'run', transfer), patch.object(download.time, 'sleep'):
                result = download.fetch(item, root, 'test')
            self.assertEqual(result['status'], 'downloaded_verified')
            self.assertEqual(target.read_bytes(), data)
            self.assertEqual(result['transferred_bytes'], len(data))
            self.assertEqual(len(result['attempts']), 2)
            self.assertEqual(result['attempts_this_invocation'], 2)
            active = json.loads(target.with_name(target.name + '.active.json').read_text())
            self.assertEqual(active['status'], 'finished')

    def test_complete_partial_verified_without_network(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            data = b'complete archive bytes'
            item, target = self.setup_file(root, data)
            target.parent.mkdir()
            part = target.with_name(target.name + '.part')
            part.write_bytes(data)
            with patch.object(download.subprocess, 'run') as run:
                result = download.fetch(item, root, 'test')
            run.assert_not_called()
            self.assertEqual(result['status'], 'verified_existing_partial')
            self.assertEqual(result['transferred_bytes'], 0)
            self.assertEqual(target.read_bytes(), data)
            self.assertFalse(part.exists())

    def test_oversized_partial_preserved_without_network(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            item, target = self.setup_file(root, b'correct')
            target.parent.mkdir()
            part = target.with_name(target.name + '.part')
            part.write_bytes(b'oversized')
            with patch.object(download.subprocess, 'run') as run:
                result = download.fetch(item, root, 'test')
            run.assert_not_called()
            self.assertEqual(result['status'], 'archive_mismatch_partial_preserved')
            self.assertEqual(part.read_bytes(), b'oversized')
            self.assertFalse(target.exists())

    def test_retry_accounting_separates_history_and_retained_bytes(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            data = b'complete archive bytes'
            item, target = self.setup_file(root, data)
            target.parent.mkdir()
            part = target.with_name(target.name + '.part')
            part.write_bytes(data[:3])
            target.with_name(target.name + '.attempts.json').write_text('[{"returncode":18,"retained_new_bytes":3}]')
            def transfer(cmd, **kwargs):
                with part.open('ab') as handle:
                    handle.write(data[3:])
                return subprocess.CompletedProcess(cmd, 0, '{"size_download":999}', '')
            with patch.object(download.subprocess, 'run', transfer):
                result = download.fetch(item, root, 'test')
            self.assertEqual(result['prior_attempt_count'], 1)
            self.assertEqual(result['attempts_this_invocation'], 1)
            self.assertEqual(len(result['attempts']), 2)
            self.assertEqual(result['transferred_bytes'], len(data) - 3)
            self.assertEqual(result['retained_new_bytes_this_invocation'], len(data) - 3)
            self.assertIn('not total network traffic', result['byte_accounting'])

    def test_retries_bounded_and_partial_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            item, target = self.setup_file(root)
            def transfer(cmd, **kwargs):
                part = Path(cmd[cmd.index('--output') + 1])
                with part.open('ab') as handle:
                    handle.write(b'x')
                return subprocess.CompletedProcess(cmd, 18, '{}', 'EOF')
            with patch.object(download.subprocess, 'run', side_effect=transfer) as run, patch.object(download.time, 'sleep'):
                result = download.fetch(item, root, 'test')
            self.assertEqual(run.call_count, 4)
            self.assertEqual(result['status'], 'download_failed_partial_preserved')
            self.assertEqual(target.with_name(target.name + '.part').read_bytes(), b'xxxx')
            self.assertEqual(result['transferred_bytes'], 4)

    def test_curl_success_without_output_is_failure(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            item, target = self.setup_file(root)
            with patch.object(download.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, '{}', '')):
                result = download.fetch(item, root, 'test')
            self.assertEqual(result['status'], 'download_failed_partial_preserved')
            self.assertFalse(target.exists())

    def test_complete_checksum_can_verify_after_curl_error(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            data = b'complete archive bytes'
            item, target = self.setup_file(root, data)
            def transfer(cmd, **kwargs):
                Path(cmd[cmd.index('--output') + 1]).write_bytes(data)
                return subprocess.CompletedProcess(cmd, 18, '{}', 'EOF after archive bytes')
            with patch.object(download.subprocess, 'run', transfer):
                result = download.fetch(item, root, 'test')
            self.assertEqual(result['status'], 'downloaded_verified')
            self.assertEqual(result['returncode'], 18)
            self.assertEqual(target.read_bytes(), data)

    def test_prior_identity_change_does_not_resume(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            item, target = self.setup_file(root)
            target.parent.mkdir()
            part = target.with_name(target.name + '.part')
            part.write_bytes(b'part')
            identity = target.with_name(target.name + '.identity.json')
            identity.write_text(json.dumps(dict(item, expected_md5='0' * 32)))
            original = identity.read_bytes()
            with patch.object(download.subprocess, 'run') as run:
                result = download.fetch(item, root, 'test')
            run.assert_not_called()
            self.assertEqual(result['status'], 'download_exception_partial_preserved')
            self.assertEqual(part.read_bytes(), b'part')
            self.assertEqual(identity.read_bytes(), original)

    def test_invalid_history_preserved_and_failure_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            item, target = self.setup_file(root)
            target.parent.mkdir()
            history = target.with_name(target.name + '.attempts.json')
            history.write_text('not json')
            with patch.object(download.subprocess, 'run') as run:
                result = download.fetch(item, root, 'test')
            run.assert_not_called()
            self.assertEqual(result['status'], 'download_exception_partial_preserved')
            self.assertEqual(history.read_text(), 'not json')

    def test_lock_conflict_does_not_overwrite_owner_state(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            item, target = self.setup_file(root)
            target.parent.mkdir()
            active = target.with_name(target.name + '.active.json')
            active.write_text('{"pid":123,"status":"running"}')
            original = active.read_bytes()
            with patch.object(download.fcntl, 'flock', side_effect=BlockingIOError), patch.object(download.subprocess, 'run') as run:
                result = download.fetch(item, root, 'test')
            run.assert_not_called()
            self.assertEqual(result['status'], 'download_in_progress')
            self.assertEqual(active.read_bytes(), original)

    def test_manifest_rejects_bad_integrity_and_paths(self):
        row = dict(sample_id='SRR1', download_urls='https://ftp.sra.ebi.ac.uk/vol/shared.fastq.gz',
                   fastq_bytes='100', fastq_md5='a' * 32)
        for key, value in (('sample_id', '../escape'), ('fastq_bytes', '0'), ('fastq_bytes', '-1'),
                           ('fastq_md5', 'z' * 32), ('download_urls', 'https://ftp.sra.ebi.ac.uk/a.fastq.gz?query=1')):
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                download.planned_items([dict(row, **{key: value})])
        self.assertEqual(download.planned_items([dict(row, fastq_md5='A' * 32)])[0]['expected_md5'], 'a' * 32)

    def test_same_basename_across_samples_keeps_distinct_paths(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            manifest = root / 'input.tsv'
            rows = [dict(sample_id=sample, download_urls=f'https://ftp.sra.ebi.ac.uk/{sample}/shared.fastq.gz',
                         fastq_bytes='100', fastq_md5='a' * 32) for sample in ('SRR1', 'SRR2')]
            with manifest.open('w') as handle:
                writer = csv.DictWriter(handle, delimiter='\t', fieldnames=list(rows[0]))
                writer.writeheader(); writer.writerows(rows)
            def fake_fetch(item, out, phase):
                return dict(item, status='verified_existing', path=str(out / item['sample_id'] / item['name']), transferred_bytes=0)
            output = root / 'out'
            with patch.object(download, 'fetch', side_effect=fake_fetch), patch.object(sys, 'argv', ['fetch', '--manifest', str(manifest), '--output', str(output)]):
                status = download.main()
            self.assertEqual(status, 0)
            with (output / 'verified_fastq.tsv').open() as handle:
                actual = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(len(actual), 2)
            for row in actual:
                self.assertIn('/' + row['sample_id'] + '/', row['path1'])
            summary = json.loads((output / 'download_summary.json').read_text())
            self.assertEqual(summary['verified_files'], 2)
            self.assertTrue(summary['complete'])

    def test_single_first_exception_still_writes_summary(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); manifest = root / 'input.tsv'; output = root / 'out'
            manifest.write_text('sample_id\tdownload_urls\tfastq_bytes\tfastq_md5\n'
                                'SRR1\thttps://ftp.sra.ebi.ac.uk/SRR1.fastq.gz\t100\t' + 'a' * 32 + '\n')
            with patch.object(download, 'fetch', side_effect=OSError('denied')), patch.object(sys, 'argv', ['fetch', '--manifest', str(manifest), '--output', str(output), '--single-first']):
                status = download.main()
            self.assertEqual(status, 1)
            summary = json.loads((output / 'download_summary.json').read_text())
            progress = json.loads((output / 'progress.json').read_text())
            self.assertFalse(summary['complete'])
            self.assertFalse(summary['retained_file_growth_accounting_complete'])
            self.assertEqual(progress['status'], 'finished')
            self.assertEqual(progress['files_finished'], 1)


if __name__ == '__main__':
    unittest.main()
