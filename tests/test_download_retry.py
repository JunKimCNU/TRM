"""Regression tests for incomplete transfers and full-size integrity failures."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

download = load('retry_fetch', ROOT/'scripts/fetch_test_panel.py')

class RecoveryTests(unittest.TestCase):
    def setUp(self):
        sleeper = patch.object(download.time, "sleep")
        sleeper.start()
        self.addCleanup(sleeper.stop)

    def fixture(self, root):
        data = b'correct archive bytes'
        item = dict(sample_id='SRR1', name='SRR1.fastq.gz',
                    url='https://ftp.sra.ebi.ac.uk/file', expected_bytes=len(data),
                    expected_md5=hashlib.md5(data).hexdigest())
        target = root/'SRR1'/'SRR1.fastq.gz'
        return item, target, target.with_name(target.name+'.part'), data

    def test_exit_zero_empty_partial_is_bounded_transport_failure_without_hashing(self):
        with tempfile.TemporaryDirectory() as directory:
            item, target, part, data = self.fixture(Path(directory))
            def transfer(command, **kwargs):
                part.touch()
                return subprocess.CompletedProcess(command, 0, '{"http_code":000}', '')
            with patch.object(download.subprocess, 'run', side_effect=transfer) as run, \
                 patch.object(download, 'digests') as digests:
                result = download.fetch(item, Path(directory), 'test')
            self.assertEqual(run.call_count, 4)
            digests.assert_not_called()
            self.assertEqual(result['status'], 'download_failed_partial_preserved')
            self.assertEqual(result['observed_bytes'], 0)
            self.assertEqual(result['missing_bytes'], len(data))
            self.assertEqual(part.read_bytes(), b'')
            self.assertFalse(target.exists())

    def test_exit_zero_truncated_transfer_resumes_without_discarding_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            item, target, part, data = self.fixture(Path(directory))
            calls = []
            def transfer(command, **kwargs):
                calls.append(command)
                if len(calls) == 1:
                    part.write_bytes(data[:7])
                else:
                    self.assertEqual(part.read_bytes(), data[:7])
                    self.assertEqual(command[command.index('--continue-at')+1], '-')
                    with part.open('ab') as handle:
                        handle.write(data[7:])
                return subprocess.CompletedProcess(command, 0, '{}', '')
            with patch.object(download.subprocess, 'run', side_effect=transfer):
                result = download.fetch(item, Path(directory), 'test')
            self.assertEqual(result['status'], 'downloaded_verified')
            self.assertEqual(result['attempts_this_invocation'], 2)
            self.assertEqual(target.read_bytes(), data)
            self.assertFalse(part.exists())

    def test_old_incomplete_mismatch_receipt_does_not_prevent_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            item, target, part, data = self.fixture(Path(directory))
            part.parent.mkdir()
            part.write_bytes(data[:7])
            old = dict(item, status='archive_mismatch_partial_preserved', observed_bytes=7)
            target.with_name(target.name+'.download.json').write_text(json.dumps(old))
            def transfer(command, **kwargs):
                self.assertEqual(part.read_bytes(), data[:7])
                with part.open('ab') as handle:
                    handle.write(data[7:])
                return subprocess.CompletedProcess(command, 0, '{}', '')
            with patch.object(download.subprocess, 'run', side_effect=transfer):
                result = download.fetch(item, Path(directory), 'test')
            self.assertEqual(result['status'], 'downloaded_verified')
            self.assertEqual(result['resumed_from_bytes'], 7)
            self.assertEqual(result['transferred_bytes'], len(data)-7)
            self.assertEqual(target.read_bytes(), data)

    def test_completed_wrong_hash_does_not_retry_or_rename(self):
        with tempfile.TemporaryDirectory() as directory:
            item, target, part, data = self.fixture(Path(directory))
            wrong = b'x'*len(data)
            def transfer(command, **kwargs):
                part.write_bytes(wrong)
                return subprocess.CompletedProcess(command, 0, '{}', '')
            with patch.object(download.subprocess, 'run', side_effect=transfer) as run:
                result = download.fetch(item, Path(directory), 'test')
            self.assertEqual(run.call_count, 1)
            self.assertEqual(result['status'], 'archive_mismatch_partial_preserved')
            self.assertEqual(part.read_bytes(), wrong)
            self.assertFalse(target.exists())

    def test_oversized_transfer_with_curl_error_is_integrity_stop(self):
        with tempfile.TemporaryDirectory() as directory:
            item, target, part, data = self.fixture(Path(directory))
            wrong = data+b'oversized'
            def transfer(command, **kwargs):
                part.write_bytes(wrong)
                return subprocess.CompletedProcess(command, 22, '{}', 'HTTP failure')
            with patch.object(download.subprocess, 'run', side_effect=transfer) as run:
                result = download.fetch(item, Path(directory), 'test')
            self.assertEqual(run.call_count, 1)
            self.assertEqual(result['status'], 'archive_mismatch_partial_preserved')
            self.assertEqual(part.read_bytes(), wrong)
            self.assertFalse(target.exists())


if __name__ == "__main__":
    unittest.main()
