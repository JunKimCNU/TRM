import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/verify_external_download.py'
spec = importlib.util.spec_from_file_location('verify_external_download', SCRIPT)
verify = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verify)


class ExternalVerificationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.incoming = self.root / 'incoming'
        self.incoming.mkdir()
        self.manifest = self.root / 'manifest.tsv'

    def row(self, name='SRR123_1.fastq.gz', payload=b'synthetic gzipped bytes', **changes):
        row = dict(file_id=name, run_accession='SRR123', archive_name=name,
                   relative_path='SRR123/' + name, url='https://example.invalid/' + name,
                   expected_bytes=str(len(payload)), expected_md5=hashlib.md5(payload).hexdigest(),
                   scientific_name='Synthetic control', task_id='task_SRR123', role='R1', representation='ena_fastq')
        row.update(changes)
        return row

    def write_manifest(self, rows):
        with self.manifest.open('w', encoding='utf-8-sig', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=verify.REQUIRED + verify.OPTIONAL, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)

    def put(self, row, data=b'synthetic gzipped bytes'):
        path = self.incoming / row['relative_path']
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        return path

    def run_verify(self, rows, **kwargs):
        self.write_manifest(rows)
        return verify.run(self.manifest, self.incoming, self.root / 'reports', **kwargs)

    def test_valid_file_and_identical_duplicate_emit_one_verified_row(self):
        row = self.row()
        source = self.put(row)
        before = source.stat()
        result = self.run_verify([row, row])
        self.assertEqual(result['counts'], dict(valid=1, missing=0, mismatch=0, error=0))
        self.assertEqual(result['exit_code'], 0)
        self.assertEqual(source.stat().st_mtime_ns, before.st_mtime_ns)
        with (self.root / 'reports/verified_files.tsv').open() as handle:
            verified = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(verified), 1)
        self.assertEqual(verified[0]['absolute_path'], str(source))
        self.assertFalse(result['project_ledgers_modified'])
        self.assertEqual(result['schema_version'], 2)
        self.assertEqual(result['program'], 'trm-external-download-verification')
        self.assertEqual(result['status'], 'complete')
        self.assertEqual(result['files'][0]['verified_file_identity'], verify.verified_file_identity(before))
        saved = json.loads((self.root / 'reports/verification.json').read_text())
        self.assertEqual(saved['files'][0]['verified_file_identity'], verify.verified_file_identity(before))

    def test_same_size_corruption_fails_md5_and_is_not_verified(self):
        row = self.row(payload=b'abcd')
        source = self.put(row, b'abce')
        result = self.run_verify([row])
        self.assertEqual(result['exit_code'], 1)
        self.assertEqual(result['files'][0]['reason'], 'md5_mismatch')
        self.assertNotIn('verified_file_identity', result['files'][0])
        self.assertEqual(source.read_bytes(), b'abce')
        self.assertEqual(len((self.root / 'reports/verified_files.tsv').read_text().splitlines()), 1)

    def test_size_mismatch_skips_hash(self):
        row = self.row(payload=b'abcd')
        self.put(row, b'a')
        result = self.run_verify([row])
        self.assertEqual(result['counts']['mismatch'], 1)
        self.assertIsNone(result['files'][0]['observed_md5'])

    def test_missing_allowed_for_partial_import(self):
        result = self.run_verify([self.row()])
        self.assertEqual(result['counts']['missing'], 1)
        self.assertEqual(result['exit_code'], 0)
        self.assertFalse(result['all_manifest_files_present_and_valid'])

    def test_require_all_makes_missing_fail(self):
        self.assertEqual(self.run_verify([self.row()], require_all=True)['exit_code'], 1)

    def test_part_sidecar_is_never_promoted(self):
        row = self.row()
        partial = self.incoming / (row['relative_path'] + '.part')
        partial.parent.mkdir()
        partial.write_bytes(b'synthetic gzipped bytes')
        result = self.run_verify([row])
        self.assertEqual(result['counts']['missing'], 1)
        self.assertTrue(partial.exists())
        self.assertFalse((self.incoming / row['relative_path']).exists())

    def test_unsafe_manifest_paths_rejected(self):
        for path in ['../x.fastq.gz', '/x.fastq.gz', 'C:/x.fastq.gz', 'C:\\x.fastq.gz',
                     'a/../x.fastq.gz', 'a//x.fastq.gz', './x.fastq.gz',
                     'quarantine/x.fastq.gz', 'failed_generations/x.fastq.gz', 'x.fastq.gz.part']:
            with self.subTest(path=path):
                row = self.row(name=path.split('/')[-1], relative_path=path)
                self.write_manifest([row])
                with self.assertRaises(verify.ValidationError):
                    verify.read_manifest(self.manifest)

    def test_conflicting_file_id_or_path_rejected(self):
        row = self.row()
        for changed in [dict(row, expected_md5='0' * 32), dict(row, file_id='different')]:
            self.write_manifest([row, changed])
            with self.assertRaises(verify.ValidationError):
                verify.read_manifest(self.manifest)

    def test_symlink_directory_escape_cannot_verify_external_file(self):
        outside = self.root / 'outside'
        outside.mkdir()
        row = self.row()
        (outside / row['archive_name']).write_bytes(b'synthetic gzipped bytes')
        (self.incoming / 'SRR123').symlink_to(outside, target_is_directory=True)
        result = self.run_verify([row])
        self.assertEqual(result['counts']['error'], 1)
        self.assertEqual(result['exit_code'], 1)

    def test_symlink_file_cannot_verify_external_file(self):
        row = self.row()
        outside = self.root / 'outside.fastq.gz'
        outside.write_bytes(b'synthetic gzipped bytes')
        path = self.incoming / row['relative_path']
        path.parent.mkdir()
        path.symlink_to(outside)
        result = self.run_verify([row])
        self.assertEqual(result['counts']['error'], 1)

    def test_output_cannot_overwrite_input(self):
        row = self.row()
        source = self.put(row)
        self.write_manifest([row])
        with self.assertRaises(verify.ValidationError):
            verify.run(self.manifest, self.incoming, self.incoming / 'reports')
        self.assertEqual(source.read_bytes(), b'synthetic gzipped bytes')

    def test_existing_report_is_never_overwritten(self):
        row = self.row()
        self.put(row)
        self.run_verify([row])
        path = self.root / 'reports/verification.json'
        original = path.read_bytes()
        with self.assertRaises(verify.ValidationError):
            verify.run(self.manifest, self.incoming, path.parent)
        self.assertEqual(path.read_bytes(), original)

    def test_replaced_file_is_not_verified_under_stale_path(self):
        row = self.row(payload=b'abcd')
        source = self.put(row, b'abcd')
        replacement = source.with_name('replacement')
        replacement.write_bytes(b'bad!')
        original_open = verify.open_regular_file
        opens = []

        def replace_before_second_open(root, relative):
            opens.append(relative)
            if len(opens) == 2:
                replacement.replace(source)
            return original_open(root, relative)

        with patch.object(verify, 'open_regular_file', side_effect=replace_before_second_open):
            result = self.run_verify([row])
        self.assertEqual(result['counts']['error'], 1)
        self.assertIn('pathname changed', result['files'][0]['reason'])
        self.assertNotIn('verified_file_identity', result['files'][0])

    def test_reports_outside_input_cannot_be_symlinks(self):
        row = self.row()
        self.put(row)
        self.write_manifest([row])
        destination = self.root / 'reports'
        destination.symlink_to(self.incoming, target_is_directory=True)
        with self.assertRaises(verify.ValidationError):
            verify.run(self.manifest, self.incoming, destination)

    def test_invalid_manifest_cli_returns_nonzero_before_reports(self):
        self.write_manifest([self.row(expected_bytes='garbage')])
        self.assertEqual(verify.main(['--manifest', str(self.manifest), '--input-root', str(self.incoming),
                                      '--output', str(self.root / 'reports')]), 2)
        self.assertFalse((self.root / 'reports').exists())


if __name__ == '__main__':
    unittest.main()
