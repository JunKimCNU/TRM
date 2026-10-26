"""Optional local MD5 receipts must bind actual files without a second hash pass."""
from __future__ import annotations

from contextlib import contextmanager
import copy
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import batch_integrity as integrity
import verify_external_download as verifier


class BatchInputVerificationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='trm-batch-receipt-')
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.incoming = self.root / 'incoming'
        self.incoming.mkdir()
        self.raw = self.incoming / 'first.fastq.gz'
        self.raw.write_bytes(b'@x\nTTAGGGTTAGGG\n+\nIIIIIIIIIIII\n')
        self.manifest = self.root / 'batch.tsv'
        self.manifest.write_text(f'sample_id\tpath1\nfirst\t{self.raw}\n')
        self.output = self.root / 'output'
        self.binary = self.root / 'scanner'
        self.binary.write_text(
            '#!/usr/bin/env python3\nimport json,pathlib,sys,time\n'
            f'root=pathlib.Path({str(self.root)!r})\n'
            'a=sys.argv\n'
            'if (root/"gate").exists():\n'
            ' (root/"ready").touch()\n'
            ' while not (root/"release").exists(): time.sleep(0.01)\n'
            'pathlib.Path(a[a.index("--output")+1]+".json").write_text(json.dumps({"summary":{"records":1,"bases":12},"candidates":[]}))\n')
        self.binary.chmod(0o755)
        self.report = self.make_report()

    def row(self, path, payload=None):
        payload = path.read_bytes() if payload is None else payload
        return {'file_id':path.name, 'run_accession':'SRR1', 'archive_name':path.name,
                'relative_path':str(path.relative_to(self.incoming)), 'url':'https://example.invalid/'+path.name,
                'expected_bytes':str(len(payload)), 'expected_md5':hashlib.md5(payload).hexdigest(),
                'scientific_name':'Synthetic control', 'task_id':'offline_fixture'}

    def make_report(self, rows=None, label='verification'):
        rows = [self.row(self.raw)] if rows is None else rows
        archive_manifest = self.root / (label + '.tsv')
        with archive_manifest.open('w', newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=verifier.REQUIRED,delimiter='\t')
            writer.writeheader();writer.writerows(rows)
        directory = self.root / label
        verifier.run(archive_manifest,self.incoming,directory,jobs=1)
        return directory / 'verification.json'

    def command(self, *options, report=True):
        args = [sys.executable,'-B',str(ROOT/'scripts/batch.py'),str(self.manifest),
                '--binary',str(self.binary),'--output-dir',str(self.output),'--jobs','1','--threads','1','--cpu-budget','1']
        if report:
            args += ['--input-verification-report',str(self.report)]
        return [*args,*options]

    def run_batch(self, *options, expected=0, report=True):
        result=subprocess.run(self.command(*options,report=report),capture_output=True,text=True,timeout=15)
        self.assertEqual(result.returncode,expected,result.stdout+result.stderr)
        return result

    def state(self):
        return json.loads((self.output/'first/state.json').read_text())

    def test_opt_in_receipt_binds_md5_and_resume_without_changing_result_bytes(self):
        self.run_batch(report=False)
        before=(self.output/'first/result.json').read_bytes()
        self.assertNotIn('input_verification',self.state()['fingerprint'])
        self.run_batch('--resume')
        self.assertEqual((self.output/'first/result.json').read_bytes(),before)
        binding=self.state()['fingerprint']['input_verification']
        self.assertEqual(binding['receipt']['sha256'],hashlib.sha256(self.report.read_bytes()).hexdigest())
        self.assertEqual(binding['files'][0]['md5'],hashlib.md5(self.raw.read_bytes()).hexdigest())
        state=(self.output/'first/state.json').read_bytes()
        self.run_batch('--resume')
        self.assertEqual((self.output/'first/state.json').read_bytes(),state)
        summary=json.loads((self.output/'summary.json').read_text())
        self.assertEqual(summary['resumed_samples'],1)
        self.assertIn('MD5',summary['input_identity_policy'])

    def test_changed_contents_with_restored_size_and_mtime_are_rejected(self):
        self.run_batch()
        info=self.raw.stat();time.sleep(0.01)
        self.raw.write_bytes(self.raw.read_bytes().replace(b'TTAGGG',b'TTAGAC'))
        os.utime(self.raw,ns=(info.st_atime_ns,info.st_mtime_ns))
        self.assertEqual((self.raw.stat().st_size,self.raw.stat().st_mtime_ns),(info.st_size,info.st_mtime_ns))
        result=self.run_batch('--resume',expected=2)
        self.assertIn('identity is stale or replaced',result.stderr)

    def test_copy_replacing_same_path_requires_reverification(self):
        duplicate=self.incoming/'copy.tmp'
        shutil.copy2(self.raw,duplicate)
        os.replace(duplicate,self.raw)
        result=self.run_batch(expected=2)
        self.assertIn('identity is stale or replaced',result.stderr)
        self.report=self.make_report(label='reverified')
        self.run_batch()

    def test_copy_to_new_path_is_not_covered_by_original_report(self):
        copied=self.incoming/'copied.fastq.gz';shutil.copy2(self.raw,copied)
        self.manifest.write_text(f'sample_id\tpath1\nfirst\t{copied}\n')
        self.assertIn('does not cover',self.run_batch(expected=2).stderr)
        self.report=self.make_report([self.row(copied)],label='copied_verified')
        self.run_batch()

    def test_only_selected_valid_rows_are_required_from_partial_report(self):
        bad=self.incoming/'bad.fastq.gz';bad.write_bytes(b'wrong')
        missing=self.incoming/'missing.fastq.gz'
        self.report=self.make_report([self.row(self.raw),self.row(bad,b'right'),self.row(missing,b'absent')],label='partial')
        document=json.loads(self.report.read_text())
        self.assertEqual(document['exit_code'],1)
        self.assertFalse(document['all_manifest_files_present_and_valid'])
        self.run_batch()

    def test_legacy_incomplete_invalid_md5_size_and_nonverified_records_rejected(self):
        original=json.loads(self.report.read_text())
        variants=[]
        def changed(edit):
            value=copy.deepcopy(original);edit(value);variants.append(value)
        changed(lambda d:d.update(schema_version=1))
        changed(lambda d:d.update(status='running'))
        changed(lambda d:d['files'][0].update(status='mismatch'))
        changed(lambda d:d['files'][0].update(expected_md5='not_an_md5',observed_md5='not_an_md5'))
        changed(lambda d:d['files'][0].update(observed_bytes=999))
        changed(lambda d:d['files'][0]['verified_file_identity'].update(size_bytes=True))
        changed(lambda d:d['files'][0].pop('verified_file_identity'))
        for index,value in enumerate(variants):
            with self.subTest(index=index):
                path=self.root/f'bad_{index}.json';path.write_text(json.dumps(value))
                with self.assertRaises(ValueError):
                    integrity.load_input_verification(path,[self.raw])
        self.assertFalse(self.output.exists())

    def test_consistent_duplicate_paths_allowed_but_conflicting_records_rejected(self):
        value=json.loads(self.report.read_text())
        duplicate=copy.deepcopy(value['files'][0]);duplicate.update(file_id='another_accession',run_accession='SRR2')
        value['files'].append(duplicate)
        path=self.root/'duplicate.json';path.write_text(json.dumps(value))
        context=integrity.load_input_verification(path,[self.raw,self.raw])
        self.assertEqual(len(context['files']),1)
        duplicate.update(expected_md5='0'*32,observed_md5='0'*32)
        path.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError,'Conflicting'):
            integrity.load_input_verification(path,[self.raw])

    def gated_batch(self):
        (self.root/'gate').touch()
        process=subprocess.Popen(self.command(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        def cleanup():
            (self.root/'release').touch()
            if process.poll() is None:process.terminate()
            process.communicate(timeout=15)
        self.addCleanup(cleanup)
        deadline=time.monotonic()+10
        while not (self.root/'ready').exists():
            if process.poll() is not None:self.fail(repr(process.communicate()))
            if time.monotonic()>deadline:self.fail('Scanner did not reach gate')
            time.sleep(0.01)
        return process

    def test_receipt_rewritten_midscan_cannot_complete(self):
        process=self.gated_batch();payload=self.report.read_bytes()
        self.report.write_bytes(payload)
        (self.root/'release').touch()
        out,err=process.communicate(timeout=15)
        self.assertEqual(process.returncode,1,out+err)
        self.assertEqual(self.state()['status'],'failed')
        self.assertIn('verification report changed',self.state()['error'])
        self.assertFalse((self.output/'first/result.json').exists())

    def test_input_changed_midscan_even_with_restored_mtime_cannot_complete(self):
        process=self.gated_batch();info=self.raw.stat()
        self.raw.write_bytes(self.raw.read_bytes().replace(b'TTAGGG',b'TTAGAC'))
        os.utime(self.raw,ns=(info.st_atime_ns,info.st_mtime_ns))
        (self.root/'release').touch()
        out,err=process.communicate(timeout=15)
        self.assertEqual(process.returncode,1,out+err)
        self.assertIn('identity is stale or replaced',self.state()['error'])

    def test_symlink_report_and_input_components_are_refused(self):
        alias=self.root/'alias.json';alias.symlink_to(self.report)
        with self.assertRaises((OSError,ValueError)):
            integrity.load_input_verification(alias,[self.raw])
        context=integrity.load_input_verification(self.report,[self.raw])
        moved=self.root/'moved';self.incoming.rename(moved)
        self.incoming.symlink_to(moved,target_is_directory=True)
        with self.assertRaises((OSError,ValueError)):
            integrity.assert_input_verification(context)

    def test_receipt_checks_never_read_archive_payload(self):
        original=integrity._open_regular_nofollow
        @contextmanager
        def guarded(path):
            with original(path) as handle:
                if Path(path)==self.raw:
                    class StatOnly:
                        def fileno(self):return handle.fileno()
                        def read(self,*args):raise AssertionError('Archive bytes were read again')
                    yield StatOnly()
                else:
                    yield handle
        with patch.object(integrity,'_open_regular_nofollow',guarded):
            context=integrity.load_input_verification(self.report,[self.raw])
            integrity.assert_input_verification(context,hash_report=True)


if __name__=='__main__':
    unittest.main()
