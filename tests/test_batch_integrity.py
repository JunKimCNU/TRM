"""Discovery snapshots, completion state and single-writer publication contracts."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest


ROOT = Path(__file__).resolve().parents[1]
BATCH = ROOT / 'scripts/batch.py'
sys.path.insert(0, str(ROOT / 'scripts'))
import batch
from batch_integrity import output_lock


class BatchIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='trm-batch-integrity-')
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.output = self.root / 'output'
        self.binary = self.root / 'scanner'
        self.manifest = self.root / 'samples.tsv'
        self.make_samples('first')
        self.write_scanner('A')

    def make_samples(self, *names):
        self.manifest.write_text('sample_id\tpath1\n' + ''.join(f'{name}\t{name}.fq\n' for name in names))
        for name in names:
            (self.root / f'{name}.fq').write_text('@x\nTTAGGGTTAGGGTTAGGGTTAGGG\n+\nIIIIIIIIIIIIIIIIIIIIIIII\n')

    def write_scanner(self, label, *, gate=False, exit_code=0):
        temporary = self.binary.with_suffix('.new')
        temporary.write_text(
            '#!/usr/bin/env python3\n'
            'import json, os, pathlib, sys, time\n'
            'args = sys.argv[1:]; source = pathlib.Path(args[1])\n'
            f'label = {label!r}\n'
            f'if {gate!r} and source.stem == "first":\n'
            '    (source.parent / "ready").write_text(str(os.getpid()))\n'
            '    deadline = time.monotonic() + 15\n'
            '    while not (source.parent / "release").exists():\n'
            '        if time.monotonic() >= deadline: raise SystemExit(12)\n'
            '        time.sleep(0.01)\n'
            f'if {exit_code}: raise SystemExit({exit_code})\n'
            'prefix = args[args.index("--output")+1]\n'
            'pathlib.Path(prefix+".json").write_text(json.dumps({"summary":{"records":1,"bases":24},"candidates":[],"actual_scanner":label}))\n')
        temporary.chmod(0o755)
        os.replace(temporary, self.binary)

    def command(self, *options):
        return [sys.executable, '-B', str(BATCH), str(self.manifest), '--binary', str(self.binary),
                '--output-dir', str(self.output), '--jobs', '1', '--threads', '1', '--cpu-budget', '1', *options]

    def run_batch(self, *options, expected=0):
        completed = subprocess.run(self.command(*options), capture_output=True, text=True, timeout=20)
        self.assertEqual(completed.returncode, expected, completed.stdout + completed.stderr)
        return completed

    def start_gated(self, *options):
        self.write_scanner('A', gate=True)
        process = subprocess.Popen(self.command(*options), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        def cleanup():
            (self.root / 'release').touch()
            if process.poll() is None:
                process.terminate()
            process.communicate(timeout=20)
        self.addCleanup(cleanup)
        deadline = time.monotonic() + 10
        while not (self.root / 'ready').exists():
            if process.poll() is not None:
                self.fail('Gated scanner exited early: ' + repr(process.communicate()))
            if time.monotonic() > deadline:
                self.fail('Scanner did not reach its gate')
            time.sleep(0.01)
        return process

    def release(self, process, expected=0):
        (self.root / 'release').touch()
        out, err = process.communicate(timeout=20)
        self.assertEqual(process.returncode, expected, out + err)

    def state(self, sample='first'):
        return json.loads((self.output / sample / 'state.json').read_text())

    def result(self, sample='first'):
        return json.loads((self.output / sample / 'result.json').read_text())

    def test_second_cohort_writer_fails_without_modifying_active_outputs(self):
        process = self.start_gated()
        paths = [self.output / 'first/state.json', self.output / 'summary.json', self.output / 'progress.json']
        before = {path: path.read_bytes() for path in paths}
        contender = self.run_batch('--resume', expected=2)
        self.assertIn('locked by another writer', contender.stderr)
        self.assertEqual(before, {path: path.read_bytes() for path in paths})
        self.release(process)
        self.assertEqual(self.state()['status'], 'complete')

    def test_native_child_keeps_namespace_locked_after_controller_is_killed(self):
        for resource_flag in ((), ('--record-resources',)):
            with self.subTest(resource_flag=resource_flag):
                self.output = self.root / ('timed' if resource_flag else 'untimed')
                (self.root / 'ready').unlink(missing_ok=True)
                (self.root / 'release').unlink(missing_ok=True)
                process = self.start_gated(*resource_flag)
                process.kill()
                process.wait(timeout=5)
                contender = self.run_batch(expected=2)
                self.assertIn('locked by another writer', contender.stderr)
                (self.root / 'release').touch()
                deadline = time.monotonic() + 10
                while True:
                    try:
                        with output_lock(self.output):
                            break
                    except ValueError:
                        if time.monotonic() >= deadline:
                            self.fail('Native child did not release inherited lock')
                        time.sleep(0.02)
                process.communicate(timeout=5)
                self.run_batch()
                self.assertEqual(self.state()['status'], 'complete')

    def test_source_binary_replacement_does_not_change_running_cohort(self):
        self.make_samples('first', 'second')
        process = self.start_gated()
        original_sha = hashlib.sha256(self.binary.read_bytes()).hexdigest()
        self.write_scanner('B')
        self.release(process)
        for name in ('first', 'second'):
            state = self.state(name)
            self.assertEqual(self.result(name)['actual_scanner'], 'A')
            self.assertEqual(state['fingerprint']['executable']['sha256'], original_sha)
            self.assertEqual(state['execution_binary']['sha256'], original_sha)
            snapshot = Path(state['command'][0])
            self.assertEqual(hashlib.sha256(snapshot.read_bytes()).hexdigest(), original_sha)
            self.assertNotEqual(snapshot, self.binary)
            self.assertEqual(snapshot.stat().st_mode & 0o222, 0)
            self.assertNotEqual(snapshot.stat().st_ino, self.binary.stat().st_ino)
        self.run_batch('--resume')
        self.assertEqual(self.result('second')['actual_scanner'], 'B')
        self.assertEqual(len(list((self.output / '.batch-binaries').glob('*/trm'))), 2)

    def test_modified_snapshot_cannot_be_published_or_report_cohort_success(self):
        process = self.start_gated()
        snapshot = Path(self.state()['execution_binary']['path'])
        snapshot.chmod(0o755)
        snapshot.write_text('#!/bin/sh\nexit 0\n')
        self.release(process, expected=1)
        self.assertEqual(self.state()['status'], 'failed')
        self.assertIn('Frozen scanner changed', self.state()['error'])
        self.assertFalse((self.output / 'first/result.json').exists())
        summary = json.loads((self.output / 'summary.json').read_text())
        self.assertFalse(summary['all_succeeded'])
        self.assertIn('Frozen scanner changed', summary['execution_error'])

    def test_failed_rerun_preserves_old_result_but_never_claims_current_success(self):
        self.run_batch()
        result = (self.output / 'first/result.json').read_bytes()
        self.write_scanner('failed', exit_code=7)
        self.run_batch(expected=1)
        self.assertEqual((self.output / 'first/result.json').read_bytes(), result)
        state = self.state()
        self.assertEqual(state['status'], 'failed')
        self.assertTrue(state['retained_previous_result_is_not_current'])
        self.assertNotIn('result_sha256', state)
        self.assertEqual(state['metrics']['exit_code'], 7)
        self.assertEqual(state['schema_version'], 2)

    def test_direct_sample_contender_does_not_change_active_state(self):
        from types import SimpleNamespace
        process = self.start_gated()
        path = self.output / 'first/state.json'
        before = path.read_bytes()
        result = batch.scan_sample({'sample_id':'first', 'paths':[self.root / 'first.fq']},
            SimpleNamespace(output_dir=self.output), {'path':str(self.binary)}, [])
        self.assertEqual(result['status'], 'failed')
        self.assertIn('locked by another writer', result['error'])
        self.assertEqual(path.read_bytes(), before)
        self.release(process)

    def test_symlink_lock_cannot_truncate_another_file(self):
        self.output.mkdir()
        victim = self.root / 'untouched'
        victim.write_text('do not change\n')
        (self.output / '.batch.lock').symlink_to(victim)
        result = self.run_batch(expected=2)
        self.assertIn('must not be a symlink', result.stderr)
        self.assertEqual(victim.read_text(), 'do not change\n')
        self.assertFalse((self.output / 'summary.json').exists())

    def test_corrupt_existing_snapshot_is_refused_without_rewriting_prior_result(self):
        self.run_batch()
        state_bytes = (self.output / 'first/state.json').read_bytes()
        snapshot = Path(self.state()['execution_binary']['path'])
        snapshot.chmod(0o755)
        snapshot.write_text('corrupt')
        result = self.run_batch('--resume', expected=2)
        self.assertIn('Existing frozen scanner is corrupt', result.stderr)
        self.assertEqual((self.output / 'first/state.json').read_bytes(), state_bytes)

    def test_existing_writable_snapshot_is_not_treated_as_immutable(self):
        self.run_batch()
        snapshot = Path(self.state()['execution_binary']['path'])
        snapshot.chmod(0o755)
        result = self.run_batch('--resume', expected=2)
        self.assertIn('not read-only', result.stderr)


if __name__ == '__main__':
    unittest.main()
