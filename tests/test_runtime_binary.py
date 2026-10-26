"""Current-runtime selection must never fall back to preserved historical scans."""
import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from archive_test_support import archive_only
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import runtime_binary as runtime
import install_runtime as installer


class RuntimeBinaryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='trm-runtime-test-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        (self.root/'Cargo.toml').write_text('[package]\nname = "trm"\nversion = "0.8.4"\n')
        (self.root/'Cargo.lock').write_text('# fixture\n')
        (self.root/'src').mkdir()
        self.source = self.make_binary('current', '0.8.4')
        self.old = self.make_binary('target/release/trm', '0.8.2')

    def make_binary(self, name, version):
        path = self.root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f'#!/bin/sh\nprintf "trm {version}\\n"\n')
        path.chmod(0o755)
        return path

    def install(self):
        return installer.install(self.source, root=self.root,
                                 expected_sha256=runtime.sha256(self.source),
                                 provenance={'source_sha256': 'recorded-at-build', 'compiler': None})

    def test_install_verifies_hash_version_and_preserves_historical_binary(self):
        original = self.old.read_bytes()
        result = self.install()
        self.assertEqual(result['version'], '0.8.4')
        self.assertEqual(runtime.resolve_binary(root=self.root, environ={}), self.root/'bin/trm')
        self.assertEqual(self.old.read_bytes(), original)
        receipt = json.loads((self.root/'bin/runtime.json').read_text())
        self.assertEqual(receipt['sha256'], runtime.sha256(self.source))
        # Move away the original results/build source; defaults stay usable.
        self.source.unlink()
        self.assertEqual(runtime.installed_identity(self.root)['version'], '0.8.4')

    def test_missing_runtime_is_actionable_and_does_not_use_existing_target(self):
        with self.assertRaisesRegex(ValueError, 'install_runtime.py'):
            runtime.resolve_binary(root=self.root, environ={})
        self.assertTrue(self.old.is_file())

    def test_explicit_then_environment_then_managed_precedence(self):
        self.install()
        self.assertEqual(runtime.resolve_binary(self.source, root=self.root, environ={'TRM_BINARY': str(self.old)}), self.source)
        self.assertEqual(runtime.resolve_binary(root=self.root, environ={'TRM_BINARY': str(self.old)}), self.old)
        self.assertEqual(runtime.resolve_binary(root=self.root, environ={}), self.root/'bin/trm')
        with self.assertRaisesRegex(ValueError, 'must not be empty'):
            runtime.resolve_binary(root=self.root, environ={'TRM_BINARY': ''})
        with self.assertRaisesRegex(ValueError, 'not an executable'):
            runtime.resolve_binary(root=self.root, environ={'TRM_BINARY': str(self.root/'missing')})

    def test_corrupt_binary_or_receipt_rejected_even_with_explicit_managed_path(self):
        self.install()
        binary = self.root/'bin/trm'
        binary.write_text(binary.read_text()+'# corrupt\n')
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            runtime.resolve_binary(binary, root=self.root, environ={})
        binary.write_bytes(self.source.read_bytes())
        (self.root/'bin/runtime.json').write_text('[]')
        with self.assertRaisesRegex(ValueError, 'Invalid runtime receipt'):
            runtime.resolve_binary(root=self.root, environ={})

    def test_version_and_source_digest_mismatch_rejected_before_install(self):
        with self.assertRaisesRegex(ValueError, 'must be executable and report trm 0.8.4'):
            installer.install(self.old, root=self.root)
        with self.assertRaisesRegex(ValueError, 'differs from'):
            installer.install(self.source, root=self.root, expected_sha256='0'*64)
        self.assertFalse((self.root/'bin').exists())
        self.install()
        (self.root/'Cargo.toml').write_text('[package]\nversion = "0.8.5"\n')
        with self.assertRaisesRegex(ValueError, 'differs from project'):
            runtime.installed_identity(self.root)

    def test_same_install_is_idempotent_and_different_install_does_not_overwrite(self):
        self.install()
        before = (self.root/'bin/runtime.json').read_bytes()
        self.assertEqual(self.install()['status'], 'already_installed')
        self.assertEqual((self.root/'bin/runtime.json').read_bytes(), before)
        self.source.write_text(self.source.read_text()+'# different build\n')
        with self.assertRaisesRegex(ValueError, 'different managed runtime'):
            self.install()
        self.assertEqual((self.root/'bin/runtime.json').read_bytes(), before)

    def test_symlink_install_directory_rejected(self):
        elsewhere = self.root/'elsewhere'; elsewhere.mkdir()
        (self.root/'bin').symlink_to(elsewhere, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'must not be a symlink'):
            self.install()
        self.assertEqual(list(elsewhere.iterdir()), [])

    def test_fresh_source_build_uses_isolated_target_without_results_directory(self):
        actual_run = subprocess.run
        build_commands = []
        def run(command, **kwargs):
            if command[0] != 'cargo':
                return actual_run(command, **kwargs)
            build_commands.append(command)
            target = Path(command[command.index('--target-dir')+1])
            self.assertNotEqual(target, self.root/'target')
            binary = target/'release/trm'; binary.parent.mkdir()
            binary.write_bytes(self.source.read_bytes()); binary.chmod(0o755)
            return subprocess.CompletedProcess(command, 0)
        original = self.old.read_bytes()
        with patch.object(installer, 'ROOT', self.root), patch.object(installer.subprocess, 'run', side_effect=run), \
                patch.object(installer.subprocess, 'check_output', return_value='fixture compiler\n'), \
                patch.object(sys, 'argv', ['install_runtime.py', '--offline', '--jobs', '2']), \
                contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(installer.main(), 0)
        self.assertEqual(len(build_commands), 1)
        self.assertIn('--locked', build_commands[0]); self.assertIn('--offline', build_commands[0])
        self.assertFalse((self.root/'results').exists())
        self.assertEqual(self.old.read_bytes(), original)
        self.assertEqual(runtime.installed_identity(self.root)['version'], '0.8.4')


class RuntimeEntrypointTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='trm-runtime-entrypoint-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.fastq = self.root/'reads.fastq'
        self.fastq.write_text('@r1\n'+'TTAGGG'*20+'\n+\n'+'I'*120+'\n')
        self.manifest = self.root/'samples.tsv'
        self.manifest.write_text('sample_id\tpath1\nfixture\treads.fastq\n')
        self.env = dict(os.environ); self.env.pop('TRM_BINARY', None)

    def run_batch(self, output, *extra, script=None):
        return subprocess.run([sys.executable, '-B', str(script or ROOT/'scripts/batch.py'),
            str(self.manifest), '--output-dir', str(output), '--threads', '1', '--jobs', '1',
            '--min-reads', '1', *extra], cwd=self.root, env=self.env,
            text=True, capture_output=True, timeout=30)

    def test_default_current_from_other_directory(self):
        current = self.root/'current-output'
        result = self.run_batch(current, '--flank-period-annotation', 'on')
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        report = json.loads((current/'fixture/result.json').read_text())
        self.assertEqual(report['version'], runtime.project_version())
        source = ROOT/'bin/trm'
        source_hash = runtime.sha256(source)
        frozen = current/'.batch-binaries'/source_hash/'trm'
        self.assertEqual(Path(report['command'][0]), frozen)
        self.assertEqual(runtime.sha256(frozen), source_hash)
        self.assertFalse(frozen.samefile(source))
        self.assertEqual(frozen.stat().st_mode & 0o222, 0)
        state = json.loads((current/'fixture/state.json').read_text())
        self.assertEqual(state['status'], 'complete')
        self.assertEqual(state['fingerprint']['executable']['path'], str(source))
        self.assertEqual(state['fingerprint']['executable']['sha256'], source_hash)
        self.assertEqual(state['execution_binary']['path'], str(frozen))
        self.assertEqual(state['execution_binary']['source_path'], str(source))
        self.assertEqual(state['execution_binary']['sha256'], source_hash)

    @archive_only
    def test_explicit_historical_override_from_other_directory(self):
        historical = ROOT/'target/release/trm'
        if not historical.is_file() or runtime.binary_version(historical) != '0.8.2':
            self.skipTest('Preserved v0.8.2 scanner is required for historical compatibility')
        previous_hash = runtime.sha256(historical)
        old = self.root/'historical-output'
        result = self.run_batch(old, '--binary', str(historical))
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        old_report = json.loads((old/'fixture/result.json').read_text())
        self.assertEqual(old_report['version'], '0.8.2')
        self.assertEqual(Path(old_report['command'][0]), old/'.batch-binaries'/previous_hash/'trm')
        self.assertEqual(runtime.sha256(Path(old_report['command'][0])), previous_hash)
        old_state = json.loads((old/'fixture/state.json').read_text())
        self.assertEqual(old_state['execution_binary']['source_path'], str(historical))
        self.assertEqual(runtime.sha256(historical), previous_hash)

    def test_fresh_wrapper_missing_runtime_fails_before_creating_output(self):
        project = self.root/'fresh-project'
        scripts = project/'scripts'; scripts.mkdir(parents=True)
        for name in ['batch.py', 'batch_common.py', 'batch_integrity.py', 'runtime_binary.py',
                     'native_diagnostic_options.py', 'quality_observation.py', 'runner_resources.py']:
            shutil.copy2(ROOT/'scripts'/name, scripts/name)
        shutil.copy2(ROOT/'Cargo.toml', project/'Cargo.toml')
        # A historical file exists; selection must still fail instead of using it.
        old = project/'target/release/trm'; old.parent.mkdir(parents=True)
        old.write_text('#!/bin/sh\nprintf "trm 0.8.2\\n"\n'); old.chmod(0o755)
        output = self.root/'should-not-exist'
        result = self.run_batch(output, script=scripts/'batch.py')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('install_runtime.py', result.stderr)
        self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
