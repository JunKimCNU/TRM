"""Twelve synthetic records validate the registered diagnostic contract."""
import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'tests/fixtures/volvariella/run_volvariella.py'
SPEC = importlib.util.spec_from_file_location('volvariella_runner_fixture', SCRIPT)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


class RunnerTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='trm-volvariella-fixture-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        binary = runtime_binary()
        paths = []
        for mate in (1, 2):
            path = self.root / f'R{mate}.fastq'
            reads = []
            for number in range(6):
                motif = module.TARGETS[number % 2]
                sequence = (motif * 30)[:150]
                if mate == 2: sequence = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
                reads.append(f'@fixture{number}/{mate}\n{sequence}\n+\n' + 'I' * len(sequence) + '\n')
            path.write_text(''.join(reads)); paths.append(path)
        prefix = self.root / 'native'
        completed = subprocess.run([str(binary), 'scan', *map(str, paths), '--output', str(prefix),
            '--threads', '1', '--paired-layout', 'fr', '--flank-period-annotation', 'on'],
            capture_output=True, text=True, timeout=20)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.source = self.root / 'native.json'; self.raw = json.loads(self.source.read_text())
        self.assessment = self.root / 'assessment.json'
        self.assessment.write_text(json.dumps({'sample_id': module.SAMPLE, 'detected_motif': module.TARGETS[1],
            'provenance': {'original_discovery_sha256': module.sha(self.source)},
            'expected_family_observations': [{'motif': 'TTAGGG', 'recovery': 'reported_positive', 'reads': 6}]}))
        self.config = {'sample_id': module.SAMPLE, 'motifs': module.TARGETS, 'max_records': 0,
            'threads': 1, 'memory_limit_mib': 512, 'downloads_enabled': False,
            'runner': module.pin(SCRIPT),
            'dependencies': [module.pin(ROOT / 'scripts' / name) for name in ('evidence_batch.py', 'expected_diagnostics.py')],
            'source': module.pin(self.source), 'assessment': module.pin(self.assessment), 'binary': module.pin(binary),
            'detector_version': self.raw['version'], 'expected_summary': {k: self.raw['summary'][k] for k in ('records', 'bases')},
            'native_config': self.raw['config'], 'junction_options': self.raw['junction_options'],
            'native_inputs': self.raw['inputs'], 'inputs': [module.fingerprint(path) for path in paths],
            'output_root': str(self.root / 'diagnostics')}
        self.config_path = self.root / 'config.json'; self.config_path.write_text(json.dumps(self.config))
        self.config_hash = module.sha(self.config_path)

    def test_full_fixture_passes_production_expected_diagnostics_join(self):
        before = module.sha(self.source)
        status = module.execute(self.config_path, self.config_hash)
        self.assertEqual(status['status'], 'complete')
        self.assertEqual(status['full_input_records'], 12)
        self.assertEqual(status['full_input_bases'], 1800)
        self.assertEqual(status['identity']['max_records'], 0)
        self.assertTrue(status['expected_diagnostics_join_accepted'])
        self.assertEqual(status['downloads_started'], 0)
        self.assertEqual(module.sha(self.source), before)
        folder = Path(self.config['output_root']) / module.SAMPLE
        evidence = json.loads((folder / 'evidence.json').read_text())
        self.assertEqual(evidence['inputs'], self.raw['inputs'])
        self.assertEqual(evidence['junction_options'], self.raw['junction_options'])
        self.assertEqual(set(evidence['target_motifs']), {module.canonical(m) for m in module.TARGETS})
        validation = json.loads((folder / 'validation.json').read_text())
        self.assertTrue(all(row['equal'] for row in validation['representative_control'].values()))
        self.assertTrue((folder / 'resource.tsv').is_file())

    def test_input_change_prevents_scan(self):
        path = Path(self.config['inputs'][0]['path']); path.write_text(path.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, 'fingerprint'):
            module.execute(self.config_path, self.config_hash)
        self.assertFalse(Path(self.config['output_root']).exists())

    def test_source_version_inputs_settings_counts_and_prefix_are_guarded_by_join(self):
        status = module.execute(self.config_path, self.config_hash)
        folder = Path(self.config['output_root']) / module.SAMPLE
        path = folder / 'evidence.json'; original = json.loads(path.read_text())
        changes = [('version', 'different'), ('inputs', []),
                   ('summary', {**original['summary'], 'records': 11}),
                   ('summary', {**original['summary'], 'bases': 1799}),
                   ('summary', {**original['summary'], 'sampled_prefix': True}),
                   ('config', {**original['config'], 'seed_bases': 30}),
                   ('junction_options', {**original['junction_options'], 'end_slack': 4})]
        for field, value in changes:
            with self.subTest(field=field, value=value):
                document = {**original, field: value}; path.write_text(json.dumps(document))
                changed = {**status, 'result_sha256': module.sha(path)}
                with self.assertRaises(ValueError): module.validate_join(self.config, self.raw, folder, changed)
        path.write_text(json.dumps(original))


if __name__ == '__main__': unittest.main()
