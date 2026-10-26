"""Native flank annotation forwarding, additive output and legacy-safe resume."""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from archive_test_support import archive_only

from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]
BATCH = ROOT/'scripts/batch.py'
SNAPSHOT = ROOT/'results/evidence_batch_scaling_20261004/before/scripts/batch.py'


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def biological_output(value, strip_annotation=False):
    """Preserve result order, examples, priors and all non-timing observations."""
    result = copy.deepcopy(value)
    result.pop('command', None)
    for key in ('elapsed_seconds', 'megabases_per_second'):
        result['summary'].pop(key, None)
    # ASLR/allocator/process startup differences are not biological results.
    result.get('memory_policy', {}).pop('startup_virtual_bytes', None)
    if strip_annotation:
        result.pop('flank_period_annotation_method', None)
        def strip(node):
            if isinstance(node, dict):
                node.pop('flank_period_annotation', None)
                for child in node.values():
                    strip(child)
            elif isinstance(node, list):
                for child in node:
                    strip(child)
        strip(result)
    return result


class BatchFlankOptionsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binary = runtime_binary()
        if not cls.binary.is_file():
            raise unittest.SkipTest('Build the current scanner or set TRM_BINARY')
        help_text = subprocess.run([str(cls.binary), '--help'], capture_output=True, text=True, check=True).stdout
        if '--flank-period-annotation' not in help_text:
            raise unittest.SkipTest('Current scanner supporting flank-period annotation required; set TRM_BINARY')

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='trm-batch-flank-')
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.output = self.root/'output'
        self.manifest = self.root/'samples.tsv'
        motif = 'TTAGGG'
        flank = 'ACGTCCGATGCTAACGTCAGTGCATCAGACGTTCGACGTA'
        sequences = [motif*12+flank, reverse_complement(motif*12+flank), motif*25, motif*12+flank[:20]]*2
        self.input = self.root/'reads.fastq'
        self.input.write_text(''.join(f'@r{i}\n{seq}\n+\n{"I"*len(seq)}\n' for i, seq in enumerate(sequences)))
        self.manifest.write_text('sample_id\tpath1\nsample\treads.fastq\n')

    def run_batch(self, *options, script=BATCH, binary=None, expected=0):
        result = subprocess.run([sys.executable, '-B', str(script), str(self.manifest),
            '--binary', str(binary or self.binary), '--output-dir', str(self.output),
            '--jobs', '1', '--threads', '1', '--cpu-budget', '1', *options],
            capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, expected, result.stdout+result.stderr)
        return result

    def load(self, name):
        return json.loads((self.output/name).read_text())

    def state(self):
        return self.load('sample/state.json')

    def result(self):
        return self.load('sample/result.json')

    def direct(self, options):
        prefix = self.root/'direct'
        result = subprocess.run([str(self.binary), 'scan', str(self.input), '--output', str(prefix), *options],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        return json.loads(prefix.with_suffix('.json').read_text())

    def test_native_equals_batch_optin_and_existing_observations_preserved(self):
        self.run_batch()
        off = self.result()
        self.assertNotIn('--flank-period-annotation', self.state()['fingerprint']['options'])
        self.assertFalse(self.load('summary.json')['flank_period_annotation'])
        self.run_batch('--flank-period-annotation', 'on', '--resume')
        state, on = self.state(), self.result()
        self.assertEqual(self.load('summary.json')['samples'][0]['status'], 'complete')
        options = state['fingerprint']['options']
        self.assertEqual(options[options.index('--flank-period-annotation')+1], 'on')
        self.assertTrue(state['flank_period_annotation'])
        self.assertTrue(state['metrics']['flank_period_annotation'])
        self.assertTrue(self.load('summary.json')['flank_period_annotation'])
        self.assertTrue(self.load('progress.json')['flank_period_annotation'])
        self.assertEqual(biological_output(on), biological_output(self.direct(options)))
        self.assertEqual(biological_output(off), biological_output(on, strip_annotation=True))
        annotations = [c['junction_evidence']['flank_period_annotation']
            for pool in ('candidates', 'rejected_candidates') for c in on[pool]]
        self.assertTrue(annotations)
        self.assertGreater(sum(a['observed_events'] for a in annotations), 0)

    def test_default_explicit_off_and_same_on_resume_preserve_artifact_bytes(self):
        self.run_batch()
        fingerprint = self.state()['fingerprint']
        artifact = (self.output/'sample/result.json').read_bytes()
        self.run_batch('--flank-period-annotation', 'off', '--resume')
        self.assertEqual(self.load('summary.json')['samples'][0]['status'], 'resumed')
        self.assertEqual(self.state()['fingerprint'], fingerprint)
        self.assertEqual((self.output/'sample/result.json').read_bytes(), artifact)
        self.run_batch('--flank-period-annotation', 'on', '--resume')
        self.assertEqual(self.load('summary.json')['samples'][0]['status'], 'complete')
        artifact = (self.output/'sample/result.json').read_bytes()
        fingerprint = self.state()['fingerprint']
        self.run_batch('--flank-period-annotation', 'on', '--resume')
        self.assertEqual(self.load('summary.json')['samples'][0]['status'], 'resumed')
        self.assertEqual(self.state()['fingerprint'], fingerprint)
        self.assertEqual((self.output/'sample/result.json').read_bytes(), artifact)
        self.run_batch('--flank-period-annotation', 'off', '--resume')
        self.assertEqual(self.load('summary.json')['samples'][0]['status'], 'complete')
        self.assertNotIn('--flank-period-annotation', self.state()['fingerprint']['options'])
        self.assertNotIn('flank_period_annotation', self.state())
        self.assertNotIn('flank_period_annotation_method', self.result())

    @unittest.skipUnless(SNAPSHOT.is_file(), 'Local pre-change batch snapshot not distributed')
    @archive_only
    def test_true_old_batch_state_resumes_without_rewriting_result_or_identity(self):
        for binary in dict.fromkeys([self.binary, (ROOT/'target/release/trm').resolve()]):
            if not binary.is_file():
                continue
            with self.subTest(binary=binary):
                self.output = self.root/binary.name/digest_path(binary)
                self.run_batch(script=SNAPSHOT, binary=binary)
                state = (self.output/'sample/state.json').read_bytes()
                artifact = (self.output/'sample/result.json').read_bytes()
                for options in [('--resume',), ('--resume', '--flank-period-annotation', 'off')]:
                    self.run_batch(*options, binary=binary)
                    self.assertEqual(self.load('summary.json')['samples'][0]['status'], 'resumed')
                    self.assertEqual((self.output/'sample/state.json').read_bytes(), state)
                    self.assertEqual((self.output/'sample/result.json').read_bytes(), artifact)

    def test_invalid_annotation_value_rejected_before_output(self):
        result = self.run_batch('--flank-period-annotation', 'yes', expected=2)
        self.assertIn('invalid choice', result.stderr)
        self.assertFalse(self.output.exists())

    @archive_only
    def test_unsupported_old_scanner_fails_instead_of_silent_unmeasured_result(self):
        old_binary = ROOT/'target/release/trm'
        if not old_binary.is_file():
            self.skipTest('Old frozen scanner is not available')
        help_text = subprocess.run([str(old_binary), '--help'], capture_output=True, text=True, check=True).stdout
        if '--flank-period-annotation' in help_text:
            self.skipTest('Default scanner already supports annotation')
        self.run_batch('--flank-period-annotation', 'on', binary=old_binary, expected=1)
        self.assertEqual(self.state()['status'], 'failed')
        self.assertTrue(self.state()['flank_period_annotation'])
        self.assertFalse(self.load('summary.json')['all_succeeded'])
        self.assertFalse((self.output/'sample/result.json').exists())


def digest_path(path):
    import hashlib
    return hashlib.sha256(str(path).encode()).hexdigest()[:12]


if __name__ == '__main__':
    unittest.main()
