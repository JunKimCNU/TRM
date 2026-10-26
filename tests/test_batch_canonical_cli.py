"""The batch wrapper must pass v0.4 exact-catalogue policy to the scanner."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]


class CanonicalBatchPolicyTest(unittest.TestCase):
    def test_canonical_is_forwarded_and_off_catalogue_family_is_retained(self):
        binary = runtime_binary()
        if not binary.is_file():
            self.skipTest('Build trm or set TRM_BINARY')
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            sequences = ['TTAGGG' * 25] * 3 + ['AAAAAAAG' * 19] * 4
            (root / 'input.fq').write_text(''.join(
                f'@r{i}\n{s}\n+\n{"I" * len(s)}\n' for i, s in enumerate(sequences)))
            (root / 'samples.tsv').write_text('sample_id\tpath1\tpath2\nsample\tinput.fq\t\n')
            result = subprocess.run([
                sys.executable, str(ROOT / 'scripts/batch.py'), str(root / 'samples.tsv'),
                '--binary', str(binary), '--output-dir', str(root / 'output'),
                '--jobs', '1', '--threads', '1', '--cpu-budget', '1',
                '--prior-policy', 'canonical'], capture_output=True, text=True, timeout=15)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            summary = json.loads((root / 'output/summary.json').read_text())
            report = json.loads(Path(summary['samples'][0]['result_json']).read_text())
            self.assertEqual(report['prior_options']['policy'], 'canonical')
            self.assertEqual(report['summary']['records'], 7)
            self.assertTrue(any(r['motif'] == 'AACCCT' for r in report['candidates']))
            self.assertTrue(any(r['motif'] == 'AAAAAAAG' for r in report['rejected_candidates']))


if __name__ == '__main__':
    unittest.main()
