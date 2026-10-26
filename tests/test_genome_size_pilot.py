"""Tests for genome-profile provenance and failure handling, without large data."""
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('genome_size_pilot', ROOT / 'scripts/run_genome_size_pilot.py')
PILOT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PILOT)


class GenomeSizePilotTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write_fit(self, status='done', sizes='10,000 bp 12,000 bp'):
        self.root.joinpath('progress.txt').write_text(f'round 0 converged\n{status}\n')
        self.root.joinpath('summary.txt').write_text(f'GenomeScope version 2.0\nGenome Haploid Length     {sizes}\n')
        self.root.joinpath('model.txt').write_text('Nonlinear regression model\n')

    def test_positive_fit_requires_terminal_success(self):
        self.write_fit()
        self.assertTrue(PILOT.model_converged(self.root))
        self.root.joinpath('progress.txt').write_text('done\nround 1 unconverged\nfail\n')
        self.assertFalse(PILOT.model_converged(self.root))

    def test_summary_file_does_not_establish_convergence(self):
        self.write_fit(status='fail', sizes='-1 bp -1 bp')
        self.assertFalse(PILOT.model_converged(self.root))

    def test_missing_or_empty_model_is_not_converged(self):
        self.write_fit()
        self.root.joinpath('model.txt').unlink()
        self.assertFalse(PILOT.model_converged(self.root))
        self.root.joinpath('model.txt').write_text('  \n')
        self.assertFalse(PILOT.model_converged(self.root))

    def test_invalid_or_nonpositive_size_is_not_converged(self):
        for sizes in ('0 bp 2 bp', '-1 bp 2 bp', 'NaN bp 2 bp', 'Inf bp Inf bp', '1e999 bp 2 bp'):
            with self.subTest(sizes=sizes):
                self.write_fit(sizes=sizes)
                self.assertFalse(PILOT.model_converged(self.root))

    def test_missing_summary_or_progress_is_not_converged(self):
        self.assertFalse(PILOT.model_converged(self.root))
        self.write_fit()
        self.root.joinpath('progress.txt').unlink()
        self.assertFalse(PILOT.model_converged(self.root))

    def test_histogram_records_weighted_tail_and_saturation(self):
        path = self.root / 'histogram'
        counts = {1: 10, 20: 5, 21: 3, 50: 2, 51: 4, 100: 1}
        path.write_text(''.join(f'{i}\t{counts.get(i, 0)}\n' for i in range(1, 101)))
        summary = PILOT.summarize_histogram(path, 100, [20, 50])
        self.assertEqual(summary['distinct_kmers'], 25)
        self.assertEqual(summary['counted_kmer_occurrences_lower_bound'], 577)
        self.assertEqual(summary['occurrences_above_fit_cutoff'], {'20': 467, '50': 304})
        self.assertAlmostEqual(summary['occurrence_fraction_above_fit_cutoff']['20'], 467 / 577)
        self.assertTrue(summary['saturated'])
        self.assertEqual(summary['kmers_at_counter_cap'], 1)

    def test_empty_mass_does_not_divide_by_zero(self):
        path = self.root / 'histogram'
        path.write_text('1\t0\n100\t0\n')
        summary = PILOT.summarize_histogram(path, 100, [20])
        self.assertFalse(summary['saturated'])
        self.assertIsNone(summary['occurrence_fraction_above_fit_cutoff']['20'])

    def test_pilot_preserves_full_histogram_and_records_actual_fit_status(self):
        files = [self.root / 'read1.fastq', self.root / 'read2.fastq']
        for file in files:
            file.write_text('@r\nACGT\n+\nIIII\n')
        originals = [file.read_bytes() for file in files]
        args = SimpleNamespace(output=self.root / 'output', tool_bin=self.root / 'bin',
                               threads=2, memory_gb=2, counter_cap=100,
                               cutoffs=[20, 50], ploidies=[1, 2])
        commands = []
        histogram = ''.join(f'{i}\t{1 if i in (1, 20, 50, 100) else 0}\n' for i in range(1, 101))

        def external(command, stem, env):
            commands.append(command)
            if Path(command[0]).name == 'kmc_tools':
                Path(command[command.index('histogram') + 1]).write_text(histogram)
            elif Path(command[0]).name == 'Rscript':
                directory = Path(command[command.index('-o') + 1])
                ploidy = int(command[command.index('-p') + 1])
                directory.joinpath('summary.txt').write_text('Genome Haploid Length 1,000 bp 2,000 bp\n')
                directory.joinpath('model.txt').write_text('Nonlinear regression model\n')
                directory.joinpath('progress.txt').write_text('done\n' if ploidy == 1 else 'fail\n')
            return {'returncode': 0}

        row = {'sample_id': 'TEST', 'path1': str(files[0]), 'path2': str(files[1])}
        with patch.object(PILOT, 'measured_run', side_effect=external):
            result = PILOT.run_one(row, 21, args, {})
        out = args.output / 'TEST' / 'k21'
        self.assertEqual((out / 'counts.histo').read_text(), histogram)
        for cutoff in args.cutoffs:
            lines = (out / f'fit_input_max{cutoff}.histo').read_text().splitlines()
            self.assertEqual(len(lines), cutoff + 1)
            self.assertEqual(lines[-1], f'{cutoff + 1}\t0')
            self.assertEqual(lines[-2], f'{cutoff}\t1')
        self.assertEqual([f['model_converged'] for f in result['fits']], [True, False, True, False])
        self.assertTrue(result['all_local_reads'])
        self.assertEqual([file.read_bytes() for file in files], originals)
        self.assertEqual((out / 'inputs.txt').read_text().splitlines(), [str(f.resolve()) for f in files])
        self.assertIn('-ci1', commands[0])
        self.assertIn('-cs100', commands[0])
        self.assertIn('-cx1000000000', commands[0])
        self.assertNotIn('-b', commands[0])
        self.assertEqual(json.loads((out / 'result.json').read_text())['sample_id'], 'TEST')


if __name__ == '__main__':
    unittest.main()
