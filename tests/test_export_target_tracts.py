"""Optional all-tract exporter: real coordinate, polarity, and failure invariants.

Build first: cargo build --release --example export_target_tracts
Override executable with TRM_EXPORT_TRACTS_BINARY. No real biological input is
used here, and passing these tests does not validate a subtelomeric mapping call.
"""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ExportTargetTractsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binary = Path(os.environ.get(
            'TRM_EXPORT_TRACTS_BINARY', ROOT / 'target/release/examples/export_target_tracts')).resolve()
        if not cls.binary.is_file():
            if 'TRM_EXPORT_TRACTS_BINARY' in os.environ:
                raise RuntimeError(f'Exporter binary missing: {cls.binary}')
            raise unittest.SkipTest('Optional exporter not built; cargo build --release --example export_target_tracts')

    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix='trm-export-target-tests-')
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(tmp.name)
        self.runs = 0

    def fasta(self, sequences, ids=None, name='input.fasta'):
        ids = ids or [f'read{i}' for i in range(len(sequences))]
        p = self.tmp / name
        p.write_text(''.join(f'>{name}\n{sequence}\n' for name, sequence in zip(ids, sequences)))
        return p

    def run_export(self, path, *options, motif='TTAGGG', success=True, output=None):
        self.runs += 1
        self.output = output or self.tmp / f'out{self.runs}'
        completed = subprocess.run(
            [str(self.binary), str(path), '--motif', motif, '--output-dir', str(self.output), *options],
            capture_output=True, text=True, timeout=20)
        if not success:
            self.assertNotEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue(completed.stderr.strip())
            self.assertFalse((self.output / 'summary.json').exists())
            return completed
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.summary = json.loads((self.output / 'summary.json').read_text())
        with (self.output / 'tracts.tsv').open() as handle:
            self.rows = list(csv.DictReader(handle, delimiter='\t'))
        lines = (self.output / 'anchors.fasta').read_text().splitlines()
        self.anchors = {lines[i][1:]: lines[i+1] for i in range(0, len(lines), 2)}
        self.assertEqual(len(lines), 2 * len(self.anchors))
        self.assertEqual(len(self.rows), self.summary['tracts'])
        self.assertEqual(len(self.anchors), self.summary['anchors'])
        self.assertEqual(self.summary['status'], 'complete')
        self.assertTrue(self.summary['validated_eof'])
        self.assertTrue(self.summary['complete_input'])
        for name, fingerprint in self.summary['outputs'].items():
            data = (self.output / name).read_bytes()
            self.assertEqual(fingerprint['size_bytes'], len(data))
            self.assertEqual(fingerprint['sha256'], hashlib.sha256(data).hexdigest())
        self.assertFalse(list(self.output.glob('*.partial')))
        return self.rows

    def assert_anchor_slice(self, row, side, sequence):
        lo, hi = int(row[side + '_anchor_start0']), int(row[side + '_anchor_end0'])
        self.assertEqual(int(row[side + '_anchor_length']), hi-lo)
        if row[side + '_anchor_status'] == 'exported':
            self.assertEqual(self.anchors[row[side + '_anchor_id']], sequence[lo:hi])
            self.assertGreaterEqual(hi-lo, 200)
        else:
            self.assertEqual(row[side + '_anchor_id'], '')
        start, end = int(row['start0']), int(row['end0'])
        self.assertTrue(hi <= start or lo >= end)

    def test_internal_tract_coordinates_and_two_real_flank_sequences(self):
        sequence = 'C'*250 + 'TTAGGG'*50 + 'C'*260
        row, = self.run_export(self.fasta([sequence]))
        self.assertEqual((int(row['start0']), int(row['end0'])), (250, 550))
        self.assertEqual((row['strand'], row['tract_length'], row['identity']), ('+', '300', '1.000000000000'))
        self.assertEqual((row['left_end_distance'], row['right_end_distance']), ('250', '260'))
        self.assert_anchor_slice(row, 'left', sequence)
        self.assert_anchor_slice(row, 'right', sequence)
        self.assertEqual(self.summary['tracts_with_anchor'], 1)

    def test_reverse_strand_is_relative_to_supplied_outward_motif(self):
        sequence = 'G'*260 + 'CCCTAA'*100 + 'G'*300
        path = self.fasta([sequence])
        row, = self.run_export(path)
        self.assertEqual((row['start0'], row['end0'], row['strand']), ('260', '860', '-'))
        self.assert_anchor_slice(row, 'left', sequence)
        self.assert_anchor_slice(row, 'right', sequence)
        same, = self.run_export(path, motif='CCCTAA')
        self.assertEqual(same['strand'], '+')

    def test_terminal_tract_and_short_seed_are_never_removed_by_500bp_or_minreads(self):
        sequence = 'TTAGGG'*4 + 'C'*240
        row, = self.run_export(self.fasta([sequence]))
        self.assertEqual((row['start0'], row['end0'], row['tract_length']), ('0', '24', '24'))
        self.assertEqual(row['left_anchor_status'], 'below_minimum')
        self.assertEqual(row['left_anchor_stop_reason'], 'read_boundary')
        self.assertEqual(row['right_anchor_status'], 'exported')
        self.assertFalse(self.summary['length_filter_applied'])
        self.assertFalse(self.summary['minimum_read_support_filter_applied'])
        self.assertEqual(self.summary['config']['seed_bases'], 24)
        self.assertEqual(self.summary['config']['min_copies'], 3)
        self.assertEqual(self.summary['config']['min_identity'], .85)

    def test_below_minimum_anchors_still_keep_the_tract_and_lengths(self):
        row, = self.run_export(self.fasta(['C'*199 + 'TTAGGG'*50 + 'C'*199]))
        self.assertEqual(row['left_anchor_length'], '199')
        self.assertEqual(row['right_anchor_length'], '199')
        self.assertEqual(self.summary['tracts_without_anchor'], 1)
        self.assertEqual(self.anchors, {})

    def test_ambiguous_bases_bound_anchors_without_skipping_to_distant_sequence(self):
        sequence = 'C'*250 + 'N' + 'C'*220 + 'TTAGGG'*50 + 'C'*250 + 'R' + 'C'*300
        row, = self.run_export(self.fasta([sequence]))
        self.assertEqual((row['left_anchor_length'], row['right_anchor_length']), ('220', '250'))
        for side in ('left', 'right'):
            self.assertEqual(row[side + '_anchor_stop_reason'], 'ambiguous_base')
            self.assert_anchor_slice(row, side, sequence)

    def test_adjacent_same_family_tracts_bound_flanks_even_below_500bp(self):
        sequence = 'TTAGGG'*40 + 'C'*300 + 'TTAGGG'*40
        first, second = self.run_export(self.fasta([sequence]))
        self.assertEqual(first['right_anchor_stop_reason'], 'same_family_tract')
        self.assertEqual(second['left_anchor_stop_reason'], 'same_family_tract')
        self.assertEqual((first['right_anchor_start0'], first['right_anchor_end0']), ('240', '540'))
        self.assertEqual((second['left_anchor_start0'], second['left_anchor_end0']), ('240', '540'))
        self.assertEqual(self.summary['max_tract_bases'], 240)
        for row in (first, second):
            for side in ('left', 'right'):
                self.assert_anchor_slice(row, side, sequence)

    def test_touching_opposite_polarity_tracts_do_not_export_repeat_as_anchor(self):
        sequence = 'TTAGGG'*50 + 'CCCTAA'*50
        rows = self.run_export(self.fasta([sequence]))
        self.assertEqual(len(rows), 2)
        self.assertEqual({r['strand'] for r in rows}, {'+', '-'})
        self.assertEqual(self.anchors, {})
        self.assertEqual(rows[0]['right_anchor_stop_reason'], 'same_family_tract')
        self.assertEqual(rows[1]['left_anchor_stop_reason'], 'same_family_tract')

    def test_symmetric_family_has_explicit_ambiguous_polarity(self):
        row, = self.run_export(self.fasta(['N'*10 + 'ACGT'*10 + 'N'*20]), motif='ACGT')
        self.assertEqual(row['strand'], 'ambiguous')
        self.assertEqual((row['start0'], row['end0']), ('10', '50'))

    def test_duplicate_input_ids_still_have_unique_safe_query_and_tract_ids(self):
        sequence = 'C'*210 + 'TTAGGG'*5 + 'C'*210
        rows = self.run_export(self.fasta([sequence, sequence], ids=['same%id/1 comment', 'same%id/1 other']))
        self.assertEqual({r['read_id'] for r in rows}, {'same%25id/1'})
        self.assertEqual([r['read_ordinal'] for r in rows], ['1', '2'])
        self.assertEqual(len({r['tract_id'] for r in rows}), 2)
        self.assertEqual(len(self.anchors), 4)
        for name in self.anchors:
            self.assertRegex(name, r'^r\d{12}_t\d{6}_[LR]$')

    def test_configured_maximum_bounds_long_flanks(self):
        sequence = 'C'*1500 + 'TTAGGG'*100 + 'C'*1500
        row, = self.run_export(self.fasta([sequence]), '--max-anchor-bp', '512')
        for side in ('left', 'right'):
            self.assertEqual(row[side + '_anchor_length'], '512')
            self.assertEqual(row[side + '_anchor_stop_reason'], 'maximum_length')
            self.assert_anchor_slice(row, side, sequence)

    def test_gzip_magic_fastq_full_eof_and_negative_records_count(self):
        positive = 'TTAGGG'*10
        negative = 'C'*80
        raw = ''.join(f'@r{i}\n{s}\n+\n'+ 'I'*len(s)+'\n' for i, s in enumerate([negative, positive, negative]))
        path = self.tmp/'compressed_without_suffix.data'
        path.write_bytes(gzip.compress(raw.encode()))
        row, = self.run_export(path)
        self.assertEqual((row['read_ordinal'], self.summary['records'], self.summary['bases']), ('2', 3, 220))
        self.assertEqual(self.summary['input']['format'], 'fastq')
        self.assertEqual(self.summary['target_reads'], 1)

    def test_valid_no_target_is_explicit_complete_zero_not_failure(self):
        rows = self.run_export(self.fasta(['C'*300, 'N'*50]))
        self.assertEqual(rows, [])
        self.assertEqual(self.summary['target_reads'], 0)
        self.assertEqual(self.summary['records'], 2)
        self.assertEqual(self.summary['max_tract_bases'], 0)

    def test_truncated_fastq_after_valid_hit_has_no_success_outputs(self):
        path = self.tmp/'broken.fastq'
        path.write_text('@good\n'+'TTAGGG'*5+'\n+\n'+'I'*30+'\n@bad\nACGT\n+\nII\n')
        self.run_export(path, success=False)
        self.assertFalse((self.output/'tracts.tsv').exists())
        self.assertFalse((self.output/'anchors.fasta').exists())
        self.assertTrue((self.output/'tracts.tsv.partial').exists())

    def test_bad_gzip_crc_has_no_success_summary(self):
        path = self.tmp/'bad.gz'
        raw = bytearray(gzip.compress(b'>read\n'+b'TTAGGG'*10+b'\n'))
        raw[-8] ^= 1
        path.write_bytes(raw)
        self.run_export(path, success=False)

    def test_existing_output_is_preserved(self):
        output = self.tmp/'existing';output.mkdir()
        sentinel = output/'marker';sentinel.write_text('keep')
        self.run_export(self.fasta(['TTAGGG'*10]), output=output, success=False)
        self.assertEqual(list(output.iterdir()), [sentinel])
        self.assertEqual(sentinel.read_text(), 'keep')

    def test_missing_hash_program_cannot_publish_success_summary(self):
        path = self.fasta(['TTAGGG'*10])
        output = self.tmp/'no-hash-program'
        completed = subprocess.run(
            [str(self.binary), str(path), '--motif', 'TTAGGG', '--output-dir', str(output)],
            capture_output=True, text=True, timeout=20, env={**os.environ, 'PATH': ''})
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn('sha256sum', completed.stderr)
        self.assertFalse((output/'summary.json').exists())

    def test_invalid_options_empty_and_invalid_sequence_are_technical_errors(self):
        path = self.fasta(['TTAGGG'*10])
        for options, motif in [((), 'TTAGGGTTAGGG'), ((), 'TTNGGG'), ((), 'AT'),
                               (('--min-anchor-bp', '199'), 'TTAGGG'),
                               (('--max-anchor-bp', '1001'), 'TTAGGG')]:
            with self.subTest(options=options, motif=motif):
                self.run_export(path, *options, motif=motif, success=False)
                self.assertFalse(self.output.exists())
        for i, data in enumerate(['', '>bad\nACG?\n']):
            path = self.tmp/f'bad{i}.fasta';path.write_text(data)
            self.run_export(path, success=False)


if __name__ == '__main__':
    unittest.main()
