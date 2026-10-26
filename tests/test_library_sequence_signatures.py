import gzip
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('sequence_signatures', ROOT / 'scripts/analyze_library_sequence_signatures.py')
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class SequenceSignatureTests(unittest.TestCase):
    def test_full_leader_orientation_and_downstream_context(self):
        motif = audit.LEADERS['SL1']
        context = b'ACGTTCGATAGCTACGATGC'
        seq = b'AA' + motif + context + b'GGG'
        forward = audit.leader_signature(seq, motif)
        reverse = audit.leader_signature(audit.revcomp(seq), motif)
        self.assertEqual(forward['anchored_5prime_records'], 1)
        self.assertEqual(forward['anchored_rc_3prime_records'], 0)
        self.assertEqual(reverse['anchored_rc_3prime_records'], 1)
        self.assertEqual(forward['contexts'], [context.decode()])
        self.assertEqual(reverse['contexts'], forward['contexts'])

    def test_internal_leader_and_insufficient_context_are_separate(self):
        motif = audit.LEADERS['SL2']
        hit = audit.leader_signature(b'AAAA' + motif + b'C' * 30, motif)
        self.assertEqual(hit['anywhere_records'], 1)
        self.assertEqual(hit['anchored_records'], 0)
        hit = audit.leader_signature(motif + b'C' * 19, motif)
        self.assertEqual(hit['anchored_records'], 1)
        self.assertEqual(hit['contexts'], [])

    def test_overlapping_ligation_and_ambiguous_windows(self):
        self.assertEqual(list(audit.positions(b'GATCGATCGATC', audit.LIGATION)), [0, 4])
        self.assertEqual(audit.valid_windows(b'ACGTACGTNACGTACGT', 8), 2)
        self.assertEqual(audit.valid_windows(b'ACGTACGN', 8), 0)

    def test_markov_probability_includes_all_seven_transitions(self):
        bases = {b: 100 for b in 'ACGT'}
        counts = {s: 25 for s in ['GA', 'AT', 'TC', 'CG']}
        self.assertAlmostEqual(audit.expected_probability(bases, counts, bases), 4 ** -8)
        self.assertAlmostEqual(audit.expected_probability(bases, counts, bases, True), 4 ** -8)
        counts['CG'] = 0
        self.assertEqual(audit.expected_probability(bases, counts, bases, True), 0)

    def test_saved_prefix_is_audited_and_edges_do_not_count_internal(self):
        sequences = [audit.LIGATION + b'C' * 30, b'C' * 8 + audit.LIGATION + b'A' * 20,
                     audit.LEADERS['SL1'] + b'ACGT' * 6, b'T' * 20 + b'G' * 20]
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / 'prefix.fastq.gz'
            p.write_bytes(gzip.compress(b''.join(b'@r\n' + seq + b'\n+\n' + b'I' * len(seq) + b'\n'
                                              for seq in sequences)))
            row = {'sample_id': 'x', 'input_slot': 'path1', 'extracted_path': str(p),
                   'extracted_sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
                   'counts': {'records': 4, 'bases': sum(map(len, sequences))},
                   'gc_fraction': 0, 'extra_identical_sequence_fraction': 0,
                   'largest_identical_sequence_count': 1}
            result = audit.inspect_prefix(row)
            self.assertEqual(result['counts']['ligation_records'], 2)
            self.assertEqual(result['counts']['ligation_internal_occurrences'], 1)
            self.assertEqual(result['counts']['any_full_leader_anchored_records'], 1)
            self.assertEqual(result['counts']['terminal_polyA_or_rc_polyT_ge20_records'], 2)
            row['extracted_sha256'] = 'incorrect'
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                audit.inspect_prefix(row)


if __name__ == '__main__':
    unittest.main()
