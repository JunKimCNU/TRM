"""Saved-query grammar contracts against independently enumerated word factors."""
import itertools
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import saved_flank_grammar as grammar

RC = str.maketrans('ACGT', 'TGCA')


def reverse_complement(sequence):
    return sequence.translate(RC)[::-1]


def finite_factors(units, maximum):
    """At most L+1 nonempty blocks cover every factor of length at most L.

    Explicitly concatenate every such word and enumerate its substrings; this
    oracle contains no graph state, recurrence, or earliest-start calculation.
    """
    factors = set()
    for blocks in itertools.product(units, repeat=maximum + 1):
        text = ''.join(blocks)
        for start in range(len(text)):
            factors.update(text[start:end] for end in
                           range(start + 1, min(len(text), start + maximum) + 1))
    return factors


def oracle(sequence, factors, cutoff):
    n = len(sequence)
    full_strands, maximal = [], []
    for strand, viewed in (('+', sequence), ('-', reverse_complement(sequence))):
        if viewed in factors:
            full_strands.append(strand)
        all_matches = [(a, b) for a in range(n) for b in range(a + 1, n + 1)
                       if viewed[a:b] in factors]
        for a, b in all_matches:
            if any(x <= a and b <= y and (a, b) != (x, y)
                   for x, y in all_matches):
                continue
            x, y = (a, b) if strand == '+' else (n-b, n-a)
            maximal.append((x, y, strand))
    selected = sorted((a, b, strand) for a, b, strand in maximal if b-a >= cutoff)
    coverage = {position for a, b, _ in selected for position in range(a, b)}
    uncovered = set(range(n)) - coverage
    outside_runs = [b-a for a in range(n) for b in range(a + 1, n + 1)
                    if set(range(a, b)) <= uncovered]
    return dict(
        full_grammar_compatible=bool(full_strands), full_grammar_strands=full_strands,
        longest_grammar_span_bp=max((b-a for a, b, _ in maximal), default=0),
        grammar_union_bp=len(coverage),
        prefix_outside_bp=min(coverage) if coverage else n,
        suffix_outside_bp=n-1-max(coverage) if coverage else n,
        max_contiguous_outside_bp=max(outside_runs, default=0),
        intervals=[dict(start0=a, end0=b, strand=strand) for a, b, strand in selected],
    )


class SavedFlankGrammarTests(unittest.TestCase):
    def test_exhaustive_short_queries_against_finite_word_enumeration(self):
        for units in (['AC', 'GTT'], ['A', 'CC'], ['AG', 'GA']):
            factors = finite_factors(units, 5)
            pure_factors = {unit: finite_factors([unit], 5) for unit in units}
            for length in range(1, 6):
                for letters in itertools.product('ACGT', repeat=length):
                    sequence = ''.join(letters)
                    cutoff = 1 + (sum(map(ord, sequence)) % 6)
                    result = grammar.analyze(sequence, units, cutoff)
                    expected = oracle(sequence, factors, cutoff)
                    for key, value in expected.items():
                        self.assertEqual(result[key], value, (units, sequence, cutoff, key))
                    compatible = [unit for unit in units if any(
                        view in pure_factors[unit] for view in
                        (sequence, reverse_complement(sequence)))]
                    self.assertEqual(result['full_compatible_units'], compatible)
                    self.assertEqual(result['full_single_unit_compatible'], bool(compatible))
                    expected_class = ('single_unit_compatible' if compatible else
                                      'mixed_only' if expected['full_grammar_compatible'] else
                                      'not_full_compatible')
                    self.assertEqual(result['compatibility_class'], expected_class)

    def test_real_mixed_query_partial_edges_and_reverse_coordinates(self):
        sequence = 'TTAGGTTAGGGTTAGGTTAGGTTAGGTTAG'
        result = grammar.analyze(sequence, ['TTAGG', 'TTAGGG'])
        self.assertEqual(len(sequence), 30)
        self.assertEqual(result['compatibility_class'], 'mixed_only')
        self.assertEqual(result['full_grammar_strands'], ['+'])
        self.assertEqual(result['grammar_union_bp'], 30)
        self.assertEqual(result['intervals'], [dict(start0=0, end0=30, strand='+')])
        self.assertEqual(result['max_contiguous_outside_bp'], 0)
        reverse = grammar.analyze(reverse_complement(sequence), ['TTAGG', 'TTAGGG'])
        self.assertEqual(reverse['full_grammar_strands'], ['-'])
        self.assertEqual(reverse['intervals'], [dict(start0=0, end0=30, strand='-')])
        # Cutting the mixed array at arbitrary phases must remain compatible.
        for a in range(6):
            for b in range(20, len(sequence)+1):
                self.assertTrue(grammar.analyze(sequence[a:b], ['TTAGG', 'TTAGGG'])
                                ['full_grammar_compatible'])

    def test_cutoff_does_not_gate_full_query_compatibility(self):
        sequence = ('TTAGGTTAGGG' * 3)[2:24]
        below = grammar.analyze(sequence, ['TTAGG', 'TTAGGG'], 20)
        above = grammar.analyze(sequence, ['TTAGG', 'TTAGGG'], 24)
        self.assertTrue(below['full_grammar_compatible'])
        self.assertTrue(above['full_grammar_compatible'])
        self.assertEqual(below['compatibility_class'], 'mixed_only')
        self.assertEqual(below['grammar_union_bp'], 22)
        self.assertEqual(above['grammar_union_bp'], 0)
        self.assertEqual(above['longest_grammar_span_bp'], 22)
        self.assertEqual(above['prefix_outside_bp'], 22)
        self.assertEqual(above['suffix_outside_bp'], 22)
        self.assertEqual(above['max_contiguous_outside_bp'], 22)

    def test_union_and_nonoverlapping_geometry_do_not_imply_nonrepeat(self):
        sequence = 'C' * 20 + 'TTAGG' * 6 + 'A' * 20
        result = grammar.analyze(sequence, ['TTAGG', 'TTAGGG'])
        self.assertEqual(result['grammar_union_bp'], 30)
        self.assertEqual(result['prefix_outside_bp'], 20)
        self.assertEqual(result['suffix_outside_bp'], 20)
        self.assertEqual(result['max_contiguous_outside_bp'], 20)
        reverse = grammar.analyze(reverse_complement(sequence), ['TTAGG', 'TTAGGG'])
        self.assertEqual(reverse['intervals'], [dict(start0=20, end0=50, strand='-')])
        # Palindromic compatibility contributes two strands, one base union.
        palindrome = grammar.analyze('AT' * 20, ['AT'], 20)
        self.assertEqual(palindrome['full_grammar_strands'], ['+', '-'])
        self.assertEqual(palindrome['grammar_union_bp'], 40)
        self.assertEqual(len(palindrome['intervals']), 2)

    def test_literal_block_phase_is_not_independently_normalized(self):
        # ACT and TGA are valid literal blocks; canonicalizing them changes
        # their junction language even though their single-unit phases match.
        original = grammar.analyze('ACTTGAACTTGA', ['ACT', 'TGA'], 1)
        changed = grammar.analyze('ACTTGAACTTGA', ['ACT', 'ATG'], 1)
        self.assertTrue(original['full_grammar_compatible'])
        self.assertFalse(changed['full_grammar_compatible'])

    def test_pure_units_short_and_maximum_length_queries(self):
        for unit in ('TTAGG', 'TTAGGG'):
            for phase in range(len(unit)):
                sequence = (unit * 20)[phase:phase+80]
                result = grammar.analyze(sequence, ['TTAGG', 'TTAGGG'])
                self.assertEqual(result['full_compatible_units'], [unit])
                self.assertEqual(result['compatibility_class'], 'single_unit_compatible')
                self.assertEqual(result['grammar_union_bp'], 80)
        result = grammar.analyze('T', ['TTAGG', 'TTAGGG'])
        self.assertTrue(result['full_single_unit_compatible'])
        self.assertEqual(result['grammar_union_bp'], 0)

    def test_invalid_input_boundaries_rejected(self):
        for sequence in ('', 'A' * 81, 'ACGN', 'acgt', b'ACGT', None, 10):
            with self.subTest(sequence=sequence), self.assertRaises(ValueError):
                grammar.analyze(sequence, ['TTAGG', 'TTAGGG'])
        for units in (None, (), [], [''], ['N'], ['ttagg'], ['A', 'A'],
                      ['A' * 65], ['A' * 64, 'C' * 64, 'G' * 64, 'T' * 64,
                                      'AC' * 32, 'CA' * 32, 'AG' * 32, 'GA' * 32, 'TC' * 32]):
            with self.subTest(units=units), self.assertRaises(ValueError):
                grammar.analyze('ACGT', units)
        for cutoff in (0, 81, True, 1.5, None):
            with self.subTest(cutoff=cutoff), self.assertRaises(ValueError):
                grammar.analyze('ACGT', ['TTAGG', 'TTAGGG'], cutoff)


if __name__ == '__main__':
    unittest.main()
