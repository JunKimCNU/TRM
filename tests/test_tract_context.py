#!/usr/bin/env python3
"""Frozen development controls and independent repeat-coverage property tests.

Synthetic results do not estimate biological sensitivity or specificity.
Run with --write-outcomes PATH to write the frozen pair-control observations.
"""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
CONTROL_DIR = ROOT / 'tests/fixtures/tract_context'
SPEC = importlib.util.spec_from_file_location('tract_context', ROOT / 'scripts/tract_context.py')
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def rc(sequence):
    return sequence.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def coverage_oracle(sequence, start, end, motif, cutoff, strand=None):
    """Mark bases of matching fixed-length windows; no maximal-run enumeration."""
    templates = [motif] if strand == '+' else [rc(motif)] if strand == '-' else [motif, rc(motif)]
    words = {''.join(template[(j + phase) % len(template)] for j in range(cutoff))
             for template in templates for phase in range(len(template))}
    covered = set()
    for i in range(start, end - cutoff + 1):
        if sequence[i:i + cutoff].upper() in words:
            covered.update(range(i, i + cutoff))
    return covered


def intervals_as_set(intervals):
    return {base for start, end in intervals for base in range(start, end)}


def control_fixtures():
    return json.loads((CONTROL_DIR / 'protocol.json').read_text())['pair_controls']


def evaluate_fixture(case):
    return M.evaluate_pair(case['sequence'], copy.deepcopy(case['left']), copy.deepcopy(case['right']),
                           case['motif'], config=copy.deepcopy(case.get('config')),
                           alignments=copy.deepcopy(case.get('alignments')),
                           expected_strand=case.get('expected_strand'))


def assert_no_promotion(testcase, item):
    """Check every explicit no-promotion flag, including nested coverage views."""
    if isinstance(item, dict):
        for key, value in item.items():
            if key in {'biological_tracts_merged', 'changes_native_classification',
                       'physical_chromosome_end_confirmed'}:
                testcase.assertIs(value, False, key)
            assert_no_promotion(testcase, value)
    elif isinstance(item, list):
        for value in item:
            assert_no_promotion(testcase, value)


class FrozenPairControls(unittest.TestCase):
    def test_protocol_unchanged(self):
        receipt = json.loads((CONTROL_DIR / 'protocol.receipt.json').read_text())
        actual = hashlib.sha256((CONTROL_DIR / 'protocol.json').read_bytes()).hexdigest()
        self.assertEqual(receipt['protocol_sha256'], actual)
        self.assertTrue(receipt['frozen_before_first_module_run'])
        protocol = json.loads((CONTROL_DIR / 'protocol.json').read_text())
        self.assertEqual(M.DEFAULT_PAIR_CONFIG, protocol['settings'])

    def test_frozen_control_outcomes(self):
        for case in control_fixtures():
            with self.subTest(control=case['id']):
                result = evaluate_fixture(case)
                self.assertEqual(result['status'], case['expected_status'])
                self.assertTrue(set(case.get('required_reason_codes', ())).issubset(result['reason_codes']))
                self.assertIs(result['changes_native_classification'], False)
                self.assertIs(result['biological_tracts_merged'], False)
                assert_no_promotion(self, result)
                self.assertIs(result['terminal_context']['physical_chromosome_end_confirmed'], False)
                if case.get('expected_alignment_status'):
                    self.assertEqual(result['alignment_evidence']['status'], case['expected_alignment_status'])
                if case.get('expected_terminal_status'):
                    self.assertEqual(result['terminal_context']['status'], case['expected_terminal_status'])

    def test_input_objects_not_mutated(self):
        for case in control_fixtures():
            case['config'] = {'max_gap_bp': 64}
            before = copy.deepcopy(case)
            M.evaluate_pair(case['sequence'], case['left'], case['right'], case['motif'],
                            config=case.get('config'), alignments=case.get('alignments'))
            self.assertEqual(case, before)

    def test_pair_chain_is_not_transitively_joined(self):
        by_id = {case['id']: case for case in control_fixtures()}
        for case_id in ['chain_ab', 'chain_bc']:
            self.assertEqual(evaluate_fixture(by_id[case_id])['status'], 'compatible_under_bounded_model')
        result = evaluate_fixture(by_id['chain_ac'])
        self.assertEqual(result['status'], 'rejected_by_bounded_model')
        self.assertIn('gap_exceeds_max_bp', result['reason_codes'])
        assert_no_promotion(self, result)

    def test_reverse_complement_preserves_bounded_sequence_outcome(self):
        # Coordinates and declared biological orientation must both be reversed.
        for case_id in ['perfect_six', 'single_substitution', 'nonrepeat_spacer', 'phase_shift_insertion']:
            case = next(x for x in control_fixtures() if x['id'] == case_id)
            original = evaluate_fixture(case)
            n = len(case['sequence'])
            def flipped(anchor):
                return dict(tract_id=anchor['tract_id'], start0=n-anchor['end0'],
                            end0=n-anchor['start0'], strand={'-': '+', '+': '-'}[anchor['strand']])
            reversed_result = M.evaluate_pair(rc(case['sequence']), flipped(case['right']),
                                              flipped(case['left']), case['motif'])
            self.assertEqual(original['status'], reversed_result['status'])
            self.assertEqual(set(original['reason_codes']), set(reversed_result['reason_codes']))

    def test_alignment_does_not_change_sequence_compatibility(self):
        cases = {case['id']: case for case in control_fixtures()}
        plain = evaluate_fixture(cases['perfect_six'])
        split = evaluate_fixture(cases['split_mapping'])
        self.assertEqual(split['status'], plain['status'])
        self.assertEqual(split['reason_codes'], plain['reason_codes'])
        self.assertFalse(split['alignment_evidence']['localization_confirmation'])

    def test_invalid_pair_metadata_and_configuration(self):
        case = next(x for x in control_fixtures() if x['id'] == 'perfect_six')
        for config in [{'min_gap_coverage': float('nan')}, {'max_gap_bp': True},
                       {'require_same_phase': 1}, {'unknown_option': 2}]:
            with self.subTest(config=config), self.assertRaises(ValueError):
                M.evaluate_pair(case['sequence'], case['left'], case['right'], case['motif'], config=config)
        invalid = dict(case['left'], tract_id='')
        with self.assertRaises(ValueError):
            M.evaluate_pair(case['sequence'], invalid, case['right'], case['motif'])
        for value in [float('nan'), float('inf'), -0.1, 1.1, True]:
            with self.subTest(identity=value), self.assertRaises(ValueError):
                M.evaluate_pair(case['sequence'], dict(case['left'], identity=value), case['right'], case['motif'])

    def test_same_query_and_chromosome_different_loci_remain_competing(self):
        case = next(x for x in control_fixtures() if x['id'] == 'perfect_six')
        first = dict(query_start0=0, query_end0=150, target='chr1', strand='+', mapq=60, tp='P',
                     target_start0=100, target_end0=250, target_length=1000)
        second = dict(first, target_start0=500, target_end0=650)
        result = M.evaluate_pair(case['sequence'], case['left'], case['right'], case['motif'],
                                 alignments=[first, second])
        evidence = result['alignment_evidence']
        self.assertEqual(evidence['distinct_primary_count'], 2)
        self.assertTrue(evidence['competing_primary_mapping'])
        self.assertEqual(evidence['status'], 'split_or_competing_primary_mapping')
        self.assertEqual(result['status'], 'compatible_under_bounded_model')
        self.assertFalse(evidence['localization_confirmation'])

    def test_exact_positioned_alignment_duplicates_count_once(self):
        case = next(x for x in control_fixtures() if x['id'] == 'perfect_six')
        positioned = dict(query_start0=0, query_end0=150, target='chr1', strand='+', mapq=60, tp='P',
                          target_start0=100, target_end0=250, target_length=1000)
        evidence = M.evaluate_pair(case['sequence'], case['left'], case['right'], case['motif'],
                                   alignments=[positioned, dict(positioned)])['alignment_evidence']
        self.assertEqual(evidence['distinct_primary_count'], 1)
        self.assertEqual(evidence['duplicate_primary_rows'], 1)
        self.assertFalse(evidence['competing_primary_mapping'])
        self.assertFalse(evidence['localization_confirmation'])

    def test_unpositioned_alignment_duplicates_retain_ambiguity(self):
        case = next(x for x in control_fixtures() if x['id'] == 'perfect_six')
        unpositioned = dict(query_start0=0, query_end0=150, target='chr1', strand='+', mapq=60, tp='P')
        evidence = M.evaluate_pair(case['sequence'], case['left'], case['right'], case['motif'],
                                   alignments=[unpositioned, dict(unpositioned)])['alignment_evidence']
        self.assertEqual(evidence['distinct_primary_count'], 2)
        self.assertEqual(evidence['duplicate_primary_rows'], 0)
        self.assertTrue(evidence['competing_primary_mapping'])
        self.assertFalse(evidence['localization_confirmation'])

    def test_partial_or_invalid_reference_coordinates_raise(self):
        case = next(x for x in control_fixtures() if x['id'] == 'perfect_six')
        base = dict(query_start0=0, query_end0=150, target='chr1', strand='+', mapq=60, tp='P')
        for coordinates in [dict(target_start0=100),
                            dict(target_start0=100, target_end0=250),
                            dict(target_start0=250, target_end0=100, target_length=1000),
                            dict(target_start0=100, target_end0=1001, target_length=1000),
                            dict(target_start0=True, target_end0=250, target_length=1000)]:
            with self.subTest(coordinates=coordinates), self.assertRaises(ValueError):
                M.evaluate_pair(case['sequence'], case['left'], case['right'], case['motif'],
                                alignments=[dict(base, **coordinates)])


class IndependentExactCoverage(unittest.TestCase):
    def test_base_index_oracle_with_variants_and_mutations(self):
        rng = random.Random(92728)
        motifs = ['TTAGGG', 'TTTAGGG', 'TTAGGC', 'TTAGAC', 'TTAGGT', 'TTAGG', 'CTCGGTTATGGG', 'ACGT']
        for motif in motifs:
            k = len(motif)
            for variant in range(9):
                source = motif if variant % 2 == 0 else rc(motif)
                seq = list('NNN' + source * 8 + 'NN')
                for _ in range(variant):
                    seq[rng.randrange(3, len(seq)-2)] = rng.choice('ACGTN')
                seq = ''.join(seq)
                start, end = (0, len(seq)) if variant % 3 else (5, len(seq)-4)
                cuts = [k, 2*k, 4*k]
                result = M.analyze_interval(seq, start, end, motif, cutoffs=cuts)
                for cutoff in cuts:
                    view = result['exact_coverage_views'][str(cutoff)]
                    for strand in [None, '+', '-']:
                        observed = view['all_strands'] if strand is None else view['by_strand'][strand]
                        expected = coverage_oracle(seq, start, end, motif, cutoff, strand)
                        with self.subTest(motif=motif, variant=variant, cutoff=cutoff, strand=strand):
                            self.assertEqual(intervals_as_set(observed['coverage_union_intervals']), expected)
                            self.assertEqual(observed['covered_bp'], len(expected))
                            self.assertEqual(observed['uncovered_bp'], end-start-len(expected))
                assert_no_promotion(self, result)

    def test_reverse_complement_exact_coverage(self):
        sequence = 'NN' + 'TTAGGG'*8 + 'CC' + 'CCCTAA'*7 + 'N'
        n = len(sequence)
        forward = M.analyze_interval(sequence, 3, n-3, 'TTAGGG')
        reverse = M.analyze_interval(rc(sequence), 3, n-3, 'TTAGGG')
        for cutoff in ['6', '12', '24']:
            for strand, opposite in [('+', '-'), ('-', '+')]:
                fwd = intervals_as_set(forward['exact_coverage_views'][cutoff]['by_strand'][strand]['coverage_union_intervals'])
                rev = intervals_as_set(reverse['exact_coverage_views'][cutoff]['by_strand'][opposite]['coverage_union_intervals'])
                self.assertEqual({n-1-i for i in fwd}, rev)

    def test_origin_and_crop_phase_invariance(self):
        seq = 'NNN' + 'TTTAGGG'*20 + 'NN'
        full = M.analyze_interval(seq, 3, 143, 'TTTAGGG')
        crop = M.analyze_interval(seq, 8, 130, 'TTTAGGG')
        phase = full['fixed_phase']['all_best_ties'][0]['absolute_template_phase']
        self.assertEqual(crop['fixed_phase']['all_best_ties'][0]['absolute_template_phase'], phase)
        shifted = M.analyze_interval('NN'+seq, 5, 145, 'TTTAGGG')
        self.assertEqual(shifted['fixed_phase']['all_best_ties'][0]['absolute_template_phase'], (phase-2)%7)
        self.assertTrue(crop['exact_segments_ge_period'][0]['left_interval_censored'])
        self.assertTrue(crop['exact_segments_ge_period'][0]['exact_continuation_left_observed'])

    def test_overlapping_known_tracts_are_counted_by_base_union(self):
        tracts = [dict(tract_id='a', start0=0, end0=40, strand='+'),
                  dict(tract_id='b', start0=30, end0=60, strand='+')]
        result = M.analyze_interval('TTAGGG'*10, 0, 60, 'TTAGGG', known_tracts=tracts)
        self.assertEqual(result['known_tract_coverage']['covered_bp'], 60)
        self.assertEqual(result['known_tract_coverage']['raw_intersection_length_sum_bp'], 70)
        self.assertEqual(result['known_tract_coverage']['double_counted_overlap_bp'], 10)
        self.assertEqual(len(result['known_tract_intersections']), 2)
        assert_no_promotion(self, result)

    def test_short_crop_empty_and_unknown_intervals(self):
        short = M.analyze_interval('TTAGGG'*5, 3, 8, 'TTAGGG')
        self.assertEqual(short['exact_segments_ge_period'], [])
        self.assertEqual(short['exact_coverage_views']['6']['all_strands']['covered_bp'], 0)
        empty = M.analyze_interval('TTAGGG', 3, 3, 'TTAGGG', quality='I'*6)
        self.assertTrue(empty['fixed_phase']['empty_interval'])
        self.assertIsNone(empty['exact_coverage_views']['6']['all_strands']['covered_fraction'])
        self.assertEqual(empty['quality']['bases'], 0)
        unknown = M.analyze_interval('N'*30, 0, 30, 'TTAGGG')
        self.assertEqual(len(unknown['fixed_phase']['all_best_ties']), 12)
        self.assertEqual(unknown['exact_coverage_views']['6']['all_strands']['covered_bp'], 0)

    def test_palindromic_strands_remain_ambiguous(self):
        result = M.analyze_interval('ACGT'*20, 0, 80, 'ACGT')
        self.assertEqual({x['strand'] for x in result['fixed_phase']['all_best_ties']}, {'+', '-'})

    def test_quality_and_lowercase_are_preserved(self):
        sequence = 'ttagggnttaggg'
        result = M.analyze_interval(sequence, 0, 13, 'ttaggg', quality='I'*6+'!'+'I'*6)
        islands = result['exact_coverage_views']['6']['all_strands']['uncovered_islands']
        self.assertEqual([(x['start0'], x['end0']) for x in islands], [(6, 7)])
        self.assertEqual(islands[0]['quality']['mean'], 0)
        self.assertEqual(result['quality']['q_below_20'], 1)

    def test_invalid_inputs_raise(self):
        for start, end in [(-1, 6), (0, 7), (5, 4), (0.0, 6), (False, 6)]:
            with self.subTest(start=start, end=end), self.assertRaises(ValueError):
                M.analyze_interval('TTAGGG', start, end, 'TTAGGG')
        for motif in ['', 'NNN', 'TTAGGGTTAGGG', 'A'*65]:
            with self.subTest(motif=motif), self.assertRaises(ValueError):
                M.analyze_interval('TTAGGG', 0, 6, motif)
        for kwargs in [dict(quality='I'), dict(quality='\x7f'*6), dict(expected_strand='?'),
                       dict(cutoffs=[0]), dict(cutoffs=[5]), dict(cutoffs=[6.0])]:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                M.analyze_interval('TTAGGG', 0, 6, 'TTAGGG', **kwargs)

    def test_invalid_nonintersecting_tract_identity_is_not_silently_ignored(self):
        for value in [float('nan'), float('inf'), -0.1, 1.1, True]:
            outside = dict(tract_id='outside', start0=90, end0=120, strand='+', identity=value)
            with self.subTest(identity=value), self.assertRaises(ValueError):
                M.analyze_interval('TTAGGG'*20, 0, 60, 'TTAGGG', known_tracts=[outside])


def write_outcomes(destination):
    observations = []
    for case in control_fixtures():
        result = evaluate_fixture(case)
        observations.append(dict(control_id=case['id'], purpose=case['purpose'],
                                 expected_status=case['expected_status'], observed_status=result['status'],
                                 status_matches=result['status'] == case['expected_status'],
                                 required_reason_codes=case.get('required_reason_codes', []),
                                 observed_reason_codes=result['reason_codes'],
                                 terminal_context=result['terminal_context'],
                                 alignment_evidence=result['alignment_evidence'], result=result))
    payload = dict(protocol_sha256=hashlib.sha256((CONTROL_DIR/'protocol.json').read_bytes()).hexdigest(),
                   module_sha256=hashlib.sha256((ROOT/'scripts/tract_context.py').read_bytes()).hexdigest(),
                   scope='Synthetic development controls; not an independent biological holdout or accuracy estimate',
                   controls=observations)
    Path(destination).write_text(json.dumps(payload, indent=2)+'\n')


if __name__ == '__main__':
    destination = None
    if '--write-outcomes' in sys.argv:
        position = sys.argv.index('--write-outcomes')
        destination = sys.argv[position+1]
        del sys.argv[position:position+2]
    program = unittest.main(exit=False)
    if destination:
        write_outcomes(destination)
    sys.exit(not program.result.wasSuccessful())
