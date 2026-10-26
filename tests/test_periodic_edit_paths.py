"""Global-path constraints, quality provenance and bounded optional evaluation."""
from copy import deepcopy
from pathlib import Path
import random
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import periodic_edit_paths as paths
from repeat_region_continuity import minimum_periodic_edit_score


class PeriodicEditPathsTests(unittest.TestCase):
    def test_exact_complete_query_and_reverse_strand(self):
        for sequence, strand in [('TTAGGGTTAGGG', '+'), ('CCCTAACCCTAA', '-')]:
            result = paths.analyze_periodic_paths(sequence, 'TTAGGG')
            self.assertEqual(result['minimum_edit_cost'], 0)
            self.assertEqual(result['all_optimal_strands'], [strand])
            self.assertEqual(result['operation_count_ranges']['match'], {'min': 12, 'max': 12})
            self.assertTrue(all(p['operation_set'] == ['match'] for p in result['query_positions']))
            self.assertTrue(all(g['max_deletions'] == 0 for g in result['deletion_boundaries']))

    def test_mandatory_substitution_keeps_actual_coordinate_and_quality(self):
        sequence = 'ACGACGATGACG'
        quality = 'I' * 7 + '!' + 'I' * 4
        result = paths.analyze_periodic_paths(sequence, 'ACG', quality, start0=10)
        changed = [r for r in result['query_positions'] if r['mandatory_edit']]
        self.assertEqual(len(changed), 1)
        self.assertEqual(changed[0]['relative_position0'], 7)
        self.assertEqual(changed[0]['absolute_position0'], 17)
        self.assertEqual(changed[0]['query_base'], 'T')
        self.assertEqual(changed[0]['template_base_set'], ['C'])
        self.assertEqual(changed[0]['operation_set'], ['substitution'])
        self.assertEqual(changed[0]['quality_character'], '!')
        self.assertEqual(changed[0]['quality_ascii'], 33)
        self.assertEqual(changed[0]['conditional_phred33'], 0)
        self.assertFalse(result['quality_encoding_verified'])
        self.assertEqual(result['operation_count_ranges']['substitution'], {'min': 1, 'max': 1})

    def test_deletion_has_no_observed_base_quality(self):
        result = paths.analyze_periodic_paths('ACGACGAGACG', 'ACG', 'I' * 11, start0=20)
        self.assertEqual(result['minimum_edit_cost'], 1)
        self.assertEqual(result['operation_count_ranges']['deletion'], {'min': 1, 'max': 1})
        gaps = [g for g in result['deletion_boundaries'] if g['max_deletions']]
        self.assertEqual(len(gaps), 1)
        self.assertEqual(gaps[0], {'relative_boundary0': 7, 'absolute_boundary0': 27,
                                  'template_base_set': ['C'], 'min_deletions': 1,
                                  'max_deletions': 1, 'quality': None})
        self.assertTrue(all(not q['possible_edit'] for q in result['query_positions']))
        self.assertTrue(all(g['quality'] is None for g in result['deletion_boundaries']))

    def test_insertion_location_can_be_ambiguous_without_mandatory_edited_base(self):
        result = paths.analyze_periodic_paths('ACGACGAACGACG', 'ACG', 'I' * 13)
        self.assertEqual(result['operation_count_ranges']['insertion'], {'min': 1, 'max': 1})
        possible = [q['relative_position0'] for q in result['query_positions'] if q['possible_edit']]
        self.assertEqual(possible, [6, 7])
        self.assertFalse(any(q['mandatory_edit'] for q in result['query_positions']))
        self.assertTrue(all(result['query_positions'][i]['operation_set'] == ['match', 'insertion'] for i in possible))

    def test_distinct_optimal_operation_types_not_one_chosen_traceback(self):
        result = paths.analyze_periodic_paths('AAC', 'AC')
        self.assertEqual(result['minimum_edit_cost'], 1)
        for op in ('substitution', 'insertion', 'deletion'):
            self.assertEqual(result['operation_count_ranges'][op], {'min': 0, 'max': 1})
        self.assertEqual(result['query_positions'][0]['operation_set'], ['match', 'substitution', 'insertion'])
        self.assertEqual(result['deletion_boundaries'][1]['template_base_set'], ['C'])
        self.assertEqual(result['deletion_boundaries'][1]['min_deletions'], 0)
        self.assertEqual(result['deletion_boundaries'][1]['max_deletions'], 1)

    def test_all_globally_tied_strands_and_template_bases(self):
        result = paths.analyze_periodic_paths('A', 'C')
        self.assertEqual(result['all_optimal_strands'], ['+', '-'])
        self.assertEqual(result['query_positions'][0]['template_base_set'], ['C', 'G'])
        self.assertEqual(result['query_positions'][0]['operation_set'], ['substitution', 'insertion'])
        self.assertTrue(result['query_positions'][0]['mandatory_edit'])
        self.assertEqual(result['operation_count_ranges']['insertion'], {'min': 0, 'max': 1})

    def test_primitive_reduction_and_input_phase_frame(self):
        first = paths.analyze_periodic_paths('CAGCAG', 'CAGCAG')
        second = paths.analyze_periodic_paths('CAGCAG', 'CAG')
        self.assertEqual(first['primitive_motif'], 'CAG')
        self.assertEqual(first['primitive_period_bp'], 3)
        self.assertEqual(first['query_positions'], second['query_positions'])
        self.assertEqual(first['orientation_results'], second['orientation_results'])
        self.assertEqual(first['orientation_results'][0]['start_phases'], [0])

    def test_empty_query_has_zero_cost_free_phase_and_no_quality_bases(self):
        result = paths.analyze_periodic_paths('', 'ACG', '', start0=50)
        self.assertTrue(result['evaluated'])
        self.assertEqual(result['minimum_edit_cost'], 0)
        self.assertIsNone(result['edit_cost_per_query_base'])
        self.assertEqual(result['query_positions'], [])
        self.assertEqual(result['all_optimal_strands'], ['+', '-'])
        for orientation in result['orientation_results']:
            self.assertEqual(orientation['start_phases'], [0, 1, 2])
            self.assertEqual(orientation['end_phases'], [0, 1, 2])
            self.assertEqual(orientation['optimal_edges'], [])
        self.assertEqual(result['deletion_boundaries'], [dict(relative_boundary0=0,
            absolute_boundary0=50, template_base_set=[], min_deletions=0, max_deletions=0, quality=None)])

    def test_quality_values_cannot_select_a_path(self):
        low = paths.analyze_periodic_paths('AAC', 'AC', '!!!')
        high = paths.analyze_periodic_paths('AAC', 'AC', 'III')
        absent = paths.analyze_periodic_paths('AAC', 'AC')
        for result in (high, absent):
            self.assertEqual(low['orientation_results'], result['orientation_results'])
            self.assertEqual(low['operation_count_ranges'], result['operation_count_ranges'])
        self.assertTrue(all(q['quality_character'] is None for q in absent['query_positions']))
        self.assertTrue(all(q['quality_ascii'] is None for q in absent['query_positions']))

    def test_malformed_arguments_and_quality_fail(self):
        cases = [(None, 'AC', None, 0), ('AC', None, None, 0), ('AC', '', None, 0),
                 ('AC', 'AC', b'II', 0), ('AC', 'AC', 'I', 0), ('AC', 'AC', ' I', 0),
                 ('AC', 'AC', 'éI', 0), ('AC', 'AC', 'II', -1), ('AC', 'AC', 'II', True)]
        for args in cases:
            with self.subTest(args=args), self.assertRaises(ValueError):
                paths.analyze_periodic_paths(*args)

    def test_ambiguous_and_oversized_inputs_are_not_evaluated_not_zero_cost(self):
        cases = [('ANC', 'AC'), ('ac', 'AC'), ('AC', 'AN'),
                 ('A' * (paths.MAX_REGION_BP + 1), 'A'), ('AC', 'A' * (paths.MAX_MOTIF_BP + 1))]
        for sequence, motif in cases:
            with self.subTest(n=len(sequence), motif=motif):
                result = paths.analyze_periodic_paths(sequence, motif)
                self.assertEqual(result['status'], 'not_evaluated')
                self.assertFalse(result['evaluated'])
                self.assertTrue(result['reason_codes'])
                self.assertIsNone(result['minimum_edit_cost'])
                self.assertIsNone(result['operation_count_ranges'])
                self.assertEqual(result['query_positions'], [])
                self.assertEqual(result['orientation_results'], [])
        with patch.object(paths, 'MAX_TOTAL_STATES', 10):
            self.assertIn('query_motif_graph_exceeds_state_bound',
                paths.analyze_periodic_paths('ACG', 'ACG')['reason_codes'])

    def test_costs_match_existing_score_only_helper_and_graph_coordinates(self):
        generator = random.Random(681)
        for _ in range(250):
            sequence = ''.join(generator.choice('ACGT') for _ in range(generator.randrange(36)))
            motif = ''.join(generator.choice('ACGT') for _ in range(generator.randrange(1, 9)))
            result = paths.analyze_periodic_paths(sequence, motif, start0=3)
            unit = result['primitive_motif']
            reverse = unit.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
            costs = [minimum_periodic_edit_score(sequence, template)['edit_cost'] for template in (unit, reverse)]
            self.assertEqual(result['minimum_edit_cost'], min(costs))
            for orientation, expected_cost in zip(result['orientation_results'], costs):
                self.assertEqual(orientation['minimum_edit_cost'], expected_cost)
                if not orientation['globally_optimal']:
                    self.assertEqual(orientation['optimal_edges'], [])
                    self.assertIsNone(orientation['operation_count_ranges'])
                for edge in orientation['optimal_edges']:
                    i, j = edge['from_query_boundary0'], edge['to_query_boundary0']
                    self.assertTrue(0 <= i <= j <= len(sequence))
                    self.assertTrue(0 <= edge['from_phase0'] < len(unit))
                    if edge['operation'] == 'deletion':
                        self.assertEqual(i, j)
                        self.assertIsNone(edge['query_base'])
                    else:
                        self.assertEqual(j, i + 1)
                        self.assertEqual(edge['query_base'], sequence[i])
            self.assertEqual(len(result['query_positions']), len(sequence))
            self.assertEqual(len(result['deletion_boundaries']), len(sequence) + 1)
            self.assertFalse(result['native_counts_changed'])
            self.assertFalse(result['classification_changed'])
            self.assertIsNone(result['biological_confirmation'])


if __name__ == '__main__':
    unittest.main()
