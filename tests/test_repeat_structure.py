"""Integer graph oracles and unchanged actual-source slices; no DNA fixtures."""
from collections import defaultdict
import hashlib
import itertools
import json
from pathlib import Path
import sys
import unittest

from archive_test_support import archive_only

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import diagnose_terminal_architecture as OLD
import repeat_structure as M


def exhaustive_paths(n, edges):
    """Enumerate every numeric edge chain; independent of endpoint dominance DP."""
    outgoing = defaultdict(list)
    for edge in edges:
        outgoing[edge[0]].append(edge)
    observations = defaultdict(list)

    def extend(origin, end, cost, minimum_ways, all_ways):
        observations[end].append((origin, cost, minimum_ways, all_ways))
        for a, b, c, minimum, total in outgoing[end]:
            extend(origin, b, cost + c, minimum_ways * minimum, all_ways * total)

    for a, b, cost, minimum, total in edges:
        extend(a, b, cost, minimum, total)
    endpoints = []
    for end in range(n + 1):
        if end not in observations:
            continue
        earliest = min(row[0] for row in observations[end])
        paths = [row for row in observations[end] if row[0] == earliest]
        minimum = min(row[1] for row in paths)
        endpoints.append({'start0': earliest, 'end0': end, 'minimum_edit_count': minimum,
                          'optimal_path_count_capped2': min(2, sum(row[2] for row in paths if row[1] == minimum)),
                          'possible_path_count_capped2': min(2, sum(row[3] for row in paths))})
    # Set-theoretic containment check, independent of the production sweep.
    maximal = [row for row in endpoints if not any(
        other['start0'] <= row['start0'] and other['end0'] >= row['end0']
        and (other['start0'], other['end0']) != (row['start0'], row['end0'])
        for other in endpoints)]
    maximal.sort(key=lambda row: (row['start0'], -row['end0'], row['minimum_edit_count']))
    return maximal, endpoints


class NumericPathEvidenceTests(unittest.TestCase):
    def test_all_1024_coordinate_graphs_with_three_cost_multiplicity_patterns(self):
        coordinates = [(a, b) for a in range(4) for b in range(a + 1, 5)]
        checked = 0
        for pattern in range(3):
            universe = [(a, b, (a + 2 * b + pattern) % 3,
                         1 + int(pattern == 1 and (a + b) % 2 == 0),
                         3 if pattern == 2 else 1 + int(pattern == 1 and (a + b) % 2 == 0))
                        for a, b in coordinates]
            for mask in range(1 << len(universe)):
                edges = [edge for i, edge in enumerate(universe) if mask & (1 << i)]
                actual = M.solve_evidence(4, edges)
                expected = exhaustive_paths(4, edges)
                self.assertEqual(actual, expected, (pattern, mask))
                old_maximal, old_endpoints = OLD.paths(4, [edge[:3] for edge in edges])
                self.assertEqual([(r['start0'], r['end0'], r['minimum_edit_count'])
                                  for r in actual[0]], old_maximal)
                self.assertEqual([(r['start0'], r['end0'], r['minimum_edit_count'])
                                  for r in actual[1]], old_endpoints)
                checked += 1
        self.assertEqual(checked, 3072)

    def test_higher_cost_alternative_is_possible_but_not_optimal(self):
        edges = [(0, 2, 0, 1, 2), (2, 4, 0, 1, 1), (0, 4, 1, 1, 1)]
        maximal, endpoints = M.solve_evidence(4, edges)
        self.assertEqual((maximal, endpoints), exhaustive_paths(4, edges))
        self.assertEqual(maximal[0]['minimum_edit_count'], 0)
        self.assertEqual(maximal[0]['optimal_path_count_capped2'], 1)
        self.assertEqual(maximal[0]['possible_path_count_capped2'], 2)

    def test_equal_cost_segmentations_and_local_assignments_are_distinct(self):
        for edges in [[(0, 2, 0, 1, 1), (2, 4, 0, 1, 1), (0, 4, 0, 1, 1)],
                      [(0, 2, 1, 2, 3), (2, 4, 0, 1, 1)]]:
            maximal, endpoints = M.solve_evidence(4, edges)
            self.assertEqual((maximal, endpoints), exhaustive_paths(4, edges))
            self.assertEqual(maximal[0]['optimal_path_count_capped2'], 2)
            self.assertEqual(maximal[0]['possible_path_count_capped2'], 2)

    def test_earliest_start_dominates_cheaper_later_start_without_free_gap(self):
        edges = [(0, 2, 5, 1, 1), (1, 2, 0, 1, 1), (2, 4, 0, 1, 1),
                 (5, 7, 0, 1, 1)]
        maximal, endpoints = M.solve_evidence(7, edges)
        self.assertEqual((maximal, endpoints), exhaustive_paths(7, edges))
        self.assertEqual([(r['start0'], r['end0']) for r in maximal], [(0, 4), (5, 7)])
        self.assertEqual(maximal[0]['minimum_edit_count'], 5)
        self.assertEqual(maximal[0]['possible_path_count_capped2'], 1)

    def test_maximal_intervals_may_cross_and_are_not_a_tiling(self):
        maximal, _ = M.solve_evidence(5, [(0, 3, 0, 1, 1), (2, 5, 0, 1, 1)])
        self.assertEqual([(r['start0'], r['end0']) for r in maximal], [(0, 3), (2, 5)])
        self.assertGreater(sum(r['end0'] - r['start0'] for r in maximal), 5)

    def test_empty_graph_and_edge_order(self):
        self.assertEqual(M.solve_evidence(6, []), ([], []))
        edges = [(0, 2, 1, 1, 2), (2, 4, 0, 2, 2), (4, 6, 1, 1, 1)]
        expected = M.solve_evidence(6, edges)
        for permutation in itertools.permutations(edges):
            self.assertEqual(M.solve_evidence(6, list(permutation)), expected)

    def test_invalid_numeric_graph_cannot_receive_evidence(self):
        bad_edges = [[(0, 0, 0, 1, 1)], [(2, 1, 0, 1, 1)], [(0, 5, 0, 1, 1)],
                     [(-1, 1, 0, 1, 1)], [(0, 1, -1, 1, 1)], [(0, 1, 0, 0, 1)],
                     [(0, 1, 0, 2, 1)], [(0, 1, 0, 1, 0)], [(0, 1, 0, True, 1)],
                     [(0, 1.0, 0, 1, 1)], [(False, 1, 0, 1, 1)],
                     [(0, 1, 0, 1, 1), (0, 1, 0, 1, 1)]]
        for edges in bad_edges:
            with self.subTest(edges=edges), self.assertRaises(ValueError):
                M.solve_evidence(4, edges)
        for n in [-1, True, 1.5]:
            with self.subTest(length=n), self.assertRaises(ValueError):
                M.solve_evidence(n, [])


class NumericCoverageCoordinateTests(unittest.TestCase):
    def test_all_243_base_classifications_form_exact_nonoverlap_partition(self):
        names = ['uncovered', 'edited_only', 'literal']
        for assignment in itertools.product(range(3), repeat=5):
            literal = [(i, i + 1) for i, kind in enumerate(assignment) if kind == 2]
            tokens = [(i, i + 1) for i, kind in enumerate(assignment) if kind]
            expected = []
            for i, kind in enumerate(assignment):
                if expected and expected[-1]['class'] == names[kind]:
                    expected[-1]['end0'] = i + 1
                else:
                    expected.append({'start0': i, 'end0': i + 1, 'class': names[kind]})
            self.assertEqual(M.coverage_segments(5, literal, tokens), expected)
            self.assertEqual(M.coverage_segments(5, list(reversed(literal)), tokens + tokens), expected)
        self.assertEqual(M.coverage_segments(0, [], []), [])

    def test_coverage_overlap_and_invalid_intervals(self):
        self.assertEqual(M.coverage_segments(7, [(1, 3), (2, 5)], [(0, 4), (3, 6)]),
                         [{'start0': 0, 'end0': 1, 'class': 'edited_only'},
                          {'start0': 1, 'end0': 5, 'class': 'literal'},
                          {'start0': 5, 'end0': 6, 'class': 'edited_only'},
                          {'start0': 6, 'end0': 7, 'class': 'uncovered'}])
        for literal, tokens in [([(0, 2)], [(0, 1)]), ([], [(1, 1)]),
                                ([], [(-1, 2)]), ([], [(0, 8)]), ([], [(False, 1)])]:
            with self.subTest(literal=literal, tokens=tokens), self.assertRaises(ValueError):
                M.coverage_segments(7, literal, tokens)

    def test_original_coordinates_and_reverse_inverse_for_every_numeric_subinterval(self):
        for offset in [0, 17, 1000]:
            for a in range(10):
                for b in range(a + 1, 11):
                    for strand in ['+', '-']:
                        x, y = M.window_coordinates(a, b, offset, offset + 10, strand)
                        self.assertEqual(y - x, b - a)
                        self.assertTrue(offset <= x < y <= offset + 10)
                        inverse = (x - offset, y - offset) if strand == '+' else (offset + 10 - y, offset + 10 - x)
                        self.assertEqual(inverse, (a, b))
        for args in [(0, 0, 0, 10, '+'), (0, 11, 0, 10, '-'),
                     (0, 1, 10, 0, '+'), (0, 1, 0, 10, '?'), (False, 1, 0, 10, '+')]:
            with self.subTest(args=args), self.assertRaises(ValueError):
                M.window_coordinates(*args)


def source_explanations(observed, units):
    """Independent direct alignment checks on unchanged observed source slices."""
    result = []
    if any(c not in b'ACGT' for c in observed):
        return result
    for index, unit in enumerate(units):
        q, k = len(observed), len(unit)
        if observed == unit:
            result.append((0, index, '=', 0))
        if k < 5:
            continue
        if q == k:
            mismatches = [i for i in range(k) if observed[i] != unit[i]]
            if len(mismatches) == 1:
                result.append((1, index, 'S', mismatches[0]))
        elif q == k + 1:
            for offset in range(q):
                if all(observed[i] == unit[i - int(i > offset)] for i in range(q) if i != offset):
                    result.append((1, index, 'I', offset))
        elif q + 1 == k:
            for offset in range(k):
                if all(observed[i - int(i > offset)] == unit[i] for i in range(k) if i != offset):
                    result.append((1, index, 'D', offset))
    return result


@archive_only
class ExistingSourceStructureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = ROOT / 'results/terminal_architecture_review_20260927/diagnosis/windows.json'
        if not path.exists():
            raise unittest.SkipTest('unchanged actual source windows unavailable')
        content = path.read_bytes()
        if hashlib.sha256(content).hexdigest() != 'f46da90b852f18b2deb73e4928a101566ef79f6b53e3ff8f40ba1eb0eaffaa24':
            raise AssertionError('frozen actual-source windows changed')
        cls.primary = [b'TTAGGG', b'TTAGTC']
        cls.hypotheses = [b'TTTAGGG']
        cls.observations = []
        for row in json.loads(content):
            if row['role'] != 'deposited_terminal_window' or row['window_bp'] != 500:
                continue
            sequence = row['sequence_outward_orientation'].encode('ascii')
            if hashlib.sha256(sequence).hexdigest() != row['oriented_slice_sha256']:
                raise AssertionError('actual oriented source slice mismatch')
            window = {'record_length': row['record_length'], 'start0': row['window_start0'],
                      'end0': row['window_end0'], 'strand': row['strand']}
            result = M.explain_window(sequence, cls.primary, cls.hypotheses, window)
            cls.observations.append((row, sequence, window, result))
        if len(cls.observations) != 6:
            raise AssertionError('require all six actual deposited 500bp end windows')

    def test_all_actual_token_explanations_and_ambiguity_flags_replay(self):
        for frozen, sequence, window, result in self.observations:
            expected = {}
            for a in range(len(sequence)):
                for size in [5, 6, 7]:
                    b = a + size
                    if b <= len(sequence):
                        choices = source_explanations(sequence[a:b], self.primary)
                        if choices:
                            expected[(a, b)] = choices
            self.assertEqual(len(result['tokens']), len(expected))
            for token in result['tokens']:
                choices = expected[(token['start0'], token['end0'])]
                actual = [(c['edit_count'], c['unit_index'], c['edit'], c['edit_offset'])
                          for c in token['explanations']]
                self.assertEqual(actual, choices)
                minimum = min(x[0] for x in choices)
                self.assertEqual(token['minimum_edit_count'], minimum)
                self.assertEqual(token['minimum_cost_variant_count'], sum(x[0] == minimum for x in choices))
                self.assertEqual(token['all_variant_count'], len(choices))
                self.assertEqual(token['unit_assignment_ambiguous'], len({x[1] for x in choices}) > 1)
                self.assertEqual(token['unit_or_operation_ambiguous'], len({(x[1], x[2]) for x in choices}) > 1)
                groups = defaultdict(list)
                for _, unit, op, offset in choices:
                    groups[(unit, op)].append(offset)
                self.assertEqual(token['edit_position_ambiguous'], any(len(x) > 1 for x in groups.values()))
            self.assertEqual(result['summary']['literal_union_bp'], frozen['literal_union_bp'])
            self.assertEqual(result['summary']['allowed_token_union_bp'], frozen['allowed_token_union_bp'])
            self.assertEqual(result['summary']['max_complete_path_bp'], frozen['max_edit_model_path_bp_before_filters'])

    def test_explanatory_seven_unit_preserves_primary_graph_and_all_insertion_offsets(self):
        observations = 0
        for _, sequence, window, result in self.observations:
            without = M.explain_window(sequence, self.primary, [], window)
            for key in ['tokens', 'coverage_segments', 'maximal_complete_paths']:
                self.assertEqual(result[key], without[key])
            expected_starts = [i for i in range(len(sequence) - 6) if sequence[i:i + 7] == self.hypotheses[0]]
            self.assertEqual([x['start0'] for x in result['literal_hypotheses']], expected_starts)
            for overlay in result['literal_hypotheses']:
                self.assertEqual(overlay['hypothesis_role'], 'separate_literal_not_in_primary_graph')
                self.assertTrue(overlay['local_unit_length_nonidentifiable'])
                edits = [(x['unit_index'], x['edit'], x['edit_offset'])
                         for x in overlay['primary_same_interval_explanations']]
                self.assertEqual(edits, [(0, 'I', 0), (0, 'I', 1), (0, 'I', 2)])
                observations += 1
        self.assertGreater(observations, 0)

    def test_actual_coverage_partition_coordinates_and_path_stops(self):
        for _, sequence, window, result in self.observations:
            n = len(sequence)
            expected = ['uncovered'] * n
            for token in result['tokens']:
                for i in range(token['start0'], token['end0']):
                    if token['minimum_edit_count'] == 0:
                        expected[i] = 'literal'
                    elif expected[i] == 'uncovered':
                        expected[i] = 'edited_only'
            actual = []
            for segment in result['coverage_segments']:
                self.assertEqual(segment['start0'], len(actual))
                actual.extend([segment['class']] * segment['span_bp'])
            self.assertEqual(actual, expected)
            for key in ['tokens', 'coverage_segments', 'maximal_complete_paths', 'literal_hypotheses', 'boundary_fragments']:
                for row in result[key]:
                    a, b = row['start0'], row['end0']
                    want = (window['end0'] - b, window['end0'] - a) if window['strand'] == '-' else (window['start0'] + a, window['start0'] + b)
                    self.assertEqual((row['original_start0'], row['original_end0']), want)
            starts = {t['start0'] for t in result['tokens']}
            for row in result['maximal_complete_paths']:
                end = row['end0']
                self.assertNotIn(end, starts)
                self.assertEqual(row['successor_status'], 'window_edge' if end == n else 'no_token_starts_at_exact_boundary')
                self.assertEqual(row['next_token_start0'], min((s for s in starts if s > end), default=None))
                self.assertEqual(row['next_base_has_token_coverage'],
                                 any(t['start0'] <= end < t['end0'] for t in result['tokens']) if end < n else None)

    def test_actual_fragments_are_proper_compatible_and_edge_scope_is_explicit(self):
        anchored_original = artificial = complete_coexistence = 0
        for _, sequence, window, result in self.observations:
            expected = set()
            for vocabulary, units in [('primary', self.primary), ('explanatory_literal', self.hypotheses)]:
                for index, unit in enumerate(units):
                    for size in range(1, len(unit)):
                        if sequence[:size] == unit[-size:]:
                            expected.add((vocabulary, index, 'left', 0, size))
                        if sequence[-size:] == unit[:size]:
                            expected.add((vocabulary, index, 'right', len(sequence) - size, len(sequence)))
            actual = {(r['vocabulary'], r['unit_index'], r['oriented_edge'], r['start0'], r['end0']) for r in result['boundary_fragments']}
            self.assertEqual(actual, expected)
            for row in result['boundary_fragments']:
                a, b = row['start0'], row['end0']
                self.assertEqual(row['observed_source_slice'].encode(), sequence[a:b])
                self.assertTrue(0 < b - a < len(row['unit']))
                self.assertEqual(row['boundary_kind'], 'artificial_crop' if row['oriented_edge'] == 'left' else 'original_record_edge')
                self.assertFalse(row['extends_detector_array'])
                self.assertFalse(row['physical_telomere_confirmed'])
                for anchor in row['adjacent_literal_anchors']:
                    self.assertEqual(sequence[anchor['start0']:anchor['end0']], anchor['literal'].encode())
                    self.assertTrue(anchor['end0'] == a or anchor['start0'] == b)
                self.assertEqual(row['anchored_original_record_fragment'], row['boundary_kind'] == 'original_record_edge' and bool(row['adjacent_literal_anchors']))
                choices = source_explanations(sequence[a:b], self.primary)
                self.assertEqual([(x['edit_count'], x['unit_index'], x['edit'], x['edit_offset']) for x in row['complete_primary_token_explanations_also_possible']], choices)
                anchored_original += row['anchored_original_record_fragment']
                artificial += row['boundary_kind'] == 'artificial_crop'
                complete_coexistence += bool(choices)
        self.assertGreater(anchored_original, 0)
        self.assertGreater(artificial, 0)
        self.assertGreater(complete_coexistence, 0)

    def test_actual_source_limits_fail_instead_of_truncating_tokens(self):
        _, sequence, window, result = max(self.observations, key=lambda x: len(x[3]['tokens']))
        for limits in [{'max_edges': 1}, {'max_explanations': 1}]:
            with self.subTest(limits=limits), self.assertRaisesRegex(ValueError, 'limit exceeded'):
                M.explain_window(sequence, self.primary, self.hypotheses, window, **limits)

    def test_invalid_api_metadata_or_units_fail_on_actual_source(self):
        _, sequence, window, _ = self.observations[0]
        for key in ['min_span', 'min_exact_anchors', 'max_edges', 'max_explanations']:
            for value in [0, -1, True, 1.5]:
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    M.explain_window(sequence, self.primary, self.hypotheses, window, **{key: value})
        for units in [[], [self.primary[0], self.primary[0]], [123], None]:
            with self.subTest(units=units), self.assertRaises(ValueError):
                M.explain_window(sequence, units, self.hypotheses, window)
        for mutation in [{'strand': '?'}, {'start0': True}, {'end0': window['end0'] + 1},
                         {'record_length': 0}, {'start0': -1}]:
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                M.explain_window(sequence, self.primary, self.hypotheses, {**window, **mutation})

    def test_original_window_crop_censoring_flags_do_not_change_paths(self):
        for _, sequence, window, result in self.observations:
            left_crop = (window['start0'] != 0) if window['strand'] == '+' else (window['end0'] != window['record_length'])
            right_crop = (window['end0'] != window['record_length']) if window['strand'] == '+' else (window['start0'] != 0)
            self.assertTrue(left_crop)
            self.assertFalse(right_crop)
            for path in result['maximal_complete_paths']:
                self.assertEqual(path['left_crop_context_may_be_missing'], path['start0'] < 7)
                self.assertFalse(path['right_crop_context_may_be_missing'])
            self.assertFalse(result['changes_detector_arrays'])
            self.assertFalse(result['changes_biological_classification'])
            self.assertTrue(result['orientation_is_requested_not_joint_strand_disambiguation'])


if __name__ == '__main__':
    unittest.main()
