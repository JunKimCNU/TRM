#!/usr/bin/env python3
"""Frozen synthetic controls and an independently enumerated finite-template oracle.

These development checks are not biological sensitivity/specificity estimates.
The protocol is frozen before the first execution of the new module.
"""
import copy
from datetime import datetime, timezone
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import random
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
CONTROL_DIR = ROOT / 'tests/fixtures/phase_bridge'
sys.path.insert(0, str(ROOT / 'scripts'))
SPEC = importlib.util.spec_from_file_location('phase_bridge', ROOT / 'scripts/phase_bridge.py')
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
ORACLE_OBSERVATIONS = []


def rc(sequence):
    return sequence.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def fixtures():
    return json.loads((CONTROL_DIR / 'protocol.json').read_text())['pair_controls']


def evaluate(case):
    return M.evaluate_pair(case['sequence'], copy.deepcopy(case['left']),
                           copy.deepcopy(case['right']), case['motif'],
                           config=copy.deepcopy(case.get('config')),
                           alignments=copy.deepcopy(case.get('alignments')),
                           expected_strand=case.get('expected_strand'))


def exhaustive_finite_template_oracle(query, motif, start_phase, end_phase,
                                      budget, max_indel_run):
    """Enumerate complete edit paths on explicit finite strings, without a DP.

    Fixed endpoint modulo allows finitely many reference lengths under the edit
    budget. Each string is materialized before exhaustive edit-path traversal.
    No cyclic states, state merging, or production alignment code is used.
    """
    minimum, best_paths = None, 0
    enumerated_terminal_paths = 0
    targets_considered = 0
    for target_length in range(max(0, len(query) - budget), len(query) + budget + 1):
        if (start_phase + target_length) % len(motif) != end_phase:
            continue
        target = ''.join(motif[(start_phase + j) % len(motif)]
                         for j in range(target_length))
        targets_considered += 1

        def visit(i, j, cost, last_indel, run):
            nonlocal minimum, best_paths, enumerated_terminal_paths
            if cost > budget:
                return
            if abs((len(query) - i) - (len(target) - j)) > budget - cost:
                return
            if i == len(query) and j == len(target):
                enumerated_terminal_paths += 1
                if minimum is None or cost < minimum:
                    minimum, best_paths = cost, 1
                elif cost == minimum:
                    best_paths += 1
                return
            if i < len(query) and j < len(target):
                visit(i + 1, j + 1, cost + (query[i] != target[j]), None, 0)
            if i < len(query):
                insertion_run = run + 1 if last_indel == 'I' else 1
                if insertion_run <= max_indel_run:
                    visit(i + 1, j, cost + 1, 'I', insertion_run)
            if j < len(target):
                deletion_run = run + 1 if last_indel == 'D' else 1
                if deletion_run <= max_indel_run:
                    visit(i, j + 1, cost + 1, 'D', deletion_run)

        visit(0, 0, 0, None, 0)
    return dict(min_edits=minimum, optimal_path_count_capped2=min(2, best_paths),
                optimal_path_ambiguous=best_paths > 1 if best_paths else None,
                terminal_paths_enumerated=enumerated_terminal_paths,
                explicit_targets_considered=targets_considered)


def oracle_cases():
    short_queries = [''.join(word) for n in range(4)
                     for word in itertools.product('AC', repeat=n)]
    rng = random.Random(20260928)
    six_queries = [''] + [''.join(rng.choice('ACGT') for _ in range(n))
                         for n in range(1, 4) for _ in range(4)]
    for motif in ['AC', 'ACG', 'TTAGGG']:
        queries = six_queries if motif == 'TTAGGG' else short_queries
        for query, start, end, budget, run in itertools.product(
                queries, range(len(motif)), range(len(motif)), range(3), [1, 2]):
            yield query, motif, start, end, budget, run


def no_promotion(test, value):
    if isinstance(value, dict):
        for key, item in value.items():
            if key in {'biological_tracts_merged', 'changes_native_classification',
                       'changes_policy_classification', 'physical_chromosome_end_confirmed',
                       'localization_confirmation'}:
                test.assertIs(item, False, key)
            no_promotion(test, item)
    elif isinstance(value, list):
        for item in value:
            no_promotion(test, item)


class FrozenPairControls(unittest.TestCase):
    def test_protocol_unchanged_and_settings_bound(self):
        receipt = json.loads((CONTROL_DIR / 'protocol.receipt.json').read_text())
        protocol = json.loads((CONTROL_DIR / 'protocol.json').read_text())
        self.assertEqual(receipt['protocol_sha256'],
                         hashlib.sha256((CONTROL_DIR / 'protocol.json').read_bytes()).hexdigest())
        self.assertTrue(receipt['frozen_before_first_module_run'])
        self.assertFalse(receipt['new_module_executed_before_freeze'])
        self.assertEqual(M.validate_config(), protocol['settings'])
        self.assertEqual(receipt['pair_controls'], len(protocol['pair_controls']))

    def test_frozen_outcomes(self):
        for case in fixtures():
            with self.subTest(control=case['id']):
                result = evaluate(case)
                self.assertEqual(result['status'], case['expected_status'])
                if 'expected_min_edits' in case:
                    self.assertEqual(result['bridge']['min_edits'], case['expected_min_edits'])
                if 'expected_path_ambiguous' in case:
                    self.assertEqual(result['bridge']['optimal_path_ambiguous'], case['expected_path_ambiguous'])
                no_promotion(self, result)

    def test_inputs_not_mutated(self):
        for case in fixtures():
            before = copy.deepcopy(case)
            evaluate(case)
            self.assertEqual(case, before)

    def test_no_transitive_join_or_physical_end_promotion(self):
        results = {case['id']: evaluate(case) for case in fixtures()
                   if case['id'] in {'chain_ab', 'chain_bc', 'chain_ac', 'internal_repeat_context'}}
        for key in ['chain_ab', 'chain_bc', 'internal_repeat_context']:
            self.assertEqual(results[key]['status'], 'compatible_under_bounded_phase_bridge')
            no_promotion(self, results[key])
        self.assertEqual(results['chain_ac']['status'], 'rejected_by_bounded_phase_bridge')

    def test_reverse_complement_outcomes(self):
        for case in fixtures():
            if case.get('alignments') or case['id'] in {'overlapping_anchors'}:
                continue
            n = len(case['sequence'])
            def flip(anchor):
                return dict(anchor, start0=n-anchor['end0'], end0=n-anchor['start0'],
                            strand={'+': '-', '-': '+'}[anchor['strand']])
            reversed_case = dict(case, sequence=rc(case['sequence']), left=flip(case['right']),
                                 right=flip(case['left']),
                                 expected_strand={'+': '-', '-': '+', None: None}[case.get('expected_strand')])
            with self.subTest(control=case['id']):
                forward_result, reverse_result = evaluate(case), evaluate(reversed_case)
                self.assertEqual(forward_result['status'], reverse_result['status'])
                if forward_result['bridge'] is not None and reverse_result['bridge'] is not None:
                    for key in ['min_edits', 'optimal_path_count_capped2', 'optimal_path_ambiguous']:
                        self.assertEqual(forward_result['bridge'][key], reverse_result['bridge'][key])

    def test_invalid_configuration(self):
        for config in [{'max_edits': True}, {'max_gap_bp': -1}, {'max_span_bp': 0},
                       {'max_edit_fraction': float('nan')}, {'max_edit_fraction': float('inf')},
                       {'max_edit_fraction': 1.1}, {'max_indel_run': -1},
                       {'min_anchor_copies': 0}, {'not_an_option': 1}]:
            with self.subTest(config=config), self.assertRaises(ValueError):
                M.validate_config(config)

    def test_bridge_path_consumes_all_gap_and_keeps_source_coordinates(self):
        for case in fixtures():
            result = evaluate(case)
            if result['status'] != 'compatible_under_bounded_phase_bridge':
                continue
            bridge = result['bridge']
            start, end = case['left']['end0'], case['right']['start0']
            query_cursor, template_cursor = start, 0
            motif = case['motif'] if case['left']['strand'] == '+' else rc(case['motif'])
            counts = dict.fromkeys('MXID', 0)
            for operation in bridge['path_ops']:
                with self.subTest(control=case['id'], operation=operation):
                    self.assertEqual(operation['query_start0'], query_cursor)
                    self.assertEqual(operation['template_start0'], template_cursor)
                    op, length = operation['operation'], operation['length_bp']
                    counts[op] += length
                    query_cursor += length if op != 'D' else 0
                    template_cursor += length if op != 'I' else 0
                    self.assertEqual(operation['query_end0'], query_cursor)
                    self.assertEqual(operation['template_end0'], template_cursor)
                    self.assertLessEqual(query_cursor, end)
                    self.assertEqual(operation['template_end_phase'],
                                     (bridge['start_phase'] + template_cursor) % len(motif))
                    observed = case['sequence'][operation['query_start0']:query_cursor]
                    expected = ''.join(motif[(bridge['start_phase'] + j) % len(motif)]
                                       for j in range(operation['template_start0'], template_cursor))
                    if op == 'M':
                        self.assertEqual(observed, expected)
                    elif op == 'X':
                        self.assertTrue(all(a != b for a, b in zip(observed, expected)))
                    if op in 'ID':
                        self.assertLessEqual(length, result['config']['max_indel_run'])
            self.assertEqual(query_cursor, end)
            self.assertEqual(bridge['gap_start0'], start)
            self.assertEqual(bridge['gap_end0'], end)
            self.assertEqual(bridge['query_consumed_bp'], end - start)
            self.assertEqual(bridge['template_consumed_bp'], template_cursor)
            self.assertEqual((bridge['start_phase'] + template_cursor) % len(motif), bridge['end_phase'])
            self.assertEqual(bridge['operation_counts'], counts)
            self.assertEqual(bridge['min_edits'], counts['X'] + counts['I'] + counts['D'])
            self.assertTrue(bridge['all_query_bases_consumed'])
            self.assertFalse(bridge['free_end_clipping_allowed'])
            self.assertFalse(bridge['free_phase_changes_allowed'])

    def test_mapping_context_stays_descriptive(self):
        case = next(case for case in fixtures() if case['id'] == 'split_primary_mapping')
        with_mapping = evaluate(case)
        without = evaluate(dict(case, alignments=None))
        self.assertEqual(with_mapping['status'], without['status'])
        self.assertEqual(with_mapping['bridge'], without['bridge'])
        self.assertEqual(with_mapping['alignment_evidence']['status'], 'split_or_competing_primary_mapping')
        self.assertFalse(with_mapping['alignment_evidence']['localization_confirmation'])

    def test_old_gap_boundary_effect_is_explicit(self):
        import tract_context
        case = next(case for case in fixtures() if case['id'] == 'short_perfect_five')
        old = tract_context.evaluate_pair(case['sequence'], case['left'], case['right'], case['motif'])
        self.assertEqual(old['status'], case['expected_old_status'])
        self.assertEqual(evaluate(case)['bridge']['min_edits'], 0)

    def test_resource_limit_is_not_sequence_rejection(self):
        case = next(case for case in fixtures() if case['id'] == 'perfect_two_copies')
        result = evaluate(dict(case, config={'max_dp_states': 1}))
        self.assertEqual(result['status'], 'not_evaluated_computation_limit')
        self.assertIsNone(result['bridge']['min_edits'])
        self.assertFalse(result['bridge']['evaluated'])
        no_promotion(self, result)

    def test_reverse_complement_equivalent_motif_has_no_resolved_direction(self):
        sequence = 'ACGT' * 15
        left = dict(tract_id='L', start0=0, end0=20, strand='+')
        right = dict(tract_id='R', start0=24, end0=44, strand='+')
        result = M.evaluate_pair(sequence, left, right, 'ACGT')
        self.assertEqual(result['status'], 'rejected_by_bounded_phase_bridge')
        self.assertIn('anchor_template_strand_ambiguous', result['reason_codes'])

    def test_invalid_coordinates_phases_and_motif_are_rejected(self):
        case = fixtures()[0]
        for change in [dict(start0=-1), dict(end0=len(case['sequence']) + 1),
                       dict(start0=True), dict(strand='?'), dict(tract_id='')]:
            with self.subTest(change=change), self.assertRaises(ValueError):
                M.evaluate_pair(case['sequence'], dict(case['left'], **change),
                                case['right'], case['motif'])
        for motif in ['', 'NNN', 'ACAC', 'A' * 65]:
            with self.subTest(motif=motif), self.assertRaises(ValueError):
                M.align_bridge('AC', motif, 0, 0, 1, 1)
        for start, end in [(-1, 0), (0, 2), (True, 0), (0, 1.0)]:
            with self.subTest(start=start, end=end), self.assertRaises(ValueError):
                M.align_bridge('AC', 'AC', start, end, 1, 1)

    def test_zero_indel_limit_disables_insertions_and_deletions(self):
        accepted = M.align_bridge('ACG', 'ACG', 0, 0, 1, 0)
        self.assertEqual(accepted['min_edits'], 0)
        rejected = M.align_bridge('ACCG', 'ACG', 0, 0, 1, 0)
        self.assertEqual(rejected['status'], 'rejected_by_bounded_phase_bridge')
        self.assertIsNone(rejected['min_edits'])


class IndependentFiniteTemplateOracle(unittest.TestCase):
    def test_exhaustive_tiny_alignments(self):
        ORACLE_OBSERVATIONS.clear()
        for query, motif, start, end, budget, run in oracle_cases():
            expected = exhaustive_finite_template_oracle(query, motif, start, end, budget, run)
            observed = M.align_bridge(query, motif, start, end, budget, run)
            parameters = dict(query=query, motif=motif, start_phase=start, end_phase=end,
                              max_edits=budget, max_indel_run=run)
            with self.subTest(**parameters):
                self.assertEqual(observed['min_edits'], expected['min_edits'])
                self.assertEqual(observed['optimal_path_count_capped2'], expected['optimal_path_count_capped2'])
                self.assertEqual(observed['optimal_path_ambiguous'], expected['optimal_path_ambiguous'])
            ORACLE_OBSERVATIONS.append(dict(**parameters, expected=expected,
                                           observed={key: observed[key] for key in
                                                     ('min_edits', 'optimal_path_count_capped2', 'optimal_path_ambiguous')}))


def write_outcomes(destination, program):
    observed = []
    for case in fixtures():
        result = evaluate(case)
        observed.append(dict(control_id=case['id'], expected_status=case['expected_status'],
                             observed_status=result['status'],
                             status_matches=result['status'] == case['expected_status'], result=result))
    result = dict(recorded_utc=datetime.now(timezone.utc).isoformat(),
                  protocol_sha256=hashlib.sha256((CONTROL_DIR / 'protocol.json').read_bytes()).hexdigest(),
                  module_sha256=hashlib.sha256((ROOT / 'scripts/phase_bridge.py').read_bytes()).hexdigest(),
                  test_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  tests_run=program.result.testsRun, tests_successful=program.result.wasSuccessful(),
                  failures=len(program.result.failures), errors=len(program.result.errors),
                  scope='Synthetic development controls; not biological sensitivity or specificity',
                  controls=observed, oracle_comparisons=len(ORACLE_OBSERVATIONS),
                  oracle_results=ORACLE_OBSERVATIONS)
    Path(destination).write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    destination = None
    if '--write-outcomes' in sys.argv:
        i = sys.argv.index('--write-outcomes')
        destination = sys.argv[i + 1]
        del sys.argv[i:i + 2]
    program = unittest.main(exit=False)
    if destination:
        write_outcomes(destination, program)
    sys.exit(not program.result.wasSuccessful())
