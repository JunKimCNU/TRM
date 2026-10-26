"""Primitive cyclic family observations remain descriptive and source-separated."""
from __future__ import annotations

from copy import deepcopy
import importlib.util
import itertools
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('family_observation_reporting', ROOT / 'scripts/family_observation_reporting.py')
reporting = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(reporting)
build = reporting.build_family_observation_report


def report(expected='TTAGGG', *, raw=None, candidates=None, observations=None, assignment=None, **kwargs):
    if raw is None:
        raw = {'summary': {'sampled_prefix': False}, 'config': {'min_k': 3, 'max_k': 32, 'seed_bases': 24, 'min_copies': 3}}
    if observations is None:
        observations = [{'motif': expected, 'reads': None, 'recovery': 'not_reported'}] if expected else []
    return build(raw, candidates or [], observations, assignment or {}, **kwargs)


def first(*args, **kwargs):
    return report(*args, **kwargs)['families'][0]


def levenshtein(left, right):
    row = list(range(len(right) + 1))
    for i, a in enumerate(left, 1):
        nxt = [i]
        for j, b in enumerate(right, 1):
            nxt.append(min(nxt[-1] + 1, row[j] + 1, row[j - 1] + (a != b)))
        row = nxt
    return row[-1]


class FamilyObservationReportingTests(unittest.TestCase):
    def test_rejected_exact_reverse_complement_power_is_observed(self):
        raw = {'candidates': [], 'rejected_candidates': [{'motif': 'CCCTAACCCTAA', 'reads': 9}]}
        row = first('TTAGGG', raw=raw)
        self.assertEqual(row['status'], 'exact_family_observed')
        self.assertEqual(row['canonical_family'], 'AACCCT')
        self.assertEqual(row['exact_original']['reads'], 9)
        self.assertEqual(row['exact_current']['read_counts'], [])
        self.assertEqual(row['exact_rows'][0]['source_bucket'], 'raw.rejected_candidates')

    def test_guillardia_single_change_is_descriptive(self):
        known = 'AAGAGAGAGAGAGAGAAGGGGGG'
        alternate = 'AAGAGAGAGAGAGAGAGGGGGGG'
        raw = {'rejected_candidates': [{'motif': alternate, 'reads': 298, 'max_tract': 149}]}
        row = first(known, raw=raw)
        nearby = row['nearby_families_all'][0]
        self.assertEqual(row['status'], 'nearby_family_observed')
        self.assertEqual(row['exact_measurement_status'], 'not_reported')
        self.assertEqual((nearby['edit_distance'], nearby['operation']), (1, 'substitution'))
        self.assertEqual(nearby['original']['reads'], 298)
        self.assertIn('not a canonical or novel', nearby['interpretation'])

    def test_cyclic_indel_is_not_linear_canonical_comparison(self):
        expected, observed = reporting._family('AAGTC'), reporting._family('AGTC')
        self.assertGreater(levenshtein(expected, observed), 1)
        row = first('AAGTC', raw={'candidates': [{'motif': 'AGTC', 'reads': 7}]})
        witness = row['nearby_families'][0]['edit_witness']
        self.assertEqual(witness['operation'], 'deletion')
        self.assertEqual(sum(a != b for a, b in zip(witness['aligned_expected_unit'], witness['aligned_observed_unit'])), 1)
        self.assertEqual(witness['aligned_observed_unit'].replace('-', ''), witness['observed_unit'])
        reverse = first('AGTC', raw={'candidates': [{'motif': 'AAGTC', 'reads': 7}]})
        self.assertEqual(reverse['nearby_families'][0]['operation'], 'insertion')

    def test_one_edit_short_units_matches_full_cyclic_dp_oracle(self):
        families = sorted({reporting._family(''.join(s)) for n in (3, 4, 5) for s in itertools.product('AC', repeat=n)
                           if len(reporting._unit(''.join(s))) >= 3})
        for expected, observed in itertools.product(families, repeat=2):
            if expected == observed:
                self.assertIsNone(reporting._nearby_witness(expected, observed))
                continue
            # Independent full dynamic programming over both cyclic cuts/strands.
            oracle = min(levenshtein(e, o) for _, _, e in reporting._rotations(expected)
                         for _, _, o in reporting._rotations(observed))
            self.assertEqual(reporting._nearby_witness(expected, observed) is not None, oracle == 1,
                             (expected, observed, oracle))

    def test_zero_absent_and_other_positive_are_different(self):
        zero = first(observations=[{'motif': 'TTAGGG', 'reads': 0, 'recovery': 'explicit_zero', 'source': 'expected_family_diagnostic'}])
        missing = first()
        other = first(raw={'candidates': [{'motif': 'ACAT', 'reads': 200}]})
        self.assertEqual(zero['exact_measurement_status'], 'explicit_zero')
        self.assertEqual(missing['exact_measurement_status'], 'not_reported')
        self.assertEqual(other['status'], 'other_repeat_families_observed')
        self.assertEqual(other['other_positive_family_count'], 1)
        self.assertNotEqual(other['status'], 'no_positive_candidates_reported')

    def test_not_reported_zero_placeholder_is_not_measurement(self):
        row = first(candidates=[{'canonical_family': 'AACCCT', 'reads': 0, 'original_discovery_reads': 0,
                                  'evidence_count_source': 'original_discovery'}],
                    observations=[{'motif': 'TTAGGG', 'reads': 0, 'recovery': 'not_reported'}])
        self.assertEqual(row['exact_measurement_status'], 'not_reported')
        self.assertEqual(row['exact_current']['read_counts'], [])

    def test_native_discovery_zero_is_reported_but_not_measured_zero(self):
        raw = {'candidates': [{'motif': 'TTAGGG', 'reads': 0}]}
        row = first(raw=raw, candidates=[{'motif': 'TTAGGG', 'reads': 0,
                                         'evidence_count_source': 'original_discovery'}])
        self.assertEqual(row['exact_measurement_status'], 'reported_zero_unverified')
        self.assertEqual(row['exact_original']['read_counts'], [0])
        self.assertFalse(row['exact_original']['explicit_zero_observed'])
        self.assertFalse(row['exact_current']['explicit_zero_observed'])
        raw['analysis_kind'] = 'targeted_evidence'
        measured = first(raw=raw)
        self.assertEqual(measured['exact_measurement_status'], 'explicit_zero')
        unspecified = first(observations=[{'motif': 'TTAGGG', 'reads': 0}])
        self.assertEqual(unspecified['exact_measurement_status'], 'reported_zero_unverified')

    def test_duplicate_conflicting_support_counts_not_hidden_by_equal_reads(self):
        raw = {'candidates': [{'motif': 'TTAGGG', 'reads': 11, 'full_repeat_reads': 2}],
               'rejected_candidates': [{'motif': 'CCCTAA', 'reads': 11, 'full_repeat_reads': 3}]}
        row = first(raw=raw)
        self.assertEqual(row['exact_measurement_status'], 'measurement_disagreement')
        self.assertEqual(row['exact_original']['reads'], 11)
        self.assertEqual(row['exact_original']['conflicting_count_fields'], ['full_repeat_reads'])

    def test_current_positive_and_diagnostic_zero_conflict(self):
        row = first(candidates=[{'canonical_family': 'AACCCT', 'reads': 17, 'original_discovery_reads': 0,
                                  'evidence_count_source': 'targeted_only'}],
                    observations=[{'motif': 'TTAGGG', 'reads': 0, 'recovery': 'explicit_zero'}])
        self.assertEqual(row['status'], 'exact_family_observed')
        self.assertEqual(row['exact_measurement_status'], 'measurement_disagreement')
        self.assertEqual(row['exact_current']['read_counts'], [0, 17])
        self.assertIsNone(row['exact_current']['reads'])
        self.assertEqual(row['exact_original']['read_counts'], [])  # fabricated historical0 is excluded

    def test_historical_positive_current_zero_remains_visible(self):
        row = first(candidates=[{'canonical_family': 'AACCCT', 'reads': 0, 'original_discovery_reads': 31,
                                  'evidence_count_source': 'measured_evidence', 'measurement_disagreement': True}],
                    observations=[{'motif': 'TTAGGG', 'reads': 0, 'recovery': 'measurement_disagreement'}])
        self.assertEqual(row['status'], 'exact_family_observed')
        self.assertTrue(row['measurement_disagreement'])
        self.assertEqual(row['exact_original']['read_counts'], [31])
        self.assertEqual(row['exact_current']['read_counts'], [0])

    def test_duplicate_alias_rows_never_sum_and_conflicts_preserved(self):
        raw = {'candidates': [{'motif': 'TTAGGG', 'reads': 11}],
               'rejected_candidates': [{'motif': 'CCCTAA', 'reads': 11}, {'motif': 'TAGGGT', 'reads': 12}]}
        row = first(raw=raw)
        self.assertEqual(len(row['exact_original']['rows']), 3)
        self.assertEqual(row['exact_original']['read_counts'], [11, 12])
        self.assertIsNone(row['exact_original']['reads'])
        self.assertEqual(row['exact_measurement_status'], 'measurement_disagreement')
        raw['rejected_candidates'].pop()
        row = first(raw=raw)
        self.assertEqual(row['exact_original']['reads'], 11)
        self.assertFalse(row['measurement_disagreement'])

    def test_all_candidates_examined_and_preview_complete_counted(self):
        motif = 'AAGTCG'
        alternatives = set()
        for i in range(len(motif)):
            for b in 'ACGT':
                candidate = motif[:i] + b + motif[i + 1:]
                if reporting._family(candidate) != reporting._family(motif):
                    alternatives.add(candidate)
        raw = {'candidates': [{'motif': 'ACAT', 'reads': 9999}] * 200,
               'rejected_candidates': [{'motif': m, 'reads': 1} for m in sorted(alternatives)]}
        row = first(motif, raw=raw)
        self.assertGreater(row['nearby_family_count'], 5)
        self.assertEqual(len(row['nearby_families']), 5)
        self.assertEqual(len(row['nearby_families_all']), row['nearby_family_count'])
        self.assertEqual(row['nearby_omitted_count'], row['nearby_family_count'] - 5)
        self.assertTrue(row['nearby_truncated'])
        self.assertTrue(all('raw.rejected_candidates' in n['buckets'] for n in row['nearby_families_all']))

    def test_nearby_survives_exact_positive_and_historical_only(self):
        row = first('AAGTC', raw={'candidates': [{'motif': 'AAGTC', 'reads': 12}]},
                    candidates=[{'motif': 'AGTC', 'reads': 0, 'original_discovery_reads': 8,
                                 'evidence_count_source': 'measured_evidence'}])
        self.assertEqual(row['status'], 'exact_family_observed')
        self.assertEqual(row['nearby_family_count'], 1)
        self.assertEqual(row['nearby_families'][0]['original']['reads'], 8)
        self.assertEqual(row['nearby_families'][0]['current']['reads'], 0)
        self.assertTrue(row['nearby_families'][0]['measurement_disagreement'])

    def test_missing_reference_is_not_a_negative_even_with_positives(self):
        value = report(None, raw={'rejected_candidates': [{'motif': 'ACAT', 'reads': 8}]},
                       assignment={'reference_status': 'missing', 'mapping_status': 'global_prior_no_curated_reference'})
        self.assertEqual(value['status'], 'no_reference')
        self.assertEqual(value['families'], [])
        self.assertEqual(value['positive_saved_family_count'], 1)
        self.assertEqual(value['reference_context']['reference_status'], 'missing')

    def test_failure_availability_and_partial_context(self):
        self.assertEqual(first(source_available=False)['status'], 'unavailable')
        self.assertEqual(first(analysis_complete=False)['status'], 'analysis_incomplete')
        partial = first(raw={'summary': {'sampled_prefix': True}})
        self.assertEqual(partial['measurement_context_status'], 'partial_input')
        self.assertEqual(partial['status'], 'analysis_incomplete')
        positive = first(raw={'candidates': [{'motif': 'TTAGGG', 'reads': 2}]}, analysis_complete=False)
        self.assertEqual(positive['status'], 'exact_family_observed')
        self.assertEqual(positive['measurement_context_status'], 'analysis_incomplete')

    def test_actual_de_novo_range_not_targeted_eligibility(self):
        motif = 'A'*32 + 'C'  # primitive33, inside fixed64 reporting budget
        value = report(motif, raw={'analysis_kind': 'targeted_evidence', 'config': {'min_k': 3, 'max_k': 32, 'seed_bases': 24, 'min_copies': 2},
                                  'candidate_evidence': [{'motif': motif, 'reads': 5}]})
        row = value['families'][0]
        self.assertTrue(row['outside_de_novo_search_range'])
        self.assertEqual(row['status'], 'exact_family_observed')
        self.assertEqual(row['nearby_comparison_status'], 'assessed')
        self.assertEqual(value['comparison_context']['native_config']['max_k'], 32)
        self.assertEqual(value['comparison_context']['nearby_max_primitive_period'], 64)

    def test_unsupported_invalid_and_small_exact_units(self):
        for invalid in ('NNAGGG', 'AG R', 'A'*64 + 'C'):
            row = first(invalid)
            self.assertEqual(row['status'], 'unsupported_expected_unit')
            self.assertEqual(row['nearby_comparison_status'], 'unsupported_period_or_invalid_unit')
        small = first('AAAA', raw={'candidate_evidence': [{'motif': 'A', 'reads': 4}], 'config': {'min_k': 1, 'max_k': 64}})
        self.assertEqual(small['status'], 'exact_family_observed')
        self.assertEqual(small['primitive_period'], 1)
        self.assertEqual(small['nearby_comparison_status'], 'unsupported_period_or_invalid_unit')

    def test_bad_counts_are_not_positive_or_zero(self):
        for bad in (-1, True, 3.5, '12'):
            row = first(raw={'candidates': [{'motif': 'TTAGGG', 'reads': bad}]})
            self.assertEqual(row['exact_measurement_status'], 'invalid_counts')
            self.assertEqual(row['exact_original']['read_counts'], [])
            self.assertNotEqual(row['status'], 'exact_family_observed')

    def test_reference_role_scope_and_no_input_mutations(self):
        raw = {'candidates': [{'motif': 'TTAGGG', 'reads': 2}], 'config': {'max_k': 64}}
        candidates = [{'motif': 'TTAGGG', 'reads': 2, 'original_discovery_reads': 2, 'candidate_status': 'unchanged'}]
        observed = [{'motif': 'TTAGGG', 'reads': 2, 'recovery': 'reported_positive', 'reference_status': 'candidate'}]
        assignment = {'review_motifs': 'TTAGGG', 'reference_scope': 'nucleomorph', 'review_reference_ids': 'x', 'reference_status': 'supported'}
        original = deepcopy((raw, candidates, observed, assignment))
        value = build(raw, candidates, observed, assignment)
        self.assertEqual((raw, candidates, observed, assignment), original)
        row = value['families'][0]
        self.assertEqual(row['reference_role'], 'review_candidate')
        self.assertEqual(row['reference_scope'], 'nucleomorph')
        self.assertEqual(row['reference_status'], 'supported')
        self.assertEqual(row['candidate_catalog_statuses'], ['candidate'])
        row['reference_rows'][0]['reads'] = 100
        self.assertEqual((raw, candidates, observed, assignment), original)


if __name__ == '__main__':
    unittest.main()
