"""Paired annotations follow count provenance without changing TRM decisions."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from trm_assessment import assess_candidate, classify, PAIR_COUNT_FIELDS


class PairedAssessmentTests(unittest.TestCase):
    def setUp(self):
        self.raw = {'motif': 'AAGCCT', 'reads': 40, 'repeat_bases': 4500,
                    'max_tract': 150, 'full_repeat_reads': 10}
        self.junction = {'full_repeat_records': 10, 'junction_reads': 20,
                         'left_forward_junction_reads': 10, 'left_reverse_junction_reads': 0,
                         'left_ambiguous_junction_reads': 0, 'right_forward_junction_reads': 0,
                         'right_reverse_junction_reads': 10, 'right_ambiguous_junction_reads': 0,
                         'internal_repeat_reads': 0}
        self.target = {**self.raw, 'junction_evidence': self.junction}
        self.context = {'measured': True, 'layout': 'fr', 'scope': 'complete_inputs',
                        'validated_pairs': 1000}
        self.pairs = dict(zip(PAIR_COUNT_FIELDS, (25, 15, 6, 4, 13, 1, 1, 4, 2)))
        self.refs = {'AAGCCT': [{'motif': 'TTAGGC', 'evidence': 'experimental',
                               'source': 'https://example.org/reference'}]}

    def assess(self, raw=None, target=None, **context):
        return assess_candidate(self.raw if raw is None else raw, target, self.refs,
                                {'AAGCCT'}, 150, **context)

    def test_complete_current_unpaired_does_not_inherit_old_pairs(self):
        raw = {**self.raw, 'paired_evidence': self.pairs}
        candidate = self.assess(raw, self.target, raw_paired_end=self.context)
        self.assertEqual(candidate['evidence_count_source'], 'measured_evidence')
        self.assertEqual(candidate['paired_evidence']['status'], 'not_measured')
        self.assertIsNone(candidate['paired_evidence']['supporting_pairs'])

    def test_native_context_required_even_when_candidate_fields_exist(self):
        target = {**self.target, 'paired_evidence': self.pairs}
        for context in (None, {'measured': False}, {**self.context, 'scope': 'unknown'}):
            with self.subTest(context=context):
                candidate = self.assess(target=target, measured_paired_end=context)
                self.assertFalse(candidate['paired_evidence']['measured'])
                self.assertIsNone(candidate['paired_evidence']['opposite_strand_pairs'])

    def test_prefix_observation_is_measured_with_explicit_partial_scope(self):
        target = {**self.target, 'paired_evidence': self.pairs}
        candidate = self.assess(target=target, measured_paired_end={
            **self.context, 'scope': 'synchronized_prefix_only'})
        self.assertTrue(candidate['paired_evidence']['measured'])
        self.assertEqual(candidate['paired_evidence']['status'], 'measured_prefix')
        self.assertEqual(candidate['paired_evidence']['scope'], 'synchronized_prefix_only')
        self.assertEqual(candidate['paired_evidence']['supporting_pairs'], 25)
        complete = self.assess(target=target, measured_paired_end=self.context)
        self.assertEqual({k: v for k, v in candidate.items() if k != 'paired_evidence'},
                         {k: v for k, v in complete.items() if k != 'paired_evidence'})

    def test_geometry_only_target_retains_original_count_and_pair_source(self):
        raw = {**self.raw, 'paired_evidence': self.pairs}
        partial = {'motif': 'AAGCCT', 'junction_evidence': self.junction,
                   'paired_evidence': {key: 0 for key in PAIR_COUNT_FIELDS}}
        candidate = self.assess(raw, partial, raw_paired_end=self.context,
                                measured_paired_end=self.context)
        self.assertEqual(candidate['paired_evidence']['supporting_pairs'], 25)
        self.assertEqual(candidate['paired_evidence']['count_source'], 'original_discovery')

    def test_targeted_counts_and_pairs_are_selected_together(self):
        raw = {**self.raw, 'reads': 80, 'paired_evidence': {**self.pairs, 'supporting_pairs': 999}}
        target = {**self.target, 'paired_evidence': self.pairs}
        candidate = self.assess(raw, target, raw_paired_end=self.context, measured_paired_end=self.context)
        self.assertEqual(candidate['reads'], 40)
        self.assertEqual(candidate['paired_evidence']['supporting_pairs'], 25)
        self.assertEqual(candidate['paired_evidence']['one_mate_only_pairs'], 10)

    def test_measured_zero_is_distinct_from_missing(self):
        zero = {**self.target, 'reads': 0, 'repeat_bases': 0, 'max_tract': 0,
                'junction_evidence': {key: 0 for key in self.junction},
                'paired_evidence': {key: 0 for key in PAIR_COUNT_FIELDS}}
        candidate = self.assess(target=zero, measured_paired_end=self.context)
        self.assertTrue(candidate['paired_evidence']['measured'])
        self.assertEqual(candidate['paired_evidence']['supporting_pairs'], 0)
        self.assertTrue(candidate['measurement_disagreement'])
        self.assertEqual(self.assess(target=self.target, measured_paired_end=self.context)
                         ['paired_evidence']['supporting_pairs'], None)

    def test_pair_orientations_do_not_change_classification_rank_or_candidate_gates(self):
        baseline = self.assess(target=self.target)
        adequacy = {'status': 'unassessable_library_or_genome'}
        for changes in ({}, {'opposite_strand_pairs': 0, 'same_strand_pairs': 15,
                             'ambiguous_strand_pairs': 0}):
            target = {**self.target, 'paired_evidence': {**self.pairs, **changes}}
            candidate = self.assess(target=target, measured_paired_end=self.context)
            self.assertEqual(classify([candidate], adequacy), classify([baseline], adequacy))
            self.assertEqual({key: value for key, value in candidate.items() if key != 'paired_evidence'},
                             {key: value for key, value in baseline.items() if key != 'paired_evidence'})

    def test_malformed_pair_counts_are_errors_not_silent_zeros(self):
        invalid = []
        for field, value in (('supporting_pairs', 26), ('opposite_strand_pairs', 15),
                             ('both_full_repeat_pairs', 5), ('r1_only_repeat_pairs', True),
                             ('same_strand_pairs', -1)):
            invalid.append({**self.pairs, field: value})
        invalid.append({key: value for key, value in self.pairs.items() if key != 'supporting_pairs'})
        for pairs in invalid:
            with self.subTest(pairs=pairs), self.assertRaises(ValueError):
                self.assess(target={**self.target, 'paired_evidence': pairs}, measured_paired_end=self.context)
        with self.assertRaises(ValueError):
            self.assess(target={**self.target, 'paired_evidence': self.pairs},
                        measured_paired_end={**self.context, 'validated_pairs': 2})

    def test_selected_source_objects_are_not_mutated(self):
        raw = {**self.raw, 'paired_evidence': self.pairs}
        target = {**self.target, 'paired_evidence': self.pairs}
        before = deepcopy((raw, target, self.context))
        self.assess(raw, target, raw_paired_end=self.context, measured_paired_end=self.context)
        self.assertEqual((raw, target, self.context), before)


if __name__ == '__main__':
    unittest.main()
