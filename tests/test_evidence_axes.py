"""Regression cases for simultaneous observation and interpretation views."""
import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from evidence_axes import build_evidence_axes as build_axes

def build_evidence_axes(document):
    return build_axes(document, global_catalogue={
        "AACCCT": [{"evidence": "experimental_terminal", "reference_status": "supported"}]})
from family_observation_reporting import build_family_observation_report


def fixture(reads=1, original=1, scope='clade'):
    candidate = dict(motif='TTAGGG', canonical_family='AACCCT', reads=reads,
        full_repeat_reads=0, reference_status='established', candidate_status='weak_or_conflicting_known_repeat',
        evidence_count_source='targeted_only', original_discovery_reads=original)
    expected = dict(motif='TTAGGG', canonical_family='AACCCT', reads=reads,
        recovery='explicit_zero' if reads == 0 else 'reported_positive')
    raw = {'candidates': [dict(motif='AACCCT', reads=original)]} if original else {}
    assignment = dict(related_motifs='TTAGGG', reference_scope=scope, reference_status='supported')
    return {'candidates': [candidate], 'analysis_complete': True,
        'family_observation': build_family_observation_report(raw, [candidate], [expected], assignment),
        'reference_agreement': {'main_trm_status': 'not_confirmed'},
        'data_adequacy': {'status': 'unassessable_library_or_genome', 'missing_assumptions': ['telomere_retention_fraction']},
        'library_context': {'limited': False, 'reasons': []}, 'externally_confirmed': False,
        'classification': 'novel_geometry_nomination'}


class EvidenceAxesTests(unittest.TestCase):
    def test_known_trace_is_visible_beside_novel_candidate_without_changing_policy(self):
        d = fixture()
        d['candidates'].append(dict(motif='AGGGTTAGGTT', canonical_family='AACCCTAACC', reads=500,
            full_repeat_reads=20, reference_status='not_in_catalogue', geometry_nomination=True))
        before = copy.deepcopy(d)
        axes = build_evidence_axes(d)
        self.assertEqual(d, before)
        self.assertEqual(axes['catalogued_positive_family_count'], 1)
        self.assertEqual(axes['alternative_candidate_count'], 1)
        self.assertIn('1 reads', axes['display']['expected_observation_axis'])
        self.assertIn('분류군 기본 기대', axes['display']['expected_observation_axis'])
        self.assertEqual(axes['main_trm_status'], 'not_confirmed')
        self.assertEqual(axes['sample_physical_terminal_status'], 'not_confirmed')

    def test_current_zero_does_not_hide_historical_positive(self):
        axes = build_evidence_axes(fixture(reads=0, original=12))
        e = axes['expected'][0]
        self.assertEqual(e['current_reads'], 0)
        self.assertEqual(e['original_reads'], 12)
        self.assertTrue(e['original_positive'])
        self.assertIn('불일치', axes['display']['expected_observation_axis'])
        self.assertIn('원래 12', axes['display']['expected_observation_axis'])

    def test_missing_is_not_zero_and_species_is_not_clade(self):
        d = fixture(scope='species')
        d['candidates'] = []
        d['family_observation'] = build_family_observation_report({}, [], [],
            dict(related_motifs='TTAGGG', reference_scope='species', reference_status='supported'))
        axes = build_evidence_axes(d)
        self.assertIsNone(axes['expected'][0]['current_reads'])
        self.assertIn('0 아님', axes['display']['expected_observation_axis'])
        self.assertIn('종별 문헌', axes['display']['expected_observation_axis'])
        zero = build_evidence_axes(fixture(reads=0, original=0))
        self.assertIn('측정 0', zero['display']['expected_observation_axis'])

    def test_library_limits_do_not_suppress_positive_and_power_is_separate(self):
        d = fixture(reads=40)
        d['library_context'] = {'limited': True, 'reasons': ['RNA_library']}
        axes = build_evidence_axes(d)
        self.assertEqual(axes['catalogued_positive_family_count'], 1)
        self.assertEqual(axes['positive_interpretation_limits'], ['RNA_library'])
        self.assertEqual(axes['negative_power_missing_assumptions'], ['telomere_retention_fraction'])

    def test_species_only_reference_is_not_global_catalogue(self):
        d = fixture()
        axes = build_axes(d, global_catalogue={})
        self.assertEqual(axes['catalogued_positive_family_count'], 0)
        self.assertTrue(axes['expected'][0]['current_positive'])

    def test_conflicting_counts_and_candidate_reference_are_not_missing_or_established(self):
        d = fixture()
        f = d['family_observation']['families'][0]
        f.update(reference_status='candidate', reference_scope='species',
                 exact_measurement_status='measurement_disagreement')
        f['exact_current'].update(reads=None, read_counts=[5, 7], positive_observed=True)
        text = build_evidence_axes(d)['display']['expected_observation_axis']
        self.assertIn('5/7', text)
        self.assertIn('문헌 후보/미확립', text)
        self.assertNotIn('현재 미보고', text)

    def test_mixed_supported_assignment_retains_review_candidate_caveat(self):
        d = fixture(scope='species')
        f = d['family_observation']['families'][0]
        f.update(reference_status='supported', reference_role='review_candidate')
        self.assertIn('문헌 후보/미확립', build_evidence_axes(d)['display']['expected_observation_axis'])

    def test_no_reference_and_malformed_counts(self):
        d = fixture()
        d['family_observation'] = {}
        d['candidates'][0]['reads'] = '0'
        axes = build_evidence_axes(d)
        self.assertEqual(axes['expected'], [])
        self.assertEqual(axes['catalogued_positive_family_count'], 0)
        self.assertIn('미지정', axes['display']['expected_observation_axis'])


if __name__ == '__main__':
    unittest.main()
