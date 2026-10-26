"""Reporting contracts: preserve evidence and separate sequence from TRM calls."""
from __future__ import annotations

import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from observation_reporting import build_observation_report


class ObservationReportingTests(unittest.TestCase):
    def candidate(self, motif, reads, score, *, original=None, source='measured_evidence'):
        return {'motif': motif, 'canonical_family': motif, 'reads': reads,
                'original_discovery_reads': reads if original is None else original,
                'ranking_score': score, 'full_repeat_reads': 0, 'repeat_bases': reads * 24,
                'max_tract_bp': 24, 'evidence_count_source': source,
                'measurement_disagreement': reads == 0 and bool(original),
                'candidate_status': 'uncatalogued_repeat'}

    def native(self, motif, reads, rank, score):
        return {'motif': motif, 'reads': reads, 'rank_score': score,
                'telomere_prior': {'discovery_rank': rank}, 'repeat_bases': reads * 24}

    def report(self, raw, candidates, representative=None, **changes):
        args = dict(source_available=True, analysis_complete=True,
                    classification='detection_inconclusive', legacy_result='미검출: 검출력/라이브러리 정보 부족',
                    agreement={'status': 'no_supported_reference', 'main_trm_status': 'not_confirmed'},
                    adequacy_status='unassessable_library_or_genome')
        args.update(changes)
        snapshot = copy.deepcopy((raw, candidates, representative, args))
        result = build_observation_report(raw, candidates, representative, **args)
        self.assertEqual(snapshot, (raw, candidates, representative, args), 'Reporting mutated decision/evidence inputs')
        return result

    def test_c24_leader_is_visible_when_reference_withheld_representative_differs(self):
        leader, selected = self.candidate('AAACCCT', 15888, 1000), self.candidate('AACCT', 315, 100)
        raw = {'candidates': [self.native('AACCT', 315, 1, 100)],
               'rejected_candidates': [self.native('AAACCCT', 15888, 0, 1000)]}
        report = self.report(raw, [leader, selected], selected,
                             classification='known_repeat_detected_boundary_unresolved', legacy_result='알려진 반복 검출·말단 경계 미확인')
        view = report['repeat_observation']
        self.assertEqual(view['native_discovery_leader']['canonical_family'], 'AAACCCT')
        self.assertEqual(view['native_discovery_leader']['reads'], 15888)
        self.assertEqual(view['native_discovery_leader']['discovery_rank_one_based'], 1)
        self.assertEqual(view['policy_score_leader']['canonical_family'], 'AAACCCT')
        self.assertEqual(report['policy_selection']['representative_family'], 'AACCT')
        self.assertEqual(report['policy_selection']['representative_policy_rank'], 2)
        self.assertTrue(report['policy_selection']['native_vs_policy_score_same_family'])
        self.assertFalse(report['policy_selection']['native_vs_representative_same_family'])
        self.assertEqual(report['decision_summary']['main_trm_status'], 'not_confirmed')

    def test_l84_positive_without_representative_is_not_displayed_as_nondetection(self):
        candidate = self.candidate('AACACCAC', 158, 1000)
        raw = {'candidates': [], 'rejected_candidates': [self.native('AACACCAC', 158, 0, 1000)]}
        report = self.report(raw, [candidate])
        self.assertEqual(report['repeat_observation']['positive_family_count'], 1)
        self.assertEqual(report['repeat_observation']['policy_score_leader']['reads'], 158)
        self.assertIsNone(report['policy_selection']['representative_motif'])
        self.assertIn('반복 서열 관측', report['decision_summary']['display_label'])
        self.assertNotIn('미검출', report['decision_summary']['display_label'])
        self.assertEqual(report['decision_summary']['classification'], 'detection_inconclusive')
        self.assertIn('미검출', report['decision_summary']['legacy_result'])

    def test_native_rank_not_score_or_array_order_drives_native_leader(self):
        # Assembly terminal ordering can outrank a numerically higher score.
        raw = {'candidates': [self.native('ACAT', 40, 1, 9000), self.native('AACCT', 30, 0, 40)]}
        report = self.report(raw, [self.candidate('ACAT', 40, 9000), self.candidate('AACCT', 30, 40)])
        self.assertEqual(report['repeat_observation']['native_discovery_leader']['canonical_family'], 'AACCT')
        self.assertEqual(report['repeat_observation']['policy_score_leader']['canonical_family'], 'ACAT')

    def test_missing_duplicate_or_conflicting_native_ranks_are_unavailable(self):
        for ranks in ((None, 0), (0, 0), (True, 1), (1, 3)):
            with self.subTest(ranks=ranks):
                raw = {'candidates': [self.native('ACAT', 40, ranks[0], 100), self.native('AACCT', 30, ranks[1], 20)]}
                report = self.report(raw, [self.candidate('ACAT', 40, 100)])
                self.assertIsNone(report['repeat_observation']['native_discovery_leader'])
                self.assertIsNotNone(report['repeat_observation']['policy_score_leader'])
        raw = {'candidates': [self.native('ACAT', 40, 0, 100)],
               'rejected_candidates': [self.native('ACAT', 41, 0, 100)]}
        self.assertIsNone(self.report(raw, [])['repeat_observation']['native_discovery_leader'])

    def test_original_positive_and_current_zero_stay_separate(self):
        raw = {'candidates': [self.native('ACAT', 40, 0, 100)]}
        report = self.report(raw, [self.candidate('ACAT', 0, 0, original=40)],
                             classification='detection_inconclusive_measurement', legacy_result='판정 보류: 탐지 경로 간 관측 불일치')
        view = report['repeat_observation']
        self.assertEqual(view['status'], 'original_positive_only')
        self.assertEqual(view['native_discovery_leader']['reads'], 40)
        self.assertIsNone(view['policy_score_leader'])
        self.assertEqual(view['original_positive_current_nonpositive'][0]['reads'], 0)
        self.assertEqual(report['decision_summary']['display_label'], '판정 보류: 탐지 경로 간 관측 불일치')

    def test_targeted_only_positive_has_no_invented_native_rank(self):
        report = self.report({'candidates': []}, [self.candidate('ACAT', 40, 100, original=0, source='targeted_only')])
        view = report['repeat_observation']
        self.assertEqual(view['status'], 'positive_observations')
        self.assertIsNone(view['native_discovery_leader'])
        self.assertEqual(view['policy_score_leader']['evidence_count_source'], 'targeted_only')

    def test_missing_no_report_and_explicit_zero_are_distinct(self):
        missing = self.report({}, [], source_available=False, analysis_complete=False,
                              classification='analysis_failed', legacy_result='분석 실패: 입력/처리 오류')
        self.assertEqual(missing['repeat_observation']['status'], 'unavailable')
        self.assertEqual(missing['repeat_observation']['measurement_status'], 'source_unavailable')
        self.assertEqual(self.report({}, [])['repeat_observation']['status'], 'no_positive_reported')
        zero = self.report({}, [self.candidate('ACAT', 0, 0)])
        self.assertEqual(zero['repeat_observation']['status'], 'reported_targets_zero')
        self.assertIsNone(zero['repeat_observation']['policy_score_leader'])

    def test_failed_partial_qc_and_library_statuses_preserve_precedence(self):
        candidate = self.candidate('ACAT', 40, 100)
        raw = {'candidates': [self.native('ACAT', 40, 0, 100)], 'summary': {'sampled_prefix': False}}
        for changes, expected in [
            ({'analysis_complete': False, 'classification': 'analysis_failed', 'legacy_result': '분석 실패'}, 'analysis_incomplete'),
            ({'adequacy_status': 'observed_qc_failure', 'classification': 'detection_inconclusive_technical', 'legacy_result': 'QC 보류'}, 'qc_limited'),
            ({'classification': 'library_context_limited', 'legacy_result': 'library 보류'}, 'complete_input_measurement'),
            ({'classification': 'taxonomic_scope_unresolved', 'legacy_result': '출처 보류'}, 'complete_input_measurement'),
        ]:
            with self.subTest(changes=changes):
                report = self.report(raw, [candidate], **changes)
                self.assertEqual(report['repeat_observation']['measurement_status'], expected)
                self.assertEqual(report['decision_summary']['display_label'], changes['legacy_result'])
                self.assertEqual(report['repeat_observation']['positive_family_count'], 1)
        raw['summary']['sampled_prefix'] = True
        report = self.report(raw, [candidate], analysis_complete=False, classification='analysis_failed', legacy_result='분석 실패')
        self.assertEqual(report['repeat_observation']['measurement_status'], 'partial_input')

    def test_reference_scope_limit_remains_headline_when_policy_selected_known_repeat(self):
        candidate = self.candidate('AACCT', 40, 100)
        report = self.report({}, [candidate], candidate, classification='known_trm_detected',
                             legacy_result='알려진 TRM 패턴 검출',
                             agreement={'status': 'library_context_limited', 'main_trm_status': 'not_confirmed'})
        self.assertEqual(report['decision_summary']['classification'], 'known_trm_detected')
        self.assertIn('라이브러리 제약', report['decision_summary']['display_label'])
        self.assertEqual(report['policy_selection']['representative_motif'], 'AACCT')

    def test_unknown_original_input_scope_is_not_labeled_full_input(self):
        for summary in ({}, {'sampled_prefix': None}):
            report = self.report({'summary': summary}, [])
            self.assertEqual(report['repeat_observation']['measurement_status'], 'input_completeness_unverified')
        full = self.report({'summary': {'sampled_prefix': False}}, [])
        self.assertEqual(full['repeat_observation']['measurement_status'], 'complete_input_measurement')

    def test_disagreement_preview_is_bounded_but_count_preserved(self):
        candidates = [self.candidate(motif, 0, 0, original=40)
                      for motif in ('ACAT', 'AACCT', 'AACCC', 'AAACCC', 'AAACCCT', 'AACCCTG')]
        report = self.report({}, candidates)
        view = report['repeat_observation']
        self.assertEqual(view['original_positive_current_nonpositive_count'], 6)
        self.assertEqual(len(view['original_positive_current_nonpositive']), 5)
        self.assertFalse(view['original_positive_current_nonpositive_preview_complete'])
        self.assertEqual(len(candidates), 6)


if __name__ == '__main__':
    unittest.main()
