"""Regression probes for conservative ambiguity and target-specific power scope."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
from assess_adequacy import REQUIRED_ASSUMPTIONS, assess_adequacy
from evidence_batch import canonical
from trm_assessment import assess_candidate, classify, expected_family_power_scope


def native_report():
    return {'config': {'mode': 'reads', 'min_k': 3, 'max_k': 32, 'min_copies': 3,
                       'seed_bases': 24, 'min_reads': 2},
            'summary': {'records': 10_000_000, 'bases': 1_500_000_000,
                        'max_record_bases': 150, 'sampled_prefix': False},
            'evidence_qc': {'records_with_ambiguity': 0, 'ambiguous_bases': 0}}


def assumptions():
    return {'genome_size_bp': 100_000_000, 'haploid_chromosome_count': 6,
            'random_wgs_assumed': True, 'telomere_retention_fraction': 0.1,
            'effective_independent_read_fraction': 0.5, 'detector_sensitivity': 0.5,
            'provenance': {key: 'Explicit synthetic hypothesis' for key in REQUIRED_ASSUMPTIONS}}


def power(report=None):
    return assess_adequacy(native_report() if report is None else report,
                           assumptions=assumptions(), include_illustrative_grid=False)


def assignment(motifs='TTAGGG'):
    return {'reference_status': 'supported', 'reference_scope': 'clade', 'related_motifs': motifs}


def known_candidate():
    raw = {'motif': 'TTAGGG', 'reads': 100, 'full_repeat_reads': 30,
           'repeat_bases': 12000, 'max_tract': 150}
    refs = {canonical('TTAGGG'): [{'motif': 'TTAGGG', 'evidence': 'experimental',
                                 'source': 'Synthetic reference'}]}
    return assess_candidate(raw, None, refs, {canonical('TTAGGG')}, 150)


class AmbiguityGuardTests(unittest.TestCase):
    def test_all_ambiguous_bases_cannot_become_clean_reads_when_record_counter_missing(self):
        report = native_report()
        report['evidence_qc'] = {'ambiguous_bases': report['summary']['bases']}
        result = power(report)
        self.assertEqual(result['status'], 'insufficient_data')
        self.assertEqual(result['read_qc']['ambiguity_record_bounds']['status'], 'bounded_from_ambiguous_bases')
        self.assertTrue(all(row['ambiguity_free_eligible_reads_lower_bound'] == 0 for row in result['power']['periods']))
        self.assertEqual(result['power']['worst_case_miss_probability_upper_bound'], 1)

    def test_missing_both_counters_withholds_power_claim_and_reports_missingness(self):
        report = native_report()
        report.pop('evidence_qc')
        result = power(report)
        self.assertEqual(result['status'], 'unassessable_library_or_genome')
        self.assertIn('ambiguity_qc_unavailable', result['reasons'])
        self.assertFalse(result['read_qc']['ambiguity_record_bounds']['bound_available'])
        self.assertEqual(result['power']['periods'][0]['effective_independent_trials_assumed_lower_bound'], 0)

    def test_measured_zero_ambiguous_bases_proves_zero_ambiguous_records(self):
        report = native_report()
        report['evidence_qc'] = {'ambiguous_bases': 0}
        result = power(report)
        self.assertEqual(result['status'], 'conditional_power_supported')
        bounds = result['read_qc']['ambiguity_record_bounds']
        self.assertEqual(bounds['upper'], 0)
        self.assertFalse(bounds['record_count_measured'])

    def test_ambiguous_bases_bound_subtracts_possible_records_not_average_coverage(self):
        report = native_report()
        report['evidence_qc'] = {'ambiguous_bases': 3}
        result = power(report)
        self.assertEqual(result['power']['periods'][0]['ambiguity_free_eligible_reads_lower_bound'], 9_999_997)

    def test_explicit_record_count_remains_authoritative_without_base_counter(self):
        report = native_report()
        report['evidence_qc'] = {'records_with_ambiguity': 10_000_000}
        result = power(report)
        self.assertEqual(result['status'], 'insufficient_data')
        self.assertEqual(result['read_qc']['ambiguity_record_bounds']['status'], 'measured_record_count')

    def test_missing_ambiguity_qc_does_not_discard_positive_known_repeat(self):
        report = native_report()
        report.pop('evidence_qc')
        result = power(report)
        self.assertEqual(classify([known_candidate()], result)[0], 'known_repeat_detected_boundary_unresolved')

    def test_assessment_preserves_input_and_never_claims_absence(self):
        report = native_report()
        original = copy.deepcopy(report)
        result = power(report)
        self.assertEqual(report, original)
        self.assertFalse(result['true_absence_inferred'])
        self.assertFalse(result['alternative_maintenance_inferred'])


class ExpectedFamilyPowerScopeTests(unittest.TestCase):
    def test_unsearched_six_base_expectation_does_not_inherit_three_base_power(self):
        report = native_report()
        report['config']['max_k'] = 3
        result = power(report)
        self.assertEqual(result['status'], 'conditional_power_supported')
        scoped = expected_family_power_scope(result, assignment())
        self.assertEqual(scoped['evaluated_periods'], [3])
        self.assertEqual(scoped['expected_families'][0]['status'], 'outside_original_scan_periods')
        decision = classify([], result, reference_assignment=assignment())
        self.assertEqual(decision[0], 'detection_inconclusive_search_scope')
        self.assertEqual(decision[2], 'resolve_expected_family_search_scope')

    def test_expected_primitive_rotations_reverse_complements_and_powers_share_scope(self):
        result = power()
        scoped = expected_family_power_scope(result, assignment('TTAGGG,TAGGGT,CCCTAA,TTAGGGTTAGGG'))
        self.assertEqual(scoped['status'], 'expected_families_covered')
        self.assertEqual(len(scoped['expected_families']), 1)
        self.assertEqual(scoped['expected_families'][0]['primitive_period'], 6)
        self.assertEqual(classify([], result, reference_assignment=assignment())[0], 'canonical_not_detected_conditional_power')

    def test_any_unmodelled_expected_family_prevents_blanket_conditional_nondetection(self):
        report = native_report()
        report['config']['max_k'] = 6
        decision = classify([], power(report), reference_assignment=assignment('TTAGGG,TTTAGGG'))
        self.assertEqual(decision[0], 'detection_inconclusive_search_scope')

    def test_scanned_period_above_32_has_no_power_model_even_when_scan_extends_to_64(self):
        report = native_report()
        report['config']['max_k'] = 64
        expected = assignment('A' * 33 + 'C')
        scoped = expected_family_power_scope(power(report), expected)
        self.assertTrue(scoped['expected_families'][0]['within_original_scan_periods'])
        self.assertFalse(scoped['expected_families'][0]['within_evaluated_power_periods'])
        self.assertEqual(classify([], power(report), reference_assignment=expected)[0], 'detection_inconclusive_search_scope')

    def test_positive_known_measurement_outside_original_search_is_still_preserved(self):
        report = native_report()
        report['config']['max_k'] = 3
        self.assertEqual(classify([known_candidate()], power(report), reference_assignment=assignment())[0],
                         'known_repeat_detected_boundary_unresolved')

    def test_no_assigned_family_has_scanned_period_label_not_general_canonical_absence(self):
        decision = classify([], power())
        self.assertEqual(decision[0], 'canonical_not_detected_conditional_power')
        self.assertIn('검사 주기', decision[1])
        self.assertEqual(decision[2], 'review_scanned_period_nondetection')

    def test_unsupported_scope_cannot_be_fabricated_from_status_only(self):
        for expected in (None, assignment()):
            with self.subTest(expected=expected):
                decision = classify([], {'status': 'conditional_power_supported'}, reference_assignment=expected)
                self.assertEqual(decision[0], 'detection_inconclusive_search_scope')

    def test_invalid_reference_unit_withholds_scope_instead_of_becoming_a_period(self):
        scoped = expected_family_power_scope(power(), assignment('TTANNN'))
        self.assertEqual(scoped['status'], 'invalid_expected_units')
        self.assertFalse(scoped['conditional_nondetection_scope_supported'])

    def test_inconsistent_scan_and_power_ranges_cannot_cover_expected_family(self):
        result = power()
        result['detector_scope']['max_k'] = 3
        self.assertFalse(expected_family_power_scope(result, assignment())['conditional_nondetection_scope_supported'])


if __name__ == '__main__':
    unittest.main()
