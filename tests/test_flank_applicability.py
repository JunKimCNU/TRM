"""Necessary-length bounds stay separate from motif observations and power."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from input_applicability import build_input_applicability

ZERO = 'no_flank_length_eligible_reads'
PRESENT = 'flank_length_eligible_reads_present_evaluation_unknown'


def report(length=151, k=27, **config):
    return dict(program='trm', version='0.8.4', analysis_kind='de_novo_discovery',
                config=dict(mode='reads', min_k=k, max_k=k, seed_bases=24, min_copies=3, **config),
                junction_options=dict(flank_bases=20, end_slack=3, flank_period_annotation=True),
                summary=dict(records=10, bases=10*length, max_record_bases=length, sampled_prefix=False),
                evidence_qc=dict(records=10, bases=10*length, max_record_bases=length,
                                 length_histogram=[dict(lower=length, upper=length, records=10)]),
                candidates=[dict(motif='AAGGAATGGGAGCAGCGGGAGCAGGCG', reads=134)],
                data_adequacy=dict(status='unassessable_library_or_genome'))


def diagnose(document, motifs=(), **kwargs):
    return build_input_applicability(document, motifs, **kwargs)['flank_period_applicability']


class FlankApplicabilityTests(unittest.TestCase):
    def test_supported_release_boundaries_in_discovery_and_targeted_modes(self):
        # Independent expected constants cover seed rounding, copy minimum,
        # candidate-period window, explicit flank length, and the 90% boundary.
        cases = [(6, 24, 3, 20, 44), (7, 24, 3, 20, 49),
                 (27, 24, 3, 20, 162), (6, 24, 5, 20, 50),
                 (6, 24, 3, 40, 64), (3, 900, 3, 20, 1001)]
        for version in ('0.8.4', '0.8.5'):
            for kind in ('de_novo_discovery', 'targeted_evidence'):
                for enabled in (True, False, None):
                    for k, seed, copies, flank, threshold in cases:
                        for length, expected in ((threshold - 1, ZERO), (threshold, PRESENT)):
                            with self.subTest(version=version, kind=kind, enabled=enabled,
                                              k=k, seed=seed, copies=copies, flank=flank,
                                              length=length):
                                raw = report(length, k)
                                raw.update(version=version, analysis_kind=kind)
                                raw['config'].update(seed_bases=seed, min_copies=copies)
                                raw['junction_options']['flank_bases'] = flank
                                if enabled is None:
                                    raw['junction_options'].pop('flank_period_annotation')
                                else:
                                    raw['junction_options']['flank_period_annotation'] = enabled
                                before = copy.deepcopy(raw)
                                result = diagnose(raw)
                                self.assertEqual(raw, before)
                                self.assertEqual(result['status'], expected)
                                self.assertEqual(result['periods'][0]['minimum_read_bp'], threshold)
                                self.assertEqual(result['model_status'], f'supported_native_v{version}')
                                self.assertEqual(result['model'],
                                                 f'native_v{version}_disjoint_seed_and_candidate_period_window')
                                self.assertEqual(result['analysis_kind'], kind)
                                self.assertEqual(result['annotation_enabled'], bool(enabled))
                                self.assertFalse(result['actual_event_measurement_inferred'])

    def test_v085_diagnostics_preserve_other_input_fields_and_evidence(self):
        raw = report(162)
        raw['candidates'][0]['junction_evidence'] = {
            'flank_period_annotation': {'observed_events': 0, 'measured_all_events': True}}
        old = build_input_applicability(raw, ['TTAGGG'])
        raw['version'] = '0.8.5'
        before = copy.deepcopy(raw)
        current = build_input_applicability(raw, ['TTAGGG'])
        self.assertEqual(raw, before)
        old_flank = old.pop('flank_period_applicability')
        new_flank = current.pop('flank_period_applicability')
        self.assertEqual(old, current)
        for name in ('native_version', 'model_status', 'model'):
            old_flank.pop(name)
            new_flank.pop(name)
        self.assertEqual(old_flank, new_flank)

    def test_future_custom_and_malformed_versions_remain_unavailable(self):
        for version in ('0.8.3', '0.8.6', '0.9.0', '0.8.4-custom', '0.8.5-custom',
                        'v0.8.5', '', None, True, 0.85, [], {}):
            with self.subTest(version=version):
                raw = report(162)
                raw['version'] = version
                result = diagnose(raw)
                self.assertEqual(result['status'], 'flank_model_unavailable')
                self.assertEqual(result['model_status'], 'unavailable')
                self.assertIsNone(result['periods'][0]['minimum_read_bp'])
                self.assertIsNone(result['annotation_enabled'])

    def test_suillus_detection_possible_but_flank_window_impossible(self):
        raw = report()
        saved = copy.deepcopy(raw)
        result = build_input_applicability(raw)
        self.assertEqual(raw, saved)
        self.assertEqual(result['periods'][0]['length_eligible_reads_lower'], 10)
        flank = result['flank_period_applicability']
        row = flank['periods'][0]
        self.assertEqual(flank['status'], ZERO)
        self.assertEqual((row['minimum_seed_bp'], row['required_window_bp'], row['minimum_read_bp']), (81,81,162))
        self.assertEqual((row['length_eligible_reads_lower'],row['length_eligible_reads_upper']), (0,0))
        self.assertFalse(flank['actual_event_measurement_inferred'])

    def test_161_162_boundary_and_mates_not_summed(self):
        for length, expected in ((151, ZERO), (161, ZERO), (162, PRESENT)):
            with self.subTest(length=length):
                raw = report(length)
                raw['paired_end'] = dict(enabled=True, pairs=5)
                self.assertEqual(diagnose(raw)['status'], expected)

    def test_actual_junction_flank_setting_can_dominate(self):
        raw = report(63, 6)
        raw['junction_options']['flank_bases'] = 40
        result = diagnose(raw)
        self.assertEqual(result['status'], ZERO)
        self.assertEqual(result['periods'][0]['minimum_read_bp'], 64)
        self.assertEqual(result['periods'][0]['required_window_bp'], 20)

    def test_full_repeat_exclusion_dominates_large_seed(self):
        for length, expected in ((1000, ZERO), (1001, PRESENT)):
            raw = report(length, 3)
            raw['config']['seed_bases'] = 900
            row = diagnose(raw)['periods'][0]
            self.assertEqual(row['minimum_read_bp'], 1001)
            self.assertEqual(row['status'], expected)

    def test_min_tract_and_representative_cap_do_not_replace_window(self):
        raw = report()
        raw['config']['min_tract'] = 500
        first = diagnose(raw)
        raw['config']['min_tract'] = 5000
        self.assertEqual(diagnose(raw), first)
        self.assertEqual(first['periods'][0]['required_window_bp'],81)

    def test_off_and_missing_known_false_do_not_change_geometry(self):
        raw = report(162)
        raw['junction_options']['flank_period_annotation'] = False
        off = diagnose(raw)
        del raw['junction_options']['flank_period_annotation']
        absent_false = diagnose(raw)
        self.assertEqual(off, absent_false)
        self.assertFalse(off['annotation_enabled'])
        self.assertEqual(off['status'], PRESENT)
        self.assertFalse(off['actual_event_measurement_inferred'])

    def test_unknown_native_or_missing_options_do_not_get_defaults(self):
        for change in ('version', 'program', 'junction_options', 'flank_bases'):
            raw = report()
            if change == 'flank_bases':
                del raw['junction_options']['flank_bases']
            else:
                del raw[change]
            with self.subTest(change=change):
                result = diagnose(raw)
                self.assertEqual(result['scope'], 'unavailable')
                self.assertIsNone(result['periods'][0]['minimum_read_bp'])
                self.assertIsNone(result['periods'][0]['length_eligible_reads_upper'])
        raw = report(); raw['version'] = '0.8.4-custom'
        self.assertEqual(diagnose(raw)['status'], 'flank_model_unavailable')

    def test_invalid_options_bool_float_and_out_of_native_range(self):
        for field, value in [('flank_bases', True), ('flank_bases',20.0), ('flank_bases',19),
                             ('flank_bases',10001), ('end_slack',True), ('end_slack',1001),
                             ('flank_period_annotation','on'), ('flank_period_annotation',None)]:
            with self.subTest(field=field,value=value):
                raw=report();raw['junction_options'][field]=value
                self.assertEqual(diagnose(raw)['status'], 'junction_options_invalid')

    def test_prefix_unverified_and_failed_results_keep_analyzed_scope(self):
        for mode in ('prefix','unspecified','failed','caller_incomplete'):
            raw=report()
            kwargs={}
            if mode=='prefix':raw['summary']['sampled_prefix']=True
            elif mode=='unspecified':del raw['summary']['sampled_prefix']
            elif mode=='failed':raw['status']='failed'
            else:kwargs['analysis_complete']=False
            result=diagnose(raw,**kwargs)
            self.assertNotEqual(result['scope'],'complete_input')
            self.assertEqual(result['periods'][0]['status'],result['status'])
            self.assertEqual(result['periods'][0]['analyzed_records_length_status'],ZERO)

    def test_invalid_qc_and_source_never_become_zero(self):
        raw=report();raw['evidence_qc']['records']=9
        result=diagnose(raw)
        self.assertEqual(result['status'],'qc_invalid')
        self.assertIsNone(result['periods'][0]['length_eligible_reads_upper'])
        for document,kwargs in [(report(),dict(source_available=False)),(None,{})]:
            result=diagnose(document,['TTAGGG'],**kwargs)
            self.assertEqual(result['status'],'source_unavailable')
            self.assertIsNone(result['families'][0]['minimum_read_bp'])

    def test_coarse_histogram_does_not_invent_exact_counts(self):
        raw=report(162)
        raw['summary'].update(records=10,bases=1610,max_record_bases=170)
        raw['evidence_qc']=dict(length_histogram=[dict(lower=150,upper=170,records=10)])
        result=diagnose(raw)
        row=result['periods'][0]
        self.assertEqual(result['status'],'flank_length_eligibility_unresolved')
        self.assertEqual(row['length_eligible_reads_lower'],0)
        self.assertGreater(row['length_eligible_reads_upper'],0)

    def test_expected_family_alias_and_range_are_separate_from_search(self):
        raw=report(151,30)
        result=diagnose(raw,['TTAGGG','CCCTAA','TTAGGGTTAGGG'])
        self.assertEqual(result['status'],ZERO)
        self.assertEqual(len(result['families']),1)
        row=result['families'][0]
        self.assertEqual(row['minimum_read_bp'],44)
        self.assertEqual(row['status'],PRESENT)
        self.assertTrue(row['outside_discovery_range'])

    def test_assembly_mode_unavailable_and_empty_reads_measured_zero_length(self):
        raw=report();raw['config']['mode']='assembly'
        self.assertEqual(diagnose(raw)['status'],'mode_not_supported')
        raw=report();raw['summary'].update(records=0,bases=0,max_record_bases=0)
        raw.pop('evidence_qc')
        self.assertEqual(diagnose(raw)['periods'][0]['length_eligible_reads_upper'],0)


if __name__=='__main__':
    unittest.main()
