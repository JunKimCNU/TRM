"""Additive joint direction/context diagnostics; no policy changes."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import unittest

from archive_test_support import archive_only

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import trm_assessment as assessment
import report_assessment as reporting

KEYS = ('left_forward', 'left_reverse', 'right_forward', 'right_reverse', 'ambiguous')
FAMILY = 'AACACTCTACACAGTGAGTAG'


def entry(sequence, **raw):
    raw = dict.fromkeys(KEYS, 0) | raw
    normalized = dict(left_forward=raw['left_forward']+raw['right_reverse'],
        left_reverse=raw['left_reverse']+raw['right_forward'], right_forward=0, right_reverse=0,
        ambiguous=raw['ambiguous'])
    return {'signature': sequence, 'sequence': sequence, 'representative_flank': sequence,
        'events': sum(raw.values()), 'quadrant_counts': raw, 'normalized_quadrant_counts': normalized}


def flank_sketch(entries):
    return {'capacity': 128, 'flank_bases': 20, 'observed_events': sum(e['events'] for e in entries),
        'distinct_lower_bound': len(entries), 'saturated': False, 'entries': entries}


def geometry(family=FAMILY, conventional=None):
    return assessment.directional_evidence(family, conventional, {'junction_reads': 23,
        'left_reverse_junction_reads': 10, 'right_forward_junction_reads': 12,
        'right_reverse_junction_reads': 1})


class DirectionalClusterTests(unittest.TestCase):
    def setUp(self):
        self.flanks = flank_sketch([entry('ACGTACGTACGTACGTACGT', left_reverse=10, right_forward=12),
                                    entry('TGCATGCATGCATGCATGCA', right_reverse=1)])

    def diagnose(self, flanks=None, chosen=None, family=FAMILY):
        return assessment.directional_flank_cluster_diagnostics(family, chosen or geometry(family), self.flanks if flanks is None else flanks)

    def test_joint_clusters_preserve_fixed_axis_and_concentration(self):
        d = self.diagnose()
        self.assertEqual(d['status'], 'measured')
        self.assertEqual(d['all_observed_clusters_exact'], 2)
        self.assertEqual(d['retained_clusters_with_consistent_events'], 1)
        self.assertEqual(d['retained_clusters_with_opposed_events'], 1)
        self.assertEqual(d['retained_consistent_events'], 22)
        self.assertEqual(d['largest_retained_consistent_cluster_event_fraction'], 1)
        self.assertEqual(d['retained_clusters_consistent_at_both_original_read_edges'], 1)
        self.assertFalse(d['changes_classification'])
        self.assertEqual(d['orientation_interpretation'], 'inferred_from_same_junction_events')

    def test_opposed_majority_in_one_cluster_not_reselected(self):
        d = self.diagnose(flank_sketch([entry('ACGTACGTACGTACGTACGT', left_reverse=1, left_forward=20)]))
        self.assertEqual(d['retained_consistent_events'], 1)
        self.assertEqual(d['retained_opposed_events'], 20)
        self.assertEqual(d['retained_clusters_with_both_directions'], 1)
        self.assertEqual(d['retained_consistent_only_clusters'], 0)

    def test_saturated_subset_is_not_population_bound(self):
        self.flanks.update(capacity=2, saturated=True, observed_events=24, distinct_lower_bound=3)
        d = self.diagnose()
        self.assertTrue(d['flank_sketch_saturated'])
        self.assertEqual(d['flank_events_not_retained'], 1)
        self.assertEqual(d['retained_cluster_count'], 2)
        self.assertEqual(d['largest_retained_consistent_cluster_event_fraction'], 1)
        for key in ('all_observed_clusters_exact', 'all_observed_clusters_with_consistent_events_exact',
                    'all_observed_clusters_with_opposed_events_exact', 'all_observed_largest_consistent_cluster_event_fraction'):
            self.assertIsNone(d[key])
        self.assertIn('not_a_population_bound', d['cluster_scope'])
        columns = reporting.directional_cluster_columns({'directional_flank_clusters': d})
        self.assertIn('보관 1 (포화', columns['consistent_flank_clusters'])
        self.assertIn('보관 22/22', columns['retained_consistent_cluster_concentration'])
        self.assertIn('23/24', columns['flank_context_completeness'])

    def test_missing_and_invalid_saturation_are_not_false(self):
        for value in (None, 0, 'False'):
            malformed = copy.deepcopy(self.flanks)
            malformed['saturated'] = value
            d = self.diagnose(malformed)
            self.assertEqual(d['status'], 'incomplete')
            self.assertIsNone(d['flank_sketch_saturated'])
            self.assertIsNone(d['retained_consistent_events'])
        del self.flanks['saturated']
        self.assertEqual(self.diagnose()['status'], 'incomplete')

    def test_missing_required_counts_and_normalized_quadrants_are_unknown(self):
        for key in ('observed_events', 'distinct_lower_bound', 'capacity', 'flank_bases'):
            malformed = copy.deepcopy(self.flanks)
            del malformed[key]
            self.assertEqual(self.diagnose(malformed)['status'], 'incomplete')
        for key in KEYS:
            malformed = copy.deepcopy(self.flanks)
            del malformed['entries'][0]['normalized_quadrant_counts'][key]
            d = self.diagnose(malformed)
            self.assertEqual(d['status'], 'incomplete')
            self.assertIsNone(d['retained_clusters_with_consistent_events'])

    def test_unsaturated_omitted_events_and_mismatching_quadrants_are_incomplete(self):
        malformed = copy.deepcopy(self.flanks)
        malformed['observed_events'] += 1
        self.assertEqual(self.diagnose(malformed)['status'], 'incomplete')
        malformed = copy.deepcopy(self.flanks)
        malformed['entries'][0]['normalized_quadrant_counts']['left_reverse'] -= 1
        self.assertEqual(self.diagnose(malformed)['status'], 'incomplete')
        malformed = copy.deepcopy(self.flanks)
        malformed['entries'][0]['quadrant_counts'].update(left_reverse=0, left_forward=10)
        d = self.diagnose(malformed)
        self.assertEqual(d['status'], 'incomplete')
        self.assertEqual(d['reason'], 'raw_and_outward_normalized_flank_quadrants_disagree')

    def test_legacy_without_raw_quadrants_keeps_direction_but_unknown_read_edges(self):
        for e in self.flanks['entries']:
            del e['quadrant_counts']
        d = self.diagnose()
        self.assertEqual(d['status'], 'measured')
        self.assertFalse(d['raw_read_edge_counts_available'])
        self.assertEqual(d['retained_consistent_events'], 22)
        self.assertIsNone(d['retained_clusters_consistent_at_both_original_read_edges'])
        self.assertTrue(all(c['consistent_left_read_edge_events'] is None for c in d['clusters_retained']))

    def test_reverse_complement_swaps_raw_edges_not_normalized_direction(self):
        original = self.diagnose()
        reverse_key = dict(left_forward='right_reverse', left_reverse='right_forward',
                           right_forward='left_reverse', right_reverse='left_forward', ambiguous='ambiguous')
        for e in self.flanks['entries']:
            q = e['quadrant_counts']
            e['quadrant_counts'] = {key: q[reverse_key[key]] for key in KEYS}
        d = self.diagnose()
        self.assertEqual(d['retained_consistent_events'], original['retained_consistent_events'])
        for before, after in zip(original['clusters_retained'], d['clusters_retained']):
            self.assertEqual(before['consistent_left_read_edge_events'], after['consistent_right_read_edge_events'])

    def test_measured_empty_and_missing_are_different(self):
        missing = assessment.directional_flank_cluster_diagnostics(FAMILY, geometry())
        self.assertEqual(missing['status'], 'unmeasured')
        self.assertIsNone(missing['retained_consistent_events'])
        d = self.diagnose(flank_sketch([]), assessment.directional_evidence(FAMILY, None, {'junction_reads': 0}))
        self.assertEqual(d['status'], 'measured')
        self.assertEqual(d['retained_consistent_events'], 0)
        self.assertEqual(d['all_observed_clusters_exact'], 0)
        self.assertIsNone(d['largest_retained_consistent_cluster_event_fraction'])
        self.assertEqual(d['orientation_interpretation'], 'not_established_no_oriented_events')
        columns = reporting.directional_cluster_columns({'directional_flank_clusters': d})
        self.assertIn('비율 없음', columns['retained_consistent_cluster_concentration'])
        self.assertEqual(reporting.validation_columns(None)['consistent_flank_clusters'], '미측정')

    def test_ambiguous_only_has_zero_consistent_count_but_no_fraction(self):
        d = self.diagnose(flank_sketch([entry('ACGTACGTACGTACGTACGT', ambiguous=3)]))
        self.assertEqual(d['retained_clusters_with_ambiguous_events'], 1)
        self.assertEqual(d['retained_consistent_events'], 0)
        self.assertIsNone(d['largest_retained_consistent_cluster_event_fraction'])

    def test_palindromic_family_does_not_invent_a_direction(self):
        d = self.diagnose(family='AT')
        self.assertEqual(d['status'], 'strand_ambiguous')
        self.assertIsNone(d['retained_clusters_with_consistent_events'])


class FrozenPrototypeRegressionTests(unittest.TestCase):
    def assert_readiness_preserved(self, before, after):
        # The optional quality collector changes this explanatory sentence,
        # not any readiness rule, count, missingness or other limitation.
        old_note = 'Base qualities are not assessed by this evidence collector; flank clustering is not a quality-aware error model.'
        new_note = 'Flank clustering is not a quality-aware error model; optional base-quality histograms are reported separately.'
        self.assertIn(old_note, before['limitations'])
        self.assertIn(new_note, after['limitations'])
        expected = copy.deepcopy(before)
        expected['limitations'] = [new_note if item == old_note else item for item in expected['limitations']]
        self.assertEqual(expected, after)

    @unittest.skipUnless((ROOT/'results/fungal_mixed_wgs_20261004/directional_context/report.json').is_file(), 'Local frozen study not distributed')
    @archive_only
    def test_six_families_match_prototype_and_old_readiness(self):
        prototype = json.loads((ROOT/'results/fungal_mixed_wgs_20261004/directional_context/report.json').read_text())
        base = ROOT/'results/overnight_validation_20261004'
        for row in prototype['families']:
            prior = row['summary']
            native = json.loads((base/f'native/{prior["sample_id"]}/result.json').read_text())
            old = json.loads((base/f'assessment/samples/{prior["sample_id"]}.json').read_text())
            raw = next(c for pool in ('candidates', 'rejected_candidates') for c in native[pool] if c['motif'] == prior['canonical_family'])
            candidate = next(c for c in old['candidates'] if c['canonical_family'] == prior['canonical_family'])
            readiness = assessment.confirmation_readiness(prior['canonical_family'], candidate['directionality'], raw['junction_evidence'], True)
            new = readiness.pop('directional_flank_clusters')
            self.assert_readiness_preserved(candidate['confirmation_readiness'], readiness)
            self.assertEqual(new['status'], 'measured')
            for key in new:
                if key in prior:
                    self.assertEqual(new[key], prior[key], (prior['sample_id'], prior['canonical_family'], key))
            for actual, expected in zip(new['clusters_retained'], row['clusters']):
                self.assertEqual(actual['representative'], expected['representative_anchor'])
                self.assertEqual(actual['retained_contexts'], expected['retained_anchor_count'])
                for key in actual.keys() & expected.keys():
                    self.assertEqual(actual[key], expected[key])

    @unittest.skipUnless((ROOT/'results/evidence_capture_reporting_20261004/before/scripts/trm_assessment.py').is_file(), 'Local baseline snapshot not distributed')
    @archive_only
    def test_all_740_saved_candidates_keep_every_existing_readiness_field(self):
        spec = importlib.util.spec_from_file_location('baseline_cluster_assessment', ROOT/'results/evidence_capture_reporting_20261004/before/scripts/trm_assessment.py')
        previous = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(previous)
        # Policy release metadata can change independently of the readiness
        # thresholds checked here. Preserve exact comparison of every rule.
        release_metadata = {'version', 'patch_scope'}
        self.assertEqual(
            {k: v for k, v in assessment.POLICY.items() if k not in release_metadata},
            {k: v for k, v in previous.POLICY.items() if k not in release_metadata},
        )
        base = ROOT/'results/overnight_validation_20261004'
        count = 0
        for run in ('SRR12816473', 'SRR18493388'):
            native = json.loads((base/f'native/{run}/result.json').read_text())
            old = json.loads((base/f'assessment/samples/{run}.json').read_text())
            candidates = {c['canonical_family']: c for c in old['candidates']}
            for pool in ('candidates', 'rejected_candidates'):
                for raw in native[pool]:
                    chosen = candidates[raw['motif']]['directionality']
                    before = previous.confirmation_readiness(raw['motif'], chosen, raw['junction_evidence'], True)
                    after = assessment.confirmation_readiness(raw['motif'], chosen, raw['junction_evidence'], True)
                    after.pop('directional_flank_clusters')
                    self.assert_readiness_preserved(before, after)
                    count += 1
        self.assertEqual(count, 740)


if __name__ == '__main__':
    unittest.main()
