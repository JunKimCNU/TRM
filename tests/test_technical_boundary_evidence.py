"""Bounded evidence ingestion: identity, integrity and conservative interpretation."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import technical_boundary_evidence as technical


class TechnicalBoundaryEvidenceTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / 'saved_diagnostic.json'
        self.source.write_text('{"complete":true}\n')
        self.path = self.root / 'evidence.json'
        self.native_sha = hashlib.sha256(b'immutable native result').hexdigest()
        self.report = {
            'sample_id': 'sampleA', 'scientific_name': 'Example species', 'tax_id': '123',
            'metadata': {'sample_accession': 'SAMN000001'},
            'provenance': {'original_discovery_sha256': self.native_sha},
            'candidates': [{'canonical_family': 'ACACGGATACGTAT', 'reads': 900}],
            'classification': 'uncatalogued_repeat',
        }
        self.document = {
            'schema_version': 1, 'scope': 'selected_junction_subset', 'thresholds_bp': [12, 16, 20],
            'source_artifacts': [{'path': self.source.name,
                                  'sha256': hashlib.sha256(self.source.read_bytes()).hexdigest()}],
            'reference_panel': [{'reference_id': 'read1', 'length_bp': 33},
                                {'reference_id': 'read2', 'length_bp': 33},
                                {'reference_id': 'other', 'length_bp': 58},
                                {'reference_id': 'short', 'length_bp': 19}],
            'samples': [{
                'sample_id': 'sampleA', 'source_result_sha256': self.native_sha,
                'scientific_name': 'Example species', 'tax_id': 123, 'biosample': 'SAMN000001',
                'families': ['ACACGGATACGTAT'],
                'availability': {'full_reads': True, 'mates': True},
                'selected_event_count': 2, 'selected_pair_count': 2, 'full_read_count': 4,
                'tiers': [
                    {'threshold_bp': 12, 'exact_match_associations': 6, 'full_reads_with_matches': 2,
                     'pairs_with_matches': 1, 'events_with_flank_overlap': 1},
                    {'threshold_bp': 16, 'exact_match_associations': 3, 'full_reads_with_matches': 2,
                     'pairs_with_matches': 1, 'events_with_flank_overlap': 1},
                    {'threshold_bp': 20, 'exact_match_associations': 0, 'full_reads_with_matches': 0,
                     'pairs_with_matches': 0, 'events_with_flank_overlap': 0}],
                'events': [
                    {'event_id': 'boundary1', 'original_pair_index_1based': 6075180, 'mate': 2,
                     'family': 'ACACGGATACGTAT', 'flank_bp': 21,
                     'maximum_exact_match_bp_overlapping_flank': 19,
                     'union_flank_overlap_by_threshold': {'12': 19, '16': 19, '20': 0},
                     'adjacent_unmatched_flank_by_threshold': {'12': 2, '16': 2, '20': 21},
                     'paired_readthrough_compatible': True},
                    {'event_id': 'boundary2', 'original_pair_index_1based': 6075181, 'mate': 1,
                     'family': 'ACACGGATACGTAT', 'flank_bp': 31,
                     'maximum_exact_match_bp_overlapping_flank': 0,
                     'union_flank_overlap_by_threshold': {'12': 0, '16': 0, '20': 0},
                     'adjacent_unmatched_flank_by_threshold': {'12': 31, '16': 31, '20': 31},
                     'paired_readthrough_compatible': False}],
            }],
        }

    @property
    def sample(self):
        return self.document['samples'][0]

    def write(self):
        self.document['payload_sha256'] = technical.evidence_payload_sha256(self.document)
        self.path.write_text(json.dumps(self.document))
        return self.path

    def load(self):
        return technical.load_evidence(self.write())

    def build(self, report=None, sample_id='sampleA'):
        return technical.build_technical_boundary_evidence(
            sample_id, self.report if report is None else report, self.load())

    def no_hit(self):
        for event in self.sample['events']:
            event['maximum_exact_match_bp_overlapping_flank'] = 0
            event['union_flank_overlap_by_threshold'] = dict.fromkeys(('12', '16', '20'), 0)
            event['adjacent_unmatched_flank_by_threshold'] = dict.fromkeys(('12', '16', '20'), event['flank_bp'])
            event['paired_readthrough_compatible'] = False
        for tier in self.sample['tiers']:
            for field in technical._TIER_FIELDS:
                tier[field] = 0

    def test_19bp_warning_is_visible_at_both_lower_tiers_without_biological_mutation(self):
        before = deepcopy(self.report)
        result = self.build()
        self.assertEqual(self.report, before)
        self.assertEqual(result['status'], 'technical_boundary_warning')
        self.assertEqual(result['flagged_event_count'], 1)
        self.assertEqual(result['paired_readthrough_compatible_event_count'], 1)
        self.assertEqual([tier['exact_match_associations'] for tier in result['tiers']], [6, 3, 0])
        self.assertEqual(result['selected_event_count'], 2)
        self.assertEqual(result['selected_pair_count'], 2)
        self.assertEqual(result['full_read_count'], 4)
        self.assertEqual(result['events'][0]['adjacent_unmatched_flank_by_threshold'], {'12': 2, '16': 2, '20': 21})
        for field in ('native_counts_changed', 'classification_changed', 'adequacy_changed', 'sequence_trimmed'):
            self.assertIs(result[field], False)
        for field in ('technical_origin_confirmed', 'physical_terminal_confirmed', 'novel_trm_confirmed'):
            self.assertIsNone(result[field])
        self.assertEqual(result['provenance']['source_result_sha256'], self.native_sha)
        self.assertEqual(result['reference_panel'][-1]['eligible_thresholds_bp'], [12, 16])

    def test_absent_evidence_and_unlisted_sample_are_unassessed_not_zero(self):
        missing = technical.build_technical_boundary_evidence('sampleA', self.report)
        unlisted = self.build(sample_id='notSelected')
        for result in (missing, unlisted):
            self.assertEqual(result['status'], 'unassessed')
            self.assertEqual(result['events'], [])
            self.assertEqual(result['tiers'], [])
            for field in (*technical._COUNT_FIELDS, 'flagged_event_count'):
                self.assertIsNone(result[field])
        self.assertEqual(unlisted['reason'], 'sample_not_in_evidence')

    def test_missing_full_reads_or_mates_is_unassessed(self):
        for missing in ('full_reads', 'mates'):
            with self.subTest(missing=missing):
                self.sample['availability'] = {'full_reads': True, 'mates': True, missing: False}
                for field in (*technical._COUNT_FIELDS, 'events', 'tiers'):
                    self.sample[field] = None
                result = self.build()
                self.assertEqual(result['reason'], missing + '_unavailable')
                self.assertIsNone(result['flagged_event_count'])
                self.assertEqual(result['events'], [])

    def test_no_hit_only_describes_selected_subset_and_limited_panel(self):
        self.no_hit()
        result = self.build()
        self.assertEqual(result['status'], 'no_exact_flank_overlap_in_selected_subset')
        self.assertEqual(result['flagged_event_count'], 0)
        self.assertEqual(result['scope'], 'selected_junction_subset')
        self.assertTrue(any('does not establish adapter-free' in line for line in result['limitations']))
        self.assertIsNone(result['technical_origin_confirmed'])

    def test_full_read_hit_outside_flank_is_not_called_no_hit(self):
        self.no_hit()
        for tier in self.sample['tiers'][:2]:
            tier.update(exact_match_associations=1, full_reads_with_matches=1, pairs_with_matches=1)
        result = self.build()
        self.assertEqual(result['status'], 'technical_reference_match_without_flank_overlap')
        self.assertEqual(result['flagged_event_count'], 0)

    def test_one_base_flank_overlap_does_not_claim_twelve_base_flank_match(self):
        event = self.sample['events'][0]
        event['maximum_exact_match_bp_overlapping_flank'] = 20
        event['union_flank_overlap_by_threshold'] = {'12': 1, '16': 1, '20': 1}
        event['adjacent_unmatched_flank_by_threshold'] = {'12': 20, '16': 20, '20': 20}
        self.sample['tiers'][2].update(exact_match_associations=1, full_reads_with_matches=1,
                                        pairs_with_matches=1, events_with_flank_overlap=1)
        result = self.build()
        self.assertEqual(result['flagged_event_count'], 1)
        self.assertIn('구간이 비반복 flank와 겹칩니다', result['summary'])
        self.assertNotIn('flank가 기술서열과 12 bp 이상', result['summary'])

    def test_wrong_native_sample_identity_or_candidate_family_is_rejected(self):
        for field, value in (('sample_id', 'wrong'), ('scientific_name', 'Other species'),
                             ('tax_id', '124'), ('candidates', [{'canonical_family': 'TTAGGG'}])):
            with self.subTest(field=field):
                report = deepcopy(self.report)
                report[field] = value
                with self.assertRaises(technical.EvidenceError):
                    self.build(report)
        report = deepcopy(self.report)
        report['provenance']['original_discovery_sha256'] = '0' * 64
        with self.assertRaisesRegex(technical.EvidenceError, 'native result SHA256 mismatch'):
            self.build(report)
        report = deepcopy(self.report)
        report['metadata']['sample_accession'] = 'SAMN999999'
        with self.assertRaisesRegex(technical.EvidenceError, 'biosample mismatch'):
            self.build(report)

    def test_missing_native_identity_and_unvalidated_mapping_are_rejected(self):
        report = deepcopy(self.report)
        del report['provenance']['original_discovery_sha256']
        with self.assertRaises(technical.EvidenceError):
            self.build(report)
        with self.assertRaisesRegex(technical.EvidenceError, 'load_evidence'):
            technical.build_technical_boundary_evidence('sampleA', self.report, {'sampleA': self.sample})

    def test_source_change_missing_source_and_payload_change_fail_explicitly(self):
        self.write()
        original = self.source.read_bytes()
        self.source.write_bytes(b'changed')
        with self.assertRaisesRegex(technical.EvidenceError, 'source artifact SHA256 mismatch'):
            technical.load_evidence(self.path)
        self.source.unlink()
        with self.assertRaisesRegex(technical.EvidenceError, 'cannot read'):
            technical.load_evidence(self.path)
        self.source.write_bytes(original)
        edited = json.loads(self.path.read_text())
        edited['samples'][0]['scientific_name'] = 'edited without new seal'
        self.path.write_text(json.dumps(edited))
        with self.assertRaisesRegex(technical.EvidenceError, 'payload seal mismatch'):
            technical.load_evidence(self.path)

    def test_dependency_paths_keep_relative_sources_and_symlinks(self):
        link = self.root / 'source_alias.json'
        link.symlink_to(self.source)
        self.document['source_artifacts'][0]['path'] = link.name
        evidence = self.load()
        self.assertEqual(evidence.source_paths, (str(link),))
        self.assertEqual(technical.source_dependencies(self.path), [str(self.path), str(link)])
        self.assertEqual(evidence.provenance['source_artifacts'][0]['path'], str(link))

    def test_duplicate_sample_event_reference_or_json_key_is_rejected(self):
        baseline = deepcopy(self.document)
        variants = []
        changed = deepcopy(baseline)
        changed['samples'].append(deepcopy(changed['samples'][0]))
        variants.append(changed)
        changed = deepcopy(baseline)
        changed['samples'][0]['events'][1]['event_id'] = 'boundary1'
        variants.append(changed)
        changed = deepcopy(baseline)
        changed['reference_panel'].append(deepcopy(changed['reference_panel'][0]))
        variants.append(changed)
        for variant in variants:
            self.document = variant
            with self.assertRaises(technical.EvidenceError):
                self.load()
        self.path.write_text('{"schema_version":1,"schema_version":1}')
        with self.assertRaises(technical.EvidenceError):
            technical.load_evidence(self.path)

    def test_invalid_counts_and_threshold_rows_are_rejected(self):
        baseline = deepcopy(self.document)
        mutations = [
            lambda s: s.update(selected_event_count=3),
            lambda s: s.update(selected_pair_count=3, full_read_count=6),
            lambda s: s.update(full_read_count=3),
            lambda s: s.update(selected_event_count=True),
            lambda s: s['tiers'][0].update(events_with_flank_overlap=0),
            lambda s: s['tiers'][0].update(exact_match_associations=1),
            lambda s: s['tiers'][1].update(exact_match_associations=7),
            lambda s: s['tiers'][0].update(pairs_with_matches=0),
            lambda s: s['tiers'][2].update(threshold_bp=16),
            lambda s: s['events'][0]['union_flank_overlap_by_threshold'].update({'20': 1}),
            lambda s: s['events'][0]['adjacent_unmatched_flank_by_threshold'].update({'12': 3}),
            lambda s: s['events'][0].update(mate=3),
            lambda s: s['events'][0].update(family='TTAGGG'),
            lambda s: s['events'][0].update(maximum_exact_match_bp_overlapping_flank=100),
            lambda s: s['availability'].update(mates=False),
        ]
        for index, mutate in enumerate(mutations):
            with self.subTest(index=index):
                self.document = deepcopy(baseline)
                mutate(self.sample)
                with self.assertRaises(technical.EvidenceError):
                    self.load()

    def test_multiple_events_on_one_pair_do_not_inflate_pair_or_read_denominator(self):
        event = deepcopy(self.sample['events'][1])
        event.update(event_id='boundary3', original_pair_index_1based=6075180, mate=1)
        self.sample['events'].append(event)
        self.sample['selected_event_count'] = 3
        result = self.build()
        self.assertEqual((result['selected_event_count'], result['selected_pair_count'], result['full_read_count']), (3, 2, 4))
        self.assertEqual(result['flagged_event_count'], 1)

    def test_empty_selected_subset_is_unassessed(self):
        self.no_hit()
        self.sample.update(events=[], selected_event_count=0, selected_pair_count=0, full_read_count=0)
        result = self.build()
        self.assertEqual(result['reason'], 'no_selected_junction_events')
        self.assertIsNone(result['flagged_event_count'])

    def test_loaded_evidence_and_returned_diagnostic_are_not_mutable_shared_state(self):
        evidence = self.load()
        sample = evidence['sampleA']
        sample['events'][0]['flank_bp'] = 999
        provenance = evidence.provenance
        provenance['source_artifacts'].clear()
        first = technical.build_technical_boundary_evidence('sampleA', self.report, evidence)
        first['events'].clear()
        second = technical.build_technical_boundary_evidence('sampleA', self.report, evidence)
        self.assertEqual(second['events'][0]['flank_bp'], 21)
        self.assertEqual(len(second['provenance']['source_artifacts']), 1)

    def test_minimal_binding_does_not_require_optional_identity_fields(self):
        for field in ('scientific_name', 'tax_id', 'biosample'):
            del self.sample[field]
        report = {'sample_id': 'sampleA', 'source_result_sha256': self.native_sha}
        self.assertEqual(self.build(report)['status'], 'technical_boundary_warning')


if __name__ == '__main__':
    unittest.main()
