"""Orthogonal summaries retain observations and every secondary nomination."""
import copy
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import summarize_interpretation as module


def fixture(sample='sample'):
    known = {'canonical_family': 'AACCCT', 'motif': 'TTAGGG', 'reads': 123,
             'reference_status': 'established', 'known_observed': True,
             'candidate_status': 'known_repeat_detected_boundary_unresolved',
             'novel_candidate': False, 'geometry_nomination': False}
    novel = {'canonical_family': 'ACAT', 'motif': 'ACAT', 'reads': 40,
             'reference_status': 'not_in_catalogue', 'known_observed': False,
             'candidate_status': 'novel_provisional_candidate',
             'novel_candidate': True, 'geometry_nomination': True}
    secondary = {**novel, 'canonical_family': 'AGAT', 'motif': 'AGAT', 'reads': 15}
    geometry = {**novel, 'canonical_family': 'ACT', 'motif': 'ACT', 'novel_candidate': False,
                'candidate_status': 'novel_geometry_nomination'}
    return {'program': 'trm-assessment', 'sample_id': sample, 'scientific_name': 'Fixture species',
        'tax_id': '999', 'classification': 'novel_provisional_candidate', 'analysis_complete': True,
        'detected_motif': 'ACAT', 'top_candidate': novel,
        'candidates': [known, novel, secondary, geometry],
        'reference_assignment': {'reference_scope': 'clade', 'reference_status': 'supported'},
        'reference_agreement': {'main_trm_status': 'not_confirmed'}, 'externally_confirmed': False,
        'family_observation': {'status': 'exact_family_observed', 'families': [
            {'canonical_family': 'AACCCT', 'expected_motifs': ['TTAGGG'], 'reference_scope': 'clade',
             'reference_status': 'supported', 'status': 'exact_family_observed',
             'exact_measurement_status': 'reported_positive',
             'exact_current': {'status': 'reported_positive', 'reads': 123, 'read_counts': [123],
                               'positive_observed': True, 'explicit_zero_observed': False},
             'exact_original': {'status': 'reported_positive', 'reads': 123}, 'nearby_family_count': 0}]},
        'data_adequacy': {'status': 'unassessable_library_or_genome',
                         'missing_assumptions': ['genome_size_bp'], 'library': {'uniform_wgs_model_blockers': ['library_selection_not_verified_uniform_genome_sampling']}},
        'metadata': {'library_selection': 'RANDOM PCR', 'whole_genome_amplification': 'unknown'},
        'library_context': {'limited': True, 'reasons': ['library_selection_not_verified_uniform_genome_sampling']}}


class SummaryTests(unittest.TestCase):
    def test_known_expected_and_secondary_provisional_axes_coexist_without_promotion(self):
        document = fixture(); before = copy.deepcopy(document)
        row, expected, queue = module.summarize_sample(document, '/fixture/sample.json', 'pin')
        self.assertEqual(row['known_positive_families'], 1)
        self.assertEqual(row['known_positive_motifs'], 'TTAGGG')
        self.assertEqual(row['known_observed_motifs'], 'TTAGGG')
        self.assertEqual(row['provisional_candidate_families'], 2)
        self.assertEqual(row['secondary_provisional_families'], 1)
        self.assertEqual(row['geometry_only_families'], 1)
        self.assertEqual(row['main_trm_status'], 'not_confirmed')
        self.assertEqual(row['reference_strength'], 'weak_clade_expectation')
        self.assertFalse(row['externally_confirmed_recorded'])
        self.assertEqual({r['canonical_family'] for r in queue if r['priority'] in {50, 60}}, {'ACAT', 'AGAT', 'ACT'})
        self.assertEqual(expected[0]['current_reads'], 123)
        self.assertEqual(document, before)

    def test_not_reported_reported_zero_and_explicit_zero_stay_distinct(self):
        for status, reads, explicit, priority in [('not_reported', None, False, 30),
                                                 ('reported_zero', 0, False, None),
                                                 ('explicit_zero', 0, True, 35)]:
            with self.subTest(status=status):
                document = fixture(); family = document['family_observation']['families'][0]
                family['exact_current'] = {'status': status, 'reads': reads, 'read_counts': [] if reads is None else [0],
                                           'positive_observed': False, 'explicit_zero_observed': explicit}
                family['exact_measurement_status'] = status
                row, expected, queue = module.summarize_sample(document, '/fixture/sample.json', 'pin')
                self.assertEqual(expected[0]['current_reads'], reads)
                self.assertEqual(expected[0]['original_reads'], 123)
                self.assertEqual(expected[0]['explicit_zero'], explicit)
                matches = [r for r in queue if r['canonical_family'] == 'AACCCT']
                self.assertEqual([r['priority'] for r in matches], [] if priority is None else [priority])

    def test_pcr_does_not_become_wga_and_disagreement_retains_both_counts(self):
        document = fixture(); family = document['family_observation']['families'][0]
        family.update(measurement_disagreement=True)
        family['exact_current'].update(status='measurement_disagreement', reads=None, read_counts=[0, 123])
        row, expected, queue = module.summarize_sample(document, '/fixture/sample.json', 'pin')
        self.assertTrue(row['ordinary_pcr_selection'])
        self.assertFalse(row['explicit_wga_blocker'])
        self.assertEqual(row['whole_genome_amplification_metadata'], 'unknown')
        self.assertEqual(expected[0]['current_read_counts'], '[0,123]')
        self.assertEqual([r['priority'] for r in queue if r['canonical_family'] == 'AACCCT'], [20])

    def test_file_outputs_deterministic_and_sources_immutable(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); assessment = root / 'assessment'; samples = assessment / 'samples'
            samples.mkdir(parents=True)
            for sample in ['z', 'a']:
                (samples / f'{sample}.json').write_text(json.dumps(fixture(sample)))
            before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in samples.iterdir()}
            first = module.summarize(assessment, root / 'first', allow_json_only=True)
            second = module.summarize(assessment, root / 'second', allow_json_only=True)
            for name in ['samples.tsv', 'expected_families.tsv', 'review_queue.tsv']:
                self.assertEqual((root / 'first' / name).read_bytes(), (root / 'second' / name).read_bytes())
            self.assertEqual(first['candidate_rows_examined'], 8)
            self.assertEqual(first['samples_with_known_and_provisional_families'], 2)
            self.assertEqual(first['expected_positive_by_reference_strength'], {'weak_clade_expectation': 2})
            self.assertEqual(first['main_trm_status_counts'], {'not_confirmed': 2})
            self.assertEqual(first['missing_power_assumption_counts'], {'genome_size_bp': 2})
            self.assertEqual(first['missing_library_protocol_samples'], 2)
            self.assertEqual(first['uniform_wgs_model_blocker_counts'], {'library_selection_not_verified_uniform_genome_sampling': 2})
            self.assertEqual(before, {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in samples.iterdir()})
            with self.assertRaises(FileExistsError): module.summarize(assessment, root / 'first', allow_json_only=True)
            with self.assertRaises(ValueError): module.summarize(assessment, assessment / 'summary', allow_json_only=True)
            with self.assertRaisesRegex(ValueError, 'authoritative samples.tsv'):
                module.summarize(assessment, root / 'implicit_fallback')

    def write_manifested_cohort(self, assessment, ids):
        samples = assessment / 'samples'; samples.mkdir(parents=True)
        for sample in ids:
            (samples / f'{sample}.json').write_text(json.dumps(fixture(sample)))
        manifest = 'sample_id\n' + ''.join(sample + '\n' for sample in ids)
        (assessment / 'samples.tsv').write_text(manifest)
        (assessment / 'metrics.json').write_text(json.dumps({'samples': len(ids), 'sample_ids': ids,
            'files': {'samples.tsv': hashlib.sha256(manifest.encode()).hexdigest()}}))
        (assessment / 'interpretation_status.json').write_text(json.dumps({'status': 'complete',
            'samples': len(ids), 'sample_ids': ids, 'completed': len(ids), 'failed': 0}))

    def test_shrunk_cohort_ignores_stale_sample_json_and_accepts_legacy_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); assessment = root / 'assessment'
            self.write_manifested_cohort(assessment, ['current'])
            (assessment / 'samples' / 'stale.json').write_text('not even valid JSON')
            summary = module.summarize(assessment, root / 'summary')
            self.assertEqual(summary['samples'], 1)
            self.assertEqual(summary['cohort_authority'], 'samples.tsv')
            self.assertEqual([Path(pin['path']).stem for pin in summary['source_pins']], ['current'])
            self.assertEqual(len(summary['cohort_control_pins']), 3)
            (assessment / 'interpretation_status.json').unlink()
            self.assertEqual(module.summarize(assessment, root / 'legacy')['samples'], 1)

    def test_incomplete_status_or_inconsistent_manifest_cannot_be_bypassed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); assessment = root / 'assessment'
            self.write_manifested_cohort(assessment, ['current'])
            status_path = assessment / 'interpretation_status.json'
            baseline_status = json.loads(status_path.read_text())
            metrics_path = assessment / 'metrics.json'
            baseline_metrics = json.loads(metrics_path.read_text())
            cases = [('status_' + status, status_path, {**baseline_status, 'status': status})
                     for status in ['preflight', 'running', 'workers_complete', 'failed']]
            cases += [('wrong_ids', status_path, {**baseline_status, 'sample_ids': ['stale']}),
                      ('not_all_completed', status_path, {**baseline_status, 'completed': 0}),
                      ('failed_sample', status_path, {**baseline_status, 'failed': 1}),
                      ('wrong_count', metrics_path, {**baseline_metrics, 'samples': 2}),
                      ('bad_hash', metrics_path, {**baseline_metrics, 'files': {'samples.tsv': 'wrong'}})]
            for name, path, value in cases:
                with self.subTest(name=name):
                    status_path.write_text(json.dumps(baseline_status))
                    metrics_path.write_text(json.dumps(baseline_metrics))
                    path.write_text(json.dumps(value))
                    with self.assertRaises(ValueError):
                        module.summarize(assessment, root / name, allow_json_only=True)

    def test_concurrent_cohort_generation_change_prevents_completed_summary(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); assessment = root / 'assessment'
            self.write_manifested_cohort(assessment, ['current'])
            original = module.summarize_sample
            for control in ['metrics.json', 'samples.tsv', 'interpretation_status.json']:
                path = assessment / control; before = path.read_bytes()
                def replace_control(*args, **kwargs):
                    result = original(*args, **kwargs)
                    path.write_bytes(before + b'\n')
                    return result
                with self.subTest(control=control), patch.object(module, 'summarize_sample', side_effect=replace_control):
                    output = root / control
                    with self.assertRaisesRegex(ValueError, 'control changed'):
                        module.summarize(assessment, output)
                    self.assertFalse((output / 'summary.json').exists())
                path.write_bytes(before)

    def test_protocol_and_model_blockers_are_independent_of_sample_classification(self):
        document = fixture()
        document.update(classification='known_trm_detected')
        document['library_context'].update(limited=False, reasons=[])
        document['metadata']['library_construction_protocol'] = 'Recorded ordinary library protocol'
        row, _, _ = module.summarize_sample(document, '/fixture/sample.json', 'pin')
        self.assertTrue(row['protocol_available'])
        self.assertFalse(row['library_limited'])
        self.assertEqual(row['uniform_wgs_model_blockers'], 'library_selection_not_verified_uniform_genome_sampling')
        self.assertFalse(row['explicit_wga_blocker'])

    def test_failed_analysis_is_reviewed_without_erasing_observations(self):
        document = fixture(); document.update(analysis_complete=False, classification='analysis_failed')
        row, expected, queue = module.summarize_sample(document, '/fixture/sample.json', 'pin')
        self.assertEqual(row['candidate_count'], 4)
        self.assertEqual(row['known_positive_families'], 1)
        self.assertIn(10, [r['priority'] for r in queue])
        document['candidates'].append(document['candidates'][0])
        with self.assertRaises(ValueError): module.summarize_sample(document, '/fixture/sample.json', 'pin')


if __name__ == '__main__': unittest.main()
