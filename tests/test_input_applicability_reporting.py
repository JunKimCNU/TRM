"""Input length diagnostics are additive, provenance-scoped report annotations.

Fixtures contain saved JSON only. No FASTQ files or detector binaries are read.
"""
from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    '_input_applicability_report_fixture', ROOT / 'tests/test_report_assessment.py')
fixture = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(fixture)


class InputApplicabilityReportingTests(unittest.TestCase):
    # Reuse the established setup without inheriting/rerunning its test methods.
    setUp = fixture.AssessmentReportTests.setUp
    command = fixture.AssessmentReportTests.command
    run_report = fixture.AssessmentReportTests.run_report
    document = fixture.AssessmentReportTests.document
    make_native_sources = fixture.AssessmentReportTests.make_native_sources

    def source_doc(self, sample='SRRknown'):
        return json.loads((self.source / sample / 'result.json').read_text())

    def store_source(self, doc, sample='SRRknown'):
        fixture.write_json(self.source / sample / 'result.json', doc)

    def set_lengths(self, bins, sample='SRRknown', prefix=False):
        """Use internally consistent small saved whole-input QC fixtures."""
        raw = self.source_doc(sample)
        records = sum(count for length, count in bins)
        bases = sum(length * count for length, count in bins)
        low, high = min(length for length, _ in bins), max(length for length, _ in bins)
        raw['summary'].update(records=records, bases=bases,
                              max_record_bases=high, sampled_prefix=prefix)
        raw['evidence_qc'].update(records=records, bases=bases,
                                 min_record_bases=low, max_record_bases=high,
                                 length_histogram=[{'lower': n, 'upper': n, 'records': c}
                                                   for n, c in bins])
        self.store_source(raw, sample)

    def expected(self, motifs, sample='SRRknown'):
        rows = fixture.read_tsv(self.assignments)
        next(r for r in rows if r['sample_id'] == sample)['related_motifs'] = ','.join(motifs)
        fixture.write_tsv(self.assignments, rows)

    def detail_rows(self, sample='SRRknown'):
        return [r for r in fixture.read_tsv(self.output / 'input_applicability.tsv')
                if r['sample_id'] == sample]

    def summary_row(self, sample='SRRknown'):
        return next(r for r in fixture.read_tsv(self.output / 'samples.tsv')
                    if r['sample_id'] == sample)

    def test_full_short_and_eligible_inputs_reach_all_report_formats(self):
        self.make_native_sources(empty=True)
        self.set_lengths([(17, 10), (18, 20), (19, 30), (20, 40)])
        self.expected(['TGTTAGGG'])
        before = {p: fixture.digest(p) for p in self.source.rglob('*.json')}
        self.run_report(native=True)
        short = self.document('SRRknown')
        eligible = self.document('SRRnone')
        view = short['input_applicability']
        other = eligible['input_applicability']
        self.assertEqual(view['status'], 'no_length_eligible_reads')
        self.assertEqual(view['scope'], 'complete_input')
        self.assertEqual(other['status'], 'length_eligible_reads_present_power_unknown')
        self.assertEqual(other['scope'], 'complete_input')
        self.assertEqual(short['schema_version'], 22)
        self.assertEqual(short['reporting_version'], '1.13.0')
        self.assertEqual(short['version'], '0.7.2')
        self.assertFalse(short['reporting_features']['input_length_applicability_changes_classification_or_adequacy'])
        self.assertFalse(short['data_adequacy']['true_absence_inferred'])
        family, = view['families']
        self.assertEqual(family['canonical_family'], fixture.canonical('TGTTAGGG'))
        self.assertEqual((family['minimum_seed_bp'], family['length_eligible_reads_lower'],
                          family['length_eligible_reads_upper']), (24, 0, 0))
        self.assertEqual(len(view['periods']), 30)
        self.assertTrue(all(x['length_eligible_reads_upper'] == 0 for x in view['periods']))
        sidecar = self.detail_rows()
        self.assertEqual(len([x for x in sidecar if x['kind'] == 'configured_period']), 30)
        expected, = [x for x in sidecar if x['kind'] == 'expected_family']
        self.assertEqual((expected['minimum_seed_bp'], expected['length_eligible_reads_upper']), ('24', '0'))
        self.assertEqual(expected['scope'], 'complete_input')
        row = self.summary_row()
        for suffix, key in [('status', 'status'), ('display', 'display_label'), ('scope', 'scope')]:
            self.assertEqual(row['input_applicability_' + suffix], view[key])
        page = (self.output / 'samples/SRRknown.html').read_text()
        index = (self.output / 'index.html').read_text()
        self.assertIn(view['display_label'], page)
        self.assertIn(view['display_label'], index)
        self.assertIn(other['display_label'], index)
        self.assertIn('TGTTAGGG', page)
        self.assertIn('../input_applicability.tsv', page)
        metrics = json.loads((self.output / 'metrics.json').read_text())
        self.assertEqual(metrics['input_applicability_status_counts']['no_length_eligible_reads'], 1)
        self.assertEqual(metrics['files']['input_applicability.tsv'], fixture.digest(self.output / 'input_applicability.tsv'))
        self.assertEqual(before, {p: fixture.digest(p) for p in self.source.rglob('*.json')})

    def test_annotation_mock_cannot_change_same_input_biological_outputs(self):
        self.make_native_sources()
        sys.path.insert(0, str(ROOT / 'scripts'))
        self.addCleanup(lambda: sys.path.remove(str(ROOT / 'scripts')))
        import report_assessment as report
        args = SimpleNamespace(source=self.source, evidence=None, catalog=self.catalog,
            assignments=self.assignments,
            taxon_references=ROOT / 'data/taxon_motif_references.tsv', taxon_reference_rows=[],
            expected_evidence=None, followup_diagnostics=None, output=self.output,
            library_evidence_rows={})
        meta = fixture.read_tsv(self.metadata)[0]
        assignment = fixture.read_tsv(self.assignments)[0]
        references = report.catalog_index(fixture.read_tsv(self.catalog))
        defaults, _ = report.read_assumptions(None)
        baseline, row0, _ = report.sample_report('SRRknown', meta, assignment, args,
                                                references, defaults, {})
        changed = copy.deepcopy(baseline['input_applicability'])
        changed.update(status='fixture_annotation_only', display_label='<script>annotation</script>',
                       recommendation='fixture explanation only')
        with mock.patch.object(report, 'build_input_applicability', return_value=changed) as builder:
            annotated, row1, _ = report.sample_report('SRRknown', meta, assignment, args,
                                                     references, defaults, {})
        self.assertEqual(builder.call_count, 1)
        self.assertEqual(annotated['input_applicability']['status'], 'fixture_annotation_only')
        for key in ['classification', 'result', 'analysis_complete', 'detected_motif',
                    'top_candidate', 'candidates', 'data_adequacy', 'reference_agreement',
                    'sequence_observation', 'expected_family_observations', 'family_observation',
                    'policy', 'next_action', 'recommendation']:
            with self.subTest(key=key):
                self.assertEqual(annotated[key], baseline[key])
        self.assertEqual(baseline['classification'], 'known_trm_detected')
        for key in ['result', 'detected_motif', 'supporting_reads', 'full_repeat_reads', 'data_adequacy']:
            self.assertEqual(row0[key], row1[key])
        rendered = ''.join(report.input_applicability_html(annotated['input_applicability']))
        self.assertIn('&lt;script&gt;annotation&lt;/script&gt;', rendered)
        self.assertNotIn('<script>annotation</script>', rendered)

    def test_prefix_is_not_promoted_to_complete_input_ineligibility(self):
        self.make_native_sources(empty=True)
        self.set_lengths([(20, 100)], prefix=True)
        self.expected(['TGTTAGGG'])
        self.run_report(native=True)
        doc = self.document('SRRknown')
        view = doc['input_applicability']
        self.assertFalse(doc['analysis_complete'])
        self.assertEqual(view['status'], 'analysis_incomplete')
        self.assertEqual(view['scope'], 'analyzed_records')
        family, = view['families']
        self.assertEqual(family['status'], 'analysis_incomplete')
        self.assertEqual(family['analyzed_records_length_status'], 'no_length_eligible_reads')
        self.assertEqual(family['length_eligible_reads_upper'], 0)
        self.assertEqual(self.summary_row()['input_applicability_status'], 'analysis_incomplete')
        self.assertTrue(all(x['input_status'] == 'analysis_incomplete' for x in self.detail_rows()))
        expected, = [x for x in self.detail_rows() if x['kind'] == 'expected_family']
        self.assertEqual(expected['analyzed_records_length_status'], 'no_length_eligible_reads')
        self.assertEqual(expected['bounds_method'], 'read_lengths_cannot_reach_seed')

    def test_source_failure_does_not_export_zero_length_eligible_count(self):
        self.make_native_sources(empty=True)
        self.run_report(native=True)
        doc = self.document('SRRfailed')
        view = doc['input_applicability']
        self.assertFalse(doc['analysis_complete'])
        self.assertEqual(view['status'], 'source_unavailable')
        self.assertEqual(view['scope'], 'unavailable')
        rows = self.detail_rows('SRRfailed')
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['kind'], 'input')
        self.assertEqual(rows[0]['length_eligible_reads_lower'], '')
        self.assertEqual(rows[0]['length_eligible_reads_upper'], '')
        self.assertIn(view['display_label'], (self.output / 'samples/SRRfailed.html').read_text())

    def test_legacy_without_histogram_retains_incomplete_scope_not_invented_full_qc(self):
        self.make_native_sources(empty=True)
        raw = self.source_doc()
        raw.pop('evidence_qc')
        self.store_source(raw)
        self.run_report(native=True)
        doc = self.document('SRRknown')
        view = doc['input_applicability']
        self.assertFalse(doc['analysis_complete'])
        self.assertEqual(view['status'], 'analysis_incomplete')
        self.assertNotEqual(view['scope'], 'complete_input')
        self.assertEqual(view['provenance']['qc_source'], 'original_summary')
        self.assertTrue(all(x['status'] == 'analysis_incomplete' for x in view['periods']))
        self.assertFalse(doc['data_adequacy']['true_absence_inferred'])

    def test_verified_targeted_qc_uses_original_discovery_config_and_source_hashes(self):
        # Base fixture has QC only in separately completed, hash-validated evidence.
        self.run_report()
        doc = self.document('SRRknown')
        view = doc['input_applicability']
        self.assertEqual(view['status'], 'length_eligible_reads_present_power_unknown')
        self.assertEqual(view['scope'], 'complete_input')
        self.assertEqual(view['provenance']['configuration_source'], 'original_discovery')
        self.assertEqual(view['provenance']['qc_source'], 'full_input_evidence')
        self.assertEqual(view['provenance']['configuration_result_sha256'],
                         fixture.digest(self.source / 'SRRknown/result.json'))
        self.assertEqual(view['provenance']['qc_result_sha256'],
                         fixture.digest(self.evidence / 'SRRknown/evidence.json'))
        self.assertEqual(doc['classification'], 'known_trm_detected')
        self.assertEqual(next(x for x in view['periods'] if x['period_bp'] == 14)['minimum_seed_bp'], 42)

    def test_family_aliases_seed_rounding_and_outside_discovery_range(self):
        self.make_native_sources(empty=True)
        self.expected(['TTGGA', 'TCCAA', 'TTGGATTGGA',
                       'TTTAGGG', 'TGTGTGGGTGTGGT'])
        raw = self.source_doc()
        raw['config']['max_k'] = 6
        self.store_source(raw)
        self.run_report(native=True)
        view = self.document('SRRknown')['input_applicability']
        self.assertEqual(len(view['periods']), 4)
        families = {x['period_bp']: x for x in view['families']}
        self.assertEqual(set(families), {5, 7, 14})
        self.assertEqual({p: x['minimum_seed_bp'] for p, x in families.items()}, {5: 25, 7: 28, 14: 42})
        self.assertEqual(set(families[5]['expected_motifs']), {'TTGGA', 'TCCAA', 'TTGGATTGGA'})
        self.assertFalse(families[5]['outside_discovery_range'])
        self.assertTrue(families[7]['outside_discovery_range'])
        self.assertTrue(families[14]['outside_discovery_range'])
        self.assertEqual(len([r for r in self.detail_rows() if r['kind'] == 'expected_family']), 3)
        self.assertIn('de novo 탐색 주기 범위 밖', (self.output / 'samples/SRRknown.html').read_text())

    def test_missing_detector_config_is_explicit_unknown_not_default_threshold(self):
        self.make_native_sources(empty=True)
        raw = self.source_doc()
        raw['config'].pop('seed_bases')
        self.store_source(raw)
        self.run_report(native=True)
        view = self.document('SRRknown')['input_applicability']
        self.assertEqual(view['status'], 'config_unavailable')
        self.assertFalse(view['periods'])
        rows = self.detail_rows()
        self.assertTrue(all(x['minimum_seed_bp'] == '' for x in rows))

    def test_invalid_histogram_is_visible_without_replacing_observation_values(self):
        self.make_native_sources(empty=True)
        raw = self.source_doc()
        raw['evidence_qc']['length_histogram'][0]['records'] -= 1
        self.store_source(raw)
        self.run_report(native=True)
        doc = self.document('SRRknown')
        view = doc['input_applicability']
        self.assertEqual(view['status'], 'qc_invalid')
        self.assertTrue(all(x['status'] == 'qc_invalid' and
                            x['length_eligible_reads_lower'] is None and
                            x['length_eligible_reads_upper'] is None
                            for x in view['periods']))
        self.assertFalse(doc['data_adequacy']['true_absence_inferred'])
        self.assertTrue(all(x['length_eligible_reads_upper'] == '' for x in self.detail_rows()))


if __name__ == '__main__':
    unittest.main()
