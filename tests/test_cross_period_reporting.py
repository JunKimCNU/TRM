"""Offline cache/report integration of the optional cross-period diagnostic.

The input is a tiny saved-evidence fixture. No FASTQ, reference search, periodic
model execution, threshold change, or biological confirmation is involved.
"""
from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'tests')]
import cross_period_evidence as evidence_module
import report_runtime
from summarize_interpretation import summarize_sample
import test_report_assessment as fixture
import test_cross_period_evidence as evidence_fixture

TABLES = ('cross_period_diagnostics.tsv', 'cross_period_families.tsv',
          'cross_period_intervals.tsv', 'cross_period_events.tsv',
          'cross_period_contexts.tsv', 'cross_period_positions.tsv')


class CrossPeriodReportingTests(unittest.TestCase):
    setUp = fixture.AssessmentReportTests.setUp
    command = fixture.AssessmentReportTests.command
    run_report = fixture.AssessmentReportTests.run_report
    document = fixture.AssessmentReportTests.document

    def save_evidence(self):
        self.payload = evidence_fixture.make_document(self.root)
        sample = self.payload['samples'][0]
        sample.update(sample_id='SRRknown', scientific_name='Fixture species', tax_id='6239',
                      biosample='SAMNfixture_SRRknown')
        # Both compared families must already occur in the native/report input.
        source_path = self.source/'SRRknown/result.json'
        target_path = self.evidence/'SRRknown/evidence.json'
        status_path = self.evidence/'SRRknown/status.json'
        source = json.loads(source_path.read_text())
        target = json.loads(target_path.read_text())
        observed = {row['motif'] for row in source['candidates']}
        for motif in sample['motifs'].values():
            family = fixture.canonical(motif)
            if family not in observed:
                candidate = deepcopy(source['candidates'][0])
                candidate['motif'] = family
                source['candidates'].append(candidate)
                row = deepcopy(target['candidate_evidence'][0])
                row['motif'] = family
                target['candidate_evidence'].append(row)
                target['target_motifs'].append(family)
                observed.add(family)
        fixture.write_json(source_path, source)
        fixture.write_json(target_path, target)
        status = json.loads(status_path.read_text())
        status['result_sha256'] = fixture.digest(target_path)
        status['identity']['source_sha256'] = fixture.digest(source_path)
        status['identity']['targets']['motifs'] = list(target['target_motifs'])
        fixture.write_json(status_path, status)
        sample['source_result_sha256'] = fixture.digest(source_path)
        self.sidecar = self.root/'cross-period.json'
        self.pin = Path(self.payload['source_artifacts'][0]['path'])
        if not self.pin.is_absolute():
            self.pin = self.root/self.pin
        self.reseal()
        return ['--cross-period-evidence', str(self.sidecar)]

    def reseal(self):
        self.payload['payload_sha256'] = evidence_module.evidence_payload_sha256(self.payload)
        fixture.write_json(self.sidecar, self.payload)

    def resources(self):
        return {row['sample_id']:row for row in fixture.read_tsv(self.output/'interpretation_resources.tsv')}

    def test_overlay_changes_only_its_namespace_and_keeps_existing_tables(self):
        extra = self.save_evidence()
        self.run_report()
        previous = self.document('SRRknown')
        candidates = (self.output/'candidates.tsv').read_bytes()
        self.run_report(extra=extra)
        current = self.document('SRRknown')
        self.assertEqual(set(previous), set(current))
        for key in previous:
            if key != 'cross_period_evidence':
                self.assertEqual(previous[key], current[key], key)
        self.assertEqual(candidates, (self.output/'candidates.tsv').read_bytes())
        self.assertEqual(summarize_sample(previous, Path('fixed.json'), '0'*64),
                         summarize_sample(current, Path('fixed.json'), '0'*64))
        diagnostic = current['cross_period_evidence']
        self.assertNotEqual(diagnostic['status'], 'unassessed')
        self.assertEqual(diagnostic['selected_interval_count'], len(self.payload['samples'][0]['intervals']))
        self.assertEqual(len(fixture.read_tsv(self.output/TABLES[2])), diagnostic['selected_interval_count'])
        self.assertEqual(len(fixture.read_tsv(self.output/TABLES[3])), diagnostic['selected_event_count'])
        self.assertEqual(len(fixture.read_tsv(self.output/TABLES[4])), diagnostic['qualified_context_count'])
        self.assertEqual(len(fixture.read_tsv(self.output/TABLES[5])), diagnostic['strict_contrast_position_count'])
        interval_rows = {r['interval_id']:r for r in fixture.read_tsv(self.output/TABLES[2])}
        for original in self.payload['samples'][0]['intervals']:
            exported = interval_rows[original['interval_id']]
            self.assertEqual(exported['deletion_quality'], '')
            self.assertEqual(json.loads(exported['candidate_operations']), original['candidate_operations'])
            self.assertEqual(json.loads(exported['comparator_operations']), original['comparator_operations'])
        position_rows = fixture.read_tsv(self.output/TABLES[5])
        for position, exported in zip(diagnostic['positions'], position_rows):
            self.assertEqual(exported['conditional_phred33'], str(position['conditional_phred33']))
            self.assertEqual(exported['quality_character'], position['quality_character'])
        self.assertFalse(current['reporting_features']['cross_period_changes_classification_or_adequacy'])
        text = (self.output/'samples/SRRknown.html').read_text()
        for phrase in ('서로 다른 반복 주기', '같은 위치', '모든 최적 편집 경로',
                       '독립적인 정확도 검증은 아닙니다', 'deletion', 'Phred+33', 'read가 덮은 길이와 다릅니다'):
            self.assertIn(phrase, text)
        metrics = json.loads((self.output/'metrics.json').read_text())
        for name in TABLES:
            self.assertIn('../'+name, text)
            self.assertEqual(metrics['files'][name], fixture.digest(self.output/name))
        row = next(r for r in fixture.read_tsv(self.output/'samples.tsv') if r['sample_id']=='SRRknown')
        self.assertEqual(row['cross_period_selected_event_count'], str(diagnostic['selected_event_count']))

    def test_omission_and_unknown_sample_are_unassessed_null(self):
        self.run_report()
        for sid in ('SRRknown', 'SRRnone', 'SRRfailed'):
            diagnostic = self.document(sid)['cross_period_evidence']
            self.assertEqual(diagnostic['status'], 'unassessed')
            self.assertTrue(all(diagnostic[key] is None for key in evidence_module.SUMMARY_FIELDS))
        for name in TABLES[1:]:
            self.assertEqual(fixture.read_tsv(self.output/name), [])
        self.run_report(extra=self.save_evidence())
        for sid in ('SRRnovel', 'SRRnone', 'SRRfailed'):
            diagnostic = self.document(sid)['cross_period_evidence']
            self.assertEqual(diagnostic['reason'], 'sample_not_in_evidence')
            self.assertTrue(all(diagnostic[key] is None for key in evidence_module.SUMMARY_FIELDS))
        rows = {r['sample_id']:r for r in fixture.read_tsv(self.output/TABLES[0])}
        self.assertEqual(rows['SRRfailed']['selected_event_count'], '')

    def test_empty_selected_sample_reports_no_observed_events_not_missing_input(self):
        extra = self.save_evidence()
        sample = self.payload['samples'][0]
        for key in ('events', 'contexts', 'intervals'):
            sample[key] = []
        self.reseal()
        self.run_report(extra=extra)
        diagnostic = self.document('SRRknown')['cross_period_evidence']
        self.assertEqual(diagnostic['selected_event_count'], 0)
        self.assertEqual(diagnostic['selected_interval_count'], 0)
        self.assertIsNone(diagnostic['inward_mate_event_count'])
        row = next(r for r in fixture.read_tsv(self.output/TABLES[0]) if r['sample_id']=='SRRknown')
        self.assertEqual(row['selected_event_count'], '0')
        self.assertEqual(row['inward_mate_event_count'], '')
        for name in TABLES[1:]:
            self.assertEqual(fixture.read_tsv(self.output/name), [])

    def test_review_text_is_escaped_without_changing_json(self):
        extra = self.save_evidence()
        hostile = '<img src=x onerror="alert(1)">'
        self.payload['samples'][0]['review_notes'] = [hostile]
        self.reseal()
        self.run_report(extra=extra)
        text = (self.output/'samples/SRRknown.html').read_text()
        self.assertNotIn(hostile, text)
        self.assertIn('&lt;img src=x onerror=&quot;alert(1)&quot;&gt;', text)
        self.assertIn(hostile, self.document('SRRknown')['cross_period_evidence']['review_notes'])

    def test_stale_native_and_identity_bindings_fail_closed(self):
        extra = self.save_evidence()
        for key, value in [('source_result_sha256', '0'*64), ('biosample', 'SAMNwrong')]:
            with self.subTest(key=key):
                old = self.payload['samples'][0][key]
                self.payload['samples'][0][key] = value
                self.reseal()
                self.run_report(success=False, extra=extra)
                self.assertEqual(json.loads((self.output/'interpretation_status.json').read_text())['status'], 'failed')
                self.payload['samples'][0][key] = old

    def test_missing_candidate_family_fails_even_with_new_native_pin(self):
        extra = self.save_evidence()
        source_path = self.source/'SRRknown/result.json'
        target_path = self.evidence/'SRRknown/evidence.json'
        status_path = self.evidence/'SRRknown/status.json'
        family = fixture.canonical(self.payload['samples'][0]['motifs']['candidate'])
        source = json.loads(source_path.read_text())
        source['candidates'] = [r for r in source['candidates'] if r['motif'] != family]
        fixture.write_json(source_path, source)
        target = json.loads(target_path.read_text())
        target['candidate_evidence'] = [r for r in target['candidate_evidence'] if r['motif'] != family]
        target['target_motifs'] = [m for m in target['target_motifs'] if m != family]
        fixture.write_json(target_path, target)
        status = json.loads(status_path.read_text())
        status['identity']['source_sha256'] = fixture.digest(source_path)
        status['result_sha256'] = fixture.digest(target_path)
        status['identity']['targets']['motifs'] = list(target['target_motifs'])
        fixture.write_json(status_path, status)
        self.payload['samples'][0]['source_result_sha256'] = fixture.digest(source_path)
        self.reseal()
        outcome = self.run_report(success=False, extra=extra)
        self.assertIn('family', outcome.stderr)

    def test_explicit_missing_and_malformed_inputs_do_not_publish_success(self):
        outcome = self.run_report(success=False,
            extra=['--cross-period-evidence', str(self.root/'missing.json')])
        self.assertIn('cannot read', outcome.stderr)
        self.assertFalse((self.output/'metrics.json').exists())
        self.sidecar = self.root/'malformed.json'
        self.sidecar.write_text('{"schema_version": NaN}')
        self.run_report(success=False, extra=['--cross-period-evidence', str(self.sidecar)])
        self.assertFalse((self.output/'metrics.json').exists())

    def test_unchanged_completed_resume_is_exact_and_changed_pin_invalidates(self):
        extra = self.save_evidence()+['--resume', '--jobs', '2']
        self.run_report(extra=extra)
        docs = {p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')}
        tables = {name:(self.output/name).read_bytes() for name in TABLES}
        self.run_report(extra=extra)
        self.assertEqual(docs, {p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')})
        self.assertEqual(tables, {name:(self.output/name).read_bytes() for name in TABLES})
        self.assertTrue(all(row['mode']=='reused' for row in self.resources().values()))
        self.pin.write_text('Changed audited-source fixture; old hash no longer binds it.\n')
        self.run_report(success=False, extra=extra)
        self.assertEqual(json.loads((self.output/'interpretation_status.json').read_text())['status'], 'failed')

    def test_changed_sealed_sidecar_recomputes_and_option_removal_clears_tables(self):
        extra = self.save_evidence()+['--resume']
        self.run_report(extra=extra)
        self.payload['samples'][0]['review_notes'].append('New annotation, no scientific promotion.')
        self.reseal()
        self.run_report(extra=extra)
        self.assertTrue(all(row['mode']=='recomputed' for row in self.resources().values()))
        self.assertIn('New annotation, no scientific promotion.',
                      self.document('SRRknown')['cross_period_evidence']['review_notes'])
        self.run_report(extra=['--resume'])
        self.assertEqual(self.document('SRRknown')['cross_period_evidence']['status'], 'unassessed')
        self.assertTrue(all(row['mode']=='recomputed' for row in self.resources().values()))
        for name in TABLES[1:]:
            self.assertEqual(fixture.read_tsv(self.output/name), [])

    def test_parallel_and_sequential_are_identical(self):
        extra = self.save_evidence()
        self.run_report(extra=extra+['--jobs', '1'])
        docs = {p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')}
        tables = {name:(self.output/name).read_bytes() for name in TABLES}
        self.run_report(extra=extra+['--jobs', '2'])
        self.assertEqual(docs, {p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')})
        self.assertEqual(tables, {name:(self.output/name).read_bytes() for name in TABLES})

    def test_output_cannot_overwrite_pinned_source(self):
        extra = self.save_evidence()
        protected = self.output/'protected-review.txt'
        protected.parent.mkdir()
        protected.write_bytes(self.pin.read_bytes())
        self.payload['source_artifacts'][0]['path'] = str(protected)
        self.reseal()
        outcome = self.run_report(success=False, extra=extra)
        self.assertIn('must not contain input provenance files', outcome.stderr)
        self.assertEqual(protected.read_bytes(), self.pin.read_bytes())

    def test_cache_receipt_requires_every_new_table(self):
        self.run_report(extra=self.save_evidence())
        receipt_path = self.output/'.report-cache/SRRknown/receipt.json'
        original = json.loads(receipt_path.read_text())
        self.assertTrue(report_runtime.verify_artifacts(original, output=self.output, sample='SRRknown'))
        bundle_path = Path(original['bundle'])
        bundle = json.loads(bundle_path.read_text())
        for missing in TABLES:
            with self.subTest(missing=missing):
                receipt = deepcopy(original)
                truncated = deepcopy(bundle)
                truncated['chunks'].remove(missing)
                fixture.write_json(bundle_path, truncated)
                receipt['artifacts'][str(bundle_path)] = report_runtime.binding(bundle_path)[0]
                receipt['artifacts'].pop(str(bundle_path.parent/missing))
                self.assertFalse(report_runtime.verify_artifacts(receipt, output=self.output, sample='SRRknown'))
        fixture.write_json(bundle_path, bundle)


if __name__ == '__main__':
    unittest.main()
