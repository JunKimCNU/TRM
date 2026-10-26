"""Report/cache integration of optional technical-boundary evidence."""
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from technical_boundary_evidence import evidence_payload_sha256
import test_report_assessment as fixture


class TechnicalBoundaryReportingTests(unittest.TestCase):
    setUp = fixture.AssessmentReportTests.setUp
    command = fixture.AssessmentReportTests.command
    run_report = fixture.AssessmentReportTests.run_report
    document = fixture.AssessmentReportTests.document

    def save_evidence(self, sample='SRRknown'):
        self.pin = self.root / 'saved-diagnostic.txt'
        self.pin.write_text('Saved selected-junction fixture; not raw reads.\n')
        self.sidecar = self.root / 'technical.json'
        event = dict(event_id='selected-event', original_pair_index_1based=7, mate=2,
            family='AAGCCT', flank_bp=21, maximum_exact_match_bp_overlapping_flank=19,
            union_flank_overlap_by_threshold={'12': 19, '16': 19, '20': 0},
            adjacent_unmatched_flank_by_threshold={'12': 2, '16': 2, '20': 21},
            paired_readthrough_compatible=True)
        tiers = [dict(threshold_bp=t, exact_match_associations=2 if t<20 else 0,
            full_reads_with_matches=2 if t<20 else 0, pairs_with_matches=1 if t<20 else 0,
            events_with_flank_overlap=1 if t<20 else 0) for t in (12, 16, 20)]
        self.payload = dict(schema_version=1, scope='selected_junction_subset',
            thresholds_bp=[12, 16, 20],
            source_artifacts=[dict(path=str(self.pin), sha256=fixture.digest(self.pin))],
            reference_panel=[dict(reference_id='test_adapter', length_bp=33)],
            samples=[dict(sample_id=sample,
                source_result_sha256=fixture.digest(self.source/sample/'result.json'),
                scientific_name='Fixture species', tax_id='6239',
                biosample='SAMNfixture_'+sample,
                families=['AAGCCT'], availability=dict(full_reads=True, mates=True),
                selected_event_count=1, selected_pair_count=1, full_read_count=2,
                tiers=tiers, events=[event])])
        self.reseal()
        return ['--technical-boundary-evidence', str(self.sidecar)]

    def reseal(self):
        self.payload['payload_sha256'] = evidence_payload_sha256(self.payload)
        fixture.write_json(self.sidecar, self.payload)

    def test_19bp_warning_reaches_html_tsv_and_preserves_all_biology(self):
        self.run_report()
        previous = self.document('SRRknown')
        self.run_report(extra=self.save_evidence())
        current = self.document('SRRknown')
        for key in previous:
            if key not in {'technical_boundary_diagnostics'}:
                self.assertEqual(previous[key], current[key], key)
        diagnostic = current['technical_boundary_diagnostics']
        self.assertEqual(diagnostic['status'], 'technical_boundary_warning')
        self.assertEqual(diagnostic['flagged_event_count'], 1)
        self.assertEqual([t['events_with_flank_overlap'] for t in diagnostic['tiers']], [1, 1, 0])
        detail = (self.output/'samples/SRRknown.html').read_text()
        self.assertIn('기술 서열과 반복 경계의 관계', detail)
        rows = fixture.read_tsv(self.output/'technical_boundary_events.tsv')
        self.assertEqual(len(rows), 3)
        self.assertEqual([r['adjacent_unmatched_flank_bp'] for r in rows], ['2','2','21'])
        self.assertEqual(fixture.read_tsv(self.output/'samples.tsv')[0]['technical_boundary_status'],
                         self.document(fixture.read_tsv(self.output/'samples.tsv')[0]['sample_id'])['technical_boundary_diagnostics']['status'])

    def test_unassessed_is_null_and_has_no_fabricated_event(self):
        self.run_report()
        for sid in ('SRRknown', 'SRRnone', 'SRRfailed'):
            diagnostic = self.document(sid)['technical_boundary_diagnostics']
            self.assertEqual(diagnostic['status'], 'unassessed')
            self.assertIsNone(diagnostic['flagged_event_count'])
            self.assertIsNone(diagnostic['selected_event_count'])
        self.assertEqual(fixture.read_tsv(self.output/'technical_boundary_events.tsv'), [])
        self.assertTrue(all(r['flagged_event_count']=='' for r in
            fixture.read_tsv(self.output/'technical_boundary_diagnostics.tsv')))

    def test_optional_subset_does_not_label_other_samples_negative(self):
        self.run_report(extra=self.save_evidence())
        for sid in ('SRRnovel', 'SRRnone', 'SRRfailed'):
            diagnostic = self.document(sid)['technical_boundary_diagnostics']
            self.assertEqual(diagnostic['status'], 'unassessed')
            self.assertIsNone(diagnostic['flagged_event_count'])

    def test_wrong_native_source_fails_report(self):
        extra=self.save_evidence()
        self.payload['samples'][0]['source_result_sha256']='0'*64
        self.reseal()
        self.run_report(success=False, extra=extra)
        self.assertEqual(json.loads((self.output/'interpretation_status.json').read_text())['status'], 'failed')

    def test_resume_keeps_evidence_and_invalidates_changed_transitive_source(self):
        extra=self.save_evidence()+['--resume', '--jobs', '2']
        self.run_report(extra=extra)
        original=self.document('SRRknown')
        self.run_report(extra=extra)
        self.assertEqual(original, self.document('SRRknown'))
        self.pin.write_text('Changed diagnostic; recorded source SHA is no longer valid.\n')
        self.run_report(success=False, extra=extra)
        self.assertEqual(json.loads((self.output/'interpretation_status.json').read_text())['status'], 'failed')

    def test_removing_optional_evidence_does_not_reuse_old_warning(self):
        self.run_report(extra=self.save_evidence()+['--resume'])
        self.run_report(extra=['--resume'])
        self.assertEqual(self.document('SRRknown')['technical_boundary_diagnostics']['status'], 'unassessed')

    def test_output_cannot_contain_transitive_source(self):
        extra=self.save_evidence()
        protected=self.output/'diagnostic.txt'
        protected.parent.mkdir(parents=True, exist_ok=True)
        protected.write_bytes(self.pin.read_bytes())
        self.payload['source_artifacts'][0]['path']=str(protected)
        self.reseal()
        outcome=self.run_report(success=False, extra=extra)
        self.assertIn('must not contain input provenance files', outcome.stderr)


if __name__=='__main__':
    unittest.main()
