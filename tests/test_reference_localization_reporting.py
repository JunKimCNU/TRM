"""Offline report/cache integration of optional saved-junction localization."""
from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'tests')]
import reference_localization_evidence as evidence_module
import report_runtime
from summarize_interpretation import summarize_sample
import test_report_assessment as fixture


TABLES = ('reference_localization_diagnostics.tsv',
          'reference_localization_families.tsv',
          'reference_localization_references.tsv',
          'reference_localization_events.tsv')
COUNTS = ('selected_event_count', 'selected_pair_count', 'assessed_event_count',
          'coordinate_supported_event_count', 'boundary_supported_event_count',
          'mate_geometry_evaluable_event_count', 'inward_mate_event_count',
          'capped_event_count', 'reference_count')


class ReferenceLocalizationReportingTests(unittest.TestCase):
    setUp = fixture.AssessmentReportTests.setUp
    command = fixture.AssessmentReportTests.command
    run_report = fixture.AssessmentReportTests.run_report
    document = fixture.AssessmentReportTests.document

    def save_evidence(self):
        self.pin = self.root / 'saved-reference-review.txt'
        self.pin.write_text('Audited selected-event fixture, no raw sequence input.\n')
        self.sidecar = self.root / 'localization.json'
        event = dict(reference_accession='GCA_000001.1', sample_id='SRRknown',
            mode='exact', family='AAGCCT', role='candidate', declared_motif='TTAGGC',
            event_id='junction1', pair_query_id='pair1',
            state='completed_independently_audited',
            retained_flank_hit_count=2, retained_context_rows=2, assessed_context_rows=2,
            assessed_flank_context_present=True,
            full_read_coordinate_supported_assessed_context_rows=1,
            raw_ungated_bridge_associations=2, qualified_bridge_associations=1,
            qualified_bridge_adjacent_reference_continuous_associations=1,
            qualified_bridge_adjacent_reference_acgt_matching_associations=1,
            same_flank_locus_mate_geometry_evaluable=True,
            event_has_same_flank_locus_inward_mate_witness=False,
            smallest_retained_same_flank_locus_inward_five_prime_separation_bp=None,
            pair_has_any_inward_witness_any_retained_locus=True,
            selected_expected_run_max_bp=120,
            selected_expected_run_max_bp_with_full_read_coordinate_support=120,
            expected_orientation_measured_run_max_bp=120,
            adjacent_expected_exact_run_ge500=False,
            selected_expected_run_reaches_reference_sequence_edge=False,
            flank_search_cap_flag=False, any_event_query_search_cap_flag=True,
            any_ambiguous_query_not_exact_assessable=False,
            fcs_overlay_assessed=None, flank_fcs_overlap=None,
            any_event_query_fcs_overlap=None,
            any_event_query_published_fcs_actions=None,
            technical_boundary_warning=True,
            query_selection_scope='conditional_orientation_stratified_saved_events',
            biological_origin_confirmed=None, physical_terminal_confirmed=None,
            main_trm_confirmed=None, novel_trm_confirmed=None)
        # Four rows are one selected event, not four observations or samples.
        events = []
        references = []
        for accession in ('GCA_000001.1', 'GCA_000002.1'):
            references.append(dict(reference_accession=accession, sample_id='SRRknown',
                repository_taxon='Fixture reference', role='assembly_context',
                sample_linkage='Different specimen', verification='fully_verified',
                linkage_limit='Assembly edge is not a physical terminus.'))
            for mode in ('exact', 'gapped'):
                row = deepcopy(event)
                row.update(reference_accession=accession, mode=mode,
                    any_ambiguous_query_not_exact_assessable=False if mode == 'exact' else None)
                events.append(row)
        self.payload = dict(schema_version=1, scope=evidence_module.SCOPE,
            source_artifacts=[dict(path=str(self.pin), sha256=fixture.digest(self.pin))],
            references=references,
            samples=[dict(sample_id='SRRknown',
                source_result_sha256=fixture.digest(self.source/'SRRknown/result.json'),
                scientific_name='Fixture species', tax_id='6239',
                biosample='SAMNfixture_SRRknown', selected_event_count=1,
                selected_pair_count=1, review_notes=['No independent physical-terminal confirmation.'],
                events=events)])
        self.reseal()
        return ['--reference-localization-evidence', str(self.sidecar)]

    def reseal(self):
        self.payload['payload_sha256'] = evidence_module.evidence_payload_sha256(self.payload)
        fixture.write_json(self.sidecar, self.payload)

    def resources(self):
        return {row['sample_id']: row for row in
                fixture.read_tsv(self.output/'interpretation_resources.tsv')}

    def test_overlay_reaches_json_html_tables_without_changing_science_or_summary(self):
        self.run_report()
        previous = self.document('SRRknown')
        candidates = (self.output/'candidates.tsv').read_bytes()
        self.run_report(extra=self.save_evidence())
        current = self.document('SRRknown')
        for key in previous:
            if key != 'reference_localization_diagnostics':
                self.assertEqual(previous[key], current[key], key)
        self.assertEqual(candidates, (self.output/'candidates.tsv').read_bytes())
        # Projection rows contain only the established scientific summary fields.
        self.assertEqual(summarize_sample(previous, Path('fixed.json'), '0'*64),
                         summarize_sample(current, Path('fixed.json'), '0'*64))
        diagnostic = current['reference_localization_diagnostics']
        self.assertEqual(diagnostic['status'], 'reviewed_subset_boundary_support')
        self.assertEqual(diagnostic['selected_event_count'], 1)
        self.assertEqual(diagnostic['selected_pair_count'], 1)
        self.assertEqual(diagnostic['reference_count'], 2)
        self.assertEqual(diagnostic['boundary_supported_event_count'], 1)
        self.assertEqual(diagnostic['mate_geometry_evaluable_event_count'], 1)
        self.assertEqual(diagnostic['inward_mate_event_count'], 0)
        self.assertEqual(diagnostic['capped_event_count'], 1)
        self.assertEqual(diagnostic['events'], self.payload['samples'][0]['events'])
        families = fixture.read_tsv(self.output/TABLES[1])
        self.assertEqual(len(families), 1)
        self.assertEqual(families[0]['selected_event_count'], '1')
        self.assertEqual(families[0]['flank_fcs_overlap_event_count'], '')
        self.assertEqual(families[0]['exact_ge500_event_count'], '0')
        self.assertEqual(len(fixture.read_tsv(self.output/TABLES[2])), 4)
        self.assertEqual(len(fixture.read_tsv(self.output/TABLES[3])), 4)
        text = (self.output/'samples/SRRknown.html').read_text()
        self.assertIn('선정 경계의 reference 위치 비교', text)
        self.assertIn('<caption>Family별', text)
        self.assertIn('행 사이 합산 금지', text)
        for name in TABLES:
            self.assertIn('../'+name, text)
        row = next(r for r in fixture.read_tsv(self.output/'samples.tsv') if r['sample_id']=='SRRknown')
        self.assertEqual(row['reference_localization_boundary_supported_event_count'], '1')

    def test_omitted_and_subset_missing_evidence_are_null_not_zero(self):
        self.run_report()
        for sid in ('SRRknown', 'SRRnone', 'SRRfailed'):
            diagnostic = self.document(sid)['reference_localization_diagnostics']
            self.assertEqual(diagnostic['status'], 'unassessed')
            self.assertTrue(all(diagnostic[k] is None for k in COUNTS))
        for name in TABLES[1:]:
            self.assertEqual(fixture.read_tsv(self.output/name), [])
        self.run_report(extra=self.save_evidence())
        for sid in ('SRRnovel', 'SRRnone', 'SRRfailed'):
            diagnostic = self.document(sid)['reference_localization_diagnostics']
            self.assertEqual(diagnostic['reason'], 'sample_not_in_evidence')
            self.assertTrue(all(diagnostic[k] is None for k in COUNTS))

    def test_coordinate_support_without_boundary_does_not_promote_or_confirm(self):
        extra = self.save_evidence()
        for event in self.payload['samples'][0]['events']:
            event.update(qualified_bridge_associations=0,
                qualified_bridge_adjacent_reference_continuous_associations=0,
                qualified_bridge_adjacent_reference_acgt_matching_associations=0)
        self.reseal()
        self.run_report(extra=extra)
        diagnostic = self.document('SRRknown')['reference_localization_diagnostics']
        self.assertEqual(diagnostic['status'], 'reviewed_subset_no_boundary_support')
        self.assertEqual(diagnostic['coordinate_supported_event_count'], 1)
        self.assertEqual(diagnostic['boundary_supported_event_count'], 0)
        self.assertEqual(diagnostic['inward_mate_event_count'], 0)
        for key in ('biological_origin_confirmed', 'physical_terminal_confirmed',
                    'main_trm_confirmed', 'novel_trm_confirmed', 'library_insert_compatibility'):
            self.assertIsNone(diagnostic[key])

    def test_reference_and_review_notes_are_html_escaped(self):
        extra = self.save_evidence()
        hostile = '<script>alert("fixture")</script>'
        self.payload['samples'][0]['review_notes'] = [hostile]
        self.payload['references'][0]['linkage_limit'] = hostile
        self.payload['references'][0]['repository_taxon'] = hostile
        self.reseal()
        self.run_report(extra=extra)
        text = (self.output/'samples/SRRknown.html').read_text()
        self.assertNotIn(hostile, text)
        self.assertIn('&lt;script&gt;alert(&quot;fixture&quot;)&lt;/script&gt;', text)
        self.assertIn(hostile, self.document('SRRknown')['reference_localization_diagnostics']['review_notes'])

    def test_unassessed_same_locus_geometry_remains_null_despite_broad_pair_witness(self):
        extra = self.save_evidence()
        for event in self.payload['samples'][0]['events']:
            event.update(same_flank_locus_mate_geometry_evaluable=False,
                event_has_same_flank_locus_inward_mate_witness=None)
        self.reseal()
        self.run_report(extra=extra)
        diagnostic = self.document('SRRknown')['reference_localization_diagnostics']
        self.assertEqual(diagnostic['mate_geometry_evaluable_event_count'], 0)
        self.assertIsNone(diagnostic['inward_mate_event_count'])
        row = next(r for r in fixture.read_tsv(self.output/TABLES[0]) if r['sample_id']=='SRRknown')
        self.assertEqual(row['mate_geometry_evaluable_event_count'], '0')
        self.assertEqual(row['inward_mate_event_count'], '')
        self.assertIn('같은 위치 mate 평가', (self.output/'samples/SRRknown.html').read_text())

    def test_stale_native_binding_fails_without_completed_report(self):
        extra = self.save_evidence()
        self.payload['samples'][0]['source_result_sha256'] = '0'*64
        self.reseal()
        self.run_report(success=False, extra=extra)
        self.assertEqual(json.loads((self.output/'interpretation_status.json').read_text())['status'], 'failed')

    def test_missing_explicit_sidecar_fails_preflight(self):
        result = self.run_report(success=False,
            extra=['--reference-localization-evidence', str(self.root/'missing.json')])
        self.assertIn('cannot read', result.stderr)
        self.assertFalse((self.output/'metrics.json').exists())

    def test_source_change_invalidates_warm_cache_and_fails_closed(self):
        extra = self.save_evidence()+['--resume', '--jobs', '2']
        self.run_report(extra=extra)
        documents = {p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')}
        self.run_report(extra=extra)
        self.assertEqual(documents, {p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')})
        self.assertTrue(all(r['mode']=='reused' for r in self.resources().values()))
        self.pin.write_text('Changed independently saved source; seal no longer binds it.\n')
        self.run_report(success=False, extra=extra)
        self.assertEqual(json.loads((self.output/'interpretation_status.json').read_text())['status'], 'failed')

    def test_option_removal_does_not_reuse_or_leave_old_event_rows(self):
        self.run_report(extra=self.save_evidence()+['--resume'])
        self.run_report(extra=['--resume'])
        self.assertEqual(self.document('SRRknown')['reference_localization_diagnostics']['status'], 'unassessed')
        self.assertTrue(all(r['mode']=='recomputed' for r in self.resources().values()))
        for name in TABLES[1:]:
            self.assertEqual(fixture.read_tsv(self.output/name), [])

    def test_sequential_and_parallel_are_identical_including_new_outputs(self):
        extra = self.save_evidence()
        self.run_report(extra=extra+['--jobs', '1'])
        before = {p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')}
        tables = {name:(self.output/name).read_bytes() for name in TABLES}
        self.run_report(extra=extra+['--jobs', '2'])
        self.assertEqual(before, {p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')})
        self.assertEqual(tables, {name:(self.output/name).read_bytes() for name in TABLES})

    def test_output_cannot_contain_transitive_small_source(self):
        extra = self.save_evidence()
        protected = self.output/'protected-review.txt'
        protected.parent.mkdir()
        protected.write_bytes(self.pin.read_bytes())
        self.payload['source_artifacts'][0]['path'] = str(protected)
        self.reseal()
        outcome = self.run_report(success=False, extra=extra)
        self.assertIn('must not contain input provenance files', outcome.stderr)
        self.assertEqual(protected.read_bytes(), self.pin.read_bytes())

    def test_required_localization_chunk_cannot_be_removed_from_sealed_receipt(self):
        self.run_report(extra=self.save_evidence())
        receipt_path = self.output/'.report-cache/SRRknown/receipt.json'
        receipt = json.loads(receipt_path.read_text())
        self.assertTrue(report_runtime.verify_artifacts(receipt, output=self.output, sample='SRRknown'))
        bundle_path = Path(receipt['bundle'])
        bundle = json.loads(bundle_path.read_text())
        missing = TABLES[3]
        bundle['chunks'].remove(missing)
        fixture.write_json(bundle_path, bundle)
        receipt['artifacts'][str(bundle_path)] = report_runtime.binding(bundle_path)[0]
        receipt['artifacts'].pop(str(bundle_path.parent/missing))
        # Even a self-consistent malicious/incomplete artifact set is rejected.
        self.assertFalse(report_runtime.verify_artifacts(receipt, output=self.output, sample='SRRknown'))


if __name__ == '__main__':
    unittest.main()
