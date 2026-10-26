"""Sealed reference evidence: binding, typed missingness, and event unions."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import reference_localization_evidence as reference

NATIVE_SHA=hashlib.sha256(b'immutable native result').hexdigest()


def fixture_event(event_id='event1',pair_id='pair1',supported=True):
    count=int(supported)
    return dict(reference_accession='GCA_example.1',sample_id='sampleA',mode='gapped',
        family='AACCT',role='candidate',declared_motif='TTAGG',event_id=event_id,
        pair_query_id=pair_id,state='completed_independently_audited',
        retained_flank_hit_count=count,retained_context_rows=count,assessed_context_rows=count,
        assessed_flank_context_present=supported,
        full_read_coordinate_supported_assessed_context_rows=count,
        raw_ungated_bridge_associations=count,qualified_bridge_associations=count,
        qualified_bridge_adjacent_reference_continuous_associations=count,
        qualified_bridge_adjacent_reference_acgt_matching_associations=count,
        same_flank_locus_mate_geometry_evaluable=False,
        event_has_same_flank_locus_inward_mate_witness=None,
        smallest_retained_same_flank_locus_inward_five_prime_separation_bp=None,
        pair_has_any_inward_witness_any_retained_locus=None,
        selected_expected_run_max_bp=34 if supported else None,
        selected_expected_run_max_bp_with_full_read_coordinate_support=34 if supported else None,
        expected_orientation_measured_run_max_bp=34 if supported else None,
        adjacent_expected_exact_run_ge500=False if supported else None,
        selected_expected_run_reaches_reference_sequence_edge=False if supported else None,
        flank_search_cap_flag=False,any_event_query_search_cap_flag=False,
        any_ambiguous_query_not_exact_assessable=None,fcs_overlay_assessed=None,
        flank_fcs_overlap=None,any_event_query_fcs_overlap=None,
        any_event_query_published_fcs_actions=None,technical_boundary_warning=False,
        query_selection_scope='conditional_orientation_stratified_saved_events',
        biological_origin_confirmed=None,physical_terminal_confirmed=None,
        main_trm_confirmed=None,novel_trm_confirmed=None)


def fixture_document(root):
    """Portable unsealed input fixture for reporting-integration tests."""
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    source=root/'source_review.json';source.write_text('{"complete":true}\n')
    return dict(schema_version=1,scope=reference.SCOPE,
        source_artifacts=[dict(path=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest())],
        references=[dict(reference_accession='GCA_example.1',sample_id='sampleA',
            repository_taxon='Example species',role='different-isolate comparator',
            sample_linkage='different specimen',verification='fully_verified',
            linkage_limit='No physical-end or specimen-origin confirmation')],
        samples=[dict(sample_id='sampleA',source_result_sha256=NATIVE_SHA,
            scientific_name='Example species',tax_id=123,biosample='SAMN000001',
            selected_event_count=2,selected_pair_count=2,review_notes=['Conditional selected events only'],
            events=[fixture_event(),fixture_event('event2','pair2',False)])])


def fixture_report():
    return dict(sample_id='sampleA',scientific_name='Example species',tax_id='123',
        metadata={'sample_accession':'SAMN000001'},
        provenance={'original_discovery_sha256':NATIVE_SHA},
        candidates=[{'canonical_family':'AACCT','reads':900}],classification='uncatalogued_repeat',
        adequacy={'status':'insufficient_information'})


class ReferenceLocalizationEvidenceTests(unittest.TestCase):
    def setUp(self):
        tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup)
        self.root=Path(tmp.name);self.path=self.root/'evidence.json'
        self.document=fixture_document(self.root);self.report=fixture_report()

    @property
    def sample(self):return self.document['samples'][0]
    @property
    def event(self):return self.sample['events'][0]
    def write(self):
        self.document['payload_sha256']=reference.evidence_payload_sha256(self.document)
        self.path.write_text(json.dumps(self.document))
        return self.path
    def load(self):return reference.load_evidence(self.write())
    def build(self,report=None,sample_id='sampleA'):
        return reference.build_reference_localization_evidence(sample_id,self.report if report is None else report,self.load())

    def test_boundary_is_diagnostic_only_and_preserves_rows(self):
        before=deepcopy(self.report);result=self.build()
        self.assertEqual(self.report,before)
        self.assertEqual(result['status'],'reviewed_subset_boundary_support')
        self.assertEqual(result['selected_event_count'],2)
        self.assertEqual(result['boundary_supported_event_count'],1)
        self.assertEqual(result['events'],self.sample['events'])
        for field in ('native_counts_changed','classification_changed','adequacy_changed'):
            self.assertIs(result[field],False)
        for field in (*reference._CONFIRMATIONS,'library_insert_compatibility'):
            self.assertIsNone(result[field])
        self.assertEqual(set(result['family_summaries'][0]),set(reference.FAMILY_SUMMARY_FIELDS))
        self.assertEqual(set(result['reference_summaries'][0]),set(reference.REFERENCE_SUMMARY_FIELDS))

    def test_reference_and_method_copies_are_event_unions(self):
        ref=deepcopy(self.document['references'][0]);ref['reference_accession']='GCA_other.1'
        self.document['references'].append(ref)
        for e in list(self.sample['events']):
            duplicate=deepcopy(e);duplicate['reference_accession']=ref['reference_accession']
            duplicate['mode']='exact';duplicate['any_ambiguous_query_not_exact_assessable']=False
            self.sample['events'].append(duplicate)
        self.event['any_event_query_search_cap_flag']=True
        result=self.build()
        self.assertEqual(len(result['events']),4)
        self.assertEqual(result['selected_event_count'],2)
        self.assertEqual(result['selected_pair_count'],2)
        self.assertEqual(result['boundary_supported_event_count'],1)
        self.assertEqual(result['capped_event_count'],1)
        self.assertEqual(result['reference_count'],2)
        self.assertEqual(len(result['family_summaries']),1)
        self.assertEqual(len(result['reference_summaries']),2)

    def test_unassessed_has_null_counts_and_empty_lists(self):
        for result in (reference.build_reference_localization_evidence('sampleA',self.report),self.build(sample_id='missing')):
            self.assertEqual(result['status'],'unassessed')
            for field in (*reference._SUMMARY_COUNTS,'reference_count'):self.assertIsNone(result[field])
            for field in ('events','family_summaries','reference_summaries','references','review_notes'):self.assertEqual(result[field],[])

    def test_zero_selected_is_unassessed_not_negative(self):
        self.sample.update(selected_event_count=0,selected_pair_count=0,events=[])
        result=self.build();self.assertEqual(result['reason'],'no_selected_junction_events')
        self.assertIsNone(result['boundary_supported_event_count'])

    def test_no_boundary_is_not_absence_and_pair_witness_is_separate(self):
        self.sample['events']=[fixture_event('event1','pair1',False),fixture_event('event2','pair2',False)]
        self.event['pair_has_any_inward_witness_any_retained_locus']=True
        result=self.build()
        self.assertEqual(result['status'],'reviewed_subset_no_boundary_support')
        self.assertIsNone(result['inward_mate_event_count'])
        self.assertEqual(result['mate_geometry_evaluable_event_count'],0)
        self.assertIsNone(result['family_summaries'][0]['exact_ge500_event_count'])
        self.assertIsNone(result['family_summaries'][0]['max_selected_expected_exact_run_bp'])
        self.assertIn('텔로미어 부재',result['recommendation'])

    def test_mate_witness_does_not_promote_boundary_or_insert_compatibility(self):
        self.sample['events']=[fixture_event('event1','pair1',False),fixture_event('event2','pair2',False)]
        self.event.update(same_flank_locus_mate_geometry_evaluable=True,
            event_has_same_flank_locus_inward_mate_witness=True,
            smallest_retained_same_flank_locus_inward_five_prime_separation_bp=13304)
        result=self.build()
        self.assertEqual(result['inward_mate_event_count'],1)
        self.assertEqual(result['mate_geometry_evaluable_event_count'],1)
        self.assertEqual(result['boundary_supported_event_count'],0)
        self.assertIsNone(result['library_insert_compatibility'])

    def test_evaluated_geometry_without_witness_is_zero_not_null(self):
        self.event.update(same_flank_locus_mate_geometry_evaluable=True,
            event_has_same_flank_locus_inward_mate_witness=False)
        result=self.build()
        self.assertEqual(result['inward_mate_event_count'],0)
        self.assertEqual(result['mate_geometry_evaluable_event_count'],1)

    def test_motif_family_equivalence_uses_primitive_rotation_and_rc(self):
        for motif in ('TTAGG','AGGTT','CCTAA','TTAGGTTAGG'):
            self.event['declared_motif']=motif
            self.sample['events'][1]['declared_motif']=motif
            self.assertEqual(self.build()['family_summaries'][0]['family'],'AACCT')
        original=deepcopy(self.document)
        for field,value in [('declared_motif','TTAGGG'),('family','TTAGG'),('family','AACCTAACCT'),
                            ('declared_motif','A'*257)]:
            self.document=deepcopy(original);self.event[field]=value
            with self.subTest(field=field,value=value),self.assertRaises(reference.EvidenceError):self.load()

    def test_fcs_mate_alternative_is_distinct_from_flank_overlap(self):
        self.event.update(fcs_overlay_assessed=True,flank_fcs_overlap=False,
            any_event_query_fcs_overlap=True,any_event_query_published_fcs_actions=['EXCLUDE'])
        row=self.build()['family_summaries'][0]
        self.assertEqual(row['flank_fcs_overlap_event_count'],0)
        self.assertEqual(row['any_query_fcs_overlap_event_count'],1)
        self.assertEqual(row['published_fcs_actions'],['EXCLUDE'])

    def test_reference_edge_and_exact_run_are_separate_from_terminal(self):
        self.event['selected_expected_run_reaches_reference_sequence_edge']=True
        row=self.build()['family_summaries'][0]
        self.assertEqual(row['reference_edge_event_count'],1)
        self.assertEqual(row['max_selected_expected_exact_run_bp'],34)
        self.assertEqual(row['exact_ge500_event_count'],0)
        self.assertIsNone(self.build()['physical_terminal_confirmed'])

    def test_exact_optional_identity_and_native_sha_binding(self):
        for field,value in [('sample_id','other'),('scientific_name','Other species'),('tax_id','124')]:
            with self.subTest(field=field):
                report=deepcopy(self.report);report[field]=value
                with self.assertRaises(reference.EvidenceError):self.build(report)
        for where in ('native','biosample','family'):
            report=deepcopy(self.report)
            if where=='native':report['provenance']['original_discovery_sha256']='0'*64
            elif where=='biosample':report['metadata']['sample_accession']='SAMN999999'
            else:report['candidates']=[{'canonical_family':'AACCCT'}]
            with self.subTest(where=where),self.assertRaises(reference.EvidenceError):self.build(report)

    def test_minimal_explicit_binding_without_optional_identities(self):
        for k in reference._IDENTITY:self.sample.pop(k)
        result=self.build({'sample_id':'sampleA','source_result_sha256':NATIVE_SHA})
        self.assertEqual(result['selected_event_count'],2)

    def test_raw_evidence_and_missing_native_binding_are_rejected(self):
        with self.assertRaisesRegex(reference.EvidenceError,'load_evidence'):
            reference.build_reference_localization_evidence('sampleA',self.report,self.document)
        report=deepcopy(self.report);report['provenance']={}
        with self.assertRaises(reference.EvidenceError):self.build(report)

    def test_missing_and_changed_source_are_explicit_errors(self):
        self.write();p=Path(self.document['source_artifacts'][0]['path'])
        p.write_text('changed')
        with self.assertRaisesRegex(reference.EvidenceError,'source artifact SHA256 mismatch'):reference.load_evidence(self.path)
        p.unlink()
        with self.assertRaisesRegex(reference.EvidenceError,'cannot read'):reference.load_evidence(self.path)

    def test_bad_seal_and_nonfinite_duplicate_json_are_rejected(self):
        self.write();data=json.loads(self.path.read_text());data['samples'][0]['review_notes']=['unsealed edit']
        self.path.write_text(json.dumps(data))
        with self.assertRaisesRegex(reference.EvidenceError,'seal mismatch'):reference.load_evidence(self.path)
        for content in ('{"schema_version":1,"schema_version":1}','{"schema_version":NaN}'):
            self.path.write_text(content)
            with self.assertRaises(reference.EvidenceError):reference.load_evidence(self.path)

    def test_wrong_schema_scope_and_metadata_are_rejected(self):
        original=deepcopy(self.document)
        changes=[lambda d:d.update(schema_version=True),lambda d:d.update(scope='whole_library'),
            lambda d:d.update(unrecognized=1),lambda d:d['references'][0].update(verification='downloaded'),
            lambda d:d['samples'][0].update(tax_id=True),lambda d:d['samples'][0].update(scientific_name='')]
        for change in changes:
            self.document=deepcopy(original);change(self.document)
            with self.subTest(change=change),self.assertRaises(reference.EvidenceError):self.load()

    def test_duplicate_reference_sample_and_event_fail(self):
        original=deepcopy(self.document)
        for field in ('references','samples','events'):
            self.document=deepcopy(original)
            target=self.sample[field] if field=='events' else self.document[field]
            target.append(deepcopy(target[0]))
            with self.subTest(field=field),self.assertRaises(reference.EvidenceError):self.load()

    def test_same_reference_can_be_reused_by_different_samples(self):
        sample=deepcopy(self.sample);sample['sample_id']='sampleB'
        for e in sample['events']:e['sample_id']='sampleB'
        ref=deepcopy(self.document['references'][0]);ref['sample_id']='sampleB'
        self.document['samples'].append(sample);self.document['references'].append(ref)
        evidence=self.load();self.assertEqual(len(evidence),2)
        result=reference.build_reference_localization_evidence('sampleA',self.report,evidence)
        self.assertEqual(len(result['references']),1)

    def test_event_reference_sample_and_cross_copy_identity_are_validated(self):
        original=deepcopy(self.document)
        for field,value in [('sample_id','sampleB'),('reference_accession','unknown'),('mode','approximate'),('state','running')]:
            self.document=deepcopy(original);self.event[field]=value
            with self.subTest(field=field),self.assertRaises(reference.EvidenceError):self.load()
        self.document=deepcopy(original)
        other=deepcopy(self.event);other.update(mode='exact',any_ambiguous_query_not_exact_assessable=False,pair_query_id='different')
        self.sample['events'].append(other)
        with self.assertRaisesRegex(reference.EvidenceError,'event identity'):self.load()

    def test_count_type_hierarchy_and_selected_denominator_are_validated(self):
        original=deepcopy(self.document)
        for field,value in [('selected_event_count',3),('selected_pair_count',True),('selected_event_count',-1)]:
            self.document=deepcopy(original);self.sample[field]=value
            with self.subTest(field=field,value=value),self.assertRaises(reference.EvidenceError):self.load()
        for field,value in [('assessed_context_rows',True),('retained_context_rows',-1),
                ('qualified_bridge_associations',2),('qualified_bridge_adjacent_reference_acgt_matching_associations',2)]:
            self.document=deepcopy(original);self.event[field]=value
            with self.subTest(field=field),self.assertRaises(reference.EvidenceError):self.load()

    def test_null_counts_runs_fcs_and_biological_promotions_fail(self):
        original=deepcopy(self.document)
        changes=[lambda d:d['samples'][0]['events'][1].update(adjacent_expected_exact_run_ge500=False),
            lambda d:d['samples'][0]['events'][0].update(expected_orientation_measured_run_max_bp=True),
            lambda d:d['samples'][0]['events'][0].update(flank_fcs_overlap=False),
            lambda d:d['samples'][0]['events'][0].update(physical_terminal_confirmed=False),
            lambda d:d['samples'][0]['events'][0].update(adjacent_expected_exact_run_ge500=True),
            lambda d:d['samples'][0]['events'][0].update(event_has_same_flank_locus_inward_mate_witness=False),
            lambda d:d['samples'][0]['events'][0].update(flank_search_cap_flag=True),
            lambda d:d['samples'][0]['events'][0].pop('technical_boundary_warning')]
        for change in changes:
            self.document=deepcopy(original);change(self.document)
            with self.subTest(change=change),self.assertRaises(reference.EvidenceError):self.load()

    def test_unknown_event_field_fails_closed(self):
        self.event['terminal_confirmed']=True
        with self.assertRaises(reference.EvidenceError):self.load()

    def test_source_dependencies_and_relative_pins(self):
        self.document['source_artifacts'][0]['path']='source_review.json'
        path=self.write();deps=reference.source_dependencies(path)
        self.assertEqual(deps,[str(path),str(self.root/'source_review.json')])
        self.assertEqual(len(reference.load_evidence(path)),1)

    def test_oversize_sidecar_or_source_and_directory_are_rejected(self):
        self.write()
        with patch.object(reference,'MAX_ARTIFACT_BYTES',10):
            with self.assertRaisesRegex(reference.EvidenceError,'small regular file'):reference.load_evidence(self.path)
        self.document['source_artifacts'][0]['path']=str(self.root)
        with self.assertRaisesRegex(reference.EvidenceError,'small regular file'):self.load()

    def test_copy_isolation(self):
        evidence=self.load();sample=evidence['sampleA'];sample['events'][0]['qualified_bridge_associations']=100
        refs=evidence.references;refs[0]['verification']='altered'
        result=reference.build_reference_localization_evidence('sampleA',self.report,evidence)
        self.assertEqual(result['boundary_supported_event_count'],1)
        self.assertEqual(result['references'][0]['verification'],'fully_verified')
        result['events'][0]['qualified_bridge_associations']=100
        self.assertEqual(evidence['sampleA']['events'][0]['qualified_bridge_associations'],1)


if __name__=='__main__':unittest.main()
