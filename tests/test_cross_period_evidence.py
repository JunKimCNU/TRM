"""Structural integrity, specimen binding and same-context reporting tests."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import cross_period_evidence as module


def make_document(root):
    """Small valid fixture; caller may change sample/native identity then reseal."""
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    source=root/'cross_period_fixture_source.json'
    source.write_text('{"audited_saved_observations":true}\n')
    fcs=dict(assessed=False,span_overlap=None,aligned_base_overlap_bp=None,published_actions=None)
    ops=lambda **kw:{k:dict(min=kw.get(k,0),max=kw.get(k,0)) for k in ('match','substitution','insertion','deletion')}
    interval=dict(interval_id='i1',pair_query_id='pair1',mate=1,start0=0,end0=6,full_read_length=12,
        sequence='TTAGGC',quality='IIIIII',candidate_cost=0,comparator_cost=1,
        candidate_operations=ops(match=6),comparator_operations=ops(match=5,insertion=1),
        candidate_ambiguous_operation_positions=0,comparator_ambiguous_operation_positions=0,deletion_quality=None,
        strict_contrasts=[dict(absolute_position0=5,query_base='C',contrast='candidate_always_match_comparator_always_observed_edit',
            quality_character='I',quality_ascii=73,conditional_phred33=40,candidate_operations=['match'],comparator_operations=['insertion'],
            candidate_template_bases=['C'],comparator_template_bases=[])])
    run=dict(complete_periods=1,contiguous_exact_bp=6,contiguous_exact_repeat_ge_500=False,start0=100,end0=106,
        family_sequence_in_reference_forward_frame='AAGCCT',matching_phase_indices_at_boundary_base=[0],
        run_reaches_context_limit=False,run_reaches_sequence_edge=False,sequence_beyond_run_unresolved=False,stop_position0=99,stop_reason='mismatch')
    context=dict(joint_id='j1',event_id='e1',reference_accession='REF1',mode='gapped',target='chr1',
        flank_hsp_id='f1',full_read_hsp_id='h1',full_read_retained_hit_index0=0,full_read_strand='+',
        full_read_query_start0=0,full_read_query_end0=12,full_read_reference_start0=100,full_read_reference_end0=112,
        boundary_read0=6,boundary_reference0=106,exact_flank_hit=None,selected_expected_run=run,
        full_hsp_observed_interval_bases=6,strict_positions_in_full_hsp=[5],mate_geometry_evaluable=False,
        same_hsp_inward_witness=None,closest_mate=None,boundary_reference_continuous=True,boundary_adjacent_acgt_matches=True,
        technical_boundary_warning=False,aligned_technical_overlap={'12':0,'16':0,'20':0},
        compartment=dict(assembly_unit='Primary Assembly',assigned_molecule='1',assigned_molecule_location_type='Chromosome',reference_sequence_bp=1000,sequence_role='assembled-molecule'),
        flank_fcs=copy.deepcopy(fcs),full_read_fcs=copy.deepcopy(fcs),mate_fcs=None,search_capped=True)
    sample=dict(sample_id='SRRknown',source_result_sha256='a'*64,scientific_name='Test species',tax_id='123',biosample='SAMN1',
        motifs={'candidate':'TTAGGC','comparator':'TTAGG'},review_notes=['Synthetic structural fixture, not biological evidence.'],
        intervals=[interval],events=[dict(event_id='e1',interval_id='i1',family='AAGCCT',role='candidate',declared_motif='TTAGGC')],contexts=[context])
    doc=dict(schema_version=module.SCHEMA_VERSION,scope=module.SCOPE,
        source_artifacts=[dict(path=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest())],samples=[sample],payload_sha256='')
    doc['payload_sha256']=module.evidence_payload_sha256(doc)
    return doc


def minimal_report(document):
    s=document['samples'][0]
    return dict(sample_id=s['sample_id'],source_result_sha256=s['source_result_sha256'],scientific_name=s['scientific_name'],
        tax_id=s['tax_id'],metadata={'sample_accession':s['biosample']},candidates=[{'canonical_family':e['family']} for e in s['events']])


class CrossPeriodEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.doc=make_document(self.root);self.path=self.root/'evidence.json'
    def write(self,doc=None,reseal=True):
        d=copy.deepcopy(self.doc if doc is None else doc)
        if reseal:d['payload_sha256']=module.evidence_payload_sha256(d)
        self.path.write_text(json.dumps(d));return module.load_evidence(self.path)
    def output(self,doc=None):
        d=self.doc if doc is None else doc
        return module.build_cross_period_evidence(d['samples'][0]['sample_id'],minimal_report(d),self.write(d))
    def mutate_bad(self,mutation):
        d=copy.deepcopy(self.doc);mutation(d['samples'][0])
        with self.assertRaises(module.EvidenceError):self.write(d)
    def test_baseline_distinguishes_missing_and_measured_zero(self):
        out=self.output()
        self.assertEqual(out['strict_contrast_position_count'],1)
        self.assertEqual(out['joint_repeat_inward_event_count'],0)
        self.assertIsNone(out['inward_mate_event_count'])
        self.assertEqual(out['selected_read_count'],1)
        self.assertEqual(out['candidate_lower_cost_interval_count'],1)
        self.assertFalse(out['quality_encoding_verified'])
        for field in module.SUMMARY_FIELDS:
            self.assertIsNone(module.build_cross_period_evidence('x',{})[field])
        missing=module.build_cross_period_evidence('other',{},self.write())
        self.assertEqual(missing['reason'],'sample_not_in_evidence')
        self.assertTrue(all(missing[k] is None for k in module.SUMMARY_FIELDS))
    def test_empty_selected_sample_has_known_zero_denominators(self):
        d=copy.deepcopy(self.doc);d['samples'][0].update(events=[],intervals=[],contexts=[])
        out=self.output(d)
        self.assertEqual(out['selected_event_count'],0);self.assertIsNone(out['inward_mate_event_count'])
        self.assertEqual(out['reason'],'no_selected_junction_events')
    def test_seal_and_source_integrity(self):
        d=copy.deepcopy(self.doc);d['samples'][0]['review_notes'].append('changed')
        with self.assertRaisesRegex(module.EvidenceError,'seal'):self.write(d,False)
        self.write();Path(self.doc['source_artifacts'][0]['path']).write_text('mutated')
        with self.assertRaisesRegex(module.EvidenceError,'SHA256'):module.load_evidence(self.path)
    def test_duplicate_and_nonfinite_json_rejected(self):
        for text in ('{"samples":[],"samples":[]}', '{"x":NaN}', '{"x":Infinity}'):
            self.path.write_text(text)
            with self.assertRaises(module.EvidenceError):module.load_evidence(self.path)
    def test_source_bounds_and_directory(self):
        with mock.patch.object(module,'MAX_ARTIFACT_BYTES',10):
            with self.assertRaises(module.EvidenceError):self.write()
        d=copy.deepcopy(self.doc);d['source_artifacts'][0]['path']=str(self.root)
        with self.assertRaises(module.EvidenceError):self.write(d)
    def test_relative_dependencies_and_duplicate_path_alias(self):
        d=copy.deepcopy(self.doc);d['source_artifacts'][0]['path']='cross_period_fixture_source.json'
        self.write(d);deps=module.source_dependencies(self.path)
        self.assertEqual(deps,[str(self.path),str(self.root/'cross_period_fixture_source.json')])
        d['source_artifacts'].append(dict(d['source_artifacts'][0],path=str(self.root/'cross_period_fixture_source.json')))
        with self.assertRaises(module.EvidenceError):self.write(d)
    def test_stat_change_during_source_read_fails(self):
        original=module.os.fstat;calls=[]
        def changing(fd):
            s=original(fd);calls.append(1)
            if len(calls)==2:
                class Changed:pass
                changed=Changed()
                for k in ('st_dev','st_ino','st_size','st_mtime_ns','st_ctime_ns','st_mode'):setattr(changed,k,getattr(s,k))
                changed.st_mtime_ns+=1
                return changed
            return s
        with mock.patch.object(module.os,'fstat',side_effect=changing):
            with self.assertRaisesRegex(module.EvidenceError,'changed while reading'):module._read_bytes(self.doc['source_artifacts'][0]['path'])
    def test_native_and_specimen_binding(self):
        evidence=self.write();report=minimal_report(self.doc)
        for field,value in [('sample_id','other'),('source_result_sha256','b'*64),('scientific_name','Other'),('tax_id','124')]:
            bad=copy.deepcopy(report);bad[field]=value
            with self.assertRaises(module.EvidenceError):module.build_cross_period_evidence('SRRknown',bad,evidence)
        bad=copy.deepcopy(report);bad['metadata']['sample_accession']='other'
        with self.assertRaises(module.EvidenceError):module.build_cross_period_evidence('SRRknown',bad,evidence)
        # Reported provenance takes precedence over an independently supplied shortcut.
        bad=copy.deepcopy(report);bad['provenance']={'original_discovery_sha256':'b'*64}
        with self.assertRaises(module.EvidenceError):module.build_cross_period_evidence('SRRknown',bad,evidence)
    def test_only_observed_event_families_need_membership(self):
        e=self.write();r=minimal_report(self.doc)
        self.assertEqual(len(r['candidates']),1)
        self.assertEqual(module.build_cross_period_evidence('SRRknown',r,e)['selected_event_count'],1)
        r['candidates']=[]
        with self.assertRaises(module.EvidenceError):module.build_cross_period_evidence('SRRknown',r,e)
        with self.assertRaises(module.EvidenceError):module.build_cross_period_evidence('SRRknown',r,{'SRRknown':self.doc['samples'][0]})
    def test_copy_isolation_and_no_input_mutation(self):
        e=self.write();r=minimal_report(self.doc);before=copy.deepcopy(r)
        out=module.build_cross_period_evidence('SRRknown',r,e)
        out['intervals'][0]['candidate_operations']['match']['min']=999
        e['SRRknown']['intervals'][0]['sequence']='X'
        e.provenance['source_artifacts'][0]['sha256']='X'
        fresh=module.build_cross_period_evidence('SRRknown',r,e)
        self.assertEqual(fresh['intervals'][0]['candidate_operations']['match']['min'],6)
        self.assertEqual(r,before);self.assertEqual(e['SRRknown']['intervals'][0]['sequence'],'TTAGGC')
    def test_bool_integer_noninteger_invalid_coordinates(self):
        for key,value in [('mate',True),('mate',3),('start0',False),('end0',0),('full_read_length',5),('candidate_cost',1.0),('candidate_cost',-1),('comparator_ambiguous_operation_positions',7)]:
            with self.subTest(key=key,value=value):self.mutate_bad(lambda s:s['intervals'][0].__setitem__(key,value))
    def test_operation_ranges_and_deletion_quality(self):
        for change in [lambda s:s['intervals'][0].__setitem__('deletion_quality',40),
                       lambda s:s['intervals'][0]['candidate_operations']['match'].__setitem__('min',7),
                       lambda s:s['intervals'][0]['comparator_operations']['insertion'].__setitem__('max',2),
                       lambda s:s['intervals'][0]['candidate_operations'].__setitem__('unknown',{'min':0,'max':0})]:self.mutate_bad(change)
    def test_ambiguous_ranges_not_combined_as_one_path(self):
        d=copy.deepcopy(self.doc);i=d['samples'][0]['intervals'][0]
        # Necessary range checks accept alternative complete paths with either
        # substitution or insertion; no independent minima are called a path.
        i['comparator_operations']['substitution']={'min':0,'max':1}
        i['comparator_operations']['insertion']={'min':0,'max':1}
        i['strict_contrasts'][0]['comparator_operations']=['substitution','insertion']
        i['strict_contrasts'][0]['comparator_template_bases']=['A']
        self.output(d)
    def test_motif_family_primitive_period(self):
        for fn in [lambda s:s['motifs'].__setitem__('candidate','TTAGGN'),
                   lambda s:s['motifs'].__setitem__('comparator','TTAGGC'),
                   lambda s:s['events'][0].__setitem__('family','TTAGGC'),
                   lambda s:s['events'][0].__setitem__('declared_motif','TTAGG'),
                   lambda s:s['events'][0].__setitem__('role','unknown')]:self.mutate_bad(fn)
    def test_strict_position_base_quality_operation_validation(self):
        for k,v in [('absolute_position0',6),('quality_ascii',True),('query_base','G'),('quality_character','H'),('conditional_phred33',39),
                    ('candidate_operations',['match','insertion']),('comparator_operations',['deletion']),('candidate_template_bases',['A']),('contrast','unknown')]:
            with self.subTest(k=k):self.mutate_bad(lambda s:s['intervals'][0]['strict_contrasts'][0].__setitem__(k,v))
        self.mutate_bad(lambda s:s['intervals'][0]['strict_contrasts'].append(copy.deepcopy(s['intervals'][0]['strict_contrasts'][0])))
    def test_duplicate_identities_and_unknown_links(self):
        self.mutate_bad(lambda s:s['events'].append(copy.deepcopy(s['events'][0])))
        self.mutate_bad(lambda s:s['events'].append(dict(s['events'][0],event_id='other')))
        self.mutate_bad(lambda s:s['intervals'].append(dict(s['intervals'][0],interval_id='other')))
        self.mutate_bad(lambda s:s['contexts'].append(dict(s['contexts'][0],joint_id='other')))
        self.mutate_bad(lambda s:s['contexts'][0].__setitem__('event_id','missing'))
    def test_geometry_null_and_actual_query_coverage(self):
        for k,v in [('same_hsp_inward_witness',False),('full_hsp_observed_interval_bases',7),('strict_positions_in_full_hsp',[4]),('strict_positions_in_full_hsp',[5,5]),('boundary_read0',2),('boundary_reference0',120)]:
            self.mutate_bad(lambda s:s['contexts'][0].__setitem__(k,v))
        d=copy.deepcopy(self.doc);d['samples'][0]['contexts'][0].update(mate_geometry_evaluable=True,same_hsp_inward_witness=False)
        self.assertEqual(self.output(d)['inward_mate_event_count'],0)
    @staticmethod
    def inward(context):
        context.update(mate_geometry_evaluable=True,same_hsp_inward_witness=True,
            closest_mate=dict(five_prime_boundary0=200,five_prime_separation_bp=100,hsp_id='m1',library_insert_compatibility=None,
                retained_hit_index0=0,strand='-',support_five_prime_boundary0=100,target='chr1'))
    def test_same_context_not_independent_maxima(self):
        d=copy.deepcopy(self.doc);s=d['samples'][0]
        other=copy.deepcopy(s['contexts'][0]);other.update(joint_id='j2',full_read_hsp_id='h2',full_read_retained_hit_index0=1,selected_expected_run=None)
        self.inward(other);s['contexts'].append(other)
        out=self.output(d)
        self.assertEqual(out['inward_mate_event_count'],1);self.assertEqual(out['joint_repeat_inward_event_count'],0)
        self.assertIsNone(out['family_summaries'][0]['max_joint_adjacent_reference_run_bp'])
        self.inward(s['contexts'][0]);out=self.output(d)
        self.assertEqual(out['joint_repeat_inward_event_count'],1)
        self.assertEqual(out['family_summaries'][0]['max_joint_adjacent_reference_run_bp'],6)
        self.assertEqual(out['selected_event_count'],1);self.assertEqual(out['qualified_context_count'],2)
    def test_mate_from_another_hsp_or_target_rejected(self):
        for key,value in [('target','chr2'),('support_five_prime_boundary0',101),('strand','+'),('five_prime_separation_bp',99),('library_insert_compatibility',True)]:
            d=copy.deepcopy(self.doc);self.inward(d['samples'][0]['contexts'][0]);d['samples'][0]['contexts'][0]['closest_mate'][key]=value
            with self.assertRaises(module.EvidenceError):self.write(d)
    def test_reference_method_copies_not_added_to_event_denominator(self):
        d=copy.deepcopy(self.doc);s=d['samples'][0]
        copy_context=copy.deepcopy(s['contexts'][0]);copy_context.update(joint_id='other',reference_accession='REF2')
        s['contexts'].append(copy_context);out=self.output(d)
        self.assertEqual(out['selected_event_count'],1);self.assertEqual(out['selected_interval_count'],1)
        self.assertEqual(out['boundary_supported_event_count'],1);self.assertEqual(out['qualified_context_count'],2)
        self.assertEqual(out['strict_contrast_position_count'],1)
    def test_adjacent_run_not_equal_to_observed_coverage(self):
        d=copy.deepcopy(self.doc);c=d['samples'][0]['contexts'][0]
        c['selected_expected_run'].update(start0=76,contiguous_exact_bp=30,complete_periods=5,stop_position0=75)
        out=self.output(d)
        self.assertEqual(out['contexts'][0]['selected_adjacent_reference_run_bp'],30)
        self.assertEqual(out['contexts'][0]['full_hsp_observed_interval_bases'],6)
    def test_exact_flank_normalized_strand_may_differ_from_read(self):
        d=copy.deepcopy(self.doc);c=d['samples'][0]['contexts'][0]
        c.update(mode='exact',flank_hsp_id=None,full_read_hsp_id=None,exact_flank_hit=dict(target='chr1',start0=106,end0=112,strand='-'))
        self.output(d)
    def test_fcs_unknown_must_not_be_false_clean_result(self):
        self.mutate_bad(lambda s:s['contexts'][0]['flank_fcs'].__setitem__('span_overlap',False))
        self.mutate_bad(lambda s:s['contexts'][0]['full_read_fcs'].__setitem__('published_actions',[]))
    def test_endpoint_touch_is_not_boundary_crossing(self):
        self.mutate_bad(lambda s:s['contexts'][0].__setitem__('full_read_query_end0',6))
        self.mutate_bad(lambda s:s['contexts'][0].__setitem__('full_read_query_start0',6))
    def test_strict_operation_template_consistency(self):
        for key,value in [('comparator_template_bases',['A']),('comparator_template_bases',['C']),('candidate_template_bases',['A','C'])]:
            self.mutate_bad(lambda s:s['intervals'][0]['strict_contrasts'][0].__setitem__(key,value))
        d=copy.deepcopy(self.doc);p=d['samples'][0]['intervals'][0]['strict_contrasts'][0]
        p.update(comparator_operations=['substitution'],comparator_template_bases=['C'])
        with self.assertRaises(module.EvidenceError):self.write(d)
    def test_unhashable_nested_values_fail_with_evidence_error(self):
        self.mutate_bad(lambda s:s['intervals'][0]['strict_contrasts'][0].__setitem__('candidate_operations',[{}]))
    def test_position_operation_requires_global_path_support(self):
        d=copy.deepcopy(self.doc);p=d['samples'][0]['intervals'][0]['strict_contrasts'][0]
        p.update(comparator_operations=['substitution'],comparator_template_bases=['A'])
        with self.assertRaisesRegex(module.EvidenceError,'complete optimal path'):self.write(d)
    def test_overlapping_physical_read_bases_must_agree(self):
        d=copy.deepcopy(self.doc);s=d['samples'][0]
        second=copy.deepcopy(s['intervals'][0]);second.update(interval_id='i2',start0=1,end0=7)
        second['strict_contrasts'][0]['absolute_position0']=6
        s['intervals'].append(second);s['events'].append(dict(s['events'][0],event_id='e2',interval_id='i2'))
        with self.assertRaisesRegex(module.EvidenceError,'overlapping'):self.write(d)
    def test_source_count_and_total_byte_bounds(self):
        with mock.patch.object(module,'MAX_SOURCE_ARTIFACTS',0):
            with self.assertRaises(module.EvidenceError):self.write()
        with mock.patch.object(module,'MAX_SOURCE_BYTES',1):
            with self.assertRaises(module.EvidenceError):self.write()
    def test_all_read_edges_and_strands_require_correct_run_side(self):
        for read_left in (True,False):
            for strand in ('+','-'):
                for mode in ('gapped','exact'):
                    with self.subTest(read_left=read_left,strand=strand,mode=mode):
                        d=copy.deepcopy(self.doc);s=d['samples'][0];i=s['intervals'][0];c=s['contexts'][0]
                        if not read_left:
                            i.update(start0=6,end0=12)
                            i['strict_contrasts'][0]['absolute_position0']=11
                            c['strict_positions_in_full_hsp']=[11]
                        c['full_read_strand']=strand
                        reference_left=read_left==(strand=='+')
                        run=c['selected_expected_run']
                        run.update(start0=100 if reference_left else 106,end0=106 if reference_left else 112,stop_position0=99 if reference_left else 112)
                        if mode=='exact':
                            c.update(mode='exact',flank_hsp_id=None,full_read_hsp_id=None,
                                exact_flank_hit=dict(target='chr1',start0=106 if reference_left else 100,end0=112 if reference_left else 106,strand='-' if strand=='+' else '+'))
                        self.output(d)
                        run.update(start0=106 if reference_left else 100,end0=112 if reference_left else 106,stop_position0=112 if reference_left else 99)
                        with self.assertRaisesRegex(module.EvidenceError,'wrong side'):self.write(d)
    def test_exact_flank_cannot_occupy_repeat_side(self):
        d=copy.deepcopy(self.doc);c=d['samples'][0]['contexts'][0]
        c.update(mode='exact',flank_hsp_id=None,full_read_hsp_id=None,exact_flank_hit=dict(target='chr1',start0=100,end0=106,strand='+'))
        with self.assertRaisesRegex(module.EvidenceError,'nonrepeat side'):self.write(d)
    def test_dependency_bounds_precede_runtime_fingerprinting(self):
        self.write()
        path=Path(self.doc['source_artifacts'][0]['path'])
        with path.open('wb') as stream:stream.truncate(module.MAX_ARTIFACT_BYTES+1)
        with self.assertRaisesRegex(module.EvidenceError,'small regular'):module.source_dependencies(self.path)
        path.unlink();path.mkdir()
        with self.assertRaisesRegex(module.EvidenceError,'small regular'):module.source_dependencies(self.path)
    def test_dependency_total_budget_guard(self):
        self.write()
        with mock.patch.object(module,'MAX_SOURCE_BYTES',1):
            with self.assertRaisesRegex(module.EvidenceError,'total source'):module.source_dependencies(self.path)
    def test_table_contracts(self):
        out=self.output()
        for name,fields in [('family_summaries',module.FAMILY_FIELDS),('intervals',module.INTERVAL_FIELDS),('events',module.EVENT_FIELDS),('contexts',module.CONTEXT_FIELDS),('positions',module.POSITION_FIELDS)]:
            self.assertTrue(out[name])
            for row in out[name]:self.assertTrue(set(fields)<=set(row),(name,set(fields)-set(row)))


if __name__=='__main__':unittest.main()
