"""Saved followup evidence cannot replace native decisions or wrong-source data."""
import copy
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from followup_diagnostics import load_followup_diagnostics


class FollowupDiagnosticsTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.directory=self.root/'followup';self.folder=self.directory/'sample'
        self.folder.mkdir(parents=True)
        self.meta={'scientific_name':'Fixture species','tax_id':'123'}
        self.source_hash='a'*64;self.raw={'version':'0.8.1',
            'inputs':[{'path':str(self.root/'never_read_raw.fastq.gz'),'bytes':999}],
            'summary':{'records':100,'bases':10000,'sampled_prefix':False},
            'candidates':[{'motif':'TTAGGG','reads':50}], 'config':{'mode':'reads'}}
        self.identity=dict(native_result_sha256=self.source_hash,native_version='0.8.1',inputs=self.raw['inputs'],
            records=100,bases=10000,**self.meta)
        self.historical={'path':self.raw['inputs'][0]['path'],'bytes':999,'sha256':'b'*64}
        self.compact={'path':str(self.root/'never_read_compact.fastq'),'size_bytes':500,'sha256':'c'*64,
            'records':2,'bases':200,'validated_eof':True}
        self.summary={'status':'complete','source_id':'source','source_records':2,'source_bases':200,
            'selected_records':2,'pairs':2,'unique_tracts':3,'unique_intervals':12,'empty_cores':1,
            'frozen_stage_counts':{'exact_boundary_blocked':1,'compatible_positive_gap':1},
            'new_model_evaluations':False,'changes_native_classification':False,'input_binding':self.compact,
            'runtime':{'wall_seconds':1.0}}
        self.table=[dict(pair_id='pair1',source_id='source',diagnostic_stage='exact_boundary_blocked',
                        native_status_changed='False',model_recomputed='False'),
                    dict(pair_id='pair2',source_id='source',diagnostic_stage='compatible_positive_gap',
                        native_status_changed='False',model_recomputed='False')]
        self.receipt={'source_identity':{'path':self.raw['inputs'][0]['path'],'bytes':999},
            'input_records':100,'input_bases':10000,'complete_input_eof':True,'records':2,'bases':200,
            'output_binding':self.compact,'source_hash_scope':'Archived checksum; no raw reread'}
        self.kind='fixed_pair_bridge';self.relation_kind='registered_full_record_extraction'
        self.status={'schema_version':1,'sample_id':'sample','status':'complete','exit_code':0,'identity':self.identity}

    def dump(self,name,value):
        p=self.root/name;p.write_text(json.dumps(value));return self.bind(p)

    @staticmethod
    def bind(path):
        return dict(path=str(path.resolve()),size_bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())

    def save(self):
        tablepath=self.root/'table.tsv'
        fields=list(self.table[0]) if self.table else ['pair_id','source_id','diagnostic_stage','native_status_changed','model_recomputed']
        with tablepath.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(self.table)
        table=self.bind(tablepath);self.summary['outputs']=[table]
        summary=self.dump('summary.json',self.summary);receipt=self.dump('receipt.json',self.receipt)
        review=self.dump('review.json',{'status':'passed','bindings':[summary,table]})
        relation_review=self.dump('relation_review.json',{'status':'passed','bindings':[receipt,self.historical]})
        self.evidence={'schema_version':1,'kind':self.kind,'diagnostic_source_id':'source',
            'artifacts':{'summary':summary,'table':table,'review':review},
            'source_relation':{'kind':self.relation_kind,'receipt':receipt,'review':relation_review,
                               'original_input_binding':self.historical},
            'data_role':'posthoc_development_diagnostics_previously_seen'}
        self.save_evidence()

    def save_evidence(self):
        p=self.folder/'evidence.json';p.write_text(json.dumps(self.evidence))
        self.status['evidence_sha256']=self.bind(p)['sha256']
        (self.folder/'status.json').write_text(json.dumps(self.status))

    def load(self,**kwargs):
        return load_followup_diagnostics('sample',self.meta,self.raw,self.source_hash,self.directory,**kwargs)

    def test_missing_not_requested_and_failed_are_not_zero_measurements(self):
        x=load_followup_diagnostics('sample',self.meta,self.raw,self.source_hash,None)
        self.assertEqual(x['status'],'not_requested');self.assertIsNone(x['measurements']);self.assertEqual(x['rows'],[])
        x=load_followup_diagnostics('absent',self.meta,self.raw,self.source_hash,self.directory)
        self.assertEqual(x['status'],'not_available');self.assertIsNone(x['measurements'])
        failure=self.dump('failure.json',{'status':'operational_failure','source_id':'source'})
        self.status.update(status='failed',exit_code=1,failure_reason='fixture input error',failure_artifact=failure)
        (self.folder/'status.json').write_text(json.dumps(self.status))
        x=self.load();self.assertEqual(x['status'],'failed');self.assertIsNone(x['measurements']);self.assertEqual(x['rows'],[])
        self.assertFalse(x['changes_native_decision'])

    def test_complete_subset_join_preserves_native_and_never_reads_sequences(self):
        self.save();before=copy.deepcopy(self.raw);x=self.load()
        self.assertEqual(x['status'],'complete');self.assertEqual(x['measurements']['registered_pairs'],2)
        self.assertEqual(self.raw,before);self.assertFalse(x['whole_sample_inference_permitted'])
        scopes={r['metric']:(r['value'],r['scope']) for r in x['rows']}
        self.assertEqual(scopes['original_records'],(100,'original_whole_input'))
        self.assertEqual(scopes['compact_records'],(2,'selected_full_records'))
        self.assertFalse(Path(self.raw['inputs'][0]['path']).exists())
        self.assertFalse(Path(self.compact['path']).exists())

    def test_empty_complete_selection_is_explicit_and_still_not_absence(self):
        self.table=[];self.summary.update(pairs=0,unique_tracts=0,unique_intervals=0,empty_cores=0,
                                         selected_records=0,frozen_stage_counts={})
        self.save();x=self.load();self.assertEqual(x['status'],'complete')
        self.assertEqual(x['measurements']['registered_pairs'],0)
        self.assertFalse(x['whole_sample_inference_permitted'])

    def test_wrong_sample_taxon_hash_version_inputs_and_totals_rejected(self):
        original=copy.deepcopy(self.status)
        mutations=[lambda:self.status.update(sample_id='other'),
            lambda:self.status['identity'].update(tax_id='999'),
            lambda:self.status['identity'].update(scientific_name='Another species'),
            lambda:self.status['identity'].update(native_result_sha256='d'*64),
            lambda:self.status['identity'].update(native_version='0.8.0'),
            lambda:self.status['identity'].update(inputs=[{'path':'other','bytes':999}]),
            lambda:self.status['identity'].update(records=99)]
        for mutate in mutations:
            self.status=copy.deepcopy(original);mutate();self.save()
            with self.assertRaises(ValueError):self.load()

    def test_partial_corrupt_and_changed_artifacts_fail_closed(self):
        with self.assertRaises(ValueError):self.load()
        self.save();(self.root/'table.tsv').write_text((self.root/'table.tsv').read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'hash or size mismatch'):self.load()
        self.save();self.status['evidence_sha256']='0'*64
        (self.folder/'status.json').write_text(json.dumps(self.status))
        with self.assertRaisesRegex(ValueError,'evidence hash'):self.load()
        with self.assertRaises(ValueError):load_followup_diagnostics('sample',self.meta,self.raw,self.source_hash,self.root/'absent')

    def test_source_relation_and_subset_denominators_are_checked(self):
        original=copy.deepcopy((self.receipt,self.summary,self.historical))
        mutations=[lambda:self.receipt['source_identity'].update(path=str(self.root/'other.fastq')),
            lambda:self.receipt.update(input_records=99),lambda:self.receipt.update(records=3),
            lambda:self.summary.update(selected_records=3),lambda:self.historical.update(bytes=1000),
            lambda:self.receipt.update(complete_input_eof=False)]
        for mutate in mutations:
            self.receipt,self.summary,self.historical=copy.deepcopy(original);mutate();self.save()
            with self.assertRaises(ValueError):self.load()

    def test_unreviewed_artifact_and_unbound_table_are_rejected(self):
        self.save();empty=self.dump('review.json',{'status':'passed','bindings':[]})
        self.evidence['artifacts']['review']=empty;self.save_evidence()
        with self.assertRaisesRegex(ValueError,'does not bind'):self.load()
        self.save();self.summary['outputs']=[]
        self.evidence['artifacts']['summary']=self.dump('summary.json',self.summary)
        self.evidence['artifacts']['review']=self.dump('review.json',{'status':'passed','bindings':[
            self.evidence['artifacts']['summary'],self.evidence['artifacts']['table']]});self.save_evidence()
        with self.assertRaisesRegex(ValueError,'summary does not bind'):self.load()

    def test_no_promotion_of_incomplete_native_analysis(self):
        self.save();x=self.load(analysis_complete=False)
        self.assertFalse(x['changes_native_decision']);self.assertFalse(x['provenance']['native_analysis_complete'])
        self.assertTrue(any('Native 분석은 완료되지' in r for r in x['recommendations']))

    def test_duplicate_rows_bad_counts_and_changed_decisions_rejected(self):
        self.table[1]['pair_id']='pair1';self.save()
        with self.assertRaisesRegex(ValueError,'unique identity'):self.load()
        self.table[1]['pair_id']='pair2';self.summary['frozen_stage_counts']['exact_boundary_blocked']=True;self.save()
        with self.assertRaisesRegex(ValueError,'nonnegative count'):self.load()
        self.summary['frozen_stage_counts']['exact_boundary_blocked']=1
        self.table[0]['native_status_changed']='True';self.save()
        with self.assertRaisesRegex(ValueError,'changed native decision'):self.load()

    def setup_mixed(self):
        self.kind='registered_mixed_windows';self.relation_kind='sequence_unique_full_record_selection'
        self.summary={'status':'complete','source':self.compact,'source_records':2,'source_bases':200,
            'selected_records':1,'historical_pairs':1,'analyses':4,
            'all_orientations_and_vocabularies_retained':True,'whole_source_biologically_analyzed':False,
            'changes_native_classification':False,'vocabularies':{'single9':['TTGGGTTAC'],'known9_11':['TTGGGTTAC','ACTTGGGTTAC']}}
        self.receipt={'source_fastq':self.raw['inputs'][0]['path'],'source_fastq_sha256':self.historical['sha256'],
            'total_full_records':100,'total_full_bases':10000,'gzip_eof_verified':True,
            'extracted_fastq':self.compact['path'],'extracted_fastq_sha256':self.compact['sha256'],
            'distinct_sequence_unique_extracted_archive_reads':2,'original_name_identity_confirmed':False,
            'caveat':'Original names were stripped; sequence association only.'}
        self.table=[]
        for name,units in self.summary['vocabularies'].items():
            for strand in ('+','-'):
                yes='True' if name=='known9_11' and strand=='+' else 'False'
                self.table.append(dict(pair_id='pair1',record_id='record1',vocabulary=name,requested_strand=strand,
                    primary_units=json.dumps(units),changes_native_classification='False',physical_telomere_confirmed='False',
                    exact_covers_original_phase_anchor_span=yes,edited_covers_original_phase_anchor_span=yes,
                    exact_full_window_path=yes,edited_full_window_path=yes))

    def test_mixed_registered_vocabularies_both_orientations_stay_separate(self):
        self.setup_mixed();self.save();x=self.load()
        self.assertEqual(x['measurements']['registered_pairs'],1);self.assertEqual(x['measurements']['analyses'],4)
        self.assertEqual(len(x['measurements']['groups']),4)
        self.assertEqual(sum(r['exact_anchor_span_paths'] for r in x['measurements']['groups']),1)
        self.assertFalse(x['provenance']['source_relation']['original_molecule_identity_confirmed'])
        self.assertTrue(any('qualifying mixed array' in r for r in x['recommendations']))

    def test_mixed_source_checksum_or_vocabulary_mismatch_rejected(self):
        self.setup_mixed();self.receipt['source_fastq_sha256']='e'*64;self.save()
        with self.assertRaisesRegex(ValueError,'historical source checksum'):self.load()
        self.receipt['source_fastq_sha256']=self.historical['sha256'];self.table[0]['primary_units']='["TTAGGG"]';self.save()
        with self.assertRaisesRegex(ValueError,'vocabulary or inference'):self.load()


if __name__=='__main__':unittest.main()
