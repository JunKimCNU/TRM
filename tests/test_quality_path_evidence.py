"""Sealed-source and scientific-boundary regression tests; no sequence input."""
import copy
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import quality_path_evidence as module


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def json_write(path,document):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(document,indent=2,ensure_ascii=False)+'\n')


def tsv_write(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader()
        writer.writerows({key:json.dumps(value,separators=(',',':')) if isinstance(value,(dict,list)) else value for key,value in row.items()} for row in rows)


def seal(document):
    document['payload_sha256']=module.payload_sha256(document)
    return document


def write_fixture(directory,sample_id='SRRknown',native_sha='a'*64,identity=None,
                  scientific_name='Example fungus',tax_id='123',biosample='SAMNexample',
                  representative_motif='TTAGGC',review_notes=None):
    directory=Path(directory)
    identity=identity or dict(scientific_name=scientific_name,tax_id=tax_id,biosample=biosample)
    document=dict(sample_id=sample_id,scientific_name=identity['scientific_name'],tax_id=identity['tax_id'],
        metadata=dict(sample_accession=identity['biosample']),provenance=dict(original_discovery_sha256=native_sha),
        native_scan_version='0.8.5',analysis_complete=True,detected_motif=representative_motif,
        classification='known_trm_detected',data_adequacy=dict(status='unassessable'),candidates=[dict(motif=representative_motif,reads=5)])
    sources=dict(prior_assessment=document,read_comparison=[],read_coverage=[],model_summary=[],interval_pair_summary=[])
    for spot in (1,2):
        for mate in (1,2):
            modeled=spot==1
            pair=f'{sample_id}__pair{spot}'
            read=f'{pair}__R{mate}' if modeled else f'{sample_id}__ordinary_prefix_spot{spot}__R{mate}'
            sources['read_comparison'].append(dict(sample_id=sample_id,query_id=read,spot_id=spot,mate=mate,bases=10,
                groups='focused' if modeled else 'ordinary',sequence_and_instrument_name_match=True,orientation='direct',
                saved_quality='constant_Q30',normalized_quality_bins='11,25,37',normalized_variable=True,
                quality_differing_positions=10,independent_biological_replication=False,original_submitted_quality_confirmed=False))
            sources['read_coverage'].append(dict(sample_id=sample_id,read_query_id=read,
                roles=['focused_oriented_mixed_boundary_pair'] if modeled else ['ordinary_unselected_prefix_control_not_biological_negative'],
                saved_intervals=1 if modeled else 0,model_quality_evaluated=modeled,
                reason='saved_intervals_annotated' if modeled else 'no_preselected_model_interval_not_a_negative_result'))
            if not modeled:continue
            iid=f'{pair}__mate{mate}__0_10'
            base=dict(sample_id=sample_id,read_query_id=read,pair_query_id=pair,spot_id=spot,mate=mate,
                interval_id=iid,start0=0,end0=10,focused=True,mixed=mate==1,event_ids=[])
            for model,motif in module.MODELS.items():
                observed={'period5':2,'period6':0,'period11':1}[model]
                deletion=2 if model=='period6' else 0
                q30=1 if model=='period5' else 0
                sources['model_summary'].append(dict(base,model=model,motif=motif,minimum_edit_cost=observed+deletion,
                    operation_count_ranges={key:dict(min=value,max=value) for key,value in dict(match=10-observed,substitution=0,insertion=observed,deletion=deletion).items()},
                    path_quality_count_ranges={key:dict(min=value,max=value) for key,value in dict(observed_edits=observed,observed_edits_Q_lt20=observed-q30,observed_edits_Q20_29=0,observed_edits_Q_ge30=q30).items()},
                    states=dict(M=10-observed,A=observed) if observed else dict(M=10),mandatory_deletion_boundaries=0,
                    possible_deletion_boundaries=3 if deletion else 0,deletion_required_but_no_mandatory_boundary=bool(deletion),
                    all_observed_positions_always_match=not observed,exclusively_gap_edits=bool(deletion),quality_used_to_optimize=False,
                    count_extrema_jointly_attainable_not_asserted=True))
            for left,right,hist in [('period5','period6',dict(AM=2,MM=8)),('period5','period11',dict(AA=1,AM=1,MM=8)),('period6','period11',dict(MA=1,MM=9))]:
                sources['interval_pair_summary'].append(dict(base,left_model=left,right_model=right,
                    left_cost=2 if left!='period11' else 1,right_cost=2 if right!='period11' else 1,bases=10,
                    state_pair_counts=hist,left_only_gap_edits=left=='period6',right_only_gap_edits=right=='period6'))
    paths=dict(prior_assessment=directory/'assessment/prior.json',read_comparison=directory/'quality/read.tsv',
        read_coverage=directory/'paths/coverage.tsv',model_summary=directory/'paths/models.tsv',interval_pair_summary=directory/'paths/pairs.tsv')
    for role,path in paths.items():
        if role in module.JSON_ROLES:json_write(path,sources[role])
        else:tsv_write(path,sources[role])
    groups=dict(assessment_completion=('assessment',('prior_assessment',)),comparison_completion=('quality',('read_comparison',)),
        path_completion=('paths',('read_coverage','model_summary','interval_pair_summary')))
    for role,(folder,children) in groups.items():
        path=directory/folder/'completion.json'
        json_write(path,dict(status='complete',artifacts_sha256={str(paths[child].relative_to(path.parent)):digest(paths[child]) for child in children}))
        paths[role]=path
    parsed={role:module._json(path.read_bytes(),path) if role in module.JSON_ROLES else module._table(path.read_bytes(),path) for role,path in paths.items()}
    payload=seal(dict(schema_version=1,scope=module.SCOPE,sample_id=sample_id,**module._derive(sample_id,parsed),review_notes=list(review_notes or []),
        source_artifacts=[dict(role=role,path=str(paths[role]),sha256=digest(paths[role])) for role in module.SOURCE_ROLES]))
    payload_path=directory/'payload.json';json_write(payload_path,payload)
    manifest=seal(dict(schema_version=1,scope=module.SCOPE,samples=[dict(sample_id=sample_id,identity=identity,original_native_sha256=native_sha,
        evidence=dict(path=payload_path.name,sha256=digest(payload_path)))]))
    manifest_path=directory/'manifest.json';json_write(manifest_path,manifest)
    return dict(manifest_path=manifest_path,payload_path=payload_path,source_path=paths['model_summary'],source_paths=paths,
        document=copy.deepcopy(document),payload=payload,manifest=manifest)


class QualityPathEvidenceTests(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory();self.addCleanup(temp.cleanup);self.directory=Path(temp.name)
        self.fx=write_fixture(self.directory)

    def load(self):return module.load_evidence(self.fx['manifest_path'])
    def build(self,loaded=None,document=None):
        return module.build_evidence('SRRknown',document or self.fx['document'],self.load() if loaded is None else loaded)
    def save_payload(self,payload):
        json_write(self.fx['payload_path'],seal(payload))
        manifest=json.loads(self.fx['manifest_path'].read_text());manifest['samples'][0]['evidence']['sha256']=digest(self.fx['payload_path'])
        json_write(self.fx['manifest_path'],seal(manifest))
    def update_source_pin(self,role,update_completion=False):
        paths=self.fx['source_paths']
        if update_completion:
            receipt_role='assessment_completion' if role=='prior_assessment' else 'comparison_completion' if role=='read_comparison' else 'path_completion'
            path=paths[receipt_role];receipt=json.loads(path.read_text())
            receipt['artifacts_sha256'][str(paths[role].relative_to(path.parent))]=digest(paths[role]);json_write(path,receipt)
        payload=json.loads(self.fx['payload_path'].read_text())
        for pin in payload['source_artifacts']:pin['sha256']=digest(paths[pin['role']])
        self.save_payload(payload)
    def mutate_source(self,role,key,value,row_index=0):
        path=self.fx['source_paths'][role];rows=list(csv.DictReader(io.StringIO(path.read_text()),delimiter='\t'))
        rows[row_index][key]=value;tsv_write(path,rows);self.update_source_pin(role,True)

    def test_valid_overlay_preserves_native_and_distinguishes_path_position_gap(self):
        before=copy.deepcopy(self.fx['document']);result=self.build()
        self.assertEqual(self.fx['document'],before);self.assertEqual(result['status'],'sources_validated')
        self.assertFalse(result['changes_classification_or_adequacy'])
        self.assertEqual(set(result['summary']),set(module.SUMMARY_FIELDS));self.assertEqual(set(result['display']),set(module.FIELDS))
        s=result['summary'];self.assertEqual((s['normalized_same_read_comparisons'],s['modeled_reads'],s['ordinary_reads_without_models']),(4,2,2))
        self.assertEqual((s['focused_reads'],s['focused_mixed_reads'],s['focused_mate_only_reads'],s['focused_pairs']),(2,1,1,1))
        self.assertEqual(s['focused_mixed_union_bases'],10);self.assertEqual(s['strict_fixed_position_contrast_intervals'],0)
        m={r['model']:r for r in result['models']}
        self.assertEqual(m['period5']['intervals_requiring_Q_ge30_observed_edit'],1)
        self.assertEqual(m['period6']['intervals_permitting_zero_observed_edit'],1)
        self.assertEqual(m['period6']['intervals_with_exclusively_gap_edits'],1)
        self.assertIsNone(m['period6']['gap_base_quality']);self.assertFalse(m['period6']['quality_used_to_optimize'])
        self.assertEqual(result['representative_motif'],'TTAGGC')
        self.assertTrue(all(s[k]=='not_confirmed' for k in ('physical_terminal_status','main_trm_status','novel_trm_status','independent_validation_status')))

    def test_absent_measurements_null_and_no_three_sample_claims(self):
        self.assertIsNone(module.load_evidence(None));self.assertEqual(module.source_dependencies(None),[])
        for sid,loaded in [('SRRknown',None),('other',self.load())]:
            result=module.build_evidence(sid,self.fx['document'],loaded)
            self.assertEqual(result['status'],'not_provided');self.assertTrue(all(v is None for v in result['summary'].values()))
            for key in ('models','read_comparisons','source_artifacts'):self.assertEqual(result[key],[])
            self.assertEqual(result['limitations'],module.ABSENT_LIMITATIONS);self.assertNotIn('Q30',' '.join(result['limitations']))

    def test_dependencies_and_staleness_even_restored_mtime(self):
        loaded=self.load();dependencies=module.source_dependencies(self.fx['manifest_path'])
        self.assertEqual(len(dependencies),10);self.assertIn(str(self.fx['source_path']),dependencies)
        path=self.fx['source_path'];st=path.stat();path.write_text(path.read_text()+'\n');os.utime(path,ns=(st.st_atime_ns,st.st_mtime_ns))
        with self.assertRaisesRegex(module.EvidenceError,'stale'):self.build(loaded)

    def test_changed_hash_missing_file_and_receipt_binding_rejected(self):
        path=self.fx['source_path'];path.write_text(path.read_text()+'\n')
        with self.assertRaisesRegex(module.EvidenceError,'hash differs'):self.load()
        self.update_source_pin('model_summary')
        with self.assertRaisesRegex(module.EvidenceError,'completed artifact seal'):self.load()
        path.unlink()
        with self.assertRaises(module.EvidenceError):self.load()

    def test_native_identity_version_and_representative_mismatch_rejected(self):
        loaded=self.load()
        changes=[lambda d:d.update(sample_id='other'),lambda d:d.update(scientific_name='other'),lambda d:d.update(tax_id='99'),
            lambda d:d['metadata'].update(sample_accession='other'),lambda d:d['provenance'].update(original_discovery_sha256='b'*64),
            lambda d:d.update(native_scan_version='0.8.6'),lambda d:d.update(analysis_complete=False),lambda d:d.update(detected_motif='TTAGGG')]
        for change in changes:
            d=copy.deepcopy(self.fx['document']);change(d)
            with self.subTest(document=d),self.assertRaises(module.EvidenceError):self.build(loaded,d)

    def test_exported_summary_or_model_cannot_override_detailed_sources(self):
        for field,key in [('summary','modeled_reads'),('models','intervals_requiring_Q_ge30_observed_edit')]:
            p=copy.deepcopy(self.fx['payload']);target=p[field][0] if field=='models' else p[field];target[key]+=1;self.save_payload(p)
            with self.subTest(field=field),self.assertRaisesRegex(module.EvidenceError,'differs from detailed'):self.load()

    def test_boolean_counts_schema_and_gap_quality_rejected(self):
        p=copy.deepcopy(self.fx['payload']);p['summary']['modeled_reads']=True;self.save_payload(p)
        with self.assertRaises(module.EvidenceError):self.load()
        for value in (0,30,False):
            p=copy.deepcopy(self.fx['payload']);p['models'][1]['gap_base_quality']=value;self.save_payload(p)
            with self.subTest(value=value),self.assertRaisesRegex(module.EvidenceError,'must be null'):self.load()
        manifest=copy.deepcopy(self.fx['manifest']);manifest['schema_version']=True;json_write(self.fx['manifest_path'],seal(manifest))
        with self.assertRaises(module.EvidenceError):self.load()

    def test_positive_confirmations_or_quality_reoptimization_rejected(self):
        for field in ('independent_biological_replication','original_submitted_quality_confirmed'):
            p=copy.deepcopy(self.fx['payload']);p['read_comparisons'][0][field]=True;self.save_payload(p)
            with self.subTest(field=field),self.assertRaises(module.EvidenceError):self.load()
        p=copy.deepcopy(self.fx['payload']);p['models'][0]['quality_used_to_optimize']=True;self.save_payload(p)
        with self.assertRaises(module.EvidenceError):self.load()

    def test_duplicate_payload_rows_samples_and_json_keys_rejected(self):
        for field in ('models','read_comparisons'):
            p=copy.deepcopy(self.fx['payload']);p[field].append(copy.deepcopy(p[field][0]));self.save_payload(p)
            with self.subTest(field=field),self.assertRaises(module.EvidenceError):self.load()
        self.save_payload(copy.deepcopy(self.fx['payload']))
        manifest=json.loads(self.fx['manifest_path'].read_text());manifest['samples'].append(copy.deepcopy(manifest['samples'][0]));json_write(self.fx['manifest_path'],seal(manifest))
        with self.assertRaisesRegex(module.EvidenceError,'duplicate sample'):self.load()
        for text in ('{"schema_version":1,"schema_version":1}','{"schema_version":NaN}'):
            self.fx['manifest_path'].write_text(text)
            with self.assertRaises(module.EvidenceError):self.load()

    def test_duplicate_source_identity_rejected_even_resealed(self):
        path=self.fx['source_paths']['model_summary'];rows=list(csv.DictReader(io.StringIO(path.read_text()),delimiter='\t'));rows.append(copy.deepcopy(rows[0]));tsv_write(path,rows)
        self.update_source_pin('model_summary',True)
        with self.assertRaisesRegex(module.EvidenceError,'duplicate source row'):self.load()

    def test_loaded_object_change_rejected_even_resealed(self):
        loaded=self.load();doc=loaded['samples']['SRRknown']['document'];doc['summary']['modeled_reads']=9
        with self.assertRaisesRegex(module.EvidenceError,'payload seal'):self.build(loaded)
        seal(doc)
        with self.assertRaisesRegex(module.EvidenceError,'changed since source validation'):self.build(loaded)

    def test_model_geometry_scope_and_optimization_checked(self):
        original=self.fx['source_paths']['model_summary'].read_text()
        for key,value in [('start0','11'),('pair_query_id','wrong'),('focused','False'),('quality_used_to_optimize','True'),('count_extrema_jointly_attainable_not_asserted','False')]:
            self.fx['source_paths']['model_summary'].write_text(original);self.mutate_source('model_summary',key,value)
            with self.subTest(key=key),self.assertRaises(module.EvidenceError):self.load()

    def test_fixed_cost_and_paired_state_arithmetic_checked(self):
        row=list(csv.DictReader(io.StringIO(self.fx['source_paths']['model_summary'].read_text()),delimiter='\t'))[0]
        q=json.loads(row['path_quality_count_ranges']);q['observed_edits']['min']=1
        self.mutate_source('model_summary','path_quality_count_ranges',q)
        with self.assertRaisesRegex(module.EvidenceError,'fixed-cost'):self.load()
        # Restore a fresh fixture, then break only the table-to-table state margin.
        fx=write_fixture(self.directory/'fresh');self.fx=fx
        self.mutate_source('interval_pair_summary','state_pair_counts',dict(MM=10))
        with self.assertRaisesRegex(module.EvidenceError,'marginals'):self.load()

    def test_same_read_identity_and_quality_semantics_checked(self):
        original=self.fx['source_paths']['read_comparison'].read_text()
        for key,value in [('sequence_and_instrument_name_match','False'),('orientation','reverse'),('normalized_variable','False'),('quality_differing_positions','11'),('normalized_quality_bins','11,11'),('original_submitted_quality_confirmed','True')]:
            self.fx['source_paths']['read_comparison'].write_text(original);self.mutate_source('read_comparison',key,value)
            with self.subTest(key=key),self.assertRaises(module.EvidenceError):self.load()

    def test_missing_columns_nested_duplicate_keys_and_nonfinite_refused(self):
        original=self.fx['source_paths']['model_summary'].read_text()
        for value in ('{"M":10,"M":10}','{"M":NaN}'):
            self.fx['source_paths']['model_summary'].write_text(original);self.mutate_source('model_summary','states',value)
            with self.subTest(value=value),self.assertRaises(module.EvidenceError):self.load()
        path=self.fx['source_paths']['model_summary'];path.write_text(original)
        rows=list(csv.DictReader(io.StringIO(path.read_text()),delimiter='\t'))
        for row in rows:row.pop('focused')
        tsv_write(path,rows);self.update_source_pin('model_summary',True)
        with self.assertRaisesRegex(module.EvidenceError,'source shape'):self.load()

    def test_manifest_payload_identity_source_roles_and_input_bounds(self):
        manifest=copy.deepcopy(self.fx['manifest']);manifest['samples'][0]['identity']['biosample']='different';json_write(self.fx['manifest_path'],seal(manifest))
        with self.assertRaisesRegex(module.EvidenceError,'manifest/payload identity'):self.load()
        json_write(self.fx['manifest_path'],self.fx['manifest'])
        with patch.object(module,'MAX_BYTES',1),self.assertRaises(module.EvidenceError):self.load()
        with patch.object(module,'MAX_TOTAL_BYTES',1),self.assertRaises(module.EvidenceError):self.load()
        p=copy.deepcopy(self.fx['payload']);p['source_artifacts'][1]=copy.deepcopy(p['source_artifacts'][0]);self.save_payload(p)
        with self.assertRaises(module.EvidenceError):self.load()
        p=copy.deepcopy(self.fx['payload']);p['source_artifacts'][1]['path']=str(self.directory/'not_read.fastq');self.save_payload(p)
        with self.assertRaisesRegex(module.EvidenceError,'format not supported'):self.load()

    def test_rows_are_copied_and_notes_preserved_for_html_escaping(self):
        p=copy.deepcopy(self.fx['payload']);p['review_notes']=['<script>escaped by renderer</script>'];self.save_payload(p)
        loaded=self.load();before=copy.deepcopy(loaded['samples']['SRRknown']['document']);result=self.build(loaded)
        result['models'][0]['gap_base_quality']=30;result['read_comparisons'][0]['normalized_quality_bins'].append(40)
        self.assertEqual(loaded['samples']['SRRknown']['document'],before)
        self.assertIn('<script>',result['display']['quality_path_limit_display'])


if __name__=='__main__':unittest.main()
