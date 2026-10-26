"""Saved same-read quality/path annotations are optional, bounded, and cache-safe."""
from copy import deepcopy
import html
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'tests')]
import quality_path_evidence as module
import report_assessment
import report_runtime
from summarize_interpretation import summarize_sample
import test_report_assessment as fixture
import test_quality_path_evidence as evidence_fixture

TABLES=('quality_path_summary.tsv','quality_path_models.tsv','quality_path_read_comparisons.tsv')


class QualityPathReportingTests(unittest.TestCase):
    command=fixture.AssessmentReportTests.command
    document=fixture.AssessmentReportTests.document

    def setUp(self):
        fixture.AssessmentReportTests.setUp(self)
        fixture.AssessmentReportTests.make_native_sources(self)
        for path in self.source.glob('*/result.json'):
            raw=json.loads(path.read_text());raw.update(program='trm',version='0.8.5')
            fixture.write_json(path,raw)

    def run_report(self,success=True,extra=()):
        return fixture.AssessmentReportTests.run_report(self,success=success,native=True,extra=extra)

    def save_evidence(self,review_notes=None):
        self.saved=evidence_fixture.write_fixture(self.root/'quality-fixture',
            sample_id='SRRknown',identity=dict(scientific_name='Fixture species',tax_id='6239',biosample='SAMNfixture_SRRknown'),
            native_sha=fixture.digest(self.source/'SRRknown/result.json'),representative_motif='TTAGGC',review_notes=review_notes)
        self.manifest=Path(self.saved['manifest_path']);self.payload_path=Path(self.saved['payload_path'])
        self.pin=Path(self.saved['source_path'])
        return ['--quality-path-evidence',str(self.manifest)]

    def modes(self):
        return {r['sample_id']:r['mode'] for r in fixture.read_tsv(self.output/'interpretation_resources.tsv')}

    def reseal_notes(self,note):
        payload=json.loads(self.payload_path.read_text());payload['review_notes'].append(note)
        payload['payload_sha256']=module.payload_sha256(payload);fixture.write_json(self.payload_path,payload)
        manifest=json.loads(self.manifest.read_text())
        manifest['samples'][0]['evidence']['sha256']=fixture.digest(self.payload_path)
        manifest['payload_sha256']=module.payload_sha256(manifest);fixture.write_json(self.manifest,manifest)

    def test_overlay_preserves_candidate_quality_and_policy_namespaces(self):
        self.run_report();before=self.document('SRRknown')
        tables={name:(self.output/name).read_bytes() for name in ('candidates.tsv','candidate_quality_observations.tsv','adequacy.tsv')}
        self.run_report(extra=self.save_evidence());after=self.document('SRRknown')
        self.assertEqual(set(before),set(after))
        for key in before:
            if key!='quality_path_evidence':self.assertEqual(before[key],after[key],key)
        for name,data in tables.items():self.assertEqual(data,(self.output/name).read_bytes(),name)
        self.assertEqual(summarize_sample(before,Path('fixed.json'),'0'*64),summarize_sample(after,Path('fixed.json'),'0'*64))
        evidence=after['quality_path_evidence'];summary=evidence['summary']
        self.assertEqual(evidence['status'],'sources_validated')
        self.assertEqual(after['detected_motif'],'TTAGGC')
        self.assertFalse(evidence['changes_classification_or_adequacy'])
        self.assertEqual(summary['physical_terminal_status'],'not_confirmed')
        self.assertEqual(summary['main_trm_status'],'not_confirmed')
        self.assertEqual(summary['novel_trm_status'],'not_confirmed')
        for model in evidence['models']:
            self.assertIsNone(model['gap_base_quality'])
            self.assertFalse(model['quality_used_to_optimize'])
        for read in evidence['read_comparisons']:
            self.assertFalse(read['independent_biological_replication'])
            self.assertFalse(read['original_submitted_quality_confirmed'])
        rendered=(self.output/'samples/SRRknown.html').read_text()
        for table in TABLES:self.assertIn('../'+table,rendered)
        for phrase in ('gap Q는 null','고정 위치 대비','품질로 최적 경로나 비용을 다시 정하지 않았고'):
            self.assertIn(phrase,rendered)
        row=next(r for r in fixture.read_tsv(self.output/TABLES[0]) if r['sample_id']=='SRRknown')
        for key in module.SUMMARY_FIELDS:self.assertEqual(row[key],str(summary[key]))
        self.assertEqual(len(fixture.read_tsv(self.output/TABLES[1])),len(evidence['models']))
        self.assertEqual(len(fixture.read_tsv(self.output/TABLES[2])),len(evidence['read_comparisons']))

    def test_absent_evidence_stays_null_in_json_and_blank_in_tsv(self):
        self.run_report()
        for sid in ('SRRknown','SRRnovel','SRRnone','SRRfailed'):
            evidence=self.document(sid)['quality_path_evidence']
            self.assertEqual(evidence['status'],'not_provided')
            self.assertTrue(all(evidence['summary'][key] is None for key in module.SUMMARY_FIELDS))
            self.assertEqual(evidence['models'],[]);self.assertEqual(evidence['read_comparisons'],[])
        self.run_report(extra=self.save_evidence())
        for sid in ('SRRnovel','SRRnone','SRRfailed'):
            self.assertEqual(self.document(sid)['quality_path_evidence']['status'],'not_provided')
        absent=next(r for r in fixture.read_tsv(self.output/TABLES[0]) if r['sample_id']=='SRRnone')
        self.assertTrue(all(absent[key]=='' for key in module.SUMMARY_FIELDS))

    def test_serial_parallel_and_resume_have_identical_artifacts(self):
        extra=self.save_evidence();self.run_report(extra=extra+['--jobs','1'])
        docs={p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')}
        tables={name:(self.output/name).read_bytes() for name in TABLES}
        self.run_report(extra=extra+['--jobs','2'])
        self.assertEqual(docs,{p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')})
        self.assertEqual(tables,{name:(self.output/name).read_bytes() for name in TABLES})
        self.run_report(extra=extra+['--jobs','2','--resume'])
        self.assertTrue(all(mode=='reused' for mode in self.modes().values()))
        self.assertEqual(docs,{p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')})

    def test_changed_source_is_rejected_before_touching_valid_warm_cache(self):
        extra=self.save_evidence()+['--jobs','2','--resume'];self.run_report(extra=extra)
        before={str(p.relative_to(self.output)):p.read_bytes() for p in self.output.rglob('*') if p.is_file()}
        self.pin.write_text(self.pin.read_text()+'\n')
        outcome=self.run_report(success=False,extra=extra)
        self.assertIn('source hash differs',outcome.stderr)
        self.assertEqual(before,{str(p.relative_to(self.output)):p.read_bytes() for p in self.output.rglob('*') if p.is_file()})
        self.assertEqual(json.loads((self.output/'interpretation_status.json').read_text())['status'],'complete')

    def test_resealed_annotation_invalidates_cache_and_removal_clears_tables(self):
        extra=self.save_evidence()+['--resume'];self.run_report(extra=extra)
        self.reseal_notes('Additional review note; observations and biological calls unchanged.')
        self.run_report(extra=extra)
        self.assertTrue(all(mode=='recomputed' for mode in self.modes().values()))
        self.assertIn('Additional review note',str(self.document('SRRknown')['quality_path_evidence']))
        self.run_report(extra=['--resume'])
        self.assertTrue(all(mode=='recomputed' for mode in self.modes().values()))
        self.assertEqual(self.document('SRRknown')['quality_path_evidence']['status'],'not_provided')
        for name in TABLES[1:]:self.assertEqual(fixture.read_tsv(self.output/name),[])

    def test_native_source_mismatch_and_missing_manifest_fail_closed(self):
        outcome=self.run_report(success=False,extra=['--quality-path-evidence',str(self.root/'missing.json')])
        self.assertIn('source unavailable',outcome.stderr)
        extra=self.save_evidence();source=self.source/'SRRknown/result.json';source.write_text(source.read_text()+'\n')
        outcome=self.run_report(success=False,extra=extra)
        self.assertIn('native SHA differs',outcome.stderr)
        self.assertEqual(json.loads((self.output/'interpretation_status.json').read_text())['status'],'failed')

    def test_cache_receipt_requires_all_new_table_chunks(self):
        self.run_report(extra=self.save_evidence())
        receipt=json.loads((self.output/'.report-cache/SRRknown/receipt.json').read_text())
        self.assertTrue(report_runtime.verify_artifacts(receipt,output=self.output,sample='SRRknown'))
        bundle_path=Path(receipt['bundle']);bundle=json.loads(bundle_path.read_text())
        for name in TABLES:
            with self.subTest(name=name):
                changed_bundle=deepcopy(bundle);changed_bundle['chunks'].remove(name);fixture.write_json(bundle_path,changed_bundle)
                changed_receipt=deepcopy(receipt)
                changed_receipt['artifacts'][str(bundle_path)]=report_runtime.binding(bundle_path)[0]
                changed_receipt['artifacts'].pop(str(bundle_path.parent/name))
                self.assertFalse(report_runtime.verify_artifacts(changed_receipt,output=self.output,sample='SRRknown'))
        fixture.write_json(bundle_path,bundle)

    def test_free_text_and_rendered_dynamic_fields_are_html_escaped(self):
        hostile='<img src=x onerror="alert(1)"> & review'
        self.run_report(extra=self.save_evidence(review_notes=[hostile]))
        evidence=self.document('SRRknown')['quality_path_evidence']
        self.assertIn(hostile,evidence['review_notes'])
        rendered=(self.output/'samples/SRRknown.html').read_text()
        self.assertNotIn(hostile,rendered);self.assertIn(html.escape(hostile,quote=True),rendered)
        altered=deepcopy(evidence);altered['models'][0]['motif']=hostile
        rendered=''.join(report_assessment.quality_path_html(altered))
        self.assertNotIn(hostile,rendered);self.assertIn(html.escape(hostile,quote=True),rendered)

if __name__=='__main__':unittest.main()
