"""Optional assembly context cannot alter biological results or cache provenance."""
from copy import deepcopy
import html
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'tests')]
import assembly_context_evidence as module
import report_assessment
import report_runtime
from summarize_interpretation import summarize_sample
import test_report_assessment as fixture
import test_assembly_context_evidence as evidence_fixture

TABLES=('assembly_context_summary.tsv','assembly_context_scaffolds.tsv',
        'assembly_context_intervals.tsv','assembly_context_events.tsv')


class AssemblyContextReportingTests(unittest.TestCase):
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
        self.saved=evidence_fixture.write_fixture(self.root/'assembly-fixture',
            sample_id='SRRknown',identity=dict(scientific_name='Fixture species',tax_id='6239',biosample='SAMNfixture_SRRknown'),
            native_sha=fixture.digest(self.source/'SRRknown/result.json'),representative_motif='TTAGGC',review_notes=review_notes)
        self.manifest=Path(self.saved['manifest_path']);self.payload_path=Path(self.saved['payload_path'])
        self.pin=Path(self.saved['source_path'])
        return ['--assembly-context-evidence',str(self.manifest)]

    def reseal_notes(self,note):
        payload=json.loads(self.payload_path.read_text());payload['review_notes'].append(note)
        payload['payload_sha256']=module.payload_sha256(payload);fixture.write_json(self.payload_path,payload)
        manifest=json.loads(self.manifest.read_text())
        manifest['samples'][0]['evidence']['sha256']=fixture.digest(self.payload_path)
        manifest['payload_sha256']=module.payload_sha256(manifest);fixture.write_json(self.manifest,manifest)

    def modes(self):
        return {r['sample_id']:r['mode'] for r in fixture.read_tsv(self.output/'interpretation_resources.tsv')}

    def test_overlay_changes_only_its_namespace_and_exports_evidence_levels(self):
        self.run_report();before=self.document('SRRknown')
        candidates=(self.output/'candidates.tsv').read_bytes()
        self.run_report(extra=self.save_evidence());after=self.document('SRRknown')
        self.assertEqual(set(before),set(after))
        for key in before:
            if key!='assembly_context_evidence':self.assertEqual(before[key],after[key],key)
        self.assertEqual(candidates,(self.output/'candidates.tsv').read_bytes())
        self.assertEqual(summarize_sample(before,Path('fixed.json'),'0'*64),summarize_sample(after,Path('fixed.json'),'0'*64))
        evidence=after['assembly_context_evidence'];summary=evidence['summary']
        self.assertEqual(evidence['status'],'sources_validated')
        self.assertEqual(after['detected_motif'],'TTAGGC')
        self.assertEqual(evidence['grammar_units'],['TTAGG','TTAGGG'])
        self.assertFalse(evidence['changes_classification_or_adequacy'])
        self.assertEqual(summary['physical_terminal_status'],'not_confirmed')
        self.assertEqual(summary['source_origin_status'],'not_confirmed')
        self.assertEqual(summary['quality_origin_status'],'unresolved_no_alternative_same_read_quality')
        self.assertEqual(summary['alternative_same_read_quality_comparisons'],0)
        self.assertTrue(all(v is False for v in evidence['confirmations'].values()))
        rendered=(self.output/'samples/SRRknown.html').read_text()
        for phrase in ('Read 관측 → assembly 확장 → 물리적 말단','assembly 정렬','정책상 대표 반복과 별개',
                       '물리적 염색체 말단 확인이 아닙니다','FCS 조치 행 0은 숙주','같은 read의 quality 비교가 아닙니다'):
            self.assertIn(phrase,rendered)
        for table in TABLES:self.assertIn('../'+table,rendered)
        exported=next(r for r in fixture.read_tsv(self.output/TABLES[0]) if r['sample_id']=='SRRknown')
        self.assertEqual(exported['read_mixed_intervals_preserved'],str(summary['read_mixed_intervals_preserved']))
        self.assertEqual(exported['supported_assembly_intervals'],str(summary['supported_assembly_intervals']))
        self.assertEqual(exported['physical_terminal_status'],'not_confirmed')
        self.assertEqual(exported['source_origin_status'],'not_confirmed')
        self.assertEqual(len(fixture.read_tsv(self.output/TABLES[1])),len(evidence['scaffolds']))
        self.assertEqual(len(fixture.read_tsv(self.output/TABLES[2])),len(evidence['intervals']))
        self.assertEqual(len(fixture.read_tsv(self.output/TABLES[3])),len(evidence['events']))

    def test_omission_and_other_samples_are_not_provided_with_null_measurements(self):
        self.run_report()
        for sid in ('SRRknown','SRRnone','SRRfailed'):
            d=self.document(sid)['assembly_context_evidence']
            self.assertEqual(d['status'],'not_provided')
            self.assertTrue(all(d['summary'][key] is None for key in module.SUMMARY_FIELDS))
            self.assertEqual(d['scaffolds'],[]);self.assertEqual(d['intervals'],[]);self.assertEqual(d['events'],[])
        self.run_report(extra=self.save_evidence())
        for sid in ('SRRnovel','SRRnone','SRRfailed'):
            self.assertEqual(self.document(sid)['assembly_context_evidence']['status'],'not_provided')
        empty=next(r for r in fixture.read_tsv(self.output/TABLES[0]) if r['sample_id']=='SRRnone')
        self.assertTrue(all(empty[key]=='' for key in module.SUMMARY_FIELDS))
        self.assertIn('관측 0건이나 TRM 부재를 뜻하지 않습니다',(self.output/'samples/SRRnone.html').read_text())

    def test_review_notes_and_dynamic_display_fields_are_html_escaped(self):
        hostile='<img src=x onerror="alert(1)"> & specimen'
        self.run_report(extra=self.save_evidence(review_notes=[hostile]))
        d=self.document('SRRknown')['assembly_context_evidence']
        self.assertIn(hostile,d['review_notes'])
        rendered=(self.output/'samples/SRRknown.html').read_text()
        self.assertNotIn(hostile,rendered);self.assertIn(html.escape(hostile,quote=True),rendered)
        # Renderer boundaries must remain safe even before validation rejects an invalid target.
        altered=deepcopy(d);altered['intervals'][0]['target']=hostile
        altered['display']['assembly_context_origin_display']=hostile
        rendered=''.join(report_assessment.assembly_context_html(altered))
        self.assertNotIn(hostile,rendered);self.assertIn(html.escape(hostile,quote=True),rendered)
        row=next(r for r in fixture.read_tsv(self.output/'samples.tsv') if r['sample_id']=='SRRknown')
        self.assertIn(hostile,row['assembly_context_limit_display'])

    def test_parallel_serial_and_unchanged_resume_are_identical(self):
        extra=self.save_evidence()
        self.run_report(extra=extra+['--jobs','1'])
        documents={p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')}
        tables={name:(self.output/name).read_bytes() for name in TABLES}
        self.run_report(extra=extra+['--jobs','2'])
        self.assertEqual(documents,{p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')})
        self.assertEqual(tables,{name:(self.output/name).read_bytes() for name in TABLES})
        self.run_report(extra=extra+['--jobs','2','--resume'])
        self.assertTrue(all(mode=='reused' for mode in self.modes().values()))
        self.assertEqual(documents,{p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')})

    def test_changed_source_cannot_hide_behind_completed_warm_cache(self):
        extra=self.save_evidence()+['--jobs','2','--resume'];self.run_report(extra=extra)
        # Validation occurs before a new report attempt. Its failed invocation
        # must not relabel or alter the preserved, previously valid artifact set.
        before={str(p.relative_to(self.output)):p.read_bytes() for p in self.output.rglob('*') if p.is_file()}
        self.pin.write_text(self.pin.read_text()+'\n')
        outcome=self.run_report(success=False,extra=extra)
        self.assertNotEqual(outcome.returncode,0)
        self.assertIn('source hash differs',outcome.stderr)
        self.assertEqual(before,{str(p.relative_to(self.output)):p.read_bytes() for p in self.output.rglob('*') if p.is_file()})
        self.assertEqual(json.loads((self.output/'interpretation_status.json').read_text())['status'],'complete')

    def test_changed_sealed_annotation_recomputes_and_option_removal_clears_rows(self):
        extra=self.save_evidence()+['--resume'];self.run_report(extra=extra)
        self.reseal_notes('New bounded annotation; no biological promotion.')
        self.run_report(extra=extra)
        self.assertTrue(all(mode=='recomputed' for mode in self.modes().values()))
        self.assertIn('New bounded annotation',str(self.document('SRRknown')['assembly_context_evidence']))
        self.run_report(extra=['--resume'])
        self.assertEqual(self.document('SRRknown')['assembly_context_evidence']['status'],'not_provided')
        self.assertTrue(all(mode=='recomputed' for mode in self.modes().values()))
        for name in TABLES[1:]:self.assertEqual(fixture.read_tsv(self.output/name),[])

    def test_changed_native_binding_and_explicit_missing_manifest_fail_closed(self):
        outcome=self.run_report(success=False,extra=['--assembly-context-evidence',str(self.root/'missing.json')])
        self.assertIn('source unavailable',outcome.stderr)
        extra=self.save_evidence()
        source=self.source/'SRRknown/result.json';source.write_text(source.read_text()+'\n')
        outcome=self.run_report(success=False,extra=extra)
        self.assertIn('native SHA differs',outcome.stderr)
        self.assertEqual(json.loads((self.output/'interpretation_status.json').read_text())['status'],'failed')

    def test_sealed_cache_requires_each_assembly_context_export(self):
        self.run_report(extra=self.save_evidence())
        receipt=json.loads((self.output/'.report-cache/SRRknown/receipt.json').read_text())
        self.assertTrue(report_runtime.verify_artifacts(receipt,output=self.output,sample='SRRknown'))
        bundle_path=Path(receipt['bundle']);bundle=json.loads(bundle_path.read_text())
        for name in TABLES:
            with self.subTest(name=name):
                modified=deepcopy(bundle);modified['chunks'].remove(name);fixture.write_json(bundle_path,modified)
                changed=deepcopy(receipt);changed['artifacts'][str(bundle_path)]=report_runtime.binding(bundle_path)[0]
                changed['artifacts'].pop(str(bundle_path.parent/name))
                self.assertFalse(report_runtime.verify_artifacts(changed,output=self.output,sample='SRRknown'))
        fixture.write_json(bundle_path,bundle)


if __name__=='__main__':unittest.main()
