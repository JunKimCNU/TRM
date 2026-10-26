"""Long-read overlays preserve native science and remain safe under report reuse."""
from copy import deepcopy
import csv
import html
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'tests')]
import long_read_context_evidence as module
import long_read_context_report as renderer
import report_runtime
from summarize_interpretation import summarize_sample
import test_report_assessment as fixture
import test_long_read_context_evidence as evidence_fixture

class LongReadContextReportingTests(unittest.TestCase):
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

    def save_evidence(self):
        saved=evidence_fixture.make_fixture(self.root/'long-read-fixture')
        self.manifest=saved['path']
        doc=json.loads(self.manifest.read_text());sample=doc['samples'][0]
        sample.update(sample_id='SRRknown',source_run='SRRknown',scientific_name='Fixture species',
            tax_id='6239',biosample='SAMNfixture_SRRknown',original_native_sha256=fixture.digest(self.source/'SRRknown/result.json'))
        sample['signature']['applicable_runs']=['SRRknown']
        self.pin=saved['reads']
        for key in ('reads','contexts','structures'):
            path=saved[key]
            with path.open(newline='') as handle:
                reader=csv.DictReader(handle,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
            for row in rows:
                row.update(sample_id='SRRknown',source_run='SRRknown')
                if 'read_id' in row:row['read_id']='SRRknown.'+row['read_id'].split('.',1)[1]
            evidence_fixture.write_tsv(path,fields,rows)
        # Keep every pinned source, including the external technical signature panel.
        module.create_manifest(self.manifest,[sample],[pin['path'] for pin in doc['source_artifacts']])
        return ['--long-read-context-evidence',str(self.manifest)]

    def modes(self):
        return {r['sample_id']:r['mode'] for r in fixture.read_tsv(self.output/'interpretation_resources.tsv')}

    def test_overlay_preserves_native_counts_candidates_policy_and_projection(self):
        self.run_report();before=self.document('SRRknown')
        unchanged={name:(self.output/name).read_bytes() for name in ('candidates.tsv','candidate_quality_observations.tsv','adequacy.tsv')}
        self.run_report(extra=self.save_evidence());after=self.document('SRRknown')
        self.assertEqual(set(before),set(after))
        for key in before:
            if key!='long_read_context_evidence':self.assertEqual(before[key],after[key],key)
        for name,data in unchanged.items():self.assertEqual(data,(self.output/name).read_bytes())
        self.assertEqual(summarize_sample(before,Path('fixed.json'),'0'*64),summarize_sample(after,Path('fixed.json'),'0'*64))
        ev=after['long_read_context_evidence']
        self.assertEqual(ev['status'],'annotated_saved_subset')
        self.assertEqual(ev['summary']['selected_reads'],3)
        self.assertIsNone(ev['physical_terminal_confirmed'])
        self.assertFalse(ev['native_counts_changed'])
        row=next(r for r in fixture.read_tsv(self.output/renderer.TABLES[0]) if r['sample_id']=='SRRknown')
        for field in module.SUMMARY_FIELDS:self.assertEqual(row[field],str(ev['summary'][field]))
        body=(self.output/'samples/SRRknown.html').read_text()
        for table in renderer.TABLES:self.assertIn('../'+table,body)

    def test_absent_evidence_and_unreviewed_context_remain_null_not_zero(self):
        self.run_report(extra=self.save_evidence())
        for sid in ('SRRnovel','SRRnone','SRRfailed'):
            ev=self.document(sid)['long_read_context_evidence']
            self.assertEqual(ev['status'],'unassessed')
            self.assertTrue(all(v is None for v in ev['summary'].values()))
        row=next(r for r in fixture.read_tsv(self.output/renderer.TABLES[0]) if r['sample_id']=='SRRnone')
        self.assertTrue(all(row[k]=='' for k in module.SUMMARY_FIELDS))
        row=next(r for r in fixture.read_tsv(self.output/renderer.TABLES[2]) if r['context_id']=='endB')
        self.assertEqual(row['structure_reviewed_read_count'],'')
        self.assertEqual(row['record_end_supported_read_count'],'')

    def test_serial_parallel_resume_and_evidence_removal(self):
        extra=self.save_evidence();self.run_report(extra=extra+['--jobs','1'])
        before={p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')}
        tables={n:(self.output/n).read_bytes() for n in renderer.TABLES}
        self.run_report(extra=extra+['--jobs','2'])
        self.assertEqual(before,{p.name:p.read_bytes() for p in (self.output/'samples').glob('*.json')})
        self.assertEqual(tables,{n:(self.output/n).read_bytes() for n in renderer.TABLES})
        self.run_report(extra=extra+['--jobs','2','--resume'])
        self.assertTrue(all(m=='reused' for m in self.modes().values()))
        doc=json.loads(self.manifest.read_text());doc['samples'][0]['notes']=['Additional provenance review']
        module.create_manifest(self.manifest,doc['samples'],[p['path'] for p in doc['source_artifacts']])
        self.run_report(extra=extra+['--resume'])
        self.assertTrue(all(m=='recomputed' for m in self.modes().values()))
        self.run_report(extra=['--resume'])
        self.assertTrue(all(m=='recomputed' for m in self.modes().values()))
        self.assertEqual(self.document('SRRknown')['long_read_context_evidence']['status'],'unassessed')
        for n in renderer.TABLES[1:]:self.assertEqual(fixture.read_tsv(self.output/n),[])

    def test_tampered_source_rejected_without_touching_warm_outputs(self):
        extra=self.save_evidence()+['--resume'];self.run_report(extra=extra)
        before={str(p.relative_to(self.output)):p.read_bytes() for p in self.output.rglob('*') if p.is_file()}
        self.pin.write_text(self.pin.read_text()+'\n')
        outcome=self.run_report(success=False,extra=extra)
        self.assertIn('SHA256 mismatch',outcome.stderr)
        self.assertEqual(before,{str(p.relative_to(self.output)):p.read_bytes() for p in self.output.rglob('*') if p.is_file()})

    def test_native_binding_mismatch_and_missing_manifest_rejected(self):
        outcome=self.run_report(success=False,extra=['--long-read-context-evidence',str(self.root/'absent.json')])
        self.assertIn('cannot read',outcome.stderr)
        extra=self.save_evidence();path=self.source/'SRRknown/result.json';path.write_text(path.read_text()+'\n')
        outcome=self.run_report(success=False,extra=extra)
        self.assertIn('SHA256 mismatch',outcome.stderr)

    def test_cache_requires_all_long_read_table_chunks(self):
        self.run_report(extra=self.save_evidence())
        receipt=json.loads((self.output/'.report-cache/SRRknown/receipt.json').read_text())
        self.assertTrue(report_runtime.verify_artifacts(receipt,output=self.output,sample='SRRknown'))
        path=Path(receipt['bundle']);bundle=json.loads(path.read_text())
        for name in renderer.TABLES:
            modified=deepcopy(bundle);modified['chunks'].remove(name);fixture.write_json(path,modified)
            altered=deepcopy(receipt);altered['artifacts'][str(path)]=report_runtime.binding(path)[0]
            altered['artifacts'].pop(str(path.parent/name))
            self.assertFalse(report_runtime.verify_artifacts(altered,output=self.output,sample='SRRknown'))
        fixture.write_json(path,bundle)

    def test_html_escapes_dynamic_text_and_disallows_script_source_links(self):
        self.save_evidence();ev=module.build_evidence('SRRknown',None,module.load_evidence(self.manifest))
        hostile='<img src=x onerror="alert(1)">'
        ev['source_run']=hostile;ev['contexts'][0]['target_record']=hostile
        ev['limitations'].append(hostile);ev['signature']['source_url']='javascript:alert(1)'
        output=''.join(renderer.render_html(ev))
        self.assertNotIn(hostile,output);self.assertIn(html.escape(hostile,quote=True),output)
        self.assertNotIn('href="javascript:',output)

if __name__=='__main__':unittest.main()
