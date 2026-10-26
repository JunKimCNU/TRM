"""Executable optional-model report contracts and immutable old-model results."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import phase_bridge as phase
import tract_context as previous


class ReportTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.base=Path(self.tmp.name);self.source=self.base/'reads.fastq.gz'
        self.seq=('TTAGGG'*20)[:97]
        content=('@r\n'+self.seq+'\n+\n'+'I'*len(self.seq)+'\n').encode()
        self.source.write_bytes(gzip.compress(content,mtime=0))
        self.left=dict(tract_id='a',start0=0,end0=48,strand='+')
        self.right=dict(tract_id='b',start0=49,end0=97,strand='+')
        self.plan=dict(schema_version=1,case_id='<script>alert(1)</script>',data_role='development_previously_seen',motif='TTAGGG',
            expected_source=dict(sha256=hashlib.sha256(self.source.read_bytes()).hexdigest(),
                size_bytes=self.source.stat().st_size,records=1,bases=len(self.seq)),
            regions=[dict(id='r:context',record_id='r',record_length=len(self.seq),start0=0,end0=len(self.seq),
                          known_tracts=[self.left,self.right])],
            pairs=[dict(id='r:pair',record_id='r',left_tract_id='a',right_tract_id='b',expected_strand='+')])
        self.output=self.base/'out';self.plan_path=self.base/'plan.json'

    def execute(self):
        self.plan_path.write_text(json.dumps(self.plan))
        return subprocess.run([sys.executable,str(ROOT/'scripts/report_phase_bridge.py'),
            '--input',str(self.source),'--plan',str(self.plan_path),'--output',str(self.output)],
            capture_output=True,text=True)

    def test_short_gap_transition_same_input_previous_unchanged_trace_and_bindings(self):
        expected=previous.evaluate_pair(self.seq,self.left,self.right,'TTAGGG',expected_strand='+')
        result=self.execute();self.assertEqual(result.returncode,0,result.stderr)
        record=json.loads((self.output/'details.jsonl').read_text())
        self.assertEqual(record['previous'],expected)
        self.assertEqual(record['previous']['status'],'rejected_by_bounded_model')
        self.assertEqual(record['same_range_previous']['status'],'rejected_by_bounded_model')
        self.assertEqual(record['same_range_previous']['config']['max_gap_bp'],1024)
        self.assertEqual(record['current']['status'],phase.SUCCESS)
        self.assertEqual(record['current']['bridge']['min_edits'],0)
        self.assertEqual(record['current']['bridge']['path_ops'][0]['query_start0'],48)
        with (self.output/'comparison.tsv').open() as f:row=next(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(row['minimum_edit_cost'],'0')
        self.assertEqual(row['edit_budget'],'1')
        summary=json.loads((self.output/'summary.json').read_text())
        self.assertTrue(summary['source']['validated_eof'])
        self.assertFalse(summary['changes_native_classification'])
        self.assertFalse(summary['biological_tracts_merged'])
        for binding in summary['outputs']:
            self.assertEqual(hashlib.sha256(Path(binding['path']).read_bytes()).hexdigest(),binding['sha256'])
        page=(self.output/'index.html').read_text()
        self.assertNotIn('<script>alert(1)</script>',page)
        self.assertIn('&lt;script&gt;alert(1)&lt;/script&gt;',page)

    def test_bad_hash_does_not_publish_success(self):
        self.plan['expected_source']['sha256']='0'*64
        result=self.execute();self.assertNotEqual(result.returncode,0)
        self.assertFalse(self.output.exists())

    def test_work_limit_is_reported_separately_from_rejection(self):
        self.plan['phase_config']={'max_dp_states':1}
        result=self.execute();self.assertEqual(result.returncode,0,result.stderr)
        current=json.loads((self.output/'details.jsonl').read_text())['current']
        self.assertEqual(current['status'],phase.COMPUTATION_LIMIT)
        self.assertIsNone(current['bridge']['min_edits'])

    def test_unknown_config_fails(self):
        self.plan['phase_config']={'unknown_option':0}
        result=self.execute();self.assertNotEqual(result.returncode,0)
        self.assertFalse(self.output.exists())

    def test_missing_selected_record_fails(self):
        self.plan['regions'][0]['record_id']='missing'
        self.plan['pairs'][0]['record_id']='missing'
        result=self.execute();self.assertNotEqual(result.returncode,0)
        self.assertFalse(self.output.exists())

    def test_existing_report_preserved(self):
        result=self.execute();self.assertEqual(result.returncode,0,result.stderr)
        before=(self.output/'summary.json').read_bytes()
        result=self.execute();self.assertNotEqual(result.returncode,0)
        self.assertEqual((self.output/'summary.json').read_bytes(),before)


if __name__=='__main__':unittest.main()
