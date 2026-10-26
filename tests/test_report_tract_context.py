"""Executable source-integrity and reporting contracts for optional diagnostics."""
from argparse import Namespace
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import report_tract_context as report


class ReaderTests(unittest.TestCase):
    def test_wrapped_records_keep_quality_and_source_positions(self):
        raw=b'@read note\nTTAGGG\nTTAGGG\n+read\nIIIIII\nHHHHHH\n'
        self.assertEqual(list(report.records_with_quality(io.BytesIO(raw),100)),
                         [('read','TTAGGGTTAGGG','IIIIIIHHHHHH')])

    def test_wrapped_fasta_and_iupac(self):
        raw=b'>a\n tta ggg\nNRY\n>b\nCCCTAA\n'
        self.assertEqual(list(report.records_with_quality(io.BytesIO(raw),100)),
                         [('a','TTAGGGNRY',None),('b','CCCTAA',None)])

    def test_truncation_mismatch_and_limits_reject(self):
        for raw,limit in [(b'@x\nTTAGGG\n+\nII\n',100),
                          (b'@x\nTTAGGG\n+y\nIIIIII\n',100),
                          (b'@x\nTTAGGG\n+\nIIIIIII\n',100),
                          (b'>x\nTTAGGG\nTTAGGG\n',10),
                          (b'>x\nTTAG!G\n',100)]:
            with self.subTest(raw=raw),self.assertRaises(ValueError):
                list(report.records_with_quality(io.BytesIO(raw),limit))


class CommandTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base=Path(self.tmp.name)
        seq='TTAGGG'*20
        self.source=self.base/'reads.fastq.gz'
        raw=('@read\n'+seq+'\n+\n'+'I'*len(seq)+'\n').encode()
        self.source.write_bytes(gzip.compress(raw,mtime=0))
        self.plan={'schema_version':1,'case_id':'<script>bad</script>',
            'data_role':'development_previously_seen','motif':'TTAGGG',
            'expected_source':{'sha256':hashlib.sha256(self.source.read_bytes()).hexdigest(),
                               'size_bytes':self.source.stat().st_size,'records':1,'bases':len(seq)},
            'regions':[{'id':'read:full','record_id':'read','record_length':len(seq),
                        'start0':0,'end0':len(seq),'expected_strand':'+',
                        'known_tracts':[{'tract_id':'left','start0':0,'end0':48,'strand':'+'},
                                        {'tract_id':'right','start0':72,'end0':120,'strand':'+'}]}],
            'pairs':[{'id':'p','record_id':'read','left_tract_id':'left','right_tract_id':'right',
                      'expected_strand':'+'}]}
        self.path=self.base/'plan.json'
        self.output=self.base/'out'

    def execute(self):
        self.path.write_text(json.dumps(self.plan))
        return subprocess.run([sys.executable,str(ROOT/'scripts/report_tract_context.py'),
            '--input',str(self.source),'--plan',str(self.path),'--output',str(self.output)],
            capture_output=True,text=True)

    def test_success_hashes_quality_metrics_no_classification_and_no_overwrite(self):
        result=self.execute()
        self.assertEqual(result.returncode,0,result.stderr)
        summary=json.loads((self.output/'summary.json').read_text())
        self.assertEqual(summary['source']['sha256'],self.plan['expected_source']['sha256'])
        self.assertTrue(summary['source']['validated_eof'])
        self.assertEqual(summary['input_bases'],120)
        self.assertFalse(summary['changes_native_classification'])
        self.assertFalse(summary['biological_tracts_merged'])
        self.assertIsInstance(summary['cutoffs'],list)
        for value in summary['outputs']:
            path=Path(value['path'])
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),value['sha256'])
        detail=json.loads((self.output/'details.jsonl').read_text())['analysis']
        self.assertTrue(detail['quality']['available'])
        self.assertEqual(detail['quality']['minimum'],40)
        page=(self.output/'index.html').read_text()
        self.assertNotIn('<script>bad</script>',page)
        self.assertIn('&lt;script&gt;bad&lt;/script&gt;',page)
        before=(self.output/'summary.json').read_bytes()
        again=self.execute()
        self.assertNotEqual(again.returncode,0)
        self.assertEqual((self.output/'summary.json').read_bytes(),before)

    def test_hash_mismatch_never_commits_success_directory(self):
        self.plan['expected_source']['sha256']='0'*64
        result=self.execute()
        self.assertNotEqual(result.returncode,0)
        self.assertFalse(self.output.exists())

    def test_missing_record_rejected(self):
        self.plan['regions'][0]['record_id']='missing'
        self.plan['pairs']=[]
        result=self.execute()
        self.assertNotEqual(result.returncode,0)
        self.assertFalse(self.output.exists())

    def test_conflicting_tract_metadata_rejected(self):
        extra=json.loads(json.dumps(self.plan['regions'][0]))
        extra['id']='another'
        extra['known_tracts'][0]['end0']=49
        self.plan['regions'].append(extra)
        result=self.execute()
        self.assertNotEqual(result.returncode,0)
        self.assertFalse(self.output.exists())

    def test_absent_mapping_not_reported_as_clean_mapping(self):
        result=self.execute()
        self.assertEqual(result.returncode,0,result.stderr)
        pair=json.loads((self.output/'pairs.json').read_text())[0]
        self.assertFalse(pair['physical_chromosome_end_confirmed'])
        self.assertEqual(pair['alignment_evidence']['status'],'not_provided')
        self.assertFalse(pair['alignment_evidence']['available'])

    def test_invalid_alignment_coordinates_rejected(self):
        self.plan['alignments']=[{'record_id':'read','query_start0':0,'query_end0':120,
            'target':'chr1','strand':'+','mapq':60,'tp':'P',
            'target_start0':400,'target_end0':520,'target_length':500}]
        result=self.execute()
        self.assertNotEqual(result.returncode,0)
        self.assertFalse(self.output.exists())

    def test_nonfinite_identity_rejected(self):
        self.plan['regions'][0]['known_tracts'][0]['identity']=float('nan')
        result=self.execute()
        self.assertNotEqual(result.returncode,0)
        self.assertFalse(self.output.exists())


if __name__=='__main__':
    unittest.main()
