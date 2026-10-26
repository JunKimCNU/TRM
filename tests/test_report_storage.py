"""Whitespace-only storage choice preserves observations and all report views."""
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'tests')]
from evidence_batch import write_json
import test_report_assessment as fixtures


class AtomicJsonStorageTests(unittest.TestCase):
    def test_default_encoding_unchanged_and_compact_roundtrip(self):
        data = {'unicode': '측정 안 됨', 'missing': None, 'zero': 0,
                'large': 2**80, 'nested': [{'positive': True, 'negative': False, 'fraction': 1/3}]}
        with tempfile.TemporaryDirectory() as folder:
            pretty, compact = [Path(folder) / name for name in ('pretty.json', 'compact.json')]
            write_json(pretty, data)
            write_json(compact, data, compact=True)
            self.assertEqual(pretty.read_text(), json.dumps(data, indent=2, ensure_ascii=False) + '\n')
            self.assertEqual(json.loads(pretty.read_bytes()), json.loads(compact.read_bytes()))
            self.assertLess(compact.stat().st_size, pretty.stat().st_size)
            self.assertTrue(compact.read_bytes().endswith(b'\n'))

    def test_failed_compact_encode_keeps_previous_file_and_cleans_temporary(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / 'result.json'
            write_json(output, {'completed': True})
            before = output.read_bytes()
            for data, error in (({'invalid': float('nan')}, ValueError), ({'invalid': object()}, TypeError)):
                with self.subTest(error=error), self.assertRaises(error):
                    write_json(output, data, compact=True)
                self.assertEqual(output.read_bytes(), before)
                self.assertEqual(list(Path(folder).iterdir()), [output])


class ReportStorageTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.AssessmentReportTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)

    def run_report(self, *args, success=True):
        proc = subprocess.run([*self.fixture.command(), *args], text=True, capture_output=True)
        self.assertEqual(proc.returncode == 0, success, proc.stdout + proc.stderr)

    def output_files(self):
        return {str(path.relative_to(self.fixture.root)): path.read_bytes()
                for root in (self.fixture.output, self.fixture.comparison)
                for path in root.rglob('*') if path.is_file()
                and not any(part.startswith('.') for part in path.relative_to(root).parts)
                and path.name not in {'interpretation_status.json', 'interpretation_resources.tsv'}}

    def checked_projection_content(self, files):
        """Verify every byte binding before comparing whitespace-independent values."""
        prefix = str(self.fixture.output.relative_to(self.fixture.root))
        manifest_name = prefix + '/summary_projections.tsv'
        metrics = json.loads(files[prefix + '/metrics.json'])
        self.assertEqual(metrics['files']['summary_projections.tsv'],
                         hashlib.sha256(files[manifest_name]).hexdigest())
        reader = csv.DictReader(io.StringIO(files[manifest_name].decode()), delimiter='\t')
        self.assertEqual(reader.fieldnames,
            ['sample_id', 'source_sha256', 'source_bytes', 'projection_sha256', 'projection_bytes'])
        rows = list(reader)
        projections = {name: json.loads(payload) for name, payload in files.items()
                       if name.startswith(prefix + '/summaries/') and name.endswith('.json')}
        self.assertEqual(len(rows), len(projections))
        self.assertEqual({row['sample_id'] for row in rows},
                         {Path(name).stem for name in projections})
        normalized = {}
        for entry in rows:
            sample = entry['sample_id']
            source_name = prefix + '/samples/' + sample + '.json'
            projection_name = prefix + '/summaries/' + sample + '.json'
            source, payload = files[source_name], files[projection_name]
            digest = hashlib.sha256(source).hexdigest()
            self.assertEqual(entry['source_sha256'], digest)
            self.assertEqual(int(entry['source_bytes']), len(source))
            self.assertEqual(entry['projection_sha256'], hashlib.sha256(payload).hexdigest())
            self.assertEqual(int(entry['projection_bytes']), len(payload))
            value = projections[projection_name]
            self.assertEqual(value['sample_id'], sample)
            self.assertEqual(value.pop('source_sha256'), digest)
            self.assertEqual(value.pop('source_bytes'), len(source))
            for row in [value['row'], *value['expected'], *value['queue']]:
                self.assertEqual(row['sample_id'], sample)
                self.assertEqual(row.pop('source_sha256'), digest)
                self.assertEqual(row['source_json'], str((self.fixture.root / source_name).resolve()))
            normalized[projection_name] = value
        return normalized, {row['sample_id']: row for row in rows}

    def test_pretty_compact_preserve_science_views_and_rebind_projection_hashes(self):
        inputs = [*self.fixture.source.rglob('*.json'), *self.fixture.evidence.rglob('*.json')]
        before = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}
        self.run_report()
        pretty = self.output_files()
        self.run_report('--json-format', 'compact')
        compact = self.output_files()
        self.assertEqual(set(pretty), set(compact))
        pretty_projections, pretty_manifest = self.checked_projection_content(pretty)
        compact_projections, compact_manifest = self.checked_projection_content(compact)
        self.assertEqual(pretty_projections, compact_projections)
        self.assertEqual(set(pretty_manifest), set(compact_manifest))
        for sample in pretty_manifest:
            old, new = pretty_manifest[sample], compact_manifest[sample]
            self.assertNotEqual(old['source_sha256'], new['source_sha256'])
            self.assertLess(int(new['source_bytes']), int(old['source_bytes']))
            self.assertNotEqual(old['projection_sha256'], new['projection_sha256'])
        for name in pretty:
            with self.subTest(file=name):
                if name.endswith('/metrics.json'):
                    a, b = json.loads(pretty[name]), json.loads(compact[name])
                    self.assertEqual(a['json_format'], 'pretty')
                    self.assertEqual(b['json_format'], 'compact')
                    for key in ('json_format', 'created_utc', 'report_wall_seconds'):
                        a.pop(key); b.pop(key)
                    # Exact hashes already checked against each manifest above.
                    self.assertNotEqual(a['files'].pop('summary_projections.tsv'),
                                        b['files'].pop('summary_projections.tsv'))
                    self.assertEqual(a, b)
                elif name in pretty_projections or name.endswith('/summary_projections.tsv'):
                    # All projection/manifest fields checked and science compared above.
                    continue
                elif name.endswith('.json'):
                    self.assertEqual(json.loads(pretty[name]), json.loads(compact[name]))
                    self.assertLess(len(compact[name]), len(pretty[name]))
                else:
                    self.assertEqual(pretty[name], compact[name])
        self.assertEqual(before, {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs})

    def test_invalid_format_rejected_before_creating_outputs(self):
        self.run_report('--json-format', 'gzip', success=False)
        self.assertFalse(self.fixture.output.exists())
        self.assertFalse(self.fixture.comparison.exists())


if __name__ == '__main__':
    unittest.main()
