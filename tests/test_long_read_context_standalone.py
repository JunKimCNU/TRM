"""Real CLI integration checks using synthetic saved annotations only."""
import csv
import hashlib
import html
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import long_read_context_evidence as evidence
from test_long_read_context_evidence import make_fixture, write_tsv


class LongReadContextStandaloneTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.fixture = make_fixture(self.directory / 'inputs', native_bound=False)
        self.sample_id = self.fixture['sample']['sample_id']
        self.output = self.directory / 'report'

    def command(self, output=None, *extra):
        return subprocess.run([
            sys.executable, '-B', str(ROOT / 'scripts/report_long_read_context.py'),
            '--evidence', str(self.fixture['path']),
            '--output', str(output or self.output), *extra,
        ], cwd=ROOT, text=True, encoding='utf-8', capture_output=True, timeout=30)

    def table(self, name):
        with (self.output / name).open(newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            self.assertEqual(len(reader.fieldnames), len(set(reader.fieldnames)))
            return list(reader)

    @staticmethod
    def snapshot(directory):
        return {str(path.relative_to(directory)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in directory.rglob('*') if path.is_file()}

    def rewrite_table(self, name, edit):
        path = self.fixture[name]
        with path.open(newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields, rows = reader.fieldnames, list(reader)
        edit(rows)
        write_tsv(path, fields, rows)

    def reseal(self):
        return evidence.create_manifest(self.fixture['path'], [self.fixture['sample']],
                                        [self.fixture[key] for key in ('reads', 'contexts', 'structures', 'signature')])

    def test_success_preserves_unprovided_policy_and_derives_small_tables(self):
        before = self.snapshot(self.directory / 'inputs')
        process = self.command()
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(json.loads(process.stdout)['status'], 'complete_evidence_only_report')
        document = json.loads((self.output / 'samples' / (self.sample_id + '.json')).read_text())
        for key in ('classification', 'detected_motif', 'data_adequacy', 'policy_assessment'):
            self.assertIsNone(document[key], key)
        self.assertIs(document['native_discovery_performed'], False)
        self.assertEqual(document['native_discovery_status'], 'not_provided')
        saved = document['long_read_context_evidence']
        self.assertEqual(saved['summary']['selected_reads'], 3)
        self.assertEqual(saved['summary']['whole_tail_signature'], 1)
        self.assertEqual(saved['summary']['partial_signature_ge20_at_read_end'], 1)
        self.assertIs(saved['native_counts_changed'], False)
        self.assertIs(saved['sequence_trimmed'], False)
        self.assertIsNone(saved['physical_terminal_confirmed'])
        self.assertIsNone(saved['main_trm_confirmed'])
        self.assertIsNone(saved['novel_trm_confirmed'])
        summary = self.table('long_read_context_summary.tsv')
        self.assertEqual(len(summary), 1)
        self.assertEqual(summary[0]['selected_reads'], '3')
        self.assertEqual(summary[0]['source_run'], 'SRR1')
        reads = self.table('long_read_context_reads.tsv')
        self.assertEqual(len(reads), 3)
        self.assertEqual({row['read_id'] for row in reads}, {'SRR1.1', 'SRR1.2', 'SRR1.3'})
        contexts = self.table('long_read_context_contexts.tsv')
        self.assertEqual(len(contexts), 2)
        unreviewed = next(row for row in contexts if row['context_id'] == 'endB')
        self.assertEqual(unreviewed['structure_reviewed_read_count'], '')
        self.assertEqual(unreviewed['record_end_supported_read_count'], '')
        json_unreviewed = next(row for row in saved['contexts'] if row['context_id'] == 'endB')
        self.assertIsNone(json_unreviewed['record_end_supported_read_count'])
        page = (self.output / 'samples' / (self.sample_id + '.html')).read_text()
        self.assertIn('미평가', page)
        self.assertIn('일반 native 탐색', page)
        completion = json.loads((self.output / 'completion.json').read_text())
        self.assertEqual(completion['native_scan_runs'], 0)
        self.assertEqual(completion['mapping_runs'], 0)
        self.assertIs(completion['classification_provided'], False)
        for artifact in completion['artifacts']:
            actual = (self.output / artifact['path']).read_bytes()
            self.assertEqual(len(actual), artifact['bytes'])
            self.assertEqual(hashlib.sha256(actual).hexdigest(), artifact['sha256'])
        self.assertEqual(self.snapshot(self.directory / 'inputs'), before)

    def test_html_escapes_taxon_reference_and_source_url(self):
        name = 'Example <script>alert("taxon")</script> & species'
        target = 'chr<script>alert("reference")</script>'
        url = 'https://example.org/primary?x=" onmouseover="alert(1)&y=<tag>'
        self.fixture['sample']['scientific_name'] = name
        self.fixture['sample']['signature']['source_url'] = url
        panel_path = self.fixture['signature']
        panel = json.loads(panel_path.read_text())
        panel['signatures'][0]['source_url'] = url
        panel_path.write_text(json.dumps(panel))
        self.rewrite_table('contexts', lambda rows: rows[0].update(target_record=target))
        self.reseal()
        process = self.command()
        self.assertEqual(process.returncode, 0, process.stderr)
        index = (self.output / 'index.html').read_text()
        page = (self.output / 'samples' / (self.sample_id + '.html')).read_text()
        self.assertIn(html.escape(name, quote=True), index)
        self.assertNotIn('<script>', index)
        self.assertIn(html.escape(target, quote=True), page)
        self.assertIn('href="' + html.escape(url, quote=True) + '"', page)
        self.assertNotIn('<script>', page)
        self.assertNotIn(' onmouseover="alert(1)', page)

    def test_existing_output_is_byte_for_byte_preserved(self):
        first = self.command()
        self.assertEqual(first.returncode, 0, first.stderr)
        before = self.snapshot(self.output)
        second = self.command()
        self.assertNotEqual(second.returncode, 0)
        self.assertIn('Output already exists', second.stderr)
        self.assertEqual(self.snapshot(self.output), before)

    def test_invalid_source_fails_before_creating_output_parent(self):
        source = self.fixture['reads']
        source.write_text(source.read_text().replace('SRR1.1', 'SRR1.999'))
        output = self.directory / 'not-created' / 'report'
        process = self.command(output)
        self.assertNotEqual(process.returncode, 0)
        self.assertIn('SHA256 mismatch', process.stderr)
        self.assertFalse(output.parent.exists())

    def test_unsafe_sample_id_cannot_be_used_as_output_path(self):
        self.fixture['sample']['sample_id'] = '../escape'
        for key in ('reads', 'contexts', 'structures'):
            self.rewrite_table(key, lambda rows: [row.update(sample_id='../escape') for row in rows])
        self.reseal()
        process = self.command()
        self.assertNotEqual(process.returncode, 0)
        self.assertIn('Unsafe sample identifier', process.stderr)
        self.assertFalse(self.output.exists())
        self.assertFalse((self.directory / 'escape.json').exists())
        self.assertFalse((self.directory / 'escape.html').exists())

    def test_output_containing_inputs_is_rejected_without_modification(self):
        # Existing input ancestors necessarily already exist. Either the stronger
        # no-overwrite gate or the explicit containment gate must refuse them.
        before = self.snapshot(self.directory)
        process = self.command(self.directory)
        self.assertNotEqual(process.returncode, 0)
        self.assertTrue('Output already exists' in process.stderr or
                        'must not contain evidence input files' in process.stderr,
                        process.stderr)
        self.assertEqual(self.snapshot(self.directory), before)
        self.assertFalse((self.directory / 'index.html').exists())

    def test_output_symlink_cannot_replace_or_modify_existing_directory(self):
        target = self.directory / 'unrelated-existing'
        target.mkdir()
        (target / 'keep.txt').write_text('preserved original data\n')
        self.output.symlink_to(target, target_is_directory=True)
        before = self.snapshot(target)
        process = self.command()
        self.assertNotEqual(process.returncode, 0)
        self.assertIn('Output already exists', process.stderr)
        self.assertIs(self.output.is_symlink(), True)
        self.assertEqual(self.snapshot(target), before)

    def test_unknown_requested_sample_fails_before_output(self):
        process = self.command(None, '--sample', 'missing-run')
        self.assertNotEqual(process.returncode, 0)
        self.assertIn('Requested sample absent', process.stderr)
        self.assertFalse(self.output.exists())


if __name__ == '__main__':
    unittest.main()
