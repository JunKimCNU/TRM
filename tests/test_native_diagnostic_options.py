"""Opt-in native diagnostics survive batch execution without changing TRM evidence."""
from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from native_diagnostic_options import add_arguments, native_diagnostic_options, validate_diagnostic_output
from quality_observation import describe_quality, describe_candidate_quality


def scientific_result(report):
    """Exclude execution/diagnostic metadata, retaining all scientific observations."""
    value = copy.deepcopy(report)
    for name in ('version', 'command', 'threads', 'targets_source', 'performance_profile', 'quality_diagnostics'):
        value.pop(name, None)
    for name in ('elapsed_seconds', 'megabases_per_second'):
        value['summary'].pop(name, None)
    value['memory_policy'].pop('startup_virtual_bytes', None)
    value['evidence_qc'].pop('quality_scores_assessed', None)
    for name in ('candidates', 'rejected_candidates', 'candidate_evidence'):
        for candidate in value.get(name, []):
            candidate.pop('quality_diagnostics', None)
    return value


class NativeDiagnosticOptionTests(unittest.TestCase):
    def test_absent_and_explicit_off_keep_legacy_option_fingerprint(self):
        parser = argparse.ArgumentParser()
        add_arguments(parser)
        for args in (SimpleNamespace(), parser.parse_args([]),
                     parser.parse_args(['--performance-profile', 'off', '--quality-diagnostics', 'off'])):
            self.assertEqual(native_diagnostic_options(args), [])
            validate_diagnostic_output({}, args)

    def test_quality_cap_requires_declared_encoding_and_bounded_integer(self):
        for quality, limit in [('off', 100), ('phred33', 0), ('phred33', -1),
                               ('phred33', 1_000_001), ('phred33', True), ('phred33', 1.5)]:
            with self.subTest(quality=quality, limit=limit), self.assertRaises(ValueError):
                native_diagnostic_options(SimpleNamespace(quality_diagnostics=quality,
                                                           quality_max_record_bases=limit))
        for limit in (1, 1_000_000):
            self.assertEqual(native_diagnostic_options(SimpleNamespace(
                quality_diagnostics='phred33', quality_max_record_bases=limit)),
                ['--quality-diagnostics', 'phred33', '--quality-max-record-bases', str(limit)])

    def test_requested_missing_diagnostics_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'performance profile is missing'):
            validate_diagnostic_output({}, SimpleNamespace(performance_profile='on'))
        with self.assertRaisesRegex(ValueError, 'quality diagnostics are missing'):
            validate_diagnostic_output({}, SimpleNamespace(quality_diagnostics='phred33'))


class NativeDiagnosticIntegrationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='trm-native-diagnostic-options-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.binary = runtime_binary()
        self.sample = 'SRRdiagnostic'
        self.manifest = self.root / 'samples.tsv'
        self.manifest.write_text('sample_id\tpath1\tpath2\nSRRdiagnostic\tR1.fastq\tR2.fastq\n')
        for mate, motif in ((1, 'TTAGGG'), (2, 'CCCTAA')):
            sequences = [motif * 25, motif * 25, motif * 16 + 'ACGT' * 13 + 'AC', 'ACGT' * 37 + 'AC']
            (self.root / f'R{mate}.fastq').write_text(''.join(
                f'@r{i}/{mate}\n{seq}\n+\n{("!" if i == 0 else "I") * len(seq)}\n'
                for i, seq in enumerate(sequences)))
        self.on = ['--performance-profile', 'on', '--quality-diagnostics', 'phred33',
                   '--quality-max-record-bases', '200']

    def command(self, script, *args):
        result = subprocess.run([sys.executable, '-B', str(ROOT / 'scripts' / script),
                                 *map(str, args)], capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def discovery(self, output, *options):
        self.command('batch.py', self.manifest, '--binary', self.binary, '--output-dir', output,
                     '--jobs', 1, '--threads', 1, '--cpu-budget', 1, '--paired-layout', 'fr', *options)
        return json.loads((output / self.sample / 'result.json').read_text())

    def evidence(self, source, output, *options):
        self.command('evidence_batch.py', '--manifest', self.manifest, '--source', source,
                     '--binary', self.binary, '--output', output, '--jobs', 1, '--threads', 1,
                     '--cpu-budget', 1, '--paired-layout', 'fr', '--catalog', ROOT / 'data/telomere_motifs.tsv',
                     *options)
        return json.loads((output / self.sample / 'evidence.json').read_text())

    def assert_measured(self, report, candidate_key):
        validate_diagnostic_output(report, SimpleNamespace(performance_profile='on',
            quality_diagnostics='phred33', quality_max_record_bases=200))
        context = describe_quality(report, 'fixture-source-sha')
        self.assertEqual(context['status'], 'measured')
        records = context['diagnostics']['input_records']
        self.assertEqual(records['quality_records'], 8)
        self.assertEqual(records['whole_read']['bases'], 1200)
        self.assertEqual(records['whole_read']['below_q20_bases'], 300)
        rows = [row for row in report[candidate_key] if row['reads']]
        self.assertTrue(rows)
        for row in rows:
            observed = describe_candidate_quality(row, context)
            self.assertIsNotNone(observed)
            self.assertEqual(observed['diagnostics']['union_tract']['bases'], row['repeat_bases'])
            self.assertFalse(observed['filtering_applied'])

    def test_discovery_default_off_resume_and_requested_diagnostics_preserve_counts(self):
        output = self.root / 'discovery'
        off = self.discovery(output)
        state_path = output / self.sample / 'state.json'
        state_before = state_path.read_bytes()
        result_before = (output / self.sample / 'result.json').read_bytes()
        self.discovery(output, '--resume', '--quality-diagnostics', 'off', '--performance-profile', 'off')
        self.assertEqual(state_path.read_bytes(), state_before)
        self.assertEqual((output / self.sample / 'result.json').read_bytes(), result_before)
        on = self.discovery(output, '--resume', *self.on)
        self.assertEqual(scientific_result(off), scientific_result(on))
        self.assert_measured(on, 'candidates')
        result_before = (output / self.sample / 'result.json').read_bytes()
        state_before = state_path.read_bytes()
        self.discovery(output, '--resume', *self.on)
        self.assertEqual(state_path.read_bytes(), state_before)
        self.assertEqual((output / self.sample / 'result.json').read_bytes(), result_before)
        with self.assertRaisesRegex(ValueError, 'retention limit differs'):
            validate_diagnostic_output(on, SimpleNamespace(quality_diagnostics='phred33',
                                                           quality_max_record_bases=201))

    def test_targeted_evidence_flags_resume_and_quality_scope(self):
        source, output = self.root / 'discovery', self.root / 'evidence'
        self.discovery(source)
        off = self.evidence(source, output)
        default_status = (output / self.sample / 'status.json').read_bytes()
        self.evidence(source, output, '--resume', '--performance-profile', 'off', '--quality-diagnostics', 'off')
        self.assertEqual((output / self.sample / 'status.json').read_bytes(), default_status)
        on = self.evidence(source, output, '--resume', *self.on)
        self.assertEqual(scientific_result(off), scientific_result(on))
        self.assert_measured(on, 'candidate_evidence')
        status_path = output / self.sample / 'status.json'
        before_status = status_path.read_bytes()
        before_result = (output / self.sample / 'evidence.json').read_bytes()
        self.evidence(source, output, '--resume', *self.on)
        self.assertEqual(status_path.read_bytes(), before_status)
        self.assertEqual((output / self.sample / 'evidence.json').read_bytes(), before_result)

    def test_cli_invalid_cap_fails_before_discovery_output_creation(self):
        output = self.root / 'invalid'
        result = subprocess.run([sys.executable, '-B', str(ROOT / 'scripts/batch.py'), str(self.manifest),
            '--binary', str(self.binary), '--output-dir', str(output), '--quality-max-record-bases', '10'],
            capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 2)
        self.assertIn('requires phred33', result.stderr)
        self.assertFalse(output.exists())

    def test_native_quality_is_retained_in_assessment_with_exact_source_scope(self):
        source, evidence = self.root / 'discovery', self.root / 'evidence'
        self.discovery(source)
        self.evidence(source, evidence)
        metadata = self.root / 'metadata.tsv'
        assignments = self.root / 'assignments.tsv'
        row = {'sample_id': self.sample, 'scientific_name': 'Fixture species', 'tax_id': '1',
               'lineage': 'Eukaryota;', 'library_strategy': 'WGS', 'library_source': 'GENOMIC',
               'library_layout': 'PAIRED', 'status': 'resolved'}
        for path, record in ((metadata, row), (assignments, {**row, 'metadata_status': 'resolved',
                'mapping_status': 'global_prior_no_curated_reference', 'related_motifs': ''})):
            with path.open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(record), delimiter='\t')
                writer.writeheader()
                writer.writerow(record)
        def assessment(label):
            output = self.root / ('assessment-' + label)
            self.command('report_assessment.py', '--source', source, '--evidence', evidence,
                         '--metadata', metadata, '--assignments', assignments, '--output', output,
                         '--comparison-output', self.root / ('comparison-' + label), '--jobs', 1)
            return output, json.loads((output / 'samples' / f'{self.sample}.json').read_text())
        _, baseline = assessment('off')
        self.evidence(source, evidence, '--resume', *self.on)
        output, report = assessment('on')
        quality = report['quality_observation']
        self.assertEqual(quality['current_measurement']['status'], 'measured')
        self.assertEqual(quality['original_discovery']['status'], 'not_measured')
        self.assertIsNone(quality['original_discovery']['diagnostics'])
        self.assertEqual(quality['current_measurement']['diagnostics']['input_records']['quality_records'], 8)
        self.assertFalse(quality['current_measurement']['changes_classification_or_adequacy'])
        self.assertEqual(report['classification'], baseline['classification'])
        self.assertEqual(report['data_adequacy']['status'], baseline['data_adequacy']['status'])
        before_candidates, after_candidates = copy.deepcopy(baseline['candidates']), copy.deepcopy(report['candidates'])
        for candidate in before_candidates + after_candidates:
            candidate.pop('quality_diagnostics', None)
        self.assertEqual(before_candidates, after_candidates)
        before_qc, after_qc = baseline['data_adequacy']['read_qc'], report['data_adequacy']['read_qc']
        self.assertFalse(before_qc['base_quality_distribution_available'])
        self.assertEqual(before_qc['base_quality_observation_status'], 'not_measured')
        self.assertTrue(after_qc['base_quality_distribution_available'])
        self.assertEqual(after_qc['base_quality_observation_status'], 'measured')
        self.assertFalse(after_qc['base_quality_used_for_power'])
        with (output / 'candidate_quality_observations.tsv').open(newline='') as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        self.assertTrue(rows)
        self.assertTrue(all(row['source_sha256'] == quality['current_measurement']['source_sha256'] for row in rows))


if __name__ == '__main__':
    unittest.main()
