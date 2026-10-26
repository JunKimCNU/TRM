"""Bounded target selection and content-aware evidence resume integration."""
from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import random
import subprocess
import sys
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch


from runtime_test_support import runtime_binary

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('evidence_batch_under_test', ROOT/'scripts/evidence_batch.py')
batch = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(batch)


def candidate(motif, reads=20, full=5, bases=1000):
    return {'motif': batch.canonical(motif), 'reads': reads, 'full_repeat_reads': full,
            'repeat_bases': bases, 'max_tract': 150}


class TargetSelectionTests(unittest.TestCase):
    def test_user_negative_examples_follow_same_generic_target_rule(self):
        motifs = ['AAAAAAAG', 'AAT', 'ACGCGGCG', 'ACAT', 'AAAT']
        report = {'candidates': [], 'rejected_candidates': [candidate(m) for m in motifs]}
        targets = batch.select_targets(report, {batch.canonical('TTAGGG')})
        self.assertEqual(set(targets['motifs']), {batch.canonical(m) for m in motifs + ['TTAGGG']})
        self.assertEqual(targets['selection']['novel_targets'], 5)
        self.assertEqual(targets['selection']['novel_candidates_not_targeted'], 0)

    def test_generic_support_and_length_requirements_apply_to_every_novel_family(self):
        rows = [candidate('ACAT', reads=4), candidate('AAAT', full=1), candidate('AT'),
                candidate('A'*32+'C'), candidate('AAAAAAAG')]
        targets = batch.select_targets({'candidates': rows}, set())
        self.assertEqual(targets['motifs'], [batch.canonical('AAAAAAAG')])
        self.assertEqual(targets['selection']['novel_screening_candidates'], 1)

    def test_cap_preserves_all_known_and_selects_strongest_novel_deterministically(self):
        known = {batch.canonical('TTAGGG'), batch.canonical('TTAGGC')}
        rows = [candidate(m, bases=b) for m, b in [('AAAT', 200), ('ACAT', 100), ('AAAAAAAG', 300)]]
        report = {'candidates': rows, 'rejected_candidates': [candidate('TTAGGG', reads=1, full=0)]}
        targets = batch.select_targets(report, known, max_novel=2)
        reverse = batch.select_targets({'candidates': list(reversed(rows))}, known, max_novel=2)
        self.assertEqual(targets, reverse)
        self.assertEqual(set(targets['motifs']), known | {batch.canonical('AAAAAAAG'), batch.canonical('AAAT')})
        self.assertEqual(targets['selection']['novel_candidates_not_targeted'], 1)
        zero = batch.select_targets(report, known, max_novel=0)
        self.assertEqual(set(zero['motifs']), known)
        self.assertEqual(zero['selection']['novel_candidates_not_targeted'], 3)

    def test_family_deduplication_and_scanner_limit(self):
        row = candidate('TTAGGG')
        duplicate = {**row, 'motif': 'CCCTAA'}
        target = batch.select_targets({'candidates': [row], 'rejected_candidates': [duplicate]}, set())
        self.assertEqual(target['motifs'], ['AACCCT'])
        rng = random.Random(97)
        known = set()
        while len(known) < 257:
            known.add(batch.canonical(''.join(rng.choices('ACGT', k=20))))
        with self.assertRaisesRegex(ValueError, '256'):
            batch.select_targets({'candidates': []}, known)


class ManifestInputTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='trm-evidence-manifest-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.manifest = self.root/'metadata'/'samples.tsv'
        self.manifest.parent.mkdir()
        self.input = self.root/'reads.fastq'
        self.input.write_text('@r\nTTAGGG\n+\nIIIIII\n')

    def normalized(self, path1, path2=''):
        return batch.normalize_manifest_inputs(
            [{'sample_id': 'sample', 'path1': path1, 'path2': path2}], self.manifest)[0]

    def test_relative_paths_use_manifest_parent_and_input_rows_are_preserved(self):
        rows = [{'sample_id': 'sample', 'path1': '../reads.fastq', 'custom': 'retained'}]
        result = batch.normalize_manifest_inputs(rows, self.manifest)
        self.assertEqual(result[0]['path1'], str(self.input.resolve()))
        self.assertEqual(result[0]['path2'], '')
        self.assertEqual(result[0]['custom'], 'retained')
        self.assertEqual(rows[0]['path1'], '../reads.fastq')
        self.assertEqual(self.normalized(str(self.input))['path1'], str(self.input.resolve()))

    def test_path1_is_required_even_when_path2_is_present(self):
        for missing in (None, '', '   '):
            with self.subTest(path1=missing), self.assertRaisesRegex(ValueError, 'path1 must not be empty'):
                self.normalized(missing, str(self.input))

    def test_duplicate_paths_symlinks_and_hardlinks_are_rejected(self):
        symlink = self.root/'symlink.fastq'
        symlink.symlink_to(self.input)
        hardlink = self.root/'hardlink.fastq'
        os.link(self.input, hardlink)
        for alias in (self.input, symlink, hardlink):
            with self.subTest(alias=alias.name), self.assertRaisesRegex(ValueError, 'same file'):
                self.normalized('../reads.fastq', str(alias))

    def test_directories_and_fifo_streams_are_rejected(self):
        stream = self.root/'input.fifo'
        os.mkfifo(stream)
        for path in (self.root, stream):
            with self.subTest(path=path), self.assertRaisesRegex(ValueError, 'regular file'):
                self.normalized(str(path))
            with self.subTest(path2=path), self.assertRaisesRegex(ValueError, 'regular file'):
                self.normalized(str(self.input), str(path))

    def test_missing_files_remain_normalized_rows_and_dash_is_not_stdin(self):
        row = self.normalized('../missing.fastq', '-')
        self.assertEqual(row['path1'], str(self.root/'missing.fastq'))
        self.assertEqual(row['path2'], str(self.manifest.parent/'-'))

    def test_invalid_manifest_fails_before_output_creation_or_scheduling(self):
        self.manifest.write_text('sample_id\tpath1\tpath2\nsample\t\t\n')
        output = self.root/'output'
        command = [sys.executable, str(ROOT/'scripts/evidence_batch.py'), '--manifest', str(self.manifest),
                   '--source', str(self.root/'source'), '--output', str(output), '--binary', '/nonexistent/trm',
                   '--jobs', '1', '--threads', '1', '--cpu-budget', '1']
        completed = subprocess.run(command, text=True, capture_output=True, timeout=10)
        self.assertEqual(completed.returncode, 2, completed.stderr)
        self.assertIn('path1 must not be empty', completed.stderr)
        self.assertFalse(output.exists())


class EvidenceResumeTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='trm-evidence-resume-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.input = self.root/'reads.fastq'
        sequences = ['TTAGGG'*25, 'CCCTAA'*25, 'ACGT'*37+'AC', 'N'*150]
        self.input.write_text(''.join(f'@r{i}\n{s}\n+\n{"I"*len(s)}\n' for i, s in enumerate(sequences)))
        self.args = SimpleNamespace(source=self.root/'source', output=self.root/'output', resume=True,
                                    threads=1, max_records=0, max_novel_targets=96)
        self.sample = 'SRRfixture'
        self.item = {'sample_id': self.sample, 'path1': str(self.input), 'path2': ''}
        self.source = self.args.source/self.sample/'result.json'
        self.source.parent.mkdir(parents=True)
        batch.write_json(self.source, {'summary': {'records': 4, 'bases': 600},
                                      'candidates': [candidate('TTAGGG', reads=2, full=2, bases=300)],
                                      'rejected_candidates': []})
        self.binary = runtime_binary()
        self.known = {batch.canonical('TTAGGG'), batch.canonical('TTAGGC')}

    def run_sample(self):
        status = batch.run_one(self.item, self.args, self.known, self.binary)
        self.assertEqual(status['status'], 'complete', status)
        return status

    def test_success_records_identity_and_preserves_source(self):
        before = self.source.read_bytes()
        status = self.run_sample()
        self.assertEqual(self.source.read_bytes(), before)
        self.assertEqual(status['identity']['source_sha256'], batch.sha256(self.source))
        self.assertEqual(status['identity']['binary_sha256'], batch.sha256(self.binary))
        output = self.args.output/self.sample/'evidence.json'
        self.assertEqual(status['result_sha256'], batch.sha256(output))
        self.assertEqual((status['records'], status['bases']), (4, 600))
        self.assertEqual(status['identity']['targets']['selection']['known_family_count'], 2)
        self.assertGreater(status['max_rss_kb'], 0)
        self.assertEqual(status['identity']['inputs'][0]['bytes'], self.input.stat().st_size)

    def test_identical_resume_does_not_start_scanner(self):
        first = self.run_sample()
        with patch.object(batch.subprocess, 'run', side_effect=AssertionError('resume must not rerun')):
            second = self.run_sample()
        self.assertEqual(first, second)

    def test_changed_source_content_invalidates_resume(self):
        first = self.run_sample()
        value = json.loads(self.source.read_text())
        value['review_note'] = 'Provenance changed without changing read totals.'
        batch.write_json(self.source, value)
        original_run = batch.subprocess.run
        with patch.object(batch.subprocess, 'run', wraps=original_run) as invoked:
            second = self.run_sample()
        self.assertEqual(invoked.call_count, 1)
        self.assertNotEqual(first['identity']['source_sha256'], second['identity']['source_sha256'])

    def test_corrupt_saved_evidence_invalidates_resume(self):
        first = self.run_sample()
        output = self.args.output/self.sample/'evidence.json'
        output.write_text(output.read_text()+'\n')
        stale_hash = batch.sha256(output)
        self.assertNotEqual(stale_hash, first['result_sha256'])
        original_run = batch.subprocess.run
        with patch.object(batch.subprocess, 'run', wraps=original_run) as invoked:
            second = self.run_sample()
        self.assertEqual(invoked.call_count, 1)
        self.assertEqual(second['result_sha256'], batch.sha256(output))
        self.assertNotEqual(second['result_sha256'], stale_hash)

    def test_changed_target_set_invalidates_resume(self):
        first = self.run_sample()
        self.known.add(batch.canonical('TTAGAC'))
        original_run = batch.subprocess.run
        with patch.object(batch.subprocess, 'run', wraps=original_run) as invoked:
            second = self.run_sample()
        self.assertEqual(invoked.call_count, 1)
        self.assertNotEqual(first['identity']['targets'], second['identity']['targets'])

    def test_original_failed_sample_does_not_open_missing_fastq(self):
        self.source.unlink()
        self.input.unlink()
        with patch.object(batch.subprocess, 'run', side_effect=AssertionError('must not run')):
            status = batch.run_one(self.item, self.args, self.known, self.binary)
        self.assertEqual(status['status'], 'source_failed')
        self.assertIn('no biological nondetection', status['reason'])

    def test_complete_scan_count_mismatch_is_failed_not_a_valid_result(self):
        source = json.loads(self.source.read_text())
        source['summary']['records'] = 5
        batch.write_json(self.source, source)
        status = batch.run_one(self.item, self.args, self.known, self.binary)
        self.assertEqual(status['status'], 'failed')
        self.assertIn('records differs', status['reason'])

    def test_prefix_scan_carries_explicit_prefix_identity(self):
        self.args.max_records = 1
        status = self.run_sample()
        self.assertEqual(status['identity']['max_records'], 1)
        evidence = json.loads((self.args.output/self.sample/'evidence.json').read_text())
        self.assertTrue(evidence['summary']['sampled_prefix'])
        self.assertEqual(evidence['summary']['records'], 1)
        self.assertEqual(status['input_bytes'], 0)
        self.assertEqual(status['input_artifact_bytes_total'], self.input.stat().st_size)
        self.assertFalse(status['compressed_bytes_processed_measured'])

    def test_direct_call_without_path1_does_not_start_subprocess(self):
        for path2 in ('', str(self.input)):
            with self.subTest(path2=path2), patch.object(batch.subprocess, 'run', side_effect=AssertionError('must not run')):
                status = batch.run_one({'sample_id': self.sample, 'path1': '', 'path2': path2},
                                       self.args, self.known, self.binary)
                self.assertEqual(status['status'], 'failed')
                self.assertIn('path1 must not be empty', status['reason'])
                self.assertEqual(json.loads((self.args.output/self.sample/'status.json').read_text()), status)

    def run_batch_from_other_directory(self, max_records=0):
        manifest = self.root/'metadata'/'samples.tsv'
        manifest.parent.mkdir(exist_ok=True)
        manifest.write_text(f'sample_id\tpath1\n{self.sample}\t../reads.fastq\n')
        elsewhere = self.root/'elsewhere'
        elsewhere.mkdir(exist_ok=True)
        command = [sys.executable, str(ROOT/'scripts/evidence_batch.py'), '--manifest', str(manifest),
                   '--source', str(self.args.source), '--catalog', str(ROOT/'data/telomere_motifs.tsv'),
                   '--binary', str(self.binary), '--output', str(self.args.output), '--jobs', '1',
                   '--threads', '1', '--max-records', str(max_records)]
        completed = subprocess.run(command, cwd=elsewhere, text=True, capture_output=True, timeout=20)
        return completed, json.loads((self.args.output/'metrics.json').read_text())

    def test_batch_from_other_cwd_resolves_manifest_relative_input(self):
        completed, metrics = self.run_batch_from_other_directory()
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(metrics['completed'], 1)
        self.assertEqual(metrics['records'], 4)
        status = json.loads((self.args.output/self.sample/'status.json').read_text())
        self.assertEqual(status['identity']['inputs'][0]['path'], str(self.input.resolve()))
        self.assertTrue(status['compressed_bytes_processed_measured'])
        self.assertGreater(metrics['compressed_MB_per_second'], 0)

    def test_prefix_batch_does_not_estimate_processed_compressed_throughput(self):
        completed, metrics = self.run_batch_from_other_directory(max_records=1)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(metrics['records'], 1)
        self.assertEqual(metrics['input_bytes'], 0)
        self.assertIsNone(metrics['compressed_MB_per_second'])

    def test_missing_file_is_recorded_as_per_sample_failure(self):
        self.input.unlink()
        completed, metrics = self.run_batch_from_other_directory()
        self.assertEqual(completed.returncode, 1, completed.stderr)
        self.assertEqual(metrics['failed'], 1)
        self.assertEqual(metrics['completed'], 0)
        status = json.loads((self.args.output/self.sample/'status.json').read_text())
        self.assertEqual(status['status'], 'failed')
        self.assertIn('reads.fastq', status['reason'])


if __name__ == '__main__':
    unittest.main()
