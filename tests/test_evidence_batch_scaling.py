"""Evidence orchestration scaling, provenance, cache scope, and failure isolation."""
from __future__ import annotations
import contextlib
import fcntl
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/evidence_batch.py'
SPEC = importlib.util.spec_from_file_location('evidence_scaling_under_test', SCRIPT)
batch = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(batch)

STUB = '''#!/usr/bin/env python3
import json,pathlib,sys,time
args=sys.argv[1:]
paths=[]
for value in args[1:]:
 if value.startswith('--'):break
 paths.append(pathlib.Path(value))
root=paths[0].parent
with (root/'launches.jsonl').open('a') as f:f.write(json.dumps(args)+'\\n')
if paths[0].name.startswith('slow'):time.sleep(0.15)
mode=(root/'mode').read_text().strip() if (root/'mode').exists() else ''
if mode=='fail':print('fixture native failure',file=sys.stderr);sys.exit(7)
if mode=='silent':sys.exit(0)
prefix=args[args.index('--output')+1]
targets=json.loads(pathlib.Path(args[args.index('--targets')+1]).read_text())['motifs']
limit=int(args[args.index('--max-records')+1]) if '--max-records' in args else 0
count=min(4,limit) if limit else 4
inputs=[{'path':str(p.resolve()),'bytes':p.stat().st_size,'modified_unix_seconds':int(p.stat().st_mtime)} for p in paths]
result={'program':'trm','version':'0.8.4','analysis_kind':'targeted_evidence',
'summary':{'records':count,'bases':26*count,'sampled_prefix':bool(limit)},'inputs':inputs,
'candidate_evidence':[{'motif':m} for m in targets],'target_motifs':targets,
'paired_end':{'measured':'--paired-layout' in args,'layout':'fr' if '--paired-layout' in args else None},
'junction_options':{'flank_period_annotation':'--flank-period-annotation' in args},'received_options':args}
if mode=='unexpected_prefix':result['summary']['sampled_prefix']=True
if mode=='wrong_inputs':result['inputs'][0]['bytes']+=1
if mode=='wrong_targets':result['candidate_evidence']=[]
if mode=='omit_paired':result.pop('paired_end')
if mode=='omit_flank':result.pop('junction_options')
pathlib.Path(prefix+'.json').write_text(json.dumps(result))
'''


class EvidenceScalingTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix='evidence-scaling-test-')
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.manifest = self.root / 'samples.tsv'
        self.source = self.root / 'source'
        self.output = self.root / 'output'
        self.catalog = self.root / 'catalog.tsv'
        self.catalog.write_text('motif\nTTAGGG\n')
        self.binary = self.root / 'scanner'
        self.binary.write_text(STUB)
        self.binary.chmod(0o755)

    def samples(self, names, paired=False):
        rows = ['sample_id\tpath1\tpath2\n']
        for index, name in enumerate(names):
            one = self.root / f'{name}.fq'
            one.write_text('A' * (index + 1))
            two = self.root / f'{name}_mate.fq'
            if paired:
                two.write_text('T' * (index + 1))
            rows.append(f'{name}\t{one.name}\t{two.name if paired else ""}\n')
            report = self.source / name / 'result.json'
            report.parent.mkdir(parents=True, exist_ok=True)
            batch.write_json(report, {'summary': {'records': 4, 'bases': 104}, 'candidates': [], 'rejected_candidates': []})
        self.manifest.write_text(''.join(rows))

    def run_batch(self, *extra, expected=0, jobs=2, timeout=30):
        command = [sys.executable, str(SCRIPT), '--manifest', str(self.manifest), '--source', str(self.source),
                   '--output', str(self.output), '--catalog', str(self.catalog), '--binary', str(self.binary),
                   '--jobs', str(jobs), '--threads', '1', '--cpu-budget', str(max(2, jobs)), *extra]
        result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return result

    def report(self, name='metrics.json', shard=None):
        out = self.output if shard is None else self.output / f'shard-{shard[0]:05d}-of-{shard[1]:05d}'
        return json.loads((out / name).read_text())

    def launches(self):
        path = self.root / 'launches.jsonl'
        return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []

    def test_bounded_futures_coalesced_checkpoints_and_manifest_order(self):
        names = [f'sample{i:03}' for i in range(50)]
        self.samples(names)
        self.run_batch('--checkpoint-seconds', '1e100', jobs=4)
        report = self.report('summary.json')
        self.assertEqual((report['max_inflight'], report['checkpoint_count']), (4, 2))
        self.assertEqual([row['sample_id'] for row in report['samples']], names)
        self.assertTrue(report['all_succeeded'])
        self.assertEqual(report['scanned_samples'], 50)
        progress = self.report('progress.json')
        self.assertEqual((progress['inflight'], progress['running_count'], progress['pending']), (0, 0, 0))
        self.assertEqual(progress['submitted'], 50)
        self.assertEqual(len(self.launches()), 50)

    def test_checkpoints_continue_during_a_long_task(self):
        self.samples(['slow_sample'])
        self.run_batch('--checkpoint-seconds', '0.025', jobs=1)
        self.assertGreater(self.report()['checkpoint_count'], 2)
        self.assertEqual(self.report()['completed'], 1)

    def test_largest_first_default_and_manifest_order_resume_do_not_change_identity(self):
        self.samples(['small', 'medium', 'large'])
        self.run_batch(jobs=1)
        self.assertEqual([Path(args[1]).name for args in self.launches()], ['large.fq', 'medium.fq', 'small.fq'])
        self.run_batch('--resume', '--schedule', 'manifest', jobs=1)
        self.assertEqual(len(self.launches()), 3)
        report = self.report()
        self.assertEqual((report['scanned_samples'], report['resumed_samples']), (0, 3))
        self.assertEqual(report['records'], 12)
        self.assertEqual(report['current_execution']['scanned']['bases'], 0)
        self.assertEqual(report['current_execution']['attempted']['user_cpu_seconds'], 0)
        self.assertEqual(report['current_execution']['resumed_artifacts']['bases'], 312)
        self.assertIsNone(report['compressed_MB_per_second'])

    def test_legacy_identity_and_status_resume_without_native_launch(self):
        self.samples(['sample'])
        self.run_batch()
        path = self.output / 'sample/status.json'
        old = json.loads(path.read_text())
        for key in ('execution_id', 'execution_version', 'subprocess_started'):
            old.pop(key, None)
        batch.write_json(path, old)
        self.assertEqual(set(old['identity']), {'source_sha256', 'inputs', 'binary_sha256', 'targets', 'threads', 'max_records'})
        self.run_batch('--resume', '--flank-period-annotation', 'off')
        self.assertEqual(len(self.launches()), 1)
        self.assertEqual(json.loads(path.read_text()), old)
        self.assertEqual(self.report()['resumed_samples'], 1)

    def test_new_native_options_are_forwarded_and_each_change_invalidates_resume(self):
        self.samples(['sample'], paired=True)
        self.run_batch()
        for extra, launches in [(['--memory-limit-mib', '256'], 2),
                                (['--memory-limit-mib', '512'], 3),
                                (['--memory-limit-mib', '512', '--paired-layout', 'fr'], 4),
                                (['--memory-limit-mib', '512', '--paired-layout', 'fr', '--flank-period-annotation', 'on'], 5)]:
            self.run_batch('--resume', *extra)
            self.assertEqual(len(self.launches()), launches)
            self.run_batch('--resume', *extra)
            self.assertEqual(len(self.launches()), launches)
        state = self.report('sample/status.json')
        options = state['identity']['native_options']
        self.assertEqual(options, ['--memory-limit-mib', '512', '--paired-layout', 'fr', '--flank-period-annotation', 'on'])
        self.assertEqual(state['command'][-len(options):], options)
        self.assertEqual(self.report()['configured_jobs_times_memory_limit_mib'], 1024)

    def test_binary_source_input_target_and_corrupt_result_changes_recompute(self):
        self.samples(['sample'])
        self.run_batch()
        self.binary.write_text(STUB + '\n# a changed executable snapshot\n')
        self.run_batch('--resume')
        self.assertTrue(list((self.output / '.binaries').glob('*/trm-evidence')))
        path = self.source / 'sample/result.json'
        source = json.loads(path.read_text()); source['note'] = 'changed'; batch.write_json(path, source)
        self.run_batch('--resume')
        with (self.root / 'sample.fq').open('a') as stream:
            stream.write('C')
        self.run_batch('--resume')
        self.catalog.write_text('motif\nTTAGGG\nTTAGGC\n')
        self.run_batch('--resume')
        (self.output / 'sample/evidence.json').write_text('{}')
        self.run_batch('--resume')
        (self.output / 'sample/status.json').write_text('{broken')
        self.run_batch('--resume')
        self.assertEqual(len(self.launches()), 7)
        self.assertEqual(self.report()['scanned_samples'], 1)

    def test_shards_cover_samples_once_and_have_independent_artifacts(self):
        names = [f'sample{i:02}' for i in range(12)]
        self.samples(names)
        collected = []
        for shard in range(4):
            self.run_batch('--shard-index', str(shard), '--shard-count', '4')
            report = self.report('summary.json', (shard, 4))
            wanted = [name for name in names if batch.sample_shard(name, 4) == shard]
            self.assertEqual([row['sample_id'] for row in report['samples']], wanted)
            collected += wanted
        self.assertEqual(sorted(collected), sorted(names))
        self.assertEqual(len(set(collected)), len(names))
        self.assertFalse((self.output / 'metrics.json').exists())
        self.assertEqual(len(self.launches()), len(names))

    def test_empty_shard_skips_unselected_missing_input_but_bad_manifest_still_fails(self):
        self.samples(['only'])
        (self.root / 'only.fq').unlink()
        empty = 1 - batch.sample_shard('only', 2)
        self.run_batch('--shard-index', str(empty), '--shard-count', '2')
        report = self.report(shard=(empty, 2))
        self.assertEqual(report['sample_count'], 0); self.assertTrue(report['all_succeeded'])
        self.assertFalse(self.launches())
        with self.manifest.open('a') as stream:
            stream.write('only\tabsent\t\n')
        self.run_batch('--shard-index', str(empty), '--shard-count', '2', expected=2)

    def test_source_failed_is_nonzero_and_quarantines_old_published_evidence(self):
        self.samples(['sample']); self.run_batch()
        saved = (self.output / 'sample/evidence.json').read_bytes()
        (self.source / 'sample/result.json').unlink()
        self.run_batch('--resume', expected=1)
        self.assertFalse(self.report()['all_succeeded']); self.assertEqual(self.report()['source_failed'], 1)
        self.assertFalse((self.output / 'sample/evidence.json').exists())
        self.assertTrue(any(path.read_bytes() == saved for path in (self.output / 'sample/attempts').glob('*/previous_published/evidence.json')))
        self.assertEqual(len(self.launches()), 1)

    def test_success_without_new_output_and_native_failure_preserve_attempts(self):
        self.samples(['sample']); self.run_batch()
        (self.root / 'mode').write_text('silent')
        self.run_batch(expected=1)
        self.assertFalse((self.output / 'sample/evidence.json').exists())
        self.assertIn('evidence.json', self.report('sample/status.json')['reason'])
        (self.root / 'mode').write_text('fail')
        self.run_batch(expected=1)
        attempts = list((self.output / 'sample/attempts').glob('attempt-*'))
        self.assertEqual(len(attempts), 3)
        self.assertTrue(any('fixture native failure' in (path / 'stderr.log').read_text() for path in attempts))
        self.assertEqual(self.report()['current_execution']['attempted']['samples'], 1)

    def test_partial_source_and_supplied_input_mismatch_never_become_complete(self):
        self.samples(['sample'])
        source = self.source / 'sample/result.json'
        value = json.loads(source.read_text()); value['summary']['sampled_prefix'] = True; batch.write_json(source, value)
        self.run_batch(expected=1)
        self.assertEqual(self.report()['source_failed'], 1)
        value['summary']['sampled_prefix'] = False
        value['inputs'] = [{'path': str(self.root / 'sample.fq'), 'bytes': 999}]; batch.write_json(source, value)
        self.run_batch(expected=1)
        self.assertIn('input bytes differs', self.report('sample/status.json')['reason'])
        self.assertFalse(self.launches())

    def test_corrupt_native_provenance_and_unrequested_prefix_fail(self):
        self.samples(['sample'])
        for mode in ('wrong_targets', 'wrong_inputs', 'unexpected_prefix'):
            (self.root / 'mode').write_text(mode)
            self.run_batch(expected=1)
            self.assertFalse(self.report()['all_succeeded'])
            self.assertFalse((self.output / 'sample/evidence.json').exists())

    def test_requested_pair_or_flank_requires_positive_native_provenance(self):
        self.samples(['sample'], paired=True)
        for mode, option in [('omit_paired', ('--paired-layout', 'fr')), ('omit_flank', ('--flank-period-annotation', 'on'))]:
            (self.root / 'mode').write_text(mode)
            self.run_batch(*option, expected=1)
            self.assertIn('provenance is missing', self.report('sample/status.json')['reason'])
            self.assertFalse((self.output / 'sample/evidence.json').exists())

    def test_prefix_bytes_unknown_and_pair_prefix_validation_precedes_output_creation(self):
        self.samples(['sample'], paired=True)
        self.run_batch('--paired-layout', 'fr', '--max-records', '1', expected=2)
        self.assertFalse(self.output.exists())
        self.run_batch('--paired-layout', 'fr', '--max-records', '2')
        report = self.report(); state = self.report('sample/status.json')
        self.assertIsNone(report['compressed_MB_per_second']); self.assertEqual(report['scanned_input_bytes'], 0)
        self.assertEqual(state['input_byte_scope'], 'prefix_processed_bytes_unknown')
        self.assertGreater(state['input_artifact_bytes_total'], 0)

    def test_invalid_cli_precedes_output_side_effects(self):
        self.samples(['sample'])
        for extra in [('--checkpoint-seconds', 'nan'), ('--checkpoint-seconds', '0'), ('--shard-index', '0'),
                      ('--shard-index', '2', '--shard-count', '2'), ('--max-records', '-1'),
                      ('--cpu-budget', '1'), ('--memory-limit-mib', '0'), ('--paired-layout', 'fr'),
                      ('--output', str(self.source / 'nested'))]:
            with self.subTest(extra=extra):
                self.run_batch(*extra, expected=2)
                self.assertFalse(self.output.exists())
        self.assertFalse(self.launches())

    def test_sample_symlink_is_rejected_without_writes_to_source(self):
        self.samples(['sample']); self.output.mkdir()
        (self.output / 'sample').symlink_to(self.source / 'sample', target_is_directory=True)
        before = {p.name: p.read_bytes() for p in (self.source / 'sample').iterdir()}
        self.run_batch(expected=1)
        self.assertEqual({p.name: p.read_bytes() for p in (self.source / 'sample').iterdir()}, before)
        receipt = self.report('.failures/sample.json')
        self.assertIn('symlink', receipt['reason'])
        self.assertFalse(self.launches())

    def test_unselected_manifest_input_is_protected_from_output_overwrite(self):
        self.samples(['selected', 'unselected'])
        self.manifest.write_text('sample_id\tpath1\nselected\tselected.fq\nunselected\toutput/metrics.json\n')
        self.run_batch('--samples', 'selected', expected=2)
        self.assertFalse(self.output.exists())
        self.assertFalse(self.launches())

    def test_output_lock_prevents_same_namespace_concurrent_execution(self):
        self.samples(['sample']); self.output.mkdir()
        with (self.output / 'runner.lock').open('a+') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            self.run_batch(expected=2)
        self.assertFalse(self.launches())

    def test_atomic_json_ignores_preexisting_part_links_and_hash_streams(self):
        victim = self.root / 'victim'; victim.write_text('preserve')
        path = self.root / 'data.json'
        path.with_name('data.json.part').symlink_to(victim)
        batch.write_json(path, {'ok': True})
        self.assertEqual(victim.read_text(), 'preserve')
        with patch.object(Path, 'read_bytes', side_effect=AssertionError('whole-file hash allocation')):
            self.assertEqual(batch.sha256(victim), hashlib.sha256(b'preserve').hexdigest())

    def test_worker_exception_is_persisted_and_other_samples_finish(self):
        self.samples(['bad', 'good'])
        args = SimpleNamespace(output=self.output, jobs=2, threads=1, execution_id='current', schedule='manifest',
             cpu_budget=2, max_records=0, memory_limit_mib=None, paired_layout=None, flank_period_annotation='off',
             shard_index=None, shard_count=None, checkpoint_seconds=1e100)
        self.output.mkdir()
        manifest = [{'sample_id': name, 'path1': str(self.root / f'{name}.fq')} for name in ('bad', 'good')]
        def execute(item, *_):
            if item['sample_id'] == 'bad':
                raise RuntimeError('injected coordinator-visible failure')
            return {'sample_id': 'good', 'status': 'complete', 'execution_id': 'current', 'records': 4, 'bases': 104}
        with patch.object(batch, 'run_one', side_effect=execute), contextlib.redirect_stdout(io.StringIO()):
            report = batch.run_scheduled(manifest, args, set(), self.binary)
        self.assertFalse(report['all_succeeded']); self.assertEqual((report['completed'], report['failed']), (1, 1))
        self.assertIn('injected', self.report('bad/status.json')['reason'])


if __name__ == '__main__':
    unittest.main()
