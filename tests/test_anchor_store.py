"""Disk-index and resource contracts using numeric metadata, never DNA fixtures."""
import copy
import itertools
import json
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import anchor_store as S
import validate_terminal_anchors as G


def alignment(kind='P', score=2000, mapq=60):
    return {'query': 'numeric_query', 'alignment_type': kind,
            'alignment_score': score, 'mapq': mapq, 'identity': .99,
            'query_coverage': 1.0, 'aligned_query_bp': 1000,
            'future_metadata': {'nested': [17, None]}}


def record(array_id='a1', source=1, ordinal=1, name='duplicate_allowed', length=1500):
    return {'array_id': array_id, 'source_index': source, 'read_ordinal': ordinal,
            'read_id': name, 'read_length': length}


class CompactAlignmentTests(unittest.TestCase):
    def test_exhaustive_numeric_row_orders_preserve_selection_and_count(self):
        alphabet = [alignment(), alignment(mapq=0), alignment('S', 1799),
                    alignment('S', 1800), alignment('I', 2100), alignment('i', -1)]
        settings = {'min_mapq': 20, 'min_identity': .85, 'min_query_coverage': .8,
                    'min_aligned_bp': 200, 'competing_score_fraction': .9}
        checked = 0
        with tempfile.TemporaryDirectory() as temporary:
            store = S.AnchorStore(Path(temporary) / 'metadata.sqlite')
            try:
                for n in range(6):
                    for indices in itertools.product(range(len(alphabet)), repeat=n):
                        rows = [copy.deepcopy(alphabet[i]) for i in indices]
                        name = 'numeric_query_' + str(checked)
                        for row in rows:
                            row['query'] = name
                            store.add_alignment(row)
                        compact = store.get(name, [])
                        self.assertEqual(len(compact), len(rows), indices)
                        self.assertEqual(G.select_alignment(compact, settings),
                                         G.select_alignment(rows, settings), indices)
                        self.assertLessEqual(list.__len__(compact), 3)
                        checked += 1
                self.assertEqual(checked, 9331)
            finally:
                store.close()

    def test_negative_zero_and_arbitrary_integer_scores_survive_json(self):
        settings = {'min_mapq': 20, 'min_identity': .85, 'min_query_coverage': .8,
                    'min_aligned_bp': 200, 'competing_score_fraction': .9}
        with tempfile.TemporaryDirectory() as temporary:
            store = S.AnchorStore(Path(temporary) / 'metadata.sqlite')
            try:
                for i, score in enumerate([-2000, 0, 2 ** 80, -(2 ** 80)]):
                    for j, competitor in enumerate([score - 1, score, score + 1]):
                        name = f'numeric_{i}_{j}'
                        rows = [alignment('I', competitor), alignment('P', score),
                                alignment('i', competitor - 1)]
                        for row in rows:
                            row['query'] = name
                            store.add_alignment(row)
                        compact = store.get(name)
                        self.assertEqual(G.select_alignment(compact, settings),
                                         G.select_alignment(rows, settings))
                        self.assertEqual(compact[0]['alignment_score'], score)
                        self.assertEqual(compact[1]['alignment_score'], competitor)
                        self.assertEqual(compact[0]['future_metadata'], rows[1]['future_metadata'])
                self.assertIsNone(store.get('absent'))
                marker = object()
                self.assertIs(store.get('absent', marker), marker)
            finally:
                store.close()

    def test_primary_quality_precedence_is_unchanged(self):
        settings = {'min_mapq': 20, 'min_identity': .85, 'min_query_coverage': .8,
                    'min_aligned_bp': 200, 'competing_score_fraction': .9}
        mutations = [{'mapq': 255}, {'identity': .84}, {'query_coverage': .79},
                     {'aligned_query_bp': 199}, {'mapq': 0, 'identity': 0}]
        with tempfile.TemporaryDirectory() as temporary:
            store = S.AnchorStore(Path(temporary) / 'metadata.sqlite')
            try:
                for i, mutation in enumerate(mutations):
                    rows = [{**alignment(), **mutation}, alignment('S', 3000)]
                    for row in rows:
                        row['query'] = str(i)
                        store.add_alignment(row)
                    self.assertEqual(G.select_alignment(store.get(str(i)), settings),
                                     G.select_alignment(rows, settings))
            finally:
                store.close()


class ExportIndexTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.store = S.AnchorStore(Path(self.temporary.name) / 'metadata.sqlite')

    def tearDown(self):
        self.store.close()
        self.temporary.cleanup()

    def test_anchor_and_array_counts_and_missing_duplicate_guards(self):
        self.store.add_anchor('a1_L', 1000)
        self.assertEqual(self.store.lengths.get('a1_L'), 1000)
        self.assertEqual(self.store.lengths.get('missing', 7), 7)
        self.store.add_array(record())
        with self.assertRaisesRegex(ValueError, 'FASTA and exported'):
            self.store.finish_export(1, 1)
        for name in ['missing', 'a1_L']:
            if name == 'a1_L':
                self.store.mark_anchor(name)
            with self.assertRaisesRegex(ValueError, 'duplicate or missing'):
                self.store.mark_anchor(name)
        with self.assertRaisesRegex(ValueError, 'duplicate FASTA'):
            self.store.add_anchor('a1_L', 1000)
        with self.assertRaisesRegex(ValueError, 'duplicate array'):
            self.store.add_array(record())
        for arrays, anchors in [(0, 1), (2, 1), (1, 0), (1, 2)]:
            with self.assertRaisesRegex(ValueError, 'counts disagree'):
                self.store.finish_export(arrays, anchors)
        self.store.finish_export(1, 1)

    def test_empty_index_is_complete_only_for_zero_counts(self):
        self.store.finish_export(0, 0)
        self.assertEqual(self.store.record_counts(), (0, 0, 0, 0))
        with self.assertRaises(ValueError):
            self.store.finish_export(1, 0)

    def test_record_identity_is_source_and_ordinal_not_read_name(self):
        rows = [record(), record('a2'), record('a3', ordinal=2),
                record('a4', source=2), record('a5', ordinal=2 ** 80)]
        for row in reversed(rows):
            self.store.add_array(row)
        self.assertEqual(self.store.record_counts(), (4, 0, 0, 0))
        for i, change in enumerate([{'name': 'changed'}, {'length': 1499}]):
            with self.assertRaisesRegex(ValueError, 'inconsistent name/length'):
                self.store.add_array(record('bad' + str(i), **change))

    def test_record_outcomes_are_unique_ordinal_unions_not_array_sums(self):
        for row in [record(), record('a2'), record('a3', ordinal=2), record('a4', source=2)]:
            self.store.add_array(row)
        def outcome(source, ordinal, joint=False, long=False, nuclear=False):
            return {'source_index': source, 'read_ordinal': ordinal,
                    'joint_terminal_context': joint, 'joint_long_model_span_context': long,
                    'joint_nuclear_terminal_context': nuclear}
        self.store.mark_outcome(outcome(1, 1, True, False, True))
        self.store.mark_outcome(outcome(1, 1, True, True, True))
        self.store.mark_outcome(outcome(1, 1))
        self.store.mark_outcome(outcome(2, 1, True))
        self.assertEqual(self.store.record_counts(), (3, 2, 1, 1))

    def test_index_cache_is_validated_and_mmap_disabled(self):
        self.assertEqual(self.store.db.execute('PRAGMA cache_size').fetchone()[0], -8192)
        self.assertEqual(self.store.db.execute('PRAGMA mmap_size').fetchone()[0], 0)
        self.assertEqual(self.store.db.execute('PRAGMA temp_store').fetchone()[0], 1)
        for bad in [0, -1, True, 1.5, '8', None]:
            with self.subTest(value=bad), self.assertRaises(ValueError):
                S.AnchorStore(Path(self.temporary.name) / 'invalid.sqlite', cache_mib=bad)


class MemoryLimitTests(unittest.TestCase):
    def test_no_request_does_not_mutate_os_limits(self):
        for pair in [(resource.RLIM_INFINITY, resource.RLIM_INFINITY), (1000, 2000)]:
            with mock.patch.object(S.resource, 'getrlimit', return_value=pair), \
                    mock.patch.object(S.resource, 'setrlimit') as setter:
                result = S.apply_memory_limit(None)
                setter.assert_not_called()
                self.assertIsNone(result['requested_memory_limit_mib'])
                self.assertEqual(result['effective_address_space_limit_bytes'],
                                 None if pair[0] == resource.RLIM_INFINITY else pair[0])

    def test_requested_limit_preserves_stricter_limits_and_hard_limit(self):
        mib = 1024 ** 2
        for soft, hard, expected in [(-1, -1, 256 * mib), (128 * mib, -1, 128 * mib),
                                     (512 * mib, 512 * mib, 256 * mib),
                                     (128 * mib, 192 * mib, 128 * mib),
                                     (-1, 128 * mib, 128 * mib)]:
            with mock.patch.object(S.resource, 'getrlimit', return_value=(soft, hard)), \
                    mock.patch.object(S.resource, 'setrlimit') as setter:
                result = S.apply_memory_limit(256)
                setter.assert_called_once_with(resource.RLIMIT_AS, (expected, hard))
                self.assertTrue(result['hard_limit_unchanged'])
                self.assertEqual(result['effective_address_space_limit_bytes'], expected)

    def test_invalid_limits_fail_before_changing_resource_state(self):
        for bad in [0, -1, True, 1.5, '256', float('inf')]:
            with mock.patch.object(S.resource, 'getrlimit') as getter, \
                    mock.patch.object(S.resource, 'setrlimit') as setter:
                with self.subTest(value=bad), self.assertRaises(ValueError):
                    S.apply_memory_limit(bad)
                getter.assert_not_called()
                setter.assert_not_called()

    def test_real_cap_in_child_blocks_numeric_allocation_without_touching_parent(self):
        before = resource.getrlimit(resource.RLIMIT_AS)
        program = '''
import json, resource
from anchor_store import apply_memory_limit
before = resource.getrlimit(resource.RLIMIT_AS)
result = apply_memory_limit(128)
try:
    memory = bytearray(256 * 1024 * 1024)
except MemoryError:
    result['allocation_failed'] = True
else:
    result['allocation_failed'] = False
result['hard_unchanged_actual'] = resource.getrlimit(resource.RLIMIT_AS)[1] == before[1]
print(json.dumps(result))
'''
        completed = subprocess.run([sys.executable, '-c', program], cwd=ROOT / 'scripts',
                                   capture_output=True, text=True, timeout=20)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        result = json.loads(completed.stdout)
        self.assertTrue(result['allocation_failed'])
        self.assertTrue(result['hard_unchanged_actual'])
        self.assertLessEqual(result['effective_address_space_limit_bytes'], 128 * 1024 ** 2)
        self.assertEqual(resource.getrlimit(resource.RLIMIT_AS), before)


if __name__ == '__main__':
    unittest.main()
