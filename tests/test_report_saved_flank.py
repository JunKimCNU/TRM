"""Standalone saved-flank report integrity, missingness and immutable-input tests."""
import copy
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import report_saved_flank as report


def entry(sequence='CAGTACGATCGTACTGACTGATTCAGACTGA', events=2):
    anchor = sequence[:20]
    packed = int(''.join(f'{"ACGT".index(base):02b}' for base in anchor), 2)
    rank = 0x6c62272e07bb014262b821756295c58d
    for byte in anchor.encode('ascii'):
        rank = ((rank ^ byte) * 0x0000000001000000000000000000013b) % (2**128)
    raw = dict(left_forward=events-1, left_reverse=0, right_forward=1,
               right_reverse=0, ambiguous=0)
    normalized = dict(left_forward=events-1, left_reverse=1, right_forward=0,
                      right_reverse=0, ambiguous=0)
    return dict(signature=f'{packed:010x}', rank=f'{rank:032x}', sequence=anchor,
                representative_flank=sequence, events=events,
                quadrant_counts=raw, normalized_quadrant_counts=normalized)


def document(entries=None, measured=True):
    if entries is None:
        entries = [entry()]
    entries = sorted(entries, key=lambda row: (row['rank'], row['signature']))
    events = sum(row['events'] for row in entries)
    sketch = dict(capacity=128, flank_bases=20, observed_events=events,
                  distinct_lower_bound=len(entries), saturated=False, entries=entries)
    return dict(program='trm', schema_version=4, version='0.8.5',
                sample_id='<script>alert(1)</script>', analysis_kind='targeted_evidence',
                summary={'sampled_prefix': False},
                candidate_evidence=[dict(motif='AACCT', period=5, reads=events,
                    status='repeat_candidate', junction_evidence=dict(junction_reads=events,
                        validation=dict(schema_version=1, measured=measured, flank_contexts=sketch)))],
                rejected_candidates=[])


class SavedFlankReporterTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / 'native.json'
        self.output = self.root / 'report'
        self.current = document()

    def build(self, doc=None, families=None, output=None):
        self.source.write_text(json.dumps(self.current if doc is None else doc), encoding='utf-8')
        return report.build_report([self.source], output or self.output,
                                   ['TTAGG', 'TTAGGG'], families or ['AACCT'])

    def read_rows(self, name, directory=None):
        with ((directory or self.output) / name).open(encoding='utf-8') as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    def parts(self, doc=None):
        candidate = (doc or self.current)['candidate_evidence'][0]
        validation = candidate['junction_evidence']['validation']
        sketch = validation['flank_contexts']
        return candidate, validation, sketch, sketch['entries']

    def test_success_pins_escaped_html_and_original_bytes_unchanged(self):
        self.source.write_text(json.dumps(self.current), encoding='utf-8')
        original = self.source.read_bytes()
        summary = report.build_report([self.source], self.output, ['TTAGG', 'TTAGGG'], ['AACCT'])
        self.assertEqual(self.source.read_bytes(), original)
        self.assertEqual(summary['source_files'][0]['sha256'], hashlib.sha256(original).hexdigest())
        self.assertFalse(summary['changes_native_counts_or_classification'])
        self.assertFalse(summary['raw_sequence_inputs_read'])
        self.assertIsNone(summary['independent_sample_count'])
        for pin in summary['outputs']:
            data = Path(pin['path']).read_bytes()
            self.assertEqual(pin['bytes'], len(data))
            self.assertEqual(pin['sha256'], hashlib.sha256(data).hexdigest())
        page = (self.output / 'index.html').read_text()
        self.assertIn('&lt;script&gt;alert(1)&lt;/script&gt;', page)
        self.assertNotIn('<script>alert(1)</script>', page)
        queries = self.read_rows('queries.tsv')
        self.assertEqual(len(queries), 1)
        self.assertEqual(queries[0]['source_anchor_events'], '2')
        self.assertEqual(queries[0]['measurement_status'], 'evaluated_saved_representative')
        self.assertNotIn('read_id', queries[0])
        self.assertEqual(queries[0]['source_json_pointer'],
                         '/candidate_evidence/0/junction_evidence/validation/flank_contexts/entries/0')
        for key in ('native_reads', 'native_junction_reads'):
            self.assertEqual(self.read_rows('families.tsv')[0][key], '2')

    def test_missing_unmeasured_zero_sketch_and_representative_are_distinct(self):
        expected = [
            ('validation_missing', 'unmeasured', '', '', 0),
            ('measured_false', 'unmeasured', '1', '', 0),
            ('sketch_missing', 'sketch_unavailable', '', '', 0),
            ('zero', 'measured_zero_retained', '0', '0', 0),
            ('representative_missing', 'representatives_incomplete', '1', '0', 1),
        ]
        for name, state, retained, assessed, query_count in expected:
            doc = document(entries=[] if name == 'zero' else None)
            candidate, validation, sketch, entries = self.parts(doc)
            if name == 'validation_missing':
                del candidate['junction_evidence']['validation']
            elif name == 'measured_false':
                validation['measured'] = False
            elif name == 'sketch_missing':
                del validation['flank_contexts']
            elif name == 'representative_missing':
                del entries[0]['representative_flank']
            output = self.root / name
            with self.subTest(name=name):
                summary = self.build(doc, output=output)
                row = self.read_rows('families.tsv', output)[0]
                self.assertEqual(row['measurement_status'], state)
                self.assertEqual(row['retained_anchors'], retained)
                self.assertEqual(row['representatives_evaluated'], assessed)
                self.assertEqual(len(self.read_rows('queries.tsv', output)), query_count)
                self.assertIsNone(summary['full_grammar_compatible_queries'])
                self.assertIsNone(summary['mixed_only_queries'])
                if name == 'representative_missing':
                    q = self.read_rows('queries.tsv', output)[0]
                    self.assertEqual(q['measurement_status'], 'representative_unavailable')
                    self.assertEqual(q['full_grammar_compatible'], '')
                    self.assertEqual(q['sequence'], '')
                    self.assertEqual(row['full_grammar_compatible_queries'], '')
                elif name == 'zero':
                    self.assertEqual(row['full_grammar_compatible_queries'], '0')

    def test_absent_candidate_is_not_measured_zero(self):
        summary = self.build(families=['AACCT', 'AACCCT'])
        rows = self.read_rows('families.tsv')
        self.assertEqual(rows[1]['measurement_status'], 'candidate_not_reported')
        self.assertEqual(rows[1]['native_reads'], '')
        self.assertEqual(rows[1]['full_grammar_compatible_queries'], '')
        self.assertEqual(summary['query_rows'], 1)

    def test_complete_mixed_query_count_is_unweighted_by_anchor_events(self):
        self.current = document([entry('TTAGGTTAGGGTTAGGTTAGGTTAGGTTAG', events=100)])
        summary = self.build()
        self.assertEqual(summary['mixed_only_queries'], 1)
        self.assertEqual(summary['full_grammar_compatible_queries'], 1)
        row = self.read_rows('queries.tsv')[0]
        self.assertEqual(row['source_anchor_events'], '100')
        self.assertEqual(row['compatibility_class'], 'mixed_only')
        self.assertEqual(row['sequence_bases'], '30')
        self.assertEqual(row['max_contiguous_outside_bp'], '0')

    def test_duplicate_sources_duplicate_families_and_no_overwrite(self):
        self.build()
        before = {p.name: p.read_bytes() for p in self.output.iterdir()}
        with self.assertRaises(ValueError):
            report.build_report([self.source], self.output, ['TTAGG', 'TTAGGG'], ['AACCT'])
        self.assertEqual({p.name: p.read_bytes() for p in self.output.iterdir()}, before)
        for paths, families in (([self.source, self.source], ['AACCT']),
                                ([self.source], ['AACCT', 'TTAGG'])):
            with self.subTest(paths=paths, families=families), self.assertRaises(ValueError):
                report.build_report(paths, self.root/'must_not_exist', ['TTAGG', 'TTAGGG'], families)
            self.assertFalse((self.root/'must_not_exist').exists())

    def test_wrong_schema_and_malformed_contexts_fail_before_output(self):
        mutations = [
            lambda d, v, s, e: d.update(schema_version=True),
            lambda d, v, s, e: d.update(version='9.9.9'),
            lambda d, v, s, e: d.update(program='foreign'),
            lambda d, v, s, e: d['summary'].update(sampled_prefix='false'),
            lambda d, v, s, e: v.update(schema_version=2),
            lambda d, v, s, e: v.update(schema_version=True),
            lambda d, v, s, e: v.update(measured=1),
            lambda d, v, s, e: s.update(capacity=127),
            lambda d, v, s, e: s.update(flank_bases=19),
            lambda d, v, s, e: s.update(observed_events=3),
            lambda d, v, s, e: s.update(observed_events=True),
            lambda d, v, s, e: s.update(distinct_lower_bound=2),
            lambda d, v, s, e: s.update(saturated=True),
            lambda d, v, s, e: s.update(entries='wrong'),
            lambda d, v, s, e: e[0].update(events=0),
            lambda d, v, s, e: e[0].update(signature='0000000000'),
            lambda d, v, s, e: e[0].update(rank='0' * 32),
            lambda d, v, s, e: e[0].update(sequence='ACGN' * 5),
            lambda d, v, s, e: e[0].update(representative_flank='ACGT' * 20),
            lambda d, v, s, e: e[0].update(representative_flank=e[0]['sequence']+'N'),
            lambda d, v, s, e: e[0].update(representative_flank=e[0]['sequence']+'A'*61),
            lambda d, v, s, e: e[0]['quadrant_counts'].update(left_forward=2),
            lambda d, v, s, e: e[0]['normalized_quadrant_counts'].update(left_forward=0, right_reverse=1),
            lambda d, v, s, e: e[0]['normalized_quadrant_counts'].update(ambiguous=-1),
        ]
        for i, mutation in enumerate(mutations):
            doc = document()
            _, validation, sketch, entries = self.parts(doc)
            mutation(doc, validation, sketch, entries)
            output = self.root / f'invalid_{i}'
            with self.subTest(index=i), self.assertRaises(ValueError):
                self.build(doc, output=output)
            self.assertFalse(output.exists())

    def test_duplicate_anchor_and_wrong_rank_order_are_rejected(self):
        valid = document([entry(), entry('GCTAACGATCGTACTGACTGATTCAGACTGA')])
        for name in ('duplicate', 'unsorted'):
            doc = copy.deepcopy(valid)
            _, _, _, entries = self.parts(doc)
            if name == 'duplicate':
                entries[1] = copy.deepcopy(entries[0])
            else:
                entries.reverse()
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.build(doc, output=self.root / name)

    def test_source_prefix_and_rejected_container_remain_visible(self):
        self.current['summary']['sampled_prefix'] = True
        self.current['rejected_candidates'] = self.current['candidate_evidence']
        self.current['candidate_evidence'] = []
        self.build()
        row = self.read_rows('families.tsv')[0]
        self.assertEqual(row['source_sampled_prefix'], 'True')
        self.assertEqual(row['candidate_container'], 'rejected_candidates')
        self.assertTrue(self.read_rows('queries.tsv')[0]['source_json_pointer'].startswith('/rejected_candidates/'))

    def test_duplicate_json_key_nonfinite_and_unsupported_grammar_fail(self):
        for text in ('{"program":"trm","program":"trm"}', '{"value":NaN}'):
            self.source.write_text(text)
            with self.subTest(text=text), self.assertRaises(ValueError):
                report.build_report([self.source], self.output, ['TTAGG', 'TTAGGG'], ['AACCT'])
            self.assertFalse(self.output.exists())
        self.source.write_text(json.dumps(self.current))
        for units, cutoff in (([], 20), (['TTAGG', 'TTAGG'], 20), (['TTAGG'], 81)):
            with self.subTest(units=units, cutoff=cutoff), self.assertRaises(ValueError):
                report.build_report([self.source], self.output, units, ['AACCT'], cutoff)
            self.assertFalse(self.output.exists())


if __name__ == '__main__':
    unittest.main()
