"""Numeric metadata/geometry and existing real-file checks; no synthetic DNA reads."""
import copy
import contextlib
import csv
import html
import io
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

from archive_test_support import archive_only
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
spec = importlib.util.spec_from_file_location('irregular_anchor_adapter', ROOT / 'scripts/validate_irregular_anchors.py')
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)


def array(side='left', strand='+', length=500):
    start = 1000 if side == 'left' else 0
    n = 1000 + length
    row = {'array_id': 'source1_read1_array1', 'source_index': 1, 'read_ordinal': 1,
           'read_id': 'same_name', 'read_length': n, 'group_id': 'published_units',
           'start0': start, 'end0': start + length, 'span_bp': length,
           'strand': strand, 'left_end_distance': start, 'right_end_distance': n - start - length,
           'model': {'matches': length - 20, 'substitutions': 10, 'insertions': 10,
                     'deletions': 2, 'model_bases': length - 8, 'exact_unit_bases': 120,
                     'identity': .9, 'optimal_parse_ambiguous': False,
                     'trace': [{'numeric_test_step': 1}],
                     'future_detector_detail': {'retain_this': True}}}
    for name, suffix in [('left', 'L'), ('right', 'R')]:
        point = row['start0'] if name == 'left' else row['end0']
        exported = name == side
        a = 0 if name == 'left' and exported else point
        b = n if name == 'right' and exported else point
        row[name + '_anchor'] = {'id': row['array_id'] + '_' + suffix if exported else '',
                                'start0': a, 'end0': b, 'length': b - a,
                                'status': 'exported' if exported else 'below_minimum',
                                'stop_reason': 'maximum_length' if exported else 'read_boundary'}
    return row


def alignment(name, strand='+', target_start=98500, qstart=0, qend=1000, mapq=60):
    span = qend - qstart
    return {'query': name, 'query_length': 1000, 'query_start': qstart, 'query_end': qend,
            'strand': strand, 'target': 'chr1', 'target_length': 100000,
            'target_start': target_start, 'target_end': target_start + span,
            'matches': span, 'block_length': span, 'mapq': mapq, 'alignment_score': span * 2,
            'alignment_type': 'P', 'identity': 1.0, 'aligned_query_bp': span,
            'query_coverage': span / 1000}


class ArrayValidationTests(unittest.TestCase):
    def validate(self, row):
        lengths = {row[side + '_anchor']['id']: row[side + '_anchor']['length']
                   for side in ['left', 'right'] if row[side + '_anchor']['id']}
        return M.validate_array(row, [{'records': 5, 'bases': 10000}], {'published_units'}, lengths)

    def test_indel_accounting_and_unknown_model_details_preserved(self):
        row = array()
        before = copy.deepcopy(row)
        self.assertEqual(len(self.validate(row)), 1)
        self.assertEqual(row, before)

    def test_source_indices_are_one_based_and_checked(self):
        for value in [0, 2, True, '1']:
            row = array(); row['source_index'] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.validate(row)

    def test_read_ordinal_cannot_exceed_source(self):
        row = array(); row['read_ordinal'] = 6
        with self.assertRaises(ValueError):
            self.validate(row)

    def test_edit_consumption_and_exact_count_rejections(self):
        for key, value in [('insertions', 11), ('model_bases', 999), ('exact_unit_bases', 999), ('matches', True)]:
            row = array(); row['model'][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.validate(row)

    def test_nonfinite_identity_is_not_a_valid_model_observation(self):
        for value in [float('nan'), float('inf'), -0.1, 1.1, True]:
            row = array(); row['model']['identity'] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.validate(row)

    def test_missing_edit_trace_is_not_completed_model_evidence(self):
        row = array(); row['model']['trace'] = []
        with self.assertRaises(ValueError):
            self.validate(row)

    def test_anchor_is_directly_adjacent_and_uniquely_named(self):
        for change in ['coordinate', 'identifier', 'status']:
            row = array()
            if change == 'coordinate':
                row['left_anchor']['end0'] -= 1
                row['left_anchor']['length'] -= 1
            elif change == 'identifier':
                row['left_anchor']['id'] = 'unrelated'
            else:
                row['left_anchor']['status'] = 'not_exported'
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.validate(row)

    def test_ambiguous_parse_is_not_ambiguous_orientation(self):
        row = array(); row['model']['optimal_parse_ambiguous'] = True
        self.validate(row)
        name = row['left_anchor']['id']
        result = M.evaluate_array(row, {name: [alignment(name)]}, M.DEFAULTS)
        self.assertTrue(result['joint_terminal_context'])
        self.assertTrue(result['model']['optimal_parse_ambiguous'])

    def test_duplicate_json_keys_fail(self):
        with self.assertRaises(ValueError):
            json.loads('{"span_bp":1,"span_bp":2}', object_pairs_hook=M.unique_object)

    def test_changed_metadata_file_fails_hash_verification(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'numeric_receipt.json'
            path.write_text('{"records":3}')
            record = {'size_bytes': path.stat().st_size, 'sha256': M.digest(path)}
            M.verify_file(path, record, 'size_bytes')
            path.write_text('{"records":4}')
            with self.assertRaises(ValueError):
                M.verify_file(path, record, 'size_bytes')


class IrregularGeometryTests(unittest.TestCase):
    def test_paf_order_preserves_primary_and_all_competing_rows(self):
        row = array(); name = row['left_anchor']['id']
        primary = alignment(name)
        secondary = {**primary, 'alignment_type': 'S', 'alignment_score': 1799,
                     'mapq': 0, 'target_start': 500, 'target_end': 1500}
        first = M.evaluate_array(row, {name: [primary, secondary]}, M.DEFAULTS)
        second = M.evaluate_array(row, {name: [secondary, primary]}, M.DEFAULTS)
        self.assertEqual(first, second)
        self.assertTrue(first['joint_terminal_context'])
        self.assertEqual(first['anchor_alignments'][0]['reported_alignments'], 2)
        # An inversion/low-MAPQ competitor still counts exactly at the AS bound.
        competitor = {**secondary, 'alignment_type': 'I', 'alignment_score': 1800}
        result = M.evaluate_array(row, {name: [competitor, primary]}, M.DEFAULTS)
        self.assertEqual(result['anchor_alignments'][0]['mapping_status'],
                         'competitive_or_split_alignment')
        self.assertFalse(result['joint_terminal_context'])

    def test_duplicate_primaries_are_not_deduplicated_or_quality_filtered(self):
        row = array(); name = row['left_anchor']['id']
        primary = alignment(name)
        for other in (copy.deepcopy(primary), {**primary, 'mapq': 0}):
            result = M.evaluate_array(row, {name: [primary, other]}, M.DEFAULTS)
            self.assertEqual(result['anchor_alignments'][0]['mapping_status'],
                             'missing_or_multiple_primary_alignments')
            self.assertEqual(result['anchor_alignments'][0]['reported_alignments'], 2)
            self.assertEqual(result['qualified_anchor_count'], 0)
            self.assertFalse(result['joint_terminal_context'])

    def test_unknown_array_metadata_and_trace_content_survive_evaluation(self):
        row = array()
        row['future_array_context'] = {'numeric_items': [17, {'nested': None}]}
        row['model']['trace'] = [{'start0': 1000, 'end0': 1500,
                                 'future_numeric_annotation': [1, 2, 3]}]
        before = copy.deepcopy(row)
        name = row['left_anchor']['id']
        result = M.evaluate_array(row, {name: [alignment(name)]}, M.DEFAULTS)
        self.assertEqual(row, before)
        for key, value in before.items():
            self.assertEqual(result[key], value)

    def test_mapping_threshold_boundary_retains_original_geometry(self):
        row = array(); name = row['left_anchor']['id']
        exact_threshold = alignment(name, target_start=98700, qstart=200,
                                    qend=1000, mapq=20)
        exact_threshold.update(matches=680, identity=.85)
        result = M.evaluate_array(row, {name: [exact_threshold]}, M.DEFAULTS)
        anchor = result['anchor_alignments'][0]
        self.assertEqual(anchor['mapping_status'], 'qualified_within_reported_alignments')
        self.assertTrue(result['joint_terminal_context'])
        self.assertEqual(anchor['mapped_junction_boundary0'], 99500)
        self.assertEqual(anchor['unaligned_anchor_junction_bp'], 0)
        self.assertEqual(anchor['query_coverage'], .8)

    def test_four_side_mapping_strand_combinations(self):
        for side in ['left', 'right']:
            for mapping_strand in ['+', '-']:
                # Expected read polarity is G for a left flank and C for right.
                row = array(side=side, strand='+' if side == 'left' else '-')
                outward_right = (side == 'left') == (mapping_strand == '+')
                start = 98500 if outward_right else 500
                name = row[side + '_anchor']['id']
                result = M.evaluate_array(row, {name: [alignment(name, mapping_strand, start)]}, M.DEFAULTS)
                with self.subTest(side=side, strand=mapping_strand):
                    self.assertTrue(result['joint_long_model_span_context'])
                    self.assertEqual(result['anchor_alignments'][0]['polarity'], 'expected')

    def test_modeled_500bp_is_not_exact_500bp(self):
        row = array(length=500)
        name = row['left_anchor']['id']
        result = M.evaluate_array(row, {name: [alignment(name)]}, M.DEFAULTS)
        self.assertTrue(result['joint_long_model_span_context'])
        self.assertEqual(result['model']['exact_unit_bases'], 120)
        self.assertNotIn('joint_long_terminal_context', result)
        self.assertNotIn('tract_ge_long_threshold', result['anchor_alignments'][0])
        self.assertFalse(result['physical_chromosome_end_confirmed'])
        self.assertTrue(result['model']['future_detector_detail']['retain_this'])

    def test_499bp_remains_short_model_context(self):
        row = array(length=499); name = row['left_anchor']['id']
        result = M.evaluate_array(row, {name: [alignment(name)]}, M.DEFAULTS)
        self.assertTrue(result['joint_terminal_context'])
        self.assertFalse(result['joint_long_model_span_context'])

    def test_ambiguous_direction_never_joint(self):
        row = array(strand='ambiguous'); name = row['left_anchor']['id']
        result = M.evaluate_array(row, {name: [alignment(name)]}, M.DEFAULTS)
        self.assertEqual(result['localization_status'], 'ambiguous_array_direction')
        self.assertFalse(result['joint_long_model_span_context'])

    def test_no_anchor_is_distinct_from_no_array(self):
        row = array()
        row['left_anchor'].update(id='', status='below_minimum', start0=1000, end0=1000, length=0)
        result = M.evaluate_array(row, {}, M.DEFAULTS)
        self.assertEqual(result['span_bp'], 500)
        self.assertEqual(result['localization_status'], 'no_exportable_adjacent_anchor')

    def test_clipping_does_not_extrapolate_boundary(self):
        row = array(); name = row['left_anchor']['id']
        result = M.evaluate_array(row, {name: [alignment(name, qend=979)]}, M.DEFAULTS)
        anchor = result['anchor_alignments'][0]
        self.assertEqual(anchor['unaligned_anchor_junction_bp'], 21)
        self.assertFalse(result['joint_terminal_context'])
        self.assertEqual(anchor['mapped_junction_boundary0'], 99479)

    def test_competitor_and_organelle_do_not_receive_terminal_credit(self):
        row = array(); name = row['left_anchor']['id']
        primary = alignment(name)
        secondary = {**primary, 'alignment_type': 'S', 'alignment_score': 1900, 'target_start': 100}
        result = M.evaluate_array(row, {name: [primary, secondary]}, M.DEFAULTS)
        self.assertEqual(result['localization_status'], 'anchor_mapping_unresolved')
        self.assertFalse(result['joint_terminal_context'])
        result = M.evaluate_array(row, {name: [primary]}, M.DEFAULTS, {'chr1': 'organelle'})
        self.assertFalse(result['joint_terminal_context'])

    def test_two_sided_context_and_disjoint_spans_are_not_merged(self):
        row = array()
        row['read_length'] = 2500
        row['right_end_distance'] = 1000
        row['right_anchor'] = {'id': row['array_id'] + '_R', 'start0': 1500, 'end0': 2500,
                               'length': 1000, 'status': 'exported', 'stop_reason': 'read_boundary'}
        mappings = {row[side + '_anchor']['id']: [alignment(row[side + '_anchor']['id'])] for side in ['left', 'right']}
        result = M.evaluate_array(row, mappings, M.DEFAULTS)
        self.assertEqual(result['localization_status'], 'two_sided_anchor_context')
        self.assertFalse(result['joint_long_model_span_context'])
        self.assertEqual(result['span_bp'], 500)


@archive_only
class ExistingActualSourceTests(unittest.TestCase):
    def test_existing_actual_export_fasta_and_paf_lengths_remain_compatible(self):
        folder = ROOT / 'results/matched_followup_20260925/dsm3286_anchors'
        receipt_path = folder / 'alignment_receipt.json'
        if not receipt_path.exists():
            self.skipTest('existing actual matched anchor artifacts unavailable')
        receipt = json.loads(receipt_path.read_text())
        queries = Path(receipt['queries']['path'])
        paf = Path(receipt['paf']['path'])
        lengths = M.fasta_lengths(queries)
        self.assertGreater(len(lengths), 0)
        checked = 0
        with paf.open() as handle:
            for line in handle:
                if not line.strip():
                    continue
                row = M.parse_paf(line)
                self.assertEqual(row['query_length'], lengths[row['query']])
                self.assertLessEqual(row['query_end'], row['query_length'])
                checked += 1
                if checked == 20:
                    break
        self.assertEqual(checked, 20)


class StreamingOutputTests(unittest.TestCase):
    def test_streamed_numeric_rows_preserve_details_late_columns_and_html(self):
        first = M.evaluate_array(array(), {}, M.DEFAULTS)
        second_input = array()
        second_input['array_id'] = 'second_array'
        second_input['read_id'] = 'escaped<&>한글'
        second_input['new_late_field'] = {'numeric': [17, 29], 'empty': None}
        name = second_input['left_anchor']['id']
        second = M.evaluate_array(second_input, {name: [alignment(name)]}, M.DEFAULTS)
        rows = [first, second]
        flat = [{k: v for k, v in row.items() if k != 'anchor_alignments'} for row in rows]
        anchors = [{'array_id': row['array_id'], 'source_index': row['source_index'],
                    'read_ordinal': row['read_ordinal'], **a}
                   for row in rows for a in row['anchor_alignments']]
        array_fields = list(dict.fromkeys(k for row in flat for k in row))
        anchor_fields = list(dict.fromkeys(k for row in anchors for k in row))
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            spool = directory / 'spool.jsonl'
            spool.write_text(''.join(json.dumps(row, ensure_ascii=False) + '\n' for row in rows))
            M.write_outputs(spool, directory, array_fields, anchor_fields)
            self.assertEqual((directory / 'details.json').read_text(),
                             json.dumps(rows, ensure_ascii=False, indent=2) + '\n')
            for filename, fields, records in [('arrays.tsv', array_fields, flat),
                                               ('anchors.tsv', anchor_fields, anchors)]:
                expected = directory / ('expected_' + filename)
                with expected.open('w', newline='') as handle:
                    writer = csv.DictWriter(handle, fields, delimiter='\t')
                    writer.writeheader()
                    for row in records:
                        writer.writerow({k: json.dumps(v, ensure_ascii=False)
                                         if isinstance(v, (dict, list)) else v for k, v in row.items()})
                self.assertEqual((directory / filename).read_bytes(), expected.read_bytes())
            with (directory / 'arrays.tsv').open() as handle:
                saved = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(saved[0]['new_late_field'], '')
            self.assertEqual(json.loads(saved[1]['new_late_field']), second_input['new_late_field'])
            expected_html_rows = ''.join('<tr>' + ''.join('<td>' + html.escape(str(row[key])) + '</td>'
                for key in ['read_id', 'span_bp', 'strand', 'localization_status',
                            'joint_long_model_span_context']) + '</tr>' for row in rows)
            self.assertEqual((directory / 'index.html').read_text(),
                             M.HTML_HEADER + expected_html_rows + '</table></html>')

    def test_empty_spool_emits_empty_outputs_not_a_biological_summary(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            spool = directory / 'spool.jsonl'
            spool.write_text('')
            M.write_outputs(spool, directory, [], [])
            self.assertEqual((directory / 'details.json').read_text(), '[]\n')
            self.assertEqual((directory / 'arrays.tsv').stat().st_size, 0)
            self.assertEqual((directory / 'anchors.tsv').stat().st_size, 0)
            self.assertFalse((directory / 'summary.json').exists())

    def test_broken_numeric_spool_does_not_silently_drop_remaining_rows(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            spool = directory / 'spool.jsonl'
            spool.write_text('{"numeric_metadata":')
            with self.assertRaises(json.JSONDecodeError):
                M.write_outputs(spool, directory, [], [])
            self.assertFalse((directory / 'summary.json').exists())


class FailurePublicationTests(unittest.TestCase):
    def test_cli_reports_numeric_allocation_or_io_failure_without_complete_result(self):
        for exception in [MemoryError('numeric allocation failure'), OSError('numeric disk failure'),
                          ValueError('numeric metadata failure')]:
            with self.subTest(error=type(exception).__name__), tempfile.TemporaryDirectory() as temporary:
                output = Path(temporary) / 'failed'
                argv = ['validate_irregular_anchors.py', '--export-dir', temporary,
                        '--assembly', str(Path(temporary) / 'unused'), '--output-dir', str(output)]
                stderr = io.StringIO()
                with mock.patch.object(sys, 'argv', argv), mock.patch.object(M, 'localize', side_effect=exception), \
                        contextlib.redirect_stderr(stderr):
                    self.assertEqual(M.main(), 1)
                failure = json.loads((output / 'failure.json').read_text())
                self.assertEqual(failure['status'], 'failed')
                self.assertFalse(failure['complete'])
                self.assertFalse(failure['biological_interpretation_available'])
                self.assertEqual(failure['error_type'], type(exception).__name__)
                self.assertEqual(json.loads(stderr.getvalue()), failure)
                self.assertEqual([p.name for p in output.iterdir()], ['failure.json'])

    def test_existing_output_is_not_modified_when_cli_fails(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            prior = output / 'summary.json'
            prior.write_text('{"numeric_existing_result":17}\n')
            before = M.digest(prior)
            argv = ['validate_irregular_anchors.py', '--export-dir', temporary,
                    '--assembly', str(output / 'unused'), '--output-dir', str(output)]
            with mock.patch.object(sys, 'argv', argv), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(M.main(), 1)
            self.assertEqual(M.digest(prior), before)
            self.assertFalse((output / 'failure.json').exists())
            self.assertEqual([p.name for p in output.iterdir()], ['summary.json'])

    def test_partial_staging_is_not_published_on_io_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / 'unpublished'
            def fail(*args):
                stage = args[5]
                (stage / 'details.json').write_text('[{"numeric":17}]\n')
                raise OSError('numeric spool failure')
            with mock.patch.object(M, '_localize', side_effect=fail), self.assertRaises(OSError):
                M.localize(temporary, str(Path(temporary) / 'unused'), None, None, None,
                           output, M.DEFAULTS)
            self.assertFalse(output.exists())
            self.assertEqual(list(Path(temporary).iterdir()), [])


if __name__ == '__main__':
    unittest.main()
