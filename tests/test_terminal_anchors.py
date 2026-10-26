import importlib.util
import copy
import json
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location('terminal_anchors', Path(__file__).resolve().parents[1] / 'scripts/validate_terminal_anchors.py')
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def tract(side='left', strand='+', length=547):
    # Left anchor precedes a right-end repeat; right anchor follows a left-end repeat.
    start = 1000 if side == 'left' else 0
    return dict(tract_id='t1', read_ordinal=1, read_id='r1', read_length=1000+length,
                tract_index=1, motif='TTAGGG', family_key='AACCCT', start0=start,
                end0=start+length, tract_length=length, identity=1.0, strand=strand,
                left_end_distance=start, right_end_distance=1000+length-start-length,
                left_anchor_id='aL' if side == 'left' else '',
                right_anchor_id='aR' if side == 'right' else '')


def aln(query='aL', strand='+', start=89000, qstart=0, qend=1000, mapq=60, score=1900, tp='P'):
    n = qend-qstart
    line = f'{query}\t1000\t{qstart}\t{qend}\t{strand}\tchr1\t90000\t{start}\t{start+n}\t{n}\t{n}\t{mapq}\ttp:A:{tp}\tAS:i:{score}\tcg:Z:{n}M'
    return M.parse_paf(line)


class TerminalAnchorTests(unittest.TestCase):
    def test_four_orientation_combinations(self):
        cases = [('left', '+', '+', 88500, 'right'),
                 ('left', '+', '-', 500, 'left'),
                 ('right', '-', '+', 500, 'left'),
                 ('right', '-', '-', 88500, 'right')]
        for side, repeat_strand, mapping_strand, start, expected_end in cases:
            with self.subTest(side=side, mapping=mapping_strand):
                t = tract(side, repeat_strand)
                name = 'aL' if side == 'left' else 'aR'
                result = M.evaluate_tract(t, {name: [aln(name, mapping_strand, start)]}, M.DEFAULTS)
                self.assertTrue(result['joint_long_terminal_context'])
                self.assertEqual(result['anchors'][0]['outward_target_end'], expected_end)
                self.assertFalse(result['physical_chromosome_end_confirmed'])

    def test_wrong_repeat_polarity(self):
        result = M.evaluate_tract(tract(strand='-'), {'aL': [aln()]}, M.DEFAULTS)
        self.assertFalse(result['joint_terminal_context'])
        self.assertEqual(result['anchors'][0]['polarity'], 'opposite')

    def test_ambiguous_repeat_strand(self):
        result = M.evaluate_tract(tract(strand='ambiguous'), {'aL': [aln()]}, M.DEFAULTS)
        self.assertFalse(result['joint_terminal_context'])

    def test_internal_assembly_anchor_is_not_terminal(self):
        result = M.evaluate_tract(tract(), {'aL': [aln(start=40000)]}, M.DEFAULTS)
        self.assertFalse(result['joint_terminal_context'])

    def test_read_end_requirement_is_separate(self):
        t = tract()
        t['right_end_distance'] = 100
        result = M.evaluate_tract(t, {'aL': [aln()]}, M.DEFAULTS)
        self.assertEqual(result['localization_status'], 'assembly_end_anchor_only')
        self.assertFalse(result['joint_long_terminal_context'])

    def test_low_or_unknown_mapq_and_competing_hits(self):
        for rows in [[aln(mapq=0)], [aln(mapq=255)], [aln(), aln(start=100, score=1850, tp='S')]]:
            result = M.evaluate_tract(tract(), {'aL': rows}, M.DEFAULTS)
            self.assertFalse(result['joint_terminal_context'])

    def test_low_identity_and_partial_alignment(self):
        low = aln()
        low['identity'] = 0.84
        for rows in [[low], [aln(qend=790)]]:
            self.assertIsNone(M.select_alignment(rows, M.DEFAULTS)[0])

    def test_query_junction_gap_both_mapping_strands(self):
        for strand in ['+', '-']:
            left = M.anchor_geometry(tract(), 'left', aln(strand=strand, start=100, qend=970), M.DEFAULTS)
            right = M.anchor_geometry(tract('right', '-'), 'right', aln(strand=strand, start=100, qstart=30), M.DEFAULTS)
            self.assertEqual(left['unaligned_anchor_junction_bp'], 30)
            self.assertEqual(right['unaligned_anchor_junction_bp'], 30)
            self.assertFalse(left['near_junction'])
            self.assertFalse(right['near_junction'])

    def test_two_sided_not_promoted(self):
        t = tract()
        t['right_anchor_id'] = 'aR'
        result = M.evaluate_tract(t, {'aL': [aln()], 'aR': [aln('aR', start=100)]}, M.DEFAULTS)
        self.assertEqual(result['localization_status'], 'two_sided_anchor_context')
        self.assertFalse(result['joint_long_terminal_context'])

    def test_short_observation_not_combined_into_long(self):
        result = M.evaluate_tract(tract(length=499), {'aL': [aln()]}, M.DEFAULTS)
        self.assertTrue(result['joint_terminal_context'])
        self.assertFalse(result['joint_long_terminal_context'])

    def test_organelle_excluded_and_annotation_optional(self):
        for compartment, expected in [('organelle', False), ('nuclear_chromosome', True)]:
            result = M.evaluate_tract(tract(), {'aL': [aln()]}, M.DEFAULTS, {'chr1': compartment})
            self.assertEqual(result['joint_nuclear_terminal_context'], expected)
        result = M.evaluate_tract(tract(), {'aL': [aln()]}, M.DEFAULTS)
        self.assertFalse(result['joint_nuclear_terminal_context'])

    def test_no_anchor_is_unresolved(self):
        t = tract()
        t['left_anchor_id'] = ''
        self.assertEqual(M.evaluate_tract(t, {}, M.DEFAULTS)['localization_status'], 'no_exportable_adjacent_anchor')

    def test_invalid_paf_and_cigar(self):
        base = 'q\t1000\t0\t1000\t+\tt\t2000\t0\t1000\t1000\t1000\t60\ttp:A:P\tAS:i:1000\tcg:Z:1000M'
        for line in [base.replace('1000M', '999M'), base.replace('cg:Z:1000M', 'cg:Z:1000N'),
                     base.replace('\t1000\t60', '\t999\t60'), base + '\ttp:A:P', 'q\t2']:
            with self.assertRaises(ValueError):
                M.parse_paf(line)
        with self.assertRaises(ValueError):
            M.parse_paf('q\t200\t0\t200\t+\tt\t200\t0\t200\t400\t400\t60\ttp:A:P\tAS:i:100\tcg:Z:200I200D')

    def test_receipt_requires_exact_source_content_and_command(self):
        with tempfile.TemporaryDirectory() as folder:
            paths = {k: Path(folder) / k for k in ['reference', 'queries', 'paf']}
            for path in paths.values():
                path.write_text('data')
            receipt = {'status': 'complete', 'exit_code': 0, 'command': ['minimap2', '-c', '--secondary=yes',
                        str(paths['reference']), str(paths['queries'])],
                       **{k: {'path': str(p), 'sha256': M.digest(p)} for k, p in paths.items()}}
            M.validate_receipt(receipt, paths)
            paths['queries'].write_text('different')
            with self.assertRaises(ValueError):
                M.validate_receipt(receipt, paths)

    def test_incomplete_or_stale_export_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for name in ['tracts.tsv', 'anchors.fasta']:
                (root/name).write_text('data')
            summary = {'status': 'complete', 'complete_input': True, 'outputs':
                       {name: {'sha256': M.digest(root/name)} for name in ['tracts.tsv', 'anchors.fasta']}}
            M.validate_export(summary, root)
            summary['complete_input'] = False
            with self.assertRaises(ValueError):
                M.validate_export(summary, root)
            summary['complete_input'] = True
            (root/'tracts.tsv').write_text('changed')
            with self.assertRaises(ValueError):
                M.validate_export(summary, root)

    def test_separate_same_read_context_never_merges_tract_lengths(self):
        source = tract(length=97)
        source['right_end_distance'] = 550
        source = M.evaluate_tract(source, {'aL': [aln()]}, M.DEFAULTS, {'chr1': 'nuclear_chromosome'})
        longer = tract(length=547)
        longer.update(tract_id='t2', start0=1099, end0=1646, left_anchor_id='',
                      left_end_distance=1099, right_end_distance=1)
        longer = M.evaluate_tract(longer, {}, M.DEFAULTS)
        context = M.same_read_context([source, longer], M.DEFAULTS)
        self.assertEqual(len(context), 1)
        self.assertEqual(context[0]['long_tract_bp'], 547)
        self.assertEqual(context[0]['intervening_read_gap_bp'], 2)
        self.assertTrue(context[0]['outward_read_order_consistent'])
        self.assertFalse(context[0]['tracts_merged'])
        self.assertFalse(longer['joint_long_terminal_context'])
        longer['read_ordinal'] = 2
        self.assertEqual(M.same_read_context([source, longer], M.DEFAULTS), [])

    def test_repeat_context_annotation_preserves_geometry_and_record_counting(self):
        # Coordinate-only diagnostics: no generated read or assembly sequences.
        first = M.evaluate_tract(tract(), {'aL': [aln()]}, M.DEFAULTS)
        second = copy.deepcopy(first)
        second['tract_id'] = 'same_read_other_tract'
        before = copy.deepcopy([first, second])
        def annotate(index, target, start, end):
            return dict(repeat_context_status='known_target_repeat_overlap',
                        known_target_repeat_overlap_bp=10, reference_span_bp=end-start,
                        audited_span_bp=end-start, unique_nonrepeat_subtelomere_confirmed=False)
        result = M.add_repeat_context([first, second], {'provenance': {'scope': 'coordinate_test'}}, annotate)
        self.assertEqual(result['joint_geometry_records_with_known_repeat_overlap'], 1)
        self.assertEqual(result['anchor_status_counts']['known_target_repeat_overlap'], 2)
        self.assertFalse(result['geometry_counts_changed'])
        for original, updated in zip(before, [first, second]):
            self.assertEqual({k: v for k, v in original.items() if k != 'anchors'},
                             {k: v for k, v in updated.items() if k != 'anchors'})
            for anchor_before, anchor_after in zip(original['anchors'], updated['anchors']):
                self.assertTrue(all(anchor_after[k] == v for k, v in anchor_before.items()))

    def test_repeat_context_unqualified_anchor_is_not_measured_zero(self):
        outcome = M.evaluate_tract(tract(), {'aL': [aln(mapq=0)]}, M.DEFAULTS)
        def must_not_run(*args):
            raise AssertionError('unqualified alignment must not receive zero overlap')
        result = M.add_repeat_context([outcome], {'provenance': {}}, must_not_run)
        self.assertEqual(outcome['anchors'][0]['repeat_context_status'], 'not_assessed_unqualified_alignment')
        self.assertNotIn('known_target_repeat_overlap_bp', outcome['anchors'][0])
        self.assertEqual(result['joint_geometry_records_with_known_repeat_overlap'], 0)


if __name__ == '__main__':
    unittest.main()
