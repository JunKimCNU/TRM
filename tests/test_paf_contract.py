"""Numeric PAF contract regressions; no generated sequence or mapper execution.

Coordinates and counts are explicit in fixtures.  Ambiguous bases still consume
coordinates, while minimap2 excludes them from PAF column 11.  These tests do not
claim that a CIGAR/count combination authenticates the underlying bases.
"""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import validate_terminal_anchors as M
import validate_irregular_anchors as IRREGULAR


def paf(cigar='100M', qspan=100, tspan=100, matches=100, block=100,
        tags=(), strand='+', score='180', tp='P', qstart=10, tstart=1000):
    """Build metadata only; expected spans/counts are deliberately independent."""
    fields = ['anchor', 200, qstart, qstart + qspan, strand, 'record', 2000,
              tstart, tstart + tspan, matches, block, 60]
    return '\t'.join(map(str, fields + [f'tp:A:{tp}', f'AS:i:{score}',
                                       f'cg:Z:{cigar}', *tags]))


def replace_column(line, index, value):
    fields = line.split('\t')
    fields[index] = str(value)
    return '\t'.join(fields)


class PafContractTests(unittest.TestCase):
    def rejects(self, lines):
        for line in lines:
            with self.subTest(paf=line), self.assertRaises(ValueError):
                M.parse_paf(line)

    def test_legacy_row_shape_and_values_are_unchanged(self):
        expected = dict(query='anchor', query_length=200, query_start=10,
                        query_end=110, strand='+', target='record',
                        target_length=2000, target_start=1000, target_end=1100,
                        matches=93, block_length=100, mapq=60,
                        alignment_score=180, alignment_type='P', identity=.93,
                        aligned_query_bp=100, query_coverage=.5)
        self.assertEqual(M.parse_paf(paf(matches=93)), expected)

    def test_zero_nn_and_consistent_optional_nm_do_not_change_row(self):
        original = M.parse_paf(paf(matches=93))
        for tags in [('nn:i:0',), ('NM:i:7',), ('nn:i:+0', 'NM:i:+7')]:
            with self.subTest(tags=tags):
                self.assertEqual(M.parse_paf(paf(matches=93, tags=tags)), original)

    def test_ambiguous_M_columns_preserve_coordinate_consumption(self):
        row = M.parse_paf(paf(matches=93, block=96, tags=('nn:i:4', 'NM:i:7')))
        self.assertEqual((row['query_start'], row['query_end']), (10, 110))
        self.assertEqual((row['target_start'], row['target_end']), (1000, 1100))
        self.assertEqual(row['aligned_query_bp'], 100)
        self.assertEqual(row['query_coverage'], .5)
        self.assertEqual(row['identity'], 93 / 96)

    def test_equals_can_include_ambiguous_equal_bases(self):
        row = M.parse_paf(paf('100=', matches=97, block=97,
                              tags=('nn:i:3', 'NM:i:3')))
        self.assertEqual(row['identity'], 1)
        self.assertEqual(row['aligned_query_bp'], 100)

    def test_X_can_include_ambiguous_mismatching_bases(self):
        row = M.parse_paf(paf('90=10X', matches=90, block=96,
                              tags=('nn:i:4', 'NM:i:10')))
        self.assertEqual(row['identity'], 90 / 96)

    def test_ambiguous_inserted_bases_are_excluded_from_block_only(self):
        row = M.parse_paf(paf('80M10I10M', qspan=100, tspan=90,
                              matches=90, block=95, tags=('nn:i:5', 'NM:i:10')))
        self.assertEqual(row['query_end'] - row['query_start'], 100)
        self.assertEqual(row['target_end'] - row['target_start'], 90)
        self.assertEqual(row['identity'], 90 / 95)

    def test_ambiguous_deleted_bases_are_excluded_from_block_only(self):
        row = M.parse_paf(paf('80M10D10M', qspan=90, tspan=100,
                              matches=90, block=95, tags=('nn:i:5', 'NM:i:10')))
        self.assertEqual(row['query_end'] - row['query_start'], 90)
        self.assertEqual(row['target_end'] - row['target_start'], 100)
        self.assertEqual(row['query_coverage'], .45)

    def test_mixed_M_equals_X_and_indels(self):
        row = M.parse_paf(paf('85M5=4X3I2D', qspan=97, tspan=96,
                              matches=86, block=94, tags=('nn:i:5', 'NM:i:13')))
        self.assertEqual(row['identity'], 86 / 94)
        self.assertEqual(row['aligned_query_bp'], 97)

    def test_mixed_M_equals_match_bounds_remain_attainable(self):
        for matches in (88, 93, 98):
            with self.subTest(matches=matches):
                row = M.parse_paf(paf('90=10M', matches=matches, block=98,
                                      tags=('nn:i:2', f'NM:i:{100 - matches}')))
                self.assertEqual(row['matches'], matches)

    def test_ambiguity_lower_bound_is_floored_at_zero(self):
        # The 12 ambiguous columns can cover all 5 '=' plus other columns.
        row = M.parse_paf(paf('5=95X', matches=0, block=88,
                              tags=('nn:i:12', 'NM:i:100')))
        self.assertEqual(row['identity'], 0)

    def test_NM_may_exceed_nonambiguous_block_length(self):
        # 10 inserted + 10 deleted columns; 5 ambiguous gap columns excluded.
        row = M.parse_paf(paf('10I10D', qspan=10, tspan=10,
                              matches=0, block=15, tags=('nn:i:5', 'NM:i:20')))
        self.assertEqual(row['identity'], 0)
        self.assertEqual(row['block_length'], 15)

    def test_nn_is_valid_without_NM(self):
        row = M.parse_paf(paf('100=', matches=97, block=97, tags=('nn:i:3',)))
        self.assertEqual(row['matches'], 97)

    def test_missing_nn_does_not_infer_hidden_ambiguous_count(self):
        self.rejects([paf('100=', matches=97, block=97),
                      paf(matches=93, block=96, tags=('NM:i:7',))])

    def test_explicit_NM_must_agree_even_without_nn(self):
        self.rejects([paf(matches=93, tags=('NM:i:6',)),
                      paf(matches=93, tags=('NM:i:8',)),
                      paf(matches=93, block=96, tags=('nn:i:4', 'NM:i:3'))])

    def test_nn_cannot_repair_query_or_reference_coordinate_mismatch(self):
        self.rejects([paf(qspan=96, matches=93, block=96, tags=('nn:i:4',)),
                      paf(tspan=96, matches=93, block=96, tags=('nn:i:4',))])

    def test_block_must_account_for_all_and_only_ambiguous_columns(self):
        self.rejects([paf(matches=90, block=95, tags=('nn:i:4',)),
                      paf(matches=90, block=97, tags=('nn:i:4',)),
                      paf(matches=0, block=1, tags=('nn:i:100',)),
                      paf(matches=0, block=0, tags=('nn:i:100',))])

    def test_match_bounds_reject_impossible_equals_X_and_indel_counts(self):
        self.rejects([paf('100=', matches=96, block=97, tags=('nn:i:3',)),
                      paf('90=10M', matches=87, block=98, tags=('nn:i:2',)),
                      paf('90=10X', matches=91, block=96, tags=('nn:i:4',)),
                      paf('80M10I10M', tspan=90, matches=91, block=95,
                          tags=('nn:i:5',)),
                      paf('10I10D', qspan=10, tspan=10, matches=1, block=20)])

    def test_both_mapping_strands_keep_forward_coordinate_fields(self):
        plus = M.parse_paf(paf(matches=93, block=96, tags=('nn:i:4',)))
        minus = M.parse_paf(paf(matches=93, block=96, tags=('nn:i:4',), strand='-'))
        self.assertEqual({**plus, 'strand': '-'}, minus)
        tract = dict(strand='+', left_end_distance=500, right_end_distance=0,
                     tract_length=550)
        for row, boundary, end in [(plus, 1100, 'right'), (minus, 1000, 'left')]:
            with self.subTest(strand=row['strand']):
                geometry = M.anchor_geometry(tract, 'left', row, M.DEFAULTS)
                self.assertEqual(geometry['mapped_junction_boundary0'], boundary)
                self.assertEqual(geometry['outward_target_end'], end)
                self.assertEqual(geometry['unaligned_anchor_junction_bp'], 90)
                self.assertEqual(geometry['aligned_query_bp'], 100)

    def test_signed_alignment_scores_and_alignment_types(self):
        for score, expected in [('+180', 180), ('-180', -180), ('0', 0)]:
            with self.subTest(score=score):
                self.assertEqual(M.parse_paf(paf(score=score))['alignment_score'], expected)
        for tp in ('P', 'S', 'I', 'i'):
            with self.subTest(tp=tp):
                self.assertEqual(M.parse_paf(paf(tp=tp))['alignment_type'], tp)

    def test_required_and_optional_consumed_tags_have_exact_types(self):
        base = paf(tags=('nn:i:0', 'NM:i:0'))
        wrong = [('cg:Z:', 'cg:i:'), ('AS:i:', 'AS:Z:'), ('tp:A:', 'tp:Z:'),
                 ('nn:i:', 'nn:f:'), ('NM:i:', 'NM:Z:')]
        self.rejects([base.replace(good, bad) for good, bad in wrong])

    def test_duplicate_tag_is_invalid_even_with_different_type(self):
        self.rejects([paf(tags=(tag,)) for tag in
                      ('cg:Z:100M', 'AS:i:180', 'tp:Z:P')] +
                     [paf(tags=('nn:i:0', 'nn:i:0')),
                      paf(tags=('NM:i:0', 'NM:Z:0')),
                      paf(tags=('zz:Z:one', 'zz:Z:two'))])

    def test_malformed_tag_names_and_missing_fields(self):
        self.rejects([paf(tags=(tag,)) for tag in
                      ('', 'nn', 'nn:i', ':i:0', 'n:i:0', 'nnn:i:0', '1n:i:0')])

    def test_bad_integer_tag_lexemes_are_not_silently_coerced(self):
        self.rejects([paf(score=value) for value in ('1.0', ' 180', '180 ', 'NaN', '+', '')])
        for tag in ('nn', 'NM'):
            self.rejects([paf(tags=(f'{tag}:i:{value}',)) for value in
                          ('-1', '0.0', ' 0', '0 ', '', 'NaN', 'True', '1_0')])

    def test_nondecimal_mandatory_integers_are_invalid(self):
        self.rejects([replace_column(paf(), index, value)
                      for index, value in [(1, '2e2'), (2, ' 10'), (3, '110 '),
                                           (7, '1_000'), (9, '100.0'), (11, '+60')]])

    def test_missing_required_tags_and_unknown_alignment_type(self):
        base = paf()
        self.rejects(['\t'.join(field for field in base.split('\t')
                                if not field.startswith(tag + ':'))
                      for tag in ('cg', 'AS', 'tp')] + [paf(tp='p'), paf(tp=''), paf(tp='PP')])

    def test_genomic_cigar_rejects_clipping_splicing_padding_and_malformed_lengths(self):
        self.rejects([paf(cigar) for cigar in
                      ('', '0M100M', '01M99M', '100m', '100Mgarbage', '100', 'M',
                       '100N', '100S', '100H', '100P', '100B', '99M1S')])

    def test_coordinate_count_and_mapq_ranges_still_apply(self):
        self.rejects([replace_column(paf(), index, value) for index, value in
                      [(0, ''), (5, ''), (4, '?'), (2, -1), (2, 110), (3, 201),
                       (7, 1100), (8, 2001), (9, -1), (9, 101), (10, 0),
                       (11, -1), (11, 256)]] + ['anchor\t200'])

    def test_unconsumed_valid_extension_tags_and_crlf_are_tolerated(self):
        base = paf()
        self.assertEqual(M.parse_paf(base), M.parse_paf(base + '\tzz:Z:one:two\tdv:f:0.07\r\n'))

    def test_irregular_adapter_uses_identical_shared_parser(self):
        self.assertIs(IRREGULAR.parse_paf, M.parse_paf)
        valid = paf('100=', matches=97, block=97, tags=('nn:i:3', 'NM:i:3'))
        self.assertEqual(IRREGULAR.parse_paf(valid), M.parse_paf(valid))
        with self.assertRaises(ValueError):
            IRREGULAR.parse_paf(valid.replace('nn:i:3', 'nn:Z:3'))

    def test_detection_and_geometry_thresholds_are_not_relaxed(self):
        self.assertEqual(M.DEFAULTS, dict(min_mapq=20, min_identity=.85,
                         min_query_coverage=.80, min_aligned_bp=200,
                         max_junction_gap=20, end_window=5000, read_end_window=50,
                         competing_score_fraction=.90, long_tract_bp=500))
        self.assertIs(IRREGULAR.DEFAULTS, M.DEFAULTS)


if __name__ == '__main__':
    unittest.main()
