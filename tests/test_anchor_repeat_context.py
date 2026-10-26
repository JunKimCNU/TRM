"""Coordinate/provenance tests only; no synthetic sequence or read fixture.

The small metadata fixture uses published DSM3286 contig/motif identifiers and
existing exact-audit interval coordinates. Mutations test numeric validation.
"""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    'anchor_repeat_context', Path(__file__).resolve().parents[1] / 'scripts/anchor_repeat_context.py')
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
SHA = '39c48a40a88e896c58e5a9528eb84f8f12c98a7572a389fec0f90c80f57d6c63'
FAMILY = 'AACCCTGACT'
CONTIG = 'CP061014.1'
LENGTH = 4100424


def fixture():
    # Existing three DSM audit rows intersect the strict geometry anchor span.
    rows = [(4100112, 4100154, 'TTAGTCAGGG'),
            (4100197, 4100233, 'TCAGGGTTAG'),
            (4100357, 4100412, 'GGGTTAGTCA')]
    hits = [dict(contig_id=CONTIG, start_0based=a, end_0based_exclusive=b,
                 tract_bp=b-a, complete_copies=(b-a)//10, partial_copy_bp=(b-a)%10,
                 matched_rotation=rotation, strand='+',
                 window_censored_left=False, window_censored_right=False)
            for a, b, rotation in rows]
    return dict(schema_version=1, method=M.METHOD, coordinate_system=M.COORDINATES,
                settings=dict(window_bp=5000, min_complete_copies=3),
                assemblies=[dict(input_sha256=SHA, expected_outward_motif='TTAGTCAGGG',
                                 full_fasta_parsed=True, reported_sequence_lengths_all_match=True,
                                 contigs=[dict(contig_id=CONTIG, length_bp=LENGTH)], hits=hits,
                                 summary=dict(contig_count=1, assembly_length_bp=LENGTH,
                                              exact_tract_count=3))])


class AnchorRepeatContextTests(unittest.TestCase):
    def load(self, audit=None, sha=SHA, family=FAMILY, lengths=None):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'audit.json'
            path.write_text(json.dumps(fixture() if audit is None else audit))
            return M.load_repeat_audit(path, sha, family, {CONTIG: LENGTH} if lengths is None else lengths)

    def test_real_coordinate_overlap_and_denominator(self):
        result = M.annotate_alignment(self.load(), CONTIG, 4100119, 4100413)
        self.assertEqual(result['reference_span_bp'], 294)
        self.assertEqual(result['known_target_repeat_overlap_bp'], 126)
        self.assertAlmostEqual(result['known_target_repeat_overlap_fraction'], 3/7)
        self.assertEqual(result['audited_span_bp'], 294)
        self.assertEqual(result['known_target_repeat_intervals'],
                         [[4100119, 4100154], [4100197, 4100233], [4100357, 4100412]])
        self.assertFalse(result['unique_nonrepeat_subtelomere_confirmed'])

    def test_union_overlapping_duplicate_and_adjacent_intervals(self):
        self.assertEqual(M.union_intervals([[30, 50], [10, 30], [15, 25], [10, 30], [60, 70]]),
                         [[10, 50], [60, 70]])

    def test_duplicate_hit_rows_do_not_double_count(self):
        audit = fixture()
        a = audit['assemblies'][0]
        a['hits'].append(copy.deepcopy(a['hits'][0]))
        a['summary']['exact_tract_count'] += 1
        result = M.annotate_alignment(self.load(audit), CONTIG, 4100119, 4100413)
        self.assertEqual(result['known_target_repeat_overlap_bp'], 126)

    def test_masking_ignores_hit_strand(self):
        for strand in ('+', '-', '+/-', '-/+'):
            audit = fixture()
            for hit in audit['assemblies'][0]['hits']:
                hit['strand'] = strand
            result = M.annotate_alignment(self.load(audit), CONTIG, 4100119, 4100413)
            self.assertEqual(result['known_target_repeat_overlap_bp'], 126)

    def test_full_scope_without_overlap_is_not_uniqueness(self):
        result = M.annotate_alignment(self.load(), CONTIG, 0, 1000)
        self.assertEqual(result['repeat_context_status'], 'no_known_overlap_in_audited_span')
        self.assertFalse(result['unique_nonrepeat_subtelomere_confirmed'])

    def test_partial_scope_without_overlap(self):
        result = M.annotate_alignment(self.load(), CONTIG, 4500, 5500)
        self.assertEqual(result['repeat_context_status'], 'partial_audit_scope_no_known_overlap')
        self.assertEqual(result['audited_span_fraction'], 0.5)

    def test_outside_scope_is_unknown_not_negative(self):
        result = M.annotate_alignment(self.load(), CONTIG, 6000, 7000)
        self.assertEqual(result['repeat_context_status'], 'outside_audit_scope')
        self.assertEqual(result['audited_span_bp'], 0)

    def test_partial_scope_with_overlap_is_explicit(self):
        result = M.annotate_alignment(self.load(), CONTIG, LENGTH-6000, LENGTH)
        self.assertEqual(result['repeat_context_status'], 'known_target_repeat_overlap')
        self.assertEqual(result['audited_scope_status'], 'partial')
        self.assertEqual(result['audited_span_bp'], 5000)
        self.assertEqual(result['known_target_repeat_overlap_bp'], 133)

    def test_short_contig_windows_are_unioned(self):
        audit = fixture()
        a = audit['assemblies'][0]
        a['contigs'][0]['length_bp'] = 8000
        a['hits'] = []
        a['summary'].update(assembly_length_bp=8000, exact_tract_count=0)
        index = self.load(audit, lengths={CONTIG: 8000})
        self.assertEqual(index['windows'][CONTIG], [[0, 8000]])
        self.assertEqual(M.annotate_alignment(index, CONTIG, 0, 8000)['audited_span_bp'], 8000)

    def test_hash_mismatch_or_ambiguous_assembly_rejected(self):
        with self.assertRaises(ValueError):
            self.load(sha='0'*64)
        audit = fixture()
        audit['assemblies'].append(copy.deepcopy(audit['assemblies'][0]))
        with self.assertRaises(ValueError):
            self.load(audit)

    def test_noncanonical_or_wrong_family_rejected(self):
        for family in ('TTAGTCAGGG', 'AACCCT', '', 'acgt'):
            with self.assertRaises(ValueError):
                self.load(family=family)

    def test_rotation_and_reverse_complement_family(self):
        # Existing outward/left motif names; no generated sequence records.
        for motif in ('TTAGTCAGGG', 'GGGTTAGTCA', 'CCCTGACTAA'):
            audit = fixture()
            audit['assemblies'][0]['expected_outward_motif'] = motif
            self.assertEqual(self.load(audit)['provenance']['family_key'], FAMILY)

    def test_length_name_duplicate_and_completeness_rejected(self):
        mutations = [lambda a: a['contigs'][0].update(length_bp=LENGTH+1),
                     lambda a: a['contigs'][0].update(contig_id='other'),
                     lambda a: a['contigs'].append(copy.deepcopy(a['contigs'][0])),
                     lambda a: a.update(full_fasta_parsed=False),
                     lambda a: a.update(reported_sequence_lengths_all_match=False)]
        for mutation in mutations:
            audit = fixture()
            mutation(audit['assemblies'][0])
            with self.assertRaises(ValueError):
                self.load(audit)

    def test_invalid_schema_method_and_window(self):
        mutations = [lambda a: a.update(schema_version=True), lambda a: a.update(schema_version=2),
                     lambda a: a.update(method='genome-wide'),
                     lambda a: a.update(coordinate_system='one-based'),
                     lambda a: a['settings'].update(window_bp=0),
                     lambda a: a['settings'].update(min_complete_copies=True)]
        for mutation in mutations:
            audit = fixture()
            mutation(audit)
            with self.assertRaises(ValueError):
                self.load(audit)

    def test_hit_bounds_window_and_length_validation(self):
        changes = [dict(start_0based=-1), dict(end_0based_exclusive=LENGTH+1),
                   dict(start_0based=6000, end_0based_exclusive=6042), dict(tract_bp=43),
                   dict(contig_id='other'), dict(start_0based=True)]
        for change in changes:
            audit = fixture()
            audit['assemblies'][0]['hits'][0].update(change)
            with self.assertRaises(ValueError):
                self.load(audit)

    def test_hit_period_copies_strand_and_censor_validation(self):
        changes = [dict(complete_copies=3), dict(partial_copy_bp=0),
                   dict(matched_rotation='TTAGGG'), dict(strand='?'),
                   dict(window_censored_right=0)]
        for change in changes:
            audit = fixture()
            audit['assemblies'][0]['hits'][0].update(change)
            with self.assertRaises(ValueError):
                self.load(audit)
        audit = fixture()
        audit['settings']['min_complete_copies'] = 6
        with self.assertRaises(ValueError):
            self.load(audit)

    def test_summary_count_disagreement_rejected(self):
        for field in ('contig_count', 'assembly_length_bp', 'exact_tract_count'):
            audit = fixture()
            audit['assemblies'][0]['summary'][field] += 1
            with self.assertRaises(ValueError):
                self.load(audit)

    def test_invalid_alignment_coordinates_rejected(self):
        index = self.load()
        for start, end in ((-1, 2), (2, 2), (3, 2), (0, LENGTH+1), (True, 3), (0, 3.0)):
            with self.assertRaises(ValueError):
                M.annotate_alignment(index, CONTIG, start, end)
        with self.assertRaises(ValueError):
            M.annotate_alignment(index, 'other', 0, 1000)

    def test_duplicate_json_key_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'audit.json'
            path.write_text('{"schema_version":1,"schema_version":1}')
            with self.assertRaises(ValueError):
                M.load_repeat_audit(path, SHA, FAMILY, {CONTIG: LENGTH})

    def test_annotation_does_not_mutate_index(self):
        index = self.load()
        before = copy.deepcopy(index)
        M.annotate_alignment(index, CONTIG, 4100119, 4100413)
        self.assertEqual(before, index)
        self.assertFalse(index['provenance']['changes_native_classification'])
        self.assertFalse(index['provenance']['changes_geometry_counts'])


if __name__ == '__main__':
    unittest.main()
