"""Reference comparison must never turn repeat abundance into a main TRM."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
from evidence_batch import canonical
from trm_assessment import (POLICY, assess_candidate, catalog_index, classify,
                            reference_agreement, reverse_complement, select_representative,
                            library_context)


def assignment(motif='CTCGGTTATGGG', scope='genus', status='supported'):
    return {'related_motifs': motif, 'reference_scope': scope, 'reference_status': status}


def candidate(motif, reads, **extra):
    return {'motif': motif, 'canonical_family': canonical(motif), 'reads': reads, **extra}


def observation(motif, reads, recovery=None, **extra):
    return {'motif': motif, 'canonical_family': canonical(motif), 'reads': reads,
            'recovery': recovery or ('reported_positive' if reads else 'explicit_zero'),
            'source': 'saved_repeat_evidence', **extra}


class ReferenceAgreementTests(unittest.TestCase):
    def test_allium_s2_both_families_preserved_and_main_trm_unassigned(self):
        expected = candidate('CTCGGTTATGGG', 16)
        top = candidate('TTTAGGG', 106)
        candidates = [top, expected]
        before = copy.deepcopy(candidates)
        result = reference_agreement(candidates, top, assignment(), [observation(expected['motif'], 16)])
        self.assertEqual(result['status'], 'observed_family_conflict')
        self.assertTrue(result['requires_conflict_review'])
        self.assertTrue(result['observed_family_conflict'])
        self.assertEqual(result['expectation_strength'], 'genus_expectation')
        self.assertEqual(result['representative_motif'], 'TTTAGGG')
        self.assertEqual(result['expected_families'][0]['reads'], 16)
        self.assertEqual(result['main_trm_status'], 'not_confirmed')
        self.assertIsNone(result['main_trm_motif'])
        self.assertEqual(candidates, before)

    def test_tetrahymena_zero_keeps_provisional_novel_observation_separate(self):
        top = candidate('AACAATATTGTAATAACCCTG', 556, novel_candidate=True,
                        novel_evidence_level='provisional_candidate')
        result = reference_agreement([top], top, assignment('TTGGGG', 'species'),
                                     [observation('TTGGGG', 0)])
        self.assertEqual(result['status'], 'expected_reference_not_observed')
        self.assertTrue(result['expected_reference_not_observed'])
        self.assertFalse(result['observed_family_conflict'])
        self.assertTrue(result['requires_conflict_review'])
        self.assertEqual(result['representative_motif'], top['motif'])
        self.assertEqual(top['novel_evidence_level'], 'provisional_candidate')
        self.assertIsNone(result['main_trm_motif'])

    def test_rotation_reverse_complement_primitive_family_agreement(self):
        for motif in ('TTAGGG', 'TAGGGT', 'CCCTAA', 'TTAGGGTTAGGG'):
            with self.subTest(motif=motif):
                top = candidate(motif, 40)
                result = reference_agreement([top], top, assignment('TTAGGG', 'species'),
                                             [observation('TTAGGG', 40)])
                self.assertEqual(result['status'], 'representative_reference_agreement')
                self.assertFalse(result['observed_family_conflict'])
                self.assertIsNone(result['main_trm_motif'])

    def test_not_reported_is_not_explicit_zero_even_with_zero_placeholder(self):
        placeholder = candidate('TTGGGG', 0)
        top = candidate('ACGTG', 80)
        result = reference_agreement([placeholder, top], top, assignment('TTGGGG', 'species'),
                                     [observation('TTGGGG', None, 'not_reported')])
        self.assertEqual(result['status'], 'expected_reference_unmeasured')
        self.assertFalse(result['expected_reference_not_observed'])
        self.assertIsNone(result['expected_families'][0]['reads'])

    def test_zero_without_explicit_measurement_is_still_unreported(self):
        result = reference_agreement([], None, assignment('TTGGGG', 'species'),
                                     [observation('TTGGGG', 0, 'not_reported')])
        self.assertEqual(result['status'], 'expected_reference_unmeasured')
        self.assertFalse(result['expected_reference_not_observed'])

    def test_trace_uses_existing_read_threshold_without_erasing_positive(self):
        for count in (1, 9, POLICY['minimum_known_reads']):
            with self.subTest(count=count):
                top = candidate('TTGGGG', count)
                result = reference_agreement([top], top, assignment('TTGGGG', 'species'),
                                             [observation('TTGGGG', count)])
                self.assertEqual(result['status'], 'expected_reference_trace' if count < 10
                                 else 'representative_reference_agreement')
                self.assertFalse(result['expected_reference_not_observed'])
                self.assertEqual(result['expected_families'][0]['reads'], count)

    def test_positive_original_zero_current_requires_measurement_review(self):
        expected = candidate('TTGGGG', 0, original_discovery_reads=12, measurement_disagreement=True)
        result = reference_agreement([expected], expected, assignment('TTGGGG', 'species'),
                                     [observation('TTGGGG', 0, 'measurement_disagreement')])
        self.assertEqual(result['status'], 'measurement_disagreement')
        self.assertFalse(result['expected_reference_not_observed'])
        self.assertEqual(result['expected_families'][0]['original_discovery_reads'], 12)

    def test_independent_zero_diagnostic_does_not_overwrite_saved_positive(self):
        expected = candidate('TTGGGG', 12)
        result = reference_agreement([expected], expected, assignment('TTGGGG', 'species'),
                                     [observation('TTGGGG', 0)])
        self.assertEqual(result['status'], 'measurement_disagreement')
        self.assertEqual(result['expected_families'][0]['reads'], 0)
        self.assertEqual(result['expected_families'][0]['candidate_reads'], 12)
        self.assertEqual(expected['reads'], 12)

    def test_original_positive_is_not_lost_when_caller_omits_disagreement_flag(self):
        expected = candidate('TTGGGG', 0, original_discovery_reads=12)
        result = reference_agreement([expected], None, assignment('TTGGGG', 'species'),
                                     [observation('TTGGGG', 0)])
        self.assertEqual(result['status'], 'measurement_disagreement')
        self.assertFalse(result['expected_reference_not_observed'])

    def test_phylogenetic_clade_default_is_not_exact_species_requirement(self):
        top = candidate('TTAGAC', 400)
        result = reference_agreement([top], top, assignment('TTAGGC', 'clade'),
                                     [observation('TTAGGC', 30)])
        self.assertEqual(result['status'], 'broad_expectation_difference')
        self.assertEqual(result['expectation_strength'], 'weak_clade_expectation')
        self.assertFalse(result['observed_family_conflict'])
        self.assertFalse(result['requires_conflict_review'])

    def test_species_and_isolate_strength_remain_distinct_from_genus(self):
        top = candidate('TTAGGC', 30)
        for scope in ('species', 'isolate', 'genus'):
            with self.subTest(scope=scope):
                result = reference_agreement([top], top, assignment('TTAGGG', scope),
                                             [observation('TTAGGG', 20)])
                self.assertEqual(result['reference_scope'], scope)
                self.assertEqual(result['status'], 'observed_family_conflict')
                self.assertEqual(result['expectation_strength'], scope +
                                 ('_expectation' if scope == 'genus' else '_reference'))

    def test_published_candidate_cannot_impose_established_conflict(self):
        top = candidate('TTAGGG', 30)
        result = reference_agreement([top], top, assignment('GGTGTTGG', 'species', 'candidate'),
                                     [observation('GGTGTTGG', 0)])
        self.assertEqual(result['status'], 'no_supported_reference')
        self.assertEqual(result['expectation_strength'], 'published_candidate_nonbinding')
        self.assertFalse(result['requires_conflict_review'])
        self.assertFalse(result['expected_reference_not_observed'])
        self.assertEqual(result['expected_families'][0]['recovery'], 'explicit_zero')

    def test_unexpected_diagnostic_target_is_nonbinding(self):
        top = candidate('CTCGGTTATGGG', 16)
        result = reference_agreement([top], top, assignment(),
                                     [observation('CTCGGTTATGGG', 16), observation('TTTAGGG', 0)])
        self.assertEqual(result['status'], 'representative_reference_agreement')
        unrelated = next(row for row in result['expected_families'] if row['motif'] == 'TTTAGGG')
        self.assertEqual(unrelated['reference_role'], 'nonbinding_comparison')

    def test_unknown_or_out_of_scope_reference_does_not_apply_broader_default(self):
        top = candidate('TTAGGG', 50)
        for status, outcome in (('unknown', 'no_supported_reference'),
                                ('out_of_scope', 'reference_out_of_scope'),
                                ('not_curated', 'no_supported_reference')):
            with self.subTest(status=status):
                result = reference_agreement([top], top, assignment('TTAGGC', 'species', status),
                                             [observation('TTAGGC', 0)])
                self.assertEqual(result['status'], outcome)
                self.assertFalse(result['requires_conflict_review'])
                self.assertFalse(result['expected_reference_not_observed'])

    def test_missing_taxonomic_scope_does_not_create_supported_species_context(self):
        result = reference_agreement([], None, assignment('TTAGGG', ''), [observation('TTAGGG', 0)])
        self.assertEqual(result['status'], 'no_supported_reference')
        self.assertTrue(result['interpretation_limited'])
        self.assertFalse(result['requires_conflict_review'])
        self.assertFalse(result['expected_reference_not_observed'])

    def test_incomplete_qc_library_and_scope_limits_take_precedence(self):
        top = candidate('TTTAGGG', 106)
        cases = [({'analysis_complete': False}, 'analysis_incomplete'),
                 ({'qc_status': 'observed_qc_failure'}, 'qc_limited'),
                 ({'library': {'limited': True}}, 'library_context_limited'),
                 ({'biological_scope': False}, 'taxonomic_scope_unresolved')]
        for options, status in cases:
            with self.subTest(status=status):
                result = reference_agreement([top], top, assignment(),
                                             [observation('CTCGGTTATGGG', 16)], **options)
                self.assertEqual(result['status'], status)
                self.assertTrue(result['interpretation_limited'])
                self.assertTrue(result['observed_family_conflict'])
                self.assertIsNone(result['main_trm_motif'])

    def test_no_reference_or_no_representative_does_not_crash_or_assign_main_trm(self):
        result = reference_agreement([], None, {}, [])
        self.assertEqual(result['status'], 'no_supported_reference')
        self.assertEqual(result['expected_families'], [])
        self.assertIsNone(result['representative_motif'])
        result = reference_agreement([], None, assignment(), [observation('CTCGGTTATGGG', 16)])
        self.assertEqual(result['status'], 'reference_observed_no_representative')
        self.assertIsNone(result['main_trm_motif'])

    def test_multiple_alternatives_must_all_be_measured_zero(self):
        context = assignment('TTAGGG,TTAGGC', 'species')
        result = reference_agreement([], None, context, [observation('TTAGGG', 0)])
        self.assertEqual(result['status'], 'expected_reference_unmeasured')
        self.assertFalse(result['expected_reference_not_observed'])
        result = reference_agreement([], None, context,
                                     [observation('TTAGGG', 0), observation('TTAGGC', 0)])
        self.assertEqual(result['status'], 'expected_reference_not_observed')
        self.assertTrue(result['expected_reference_not_observed'])

    def test_aliases_collapse_without_double_counting(self):
        expected = 'CTCGGTTATGGG'
        rc = reverse_complement(expected)
        top = candidate(expected, 16)
        result = reference_agreement([top], top, assignment(expected + ',' + rc),
                                     [observation(expected, 16), observation(rc, 16)])
        self.assertEqual(len(result['expected_families']), 1)
        self.assertEqual(result['expected_families'][0]['reads'], 16)
        self.assertEqual(result['status'], 'representative_reference_agreement')

    def test_conflicting_alias_measurements_never_select_larger_count(self):
        result = reference_agreement([], None, assignment('TTAGGG', 'species'),
                                     [observation('TTAGGG', 16), observation('CCCTAA', 0)])
        self.assertEqual(result['status'], 'measurement_disagreement')
        self.assertIsNone(result['expected_families'][0]['reads'])
        self.assertEqual(len(result['expected_families'][0]['observations']), 2)

    def test_supporting_reads_alias_and_invalid_explicit_zero(self):
        result = reference_agreement([], None, assignment('TTAGGG', 'species'),
                                     [{'motif': 'TTAGGG', 'supporting_reads': 7, 'recovery': 'reported_positive'}])
        self.assertEqual(result['status'], 'expected_reference_trace')
        with self.assertRaises(ValueError):
            reference_agreement([], None, assignment('TTAGGG', 'species'),
                                [observation('TTAGGG', None, 'explicit_zero')])

    def test_orthogonal_annotation_preserves_classification_selection_and_pair_counts(self):
        refs = catalog_index([{'motif': 'TTTAGGG', 'evidence': 'experimental', 'source': 'fixture'},
                              {'motif': 'CTCGGTTATGGG', 'evidence': 'experimental', 'source': 'fixture'}])
        candidates = [assess_candidate({'motif': motif, 'reads': reads, 'repeat_bases': reads * 120,
                                       'full_repeat_reads': reads - 3, 'max_tract': 150},
                                      None, refs, {canonical('CTCGGTTATGGG')}, 150)
                      for motif, reads in (('TTTAGGG', 106), ('CTCGGTTATGGG', 16))]
        context = assignment()
        adequacy = {'status': 'unassessable_library_or_genome'}
        classification = classify(candidates, adequacy, reference_assignment=context)
        top = select_representative(candidates, classification, context)
        before = copy.deepcopy(candidates)
        result = reference_agreement(candidates, top, context, [observation('CTCGGTTATGGG', 16)])
        self.assertEqual(result['status'], 'observed_family_conflict')
        self.assertEqual(candidates, before)
        self.assertEqual(classify(candidates, adequacy, reference_assignment=context), classification)
        self.assertIs(select_representative(candidates, classification, context), top)

    def test_reviewed_taxonomic_identity_conflict_preserves_registered_taxon_and_not_rna(self):
        metadata = {'scientific_name': 'Tetrahymena thermophila', 'tax_id': '5911',
                    'library_strategy': 'WGS', 'library_source': 'GENOMIC'}
        reviewed = {'library_review_code': 'taxonomic_identity_conflict',
                    'detailed_evidence_path': '/fixture/identity_review.json',
                    'evidence_sha256': 'a' * 64}
        before = copy.deepcopy((metadata, reviewed))
        library = library_context(metadata, reviewed)
        self.assertTrue(library['limited'])
        self.assertTrue(library['identity_conflict_supported'])
        self.assertEqual(library['reasons'], ['taxonomic_identity_conflict'])
        self.assertFalse(library['RNA_component_supported'])
        self.assertFalse(library['metadata_contradiction_proves_RNA'])
        self.assertEqual(library['evidence'][0]['source'], reviewed['detailed_evidence_path'])
        self.assertEqual(library['evidence'][0]['sha256'], reviewed['evidence_sha256'])
        self.assertIn('no automatic taxon reassignment', library['interpretation'])
        self.assertNotIn('assigned_tax_id', library)
        self.assertEqual((metadata, reviewed), before)

    def test_unreviewed_low_mapping_or_metadata_code_does_not_establish_identity_conflict(self):
        library = library_context({'library_strategy': 'WGS', 'library_source': 'GENOMIC',
                                   'library_review_code': 'taxonomic_identity_conflict',
                                   'mapping_fraction': 0.00002})
        self.assertFalse(library['identity_conflict_supported'])
        self.assertFalse(library['limited'])
        self.assertFalse(library['RNA_component_supported'])
        library = library_context({}, {'library_review_code': 'metadata_conflict_unresolved'})
        self.assertTrue(library['limited'])
        self.assertFalse(library['identity_conflict_supported'])

    def test_reviewed_identity_scope_blocks_nomination_promotion_without_changing_scores(self):
        raw = {'motif': 'ACAT', 'reads': 120, 'repeat_bases': 14400,
               'full_repeat_reads': 60, 'max_tract': 150}
        targeted = {**raw, 'junction_evidence': {'full_repeat_records': 60, 'junction_reads': 24,
                    'left_forward_junction_reads': 12, 'right_reverse_junction_reads': 12,
                    'left_reverse_junction_reads': 0, 'right_forward_junction_reads': 0}}
        observed = assess_candidate(raw, targeted, {}, set(), 150)
        candidates = [observed]
        before = copy.deepcopy(candidates)
        adequacy = {'status': 'unassessable_library_or_genome'}
        self.assertEqual(classify(candidates, adequacy)[0], 'novel_geometry_nomination')
        library = library_context({}, {'library_review_code': 'taxonomic_identity_conflict'})
        classification = classify(candidates, adequacy, biological_scope=False, library=library)
        self.assertEqual(classification[0], 'taxonomic_scope_unresolved')
        top = select_representative(candidates, classification)
        agreement = reference_agreement(candidates, top, assignment('TTGGGG', 'species'),
                                         [observation('TTGGGG', 0)], biological_scope=False,
                                         library=library)
        self.assertEqual(agreement['status'], 'taxonomic_scope_unresolved')
        self.assertTrue(agreement['interpretation_limited'])
        self.assertEqual(agreement['expected_families'][0]['recovery'], 'explicit_zero')
        self.assertIsNone(agreement['main_trm_motif'])
        self.assertIs(top, observed)
        self.assertEqual(candidates, before)


if __name__ == '__main__':
    unittest.main()
