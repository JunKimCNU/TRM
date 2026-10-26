"""Clade-relative variation is orthogonal to catalogue status and TRM claims."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
from clade_variant_reporting import build_clade_variant_report
from evidence_batch import canonical
from prepare_taxon_priors import load_references


ROOT = Path(__file__).parents[1]


def metadata(**extra):
    row = {'sample_id': 'synthetic', 'status': 'resolved', 'scientific_name': 'Synthetic nematode',
           'tax_id': '99999999', 'lineage': 'Eukaryota;Metazoa;Nematoda',
           'lineage_tax_ids': '2759;33208;6231', 'library_source': 'GENOMIC'}
    row.update(extra)
    return row


def candidate(motif='TTAGGG', reads=30, **extra):
    row = {'motif': motif, 'reads': reads, 'reference_status': 'established',
           'known_detected': True, 'candidate_status': 'known_repeat_supported',
           'novel_candidate': False, 'geometry_nomination': False}
    row.update(extra)
    return row


def assignment(motifs='TTAGGC', scope='clade', **extra):
    row = {'reference_scope': scope, 'reference_status': 'supported', 'related_motifs': motifs,
           'repeat_architecture': 'exact_short_unit', 'reference_ids': 'synthetic_reference'}
    row.update(extra)
    return row


def catalogue(*motifs):
    return {canonical(motif): [{'motif': motif, 'evidence': 'experimental'}] for motif in motifs}


class CladeVariantReportingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.references = load_references(ROOT / 'data/taxon_motif_references.tsv')

    def build(self, candidates=None, **options):
        return build_clade_variant_report(
            [candidate()] if candidates is None else candidates,
            options.pop('metadata', metadata()), options.pop('assignment', assignment()),
            options.pop('taxon_reference_rows', self.references),
            global_catalogue=options.pop('global_catalogue', catalogue('TTAGGG', 'TTAGGC', 'TTAGAC')),
            **options)

    def test_globally_known_animal_motif_remains_nematode_clade_variant(self):
        report = self.build()
        self.assertEqual(report['status'], 'clade_variant_observed')
        row = report['families'][0]
        self.assertEqual(row['global_catalogue_status'], 'globally_established_family')
        self.assertEqual(row['clade_relationship'], 'differs_from_clade_default')
        self.assertEqual(report['default_families'], [canonical('TTAGGC')])
        self.assertTrue(report['review_required'])
        self.assertEqual(row['next_action'], 'validate_clade_variant_terminal_location')
        self.assertFalse(report['biological_novelty_claimed'])
        self.assertFalse(report['literature_first_report_claimed'])

    def test_rotation_reverse_complement_and_tandem_powers_are_same_default(self):
        for motif in ('TTAGGC', 'TAGGCT', 'GCCTAA', 'TTAGGCTTAGGC'):
            with self.subTest(motif=motif):
                report = self.build([candidate(motif)], representative={'motif': 'GCCTAA'})
                self.assertEqual(report['status'], 'clade_canonical_observed')
                self.assertTrue(report['families'][0]['is_representative'])
                self.assertFalse(report['review_required'])

    def test_published_species_exception_is_not_new_biology(self):
        report = self.build([candidate('TTAGAC')], assignment=assignment('TTAGAC', 'species'))
        row = report['families'][0]
        self.assertEqual(row['status'], 'published_clade_exception_observed')
        self.assertEqual(row['scoped_reference_relationship'], 'published_species_or_isolate_exception')
        self.assertTrue(row['published_scoped_exception'])
        self.assertEqual(report['published_exception_family_count'], 1)
        self.assertEqual(row['next_action'], 'validate_published_clade_exception_in_sample')
        self.assertFalse(row['biological_novelty_claimed'])

    def test_isolate_exception_and_genus_expectation_have_different_scope(self):
        isolated = self.build(assignment=assignment('TTAGGG', 'isolate'))['families'][0]
        genus = self.build(assignment=assignment('TTAGGG', 'genus'))['families'][0]
        self.assertTrue(isolated['published_scoped_exception'])
        self.assertFalse(genus['published_scoped_exception'])
        self.assertEqual(genus['scoped_reference_relationship'], 'supported_genus_expectation')

    def test_unknown_taxonomy_and_identity_conflict_cannot_be_inferred_from_repeat(self):
        unresolved = self.build(metadata=metadata(status='unresolved'))
        absent_lineage = self.build(metadata=metadata(lineage=''))
        conflict = self.build(library={'identity_conflict_supported': True, 'limited': True})
        for report in (unresolved, absent_lineage, conflict):
            self.assertEqual(report['status'], 'taxonomic_scope_unresolved')
            self.assertFalse(report['comparison_permitted'])
            self.assertFalse(report['review_required'])

    def test_unknown_dipteran_context_blocks_insect_and_animal_default(self):
        report = self.build(metadata=metadata(lineage='Eukaryota;Metazoa;Arthropoda;Insecta;Diptera',
                                             lineage_tax_ids='2759;33208;6656;50557;7147'))
        self.assertEqual(report['status'], 'nearest_clade_has_no_supported_default')
        self.assertEqual(report['clade_context']['reference_status'], 'unknown')
        self.assertEqual(report['default_families'], [])
        self.assertFalse(report['review_required'])

    def test_unknown_fungi_do_not_inherit_a_universal_fungal_motif(self):
        report = self.build(metadata=metadata(lineage='Eukaryota;Fungi', lineage_tax_ids='2759;4751'))
        self.assertEqual(report['status'], 'nearest_clade_has_no_supported_default')
        self.assertFalse(report['variant_observed'])

    def test_more_specific_unknown_candidate_and_complex_contexts_block_interpretation(self):
        cases = [assignment('', 'species', reference_status=status)
                 for status in ('unknown', 'out_of_scope', 'candidate')]
        cases.append(assignment('TTAGGG', 'species', repeat_architecture='exact_short_unit,irregular'))
        for scoped in cases:
            with self.subTest(scoped=scoped):
                report = self.build(assignment=scoped)
                self.assertEqual(report['status'], 'specific_reference_blocks_clade_interpretation')
                self.assertFalse(report['comparison_permitted'])

    def test_library_limited_observations_remain_visible_without_actionable_attribution(self):
        report = self.build(library={'limited': True, 'reasons': ['nonuniform_sampling']})
        self.assertEqual(report['status'], 'library_context_limited')
        self.assertEqual(len(report['families']), 1)
        self.assertEqual(report['families'][0]['reads'], 30)
        self.assertTrue(report['families'][0]['comparison_limited_by_library'])
        self.assertIsNone(report['variant_observed'])
        self.assertFalse(report['review_required'])

    def test_known_positive_does_not_require_geometry_to_be_listed_or_reviewed(self):
        report = self.build([candidate(geometry={'status': 'insufficient_junctions'})])
        row = report['families'][0]
        self.assertTrue(row['known_repeat_pattern_supported'])
        self.assertTrue(row['review_required'])
        self.assertFalse(row['terminal_location_confirmed'])

    def test_trace_and_sequence_only_observation_do_not_enter_actionable_queue(self):
        trace = self.build([candidate(reads=2, known_detected=False)])
        sequence_only = self.build([candidate(reads=30, known_detected=False)])
        self.assertEqual(trace['status'], 'clade_variant_trace_only')
        self.assertEqual(trace['variant_trace_family_count'], 1)
        self.assertTrue(trace['families'][0]['sequence_only_observation'])
        self.assertFalse(trace['review_required'])
        self.assertEqual(sequence_only['status'], 'clade_variant_observed')
        self.assertFalse(sequence_only['families'][0]['known_repeat_pattern_supported'])
        self.assertTrue(sequence_only['families'][0]['sequence_only_observation'])
        self.assertFalse(sequence_only['review_required'])

    def test_secondary_variant_is_not_hidden_when_canonical_is_representative(self):
        report = self.build([candidate('TTAGGC', reads=1000), candidate()], representative={'motif': 'TTAGGC'})
        self.assertEqual(report['status'], 'clade_canonical_and_variant_observed')
        self.assertEqual(len(report['families']), 2)
        self.assertEqual(len(report['review_items']), 1)
        self.assertEqual(report['review_items'][0]['canonical_family'], canonical('TTAGGG'))

    def test_sample_specific_extension_does_not_fabricate_global_catalogue_membership(self):
        report = self.build([candidate('TTAGAC')], global_catalogue=catalogue('TTAGGG'),
                            assignment=assignment('TTAGAC', 'species'))
        row = report['families'][0]
        self.assertEqual(row['global_catalogue_status'], 'not_in_global_catalogue')
        self.assertEqual(row['sample_reference_status'], 'established')
        self.assertTrue(row['published_scoped_exception'])

    def test_uncatalogued_existing_nomination_is_visible_without_promoting_other_microsatellites(self):
        nominated = candidate('TTAGTC', reference_status='unreferenced', known_detected=False, novel_candidate=True)
        tail = candidate('AAAAT', reads=10000, reference_status='unreferenced', known_detected=False)
        report = self.build([nominated, tail])
        self.assertEqual(len(report['families']), 1)
        self.assertEqual(report['families'][0]['canonical_family'], canonical('TTAGTC'))
        self.assertTrue(report['review_required'])
        self.assertEqual(report['unqualified_observed_families_not_copied'], 1)
        self.assertFalse(report['changes_candidate_counts_or_ranking'])

    def test_measurement_disagreement_is_review_not_positive_variant(self):
        report = self.build([candidate(reads=0, original_discovery_reads=100, measurement_disagreement=True,
                                     known_detected=False)])
        self.assertEqual(report['status'], 'measurement_disagreement')
        self.assertFalse(report['variant_observed'])
        self.assertEqual(report['variant_observed_family_count'], 0)
        row = report['families'][0]
        self.assertEqual((row['reads'], row['original_discovery_reads']), (0, 100))
        self.assertEqual(row['next_action'], 'resolve_variant_measurement_disagreement')

    def test_supplied_taxids_are_not_overridden_with_matching_lineage_names(self):
        report = self.build(metadata=metadata(lineage_tax_ids='2759;99999998'))
        self.assertEqual(report['status'], 'no_curated_clade_default')
        self.assertFalse(report['comparison_permitted'])

    def test_observation_does_not_alter_inputs_or_existing_decisions(self):
        inputs = [[candidate(), candidate('TTAGGC')], metadata(), assignment(), copy.deepcopy(self.references)]
        before = copy.deepcopy(inputs)
        report = build_clade_variant_report(*inputs, global_catalogue=catalogue('TTAGGG'))
        self.assertEqual(inputs, before)
        self.assertFalse(report['changes_classification'])
        self.assertFalse(report['changes_adequacy'])

    def test_incomplete_analysis_does_not_create_actionable_review(self):
        report = self.build(analysis_complete=False)
        self.assertEqual(report['status'], 'analysis_incomplete')
        self.assertEqual(report['families'], [])
        self.assertFalse(report['review_required'])

    def test_duplicate_measurements_are_never_summed_into_stronger_evidence(self):
        with self.assertRaisesRegex(ValueError, 'one assessed row per primitive family'):
            self.build([candidate(), candidate('CCCTAA')])


if __name__ == '__main__':
    unittest.main()
