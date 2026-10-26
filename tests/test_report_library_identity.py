"""Guard interpretation boundaries in the one-off library review renderer."""
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('identity_report', ROOT / 'scripts/report_library_identity.py')
report = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(report)


def sequence(role='RNA_PolyA_metadata_conflict', sid='ERR13319074'):
    return {'sample_id': sid, 'role': role, 'full_SL_anchored_records': 0,
            'ligation_read_fraction': 0.002, 'ligation_internal_enrichment_vs_iid': 0.9,
            'ligation_internal_enrichment_vs_markov1': 0.7}


class LibraryIdentityReportTests(unittest.TestCase):
    def test_rna_metadata_and_leader_counts_do_not_alone_assign_identity(self):
        for sid in ('ERR13319074', 'ERR13319101', 'ERR13623718', 'ERR13623736'):
            for hits in (0, 1000):
                row = sequence(sid=sid)
                row['full_SL_anchored_records'] = hits
                self.assertEqual(report.interpretation(row)['library_review_code'],
                                 'metadata_conflict_unresolved')

    def test_reviewed_ju765_alignment_supports_component_without_purity(self):
        alignment = {'counts': {'strong_spliced_primary_records': 25,
                                'high_mapq_primary_mapped_records': 100},
                     'high_mapq_strong_spliced_fraction': 0.25,
                     'high_mapq_aligned_reference_ACGT_fractions': {'exonic': 0.92}}
        review = {'review_label': 'RNA_supported', 'library_purity_estimated': False}
        row = sequence(sid='ERR13623736')
        verdict = report.interpretation(row, alignment, review)
        self.assertEqual(verdict['library_review_code'], 'RNA_component_supported')
        self.assertIn('혼합비', verdict['next_step'])
        self.assertIn('미확정', verdict['next_step'])
        alignment['high_mapq_strong_spliced_fraction'] = 0.50
        with self.assertRaisesRegex(ValueError, 'denominator'):
            report.interpretation(row, alignment, review)
        alignment['high_mapq_strong_spliced_fraction'] = 0.25
        review['library_purity_estimated'] = True
        with self.assertRaisesRegex(ValueError, 'purity'):
            report.interpretation(row, alignment, review)

    def test_hic_annotation_without_enrichment_stays_annotation_only(self):
        row = sequence(role='HiC_protocol_annotation')
        self.assertEqual(report.interpretation(row)['library_review_code'], 'HiC_annotation_only')
        row['ligation_internal_enrichment_vs_markov1'] = 15
        row['ligation_internal_enrichment_vs_iid'] = 20
        verdict = report.interpretation(row)
        self.assertEqual(verdict['library_review_code'], 'HiC_supported')
        self.assertIn('균일 WGS가 아니', verdict['TRM_interpretation_impact'])

    def test_genomic_context_never_claims_purity_or_adequacy(self):
        verdict = report.interpretation(sequence(role='ordinary_genomic_comparison'))
        self.assertEqual(verdict['library_review_code'], 'WGS_context_compatible')
        self.assertIn('증명되지는 않습니다', verdict['TRM_interpretation_impact'])


if __name__ == '__main__':
    unittest.main()
