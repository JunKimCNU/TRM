"""New quality diagnostics remain observations, with missingness and source binding."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from quality_observation import describe_quality, describe_candidate_quality


def histogram(bases, score=40):
    bins = [0] * 94
    bins[score] = bases
    return {'phred33_counts': bins, 'bases': bases,
            'below_q20_bases': bases if score < 20 else 0,
            'below_q30_bases': bases if score < 30 else 0}


def fixture():
    records = {'observed_records': 2, 'quality_records': 1,
               'fasta_records_without_quality': 1, 'fastq_records_above_retention_limit': 0,
               'whole_read': histogram(150)}
    report = {'summary': {'records': 2, 'bases': 300}, 'quality_diagnostics': {
        'schema_version': 1, 'encoding': 'phred33',
        'encoding_provenance': 'explicit_user_declaration_not_inferred',
        'input_scope': 'complete_inputs', 'input_records': records,
        'filtering_applied': False, 'variant_discriminating_bases_assessed': False}}
    candidate = {'reads': 2, 'repeat_bases': 180, 'quality_diagnostics': {
        'supporting_records': records, 'union_tract_bases': 180,
        'union_tract_bases_without_quality': 90, 'union_tract': histogram(90)}}
    return report, candidate


class QualityObservationTests(unittest.TestCase):
    def test_missing_old_or_targeted_quality_never_borrows_other_measurements(self):
        report, candidate = fixture()
        absent = describe_quality({}, 'old')
        self.assertEqual(absent['status'], 'not_measured')
        self.assertIsNone(describe_candidate_quality(candidate, absent))
        context = describe_quality(report, 'new')
        self.assertIsNone(describe_candidate_quality({'reads': 2}, context))

    def test_partial_quality_preserves_denominators_and_does_not_filter(self):
        report, candidate = fixture()
        before = copy.deepcopy((report, candidate))
        context = describe_quality(report, 'source')
        value = describe_candidate_quality(candidate, context)
        self.assertEqual(context['status'], 'partially_measured')
        self.assertEqual(value['diagnostics']['union_tract']['bases'], 90)
        self.assertEqual(value['diagnostics']['union_tract_bases_without_quality'], 90)
        self.assertFalse(value['filtering_applied'])
        self.assertFalse(value['variant_discriminating_bases_assessed'])
        self.assertEqual((report, candidate), before)

    def test_conflicting_histogram_counts_encoding_or_tracts_are_rejected(self):
        report, candidate = fixture()
        bad = copy.deepcopy(report)
        bad['quality_diagnostics']['input_records']['whole_read']['phred33_counts'][0] += 1
        with self.assertRaisesRegex(ValueError, 'histogram totals'):
            describe_quality(bad, 'source')
        bad = copy.deepcopy(report)
        bad['quality_diagnostics']['encoding'] = 'inferred'
        with self.assertRaisesRegex(ValueError, 'contract'):
            describe_quality(bad, 'source')
        candidate['repeat_bases'] += 1
        with self.assertRaisesRegex(ValueError, 'tract counts'):
            describe_candidate_quality(candidate, describe_quality(report, 'source'))

    def test_complete_quality_must_cover_every_analyzed_base(self):
        report, _ = fixture()
        records = report['quality_diagnostics']['input_records']
        records.update(quality_records=2, fasta_records_without_quality=0)
        for bases in (2, 150, 299):
            with self.subTest(bases=bases):
                records['whole_read'] = histogram(bases)
                with self.assertRaisesRegex(ValueError, 'Complete quality bases'):
                    describe_quality(report, 'source')
        records['whole_read'] = histogram(300)
        self.assertEqual(describe_quality(report, 'source')['status'], 'measured')

    def test_quality_bearing_records_require_nonempty_whole_read_counts(self):
        for bases in (0, 1):
            with self.subTest(bases=bases):
                report, candidate = fixture()
                context = describe_quality(report, 'source')
                records = report['quality_diagnostics']['input_records']
                records.update(quality_records=2, fasta_records_without_quality=0,
                               whole_read=histogram(bases))
                with self.assertRaisesRegex(ValueError, 'fewer than quality-bearing'):
                    describe_quality(report, 'source')
                with self.assertRaisesRegex(ValueError, 'fewer than quality-bearing'):
                    describe_candidate_quality(candidate, context)

    def test_histogram_totals_reject_booleans_and_floats_even_when_equal(self):
        for key in ('bases', 'below_q20_bases', 'below_q30_bases'):
            for replacement in (False, 0.0):
                with self.subTest(key=key, replacement=repr(replacement)):
                    report, candidate = fixture()
                    context = describe_quality(report, 'source')
                    for records in (report['quality_diagnostics']['input_records'],
                                    candidate['quality_diagnostics']['supporting_records']):
                        records.update(quality_records=0, fasta_records_without_quality=2,
                                       whole_read=histogram(0))
                        records['whole_read'][key] = replacement
                    with self.assertRaisesRegex(ValueError, 'Invalid quality count'):
                        describe_quality(report, 'source')
                    with self.assertRaisesRegex(ValueError, 'Invalid quality count'):
                        describe_candidate_quality(candidate, context)

    def test_sequence_summary_counts_are_required_nonnegative_integers(self):
        for key in ('records', 'bases'):
            for bad in (None, False, 2.0, -1):
                with self.subTest(key=key, bad=repr(bad)):
                    report, _ = fixture()
                    report['summary'][key] = bad
                    with self.assertRaisesRegex(ValueError, 'Invalid quality count'):
                        describe_quality(report, 'source')

    def test_unavailable_quality_and_zero_hit_candidates_remain_valid(self):
        for missing_key in ('fasta_records_without_quality', 'fastq_records_above_retention_limit'):
            with self.subTest(missing_key=missing_key):
                report, _ = fixture()
                records = report['quality_diagnostics']['input_records']
                records.update(quality_records=0, fasta_records_without_quality=0,
                               fastq_records_above_retention_limit=0, whole_read=histogram(0))
                records[missing_key] = 2
                context = describe_quality(report, 'source')
                self.assertEqual(context['status'], 'unavailable_in_analyzed_records')
                empty = dict(records, observed_records=0)
                empty[missing_key] = 0
                candidate = {'reads': 0, 'repeat_bases': 0, 'quality_diagnostics': {
                    'supporting_records': empty, 'union_tract_bases': 0,
                    'union_tract_bases_without_quality': 0, 'union_tract': histogram(0)}}
                self.assertEqual(describe_candidate_quality(candidate, context)['status'], 'unavailable')


if __name__ == '__main__':
    unittest.main()
