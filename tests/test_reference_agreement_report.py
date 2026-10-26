"""Report joins preserve native observations while exposing reference conflicts."""
import copy
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parent))
import test_report_assessment as fixtures


class ReferenceAgreementReportTests(unittest.TestCase):
    def setUp(self):
        # Compose the saved-artifact fixture; do not inherit and rerun its suite.
        self.fixture = fixtures.AssessmentReportTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.sample = 'SRRknown'
        self.source = self.fixture.source / self.sample / 'result.json'
        self.diagnostics = self.fixture.root / 'expected_diagnostics'
        self.meta = next(row for row in fixtures.read_tsv(self.fixture.metadata)
                         if row['sample_id'] == self.sample)
        self.assignment = next(row for row in fixtures.read_tsv(self.fixture.assignments)
                               if row['sample_id'] == self.sample)
        source = json.loads(self.source.read_text())
        evidence = json.loads((self.fixture.evidence / self.sample / 'evidence.json').read_text())
        source.update(version='0.8.0', inputs=[{'path': '/fixture/reads.fastq.gz', 'bytes': 123,
                                              'modified_unix_seconds': 1}],
                      evidence_qc=evidence['evidence_qc'], junction_options={'end_slack': 3, 'flank_bases': 20})
        source['config'].update(min_identity=.85, min_tract=500, end_window=50)
        self.raw = source
        self.template = evidence['candidate_evidence'][0]

    def configure(self, expected, candidates, *, scope='species', reference_status='supported'):
        self.assignment.update(related_motifs=expected, reference_scope=scope,
                               reference_status=reference_status, repeat_architecture='exact_short_unit')
        fixtures.write_tsv(self.fixture.metadata, [self.meta])
        fixtures.write_tsv(self.fixture.assignments, [self.assignment])
        fixtures.write_tsv(self.fixture.library, [row for row in fixtures.read_tsv(self.fixture.library)
                                                  if row['sample_id'] == self.sample])
        motifs = list(dict.fromkeys([m for m in expected.split(',') if m] + ['TTTAGGG']))
        fixtures.write_tsv(self.fixture.catalog, [{'motif': motif, 'lineage': 'Fixture species',
                               'evidence': 'experimental', 'source': 'https://example.test/reference'}
                              for motif in motifs])
        self.raw['candidates'] = []
        for motif, reads, full in candidates:
            row = copy.deepcopy(self.template)
            row.update(motif=fixtures.canonical(motif), reads=reads, full_repeat_reads=full,
                       repeat_bases=reads * 120, max_tract=150)
            row['junction_evidence'].update(full_repeat_records=full, junction_reads=6,
                left_forward_junction_reads=3, right_reverse_junction_reads=3,
                left_reverse_junction_reads=0, right_forward_junction_reads=0)
            self.raw['candidates'].append(row)
        fixtures.write_json(self.source, self.raw)

    def save_diagnostic(self, measured):
        rows = []
        for motif, reads, full in measured:
            native = next((row for row in self.raw['candidates']
                           if row['motif'] == fixtures.canonical(motif)), None)
            rows.append({'motif': fixtures.canonical(motif), 'reads': reads, 'full_repeat_reads': full,
                         'repeat_bases': native['repeat_bases'] if native else reads * 120,
                         'max_tract': native['max_tract'] if native else 150 if reads else 0})
        document = {**copy.deepcopy(self.raw), 'analysis_kind': 'targeted_evidence',
                    'candidate_evidence': rows}
        path = self.diagnostics / self.sample / 'evidence.json'
        fixtures.write_json(path, document)
        fixtures.write_json(path.parent / 'status.json', {
            'sample_id': self.sample, 'status': 'complete', 'exit_code': 0,
            'identity': {'source_sha256': fixtures.digest(self.source), 'motifs': [row[0] for row in measured]},
            'result_sha256': fixtures.digest(path)})

    def report(self, diagnostic=False, *, success=True, extra=()):
        extra = (('--expected-evidence', str(self.diagnostics)) if diagnostic else ()) + tuple(extra)
        completed = self.fixture.run_report(native=True, extra=extra, success=success)
        return self.fixture.document(self.sample) if success else completed

    def test_s2_conflict_headline_does_not_change_native_counts_ranking_or_adequacy(self):
        self.configure('CTCGGTTATGGG', [('TTTAGGG', 106, 60), ('CTCGGTTATGGG', 16, 6)], scope='genus')
        source_hash = fixtures.digest(self.source)
        baseline = self.report()
        self.save_diagnostic([('CTCGGTTATGGG', 16, 6), ('TTTAGGG', 106, 60)])
        document = self.report(diagnostic=True)
        self.assertEqual(document['reference_agreement']['status'], 'observed_family_conflict')
        self.assertEqual(document['detected_motif'], 'TTTAGGG')
        self.assertEqual(document['expected_family_observations'][0]['reads'], 16)
        self.assertEqual(document['top_candidate']['reads'], 106)
        self.assertEqual(document['candidates'], baseline['candidates'])
        self.assertEqual(document['classification'], baseline['classification'])
        self.assertEqual(document['data_adequacy'], baseline['data_adequacy'])
        self.assertEqual(fixtures.digest(self.source), source_hash)
        self.assertEqual(document['next_action'], 'review_reference_and_representative')
        self.assertIn('main TRM 해석 보류', document['trm_interpretation'])
        html = (self.fixture.output / 'samples' / f'{self.sample}.html').read_text()
        # Observation reporting separates the primary decision, reference conflict,
        # and policy representative; the conflict must remain visible in its row.
        self.assertIn('<p class="result">' + document['decision_summary']['display_label'], html)
        self.assertIn('<p>계통 비교: ' + document['trm_interpretation'] + '</p>', html)
        self.assertIn('정책상 대표 반복: <b>TTTAGGG</b>', html)
        self.assertIn('Main TRM: <b>' + document['decision_summary']['main_trm_display'], html)
        summary = fixtures.read_tsv(self.fixture.output / 'samples.tsv')[0]
        self.assertEqual(summary['reference_agreement_status'], 'observed_family_conflict')
        self.assertEqual(summary['supporting_reads'], '106')

    def test_tetrahymena_expected_zero_is_provenance_join_only(self):
        self.configure('TTGGGG', [('ACAT', 40, 10)])
        baseline = self.report()
        self.assertEqual(baseline['reference_agreement']['status'], 'expected_reference_unmeasured')
        self.assertIsNone(baseline['expected_family_observations'][0]['reads'])
        self.save_diagnostic([('TTGGGG', 0, 0)])
        document = self.report(diagnostic=True)
        self.assertEqual(document['reference_agreement']['status'], 'expected_reference_not_observed')
        self.assertEqual(document['expected_family_observations'][0]['recovery'], 'explicit_zero')
        self.assertEqual(document['expected_family_observations'][0]['reads'], 0)
        self.assertEqual(document['candidates'], baseline['candidates'])
        self.assertEqual(document['classification'], baseline['classification'])
        self.assertEqual(document['data_adequacy'], baseline['data_adequacy'])
        self.assertNotIn(fixtures.canonical('TTGGGG'), [row['canonical_family'] for row in document['candidates']])
        self.assertIsNone(document['reference_agreement']['main_trm_motif'])
        provenance = document['provenance']['expected_family_diagnostic']
        self.assertEqual(provenance['source_sha256'], fixtures.digest(self.source))
        self.assertEqual(provenance['result_sha256'], fixtures.digest(provenance['result_file']))
        self.assertEqual(provenance['used_expected_families'], [fixtures.canonical('TTGGGG')])

    def test_unrelated_target_cannot_create_expected_zero_or_new_candidate(self):
        self.configure('', [('ACAT', 40, 10)], reference_status='not_curated')
        baseline = self.report()
        self.save_diagnostic([('TTGGGG', 0, 0)])
        document = self.report(diagnostic=True)
        self.assertEqual(document['expected_family_observations'], [])
        self.assertEqual(document['reference_agreement']['status'], 'no_supported_reference')
        self.assertFalse(document['reference_agreement']['expected_reference_not_observed'])
        self.assertEqual(document['candidates'], baseline['candidates'])
        self.assertEqual(document['provenance']['expected_family_diagnostic']['used_expected_families'], [])

    def test_stale_diagnostic_source_is_rejected_by_report(self):
        self.configure('TTGGGG', [('ACAT', 40, 10)])
        self.save_diagnostic([('TTGGGG', 0, 0)])
        path = self.diagnostics / self.sample / 'status.json'
        status = json.loads(path.read_text())
        status['identity']['source_sha256'] = 'stale-source'
        fixtures.write_json(path, status)
        completed = self.report(diagnostic=True, success=False)
        self.assertIn('source hash differs', completed.stderr)

    def test_library_limited_known_repeat_recommendation_preserves_limit(self):
        self.configure('TTTAGGG', [('TTTAGGG', 106, 60)])
        self.raw['candidates'][0]['junction_evidence'].update(
            junction_reads=24, left_forward_junction_reads=12, right_reverse_junction_reads=12)
        fixtures.write_json(self.source, self.raw)
        library = fixtures.read_tsv(self.fixture.library)
        library[0]['library_strategy'] = 'RNA-Seq'
        fixtures.write_tsv(self.fixture.library, library)
        document = self.report()
        self.assertEqual(document['reference_agreement']['status'], 'library_context_limited')
        self.assertEqual(document['observation_next_action'], 'canonical_screen_complete')
        self.assertEqual(document['next_action'], 'resolve_library_context')
        self.assertNotIn('새 motif 검증 대상에서는 제외', document['recommendation'])

    def save_identity_review(self):
        path = self.fixture.root / 'identity_review.json'
        fixtures.write_json(path, {'sample_id': self.sample,
                                  'finding': 'Fixture positive reference-consistency conflict; no taxon reassignment.'})
        row = {'sample_id': self.sample, 'scientific_name': self.meta['scientific_name'],
               'tax_id': self.meta['tax_id'], 'library_review_code': 'taxonomic_identity_conflict',
               'detailed_evidence_path': path.name, 'evidence_sha256': fixtures.digest(path),
               'source_result_sha256': fixtures.digest(self.source),
               'previous_full_read_records': self.raw['summary']['records']}
        table = self.fixture.root / 'identity_review.tsv'
        fixtures.write_tsv(table, [row])
        return table, row

    def test_reviewed_identity_conflict_limits_attribution_but_preserves_metadata_and_observations(self):
        self.configure('TTGGGG', [('ACAT', 40, 10)])
        self.raw['candidates'][0]['junction_evidence'].update(
            junction_reads=24, left_forward_junction_reads=12, right_reverse_junction_reads=12)
        fixtures.write_json(self.source, self.raw)
        self.save_diagnostic([('TTGGGG', 0, 0)])
        baseline = self.report(diagnostic=True)
        self.assertEqual(baseline['classification'], 'novel_geometry_nomination')
        table, row = self.save_identity_review()
        document = self.report(diagnostic=True, extra=('--library-evidence', str(table)))
        self.assertEqual(document['classification'], 'taxonomic_scope_unresolved')
        self.assertEqual(document['reference_agreement']['status'], 'taxonomic_scope_unresolved')
        self.assertTrue(document['library_context']['identity_conflict_supported'])
        self.assertFalse(document['library_context']['RNA_component_supported'])
        self.assertFalse(document['library_context']['metadata_contradiction_proves_RNA'])
        self.assertEqual(document['metadata'], baseline['metadata'])
        self.assertEqual(document['tax_id'], self.meta['tax_id'])
        self.assertEqual(document['scientific_name'], self.meta['scientific_name'])
        self.assertEqual(document['reference_assignment'], baseline['reference_assignment'])
        self.assertEqual(document['candidates'], baseline['candidates'])
        self.assertEqual(document['top_candidate'], baseline['top_candidate'])
        self.assertEqual(document['data_adequacy'], baseline['data_adequacy'])
        self.assertEqual(document['expected_family_observations'][0]['recovery'], 'explicit_zero')
        self.assertEqual(document['reviewed_library_evidence']['source_result_sha256'], row['source_result_sha256'])
        self.assertEqual(document['next_action'], 'check_host_or_mixed_sample')
        self.assertIsNone(document['reference_agreement']['main_trm_motif'])

    def test_identity_review_stale_source_hash_is_rejected(self):
        self.configure('TTGGGG', [('ACAT', 40, 10)])
        table, row = self.save_identity_review()
        fixtures.write_tsv(table, [{**row, 'source_result_sha256': '0' * 64}])
        completed = self.report(success=False, extra=('--library-evidence', str(table)))
        self.assertIn('reviewed identity conflict source result hash differs', completed.stderr)

    def test_identity_review_requires_exact_metadata_and_evidence_hash(self):
        self.configure('TTGGGG', [('ACAT', 40, 10)])
        table, row = self.save_identity_review()
        for change, error in (({'tax_id': '999'}, 'tax_id differs'),
                              ({'scientific_name': 'Different taxon'}, 'scientific_name differs'),
                              ({'evidence_sha256': '0' * 64}, 'provenance hash differs'),
                              ({'source_result_sha256': ''}, 'requires source_result_sha256'),
                              ({'evidence_sha256': ''}, 'requires evidence_sha256')):
            with self.subTest(change=change):
                fixtures.write_tsv(table, [{**row, **change}])
                completed = self.report(success=False, extra=('--library-evidence', str(table)))
                self.assertIn(error, completed.stderr)

    def test_rejected_nearby_family_is_reported_without_inventing_expected_zero(self):
        expected = 'AAGAGAGAGAGAGAGAAGGGGGG'
        observed = 'AAGAGAGAGAGAGAGAGGGGGGG'
        self.configure(expected, [(observed, 68, 17)])
        self.raw['rejected_candidates'] = self.raw.pop('candidates')
        self.raw['candidates'] = []
        fixtures.write_json(self.source, self.raw)
        before_hash = fixtures.digest(self.source)
        document = self.report()
        view = document['family_observation']
        self.assertEqual(view['status'], 'nearby_family_observed')
        family = view['families'][0]
        self.assertEqual(family['exact_measurement_status'], 'not_reported')
        neighbor = next(row for row in family['nearby_families_all']
                        if row['canonical_family'] == fixtures.canonical(observed))
        self.assertEqual(neighbor['edit_distance'], 1)
        self.assertEqual(neighbor['original']['reads'], 68)
        self.assertFalse(document['reporting_features']['family_observation_changes_classification_or_adequacy'])
        self.assertEqual(document['reference_agreement']['main_trm_status'], 'not_confirmed')
        self.assertEqual(fixtures.digest(self.source), before_hash)
        table = fixtures.read_tsv(self.fixture.output / 'expected_family_observations.tsv')[0]
        self.assertEqual(table['current_reads'], '')
        self.assertEqual(table['original_reads'], '')
        neighbors = fixtures.read_tsv(self.fixture.output / 'nearby_families.tsv')
        self.assertTrue(any(row['observed_family'] == fixtures.canonical(observed)
                            and row['original_reads'] == '68' for row in neighbors))
        html = (self.fixture.output / 'samples' / f'{self.sample}.html').read_text()
        self.assertIn('기대 family와 근접 반복 관측', html)
        self.assertIn('nearby_families.tsv', html)

    def test_explicit_expected_zero_keeps_nearby_observation_and_policy(self):
        self.configure('TTAGGG', [('TTTAGGG', 40, 10)])
        baseline = self.report()
        self.save_diagnostic([('TTAGGG', 0, 0)])
        document = self.report(diagnostic=True)
        self.assertEqual(document['family_observation']['status'], 'nearby_family_observed')
        family = document['family_observation']['families'][0]
        self.assertEqual(family['exact_measurement_status'], 'explicit_zero')
        self.assertEqual(family['exact_current']['reads'], 0)
        for key in ('classification', 'top_candidate', 'candidates', 'data_adequacy'):
            self.assertEqual(document[key], baseline[key])
        self.assertEqual(document['schema_version'], 22)
        self.assertEqual(document['reporting_version'], '1.13.0')

    def test_unrelated_positive_candidates_are_not_reported_as_no_repeats(self):
        self.configure('TTAGGG', [('ACAT', 40, 10)])
        document = self.report()
        self.assertEqual(document['family_observation']['status'], 'other_repeat_families_observed')
        self.assertEqual(document['family_observation']['families'][0]['exact_measurement_status'], 'not_reported')
        self.assertEqual(document['family_observation']['families'][0]['nearby_family_count'], 0)


if __name__ == '__main__':
    unittest.main()
