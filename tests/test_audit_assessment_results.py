"""Independent artifact checker catches altered outputs and stale evidence."""
import importlib.util
import copy
import json
from pathlib import Path
import subprocess
import sys
import unittest

import test_report_assessment as fixtures


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('assessment_artifact_audit', ROOT/'scripts/audit_assessment_results.py')
audit_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit_module)


class AssessmentArtifactAuditTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.AssessmentReportTests('test_end_to_end_known_novel_unassessable_and_failed_are_distinct')
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        f = self.fixture
        binary = f.evidence/'trm-evidence'
        binary.write_bytes(b'Fixture binary artifact; never executed')
        original_hashes = []
        for source in sorted(f.source.glob('*/result.json')):
            sid = source.parent.name
            raw = json.loads(source.read_text())
            original = f.root/'original'/sid/'result.json'
            fixtures.write_json(original, raw)
            state = {'status': 'complete', 'result_sha256': fixtures.digest(original)}
            fixtures.write_json(original.parent/'state.json', state)
            fixtures.write_json(source.parent/'state.json', {'status': 'complete', 'result_sha256': fixtures.digest(source)})
            fixtures.write_json(source.parent/'provenance.json', {
                'source_result_json': str(original), 'original_result_sha256': fixtures.digest(original),
                'source_state': state, 'result_sha256': fixtures.digest(source)})
            status_path = f.evidence/sid/'status.json'
            status = json.loads(status_path.read_text())
            status['identity']['binary_sha256'] = fixtures.digest(binary)
            status.update(records=raw['summary']['records'], bases=raw['summary']['bases'])
            fixtures.write_json(status_path, status)
            original_hashes.append({'sample_id': sid, 'path': str(source), 'sha256': fixtures.digest(source)})
        decision = f.root/'decision.json'
        fixtures.write_json(decision, {'input_result_sha256': original_hashes})
        pilot = f.root/'pilot.json'
        fixtures.write_json(pilot, {'samples': [{'sample_id': 'SRRknown', 'detected_motif': 'TTAGGC'}]})
        fixtures.write_json(f.evidence/'metrics.json', {'sample_count': 4, 'completed': 3,
            'source_failed': 1, 'failed': 0, 'records': 30_000_000, 'bases': 4_500_000_000})
        f.run_report()
        self.args = audit_module.parse_args([
            '--assessment', str(f.output), '--comparison', str(f.comparison), '--source', str(f.source),
            '--evidence', str(f.evidence), '--metadata', str(f.metadata), '--assignments', str(f.assignments),
            '--catalog', str(f.catalog), '--taxon-references', str(ROOT/'data/taxon_motif_references.tsv'),
            '--decision-audit', str(decision), '--pilot-check', str(pilot),
            '--expected-samples', '4', '--expected-source-failures', '1',
            '--expected-pilot-samples', '1', '--expected-mhapla-samples', '0'])

    def failures(self):
        return [row for row in audit_module.run_audit(self.args)['checks'] if row['status'] == 'failed']

    def evidence_update(self, sample='SRRknown', reads=40, full=3):
        f = self.fixture
        path = f.evidence/sample/'evidence.json'
        evidence = json.loads(path.read_text())
        candidate = evidence['candidate_evidence'][0]
        candidate['reads'] = reads
        if reads == 0:
            candidate.update(repeat_bases=0, max_tract=0)
        candidate['junction_evidence'].update(full_repeat_records=full, junction_reads=0,
            left_forward_junction_reads=0, left_reverse_junction_reads=0,
            right_forward_junction_reads=0, right_reverse_junction_reads=0)
        fixtures.write_json(path, evidence)
        status_path = path.parent/'status.json'
        status = json.loads(status_path.read_text())
        status['result_sha256'] = fixtures.digest(path)
        fixtures.write_json(status_path, status)
        f.run_report()
        return f.document(sample)

    def observation_check(self, document, sample='SRRknown'):
        f = self.fixture
        return audit_module.check_observation_contract(document,
            json.loads((f.source/sample/'result.json').read_text()),
            json.loads((f.evidence/sample/'evidence.json').read_text()))

    def test_complete_fixture_and_independent_canonicalization_pass(self):
        self.assertEqual(self.failures(), [])
        self.assertEqual(audit_module.canonical('TTAGACTTAGAC'), 'AAGTCT')
        self.assertEqual(audit_module.canonical('TTAGGT'), 'AAACCT')

    def test_missing_direction_is_blank_but_measured_zero_is_not(self):
        f = self.fixture
        path = f.evidence/'SRRknown/evidence.json'
        evidence = json.loads(path.read_text())
        zero_geometry = {key: 0 for key in evidence['candidate_evidence'][0]['junction_evidence']}
        zero_geometry['full_repeat_records'] = 10
        for measured in (False, True):
            with self.subTest(measured=measured):
                candidate = evidence['candidate_evidence'][0]
                if measured:
                    candidate['junction_evidence'] = zero_geometry
                else:
                    candidate.pop('junction_evidence')
                fixtures.write_json(path, evidence)
                status_path = path.parent/'status.json'
                status = json.loads(status_path.read_text())
                status['result_sha256'] = fixtures.digest(path)
                fixtures.write_json(status_path, status)
                f.run_report()
                self.assertEqual(self.failures(), [])
                table = f.output/'candidates.tsv'
                rows = fixtures.read_tsv(table)
                row = next(row for row in rows if row['sample_id'] == 'SRRknown')
                self.assertEqual(row['left_consistent_junctions'], '0' if measured else '')
                row['left_consistent_junctions'] = '' if measured else '0'
                fixtures.write_tsv(table, rows)
                self.assertTrue(any('direction AAGCCT left_consistent' in str(failure)
                                    for failure in self.failures()))

    def test_mixed_evidence_view_checks_each_original_binary(self):
        f = self.fixture
        other = f.root/'other-evidence'
        other.mkdir()
        binary = other/'trm-evidence'
        binary.write_bytes(b'Another collector fixture; never executed')
        selected = f.evidence/'SRRknown'
        selected.rename(other/'SRRknown')
        selected.symlink_to(other/'SRRknown', target_is_directory=True)
        status_path = selected/'status.json'
        status = json.loads(status_path.read_text())
        status['identity']['binary_sha256'] = fixtures.digest(binary)
        fixtures.write_json(status_path, status)
        f.run_report()
        self.assertEqual(self.failures(), [])
        binary.write_bytes(b'Tampered collector fixture; never executed')
        failures = self.failures()
        self.assertTrue(any(failure.get('sample_id') == 'SRRknown'
                            and 'frozen evidence executable hash' in str(failure)
                            for failure in failures))

    def test_observation_only_call_passes_without_relaxing_pattern_threshold(self):
        document = self.evidence_update(full=3)
        self.assertEqual(document['classification'], 'known_repeat_observed_pattern_unresolved')
        self.observation_check(document)
        tampered = copy.deepcopy(document)
        tampered['classification'] = 'known_trm_detected'
        with self.assertRaisesRegex(AssertionError, 'known sample call'):
            self.observation_check(tampered)
        tampered = copy.deepcopy(document)
        tampered['candidates'][0]['known_detected'] = True
        with self.assertRaisesRegex(AssertionError, 'unchanged gates'):
            self.observation_check(tampered)

    def test_expected_missing_family_has_none_counts_and_is_not_explicit_zero(self):
        f = self.fixture
        rows = fixtures.read_tsv(f.assignments)
        rows[0]['related_motifs'] = 'TTAGGC,TTAGGG'
        fixtures.write_tsv(f.assignments, rows)
        f.run_report()
        document = f.document('SRRknown')
        self.observation_check(document)
        missing = next(row for row in document['expected_family_observations'] if row['motif'] == 'TTAGGG')
        self.assertEqual(missing['recovery'], 'not_reported')
        self.assertIsNone(missing['reads'])
        missing['reads'] = 0
        with self.assertRaisesRegex(AssertionError, 'not-reported expected count'):
            self.observation_check(document)
        missing['recovery'] = 'explicit_zero'
        missing['full_repeat_reads'] = 0
        with self.assertRaisesRegex(AssertionError, 'unsaved expected family recovery'):
            self.observation_check(document)

    def test_positive_original_zero_measured_cannot_be_hidden_as_explicit_zero(self):
        document = self.evidence_update(reads=0, full=0)
        self.assertEqual(document['classification'], 'detection_inconclusive_measurement')
        self.observation_check(document)
        expected = document['expected_family_observations'][0]
        self.assertEqual(expected['recovery'], 'measurement_disagreement')
        expected['recovery'] = 'explicit_zero'
        with self.assertRaisesRegex(AssertionError, 'hides a positive original'):
            self.observation_check(document)

    def test_explicit_measured_zero_remains_distinct_from_unreported_family(self):
        f = self.fixture
        source_path = f.source/'SRRknown/result.json'
        raw = json.loads(source_path.read_text())
        raw['candidates'][0].update(reads=0, repeat_bases=0, full_repeat_reads=0, max_tract=0)
        fixtures.write_json(source_path, raw)
        state_path = source_path.parent/'state.json'
        state = json.loads(state_path.read_text())
        state['result_sha256'] = fixtures.digest(source_path)
        fixtures.write_json(state_path, state)
        status_path = f.evidence/'SRRknown/status.json'
        status = json.loads(status_path.read_text())
        status['identity']['source_sha256'] = fixtures.digest(source_path)
        fixtures.write_json(status_path, status)
        document = self.evidence_update(reads=0, full=0)
        self.observation_check(document)
        expected = document['expected_family_observations'][0]
        self.assertEqual(expected['recovery'], 'explicit_zero')
        self.assertEqual(expected['reads'], 0)
        expected['recovery'] = 'not_reported'
        expected['reads'] = expected['full_repeat_reads'] = None
        with self.assertRaisesRegex(AssertionError, 'saved measured zero family recovery'):
            self.observation_check(document)

    def test_trace_classification_preserves_actual_small_count(self):
        f = self.fixture
        rows = fixtures.read_tsv(f.assignments)
        rows[0]['reference_scope'] = 'species'
        fixtures.write_tsv(f.assignments, rows)
        document = self.evidence_update(reads=3, full=0)
        self.assertEqual(document['classification'], 'published_repeat_trace_observed')
        self.observation_check(document)
        document['classification'] = 'known_repeat_observed_pattern_unresolved'
        with self.assertRaisesRegex(AssertionError, 'observation count gate'):
            self.observation_check(document)

    def test_saved_sequence_counts_must_match_displayed_representative(self):
        document = self.evidence_update()
        candidate = document['candidates'][0]
        candidate['reads'] += 1
        candidate['sequence_observation']['current_reads'] += 1
        document['top_candidate'] = candidate
        document['sequence_observation'] = candidate['sequence_observation']
        with self.assertRaisesRegex(AssertionError, 'saved selected reads'):
            self.observation_check(document)

    def test_top_candidate_cannot_be_an_unrelated_family_or_stale_count_copy(self):
        document = self.evidence_update()
        document['top_candidate']['reads'] += 1
        with self.assertRaisesRegex(AssertionError, 'representative differs'):
            self.observation_check(document)

    def test_candidate_reference_with_experimental_token_stays_candidate(self):
        f = self.fixture
        rows = fixtures.read_tsv(f.catalog)
        rows[0]['reference_status'] = 'candidate'
        fixtures.write_tsv(f.catalog, rows)
        document = self.evidence_update()
        self.assertEqual(document['classification'], 'published_candidate_observed_pattern_unresolved')
        self.observation_check(document)
        taxon_rows = audit_module.index_rows(audit_module.read_tsv(self.args.taxon_references), 'reference_id')
        audit_module.check_reference_evidence(document, rows, taxon_rows)
        document['candidates'][0]['reference_status'] = 'established'
        with self.assertRaisesRegex(AssertionError, 'reference certainty'):
            audit_module.check_reference_evidence(document, rows, taxon_rows)

    def test_candidate_taxon_reference_does_not_inherit_established_evidence_token(self):
        document = self.evidence_update()
        candidate = document['candidates'][0]
        candidate['reference_status'] = 'published_candidate'
        document['reference_assignment'].update(reference_status='candidate', related_motifs='',
            review_reference_ids='candidate_fixture', review_motifs='TTAGGC')
        row = {'reference_id': 'candidate_fixture', 'related_motifs': 'TTAGGC',
               'evidence': 'experimental', 'reference_status': 'candidate',
               'repeat_architecture': 'exact_short_unit', 'source': 'https://example.org/reference'}
        audit_module.check_reference_evidence(document, [], {'candidate_fixture': row})
        candidate['reference_status'] = 'established'
        with self.assertRaisesRegex(AssertionError, 'reference certainty'):
            audit_module.check_reference_evidence(document, [], {'candidate_fixture': row})

    def test_unrelated_novel_representative_cannot_hide_observed_exact_species_reference(self):
        document = self.evidence_update()
        novel = copy.deepcopy(self.fixture.document('SRRnovel')['top_candidate'])
        novel['expected_family_match'] = False
        document['candidates'].append(novel)
        document['reference_assignment']['reference_scope'] = 'species'
        document.update(classification='novel_geometry_nomination', top_candidate=novel,
                        detected_motif=novel['motif'], sequence_observation=novel['sequence_observation'])
        with self.assertRaisesRegex(AssertionError, 'applicable observed reference precedes'):
            audit_module.check_observation_contract(document)

    def test_observation_state_does_not_bypass_library_limit(self):
        document = self.evidence_update()
        document['library_context']['limited'] = True
        with self.assertRaisesRegex(AssertionError, 'library block'):
            self.observation_check(document)

    def test_invented_both_end_status_without_enough_junctions_is_rejected(self):
        document = self.fixture.document('SRRknown')
        document['candidates'][0]['directionality']['left_consistent'] = 1
        with self.assertRaisesRegex(AssertionError, 'both-end junction counts'):
            self.observation_check(document)

    def test_legacy_schema_does_not_require_new_observation_fields(self):
        document = self.fixture.document('SRRknown')
        document.update(schema_version=5, version='0.5.0')
        document.pop('sequence_observation')
        document.pop('expected_family_observations')
        for candidate in document['candidates']:
            for field in ('sequence_observation', 'known_observed', 'reference_observed'):
                candidate.pop(field)
        self.assertIn('legacy schema', self.observation_check(document))
        audit_module.check_interpretation(document, {'absence_inferred': 'False', 'alt_inferred': 'False'})

    def test_stale_prior_raw_and_evidence_hashes_are_detected(self):
        f = self.fixture
        for path in (f.source/'SRRknown/result.json', f.root/'original/SRRknown/result.json',
                     f.evidence/'SRRknown/evidence.json'):
            with self.subTest(path=path):
                before = path.read_bytes()
                path.write_bytes(before+b'\n')
                failures = self.failures()
                self.assertTrue(any(row['check'] == 'original_and_evidence_provenance' for row in failures), failures)
                path.write_bytes(before)

    def test_table_count_or_sample_link_tampering_is_detected(self):
        f = self.fixture
        path = f.output/'samples/SRRknown.html'
        path.write_text(path.read_text().replace('SRRknown.json', 'missing.json'))
        metrics_path = f.output/'metrics.json'
        metrics = json.loads(metrics_path.read_text())
        metrics['samples'] += 1
        fixtures.write_json(metrics_path, metrics)
        names = {row['check'] for row in self.failures()}
        self.assertIn('sample_html_links', names)
        self.assertIn('aggregate_metrics_and_source_failures', names)

    def test_missing_assumptions_cannot_be_hidden_by_adequacy_upgrade(self):
        path = self.fixture.output/'samples/SRRnone.json'
        document = json.loads(path.read_text())
        document['data_adequacy']['status'] = 'conditional_power_supported'
        document['data_adequacy']['missing_assumptions'] = []
        fixtures.write_json(path, document)
        failures = self.failures()
        self.assertTrue(any(row['check'] == 'no_absence_or_unjustified_power_upgrade' for row in failures), failures)

    def test_failed_source_cannot_be_relabelled_biological_negative(self):
        path = self.fixture.output/'samples/SRRfailed.json'
        document = json.loads(path.read_text())
        document['classification'] = 'true_absence'
        document['data_adequacy']['true_absence_inferred'] = True
        fixtures.write_json(path, document)
        names = {row['check'] for row in self.failures()}
        self.assertIn('original_and_evidence_provenance', names)
        self.assertIn('no_absence_or_unjustified_power_upgrade', names)

    def test_missing_required_artifact_emits_json_and_nonzero_exit(self):
        output = self.fixture.root/'failed-audit.json'
        run = subprocess.run([sys.executable, str(ROOT/'scripts/audit_assessment_results.py'),
            '--assessment', str(self.fixture.root/'missing'), '--output', str(output)],
            text=True, capture_output=True, timeout=10)
        self.assertEqual(run.returncode, 1, run.stderr)
        result = json.loads(run.stdout)
        self.assertEqual(result, json.loads(output.read_text()))
        self.assertFalse(result['passed'])
        self.assertEqual(result['checks'][0]['check'], 'required_artifacts')


if __name__ == '__main__':
    unittest.main()
