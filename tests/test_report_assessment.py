"""End-to-end saved-artifact interpretation without reading sequencing files."""
from __future__ import annotations

import csv
import hashlib
from html.parser import HTMLParser
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from urllib.parse import urlsplit, unquote


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/report_assessment.py'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2))


def write_tsv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def canonical(sequence):
    reverse = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    return min(value[position:]+value[:position] for value in (sequence, reverse)
               for position in range(len(sequence)))


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = []

    def handle_starttag(self, tag, attributes):
        if tag == 'a':
            self.hrefs.extend(value for key, value in attributes if key == 'href')


class AssessmentReportTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='trm-assessment-report-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root/'source'
        self.evidence = self.root/'evidence'
        self.output = self.root/'current'
        self.comparison = self.root/'comparison'
        self.source.mkdir()
        self.evidence.mkdir()
        self.catalog = self.root/'catalog.tsv'
        write_tsv(self.catalog, [{'motif': 'TTAGGC', 'lineage': 'Eukaryota; Nematoda',
                                  'evidence': 'experimental', 'source': 'https://example.org/reference'}])
        self.metadata = self.root/'metadata.tsv'
        self.assignments = self.root/'assignments.tsv'
        self.library = self.root/'library.tsv'
        self.duplicates = self.root/'duplicates.json'
        meta, assigned, library = [], [], []
        for sid in ('SRRknown', 'SRRnovel', 'SRRnone', 'SRRfailed'):
            meta.append({'sample_id': sid, 'scientific_name': 'Fixture species', 'tax_id': '6239',
                         'status': 'resolved', 'lineage': 'Eukaryota; Metazoa; Nematoda;',
                         'sample_accession': 'SAMNfixture_'+sid, 'study_accession': 'PRJfixture',
                         'library_strategy': 'WGS', 'library_source': 'GENOMIC', 'library_layout': 'PAIRED'})
            assigned.append({'sample_id': sid, 'scientific_name': 'Fixture species', 'tax_id': '6239',
                             'metadata_status': 'resolved', 'lineage': 'Eukaryota; Metazoa; Nematoda;',
                             'mapping_status': 'exact_species_reference' if sid == 'SRRknown' else 'global_prior_no_curated_reference',
                             'related_motifs': 'TTAGGC' if sid == 'SRRknown' else '',
                             'prior_source': 'https://example.org/reference' if sid == 'SRRknown' else ''})
            library.append({'sample_id': sid, 'status': 'complete', 'error': 'fixture library-note',
                            'scientific_name': 'must not override taxonomy', 'tax_id': '99999',
                            'library_selection': 'RANDOM', 'library_strategy': 'WGS',
                            'source_url': 'https://example.org/library-context'})
            if sid == 'SRRfailed':
                write_json(self.evidence/sid/'status.json', {'sample_id': sid, 'status': 'source_failed',
                                                          'reason': 'Original FASTQ was truncated'})
                continue
            family = 'AAGCCT' if sid == 'SRRknown' else 'ACAT'
            candidate = {'motif': family, 'reads': 40, 'repeat_bases': 4500, 'max_tract': 150,
                         'full_repeat_reads': 10, 'long_tract_reads': 0,
                         'telomere_prior': {'matched_reference': 'TTAGGG'},
                         'support': {'status': 'insufficient_background_resolution'}}
            if sid == 'SRRnone':
                candidate.update(reads=2, repeat_bases=48, max_tract=24, full_repeat_reads=0)
            summary = {'records': 10_000_000, 'bases': 1_500_000_000, 'max_record_bases': 150,
                       'sampled_prefix': False, 'call': 'candidates_only', 'elapsed_seconds': 1.0}
            raw = {'schema_version': 3, 'config': {'mode': 'reads', 'min_k': 3, 'max_k': 32,
                    'seed_bases': 24, 'min_copies': 3, 'min_reads': 2},
                   'summary': summary, 'candidates': [candidate], 'rejected_candidates': []}
            source_path = self.source/sid/'result.json'
            write_json(source_path, raw)
            geometry = {'full_repeat_records': candidate['full_repeat_reads'], 'junction_reads': 20,
                        'left_forward_junction_reads': 10 if sid == 'SRRknown' else 0,
                        'right_reverse_junction_reads': 10 if sid == 'SRRknown' else 0,
                        'left_reverse_junction_reads': 10 if sid == 'SRRnovel' else 0,
                        'right_forward_junction_reads': 10 if sid == 'SRRnovel' else 0,
                        'left_ambiguous_junction_reads': 0, 'right_ambiguous_junction_reads': 0,
                        'internal_repeat_reads': 0}
            if sid == 'SRRnone':
                geometry['junction_reads'] = 0
            target_candidate = {**candidate, 'junction_evidence': geometry}
            targeted = {'schema_version': 4, 'analysis_kind': 'targeted_evidence', 'summary': summary,
                        'target_motifs': [family], 'candidate_evidence': [target_candidate],
                        'evidence_qc': {'records': 10_000_000, 'bases': 1_500_000_000,
                                        'min_record_bases': 150, 'max_record_bases': 150,
                                        'ambiguous_bases': 0, 'records_with_ambiguity': 0,
                                        'quality_scores_assessed': False, 'duplicates_assessed': False,
                                        'length_histogram': [{'lower': 150, 'upper': 150, 'records': 10_000_000}]}}
            evidence_path = self.evidence/sid/'evidence.json'
            write_json(evidence_path, targeted)
            targets = {'motifs': [family], 'selection': {'novel_target_limit': 96, 'novel_candidates_not_targeted': 0}}
            write_json(self.evidence/sid/'targets.json', targets)
            write_json(self.evidence/sid/'status.json', {'sample_id': sid, 'status': 'complete',
                'result_sha256': digest(evidence_path), 'result_json': str(evidence_path.resolve()),
                'identity': {'source_sha256': digest(source_path), 'targets': targets,
                             'binary_sha256': '0'*64, 'inputs': [], 'max_records': 0, 'threads': 1}})
        write_tsv(self.metadata, meta)
        write_tsv(self.assignments, assigned)
        write_tsv(self.library, library)
        # Uses the actual archive audit schema, not guessed aliases.
        write_json(self.duplicates, {'archive_duplicate_groups': [
            {'group_id': 'archive_duplicate_fixture', 'accessions': ['SRRknown', 'SRRnovel'],
             'conclusion': 'archive_duplicate_deposit', 'identical_all_aggregate_candidate_evidence': True}]})

    def command(self):
        return [sys.executable, str(SCRIPT), '--source', str(self.source), '--evidence', str(self.evidence),
                '--metadata', str(self.metadata), '--assignments', str(self.assignments),
                '--catalog', str(self.catalog), '--library-context', str(self.library),
                '--duplicate-audit', str(self.duplicates), '--output', str(self.output),
                '--comparison-output', str(self.comparison)]

    def run_report(self, success=True, native=False, extra=()):
        command = self.command()
        command.extend(extra)
        if native:
            position = command.index('--evidence')
            del command[position:position+2]
        completed = subprocess.run(command, text=True, capture_output=True, timeout=20)
        if success:
            self.assertEqual(completed.returncode, 0, completed.stderr)
        else:
            self.assertNotEqual(completed.returncode, 0, 'Invalid artifact/report inputs unexpectedly accepted')
        return completed

    def document(self, sample):
        return json.loads((self.output/'samples'/f'{sample}.json').read_text())

    def test_source_and_metadata_must_be_explicit(self):
        for option in ('--source', '--metadata'):
            with self.subTest(option=option):
                command = self.command()
                position = command.index(option)
                del command[position:position+2]
                completed = subprocess.run(command, capture_output=True, text=True, timeout=20)
                self.assertEqual(completed.returncode, 2)
                self.assertIn('required', completed.stderr)
                self.assertIn(option, completed.stderr)
                self.assertFalse(self.output.exists())

    def test_duplicate_audit_is_not_loaded_from_historical_default_path(self):
        historical = self.root/'results/motif_review_20260920/duplicate_audit.json'
        write_json(historical, {'archive_duplicate_groups': [
            {'accessions': ['SRRknown', 'SRRnovel'], 'conclusion': 'unrelated historical fixture'}]})
        command = self.command()
        position = command.index('--duplicate-audit')
        del command[position:position+2]
        command += ['--taxon-references', str(ROOT/'data/taxon_motif_references.tsv')]
        completed = subprocess.run(command, cwd=self.root, capture_output=True, text=True, timeout=20)
        self.assertEqual(completed.returncode, 0, completed.stdout+completed.stderr)
        for sample in ('SRRknown', 'SRRnovel'):
            self.assertEqual(self.document(sample)['independent_validation']['same_archive_read_content_group'], [])

    def test_end_to_end_known_novel_unassessable_and_failed_are_distinct(self):
        self.run_report()
        known, novel, none, failed = (self.document(sid) for sid in ('SRRknown', 'SRRnovel', 'SRRnone', 'SRRfailed'))
        self.assertEqual(known['classification'], 'known_trm_detected')
        self.assertEqual(known['detected_motif'], 'TTAGGC')
        self.assertEqual(known['top_candidate']['canonical_family'], 'AAGCCT')
        self.assertEqual(known['top_candidate']['background_status'], 'insufficient_background_resolution')
        self.assertEqual(novel['classification'], 'novel_geometry_nomination')
        self.assertEqual(novel['next_action'], 'complete_candidate_evidence')
        self.assertTrue(novel['top_candidate']['geometry_nomination'])
        self.assertFalse(novel['top_candidate']['novel_candidate'])
        self.assertIsNone(novel['top_candidate']['confirmation_readiness']['distinct_junction_sequences_exact'])
        self.assertEqual(novel['top_candidate']['canonical_family'], 'ACAT')
        self.assertEqual(none['classification'], 'detection_inconclusive')
        self.assertEqual(none['repeat_observation']['policy_score_leader']['reads'], 2)
        self.assertIsNone(none['policy_selection']['representative_motif'])
        self.assertIn('반복 서열 관측', none['decision_summary']['display_label'])
        self.assertIn('미검출', none['result'])  # Compatibility field is unchanged.
        self.assertNotIn('미검출', none['decision_summary']['display_label'])
        sample_html = (self.output/'samples/SRRnone.html').read_text()
        self.assertIn('현재 양성 반복 관측', sample_html)
        self.assertIn('정책상 대표 반복', sample_html)
        none_row = next(r for r in read_tsv(self.output/'samples.tsv') if r['sample_id'] == 'SRRnone')
        self.assertEqual(none_row['policy_score_leader_reads'], '2')
        self.assertIn('반복 서열 관측', none_row['decision_display'])
        self.assertEqual(none['data_adequacy']['status'], 'unassessable_library_or_genome')
        self.assertFalse(none['data_adequacy']['true_absence_inferred'])
        self.assertEqual(failed['classification'], 'analysis_failed')
        self.assertEqual(failed['data_adequacy']['status'], 'observed_qc_failure')
        self.assertIsNone(failed['detected_motif'])
        self.assertEqual(len(read_tsv(self.output/'samples.tsv')), 4)
        self.assertGreaterEqual(len(read_tsv(self.output/'candidates.tsv')), 3)

    def test_canonical_display_anchor_and_comparison_outputs_are_separate(self):
        self.run_report()
        comparisons = {row['sample_id']: row for row in read_tsv(self.comparison/'comparison.tsv')}
        row = comparisons['SRRknown']
        self.assertEqual(row['previous_canonical_key'], 'AAGCCT')
        self.assertEqual(row['previous_reference_anchor'], 'TTAGGG')
        self.assertEqual(row['previous_motif_display'], 'TTAGGC')
        self.assertEqual(row['updated_motif_display'], 'TTAGGC')
        self.assertEqual(row['motif_family_changed'], 'False')
        self.assertTrue((self.output/'index.html').is_file())
        self.assertTrue((self.comparison/'index.html').is_file())
        self.assertFalse((self.output/'comparison.tsv').exists())
        self.assertNotEqual((self.output/'index.html').read_text(), (self.comparison/'index.html').read_text())

    def test_missing_and_zero_direction_counts_survive_json_tsv_and_html(self):
        self.run_report()
        measured_positive = self.document('SRRnovel')['top_candidate']
        evidence_path = self.evidence/'SRRknown/evidence.json'
        target = json.loads(evidence_path.read_text())
        candidate = target['candidate_evidence'][0]
        zero_geometry = {key: 0 for key in candidate['junction_evidence']}
        zero_geometry['full_repeat_records'] = candidate['full_repeat_reads']
        for measured in (False, True):
            with self.subTest(measured=measured):
                if measured:
                    candidate['junction_evidence'] = zero_geometry
                else:
                    candidate.pop('junction_evidence')
                write_json(evidence_path, target)
                status_path = self.evidence/'SRRknown/status.json'
                status = json.loads(status_path.read_text())
                status['result_sha256'] = digest(evidence_path)
                write_json(status_path, status)
                self.run_report()
                document = self.document('SRRknown')
                self.assertEqual(document['classification'], 'known_repeat_detected_boundary_unresolved')
                direction = document['top_candidate']['directionality']
                self.assertEqual(direction['status'], 'no_junction_support' if measured else 'not_measured')
                for key in ('left_consistent', 'right_consistent', 'opposite_junctions'):
                    self.assertEqual(direction[key], 0 if measured else None)
                for table in ('samples.tsv', 'candidates.tsv'):
                    row = next(row for row in read_tsv(self.output/table) if row['sample_id'] == 'SRRknown')
                    for key in ('left_consistent_junctions', 'right_consistent_junctions', 'opposite_junctions'):
                        self.assertEqual(row[key], '0' if measured else '')
                page = (self.output/'samples/SRRknown.html').read_text()
                expected = '0' if measured else '미측정'
                self.assertIn('반대 방향 경계 event: '+expected+'.', page)
                self.assertIn('3′ <b>'+expected+'</b></div>', page)
                # Another family's positive evidence, rank and nomination do not change.
                self.assertEqual(self.document('SRRnovel')['top_candidate'], measured_positive)

    def test_original_hashes_preserved_and_metadata_join_does_not_clobber_taxonomy(self):
        paths = [*self.source.rglob('*.json'), *self.evidence.rglob('*.json')]
        before = {path: digest(path) for path in paths}
        self.run_report()
        self.assertEqual(before, {path: digest(path) for path in paths})
        doc = self.document('SRRknown')
        self.assertEqual(doc['provenance']['original_discovery_sha256'], before[self.source/'SRRknown/result.json'])
        self.assertEqual(doc['provenance']['targeted_evidence_sha256'], before[self.evidence/'SRRknown/evidence.json'])
        self.assertEqual(doc['metadata']['status'], 'resolved')
        self.assertEqual(doc['metadata']['tax_id'], '6239')
        self.assertEqual(doc['metadata']['scientific_name'], 'Fixture species')
        self.assertEqual(doc['metadata']['library_context_status'], 'complete')
        self.assertEqual(doc['metadata']['library_selection'], 'RANDOM')
        self.assertEqual(doc['independent_validation']['same_archive_read_content_group'], ['SRRknown', 'SRRnovel'])
        self.assertTrue(doc['independent_validation']['search_is_live_link_not_a_verified_alternative_dataset'])

    def test_stale_original_or_targeted_content_hash_is_rejected(self):
        for kind, filename in (('source', 'result.json'), ('evidence', 'evidence.json')):
            with self.subTest(kind=kind):
                base = self.source if kind == 'source' else self.evidence
                path = base/'SRRknown'/filename
                old = path.read_bytes()
                path.write_bytes(old + b'\n')
                self.run_report(success=False)
                path.write_bytes(old)

    def test_completed_evidence_counts_and_prefix_must_match_original(self):
        path = self.evidence/'SRRknown/evidence.json'
        status_path = self.evidence/'SRRknown/status.json'
        original = json.loads(path.read_text())
        for key, value in (('records', 999), ('bases', 999), ('sampled_prefix', True)):
            with self.subTest(key=key):
                evidence = json.loads(json.dumps(original))
                evidence['summary'][key] = value
                write_json(path, evidence)
                status = json.loads(status_path.read_text())
                status['result_sha256'] = digest(path)
                write_json(status_path, status)
                self.run_report(success=False)

    def test_duplicate_metadata_or_assignment_sample_is_rejected(self):
        for path in (self.metadata, self.assignments):
            with self.subTest(path=path.name):
                original = path.read_bytes()
                rows = read_tsv(path)
                write_tsv(path, rows + [rows[0]])
                self.run_report(success=False)
                path.write_bytes(original)

    def test_html_escapes_metadata_text(self):
        rows = read_tsv(self.metadata)
        rows[0]['scientific_name'] = '<script>alert(1)</script>'
        write_tsv(self.metadata, rows)
        # Keep assignment identity consistent where strict join validation applies.
        assigned = read_tsv(self.assignments)
        assigned[0]['scientific_name'] = rows[0]['scientific_name']
        write_tsv(self.assignments, assigned)
        self.run_report()
        rendered = (self.output/'index.html').read_text()
        self.assertIn('&lt;script&gt;alert(1)&lt;/script&gt;', rendered)
        self.assertNotIn('<script>alert(1)</script>', rendered)

    def test_comparison_directory_must_not_alias_original_or_evidence(self):
        for protected in (self.source, self.evidence):
            with self.subTest(protected=protected.name):
                self.comparison = protected
                self.run_report(success=False)

    def make_native_sources(self, empty=False):
        for source in self.source.glob('*/result.json'):
            raw = json.loads(source.read_text())
            targeted = json.loads((self.evidence/source.parent.name/'evidence.json').read_text())
            raw['schema_version'] = 4
            raw['analysis_kind'] = 'de_novo_discovery'
            raw['version'] = '0.4.0'
            raw['evidence_qc'] = targeted['evidence_qc']
            raw['candidates'] = [] if empty else targeted['candidate_evidence']
            raw['rejected_candidates'] = []
            if empty:
                raw['summary']['call'] = 'no_call'
                raw['summary']['reported_candidates'] = 0
            write_json(source, raw)

    def test_native_v040_scan_reuses_geometry_without_secondary_evidence(self):
        self.make_native_sources()
        before = {path: digest(path) for path in self.source.rglob('*.json')}
        self.run_report(native=True)
        known = self.document('SRRknown')
        self.assertEqual(known['classification'], 'known_trm_detected')
        self.assertEqual(known['detected_motif'], 'TTAGGC')
        self.assertTrue(known['analysis_complete'])
        self.assertEqual(known['provenance']['targeted_evidence_result'],
                         str((self.source/'SRRknown/result.json').resolve()))
        self.assertEqual(known['provenance']['original_discovery_sha256'],
                         known['provenance']['targeted_evidence_sha256'])
        self.assertEqual(self.document('SRRnovel')['classification'], 'novel_geometry_nomination')
        self.assertEqual(before, {path: digest(path) for path in self.source.rglob('*.json')})

    def test_native_zero_candidate_scan_is_complete_but_does_not_establish_absence(self):
        self.make_native_sources(empty=True)
        self.run_report(native=True)
        for sample in ('SRRknown', 'SRRnovel', 'SRRnone'):
            with self.subTest(sample=sample):
                doc = self.document(sample)
                self.assertTrue(doc['analysis_complete'])
                self.assertEqual(doc['classification'], 'detection_inconclusive')
                self.assertIsNone(doc['detected_motif'])
                self.assertFalse(doc['data_adequacy']['true_absence_inferred'])

    @staticmethod
    def pair_fixture():
        return {'supporting_pairs': 25, 'both_mates_repeat_pairs': 15,
                'r1_only_repeat_pairs': 6, 'r2_only_repeat_pairs': 4,
                'opposite_strand_pairs': 13, 'same_strand_pairs': 1,
                'ambiguous_strand_pairs': 1, 'both_full_repeat_pairs': 4,
                'exactly_one_full_repeat_pairs': 2}

    @staticmethod
    def pair_context():
        return {'measured': True, 'layout': 'fr', 'scope': 'complete_inputs',
                'validated_pairs': 5_000_000}

    def test_native_pair_annotation_reaches_html_and_tsv_without_changing_decisions(self):
        self.make_native_sources()
        self.run_report(native=True)
        before = {sid: self.document(sid) for sid in ('SRRknown', 'SRRnovel')}
        self.assertNotIn('Paired-end 일관성', (self.output/'index.html').read_text())
        for sid in before:
            path = self.source/sid/'result.json'
            raw = json.loads(path.read_text())
            raw['paired_end'] = self.pair_context()
            raw['candidates'][0]['paired_evidence'] = self.pair_fixture()
            write_json(path, raw)
        self.run_report(native=True)
        for sid, previous in before.items():
            current = self.document(sid)
            self.assertEqual(current['classification'], previous['classification'])
            self.assertEqual(current['data_adequacy'], previous['data_adequacy'])
            self.assertEqual(current['top_candidate']['ranking_score'], previous['top_candidate']['ranking_score'])
            self.assertEqual(current['top_candidate']['novel_candidate'], previous['top_candidate']['novel_candidate'])
            self.assertEqual(current['top_candidate']['paired_evidence']['supporting_pairs'], 25)
            rendered = (self.output/f'samples/{sid}.html').read_text()
            self.assertIn('Paired-end 일관성', rendered)
            self.assertIn('역상보 13', rendered)
            self.assertIn('한쪽만 10', rendered)
            self.assertIn('novel 후보 승격', rendered)
        rows = {row['sample_id']: row for row in read_tsv(self.output/'candidates.tsv')}
        self.assertEqual(rows['SRRknown']['pair_supporting_pairs'], '25')
        self.assertEqual(rows['SRRknown']['pair_one_mate_only_pairs'], '10')
        self.assertEqual(rows['SRRnone']['pair_measurement_status'], 'not_measured')
        self.assertEqual(rows['SRRnone']['pair_supporting_pairs'], '')
        self.assertIn('Paired-end 일관성', (self.output/'index.html').read_text())

    def test_current_unpaired_target_does_not_display_old_paired_discovery(self):
        source_path = self.source/'SRRknown/result.json'
        raw = json.loads(source_path.read_text())
        raw['paired_end'] = self.pair_context()
        raw['candidates'][0]['paired_evidence'] = self.pair_fixture()
        write_json(source_path, raw)
        status_path = self.evidence/'SRRknown/status.json'
        status = json.loads(status_path.read_text())
        status['identity']['source_sha256'] = digest(source_path)
        write_json(status_path, status)
        self.run_report()
        paired = self.document('SRRknown')['top_candidate']['paired_evidence']
        self.assertEqual(paired['status'], 'not_measured')
        self.assertEqual(paired['count_source'], 'measured_evidence')
        self.assertIsNone(paired['supporting_pairs'])
        self.assertNotIn('Paired-end 일관성', (self.output/'samples/SRRknown.html').read_text())

    def test_complete_paired_target_is_displayed_from_current_measurement(self):
        evidence_path = self.evidence/'SRRknown/evidence.json'
        targeted = json.loads(evidence_path.read_text())
        targeted['paired_end'] = self.pair_context()
        targeted['candidate_evidence'][0]['paired_evidence'] = self.pair_fixture()
        write_json(evidence_path, targeted)
        status_path = self.evidence/'SRRknown/status.json'
        status = json.loads(status_path.read_text())
        status['result_sha256'] = digest(evidence_path)
        write_json(status_path, status)
        self.run_report()
        paired = self.document('SRRknown')['top_candidate']['paired_evidence']
        self.assertTrue(paired['measured'])
        self.assertEqual(paired['supporting_pairs'], 25)
        self.assertEqual(paired['count_source'], 'measured_evidence')

    def test_actual_native_paired_prefix_stays_measured_and_is_labeled_partial(self):
        from runtime_test_support import runtime_binary
        binary = runtime_binary()
        paths = []
        for mate, motif in ((1, 'TTAGGC'), (2, 'GCCTAA')):
            path = self.root/f'prefix_R{mate}.fq'
            sequence = motif * 25
            path.write_text(''.join(f'@pair{i}/{mate}\n{sequence}\n+\n'+
                                    'I'*len(sequence)+'\n' for i in range(3)))
            paths.append(path)
        run = subprocess.run([str(binary), 'scan', *map(str, paths), '--paired-layout', 'fr',
                              '--max-records', '4', '--threads', '1', '--output',
                              str(self.source/'SRRknown/result')], capture_output=True, text=True, timeout=20)
        self.assertEqual(run.returncode, 0, run.stdout+run.stderr)
        raw = json.loads((self.source/'SRRknown/result.json').read_text())
        self.assertEqual(raw['summary']['records'], 4)
        self.assertTrue(raw['summary']['sampled_prefix'])
        self.run_report(native=True)
        document = self.document('SRRknown')
        self.assertFalse(document['analysis_complete'])
        self.assertEqual(document['classification'], 'analysis_failed')
        self.assertEqual(document['data_adequacy']['status'], 'observed_qc_failure')
        candidate = next(c for c in document['candidates'] if c['canonical_family'] == 'AAGCCT')
        paired = candidate['paired_evidence']
        self.assertTrue(paired['measured'])
        self.assertEqual(paired['status'], 'measured_prefix')
        self.assertEqual(paired['scope'], 'synchronized_prefix_only')
        self.assertEqual(paired['supporting_pairs'], 2)
        self.assertIn('앞부분 pair만 검사', (self.output/'samples/SRRknown.html').read_text())
        row = next(r for r in read_tsv(self.output/'candidates.tsv') if r['sample_id'] == 'SRRknown')
        self.assertEqual(row['pair_measurement_status'], 'measured_prefix')
        self.assertEqual(row['pair_scope'], 'synchronized_prefix_only')

    def test_native_v050_provisional_keeps_external_confirmation_outstanding(self):
        self.make_native_sources()
        path = self.source/'SRRnovel/result.json'
        raw = json.loads(path.read_text())
        raw.update(version='0.5.0', schema_version=5)
        q = {'left_forward': 0, 'left_reverse': 1, 'right_forward': 0,
             'right_reverse': 0, 'ambiguous': 0}
        raw['candidates'][0]['junction_evidence']['validation'] = {
            'schema_version': 1, 'measured': True, 'clean_two_flank_repeat_records': 0,
            'junction_read_signatures': {'capacity': 128, 'observed_records': 20,
                'distinct_lower_bound': 20, 'saturated': False,
                'entries': [{'signature': f'{i:032x}', 'records': 1,
                             'normalized_quadrant_counts': q} for i in range(20)]},
            'flank_contexts': {'capacity': 128, 'flank_bases': 20, 'observed_events': 20,
                'distinct_lower_bound': 2, 'saturated': False,
                'entries': [{'signature': str(i), 'sequence': seq, 'representative_flank': seq,
                             'events': 10, 'normalized_quadrant_counts': {**q, 'left_reverse': 10}}
                            for i, seq in enumerate(['ACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCA'])]}}
        write_json(path, raw)
        self.run_report(native=True)
        doc = self.document('SRRnovel')
        top = doc['top_candidate']
        self.assertEqual(doc['classification'], 'novel_provisional_candidate')
        self.assertTrue(top['confirmation_readiness']['provisional_ready'])
        self.assertEqual(top['confirmation_readiness']['missing_evidence'], [])
        self.assertFalse(top['externally_confirmed'])
        self.assertIn('chromosome_terminal_localization_and_orientation', top['external_confirmation']['missing_evidence'])
        summary = {row['sample_id']: row for row in read_tsv(self.output/'samples.tsv')}['SRRnovel']
        self.assertIn('염색체 말단', summary['confirmation_missing'])
        self.assertEqual(summary['distinct_junction_sequences'], '20')
        self.assertEqual(doc['schema_version'], 22)
        self.assertEqual(doc['version'], '0.7.2')
        self.assertEqual(doc['reporting_version'], '1.13.0')
        self.assertEqual(summary['consistent_flank_clusters'], '2 (관측 전체)')
        self.assertEqual(summary['retained_consistent_cluster_concentration'], '보관 10/20 (50.0%)')
        annotation = top['confirmation_readiness']['directional_flank_clusters']
        self.assertFalse(annotation['raw_read_edge_counts_available'])
        self.assertIsNone(annotation['retained_clusters_consistent_at_both_original_read_edges'])
        detail = (self.output/summary['details']).read_text()
        self.assertIn('방향과 flank 군집의 결합 관측', detail)
        self.assertIn('보관 10/20 (50.0%)', detail)
        self.assertIn('전체 군집 수의 하한도 아닙니다', detail)

    def test_all_native_raw_observations_preserved_even_below_old_reporting_threshold(self):
        self.make_native_sources()
        path = self.source/'SRRnovel/result.json'
        raw = json.loads(path.read_text())
        tiny = {'motif': 'AAAT', 'reads': 1, 'repeat_bases': 28, 'max_tract': 28,
                'full_repeat_reads': 0, 'junction_evidence': {'full_repeat_records': 0, 'junction_reads': 0}}
        raw['rejected_candidates'].append(tiny)
        write_json(path, raw)
        self.run_report(native=True)
        tiny_out = next(c for c in self.document('SRRnovel')['candidates'] if c['canonical_family'] == canonical('AAAT'))
        self.assertEqual(tiny_out['reads'], 1)
        self.assertEqual(tiny_out['confirmation_readiness']['status'], 'unmeasured')
        self.assertIsNone(tiny_out['confirmation_readiness']['clean_two_flank_repeat_records'])

    def test_genuine_previous_assessment_comparison_never_substitutes_raw_scan(self):
        previous = self.root/'previous'
        old = {'program': 'trm-assessment', 'version': '0.4.0', 'sample_id': 'SRRnovel',
               'tax_id': '6239', 'classification': 'novel_trm_candidate', 'result': '이전 후보',
               'detected_motif': 'ACAT', 'top_candidate': {'canonical_family': 'ACAT'},
               'candidates': [{'novel_candidate': True}]}
        write_json(previous/'samples/SRRnovel.json', old)
        before = digest(previous/'samples/SRRnovel.json')
        self.run_report(extra=['--previous-assessment', str(previous)])
        rows = {row['sample_id']: row for row in read_tsv(self.comparison/'previous_assessment_comparison.tsv')}
        self.assertEqual(rows['SRRnovel']['previous_classification'], 'novel_trm_candidate')
        self.assertEqual(rows['SRRnovel']['updated_classification'], 'novel_geometry_nomination')
        self.assertEqual(rows['SRRnovel']['previous_geometry_nominations'], '1')
        self.assertEqual(rows['SRRknown']['comparison_status'], 'previous_assessment_missing')
        self.assertEqual(before, digest(previous/'samples/SRRnovel.json'))
        self.assertIn('previous_assessment_comparison.html', (self.comparison/'index.html').read_text())
        write_json(previous/'samples/SRRnovel.json', json.loads((self.source/'SRRnovel/result.json').read_text()))
        self.run_report(success=False, extra=['--previous-assessment', str(previous)])

    def test_exact_alternative_architecture_context_precedes_novel_promotion(self):
        rows = read_tsv(self.assignments)
        rows[1].update(reference_status='out_of_scope', reference_scope='species',
                       repeat_architecture='irregular', mapping_status='exact_species_reference')
        write_tsv(self.assignments, rows)
        self.run_report()
        doc = self.document('SRRnovel')
        self.assertEqual(doc['classification'], 'known_alternative_architecture')
        self.assertTrue(doc['top_candidate']['geometry_nomination'])

    def test_reviewed_library_evidence_identity_provenance_and_context_guard(self):
        table, evidence = self.root/'review.tsv', self.root/'review.json'
        write_json(evidence, {'sample_id': 'SRRnovel', 'findings': 'fixture only'})
        row = {'sample_id': 'SRRnovel', 'scientific_name': 'Fixture species', 'tax_id': '6239',
               'library_review_code': 'metadata_conflict_unresolved', 'detailed_evidence_path': evidence.name,
               'evidence_sha256': digest(evidence)}
        write_tsv(table, [row])
        self.run_report(extra=['--library-evidence', str(table)])
        doc = self.document('SRRnovel')
        self.assertEqual(doc['classification'], 'library_context_limited')
        self.assertTrue(doc['top_candidate']['geometry_nomination'])
        self.assertFalse(doc['library_context']['RNA_component_supported'])
        self.assertFalse(doc['library_context']['metadata_contradiction_proves_RNA'])
        self.assertEqual(doc['reviewed_library_evidence']['evidence_sha256'], digest(evidence))
        for field, value in [('tax_id', '0'), ('sample_id', 'SRRunknown'), ('evidence_sha256', '0'*64)]:
            with self.subTest(field=field):
                write_tsv(table, [{**row, field: value}])
                self.run_report(success=False, extra=['--library-evidence', str(table)])

    def test_legacy_without_secondary_or_native_evidence_is_explicit_failure(self):
        self.run_report(native=True)
        doc = self.document('SRRknown')
        self.assertFalse(doc['analysis_complete'])
        self.assertEqual(doc['classification'], 'analysis_failed')
        self.assertIsNone(doc['detected_motif'])

    def test_surrounding_taxonomy_whitespace_does_not_change_identity_or_scope(self):
        rows = read_tsv(self.metadata)
        for row in rows:
            for field in ('scientific_name', 'tax_id', 'lineage'):
                row[field] = '  '+row[field]+'  '
        write_tsv(self.metadata, rows)
        self.run_report()
        known = self.document('SRRknown')
        self.assertEqual(known['classification'], 'known_trm_detected')
        self.assertEqual(known['tax_id'].strip(), '6239')

    def test_individual_html_links_exist_and_metadata_is_escaped(self):
        injected = '<img src=x onerror="alert(1)"> & fixture'
        for path in (self.metadata, self.assignments):
            rows = read_tsv(path)
            rows[0]['scientific_name'] = injected
            write_tsv(path, rows)
        self.run_report()
        rows = read_tsv(self.output/'samples.tsv')
        for row in rows:
            page = self.output/row['details']
            self.assertEqual(page.suffix, '.html')
            self.assertTrue(page.is_file(), page)
            parser = Links()
            parser.feed(page.read_text())
            self.assertIn(row['sample_id']+'.json', parser.hrefs)
            self.assertIn('../index.html', parser.hrefs)
            for href in parser.hrefs:
                parsed = urlsplit(href)
                if not parsed.scheme and parsed.path:
                    self.assertTrue((page.parent/unquote(parsed.path)).is_file(), href)
        page = (self.output/'samples/SRRknown.html').read_text()
        self.assertIn('&lt;img src=x onerror=&quot;alert(1)&quot;&gt; &amp; fixture', page)
        self.assertNotIn('<img src=x', page)

    def test_sourced_library_override_preserves_archive_and_literature_provenance(self):
        rows = read_tsv(self.metadata)
        rows[0]['whole_genome_amplification'] = 'false'
        rows[0]['amplification_method'] = 'archive method unspecified'
        write_tsv(self.metadata, rows)
        overrides = self.root/'library_overrides.tsv'
        override = {'sample_id': 'SRRknown', 'whole_genome_amplification': 'true',
                    'amplification_method': 'MDA', 'library_evidence_source': 'https://example.org/protocol',
                    'notes': 'Fixture literature annotation; not an actual sample protocol.'}
        write_tsv(overrides, [override])
        original = self.metadata.read_bytes()
        self.run_report(extra=['--library-overrides', str(overrides)])
        doc = self.document('SRRknown')
        self.assertEqual(doc['metadata']['whole_genome_amplification'], 'true')
        self.assertEqual(doc['metadata']['amplification_method'], 'MDA')
        self.assertEqual(doc['metadata']['archive_whole_genome_amplification'], 'false')
        self.assertEqual(doc['metadata']['archive_amplification_method'], 'archive method unspecified')
        self.assertEqual(doc['metadata']['library_literature_annotation'], override)
        self.assertEqual(doc['metadata']['status'], 'resolved')
        self.assertEqual(self.metadata.read_bytes(), original)
        self.assertIn('explicit_whole_genome_amplification',
                      doc['data_adequacy']['library']['uniform_wgs_model_blockers'])

    def test_mhapla_supported_exact_reference_is_established_not_novel(self):
        motif = 'CCCAAGGTTTAAAAGG'
        family = canonical(motif)
        for path in (self.metadata, self.assignments):
            rows = read_tsv(path)
            rows[0].update(scientific_name='Meloidogyne hapla', tax_id='6305',
                           lineage='Eukaryota; Metazoa; Nematoda; Meloidogyne;')
            if path == self.assignments:
                rows[0].update(related_motifs=motif, reference_status='supported',
                               reference_ids='mhapla_fixture', reference_evidence='experimental_terminal',
                               repeat_architecture='exact_short_unit', prior_source='https://example.org/mhapla')
            write_tsv(path, rows)
        taxon_references = self.root/'taxon_references.tsv'
        write_tsv(taxon_references, [{'reference_id': 'mhapla_fixture', 'scientific_name': 'Meloidogyne hapla',
            'reference_taxid': '6305', 'scope': 'species', 'related_motifs': motif,
            'evidence': 'experimental_terminal', 'source': 'https://example.org/mhapla',
            'notes': 'Fixture reference context, no new biological claim.', 'reference_status': 'supported',
            'repeat_architecture': 'exact_short_unit', 'sample_ids': ''}])
        source_path = self.source/'SRRknown/result.json'
        raw = json.loads(source_path.read_text())
        raw['candidates'][0]['motif'] = family
        write_json(source_path, raw)
        evidence_path = self.evidence/'SRRknown/evidence.json'
        evidence = json.loads(evidence_path.read_text())
        evidence['target_motifs'] = [family]
        evidence['candidate_evidence'][0]['motif'] = family
        q = evidence['candidate_evidence'][0]['junction_evidence']
        forward = motif in {family[p:]+family[:p] for p in range(len(family))}
        q.update(left_forward_junction_reads=0 if forward else 10,
                 right_reverse_junction_reads=0 if forward else 10,
                 left_reverse_junction_reads=10 if forward else 0,
                 right_forward_junction_reads=10 if forward else 0)
        write_json(evidence_path, evidence)
        status_path = self.evidence/'SRRknown/status.json'
        status = json.loads(status_path.read_text())
        status['identity']['source_sha256'] = digest(source_path)
        status['identity']['targets']['motifs'] = [family]
        status['result_sha256'] = digest(evidence_path)
        write_json(status_path, status)
        self.run_report(extra=['--taxon-references', str(taxon_references)])
        doc = self.document('SRRknown')
        self.assertEqual(doc['classification'], 'known_trm_detected')
        self.assertEqual(doc['detected_motif'], motif)
        self.assertEqual(doc['top_candidate']['reference_status'], 'established')
        self.assertFalse(doc['top_candidate']['novel_candidate'])
        self.assertIn('https://example.org/mhapla', doc['top_candidate']['reference_sources'])

    def test_template_candidate_observation_survives_weak_full_read_support(self):
        rows = read_tsv(self.assignments)
        row = next(r for r in rows if r['sample_id'] == 'SRRnovel')
        row.update(reference_status='candidate', reference_scope='species', review_reference_ids='template_fixture',
                   review_motifs='ACAT', repeat_architecture='exact_short_unit')
        write_tsv(self.assignments, rows)
        refs = self.root/'candidate_references.tsv'
        write_tsv(refs, [{'reference_id': 'template_fixture', 'reference_status': 'candidate',
            'related_motifs': 'ACAT', 'repeat_architecture': 'exact_short_unit', 'scope': 'species',
            'scientific_name': 'Fixture species', 'reference_taxid': '6239',
            'evidence': 'telomerase_template', 'source': 'https://example.org/template'}])
        source = self.source/'SRRnovel/result.json'
        raw = json.loads(source.read_text())
        raw['candidates'][0]['full_repeat_reads'] = 3
        write_json(source, raw)
        evidence_path = self.evidence/'SRRnovel/evidence.json'
        evidence = json.loads(evidence_path.read_text())
        evidence['candidate_evidence'][0]['full_repeat_reads'] = 3
        evidence['candidate_evidence'][0]['junction_evidence']['full_repeat_records'] = 3
        write_json(evidence_path, evidence)
        status_path = self.evidence/'SRRnovel/status.json'
        status = json.loads(status_path.read_text())
        status['identity']['source_sha256'] = digest(source)
        status['result_sha256'] = digest(evidence_path)
        write_json(status_path, status)
        before = digest(source)
        self.run_report(extra=['--taxon-references', str(refs)])
        doc = self.document('SRRnovel')
        self.assertEqual(doc['classification'], 'published_candidate_observed_pattern_unresolved')
        self.assertEqual(doc['expected_motifs'], ['ACAT'])
        self.assertEqual(doc['top_candidate']['reference_status'], 'published_candidate')
        self.assertFalse(doc['top_candidate']['known_detected'])
        self.assertFalse(doc['top_candidate']['novel_candidate'])
        self.assertEqual(doc['sequence_observation']['status'], 'sequence_observed')
        self.assertEqual(doc['expected_family_observations'][0]['reads'], 40)
        self.assertEqual(doc['expected_family_observations'][0]['full_repeat_reads'], 3)
        self.assertEqual(doc['version'], '0.7.2')
        self.assertEqual(digest(source), before)
        summary = next(r for r in read_tsv(self.output/'samples.tsv') if r['sample_id'] == 'SRRnovel')
        self.assertIn('40 reads / full 3', summary['expected_repeat_observation'])
        self.assertNotIn('미검출', summary['result'])

    def test_missing_expected_family_is_not_displayed_as_explicit_zero(self):
        rows = read_tsv(self.assignments)
        row = next(r for r in rows if r['sample_id'] == 'SRRnone')
        row.update(related_motifs='TTAGGC', reference_status='supported', reference_scope='species')
        write_tsv(self.assignments, rows)
        self.run_report()
        doc = self.document('SRRnone')
        observed = doc['expected_family_observations'][0]
        self.assertEqual(observed['recovery'], 'not_reported')
        self.assertIsNone(observed['reads'])
        self.assertIsNone(observed['full_repeat_reads'])
        self.assertNotIn(canonical('TTAGGC'), {c['canonical_family'] for c in doc['candidates']})
        self.assertNotIn(canonical('TTAGGC'), {c['canonical_family'] for c in read_tsv(self.output/'candidates.tsv')
                                             if c['sample_id'] == 'SRRnone'})
        self.assertIn('명시적 0이 아님', (self.output/'samples/SRRnone.html').read_text())


class FlankApplicabilityReportingTests(unittest.TestCase):
    """Saved JSON integration only: no production FASTQ or new detector pass."""
    setUp = AssessmentReportTests.setUp
    command = AssessmentReportTests.command
    run_report = AssessmentReportTests.run_report
    document = AssessmentReportTests.document
    make_native_sources = AssessmentReportTests.make_native_sources

    def supported_source(self, sample='SRRknown', **options):
        path = self.source/sample/'result.json'
        raw = json.loads(path.read_text())
        raw.update(program='trm', version='0.8.4',
                   junction_options={'flank_bases': 20, 'end_slack': 3,
                                     'flank_period_annotation': True, **options},
                   flank_period_annotation_method={'schema_version': 1, 'scope': 'fixture all-event annotation'})
        write_json(path, raw)
        return raw

    def flank(self, sample='SRRknown'):
        return self.document(sample)['input_applicability']['flank_period_applicability']

    def detail_rows(self, sample='SRRknown'):
        return [r for r in read_tsv(self.output/'flank_period_applicability.tsv') if r['sample_id'] == sample]

    def test_period_and_expected_family_lengths_reach_all_formats_without_biology_changes(self):
        self.make_native_sources()
        assignments = read_tsv(self.assignments)
        assignments[0]['related_motifs'] += ','+'A'*26+'C'
        write_tsv(self.assignments, assignments)
        self.run_report(native=True)
        previous = self.document('SRRknown')
        raw = self.supported_source()
        original_hashes = {p: digest(p) for p in self.source.rglob('*.json')}
        self.run_report(native=True)
        current = self.document('SRRknown')
        for key in ('classification', 'result', 'detected_motif', 'top_candidate', 'candidates',
                    'data_adequacy', 'reference_agreement', 'sequence_observation',
                    'expected_family_observations', 'family_observation', 'policy_selection',
                    'decision_summary', 'next_action', 'recommendation'):
            with self.subTest(key=key):
                self.assertEqual(previous[key], current[key])
        for key in previous['input_applicability']:
            if key not in ('provenance', 'flank_period_applicability'):
                self.assertEqual(previous['input_applicability'][key], current['input_applicability'][key])
        flank = self.flank()
        self.assertTrue(flank['annotation_enabled'])
        self.assertEqual(flank['model_status'], 'supported_native_v0.8.4')
        long_period = next(r for r in flank['periods'] if r['period_bp'] == 27)
        self.assertEqual([long_period[k] for k in ('minimum_seed_bp', 'required_window_bp',
                          'junction_flank_bp', 'minimum_read_bp')], [81, 81, 20, 162])
        self.assertEqual(long_period['length_eligible_reads_upper'], 0)
        self.assertEqual(long_period['status'], 'no_flank_length_eligible_reads')
        expected = next(r for r in flank['families'] if r['period_bp'] == 27)
        self.assertEqual(expected['minimum_read_bp'], 162)
        self.assertEqual(expected['length_eligible_reads_upper'], 0)
        provenance = current['input_applicability']['provenance']
        self.assertEqual(provenance['native_model_metadata'], {k: raw[k] for k in (
            'program', 'version', 'analysis_kind', 'junction_options', 'flank_period_annotation_method')})
        rows = self.detail_rows()
        self.assertEqual(sum(r['kind'] == 'configured_period' for r in rows), 30)
        self.assertEqual(sum(r['kind'] == 'expected_family' for r in rows), 2)
        exported = next(r for r in rows if r['kind'] == 'configured_period' and r['period_bp'] == '27')
        self.assertEqual((exported['minimum_read_bp'], exported['length_eligible_reads_upper']), ('162', '0'))
        summary = next(r for r in read_tsv(self.output/'samples.tsv') if r['sample_id'] == 'SRRknown')
        for suffix, key in (('status', 'status'), ('display', 'display_label'), ('scope', 'scope')):
            self.assertEqual(summary['flank_period_applicability_'+suffix], flank[key])
        page = (self.output/'samples/SRRknown.html').read_text()
        self.assertIn('설정된 탐색 주기별 필요조건', page)
        self.assertIn('기대 family별 필요조건', page)
        self.assertIn('../flank_period_applicability.tsv', page)
        self.assertIn('실제 junction·flank 관측 계수가 아닙니다', page)
        self.assertIn('Novel 후보의 primitive 단위 길이', page)
        self.assertIn(flank['display_label'], (self.output/'index.html').read_text())
        metrics = json.loads((self.output/'metrics.json').read_text())
        self.assertEqual(metrics['schema_version'], 22)
        self.assertEqual(metrics['flank_period_applicability_status_counts'][flank['status']], 1)
        self.assertEqual(metrics['files']['flank_period_applicability.tsv'],
                         digest(self.output/'flank_period_applicability.tsv'))
        self.assertEqual(original_hashes, {p: digest(p) for p in self.source.rglob('*.json')})

    def test_legacy_and_missing_sources_export_unmeasured_not_zero(self):
        self.make_native_sources()
        self.run_report(native=True)
        flank = self.flank()
        self.assertEqual(flank['status'], 'flank_model_unavailable')
        self.assertIsNone(flank['annotation_enabled'])
        self.assertTrue(all(r['minimum_read_bp'] is None and r['length_eligible_reads_upper'] is None
                            for r in flank['periods']+flank['families']))
        self.assertTrue(all(r['minimum_read_bp'] == '' and r['length_eligible_reads_upper'] == ''
                            for r in self.detail_rows()))
        self.assertNotIn('junction_options', self.document('SRRknown')['input_applicability']['provenance']['native_model_metadata'])
        failed, = self.detail_rows('SRRfailed')
        self.assertEqual(failed['kind'], 'input')
        self.assertEqual(failed['status'], 'source_unavailable')
        self.assertEqual(failed['length_eligible_reads_upper'], '')
        self.assertIn('미평가', (self.output/'samples/SRRknown.html').read_text())

    def test_prefix_retains_bounds_without_full_input_conclusion(self):
        self.make_native_sources()
        raw = self.supported_source()
        raw['summary']['sampled_prefix'] = True
        write_json(self.source/'SRRknown/result.json', raw)
        self.run_report(native=True)
        flank = self.flank()
        self.assertEqual((flank['status'], flank['scope']), ('analysis_incomplete', 'analyzed_records'))
        period = next(r for r in flank['periods'] if r['period_bp'] == 27)
        self.assertEqual(period['length_eligible_reads_upper'], 0)
        self.assertEqual(period['analyzed_records_length_status'], 'no_flank_length_eligible_reads')
        self.assertEqual(period['status'], 'analysis_incomplete')
        self.assertFalse(self.document('SRRknown')['analysis_complete'])
        self.assertTrue(all(r['status'] == 'analysis_incomplete' and r['scope'] == 'analyzed_records'
                            for r in self.detail_rows()))
        self.assertIn('분석된 records만', (self.output/'samples/SRRknown.html').read_text())

    def test_disabled_missing_and_invalid_options_keep_their_distinct_meaning(self):
        self.make_native_sources()
        raw = self.supported_source(flank_period_annotation=False)
        for options, expected_status, enabled in (
                (raw['junction_options'], 'flank_length_eligible_reads_present_evaluation_unknown', False),
                (None, 'junction_options_unavailable', None),
                ({'flank_bases': True, 'end_slack': 3}, 'junction_options_invalid', None)):
            with self.subTest(options=options):
                changed = {**raw, 'junction_options': options}
                write_json(self.source/'SRRknown/result.json', changed)
                self.run_report(native=True)
                self.assertEqual(self.flank()['status'], expected_status)
                self.assertIs(self.flank()['annotation_enabled'], enabled)
                if enabled is None:
                    self.assertTrue(all(r['length_eligible_reads_upper'] is None for r in self.flank()['periods']))

    def test_original_model_options_are_not_replaced_by_differently_configured_targeted_scan(self):
        raw = self.supported_source(flank_bases=70)
        status_path = self.evidence/'SRRknown/status.json'
        status = json.loads(status_path.read_text())
        status['identity']['source_sha256'] = digest(self.source/'SRRknown/result.json')
        evidence_path = self.evidence/'SRRknown/evidence.json'
        targeted = json.loads(evidence_path.read_text())
        targeted.update(program='trm', version='99.0.0',
                        junction_options={'flank_bases': 20, 'end_slack': 3, 'flank_period_annotation': False})
        write_json(evidence_path, targeted)
        status['result_sha256'] = digest(evidence_path)
        write_json(status_path, status)
        self.run_report()
        family, = self.flank()['families']
        self.assertEqual((family['minimum_seed_bp'], family['junction_flank_bp'], family['minimum_read_bp']),
                         (24, 70, 94))
        self.assertTrue(self.flank()['annotation_enabled'])
        provenance = self.document('SRRknown')['input_applicability']['provenance']
        self.assertEqual(provenance['native_model_metadata']['junction_options'], raw['junction_options'])
        self.assertEqual(provenance['configuration_result_sha256'], status['identity']['source_sha256'])
        self.assertEqual(provenance['qc_source'], 'full_input_evidence')
        self.assertEqual(provenance['qc_result_sha256'], status['result_sha256'])


class SampleReferenceIndexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('assessment_report_reference_contract', SCRIPT)
        cls.module = importlib.util.module_from_spec(spec)
        sys.path.insert(0, str(ROOT/'scripts'))
        try:
            spec.loader.exec_module(cls.module)
        finally:
            sys.path.pop(0)

    def fixture(self):
        motif = 'CCCAAGGTTTAAAAGG'
        assignment = {'scientific_name': 'Meloidogyne hapla', 'tax_id': '6305',
                      'mapping_status': 'exact_species_reference', 'reference_status': 'supported',
                      'reference_ids': 'mhapla_fixture', 'related_motifs': motif,
                      'repeat_architecture': 'exact_short_unit'}
        reference = {'reference_id': 'mhapla_fixture', 'reference_status': 'supported',
                     'repeat_architecture': 'exact_short_unit', 'related_motifs': motif,
                     'evidence': 'experimental_terminal', 'source': 'https://example.org/mhapla',
                     'scope': 'species', 'reference_taxid': '6305', 'scientific_name': 'Meloidogyne hapla'}
        return assignment, reference

    def test_supported_exact_assignment_adds_reference_without_mutating_global_catalogue(self):
        assignment, row = self.fixture()
        original = {'AAGCCT': [{'motif': 'TTAGGC', 'evidence': 'experimental',
                               'source': 'https://example.org/global'}]}
        before = json.loads(json.dumps(original))
        merged = self.module.sample_reference_index(original, assignment, [row])
        family = canonical(row['related_motifs'])
        self.assertEqual(merged[family][0]['motif'], row['related_motifs'])
        self.assertEqual(merged[family][0]['reference_id'], 'mhapla_fixture')
        self.assertEqual(original, before)
        self.assertNotIn(family, original)

    def test_candidate_unknown_and_out_of_scope_statuses_never_become_established(self):
        for status in ('candidate', 'unknown', 'out_of_scope', 'mixed', ''):
            for side in ('assignment', 'reference'):
                with self.subTest(status=status, side=side):
                    assignment, row = self.fixture()
                    (assignment if side == 'assignment' else row)['reference_status'] = status
                    self.assertEqual(self.module.sample_reference_index({}, assignment, [row]), {})

    def test_already_applicable_supported_clade_assignment_is_sample_local(self):
        assignment, row = self.fixture()
        # Taxonomic applicability is validated upstream; reporting validates the join.
        assignment.update(mapping_status='curated_clade_reference', reference_scope='clade')
        row.update(scope='clade', evidence='clade_consensus', scientific_name='Fixture clade',
                   reference_taxid='999999')
        merged = self.module.sample_reference_index({}, assignment, [row])
        self.assertEqual(merged[canonical(row['related_motifs'])][0]['scope'], 'clade')

    def test_supported_mixed_short_units_join_as_documented_component_families(self):
        assignment, row = self.fixture()
        motifs = 'TTAGGC,TTAGGCTTAGAC'
        assignment.update(related_motifs=motifs, repeat_architecture='mixed_short_units')
        row.update(related_motifs=motifs, repeat_architecture='mixed_short_units')
        merged = self.module.sample_reference_index({}, assignment, [row])
        self.assertEqual(set(merged), {canonical(motif) for motif in motifs.split(',')})
        for references in merged.values():
            self.assertEqual(references[0]['repeat_architecture'], 'mixed_short_units')

    def test_reference_id_motif_evidence_and_architecture_all_must_qualify(self):
        changes = ({'reference_id': 'unrelated_reference'}, {'related_motifs': 'TTAGGC'},
                   {'evidence': 'comparative'}, {'evidence': 'computational'},
                   {'repeat_architecture': 'irregular'}, {'repeat_architecture': 'long_complex_family'},
                   {'related_motifs': 'A'*32+'C'})
        for change in changes:
            with self.subTest(change=change):
                assignment, row = self.fixture()
                row.update(change)
                self.assertEqual(self.module.sample_reference_index({}, assignment, [row]), {})

    def test_global_conventional_orientation_remains_first_when_context_adds_same_family(self):
        assignment, row = self.fixture()
        motif = row['related_motifs']
        family = canonical(motif)
        conventional = motif.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        global_rows = {family: [{'motif': conventional, 'evidence': 'experimental',
                                'source': 'https://example.org/existing-global'}]}
        merged = self.module.sample_reference_index(global_rows, assignment, [row])
        self.assertEqual(merged[family][0]['motif'], conventional)
        self.assertEqual(merged[family][1]['motif'], motif)
        self.assertEqual(len(global_rows[family]), 1)

    def test_review_id_and_motif_join_candidate_without_establishing_it(self):
        assignment, row = self.fixture()
        assignment.update(reference_status='candidate', reference_ids='', related_motifs='',
                          review_reference_ids=row['reference_id'], review_motifs=row['related_motifs'])
        row.update(reference_status='candidate', evidence='telomerase_template')
        merged = self.module.sample_reference_index({}, assignment, [row])
        family = canonical(row['related_motifs'])
        self.assertEqual(merged[family][0]['reference_status'], 'candidate')
        raw = {'motif': family, 'reads': 40, 'repeat_bases': 4000, 'max_tract': 151, 'full_repeat_reads': 10}
        candidate = self.module.assess_candidate(raw, None, merged, {family}, 151)
        self.assertEqual(candidate['reference_status'], 'published_candidate')
        self.assertFalse(candidate['known_detected'])
        for change in ({'review_reference_ids': 'wrong'}, {'review_motifs': 'TTAGGC'}, {'reference_status': 'unknown'}):
            with self.subTest(change=change):
                self.assertEqual(self.module.sample_reference_index({}, {**assignment, **change}, [row]), {})


if __name__ == '__main__':
    unittest.main()
