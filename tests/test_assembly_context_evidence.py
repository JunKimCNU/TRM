"""Small sealed fixtures exercise optional assembly annotations without sequence IO."""
import copy
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import assembly_context_evidence as module


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def json_write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def tsv_write(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def seal(value):
    value['payload_sha256'] = module.payload_sha256(value)
    return value


def write_fixture(directory, sample_id='SRRknown', identity=None, native_sha=None,
                  representative_motif='TTAGGC', review_notes=None):
    """Return paths plus a matching current report; all fixture inputs are small."""
    directory = Path(directory).absolute()
    identity = identity or dict(scientific_name='Fixture species', tax_id='6239',
                                biosample='SAMNfixture_' + sample_id)
    native_sha = native_sha or 'a' * 64
    document = dict(sample_id=sample_id, scientific_name=identity['scientific_name'],
        tax_id=identity['tax_id'], metadata=dict(sample_accession=identity['biosample']),
        provenance=dict(original_discovery_sha256=native_sha), native_scan_version='0.8.5',
        analysis_complete=True, detected_motif=representative_motif,
        classification='known_trm_detected', data_adequacy=dict(status='unassessable'),
        candidates=[dict(motif=representative_motif, reads=5)])
    group_paths = {
        'prior_assessment': directory/'assessment'/'prior.json',
        'read_summary': directory/'matched'/'read.tsv',
        'oriented_contexts': directory/'matched'/'oriented.tsv',
        'scaffold_annotations': directory/'origin'/'scaffolds.tsv',
        'fcs_review': directory/'origin'/'fcs.tsv',
        'grammar_intervals': directory/'origin'/'grammar.tsv',
        'event_extensions': directory/'origin'/'extensions.tsv',
        'quality_summary': directory/'origin'/'quality.tsv',
    }
    sources = {'prior_assessment': document,
        'read_summary': [dict(sample_id=sample_id, method='gapped', strict_mixed_intervals=4,
            identical_observed_mixed_intervals=2, oriented_same_hsp_inward_mate_join_events=1,
            oriented_same_hsp_inward_mate_join_pairs=1)],
        'quality_summary': [dict(sample_id=sample_id, quality_origin='unresolved',
            alternative_same_read_quality_compared=False, original_submitted_http403=2,
            normalized_http206=1, all_observed_bases_q30=True, local_head_complete_reads=10,
            local_head_bases=1500, normalized_decoded_fastq_records=0)],
        'fcs_review': [dict(sample_id=sample_id, report_present=True, report_data_rows=0,
            selected_scaffold_action_rows=0, positive_scaffold_taxonomy_report_available=False)],
        'scaffold_annotations': [dict(sample_id=sample_id, assembly_accession='GCA_fixture.1',
            biosample=identity['biosample'], scaffold_accession_version='SCAF1.1', length_bp=100,
            genbank_source_organism_label=identity['scientific_name'], genbank_taxon_dbxref=identity['tax_id'],
            assembly_sequence_role='unplaced-scaffold', assembly_unit='Primary Assembly',
            submitter_seqid_voucher_prefix_agrees=True, fcs_report_status='present_no_action_rows',
            fcs_scaffold_action='not_listed_in_present_assembly_report',
            host_origin_independently_confirmed=False, physical_terminal_confirmed=False,
            independent_biological_replication=False)],
        'grammar_intervals': [dict(sample_id=sample_id, assembly_accession='GCA_fixture.1',
            target='SCAF1.1', grammar_interval_id='interval1', strand='+', start0=10, end0=100,
            length_bp=90, grammar_class='mixed_grammar_neither_pure_phase_compatible',
            left_stop='incompatible_acgt', right_stop='scaffold_boundary')],
        'event_extensions': [dict(sample_id=sample_id, event_id='event1', fullread_hsp_id='hsp1',
            pair_query_id='pair1', target='SCAF1.1', g_strand_in_reference='+',
            saved_interval_start0=20, saved_interval_end0=40, saved_interval_bp=20,
            grammar_interval_id='interval1', extended_start0=10, extended_end0=100,
            extended_exact_grammar_bp=90, added_left_bp=10, added_right_bp=60,
            outward_residual_bp=0, outward_stop='scaffold_boundary', scaffold_length_bp=100,
            physical_terminal_confirmed=False, independent_evidence=False)],
        'oriented_contexts': [dict(sample_id=sample_id, event_id='event1', fullread_hsp_id='hsp1',
            pair_query_id='pair1', read_query_id='read1', target='SCAF1.1',
            mixed_interval_start0=20, mixed_interval_end0=40, scaffold_length_bp=100,
            repeat_outward_side='right', inward_mate_witness_hsp_ids='mate_hsp1',
            retained_qualifying_fullread_hsps=3, retained_qualifying_fullread_target_scaffolds=2,
            known_repeat_overlap_in_exterior_bp=5, clean_exterior_confirmed=False,
            fullread_retention_truncated=False, fullread_target_search_cap_possibly_reached=True,
            fullread_per_target_search_cap_possibly_reached=False)],
    }
    for role, path in group_paths.items():
        if role == 'prior_assessment':
            json_write(path, sources[role])
        else:
            tsv_write(path, sources[role])
    paths = dict(group_paths)
    groups = {'assessment_completion': ('assessment', ['prior_assessment']),
        'matched_completion': ('matched', ['read_summary', 'oriented_contexts']),
        'origin_completion': ('origin', ['scaffold_annotations', 'fcs_review', 'grammar_intervals',
                                        'event_extensions', 'quality_summary'])}
    for role, (folder, children) in groups.items():
        path = directory/folder/'completion.json'
        receipt = dict(status='complete', artifacts_sha256={
            str(paths[child].relative_to(path.parent)): digest(paths[child]) for child in children})
        json_write(path, receipt)
        paths[role] = path
    parsed = {role: json.loads(path.read_text()) if role in module.JSON_ROLES else
        list(csv.DictReader(io.StringIO(path.read_text()), delimiter='\t')) for role, path in paths.items()}
    payload = dict(schema_version=1, scope=module.SCOPE, sample_id=sample_id,
        **module._derive(sample_id, parsed), grammar_units=['TTAGG', 'TTAGGG'],
        same_source_assembly=True, confirmations={key: False for key in module.CONFIRMATIONS},
        review_notes=list(review_notes or []), source_artifacts=[dict(role=role,
        path=str(paths[role]), sha256=digest(paths[role])) for role in module.SOURCE_ROLES])
    seal(payload)
    payload_path = directory/'payload.json'
    json_write(payload_path, payload)
    manifest = seal(dict(schema_version=1, scope=module.SCOPE, samples=[dict(sample_id=sample_id,
        identity=identity, original_native_sha256=native_sha,
        evidence=dict(path=payload_path.name, sha256=digest(payload_path)))]))
    manifest_path = directory/'manifest.json'
    json_write(manifest_path, manifest)
    return dict(manifest_path=manifest_path, payload_path=payload_path,
                source_path=paths['read_summary'], source_paths=paths,
                document=copy.deepcopy(document), payload=payload, manifest=manifest)


class AssemblyContextEvidenceTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.directory = Path(temp.name)
        self.fx = write_fixture(self.directory)

    def load(self):
        return module.load_evidence(self.fx['manifest_path'])

    def build(self, loaded=None, document=None):
        return module.build_evidence('SRRknown', document or self.fx['document'],
                                     self.load() if loaded is None else loaded)

    def save_payload(self, payload):
        json_write(self.fx['payload_path'], seal(payload))
        manifest = json.loads(self.fx['manifest_path'].read_text())
        manifest['samples'][0]['evidence']['sha256'] = digest(self.fx['payload_path'])
        json_write(self.fx['manifest_path'], seal(manifest))

    def save_manifest(self, manifest):
        json_write(self.fx['manifest_path'], seal(manifest))

    def update_source_pin(self, role, update_completion=False):
        paths = self.fx['source_paths']
        if update_completion:
            receipt_role = ('assessment_completion' if role == 'prior_assessment' else
                            'matched_completion' if role in ('read_summary', 'oriented_contexts') else
                            'origin_completion')
            receipt_path = paths[receipt_role]
            receipt = json.loads(receipt_path.read_text())
            receipt['artifacts_sha256'][str(paths[role].relative_to(receipt_path.parent))] = digest(paths[role])
            json_write(receipt_path, receipt)
        payload = json.loads(self.fx['payload_path'].read_text())
        for pin in payload['source_artifacts']:
            pin['sha256'] = digest(paths[pin['role']])
        self.save_payload(payload)

    def test_valid_sealed_overlay_and_distinct_evidence_levels(self):
        before = copy.deepcopy(self.fx['document'])
        result = self.build()
        self.assertEqual(self.fx['document'], before)
        self.assertEqual(result['status'], 'sources_validated')
        self.assertFalse(result['changes_classification_or_adequacy'])
        self.assertEqual(result['summary']['read_mixed_intervals_selected'], 4)
        self.assertEqual(result['summary']['read_mixed_intervals_preserved'], 2)
        self.assertEqual(result['summary']['supported_assembly_intervals'], 1)
        self.assertEqual(result['summary']['outward_edge_events'], 1)
        self.assertEqual(result['summary']['max_supported_assembly_bp'], 90)
        self.assertEqual(result['summary']['physical_terminal_status'], 'not_confirmed')
        self.assertEqual(result['representative_motif'], 'TTAGGC')
        self.assertEqual(result['grammar_units'], ['TTAGG', 'TTAGGG'])
        self.assertIn('assembly 정렬에서 구간 전체 동일 보존', result['display']['assembly_context_read_display'])
        self.assertTrue(result['events'][0]['search_capped'])
        self.assertFalse(result['events'][0]['clean_exterior_confirmed'])
        self.assertTrue(all(value is False for value in result['confirmations'].values()))
        self.assertEqual(set(result['summary']), set(module.SUMMARY_FIELDS))
        self.assertEqual(set(result['display']), set(module.FIELDS))

    def test_absent_is_unmeasured_without_three_sample_claims(self):
        for sample, loaded in [('SRRknown', None), ('another_run', self.load())]:
            result = module.build_evidence(sample, self.fx['document'], loaded)
            self.assertEqual(result['status'], 'not_provided')
            self.assertTrue(all(value is None for value in result['summary'].values()))
            for key in ('scaffolds', 'intervals', 'events', 'source_artifacts'):
                self.assertEqual(result[key], [])
            self.assertEqual(result['limitations'], module.ABSENT_LIMITATIONS)
            self.assertNotIn('Q30', ' '.join(result['limitations']))
            self.assertFalse(result['changes_classification_or_adequacy'])
        self.assertIsNone(module.load_evidence(None))
        self.assertEqual(module.source_dependencies(None), [])

    def test_identity_native_version_and_representative_refuse_mismatch(self):
        mutations = [lambda d: d.update(sample_id='other'), lambda d: d.update(scientific_name='other'),
            lambda d: d.update(tax_id='999'), lambda d: d['metadata'].update(sample_accession='other'),
            lambda d: d['provenance'].update(original_discovery_sha256='b'*64),
            lambda d: d.update(native_scan_version='0.8.6'), lambda d: d.update(analysis_complete=False),
            lambda d: d.update(detected_motif='TTAGGG')]
        loaded = self.load()
        for mutation in mutations:
            document = copy.deepcopy(self.fx['document']); mutation(document)
            with self.subTest(document=document), self.assertRaises(module.EvidenceError):
                self.build(loaded, document)

    def test_all_dependencies_are_registered_and_stale_rejected(self):
        loaded = self.load()
        paths = module.source_dependencies(self.fx['manifest_path'])
        self.assertEqual(len(paths), 13)
        self.assertIn(str(self.fx['source_path']), paths)
        source = self.fx['source_path']; st = source.stat()
        original = source.read_text()
        source.write_text(original.replace('gapped', 'gapPed'))
        os.utime(source, ns=(st.st_atime_ns, st.st_mtime_ns))
        with self.assertRaisesRegex(module.EvidenceError, 'stale'):
            self.build(loaded)

    def test_missing_source_refused(self):
        self.fx['source_path'].unlink()
        with self.assertRaises(module.EvidenceError): self.load()

    def test_changed_detailed_source_digest_refused(self):
        self.fx['source_path'].write_text(self.fx['source_path'].read_text() + '\n')
        with self.assertRaisesRegex(module.EvidenceError, 'hash differs'): self.load()

    def test_resealed_source_requires_prior_completion_binding(self):
        self.fx['source_path'].write_text(self.fx['source_path'].read_text().replace('\t4\t2\t', '\t5\t2\t'))
        self.update_source_pin('read_summary')
        with self.assertRaisesRegex(module.EvidenceError, 'completed artifact seal'): self.load()

    def test_resealed_counts_cannot_diverge_from_detailed_source(self):
        payload = copy.deepcopy(self.fx['payload'])
        payload['measurements']['read_mixed_intervals_selected'] = 5
        self.save_payload(payload)
        with self.assertRaisesRegex(module.EvidenceError, 'differs from detailed'): self.load()

    def test_resealed_manifest_and_payload_identity_must_agree(self):
        manifest = copy.deepcopy(self.fx['manifest'])
        manifest['samples'][0]['identity']['biosample'] = 'different'
        self.save_manifest(manifest)
        with self.assertRaisesRegex(module.EvidenceError, 'manifest/payload identity'): self.load()

    def test_false_is_not_an_integer_count_or_schema(self):
        original = copy.deepcopy(self.fx['payload'])
        for location, key in [('measurements', 'original_fastq_http403'), ('events', 'saved_interval_start0'), ('scaffolds', 'length_bp')]:
            payload = copy.deepcopy(original)
            target = payload[location][0] if isinstance(payload[location], list) else payload[location]
            target[key] = False
            self.save_payload(payload)
            with self.subTest(key=key), self.assertRaises(module.EvidenceError): self.load()
        manifest = copy.deepcopy(self.fx['manifest']); manifest['schema_version'] = True
        self.save_manifest(manifest)
        with self.assertRaises(module.EvidenceError): self.load()

    def test_zero_is_not_confirmation_boolean(self):
        payload = copy.deepcopy(self.fx['payload']); payload['confirmations']['novel_trm'] = 0
        self.save_payload(payload)
        with self.assertRaises(module.EvidenceError): self.load()

    def test_positive_physical_or_independence_claim_is_rejected(self):
        original = copy.deepcopy(self.fx['payload'])
        for key in module.CONFIRMATIONS:
            payload = copy.deepcopy(original); payload['confirmations'][key] = True
            self.save_payload(payload)
            with self.subTest(key=key), self.assertRaises(module.EvidenceError): self.load()

    def test_duplicate_rows_and_samples_rejected(self):
        original = copy.deepcopy(self.fx['payload'])
        for field in ('scaffolds', 'intervals', 'events'):
            payload = copy.deepcopy(original); payload[field].append(copy.deepcopy(payload[field][0]))
            self.save_payload(payload)
            with self.subTest(field=field), self.assertRaises(module.EvidenceError): self.load()
        self.save_payload(original)
        manifest = copy.deepcopy(self.fx['manifest']); manifest['samples'].append(copy.deepcopy(manifest['samples'][0]))
        self.save_manifest(manifest)
        with self.assertRaisesRegex(module.EvidenceError, 'duplicate sample'): self.load()

    def test_duplicate_json_keys_and_nonfinite_numbers_rejected(self):
        path = self.fx['manifest_path']
        for text in ('{"schema_version":1,"schema_version":1}', '{"schema_version":NaN}'):
            path.write_text(text)
            with self.subTest(text=text), self.assertRaises(module.EvidenceError): self.load()

    def test_mutated_loaded_payload_refused(self):
        loaded = self.load()
        loaded['samples']['SRRknown']['document']['measurements']['read_mixed_intervals_selected'] = 8
        with self.assertRaisesRegex(module.EvidenceError, 'payload seal'): self.build(loaded)

    def test_geometry_orientation_and_denominators_checked(self):
        original = copy.deepcopy(self.fx['payload'])
        mutations = [lambda p: p['events'][0].update(saved_interval_start0=0),
            lambda p: p['events'][0].update(added_right_bp=61),
            lambda p: p['intervals'][0].update(strand='-'),
            lambda p: p['intervals'][0].update(outward_residual_bp=2),
            lambda p: p['measurements'].update(read_mixed_intervals_preserved=0),
            lambda p: p['events'][0].update(same_hsp_inward_mate_witness=False),
            lambda p: p['events'][0].update(boundary_outward_direction_consistent=False),
            lambda p: p['events'][0].update(scaffold_length_bp=200)]
        for mutation in mutations:
            payload = copy.deepcopy(original); mutation(payload); self.save_payload(payload)
            with self.subTest(payload=payload), self.assertRaises(module.EvidenceError): self.load()

    def test_source_coordinates_and_same_hsp_binding_are_checked(self):
        path = self.fx['source_paths']['oriented_contexts']
        rows = list(csv.DictReader(io.StringIO(path.read_text()), delimiter='\t'))
        rows[0]['mixed_interval_start0'] = '21'; tsv_write(path, rows)
        self.update_source_pin('oriented_contexts', update_completion=True)
        with self.assertRaisesRegex(module.EvidenceError, 'source coordinates'): self.load()

    def test_malformed_source_columns_get_controlled_error(self):
        path = self.fx['source_paths']['quality_summary']
        rows = list(csv.DictReader(io.StringIO(path.read_text()), delimiter='\t')); rows[0].pop('quality_origin')
        tsv_write(path, rows); self.update_source_pin('quality_summary', update_completion=True)
        with self.assertRaisesRegex(module.EvidenceError, 'source shape'): self.load()

    def test_q30_display_requires_observed_positive_prefix_and_no_alternate_fastq(self):
        path = self.fx['source_paths']['quality_summary']
        original = list(csv.DictReader(io.StringIO(path.read_text()), delimiter='\t'))
        for key, value in [('all_observed_bases_q30', 'False'), ('local_head_complete_reads', '0'),
                           ('local_head_bases', '0'), ('normalized_decoded_fastq_records', '1')]:
            rows = copy.deepcopy(original); rows[0][key] = value
            tsv_write(path, rows); self.update_source_pin('quality_summary', update_completion=True)
            with self.subTest(key=key), self.assertRaisesRegex(module.EvidenceError, 'Q30 prefix'): self.load()

    def test_bounds_and_extension_allowlist_prevent_large_input_reads(self):
        with patch.object(module, 'MAX_TOTAL_BYTES', 1), self.assertRaises(module.EvidenceError): self.load()
        with patch.object(module, 'MAX_BYTES', 1), self.assertRaises(module.EvidenceError): self.load()
        payload = copy.deepcopy(self.fx['payload'])
        pin = next(p for p in payload['source_artifacts'] if p['role'] == 'read_summary')
        pin['path'] = str(self.directory/'not_read.fasta'); self.save_payload(payload)
        with self.assertRaisesRegex(module.EvidenceError, 'format not supported'): self.load()

    def test_duplicate_source_roles_and_paths_refused(self):
        payload = copy.deepcopy(self.fx['payload'])
        payload['source_artifacts'][1] = copy.deepcopy(payload['source_artifacts'][0])
        self.save_payload(payload)
        with self.assertRaises(module.EvidenceError): self.load()

    def test_html_sensitive_review_note_is_preserved_for_renderer(self):
        payload = copy.deepcopy(self.fx['payload'])
        payload['review_notes'] = ['<script>must be escaped by HTML renderer</script>']
        self.save_payload(payload)
        result = self.build()
        self.assertIn('<script>', result['display']['assembly_context_limit_display'])
        self.assertIn('<script>', result['review_notes'][0])

    def test_result_rows_are_copied_without_changing_loaded_payload(self):
        loaded = self.load(); before = copy.deepcopy(loaded['samples']['SRRknown']['document'])
        result = self.build(loaded); result['events'][0]['search_capped'] = False
        result['confirmations']['physical_terminal'] = True
        self.assertEqual(loaded['samples']['SRRknown']['document'], before)


if __name__ == '__main__':
    unittest.main()
