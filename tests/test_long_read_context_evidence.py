"""Integrity, coordinate and identity gates for optional long-read overlays."""
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import long_read_context_evidence as evidence


NATIVE_SHA = hashlib.sha256(b'original immutable native discovery').hexdigest()


def write_tsv(path, fields, rows):
    with Path(path).open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_fixture(directory, *, native_bound=True):
    """Small reusable fixture for CLI/cache tests; no biological result invented."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    identity = {'sample_id': 'SRR1', 'source_run': 'SRR1'}
    reads = []
    for i, (context, length, start, end, direction, maximum, tail, terminal) in enumerate((
        ('endA', 1000, 900, 974, 'increasing', 26, 26, 26),
        ('endB', 1000, 25, 100, 'decreasing', 21, 25, 21),
        ('endA', 1000, 900, 960, 'increasing', 0, 40, 0),
    ), 1):
        reads.append({**identity, 'read_key': 'r' + str(i), 'read_id': 'SRR1.' + str(i),
                      'context_id': context, 'tract_id': 'tract' + str(i),
                      'read_length': length, 'tract_start0': start, 'tract_end0': end,
                      'query_outward_direction': direction, 'outward_tail_bp': tail,
                      'native_tract_bp': end - start, 'HSP_outward_clip_bp': tail,
                      'outward_tail_maximum_exact_match_bp': maximum,
                      **{'outward_tail_matched_bp_ge' + str(t): maximum if maximum >= t else 0
                         for t in evidence.THRESHOLDS},
                      'outward_tail_full26_anywhere': maximum == 26,
                      'outward_tail_full26_at_read_end': terminal == 26,
                      'terminal20_match_bp': terminal,
                      'boundary_crossing_exact12_hits': 1 if i == 2 else 0,
                      'boundary_crossing_max_native_overlap_bp': 2 if i == 2 else 0})
    contexts = [{**identity, 'context_id': 'endA', 'target_record': 'contigA',
                 'reference_side': 'right', 'reference_length': 2000, 'linked_read_count': 2,
                 'selected_boundary_start0': 1900, 'selected_boundary_end0': 2000,
                 'selected_boundary_same_hsp_count': 2, 'selected_boundary_within100bp_count': 2,
                 'local_structure_reviewed': True, 'structure_reviewed_read_count': 1,
                 'distal_boundary_supported_read_count': 1, 'record_end_supported_read_count': 1,
                 'other_native_tract_supported_read_count': 1},
                {**identity, 'context_id': 'endB', 'target_record': 'contigB',
                 'reference_side': 'left', 'reference_length': 1000, 'linked_read_count': 1,
                 'selected_boundary_start0': 0, 'selected_boundary_end0': 100,
                 'selected_boundary_same_hsp_count': 1, 'selected_boundary_within100bp_count': 0,
                 'local_structure_reviewed': False, 'structure_reviewed_read_count': None,
                 'distal_boundary_supported_read_count': None, 'record_end_supported_read_count': None,
                 'other_native_tract_supported_read_count': None}]
    structures = [{**identity, 'read_key': 'r1', 'context_id': 'endA', 'alignment_id': 'HSP1',
                   'mapping_strand': '+', 'query_length': 1000, 'target_start0': 1026,
                   'target_end0': 2000, 'target_record_end': 2000, 'unmapped_reference_suffix_bp': 0,
                   'selected_native_tract_id': 'tract1', 'selected_native_start0': 900,
                   'selected_native_end0': 974, 'distal_boundary_reached': True,
                   'outermost_array_reached': True, 'reference_record_last_base_aligned': True,
                   'other_native_tract_count': 1, 'other_native_tract_distal_boundary_support_count': 1,
                   'selected_native_tract_distal_boundary_support': False}]
    paths = [directory / name for name in ('reads.tsv', 'contexts.tsv', 'structures.tsv')]
    for path, fields, rows in zip(paths, (evidence.READ_FIELDS, evidence.CONTEXT_FIELDS, evidence.STRUCTURE_FIELDS),
                                 (reads, contexts, structures)):
        write_tsv(path, fields, rows)
    sample = {**identity, 'scientific_name': 'Example species', 'tax_id': 123,
              'biosample': 'SAMN001', 'original_native_sha256': NATIVE_SHA if native_bound else None,
              'signature': {'signature_id': 'external_pcr_26', 'sequence': 'AAGCAGTGGTATCAACGCAGAGTACT',
                            'length_bp': 26, 'source_url': 'https://example.org/primary-protocol.pdf',
                            'applicable_runs': ['SRR1'], 'known_signature': True,
                            'sequence_derived_from_reads': False},
              'reads_source': str(paths[0]), 'contexts_source': str(paths[1]), 'structures_source': str(paths[2])}
    panel_path = directory / 'signature_panel.json'
    panel_path.write_text(json.dumps({'signatures': [sample['signature']]}))
    sample['signature_source'] = str(panel_path)
    manifest = directory / 'evidence.json'
    evidence.create_manifest(manifest, [sample], [*paths, panel_path])
    native = {'sample_id': 'SRR1', 'scientific_name': 'Example species', 'tax_id': '123',
              'analysis_complete': True, 'metadata': {'sample_accession': 'SAMN001'},
              'provenance': {'original_discovery_sha256': NATIVE_SHA},
              'candidates': [{'motif': 'TTGGG', 'reads': 10, 'rank': 1}],
              'classification': 'unchanged', 'adequacy': {'status': 'unchanged'}}
    return {'path': manifest, 'native': native, 'sample': sample,
            'reads': paths[0], 'contexts': paths[1], 'structures': paths[2], 'signature': panel_path}


class LongReadContextEvidenceTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.fixture = make_fixture(temporary.name)
        self.path = self.fixture['path']

    def load(self):
        return evidence.load_evidence(self.path)

    def modify_manifest(self, change):
        document = json.loads(self.path.read_text())
        change(document)
        document['payload_sha256'] = evidence.payload_sha256(document)
        self.path.write_text(json.dumps(document))

    def modify_table(self, name, change, reseal=True):
        path = self.fixture[name]
        with path.open(newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields, rows = reader.fieldnames, list(reader)
        change(rows)
        write_tsv(path, fields, rows)
        if reseal:
            def pins(document):
                for pin in document['source_artifacts']:
                    if pin['path'] == str(path):
                        pin['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
            self.modify_manifest(pins)

    def rejects(self, substring=None):
        with self.assertRaises(evidence.EvidenceError) as caught:
            self.load()
        if substring:
            self.assertIn(substring, str(caught.exception))

    def test_selected_counts_are_derived_and_no_native_fields_mutated(self):
        original = deepcopy(self.fixture['native'])
        result = evidence.build_evidence('SRR1', self.fixture['native'], self.load())
        self.assertEqual(self.fixture['native'], original)
        self.assertEqual(result['summary']['selected_reads'], 3)
        self.assertEqual(result['summary']['whole_tail_signature'], 1)
        self.assertEqual(result['summary']['partial_signature_ge20_at_read_end'], 1)
        self.assertEqual(result['summary']['structure_reviewed_reads'], 1)
        self.assertFalse(result['changes_classification_or_adequacy'])
        self.assertIsNone(result['physical_terminal_confirmed'])

    def test_missing_evidence_is_unassessed_with_nulls(self):
        result = evidence.build_evidence('SRR1', self.fixture['native'])
        self.assertEqual(result['status'], 'unassessed')
        self.assertTrue(all(value is None for value in result['summary'].values()))
        self.assertEqual(result['reads'], [])
        self.assertEqual(result['contexts'], [])
        self.assertEqual(evidence.load_evidence(None), {})
        self.assertEqual(evidence.source_dependencies(None), [])

    def test_unlisted_sample_remains_unassessed(self):
        result = evidence.build_evidence('sampleB', self.fixture['native'], self.load())
        self.assertEqual(result['reason'], 'sample_not_in_evidence')
        self.assertIsNone(result['summary']['selected_reads'])

    def test_unbound_source_works_only_standalone(self):
        self.modify_manifest(lambda d: d['samples'][0].update(original_native_sha256=None))
        loaded = self.load()
        self.assertEqual(evidence.build_evidence('SRR1', None, loaded)['status'], 'annotated_saved_subset')
        with self.assertRaisesRegex(evidence.EvidenceError, 'standalone only'):
            evidence.build_evidence('SRR1', self.fixture['native'], loaded)

    def test_incomplete_native_document_rejected(self):
        self.fixture['native']['analysis_complete'] = False
        with self.assertRaisesRegex(evidence.EvidenceError, 'analysis_complete'):
            evidence.build_evidence('SRR1', self.fixture['native'], self.load())

    def test_native_hash_mismatch_rejected(self):
        self.fixture['native']['provenance']['original_discovery_sha256'] = '0' * 64
        with self.assertRaisesRegex(evidence.EvidenceError, 'SHA256 mismatch'):
            evidence.build_evidence('SRR1', self.fixture['native'], self.load())

    def test_native_sample_taxon_name_and_biosample_must_all_match(self):
        for key in ('sample_id', 'scientific_name', 'tax_id', 'biosample'):
            with self.subTest(key=key):
                document = deepcopy(self.fixture['native'])
                if key == 'biosample':
                    document['metadata']['sample_accession'] = 'SAMN_OTHER'
                else:
                    document[key] = 'other'
                with self.assertRaisesRegex(evidence.EvidenceError, 'identity mismatch'):
                    evidence.build_evidence('SRR1', document, self.load())

    def test_unreviewed_context_remains_null(self):
        result = evidence.build_evidence('SRR1', None, self.load())
        context = result['contexts'][1]
        self.assertIs(context['local_structure_reviewed'], False)
        self.assertIsNone(context['structure_reviewed_read_count'])
        self.assertIsNone(context['record_end_supported_read_count'])

    def test_resealed_summary_fabrication_rejected(self):
        self.modify_manifest(lambda d: d['samples'][0]['summary'].update(selected_reads=4))
        self.rejects('disagrees')

    def test_summary_boolean_is_not_a_count(self):
        self.modify_manifest(lambda d: d['samples'][0]['summary'].update(whole_tail_signature=True))
        self.rejects('integer')

    def test_current_dependency_lookup_rejects_source_tampering(self):
        evidence.source_dependencies(self.path)  # warm cache has already loaded the sidecar
        self.fixture['reads'].write_text(self.fixture['reads'].read_text().replace('SRR1.1', 'SRR1.9'))
        with self.assertRaisesRegex(evidence.EvidenceError, 'source SHA256 mismatch'):
            evidence.source_dependencies(self.path)

    def test_manifest_tampering_rejected_without_reseal(self):
        self.path.write_text(self.path.read_text().replace('Example species', 'Other species'))
        self.rejects('payload SHA256 mismatch')

    def test_resealed_invalid_tract_coordinates_rejected(self):
        self.modify_table('reads', lambda rows: rows[0].update(tract_end0='1001'))
        self.rejects('coordinates')

    def test_resealed_tail_length_mismatch_rejected(self):
        self.modify_table('reads', lambda rows: rows[0].update(outward_tail_bp='25'))
        self.rejects('tail')

    def test_resealed_reverse_complement_direction_mismatch_rejected(self):
        self.modify_table('structures', lambda rows: rows[0].update(mapping_strand='-'))
        self.rejects('reverse-complement')

    def test_thresholds_must_be_nested(self):
        self.modify_table('reads', lambda rows: rows[0].update(outward_tail_matched_bp_ge16='27'))
        self.rejects('nested')

    def test_full_signature_cannot_exceed_partial_length(self):
        self.modify_table('reads', lambda rows: rows[1].update(outward_tail_full26_at_read_end='True'))
        self.rejects('terminal signature')

    def test_short_terminal_match_not_implicitly_assessed(self):
        self.modify_table('reads', lambda rows: rows[1].update(terminal20_match_bp='19'))
        self.rejects('threshold')

    def test_duplicate_read_rejected_even_if_resealed(self):
        self.modify_table('reads', lambda rows: rows.append(deepcopy(rows[0])))
        self.rejects('duplicate')

    def test_wrong_run_identity_rejected(self):
        self.modify_table('reads', lambda rows: rows[0].update(source_run='SRR2'))
        self.rejects('identity mismatch')

    def test_read_accession_cannot_cross_run(self):
        self.modify_table('reads', lambda rows: rows[0].update(read_id='SRR2.1'))
        self.rejects('source_run')

    def test_signature_has_to_be_applicable_to_this_run(self):
        self.modify_manifest(lambda d: d['samples'][0]['signature'].update(applicable_runs=['SRR2']))
        self.rejects('not applicable')

    def test_signature_length_and_sequence_disagree_rejected(self):
        self.modify_manifest(lambda d: d['samples'][0]['signature'].update(length_bp=27))
        self.rejects('signature length')

    def test_read_derived_signature_rejected(self):
        self.modify_manifest(lambda d: d['samples'][0]['signature'].update(sequence_derived_from_reads=True))
        self.rejects('derived')

    def test_non_external_signature_url_rejected(self):
        self.modify_manifest(lambda d: d['samples'][0]['signature'].update(source_url='file:///tmp/adapter'))
        self.rejects('HTTPS')

    def test_resealed_context_denominator_rejected(self):
        self.modify_table('contexts', lambda rows: rows[0].update(linked_read_count='3'))
        self.rejects('does not match')

    def test_unreviewed_zero_is_rejected_not_silently_interpreted(self):
        self.modify_table('contexts', lambda rows: rows[1].update(structure_reviewed_read_count='0'))
        self.rejects('must be null')

    def test_resealed_structure_read_not_in_subset_rejected(self):
        self.modify_table('structures', lambda rows: rows[0].update(read_key='r999'))
        self.rejects('unknown structure')

    def test_resealed_structure_coordinate_denominator_rejected(self):
        self.modify_table('structures', lambda rows: rows[0].update(target_record_end='1999'))
        self.rejects('coordinates')

    def test_other_tract_support_cannot_exceed_other_tract_count(self):
        self.modify_table('structures', lambda rows: rows[0].update(other_native_tract_distal_boundary_support_count='2'))
        self.rejects('inconsistent')

    def test_record_end_support_requires_alignment_to_last_base(self):
        self.modify_table('structures', lambda rows: rows[0].update(target_end0='1999', unmapped_reference_suffix_bp='1'))
        self.rejects('inconsistent')

    def test_duplicate_json_key_rejected(self):
        text = self.path.read_text()
        self.path.write_text(text.replace('"schema_version": 1,', '"schema_version": 1, "schema_version": 1,', 1))
        self.rejects('malformed JSON')

    def test_source_pin_required_for_tables(self):
        self.modify_manifest(lambda d: d['source_artifacts'].pop(0))
        self.rejects('no source pin')

    def test_duplicate_source_pin_rejected(self):
        self.modify_manifest(lambda d: d['source_artifacts'].append(deepcopy(d['source_artifacts'][0])))
        self.rejects('duplicate source')

    def test_no_structure_means_null_not_zero_summary(self):
        def change(rows):
            for row in rows:
                row['local_structure_reviewed'] = 'False'
                for field in ('structure_reviewed_read_count', 'distal_boundary_supported_read_count',
                              'record_end_supported_read_count', 'other_native_tract_supported_read_count'):
                    row[field] = ''
        self.modify_table('contexts', change)
        def clear(document):
            sample = document['samples'][0]
            sample['structures_source'] = None
            for field in ('structure_reviewed_reads', 'other_tract_support_reads',
                          'distal_boundary_support_reads', 'record_end_support_reads'):
                sample['summary'][field] = None
        self.modify_manifest(clear)
        result = evidence.build_evidence('SRR1', None, self.load())
        self.assertIsNone(result['summary']['distal_boundary_support_reads'])
        self.assertEqual(result['structures'], [])

    def test_loaded_summary_mutation_is_rejected(self):
        loaded = self.load()
        loaded['SRR1']['summary']['selected_reads'] += 1
        with self.assertRaisesRegex(evidence.EvidenceError, 'changed after validation'):
            evidence.build_evidence('SRR1', None, loaded)

    def test_source_mutation_between_load_and_build_rejected(self):
        loaded = self.load()
        self.fixture['reads'].write_text(self.fixture['reads'].read_text().replace('SRR1.1', 'SRR1.9'))
        with self.assertRaisesRegex(evidence.EvidenceError, 'source SHA256 mismatch'):
            evidence.build_evidence('SRR1', None, loaded)

    def test_alias_cannot_bind_another_native_run(self):
        for name in ('reads', 'contexts', 'structures'):
            self.modify_table(name, lambda rows: [row.update(sample_id='OTHER_RUN') for row in rows])
        self.modify_manifest(lambda d: d['samples'][0].update(sample_id='OTHER_RUN'))
        native = deepcopy(self.fixture['native'])
        native['sample_id'] = 'OTHER_RUN'
        with self.assertRaisesRegex(evidence.EvidenceError, 'sample_id must equal source_run'):
            evidence.build_evidence('OTHER_RUN', native, self.load())

    def test_resealed_signature_sequence_substitution_rejected(self):
        self.modify_manifest(lambda d: d['samples'][0]['signature'].update(sequence='A' * 26))
        self.rejects('differs from pinned source panel')

    def test_signature_source_url_substitution_rejected(self):
        self.modify_manifest(lambda d: d['samples'][0]['signature'].update(source_url='https://other.example/paper'))
        self.rejects('differs from pinned source panel')

    def test_signature_reverse_complement_mismatch_rejected(self):
        self.modify_manifest(lambda d: d['samples'][0]['signature'].update(reverse_complement='A' * 26))
        self.rejects('reverse complement')

    def test_reviewed_structure_cannot_claim_an_empty_review(self):
        self.modify_table('contexts', lambda rows: rows[0].update(structure_reviewed_read_count='0'))
        self.rejects('at least one reviewed read')

    def test_bounded_source_read_refuses_oversize_artifact(self):
        with self.fixture['reads'].open('wb') as handle:
            handle.truncate(evidence.MAX_ARTIFACT_BYTES + 1)
        self.rejects('small regular file')


if __name__ == '__main__':
    unittest.main()
