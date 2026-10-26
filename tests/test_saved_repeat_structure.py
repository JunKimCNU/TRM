import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path[:0] = [str(Path(__file__).resolve().parents[1]/'scripts'), str(Path(__file__).parent)]
import saved_repeat_structure as module
from report_mixed_arrays import GROUP_COUNTS
from test_quality_observation import histogram


def fixture():
    raw = {'program': 'trm', 'version': '0.8.5', 'analysis_kind': 'de_novo_discovery',
        'summary': {'records': 2, 'bases': 300, 'sampled_prefix': False},
        'inputs': [{'path': '/test/input.fastq.gz', 'bytes': 100, 'modified_unix_seconds': 123}]}
    quality = copy.deepcopy(raw)
    quality['analysis_kind'] = 'targeted_evidence'
    quality['quality_diagnostics'] = {
        'schema_version': 1, 'encoding': 'phred33', 'encoding_provenance': 'explicit_user_declaration_not_inferred',
        'input_scope': 'complete_inputs', 'filtering_applied': False, 'variant_discriminating_bases_assessed': False,
        'input_records': {'observed_records': 2, 'quality_records': 2, 'fasta_records_without_quality': 0,
            'fastq_records_above_retention_limit': 0, 'whole_read': histogram(300, 30)}}
    mixed = copy.deepcopy(raw)
    mixed['analysis_kind'] = 'targeted_mixed_arrays'
    group = dict.fromkeys(GROUP_COUNTS, 0)
    group.update(name='components', units=['TTAGG', 'TTAGGG'], records=2, mixed_records=1,
        nonqualifying_records=1, total_mixed_arrays=1, full_mixed_records=1, entire_mixed_records=1,
        repeat_bases=150, max_array_bp=150, max_any_array_bp=150, examples=[])
    mixed['mixed_arrays'] = {'mode': 'reads', 'settings': {'exact_only': True, 'full_fraction': .9}, 'groups': [group]}
    doc = {'sample_id': 'RUN', 'scientific_name': 'Species test', 'tax_id': '10',
        'metadata': {'sample_accession': 'BIO'}, 'analysis_complete': True,
        'provenance': {'original_discovery_sha256': 'a'*64}, 'classification': 'known_trm_detected',
        'candidates': [{'motif': 'TTAGG', 'reads': 2}], 'data_adequacy': {'status': 'unassessable'}}
    return raw, quality, mixed, doc


class SavedStructureTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.path = self.folder/'manifest.json'
        self.raw, self.quality, self.mixed, self.doc = fixture()

    def write(self, quality=True, mixed=True):
        row = {'sample_id': 'RUN', 'identity': {'scientific_name': 'Species test', 'tax_id': '10', 'biosample': 'BIO'},
            'original_native_sha256': 'a'*64, 'quality_result': None, 'mixed_result': None}
        for kind, flag in [('quality', quality), ('mixed', mixed)]:
            if flag:
                p = self.folder/(kind+'.json');p.write_text(json.dumps(getattr(self, kind)))
                row[kind+'_result'] = {'path': p.name, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
        value = {'schema_version': 1, 'samples': [row]}
        self.save_manifest(value)
        return value

    def save_manifest(self, value):
        value['payload_sha256'] = module.payload_sha256(value)
        self.path.write_text(json.dumps(value))

    def build(self):
        return module.build_evidence('RUN', self.doc, self.raw, module.load_evidence(self.path))

    def test_constant_quality_and_mixed_counts_without_biological_mutation(self):
        self.write(); before = copy.deepcopy((self.raw, self.doc))
        result = self.build()
        self.assertTrue(result['quality']['constant_across_complete_input'])
        self.assertEqual(result['quality']['quality_origin'], 'unresolved')
        self.assertEqual(result['mixed']['groups'][0]['full_mixed_records'], 1)
        self.assertFalse(result['changes_classification_or_adequacy'])
        self.assertEqual((self.raw, self.doc), before)
        self.assertEqual(len(module.source_dependencies(self.path)), 3)

    def test_any_constant_value_is_flagged_not_q30_blacklist(self):
        self.quality['quality_diagnostics']['input_records']['whole_read'] = histogram(300, 7)
        self.write();result = self.build()
        self.assertIn('Q7', result['display']['saved_quality_display'])
        self.assertTrue(result['quality']['constant_across_complete_input'])

    def test_binned_qualities_are_not_calibrated_variant_confidence(self):
        h = self.quality['quality_diagnostics']['input_records']['whole_read']
        h['phred33_counts'][30] = 100;h['phred33_counts'][40] = 200
        self.write();quality = self.build()['quality']
        self.assertEqual(quality['distribution'], 'multiple_values')
        self.assertFalse(quality['quality_calibration_verified'])

    def test_partial_constant_does_not_claim_whole_input_constancy(self):
        r = self.quality['quality_diagnostics']['input_records']
        r.update(quality_records=1, fasta_records_without_quality=1, whole_read=histogram(150, 30))
        self.write();result = self.build()
        self.assertFalse(result['quality']['constant_across_complete_input'])
        self.assertIn('quality 보유 부분', result['quality']['display'])

    def test_missing_sources_and_quality_are_unmeasured_not_zero(self):
        result = module.build_evidence('RUN', self.doc, self.raw, None)
        self.assertIsNone(result['quality']);self.assertIsNone(result['mixed'])
        self.quality.pop('quality_diagnostics');self.write(mixed=False)
        result = self.build()
        self.assertEqual(result['quality']['distribution'], 'unmeasured')
        self.assertIsNone(result['quality']['quality_bases'])
        self.assertIsNone(result['mixed'])

    def test_payload_or_source_tampering_is_rejected(self):
        value = self.write();value['samples'][0]['sample_id'] = 'OTHER'
        self.path.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, 'payload'):module.load_evidence(self.path)
        self.write();(self.folder/'quality.json').write_text('{}')
        with self.assertRaisesRegex(ValueError, 'source hash'):module.load_evidence(self.path)

    def test_duplicate_ids_and_json_keys_are_rejected(self):
        value = self.write();value['samples'].append(copy.deepcopy(value['samples'][0]));self.save_manifest(value)
        with self.assertRaisesRegex(ValueError, 'duplicate'):module.load_evidence(self.path)
        self.path.write_text('{"schema_version":1,"schema_version":1}')
        with self.assertRaisesRegex(ValueError, 'malformed'):module.load_evidence(self.path)

    def test_stale_hash_wrong_biosample_and_incomplete_analysis_rejected(self):
        self.write()
        for field in ('hash', 'biosample', 'incomplete'):
            with self.subTest(field=field):
                doc = copy.deepcopy(self.doc)
                if field == 'hash':doc['provenance']['original_discovery_sha256'] = 'b'*64
                if field == 'biosample':doc['metadata']['sample_accession'] = 'OTHER'
                if field == 'incomplete':doc['analysis_complete'] = False
                with self.assertRaises(ValueError):module.build_evidence('RUN', doc, self.raw, module.load_evidence(self.path))

    def test_partial_input_wrong_version_input_size_or_counts_rejected(self):
        original = copy.deepcopy(self.quality)
        for field in ('prefix', 'version', 'input', 'counts', 'bool_count', 'quality_scope'):
            self.quality = copy.deepcopy(original)
            with self.subTest(field=field):
                if field == 'prefix':self.quality['summary']['sampled_prefix'] = True
                if field == 'version':self.quality['version'] = '0.8.4'
                if field == 'input':self.quality['inputs'][0]['bytes'] += 1
                if field == 'counts':self.quality['summary']['records'] += 1
                if field == 'bool_count':self.quality['summary']['records'] = True
                if field == 'quality_scope':self.quality['quality_diagnostics']['input_scope'] = 'analyzed_prefix'
                self.write()
                with self.assertRaises(ValueError):self.build()

    def test_inconsistent_quality_and_mixed_counts_rejected(self):
        self.quality['quality_diagnostics']['input_records']['whole_read']['bases'] += 1
        self.write()
        with self.assertRaisesRegex(ValueError, 'histogram'):self.build()
        self.quality = fixture()[1]
        self.mixed['mixed_arrays']['groups'][0]['full_mixed_records'] = 3
        self.write()
        with self.assertRaisesRegex(ValueError, 'subset counts'):self.build()

    def test_matching_malformed_input_provenance_is_rejected(self):
        original = copy.deepcopy((self.raw, self.quality))
        for key, value in [('bytes', -1), ('bytes', True), ('modified_unix_seconds', -1),
                           ('modified_unix_seconds', True), ('modified_unix_seconds', '123'),
                           ('path', ''), ('path', None)]:
            with self.subTest(key=key, value=value):
                self.raw, self.quality = copy.deepcopy(original)
                for report in (self.raw, self.quality):report['inputs'][0][key] = value
                self.write(mixed=False)
                with self.assertRaisesRegex(ValueError, 'provenance'):self.build()
        for value in (None, [], {}):
            self.raw, self.quality = copy.deepcopy(original)
            for report in (self.raw, self.quality):report['inputs'] = value
            self.write(mixed=False)
            with self.assertRaisesRegex(ValueError, 'provenance'):self.build()

    def test_unreviewed_matching_native_version_is_rejected(self):
        self.raw['version'] = self.quality['version'] = '0.9.0'
        self.write(mixed=False)
        with self.assertRaisesRegex(ValueError, 'version'):self.build()


if __name__ == '__main__':unittest.main()
