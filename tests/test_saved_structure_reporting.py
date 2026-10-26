"""Saved diagnostics reach exports/cache without changing biological results."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest

sys.path[:0] = [str(Path(__file__).resolve().parents[1]/'scripts'), str(Path(__file__).parent)]
import test_report_assessment as fixtures
from test_saved_repeat_structure import fixture
from test_quality_observation import histogram
from saved_repeat_structure import payload_sha256


class SavedStructureReportTests(unittest.TestCase):
    def test_overlay_export_scope_preservation_and_stale_cache_refusal(self):
        f = fixtures.AssessmentReportTests();f.setUp();self.addCleanup(f.doCleanups)
        f.make_native_sources()
        source = f.source/'SRRknown/result.json'
        raw = json.loads(source.read_text())
        raw.update(program='trm', version='0.8.5', inputs=fixture()[0]['inputs'])
        fixtures.write_json(source, raw)
        f.run_report(native=True)
        before = f.document('SRRknown')
        quality = copy.deepcopy(raw)
        quality['analysis_kind'] = 'targeted_evidence'
        quality['quality_diagnostics'] = fixture()[1]['quality_diagnostics']
        count, bases = raw['summary']['records'], raw['summary']['bases']
        quality['quality_diagnostics']['input_records'].update(observed_records=count,
            quality_records=count, whole_read=histogram(bases, 30))
        path = f.root/'quality.json';fixtures.write_json(path, quality)
        row = {'sample_id': 'SRRknown', 'identity': {'scientific_name': 'Fixture species',
            'tax_id': '6239', 'biosample': 'SAMNfixture_SRRknown'},
            'original_native_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
            'quality_result': {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()},
            'mixed_result': None}
        manifest = {'schema_version': 1, 'samples': [row]}
        manifest['payload_sha256'] = payload_sha256(manifest)
        mp = f.root/'structure.json';fixtures.write_json(mp, manifest)
        args = ('--saved-structure-evidence', str(mp), '--resume')
        f.run_report(native=True, extra=args)
        after = f.document('SRRknown')
        for key in ('candidates','classification','top_candidate','expected_family_observations',
                    'data_adequacy','reference_agreement','quality_observation'):
            self.assertEqual(before[key], after[key], key)
        self.assertTrue(after['saved_repeat_structure']['quality']['constant_across_complete_input'])
        self.assertIsNone(after['saved_repeat_structure']['mixed'])
        self.assertEqual(f.document('SRRnovel')['saved_repeat_structure']['status'], 'not_provided')
        for name in ('samples.tsv', 'saved_repeat_structure.tsv', 'index.html', 'samples/SRRknown.html'):
            self.assertIn('Q30', (f.output/name).read_text())
        receipt = json.loads((f.output/'.report-cache/SRRknown/receipt.json').read_text())
        self.assertIn(str(path), json.dumps(receipt))
        # Valid saved artifacts may be resumed; altered source bytes cannot hide
        # behind a completed cache receipt.
        f.run_report(native=True, extra=args)
        path.write_text(path.read_text()+'\n')
        result = f.run_report(success=False, native=True, extra=args)
        self.assertIn('source hash differs', result.stderr)


if __name__ == '__main__':unittest.main()
