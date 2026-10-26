"""Cache dependencies include transitive policy code, not acquisition code."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'tests')]
from report_code_dependencies import local_dependency_paths
import test_report_assessment as fixtures


class ReportDependencyScopeTests(unittest.TestCase):
    def test_transitive_local_branch_and_literal_imports_are_tracked(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'main.py').write_text('import policy\nif False:\n from helper import value\n')
            (root / 'policy.py').write_text('import main\n')
            (root / 'helper.py').write_text('import importlib\nimportlib.import_module("late")\n')
            (root / 'late.py').write_text('value=1\n')
            (root / 'download.py').write_text('raise RuntimeError("unrelated")\n')
            self.assertEqual({p.name for p in local_dependency_paths(root, ('main.py',))},
                             {'main.py', 'policy.py', 'helper.py', 'late.py'})
            (root / 'late.py').write_text('__import__(module_name)\n')
            with self.assertRaisesRegex(ValueError, 'explicit literal'):
                local_dependency_paths(root, ('main.py',))

    def test_actual_report_closure_tracks_policy_and_excludes_downloads(self):
        names = {p.name for p in local_dependency_paths(ROOT / 'scripts')}
        self.assertTrue({'trm_assessment.py', 'assess_adequacy.py', 'native_result_integrity.py',
                         'summary_projection.py', 'report_code_dependencies.py',
                         'reference_localization_evidence.py', 'cross_period_evidence.py'} <= names)
        self.assertNotIn('verify_external_download.py', names)
        self.assertNotIn('batch.py', names)

    def test_duplicate_group_change_only_recomputes_affected_samples(self):
        fixture = fixtures.AssessmentReportTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        fixture.make_native_sources()
        fixture.run_report(native=True)
        audit = json.loads(fixture.duplicates.read_text())
        audit['archive_duplicate_groups'].append({'accessions': ['SRRnone', 'outside_cohort']})
        fixtures.write_json(fixture.duplicates, audit)
        fixture.run_report(native=True, extra=('--resume',))
        rows = {r['sample_id']: r for r in fixtures.read_tsv(fixture.output / 'interpretation_resources.tsv')}
        self.assertEqual(rows['SRRnone']['mode'], 'recomputed')
        self.assertEqual(rows['SRRknown']['mode'], 'reused')
        for sample in ('SRRknown', 'SRRnone'):
            receipt = json.loads((fixture.output / '.report-cache' / sample / 'receipt.json').read_text())
            task = json.loads((Path(receipt['bundle']).parent / 'task.json').read_text())
            self.assertEqual(set(task['duplicate_groups']), {sample})
            self.assertLessEqual(set(task['args']['library_evidence_rows']), {sample})


if __name__ == '__main__':
    unittest.main()
