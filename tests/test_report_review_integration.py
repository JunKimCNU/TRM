"""Saved-native report regressions for guarded power, clade variants and projections.

All mutable input/output artifacts are synthetic and live in TemporaryDirectory.
No FASTQ, network transfer, production result, or source-code mutation is used.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from assess_adequacy import REQUIRED_ASSUMPTIONS
from evidence_batch import canonical
from report_code_dependencies import local_dependency_paths
from summarize_interpretation import summarize


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2), encoding='utf-8')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def native_candidate(motif, reads):
    family = canonical(motif)
    forward = any(family[n:] + family[:n] == motif for n in range(len(family)))
    geometry = {'full_repeat_records': 10, 'junction_reads': 20,
                'left_forward_junction_reads': 0 if forward else 10,
                'left_reverse_junction_reads': 10 if forward else 0,
                'right_forward_junction_reads': 10 if forward else 0,
                'right_reverse_junction_reads': 0 if forward else 10,
                'left_ambiguous_junction_reads': 0, 'right_ambiguous_junction_reads': 0,
                'internal_repeat_reads': 0}
    return {'motif': family, 'reads': reads, 'repeat_bases': reads * 100,
            'max_tract': 150, 'full_repeat_reads': 10, 'long_tract_reads': 0,
            'junction_evidence': geometry}


class ReportReviewIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='trm-review-integration-')
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = Path(cls.temporary.name)
        cls.source, cls.evidence, cls.output = (cls.root / name for name in ('source', 'evidence', 'assessment'))
        cls.metadata, cls.assignments = cls.root / 'metadata.tsv', cls.root / 'assignments.tsv'
        cls.catalog, cls.assumptions = cls.root / 'catalog.tsv', cls.root / 'assumptions.json'
        write_tsv(cls.catalog, [{'motif': motif, 'lineage': 'Eukaryota', 'evidence': 'experimental',
                                'source': 'https://example.org/synthetic-control'} for motif in ('TTAGGC', 'TTAGGG')])
        explicit = {'genome_size_bp': 100_000_000, 'haploid_chromosome_count': 6,
                    'random_wgs_assumed': True, 'telomere_retention_fraction': 0.1,
                    'effective_independent_read_fraction': 0.5, 'detector_sensitivity': 0.5,
                    'provenance': {key: 'Explicit synthetic hypothesis' for key in REQUIRED_ASSUMPTIONS}}
        write_json(cls.assumptions, {'defaults': explicit})
        meta, assignments = [], []
        for sample in ('SRRcanonical', 'SRRexception', 'SRRperiod', 'SRRfailed'):
            item = {'sample_id': sample, 'status': 'resolved', 'scientific_name': 'Synthetic nematode',
                    'tax_id': '99999999', 'lineage': 'Eukaryota;Metazoa;Nematoda',
                    'lineage_tax_ids': '2759;33208;6231', 'library_strategy': 'WGS',
                    'library_source': 'GENOMIC', 'library_selection': 'RANDOM',
                    'library_layout': 'PAIRED'}
            if sample == 'SRRexception':
                item.update(scientific_name='Panagrellus redivivus', tax_id='6233')
            scoped = {'sample_id': sample, 'scientific_name': item['scientific_name'], 'tax_id': item['tax_id'],
                      'lineage': item['lineage'], 'metadata_status': 'resolved',
                      'reference_scope': 'clade', 'reference_status': 'supported',
                      'repeat_architecture': 'exact_short_unit', 'related_motifs': 'TTAGGC',
                      'reference_ids': 'nematoda_default', 'review_reference_ids': '',
                      'mapping_status': 'clade_default_reference', 'prior_source': 'https://example.org/synthetic-control'}
            if sample == 'SRRexception':
                scoped.update(reference_scope='species', related_motifs='TTAGAC',
                              reference_ids='predivivus_lim2023_published', mapping_status='exact_species_reference')
            rows = ([native_candidate('TTAGGC', 1000), native_candidate('TTAGGG', 40)]
                    if sample == 'SRRcanonical' else [native_candidate('TTAGAC', 40)]
                    if sample == 'SRRexception' else [] if sample == 'SRRperiod' else [native_candidate('TTAGGG', 40)])
            summary = {'records': 10_000_000, 'bases': 1_500_000_000, 'max_record_bases': 150,
                       'sampled_prefix': False, 'call': 'candidates_only', 'elapsed_seconds': 1.0}
            qc = {'records': summary['records'], 'bases': summary['bases'], 'min_record_bases': 150,
                  'max_record_bases': 150, 'ambiguous_bases': 0, 'records_with_ambiguity': 0,
                  'quality_scores_assessed': False, 'duplicates_assessed': False,
                  'length_histogram': [{'lower': 150, 'upper': 150, 'records': summary['records']}]}
            raw = {'schema_version': 4, 'config': {'mode': 'reads', 'min_k': 3,
                   'max_k': 3 if sample == 'SRRperiod' else 32, 'seed_bases': 24, 'min_copies': 3, 'min_reads': 2},
                   'summary': summary, 'evidence_qc': qc, 'candidates': rows, 'rejected_candidates': []}
            source = cls.source / sample / 'result.json'
            write_json(source, raw)
            write_json(source.with_name('state.json'), {'sample_id': sample,
                'status': 'failed' if sample == 'SRRfailed' else 'complete', 'result_sha256': digest(source)})
            targeted = {'schema_version': 4, 'analysis_kind': 'targeted_evidence', 'summary': summary,
                        'candidate_evidence': rows, 'evidence_qc': qc, 'target_motifs': [row['motif'] for row in rows]}
            evidence = cls.evidence / sample / 'evidence.json'
            write_json(evidence, targeted)
            write_json(evidence.with_name('status.json'), {'sample_id': sample, 'status': 'complete',
                'result_sha256': digest(evidence), 'identity': {'source_sha256': digest(source), 'targets': {}}})
            meta.append(item)
            assignments.append(scoped)
        write_tsv(cls.metadata, meta)
        write_tsv(cls.assignments, assignments)
        cls.command = [sys.executable, str(ROOT / 'scripts/report_assessment.py'), '--source', str(cls.source),
                       '--metadata', str(cls.metadata), '--assignments', str(cls.assignments),
                       '--catalog', str(cls.catalog), '--taxon-references', str(ROOT / 'data/taxon_motif_references.tsv'),
                       '--assumptions', str(cls.assumptions), '--output', str(cls.output),
                       '--comparison-output', str(cls.root / 'comparison')]
        process = subprocess.run(cls.command, text=True, capture_output=True, timeout=60)
        if process.returncode:
            raise AssertionError(process.stdout + process.stderr)

    def document(self, sample):
        return json.loads((self.output / 'samples' / (sample + '.json')).read_text())

    def copied_report(self, name):
        target = self.root / name
        target.mkdir()
        for file in ('interpretation_status.json', 'metrics.json', 'samples.tsv', 'summary_projections.tsv'):
            shutil.copyfile(self.output / file, target / file)
        shutil.copytree(self.output / 'samples', target / 'samples')
        shutil.copytree(self.output / 'summaries', target / 'summaries')
        return target

    def test_canonical_completion_keeps_secondary_globally_known_clade_variant_action(self):
        report = self.document('SRRcanonical')
        self.assertEqual(report['classification'], 'known_trm_detected')
        self.assertEqual(report['detected_motif'], 'TTAGGC')
        self.assertEqual(report['next_action'], 'review_clade_variant')
        clade = report['clade_variant_observation']
        self.assertEqual(clade['status'], 'clade_canonical_and_variant_observed')
        self.assertEqual(clade['variant_observed_family_count'], 1)
        self.assertEqual(clade['canonical_observed_family_count'], 1)
        self.assertEqual(clade['review_items'][0]['global_catalogue_status'], 'globally_established_family')
        self.assertFalse(clade['biological_novelty_claimed'])
        html = (self.output / 'samples/SRRcanonical.html').read_text()
        self.assertIn('clade_variant_observations.tsv', html)
        self.assertIn('분류군 variant', html)
        rows = [row for row in read_tsv(self.output / 'clade_variant_observations.tsv')
                if row['sample_id'] == 'SRRcanonical']
        self.assertEqual(len(rows), 2)
        self.assertEqual(sum(row['review_required'] == 'True' for row in rows), 1)

    def test_species_catalogue_extension_preserves_published_exception_axis(self):
        report = self.document('SRRexception')
        self.assertEqual(report['classification'], 'known_trm_detected')
        row = report['clade_variant_observation']['families'][0]
        self.assertEqual(row['global_catalogue_status'], 'not_in_global_catalogue')
        self.assertEqual(row['sample_reference_status'], 'established')
        self.assertEqual(row['status'], 'published_clade_exception_observed')
        self.assertTrue(row['published_scoped_exception'])
        self.assertFalse(row['biological_novelty_claimed'])

    def test_explicit_three_base_power_does_not_classify_six_base_expected_absence(self):
        report = self.document('SRRperiod')
        self.assertEqual(report['data_adequacy']['status'], 'conditional_power_supported')
        self.assertEqual(report['classification'], 'detection_inconclusive_search_scope')
        # Reference-agreement guidance has independent display precedence;
        # the guarded classifier action remains explicit on its own axis.
        self.assertEqual(report['observation_next_action'], 'resolve_expected_family_search_scope')
        self.assertEqual(report['next_action'], 'measure_expected_families_before_interpreting_nondetection')
        self.assertEqual(report['expected_family_power_scope']['status'], 'expected_family_outside_power_scope')
        self.assertFalse(report['data_adequacy']['true_absence_inferred'])

    def test_failed_native_state_cannot_be_rescued_by_complete_targeted_receipt(self):
        command = self.command.copy()
        command[command.index('--output') + 1] = str(self.root / 'targeted-assessment')
        command[command.index('--comparison-output') + 1] = str(self.root / 'targeted-comparison')
        command += ['--evidence', str(self.evidence)]
        process = subprocess.run(command, text=True, capture_output=True, timeout=60)
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
        report = json.loads((self.root / 'targeted-assessment/samples/SRRfailed.json').read_text())
        self.assertEqual(report['classification'], 'analysis_failed')
        self.assertFalse(report['analysis_complete'])
        self.assertEqual(report['data_adequacy']['status'], 'observed_qc_failure')
        self.assertEqual(report['provenance']['native_result_validation']['status'], 'source_not_complete')
        self.assertFalse(report['clade_variant_observation']['review_required'])
        self.assertIsNone(report['top_candidate'])

    def test_projection_and_full_json_tables_equal_including_secondary_review(self):
        projection = self.root / 'projection-summary'
        independent = self.root / 'independent-summary'
        a = summarize(self.output, projection)
        b = summarize(self.output, independent, use_projections=False)
        self.assertEqual(a['producer_projection_samples'], 4)
        self.assertEqual(b['independently_recalculated_samples'], 4)
        for name in ('samples.tsv', 'expected_families.tsv', 'review_queue.tsv'):
            self.assertEqual((projection / name).read_bytes(), (independent / name).read_bytes(), name)
        queue = read_tsv(projection / 'review_queue.tsv')
        variant = [row for row in queue if row['sample_id'] == 'SRRcanonical'
                   and row['reason'] == 'clade_variant_requires_scoped_validation']
        exception = [row for row in queue if row['sample_id'] == 'SRRexception'
                     and row['reason'] == 'published_clade_exception_requires_sample_validation']
        self.assertEqual(len(variant), 1)
        self.assertEqual(variant[0]['is_representative'], 'False')
        self.assertEqual(len(exception), 1)

    def test_projected_and_full_json_refuse_running_state_even_with_old_complete_artifacts(self):
        copied = self.copied_report('running-report')
        state = json.loads((copied / 'interpretation_status.json').read_text())
        state['status'] = 'running'
        write_json(copied / 'interpretation_status.json', state)
        for mode in (True, False):
            with self.subTest(projections=mode):
                with self.assertRaisesRegex(ValueError, 'status is not complete'):
                    summarize(copied, self.root / ('running-summary-' + str(mode)), use_projections=mode)

    def test_projection_checks_complete_source_bytes_not_only_summary_sidecar(self):
        copied = self.copied_report('tampered-source')
        source = copied / 'samples/SRRcanonical.json'
        content = source.read_bytes()
        self.assertIn(b'"reads": 1000', content)
        source.write_bytes(content.replace(b'"reads": 1000', b'"reads": 1001', 1))
        with self.assertRaisesRegex(ValueError, 'Assessment bytes changed'):
            summarize(copied, self.root / 'tampered-source-summary')

    def test_projection_refuses_modified_projection_and_manifest(self):
        for what in ('projection', 'manifest'):
            with self.subTest(what=what):
                copied = self.copied_report('tampered-' + what)
                path = copied / ('summaries/SRRcanonical.json' if what == 'projection' else 'summary_projections.tsv')
                path.write_bytes(path.read_bytes() + b' ')
                with self.assertRaisesRegex(ValueError, 'projection.*(differ|checksum)'):
                    summarize(copied, self.root / ('tampered-' + what + '-summary'))

    def test_transitive_dependency_closure_contains_policy_and_projection_but_not_downloaders(self):
        names = {path.name for path in local_dependency_paths(ROOT / 'scripts')}
        self.assertTrue({'report_assessment.py', 'report_runtime.py', 'trm_assessment.py', 'assess_adequacy.py',
                         'clade_variant_reporting.py', 'prepare_taxon_priors.py', 'reference_evidence.py',
                         'native_result_integrity.py', 'summary_projection.py', 'summarize_interpretation.py'} <= names)
        self.assertFalse(any('download' in name.lower() or 'windows' in name.lower() for name in names))

    def test_dependency_closure_tracks_function_local_imports_and_literal_dynamic_imports(self):
        directory = self.root / 'dependency-fixture'
        directory.mkdir()
        (directory / 'report_assessment.py').write_text('def branch():\n    import policy\n')
        (directory / 'report_runtime.py').write_text('import importlib\nimportlib.import_module("adapter")\n')
        (directory / 'policy.py').write_text('from biology import evaluate\n')
        (directory / 'biology.py').write_text('def evaluate():\n    return None\n')
        (directory / 'adapter.py').write_text('import policy\n')
        (directory / 'unused_download.py').write_text('raise RuntimeError("must not import")\n')
        names = {path.name for path in local_dependency_paths(directory)}
        self.assertEqual(names, {'report_assessment.py', 'report_runtime.py', 'policy.py', 'biology.py', 'adapter.py'})
        (directory / 'adapter.py').write_text('import importlib\nimportlib.import_module(variable)\n')
        with self.assertRaisesRegex(ValueError, 'explicit literal dependency'):
            local_dependency_paths(directory)


if __name__ == '__main__':
    unittest.main()
