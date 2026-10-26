"""Taxonomy/evidence contracts over tiny saved fixtures, never production reads."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(Path(__file__).parent))

import prepare_taxon_priors as priors
import report_assessment as reporting
from reference_evidence import SUPPORTED_REFERENCE_EVIDENCE
from trm_assessment import SUPPORTED_EVIDENCE, assess_candidate, resolved_eukaryotic_scope
import test_report_assessment as fixtures


class TaxonomyScopeTests(unittest.TestCase):
    def test_only_resolved_positive_eukaryotic_metadata_permits_attribution(self):
        metadata = {'status': 'resolved', 'scientific_name': 'Fixture organism',
                    'tax_id': '6239', 'lineage': 'Eukaryota; Nematoda;',
                    'library_source': 'GENOMIC'}
        self.assertTrue(resolved_eukaryotic_scope(metadata))
        for change in ({'status': ''}, {'status': 'lookup_failed'}, {'lineage': ''},
                       {'tax_id': ''}, {'scientific_name': ''}, {'tax_id': '0'},
                       {'lineage': 'Bacteria;'}, {'lineage': 'Eukaryota; Bacteria;'},
                       {'library_source': 'METAGENOMIC'}, {'metagenome': True},
                       {'metagenome': 'yes'}, {'metagenome': '1'}):
            with self.subTest(change=change):
                self.assertFalse(resolved_eukaryotic_scope({**metadata, **change}))
        self.assertTrue(resolved_eukaryotic_scope({**metadata, 'tax_id': '2759', 'lineage': ''}))


class BiologyReportContracts(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.AssessmentReportTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)

    def test_failed_or_missing_taxonomy_retains_observations_without_positive_attribution(self):
        fixture = self.fixture
        fixture.run_report()
        baseline = fixture.document('SRRknown')
        self.assertEqual(baseline['classification'], 'known_trm_detected')
        original_hash = fixtures.digest(fixture.source / 'SRRknown/result.json')
        cases = [('resolved', ''), ('lookup_failed', ''), ('', ''),
                 ('lookup_failed', 'Eukaryota; Metazoa; Nematoda;')]
        for index, (status, lineage) in enumerate(cases):
            with self.subTest(status=status, lineage=lineage):
                for path in (fixture.metadata, fixture.assignments):
                    rows = fixtures.read_tsv(path)
                    row = next(r for r in rows if r['sample_id'] == 'SRRknown')
                    row['lineage'] = lineage
                    row['status' if path == fixture.metadata else 'metadata_status'] = status
                    fixtures.write_tsv(path, rows)
                fixture.output = fixture.root / f'taxonomy-{index}'
                fixture.comparison = fixture.root / f'comparison-{index}'
                fixture.run_report()
                document = fixture.document('SRRknown')
                self.assertEqual(document['classification'], 'taxonomic_scope_unresolved')
                self.assertEqual(document['reference_agreement']['status'], 'taxonomic_scope_unresolved')
                self.assertEqual(document['candidates'], baseline['candidates'])
                self.assertEqual(document['detected_motif'], baseline['detected_motif'])
                self.assertEqual(document['reference_agreement']['main_trm_status'], 'not_confirmed')
                self.assertFalse(document['externally_confirmed'])
                self.assertEqual(fixtures.digest(fixture.source / 'SRRknown/result.json'), original_hash)

    def test_actual_panagrellus_reference_survives_prepare_and_report_roundtrip(self):
        fixture = self.fixture
        references = ROOT / 'data/taxon_motif_references.tsv'
        actual = next(row for row in priors.load_references(references)
                      if row['reference_id'] == 'predivivus_lim2023_published')
        self.assertEqual(actual['evidence'], 'published_species_reference')
        metadata = fixtures.read_tsv(fixture.metadata)
        next(row for row in metadata if row['sample_id'] == 'SRRknown').update(
            scientific_name='Panagrellus redivivus', tax_id='6233',
            lineage='Eukaryota; Metazoa; Nematoda; Panagrellus;')
        fixtures.write_tsv(fixture.metadata, metadata)
        prepared = fixture.root / 'prepared'
        command = [sys.executable, str(ROOT / 'scripts/prepare_taxon_priors.py'),
                   '--metadata', str(fixture.metadata), '--references', str(references),
                   '--output-dir', str(prepared)]
        completed = subprocess.run(command, capture_output=True, text=True, timeout=10)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        fixture.assignments = prepared / 'assignments.tsv'
        assigned = next(row for row in fixtures.read_tsv(fixture.assignments)
                        if row['sample_id'] == 'SRRknown')
        self.assertEqual(assigned['reference_ids'], actual['reference_id'])
        self.assertEqual(assigned['related_motifs'], 'TTAGAC')
        family = fixtures.canonical('TTAGAC')
        source = fixture.source / 'SRRknown/result.json'
        evidence = fixture.evidence / 'SRRknown/evidence.json'
        raw = json.loads(source.read_text())
        raw['candidates'][0]['motif'] = family
        measured = json.loads(evidence.read_text())
        measured['target_motifs'] = [family]
        measured['candidate_evidence'][0]['motif'] = family
        fixtures.write_json(source, raw)
        fixtures.write_json(evidence, measured)
        status_path = fixture.evidence / 'SRRknown/status.json'
        status = json.loads(status_path.read_text())
        status['identity']['source_sha256'] = fixtures.digest(source)
        status['result_sha256'] = fixtures.digest(evidence)
        fixtures.write_json(status_path, status)
        before = fixtures.digest(source)
        fixture.run_report(extra=['--taxon-references', str(references)])
        document = fixture.document('SRRknown')
        candidate = next(row for row in document['candidates'] if row['canonical_family'] == family)
        self.assertEqual(candidate['reference_status'], 'established')
        self.assertIn(actual['evidence'], candidate['reference_evidence'])
        self.assertIn(actual['source'], candidate['reference_sources'])
        self.assertEqual(candidate['reads'], 40)
        self.assertEqual(document['reference_assignment']['related_motifs'], 'TTAGAC')
        self.assertEqual(document['reference_agreement']['main_trm_status'], 'not_confirmed')
        self.assertFalse(candidate['novel_candidate'])
        self.assertFalse(document['externally_confirmed'])
        self.assertEqual(fixtures.digest(source), before)


class ReferenceEvidenceContracts(unittest.TestCase):
    def setUp(self):
        self.reference = next(row for row in priors.load_references(ROOT / 'data/taxon_motif_references.tsv')
                              if row['reference_id'] == 'predivivus_lim2023_published')
        self.metadata = {'sample_id': 'fixture', 'status': 'resolved',
                         'scientific_name': 'Panagrellus redivivus', 'tax_id': '6233',
                         'lineage': 'Eukaryota; Nematoda; Panagrellus;'}

    def test_every_supported_evidence_type_roundtrips_without_changing_counts(self):
        self.assertEqual(priors.DIRECT_EVIDENCE, SUPPORTED_EVIDENCE)
        for evidence in SUPPORTED_REFERENCE_EVIDENCE:
            with self.subTest(evidence=evidence):
                reference = {**self.reference, 'evidence': evidence}
                with tempfile.TemporaryDirectory(prefix='trm-reference-contract-') as temporary:
                    table = Path(temporary) / 'references.tsv'
                    serialized = {key: value for key, value in reference.items()
                                  if key not in {'motifs', 'sample_ids'}}
                    priors.write_table(table, [serialized], list(serialized))
                    reference = priors.load_references(table)[0]
                assignments, emitted = priors.assign([self.metadata], [reference])
                self.assertEqual(emitted[0]['related_motifs'], 'TTAGAC')
                index = reporting.sample_reference_index({}, assignments[0], [reference])
                family = fixtures.canonical('TTAGAC')
                self.assertEqual(index[family][0]['evidence'], evidence)
                raw = {'motif': family, 'reads': 40, 'full_repeat_reads': 10,
                       'repeat_bases': 4500, 'max_tract': 150}
                candidate = assess_candidate(raw, None, index, {family}, 150)
                self.assertEqual(candidate['reference_status'], 'established')
                self.assertEqual(candidate['reads'], 40)
                self.assertFalse(candidate['novel_candidate'])

    def test_supported_token_never_promotes_candidate_unknown_or_out_of_scope(self):
        for status in ('candidate', 'unknown', 'out_of_scope'):
            with self.subTest(status=status):
                reference = {**self.reference, 'reference_status': status}
                assignments, emitted = priors.assign([self.metadata], [reference])
                self.assertEqual(emitted, [])
                self.assertEqual(assignments[0]['reference_status'], status)
                index = reporting.sample_reference_index({}, assignments[0], [reference])
                if status == 'candidate':
                    family = fixtures.canonical('TTAGAC')
                    raw = {'motif': family, 'reads': 40, 'full_repeat_reads': 10,
                           'repeat_bases': 4500, 'max_tract': 150}
                    candidate = assess_candidate(raw, None, index, set(), 150)
                    self.assertEqual(candidate['reference_status'], 'published_candidate')
                    self.assertFalse(candidate['known_detected'])
                    self.assertFalse(candidate['novel_candidate'])
                else:
                    self.assertEqual(index, {})


if __name__ == '__main__':
    unittest.main()
