"""Descriptive unit diagnostic and complete-source CLI contract tests."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
import repeat_unit_specificity as m
import tract_context
from repeat_region_continuity import minimum_periodic_edit_score


class SpecificityTests(unittest.TestCase):
    def test_exact_canonical_and_known_variants(self):
        for unit in ['TTAGGG', 'TTTAGGG', 'TTAGG', 'TTAGGC', 'TTAGAC', 'TTAGGT']:
            with self.subTest(unit=unit):
                x = m.evaluate_unit_specificity(unit * 30, unit)
                self.assertEqual(x['status'], 'unit_specific_evidence')
                self.assertEqual(x['cyclic_edit']['best_edit_cost'], 0)
                self.assertEqual(x['model_comparison']['fixed_normalized_gain'], 1)
                self.assertFalse(x['biological_repeat_unit_confirmed'])

    def test_exact_composition_biased_units_not_blacklisted(self):
        for unit in ['A' * 17 + 'G', 'A' * 14 + 'G', 'A' * 63 + 'G', 'AAAAAAAG', 'AAAT']:
            x = m.evaluate_unit_specificity(unit * 8, unit)
            self.assertEqual(x['status'], 'unit_specific_evidence', unit)
            for phase in x['fixed_phase']['all_best_ties']:
                for counts in phase['discriminating_bases']:
                    self.assertEqual(counts['precision'], 1)
                    self.assertEqual(counts['recall'], 1)

    def test_pure_a_does_not_support_a14g_despite_high_identity(self):
        x = m.evaluate_unit_specificity('A' * 1500, 'A' * 14 + 'G')
        self.assertEqual(x['status'], 'composition_explains_fit')
        self.assertGreater(x['fixed_phase']['all_best_ties'][0]['identity'], .93)
        self.assertIsNone(x['model_comparison']['edit_normalized_gain'])
        for phase in x['fixed_phase']['all_best_ties']:
            self.assertEqual(phase['discriminating_bases'][0]['precision'], 0)
            self.assertIsNone(phase['discriminating_bases'][0]['recall'])

    def test_homopolymer_unit_composition_only_with_null_gains(self):
        x = m.evaluate_unit_specificity('A' * 500, 'AAAA')
        self.assertEqual(x['primitive_motif'], 'A')
        self.assertTrue(x['supplied_motif_reduced'])
        self.assertEqual(x['status'], 'composition_explains_fit')
        self.assertEqual(x['cyclic_edit']['best_edit_cost'], 0)
        self.assertIsNone(x['model_comparison']['fixed_normalized_gain'])

    def test_shuffled_biased_composition_not_unit_specific(self):
        seq = list(('A' * 17 + 'G') * 40)
        random.Random(9401).shuffle(seq)
        x = m.evaluate_unit_specificity(''.join(seq), 'A' * 17 + 'G')
        self.assertNotEqual(x['status'], 'unit_specific_evidence')
        self.assertEqual(x['majority_baseline']['base_counts']['G'], 40)

    def test_indels_preserve_model_disagreement(self):
        source = 'TTAGGG' * 60
        seq = source[:180] + 'C' + source[180:]
        x = m.evaluate_unit_specificity(seq, 'TTAGGG')
        self.assertEqual(x['cyclic_edit']['best_edit_cost'], 1)
        self.assertEqual(x['status'], 'ambiguous')
        self.assertTrue(x['model_comparison']['operational_model_disagreement'])
        self.assertIn('fixed_phase_and_indel_tolerant_operational_support_disagree', x['reason_codes'])

    def test_no_free_phase_reset_or_query_clipping(self):
        unit = 'ACGT'
        seq = unit * 15 + 'A' * 25 + unit * 15
        x = m.evaluate_unit_specificity(seq, unit)
        self.assertEqual(x['cyclic_edit']['best_edit_cost'], minimum_periodic_edit_score(seq, unit)['edit_cost'])
        self.assertGreater(x['cyclic_edit']['best_edit_cost'], 0)
        self.assertEqual(x['query_bp'], len(seq))

    def test_reverse_complement_and_rotations(self):
        unit = 'A' * 17 + 'G'
        seq = unit * 10 + 'AAG'
        baseline = m.evaluate_unit_specificity(seq, unit)
        for shift in range(len(unit)):
            rotated = unit[shift:] + unit[:shift]
            for s in [seq, tract_context.reverse_complement(seq)]:
                x = m.evaluate_unit_specificity(s, rotated, start0=115)
                self.assertEqual(x['status'], baseline['status'])
                self.assertEqual(x['cyclic_edit']['best_edit_cost'], baseline['cyclic_edit']['best_edit_cost'])
                self.assertEqual(x['fixed_phase']['all_best_ties'][0]['matches_bp'], baseline['fixed_phase']['all_best_ties'][0]['matches_bp'])
                self.assertEqual(x['start0'], 115)

    def test_all_phase_and_majority_ties_preserved(self):
        x = m.evaluate_unit_specificity('ACGT' * 10, 'ACGT')
        self.assertEqual(x['fixed_phase']['all_best_tie_count'], 2)
        self.assertEqual(x['fixed_phase']['all_best_distinct_phase_templates'], 1)
        self.assertEqual(x['majority_baseline']['all_tied_majority_bases'], list('ACGT'))
        for row in x['fixed_phase']['all_best_ties']:
            self.assertEqual(len(row['discriminating_bases']), 4)
        x = m.evaluate_unit_specificity('A' * 30, 'AAG')
        self.assertEqual(x['fixed_phase']['all_best_tie_count'], 3)

    def test_finite_template_composition_null_is_phase_specific(self):
        x = m.evaluate_unit_specificity('AAGGA', 'AAGG')
        best = x['fixed_phase']['all_best_ties'][0]
        self.assertAlmostEqual(best['composition_null_expected_matches'], 13 / 5)
        self.assertAlmostEqual(best['composition_null_expected_identity'], 13 / 25)
        self.assertIn('Not corrected', x['fixed_phase']['composition_null_scope'])

    def test_three_copy_and_discriminating_event_scope(self):
        x = m.evaluate_unit_specificity('TTAGGG' * 2, 'TTAGGG')
        self.assertEqual(x['status'], 'ambiguous')
        self.assertIn('insufficient_complete_unit_copies_for_operational_support', x['reason_codes'])
        x = m.evaluate_unit_specificity(('A' * 17 + 'G') * 3, 'A' * 17 + 'G')
        self.assertEqual(x['status'], 'unit_specific_evidence')

    def test_invalid_symbols_and_empty_not_negative_evidence(self):
        for sequence, unit in [('', 'A'), ('AN', 'A'), ('RYS', 'AC'), ('ACGT', ''), ('ACGT', 'N'), ('ACGT', 'A' * 257)]:
            x = m.evaluate_unit_specificity(sequence, unit)
            self.assertEqual(x['status'], 'not_evaluated')
            self.assertFalse(x['evaluated'])
            self.assertIsNone(x['fixed_phase'])

    def test_lowercase_explicit_normalization_and_power_reduction(self):
        x = m.evaluate_unit_specificity('ttaggg' * 12, 'ttagggttaggg')
        self.assertEqual(x['primitive_motif'], 'TTAGGG')
        self.assertEqual(x['supplied_motif'], 'ttagggttaggg')
        self.assertEqual(x['status'], 'unit_specific_evidence')

    def test_resource_limits_do_not_score_prefix(self):
        x = m.evaluate_unit_specificity('TTAGGG' * 2000, 'TTAGGG')
        self.assertEqual(x['status'], 'not_evaluated')
        self.assertEqual(x['query_bp'], 12000)
        self.assertIn('interval_exceeds_registered_length_bound', x['reason_codes'])
        x = m.evaluate_unit_specificity('TTAGGG' * 10, 'TTAGGG', config={'max_dp_cells': 719})
        self.assertEqual(x['status'], 'not_evaluated')
        self.assertIn('interval_exceeds_registered_dp_work_bound', x['reason_codes'])

    def test_config_cannot_fit_thresholds_to_controls(self):
        for c in [{'minimum_fixed_identity': .85}, {'minimum_complete_unit_copies': 2}, {'max_region_bp': True}, {'max_region_bp': 100001}, {'unknown': 1}]:
            with self.assertRaises(ValueError):
                m.evaluate_unit_specificity('ACGT' * 10, 'ACGT', config=c)
        self.assertEqual(m.validate_config({'max_region_bp': 100000})['max_region_bp'], 100000)

    def test_residue_generalization_matches_existing_public_oracle(self):
        rng = random.Random(551203)
        for k in range(1, 65):
            unit = ''.join(rng.choices('ACGT', k=k))
            unit = m.primitive_unit(unit)
            for length in [1, len(unit) + 1, 3 * len(unit) + 2]:
                seq = ''.join(rng.choices('ACGT', k=length))
                a = tract_context.fixed_phase_scores(seq, 17, unit)
                b = m._residue_fixed_phase_scores(seq, 17, unit)
                self.assertEqual(a['scores'], b['scores'])
                self.assertEqual(a['all_best_ties'], b['all_best_ties'])

    def test_long_supplied_units_supported_with_explicit_complete_query(self):
        for size in [65, 256]:
            unit = 'A' * (size - 1) + 'G'
            x = m.evaluate_unit_specificity(unit * 3, unit)
            self.assertEqual(x['primitive_period_bp'], size)
            self.assertEqual(x['status'], 'unit_specific_evidence')
            self.assertEqual(len(x['fixed_phase']['scores']), 2 * size)

    def test_homopolymer_edit_baseline_matches_existing_dp(self):
        rng = random.Random(44701)
        for _ in range(20):
            seq = ''.join(rng.choices('ACGT', k=73))
            x = m.evaluate_unit_specificity(seq, 'TTAGGG')
            for base, cost in x['majority_baseline']['all_base_homopolymer_costs'].items():
                self.assertEqual(cost, minimum_periodic_edit_score(seq, base)['edit_cost'])

    def test_json_has_no_nan_or_unreported_positive_claim(self):
        x = m.evaluate_unit_specificity('A' * 100, 'A' * 14 + 'G')
        json.dumps(x, allow_nan=False)
        self.assertFalse(x['changes_native_classification'])
        self.assertFalse(x['biological_absence_or_novelty_inference_permitted'])


class CLITests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.folder = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def prepare(self, content=None, gzip_source=False):
        data = content or ('>other comment\nACGT\n>selected original header\n' + 'TTAGGG' * 30 + '\n')
        path = self.folder / ('reads.fa.gz' if gzip_source else 'reads.fa')
        path.write_bytes(gzip.compress(data.encode()) if gzip_source else data.encode())
        plan = {'schema_version': 1, 'case_id': 'tiny_test', 'data_role': 'development_previously_seen',
                'source': {'path': path.name},
                'expected_source': {'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'size_bytes': path.stat().st_size, 'records': 2, 'bases': 184},
                'regions': [{'id': 'test', 'record_id': 'selected', 'record_length': 180, 'start0': 0, 'end0': 180, 'strand': '+', 'motif': 'TTAGGG'}]}
        return path, plan

    def invoke(self, plan, extra=()):
        path = self.folder / 'plan.json'
        path.write_text(json.dumps(plan))
        output = self.folder / 'output'
        result = subprocess.run([sys.executable, str(SCRIPTS / 'assess_repeat_unit_specificity.py'), '--plan', str(path), '--output', str(output), *extra], text=True, capture_output=True)
        return result, output

    def test_actual_fasta_complete_source_provenance(self):
        source, plan = self.prepare()
        run, output = self.invoke(plan, ['--memory-limit-mib', '256'])
        self.assertEqual(run.returncode, 0, run.stderr)
        summary = json.loads((output / 'summary.json').read_text())
        self.assertEqual((summary['input_records'], summary['input_bases']), (2, 184))
        self.assertEqual(summary['source']['sha256'], plan['expected_source']['sha256'])
        row = json.loads((output / 'details.jsonl').read_text())
        self.assertEqual(row['source_record_ordinal'], 2)
        self.assertEqual(row['analysis']['status'], 'unit_specific_evidence')
        self.assertEqual(row['full_record_normalized_sequence_sha256'], hashlib.sha256(('TTAGGG' * 30).encode()).hexdigest())
        self.assertFalse(row['full_header_preserved'])
        self.assertEqual(row['registered_strand'], '+')
        self.assertIn('both motif strands', row['registered_strand_role'])
        self.assertTrue(summary['registered_strand_is_metadata_only'])
        self.assertFalse(summary['changes_native_classification'])

    def test_wrapped_fastq_quality_and_original_coordinates(self):
        seq = 'CCC' + 'TTAGGG' * 30 + 'AAA'
        fq = '@other\nACGT\n+\n!!!!\n@selected comment\n' + seq[:80] + '\n' + seq[80:] + '\n+selected\n' + 'I' * 90 + '\n' + 'I' * (len(seq) - 90) + '\n'
        source, plan = self.prepare(fq, gzip_source=True)
        plan['expected_source']['bases'] = 190
        plan['regions'][0].update(record_length=186, start0=3, end0=183)
        run, output = self.invoke(plan)
        self.assertEqual(run.returncode, 0, run.stderr)
        row = json.loads((output / 'details.jsonl').read_text())
        self.assertEqual((row['analysis']['start0'], row['analysis']['end0']), (3, 183))
        self.assertEqual(row['quality']['bases'], 180)
        self.assertEqual(row['quality']['mean'], 40)

    def test_bad_source_hash_fails_without_published_output(self):
        source, plan = self.prepare()
        plan['expected_source']['sha256'] = '0' * 64
        run, output = self.invoke(plan)
        self.assertNotEqual(run.returncode, 0)
        self.assertFalse(output.exists())

    def test_bad_record_counts_fail_even_after_selected_hit(self):
        source, plan = self.prepare()
        plan['expected_source']['records'] = 1
        run, output = self.invoke(plan)
        self.assertNotEqual(run.returncode, 0)
        self.assertFalse(output.exists())

    def test_bad_gzip_crc_rejected_even_with_rebound_sha(self):
        source, plan = self.prepare(gzip_source=True)
        raw = bytearray(source.read_bytes())
        raw[-8] ^= 1
        source.write_bytes(raw)
        plan['expected_source']['sha256'] = hashlib.sha256(raw).hexdigest()
        run, output = self.invoke(plan)
        self.assertNotEqual(run.returncode, 0)
        self.assertFalse(output.exists())

    def test_selected_duplicate_or_missing_fails(self):
        for duplicate in [True, False]:
            with self.subTest(duplicate=duplicate):
                source, plan = self.prepare(('>selected\n' + 'TTAGGG' * 30 + '\n') * 2 if duplicate else '>other\nACGT\n>absent\n' + 'TTAGGG' * 30 + '\n')
                if duplicate:
                    plan['expected_source']['bases'] = 360
                run, output = self.invoke(plan)
                self.assertNotEqual(run.returncode, 0)
                self.assertFalse(output.exists())

    def test_not_evaluated_interval_retained_in_output(self):
        source, plan = self.prepare()
        plan['config'] = {'max_region_bp': 100}
        run, output = self.invoke(plan)
        self.assertEqual(run.returncode, 0, run.stderr)
        row = json.loads((output / 'details.jsonl').read_text())
        self.assertEqual(row['analysis']['status'], 'not_evaluated')
        self.assertEqual(row['analysis']['query_bp'], 180)

    def test_total_work_limit_fails_before_source_scoring(self):
        source, plan = self.prepare()
        run, output = self.invoke(plan, ['--max-total-dp-cells', '10'])
        self.assertNotEqual(run.returncode, 0)
        self.assertFalse(output.exists())

    def test_no_output_overwrite(self):
        source, plan = self.prepare()
        run, output = self.invoke(plan)
        self.assertEqual(run.returncode, 0, run.stderr)
        prior = (output / 'summary.json').read_bytes()
        run, output = self.invoke(plan)
        self.assertNotEqual(run.returncode, 0)
        self.assertEqual((output / 'summary.json').read_bytes(), prior)


if __name__ == '__main__':
    unittest.main()
