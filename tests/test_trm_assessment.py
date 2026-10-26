"""Semantic tests of evidence interpretation, independent of scanner tests."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))
import trm_assessment as interpretation
from evidence_batch import canonical


def reference(motif="TTAGGG", evidence="experimental"):
    return {"motif": motif, "evidence": evidence, "lineage": "example species",
            "source": "https://example.test/direct-motif-evidence"}


def raw(motif="TTAGGG", reads=100, full=60, repeat_bases=12000, span=150):
    return {"motif": motif, "reads": reads, "full_repeat_reads": full,
            "repeat_bases": repeat_bases, "max_tract": span}


def targeted(motif="TTAGGG", lf=12, lr=0, rf=0, rr=12, junction_reads=None,
             reads=100, full=60, repeat_bases=12000, span=150):
    if junction_reads is None:
        junction_reads = lf + lr + rf + rr
    return {"motif": canonical(motif), "reads": reads, "repeat_bases": repeat_bases,
            "max_tract": span, "junction_evidence": {
                "left_forward_junction_reads": lf, "left_reverse_junction_reads": lr,
                "right_forward_junction_reads": rf, "right_reverse_junction_reads": rr,
                "junction_reads": junction_reads, "full_repeat_records": full}}


def assess(candidate=None, target=None, catalogue=None, expected=None):
    return interpretation.assess_candidate(candidate or raw(), target,
        interpretation.catalog_index(catalogue if catalogue is not None else [reference()]),
        set(expected or []), 150)


def classify(candidates, status="unassessable_library_or_genome", **options):
    adequacy = {"status": status}
    if status == 'conditional_power_supported':
        # A conditional claim must identify the period models that were run.
        adequacy.update(power={'periods': [{'period': k} for k in range(3, 33)]},
                        detector_scope={'min_k': 3, 'max_k': 32})
    return interpretation.classify(candidates, adequacy, **options)[0]


def validation(target, unique=24, repeats=1, internal=0):
    """Synthetic evidence contract; each signature votes once regardless of PCR count."""
    q = lambda n: {'left_forward': n, 'left_reverse': 0, 'right_forward': 0,
                   'right_reverse': 0, 'ambiguous': 0}
    observed = unique * repeats
    target['junction_evidence'].update(junction_reads=observed,
        left_forward_junction_reads=observed//2, right_reverse_junction_reads=observed-observed//2,
        left_reverse_junction_reads=0, right_forward_junction_reads=0)
    target['reads'] = max(target['reads'], observed + target['junction_evidence']['full_repeat_records'])
    flanks = ['ACGTACGTACGTACGTACGT', 'TGCATGCATGCATGCATGCA']
    target['junction_evidence']['validation'] = {
        'schema_version': 1, 'measured': True, 'clean_two_flank_repeat_records': internal,
        'junction_read_signatures': {'capacity': 128, 'observed_records': observed,
            'distinct_lower_bound': unique, 'saturated': False,
            'entries': [{'signature': f'{i:032x}', 'records': repeats, 'normalized_quadrant_counts': q(repeats)}
                        for i in range(unique)]},
        'flank_contexts': {'capacity': 128, 'observed_events': observed,
            'distinct_lower_bound': 2, 'saturated': False, 'flank_bases': 20,
            'entries': [{'signature': str(i), 'sequence': seq, 'representative_flank': seq,
                         'events': observed//2 if i == 0 else observed-observed//2,
                         'normalized_quadrant_counts': q(observed//2 if i == 0 else observed-observed//2)}
                        for i, seq in enumerate(flanks)]}}
    return target


class InterpretationTests(unittest.TestCase):
    def test_published_g_strand_requires_left_c_rich_right_g_rich(self):
        candidate = assess(target=targeted())
        self.assertEqual(candidate["canonical_family"], "AACCCT")
        self.assertEqual(candidate["motif"], "TTAGGG")
        self.assertEqual(candidate["directionality"]["status"], "consistent_both_ends")
        self.assertEqual(candidate["directionality"]["left_consistent"], 12)
        self.assertEqual(candidate["directionality"]["right_consistent"], 12)
        self.assertEqual(candidate["candidate_status"], "known_trm_pattern_supported")
        self.assertEqual(classify([candidate]), "known_trm_detected")

    def test_opposite_quadrants_are_conflict_for_published_orientation(self):
        candidate = assess(target=targeted(lf=0, lr=12, rf=12, rr=0))
        self.assertEqual(candidate["directionality"]["status"], "orientation_conflict")
        self.assertFalse(candidate["known_detected"])
        self.assertEqual(candidate["candidate_status"], "weak_or_conflicting_known_repeat")

    def test_motif_rotation_reverse_complement_and_repeated_unit_invariance(self):
        expected = assess(target=targeted())
        for motif in ("TTAGGG", "TAGGGT", "CCCTAA", "TTAGGGTTAGGG"):
            candidate = assess(raw(motif), targeted(motif))
            self.assertEqual(candidate, expected)

    def test_reverse_complementing_reads_swaps_consistent_ends_not_conclusion(self):
        a = assess(target=targeted(lf=8, lr=1, rf=2, rr=12))
        b = assess(target=targeted(lf=12, lr=2, rf=1, rr=8))
        self.assertEqual(a["candidate_status"], b["candidate_status"])
        self.assertEqual(a["directionality"]["consistent_fraction"], b["directionality"]["consistent_fraction"])
        self.assertEqual(a["directionality"]["left_consistent"], b["directionality"]["right_consistent"])

    def test_full_repeat_only_known_signal_remains_boundary_unresolved(self):
        candidate = assess(target=targeted(lf=0, rr=0))
        self.assertTrue(candidate["known_detected"])
        self.assertEqual(candidate["candidate_status"], "known_repeat_detected_boundary_unresolved")
        self.assertEqual(candidate["directionality"]["status"], "no_junction_support")
        self.assertEqual(classify([candidate]), "known_repeat_detected_boundary_unresolved")
        self.assertEqual(interpretation.classify([candidate], {"status": "unassessable_library_or_genome"})[2], "review_boundary_support")

    def test_unmeasured_geometry_is_not_a_negative_test(self):
        candidate = assess()
        self.assertEqual(candidate["directionality"]["status"], "not_measured")
        self.assertIsNone(candidate["unique_junction_reads"])
        self.assertTrue(any("not a negative test" in note for note in candidate["notes"]))
        self.assertEqual(classify([candidate]), "known_repeat_detected_boundary_unresolved")

    def test_unmeasured_direction_counts_remain_distinct_from_measured_zero(self):
        missing = assess()
        zero = assess(target=targeted(lf=0, rr=0))
        for field in ('consistent_junctions', 'opposite_junctions', 'left_consistent',
                      'right_consistent', 'oriented_junction_events', 'unique_junction_reads'):
            with self.subTest(field=field):
                self.assertIsNone(missing['directionality'][field])
                self.assertEqual(zero['directionality'][field], 0)
        self.assertEqual(missing['directionality']['status'], 'not_measured')
        self.assertEqual(zero['directionality']['status'], 'no_junction_support')
        for candidate in (missing, zero):
            self.assertEqual(candidate['ranking_score'], 12000)
            self.assertFalse(candidate['geometry_nomination'])
            self.assertFalse(candidate['novel_candidate'])
            self.assertEqual(classify([candidate]), 'known_repeat_detected_boundary_unresolved')

    def test_unknown_full_repeat_only_does_not_become_novel_trm(self):
        candidate = assess(raw("AAT"), targeted("AAT", lf=0, rr=0), catalogue=[])
        self.assertEqual(candidate["candidate_status"], "uncatalogued_repeat")
        self.assertFalse(candidate["novel_candidate"])

    def test_unknown_family_can_qualify_from_geometry_without_identity_blacklist(self):
        candidate = assess(raw("AAT"), targeted("AAT"), catalogue=[])
        self.assertTrue(candidate["geometry_nomination"])
        self.assertFalse(candidate["novel_candidate"])
        self.assertEqual(candidate["candidate_status"], "novel_geometry_nomination")
        self.assertEqual(candidate["directionality"]["orientation_basis"], "inferred_from_junctions_not_a_published_orientation")
        self.assertEqual(classify([candidate]), "novel_geometry_nomination")
        self.assertEqual(candidate['confirmation_readiness']['status'], 'unmeasured')
        self.assertIsNone(candidate['confirmation_readiness']['clean_two_flank_repeat_records'])

    def test_palindromic_family_cannot_supply_orientation(self):
        candidate = assess(raw("ACGT"), targeted("ACGT"), catalogue=[])
        self.assertEqual(candidate["directionality"]["status"], "strand_ambiguous")
        self.assertFalse(candidate["novel_candidate"])

    def test_two_ends_on_five_records_are_not_ten_junction_reads(self):
        candidate = assess(raw("AAC"), targeted("AAC", lf=5, rr=5, junction_reads=5), catalogue=[])
        self.assertEqual(candidate["directionality"]["oriented_junction_events"], 10)
        self.assertEqual(candidate["unique_junction_reads"], 5)
        self.assertFalse(candidate["novel_candidate"])
        self.assertEqual(candidate["full_to_junction_ratio"], 12)

    def test_missing_unique_junction_count_cannot_pass_novel_read_gate(self):
        target = targeted("AAC")
        del target["junction_evidence"]["junction_reads"]
        candidate = assess(raw("AAC"), target, catalogue=[])
        self.assertEqual(candidate["directionality"]["status"], "junction_record_count_unavailable")
        self.assertFalse(candidate["novel_candidate"])

    def test_one_end_only_does_not_establish_bidirectional_boundary_pattern(self):
        candidate = assess(raw("AAC"), targeted("AAC", lf=24, rr=0), catalogue=[])
        self.assertEqual(candidate["directionality"]["status"], "one_end_only")
        self.assertFalse(candidate["novel_candidate"])

    def test_orientation_wilson_quantity_is_descriptive_not_pvalue_or_fdr(self):
        candidate = assess(target=targeted())
        direction = candidate["directionality"]
        self.assertIn("not p-values or FDR", direction["interpretation"])
        self.assertIn("events can overlap", direction["interpretation"])
        self.assertNotIn("p_value", direction)
        self.assertNotIn("fdr", direction)

    def test_comparative_reference_is_published_candidate_not_novel_or_established(self):
        candidate = assess(raw("TTATTGGG"), targeted("TTATTGGG"),
                           catalogue=[reference("TTATTGGG", "comparative")])
        self.assertEqual(candidate["reference_status"], "published_candidate")
        self.assertEqual(candidate["candidate_status"], "published_candidate_pattern_supported")
        self.assertFalse(candidate["novel_candidate"])
        self.assertFalse(candidate["known_detected"])
        self.assertEqual(classify([candidate]), "published_candidate_supported")

    def test_nearby_reference_anchor_never_relabels_actual_motif(self):
        observed = raw("TTAGAC")
        observed["telomere_prior"] = {"matched_reference": "TTAGGG", "eligible": True}
        candidate = assess(observed, targeted("TTAGAC"), catalogue=[reference("TTAGGG")])
        self.assertEqual(candidate["reference_status"], "not_in_catalogue")
        self.assertEqual(canonical(candidate["motif"]), canonical("TTAGAC"))
        self.assertEqual(candidate["published_motifs"], [])
        self.assertNotEqual(candidate["motif"], "TTAGGG")

    def test_credible_known_boundary_uncertainty_precedes_novel_primary_without_erasure(self):
        known = assess(raw(reads=10, full=5, repeat_bases=750), targeted(lf=0, rr=0, reads=10, full=5, repeat_bases=750))
        novel = assess(raw("AAC", reads=10000, full=5000, repeat_bases=750000),
                       targeted("AAC", reads=10000, full=5000, repeat_bases=750000), catalogue=[])
        self.assertGreater(novel["ranking_score"], known["ranking_score"])
        self.assertEqual(classify([known, novel]), "known_repeat_detected_boundary_unresolved")
        self.assertTrue(novel["geometry_nomination"])
        self.assertFalse(novel["novel_candidate"])

    def test_minor_known_pattern_does_not_override_dominant_known_family(self):
        full_only = assess(raw(reads=10000, full=5000, repeat_bases=750000),
                           targeted(lf=0, rr=0, reads=10000, full=5000, repeat_bases=750000))
        known_pattern = assess(raw("TTAGGC", reads=50, full=10, repeat_bases=6000),
                               targeted("TTAGGC", reads=50, full=10, repeat_bases=6000), catalogue=[reference("TTAGGC")])
        self.assertGreater(full_only["ranking_score"], known_pattern["ranking_score"])
        self.assertEqual(classify([full_only, known_pattern]), "known_repeat_detected_boundary_unresolved")
        self.assertEqual(classify([known_pattern, full_only]), "known_repeat_detected_boundary_unresolved")

    def test_partial_known_directionality_is_not_opposite_conflict(self):
        candidate = assess(raw("TTAGGT", reads=13739, full=9080, repeat_bases=1581657, span=125),
            targeted("TTAGGT", lf=387, lr=270, rf=21, rr=198, junction_reads=860,
                     reads=13739, full=9080, repeat_bases=1581657, span=125),
            catalogue=[reference("TTAGGT", "published_terminal_long_read")])
        self.assertAlmostEqual(candidate["directionality"]["consistent_fraction"], 585/876)
        self.assertEqual(candidate["directionality"]["status"], "partial_directional_support")
        self.assertTrue(candidate["known_detected"])
        self.assertEqual(candidate["candidate_status"], "known_repeat_detected_boundary_unresolved")
        self.assertEqual(classify([candidate]), "known_repeat_detected_boundary_unresolved")

    def test_balanced_junctions_are_mixed_not_an_opposite_orientation_conflict(self):
        candidate = assess(target=targeted(lf=10, lr=10, rf=10, rr=10))
        self.assertEqual(candidate["directionality"]["status"], "mixed_directionality")
        self.assertTrue(candidate["known_detected"])
        self.assertEqual(candidate["candidate_status"], "known_repeat_detected_boundary_unresolved")

    def test_expected_motif_bonus_is_bounded_and_never_changes_observed_identity(self):
        expected = assess(raw(reads=10, full=5, repeat_bases=100), expected={canonical("TTAGGG")})
        variant = assess(raw("TTAGAC", reads=100, full=60, repeat_bases=1000), catalogue=[])
        self.assertEqual(expected["ranking_score"], 125)
        self.assertGreater(variant["ranking_score"], expected["ranking_score"])
        self.assertEqual(canonical(variant["motif"]), canonical("TTAGAC"))

    def test_targeted_positive_not_in_original_discovery_is_retained(self):
        candidate = assess(raw(reads=0, full=0, repeat_bases=0, span=0), targeted())
        self.assertEqual(candidate["evidence_count_source"], "targeted_only")
        self.assertEqual(candidate["original_discovery_reads"], 0)
        self.assertEqual(candidate["reads"], 100)
        self.assertEqual(candidate["full_repeat_reads"], 60)
        self.assertEqual(candidate["repeat_bases"], 12000)
        self.assertEqual(candidate["max_tract_bp"], 150)
        self.assertEqual(candidate["candidate_status"], "known_trm_pattern_supported")

    def test_complete_measured_counts_are_used_together_and_original_values_preserved(self):
        candidate = assess(raw(reads=25, full=12, repeat_bases=2250, span=100),
                           targeted(reads=30, full=15, repeat_bases=3600, span=150, lf=0, rr=0))
        self.assertEqual(candidate["reads"], 30)
        self.assertEqual(candidate["full_repeat_reads"], 15)
        self.assertEqual(candidate["repeat_bases"], 3600)
        self.assertEqual(candidate["max_tract_bp"], 150)
        self.assertEqual(candidate["targeted_reads"], 30)
        self.assertEqual(candidate["targeted_full_repeat_reads"], 15)
        self.assertEqual(candidate["original_discovery_reads"], 25)
        self.assertEqual(candidate["original_discovery_full_repeat_reads"], 12)
        self.assertEqual(candidate["original_discovery_repeat_bases"], 2250)
        self.assertEqual(candidate["original_discovery_max_tract_bp"], 100)
        self.assertEqual(candidate["evidence_count_source"], "measured_evidence")
        self.assertTrue(candidate["measured_count_set_complete"])
        self.assertTrue(any("counts differ" in note for note in candidate["notes"]))

    def test_zero_vs_one_legacy_read_cannot_suppress_the_same_complete_measured_signal(self):
        result = []
        for old_reads in (0, 1):
            candidate = assess(raw(reads=old_reads, full=0, repeat_bases=24*old_reads, span=24 if old_reads else 0),
                               targeted(reads=10, full=5, repeat_bases=1500, lf=0, rr=0))
            self.assertEqual((candidate["reads"], candidate["full_repeat_reads"]), (10, 5))
            self.assertEqual(candidate["candidate_status"], "known_repeat_detected_boundary_unresolved")
            result.append(classify([candidate]))
        self.assertEqual(result[0], result[1])

    def test_lower_complete_measured_counts_are_not_replaced_by_legacy_maxima(self):
        candidate = assess(raw(reads=100, full=60, repeat_bases=12000, span=150),
                           targeted(reads=8, full=3, repeat_bases=480, span=90, lf=0, rr=0))
        self.assertEqual((candidate["reads"], candidate["full_repeat_reads"], candidate["repeat_bases"], candidate["max_tract_bp"]),
                         (8, 3, 480, 90))
        self.assertFalse(candidate["known_detected"])
        self.assertEqual(candidate["original_discovery_reads"], 100)

    def test_geometry_only_or_incomplete_count_sets_fall_back_without_mixing(self):
        full_target = targeted(reads=30, full=15, repeat_bases=3600, span=150)
        geometry_only = {"motif": full_target["motif"], "junction_evidence": full_target["junction_evidence"]}
        partials = [geometry_only]
        for omitted in ("reads", "repeat_bases", "max_tract"):
            partial = copy.deepcopy(full_target)
            del partial[omitted]
            partials.append(partial)
        partial = copy.deepcopy(full_target)
        del partial["junction_evidence"]["full_repeat_records"]
        partials.append(partial)
        for target in partials:
            with self.subTest(target=target):
                candidate = assess(raw(reads=25, full=12, repeat_bases=2250, span=100), target)
                self.assertEqual((candidate["reads"], candidate["full_repeat_reads"], candidate["repeat_bases"], candidate["max_tract_bp"]),
                                 (25, 12, 2250, 100))
                self.assertEqual(candidate["evidence_count_source"], "original_discovery")
                self.assertFalse(candidate["measured_count_set_complete"])
                self.assertEqual(candidate["directionality"]["status"], "consistent_both_ends")

    def test_zero_measured_positive_original_is_explicit_measurement_disagreement(self):
        candidate = assess(raw(reads=25, full=12, repeat_bases=2250, span=100),
                           targeted(reads=0, full=0, repeat_bases=0, span=0, lf=0, rr=0))
        self.assertEqual(candidate["reads"], 0)
        self.assertEqual(candidate["original_discovery_reads"], 25)
        self.assertEqual(candidate["original_discovery_full_repeat_reads"], 12)
        self.assertTrue(candidate["measurement_disagreement"])
        self.assertEqual(candidate["measurement_disagreement_reason"], "positive_original_zero_measured")
        self.assertEqual(candidate["evidence_count_source"], "measured_evidence")
        self.assertTrue(any("Neither result proves motif absence or ALT" in note for note in candidate["notes"]))

    def test_invalid_complete_measured_counts_are_rejected(self):
        for target in (targeted(reads=-1), targeted(reads=1.5), targeted(reads=True), targeted(reads=5, full=6)):
            with self.subTest(target=target), self.assertRaisesRegex(ValueError, 'Measured'):
                assess(raw(), target)

    def test_unresolved_measurement_disagreement_blocks_conditional_nondetection(self):
        disputed = assess(raw(), targeted(reads=0, full=0, repeat_bases=0, span=0,
                                          lf=0, lr=0, rf=0, rr=0))
        self.assertEqual(classify([disputed], 'conditional_power_supported'),
                         'detection_inconclusive_measurement')
        positive = assess(raw('TTAGGC'), targeted('TTAGGC'), catalogue=[reference('TTAGGC')])
        self.assertEqual(classify([disputed, positive], 'conditional_power_supported'),
                         'known_trm_detected')

    def test_long_observed_tract_is_not_denied_based_on_short_mean_length(self):
        candidate = assess(raw(span=600), targeted(span=600))
        self.assertFalse(any("no observed contiguous" in note for note in candidate["notes"]))

    def test_adequacy_uncertainty_categories_do_not_claim_true_absence(self):
        states = {"observed_qc_failure": "detection_inconclusive_technical",
                  "insufficient_data": "detection_inconclusive_depth",
                  "unassessable_library_or_genome": "detection_inconclusive",
                  "conditional_power_supported": "canonical_not_detected_conditional_power"}
        for adequacy, expected in states.items():
            with self.subTest(adequacy=adequacy):
                self.assertEqual(classify([], adequacy), expected)
                self.assertNotIn("absence", expected)
                self.assertNotIn("alt", expected)

    def test_invalid_qc_precedes_positive_call_but_unassessable_library_does_not(self):
        known = assess(target=targeted())
        self.assertEqual(classify([known], "observed_qc_failure"), "detection_inconclusive_technical")
        self.assertEqual(classify([known], "unassessable_library_or_genome"), "known_trm_detected")
        self.assertEqual(classify([known], analysis_complete=False), "analysis_failed")
        self.assertEqual(classify([known], biological_scope=False), "taxonomic_scope_unresolved")

    def test_inconsistent_family_or_junction_counts_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "family"):
            assess(raw(), targeted("TTAGGC"))
        with self.assertRaisesRegex(ValueError, "junction"):
            assess(raw(), targeted(lf=12, rr=0, junction_reads=5))

    def test_assessment_does_not_mutate_source_evidence(self):
        candidate, target = raw(), targeted()
        originals = copy.deepcopy((candidate, target))
        assess(candidate, target)
        self.assertEqual((candidate, target), originals)


class ObservationStageTests(unittest.TestCase):
    assignment = {'reference_scope': 'species', 'reference_status': 'supported',
                  'repeat_architecture': 'exact_short_unit'}

    def observed(self, reads=100, full=3, evidence='experimental', expected=True,
                 target_options=None):
        return assess(raw(reads=reads, full=full),
                      targeted(reads=reads, full=full, lf=0, rr=0, **(target_options or {})),
                      catalogue=[reference(evidence=evidence)],
                      expected={canonical('TTAGGG')} if expected else set())

    def novel(self):
        return assess(raw('AAC', repeat_bases=900000),
                      validation(targeted('AAC', repeat_bases=900000)), catalogue=[])

    def test_observation_is_separate_from_full_repeat_pattern_support(self):
        for full in (0, 3):
            with self.subTest(full=full):
                candidate = self.observed(full=full)
                self.assertEqual(candidate['sequence_observation']['status'], 'sequence_observed')
                self.assertTrue(candidate['known_observed'])
                self.assertTrue(candidate['reference_observed'])
                self.assertFalse(candidate['known_detected'])
                self.assertEqual(candidate['full_repeat_reads'], full)
                self.assertEqual(classify([candidate]), 'known_repeat_observed_pattern_unresolved')
                self.assertEqual(interpretation.POLICY['minimum_full_repeat_reads'], 5)

    def test_no_positive_counts_do_not_become_observation_from_a_reference(self):
        candidate = assess(raw(reads=0, full=0, repeat_bases=0, span=0),
                           targeted(reads=0, full=0, repeat_bases=0, span=0, lf=0, rr=0),
                           expected={canonical('TTAGGG')})
        self.assertFalse(candidate['sequence_observation']['observed'])
        self.assertFalse(candidate['reference_observed'])
        self.assertEqual(candidate['sequence_observation']['status'], 'no_sequence_observed')
        code = classify([candidate], 'conditional_power_supported', reference_assignment=self.assignment)
        self.assertEqual(code, 'canonical_not_detected_conditional_power')
        self.assertIsNone(interpretation.select_representative([candidate], code, self.assignment))

    def test_expected_observed_reference_precedes_unrelated_novel_without_erasure(self):
        expected, novel = self.observed(), self.novel()
        self.assertGreater(novel['ranking_score'], expected['ranking_score'])
        self.assertTrue(novel['novel_candidate'])
        originals = copy.deepcopy([expected, novel])
        for candidates in ([expected, novel], [novel, expected]):
            result = interpretation.classify(candidates, {'status': 'conditional_power_supported'},
                                             reference_assignment=self.assignment)
            self.assertEqual(result[0], 'known_repeat_observed_pattern_unresolved')
            self.assertEqual(interpretation.select_representative(candidates, result, self.assignment), expected)
        self.assertEqual([expected, novel], originals)

    def test_isolate_reference_has_the_same_observation_priority(self):
        expected, novel = self.observed(), self.novel()
        assignment = dict(self.assignment, reference_scope='isolate')
        self.assertEqual(classify([expected, novel], reference_assignment=assignment),
                         'known_repeat_observed_pattern_unresolved')

    def test_clade_expectation_does_not_receive_exact_species_priority(self):
        expected, novel = self.observed(), self.novel()
        assignment = dict(self.assignment, reference_scope='clade')
        self.assertEqual(classify([expected, novel], reference_assignment=assignment),
                         'novel_provisional_candidate')

    def test_global_weak_reference_does_not_suppress_a_strong_novel_candidate(self):
        for reads in (1, 100):
            with self.subTest(reads=reads):
                incidental, novel = self.observed(reads=reads, full=0, expected=False), self.novel()
                code = classify([incidental, novel])
                self.assertEqual(code, 'novel_provisional_candidate')
                self.assertEqual(interpretation.select_representative([incidental, novel], code), novel)
                self.assertEqual(incidental['sequence_observation']['observed'], True)

    def test_template_only_reference_is_observed_without_becoming_established(self):
        candidate = self.observed(evidence='telomerase_template')
        self.assertEqual(candidate['reference_status'], 'published_candidate')
        self.assertTrue(candidate['reference_observed'])
        self.assertFalse(candidate['known_observed'])
        self.assertFalse(candidate['known_detected'])
        result = interpretation.classify([candidate, self.novel()], {'status': 'conditional_power_supported'},
                                         reference_assignment=dict(self.assignment, reference_status='candidate'))
        self.assertEqual(result[0], 'published_candidate_observed_pattern_unresolved')
        self.assertEqual(result[2], 'validate_published_candidate')
        self.assertEqual(interpretation.select_representative([candidate, self.novel()], result, self.assignment), candidate)

    def test_candidate_reference_status_cannot_be_promoted_by_evidence_token(self):
        ref = dict(reference(), reference_status='candidate')
        candidate = assess(target=targeted(), catalogue=[ref])
        self.assertEqual(candidate['reference_status'], 'published_candidate')
        self.assertFalse(candidate['known_detected'])
        self.assertFalse(candidate['known_observed'])
        self.assertTrue(candidate['reference_observed'])
        self.assertEqual(classify([candidate]), 'published_candidate_supported')

    def test_observed_expected_orientation_conflict_never_becomes_pattern_support(self):
        candidate = assess(target=targeted(lf=0, lr=12, rf=12, rr=0), expected={canonical('TTAGGG')})
        self.assertEqual(candidate['directionality']['status'], 'orientation_conflict')
        self.assertFalse(candidate['known_detected'])
        self.assertTrue(candidate['known_observed'])
        self.assertEqual(classify([candidate, self.novel()], reference_assignment=self.assignment),
                         'known_repeat_observed_pattern_unresolved')

    def test_library_and_qc_limits_preserve_separate_sequence_observation(self):
        expected, novel = self.observed(), self.novel()
        for status, kwargs, expected_code in (
                ('conditional_power_supported', {'library': {'limited': True}}, 'library_context_limited'),
                ('observed_qc_failure', {}, 'detection_inconclusive_technical'),
                ('conditional_power_supported', {'biological_scope': False}, 'taxonomic_scope_unresolved')):
            with self.subTest(code=expected_code):
                result = interpretation.classify([novel, expected], {'status': status},
                    reference_assignment=self.assignment, **kwargs)
                self.assertEqual(result[0], expected_code)
                self.assertEqual(interpretation.select_representative([novel, expected], result, self.assignment), expected)
                self.assertEqual(expected['sequence_observation']['status'], 'sequence_observed')

    def test_expected_trace_is_preserved_before_nondetection_but_does_not_override_novel(self):
        trace = self.observed(reads=3, full=0)
        self.assertEqual(trace['sequence_observation']['status'], 'trace_sequence_observed')
        self.assertFalse(trace['sequence_observation']['representative_eligible'])
        code = classify([trace], 'conditional_power_supported', reference_assignment=self.assignment)
        self.assertEqual(code, 'published_repeat_trace_observed')
        self.assertEqual(interpretation.select_representative([trace], code, self.assignment), trace)
        self.assertEqual(classify([trace, self.novel()], reference_assignment=self.assignment),
                         'novel_provisional_candidate')

    def test_weak_current_counts_are_not_promoted_by_larger_historical_counts(self):
        candidate = assess(raw(reads=100, full=60),
                           targeted(reads=3, full=0, lf=0, rr=0), expected={canonical('TTAGGG')})
        observation = candidate['sequence_observation']
        self.assertEqual(observation['current_reads'], 3)
        self.assertEqual(observation['original_discovery_reads'], 100)
        self.assertEqual(observation['status'], 'trace_sequence_observed')
        self.assertFalse(observation['representative_eligible'])
        self.assertEqual(classify([candidate], reference_assignment=self.assignment), 'published_repeat_trace_observed')

    def test_measurement_disagreement_retains_observation_without_positive_promotion(self):
        candidate = assess(raw(reads=100, full=60),
                           targeted(reads=0, full=0, repeat_bases=0, span=0, lf=0, rr=0),
                           expected={canonical('TTAGGG')})
        observation = candidate['sequence_observation']
        self.assertEqual(observation['status'], 'measurement_disagreement')
        self.assertTrue(observation['observed'])
        self.assertFalse(observation['representative_eligible'])
        self.assertEqual(observation['current_reads'], 0)
        code = classify([candidate], 'conditional_power_supported', reference_assignment=self.assignment)
        self.assertEqual(code, 'detection_inconclusive_measurement')
        self.assertEqual(interpretation.select_representative([candidate], code, self.assignment), candidate)

    def test_known_alternative_architecture_is_not_overwritten_by_incidental_known_pattern(self):
        candidate = assess(target=targeted())
        assignment = {'reference_scope': 'species', 'reference_status': 'out_of_scope',
                      'repeat_architecture': 'irregular'}
        code = classify([candidate], reference_assignment=assignment)
        self.assertEqual(code, 'known_alternative_architecture')
        self.assertEqual(interpretation.select_representative([candidate], code, assignment), candidate)
        self.assertTrue(candidate['known_detected'])

    def test_representative_uses_same_dominant_known_family_as_classification(self):
        full_only = assess(raw(reads=10000, full=5000, repeat_bases=750000),
                           targeted(lf=0, rr=0, reads=10000, full=5000, repeat_bases=750000))
        known_pattern = assess(raw('TTAGGC'), targeted('TTAGGC'), catalogue=[reference('TTAGGC')])
        candidates = [known_pattern, full_only]
        code = classify(candidates)
        self.assertEqual(code, 'known_repeat_detected_boundary_unresolved')
        self.assertEqual(interpretation.select_representative(candidates, code), full_only)


class ProvisionalEvidenceTests(unittest.TestCase):
    def candidate(self, target=None):
        return assess(raw('AAT'), target or validation(targeted('AAT')), catalogue=[])

    def test_complete_unique_evidence_promotes_only_to_provisional_without_genome_size(self):
        candidate = self.candidate()
        self.assertTrue(candidate['geometry_nomination'])
        self.assertTrue(candidate['novel_candidate'])
        self.assertEqual(candidate['novel_evidence_level'], 'provisional_candidate')
        self.assertEqual(classify([candidate]), 'novel_provisional_candidate')
        self.assertFalse(candidate['externally_confirmed'])
        self.assertTrue(candidate['confirmation_readiness']['provisional_ready'])

    def test_pcr_multiplicity_does_not_create_distinct_sequence_support(self):
        candidate = self.candidate(validation(targeted('AAT'), unique=2, repeats=100))
        self.assertTrue(candidate['geometry_nomination'])
        self.assertFalse(candidate['novel_candidate'])
        ready = candidate['confirmation_readiness']
        self.assertEqual(ready['distinct_junction_sequences_exact'], 2)
        self.assertFalse(ready['checks']['minimum_distinct_junction_sequences'])

    def test_duplicate_skew_cannot_override_distinct_sequence_direction(self):
        target = validation(targeted('AAT'), unique=24)
        entries = target['junction_evidence']['validation']['junction_read_signatures']['entries']
        for entry in entries[4:]:
            entry['normalized_quadrant_counts'].update(left_forward=0, left_reverse=1)
        candidate = self.candidate(target)
        self.assertTrue(candidate['geometry_nomination'])
        self.assertFalse(candidate['novel_candidate'])
        self.assertAlmostEqual(candidate['confirmation_readiness']['deduplicated_orientation']['consistent_fraction'], 4/24)

    def test_conflicting_quadrants_on_one_signature_do_not_get_two_votes(self):
        target = validation(targeted('AAT'))
        for entry in target['junction_evidence']['validation']['junction_read_signatures']['entries'][:10]:
            entry['normalized_quadrant_counts']['left_reverse'] = 1
        candidate = self.candidate(target)
        dedup = candidate['confirmation_readiness']['deduplicated_orientation']
        self.assertEqual(dedup['ambiguous_or_conflicting_unique_retained'], 10)
        self.assertEqual(dedup['consistent_unique_retained'], 14)
        self.assertFalse(candidate['novel_candidate'])

    def test_reverse_complement_read_quadrants_collapse_to_same_direction_class(self):
        a = validation(targeted('AAT'))
        b = copy.deepcopy(a)
        for entry in b['junction_evidence']['validation']['junction_read_signatures']['entries']:
            entry['normalized_quadrant_counts'].update(left_forward=0, right_reverse=1)
        self.assertEqual(self.candidate(a)['confirmation_readiness']['deduplicated_orientation'],
                         self.candidate(b)['confirmation_readiness']['deduplicated_orientation'])

    def test_near_identical_flanks_merge_into_one_error_tolerant_cluster(self):
        target = validation(targeted('AAT'))
        entries = target['junction_evidence']['validation']['flank_contexts']['entries']
        entries[1]['sequence'] = 'TCGTACGTACGTACGTACGA'
        entries[1]['representative_flank'] = entries[1]['sequence']
        candidate = self.candidate(target)
        self.assertEqual(candidate['confirmation_readiness']['physical_flank_cluster_count'], 1)
        self.assertFalse(candidate['novel_candidate'])

    def test_boundary_offset_representatives_merge_despite_different_20mers(self):
        rep = 'ACGTCGATCGTAGCATGCTAGCAGTTCGATGCCGATCGATCTGACGATCTAGCATGCGTACGAT'
        a, b = rep[:55], rep[7:62]
        self.assertGreater(sum(x != y for x, y in zip(a[:20], b[:20])), 2)
        groups = interpretation.flank_clusters([a[:20], b[:20]], {a[:20]: a, b[:20]: b})
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups, interpretation.flank_clusters([b[:20], a[:20]], {a[:20]: a, b[:20]: b}))

    def test_internal_array_counterevidence_holds_promotion_without_erasing_nomination(self):
        candidate = self.candidate(validation(targeted('AAT'), internal=1))
        self.assertTrue(candidate['geometry_nomination'])
        self.assertFalse(candidate['novel_candidate'])
        self.assertIn('coexist', candidate['confirmation_readiness']['internal_counterevidence_interpretation'])
        self.assertEqual(classify([candidate]), 'novel_geometry_nomination')

    def test_saturated_reads_use_omitted_records_not_subset_fraction_or_cardinality_bound(self):
        target = validation(targeted('AAT'), unique=12)
        sketch = target['junction_evidence']['validation']['junction_read_signatures']
        sketch.update(capacity=10, saturated=True, distinct_lower_bound=11, entries=sketch['entries'][:10])
        candidate = self.candidate(target)
        d = candidate['confirmation_readiness']['deduplicated_orientation']
        self.assertIsNone(d['consistent_fraction'])
        self.assertEqual(d['omitted_records_upper_bound_on_unseen_unique_sequences'], 2)
        self.assertAlmostEqual(d['consistent_fraction_lower_bound'], 10/12)
        self.assertTrue(candidate['novel_candidate'])
        # The same retained subset cannot support a denominator containing 30 unseen records.
        target = validation(targeted('AAT'), unique=40)
        sketch = target['junction_evidence']['validation']['junction_read_signatures']
        sketch.update(capacity=10, saturated=True, distinct_lower_bound=11, entries=sketch['entries'][:10])
        candidate = self.candidate(target)
        self.assertFalse(candidate['novel_candidate'])
        self.assertEqual(candidate['confirmation_readiness']['deduplicated_orientation']['status'], 'unresolved_saturation')

    def test_saturated_flank_clusters_are_not_asserted_as_global_lower_bound(self):
        target = validation(targeted('AAT'))
        sketch = target['junction_evidence']['validation']['flank_contexts']
        sketch.update(capacity=2, saturated=True, distinct_lower_bound=3, observed_events=25)
        candidate = self.candidate(target)
        ready = candidate['confirmation_readiness']
        self.assertEqual(ready['distinct_physical_flanks_lower_bound'], 3)
        self.assertEqual(ready['retained_flank_cluster_count'], 2)
        self.assertIsNone(ready['physical_flank_cluster_count'])
        self.assertFalse(candidate['novel_candidate'])

    def test_missing_partial_validation_is_incomplete_not_zero(self):
        target = validation(targeted('AAT'))
        del target['junction_evidence']['validation']['clean_two_flank_repeat_records']
        ready = self.candidate(target)['confirmation_readiness']
        self.assertEqual(ready['status'], 'incomplete')
        self.assertIsNone(ready['clean_two_flank_repeat_records'])

    def test_invalid_sketch_count_contract_is_rejected(self):
        target = validation(targeted('AAT'))
        target['junction_evidence']['validation']['junction_read_signatures']['observed_records'] += 1
        with self.assertRaisesRegex(ValueError, 'Unsaturated'):
            self.candidate(target)

    def test_exact_species_architecture_precedes_novel_but_does_not_erase_it(self):
        candidate = self.candidate()
        assignment = {'reference_scope': 'species', 'reference_status': 'out_of_scope', 'repeat_architecture': 'irregular'}
        self.assertEqual(classify([candidate], reference_assignment=assignment), 'known_alternative_architecture')
        self.assertTrue(candidate['novel_candidate'])
        assignment['reference_scope'] = 'clade'
        self.assertEqual(classify([candidate], reference_assignment=assignment), 'novel_provisional_candidate')

    def test_library_context_limits_novelty_without_claiming_RNA_from_metadata(self):
        library = interpretation.library_context({'library_strategy': 'WGS', 'library_source': 'GENOMIC',
                                                 'library_construction_protocol': 'RNA PolyA'})
        self.assertTrue(library['limited'])
        self.assertFalse(library['RNA_component_supported'])
        self.assertFalse(library['metadata_contradiction_proves_RNA'])
        self.assertEqual(classify([self.candidate()], library=library), 'library_context_limited')
        self.assertEqual(classify([], 'conditional_power_supported', library=library), 'library_context_limited')
        self.assertEqual(classify([assess(target=targeted())], library=library), 'known_trm_detected')

    def test_full_to_junction_ratio_warning_does_not_impose_universal_length_rule(self):
        target = validation(targeted('AAT', full=5, span=75))
        candidate = self.candidate(target)
        self.assertTrue(candidate['long_array_full_to_junction_warning'])
        self.assertTrue(candidate['novel_candidate'])
        self.assertFalse(candidate['externally_confirmed'])


if __name__ == "__main__":
    unittest.main()
