"""Orthogonal clade-default/variant observations, never a new TRM claim.

Global catalogue membership, a published species/isolate exception, and a
sequence difference from a weak clade default are separate axes. This module
does not alter candidate counts, ranking, biological classification or power.
"""
from __future__ import annotations

from copy import deepcopy

from evidence_batch import canonical
from prepare_taxon_priors import assign, motif_list, status_of
from reference_evidence import SUPPORTED_REFERENCE_EVIDENCE
from trm_assessment import resolved_eukaryotic_scope, sequence_observation


def _families(value):
    return {canonical(motif) for motif in motif_list(value)} if value else set()


def catalogue_status(family, catalogue):
    if catalogue is None:
        return 'not_provided'
    rows = catalogue.get(family, [])
    if not rows:
        return 'not_in_global_catalogue'
    established = any(row.get('evidence') in SUPPORTED_REFERENCE_EVIDENCE
                      and row.get('reference_status') in (None, '', 'supported', 'established') for row in rows)
    return 'globally_established_family' if established else 'globally_catalogued_unestablished_family'


def build_clade_variant_report(candidates, metadata, assignment, taxon_reference_rows, *,
                              global_catalogue=None, representative=None,
                              analysis_complete=True, library=None):
    """Compare relevant observed families with the nearest curated clade context.

    Uses the existing taxonomy-matching implementation on clade rows only.
    More-specific unknown/alternative contexts remain a block on interpretation;
    supported species/isolate exceptions remain known published exceptions.
    Known/reference traces and existing novel nominations remain visible. The
    enormous unqualified microsatellite tail is counted, not copied: its full
    observations remain in the input candidate table.
    """
    metadata, assignment = metadata or {}, assignment or {}
    library = library or {}
    result = {
        'schema_version': 1, 'status': 'not_assessed', 'comparison_permitted': False,
        'clade_context': None, 'default_families': [], 'families': [], 'review_items': [],
        'review_required': False, 'variant_observed': False, 'canonical_observed': False,
        'variant_observed_family_count': 0, 'canonical_observed_family_count': 0,
        'variant_trace_family_count': 0, 'published_exception_family_count': 0,
        'unqualified_observed_families_not_copied': 0,
        'scoped_reference': {key: deepcopy(assignment.get(key)) for key in
                             ('reference_scope', 'reference_status', 'related_motifs',
                              'reference_ids', 'review_reference_ids', 'repeat_architecture')},
        'changes_classification': False, 'changes_candidate_counts_or_ranking': False,
        'changes_adequacy': False, 'biological_novelty_claimed': False,
        'literature_first_report_claimed': False, 'terminal_location_confirmed': False,
        'scope': 'Observed primitive repeat families compared with a weak nearest-clade default; this is neither taxon assignment nor a new telomeric-repeat confirmation.',
        'limitations': [
            'A globally catalogued family can differ from a clade default; these are separate descriptions.',
            'Published species/isolate exceptions remain published exceptions, not novel discoveries.',
            'A clade default is not species ground truth, a phylogenetic reconstruction, or a universal motif.',
            'Short-read sequence and boundary observations do not establish physical chromosome-terminal location.',
            'Unqualified uncatalogued repeats remain in the full candidate table; omission here is not rejection.',
            'Library/taxonomic limits suppress actionable clade interpretation without deleting sequence observations.',
        ]}
    if not analysis_complete:
        result['status'] = 'analysis_incomplete'
        return result
    if not resolved_eukaryotic_scope(metadata) or library.get('identity_conflict_supported'):
        result['status'] = 'taxonomic_scope_unresolved'
        return result
    # Isolate and species assignment should never be replaced with an inferred
    # clade motif. A separate clade-only lookup is a descriptive comparison.
    clade_rows = []
    for original in taxon_reference_rows or []:
        if original.get('scope') != 'clade':
            continue
        row = dict(original)
        row['reference_status'] = status_of(row)
        row['repeat_architecture'] = row.get('repeat_architecture') or 'exact_short_unit'
        row['motifs'] = motif_list(row['related_motifs']) if row.get('related_motifs') else []
        clade_rows.append(row)
    sample = dict(metadata)
    sample.setdefault('sample_id', str(metadata.get('run_accession') or 'sample'))
    # Match the report's validated, whitespace-insensitive taxonomy identity.
    # Preserve the original metadata/provenance outside this local lookup.
    for key in ('sample_id', 'status', 'scientific_name', 'tax_id', 'lineage'):
        if isinstance(sample.get(key), str):
            sample[key] = sample[key].strip()
    clade_assignment = assign([sample], clade_rows)[0][0]
    result['clade_context'] = {key: deepcopy(clade_assignment.get(key)) for key in
                               ('reference_status', 'reference_scope', 'mapping_status',
                                'mapping_basis', 'related_motifs', 'reference_ids', 'reference_evidence',
                                'prior_source', 'review_reference_ids', 'suppressed_reference_ids',
                                'repeat_architecture', 'notes')}
    selected_ids = set((clade_assignment.get('reference_ids') or '').split(','))
    result['clade_context']['taxa'] = [
        {'scientific_name': row['scientific_name'], 'tax_id': row['reference_taxid'],
         'reference_id': row['reference_id']} for row in clade_rows if row['reference_id'] in selected_ids]
    defaults = _families(clade_assignment.get('related_motifs'))
    result['default_families'] = sorted(defaults)
    if clade_assignment.get('reference_status') != 'supported' or not defaults:
        result['status'] = ('no_curated_clade_default' if clade_assignment.get('reference_status') == 'not_curated'
                            else 'nearest_clade_has_no_supported_default')
        return result
    narrower = assignment.get('reference_scope') in {'species', 'isolate', 'genus'}
    blocked_specific = narrower and (
        assignment.get('reference_status') in {'candidate', 'unknown', 'out_of_scope'}
        or bool(set((assignment.get('repeat_architecture') or '').split(','))
                & {'irregular', 'retrotransposon', 'long_complex_family', 'diverse'}))
    if blocked_specific:
        result['status'] = 'specific_reference_blocks_clade_interpretation'
        return result
    supported_scoped = (narrower and assignment.get('reference_status') == 'supported')
    scoped_families = _families(assignment.get('related_motifs')) if supported_scoped else set()
    representative_family = canonical(representative['motif']) if representative and representative.get('motif') else None
    permitted = not library.get('limited')
    result['comparison_permitted'] = permitted
    rows, reviews, seen = [], [], set()
    for candidate in candidates or []:
        observation = sequence_observation(candidate)
        if not observation['observed']:
            continue
        family = canonical(candidate.get('canonical_family') or candidate['motif'])
        if family in seen:
            raise ValueError('Clade comparison requires one assessed row per primitive family; do not combine duplicate measurements')
        seen.add(family)
        reference_observation = candidate.get('reference_status') in {'established', 'published_candidate'}
        nominated = bool(candidate.get('novel_candidate') or candidate.get('geometry_nomination')
                         or candidate.get('candidate_status') == 'published_candidate_pattern_supported')
        if family not in defaults and not reference_observation and not nominated:
            result['unqualified_observed_families_not_copied'] += 1
            continue
        different = family not in defaults
        positive = observation['current_reads'] > 0 and not observation['measurement_disagreement']
        enough = observation['representative_eligible']
        scoped_match = bool(supported_scoped and family in scoped_families)
        published_exception = scoped_match and different and assignment.get('reference_scope') in {'species', 'isolate'}
        if observation['measurement_disagreement']:
            state = 'measurement_disagreement'
        elif not positive:
            state = 'historical_observation_only'
        elif different and not enough:
            state = 'clade_variant_trace'
        elif different and published_exception:
            state = 'published_clade_exception_observed'
        elif different:
            state = 'clade_variant_observed'
        else:
            state = 'clade_canonical_observed' if enough else 'clade_canonical_trace'
        if published_exception:
            relationship = 'published_species_or_isolate_exception'
        elif scoped_match:
            relationship = 'supported_genus_expectation' if assignment.get('reference_scope') == 'genus' else 'supported_scoped_reference'
        elif supported_scoped:
            relationship = 'differs_from_supported_scoped_reference'
        else:
            relationship = 'no_narrower_supported_reference'
        pattern = bool(candidate.get('known_detected') or nominated)
        review = permitted and different and (
            observation['measurement_disagreement'] or (enough and pattern))
        action = ('resolve_variant_measurement_disagreement' if observation['measurement_disagreement'] else
                  'validate_published_clade_exception_in_sample' if published_exception else
                  'validate_clade_variant_terminal_location' if pattern else
                  'review_variant_sequence_without_terminal_claim')
        item = {
            'canonical_family': family, 'motif': candidate.get('motif') or family,
            'reads': observation['current_reads'], 'original_discovery_reads': observation['original_discovery_reads'],
            'sequence_observation_status': observation['status'],
            'candidate_status': candidate.get('candidate_status'),
            'novel_evidence_level': candidate.get('novel_evidence_level'),
            'known_repeat_pattern_supported': bool(candidate.get('known_detected')),
            'existing_candidate_nomination': nominated,
            'sequence_only_observation': positive and not pattern,
            'global_catalogue_status': catalogue_status(family, global_catalogue),
            'sample_reference_status': candidate.get('reference_status'),
            'clade_relationship': 'differs_from_clade_default' if different else 'matches_clade_default',
            'scoped_reference_relationship': relationship,
            'published_scoped_exception': published_exception,
            'status': state, 'is_representative': family == representative_family,
            'comparison_limited_by_library': not permitted,
            'review_required': review,
            'review_reason': ('measurement_disagreement' if observation['measurement_disagreement'] else
                              'nondefault_clade_family_with_existing_pattern_or_nomination') if review else None,
            'next_action': action if review else None,
            'validation_needed': (['resolve_discordant_measurements'] if observation['measurement_disagreement'] else
                                  ['chromosome_terminal_localization_and_orientation',
                                   'independent_biological_sample_or_library_if_claiming_new_biology',
                                   'scoped_literature_review_before_any_novelty_claim']) if review else [],
            'biological_novelty_claimed': False, 'terminal_location_confirmed': False,
        }
        rows.append(item)
        if review:
            reviews.append({key: deepcopy(item[key]) for key in
                            ('canonical_family', 'motif', 'status', 'global_catalogue_status',
                             'scoped_reference_relationship', 'published_scoped_exception',
                             'review_reason', 'next_action', 'validation_needed')})
    rows.sort(key=lambda row: row['canonical_family'])
    reviews.sort(key=lambda row: row['canonical_family'])
    variants = [row for row in rows if row['status'] in {'clade_variant_observed', 'published_clade_exception_observed'}]
    canonicals = [row for row in rows if row['status'] == 'clade_canonical_observed']
    traces = [row for row in rows if row['status'] == 'clade_variant_trace']
    if not permitted:
        status = 'library_context_limited'
    elif variants and canonicals:
        status = 'clade_canonical_and_variant_observed'
    elif variants:
        status = 'clade_variant_observed'
    elif canonicals:
        status = 'clade_canonical_observed'
    elif traces:
        status = 'clade_variant_trace_only'
    elif any(row['status'] == 'measurement_disagreement' for row in rows):
        status = 'measurement_disagreement'
    else:
        status = 'no_relevant_family_observation'
    result.update(status=status, families=rows, review_items=reviews, review_required=bool(reviews),
                  variant_observed=bool(variants) if permitted else None,
                  canonical_observed=bool(canonicals) if permitted else None,
                  variant_observed_family_count=len(variants), canonical_observed_family_count=len(canonicals),
                  variant_trace_family_count=len(traces),
                  published_exception_family_count=sum(row['published_scoped_exception'] and row['status'] == 'published_clade_exception_observed' for row in rows))
    return result
