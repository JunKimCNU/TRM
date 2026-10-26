"""Transparent short-read TRM interpretation; no species inferred from motifs.

Thresholds below are screening rules, not calibrated probabilities or FDR.
Background-control uncertainty remains an annotation rather than a universal
exclusion rule. Unknown repeat families require direct junction evidence.
"""
from __future__ import annotations

import math
import re
from evidence_batch import canonical
from reference_evidence import SUPPORTED_REFERENCE_EVIDENCE

SUPPORTED_EVIDENCE = SUPPORTED_REFERENCE_EVIDENCE
POLICY = {
    'version': '0.7.2', 'minimum_known_reads': 10, 'minimum_full_repeat_reads': 5,
    'minimum_novel_junction_reads': 10, 'minimum_novel_each_end': 3,
    'minimum_consistent_fraction': 0.8, 'minimum_novel_wilson_lower': 0.6,
    'minimum_partial_consistent_fraction': 0.6,
    'maximum_expected_family_multiplier': 1.25,
    'minimum_distinct_junction_sequences': 10,
    'minimum_physical_flank_clusters': 2,
    'flank_cluster_hamming_distance': 2,
    'flank_cluster_max_offset': 16,
    'flank_cluster_min_overlap': 30,
    'minimum_deduplicated_consistency': 0.8,
    'clean_two_flank_review_trigger_records': 1,
    'novel_evidence_levels': ['uncatalogued_repeat', 'geometry_nomination', 'provisional_candidate', 'externally_confirmed'],
    'external_confirmation': 'Never assigned by this short-read pipeline; independent chromosome-terminal evidence is required.',
    'policy_design': 'v0.5 pattern thresholds and v0.6 sequence-observation, classification and representative-selection rules are unchanged. v0.7 adds an orthogonal comparison of the observed representative with applicable literature expectations; it never assigns the main TRM, forces an expected family, removes other candidates or upgrades short-read evidence to terminal confirmation.',
    'observation_rule': 'At least 10 current selected-source reads support a representative sequence observation, independently of full-repeat and boundary pattern support. Smaller positive observations remain traces. An applicable exact species/isolate published family may precede an unrelated novel nomination, without erasing alternative candidates or establishing chromosome-terminal location.',
    'novelty_scope': 'Absent from this version of the curated catalogue, not a claim of literature-wide novelty.',
    'calibration': 'Transparent screening thresholds; independent biological validation still required.',
    'count_source_rule': 'Use a complete current measured set of reads, full-repeat records, repeat bases and max tract together. Geometry-only evidence falls back to the complete original count set. Preserve original counts and flag positive-original/zero-measured disagreements; never choose per-field maxima.',
    'reference_agreement_rule': 'Compare exact primitive repeat families after rotation/reverse-complement normalization. Retain species/isolate, genus and weak clade scopes; published candidates and unknown/out-of-scope contexts impose no established-reference requirement. Missing reports, explicit measured zeros, traces and measurement disagreements remain distinct. A mismatch is a review annotation, not proof of contamination, biological novelty, telomere absence or a different main TRM.',
    'patch_scope': 'v0.7.2 requires an explicit evaluated period scope before conditional nondetection and covers every assigned expected primitive family within that model. Resolved eukaryotic metadata remains required for taxonomic attribution. Positive observations, pattern thresholds and novel-candidate promotion rules are unchanged.',
}


def _count(value, label):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f'{label} must be a nonnegative integer')
    return value


def flank_clusters(sequences, representatives=None):
    """Deterministic error/jitter-tolerant components; operational, not loci.

    Transitive links deliberately merge close anchors conservatively. With a
    saturated sketch, unseen bridging anchors could merge retained components;
    the retained component count is therefore NOT a lower bound on full-pool
    components. Different offsets at one locus can also remain separate.
    """
    seqs = sorted(set(sequences))
    if any(len(s) != 20 or set(s) - set('ACGT') for s in seqs):
        raise ValueError('Physical flank anchors must contain exactly 20 ACGT bases')
    reps = representatives or {s: s for s in seqs}
    for sequence in seqs:
        rep = reps.get(sequence, '')
        if not 20 <= len(rep) <= 80 or not rep.startswith(sequence) or set(rep) - set('ACGT'):
            raise ValueError('Representative flank must start with its normalized 20-bp anchor and contain 20..80 ACGT bases')
    # Two bits per validated A/C/G/T base. Folding each XOR lane with its
    # neighbor counts differing BASES, including substitutions differing in
    # both bits. This is the same Hamming/overlap rule as character comparison.
    encoding = str.maketrans('ACGT', '0123')
    lengths = [len(reps[s]) for s in seqs]
    packed = [int(reps[s].translate(encoding), 4) for s in seqs]
    anchors = [value >> (2 * (length - 20)) for value, length in zip(packed, lengths)]
    lane_bits = int('55' * 20, 16)  # one low bit per base, up to 80 bases
    anchor_mask = lane_bits & ((1 << 40) - 1)
    distance = POLICY['flank_cluster_hamming_distance']
    overlap = POLICY['flank_cluster_min_overlap']
    max_offset = POLICY['flank_cluster_max_offset']
    alignments = {}

    def close(i, j):
        difference = anchors[i] ^ anchors[j]
        if ((difference | (difference >> 1)) & anchor_mask).bit_count() <= distance:
            return True
        left_length, right_length = lengths[i], lengths[j]
        if min(left_length, right_length) < overlap:
            return False
        key = (left_length, right_length)
        offsets = alignments.get(key)
        if offsets is None:
            offsets = []
            low = max(-max_offset, overlap - right_length)
            high = min(max_offset, left_length - overlap)
            for offset in range(low, high + 1):
                a, b = max(offset, 0), max(-offset, 0)
                length = min(left_length - a, right_length - b)
                # Right shift puts the compared substring at the low end;
                # the lane mask removes bases outside that overlap.
                offsets.append((2 * (left_length - a - length),
                                2 * (right_length - b - length),
                                lane_bits & ((1 << (2 * length)) - 1)))
            alignments[key] = offsets
        left, right = packed[i], packed[j]
        for left_shift, right_shift, mask in offsets:
            difference = (left >> left_shift) ^ (right >> right_shift)
            if ((difference | (difference >> 1)) & mask).bit_count() <= distance:
                return True
        return False
    parents = list(range(len(seqs)))
    def find(i):
        while parents[i] != i:
            parents[i] = parents[parents[i]]
            i = parents[i]
        return i
    for i in range(len(seqs)):
        for j in range(i):
            if close(i, j):
                left, right = find(i), find(j)
                parents[max(left, right)] = min(left, right)
    groups = {}
    for i, sequence in enumerate(seqs):
        groups.setdefault(find(i), []).append(sequence)
    return [groups[key] for key in sorted(groups)]


def directional_flank_cluster_diagnostics(family, geometry, flanks=None, clusters=None):
    """Joint context/direction annotation; never a classification gate.

    The geometry axis stays fixed before clustering. A saturated retained
    component count is not a population lower bound: omitted anchors may bridge
    retained components. Missing new metadata must not turn into zero support.
    """
    metric_names = (
        'flank_events_observed', 'flank_events_retained', 'flank_events_not_retained',
        'retained_anchors', 'retained_cluster_count', 'retained_clusters_with_consistent_events',
        'retained_clusters_with_opposed_events', 'retained_clusters_with_both_directions',
        'retained_consistent_only_clusters', 'retained_opposed_only_clusters',
        'retained_clusters_with_ambiguous_events', 'retained_clusters_consistent_at_both_original_read_edges',
        'retained_consistent_events', 'retained_opposed_events', 'retained_ambiguous_events',
        'largest_retained_cluster_events', 'largest_retained_consistent_cluster_events',
        'largest_retained_cluster_event_fraction', 'largest_retained_consistent_cluster_event_fraction',
        'all_observed_clusters_exact', 'all_observed_clusters_with_consistent_events_exact',
        'all_observed_clusters_with_opposed_events_exact',
        'all_observed_largest_consistent_cluster_event_fraction')
    result = {'schema_version': 1, 'status': 'unmeasured', 'reason': 'complete_flank_validation_unavailable',
        'orientation_basis': geometry.get('orientation_basis'), 'display_motif': geometry.get('display_motif'),
        'orientation_interpretation': 'unmeasured', 'flank_sketch_saturated': None,
        'raw_read_edge_counts_available': None, 'cluster_scope': 'unmeasured', 'clusters_retained': None,
        'changes_classification': False,
        'interpretation': 'Operational context groups and event concentration, not independent molecules, loci or chromosome ends. Keep the original published/inferred axis fixed. Saturated retained groups are neither exact population counts nor lower bounds; retained proportions are not population estimates. Read signatures cannot be linked to context clusters from this stored schema.',
        **dict.fromkeys(metric_names)}
    if not isinstance(flanks, dict):
        return result
    result.update(status='incomplete', reason='flank_direction_metadata_incomplete_or_inconsistent')
    keys = ('left_forward', 'left_reverse', 'right_forward', 'right_reverse', 'ambiguous')
    def valid_count(value):
        return isinstance(value, int) and not isinstance(value, bool) and value >= 0
    def counts(q):
        return isinstance(q, dict) and all(valid_count(q.get(key)) for key in keys)
    if not isinstance(flanks.get('saturated'), bool):
        return result
    if not all(valid_count(flanks.get(key)) for key in ('observed_events', 'capacity', 'distinct_lower_bound')):
        return result
    entries = flanks.get('entries')
    if not isinstance(entries, list) or not flanks['capacity'] or len(entries) > flanks['capacity']:
        return result
    if any(not isinstance(entry, dict) or not valid_count(entry.get('events')) or not entry['events']
           or not counts(entry.get('normalized_quadrant_counts')) for entry in entries):
        return result
    if any(sum(entry['normalized_quadrant_counts'][key] for key in keys) != entry['events'] for entry in entries):
        return result
    retained = sum(entry['events'] for entry in entries)
    observed, bound, saturated = flanks['observed_events'], flanks['distinct_lower_bound'], flanks['saturated']
    if retained > observed or bound < len(entries) or bound > observed:
        return result
    if not saturated and (bound != len(entries) or observed != retained):
        return result
    if saturated and (len(entries) != flanks['capacity'] or bound <= len(entries) or observed-retained < bound-len(entries)):
        return result
    if flanks.get('flank_bases') != 20:
        return result
    try:
        keyed = {entry['sequence']: entry for entry in entries}
        if len(keyed) != len(entries):
            return result
        if clusters is None:
            clusters = flank_clusters(list(keyed), {key: entry['representative_flank'] for key, entry in keyed.items()})
    except (KeyError, TypeError, ValueError):
        return result
    display = geometry.get('display_motif')
    if reverse_complement(family) in rotations(family):
        result.update(status='strand_ambiguous', reason='palindromic_family_has_no_directional_axis')
        return result
    if display not in rotations(family) | rotations(reverse_complement(family)):
        result.update(reason='saved_geometry_orientation_unavailable')
        return result
    chosen_a = display in rotations(family)
    consistent_keys = ('left_reverse', 'right_forward') if chosen_a else ('left_forward', 'right_reverse')
    opposed_keys = ('left_forward', 'right_reverse') if chosen_a else ('left_reverse', 'right_forward')
    def direction(q):
        return sum(q[key] for key in consistent_keys), sum(q[key] for key in opposed_keys), q['ambiguous']
    raw_complete = all('quadrant_counts' in entry for entry in entries)
    for entry in entries:
        raw = entry.get('quadrant_counts')
        if raw is None:
            raw_complete = False
            continue
        if not counts(raw) or sum(raw[key] for key in keys) != entry['events']:
            return result
        normalized = entry['normalized_quadrant_counts']
        expected = dict(left_forward=raw['left_forward']+raw['right_reverse'],
                        left_reverse=raw['left_reverse']+raw['right_forward'],
                        right_forward=0, right_reverse=0, ambiguous=raw['ambiguous'])
        if normalized != expected:
            result.update(reason='raw_and_outward_normalized_flank_quadrants_disagree')
            return result
    detailed = []
    for group in clusters:
        members = [keyed[key] for key in group]
        q = {key: sum(entry['normalized_quadrant_counts'][key] for entry in members) for key in keys}
        consistent, opposed, ambiguous = direction(q)
        raw = {key: sum(entry['quadrant_counts'][key] for entry in members) for key in keys} if raw_complete else None
        detailed.append({'representative': group[0], 'retained_contexts': len(group),
            'events': sum(entry['events'] for entry in members), 'consistent_events': consistent,
            'opposed_events': opposed, 'ambiguous_events': ambiguous,
            'consistent_left_read_edge_events': raw[consistent_keys[0]] if raw is not None else None,
            'consistent_right_read_edge_events': raw[consistent_keys[1]] if raw is not None else None})
    consistent = sum(group['consistent_events'] for group in detailed)
    largest = max((group['events'] for group in detailed), default=0)
    largest_consistent = max((group['consistent_events'] for group in detailed), default=0)
    consistent_clusters = sum(group['consistent_events'] > 0 for group in detailed)
    opposed_clusters = sum(group['opposed_events'] > 0 for group in detailed)
    fraction = largest_consistent / consistent if consistent else None
    for group in detailed:
        group['fraction_of_retained_flank_events'] = group['events']/retained if retained else None
        group['fraction_of_retained_consistent_events'] = group['consistent_events']/consistent if consistent else None
    result.update(status='measured', reason=None, flank_sketch_saturated=saturated,
        raw_read_edge_counts_available=raw_complete,
        orientation_interpretation=('not_established_no_oriented_events' if geometry.get('oriented_junction_events') == 0 else
            'published_orientation_in_saved_catalogue' if geometry.get('orientation_basis') == 'published_motif_orientation' else
            'inferred_from_same_junction_events'),
        cluster_scope='retained_sketch_only_not_a_population_bound' if saturated else 'all_observed_usable_flank_anchors',
        flank_events_observed=observed, flank_events_retained=retained, flank_events_not_retained=observed-retained,
        retained_anchors=len(entries), retained_cluster_count=len(detailed),
        retained_clusters_with_consistent_events=consistent_clusters,
        retained_clusters_with_opposed_events=opposed_clusters,
        retained_clusters_with_both_directions=sum(group['consistent_events'] > 0 and group['opposed_events'] > 0 for group in detailed),
        retained_consistent_only_clusters=sum(group['consistent_events'] > 0 and group['opposed_events'] == group['ambiguous_events'] == 0 for group in detailed),
        retained_opposed_only_clusters=sum(group['opposed_events'] > 0 and group['consistent_events'] == group['ambiguous_events'] == 0 for group in detailed),
        retained_clusters_with_ambiguous_events=sum(group['ambiguous_events'] > 0 for group in detailed),
        retained_clusters_consistent_at_both_original_read_edges=(sum(group['consistent_left_read_edge_events'] > 0 and group['consistent_right_read_edge_events'] > 0 for group in detailed) if raw_complete else None),
        retained_consistent_events=consistent, retained_opposed_events=sum(group['opposed_events'] for group in detailed),
        retained_ambiguous_events=sum(group['ambiguous_events'] for group in detailed),
        largest_retained_cluster_events=largest, largest_retained_consistent_cluster_events=largest_consistent,
        largest_retained_cluster_event_fraction=largest/retained if retained else None,
        largest_retained_consistent_cluster_event_fraction=fraction,
        all_observed_clusters_exact=len(detailed) if not saturated else None,
        all_observed_clusters_with_consistent_events_exact=consistent_clusters if not saturated else None,
        all_observed_clusters_with_opposed_events_exact=opposed_clusters if not saturated else None,
        all_observed_largest_consistent_cluster_event_fraction=fraction if not saturated else None,
        clusters_retained=detailed)
    return result


def confirmation_readiness(family, geometry, junction, measured_complete):
    """Audit distinct sequence/context evidence without inventing independence."""
    result = {
        'status': 'unmeasured', 'provisional_ready': False,
        'validation_measured': False, 'checks': {}, 'missing_evidence': [],
        'distinct_junction_sequences_lower_bound': None,
        'distinct_junction_sequences_exact': None,
        'distinct_physical_flanks_lower_bound': None,
        'physical_flank_cluster_count': None, 'retained_flank_cluster_count': None,
        'clean_two_flank_repeat_records': None,
        'deduplicated_orientation': {'status': 'unmeasured', 'consistent_fraction': None},
        'directional_flank_clusters': directional_flank_cluster_diagnostics(family, geometry),
        'flank_cluster_rule': 'Connected components: outward-normalized 20-bp anchor Hamming<=2 OR representative-flank overlap>=30 bp at offsets -16..16 with <=2 mismatches. Operational error/jitter-tolerant groups, not loci. Omitted anchors may bridge groups.',
        'limitations': ['Distinct read sequences are not independent molecules; identical sequences and reverse complements are collapsed, but PCR copies with errors and overlapping fragments can remain.',
                       'Flank clustering is not a quality-aware error model; optional base-quality histograms are reported separately.',
                       'Distinct flank contexts/clusters are not chromosome counts or mapped subtelomeric loci.',
                       'Sequence hashes have theoretical collision risk; saturated sketches give bounds, not unbiased population ratios.',
                       'No short-read evidence level establishes chromosome-terminal location or external confirmation.'],
    }
    validation = junction.get('validation') if isinstance(junction, dict) else None
    required = ('schema_version', 'measured', 'clean_two_flank_repeat_records', 'junction_read_signatures', 'flank_contexts')
    if not measured_complete or not isinstance(validation, dict) or not validation.get('measured'):
        result['missing_evidence'] = ['complete_v050_distinct_sequence_and_flank_measurement']
        return result
    missing = [key for key in required if key not in validation]
    if missing or validation.get('schema_version') != 1:
        result.update(status='incomplete', missing_evidence=['validation_schema_or_fields_incomplete'])
        return result
    sketches = {}
    for name, observed_key, count_key in (('junction_read_signatures', 'observed_records', 'records'),
                                         ('flank_contexts', 'observed_events', 'events')):
        sketch = validation[name]
        needed = ('capacity', 'distinct_lower_bound', 'saturated', 'entries', observed_key)
        if not isinstance(sketch, dict) or any(key not in sketch for key in needed):
            result.update(status='incomplete', missing_evidence=[name+'_fields_incomplete'])
            return result
        cap = _count(sketch['capacity'], name+'.capacity')
        observed = _count(sketch[observed_key], name+'.'+observed_key)
        bound = _count(sketch['distinct_lower_bound'], name+'.distinct_lower_bound')
        entries = sketch['entries']
        if not cap or not isinstance(sketch['saturated'], bool) or not isinstance(entries, list) or len(entries) > cap:
            raise ValueError(f'Invalid {name} sketch shape')
        signatures, total = set(), 0
        for entry in entries:
            signature = entry.get('signature')
            if not isinstance(signature, str) or not signature or signature in signatures:
                raise ValueError(f'Invalid or repeated {name} signature')
            signatures.add(signature)
            n = _count(entry.get(count_key), name+'.'+count_key)
            if not n:
                raise ValueError(f'{name} retained entry has zero observations')
            total += n
            q = entry.get('normalized_quadrant_counts')
            if not isinstance(q, dict) or any(key not in q for key in ('left_forward', 'left_reverse', 'right_forward', 'right_reverse', 'ambiguous')):
                result.update(status='incomplete', missing_evidence=[name+'_normalized_orientation_incomplete'])
                return result
            for key, value in q.items():
                _count(value, name+'.normalized_quadrant_counts.'+key)
        if total > observed or bound < len(entries) or bound > observed:
            raise ValueError(f'{name} counts contradict retained entries')
        if not sketch['saturated'] and (bound != len(entries) or total != observed):
            raise ValueError(f'Unsaturated {name} must account for all observations and distinct entries')
        if sketch['saturated'] and (len(entries) != cap or bound <= len(entries) or observed - total < bound - len(entries)):
            raise ValueError(f'Saturated {name} lacks evidence of omitted entries')
        sketches[name] = (sketch, total)
    reads, retained_records = sketches['junction_read_signatures']
    flanks, _ = sketches['flank_contexts']
    if flanks.get('flank_bases') != 20:
        result.update(status='incomplete', missing_evidence=['physical_flank_length_not_20'])
        return result
    if reads['observed_records'] != junction.get('junction_reads'):
        raise ValueError('Read signature observations differ from junction_reads')
    internal = _count(validation['clean_two_flank_repeat_records'], 'clean_two_flank_repeat_records')
    chosen_a = geometry.get('display_motif') in rotations(family)
    consistent = opposite = ambiguous = 0
    for entry in reads['entries']:
        q = entry['normalized_quadrant_counts']
        a = bool(q['left_reverse'] or q['right_forward'])
        b = bool(q['left_forward'] or q['right_reverse'])
        if q['ambiguous'] or a == b:
            ambiguous += 1
        elif a == chosen_a:
            consistent += 1
        else:
            opposite += 1
    retained = len(reads['entries'])
    omitted_records = reads['observed_records'] - retained_records
    maximum_distinct = retained + omitted_records
    lower = consistent / maximum_distinct if maximum_distinct else None
    upper = (consistent + omitted_records) / maximum_distinct if maximum_distinct else None
    exact_fraction = consistent / retained if retained and not reads['saturated'] else None
    support = (lower is not None and lower >= POLICY['minimum_deduplicated_consistency']
               and geometry.get('status') != 'strand_ambiguous')
    dedup = {'status': 'supported_by_conservative_bound' if support and reads['saturated'] else (
                       'supported' if support else 'unresolved_saturation' if reads['saturated'] else 'not_supported'),
             'chosen_orientation_basis': 'Same direction as raw published/inferred geometry; not reselected after deduplication.',
             'consistent_unique_retained': consistent, 'opposite_unique_retained': opposite,
             'ambiguous_or_conflicting_unique_retained': ambiguous,
             'retained_unique_signatures': retained, 'saturated': reads['saturated'],
             'omitted_records_upper_bound_on_unseen_unique_sequences': omitted_records,
             'total_unique_sequences_upper_bound': maximum_distinct,
             'consistent_fraction': exact_fraction,
             'consistent_fraction_lower_bound': lower, 'consistent_fraction_upper_bound': upper,
             'interpretation': 'One vote per RC-collapsed read signature. Conflicting or ambiguous signatures count against consistency. Bounds assign all omitted records to distinct adverse/favorable sequences; no unbiased sketch ratio, independence, p-value or FDR claim.'}
    if any('representative_flank' not in entry for entry in flanks['entries']):
        result.update(status='incomplete', missing_evidence=['flank_representatives_unmeasured'])
        return result
    clusters = flank_clusters([entry.get('sequence', '') for entry in flanks['entries']],
                              {entry['sequence']: entry['representative_flank'] for entry in flanks['entries']})
    if len({entry['sequence'] for entry in flanks['entries']}) != len(flanks['entries']):
        raise ValueError('Duplicate physical flank sequence keys')
    flank_count = len(clusters) if not flanks['saturated'] else None
    checks = {
        'complete_validation_measurement': True,
        'minimum_distinct_junction_sequences': reads['distinct_lower_bound'] >= POLICY['minimum_distinct_junction_sequences'],
        'deduplicated_orientation_consistent': support if not reads['saturated'] or support else None,
        'minimum_error_tolerant_flank_clusters': flank_count >= POLICY['minimum_physical_flank_clusters'] if flank_count is not None else None,
        'no_clean_two_flank_counterevidence': internal < POLICY['clean_two_flank_review_trigger_records'],
    }
    result.update(status='measured', validation_measured=True, checks=checks,
        missing_evidence=[key for key, value in checks.items() if value is not True],
        provisional_ready=all(value is True for value in checks.values()),
        distinct_junction_sequences_lower_bound=reads['distinct_lower_bound'],
        distinct_junction_sequences_exact=reads['distinct_lower_bound'] if not reads['saturated'] else None,
        distinct_physical_flanks_lower_bound=flanks['distinct_lower_bound'],
        physical_flank_cluster_count=flank_count, retained_flank_cluster_count=len(clusters),
        flank_context_sketch_saturated=flanks['saturated'],
        directional_flank_clusters=directional_flank_cluster_diagnostics(family, geometry, flanks, clusters),
        clean_two_flank_repeat_records=internal, deduplicated_orientation=dedup,
        internal_counterevidence_interpretation='A tract with clean flanks on both sides occurs within a read. This requires review before stronger nomination; interstitial copies can coexist with genuine telomeres, so it does not disprove terminal copies.',
        flank_clusters_retained=[{'representative': group[0], 'retained_contexts': len(group)} for group in clusters])
    return result


def resolved_eukaryotic_scope(metadata):
    """Require positive metadata evidence before attributing eukaryotic TRMs.

    Lookup failure or a missing lineage cannot increase confidence. This gate
    limits interpretation only; it never filters observed repeat families or
    infers taxonomy from their sequence. An explicit Eukaryota TaxID can stand
    alone because that root taxon need not contain itself in its lineage.
    """
    if str(metadata.get('status', '')).strip() != 'resolved':
        return False
    taxid = str(metadata.get('tax_id', '')).strip()
    if not re.fullmatch(r'[1-9][0-9]*', taxid) or not str(metadata.get('scientific_name', '')).strip():
        return False
    lineage = {token.strip() for token in str(metadata.get('lineage', '')).split(';')}
    mixed = (str(metadata.get('metagenome', '')).strip().lower() in {'true', 'yes', '1'}
             or str(metadata.get('library_source', '')).strip().upper() == 'METAGENOMIC')
    return not mixed and not lineage.intersection({'Bacteria', 'Archaea'}) and (
        'Eukaryota' in lineage or taxid == '2759')


def library_context(metadata, reviewed=None, adequacy=None):
    """Reviewed identity/representation limits never silently reassign a taxon."""
    evidence, reasons = [], []
    reviewed = reviewed or {}
    code = reviewed.get('library_review_code')
    identity_conflict = code == 'taxonomic_identity_conflict'
    if code in {'RNA_component_supported', 'metadata_conflict_unresolved', 'HiC_supported',
                'HiC_annotation_only', 'taxonomic_identity_conflict'}:
        reasons.append(code)
        evidence.append({'type': 'reviewed_library_evidence', 'code': code,
                         'source': reviewed.get('detailed_evidence_path'), 'sha256': reviewed.get('evidence_sha256')})
    for field in ('library_strategy', 'library_source', 'library_selection', 'library_construction_protocol',
                  'library_construction', 'library_description'):
        text = str(metadata.get(field, ''))
        if re.search(r'\bRNA\b|RNA[- ]?seq|\bcDNA\b|poly[ -]?A|Hi[- ]?C|chromatin.*conformation', text, re.I):
            reasons.append('rna_or_hic_metadata_requires_review')
            evidence.append({'type': 'archive_metadata', 'field': field, 'value': text})
    blockers = (adequacy or {}).get('library', {}).get('uniform_wgs_model_blockers', [])
    if blockers:
        reasons.extend(blockers)
    return {'limited': bool(reasons), 'reasons': sorted(set(reasons)), 'evidence': evidence,
            'reviewed_code': code, 'RNA_component_supported': code == 'RNA_component_supported',
            'identity_conflict_supported': identity_conflict,
            'metadata_contradiction_proves_RNA': False,
            'interpretation': (
                'Reviewed sequence evidence supports a conflict between the registered taxon and observed reference consistency. '
                'Registered-taxonomy TRM interpretation requires review; no automatic taxon reassignment, RNA identity, '
                'sample purity or mixture proportion is established. Observed repeats remain available.'
                if identity_conflict else
                'Metadata contradictions or nonuniform sampling limit genomic novelty/nondetection interpretation; '
                'metadata alone does not establish RNA identity, purity, or mixture proportion.')}


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def rotations(sequence):
    return {sequence[i:] + sequence[:i] for i in range(len(sequence))}


def wilson_lower(success, total):
    if not total:
        return 0.0
    p = success / total
    z2 = 1.96 ** 2
    return max(0.0, (p + z2/(2*total) - 1.96*math.sqrt(p*(1-p)/total + z2/(4*total**2))) / (1+z2/total))


def catalog_index(rows):
    indexed = {}
    for row in rows:
        indexed.setdefault(canonical(row['motif']), []).append(row)
    return indexed


def directional_evidence(family, conventional, junction):
    """Canonical F/R quadrants transformed to a conventional G-strand motif.

    Expected: reverse complement at the read's left edge, conventional motif at
    its right edge. Unknown orientation is selected from data and labelled so.
    A palindromic family cannot give orientation evidence.
    """
    if junction is None:
        return {'status': 'not_measured', 'consistent_junctions': None,
                'opposite_junctions': None, 'left_consistent': None, 'right_consistent': None,
                'consistent_fraction': None, 'wilson_lower_descriptive': None,
                'orientation_basis': 'not_measured', 'display_motif': conventional or family,
                'unique_junction_reads': None, 'oriented_junction_events': None}
    q = {key: int(junction.get(key, 0)) for key in (
        'left_forward_junction_reads', 'left_reverse_junction_reads',
        'right_forward_junction_reads', 'right_reverse_junction_reads')}
    unique_junctions = junction.get('junction_reads')
    if unique_junctions is not None:
        unique_junctions = int(unique_junctions)
        if unique_junctions < 0 or any(count < 0 or count > unique_junctions for count in q.values()):
            raise ValueError('Invalid junction record/quadrant counts')
    elif any(count < 0 for count in q.values()):
        raise ValueError('Junction quadrant counts must be nonnegative')
    if reverse_complement(family) in rotations(family):
        return {'status': 'strand_ambiguous', 'consistent_junctions': 0,
                'opposite_junctions': 0, 'left_consistent': 0, 'right_consistent': 0,
                'consistent_fraction': None, 'wilson_lower_descriptive': None,
                'orientation_basis': 'palindromic_family', 'display_motif': conventional or family,
                'canonical_quadrants': q, 'unique_junction_reads': unique_junctions,
                'oriented_junction_events': 0}
    a = (q['left_reverse_junction_reads'], q['right_forward_junction_reads'])
    b = (q['left_forward_junction_reads'], q['right_reverse_junction_reads'])
    if conventional:
        forward = conventional in rotations(family)
        chosen, opposite = (a, b) if forward else (b, a)
        basis = 'published_motif_orientation'
        display = conventional
    else:
        forward = sum(a) >= sum(b)
        chosen, opposite = (a, b) if forward else (b, a)
        basis = 'inferred_from_junctions_not_a_published_orientation'
        display = family if forward else min(rotations(reverse_complement(family)))
    consistent, discordant = sum(chosen), sum(opposite)
    total = consistent + discordant
    fraction = consistent/total if total else None
    lower = wilson_lower(consistent, total) if total else None
    status = 'no_junction_support'
    if total:
        status = 'limited_junction_support'
    if total and unique_junctions is None:
        status = 'junction_record_count_unavailable'
    if total >= POLICY['minimum_novel_junction_reads'] and unique_junctions is not None and unique_junctions >= POLICY['minimum_novel_junction_reads']:
        opposite_lower = wilson_lower(discordant, total)
        if fraction >= POLICY['minimum_consistent_fraction']:
            status = 'one_end_only'
            if min(chosen) >= POLICY['minimum_novel_each_end'] and lower >= POLICY['minimum_novel_wilson_lower']:
                status = 'consistent_both_ends'
        elif (1 - fraction >= POLICY['minimum_consistent_fraction']
              and opposite_lower >= POLICY['minimum_novel_wilson_lower']):
            status = 'orientation_conflict'
        elif (fraction >= POLICY['minimum_partial_consistent_fraction'] and lower > 0.5
              and min(chosen) >= POLICY['minimum_novel_each_end']):
            status = 'partial_directional_support'
        else:
            status = 'mixed_directionality'
    return {'status': status, 'consistent_junctions': consistent,
            'opposite_junctions': discordant, 'left_consistent': chosen[0],
            'right_consistent': chosen[1], 'consistent_fraction': fraction,
            'wilson_lower_descriptive': lower, 'orientation_basis': basis,
            'display_motif': display, 'canonical_quadrants': q,
            'unique_junction_reads': unique_junctions, 'oriented_junction_events': total,
            'interpretation': 'Read-boundary geometry; not physical chromosome-end localization. Quadrant events can overlap on one record; reads and mates can be dependent. Wilson bounds are descriptive screening quantities, not p-values or FDR.'}


def sequence_observation(candidate):
    """Describe sequence counts without requiring a telomere-pattern call.

    Current and historical counts remain separate. In particular, an original
    positive/current-zero disagreement must not acquire a positive current
    count or be reduced to an unqualified negative result.
    """
    reads = _count(candidate.get('reads', 0), 'reads')
    original = _count(candidate.get('original_discovery_reads', reads),
                      'original_discovery_reads')
    disagreement = bool(candidate.get('measurement_disagreement'))
    minimum = POLICY['minimum_known_reads']
    eligible = reads >= minimum and not disagreement
    if disagreement:
        status = 'measurement_disagreement'
        interpretation = ('Positive original evidence and zero current measured reads conflict. '
                          'Retain both count sets; resolve the discrepancy before a positive pattern or nondetection interpretation.')
    elif eligible:
        status = 'sequence_observed'
        interpretation = ('The repeat family is observed at the representative-read threshold. '
                          'This is sequence evidence, not chromosome-terminal localization or a supported telomere pattern.')
    elif reads:
        status = 'trace_sequence_observed'
        interpretation = ('Positive repeat-family evidence is below the representative-read threshold. '
                          'Retain the trace; it does not by itself replace a stronger candidate or establish a telomeric repeat.')
    else:
        status = 'no_sequence_observed'
        interpretation = ('No positive sequence observation is reported in the selected evidence; an unreported family is not an explicit measured zero. '
                          'Scanner sensitivity, input completeness and library representation still limit nondetection interpretation.')
    return {'status': status, 'observed': bool(reads or original),
            'current_reads': reads, 'original_discovery_reads': original,
            'count_source': candidate.get('evidence_count_source', 'original_discovery'),
            'minimum_representative_reads': minimum, 'representative_eligible': eligible,
            'measurement_disagreement': disagreement, 'interpretation': interpretation}


PAIR_COUNT_FIELDS = ('supporting_pairs', 'both_mates_repeat_pairs', 'r1_only_repeat_pairs',
                     'r2_only_repeat_pairs', 'opposite_strand_pairs', 'same_strand_pairs',
                     'ambiguous_strand_pairs', 'both_full_repeat_pairs',
                     'exactly_one_full_repeat_pairs')


def paired_annotation(counts, context, count_source, reads, full_reads):
    """Describe only the paired measurement attached to the selected count set.

    An old paired discovery cannot fill missing paired evidence in a newer
    complete unpaired measurement. Native top-level provenance must declare an
    FR measurement with explicit full-input or synchronized-prefix scope;
    row-shaped legacy metadata alone is insufficient.
    This annotation is never consumed by ranking, classification or adequacy.
    """
    result = {'status': 'not_measured', 'measured': False, 'count_source': count_source,
              'layout': None, 'scope': None, **dict.fromkeys(PAIR_COUNT_FIELDS), 'one_mate_only_pairs': None,
              'interpretation': 'Pair orientation is descriptive consistency evidence only. '
                  'One-mate support is neutral; opposite orientation also occurs in internal repeats. '
                  'Counts do not establish independent molecules, terminal location, novelty or data adequacy.'}
    if (not isinstance(context, dict) or context.get('measured') is not True
            or context.get('layout') != 'fr'
            or context.get('scope') not in ('complete_inputs', 'synchronized_prefix_only')):
        return result
    native = counts.get('paired_evidence')
    if native is None:
        return result
    if not isinstance(native, dict) or any(key not in native for key in PAIR_COUNT_FIELDS):
        raise ValueError('Paired evidence count fields are incomplete')
    values = {key: _count(native[key], 'paired_evidence.'+key) for key in PAIR_COUNT_FIELDS}
    total_pairs = _count(context.get('validated_pairs'), 'paired_end.validated_pairs')
    supporting = values['supporting_pairs']
    both = values['both_mates_repeat_pairs']
    one_mate = values['r1_only_repeat_pairs'] + values['r2_only_repeat_pairs']
    if (supporting != both + one_mate or both != sum(values[key] for key in
            ('opposite_strand_pairs', 'same_strand_pairs', 'ambiguous_strand_pairs'))
            or supporting > total_pairs or reads != supporting + both
            or values['both_full_repeat_pairs'] > both
            or values['both_full_repeat_pairs'] + values['exactly_one_full_repeat_pairs'] > supporting
            or full_reads != 2 * values['both_full_repeat_pairs'] + values['exactly_one_full_repeat_pairs']):
        raise ValueError('Paired evidence counts conflict with the selected read count set')
    result.update(values, status=('measured_prefix' if context['scope'] == 'synchronized_prefix_only' else 'measured'),
                  measured=True, layout='fr', scope=context['scope'],
                  one_mate_only_pairs=one_mate)
    return result


def assess_candidate(raw, targeted, references, expected_families, mean_length, *,
                     raw_paired_end=None, measured_paired_end=None):
    family = canonical(raw['motif'])
    if targeted and 'motif' in targeted and canonical(targeted['motif']) != family:
        raise ValueError('Targeted evidence family does not match original candidate family')
    refs = references.get(family, [])
    established = [r for r in refs if r.get('evidence') in SUPPORTED_EVIDENCE
                   and r.get('reference_status') != 'candidate']
    # Every displayed reference is an exact family match, never a nearby anchor.
    conventional = (established or refs or [{}])[0].get('motif')
    geometry = directional_evidence(family, conventional,
                                   targeted.get('junction_evidence') if targeted else None)
    legacy_reads = int(raw.get('reads', 0))
    measured_complete = (targeted is not None
                         and all(key in targeted for key in ('reads', 'repeat_bases', 'max_tract'))
                         and isinstance(targeted.get('junction_evidence'), dict)
                         and 'full_repeat_records' in targeted['junction_evidence'])
    counts = targeted if measured_complete else raw
    full_value = (counts['junction_evidence']['full_repeat_records'] if measured_complete
                  else counts.get('full_repeat_reads', 0))
    if measured_complete:
        values = [counts[key] for key in ('reads', 'repeat_bases', 'max_tract')] + [full_value]
        if any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in values):
            raise ValueError('Measured evidence counts must be nonnegative integers')
        if full_value > counts['reads']:
            raise ValueError('Measured full-repeat records exceed measured reads')
    full = int(full_value)
    reads = int(counts.get('reads', 0))
    span = int(counts.get('max_tract', 0))
    repeated_bases = int(counts.get('repeat_bases', 0))
    targeted_only = measured_complete and legacy_reads == 0 and reads > 0
    measurement_disagreement = measured_complete and legacy_reads > 0 and reads == 0
    count_source = 'targeted_only' if targeted_only else ('measured_evidence' if measured_complete else 'original_discovery')
    pairs = paired_annotation(counts, measured_paired_end if measured_complete else raw_paired_end,
                              count_source, reads, full)
    consistent = geometry['consistent_junctions']
    opposite = geometry['opposite_junctions']
    total_junction = geometry['unique_junction_reads']
    full_and_junction = full >= POLICY['minimum_full_repeat_reads'] and geometry['status'] == 'consistent_both_ends'
    known_signal = bool(established) and reads >= POLICY['minimum_known_reads'] and full >= POLICY['minimum_full_repeat_reads']
    # Strong opposite orientation conflicts with a conventional telomere model.
    if geometry['status'] == 'orientation_conflict':
        known_signal = False
    geometry_nomination = (not refs and full_and_junction and total_junction is not None
                    and total_junction >= POLICY['minimum_novel_junction_reads']
                    and 3 <= len(family) <= 32)
    readiness = confirmation_readiness(family, geometry,
        targeted.get('junction_evidence') if targeted else None, measured_complete)
    novel_signal = geometry_nomination and readiness['provisional_ready'] and not measurement_disagreement
    if established:
        status = 'known_trm_pattern_supported' if known_signal and full_and_junction else (
            'known_repeat_detected_boundary_unresolved' if known_signal else 'weak_or_conflicting_known_repeat')
    elif refs:
        status = 'published_candidate_pattern_supported' if full_and_junction else 'published_candidate_unresolved'
    else:
        status = ('novel_provisional_candidate' if novel_signal else
                  'novel_geometry_nomination' if geometry_nomination else 'uncatalogued_repeat')
    notes = []
    if span < 500:
        notes.append('The selected measurement has no observed contiguous >=500-bp tract. Array length remains unresolved; 500 bp is not a universal biological requirement or a candidate promotion gate.')
    if full and total_junction and full < total_junction:
        notes.append('Boundary records outnumber full-repeat records; differs from the simple long-array random-fragment expectation, but library/tract lengths can alter this ratio.')
    if geometry['status'] == 'not_measured':
        notes.append('New junction evidence not measured for this family; absence of support is not a negative test.')
    if measured_complete and (reads, full, repeated_bases, span) != (
            legacy_reads, int(raw.get('full_repeat_reads', 0)), int(raw.get('repeat_bases', 0)), int(raw.get('max_tract', 0))):
        notes.append('Current measured and original de novo counts differ; retain both measurements and inspect detector sensitivity.')
    if targeted is not None and not measured_complete:
        notes.append('A complete measured count set was not supplied; all displayed counts retain original discovery values while available geometry is annotated.')
    if targeted_only:
        notes.append('This family was not reported by original discovery; displayed counts come from positive targeted evidence and must not be described as nondetection.')
    if measurement_disagreement:
        notes.append('Zero current measured reads conflict with positive original discovery. The original positive observation is retained; resolve this disagreement before interpreting nondetection. Neither result proves motif absence or ALT.')
    if raw.get('support', {}).get('status') == 'insufficient_background_resolution':
        notes.append('Single-control depth is insufficient for the old enrichment gate; this is not motif absence.')
    if geometry_nomination and not novel_signal:
        notes.append('Legacy read-boundary nomination is retained. Stronger provisional evidence is incomplete or has counterevidence; this is not a negative test of terminal copies.')
    if readiness.get('clean_two_flank_repeat_records'):
        notes.append('Clean flanks on both sides of a tract require internal-array review. These are within-read observations; internal chromosomal location is not established and terminal copies can coexist.')
    # The abundance score keeps its bounded expectation contribution. Reporting
    # an observed applicable reference separately is not a change to this score.
    score = repeated_bases * (POLICY['maximum_expected_family_multiplier'] if family in expected_families else 1.0)
    if geometry['status'] == 'consistent_both_ends':
        score *= 2.0
    elif geometry['status'] == 'orientation_conflict':
        score *= 0.1
    result = {
        'motif': geometry['display_motif'], 'canonical_family': family,
        'reference_status': 'established' if established else ('published_candidate' if refs else 'not_in_catalogue'),
        'published_motifs': sorted({r['motif'] for r in refs}),
        'reference_sources': sorted({r['source'] for r in refs}),
        'reference_evidence': sorted({r['evidence'] for r in refs}),
        'expected_family_match': family in expected_families if expected_families else None,
        'candidate_status': status, 'reads': reads, 'repeat_bases': repeated_bases,
        'full_repeat_reads': full, 'max_tract_bp': span, 'directionality': geometry,
        'unique_junction_reads': total_junction,
        'full_to_junction_ratio': full/total_junction if total_junction else None,
        'targeted_reads': targeted.get('reads') if targeted else None,
        'targeted_full_repeat_reads': targeted.get('junction_evidence', {}).get('full_repeat_records') if targeted else None,
        'targeted_repeat_bases': targeted.get('repeat_bases') if targeted else None,
        'targeted_max_tract_bp': targeted.get('max_tract') if targeted else None,
        'original_discovery_reads': legacy_reads,
        'original_discovery_full_repeat_reads': int(raw.get('full_repeat_reads', 0)),
        'original_discovery_repeat_bases': int(raw.get('repeat_bases', 0)),
        'original_discovery_max_tract_bp': int(raw.get('max_tract', 0)),
        'evidence_count_source': count_source,
        'paired_evidence': pairs,
        'measured_count_set_complete': measured_complete,
        'measurement_disagreement': measurement_disagreement,
        'measurement_disagreement_reason': 'positive_original_zero_measured' if measurement_disagreement else None,
        'background_status': raw.get('support', {}).get('status', 'not_assessed'),
        'ranking_score': score, 'novel_candidate': novel_signal,
        'geometry_nomination': geometry_nomination,
        'legacy_novel_candidate': geometry_nomination,
        'novel_evidence_level': ('known_reference_observation' if established else 'published_candidate' if refs else
                                 'provisional_candidate' if novel_signal else 'geometry_nomination' if geometry_nomination else 'uncatalogued_repeat'),
        'confirmation_readiness': readiness,
        'externally_confirmed': False,
        'external_confirmation': {
            'status': 'not_established_by_short_reads',
            'missing_evidence': ['chromosome_terminal_localization_and_orientation',
                                 'independent_biological_sample_or_library_validation',
                                 'literature_review_beyond_this_catalogue'],
            'scope': 'Required before confirming a catalogue-external biological TRM; internal screening readiness is not external confirmation.'},
        'long_array_full_to_junction_warning': bool(full and total_junction and full < total_junction),
        'known_detected': known_signal, 'notes': notes,
    }
    observation = sequence_observation(result)
    result.update(sequence_observation=observation,
                  known_observed=bool(established) and observation['representative_eligible'],
                  reference_observed=bool(refs) and observation['representative_eligible'])
    return result


def _ranked(candidates):
    return sorted(candidates, key=lambda candidate: (-candidate['ranking_score'],
                                                     candidate['canonical_family']))


def _expected_reference_observations(candidates, assignment, include_trace=False):
    """Require an applicable exact taxon reference and actual sequence counts."""
    if (assignment or {}).get('reference_scope') not in {'species', 'isolate'}:
        return []
    observed = []
    for candidate in candidates:
        if (candidate.get('expected_family_match') is not True
                or candidate.get('reference_status') not in {'established', 'published_candidate'}):
            continue
        observation = sequence_observation(candidate)
        if observation['representative_eligible'] or (
                include_trace and observation['status'] == 'trace_sequence_observed'):
            observed.append(candidate)
    return _ranked(observed)


def _alternative_architecture(assignment):
    assignment = assignment or {}
    return (assignment.get('reference_scope') in {'species', 'isolate'}
            and assignment.get('reference_status') == 'out_of_scope'
            and assignment.get('repeat_architecture') in {'irregular', 'retrotransposon', 'long_complex_family'})


def reference_agreement(candidates, representative, assignment, expected_observations, *,
                        analysis_complete=True, biological_scope=True, library=None,
                        qc_status=None):
    """Annotate reference agreement without changing candidate selection.

    ``expected_observations`` is the report's explicit measurement table. A
    missing family is never converted to a measured zero. Its selected counts
    remain separate from discovery counts, including when an independent
    targeted measurement disagrees. The only read threshold used here is the
    existing representative-observation threshold; it is not a new filter.
    """
    assignment = assignment or {}
    representative = representative or {}
    scope = assignment.get('reference_scope') or None
    reference_status = assignment.get('reference_status') or 'not_curated'
    strength = ({'species': 'species_reference', 'isolate': 'isolate_reference',
                 'genus': 'genus_expectation', 'clade': 'weak_clade_expectation'}.get(scope,
                 'unspecified_scope') if reference_status == 'supported' else
                'published_candidate_nonbinding' if reference_status == 'candidate' else
                'out_of_scope' if reference_status == 'out_of_scope' else 'unknown_reference')

    def family_of(row):
        motif = row.get('canonical_family') or row.get('motif')
        if not isinstance(motif, str) or not motif or set(motif.upper()) - set('ACGT'):
            raise ValueError('Reference agreement requires a nonempty ACGT repeat family')
        return canonical(motif)

    by_family = {family_of(candidate): candidate for candidate in candidates}
    related = [m.strip() for m in assignment.get('related_motifs', '').split(',') if m.strip()]
    expected = {}
    for motif in related:
        expected.setdefault(family_of({'motif': motif}), {'motif': motif, 'observations': []})
    binding_families = (set(expected) if reference_status == 'supported'
                        and scope in ('species', 'isolate', 'genus', 'clade') else set())
    for observation in expected_observations or []:
        family = family_of(observation)
        expected.setdefault(family, {'motif': observation.get('motif', family), 'observations': []})[
            'observations'].append(observation)

    families = []
    for family, group in sorted(expected.items()):
        candidate = by_family.get(family)
        candidate_reads = _count(candidate.get('reads', 0), 'candidate.reads') if candidate else None
        original_reads = (_count(candidate.get('original_discovery_reads', candidate_reads),
                                 'candidate.original_discovery_reads') if candidate else None)
        disagreement = bool(candidate and candidate.get('measurement_disagreement'))
        rows, measured_counts, sources = group['observations'], set(), set()
        explicit_zero = False
        for row in rows:
            recovery = row.get('recovery', 'not_reported')
            reads = row.get('reads', row.get('supporting_reads'))
            if reads is not None:
                measured_counts.add(_count(reads, 'expected_observation.reads'))
            if recovery == 'explicit_zero':
                if reads != 0:
                    raise ValueError('An explicit zero requires a measured read count of zero')
                explicit_zero = True
            if recovery == 'measurement_disagreement':
                disagreement = True
            if row.get('source'):
                sources.add(str(row['source']))
        # Alias rows must not choose a favorable count by order or maximum.
        disagreement |= len(measured_counts) > 1
        reads = next(iter(measured_counts)) if len(measured_counts) == 1 else None
        if not rows and candidate_reads:
            reads = candidate_reads
            sources.add('selected_candidate_count_set')
        # A positive saved discovery and a targeted zero are contradictory,
        # even if the caller has not yet rebuilt its candidate assessment.
        disagreement |= bool(explicit_zero and (candidate_reads or original_reads))
        if disagreement:
            recovery = 'measurement_disagreement'
        elif reads is not None and reads >= POLICY['minimum_known_reads']:
            recovery = 'sequence_observed'
        elif reads:
            recovery = 'trace_sequence_observed'
        elif explicit_zero:
            recovery = 'explicit_zero'
        else:
            recovery = 'not_reported'
        families.append({
            'motif': group['motif'], 'canonical_family': family,
            'recovery': recovery, 'reads': reads,
            'candidate_reads': candidate_reads,
            'original_discovery_reads': original_reads,
            'source': ';'.join(sorted(sources)) or 'family_not_in_saved_report',
            'reference_role': 'established_expectation' if family in binding_families else 'nonbinding_comparison',
            'measurement_disagreement': disagreement,
            'observations': [dict(row) for row in rows],
        })
    applicable = [row for row in families if row['canonical_family'] in binding_families]
    representative_family = family_of(representative) if representative else None
    representative_reads = _count(representative.get('reads', 0), 'representative.reads') if representative else 0
    representative_observed = bool(representative_reads and not representative.get('measurement_disagreement'))
    differs = bool(representative_observed and applicable and representative_family not in binding_families)
    positive = [row for row in applicable if row['recovery'] in ('sequence_observed', 'trace_sequence_observed')]
    observed = [row for row in applicable if row['recovery'] == 'sequence_observed']
    all_zero = bool(applicable) and all(row['recovery'] == 'explicit_zero' for row in applicable)
    mismatch = bool(differs and positive and scope in ('species', 'isolate', 'genus'))
    measurement_disagreement = (bool(representative.get('measurement_disagreement'))
                               or any(row['measurement_disagreement'] for row in applicable))
    requires_review = bool(scope in ('species', 'isolate', 'genus') and
                           (mismatch or all_zero or measurement_disagreement))
    limitations = ['Sequence agreement or abundance does not establish the sample main TRM or physical chromosome-terminal location.',
                   'Explicit zero means zero in the specified measured input and detector; it does not prove biological absence or sufficient detection power.',
                   'This comparison does not alter candidate counts, ranking, classification, selection, or data adequacy; all candidates remain available.']
    if scope == 'genus':
        limitations.append('A genus expectation is broader than a species/isolate result; species exceptions and sample provenance require review.')
    elif scope == 'clade':
        limitations.append('A clade default is weak and allows phylogenetic exceptions; disagreement is not an established species-reference conflict.')
    reasons = []
    if differs:
        reasons.append('representative_family_differs_from_applicable_expectation')
    if mismatch:
        reasons.append('both_expected_and_different_representative_families_are_observed')
    if all_zero:
        reasons.append('all_applicable_expected_families_have_explicit_measured_zero')
    if measurement_disagreement:
        reasons.append('expected_or_representative_measurements_disagree')
    if any(row['recovery'] == 'not_reported' for row in applicable):
        reasons.append('one_or_more_expected_families_unreported_not_measured_zero')
    if positive and not observed:
        reasons.append('expected_evidence_is_below_existing_representative_read_threshold')

    limited = True
    if not analysis_complete:
        status, action = 'analysis_incomplete', 'repair_or_complete_analysis'
    elif qc_status == 'observed_qc_failure':
        status, action = 'qc_limited', 'resolve_qc_before_reference_interpretation'
    elif not biological_scope:
        status, action = 'taxonomic_scope_unresolved', 'check_host_or_mixed_sample'
    elif (library or {}).get('limited'):
        status, action = 'library_context_limited', 'resolve_library_context'
    elif reference_status == 'out_of_scope':
        status, action = 'reference_out_of_scope', 'review_known_architecture'
    elif not applicable or scope not in ('species', 'isolate', 'genus', 'clade'):
        status, action = 'no_supported_reference', 'review_literature_and_terminal_evidence'
    else:
        limited = False
        if measurement_disagreement:
            status, action = 'measurement_disagreement', 'resolve_measurement_disagreement'
        elif mismatch and observed:
            status, action = 'observed_family_conflict', 'review_reference_and_representative'
        elif all_zero:
            status, action = 'expected_reference_not_observed', 'review_expected_family_nondetection'
        elif differs and scope == 'clade':
            status, action = 'broad_expectation_difference', 'review_clade_exception_without_forcing_default'
        elif observed and representative_observed and representative_family in binding_families:
            status, action = 'representative_reference_agreement', 'retain_observation_terminal_confirmation_separate'
        elif observed:
            status, action = 'reference_observed_no_representative', 'review_representative_observation'
        elif positive:
            status, action = 'expected_reference_trace', 'review_trace_observation'
        else:
            status, action = 'expected_reference_unmeasured', 'measure_expected_families_before_interpreting_nondetection'
    if limited:
        reasons.insert(0, status)
    if reference_status == 'candidate':
        reasons.append('published_candidate_is_not_an_established_reference_requirement')
    return {
        'schema_version': 1, 'status': status, 'reasons': reasons,
        'representative_motif': representative.get('motif'), 'representative_family': representative_family,
        'reference_scope': scope, 'reference_status': reference_status,
        'expectation_strength': strength, 'expected_families': families,
        'representative_differs_from_expectation': differs,
        'observed_family_conflict': mismatch,
        'expected_reference_not_observed': all_zero,
        'requires_conflict_review': requires_review,
        'interpretation_limited': limited, 'limitations': limitations,
        'main_trm_status': 'not_confirmed', 'main_trm_motif': None, 'next_action': action,
    }


def _observation_classification(candidate):
    if candidate['candidate_status'] == 'published_candidate_pattern_supported':
        return ('published_candidate_supported', '문헌 후보와 일치: 독립 검증 필요', 'validate_independent_biosample')
    if candidate['reference_status'] == 'established':
        return ('known_repeat_observed_pattern_unresolved', '알려진 반복 서열 관측·TRM 패턴 미확인',
                'review_observed_repeat')
    return ('published_candidate_observed_pattern_unresolved', '문헌 후보 서열 관측·TRM 패턴 미확인',
            'validate_published_candidate')


def select_representative(candidates, classification, assignment=None):
    """Choose the displayed sequence while preserving every alternative.

    ``classification`` is a status string or the legacy classify() three-tuple.
    A representative under a library/QC/architecture-limited status is only an
    observed sequence; it does not override the limiting sample conclusion.
    """
    code = classification[0] if isinstance(classification, (tuple, list)) else classification
    ordered = _ranked(candidates)
    expected = _expected_reference_observations(ordered, assignment)
    pools = {
        'known_trm_detected': [c for c in ordered if c['known_detected']],
        'known_repeat_detected_boundary_unresolved': [c for c in ordered if c['known_detected']],
        'novel_provisional_candidate': [c for c in ordered if c['novel_candidate']],
        'novel_geometry_nomination': [c for c in ordered if c.get('geometry_nomination')],
        'published_candidate_supported': [c for c in ordered if c['candidate_status'] == 'published_candidate_pattern_supported'],
        'known_repeat_observed_pattern_unresolved': [c for c in ordered if c['reference_status'] == 'established' and sequence_observation(c)['representative_eligible']],
        'published_candidate_observed_pattern_unresolved': [c for c in ordered if c['reference_status'] == 'published_candidate' and sequence_observation(c)['representative_eligible']],
        'published_repeat_trace_observed': _expected_reference_observations(ordered, assignment, include_trace=True),
        'detection_inconclusive_measurement': [c for c in ordered if c.get('measurement_disagreement')],
    }
    if code in pools:
        pool = pools[code]
        if code in {'known_repeat_observed_pattern_unresolved',
                    'published_candidate_observed_pattern_unresolved', 'published_candidate_supported'}:
            preferred = [c for c in expected if c in pool]
            return (preferred or pool or [None])[0]
        return (pool or [None])[0]
    if code in {'library_context_limited', 'detection_inconclusive_technical',
                'taxonomic_scope_unresolved', 'known_alternative_architecture'}:
        # Limited contexts retain the strongest reference observation first;
        # its identity is displayed separately from the sample-level decision.
        for pool in (expected, [c for c in ordered if c['known_detected']],
                     [c for c in ordered if c['novel_candidate']],
                     [c for c in ordered if c.get('geometry_nomination')],
                     [c for c in ordered if c['reference_status'] != 'not_in_catalogue'
                      and sequence_observation(c)['representative_eligible']],
                     [c for c in ordered if sequence_observation(c)['observed']]):
            if pool:
                return pool[0]
    return None


def expected_family_power_scope(adequacy, assignment=None):
    """Bind a conditional nondetection statement to the periods actually modelled.

    This is not a negative measurement of a target. Registered targeted scans
    can observe motifs outside the original de novo range, but they do not
    extend the original scanner's conditional power model by themselves.
    """
    assignment = assignment or {}
    value = assignment.get('related_motifs') or ''
    motifs = [item.strip() for item in value.split(',') if item.strip()] if isinstance(value, str) else list(value)
    grouped, invalid = {}, []
    for motif in motifs:
        if not isinstance(motif, str) or not motif or set(motif.upper()) - set('ACGT'):
            invalid.append(motif)
            continue
        grouped.setdefault(canonical(motif), set()).add(motif)
    power = adequacy.get('power') or {}
    period_rows = power.get('periods')
    available = (isinstance(period_rows, list) and bool(period_rows)
                 and all(isinstance(row, dict) and isinstance(row.get('period'), int)
                         and not isinstance(row['period'], bool) and 1 <= row['period'] <= 32
                         for row in period_rows))
    periods = sorted({row['period'] for row in period_rows}) if available else []
    scope = adequacy.get('detector_scope') or {}
    minimum, maximum = scope.get('min_k'), scope.get('max_k')
    scan_known = (isinstance(minimum, int) and not isinstance(minimum, bool)
                  and isinstance(maximum, int) and not isinstance(maximum, bool)
                  and 1 <= minimum <= maximum)
    families = []
    for family, spellings in sorted(grouped.items()):
        period = len(family)
        in_scan = minimum <= period <= maximum if scan_known else None
        covered = period in periods if available else None
        families.append({'canonical_family': family, 'expected_motifs': sorted(spellings),
                         'primitive_period': period, 'within_original_scan_periods': in_scan,
                         'within_evaluated_power_periods': covered,
                         'status': ('power_scope_unavailable' if not available else
                                    'outside_original_scan_periods' if in_scan is False else
                                    'outside_evaluated_power_periods' if not covered else
                                    'covered_by_evaluated_period_model')})
    if invalid:
        status = 'invalid_expected_units'
    elif not available:
        status = 'power_scope_unavailable'
    elif any(row['within_evaluated_power_periods'] is not True
             or row['within_original_scan_periods'] is False for row in families):
        status = 'expected_family_outside_power_scope'
    elif families:
        status = 'expected_families_covered'
    else:
        status = 'no_expected_family_assigned'
    return {'schema_version': 1, 'status': status, 'evaluated_periods': periods,
            'original_scan_period_range': {'min_k': minimum, 'max_k': maximum} if scan_known else None,
            'expected_families': families, 'invalid_expected_units': invalid,
            'expected_families_modelled': status == 'expected_families_covered' if families or invalid else None,
            'conditional_nondetection_scope_supported': status in {'expected_families_covered', 'no_expected_family_assigned'},
            'scope': 'Only evaluated raw-discovery period models; not evidence of motif absence or a validation of unmodelled targeted scans.',
            'targeted_measurement_extends_power_scope': False}


def classify(candidates, adequacy, analysis_complete=True, biological_scope=True,
             reference_assignment=None, library=None):
    if not analysis_complete:
        return ('analysis_failed', '분석 실패: 입력/처리 오류', 'repair_or_replace_input')
    if adequacy['status'] == 'observed_qc_failure':
        return ('detection_inconclusive_technical', '판정 보류: 데이터/QC 불일치', 'replace_or_check_library')
    if not biological_scope:
        return ('taxonomic_scope_unresolved', '진핵생물 TRM 판정 범위 확인 필요', 'check_host_or_mixed_sample')
    assignment = reference_assignment or {}
    if _alternative_architecture(assignment):
        return ('known_alternative_architecture', '문헌상 비정형 말단 구조: 관측 반복은 별도 후보로 보존', 'review_known_architecture')
    known = [c for c in candidates if c['known_detected']]
    if known:
        best_known = min(known, key=lambda candidate: (-candidate['ranking_score'], candidate['canonical_family']))
        if best_known['candidate_status'] == 'known_trm_pattern_supported':
            return ('known_trm_detected', '알려진 TRM 패턴 검출', 'canonical_screen_complete')
        return ('known_repeat_detected_boundary_unresolved', '알려진 반복 검출·말단 경계 미확인', 'review_boundary_support')
    if (library or {}).get('limited'):
        return ('library_context_limited', '라이브러리 제약: 관측 반복의 genomic TRM 해석 보류', 'resolve_library_context')
    expected = _expected_reference_observations(candidates, assignment)
    if expected:
        return _observation_classification(expected[0])
    if any(c['novel_candidate'] for c in candidates):
        return ('novel_provisional_candidate', '미등록 반복 잠정 후보: 중복·flank 근거 확보, 말단 검증 필요', 'validate_terminal_location')
    if any(c.get('geometry_nomination') for c in candidates):
        return ('novel_geometry_nomination', '미등록 반복 경계 후보: 추가 검증 근거 부족', 'complete_candidate_evidence')
    if any(c['candidate_status'] == 'published_candidate_pattern_supported' for c in candidates):
        return ('published_candidate_supported', '문헌 후보와 일치: 독립 검증 필요', 'validate_independent_biosample')
    references_observed = _ranked([c for c in candidates
        if c['reference_status'] in {'established', 'published_candidate'}
        and sequence_observation(c)['representative_eligible']])
    if references_observed:
        return _observation_classification(references_observed[0])
    if any(c.get('measurement_disagreement') for c in candidates):
        return ('detection_inconclusive_measurement', '판정 보류: 탐지 경로 간 관측 불일치',
                'review_measurement_disagreement')
    if _expected_reference_observations(candidates, assignment, include_trace=True):
        return ('published_repeat_trace_observed', '예상 문헌 반복의 소수 서열 관측·판정 보류',
                'review_trace_observation')
    status = adequacy['status']
    if status == 'insufficient_data':
        return ('detection_inconclusive_depth', '미검출: 검출력 부족', 'obtain_more_suitable_wgs')
    if status == 'conditional_power_supported':
        power_scope = expected_family_power_scope(adequacy, assignment)
        if not power_scope['conditional_nondetection_scope_supported']:
            return ('detection_inconclusive_search_scope', '미검출 판단 보류: 예상 TRM의 탐색·검출력 범위 미확인',
                    'resolve_expected_family_search_scope')
        if power_scope['status'] == 'no_expected_family_assigned':
            return ('canonical_not_detected_conditional_power', '조건부 검출력 지지·검사 주기 내 반복 미검출',
                    'review_scanned_period_nondetection')
        return ('canonical_not_detected_conditional_power', '조건부 검출력 지지·모델 범위 내 예상 TRM 미검출', 'validate_alternative_telomere_hypothesis')
    return ('detection_inconclusive', '미검출: 검출력/라이브러리 정보 부족', 'resolve_adequacy_then_validate')
