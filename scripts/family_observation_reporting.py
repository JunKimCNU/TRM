"""Describe saved repeat families without changing detector or biological policy.

Exact identity uses primitive units and cyclic/reverse-complement equivalence.
The optional one-edit neighborhood is deliberately limited to primitive periods
3..64, independently of the de novo scan's configured min_k/max_k. All saved
rows remain traceable; overlapping sources are never added together.
"""
from __future__ import annotations

from copy import deepcopy

NEARBY_MIN_PERIOD = 3
NEARBY_MAX_PERIOD = 64
NEARBY_DISPLAY_LIMIT = 5
_COMPLEMENT = str.maketrans('ACGT', 'TGCA')
_COUNT_SOURCES = {'measured_evidence', 'targeted_only', 'expected_family_diagnostic'}
_COUNT_FIELDS = ('full_repeat_reads', 'repeat_bases', 'max_tract', 'max_tract_bp')


def _count(value):
    return value if isinstance(value, int) and not isinstance(value, bool) and value >= 0 else None


def _unit(value):
    if not isinstance(value, str):
        return None
    value = value.strip().upper()
    if not value or set(value) - set('ACGT'):
        return None
    # Prefix function finds the primitive period in linear time, including powers.
    prefix = [0] * len(value)
    for i in range(1, len(value)):
        j = prefix[i - 1]
        while j and value[i] != value[j]:
            j = prefix[j - 1]
        if value[i] == value[j]:
            j += 1
        prefix[i] = j
    period = len(value) - prefix[-1]
    return value[:period] if len(value) % period == 0 else value


def _rotations(value):
    for strand, sequence in (('+', value), ('-', value.translate(_COMPLEMENT)[::-1])):
        for offset in range(len(value)):
            yield strand, offset, sequence[offset:] + sequence[:offset]


def _family(value):
    unit = _unit(value)
    return min(rotation for _, _, rotation in _rotations(unit)) if unit else None


def _linear_one_edit(expected, observed):
    """Return one linear edit's witness, or None; O(length) and no DP matrix."""
    if abs(len(expected) - len(observed)) > 1:
        return None
    i = 0
    while i < min(len(expected), len(observed)) and expected[i] == observed[i]:
        i += 1
    if len(expected) == len(observed):
        if i == len(expected) or expected[i + 1:] != observed[i + 1:]:
            return None
        operation = 'substitution'
        aligned_expected, aligned_observed = expected, observed
    elif len(observed) == len(expected) + 1:
        if expected[i:] != observed[i + 1:]:
            return None
        operation = 'insertion'
        aligned_expected, aligned_observed = expected[:i] + '-' + expected[i:], observed
    else:
        if expected[i + 1:] != observed[i:]:
            return None
        operation = 'deletion'
        aligned_expected, aligned_observed = expected, observed[:i] + '-' + observed[i:]
    return {'operation': operation, 'edit_position0': i,
            'aligned_expected_unit': aligned_expected, 'aligned_observed_unit': aligned_observed}


def _nearby_witness(expected, observed):
    if expected == observed or abs(len(expected) - len(observed)) > 1:
        return None
    if not all(NEARBY_MIN_PERIOD <= len(s) <= NEARBY_MAX_PERIOD for s in (expected, observed)):
        return None
    # With at most one edit, fixing the expected cyclic cut loses no solution:
    # the corresponding observed cut can always be rotated to that boundary.
    # A deletion at the cut is represented at position0; insertion can be last.
    for strand, offset, rotated in _rotations(observed):
        witness = _linear_one_edit(expected, rotated)
        if witness:
            return {**witness, 'edit_distance': 1, 'expected_unit': expected,
                    'observed_unit': rotated, 'expected_strand': '+', 'expected_offset0': 0,
                    'observed_strand': strand, 'observed_offset0': offset,
                    'orientation_basis': 'Primitive canonical family strings, not biological chromosome strands.'}
    return None


def _observation_row(row, bucket, index, scope, *, count_field='reads', reported=True,
                     zero_measured=False):
    motif = row.get('canonical_family') or row.get('motif')
    family = _family(motif)
    value = row.get(count_field)
    numeric = _count(value)
    supplied_family = _family(row.get('canonical_family'))
    motif_family = _family(row.get('motif'))
    identity_conflict = bool(supplied_family and motif_family and supplied_family != motif_family)
    return {'canonical_family': family, 'observed_motif': row.get('motif', motif),
            'primitive_unit': _unit(motif), 'source_bucket': bucket, 'source_index': index,
            'source_locator': f'{bucket}[{index}].{count_field}', 'count_scope': scope,
            'count_field': count_field, 'reads': numeric, 'reported_value': deepcopy(value),
            'count_observed': reported and numeric is not None,
            'explicit_zero_measurement': reported and numeric == 0 and zero_measured,
            'invalid_read_count': value is not None and numeric is None,
            'identity_conflict': identity_conflict,
            'measurement_disagreement': bool(row.get('measurement_disagreement')) or row.get('recovery') == 'measurement_disagreement',
            'source': deepcopy(row.get('source')), 'recovery': row.get('recovery'),
            'evidence_count_source': row.get('evidence_count_source'),
            'candidate_status': row.get('candidate_status'),
            **{key: deepcopy(row.get(key)) if count_field == 'reads' else None for key in _COUNT_FIELDS}}


def _view(rows):
    values = sorted({r['reads'] for r in rows if r['count_observed']})
    differing_fields = []
    for field in ('reads', 'full_repeat_reads', 'repeat_bases', 'max_tract_bp'):
        observed = set()
        for row in rows:
            value = row.get('max_tract_bp') if field == 'max_tract_bp' else row.get(field)
            if field == 'max_tract_bp' and value is None:
                value = row.get('max_tract')
            if row['count_observed'] and _count(value) is not None:
                observed.add(value)
        if len(observed) > 1:
            differing_fields.append(field)
    conflict = bool(differing_fields) or any(r['measurement_disagreement'] or r['identity_conflict'] for r in rows)
    explicit_zero = any(r['explicit_zero_measurement'] for r in rows)
    status = ('measurement_disagreement' if conflict else 'reported_positive' if any(values)
              else 'explicit_zero' if values == [0] and explicit_zero
              else 'reported_zero_unverified' if values == [0]
              else 'invalid_counts' if any(r['invalid_read_count'] for r in rows)
              else 'not_reported')
    return {'status': status, 'reads': values[0] if len(values) == 1 else None,
            'read_counts': values, 'positive_observed': any(v > 0 for v in values),
            'explicit_zero_observed': explicit_zero, 'reported_zero_observed': 0 in values,
            'conflicting_count_fields': differing_fields, 'rows': rows}


def _group(family, rows):
    current = _view([r for r in rows if r['count_scope'] == 'current'])
    original = _view([r for r in rows if r['count_scope'] == 'original'])
    cross_difference = (bool(current['read_counts'] and original['read_counts'])
                        and current['read_counts'] != original['read_counts'])
    disagreement = (current['status'] == 'measurement_disagreement'
                    or original['status'] == 'measurement_disagreement' or cross_difference)
    return {'canonical_family': family, 'primitive_period': len(family),
            'current': current, 'original': original, 'rows': rows,
            'positive_observed': current['positive_observed'] or original['positive_observed'],
            'measurement_disagreement': disagreement,
            'counts_differ_between_current_and_original': cross_difference,
            'sources': [r['source_locator'] for r in rows],
            'buckets': sorted({r['source_bucket'] for r in rows})}


def _motifs(value):
    if isinstance(value, str):
        return [m.strip() for m in value.split(',') if m.strip()]
    return [m for m in value or [] if isinstance(m, str) and m.strip()]


def build_family_observation_report(raw, candidates, expected_observations, assignment, *,
                                    source_available=True, analysis_complete=True):
    """Add an observation view; do not mutate input, rank, or classify a TRM.

    Expected observations are authoritative diagnostic rows. related_motifs in
    the assignment also supplies expectations; unjoined review_motifs do not.
    Current assessed zero placeholders without measured provenance are not zeros.
    Exact matches remain visible even outside the actual de novo period range;
    registered targeted evidence can measure such families separately.
    """
    raw, assignment = raw or {}, assignment or {}
    candidates, expected_observations = candidates or [], expected_observations or []
    by_family, invalid_rows = {}, []

    def add(row):
        if row['canonical_family'] is None:
            invalid_rows.append(row)
        else:
            by_family.setdefault(row['canonical_family'], []).append(row)

    for bucket in ('candidates', 'rejected_candidates', 'candidate_evidence'):
        for index, row in enumerate(raw.get(bucket, []) or []):
            add(_observation_row(row, f'raw.{bucket}', index, 'original',
                                 zero_measured=raw.get('analysis_kind') == 'targeted_evidence'))
    raw_families = set(by_family)
    for index, row in enumerate(candidates):
        family = _family(row.get('canonical_family') or row.get('motif'))
        count = _count(row.get('reads'))
        reported = (bool(count) or row.get('evidence_count_source') in _COUNT_SOURCES or family in raw_families)
        add(_observation_row(row, 'assessed_candidates', index, 'current', reported=reported,
                             zero_measured=row.get('evidence_count_source') in _COUNT_SOURCES))
        historical = _count(row.get('original_discovery_reads'))
        if historical and historical > 0 or row.get('original_discovery_reported') is True:
            add(_observation_row(row, 'assessed_candidates', index, 'original', count_field='original_discovery_reads'))
    for index, row in enumerate(expected_observations):
        # not_reported explicitly overrides numeric placeholder0, not positive evidence.
        reported = row.get('recovery') != 'not_reported' or bool(_count(row.get('reads')))
        add(_observation_row(row, 'expected_observations', index, 'current', reported=reported,
                             zero_measured=(row.get('recovery') == 'explicit_zero'
                                            or row.get('source') in _COUNT_SOURCES
                                            or row.get('evidence_count_source') in _COUNT_SOURCES)))
        if _count(row.get('original_discovery_reads')) or row.get('original_discovery_reported') is True:
            add(_observation_row(row, 'expected_observations', index, 'original', count_field='original_discovery_reads'))
    groups = {family: _group(family, rows) for family, rows in by_family.items()}
    positives = {family: group for family, group in groups.items() if group['positive_observed']}

    expected = {}
    for row in expected_observations:
        motif = row.get('motif') or row.get('canonical_family')
        family = _family(row.get('canonical_family') or motif)
        key = family if family is not None else ('invalid', str(motif))
        item = expected.setdefault(key, {'canonical_family': family, 'expected_motifs': [], 'reference_rows': []})
        if motif not in item['expected_motifs']:
            item['expected_motifs'].append(motif)
        item['reference_rows'].append(deepcopy(row))
    for motif in _motifs(assignment.get('related_motifs')):
        family = _family(motif)
        key = family if family is not None else ('invalid', motif)
        item = expected.setdefault(key, {'canonical_family': family, 'expected_motifs': [], 'reference_rows': []})
        if motif not in item['expected_motifs']:
            item['expected_motifs'].append(motif)
    config = raw.get('config', {}) or {}
    scan_min, scan_max = _count(config.get('min_k')), _count(config.get('max_k'))
    summary = raw.get('summary', {}) or {}
    context_status = ('source_unavailable' if not source_available else 'partial_input' if summary.get('sampled_prefix') is True
                      else 'analysis_incomplete' if not analysis_complete else 'complete_input' if summary.get('sampled_prefix') is False
                      else 'input_completeness_unverified')
    review_families = {_family(m) for m in _motifs(assignment.get('review_motifs'))}
    related_families = {_family(m) for m in _motifs(assignment.get('related_motifs'))}
    families = []
    for item in expected.values():
        family = item['canonical_family']
        exact = groups.get(family) if family else None
        exact = exact or _group(family or '', [])
        nearby = []
        if family:
            for observed, group in positives.items():
                witness = _nearby_witness(family, observed)
                if witness:
                    nearby.append({**group, 'edit_distance': 1, 'operation': witness['operation'],
                                   'edit_witness': witness,
                                   'interpretation': 'One edit between primitive cyclic units; descriptive similarity only, not a canonical or novel TRM assignment.'})
        nearby.sort(key=lambda r: (r['edit_distance'], r['canonical_family']))
        nearby_set = {r['canonical_family'] for r in nearby}
        other_families = sorted(set(positives) - nearby_set - {family})
        valid_unit = bool(family)
        comparison_supported = valid_unit and NEARBY_MIN_PERIOD <= len(family) <= NEARBY_MAX_PERIOD
        if exact['positive_observed']:
            status, label = 'exact_family_observed', '예상 primitive family와 정확 일치하는 반복 관측'
        elif nearby:
            status, label = 'nearby_family_observed', '예상 family에서 한 염기 차이인 반복 관측 · 생물학적 판정 아님'
        elif not valid_unit or len(family) > NEARBY_MAX_PERIOD:
            status, label = 'unsupported_expected_unit', '예상 반복 단위가 비교 범위 밖 · 미검출 판정 아님'
        elif other_families:
            status, label = 'other_repeat_families_observed', '다른 반복 family 관측 · 예상 family와 구분'
        elif not source_available:
            status, label = 'unavailable', '원본 결과 없음 · 관측 평가 불가'
        elif not analysis_complete or summary.get('sampled_prefix') is True:
            status, label = 'analysis_incomplete', '분석 미완료/일부 입력 · 관측 평가 제한'
        else:
            status, label = 'no_positive_candidates_reported', '양성 후보 미보고 · 반복서열의 생물학적 부재 아님'
        exact_values = exact['current']['read_counts'] + exact['original']['read_counts']
        measurement = ('measurement_disagreement' if exact['measurement_disagreement']
                       else 'reported_positive' if exact['positive_observed']
                       else 'explicit_zero' if exact['current']['explicit_zero_observed'] or exact['original']['explicit_zero_observed']
                       else 'reported_zero_unverified' if 0 in exact_values
                       else 'invalid_counts' if 'invalid_counts' in (exact['current']['status'], exact['original']['status'])
                       else 'not_reported')
        if exact['measurement_disagreement']:
            label += ' · 관측 계수 불일치 보존'
        rows = item['reference_rows']
        # assess_candidate.reference_status is catalogue membership, whereas
        # assignment.reference_status concerns the scoped species/clade claim.
        # They are different evidence axes and must not override each other.
        reference_status = assignment.get('reference_status')
        catalogue_statuses = sorted({str(r['reference_status']) for r in rows if r.get('reference_status')})
        reference_scope = next((r.get('reference_scope') for r in rows if r.get('reference_scope')), assignment.get('reference_scope'))
        role = next((r.get('reference_role') for r in rows if r.get('reference_role')), None)
        role = role or ('review_candidate' if family in review_families and family not in related_families else 'expected_reference')
        outside = (bool(len(family) < scan_min or len(family) > scan_max)
                   if family and scan_min is not None and scan_max is not None else None)
        families.append({**item, 'primitive_period': len(family) if family else None,
                         'status': status, 'display_label': label, 'reference_status': reference_status,
                         'candidate_catalog_statuses': catalogue_statuses,
                         'reference_scope': reference_scope, 'reference_role': role,
                         'exact_measurement_status': measurement, 'measurement_context_status': context_status,
                         'exact_current': exact['current'], 'exact_original': exact['original'],
                         'exact_rows': exact['rows'], 'measurement_disagreement': exact['measurement_disagreement'],
                         'nearby_families': nearby[:NEARBY_DISPLAY_LIMIT], 'nearby_families_all': nearby,
                         'nearby_family_count': len(nearby), 'nearby_displayed_count': min(len(nearby), NEARBY_DISPLAY_LIMIT),
                         'nearby_omitted_count': max(0, len(nearby) - NEARBY_DISPLAY_LIMIT),
                         'nearby_truncated': len(nearby) > NEARBY_DISPLAY_LIMIT,
                         'nearby_order': 'edit_distance_then_primitive_canonical_family',
                         'nearby_comparison_status': 'assessed' if comparison_supported else 'unsupported_period_or_invalid_unit',
                         'other_positive_families': other_families, 'other_positive_family_count': len(other_families),
                         'outside_de_novo_search_range': outside,
                         'expected_outside_native_period_range': outside,
                         'targeted_registration_may_measure_outside_de_novo_range': True})
    families.sort(key=lambda r: (r['canonical_family'] or '', str(r['expected_motifs'])))
    statuses = {r['status'] for r in families}
    if not families:
        status, label = 'no_reference', '문헌 예상 family 미지정 · 미검출 판정 아님'
    elif len(statuses) == 1:
        status, label = families[0]['status'], families[0]['display_label']
    else:
        status, label = 'mixed_reference_observations', '예상 family별 관측 결과가 다름 · 개별 행 확인'
    return {'families': families, 'status': status, 'display_label': label,
            'reference_family_count': len(families), 'positive_saved_family_count': len(positives),
            'positive_saved_families': sorted(positives), 'invalid_observation_rows': invalid_rows,
            'measurement_context_status': context_status,
            'reference_context': {k: deepcopy(assignment.get(k)) for k in ('reference_status', 'reference_scope', 'mapping_status', 'repeat_architecture', 'reference_ids', 'review_reference_ids')},
            'comparison_context': {
                'exact_identity': 'Primitive unit plus all rotations and reverse complement; observations remain separate by source and epoch.',
                'nearby_metric': 'At most one substitution/insertion/deletion between primitive units over rotations and reverse complements; exact identities excluded from nearby list.',
                'nearby_min_primitive_period': NEARBY_MIN_PERIOD, 'nearby_max_primitive_period': NEARBY_MAX_PERIOD,
                'nearby_display_limit': NEARBY_DISPLAY_LIMIT,
                'native_config': {k: deepcopy(config.get(k)) for k in ('mode', 'min_k', 'max_k', 'seed_bases', 'min_copies', 'min_identity', 'min_reads', 'min_tract')},
                'native_analysis_kind': raw.get('analysis_kind'),
                'scan_range_interpretation': 'min_k/max_k describe de novo search, not registered targeted-measurement eligibility. The fixed reporting comparison budget is separate.',
                'search_scope': 'All saved raw primary/rejected/targeted rows, all current assessed rows, and positive historical observations; no top-N input filtering.',
                'counts': 'Never sum duplicate, current, historical, or diagnostic rows. reads is null when distinct values occur within a scope; read_counts and source rows retain them.'},
            'observation_scope': 'Descriptive sequence-family observations only. No canonical/novel promotion, terminal-location claim, adequacy inference, candidate ranking or classification change.'}
