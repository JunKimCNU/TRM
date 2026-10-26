"""Compact simultaneous evidence views; no detector or policy decisions here."""
from __future__ import annotations

from clade_variant_reporting import catalogue_status


AXIS_FIELDS = ['expected_observation_axis', 'catalogued_observation_axis',
               'alternative_candidate_axis', 'terminal_confirmation_axis',
               'positive_interpretation_axis', 'negative_power_axis']


def build_evidence_axes(document, *, global_catalogue=None):
    """Project existing evidence without choosing a winner or merging counts.

    The family observation module owns zero/missing/conflict semantics. Original
    and current measurements remain separate, including a current zero after an
    original positive. Catalogue identity is never a taxon canonical assignment.
    """
    candidates = document.get('candidates', [])
    by_family = {c['canonical_family']: c for c in candidates}
    expected = []
    for family in document.get('family_observation', {}).get('families', []):
        current, original = family.get('exact_current') or {}, family.get('exact_original') or {}
        candidate = by_family.get(family['canonical_family'], {})
        expected.append({
            'canonical_family': family['canonical_family'],
            'expected_motifs': family.get('expected_motifs', []),
            'reference_scope': family.get('reference_scope'),
            'reference_status': family.get('reference_status'),
            'reference_role': family.get('reference_role'),
            'observation_status': family.get('exact_measurement_status'),
            'current_reads': current.get('reads'), 'original_reads': original.get('reads'),
            'current_read_counts': current.get('read_counts', []),
            'original_read_counts': original.get('read_counts', []),
            'current_positive': current.get('positive_observed', False),
            'original_positive': original.get('positive_observed', False),
            'assessed_candidate_full_repeat_reads': candidate.get('full_repeat_reads'),
            'candidate_status': candidate.get('candidate_status'),
        })
    def positive(row):
        value = row.get('reads')
        return isinstance(value, int) and not isinstance(value, bool) and value > 0

    def compact(row):
        return {key: row.get(key) for key in ('motif', 'canonical_family', 'reads',
            'full_repeat_reads', 'reference_status', 'candidate_status',
            'evidence_count_source', 'measurement_disagreement')}

    # Existing policy ordering is retained only for bounded previews. Nothing
    # outside the preview is discarded from the assessment candidate table.
    registered = [c for c in candidates if positive(c) and
                  catalogue_status(c['canonical_family'], global_catalogue) == 'globally_established_family']
    alternatives = [c for c in candidates if positive(c) and
                    (c.get('novel_candidate') is True or c.get('geometry_nomination') is True)]
    library = document.get('library_context') or {}
    adequacy = document.get('data_adequacy') or {}
    agreement = document.get('reference_agreement') or {}
    labels = {'reported_positive': '관측', 'explicit_zero': '측정 0',
              'not_reported': '미보고 (0 아님)', 'measurement_disagreement': '측정 간 불일치',
              'reported_zero_unverified': '0 보고·측정 출처 미확인'}

    def expected_label(row):
        scope = {'clade': '분류군 기본 기대', 'species': '종별 문헌'}.get(
            row['reference_scope'], row['reference_scope'] or '참조 범위 미지정')
        if (row['reference_status'] not in ('supported', 'established')
                or row['reference_role'] == 'review_candidate'):
            scope += '·문헌 후보/미확립'
        current = ('/'.join(map(str, row['current_read_counts'])) + ' (계수 불일치)'
                   if len(row['current_read_counts']) > 1 else
                   '미보고' if row['current_reads'] is None else str(row['current_reads']))
        original = (f"; 원래 {row['original_reads']}" if row['original_reads'] is not None and
                    row['original_reads'] != row['current_reads'] else '')
        return (f"{'/'.join(row['expected_motifs']) or row['canonical_family']} [{scope}]: "
                f"{labels.get(row['observation_status'], row['observation_status'])}, "
                f"현재 {current} reads{original}")

    def preview_label(rows, empty):
        items = [f"{c['motif']} {c['reads']} reads (full {c.get('full_repeat_reads')})" for c in rows[:5]]
        if len(rows) > 5:
            items.append(f'외 {len(rows)-5} family; 전체 후보표 참조')
        return '; '.join(items) if items else empty

    positive_limits = list(library.get('reasons') or [])
    if not document.get('analysis_complete'):
        positive_limits.append('analysis_incomplete')
    if agreement.get('status') in {'qc_limited', 'taxonomic_scope_unresolved', 'reference_out_of_scope'}:
        positive_limits.append(agreement['status'])
    main_status = agreement.get('main_trm_status', 'not_confirmed')
    physical_status = 'externally_confirmed' if document.get('externally_confirmed') is True else 'not_confirmed'
    displays = {
        'expected_observation_axis': '; '.join(map(expected_label, expected)) or '종·분류군 기대 motif 미지정',
        'catalogued_observation_axis': preview_label(registered, '확립된 catalogue family의 현재 양성 관측 없음'),
        'alternative_candidate_axis': preview_label(alternatives, '현재 정책의 잠정·구조 후보 없음'),
        'terminal_confirmation_axis': f'Main TRM: {main_status}; 이 표본의 물리적 말단: {physical_status}',
        'positive_interpretation_axis': '; '.join(positive_limits) or '등록된 추가 제약 없음 · 말단 귀속은 별도 검증',
        'negative_power_axis': adequacy.get('status', 'unassessed') + ' · 반복 관측 여부와 별도',
    }
    if global_catalogue is None:
        displays['catalogued_observation_axis'] = '전역 catalogue 미제공·평가 불가'
    return {'schema_version': 1, 'expected': expected,
            'catalogued_positive_family_count': len(registered),
            'catalogued_positive_preview': [compact(c) for c in registered[:5]],
            'alternative_candidate_count': len(alternatives),
            'alternative_candidate_preview': [compact(c) for c in alternatives[:5]],
            'preview_limit': 5, 'preview_order': 'existing_candidate_policy_order',
            'main_trm_status': main_status, 'sample_physical_terminal_status': physical_status,
            'positive_interpretation_limits': positive_limits,
            'negative_power_status': adequacy.get('status', 'unassessed'),
            'negative_power_missing_assumptions': adequacy.get('missing_assumptions', []),
            'display': displays, 'changes_classification_or_adequacy': False,
            'interpretation': 'Parallel, nonexclusive evidence axes. A global catalogue match is not a taxon canonical or sample terminal assignment. Read counts are not independent molecules. Candidate previews do not replace the full candidate table.'}
