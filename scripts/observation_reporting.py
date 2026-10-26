"""Additive observation display; never selects candidates or changes TRM policy."""
from __future__ import annotations

from evidence_batch import canonical

REPORTING_VERSION = '1.13.0'


def _positive(value):
    return isinstance(value, int) and not isinstance(value, bool) and value > 0


def _observation(candidate):
    """Keep current and original counts separate, including current zero."""
    return {key: candidate.get(key) for key in (
        'motif', 'canonical_family', 'reads', 'full_repeat_reads', 'repeat_bases',
        'max_tract_bp', 'original_discovery_reads', 'evidence_count_source',
        'measurement_disagreement', 'ranking_score', 'candidate_status')}


def build_observation_report(raw, candidates, representative, *, source_available,
                             analysis_complete, classification, legacy_result,
                             agreement, adequacy_status):
    """Views over saved evidence; absence of an emitted family is never zero.

    Native ranks are only recovered from an unambiguous complete set of stored
    ranks over the raw candidates/rejected union. Targeted-only positives are
    shown in the policy observation table without inventing native ranks.
    """
    native_rows = {}
    duplicate_conflict = False
    for candidate in raw.get('candidates', []) + raw.get('rejected_candidates', []):
        family = canonical(candidate['motif'])
        if family in native_rows and native_rows[family] != candidate:
            duplicate_conflict = True
        native_rows.setdefault(family, candidate)
    ranks = [c.get('telomere_prior', {}).get('discovery_rank') for c in native_rows.values()]
    ranks_valid = (bool(ranks) and not duplicate_conflict
                   and all(isinstance(r, int) and not isinstance(r, bool) and r >= 0 for r in ranks)
                   and set(ranks) == set(range(len(native_rows))))
    native_positive = [(family, c) for family, c in native_rows.items() if _positive(c.get('reads'))]
    native_leader = None
    if ranks_valid and native_positive:
        family, row = min(native_positive, key=lambda item: item[1]['telomere_prior']['discovery_rank'])
        rank = row['telomere_prior']['discovery_rank']
        native_leader = {
            'canonical_family': family, 'motif': row['motif'],
            'native_display_annotation': row.get('telomere_prior', {}).get('observed_motif_display'),
            'reads': row.get('reads'), 'repeat_bases': row.get('repeat_bases'),
            'full_repeat_reads': row.get('full_repeat_reads'), 'max_tract_bp': row.get('max_tract'),
            'rank_score': row.get('rank_score'), 'discovery_rank_zero_based': rank,
            'discovery_rank_one_based': rank + 1, 'count_source': 'original_native_discovery',
        }
    positive = [c for c in candidates if _positive(c.get('reads'))]
    historical = [c for c in candidates if _positive(c.get('original_discovery_reads'))]
    historical_nonpositive = [c for c in historical if not _positive(c.get('reads'))]
    ordered = sorted(positive, key=lambda c: (-c['ranking_score'], c['canonical_family']))
    policy_leader = _observation(ordered[0]) if ordered else None
    observed = bool(positive or native_positive or historical)
    sampled = raw.get('summary', {}).get('sampled_prefix')
    if not source_available:
        measurement_status = 'source_unavailable'
        measurement_label = '원본 결과 없음 · 관측 평가 불가'
    elif sampled is True:
        measurement_status = 'partial_input'
        measurement_label = '입력 일부만 분석 · 전체 입력 판정 보류'
    elif not analysis_complete:
        measurement_status = 'analysis_incomplete'
        measurement_label = '분석 미완료 · 남아 있는 관측은 보존'
    elif adequacy_status == 'observed_qc_failure':
        measurement_status = 'qc_limited'
        measurement_label = '측정 완료 · QC 불일치로 해석 보류'
    elif sampled is not False:
        measurement_status = 'input_completeness_unverified'
        measurement_label = '입력 전체 분석 여부 미확인 · 전체 입력으로 간주 불가'
    else:
        measurement_status = 'complete_input_measurement'
        measurement_label = '전체 입력 측정 완료 · 검출력 충분 여부는 별도'
    if positive:
        observation_status, observation_label = 'positive_observations', '반복 서열 관측'
    elif native_positive or historical:
        observation_status, observation_label = 'original_positive_only', '원래 양성 관측 보존 · 현재 측정 상태 확인 필요'
    elif not source_available:
        observation_status, observation_label = 'unavailable', '관측 평가 불가'
    elif candidates and all(c.get('reads') == 0 and c.get('evidence_count_source') in
                            {'measured_evidence', 'targeted_only'} for c in candidates):
        observation_status, observation_label = 'reported_targets_zero', '보고된 표적에서 0건 · 다른 family의 부재는 아님'
    else:
        observation_status, observation_label = 'no_positive_reported', '양성 반복 미보고 · 명시적 전체 0건은 아님'
    # Only replace ambiguous presentation text. The policy classification and
    # legacy result remain untouched, including all failure/library precedence.
    display_label = legacy_result
    inconclusive = {
        'detection_inconclusive': '검출력·라이브러리 정보 부족',
        'detection_inconclusive_depth': '조건부 검출력 부족',
        'detection_inconclusive_search_scope': '기대 family의 탐색·검출력 범위 부족',
    }
    if classification in inconclusive:
        prefix = '반복 서열 관측' if observed else '양성 반복 미보고'
        display_label = prefix + ' · TRM 판정 보류 (' + inconclusive[classification] + ')'
    interpretation_limits = {
        'analysis_incomplete': '분석 미완료 · TRM 해석 보류',
        'qc_limited': 'QC 제약 · TRM 해석 보류',
        'taxonomic_scope_unresolved': '생물학적 출처·분류 범위 확인 필요 · TRM 귀속 보류',
        'library_context_limited': '라이브러리 제약 · genomic TRM 해석 보류',
        'reference_out_of_scope': '문헌상 비정형 말단 구조 · 단일 short TRM 해석 범위 밖',
    }
    display_label = interpretation_limits.get(agreement['status'], display_label)
    policy_rank = next((i + 1 for i, c in enumerate(candidates)
                        if representative and c['canonical_family'] == representative['canonical_family']), None)
    def same_family(left, right):
        return left['canonical_family'] == right['canonical_family'] if left and right else None
    return {
        'repeat_observation': {
            'status': observation_status, 'display_label': observation_label,
            'measurement_status': measurement_status, 'measurement_label': measurement_label,
            'source_available': source_available, 'analysis_complete': analysis_complete,
            'positive_family_count': len(positive), 'original_positive_family_count': len(native_positive),
            'native_discovery_leader': native_leader,
            'native_ranking_status': 'stored_ranks_available' if ranks_valid else 'stored_ranks_unavailable_or_ambiguous',
            'native_ranking_basis': 'Stored native discovery_rank before reference annotation; scanner score/assembly terminal ordering over emitted families, including prior-rejected families. Not an unbiased abundance rank or TRM assignment.',
            'native_rank_scope': 'Original reported families after detector/seed/min_reads thresholds; unreported families are not measured zeros.',
            'policy_score_leader': policy_leader,
            'policy_ranking_basis': 'Python ranking_score over positive current assessment counts: repeat bases, expected-family multiplier and boundary-direction factors. Reference-dependent; not the native discovery order.',
            'top_positive_observations': [_observation(c) for c in ordered[:5]],
            'top_positive_observations_limit': 5,
            'top_positive_observations_order': 'policy_ranking_score_then_canonical_family',
            'original_positive_current_nonpositive': [_observation(c) for c in historical_nonpositive[:5]],
            'original_positive_current_nonpositive_count': len(historical_nonpositive),
            'original_positive_current_nonpositive_preview_limit': 5,
            'original_positive_current_nonpositive_preview_complete': len(historical_nonpositive) <= 5,
            'interpretation': 'Positive repeats may be satellites or internal arrays. No displayed leader assigns a main TRM or establishes physical terminal location. Counts are reads/records, not independent molecules.',
        },
        'policy_selection': {
            'representative_motif': representative.get('motif') if representative else None,
            'representative_family': representative.get('canonical_family') if representative else None,
            'representative_policy_rank': policy_rank,
            'native_vs_policy_score_same_family': same_family(native_leader, policy_leader),
            'native_vs_representative_same_family': same_family(native_leader, representative),
            'policy_score_vs_representative_same_family': same_family(policy_leader, representative),
            'selection_status': 'representative_selected' if representative else 'no_policy_representative',
            'basis': 'Existing classification-specific select_representative policy; may differ from both native discovery leader and Python score leader. No representative does not imply no observed repeat.',
        },
        'decision_summary': {
            'display_label': display_label, 'classification': classification,
            'legacy_result': legacy_result, 'reference_agreement_status': agreement['status'],
            'main_trm_status': agreement['main_trm_status'],
            'main_trm_display': '미확정 · 물리적 말단 검증 필요' if agreement['main_trm_status'] == 'not_confirmed' else agreement['main_trm_status'],
            'main_trm_motif': agreement.get('main_trm_motif'),
            'data_adequacy_status': adequacy_status,
        },
    }
