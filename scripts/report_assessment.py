#!/usr/bin/env python3
"""Create separate current-version interpretation and previous-version comparison.

Join archived de novo evidence, full-input targeted geometry, curated taxonomy,
and explicit adequacy assumptions. Does not edit any input result.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import html
import json
from pathlib import Path
import time
from urllib.parse import urlencode

from assess_adequacy import assess_adequacy
from evidence_batch import canonical, index_samples, read_tsv, sha256, write_json
from expected_diagnostics import supplement_expected_observations
from followup_diagnostics import load_followup_diagnostics
from observation_reporting import REPORTING_VERSION, build_observation_report
from evidence_axes import AXIS_FIELDS, build_evidence_axes
from family_observation_reporting import build_family_observation_report
from input_applicability import build_input_applicability
from technical_boundary_evidence import load_evidence, build_technical_boundary_evidence
import reference_localization_evidence as localization_evidence
import cross_period_evidence
import saved_repeat_structure
import assembly_context_evidence
import quality_path_evidence
import long_read_context_evidence
import long_read_context_report
from native_result_integrity import validate_native_result
from clade_variant_reporting import build_clade_variant_report
from quality_observation import describe_quality, describe_candidate_quality
from trm_assessment import (POLICY, SUPPORTED_EVIDENCE, assess_candidate, catalog_index, classify,
                            reverse_complement, library_context, select_representative,
                            sequence_observation, PAIR_COUNT_FIELDS, reference_agreement,
                            resolved_eukaryotic_scope, expected_family_power_scope)

SAMPLE_FIELDS = ['sample_id', 'scientific_name', 'result', 'detected_motif', 'expected_motifs',
                 'sequence_observation', 'expected_repeat_observation',
                 'reference_status', 'repeat_architecture', 'reference_context',
                 'evidence_level', 'supporting_reads', 'full_repeat_reads', 'left_consistent_junctions',
                 'right_consistent_junctions', 'opposite_junctions', 'direction_consistency',
                 'data_adequacy', 'records', 'bases', 'unmeasured_novel_families',
                 'library_context', 'why', 'confirmation_missing', 'geometry_nominations',
                 'provisional_candidates', 'distinct_junction_sequences', 'physical_flank_clusters',
                 'deduplicated_orientation', 'clean_two_flank_records',
                 'next_action', 'sra_search', 'details']
SAMPLE_FIELDS += ['paired_consistency']
SAMPLE_FIELDS += AXIS_FIELDS
SAMPLE_FIELDS += list(saved_repeat_structure.FIELDS)
SAMPLE_FIELDS += list(assembly_context_evidence.FIELDS)
SAMPLE_FIELDS += list(quality_path_evidence.FIELDS)
SAMPLE_FIELDS += list(long_read_context_evidence.FIELDS)
SAMPLE_FIELDS += ['trm_interpretation', 'reference_agreement_status', 'reference_conflict_review',
                  'main_trm_status', 'reference_review_recommendation']
SAMPLE_FIELDS += ['observation_status', 'measurement_status', 'decision_display',
                  'native_discovery_leader', 'native_discovery_reads',
                  'policy_score_leader', 'policy_score_leader_reads',
                  'policy_representative', 'positive_repeat_families', 'main_trm_display',
                  'native_vs_policy_score_same_family', 'native_vs_representative_same_family',
                  'policy_score_vs_representative_same_family', 'reporting_version']
SAMPLE_FIELDS += ['followup_diagnostic_status', 'followup_diagnostic_display',
                  'followup_diagnostic_scope', 'followup_diagnostic_summary',
                  'followup_diagnostic_recommendations']
SAMPLE_FIELDS += ['family_observation_status', 'family_observation_display',
                  'family_observation_summary']
SAMPLE_FIELDS += ['input_applicability_status', 'input_applicability_display',
                 'input_applicability_scope', 'input_applicability_recommendation']
SAMPLE_FIELDS += ['flank_period_applicability_status', 'flank_period_applicability_display',
                 'flank_period_applicability_scope', 'flank_period_applicability_recommendation']
SAMPLE_FIELDS += ['technical_boundary_status', 'technical_boundary_display',
                 'technical_boundary_scope', 'technical_boundary_summary',
                 'technical_boundary_recommendation', 'technical_boundary_selected_events',
                 'technical_boundary_flagged_events']
LOCALIZATION_COUNTS = ('selected_event_count', 'selected_pair_count', 'assessed_event_count',
    'coordinate_supported_event_count', 'boundary_supported_event_count',
    'mate_geometry_evaluable_event_count', 'inward_mate_event_count',
    'capped_event_count', 'reference_count')
SAMPLE_FIELDS += ['reference_localization_status', 'reference_localization_display',
                 'reference_localization_scope', 'reference_localization_summary']
SAMPLE_FIELDS += ['reference_localization_' + key for key in LOCALIZATION_COUNTS]
REFERENCE_LOCALIZATION_FIELDS = ['sample_id', 'scientific_name', 'status', 'reason',
    'display_label', 'scope', 'summary', 'recommendation', *LOCALIZATION_COUNTS]
REFERENCE_LOCALIZATION_TABLES = ('reference_localization_diagnostics.tsv',
    'reference_localization_families.tsv', 'reference_localization_references.tsv',
    'reference_localization_events.tsv')
CROSS_PERIOD_COUNTS = cross_period_evidence.SUMMARY_FIELDS
SAMPLE_FIELDS += ['cross_period_status', 'cross_period_display',
                 'cross_period_scope', 'cross_period_summary']
SAMPLE_FIELDS += ['cross_period_' + key for key in CROSS_PERIOD_COUNTS]
CROSS_PERIOD_FIELDS = ['sample_id', 'scientific_name', 'status', 'reason',
    'display_label', 'scope', 'summary', 'recommendation', *CROSS_PERIOD_COUNTS]
CROSS_PERIOD_TABLES = ('cross_period_diagnostics.tsv', 'cross_period_families.tsv',
    'cross_period_intervals.tsv', 'cross_period_events.tsv',
    'cross_period_contexts.tsv', 'cross_period_positions.tsv')
ASSEMBLY_CONTEXT_FIELDS = ('sample_id', 'scientific_name', 'status', 'scope',
                          *assembly_context_evidence.SUMMARY_FIELDS)
ASSEMBLY_CONTEXT_TABLES = ('assembly_context_summary.tsv', 'assembly_context_scaffolds.tsv',
                           'assembly_context_intervals.tsv', 'assembly_context_events.tsv')
QUALITY_PATH_FIELDS = ('sample_id', 'scientific_name', 'status', 'scope',
                       *quality_path_evidence.SUMMARY_FIELDS)
QUALITY_PATH_TABLES = ('quality_path_summary.tsv', 'quality_path_models.tsv',
                       'quality_path_read_comparisons.tsv')
LONG_READ_CONTEXT_TABLES = long_read_context_report.TABLES
TECHNICAL_BOUNDARY_FIELDS = ['sample_id', 'scientific_name', 'status', 'scope',
    'selected_event_count', 'selected_pair_count', 'full_read_count', 'flagged_event_count',
    'threshold_bp', 'exact_match_associations', 'full_reads_with_matches',
    'pairs_with_matches', 'events_with_flank_overlap']
TECHNICAL_EVENT_FIELDS = ['sample_id', 'event_id', 'original_pair_index_1based', 'mate',
    'family', 'flank_bp', 'maximum_exact_match_bp_overlapping_flank',
    'paired_readthrough_compatible', 'threshold_bp', 'union_flank_overlap_bp',
    'adjacent_unmatched_flank_bp']
INPUT_APPLICABILITY_FIELDS = ['sample_id', 'scientific_name', 'kind', 'scope',
    'input_status', 'status', 'display_label', 'canonical_family', 'expected_motifs',
    'period_bp', 'minimum_seed_bp', 'length_eligible_reads_lower',
    'length_eligible_reads_upper', 'outside_discovery_range',
    'analyzed_records_length_status', 'bounds_method']
FLANK_APPLICABILITY_FIELDS = ['sample_id', 'scientific_name', 'kind', 'scope',
    'input_status', 'status', 'display_label', 'annotation_enabled', 'model_status',
    'canonical_family', 'expected_motifs', 'period_bp', 'minimum_seed_bp',
    'required_window_bp', 'junction_flank_bp', 'minimum_read_bp',
    'length_eligible_reads_lower', 'length_eligible_reads_upper', 'outside_discovery_range',
    'analyzed_records_length_status', 'bounds_method']
DIRECTIONAL_CLUSTER_FIELDS = ['directional_cluster_measurement', 'consistent_flank_clusters',
    'opposed_flank_clusters', 'retained_consistent_cluster_concentration',
    'flank_context_completeness', 'flank_orientation_basis']
SAMPLE_FIELDS += DIRECTIONAL_CLUSTER_FIELDS
SAMPLE_FIELDS += ['clade_variant_status', 'clade_variant_display', 'clade_variant_families',
                 'clade_variant_review_required', 'published_clade_exception_families']
SAMPLE_FIELDS += ['quality_observation_status', 'quality_observation_display']
QUALITY_FIELDS = ['sample_id', 'canonical_family', 'motif', 'status', 'encoding',
                  'supporting_records', 'quality_records', 'union_tract_bases',
                  'measured_tract_quality_bases', 'tract_bases_without_quality',
                  'tract_bases_below_q20', 'tract_bases_below_q30', 'source_sha256',
                  'variant_discriminating_bases_assessed', 'filtering_applied']
CLADE_VARIANT_FIELDS = ['sample_id', 'scientific_name', 'canonical_family', 'motif',
    'reads', 'original_discovery_reads', 'status', 'global_catalogue_status',
    'clade_relationship', 'scoped_reference_relationship', 'published_scoped_exception',
    'known_repeat_pattern_supported', 'existing_candidate_nomination',
    'sequence_only_observation', 'comparison_limited_by_library', 'is_representative',
    'review_required', 'review_reason', 'next_action', 'validation_needed',
    'biological_novelty_claimed', 'terminal_location_confirmed']
FAMILY_FIELDS = ['sample_id', 'scientific_name', 'expected_family', 'expected_motifs',
                 'reference_status', 'reference_scope', 'candidate_catalog_statuses',
                 'status', 'display_label',
                 'exact_measurement_status', 'current_status', 'original_status',
                 'current_reads', 'original_reads',
                 'current_read_counts', 'original_read_counts', 'measurement_disagreement',
                 'nearby_family_count', 'nearby_displayed_count', 'nearby_truncated',
                 'other_positive_family_count', 'outside_de_novo_search_range',
                 'measurement_context_status', 'nearby_comparison_status']
NEARBY_FIELDS = ['sample_id', 'scientific_name', 'expected_family', 'observed_family',
                 'edit_distance', 'operation', 'current_reads', 'original_reads',
                 'current_read_counts', 'original_read_counts', 'measurement_disagreement',
                 'sources', 'buckets', 'edit_witness']
FOLLOWUP_FIELDS = ['sample_id', 'scientific_name', 'status', 'display',
                   'metric', 'label', 'value', 'unit', 'scope']
COMPARISON_FIELDS = ['sample_id', 'scientific_name', 'previous_call', 'previous_canonical_key',
                     'previous_reference_anchor', 'previous_motif_display', 'updated_result',
                     'updated_motif_display', 'updated_canonical_key', 'motif_family_changed',
                     'interpretation_changed', 'reason', 'previous_result', 'updated_result_file']
CANDIDATE_FIELDS = ['sample_id', 'motif', 'canonical_family', 'reference_status', 'expected_family_match',
                    'candidate_status', 'reads', 'full_repeat_reads', 'repeat_bases', 'max_tract_bp',
                    'left_consistent_junctions', 'right_consistent_junctions', 'opposite_junctions',
                    'direction_consistency', 'direction_status', 'orientation_basis',
                    'full_to_junction_ratio', 'background_status', 'reference_sources', 'notes']
CANDIDATE_FIELDS += ['novel_evidence_level', 'geometry_nomination', 'novel_candidate',
                     'confirmation_status', 'confirmation_missing', 'distinct_junction_sequences',
                     'physical_flank_clusters', 'deduplicated_orientation', 'clean_two_flank_records',
                     'externally_confirmed']
CANDIDATE_FIELDS += ['pair_measurement_status', 'pair_count_source', 'pair_layout', 'pair_scope']
CANDIDATE_FIELDS += ['pair_'+key for key in PAIR_COUNT_FIELDS] + ['pair_one_mate_only_pairs']
CANDIDATE_FIELDS += DIRECTIONAL_CLUSTER_FIELDS
PREVIOUS_ASSESSMENT_FIELDS = ['sample_id', 'scientific_name', 'comparison_status', 'previous_version',
    'updated_version', 'previous_classification', 'updated_classification', 'previous_result',
    'updated_result', 'previous_motif', 'updated_motif', 'previous_canonical_key', 'updated_canonical_key',
    'motif_family_changed', 'classification_changed', 'previous_geometry_nominations',
    'updated_geometry_nominations', 'updated_provisional_candidates', 'confirmation_missing',
    'reason', 'previous_assessment_file', 'previous_assessment_sha256', 'updated_result_file']
PREVIOUS_ASSESSMENT_FIELDS += ['previous_reference_agreement', 'updated_reference_agreement',
                              'reference_interpretation_added', 'reference_interpretation_changed']
PREVIOUS_ASSESSMENT_FIELDS += ['previous_schema_version', 'updated_schema_version',
                              'previous_followup_status', 'updated_followup_status',
                              'followup_diagnostics_added']
DISPLAY_NAMES = {
    'quality_path_status': '같은 read 품질·최적경로 연결 상태',
    'quality_path_read_display': '선정 read의 대체 품질 비교',
    'quality_path_scope_display': '집중 혼합 구간·read·pair 분모',
    'quality_path_model_display': '최적경로 품질 제약·모호성',
    'quality_path_limit_display': '품질·gap·물리적 말단 한계',
    'assembly_context_status': '저장 assembly 문맥 연결 상태',
    'assembly_context_read_display': '선정 read의 혼합 구간 보존',
    'assembly_context_extension_display': '같은 표본 assembly 문법 확장',
    'assembly_context_origin_display': 'Scaffold 기원·구획',
    'assembly_context_quality_display': '원 quality 접근·비교',
    'assembly_context_limit_display': '말단·main/novel 확인 한계',
    'cross_period_display': '서로 다른 반복 주기와 같은 위치 근거',
    'cross_period_scope': '반복 주기 진단 범위',
    'cross_period_summary': '반복 주기 진단 요약',
    'reference_localization_display': '선정 경계의 reference 위치 비교',
    'reference_localization_scope': 'Reference 비교 범위',
    'reference_localization_summary': 'Reference 비교 요약',
    'quality_observation_display': '염기 quality 관측',
    'clade_variant_display': '분류군 canonical과 다른 반복',
    'clade_variant_families': '분류군 variant 관측 family 수',
    'clade_variant_review_required': '분류군 variant 후속 검토',
    'technical_boundary_display': '기술 서열과 경계의 관계',
    'technical_boundary_scope': '기술 서열 진단 범위',
    'technical_boundary_summary': '기술 서열 진단 요약',
    'technical_boundary_recommendation': '경계 검증 권고',
    'input_applicability_status': '입력 길이 적합성 코드',
    'input_applicability_display': '현재 탐색 조건의 입력 길이 적합성',
    'input_applicability_scope': '입력 길이 평가 범위',
    'input_applicability_recommendation': '입력 길이 확인 권고',
    'flank_period_applicability_status': 'Flank 주기 평가 길이 코드',
    'flank_period_applicability_display': 'Flank 주기 평가의 길이 필요조건',
    'flank_period_applicability_scope': 'Flank 길이 진단 범위',
    'flank_period_applicability_recommendation': 'Flank 길이 확인 권고',
    'directional_cluster_measurement': '방향별 flank 군집 측정',
    'consistent_flank_clusters': '일치 방향 사건이 있는 flank 군집',
    'opposed_flank_clusters': '반대 방향 사건이 있는 flank 군집',
    'retained_consistent_cluster_concentration': '보관 일치 사건의 최대 군집 집중도',
    'flank_context_completeness': 'Flank 사건 보관 범위',
    'flank_orientation_basis': 'Flank 방향 기준',
    'saved_structure_status': '저장된 구조·quality 연결',
    'saved_quality_display': '전체 입력 quality 분포의 추가 관측',
    'saved_mixed_display': '혼합 배열의 추가 관측',
    'expected_observation_axis': '기대 motif 관측·참조 범위',
    'catalogued_observation_axis': '전역 등록 motif 관측 (종 canonical 지정 아님)',
    'alternative_candidate_axis': '함께 관측된 잠정·구조 후보',
    'terminal_confirmation_axis': 'Main TRM·표본 말단 확인',
    'positive_interpretation_axis': '양성 관측의 귀속 제약',
    'negative_power_axis': '미검출 해석의 검출력',
    'family_observation_status': '기대 family 관측 상태',
    'family_observation_display': '기대·근접 family 관측',
    'family_observation_summary': '기대 family별 관측 요약',
    'followup_diagnostic_display': '선택 구간 후속 진단',
    'followup_diagnostic_summary': '후속 진단 요약',
    'followup_diagnostic_scope': '후속 진단 적용 범위',
    'followup_diagnostic_recommendations': '후속 진단 확인 권고',
    'previous_schema_version': '이전 보고서 schema', 'updated_schema_version': '새 보고서 schema',
    'previous_followup_status': '이전 후속 진단 상태', 'updated_followup_status': '새 후속 진단 상태',
    'followup_diagnostics_added': '선택 구간 진단 추가',
    'observation_status': '반복 서열 관측', 'measurement_status': '입력 측정 상태',
    'decision_display': 'TRM 판정', 'native_discovery_leader': 'Native 탐색 우선 반복 (원래 관측)',
    'native_discovery_reads': 'Native 우선 반복 reads (원래 관측)',
    'policy_score_leader': '해석 점수 1위 반복 (현재 계수)',
    'policy_score_leader_reads': '해석 점수 1위 reads (현재 계수)',
    'policy_representative': '정책상 대표 반복', 'positive_repeat_families': '현재 양성 반복 family 수',
    'main_trm_display': 'Main TRM 확정 상태',
    'sample_id': '샘플', 'scientific_name': '종 이름', 'result': '관측 근거 분류',
    'detected_motif': '대표 관측 motif', 'expected_motifs': '문헌상 기대 motif',
    'sequence_observation': '대표 서열 관측', 'expected_repeat_observation': '문헌 예상 서열 관측',
    'reference_status': '분류군 문헌 상태', 'repeat_architecture': '문헌상 반복 구조',
    'reference_context': '분류군 문헌 해석', 'unmeasured_novel_families': '경계 미측정 미등록 families',
    'evidence_level': '관측 근거', 'supporting_reads': '지지 reads',
    'full_repeat_reads': '반복으로 덮인 reads', 'left_consistent_junctions': '5′ 방향 일치 경계',
    'right_consistent_junctions': '3′ 방향 일치 경계', 'opposite_junctions': '반대 방향 경계',
    'direction_consistency': '방향 일치율', 'data_adequacy': '미검출 해석을 위한 데이터 평가',
    'records': '전체 reads', 'bases': '전체 bases', 'next_action': '다음 확인 단계',
    'sra_search': '독립 자료 검색', 'details': '상세 결과',
    'previous_call': '이전 판정', 'previous_canonical_key': '이전 내부 motif 표기',
    'previous_reference_anchor': '이전 비교 기준', 'previous_motif_display': '이전 관측 motif',
    'updated_result': '새 판정', 'updated_motif_display': '새 관측 motif',
    'updated_canonical_key': '새 내부 motif 표기', 'motif_family_changed': 'Motif family 변경',
    'interpretation_changed': '해석 체계 변경', 'reason': '변경 설명',
    'previous_result': '이전 파일', 'updated_result_file': '새 파일',
    'library_context': '시료·라이브러리 해석', 'why': '판정 이유', 'confirmation_missing': '추가 확인 필요',
    'geometry_nominations': '기존 경계 후보 수', 'provisional_candidates': '추가 근거 잠정 후보 수',
    'distinct_junction_sequences': '서로 다른 경계 read 서열', 'physical_flank_clusters': '20-bp flank 군집',
    'deduplicated_orientation': '서열 중복 보정 방향성', 'clean_two_flank_records': '양쪽 clean flank reads',
    'previous_classification': '이전 assessed 판정', 'updated_classification': '새 assessed 판정',
    'previous_version': '이전 버전', 'updated_version': '새 버전', 'classification_changed': '판정 변경',
    'previous_motif': '이전 관측 motif', 'updated_motif': '새 관측 motif',
    'paired_consistency': 'Paired-end 일관성 (보조 관측)',
    'trm_interpretation': 'TRM 해석 · 계통 기대와 관측의 관계',
    'reference_agreement_status': '계통 비교 상태', 'reference_conflict_review': '기대·관측 충돌 검토 필요',
    'main_trm_status': 'Main TRM 확정 상태', 'reference_review_recommendation': '계통 비교 후 우선 확인',
    'previous_reference_agreement': '이전 계통 비교', 'updated_reference_agreement': '새 계통 비교',
    'reference_interpretation_added': '계통 비교 새로 제공',
    'reference_interpretation_changed': '계통 비교 변경',
}
DISPLAY_NAMES.update({'long_read_context_status': '긴 read 추가 근거 상태',
                      'long_read_context_source_run': '긴 read 근거 run',
                      **{'long_read_context_' + key: value
                         for key, value in long_read_context_report.LABELS.items()}})
REFERENCE_AGREEMENT_LABELS = {
    'analysis_incomplete': '입력 분석 미완료 · TRM 해석 보류',
    'qc_limited': 'QC 제약 · TRM 해석 보류',
    'taxonomic_scope_unresolved': '생물학적 출처·분류 범위 확인 필요',
    'library_context_limited': 'Library 제약 · genomic TRM 해석 보류',
    'reference_out_of_scope': '문헌상 비정형 구조 · 단일 short TRM 해석 범위 밖',
    'no_supported_reference': '확립된 계통 기대 미지정 · 관측 후보의 외부 검증 필요',
    'measurement_disagreement': '관측 경로 간 계수 불일치 · 재검토 필요',
    'observed_family_conflict': '계통 기대와 대표 관측 반복 불일치 · main TRM 해석 보류',
    'broad_expectation_difference': '상위 분류군 기대와 차이 · 종별 예외·출처 검토 필요',
    'expected_reference_not_observed': '기대 family 표적 미검출 · 원인 확인 필요',
    'expected_reference_trace': '기대 family 소수 관측 · main TRM 해석 보류',
    'representative_reference_agreement': '계통 기대와 대표 관측 family 일치 · 말단 위치는 별도 검증',
    'reference_observed_no_representative': '기대 family 관측 · 대표 후보 및 말단 근거 확인 필요',
    'expected_reference_unmeasured': '기대 family 미보고 · 명시적 0으로 해석 불가',
}


def reference_review_text(agreement):
    status = agreement['status']
    limited = {
        'analysis_incomplete': '입력·처리 오류를 해결하고 전체 분석을 완료하세요. 남은 관측이나 표적 진단으로 분석 실패를 덮지 않습니다.',
        'qc_limited': 'QC 불일치의 원인을 해결한 뒤 TRM을 해석하세요. 관측 반복은 보존하지만 기술적 실패를 미검출로 해석하지 않습니다.',
        'taxonomic_scope_unresolved': 'Host·혼합 시료·등록 분류군과 관측 서열의 출처를 확인하세요. 반복 관측만으로 특정 종의 TRM을 지정하지 않습니다.',
        'library_context_limited': 'Library 종류·표기 충돌·genomic 대표성의 제약을 확인하세요. 알려진 반복이 관측돼도 genomic main TRM·신규성 해석을 완료한 것으로 처리하지 않습니다.',
        'reference_out_of_scope': '문헌상 비정형·긴·혼합 말단 구조를 별도로 확인하세요. 현재 짧은 반복의 관측만으로 그 구조를 대체하거나 main TRM을 지정하지 않습니다.',
    }
    if status in limited:
        return limited[status]
    if status == 'observed_family_conflict':
        return ('기대 family와 대표 반복을 모두 보존했습니다. 반복량이나 pair 방향만으로 main TRM을 선택하지 말고, '
                'library 간 재현성·분류학적 출처·long-read/assembly 말단 위치를 비교하세요.')
    if status == 'broad_expectation_difference':
        return '상위 분류군 motif에는 종별 예외가 있습니다. 관측 반복을 제외하거나 기대값으로 바꾸지 말고 더 가까운 계통·말단 근거를 확인하세요.'
    if status in {'expected_reference_not_observed', 'expected_reference_unmeasured', 'expected_reference_trace'}:
        return ('미보고·표적 0건·소수 관측을 구분하세요. 현재 입력의 분류학적 정체성·library 대표성·검출 조건을 확인하고 '
                '독립 genomic WGS와 말단 자료로 재검증하세요. 미등록 후보는 보존하되 신규 TRM·부재를 확정하지 않습니다.')
    if status == 'measurement_disagreement':
        return '원래 관측과 표적 재측정의 모든 계수를 보존했습니다. 입력·검출 설정·탐지 경로를 비교해 불일치 원인을 확인하세요.'
    return ''
ADEQUACY_LABELS = {'observed_qc_failure': '입력/처리/QC 오류 확인',
    'insufficient_data': '데이터 또는 조건부 검출력 부족',
    'unassessable_library_or_genome': 'Genome/library 정보 부족 — 실패 배제 불가',
    'conditional_power_supported': '명시한 가정 아래 검출력 충분'}
EVIDENCE_LABELS = {'known_trm_pattern_supported': '반복 reads + 양쪽 경계 방향성',
    'known_repeat_detected_boundary_unresolved': '알려진 반복 검출·경계 근거 미확정',
    'novel_geometry_nomination': '미등록 반복·경계 후보 (추가 검증 근거 부족)',
    'novel_provisional_candidate': '미등록 반복 잠정 후보 (서열 중복·flank 근거, 말단 미확인)',
    'uncatalogued_repeat': '미등록 반복 관측·경계 후보 기준 미충족',
    'published_candidate_pattern_supported': '문헌 후보 + 양쪽 경계 방향성',
    'no_supported_call': '지지되는 판정 없음',
    'weak_or_conflicting_known_repeat': '알려진 반복·지지 부족/방향 충돌',
    'published_candidate_unresolved': '문헌 후보·경계 지지 부족',
    'unresolved_repeat': '일반 반복·TRM 후보 기준 미충족'}
ARCHITECTURE_LABELS = {'exact_short_unit': '단순 짧은 반복', 'long_complex_family': '긴·복합 반복',
    'mixed_short_units': '여러 단위의 혼합 반복', 'diverse': '다양함·단일 motif 미지정',
    'irregular': '불규칙 반복', 'retrotransposon': '레트로트랜스포존 기반'}
CHECK_LABELS = {
    'complete_v050_distinct_sequence_and_flank_measurement': 'v0.5 경계 서열·flank 측정 필요',
    'minimum_distinct_junction_sequences': '서로 다른 경계 read 서열 10개 이상 필요',
    'deduplicated_orientation_consistent': '서열 중복 보정 후 방향성 확인 필요',
    'minimum_error_tolerant_flank_clusters': '완전한 flank 표본에서 2개 이상 서열 군집 확인 필요',
    'no_clean_two_flank_counterevidence': '양쪽 clean flank가 있는 반복의 위치 확인 필요',
    'chromosome_terminal_localization_and_orientation': '염색체 말단 위치·방향의 외부 검증 필요',
    'independent_biological_sample_or_library_validation': '독립 생물학적 시료·라이브러리 검증 필요',
    'literature_review_beyond_this_catalogue': '현재 catalogue 외 문헌의 novelty 검토 필요',
}
OBSERVATION_LABELS = {
    'sequence_observed': '반복 서열 관측 (말단 패턴 판정과 별개)',
    'trace_sequence_observed': '소수 지지 서열 관측',
    'no_sequence_observed': '보고된 지지 서열 없음',
    'measurement_disagreement': '탐지 경로 사이 관측 불일치',
    'not_reported': '기본 탐색에서 미보고 (명시적 0이 아님)',
}


def directional_cluster_columns(readiness):
    diagnostic = readiness.get('directional_flank_clusters') or {}
    status = diagnostic.get('status', 'unmeasured')
    label = {'measured': '측정됨', 'unmeasured': '미측정', 'incomplete': '자료 불완전·미확정',
             'strand_ambiguous': '방향 구분 불가'}.get(status, '미확정')
    result = dict.fromkeys(DIRECTIONAL_CLUSTER_FIELDS, label)
    result['directional_cluster_measurement'] = label
    if status != 'measured':
        return result
    saturated = diagnostic.get('flank_sketch_saturated')
    for column, side in [('consistent_flank_clusters', 'consistent'), ('opposed_flank_clusters', 'opposed')]:
        exact = diagnostic.get(f'all_observed_clusters_with_{side}_events_exact')
        retained = diagnostic.get(f'retained_clusters_with_{side}_events')
        result[column] = (f'{exact} (관측 전체)' if exact is not None else
                          f'보관 {retained} (포화·전체 미확정)' if saturated is True and retained is not None else '미확정')
    ratio = diagnostic.get('largest_retained_consistent_cluster_event_fraction')
    numerator = diagnostic.get('largest_retained_consistent_cluster_events')
    denominator = diagnostic.get('retained_consistent_events')
    result['retained_consistent_cluster_concentration'] = (
        f'보관 {numerator}/{denominator} ({ratio:.1%})' if ratio is not None else
        '보관 일치 사건 0·비율 없음' if denominator == 0 else '미확정')
    retained, observed = diagnostic.get('flank_events_retained'), diagnostic.get('flank_events_observed')
    result['flank_context_completeness'] = (f'{retained}/{observed} 사건 보관 ('+
        ('포화·전체 분포 미확정' if saturated is True else '관측 전체' if saturated is False else '보관 범위 미확정')+')'
        if retained is not None and observed is not None else '미확정')
    result['flank_orientation_basis'] = {
        'published_orientation_in_saved_catalogue': '문헌 방향 (저장 catalogue)',
        'inferred_from_same_junction_events': '같은 경계 관측에서 추정·G-strand 미확립',
        'not_established_no_oriented_events': '방향 사건 없음·방향 미확립',
    }.get(diagnostic.get('orientation_interpretation'), '방향 근거 미확정')
    return result


def validation_columns(candidate):
    readiness = (candidate or {}).get('confirmation_readiness', {})
    exact = readiness.get('distinct_junction_sequences_exact')
    lower = readiness.get('distinct_junction_sequences_lower_bound')
    count = str(exact) if exact is not None else f'>={lower} (상한 미확정)' if lower is not None else '미측정'
    clusters = readiness.get('physical_flank_cluster_count')
    cluster_label = str(clusters) if clusters is not None else ('sketch 포화·전체 군집 수 미확정' if readiness.get('flank_context_sketch_saturated') else '미측정')
    d = readiness.get('deduplicated_orientation', {})
    ratio = d.get('consistent_fraction')
    bound = d.get('consistent_fraction_lower_bound')
    direction = f'{ratio:.1%} (고유 서열, 독립 분자 아님)' if ratio is not None else (
        f'보수적 하한 {bound:.1%} (sketch 포화)' if bound is not None else '미측정')
    missing = list(readiness.get('missing_evidence', []))
    if (candidate or {}).get('reference_status') == 'not_in_catalogue':
        missing.extend(candidate.get('external_confirmation', {}).get('missing_evidence', []))
    return {'confirmation_status': readiness.get('status', 'unmeasured'),
            'confirmation_missing': '; '.join(CHECK_LABELS.get(code, code) for code in missing),
            'distinct_junction_sequences': count, 'physical_flank_clusters': cluster_label,
            'deduplicated_orientation': direction,
            'clean_two_flank_records': readiness.get('clean_two_flank_repeat_records'),
            **directional_cluster_columns(readiness)}


def paired_columns(candidate):
    pairs = (candidate or {}).get('paired_evidence') or {}
    return {'pair_measurement_status': pairs.get('status', 'not_measured'),
            'pair_count_source': pairs.get('count_source'), 'pair_layout': pairs.get('layout'),
            'pair_scope': pairs.get('scope'),
            **{'pair_'+key: pairs.get(key) for key in (*PAIR_COUNT_FIELDS, 'one_mate_only_pairs')}}


def paired_summary(candidate):
    pairs = (candidate or {}).get('paired_evidence') or {}
    if not pairs.get('measured'):
        return '미측정'
    scope = '앞부분 pair만 검사 · ' if pairs.get('scope') == 'synchronized_prefix_only' else ''
    return (scope+f"지지 {pairs['supporting_pairs']:,} · 양쪽 {pairs['both_mates_repeat_pairs']:,} "
            f"(역상보 {pairs['opposite_strand_pairs']:,} / 같은 방향 {pairs['same_strand_pairs']:,} / "
            f"방향 모호 {pairs['ambiguous_strand_pairs']:,}) · 한쪽만 {pairs['one_mate_only_pairs']:,}")


def load_library_evidence(path, metadata):
    """Read a separately reviewed table; validate identities and file provenance."""
    if not path:
        return {}
    accepted = {'RNA_component_supported', 'metadata_conflict_unresolved', 'HiC_supported',
                'HiC_annotation_only', 'WGS_context_compatible', 'taxonomic_identity_conflict'}
    result = index_samples(read_tsv(path), 'reviewed library evidence')
    for sample, row in result.items():
        if sample not in metadata:
            raise ValueError(f'Library evidence sample absent from metadata: {sample}')
        for field in ('scientific_name', 'tax_id'):
            if not row.get(field) or str(row[field]).strip() != str(metadata[sample].get(field, '')).strip():
                raise ValueError(f'{sample}: library evidence {field} differs from metadata')
        if row.get('library_review_code') not in accepted:
            raise ValueError(f'{sample}: unrecognized library evidence review code')
        if row.get('library_review_code') == 'taxonomic_identity_conflict':
            for field in ('source_result_sha256', 'evidence_sha256'):
                value = row.get(field, '')
                if len(value) != 64 or any(c not in '0123456789abcdef' for c in value):
                    raise ValueError(f'{sample}: reviewed identity conflict requires {field}')
        evidence = row.get('evidence_path') or row.get('detailed_evidence_path')
        if not evidence:
            raise ValueError(f'{sample}: library evidence requires a local provenance file')
        source = Path(evidence)
        if not source.is_absolute():
            source = Path(path).resolve().parent/source
        if not source.is_file():
            raise ValueError(f'{sample}: library evidence provenance file missing')
        digest = sha256(source)
        if row.get('evidence_sha256') and row['evidence_sha256'] != digest:
            raise ValueError(f'{sample}: library evidence provenance hash differs')
        row.update(detailed_evidence_path=str(source.resolve()), evidence_sha256=digest,
                   library_evidence_table=str(Path(path).resolve()), library_evidence_table_sha256=sha256(path))
        if row.get('previous_assessment_path') or row.get('previous_assessment_sha256'):
            previous = Path(row.get('previous_assessment_path', ''))
            if not previous.is_absolute():
                previous = Path(path).resolve().parent/previous
            if not previous.is_file() or sha256(previous) != row.get('previous_assessment_sha256'):
                raise ValueError(f'{sample}: library audit previous assessment hash differs')
            doc = json.loads(previous.read_text())
            if doc.get('sample_id') != sample or str(doc.get('tax_id')) != str(metadata[sample].get('tax_id')):
                raise ValueError(f'{sample}: library audit source identity differs')
    return result


def compare_previous_assessment(document, directory, output):
    sample = document['sample_id']
    previous_path = Path(directory)/'samples'/f'{sample}.json'
    if not previous_path.is_file():
        previous_path = Path(directory)/f'{sample}.json'
    previous = json.loads(previous_path.read_text()) if previous_path.is_file() else None
    if previous is not None:
        if previous.get('program') != 'trm-assessment' or previous.get('sample_id') != sample:
            raise ValueError(f'{sample}: --previous-assessment must contain matching assessed JSON, not raw scans')
        if str(previous.get('tax_id')) != str(document.get('tax_id')):
            raise ValueError(f'{sample}: previous assessment taxon differs')
    old = previous or {}
    old_key = (old.get('top_candidate') or {}).get('canonical_family', '')
    new_key = (document.get('top_candidate') or {}).get('canonical_family', '')
    legacy_geometry = sum(c.get('geometry_nomination', c.get('novel_candidate', False)) for c in old.get('candidates', []))
    columns = validation_columns(document.get('top_candidate'))
    return {'sample_id': sample, 'scientific_name': document['scientific_name'],
        'comparison_status': 'compared_assessed_json' if previous is not None else 'previous_assessment_missing',
        'previous_version': old.get('version'), 'updated_version': document['version'],
        'previous_schema_version': old.get('schema_version'), 'updated_schema_version': document['schema_version'],
        'previous_followup_status': old.get('followup_diagnostics', {}).get('status', 'not_reported'),
        'updated_followup_status': document.get('followup_diagnostics', {}).get('status', 'not_reported'),
        'followup_diagnostics_added': (previous is not None
            and old.get('followup_diagnostics', {}).get('status') != 'complete'
            and document.get('followup_diagnostics', {}).get('status') == 'complete'),
        'previous_classification': old.get('classification'), 'updated_classification': document['classification'],
        'previous_result': old.get('result'), 'updated_result': document['result'],
        'previous_motif': old.get('detected_motif'), 'updated_motif': document.get('detected_motif'),
        'previous_canonical_key': old_key, 'updated_canonical_key': new_key,
        'motif_family_changed': old_key != new_key if previous is not None else None,
        'classification_changed': old.get('classification') != document['classification'] if previous is not None else None,
        'previous_reference_agreement': old.get('reference_agreement', {}).get('status'),
        'updated_reference_agreement': document.get('reference_agreement', {}).get('status'),
        'reference_interpretation_added': previous is not None and 'reference_agreement' not in old,
        'reference_interpretation_changed': (
            old['reference_agreement'].get('status') != document.get('reference_agreement', {}).get('status')
            if 'reference_agreement' in old else None),
        'previous_geometry_nominations': legacy_geometry if previous is not None else None,
        'updated_geometry_nominations': sum(c['geometry_nomination'] for c in document['candidates']),
        'updated_provisional_candidates': sum(c['novel_candidate'] for c in document['candidates']),
        'confirmation_missing': columns['confirmation_missing'],
        'reason': '이전 assessed 판정과 새 증거 단계·분류군/라이브러리 조건을 비교합니다. 경계 후보는 보존하며 말단 확인은 별도입니다.' if previous is not None else '이 샘플의 이전 assessed JSON이 없어 비교하지 않았습니다.',
        'previous_assessment_file': str(previous_path.resolve()) if previous is not None else None,
        'previous_assessment_sha256': sha256(previous_path) if previous is not None else None,
        'updated_result_file': str((Path(output)/'samples'/f'{sample}.json').resolve())}


def write_tsv(path, rows, fields):
    import csv
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow({key: '' if row.get(key) is None else row.get(key, '') for key in fields})


def write_html(path, rows, fields, title, explanation, links=()):
    """Portable, searchable table with literal text escaped; no external assets."""
    headers = ''.join('<th>'+html.escape(DISPLAY_NAMES.get(field, field))+'</th>' for field in fields)
    body = []
    for row in rows:
        cells = []
        for key in fields:
            value = '' if row.get(key) is None else str(row.get(key, ''))
            if key == 'repeat_architecture':
                value = ARCHITECTURE_LABELS.get(value, value)
            if key == 'sra_search' and value.startswith('https://www.ncbi.nlm.nih.gov/sra/'):
                cell = '<a href="'+html.escape(value, quote=True)+'">SRA 검색</a>'
            elif key == 'details':
                cell = '<a href="'+html.escape(value, quote=True)+'">상세 보기</a>'
            else:
                cell = html.escape(value)
            cells.append('<td>'+cell+'</td>')
        body.append('<tr>'+''.join(cells)+'</tr>')
    navigation = ' · '.join('<a href="'+html.escape(target, quote=True)+'">'+html.escape(label)+'</a>' for label, target in links)
    text = ('<!doctype html><html lang="ko"><meta charset="utf-8"><title>'+html.escape(title)+
            '</title><style>body{font:15px system-ui;margin:24px;color:#182b3a}h1{font-size:24px}'
            'input{padding:9px;width:360px;max-width:80vw;margin:12px 0}table{border-collapse:collapse}'
            'th,td{padding:8px 10px;border-bottom:1px solid #ddd;vertical-align:top;text-align:left}'
            'th{position:sticky;top:0;background:#eaf1f5;cursor:pointer;white-space:nowrap}'
            'tbody tr:hover{background:#f3f8fa}p{max-width:1050px;line-height:1.6}a{color:#0865a5}'
            '.scroll{overflow:auto;max-height:75vh}small{color:#526675}</style><h1>'+html.escape(title)+
            '</h1><p>'+html.escape(explanation)+'</p><p>'+navigation+'</p><input id="query" placeholder="종 이름, accession, 결과로 검색">'
            '<small id="count"></small><div class="scroll"><table><thead><tr>'+headers+
            '</tr></thead><tbody>'+''.join(body)+'</tbody></table></div><script>'
            'const input=document.getElementById("query"),count=document.getElementById("count"),'
            'tbody=document.querySelector("tbody");function filter(){let n=0;for(const r of tbody.rows){'
            'r.hidden=!r.textContent.toLowerCase().includes(input.value.toLowerCase());if(!r.hidden)n++;}'
            'count.textContent=n+" samples"}input.oninput=filter;filter();'
            'document.querySelectorAll("th").forEach((h,i)=>{let asc=true;h.onclick=()=>{'
            'const rows=[...tbody.rows];rows.sort((a,b)=>{const x=a.cells[i].textContent,y=b.cells[i].textContent;'
            'return (asc?1:-1)*(x!==""&&y!==""&&Number.isFinite(+x)&&Number.isFinite(+y)?+x-+y:x.localeCompare(y))});'
            'rows.forEach(r=>tbody.appendChild(r));asc=!asc}});</script></html>')
    Path(path).write_text(text)


def followup_html(diagnostic):
    """Render supplemental diagnostics without changing the decision narrative."""
    if not diagnostic:
        return []
    esc = lambda value: html.escape(str(value), quote=True)
    content = ['<h2>선택 구간의 후속 구조 진단</h2>',
               '<p><b>'+esc(diagnostic['display'])+'</b></p>',
               '<p>'+esc(diagnostic['summary'])+'</p>']
    if diagnostic['status'] in ('not_requested', 'not_available'):
        return content
    if diagnostic['status'] == 'failed':
        reason = (diagnostic.get('provenance') or {}).get('failure_reason')
        if reason:
            content.append('<p>기록된 실패 사유: '+esc(reason)+'</p>')
    if diagnostic.get('scope'):
        content.append('<p>적용 범위: '+esc(diagnostic['scope'])+'</p>')
    if diagnostic.get('rows'):
        content.append('<table><tr><th>측정 항목</th><th>값</th><th>단위</th><th>분모·적용 범위</th></tr>')
        for row in diagnostic['rows']:
            value = row['value']
            display = '미측정' if value is None else f'{value:,}' if type(value) is int else str(value)
            content.append('<tr><td>'+esc(row['label'])+'</td><td>'+esc(display)+
                           '</td><td>'+esc(row['unit'])+'</td><td>'+esc(row['scope'])+'</td></tr>')
        content.append('</table>')
    if diagnostic.get('recommendations'):
        content.append('<p>이 진단에 따른 확인 권고:</p><ul>')
        content.extend('<li>'+esc(value)+'</li>' for value in diagnostic['recommendations'])
        content.append('</ul>')
    content.append('<p class="note">선택된 구간의 구조 설명입니다. 전체 표본의 발생률·독립 분자 수·'
                   '검출력을 추정하지 않으며, 기존 TRM 판정·대표 motif·데이터 충분성 평가를 바꾸지 않습니다. '
                   '진단 미제공·진단 실행 실패는 반복 미검출이나 원래 분석의 실패와 구분합니다. '
                   '출처와 검사 범위는 <a href="../followup_diagnostics.tsv">진단 TSV</a> 및 이 표본의 JSON에 보존합니다.</p>')
    return content


def _family_count_text(measurement):
    measurement = measurement or {}
    reads = measurement.get('reads')
    if reads is not None:
        return (f'{reads:,} (표적 0건 검증 미확인)'
                if measurement.get('status') == 'reported_zero_unverified' else f'{reads:,}')
    counts = measurement.get('read_counts', [])
    if counts:
        return '상이한 계수: '+', '.join(str(value) for value in counts)
    return '미보고/미측정'


def family_observation_tables(sample, scientific_name, report):
    """Full reference/neighbor tables; never sum counts across observations."""
    families, neighbors = [], []
    for family in report.get('families', []):
        current, original = family.get('exact_current', {}), family.get('exact_original', {})
        families.append({
            'sample_id': sample, 'scientific_name': scientific_name,
            'expected_family': family['canonical_family'],
            'expected_motifs': ','.join(str(m) for m in family.get('expected_motifs', []) if m is not None),
            'candidate_catalog_statuses': json.dumps(family.get('candidate_catalog_statuses', [])),
            **{key: family.get(key) for key in ('reference_status', 'reference_scope',
                'status', 'display_label', 'exact_measurement_status', 'measurement_disagreement',
                'nearby_family_count', 'nearby_truncated', 'other_positive_family_count',
                'outside_de_novo_search_range', 'measurement_context_status',
                'nearby_comparison_status')},
            'current_reads': current.get('reads'), 'original_reads': original.get('reads'),
            'current_status': current.get('status'), 'original_status': original.get('status'),
            'current_read_counts': json.dumps(current.get('read_counts', [])),
            'original_read_counts': json.dumps(original.get('read_counts', [])),
            'nearby_displayed_count': len(family.get('nearby_families', [])),
        })
        for neighbor in family.get('nearby_families_all', []):
            current, original = neighbor.get('current', {}), neighbor.get('original', {})
            neighbors.append({
                'sample_id': sample, 'scientific_name': scientific_name,
                'expected_family': family['canonical_family'],
                'observed_family': neighbor['canonical_family'],
                **{key: neighbor.get(key) for key in ('edit_distance', 'operation', 'measurement_disagreement')},
                'current_reads': current.get('reads'), 'original_reads': original.get('reads'),
                'current_read_counts': json.dumps(current.get('read_counts', [])),
                'original_read_counts': json.dumps(original.get('read_counts', [])),
                **{key: json.dumps(neighbor.get(key), ensure_ascii=False, sort_keys=True)
                   for key in ('sources', 'buckets', 'edit_witness')},
            })
    return families, neighbors


def family_observation_html(report):
    if not report:
        return []
    esc = lambda value: html.escape(str(value), quote=True)
    contexts = {'source_unavailable': '원본 결과 없음', 'partial_input': '일부 입력 분석',
                'analysis_incomplete': '분석 미완료', 'complete_input': '전체 입력 결과',
                'input_completeness_unverified': '입력 전체 분석 여부 미확인'}
    measurements = {'reported_positive': '양성 관측', 'explicit_zero': '표적 측정 0건',
                    'not_reported': '미보고 · 0건 아님',
                    'measurement_disagreement': '측정 간 불일치',
                    'reported_zero_unverified': '저장값 0 · 표적 측정 출처 미확인',
                    'invalid_counts': '유효하지 않은 계수'}
    context = report.get('measurement_context_status', 'input_completeness_unverified')
    content = ['<h2>기대 family와 근접 반복 관측</h2>',
               '<p><b>'+esc(report['display_label'])+'</b> · '+
               esc(contexts.get(context, context))+'</p>']
    if report.get('families'):
        content.append('<table><tr><th>기대 family</th><th>관측 상태</th><th>정확 family 측정 상태</th>'
                       '<th>현재 reads</th><th>원래 reads</th><th>한 염기 차이 family 수</th></tr>')
        for family in report['families']:
            content.append('<tr><td>'+esc(family['canonical_family'])+'</td><td>'+esc(family['display_label'])+
                '</td><td>'+esc(measurements.get(family['exact_measurement_status'], family['exact_measurement_status']))+'</td><td>'+
                esc(_family_count_text(family.get('exact_current')))+'</td><td>'+
                esc(_family_count_text(family.get('exact_original')))+'</td><td>'+
                esc(family['nearby_family_count'])+'</td></tr>')
        content.append('</table>')
        for family in report['families']:
            if not family.get('nearby_families'):
                continue
            content.append('<p>기대 '+esc(family['canonical_family'])+'의 근접 family: 전체 '+
                           esc(family['nearby_family_count'])+'개 중 '+
                           esc(len(family['nearby_families']))+'개 표시</p><ul>')
            for neighbor in family['nearby_families']:
                content.append('<li>'+esc(neighbor['canonical_family'])+' · '+esc(neighbor['operation'])+
                    ' '+esc(neighbor['edit_distance'])+'회 · 현재 '+
                    esc(_family_count_text(neighbor.get('current')))+' / 원래 '+
                    esc(_family_count_text(neighbor.get('original')))+' reads'+
                    (' · 측정 불일치' if neighbor.get('measurement_disagreement') else '')+'</li>')
            content.append('</ul>')
    content.append('<p class="note">전체 저장 primary·rejected 후보와 현재 관측을 비교합니다. '
        'Primitive 단위의 회전·역상보를 같은 family로 묶고, 한 염기 치환·삽입·결실을 근접 관측으로 표시합니다. '
        '근접 비교 범위는 primitive 3–64 bp이며, 범위 밖 단위의 정확 관측은 보존합니다. '
        '근접성은 정확 일치나 TRM 신뢰도 점수가 아니며 canonical·novel 판정을 바꾸지 않습니다. '
        '미보고는 측정된 0이나 생물학적 부재가 아닙니다. 부분 입력·분석 실패·library 제약은 그대로 적용됩니다. '
        '<a href="../expected_family_observations.tsv">기대 family 전체 표</a> · '
        '<a href="../nearby_families.tsv">근접 family 전체 표</a>. '
        '출처별 계수·측정 불일치·실제 탐색 범위는 JSON에 보존합니다.</p>')
    return content


def input_applicability_rows(document):
    """Export separate diagnostic rows without editing observed family tables."""
    report = document.get('input_applicability') or {}
    result = []
    for kind, values in (('configured_period', report.get('periods', [])),
                         ('expected_family', report.get('families', []))):
        for value in values:
            result.append({'sample_id': document['sample_id'],
                'scientific_name': document['scientific_name'], 'kind': kind,
                'scope': report.get('scope'), 'input_status': report.get('status'),
                **{key: value.get(key) for key in INPUT_APPLICABILITY_FIELDS
                   if key not in ('sample_id', 'scientific_name', 'kind', 'scope', 'input_status')},
                'expected_motifs': ','.join(value.get('expected_motifs', []))})
    # Missing or invalid QC must remain visible in the per-input export.
    if not result:
        result.append({'sample_id': document['sample_id'], 'scientific_name': document['scientific_name'],
                       'kind': 'input', 'scope': report.get('scope'),
                       'input_status': report.get('status'), 'status': report.get('status'),
                       'display_label': report.get('display_label')})
    return result


def input_applicability_html(report):
    if not report:
        return []
    esc = lambda value: html.escape(str(value), quote=True)
    scopes = {'complete_input': '완료된 전체 입력', 'analyzed_records': '분석된 records만',
              'unverified_input': '입력 전체 여부 미확인', 'unavailable': '평가 범위 미확인'}
    content = ['<h2>입력 길이와 탐색 조건</h2>',
               '<p><b>'+esc(report['display_label'])+'</b> · '+
               esc(scopes.get(report['scope'], report['scope']))+'</p>',
               '<p>'+esc(report.get('recommendation', ''))+'</p>']
    if report.get('families'):
        content.append('<table><tr><th>기대 motif (회전·RC 통합)</th><th>단위 bp</th>'
                       '<th>최소 seed bp</th><th>길이 적격 reads 하한–상한</th><th>길이 진단</th></tr>')
        for family in report['families']:
            bounds = ('미평가' if family.get('length_eligible_reads_lower') is None else
                      f'{family["length_eligible_reads_lower"]:,}–{family["length_eligible_reads_upper"]:,}')
            content.append('<tr><td>'+esc(', '.join(family.get('expected_motifs', [])))+
                '</td><td>'+esc(family.get('period_bp', '미평가'))+'</td><td>'+esc(
                    family['minimum_seed_bp'] if family.get('minimum_seed_bp') is not None else '미평가')+
                '</td><td>'+esc(bounds)+'</td><td>'+esc(family.get('display_label'))+
                (' · de novo 탐색 주기 범위 밖' if family.get('outside_discovery_range') else '')+'</td></tr>')
        content.append('</table>')
    content.append('<p class="note">최소 seed 길이 = k × max(min_copies, ceil(seed_bases/k)). '
        '현재 입력의 길이와 원래 discovery 설정만 비교합니다. Paired mates의 길이는 합치지 않습니다. '
        '길이 조건 충족은 텔로미어의 충분한 표집·서열 품질·library 보존·검출력을 보장하지 않습니다. '
        '범위 밖 기대 family도 별도 표적 검사로 관측될 수 있습니다. '
        '기존 관측 0·미보고·후보 순위·판정은 바꾸지 않으며 부재·ALT·novel의 근거로 사용하지 않습니다. '
        '<a href="../input_applicability.tsv">탐색 주기·기대 family별 길이 진단 TSV</a></p>')
    return content


def flank_applicability_rows(document):
    """Keep necessary-length bounds distinct from measured flank event counts."""
    report = (document.get('input_applicability') or {}).get('flank_period_applicability') or {}
    common = {'sample_id': document['sample_id'], 'scientific_name': document['scientific_name'],
              'scope': report.get('scope'), 'input_status': report.get('status'),
              'annotation_enabled': report.get('annotation_enabled'),
              'model_status': report.get('model_status')}
    result = []
    for kind, values in (('configured_period', report.get('periods', [])),
                         ('expected_family', report.get('families', []))):
        for value in values:
            result.append({**common, 'kind': kind,
                **{key: value.get(key) for key in FLANK_APPLICABILITY_FIELDS
                   if key not in (*common, 'kind')},
                'expected_motifs': ','.join(value.get('expected_motifs', []))})
    if not result:
        result.append({**common, 'kind': 'input', 'status': report.get('status'),
                       'display_label': report.get('display_label')})
    return result


def flank_applicability_html(report):
    if not report:
        return []
    esc = lambda value: html.escape(str(value), quote=True)
    show = lambda value: esc('미평가' if value is None else value)
    scopes = {'complete_input': '완료된 전체 입력', 'analyzed_records': '분석된 records만',
              'unverified_input': '입력 전체 여부 미확인', 'unavailable': '평가 범위 미확인'}
    enabled = report.get('annotation_enabled')
    enabled_label = '설정됨' if enabled is True else '설정되지 않음' if enabled is False else '미확인'
    content = ['<h2>Flank 주기 평가의 길이 필요조건</h2>',
               '<p><b>'+esc(report.get('display_label'))+'</b> · '+
               esc(scopes.get(report.get('scope'), report.get('scope')))+'</p>',
               '<p>원래 scan의 flank 주기 annotation: '+esc(enabled_label)+
               ' · 길이 모델: '+show(report.get('model_status'))+'</p>',
               '<p>'+esc(report.get('recommendation', ''))+'</p>']
    for title, values, family in (
            ('기대 family별 필요조건', report.get('families', []), True),
            ('설정된 탐색 주기별 필요조건', report.get('periods', []), False)):
        if not values:
            continue
        content.append('<h3>'+title+'</h3><table><tr>'+(
            '<th>기대 motif (회전·RC 통합)</th>' if family else '')+
            '<th>단위 bp</th><th>최소 seed bp</th><th>주기 평가 창 bp</th>'
            '<th>Junction flank bp</th><th>최소 read bp</th>'
            '<th>길이 필요조건 충족 reads 하한–상한</th><th>길이 진단</th></tr>')
        for value in values:
            low, high = value.get('length_eligible_reads_lower'), value.get('length_eligible_reads_upper')
            bounds = '미평가' if low is None or high is None else f'{low:,}–{high:,}'
            content.append('<tr>'+('<td>'+esc(', '.join(value.get('expected_motifs', [])))+'</td>'
                if family else '')+''.join('<td>'+show(value.get(key))+'</td>' for key in (
                    'period_bp', 'minimum_seed_bp', 'required_window_bp', 'junction_flank_bp',
                    'minimum_read_bp'))+'<td>'+esc(bounds)+'</td><td>'+show(value.get('display_label'))+
                (' · de novo 탐색 주기 범위 밖' if value.get('outside_discovery_range') else '')+'</td></tr>')
        content.append('</table>')
    content.append('<p class="note">최소 seed d = k × max(min_copies, ceil(seed_bases/k)), '
        '주기 평가 창 W = max(20, 3k), junction flank 설정 F일 때 최소 read 길이는 '
        'max(d + max(W, F), floor(10d/9) + 1)입니다. '
        '실제 tract가 더 길거나 flank에 모호한 염기·다른 같은-family tract가 있으면 추가 길이가 필요할 수 있습니다. '
        '이 수치는 원래 설정과 read 길이의 필요조건이며 실제 junction·flank 관측 계수가 아닙니다. '
        'Annotation이 설정되지 않았거나 실제 사건이 미보고이면 측정된 0으로 해석하지 않습니다. '
        'Novel 후보의 primitive 단위 길이는 설정된 주기 표에서 확인하며, 범위 밖 후보는 그 표로 평가하지 않습니다. '
        'Paired mates 길이는 합치지 않습니다. 조건 충족도 clean flank·말단 위치·독립 분자·검출력을 보장하지 않으며 '
        '후보 계수·순위·판정·adequacy를 바꾸지 않습니다. '
        '<a href="../flank_period_applicability.tsv">주기·기대 family별 flank 길이 진단 TSV</a></p>')
    return content


def technical_boundary_rows(document):
    diagnostic = document['technical_boundary_diagnostics']
    base = {key: diagnostic.get(key) for key in (
        'status', 'scope', 'selected_event_count', 'selected_pair_count',
        'full_read_count', 'flagged_event_count')}
    base.update(sample_id=document['sample_id'], scientific_name=document.get('scientific_name'))
    for tier in diagnostic.get('tiers') or [{}]:
        yield {**base, **tier}


def technical_event_rows(document):
    for event in document['technical_boundary_diagnostics'].get('events', []):
        for threshold, overlap in event.get('union_flank_overlap_by_threshold', {}).items():
            yield {**event, 'sample_id': document['sample_id'], 'threshold_bp': int(threshold),
                   'union_flank_overlap_bp': overlap,
                   'adjacent_unmatched_flank_bp': event['adjacent_unmatched_flank_by_threshold'][threshold]}


def technical_boundary_html(diagnostic):
    if not diagnostic:
        return []
    esc = lambda value: html.escape(str(value), quote=True)
    content = ['<h2>기술 서열과 반복 경계의 관계</h2>',
        '<p><b>'+esc(diagnostic['display_label'])+'</b></p>',
        '<p>'+esc(diagnostic['summary'])+'</p>',
        '<p>'+esc(diagnostic['recommendation'])+'</p>']
    if diagnostic.get('tiers'):
        content.append('<table><tr><th>Exact 일치 최소 길이</th><th>일치한 reads</th>'
                       '<th>일치한 pairs</th><th>Flank와 겹친 사건</th></tr>')
        for tier in diagnostic['tiers']:
            content.append('<tr>'+''.join('<td>'+esc(tier[key])+'</td>' for key in (
                'threshold_bp', 'full_reads_with_matches', 'pairs_with_matches',
                'events_with_flank_overlap'))+'</tr>')
        content.append('</table>')
    content.append('<p class="note">선정된 경계의 제한된 exact-match 패널 검사입니다. '
        '기준별 계수는 서로 겹치므로 합산하지 않습니다. 미일치는 adapter-free·유전체 기원·'
        '고유 subtelomere를 증명하지 않으며, 경고도 motif 자동 제외나 kit 확정이 아닙니다. '
        'Native 관측·판정·미검출 검출력은 유지합니다. '
        '<a href="../technical_boundary_diagnostics.tsv">기준별 요약 TSV</a> · '
        '<a href="../technical_boundary_events.tsv">사건별 경계 TSV</a></p>')
    return content


def reference_localization_rows(document, key=None):
    """Serialize reviewed rows only; never recalculate alignments or decisions."""
    diagnostic = document['reference_localization_diagnostics']
    rows = [diagnostic] if key is None else diagnostic.get(key, [])
    for row in rows:
        values = {**row, 'sample_id': document['sample_id'],
                  'scientific_name': document.get('scientific_name')}
        yield {name: json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False)
               if isinstance(value, (list, dict)) else value for name, value in values.items()}


def cross_period_rows(document, key=None):
    """Export sealed saved-interval diagnostics without rescoring native results."""
    diagnostic = document['cross_period_evidence']
    rows = [diagnostic] if key is None else diagnostic.get(key, [])
    for row in rows:
        values = {**row, 'sample_id': document['sample_id'],
                  'scientific_name': document.get('scientific_name')}
        yield {name: json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False)
               if isinstance(value, (list, dict)) else value for name, value in values.items()}


def assembly_context_rows(document, key=None):
    """Export saved annotations, with absent measurements left as empty cells."""
    evidence = document['assembly_context_evidence']
    rows = ([dict(status=evidence['status'], scope=evidence['scope'], **evidence['summary'])]
            if key is None else evidence.get(key, []))
    for row in rows:
        values = {**row, 'sample_id': document['sample_id'],
                  'scientific_name': document.get('scientific_name')}
        yield {name: json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False)
               if isinstance(value, (list, dict)) else value for name, value in values.items()}


def quality_path_rows(document, key=None):
    """Export the separate saved-quality diagnostic; missing is null, not zero."""
    evidence = document['quality_path_evidence']
    rows = ([dict(status=evidence['status'], scope=evidence['scope'], **evidence['summary'])]
            if key is None else evidence.get(key, []))
    for row in rows:
        values = {**row, 'sample_id': document['sample_id'],
                  'scientific_name': document.get('scientific_name')}
        yield {name: json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False)
               if isinstance(value, (list, dict)) else value for name, value in values.items()}


def quality_path_html(evidence):
    if not evidence:
        return []
    def esc(value):
        return html.escape('미평가' if value is None else str(value), quote=True)
    content = ['<h2>같은 read의 대체 품질과 최적경로 제약</h2>']
    if evidence['status'] == 'not_provided':
        content.append('<p>선택적 품질·경로 근거 미제공. 추가 측정은 미평가이며 관측 0건이나 TRM 부재를 뜻하지 않습니다.</p>')
    else:
        content.append('<table><tr><th>검토 범위</th><th>저장된 근거</th></tr>')
        for key in quality_path_evidence.FIELDS[1:]:
            content.append('<tr><th>' + esc(DISPLAY_NAMES[key]) + '</th><td>' +
                           esc(evidence['display'][key]) + '</td></tr>')
        content.append('</table><p>아래 수치는 미리 선정한 집중 혼합 구간에 한정합니다. '
            '고정 위치 대비는 한 모델의 모든 최적 경로에서 match(M), 다른 모델의 모든 최적 경로에서 '
            '관측 염기 edit(E)인 경우입니다. 경로마다 match/edit가 달라지는 A를 확정 위치로 바꾸지 않습니다. '
            '고정 위치 대비 0건과 경로 전체에서 높은 Q의 관측 edit가 필요하다는 제약은 서로 다릅니다.</p>')
        columns = (('motif', '주기 모델'), ('focused_mixed_intervals', '집중 혼합 구간'),
            ('intervals_requiring_Q_ge30_observed_edit', '모든 최적 경로에 Q≥30 관측 edit 필요'),
            ('intervals_permitting_zero_observed_edit', '관측 edit 0인 최적 경로 있음'),
            ('intervals_with_exclusively_gap_edits', '모든 최적 경로가 gap만으로 설명'),
            ('deletion_required_without_fixed_mandatory_boundary', 'Deletion 필요·고정 boundary 없음'),
            ('strict_position_contrast_intervals', '엄격 고정 위치 대비 구간'))
        content.append('<table><tr>' + ''.join('<th>'+esc(label)+'</th>' for _,label in columns) + '</tr>')
        for row in evidence['models']:
            content.append('<tr>' + ''.join('<td>'+esc(row.get(key))+'</td>' for key,_ in columns) + '</tr>')
        content.append('</table><p class="note">각 칸은 구간 수이며 독립 분자 수가 아닙니다. '
            '관측 edit는 모델상의 insertion/substitution이며 실제 sequencing error로 확정한 것이 아닙니다. '
            '관측 염기가 없는 deletion의 gap Q는 null입니다. 인접 염기의 Q를 gap 신뢰도로 대체하지 않습니다. '
            '품질로 최적 경로나 비용을 다시 정하지 않았고, 높은 Q만으로 variant 정확도·물리적 말단을 확정하지 않습니다. '
            '11bp 고정 교대 모델은 임의의 5/6bp 혼합 문법 전체와 다릅니다.</p>')
        content.append('<p>' + ' '.join(esc(note) for note in evidence['limitations']) + '</p>')
    content.append('<p>' + ' · '.join('<a href="../'+name+'">'+label+'</a>' for name,label in
        zip(QUALITY_PATH_TABLES, ('품질·분모 요약 TSV', '집중 모델별 TSV', '같은 read 품질 비교 TSV'))) + '</p>')
    return content


def assembly_context_html(evidence, quality_path=None):
    if not evidence:
        return []
    def esc(value):
        return html.escape('미평가' if value is None else str(value), quote=True)
    content = ['<h2>Read 관측 → assembly 확장 → 물리적 말단</h2>']
    if quality_path and quality_path['status'] == 'sources_validated':
        content.append('<p class="note">이 assembly 문맥의 품질 접근·비교 수치는 이전 검토 시점의 기록입니다. '
            '당시 대체 품질 비교 0건을 현재 미확보로 해석하지 마세요. 후속 same-read 비교와 최적경로 품질은 '
            '아래 별도 항목에 연결했습니다. 원 source 계수와 provenance는 보존합니다.</p>')
    if evidence['status'] == 'not_provided':
        content.append('<p>저장 assembly 문맥 근거 미제공. 관측 0건이나 TRM 부재를 뜻하지 않습니다.</p>')
    else:
        content.append('<table><tr><th>근거 층위</th><th>검토 결과</th></tr>')
        for key in assembly_context_evidence.FIELDS[1:]:
            content.append('<tr><th>' + esc(DISPLAY_NAMES[key]) + '</th><td>' +
                           esc(evidence['display'][key]) + '</td></tr>')
        content.append('</table><p>검토한 구조는 TTAGG·TTAGGG 혼합 배열입니다. '
            '정책상 대표 반복과 별개이며 B. rubropictus의 대표 후보를 대체하지 않습니다. '
            '선정 read 구간, assembly 구간, 사건, HSP, pair는 다른 분모입니다. '
            '같은 표본의 short-read assembly이므로 독립 생물학적 검증이 아닙니다.</p>')
        content.append('<table><tr><th>Assembly</th><th>Scaffold</th><th>0-based 시작</th>'
            '<th>끝(exclusive)</th><th>혼합 문법 길이(bp)</th><th>방향쪽 scaffold 끝 도달</th></tr>')
        for row in evidence['intervals']:
            content.append('<tr>' + ''.join('<td>' + esc(row.get(key)) + '</td>' for key in
                ('assembly_accession', 'target', 'start0', 'end0', 'length_bp', 'outward_scaffold_edge')) + '</tr>')
        content.append('</table><p class="note">Scaffold 끝 도달은 물리적 염색체 말단 확인이 아닙니다. '
            'FCS 조치 행 0은 숙주·핵 유래의 확정이 아니며, 고정 Q30의 발생 단계와 오류율은 미확정입니다. '
            '원 FASTQ 접근 실패와 binary header 접근은 같은 read의 quality 비교가 아닙니다. '
            '기존 후보·계수·판정·미검출 검출력은 변경하지 않습니다.</p>')
    content.append('<p>' + ' · '.join('<a href="../' + name + '">' + label + '</a>' for name, label in
        zip(ASSEMBLY_CONTEXT_TABLES, ('층위별 요약 TSV', 'Scaffold 출처 TSV', 'Assembly 구간 TSV', '사건·HSP TSV'))) + '</p>')
    return content


def cross_period_html(diagnostic):
    if not diagnostic:
        return []
    def esc(value):
        return html.escape('미평가' if value is None else str(value), quote=True)
    content = ['<h2>서로 다른 반복 주기와 같은 위치의 근거</h2>',
               '<p><b>' + esc(diagnostic['display_label']) + '</b></p>',
               '<p>' + esc(diagnostic['summary']) + '</p>',
               '<p>범위: ' + esc(diagnostic['scope']) + '</p>']
    if diagnostic.get('reason'):
        content.append('<p>상태 사유: ' + esc(diagnostic['reason']) + '</p>')
    links = '<p>' + ' · '.join('<a href="../' + name + '">' + label + '</a>' for name, label in
        zip(CROSS_PERIOD_TABLES, ('진단 요약 TSV', 'Family TSV', '원래 구간·모델 TSV',
            '사건 TSV', '같은 위치의 근거 TSV', '실제 관측 위치·품질 TSV'))) + '</p>'
    if diagnostic.get('selected_interval_count') is None:
        content.append('<p class="note">선정 구간 진단을 평가하지 않았습니다. 미평가는 관측 0건이나 '
                       'TRM 부재를 뜻하지 않으며 기존 native 관측·판정·adequacy는 유지합니다.</p>')
        content.extend(['<p>' + esc(diagnostic['recommendation']) + '</p>', links])
        return content
    content.append('<p>선정 범위: 사건 ' + esc(diagnostic['selected_event_count']) +
        '건 · pair ' + esc(diagnostic['selected_pair_count']) +
        '개 · 구간 포함 read ' + esc(diagnostic['selected_read_count']) +
        '개 · 원래 native 구간 ' + esc(diagnostic['selected_interval_count']) +
        '개. <b>선택 family와의 모델 일관성 (선택 조건의 영향; 독립 검증 아님)</b>: '
        '후보 모델 비용이 낮은 구간 ' + esc(diagnostic['candidate_lower_cost_interval_count']) +
        '개 · 비교 모델 비용이 낮은 구간 ' + esc(diagnostic['comparator_lower_cost_interval_count']) +
        '개 · 동일 비용 ' + esc(diagnostic['equal_cost_interval_count']) +
        '개. 엄격한 관측 염기 대비 ' + esc(diagnostic['strict_contrast_position_count']) +
        '개는 motif 정확도나 novel 확정 점수가 아닙니다.</p>')
    labels = [('declared_motif', 'Motif'), ('role', '비교 역할'),
              ('selected_event_count', '선정 사건'),
              ('selected_interval_count', '원래 native 구간'),
              ('boundary_supported_event_count', '같은 위치 경계 연결'),
              ('mate_geometry_evaluable_event_count', '같은 위치 mate 평가'),
              ('inward_mate_event_count', '같은 위치 inward mate'),
              ('joint_repeat_inward_event_count', '인접 반복 + inward mate'),
              ('joint_continuous_matching_repeat_inward_event_count', '연속 일치 반복 + inward mate'),
              ('max_joint_adjacent_reference_run_bp', '함께 확인된 reference 반복 최대 bp')]
    rows = diagnostic.get('family_summaries') or []
    if rows:
        content.append('<table><caption>Family별 고유 사건 — reference·방법·HSP 간 합집합</caption>'
            '<thead><tr>' + ''.join('<th scope="col">' + esc(label) + '</th>' for _, label in labels) +
            '</tr></thead><tbody>')
        for row in rows:
            content.append('<tr>' + ''.join('<td>' + esc(row.get(key)) + '</td>' for key, _ in labels) + '</tr>')
        content.append('</tbody></table>')
    content.append('<p class="note">모델은 선정된 원래 native 구간 전체에서 양 가닥·모든 시작 위상과 '
        '모든 최적 편집 경로를 비교합니다. 선택된 family의 비용이 더 낮아도 독립적인 정확도 검증은 아닙니다. '
        '삽입·삭제 위치와 위상은 여러 최적 해를 가질 수 있고, 각 연산 수의 범위를 동시에 달성한다고 가정하지 않습니다. '
        '엄격한 관측 위치 대비가 없다고 두 motif가 동등하거나 후보가 부재하다고 판단하지 않습니다.</p>')
    content.append('<p class="note">Quality는 실제 관측 염기에만 연결한 조건부 Phred+33이며 encoding·calibration·authenticity는 '
        '미검증입니다. Template deletion에는 관측 염기가 없어 quality가 null입니다. '
        '모델의 gap과 reference 정렬의 gap을 같은 사건으로 간주하지 않습니다. 평균 quality 또는 높은 quality만으로 variant 정확성을 확정하지 않습니다.</p>')
    content.append('<p class="note">같은 위치의 근거는 동일 사건·flank context·full-read HSP에서 결합합니다. '
        '인접 reference 반복 길이는 read가 덮은 길이와 다릅니다. 0건은 관측된 양성 사건이 없다는 뜻이며, '
        'mate 평가 분모가 0인 경우 inward 수는 미평가입니다. Reference·방법·HSP 또는 mate를 독립 표본으로 합산하지 않습니다. '
        '검색 상한·표본 연결·구획·기술 서열·FCS 경고를 유지하고, 물리적 말단·main/novel TRM·생물학적 기원을 확정하지 않습니다. '
        'Native 관측·판정·adequacy는 유지합니다.</p>')
    for key, title in [('review_notes', '개별 검토 주석'), ('limitations', '해석의 한계')]:
        if diagnostic.get(key):
            content.append('<details><summary>' + title + '</summary><ul>')
            content.extend('<li>' + esc(note) + '</li>' for note in diagnostic[key])
            content.append('</ul></details>')
    content.append('<p>' + esc(diagnostic['recommendation']) + '</p>')
    content.append(links)
    return content


def reference_localization_html(diagnostic):
    if not diagnostic:
        return []
    def esc(value):
        return html.escape('미평가' if value is None else str(value), quote=True)
    content = ['<h2>선정 경계의 reference 위치 비교</h2>',
               '<p><b>' + esc(diagnostic['display_label']) + '</b></p>',
               '<p>' + esc(diagnostic['summary']) + '</p>',
               '<p>범위: ' + esc(diagnostic['scope']) + '</p>']
    if diagnostic.get('reason'):
        content.append('<p>상태 사유: ' + esc(diagnostic['reason']) + '</p>')
    labels = [('declared_motif', 'Motif'), ('role', '비교 역할'),
              ('selected_event_count', '선정 사건'), ('assessed_event_count', '위치 문맥 평가'),
              ('coordinate_supported_event_count', 'Full-read 좌표 일치'),
              ('boundary_supported_event_count', '실제 경계 횡단'),
              ('mate_geometry_evaluable_event_count', '같은 위치 mate 평가'),
              ('inward_mate_event_count', '같은 위치 inward mate'),
              ('capped_event_count', '검색 상한 경고 사건')]
    def table(title, rows, columns):
        if not rows:
            return
        content.append('<table><caption>' + esc(title) + '</caption><thead><tr>' +
                       ''.join('<th scope="col">' + esc(label) + '</th>' for _, label in columns) +
                       '</tr></thead><tbody>')
        for row in rows:
            content.append('<tr>' + ''.join('<td>' + esc(row.get(key)) + '</td>' for key, _ in columns) + '</tr>')
        content.append('</tbody></table>')
    table('Family별 고유 선정 사건 — reference·방법 간 합집합', diagnostic.get('family_summaries'), labels)
    if diagnostic.get('reference_summaries'):
        content.append('<details><summary>Reference·방법별 비교 펼치기</summary>')
        table('Reference·방법별 계수 — 행 사이 합산 금지', diagnostic['reference_summaries'],
              [('reference_accession', 'Reference'), ('mode', '방법'), *labels])
        content.append('</details>')
    if diagnostic.get('references'):
        content.append('<details><summary>Reference와 원래 표본의 관계</summary><ul>')
        for reference in diagnostic['references']:
            content.append('<li>' + ' · '.join(esc(reference.get(key)) for key in
                ('reference_accession', 'repository_taxon', 'sample_linkage', 'linkage_limit')) + '</li>')
        content.append('</ul></details>')
    for key, title in [('review_notes', '개별 검토 주석'), ('limitations', '해석의 한계')]:
        if diagnostic.get(key):
            content.append('<p>' + title + '</p><ul>')
            content.extend('<li>' + esc(note) + '</li>' for note in diagnostic[key])
            content.append('</ul>')
    content.append('<p>' + esc(diagnostic['recommendation']) + '</p>')
    content.append('<p class="note">선정 사건의 좌표 일치·실제 경계 횡단·같은 위치 inward mate는 서로 다른 근거입니다. '
        'Mate 방향만으로 library insert 적합성·독립 분자·구획 기원·물리적 말단·main/novel TRM을 확정하지 않습니다. '
        '원래 후보·계수·순위·판정·검출력은 유지합니다. 미평가는 관측된 0건과 다릅니다. '
        '<a href="../reference_localization_diagnostics.tsv">진단 요약 TSV</a> · '
        '<a href="../reference_localization_families.tsv">Family 합집합 TSV</a> · '
        '<a href="../reference_localization_references.tsv">Reference·방법 TSV</a> · '
        '<a href="../reference_localization_events.tsv">사건·reference·방법 TSV</a></p>')
    return content


def clade_variant_html(report):
    esc = lambda value: html.escape(str(value), quote=True)
    if not report:
        return []
    content = ['<h2>전역 등록 여부와 분류군 variant</h2>',
               '<p>' + esc(report['display_label']) + '. 전역 catalogue 등록과 분류군 기본 motif 차이를 별도로 비교합니다. '
               '새로운 발견·물리적 말단·main TRM 확정은 아닙니다.</p>']
    if report.get('default_families'):
        content.append('<p>분류군 기본 family: ' + esc(', '.join(report['default_families'])) + '</p>')
    content.append('<p><a href="../clade_variant_observations.tsv">전체 family 비교 TSV</a></p>')
    if report.get('families'):
        content.append('<table><tr><th>Motif</th><th>Reads</th><th>전역 catalogue</th><th>분류군 비교</th><th>종/표본 문헌 관계</th><th>후속 확인</th></tr>')
        for row in report['families'][:100]:
            content.append('<tr>' + ''.join('<td>' + esc(row.get(key, '')) + '</td>' for key in
                           ('motif', 'reads', 'global_catalogue_status', 'status', 'scoped_reference_relationship', 'next_action')) + '</tr>')
        content.append('</table><p>화면에는 최대 100개를 표시하며 전체 비교는 TSV/JSON에 보존합니다.</p>')
    return content


def write_sample_html(path, document):
    esc = lambda value: html.escape(str(value), quote=True)
    top = document['top_candidate']
    adeq = document['data_adequacy']
    qc = adeq.get('read_qc') or {}
    observed = document['repeat_observation']
    native = observed['native_discovery_leader']
    policy_leader = observed['policy_score_leader']
    content = [f'<h1>{esc(document["scientific_name"])} · {esc(document["sample_id"])}</h1>',
               f'<p class="result">{esc(document["decision_summary"]["display_label"])}</p>',
               f'<p><b>{esc(observed["measurement_label"])}</b> · {esc(observed["display_label"])}</p>',
               f'<p>Native 탐색 우선 반복: <b>{esc(native["canonical_family"]) if native else "순위 미제공"}</b> '
               f'({format(native["reads"], ",") + " reads · 원래 관측" if native else "저장된 native 순위로 지정 불가"}) · '
               f'해석 점수 1위 반복: <b>{esc(policy_leader["motif"]) if policy_leader else "현재 양성 관측 없음"}</b> '
               f'({format(policy_leader["reads"], ",") + " reads · 현재 계수" if policy_leader else "미지정"})</p>',
               f'<p>정책상 대표 반복: <b>{esc(top["motif"]) if top else "미지정"}</b> · '
               f'Main TRM: <b>{esc(document["decision_summary"]["main_trm_display"])}</b></p>',
               f'<p>계통 비교: {esc(document.get("trm_interpretation", ""))}</p>',
               f'<p>{esc(document["recommendation"])}</p>',
               '<p class="note">Native 탐색 순위는 참조 annotation 전의 저장된 순위이며 검출·출력 기준을 통과한 반복만 비교합니다. '
               '원래 native 계수와 현재 해석 계수는 다른 측정 경로일 수 있습니다. '
               '해석 점수에는 기대 family·경계 방향의 가중치가 포함됩니다. 정책상 대표는 분류별 선택 규칙을 따릅니다. '
               '어느 1위도 main TRM 확정을 뜻하지 않으며 일반 위성·내부 반복일 수 있습니다. 대표 미지정도 서열 부재가 아닙니다.</p>']
    axes = document.get('evidence_axes', {}).get('display', {})
    content[1:1] = ['<h2>동시 관측·확인 상태</h2>',
        '<table>' + ''.join('<tr><th>' + esc(DISPLAY_NAMES.get(key, key)) + '</th><td>' +
        esc(axes.get(key, '미평가')) + '</td></tr>' for key in AXIS_FIELDS) + '</table>',
        '<p>각 항목은 동시에 성립할 수 있습니다. 전역 motif 등록, 이 종의 기대, 후보 판정, 물리적 말단 확인을 구분합니다.</p>']
    structure = document.get('saved_repeat_structure') or {}
    if structure.get('status') != 'not_provided':
        content += ['<h2>저장한 전체 입력 구조·quality 관측</h2>',
            '<p>' + esc(structure.get('display', {}).get('saved_quality_display', '미평가')) + '</p>',
            '<p>' + esc(structure.get('display', {}).get('saved_mixed_display', '미평가')) + '</p>',
            '<p>별도 원자료의 계수를 보존한 추가 관측입니다. 전체 입력 검사와 선정 구간 비교, 단순 motif와 혼합 배열을 구분합니다. '
            '고정 quality의 원인은 미확정이며, 여러 quality 값이 있어도 실제 오류율을 보정한 것은 아닙니다. '
            '혼합·성분 read 수를 합산하지 않으며 물리적 말단과 main TRM은 별도로 검증해야 합니다.</p>']
    content.extend(assembly_context_html(document.get('assembly_context_evidence'), document.get('quality_path_evidence')))
    content.extend(quality_path_html(document.get('quality_path_evidence')))
    content.extend(long_read_context_report.render_html(document.get('long_read_context_evidence')))
    content.extend(clade_variant_html(document.get('clade_variant_observation')))
    quality = (document.get('quality_observation') or {}).get('current_measurement') or {}
    content.append('<h2>염기 quality 관측</h2><p>' + esc(quality.get('status', 'not_measured')) +
        ' · 후보 구분 염기의 정확도·분자 독립성·말단 회수율은 별도 검증 대상입니다. '
        '<a href="../candidate_quality_observations.tsv">측정된 후보 quality TSV</a></p>')
    content.extend(technical_boundary_html(document.get('technical_boundary_diagnostics')))
    content.extend(reference_localization_html(document.get('reference_localization_diagnostics')))
    content.extend(cross_period_html(document.get('cross_period_evidence')))
    content.extend(input_applicability_html(document.get('input_applicability')))
    content.extend(flank_applicability_html(
        (document.get('input_applicability') or {}).get('flank_period_applicability')))
    content.extend(family_observation_html(document.get('family_observation')))
    comparison_labels = []
    for field, label in [('native_vs_policy_score_same_family', 'Native 우선 ↔ 해석 점수 1위'),
                         ('native_vs_representative_same_family', 'Native 우선 ↔ 정책상 대표')]:
        same = document['policy_selection'][field]
        if same is not None:
            comparison_labels.append(label + ': ' + ('같은 family' if same else '다른 family'))
    if comparison_labels:
        content.append('<p>'+esc(' · '.join(comparison_labels))+
                       '. Motif 회전·역상보 표기가 달라도 canonical family가 같으면 같은 반복입니다.</p>')
    if observed['top_positive_observations']:
        content += ['<h2>현재 양성 반복 관측 · 해석 점수순 상위 5개</h2>',
                    '<table><tr><th>Motif / family</th><th>현재 reads</th><th>원래 discovery reads</th><th>현재 계수 출처</th></tr>']
        for row in observed['top_positive_observations']:
            content.append('<tr><td>'+esc(row['motif'])+' / '+esc(row['canonical_family'])+'</td><td>'+esc(row['reads'])+
                           '</td><td>'+esc(row['original_discovery_reads'])+'</td><td>'+esc(row['evidence_count_source'])+'</td></tr>')
        content.append('</table><p>전체 후보와 순위·계수는 아래 표와 JSON/TSV에 보존됩니다. 원래 양성·현재 0건 충돌도 JSON에 별도 보존됩니다.</p>')
    agreement = document.get('reference_agreement')
    if agreement:
        content.append('<p class="note">Main TRM의 물리적 말단 위치는 현재 short-read 분석으로 확정하지 않습니다. '
                       '대표 후보를 기대 motif로 강제 교체하거나, 기대와 다른 관측을 제거하지 않습니다.</p>')
    expected = ', '.join(document['expected_motifs']) or '확립된 단일 short motif 지정 없음'
    assignment = document['reference_assignment']
    content += ['<h2>분류군 문헌</h2>', f'<p>기대 motif: <b>{esc(expected)}</b> · '
                f'문헌 상태: {esc(assignment.get("reference_status", "미지정"))} · '
                f'반복 구조: {esc(ARCHITECTURE_LABELS.get(assignment.get("repeat_architecture"), "미지정"))}</p>',
                f'<p>{esc(document["reference_context"])}</p>']
    expected_rows = document.get('expected_family_observations', [])
    if expected_rows:
        content.append('<table><tr><th>문헌 예상 family</th><th>관측 reads</th><th>Full-repeat reads</th><th>관측 해석</th></tr>')
        for row in expected_rows:
            explanation = ('미보고: 출력/seed 기준 아래일 수 있으며 명시적 0이 아님' if row['recovery'] == 'not_reported' else
                           '탐지 경로 간 관측 불일치: 이전 양성 관측 보존' if row['recovery'] == 'measurement_disagreement' else
                           'Targeted 검사에서 명시적 0; 생물학적 부재를 뜻하지 않음' if row['recovery'] == 'explicit_zero' else
                           EVIDENCE_LABELS.get(row['candidate_status'], row['candidate_status']) or '반복 서열 관측')
            content.append('<tr><td>'+esc(row['motif'])+'</td><td>'+esc(row['reads'] if row['reads'] is not None else '미보고')+
                           '</td><td>'+esc(row['full_repeat_reads'] if row['full_repeat_reads'] is not None else '미보고')+
                           '</td><td>'+esc(explanation)+'</td></tr>')
        content.append('</table>')
    diagnostic = document.get('provenance', {}).get('expected_family_diagnostic')
    if diagnostic:
        content.append('<p>문헌 기대 family는 같은 전체 입력의 별도 표적 측정으로 보완했습니다. '
                       '원래 후보 순위·경계·pair 계수·검출력 평가는 유지합니다. '
                       '<a href="'+esc(diagnostic['result_file'])+'">표적 결과</a> · '
                       '<a href="'+esc(diagnostic['status_file'])+'">완료·출처 기록</a></p>')
    library = document.get('library_context', {})
    content += ['<h2>시료 정체성·라이브러리와 판정 범위</h2>',
                '<p>'+esc('해석 제약: '+', '.join(library.get('reasons', [])) if library.get('limited') else
                         '확인된 library 제약은 없지만 WGS 표기만으로 순도·말단 회수율·독립 분자를 보장하지 않습니다.')+'</p>',
                '<p class="note">Metadata 충돌은 RNA의 직접 증거가 아닙니다. 알려진 motif 검출은 관측으로 보존하며, '
                '제약이 있는 자료의 미등록 후보·미검출에서 genomic novelty나 부재를 추론하지 않습니다.</p>']
    if library.get('identity_conflict_supported'):
        reviewed = document.get('reviewed_library_evidence') or {}
        content.append('<p class="note"><b>검토된 등록 분류군–서열 불일치:</b> archive의 종 이름은 보존했지만 '
                       '관측 반복을 그 종의 TRM으로 귀속하지 않습니다. 이 표시는 RNA 판정이나 자동 종 이름 교정이 아닙니다. '
                       '<a href="'+esc(reviewed.get('detailed_evidence_path', ''))+'">서열·공개 metadata 검토 근거</a></p>')
    published = document['metadata'].get('published_sample_context')
    if published:
        content.append('<p>같은 accession을 분석한 기존 논문: '
            f'putative TRM <b>{esc(published.get("published_putative_trm", ""))}</b> · '
            f'{esc(published.get("published_classification", ""))}. '
            '<a href="'+esc(published['source'])+'">출처</a></p>')
        content.append('<p class="note">'+esc(published.get('interpretation', ''))+
            ' 같은 입력의 재분석은 독립 생물학적 검증이 아닙니다. 논문의 ND 표기도 TRM 부재 확정이 아닙니다.</p>')
    if top:
        direction = top['directionality']
        fraction = direction.get('consistent_fraction')
        consistency = f'{fraction:.1%}' if fraction is not None else '측정 근거 없음'
        junction_count = top.get('unique_junction_reads')
        junction_display = f'{junction_count:,}' if junction_count is not None else '미측정'
        direction_counts = {key: f'{direction[key]:,}' if direction[key] is not None else '미측정'
                            for key in ('left_consistent', 'right_consistent', 'opposite_junctions')}
        content += ['<h2>관측 근거</h2>',
            f'<p>지지 reads <b>{top["reads"]:,}</b> · 반복으로 덮인 reads <b>{top["full_repeat_reads"]:,}</b> · '
            f'경계가 있는 read records <b>{junction_display}</b> (리드당 1회 집계)</p>',
            f'<div class="geometry">5′ {esc(reverse_complement(top["motif"]))} → flank '
            f'<b>{direction_counts["left_consistent"]}</b> &nbsp; | &nbsp; '
            f'flank → {esc(top["motif"])} 3′ <b>{direction_counts["right_consistent"]}</b></div>',
            f'<p>예상 방향 일치율: <b>{consistency}</b> · 반대 방향 경계 event: {direction_counts["opposite_junctions"]}. '
            '동일 read의 양쪽 경계는 event에 각각 들어갈 수 있습니다.</p>',
            '<p class="note">반복이 read의 ≥90%를 덮으면 full-repeat로 집계하고 경계에서 제외합니다. '
            'Flank는 복잡도 검사를 통과한 서열이며 유일한 subtelomeric locus로 확인한 것은 아닙니다. '
            '왼쪽/오른쪽 read 경계는 같은 물리적 locus의 두 strand에서도 나올 수 있습니다. '
            '이 개수는 PCR 중복을 제거한 독립 분자 수가 아닙니다. '
            '이 도식은 예상 read 경계 방향이며 서로 다른 염색체 말단의 직접 관측을 뜻하지 않습니다.</p>']
        if top.get('paired_evidence', {}).get('measured'):
            content += ['<h3>Paired-end 일관성 · 보조 관측</h3>',
                '<p>'+esc(paired_summary(top))+'</p>',
                '<p class="note">현재 read 계수와 같은 측정에서 얻은 FR mate 쌍입니다. '
                '같은 family가 양쪽 mate에서 검출될 때 역상보 방향인지 확인합니다. '
                '한쪽 mate에만 반복이 있어도 정상일 수 있으며 제외 사유가 아닙니다. '
                '내부 반복도 역상보 pair를 만들 수 있습니다. 이 항목은 TRM·novel 후보 승격, '
                '검출력 판정이나 독립 분자 수 추정에 사용하지 않습니다.</p>']
        columns = validation_columns(top)
        readiness = top['confirmation_readiness']
        internal = columns['clean_two_flank_records']
        content += ['<h2>후보 단계와 추가 검증 근거</h2>',
            '<p><b>'+esc(EVIDENCE_LABELS.get(top['candidate_status'], top['candidate_status']))+'</b></p>',
            '<p>서로 다른 경계 read 서열: <b>'+esc(columns['distinct_junction_sequences'])+'</b> · '
            '20-bp flank 군집: <b>'+esc(columns['physical_flank_clusters'])+'</b> · '
            '중복 보정 방향성: <b>'+esc(columns['deduplicated_orientation'])+'</b></p>',
            '<p>양쪽 clean flank를 가진 반복 read: <b>'+esc(internal if internal is not None else '미측정')+'</b>. '
            '이런 관측은 내부 배열 여부를 검토할 근거입니다. Genuine telomere와 interstitial copy가 함께 존재할 수 있어 motif의 부정 판정은 아닙니다.</p>',
            '<p>추가 확인: '+esc(columns['confirmation_missing'] or '측정된 잠정 후보 조건 충족 여부와 별개로, 독립 자료 및 물리적 염색체 말단 확인 필요')+'</p>',
            '<p class="note">동일 read와 역상보 서열을 합쳐 세었으며 독립 분자 수는 아닙니다. '
            'Flank 군집은 20-bp Hamming ≤2 또는 확장 flank의 ±16-bp offset·≥30-bp 겹침·≤2 mismatch 연결 성분입니다. 독립 locus·염색체 수가 아니고 서로 다른 offset도 영향을 줍니다. '
            '포화 sketch의 고유 수는 하한이며, 관측 subset 비율을 전체 비율로 사용하지 않습니다. '
            'Flank sketch가 포화하면 미관측 서열이 군집을 연결할 수 있어 전체 군집 수를 미확정으로 둡니다. '
            'Flank 군집은 quality 기반 오류 모델이 아니며, 선택적으로 측정한 quality 분포는 별도 항목에서 확인합니다.</p>',
            '<h3>방향과 flank 군집의 결합 관측</h3>',
            '<p>일치 방향 사건이 있는 군집: <b>'+esc(columns['consistent_flank_clusters'])+'</b> · '
            '반대 방향 사건이 있는 군집: <b>'+esc(columns['opposed_flank_clusters'])+'</b></p>',
            '<p>보관 일치 사건의 최대 군집 집중도: <b>'+esc(columns['retained_consistent_cluster_concentration'])+'</b> · '
            +esc(columns['flank_context_completeness'])+'</p>',
            '<p>방향 기준: '+esc(columns['flank_orientation_basis'])+'</p>',
            '<p class="note">전체 관측에서 정한 방향을 고정하고 같은 방향의 사건이 어느 군집에 있는지 함께 봅니다. '
            '군집마다 다수 방향을 다시 고르지 않으며 이 진단은 기존 판정을 바꾸지 않습니다. '
            '포화 sketch의 보관 군집 수는 전체 군집 수의 하한도 아닙니다. 빠진 anchor가 군집을 연결할 수 있습니다. '
            '집중도는 보관 사건의 기술통계이며 독립 분자·고유 locus·염색체 말단 수가 아닙니다. '
            '관측 방향에서 추정한 반복의 생물학적 G-strand는 별도 확인이 필요합니다.</p>',
            '<p class="note">미등록 반복 → 경계 후보 → 추가 근거가 있는 잠정 후보의 순서입니다. '
            '외부 말단 검증 완료 단계는 short-read 검사로 부여하지 않습니다. 500 bp는 보편적인 telomere 길이 조건이나 후보 승격 문턱이 아닙니다.</p>']
        if top.get('long_array_full_to_junction_warning'):
            content.append('<p class="note">Full-repeat보다 junction reads가 많습니다. 단순 장배열·무작위 fragment 모델과의 차이를 검토하되 read 길이·array 길이·library 효과가 있으므로 일률적 제외 조건으로 사용하지 않습니다.</p>')
    content += ['<h2>데이터와 미검출 해석</h2>', f'<p><b>{esc(ADEQUACY_LABELS[adeq["status"]])}</b></p>',
                f'<p>전체 reads: {qc.get("records",0):,} · 전체 bases: {qc.get("bases",0):,} · '
                f'평균 read 길이: {qc.get("mean_record_bases") or 0:.1f} bp</p>']
    blockers = adeq.get('library', {}).get('uniform_wgs_model_blockers', [])
    if blockers:
        content.append('<p>모델 적용 전 확인할 library 조건: '+esc(', '.join(blockers))+'</p>')
    if adeq.get('missing_assumptions'):
        content.append('<p>Genome size·염색체 수·무작위 WGS 가정·말단 회수율·유효 독립 read 비율·'
                       '탐지 민감도 중 필요한 정보가 빠져 있습니다. 많은 read 수와 WGS 표기만으로 '
                       '기술적 미검출을 배제하지 않습니다.</p>')
    show_pairs = any(c.get('paired_evidence', {}).get('measured') for c in document['candidates'])
    content += ['<h2>관측 후보</h2>', '<p>주 판정과 보조 후보를 구분하세요. 미등록 후보는 이 catalogue에 '
                '없는 family라는 뜻이며 문헌 전체에서 새롭다고 확정한 것이 아닙니다.</p>',
                '<table><tr><th>Motif</th><th>근거 단계</th><th>Reads</th><th>Full-repeat</th><th>원래 방향 일치</th><th>서로 다른 경계 서열</th><th>Flank 군집</th><th>일치 방향 flank 군집</th><th>보관 일치 사건의 최대 군집 집중도</th><th>추가 확인</th>'+
                ('<th>Paired-end 일관성 · 보조 관측</th>' if show_pairs else '')+'</tr>']
    for candidate in document['candidates'][:25]:
        ratio = candidate['directionality'].get('consistent_fraction')
        columns = validation_columns(candidate)
        content.append('<tr><td>'+esc(candidate['motif'])+'</td><td>'+esc(EVIDENCE_LABELS.get(candidate['candidate_status'], candidate['candidate_status']))+
                       f'</td><td>{candidate["reads"]:,}</td><td>{candidate["full_repeat_reads"]:,}</td><td>'+
                       (f'{ratio:.1%}' if ratio is not None else '미측정')+'</td><td>'+esc(columns['distinct_junction_sequences'])+
                       '</td><td>'+esc(columns['physical_flank_clusters'])+'</td><td>'+esc(columns['consistent_flank_clusters'])+
                       '</td><td>'+esc(columns['retained_consistent_cluster_concentration'])+'</td><td>'+esc(columns['confirmation_missing'])+'</td>'+
                       ('<td>'+esc(paired_summary(candidate))+'</td>' if show_pairs else '')+'</tr>')
    content.append('</table>')
    if show_pairs:
        content.append('<p class="note">Pair 값은 같은 현재 측정의 read 계수에만 연결합니다. 미측정은 0이 아니며, '
                       'pair 방향은 일관성 보조 관측으로만 사용합니다. 말단 위치·신규성·데이터 충분성을 확정하지 않습니다.</p>')
    if len(document['candidates']) > 25:
        content.append(f'<p>표에는 점수순 25개를 표시했습니다. 전체 {len(document["candidates"]):,}개 관측 후보는 '
                       '<a href="../candidates.tsv">후보 TSV</a>와 전체 JSON에 보존되어 있습니다.</p>')
    selection = document['target_selection']
    unmeasured = selection.get('novel_candidates_not_targeted', 0)
    if unmeasured:
        content.append(f'<p class="note">전체 입력 read를 읽었지만 기존 결과 업그레이드의 경계 측정은 알려진 family와 '
                       f'우선순위 후보 {selection.get("novel_targets", 0):,}개에 한정됐습니다. '
                       f'선정 기준을 통과한 나머지 미등록 family {unmeasured:,}개는 경계를 측정하지 않았습니다. '
                       '이들의 미검출 또는 방향 불일치를 뜻하지 않습니다. 원래 반복 관측값은 보존되어 있습니다.</p>')
    scenarios = adeq.get('illustrative_grid', {}).get('scenarios', [])
    if scenarios:
        content += ['<details><summary>가상 조건별 검출력 보기 — 실제 샘플 추정값이 아닙니다</summary>',
                    '<p>모든 2 × 염색체 수 말단에 같은 motif의 500 bp 배열·유효 독립 비율 0.5·탐지 민감도 0.5를 가정했습니다. '
                    '이 값들은 실제 sample에서 측정되지 않았으며 충분성 판정을 올리지 않습니다. '
                    '확률은 raw detector 임계값 미도달에 관한 조건부 값이며 부재 확률이 아닙니다.</p>',
                    '<table><tr><th>Genome</th><th>Haploid 염색체 수</th><th>말단 회수율</th><th>조건부 미검출 상한</th></tr>']
        for scenario in scenarios:
            miss = scenario['worst_case_miss_probability_upper_bound']
            formatted = '부동소수점 표시 한계 미만' if miss == 0 else f'{miss:.3g}'
            content.append(f'<tr><td>{scenario["genome_size_bp"]/1e6:g} Mb</td><td>{scenario["haploid_chromosome_count"]}</td>'
                           f'<td>{scenario["telomere_retention_fraction"]:.0%}</td><td>{formatted}</td></tr>')
        content.append('</table></details>')
    content.extend(followup_html(document.get('followup_diagnostics')))
    validation = document['independent_validation']
    content += ['<h2>후속 확인</h2>', f'<p>{esc(validation["download_advice"])}</p>',
        f'<p>현재 BioSample: {esc(validation.get("current_biosample") or "미확인")}</p>',
        '<p><a href="'+esc(validation['sra_search_url'])+'">NCBI SRA에서 해당 분류군 WGS 검색</a></p>',
        '<p>다른 BioSample·분리주·라이브러리를 우선하고 ReadHash/RunHash 및 연구 출처를 확인하세요. '
        '검색 링크는 독립성이 확인된 대체 데이터의 추천 목록이 아닙니다.</p>',
        '<p class="note">TRM 부재, >32 bp 반복, 불규칙 말단 서열, library 문제는 short-read 미검출 하나로 '
        '확정 구분할 수 없습니다. ALT 등 유지 기전도 별도 확인이 필요합니다.</p>',
        f'<p><a href="{esc(document["sample_id"])}.json">전체 JSON·출처·계산값</a> · <a href="../index.html">전체 샘플</a></p>']
    duplicates = validation['same_archive_read_content_group']
    if duplicates:
        content.append('<p class="note">Archive ReadHash/RunHash 중복 등록군: '+esc(', '.join(duplicates))+
                       '. 서로를 독립 검증 데이터로 세지 않습니다.</p>')
    Path(path).write_text('<!doctype html><html lang="ko"><meta charset="utf-8"><title>'+esc(document['sample_id'])+
        '</title><style>body{font:16px system-ui;color:#193040;max-width:1100px;margin:35px auto;padding:0 20px;line-height:1.65}'
        'h1{font-size:25px}h2{font-size:20px;margin-top:30px}.result{background:#eaf4f7;padding:16px;font-size:20px}'
        '.geometry{padding:12px;background:#f3f6f8;font-family:monospace}.note{color:#526775;font-size:14px}'
        'table{border-collapse:collapse}td,th{padding:7px 12px;border-bottom:1px solid #ddd;text-align:left}'
        'a{color:#075f97}details{margin:25px 0}summary{cursor:pointer;font-weight:600}</style>'+''.join(content)+'</html>')


def read_assumptions(path):
    if not path:
        return {}, {}
    value = json.loads(Path(path).read_text())
    return value.get('defaults', {}), value.get('samples', {})


def sra_link(meta):
    taxid = str(meta.get('tax_id', '')).strip()
    query = f'txid{taxid}[Organism:exp] AND "strategy wgs"[Properties]' if taxid.isdigit() else (
        '"'+meta.get('scientific_name', '')+'"[Organism] AND "strategy wgs"[Properties]')
    return 'https://www.ncbi.nlm.nih.gov/sra/?'+urlencode({'term': query})


def sample_reference_index(references, assignment, taxon_rows):
    """Add exact published taxon references without expanding global canon.

    Match the curated reference ID and motif together; an expected motif alone
    never establishes a published reference or its orientation.
    Applicability is established by prepare_taxon_priors; this join checks
    IDs/families/status and does not re-infer taxonomic scope. Additions remain
    sample-local and do not change the global catalogue. Candidate/template
    rows retain their candidate status and never become established evidence.
    """
    result = {family: list(rows) for family, rows in references.items()}
    if assignment.get('reference_status') not in ('supported', 'candidate'):
        return result
    supported_ids = set(assignment.get('reference_ids', '').split(','))
    candidate_ids = set(assignment.get('review_reference_ids', '').split(','))
    supported_families = {canonical(m) for m in assignment.get('related_motifs', '').split(',') if m}
    candidate_families = {canonical(m) for m in assignment.get('review_motifs', '').split(',') if m}
    for row in taxon_rows:
        supported = (assignment.get('reference_status') == 'supported'
                     and row.get('reference_id') in supported_ids
                     and row.get('reference_status') == 'supported'
                     and row.get('evidence') in SUPPORTED_EVIDENCE)
        candidate = (row.get('reference_id') in candidate_ids
                     and row.get('reference_status') == 'candidate'
                     and bool(row.get('evidence')) and row.get('source', '').startswith(('https://', 'http://')))
        if (not (supported or candidate)
                or row.get('repeat_architecture') not in ('exact_short_unit', 'mixed_short_units')):
            continue
        for motif in row.get('related_motifs', '').split(','):
            if not 3 <= len(motif) <= 32 or set(motif) - set('ACGT'):
                continue
            family = canonical(motif)
            if family not in (supported_families if supported else candidate_families):
                continue
            reference = {'motif': motif, 'evidence': row['evidence'], 'source': row['source'],
                         'reference_id': row['reference_id'], 'scope': row['scope'],
                         'reference_status': row['reference_status'],
                         'repeat_architecture': row['repeat_architecture'],
                         'scientific_name': row['scientific_name']}
            result.setdefault(family, []).append(reference)
    return result


def reference_context(assignment):
    status = assignment.get('reference_status', '')
    architecture = assignment.get('repeat_architecture', '')
    if status == 'out_of_scope' or architecture in ('long_complex_family', 'irregular', 'retrotransposon'):
        return '문헌상 긴·복합·불규칙 말단 구조: 검출된 짧은 반복만으로 이 종의 말단 motif를 지정하지 않습니다.'
    if status == 'candidate':
        return '종/근연 분류군의 문헌 후보입니다. Template·서열 추정은 말단 위치의 실험적 확정과 구분하며, 현재 표본의 관측과 비교합니다.'
    if status in ('unknown', 'mixed'):
        return '종/근연 분류군의 문헌은 단일 확립 motif를 지정하지 않습니다. 상위 분류군 기본값을 강제하지 않습니다.'
    if architecture == 'mixed_short_units':
        return '문헌상 여러 단위가 섞인 말단 반복입니다. 개별 component 검출로 전체 배열 구조를 복원했다고 해석하지 않습니다.'
    if assignment.get('reference_scope') == 'clade':
        return '분류군의 대표 motif이며 해당 종에서 검증된 말단 서열이라는 뜻은 아닙니다.'
    return '문헌의 정확한 motif family를 참조했습니다. 현재 샘플의 말단 위치는 별도 검증 대상입니다.'


def sample_report(sample, meta, assignment, args, references, assumptions, duplicate_groups):
    global_references = references
    references = sample_reference_index(references, assignment, getattr(args, 'taxon_reference_rows', []))
    source = args.source/sample/'result.json'
    raw = json.loads(source.read_text()) if source.is_file() else {'summary': {}, 'candidates': []}
    source_hash = sha256(source) if source.is_file() else None
    native_integrity = validate_native_result(source, sample, source_hash)
    if args.evidence is None:
        # New v0.4 scans already collected these measurements in their one
        # de novo pass; never require a redundant FASTQ pass for a fresh scan.
        evidence_path = source
        all_rows = {canonical(c['motif']): c for c in raw.get('candidates', [])+raw.get('rejected_candidates', [])}
        complete = (native_integrity['analysis_eligible'] and 'evidence_qc' in raw
                    and raw.get('summary', {}).get('sampled_prefix') is False
                    and all('junction_evidence' in c for c in all_rows.values()))
        scan_status = {'status': 'complete' if complete else 'missing_native_evidence',
                       'operation': 'reuse_native_scan',
                       'reason': '' if complete else 'Full native scan evidence is required; legacy results need --evidence.'}
        targeted = {'summary': raw.get('summary', {}), 'evidence_qc': raw.get('evidence_qc'),
                    'candidate_evidence': list(all_rows.values()),
                    'paired_end': raw.get('paired_end'),
                    'quality_diagnostics': raw.get('quality_diagnostics')} if complete else {}
    else:
        evidence_path = args.evidence/sample/'evidence.json'
        status_path = args.evidence/sample/'status.json'
        scan_status = json.loads(status_path.read_text()) if status_path.is_file() else {'status': 'not_run'}
        complete = native_integrity['analysis_eligible'] and scan_status.get('status') == 'complete'
        targeted = json.loads(evidence_path.read_text()) if complete else {}
    if not native_integrity['analysis_eligible']:
        scan_status = {'status': native_integrity['status'], 'operation': 'validate_native_source',
                       'reason': native_integrity['reason']}
    evidence_hash = (source_hash if args.evidence is None else sha256(evidence_path)) if complete else None
    discovery_quality = describe_quality(raw, source_hash)
    measured_quality = describe_quality(targeted, evidence_hash)
    if complete and args.evidence is not None:
        if scan_status.get('identity', {}).get('source_sha256') != source_hash:
            raise ValueError(f'{sample}: evidence source hash differs from original discovery result')
        if scan_status.get('result_sha256') != evidence_hash:
            raise ValueError(f'{sample}: evidence result hash differs from recorded completed scan')
    for field in ('scientific_name', 'tax_id', 'lineage'):
        if str(assignment.get(field, '')).strip() != str(meta.get(field, '')).strip():
            raise ValueError(f'{sample}: reference assignment {field} differs from metadata')
    expected_motifs = [m for m in assignment.get('related_motifs', '').split(',') if m]
    # A reviewed candidate is an expectation for comparison, never an emitted
    # established native prior. Require the same ID/family join as above.
    review_ids = set(assignment.get('review_reference_ids', '').split(','))
    for motif in assignment.get('review_motifs', '').split(','):
        if motif and any(r.get('reference_id') in review_ids and r.get('reference_status') == 'candidate'
                         for r in references.get(canonical(motif), [])):
            expected_motifs.append(motif)
    expected_motifs = list(dict.fromkeys(expected_motifs))
    expected_families = {canonical(m) for m in expected_motifs}
    # Verify geometry and original scan describe the same number of sequences.
    if complete:
        for key in ('records', 'bases'):
            if raw['summary'][key] != targeted['summary'][key]:
                raise ValueError(f'{sample}: {key} mismatch')
        if targeted['summary'].get('sampled_prefix'):
            raise ValueError(f'{sample}: prefix-only geometry is not a complete analysis')
    power_input = {'summary': raw.get('summary', {}), 'config': raw.get('config', {})}
    # The flank-length model is version/option gated. Preserve absence instead
    # of inventing current defaults for legacy or separately targeted scans.
    native_model_metadata = {key: raw[key] for key in (
        'program', 'version', 'analysis_kind', 'junction_options', 'flank_period_annotation_method') if key in raw}
    power_input.update(native_model_metadata)
    if targeted.get('evidence_qc'):
        power_input['evidence_qc'] = targeted['evidence_qc']
    if targeted.get('quality_diagnostics') is not None:
        power_input['quality_diagnostics'] = targeted['quality_diagnostics']
    if not complete:
        power_input['analysis_status'] = 'failed'
    adequacy = assess_adequacy(power_input, meta, assumptions)
    if not complete:
        adequacy.update(status='observed_qc_failure', analysis_failure=scan_status.get('reason', 'Evidence scan incomplete'))
    raw_rows = {canonical(c['motif']): c for c in raw.get('candidates', [])+raw.get('rejected_candidates', [])}
    evidence_rows = {c['motif']: c for c in targeted.get('candidate_evidence', [])}
    # Keep all original observed candidates in a traceable archived source.
    # Interpret every newly measured target, all known families, and any expected family.
    families = set(raw_rows) | set(evidence_rows) | set(references) | expected_families
    mean_length = raw['summary'].get('bases', 0)/max(raw['summary'].get('records', 0), 1)
    candidates = []
    for family in sorted(families):
        original = raw_rows.get(family, {'motif': family, 'reads': 0, 'repeat_bases': 0,
                                      'max_tract': 0, 'full_repeat_reads': 0})
        candidate = assess_candidate(original, evidence_rows.get(family), references, expected_families, mean_length,
                                     raw_paired_end=raw.get('paired_end'),
                                     measured_paired_end=targeted.get('paired_end'))
        measured_counts = candidate['evidence_count_source'] in {'measured_evidence', 'targeted_only'}
        quality = describe_candidate_quality(evidence_rows.get(family, {}) if measured_counts else original,
                                             measured_quality if measured_counts else discovery_quality)
        if quality is not None:
            candidate['quality_diagnostics'] = quality
        # Expected but unreported families belong in the expectation table,
        # not a fabricated zero-count row in the observed candidate table.
        # Explicit measured zeros and positive historical evidence are retained.
        if (candidate['reads'] or candidate['original_discovery_reads']
                or (family in expected_families and (family in raw_rows or family in evidence_rows))):
            candidates.append(candidate)
    candidates.sort(key=lambda c: (-c['ranking_score'], c['canonical_family']))
    scope = resolved_eukaryotic_scope(meta)
    reviewed = getattr(args, 'library_evidence_rows', {}).get(sample)
    if reviewed and reviewed.get('previous_full_read_records') and raw.get('summary', {}).get('records') is not None:
        if int(reviewed['previous_full_read_records']) != raw['summary']['records']:
            raise ValueError(f'{sample}: reviewed library evidence concerns different full-input read totals')
    if reviewed and reviewed.get('library_review_code') == 'taxonomic_identity_conflict':
        # A separately reviewed positive identity conflict limits attribution to
        # the archive taxon. Preserve that taxon and every observed family.
        if reviewed.get('source_result_sha256') != source_hash:
            raise ValueError(f'{sample}: reviewed identity conflict source result hash differs')
        scope = False
    library = library_context(meta, reviewed, adequacy)
    code, label, action = classify(candidates, adequacy, complete, scope, assignment, library)
    top = select_representative(candidates, code, assignment)
    observation = sequence_observation(top) if top else {
        'status': 'no_sequence_observed', 'observed': False, 'current_reads': None,
        'interpretation': 'No representative supported sequence is reported; this is not a complete zero-hit measurement or biological absence.'}
    expected_observations = []
    by_family = {c['canonical_family']: c for c in candidates}
    for motif in expected_motifs:
        family = canonical(motif)
        candidate = by_family.get(family)
        saved = family in raw_rows or family in evidence_rows
        recovery = ('reported_positive' if candidate and candidate['reads'] > 0 else
                    'measurement_disagreement' if candidate and candidate.get('measurement_disagreement') else
                    'explicit_zero' if saved and candidate and candidate.get('evidence_count_source') in ('measured_evidence', 'targeted_only')
                    else 'not_reported')
        expected_observations.append({
            'motif': motif, 'canonical_family': family,
            'recovery': recovery,
            'reads': candidate['reads'] if recovery != 'not_reported' else None,
            'full_repeat_reads': candidate['full_repeat_reads'] if recovery != 'not_reported' else None,
            'candidate_status': candidate['candidate_status'] if candidate else None,
            'reference_status': candidate['reference_status'] if candidate else None,
            'source': 'saved_repeat_evidence' if saved else 'family_not_in_saved_report',
            'interpretation': 'Sequence observation does not establish a physical terminal location. Missing native families may fall below reporting/seed thresholds; missing is not explicit zero.'})
    expected_observations, expected_diagnostic = supplement_expected_observations(
        sample, expected_observations, getattr(args, 'expected_evidence', None), raw,
        source_hash, analysis_complete=complete)
    agreement = reference_agreement(candidates, top, assignment, expected_observations,
                                    analysis_complete=complete, biological_scope=scope,
                                    library=library, qc_status=adequacy['status'])
    agreement_label = REFERENCE_AGREEMENT_LABELS.get(agreement['status'], agreement['status'])
    review_recommendation = reference_review_text(agreement)
    expected_observation_text = '; '.join(
        f'{o["motif"]}: {o["reads"]:,} reads / full {o["full_repeat_reads"]:,}' if o['reads'] is not None
        else f'{o["motif"]}: 미보고 (명시적 0 아님)' for o in expected_observations) or '문헌 예상 family 미지정'
    recommendations = {
        'canonical_screen_complete': '알려진 motif 검출. 경계 방향성 수준은 별도 확인하며, 새 motif 검증 대상에서는 제외합니다.',
        'review_boundary_support': '알려진 반복은 검출됐지만 말단 경계 근거는 부족합니다. 미검출/새 TRM으로 분류하지 않으며 말단 확인이 필요하면 독립 자료를 검토하세요.',
        'review_observed_repeat': '알려진 반복 서열은 관측됐습니다. Full-repeat·경계 방향·내부 배열 가능성을 검토하고 말단 패턴 판정과 구분하세요. 다른 반복 후보는 보존합니다.',
        'validate_published_candidate': '문헌 후보 서열이 관측됐습니다. Template/계산 추정은 실험적 말단 검증과 구분하며, 독립 long-read·assembly로 위치와 방향을 확인하세요.',
        'review_trace_observation': '문헌 family의 소수 지지 read가 있습니다. 무검출로 처리하지 않으며 오류·오염·낮은 회수율을 검토하고 독립 자료로 확인하세요.',
        'validate_independent_biosample': '다른 BioSample·분리주·라이브러리의 WGS 또는 long-read로 재검증하세요. 같은 ReadHash의 재등록 자료는 독립 검증으로 세지 않습니다.',
        'validate_alternative_telomere_hypothesis': '명시한 가정 아래 단순 short TRM의 미검출 가능성이 낮습니다. 독립 WGS와 long-read/assembly로 미등록·긴·불규칙 말단 서열 가설을 확인하세요.',
        'resolve_expected_family_search_scope': '기대 family의 주기가 원래 탐색 또는 검출력 모델 범위에 포함되지 않았습니다. 해당 family를 평가한 설정·표적 측정과 검출력 근거를 확인하세요. 현재 결과로 canonical TRM 미검출을 주장하지 않습니다.',
        'review_scanned_period_nondetection': '명시한 가정과 실제 검사한 주기에 한정한 미검출입니다. 종별 기대 family를 확인하고, 검사하지 않은 주기·불규칙 배열 및 라이브러리 회수율을 별도로 검토하세요.',
        'resolve_adequacy_then_validate': 'Genome size·염색체 수·library prep·말단 회수율·유효 독립 read 수가 부족해 기술적 미검출을 배제할 수 없습니다. 메타데이터 확인 또는 독립 WGS 검증이 필요합니다.',
        'replace_or_check_library': '데이터/라이브러리 제약을 확인하고 적합한 독립 genomic WGS로 다시 분석하세요.',
        'review_measurement_disagreement': '이전 양성 관측과 새 경로의 0건이 충돌합니다. 탐지 설정·민감도 차이를 확인한 뒤 미검출을 해석하세요.',
        'obtain_more_suitable_wgs': '명시한 검출 모델에서 데이터가 부족합니다. 적합한 WGS를 추가 확보하세요.',
        'repair_or_replace_input': 'FASTQ 입력/처리 오류를 해결한 뒤 재분석하세요. 생물학적 미검출로 해석하지 않습니다.',
        'check_host_or_mixed_sample': '등록 분류군과 host/혼합 시료 관계를 확인하세요. 관측 반복의 생물학적 출처는 미확인입니다.',
        'review_known_architecture': '해당 종은 문헌상 불규칙·긴·레트로트랜스포존 기반 말단 구조입니다. 짧은 반복 후보는 보존하지만 종의 TRM으로 지정하지 않습니다. 알려진 구조와 후보의 물리적 위치를 독립 자료로 검토하세요.',
        'resolve_library_context': 'RNA/Hi-C 또는 library 표기·representation 제약을 확인하세요. Metadata 충돌만으로 RNA를 확정하지 않으며, 관측 반복은 보존하고 genomic novelty·부재 해석은 보류합니다.',
        'complete_candidate_evidence': '경계 후보는 관측됐습니다. 중복 보정 서열·flank 다양성과 양쪽 비반복 flank 근거를 확인하고, 독립 genomic 자료 및 말단 위치 증거를 확보하세요.',
        'validate_terminal_location': '추가 서열·flank 기준을 통과한 잠정 후보입니다. 독립 genomic 자료와 long-read/assembly의 물리적 말단 위치로 확인하세요. Short reads만으로 새 TRM을 확정하지 않습니다.',
    }
    targets = scan_status.get('identity', {}).get('targets', {})
    old_metrics_path = args.source/sample/'original_scan_metrics.json'
    if not old_metrics_path.is_file():
        old_metrics_path = args.source/sample/'metrics.json'
    old_metrics = json.loads(old_metrics_path.read_text()) if old_metrics_path.is_file() else {}
    execution = {'sample_id': sample, 'analysis_complete': complete,
        'evidence_operation': scan_status.get('operation', 'targeted_full_input_evidence'),
        'original_scan_wall_seconds': old_metrics.get('wall_seconds'),
        'original_scan_max_rss_kb': old_metrics.get('max_rss_kb'),
        'added_evidence_wall_seconds': scan_status.get('wall_seconds'),
        'added_evidence_user_cpu_seconds': scan_status.get('user_cpu_seconds'),
        'added_evidence_system_cpu_seconds': scan_status.get('system_cpu_seconds'),
        'added_evidence_max_rss_kb': scan_status.get('max_rss_kb'),
        'added_evidence_threads': scan_status.get('identity', {}).get('threads'),
        'input_compressed_bytes': scan_status.get('input_bytes'),
        'records': raw['summary'].get('records'), 'bases': raw['summary'].get('bases'),
        'comparison_limit': 'Original all-period discovery and added targeted evidence perform different work; concurrent load/threads differ. These are measured stages, not a speedup benchmark.'}
    document = {
        'schema_version': 22, 'program': 'trm-assessment', 'version': POLICY['version'],
        'reporting_version': REPORTING_VERSION,
        'native_scan_version': raw.get('version'),
        'reporting_features': {'paired_consistency_annotation': 1,
                              'directional_flank_cluster_annotation': 1,
                              'directional_flank_clusters_change_classification_or_adequacy': False,
                              'pair_evidence_changes_classification_or_adequacy': False,
                              'reference_agreement': 1,
                              'expected_diagnostic_table_preserves_native_candidates': True,
                              'observation_decision_separation': 1,
                              'observation_view_changes_classification_or_adequacy': False},
        'sample_id': sample, 'scientific_name': meta.get('scientific_name'), 'tax_id': meta.get('tax_id'),
        'classification': code, 'result': label, 'analysis_complete': complete,
        'detected_motif': top['motif'] if top else None, 'top_candidate': top,
        'sequence_observation': observation, 'expected_family_observations': expected_observations,
        'reference_agreement': agreement, 'trm_interpretation': agreement_label,
        'expected_motifs': expected_motifs, 'reference_assignment': assignment,
        'reference_context': reference_context(assignment),
        'library_context': library, 'reviewed_library_evidence': reviewed,
        'motif_assignment_scope': 'Observed repeat family, not assignment of a species terminal motif; inspect primary context and evidence level.',
        'externally_confirmed': False,
        'data_adequacy': adequacy, 'candidates': candidates,
        'quality_observation': {'current_measurement': measured_quality,
                                'original_discovery': discovery_quality,
                                'changes_classification_or_adequacy': False},
        'native_performance_profile': raw.get('performance_profile'),
        'expected_family_power_scope': expected_family_power_scope(adequacy, assignment),
        'execution': execution,
        'next_action': agreement['next_action'] if review_recommendation else action,
        'recommendation': review_recommendation or recommendations[action],
        'observation_next_action': action, 'observation_recommendation': recommendations[action],
        'independent_validation': {'sra_search_url': sra_link(meta),
            'same_archive_read_content_group': duplicate_groups.get(sample, []),
            'current_biosample': meta.get('sample_accession'),
            'independence_rule': 'Prefer another biological specimen/BioSample and separate library; verify ReadHash/RunHash and provenance. Another accession alone is insufficient.',
            'download_advice': 'SRA 검색 결과에서 독립 genomic WGS 또는 long-read를 골라 다운로드한 후 같은 분석을 수행하세요.',
            'search_is_live_link_not_a_verified_alternative_dataset': True},
        'unresolved_alternatives': [] if code == 'known_trm_detected' else [
            'Technical nondetection including inadequate representation of chromosome ends',
            'Uncatalogued short tandem repeat', 'Long (>32-bp), mixed or irregular telomeric sequence',
            'Absence of a simple short tandem TRM / other telomere architecture; maintenance mechanism unresolved'],
        'limitations': ['Short-read boundary orientation does not prove physical chromosome-terminal location.',
                       'Both read-edge orientations can arise from opposite strands of one physical junction; they do not establish distinct chromosome ends.',
                       'No result establishes ALT or telomere absence.',
                       'Novel-family screening and inferred orientation are not calibrated for multiple-family false discovery rates.',
                       'RC-collapsed distinct sequences are not independent molecules; optional quality histograms do not assess variant-discriminating bases, independent molecules or PCR effects.',
                       'Flank clusters merge Hamming<=2 anchors or >=30-bp representative overlaps at offsets +/-16 with <=2 mismatches; these are operational groups, not loci. Saturated sketches cannot establish full cluster count.',
                       'Nonrepeat flanks are not assumed unique; karyotype is not inferred from junction count.',
                       'Canonical means an exact published family in the global catalogue or a supported exact taxon reference; it is not terminal localization in this sample.'],
        'target_selection': targets.get('selection', {}) if args.evidence is not None else {
            'mode': 'native_full_discovery_evidence', 'novel_candidates_not_targeted': 0,
            'rule': 'All raw observed families retained and assessed; native geometry reused. Missing validation in older scans remains unmeasured, not zero.'}, 'policy': POLICY,
        'input_summary': raw.get('summary', {}), 'metadata': meta,
        'provenance': {'original_discovery_result': str(source.resolve()),
            'original_discovery_sha256': source_hash,
            'native_result_validation': native_integrity,
            'targeted_evidence_result': str(evidence_path.resolve()) if complete else None,
            'targeted_evidence_sha256': evidence_hash,
            'targeted_scan_status': scan_status.get('status'),
            'all_original_repeat_families_retained_in_source': len(raw_rows),
            'catalogue_sha256': sha256(args.catalog), 'assignment_sha256': sha256(args.assignments),
            'taxon_reference_sha256': sha256(args.taxon_references) if getattr(args, 'taxon_references', None) else None,
            'expected_family_diagnostic': expected_diagnostic,
            'library_evidence_table_sha256': reviewed.get('library_evidence_table_sha256') if reviewed else None}}
    document.update(build_observation_report(raw, candidates, top, source_available=source.is_file(),
        analysis_complete=complete, classification=code, legacy_result=label,
        agreement=agreement, adequacy_status=adequacy['status']))
    document['family_observation'] = build_family_observation_report(
        raw, candidates, expected_observations, assignment,
        source_available=source.is_file(), analysis_complete=complete)
    document['reporting_features'].update(
        exact_and_nearby_family_observation=1,
        family_observation_changes_classification_or_adequacy=False)
    document['evidence_axes'] = build_evidence_axes(document, global_catalogue=global_references)
    document['reporting_features'].update(parallel_evidence_axes=1,
        parallel_evidence_axes_change_classification_or_adequacy=False)
    structure = saved_repeat_structure.build_evidence(sample, document, raw,
        saved_repeat_structure.load_evidence(getattr(args, 'saved_structure_evidence', None)))
    document['saved_repeat_structure'] = structure
    document['reporting_features'].update(saved_repeat_structure=1,
        saved_repeat_structure_changes_classification_or_adequacy=False)
    document['assembly_context_evidence'] = assembly_context_evidence.build_evidence(
        sample, document, assembly_context_evidence.load_evidence(
            getattr(args, 'assembly_context_evidence', None)))
    document['reporting_features'].update(assembly_context_evidence=1,
        assembly_context_changes_classification_or_adequacy=False)
    document['quality_path_evidence'] = quality_path_evidence.build_evidence(
        sample, document, quality_path_evidence.load_evidence(
            getattr(args, 'quality_path_evidence', None)))
    document['reporting_features'].update(quality_path_evidence=1,
        quality_path_changes_classification_or_adequacy=False)
    document['long_read_context_evidence'] = long_read_context_evidence.build_evidence(
        sample, document, long_read_context_evidence.load_evidence(
            getattr(args, 'long_read_context_evidence', None)))
    document['reporting_features'].update(long_read_context_evidence=1,
        long_read_context_changes_classification_or_adequacy=False)
    # Length is a necessary condition, separate from observations and power.
    applicability = build_input_applicability(power_input, expected_motifs,
        analysis_complete=complete, source_available=source.is_file())
    applicability['provenance'] = {
        'configuration_source': 'original_discovery', 'configuration_result_sha256': source_hash,
        'qc_source': 'full_input_evidence' if targeted.get('evidence_qc') else 'original_summary',
        'qc_result_sha256': evidence_hash if targeted.get('evidence_qc') else source_hash,
        'native_model_metadata': native_model_metadata,
        'applies_to': 'original_discovery_settings; does not reinterpret differently configured targeted measurements'}
    document['input_applicability'] = applicability
    document['reporting_features'].update(input_length_applicability=1,
        input_length_applicability_changes_classification_or_adequacy=False,
        flank_period_length_applicability=1,
        flank_period_length_applicability_changes_classification_or_adequacy=False)
    # Add explanatory evidence only after all scientific decisions have been made.
    technical_path = getattr(args, 'technical_boundary_evidence', None)
    technical = build_technical_boundary_evidence(sample, document,
        load_evidence(technical_path) if technical_path else None)
    document['technical_boundary_diagnostics'] = technical
    document['reporting_features'].update(technical_boundary_diagnostics=1,
        technical_boundary_changes_classification_or_adequacy=False)
    localization_path = getattr(args, 'reference_localization_evidence', None)
    localization = localization_evidence.build_reference_localization_evidence(sample, document,
        localization_evidence.load_evidence(localization_path) if localization_path else None)
    document['reference_localization_diagnostics'] = localization
    document['reporting_features'].update(reference_localization_diagnostics=1,
        reference_localization_changes_classification_or_adequacy=False)
    cross_period_path = getattr(args, 'cross_period_evidence', None)
    periodic = cross_period_evidence.build_cross_period_evidence(sample, document,
        cross_period_evidence.load_evidence(cross_period_path) if cross_period_path else None)
    document['cross_period_evidence'] = periodic
    document['reporting_features'].update(cross_period_diagnostics=1,
        cross_period_changes_classification_or_adequacy=False)
    followup = load_followup_diagnostics(sample, meta, raw, source_hash,
        getattr(args, 'followup_diagnostics', None), analysis_complete=complete)
    document['followup_diagnostics'] = followup
    document['reporting_features'].update(
        selected_interval_followup=1, followup_changes_classification_or_adequacy=False)
    clade = build_clade_variant_report(candidates, meta, assignment,
        getattr(args, 'taxon_reference_rows', []), global_catalogue=global_references,
        representative=top, analysis_complete=complete, library=library)
    clade_labels = {
        'clade_canonical_and_variant_observed': '분류군 기본 반복과 다른 variant가 함께 관측됨',
        'clade_variant_observed': '분류군 기본 반복과 다른 variant 관측',
        'clade_canonical_observed': '분류군 기본 반복 관측',
        'clade_variant_trace_only': '분류군 variant의 소수 서열 관측',
        'library_context_limited': '라이브러리 제약으로 분류군 해석 보류',
        'taxonomic_scope_unresolved': '분류군 정보 미확정',
        'analysis_incomplete': '분석 미완료로 분류군 비교 보류',
        'specific_reference_blocks_clade_interpretation': '종/하위 분류군의 문헌 구조에 따라 기본 반복 비교 보류',
        'nearest_clade_has_no_supported_default': '해당 분류군의 단일 기본 반복 근거 없음',
        'no_curated_clade_default': '비교할 분류군 기본 반복 미등록',
        'measurement_disagreement': '분류군 비교 전 측정 불일치 확인 필요',
        'no_relevant_family_observation': '분류군 비교에 사용할 반복 관측 미보고'}
    clade['display_label'] = clade_labels.get(clade['status'], clade['status'])
    document['clade_variant_observation'] = clade
    document['reporting_features']['clade_variant_observation'] = 1
    if clade['review_required'] and action == 'canonical_screen_complete':
        document['next_action'] = 'review_clade_variant'
        document['recommendation'] = ('전역적으로 알려진 반복이지만 분류군 기본 motif와 다른 family가 관측됐습니다. '
            '종별로 이미 알려진 예외인지 문헌을 확인하고, 새로운 생물학적 주장을 하려면 물리적 말단 위치와 독립 표본을 검증하세요.')
    direction = top['directionality'] if top else {}
    row = {'sample_id': sample, 'scientific_name': meta.get('scientific_name', ''),
           'result': label, 'detected_motif': top['motif'] if top else '',
           'sequence_observation': OBSERVATION_LABELS.get(observation['status'], observation['status']),
           'expected_repeat_observation': expected_observation_text,
           'expected_motifs': ','.join(expected_motifs),
           'reference_status': assignment.get('reference_status', ''),
           'repeat_architecture': assignment.get('repeat_architecture', ''),
           'reference_context': document['reference_context'],
           'library_context': ('등록 분류군–서열 불일치: TRM의 종별 귀속 보류' if library.get('identity_conflict_supported') else
                               '라이브러리 해석 제약: '+', '.join(library['reasons']) if library['limited'] else
                               '확인된 해석 제약 없음·순도/충분성 미확정'),
           'why': (document['reference_context'] if code == 'known_alternative_architecture' else
                   '라이브러리 조건 확인 필요: '+', '.join(library['reasons']) if code == 'library_context_limited' else
                   EVIDENCE_LABELS.get(top['candidate_status'], top['candidate_status']) if top else label),
           'geometry_nominations': sum(c['geometry_nomination'] for c in candidates),
           'provisional_candidates': sum(c['novel_candidate'] for c in candidates),
           **validation_columns(top),
           'evidence_level': EVIDENCE_LABELS.get(top['candidate_status'], top['candidate_status']) if top else EVIDENCE_LABELS['no_supported_call'],
           'supporting_reads': top['reads'] if top else '',
           'full_repeat_reads': top['full_repeat_reads'] if top else '',
           'left_consistent_junctions': direction.get('left_consistent'),
           'right_consistent_junctions': direction.get('right_consistent'),
           'opposite_junctions': direction.get('opposite_junctions'),
           'direction_consistency': round(direction['consistent_fraction'],4) if direction.get('consistent_fraction') is not None else None,
           'data_adequacy': ADEQUACY_LABELS[adequacy['status']], 'records': raw['summary'].get('records'),
           'bases': raw['summary'].get('bases'), 'next_action': recommendations[action],
           'unmeasured_novel_families': document['target_selection'].get('novel_candidates_not_targeted', 0),
           'paired_consistency': paired_summary(top),
           'trm_interpretation': agreement_label, 'reference_agreement_status': agreement['status'],
           'reference_conflict_review': agreement['requires_conflict_review'],
           'main_trm_status': agreement['main_trm_status'],
           'reference_review_recommendation': review_recommendation,
           'sra_search': sra_link(meta), 'details': f'samples/{sample}.html'}
    row['next_action'] = document['recommendation']
    row.update(quality_observation_status=measured_quality['status'],
        quality_observation_display={'not_measured': 'Quality 미측정', 'measured': '명시한 Phred+33 quality 관측',
            'partially_measured': '일부 read의 quality 관측·결측 별도 보존',
            'unavailable_in_analyzed_records': '입력에 측정 가능한 quality 없음'}[measured_quality['status']])
    row.update(clade_variant_status=clade['status'], clade_variant_display=clade['display_label'],
        clade_variant_families=clade['variant_observed_family_count'],
        clade_variant_review_required=clade['review_required'],
        published_clade_exception_families=clade['published_exception_family_count'])
    row.update(technical_boundary_status=technical['status'],
        technical_boundary_display=technical['display_label'],
        technical_boundary_scope=technical['scope'],
        technical_boundary_summary=technical['summary'],
        technical_boundary_recommendation=technical['recommendation'],
        technical_boundary_selected_events=technical['selected_event_count'],
        technical_boundary_flagged_events=technical['flagged_event_count'])
    row.update(reference_localization_status=localization['status'],
        reference_localization_display=localization['display_label'],
        reference_localization_scope=localization['scope'],
        reference_localization_summary=localization['summary'])
    row.update({'reference_localization_' + key: localization[key] for key in LOCALIZATION_COUNTS})
    row.update(cross_period_status=periodic['status'],
        cross_period_display=periodic['display_label'], cross_period_scope=periodic['scope'],
        cross_period_summary=periodic['summary'])
    row.update({'cross_period_' + key: periodic[key] for key in CROSS_PERIOD_COUNTS})
    row.update(input_applicability_status=applicability['status'],
        input_applicability_display=applicability['display_label'],
        input_applicability_scope=applicability['scope'],
        input_applicability_recommendation=applicability['recommendation'])
    flank_applicability = applicability['flank_period_applicability']
    row.update(flank_period_applicability_status=flank_applicability['status'],
        flank_period_applicability_display=flank_applicability['display_label'],
        flank_period_applicability_scope=flank_applicability['scope'],
        flank_period_applicability_recommendation=flank_applicability['recommendation'])
    row.update(document['evidence_axes']['display'])
    row.update(document['saved_repeat_structure']['display'])
    row.update(document['assembly_context_evidence']['display'])
    row.update(document['quality_path_evidence']['display'])
    row.update(document['long_read_context_evidence']['display'])
    family_view = document['family_observation']
    row.update(family_observation_status=family_view['status'],
        family_observation_display=family_view['display_label'],
        family_observation_summary='; '.join(
            (family['canonical_family'] or ','.join(str(m) for m in family.get('expected_motifs', [])
                if m is not None) or '해석 불가 단위')+': '+family['display_label']
            for family in family_view.get('families', [])))
    row.update(followup_diagnostic_status=followup['status'],
        followup_diagnostic_display=followup['display'],
        followup_diagnostic_scope=followup['scope'],
        followup_diagnostic_summary=followup['summary'],
        followup_diagnostic_recommendations='; '.join(followup['recommendations']))
    observed = document['repeat_observation']
    native_leader = observed['native_discovery_leader'] or {}
    policy_leader = observed['policy_score_leader'] or {}
    row.update(observation_status=observed['display_label'], measurement_status=observed['measurement_label'],
        decision_display=document['decision_summary']['display_label'],
        native_discovery_leader=native_leader.get('canonical_family'), native_discovery_reads=native_leader.get('reads'),
        policy_score_leader=policy_leader.get('motif'), policy_score_leader_reads=policy_leader.get('reads'),
        policy_representative=top['motif'] if top else '미지정',
        main_trm_display=document['decision_summary']['main_trm_display'],
        positive_repeat_families=observed['positive_family_count'], reporting_version=REPORTING_VERSION)
    for field in ('native_vs_policy_score_same_family', 'native_vs_representative_same_family',
                  'policy_score_vs_representative_same_family'):
        row[field] = document['policy_selection'][field]
    execution['expected_diagnostic_wall_seconds'] = expected_diagnostic.get('wall_seconds') if expected_diagnostic else None
    execution['expected_diagnostic_max_rss_kb'] = expected_diagnostic.get('max_rss_kb') if expected_diagnostic else None
    previous = raw.get('candidates', [{}])[0] if raw.get('candidates') else {}
    previous_key = previous.get('motif', '')
    prior_refs = references.get(previous_key, [])
    previous_display = prior_refs[0]['motif'] if prior_refs else previous_key
    updated_key = top['canonical_family'] if top else ''
    changed = previous_key != updated_key
    comparison = {'sample_id': sample, 'scientific_name': row['scientific_name'],
        'previous_call': raw.get('summary', {}).get('call', 'source_failed'),
        'previous_canonical_key': previous_key,
        'previous_reference_anchor': previous.get('telomere_prior', {}).get('matched_reference', ''),
        'previous_motif_display': previous_display, 'updated_result': label,
        'updated_motif_display': row['detected_motif'], 'updated_canonical_key': updated_key,
        'motif_family_changed': changed if source.is_file() and complete else None,
        'interpretation_changed': True,
        'reason': ('Input/source analysis incomplete; no biological motif comparison' if not complete else
                   'Neither result assigns a supported primary; nondetection remains unresolved' if not previous_key and not updated_key else
                   'Updated evidence requirements / control treated as annotation / direct reference coverage' if changed else
                   'Same family; conventional motif display, exact references and junction/adequacy interpretation added'),
        'previous_result': str(source.resolve()), 'updated_result_file': str((args.output/f'samples/{sample}.json').resolve())}
    return document, row, comparison


def sample_artifacts(document, row, comparison, args, directory):
    """Materialize one complete sample; detailed candidate rows never cross IPC."""
    from report_runtime import atomic_json
    sample = document['sample_id']
    meta = document['metadata']
    directory = Path(directory)
    write_json(directory/'assessment.json', document, compact=args.json_format == 'compact')
    from summarize_interpretation import summarize_sample, VERSION as SUMMARY_VERSION
    from summary_projection import write_projection
    write_projection(document, directory/'assessment.json',
                     (args.output/'samples'/f'{sample}.json').resolve(), directory/'summary_projection.json',
                     summarize_sample=summarize_sample, version=SUMMARY_VERSION)
    write_sample_html(directory/'assessment.html', document)
    chunks = []
    def table(name, rows, fields):
        write_tsv(directory/name, rows, fields)
        chunks.append(name)
    table('input_applicability.tsv', input_applicability_rows(document), INPUT_APPLICABILITY_FIELDS)
    table('flank_period_applicability.tsv', flank_applicability_rows(document), FLANK_APPLICABILITY_FIELDS)
    table('technical_boundary_diagnostics.tsv', technical_boundary_rows(document), TECHNICAL_BOUNDARY_FIELDS)
    table('technical_boundary_events.tsv', technical_event_rows(document), TECHNICAL_EVENT_FIELDS)
    table('reference_localization_diagnostics.tsv', reference_localization_rows(document), REFERENCE_LOCALIZATION_FIELDS)
    for name, key, fields in (
        ('reference_localization_families.tsv', 'family_summaries', localization_evidence.FAMILY_SUMMARY_FIELDS),
        ('reference_localization_references.tsv', 'reference_summaries', localization_evidence.REFERENCE_SUMMARY_FIELDS),
        ('reference_localization_events.tsv', 'events', localization_evidence.EVENT_FIELDS)):
        table(name, reference_localization_rows(document, key),
              list(dict.fromkeys(['sample_id', 'scientific_name', *fields])))
    table('cross_period_diagnostics.tsv', cross_period_rows(document), CROSS_PERIOD_FIELDS)
    for name, key, fields in (
        ('cross_period_families.tsv', 'family_summaries', cross_period_evidence.FAMILY_FIELDS),
        ('cross_period_intervals.tsv', 'intervals', cross_period_evidence.INTERVAL_FIELDS),
        ('cross_period_events.tsv', 'events', cross_period_evidence.EVENT_FIELDS),
        ('cross_period_contexts.tsv', 'contexts', cross_period_evidence.CONTEXT_FIELDS),
        ('cross_period_positions.tsv', 'positions', cross_period_evidence.POSITION_FIELDS)):
        table(name, cross_period_rows(document, key),
              list(dict.fromkeys(['sample_id', 'scientific_name', *fields])))
    table('assembly_context_summary.tsv', assembly_context_rows(document), ASSEMBLY_CONTEXT_FIELDS)
    for name, key, fields in (
        ('assembly_context_scaffolds.tsv', 'scaffolds', assembly_context_evidence.SCAFFOLD_FIELDS),
        ('assembly_context_intervals.tsv', 'intervals', assembly_context_evidence.INTERVAL_FIELDS),
        ('assembly_context_events.tsv', 'events', assembly_context_evidence.EVENT_FIELDS)):
        table(name, assembly_context_rows(document, key),
              list(dict.fromkeys(['sample_id', 'scientific_name', *fields])))
    table('quality_path_summary.tsv', quality_path_rows(document), QUALITY_PATH_FIELDS)
    for name, fields, rows in long_read_context_report.tables(document['long_read_context_evidence'],
            sample, document.get('scientific_name')):
        table(name, rows, fields)
    for name, key, fields in (
        ('quality_path_models.tsv', 'models', quality_path_evidence.MODEL_FIELDS),
        ('quality_path_read_comparisons.tsv', 'read_comparisons', quality_path_evidence.READ_FIELDS)):
        table(name, quality_path_rows(document, key),
              list(dict.fromkeys(['sample_id', 'scientific_name', *fields])))
    expected_rows, nearby_rows = family_observation_tables(
        sample, meta.get('scientific_name', ''), document['family_observation'])
    table('expected_family_observations.tsv', expected_rows, FAMILY_FIELDS)
    table('nearby_families.tsv', nearby_rows, NEARBY_FIELDS)
    followup = document['followup_diagnostics']
    followup_rows = ({'sample_id': sample, 'scientific_name': meta.get('scientific_name', ''),
        'status': followup['status'], 'display': followup['display'], **metric}
        for metric in followup.get('rows') or [{'metric': 'diagnostic_status',
            'label': followup['display'], 'value': None, 'unit': '', 'scope': followup['scope']}])
    table('followup_diagnostics.tsv', followup_rows, FOLLOWUP_FIELDS)
    table('runtime.tsv', [document['execution']], list(document['execution']))
    adequacy = document['data_adequacy']
    adequacy_row = {'sample_id': sample, 'status': adequacy['status'],
        'description': ADEQUACY_LABELS[adequacy['status']],
        'reasons': ';'.join(adequacy.get('reasons', [])),
        'missing_assumptions': ';'.join(adequacy.get('missing_assumptions', [])),
        'library_model_blockers': ';'.join(adequacy.get('library', {}).get('uniform_wgs_model_blockers', [])),
        'wga_mda_status': meta.get('wga_mda_status', 'not_queried'),
        'protocol_available': bool(meta.get('library_construction_protocol')),
        'conditional_power_calculation_available': adequacy.get('power') is not None,
        'absence_inferred': False, 'alt_inferred': False}
    table('adequacy.tsv', [adequacy_row], list(adequacy_row))
    scenarios = [{'sample_id': sample, 'illustrative_only': True, **scenario}
        for scenario in adequacy.get('illustrative_grid', {}).get('scenarios', [])]
    if scenarios:
        table('illustrative_power_scenarios.tsv', scenarios, list(scenarios[0]))
    def projected_candidates():
        for candidate in document['candidates']:
            d = candidate['directionality']
            candidate_row = {**candidate, 'sample_id': sample,
                **validation_columns(candidate), **paired_columns(candidate),
                'left_consistent_junctions': d['left_consistent'], 'right_consistent_junctions': d['right_consistent'],
                'opposite_junctions': d['opposite_junctions'], 'direction_consistency': d['consistent_fraction'],
                'direction_status': d['status'], 'orientation_basis': d['orientation_basis'],
                'reference_sources': ';'.join(candidate['reference_sources']), 'notes': ' '.join(candidate['notes'])}
            yield {key: candidate_row.get(key) for key in CANDIDATE_FIELDS}
    table('candidates.tsv', projected_candidates(), CANDIDATE_FIELDS)
    table('clade_variant_observations.tsv', (
        {key: ({'sample_id': sample, 'scientific_name': meta.get('scientific_name', ''), **item}).get(key)
         for key in CLADE_VARIANT_FIELDS}
        for item in document['clade_variant_observation']['families']), CLADE_VARIANT_FIELDS)
    def quality_rows():
        for candidate in document['candidates']:
            observed_quality = candidate.get('quality_diagnostics')
            if observed_quality is None:
                continue
            diagnostic = observed_quality['diagnostics']
            yield {'sample_id': sample, 'canonical_family': candidate['canonical_family'],
                'motif': candidate['motif'], 'status': observed_quality['status'],
                'encoding': observed_quality['encoding'],
                'supporting_records': diagnostic['supporting_records']['observed_records'],
                'quality_records': diagnostic['supporting_records']['quality_records'],
                'union_tract_bases': diagnostic['union_tract_bases'],
                'measured_tract_quality_bases': diagnostic['union_tract']['bases'],
                'tract_bases_without_quality': diagnostic['union_tract_bases_without_quality'],
                'tract_bases_below_q20': diagnostic['union_tract']['below_q20_bases'],
                'tract_bases_below_q30': diagnostic['union_tract']['below_q30_bases'],
                'source_sha256': observed_quality['source_sha256'],
                'variant_discriminating_bases_assessed': False, 'filtering_applied': False}
    table('candidate_quality_observations.tsv', quality_rows(), QUALITY_FIELDS)
    bundle = {'sample_id': sample, 'summary': row, 'comparison': comparison,
        'assessed_comparison': compare_previous_assessment(document, args.previous_assessment, args.output)
            if args.previous_assessment else None,
        'classification': document['classification'], 'followup_status': followup['status'],
        'family_status': document['family_observation']['status'], 'adequacy_status': adequacy['status'],
        'applicability_status': document['input_applicability']['status'],
        'flank_applicability_status': document['input_applicability']['flank_period_applicability']['status'],
        'technical_boundary_status': document['technical_boundary_diagnostics']['status'],
        'reference_localization_status': document['reference_localization_diagnostics']['status'],
        'cross_period_status': document['cross_period_evidence']['status'],
        'assembly_context_status': document['assembly_context_evidence']['status'],
        'quality_path_status': document['quality_path_evidence']['status'],
        'long_read_context_status': document['long_read_context_evidence']['status'],
        'chunks': chunks,
        'counts': {'candidate_rows': len(document['candidates']), 'family_rows': len(expected_rows),
                   'nearby_rows': len(nearby_rows),
                   'novel_candidate_rows': sum(c['novel_candidate'] for c in document['candidates']),
                   'geometry_nomination_rows': sum(c['geometry_nomination'] for c in document['candidates'])}}
    atomic_json(directory/'bundle.json', bundle)
    return bundle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True, help='Sample directories containing the native result.json files')
    parser.add_argument('--evidence', type=Path, help='Optional complete targeted measurements; omit when the native scan already includes junction evidence and QC')
    parser.add_argument('--expected-evidence', type=Path,
                        help='Optional full-input targeted diagnostics for the expectation table only; preserves native candidates, ranking and adequacy')
    parser.add_argument('--followup-diagnostics', type=Path,
                        help='Optional reviewed selected-interval diagnostics, joined by exact native input provenance; adds explanations without changing TRM decisions or adequacy')
    parser.add_argument('--technical-boundary-evidence', type=Path,
                        help='Optional sealed selected-junction technical-sequence evidence JSON; adds warnings without filtering motifs or changing TRM decisions')
    parser.add_argument('--reference-localization-evidence', type=Path,
                        help='Optional sealed saved-junction reference review; displays localization evidence without changing native counts, TRM decisions or adequacy')
    parser.add_argument('--saved-structure-evidence', type=Path,
        help='Optional sealed manifest linking saved full-input quality/mixed native results; annotation only')
    parser.add_argument('--quality-path-evidence', type=Path,
        help='Optional sealed same-read quality and saved all-optimal-path review; annotation only, no native count/classification or error-calibration changes')
    parser.add_argument('--long-read-context-evidence', type=Path,
        help='Optional sealed same-run long-read technical/structure annotations with native SHA binding; targeted studies without a native result use report_long_read_context.py')
    parser.add_argument('--assembly-context-evidence', type=Path,
        help='Optional sealed saved assembly-context review; separates read observations, assembly extension and physical-terminal confirmation without changing TRM decisions')
    parser.add_argument('--cross-period-evidence', type=Path,
                        help='Optional sealed saved-interval cross-period diagnostic JSON; preserves native counts, classification and adequacy')
    parser.add_argument('--metadata', type=Path, required=True, help='Metadata TSV for the samples in this analysis')
    parser.add_argument('--library-context', type=Path)
    parser.add_argument('--library-evidence', type=Path, help='Reviewed library identity TSV with sample/taxon identities and local evidence provenance')
    parser.add_argument('--previous-assessment', type=Path, help='Prior assessed JSON directory (samples/ID.json or ID.json); emits genuine version comparison separately from raw comparison')
    parser.add_argument('--library-overrides', type=Path, help='Separately sourced literature annotations; archive metadata retained')
    parser.add_argument('--published-context', type=Path, help='Prior published analysis of exact accessions; annotation only')
    parser.add_argument('--assignments', type=Path, required=True)
    parser.add_argument('--catalog', type=Path, default=Path('data/telomere_motifs.tsv'))
    parser.add_argument('--taxon-references', type=Path, default=Path('data/taxon_motif_references.tsv'))
    parser.add_argument('--assumptions', type=Path)
    parser.add_argument('--duplicate-audit', type=Path, help='Optional explicit archive-duplicate audit for these samples; no historical local audit is loaded by default')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--comparison-output', type=Path, required=True)
    parser.add_argument('--json-format', choices=('pretty', 'compact'), default='pretty',
                        help='JSON whitespace only; compact retains every field and uses faster one-shot encoding. Default: pretty')
    parser.add_argument('--jobs', type=int, default=1, help='Concurrent sample interpretation processes; default 1')
    parser.add_argument('--resume', action='store_true', help='Reuse only completed sample artifacts with matching inputs and checksums')
    parser.add_argument('--memory-limit-mib', type=int, help='Hard address-space limit per interpretation worker; unset means inherited OS limit')
    parser.add_argument('--memory-budget-mib', type=int, help='Require jobs times per-worker cap to fit this worker budget; coordinator memory is additional')
    args = parser.parse_args()
    from report_runtime import validate_limits
    try:
        validate_limits(args)
    except ValueError as error:
        parser.error(str(error))
    def overlap(left, right):
        left, right = left.resolve(), right.resolve()
        return left == right or left in right.parents or right in left.parents
    outputs = (args.output, args.comparison_output)
    input_directories = ([args.source] + ([args.evidence] if args.evidence else [])
                         + ([args.expected_evidence] if args.expected_evidence else [])
                         + ([args.followup_diagnostics] if args.followup_diagnostics else [])
                         + ([args.previous_assessment] if args.previous_assessment else []))
    if overlap(*outputs) or any(overlap(output, source) for output in outputs for source in input_directories):
        parser.error('Current, comparison and input directories must be separate and not nested')
    input_files = [getattr(args, name) for name in ('metadata', 'assignments', 'catalog',
        'taxon_references', 'library_context', 'library_evidence', 'library_overrides',
        'published_context', 'assumptions', 'duplicate_audit') if getattr(args, name)]
    if args.technical_boundary_evidence:
        from technical_boundary_evidence import source_dependencies
        input_files.extend(Path(path) for path in source_dependencies(args.technical_boundary_evidence))
    if args.reference_localization_evidence:
        input_files.extend(Path(path) for path in localization_evidence.source_dependencies(args.reference_localization_evidence))
    if args.saved_structure_evidence:
        input_files.extend(Path(path) for path in saved_repeat_structure.source_dependencies(args.saved_structure_evidence))
    if args.assembly_context_evidence:
        input_files.extend(Path(path) for path in assembly_context_evidence.source_dependencies(args.assembly_context_evidence))
    if args.quality_path_evidence:
        input_files.extend(Path(path) for path in quality_path_evidence.source_dependencies(args.quality_path_evidence))
    if args.long_read_context_evidence:
        input_files.extend(Path(path) for path in long_read_context_evidence.source_dependencies(args.long_read_context_evidence))
    if args.cross_period_evidence:
        input_files.extend(Path(path) for path in cross_period_evidence.source_dependencies(args.cross_period_evidence))
    if any(output.resolve() == source.resolve() or output.resolve() in source.resolve().parents
           for output in outputs for source in input_files):
        parser.error('Output directories must not contain input provenance files')
    from report_runtime import output_locks, atomic_json, utc
    with output_locks(args.output, args.comparison_output) as lock_fds:
        import uuid
        args._runtime_lock_fds = lock_fds
        args._execution_id = uuid.uuid4().hex
        execution_directory = args.output/'.report-executions'/args._execution_id
        invocation = {'execution_id': args._execution_id, 'status': 'preflight', 'started_utc': utc()}
        atomic_json(args.output/'interpretation_status.json', invocation)
        atomic_json(execution_directory/'status.json', invocation)
        atomic_json(execution_directory/'invocation.json', {
            key: str(value.resolve()) if isinstance(value, Path) else value
            for key, value in vars(args).items() if not key.startswith('_')})
        try:
            execute_report(args, parser)
        except BaseException as error:
            status_path = args.output/'interpretation_status.json'
            status = json.loads(status_path.read_text()) if status_path.is_file() else invocation
            status.update(status='failed', finished_utc=utc(), error=str(error),
                note='Incomplete aggregation is not a completed report. Prior metrics/index may describe an earlier successful execution; verified sample caches may be resumed.')
            atomic_json(status_path, status)
            atomic_json(execution_directory/'status.json', status)
            raise


def execute_report(args, parser):
    args.output.mkdir(parents=True, exist_ok=True)
    args.comparison_output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    from report_runtime import prepare_globals
    args.runtime_globals = prepare_globals(args)
    metadata = index_samples(read_tsv(args.metadata), 'metadata')
    assignments = index_samples(read_tsv(args.assignments), 'assignments')
    if metadata.keys() != assignments.keys():
        raise ValueError('Metadata and assignment sample sets differ')
    if args.library_context:
        for row in read_tsv(args.library_context):
            if row['sample_id'] in metadata:
                for key, value in row.items():
                    if not value or key in {'sample_id', 'scientific_name', 'tax_id', 'lineage', 'sample_accession', 'study_accession'}:
                        continue
                    output_key = 'library_context_'+key if key in {'status', 'error', 'source_url', 'response_sha256'} else key
                    metadata[row['sample_id']][output_key] = value
    if args.library_overrides:
        for sample, row in index_samples(read_tsv(args.library_overrides), 'library literature annotations').items():
            if sample not in metadata:
                raise ValueError(f'Library literature annotation sample absent from metadata: {sample}')
            if not row.get('library_evidence_source', '').startswith('http'):
                raise ValueError(f'{sample}: library literature annotation requires source URL')
            metadata[sample]['library_literature_annotation'] = row
            for key in ('whole_genome_amplification', 'amplification_method'):
                if row.get(key):
                    metadata[sample]['archive_'+key] = metadata[sample].get(key, '')
                    metadata[sample][key] = row[key]
    if args.published_context:
        for sample, row in index_samples(read_tsv(args.published_context), 'published sample context').items():
            if sample not in metadata:
                raise ValueError(f'Published context sample absent from metadata: {sample}')
            if not row.get('source', '').startswith(('https://', 'http://')):
                raise ValueError(f'{sample}: published context requires source URL')
            metadata[sample]['published_sample_context'] = row
    try:
        args.library_evidence_rows = load_library_evidence(args.library_evidence, metadata)
    except (ValueError, OSError) as error:
        raise ValueError(str(error))
    references = catalog_index(read_tsv(args.catalog))
    args.taxon_reference_rows = read_tsv(args.taxon_references)
    defaults, per_sample = read_assumptions(args.assumptions)
    duplicates = {}
    if args.duplicate_audit and args.duplicate_audit.is_file():
        audit = json.loads(args.duplicate_audit.read_text())
        for group in audit.get('archive_duplicate_groups', []):
            runs = group.get('accessions', group.get('runs', group.get('sample_ids', [])))
            if runs and isinstance(runs[0], dict):
                runs = [r.get('sample_id', r.get('run_accession')) for r in runs]
            for run in runs:
                duplicates[run] = runs
    from report_runtime import run_samples, merge_chunks, read_bundle, validate_receipts, assert_unchanged
    receipts, runtime_status, execution_directory = run_samples(args, metadata, assignments,
        references, defaults, per_sample, duplicates, args.runtime_globals)
    summaries, comparisons, assessed_comparisons = [], [], []
    codes, followup_counts, family_counts, adequacy_counts, applicability_counts, flank_applicability_counts, technical_counts, localization_counts, cross_period_counts = (
        Counter() for _ in range(9))
    table_counts = Counter()
    assembly_context_counts = Counter()
    quality_path_counts = Counter()
    long_read_context_counts = Counter()
    for receipt in receipts:
        bundle = read_bundle(receipt)
        summaries.append(bundle['summary'])
        comparisons.append(bundle['comparison'])
        if bundle['assessed_comparison'] is not None:
            assessed_comparisons.append(bundle['assessed_comparison'])
        codes[bundle['classification']] += 1
        followup_counts[bundle['followup_status']] += 1
        family_counts[bundle['family_status']] += 1
        adequacy_counts[bundle['adequacy_status']] += 1
        applicability_counts[bundle['applicability_status']] += 1
        flank_applicability_counts[bundle['flank_applicability_status']] += 1
        technical_counts[bundle['technical_boundary_status']] += 1
        localization_counts[bundle['reference_localization_status']] += 1
        cross_period_counts[bundle['cross_period_status']] += 1
        assembly_context_counts[bundle['assembly_context_status']] += 1
        quality_path_counts[bundle['quality_path_status']] += 1
        long_read_context_counts[bundle['long_read_context_status']] += 1
        table_counts.update(bundle['counts'])
    write_tsv(args.output/'samples.tsv', summaries, SAMPLE_FIELDS)
    write_tsv(args.output/'saved_repeat_structure.tsv', summaries,
              ['sample_id', 'scientific_name', *saved_repeat_structure.FIELDS, 'details'])
    write_tsv(args.output/'evidence_axes.tsv', summaries,
              ['sample_id', 'scientific_name', *AXIS_FIELDS, 'details'])
    from summary_projection import write_manifest as write_projection_manifest
    write_projection_manifest(receipts, args.output)
    for name in ('input_applicability.tsv', 'flank_period_applicability.tsv', 'followup_diagnostics.tsv', 'candidates.tsv',
                 'technical_boundary_diagnostics.tsv', 'technical_boundary_events.tsv',
                 *REFERENCE_LOCALIZATION_TABLES, *CROSS_PERIOD_TABLES, *ASSEMBLY_CONTEXT_TABLES, *QUALITY_PATH_TABLES, *LONG_READ_CONTEXT_TABLES,
                 'clade_variant_observations.tsv',
                 'candidate_quality_observations.tsv',
                 'expected_family_observations.tsv', 'nearby_families.tsv', 'runtime.tsv',
                 'adequacy.tsv', 'illustrative_power_scenarios.tsv'):
        merge_chunks(receipts, name, args.output / name)
    write_tsv(args.comparison_output/'comparison.tsv', comparisons, COMPARISON_FIELDS)
    if assessed_comparisons:
        write_tsv(args.comparison_output/'previous_assessment_comparison.tsv', assessed_comparisons, PREVIOUS_ASSESSMENT_FIELDS)
        write_html(args.comparison_output/'previous_assessment_comparison.html', assessed_comparisons,
            ['sample_id', 'scientific_name', 'previous_version', 'updated_version',
             'previous_classification', 'updated_classification', 'previous_motif', 'updated_motif',
             'previous_reference_agreement', 'updated_reference_agreement', 'reference_interpretation_added',
             'previous_schema_version', 'updated_schema_version',
             'previous_followup_status', 'updated_followup_status', 'followup_diagnostics_added',
             'classification_changed', 'motif_family_changed', 'confirmation_missing', 'reason'],
            '이전 assessed 결과 → v'+POLICY['version']+' assessed 결과',
            '이전 assessment JSON의 실제 classification·motif와 비교합니다. 원래 경계 후보는 보존하고 추가 검증 근거·구조/라이브러리 제약을 구분합니다. Short reads만으로 말단 위치를 확정하지 않습니다.')
    overview_fields = ['sample_id','scientific_name', *AXIS_FIELDS, *saved_repeat_structure.FIELDS,
                       *assembly_context_evidence.FIELDS, *quality_path_evidence.FIELDS, *long_read_context_evidence.FIELDS, 'measurement_status','observation_status',
                       'native_discovery_leader','native_discovery_reads','policy_score_leader','policy_score_leader_reads',
                       'policy_representative','decision_display','main_trm_display','trm_interpretation',
                       'expected_repeat_observation','family_observation_display',
                       'family_observation_summary','library_context','data_adequacy',
                       'clade_variant_display','clade_variant_review_required',
                       'quality_observation_display',
                       'technical_boundary_display','technical_boundary_summary',
                       'reference_localization_display','reference_localization_summary',
                       'cross_period_display','cross_period_summary',
                       'input_applicability_display','input_applicability_scope',
                       'flank_period_applicability_display','flank_period_applicability_scope','details']
    if any(row['paired_consistency'] != '미측정' for row in summaries):
        overview_fields.insert(-1, 'paired_consistency')
    if args.followup_diagnostics:
        overview_fields[-1:-1] = ['followup_diagnostic_display', 'followup_diagnostic_summary']
    write_html(args.output/'index.html', summaries, overview_fields, 'TRM 해석 v'+POLICY['version']+' 결과',
        '측정 상태·반복 관측·정책상 대표·TRM 판정·main TRM 확정을 구분합니다. Native 우선 반복은 참조 annotation 전의 '
        '저장된 탐색 순위이며 현재 해석 점수와 다릅니다. 높은 순위도 내부·위성 반복일 수 있습니다. '
        '대표 미지정은 서열 부재가 아니며, 미보고와 명시적 0건을 구분합니다. 호환용 result/detected_motif 필드는 TSV/JSON에 보존합니다.')
    write_html(args.comparison_output/'index.html', comparisons, COMPARISON_FIELDS, 'Raw scan → TRM 해석 v'+POLICY['version']+' 비교',
        '이 표는 raw scan과 assessment를 비교하며 이전 assessed 버전 비교가 아닙니다. --previous-assessment를 지정하면 previous_assessment_comparison.tsv/html에 실제 이전 판정과의 비교가 별도로 생성됩니다.',
        links=[('실제 이전 assessed 판정 비교', 'previous_assessment_comparison.html'),
               ('이전 판정 비교 TSV', 'previous_assessment_comparison.tsv')] if assessed_comparisons else [])
    metrics = {'version': POLICY['version'], 'reporting_version': REPORTING_VERSION,
               'json_format': args.json_format,
               'schema_version': 22, 'followup_diagnostic_status_counts': dict(followup_counts),
               'clade_variant_status_counts': dict(Counter(row['clade_variant_status'] for row in summaries)),
               'clade_variant_review_samples': sum(row['clade_variant_review_required'] for row in summaries),
               'technical_boundary_status_counts': dict(technical_counts),
               'reference_localization_status_counts': dict(localization_counts),
               'cross_period_status_counts': dict(cross_period_counts),
               'input_applicability_status_counts': dict(applicability_counts),
               'flank_period_applicability_status_counts': dict(flank_applicability_counts),
               'family_observation_status_counts': dict(family_counts),
               'expected_family_rows': table_counts['family_rows'], 'nearby_family_rows': table_counts['nearby_rows'],
               'created_utc': datetime.now(timezone.utc).isoformat(),
               'samples': len(summaries), 'sample_ids': sorted(metadata), 'candidate_rows': table_counts['candidate_rows'], 'classifications': dict(codes),
               'data_adequacy_counts': dict(adequacy_counts),
               'novel_candidate_rows_including_secondary': table_counts['novel_candidate_rows'],
               'geometry_nomination_rows_including_secondary': table_counts['geometry_nomination_rows'],
               'previous_assessment_comparisons': len(assessed_comparisons),
               'eligible_novel_families_without_new_junction_measurement': sum(r['unmeasured_novel_families'] for r in summaries),
               'top_family_changes': sum(c['motif_family_changed'] is True for c in comparisons),
               'report_wall_seconds': time.monotonic()-started,
               'files': {'samples.tsv': sha256(args.output/'samples.tsv'),
                         'evidence_axes.tsv': sha256(args.output/'evidence_axes.tsv'),
                         'saved_repeat_structure.tsv': sha256(args.output/'saved_repeat_structure.tsv'),
                         'summary_projections.tsv': sha256(args.output/'summary_projections.tsv'),
                         'clade_variant_observations.tsv': sha256(args.output/'clade_variant_observations.tsv'),
                         'candidate_quality_observations.tsv': sha256(args.output/'candidate_quality_observations.tsv'),
                         'input_applicability.tsv': sha256(args.output/'input_applicability.tsv'),
                         'flank_period_applicability.tsv': sha256(args.output/'flank_period_applicability.tsv'),
                         'technical_boundary_diagnostics.tsv': sha256(args.output/'technical_boundary_diagnostics.tsv'),
                         'technical_boundary_events.tsv': sha256(args.output/'technical_boundary_events.tsv'),
                         'candidates.tsv': sha256(args.output/'candidates.tsv'),
                         'followup_diagnostics.tsv': sha256(args.output/'followup_diagnostics.tsv'),
                         'expected_family_observations.tsv': sha256(args.output/'expected_family_observations.tsv'),
                         'nearby_families.tsv': sha256(args.output/'nearby_families.tsv'),
                         'comparison.tsv': sha256(args.comparison_output/'comparison.tsv')}}
    metrics['assembly_context_status_counts'] = dict(assembly_context_counts)
    metrics['quality_path_status_counts'] = dict(quality_path_counts)
    metrics['long_read_context_status_counts'] = dict(long_read_context_counts)
    metrics['files'].update({name: sha256(args.output/name) for name in (
        *REFERENCE_LOCALIZATION_TABLES, *CROSS_PERIOD_TABLES, *ASSEMBLY_CONTEXT_TABLES, *QUALITY_PATH_TABLES, *LONG_READ_CONTEXT_TABLES)})
    if assessed_comparisons:
        metrics['files']['previous_assessment_comparison.tsv'] = sha256(args.comparison_output/'previous_assessment_comparison.tsv')
    validate_receipts(receipts)
    assert_unchanged(args.runtime_globals['pins'])
    assert_unchanged(args.runtime_globals['code_pins'])
    write_json(args.output/'metrics.json', metrics)
    from report_runtime import atomic_json, utc
    runtime_status.update(status='complete', finished_utc=utc(), report_wall_seconds=time.monotonic()-started)
    atomic_json(args.output/'interpretation_status.json', runtime_status)
    atomic_json(execution_directory/'status.json', runtime_status)
    print(json.dumps(metrics, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 3 and sys.argv[1] == '--_worker-task':
        from report_runtime import worker
        worker(Path(sys.argv[2]), sys.modules[__name__])
    else:
        main()
