"""Length-only opportunity for the supported native candidate-period annotation.

These bounds count reads of sufficient length, not observed junctions, successful
annotations, independent molecules, or telomere detection power.
"""
from assess_adequacy import _read_qc, eligible_read_bounds


# These releases share the reviewed seed, disjoint flank and 90% full-tract
# contract. Accept releases explicitly; a future/custom build needs review.
SUPPORTED_NATIVE_VERSIONS = ('0.8.4', '0.8.5')


LABELS = {
    'no_flank_length_eligible_reads': '반복과 flank 검사 창을 함께 담을 길이의 리드 없음',
    'flank_length_eligible_reads_present_evaluation_unknown': '반복·flank 길이 필요조건 충족 리드 있음 · 실제 평가 여부 별도',
    'flank_length_eligibility_unresolved': '반복·flank 길이 필요조건 충족 여부 불명',
    'flank_model_unavailable': 'native flank 길이 모델의 적용 근거 부족',
    'junction_options_unavailable': '실제 junction 설정 정보 부족',
    'junction_options_invalid': '실제 junction 설정이 유효하지 않음',
}
LENGTH_STATES = tuple(LABELS)[:3]
LIMITATIONS = [
    '길이 bounds는 반복과 flank를 한 read에 담기 위한 필요조건만 평가하며 실제 경계 사건이나 평가 성공 수가 아닙니다.',
    '주석 on/off는 설정이며 candidate별 실제 주석의 유무·measured_all_events·사건 수와 구분합니다.',
    'ACGT 연속성, 실제 seed, 염기 품질, 다른 tract, flank 복잡도와 read 끝 위치는 검사하지 않습니다.',
    'Full-repeat 제외 조건은 최소 seed 길이에만 적용한 하한이며 실제 확장 tract가 길면 추가로 제외될 수 있습니다.',
    'Paired mate 길이는 합치지 않으며 주기별·family별 적격 read 수는 겹치므로 합산할 수 없습니다.',
    '길이 부족 또는 충족은 motif 미검출, TRM 부재·ALT·novel, 물리적 말단 또는 검출력 충분성의 판정이 아닙니다.',
    '대표 flank의80 bp 저장 상한은 모든 실제 사건에 적용하는 주석 창의 길이 한계가 아닙니다.',
    'Targeted 결과의 config min_k/max_k는 실제 target 목록이 아닙니다. 범위 밖 family도 표적 검사될 수 있습니다.',
]


def _status(lower, upper):
    if upper == 0:
        return LENGTH_STATES[0]
    return LENGTH_STATES[1] if lower > 0 else LENGTH_STATES[2]


def build_flank_applicability(report, seed_result):
    """Reuse validated seed/QC scope, adding a separate non-decision diagnostic.

    The native contract is limited to the reviewed trm v0.8.4/v0.8.5 releases.
    Unknown versions/options never acquire guessed defaults. This function does
    not inspect or mutate candidates, their annotations, or the seed diagnostic.
    """
    labels = dict(LABELS)
    # Propagate the existing input diagnostic's source/QC/completeness reasons.
    labels[seed_result['status']] = seed_result['display_label']
    result = dict(schema_version=1, status=None, display_label=None, scope='unavailable',
                  reason=None, recommendation=None, model_status='unavailable',
                  model='reviewed_native_disjoint_seed_and_candidate_period_window',
                  native_program=report.get('program') if isinstance(report, dict) else None,
                  native_version=report.get('version') if isinstance(report, dict) else None,
                  analysis_kind=report.get('analysis_kind') if isinstance(report, dict) else None,
                  annotation_enabled=None, actual_event_measurement_inferred=False,
                  junction_options=None, periods=[], families=[],
                  records=seed_result.get('records'), bases=seed_result.get('bases'),
                  max_record_bases=seed_result.get('max_record_bases'),
                  formula='max(d + max(W,F), floor(10*d/9)+1); d=k*max(min_copies,ceil(seed_bases/k)); W=max(20,3*k)',
                  limitations=list(LIMITATIONS))
    for name in ('periods', 'families'):
        for source in seed_result[name]:
            row = {key: source[key] for key in ('kind', 'period_bp', 'canonical_family',
                   'expected_motifs', 'outside_discovery_range') if key in source}
            row.update(minimum_seed_bp=source.get('minimum_seed_bp'), required_window_bp=None,
                       junction_flank_bp=None, full_repeat_minimum_read_bp=None,
                       minimum_read_bp=None, length_eligible_reads_lower=None,
                       length_eligible_reads_upper=None, bounds_method=None,
                       analyzed_records_length_status=None)
            result[name].append(row)

    def finish(status, scope, reason, recommendation):
        result.update(status=status, display_label=labels[status], scope=scope,
                      reason=reason, recommendation=recommendation)
        for row in result['periods'] + result['families']:
            effective = row['analyzed_records_length_status'] if scope == 'complete_input' else status
            row.update(scope=scope, status=effective, display_label=labels[effective])
        return result

    def unavailable(status, reason):
        return finish(status, 'unavailable', reason,
                      '원래 native 버전·junction 설정과 검증된 전체 입력의 길이 QC를 확인하세요. 미평가는0이 아닙니다.')

    if seed_result['scope'] == 'unavailable':
        return unavailable(seed_result['status'], seed_result['reason'])
    if (not isinstance(report, dict) or report.get('program') != 'trm' or
            report.get('version') not in SUPPORTED_NATIVE_VERSIONS):
        return unavailable('flank_model_unavailable', '검토된 trm v0.8.4/v0.8.5의 seed·flank·full-repeat 계약을 확인할 수 없습니다.')
    version = report['version']
    result['model'] = f'native_v{version}_disjoint_seed_and_candidate_period_window'
    options = report.get('junction_options')
    if options is None:
        return unavailable('junction_options_unavailable', 'junction_options가 없어 실제 flank 길이를 확인할 수 없습니다.')
    if not isinstance(options, dict):
        return unavailable('junction_options_invalid', 'junction_options는 객체여야 합니다.')
    if not all(key in options for key in ('flank_bases', 'end_slack')):
        return unavailable('junction_options_unavailable', 'flank_bases 또는 end_slack가 누락됐습니다.')
    flank, slack = options['flank_bases'], options['end_slack']
    enabled = options.get('flank_period_annotation', False)  # Both reviewed releases omit only false.
    if (type(flank) is not int or not 20 <= flank <= 10_000 or
            type(slack) is not int or not 0 <= slack <= 1000 or type(enabled) is not bool):
        return unavailable('junction_options_invalid', '실제 flank_bases/end_slack/on-off 설정의 형식 또는 범위가 유효하지 않습니다.')
    result.update(model_status=f'supported_native_v{version}', annotation_enabled=enabled,
                  junction_options=dict(flank_bases=flank, end_slack=slack, flank_period_annotation=enabled))
    # Seed diagnostics already validate QC (including bool/float lookalikes) and
    # keep analyzed-prefix bounds separate from whole-input conclusions.
    qc = _read_qc(report)
    for row in result['periods'] + result['families']:
        seed = row['minimum_seed_bp']
        window = max(20, 3 * row['period_bp'])
        full_minimum = (10 * seed) // 9 + 1
        minimum = max(seed + max(window, flank), full_minimum)
        bounds = eligible_read_bounds(qc['records'], qc['bases'], qc['max_record_bases'],
                                     minimum, qc['length_histogram'])
        row.update(required_window_bp=window, junction_flank_bp=flank,
                   full_repeat_minimum_read_bp=full_minimum, minimum_read_bp=minimum,
                   length_eligible_reads_lower=bounds['lower'], length_eligible_reads_upper=bounds['upper'],
                   bounds_method=bounds['method'],
                   analyzed_records_length_status=_status(bounds['lower'], bounds['upper']))
    if seed_result['scope'] != 'complete_input':
        return finish(seed_result['status'], seed_result['scope'], seed_result['reason'],
                      '현재 bounds는 보고된 records에만 해당합니다. 전체 입력 처리 완료와 길이 QC를 확인하세요.')
    states = [r['analyzed_records_length_status'] for r in result['periods']]
    if all(state == LENGTH_STATES[0] for state in states):
        return finish(LENGTH_STATES[0], 'complete_input',
                      '모든 설정 주기에서 seed와 flank 검사를 함께 담을 길이의 read 수 상한이0입니다.',
                      '더 긴 개별 read가 필요합니다. 기존 motif 관측과 native 주석의 실제 측정 상태는 별도로 유지하세요.')
    if LENGTH_STATES[1] in states:
        return finish(LENGTH_STATES[1], 'complete_input',
                      '일부 설정 주기에서 길이 필요조건을 충족합니다. 모든 주기나 실제 flank 평가 성공을 뜻하지 않습니다.',
                      '후보 period별 행과 실제 native flank 주석을 함께 확인하세요. 길이 충족만으로 말단·검출력을 확정하지 마세요.')
    return finish(LENGTH_STATES[2], 'complete_input',
                  '보존된 길이 분포 bounds만으로 반복과 flank를 함께 담을 read의 존재를 확정할 수 없습니다.',
                  '세밀한 read 길이 QC 또는 더 긴 입력을 확보하기 전에는 평가 성공·실패를 단정하지 마세요.')
