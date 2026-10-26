"""Length-only applicability of the recorded exact-seed discovery settings.

This diagnostic reads saved JSON only. It neither changes candidate selection nor
estimates sensitivity, telomere representation, or biological absence.
"""

from assess_adequacy import _read_qc, eligible_read_bounds, required_seed_bases
from reference_seed_opportunity import orbit, primitive


LABELS = {
    "no_length_eligible_reads": "최소 seed 길이를 만족하는 리드 없음",
    "length_eligible_reads_present_power_unknown": "길이 조건을 만족하는 리드 있음 · 검출력 미평가",
    "length_eligibility_unresolved": "길이 조건 충족 여부를 확정할 수 없음",
    "analysis_incomplete": "분석 미완료 · 처리된 리드만 평가",
    "input_completeness_unverified": "전체 입력 처리 여부 미확인",
    "qc_invalid": "리드 길이 QC 정보 누락 또는 불일치",
    "source_unavailable": "분석 원본 결과 없음",
    "mode_not_supported": "리드 길이 진단이 지원하지 않는 분석 모드",
    "config_unavailable": "실제 검색 설정 정보 부족",
    "config_invalid": "실제 검색 설정이 유효하지 않음",
}
LIMITATIONS = (
    "길이 조건은 필요조건이며 정확한 반복 seed, 염기 품질, 텔로미어 포함 여부를 확인하지 않습니다.",
    "길이 조건을 만족하는 리드가 있어도 검출력, 유전체 coverage, 독립 분자 수 또는 텔로미어 보존을 보장하지 않습니다.",
    "주기별·family별 리드 수는 서로 겹치므로 합산할 수 없으며 paired read도 개별 리드로 계산합니다.",
    "Family 계산은 primitive unit과 회전·역상보를 통합한 길이 진단이며 불규칙·혼합 배열을 고정 tandem repeat으로 증명하지 않습니다.",
    "검색 범위 밖 family의 길이 계산은 기록된 seed 설정을 적용한 참고값이며 별도 targeted 검색의 실행 여부나 가능 여부를 판정하지 않습니다.",
    "관측된 입력과 설정에서의 길이 부적합은 생물학적 TRM 부재, 신규 TRM 또는 ALT의 근거가 아닙니다.",
)
CONFIG_FIELDS = ("mode", "min_k", "max_k", "seed_bases", "min_copies")
FAILED_STATES = {"failed", "error", "truncated", "incomplete", "running", "pending"}


def _expected_families(expected_motifs):
    """Share the reference helper's primitive/rotation/RC convention, not seeds."""
    if isinstance(expected_motifs, str):
        expected_motifs = (expected_motifs,)
    grouped, invalid = {}, []
    try:
        iterator = iter(expected_motifs)
    except TypeError:
        return [], [{"input_index": None, "reason": "expected_motifs must be an iterable of DNA strings"}]
    for index, motif in enumerate(iterator):
        if not isinstance(motif, str) or not 1 <= len(motif) <= 64 or not set(motif.upper()) <= set("ACGT"):
            invalid.append({"input_index": index,
                            "expected_motif": motif if isinstance(motif, (str, int, float, bool, type(None))) else None,
                            "reason": "expected motif requires 1..64 literal A/C/G/T bases"})
            continue
        motif = motif.upper()
        unit = primitive(motif.encode("ascii"))
        family = orbit(unit)[0].decode("ascii")
        grouped.setdefault(family, set()).add(motif)
    return [{"canonical_family": family, "expected_motifs": sorted(motifs),
             "period_bp": len(family)} for family, motifs in sorted(grouped.items())], invalid


def _length_status(lower, upper):
    if upper == 0:
        return "no_length_eligible_reads"
    if lower > 0:
        return "length_eligible_reads_present_power_unknown"
    return "length_eligibility_unresolved"


def _finish(result, status, scope, reason, recommendation):
    result.update(status=status, display_label=LABELS[status], scope=scope,
                  reason=reason, recommendation=recommendation)
    for row in result["periods"] + result["families"]:
        row["scope"] = scope
        # Computed bounds may remain useful for an incomplete analyzed prefix,
        # but its displayed status must never assert whole-input suitability.
        effective = (row["analyzed_records_length_status"] if scope == "complete_input"
                     else status)
        row.update(status=effective, display_label=LABELS[effective])
    return result


def _build_seed_applicability(report, expected_motifs=(), *, analysis_complete=True,
                             source_available=True):
    """Describe necessary read-length conditions without a biological decision.

    ``complete_input`` requires both the caller's completion gate and an explicit
    false ``summary.sampled_prefix``. Legacy summaries without a histogram can
    supply conservative aggregate bounds. Missing settings never use defaults.
    Invalid expected references are listed separately and do not invalidate an
    otherwise usable diagnostic for the configured discovery periods.
    """
    families, invalid = _expected_families(expected_motifs)
    result = {"schema_version": 1, "kind": "input", "status": None,
              "display_label": None, "scope": "unavailable", "reason": None,
              "recommendation": None, "records": None, "bases": None,
              "max_record_bases": None, "config": None, "periods": [],
              "families": [], "invalid_expected_motifs": invalid,
              "limitations": list(LIMITATIONS)}
    for family in families:
        result["families"].append({**family, "kind": "expected_family",
                                   "outside_discovery_range": None,
                                   "minimum_seed_bp": None,
                                   "length_eligible_reads_lower": None,
                                   "length_eligible_reads_upper": None,
                                   "bounds_method": None,
                                   "analyzed_records_length_status": None})

    def unavailable(status, reason):
        return _finish(result, status, "unavailable", reason,
                       "검증된 전체 입력 결과와 실제 reads 모드 검색 설정·길이 QC를 확보한 뒤 다시 평가하세요.")

    if source_available is not True or not isinstance(report, dict):
        return unavailable("source_unavailable", "길이 진단에 사용할 분석 원본 결과가 없습니다.")
    if "config" not in report:
        return unavailable("config_unavailable", "config가 없어 실제 seed 조건을 알 수 없습니다.")
    config = report["config"]
    if not isinstance(config, dict):
        return unavailable("config_invalid", "config는 객체여야 합니다.")
    missing = [key for key in CONFIG_FIELDS if key not in config]
    if missing:
        return unavailable("config_unavailable", "필수 검색 설정 누락: " + ", ".join(missing))
    if not isinstance(config["mode"], str):
        return unavailable("config_invalid", "config.mode는 문자열이어야 합니다.")
    if config["mode"] != "reads":
        return unavailable("mode_not_supported", "이 진단은 reads 모드의 개별 리드 길이에만 적용합니다.")
    limits = {"min_k": (1, 64), "max_k": (1, 64),
              "seed_bases": (1, 100_000_000), "min_copies": (2, 1_000_000)}
    for key, (lower, upper) in limits.items():
        value = config[key]
        if isinstance(value, bool) or not isinstance(value, int) or not lower <= value <= upper:
            return unavailable("config_invalid", f"config.{key}는 {lower}..{upper} 범위의 정수여야 합니다.")
    if config["min_k"] > config["max_k"]:
        return unavailable("config_invalid", "config.min_k가 max_k보다 큽니다.")
    result["config"] = {key: config[key] for key in CONFIG_FIELDS}
    for period in range(config["min_k"], config["max_k"] + 1):
        result["periods"].append({"kind": "configured_period", "period_bp": period,
                                  "minimum_seed_bp": None,
                                  "length_eligible_reads_lower": None,
                                  "length_eligible_reads_upper": None,
                                  "bounds_method": None,
                                  "analyzed_records_length_status": None})
    for row in result["periods"] + result["families"]:
        row["minimum_seed_bp"] = required_seed_bases(row["period_bp"], config["seed_bases"], config["min_copies"])
    for row in result["families"]:
        row["outside_discovery_range"] = not config["min_k"] <= row["period_bp"] <= config["max_k"]
    try:
        qc = _read_qc(report)
        # The shared parser compares duplicate totals by equality; disallow
        # JSON bool/float lookalikes as well as contradictory numeric totals.
        for key in ("records", "bases", "max_record_bases"):
            if key in report.get("evidence_qc", {}):
                value = report["evidence_qc"][key]
                if isinstance(value, bool) or not isinstance(value, int):
                    raise ValueError(f"evidence_qc.{key} must be an integer")
        for row in result["periods"] + result["families"]:
            bound = eligible_read_bounds(qc["records"], qc["bases"], qc["max_record_bases"],
                                         row["minimum_seed_bp"], qc["length_histogram"])
            if not 0 <= bound["lower"] <= bound["upper"] <= qc["records"]:
                raise ValueError("inconsistent eligible-read bounds")
            row.update(length_eligible_reads_lower=bound["lower"],
                       length_eligible_reads_upper=bound["upper"], bounds_method=bound["method"],
                       analyzed_records_length_status=_length_status(bound["lower"], bound["upper"]))
    except (ValueError, AttributeError, TypeError, OverflowError) as error:
        for row in result["periods"] + result["families"]:
            row.update(length_eligible_reads_lower=None, length_eligible_reads_upper=None,
                       bounds_method=None, analyzed_records_length_status=None)
        return unavailable("qc_invalid", "길이 QC를 검증할 수 없습니다: " + str(error))
    for key in ("records", "bases", "max_record_bases"):
        result[key] = qc[key]
    result["length_distribution_available"] = qc["length_distribution_available"]
    result["sampled_prefix"] = qc["sampled_prefix"]
    incomplete = (analysis_complete is not True or qc["sampled_prefix"] is True
                  or any(str(report.get(key, "")).lower() in FAILED_STATES
                         for key in ("status", "state", "analysis_status")))
    if incomplete:
        return _finish(result, "analysis_incomplete", "analyzed_records",
                       "분석 미완료 또는 prefix 처리 결과입니다. 길이 bounds는 처리된 리드에만 해당합니다.",
                       "전체 입력을 검증하며 분석을 완료한 뒤 미검출 결과를 해석하세요.")
    if qc["sampled_prefix"] is not False:
        return _finish(result, "input_completeness_unverified", "unverified_input",
                       "sampled_prefix=False가 명시되지 않아 전체 입력 처리 여부를 확인할 수 없습니다.",
                       "전체 입력 처리 완료 여부를 검증하세요. 현재 길이 bounds는 보고된 리드에만 해당합니다.")
    statuses = [row["analyzed_records_length_status"] for row in result["periods"]]
    if all(status == "no_length_eligible_reads" for status in statuses):
        status = "no_length_eligible_reads"
        reason = "설정된 모든 검색 주기의 최소 seed 길이를 만족하는 리드 수 상한이 0입니다."
        recommendation = "현재 입력 길이와 seed 조건의 부적합을 확인하고 더 긴 리드 또는 명시적으로 재검토한 검색 조건을 사용하세요."
    elif "length_eligible_reads_present_power_unknown" in statuses:
        status = "length_eligible_reads_present_power_unknown"
        reason = "적어도 한 검색 주기에서 최소 seed 길이를 만족하는 리드가 있습니다. 모든 주기의 분석 가능성을 뜻하지 않습니다."
        recommendation = "예상 family별 길이 조건과 별도의 데이터 적정성·관측 근거를 함께 확인하세요."
    else:
        status = "length_eligibility_unresolved"
        reason = "집계값·구간형 길이 분포만으로 최소 seed 길이를 만족하는 리드가 있는지 확정할 수 없습니다."
        recommendation = "더 세밀한 길이 QC가 확보되기 전까지 길이 부적합이나 검출력 충분을 단정하지 마세요."
    return _finish(result, status, "complete_input", reason, recommendation)


def build_input_applicability(report, expected_motifs=(), *, analysis_complete=True,
                              source_available=True):
    """Preserve seed diagnostics and add a separate flank length opportunity."""
    from flank_applicability import build_flank_applicability

    result = _build_seed_applicability(report, expected_motifs,
                                     analysis_complete=analysis_complete,
                                     source_available=source_available)
    flank = build_flank_applicability(report, result)
    result['schema_version'] = 2
    result['flank_period_applicability'] = flank
    return result
