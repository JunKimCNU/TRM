#!/usr/bin/env python3
"""Render a separate genome-size feasibility report without changing TRM calls.

Only saved metadata, KMC diagnostics and GenomeScope output files are read.
An accepted fit is an exploratory model result, never a biological validation.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import html
import io
import json
import math
from pathlib import Path
import re


NUMBER = r"[-+]?(?:\d[\d,]*\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
ESTIMATE_FIELDS = (
    "sample_id", "scientific_name", "k", "ploidy_assumption", "max_kmer_coverage",
    "fit_status", "fit_guard_reasons", "progress_final_token", "execution_returncode",
    "genome_haploid_length_min_bp", "genome_haploid_length_max_bp",
    "raw_total_bases", "nominal_base_coverage_min", "nominal_base_coverage_max",
    "model_kmercov", "model_length", "all_kmer_fit_score_percent", "full_model_fit_score_percent",
    "read_error_rate_min_percent", "read_error_rate_max_percent",
    "genome_repeat_length_min_bp", "genome_repeat_length_max_bp",
    "counter_cap", "kmers_at_counter_cap", "counter_saturated",
    "counted_kmer_occurrences_lower_bound", "occurrences_above_fit_cutoff",
    "counted_occurrence_fraction_above_fit_cutoff", "kmc_total_occurrences",
    "exact_tail_guard_passed", "exact_occurrences_above_fit_cutoff", "exact_occurrence_fraction_above_fit_cutoff",
    "counter_clipped_occurrences",
    "all_local_reads", "canonical_kmers",
    "quality_trimming", "adapter_trimming", "summary_path", "model_path",
    "progress_path", "linear_plot_path", "log_plot_path", "result_path",
)
RUNTIME_FIELDS = (
    "sample_id", "k", "stage", "ploidy_assumption", "max_kmer_coverage",
    "returncode", "wall_seconds", "user_cpu_seconds", "system_cpu_seconds",
    "max_rss_kb", "execution_path",
)
PROTOCOL_FIELDS = ("sample_id", "scientific_name", "previous_trm_classification", "previous_detected_motif",
                   "library_strategy", "library_source", "library_construction_protocol", "protocol_conflict_flags",
                   "recommendation", "library_source_url", "previous_assessment_json")
STATUS_LABELS = {
    "candidate_for_spectrum_pilot": "spectrum 분석 후보",
    "review_protocol_conflict": "WGS / protocol 불일치",
    "review_library_structure": "library 구조 검토",
    "review_population_model": "집단 혼합 모델 검토",
    "review_non_genomic_wgs": "genomic WGS 여부 검토",
    "outside_eukaryotic_trm_scope": "진핵생물 분석 범위 밖",
    "nonstandard_platform": "platform 오류 모델 검토",
    "review_amplification": "증폭 편향 검토",
    "review_input_qc": "입력 QC 검토",
    "failed_input": "입력 분석 실패",
}


def read_json(path: Path, default=None):
    return json.loads(path.read_text()) if path.is_file() else default


def summary_interval(text: str, property_name: str, *, ordered=True):
    """Parse exactly two finite numeric values on a named GenomeScope row."""
    lines = [line[len(property_name):].strip() for line in text.splitlines()
             if line.startswith(property_name) and
             (len(line) == len(property_name) or line[len(property_name)].isspace())]
    if len(lines) != 1:
        return None
    tokens = re.fullmatch(rf"({NUMBER})\s*(?:bp|%)?\s+({NUMBER})\s*(?:bp|%)?", lines[0])
    if not tokens:
        return None
    values = tuple(float(value.replace(",", "")) for value in tokens.groups())
    return values if all(math.isfinite(value) for value in values) and (not ordered or values[0] <= values[1]) else None


def model_parameter(text: str, name: str):
    lines = [line for line in text.splitlines() if re.match(rf"^\s*{re.escape(name)}\s+", line)]
    if len(lines) != 1:
        return None
    tokens = lines[0].split()
    try:
        value = float(tokens[1])
    except (ValueError, IndexError):
        return None
    return value if math.isfinite(value) else None


def inspect_fit(directory: Path):
    """Require completed execution, a success marker and positive size/coverage.

    The installed R package reports `done` on success; other GenomeScope2
    distributions use `pass`. Both still require valid fitted parameters.
    A summary file alone is not a successful fit.
    """
    execution = read_json(directory / "fit.execution.json")
    summary = (directory / "summary.txt").read_text() if (directory / "summary.txt").is_file() else ""
    model = (directory / "model.txt").read_text() if (directory / "model.txt").is_file() else ""
    progress = (directory / "progress.txt").read_text().split() if (directory / "progress.txt").is_file() else []
    marker = progress[-1] if progress else ""
    length = summary_interval(summary, "Genome Haploid Length")
    params = {key: model_parameter(model, key) for key in ("length", "kmercov")}
    reasons = []
    if execution is None:
        reasons.append("execution_not_finished_or_not_recorded")
    elif execution.get("returncode") != 0:
        reasons.append("nonzero_exit")
    if marker not in {"pass", "done"}:
        reasons.append("missing_success_marker")
    if length is None or length[0] <= 0:
        reasons.append("nonpositive_or_missing_genome_length")
    for key, value in params.items():
        if value is None or value <= 0:
            reasons.append(f"nonpositive_or_missing_model_{key}")
    status = "pending_or_unrecorded" if execution is None else ("rejected_model_output" if reasons else "exploratory_fit")
    values = {"fit_status": status, "fit_guard_reasons": ";".join(reasons),
              "progress_final_token": marker,
              "execution_returncode": execution.get("returncode") if execution else None,
              "model_kmercov": params["kmercov"], "model_length": params["length"]}
    for source, target in (("Genome Haploid Length", "genome_haploid_length"),
                           ("Genome Repeat Length", "genome_repeat_length"),
                           ("Read Error Rate", "read_error_rate")):
        interval = summary_interval(summary, source)
        unit = "bp" if "Length" in source else "percent"
        for index, bound in enumerate(("min", "max")):
            values[f"{target}_{bound}_{unit}"] = interval[index] if interval else None
    scores = summary_interval(summary, "Model Fit", ordered=False)
    values["all_kmer_fit_score_percent"] = scores[0] if scores else None
    values["full_model_fit_score_percent"] = scores[1] if scores else None
    # Preserve rejected values for diagnosis; consumers must check fit_status.
    return values


def relative(path: Path, root: Path):
    return str(path.resolve().relative_to(root.resolve()))


def exact_tail_audit(histogram, stats, cutoff):
    """Recover exact excluded occurrence mass from KMC totals when auditable.

    Counter saturation affects bins above these cutoffs. The retained bins are
    exact only when cutoff < counter cap and no KMC threshold filters fired.
    """
    total = stats.get("#Total no. of k-mers")
    observed = histogram.get("counted_kmer_occurrences_lower_bound")
    tail = histogram.get("occurrences_above_fit_cutoff", {}).get(str(cutoff))
    guard = (isinstance(total, int) and total > 0 and isinstance(observed, int) and
             isinstance(tail, int) and 0 <= tail <= observed <= total and
             stats.get("#k-mers_below_min_threshold") == 0 and
             stats.get("#k-mers_above_max_threshold") == 0 and
             stats.get("#Unique_counted_k-mers") == histogram.get("distinct_kmers") and
             isinstance(histogram.get("counter_cap"), int) and cutoff < histogram["counter_cap"])
    excluded = total - (observed - tail) if guard else None
    return {"kmc_total_occurrences": total, "exact_tail_guard_passed": guard,
            "exact_occurrences_above_fit_cutoff": excluded,
            "exact_occurrence_fraction_above_fit_cutoff": excluded / total if guard else None,
            "counter_clipped_occurrences": total - observed if guard else None}


def collect_pilot(root: Path, samples: list[dict]):
    pilot = root / "pilot"
    configuration = read_json(pilot / "run_configuration.json", {})
    arguments = configuration.get("arguments", {})
    names = {row["sample_id"]: row["scientific_name"] for row in samples}
    bases = {row["sample_id"]: int(row["bases"]) if row.get("bases") else None for row in samples}
    estimates, runtimes, artifacts = [], [], set()
    for sid in arguments.get("samples", []):
        if sid not in names or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sid):
            raise ValueError(f"Unknown or unsafe pilot sample ID: {sid}")
        for k in arguments.get("k", []):
            directory = pilot / sid / f"k{k}"
            result_path = directory / "result.json"
            result = read_json(result_path, {})
            if result and (result.get("sample_id") != sid or result.get("k") != k):
                raise ValueError(f"Pilot identity mismatch: {result_path}")
            hist = result.get("histogram", {})
            stats_path = directory / "kmc_summary.json"
            kmc_stats = read_json(stats_path, {}).get("Stats", {})
            if stats_path.is_file():
                artifacts.add(stats_path)
            if result:
                artifacts.add(result_path)
            for cutoff in arguments.get("cutoffs", []):
                for ploidy in arguments.get("ploidies", []):
                    out = directory / f"p{ploidy}_max{cutoff}"
                    row = {"sample_id": sid, "scientific_name": names[sid], "k": k,
                           "ploidy_assumption": ploidy, "max_kmer_coverage": cutoff,
                           **inspect_fit(out), "counter_cap": hist.get("counter_cap", arguments.get("counter_cap")),
                           "kmers_at_counter_cap": hist.get("kmers_at_counter_cap"),
                           "counter_saturated": hist.get("saturated"),
                           "counted_kmer_occurrences_lower_bound": hist.get("counted_kmer_occurrences_lower_bound"),
                           "occurrences_above_fit_cutoff": hist.get("occurrences_above_fit_cutoff", {}).get(str(cutoff)),
                           "counted_occurrence_fraction_above_fit_cutoff": hist.get("occurrence_fraction_above_fit_cutoff", {}).get(str(cutoff)),
                           **exact_tail_audit(hist, kmc_stats, cutoff),
                           "all_local_reads": result.get("all_local_reads"),
                           "canonical_kmers": result.get("canonical"),
                           "quality_trimming": result.get("quality_trimming"),
                           "adapter_trimming": result.get("adapter_trimming"),
                           "result_path": relative(result_path, root) if result else ""}
                    row["raw_total_bases"] = bases[sid]
                    valid = row["fit_status"] == "exploratory_fit" and bases[sid] is not None
                    row["nominal_base_coverage_min"] = bases[sid] / row["genome_haploid_length_max_bp"] if valid else None
                    row["nominal_base_coverage_max"] = bases[sid] / row["genome_haploid_length_min_bp"] if valid else None
                    for field, filename in (("summary_path", "summary.txt"), ("model_path", "model.txt"),
                                            ("progress_path", "progress.txt"), ("linear_plot_path", "linear_plot.png"),
                                            ("log_plot_path", "log_plot.png")):
                        path = out / filename
                        row[field] = relative(path, root) if path.is_file() else ""
                        if path.is_file() and path.suffix == ".txt":
                            artifacts.add(path)
                    estimates.append(row)
            for path in directory.glob("**/*.execution.json"):
                execution = read_json(path)
                parent_match = re.fullmatch(r"p(\d+)_max(\d+)", path.parent.name)
                runtimes.append({"sample_id": sid, "k": k, "stage": path.name.split(".")[0],
                                 "ploidy_assumption": int(parent_match[1]) if parent_match else None,
                                 "max_kmer_coverage": int(parent_match[2]) if parent_match else None,
                                 **{key: execution.get(key) for key in RUNTIME_FIELDS[5:-1]},
                                 "execution_path": relative(path, root)})
                artifacts.add(path)
    runtimes.sort(key=lambda row: (row["sample_id"], row["k"], row["stage"],
                                  row["max_kmer_coverage"] or 0, row["ploidy_assumption"] or 0))
    final = read_json(pilot / "results.json")
    if final is not None:
        expected = {(sid, k) for sid in arguments.get("samples", []) for k in arguments.get("k", [])}
        identities = [(item.get("sample_id"), item.get("k")) for item in final.get("results", [])]
        if len(identities) != len(expected) or set(identities) != expected:
            raise ValueError("Final pilot results do not match configured sample / k combinations")
        if any(row["fit_status"] == "pending_or_unrecorded" for row in estimates):
            raise ValueError("Final pilot results exist but a configured fit has no completed execution record")
    for path in (pilot / "results.json", pilot / "run_configuration.json"):
        if path.is_file():
            artifacts.add(path)
    return estimates, runtimes, configuration, final, artifacts


def sample_model_ranges(estimates):
    ranges = []
    for sid in dict.fromkeys(row["sample_id"] for row in estimates):
        rows = [row for row in estimates if row["sample_id"] == sid]
        valid = [row for row in rows if row["fit_status"] == "exploratory_fit"]
        low = min((r["genome_haploid_length_min_bp"] for r in valid), default=None)
        high = max((r["genome_haploid_length_max_bp"] for r in valid), default=None)
        ranges.append({"sample_id": sid, "scientific_name": rows[0]["scientific_name"],
                       "planned_models": len(rows), "exploratory_models": len(valid),
                       "minimum_over_models_bp": low, "maximum_over_models_bp": high,
                       "maximum_to_minimum_ratio": high / low if low else None,
                       "minimum_all_kmer_fit_score_percent": min((r["all_kmer_fit_score_percent"] for r in valid if r["all_kmer_fit_score_percent"] is not None), default=None),
                       "review_required": True,
                       "interpretation": "Range across explicit model assumptions, not a biological confidence interval or a chosen genome-size estimate."})
    return ranges


def join_previous_observations(samples):
    """Expose protocol conflicts without silently changing frozen TRM outputs."""
    artifacts, review = set(), []
    conflict_codes = {"rna_protocol_conflicts_with_wgs_label", "hic_protocol_conflicts_with_wgs_label"}
    for row in samples:
        path = Path(row["assessment_json"])
        data = read_json(path)
        if data is None or data.get("sample_id") != row["sample_id"]:
            raise ValueError(f"Missing or mismatched previous assessment: {path}")
        if hashlib.sha256(path.read_bytes()).hexdigest() != row["assessment_sha256"]:
            raise ValueError(f"Previous assessment changed after metadata audit: {path}")
        row["previous_trm_classification"] = data.get("classification", "")
        row["previous_detected_motif"] = data.get("detected_motif") or ""
        artifacts.add(path)
        conflicts = conflict_codes.intersection(row.get("all_flags", "").split(";"))
        if conflicts:
            review.append({key: row.get(key, "") for key in PROTOCOL_FIELDS})
            review[-1].update({"protocol_conflict_flags": ";".join(sorted(conflicts)),
                               "previous_assessment_json": str(path),
                               "recommendation": "Resolve the WGS/protocol metadata conflict and verify ordinary genomic WGS before interpreting prior repeat observations as telomeric or using a genome-size/absence model; obtain independent conventional WGS if unresolved."})
    return artifacts, review


def number(value, digits=0, scale=1):
    if value in (None, ""):
        return "미확인"
    return f"{float(value) / scale:,.{digits}f}"


def e(value):
    return html.escape(str(value), quote=True)


def link(path, label):
    return f'<a href="{e(path)}">{e(label)}</a>' if path else ""


def render_html(samples, metadata, audit, estimates, runtimes, summary):
    cards = []
    for row in summary["pilot_sample_model_ranges"]:
        low, high = row["minimum_over_models_bp"], row["maximum_over_models_bp"]
        extent = f'{number(low, 2, 1e6)}–{number(high, 2, 1e6)} Mb' if low else "완료된 유효 모델 없음"
        cards.append(f'<article><b>{e(row["sample_id"])}</b><br><em>{e(row["scientific_name"])}</em>'
                     f'<div class="big">{extent}</div><span>{row["exploratory_models"]}/{row["planned_models"]} 모델의 계산 조건 통과</span>'
                     f'<br>최대 / 최소 {number(row["maximum_to_minimum_ratio"], 2)}배'
                     f'<br>all-kmer fit 최저 {number(row["minimum_all_kmer_fit_score_percent"], 2)}% · 해석 검토 필요</article>')
    model_rows = []
    for row in estimates:
        accepted = row["fit_status"] == "exploratory_fit"
        length = (f'{number(row["genome_haploid_length_min_bp"], 2, 1e6)}–{number(row["genome_haploid_length_max_bp"], 2, 1e6)}'
                  if accepted else "—")
        fraction = row["exact_occurrence_fraction_above_fit_cutoff"]
        tail = number(float(fraction) * 100, 3) + "%" if fraction is not None else "정확한 비율 미확인"
        saturation = ("없음" if row["counter_saturated"] is False else
                      (f'있음 ({number(row["kmers_at_counter_cap"])})' if row["counter_saturated"] else "미확인"))
        status = {"exploratory_fit": "탐색적 추정", "rejected_model_output": "모델 결과 거부",
                  "pending_or_unrecorded": "미완료 / 기록 대기"}[row["fit_status"]]
        plots = " · ".join(filter(None, (link(row["linear_plot_path"], "linear"), link(row["log_plot_path"], "log"),
                                             link(row["summary_path"], "summary"), link(row["model_path"], "model"), link(row["progress_path"], "progress"))))
        model_rows.append(f'<tr><td>{e(row["sample_id"])}</td><td>{row["k"]}</td><td>{row["ploidy_assumption"]}</td>'
                          f'<td>{number(row["max_kmer_coverage"])}</td><td title="{e(row["fit_guard_reasons"])}">{status}</td>'
                          f'<td>{length}</td><td>{number(row["model_kmercov"], 2)}</td>'
                          f'<td>{number(row["nominal_base_coverage_min"], 1)}–{number(row["nominal_base_coverage_max"], 1)}×</td>'
                          f'<td>{number(row["all_kmer_fit_score_percent"], 2)} / {number(row["full_model_fit_score_percent"], 2)}%</td>'
                          f'<td>{tail}</td><td>{saturation}</td><td>{plots}</td></tr>')
    sample_rows = []
    for row in samples:
        flags = row.get("all_flags", "")
        try:
            reasons = json.loads(row.get("review_reasons") or "[]")
        except ValueError:
            reasons = []
        details = "".join(f'<li><code>{e(item.get("code", ""))}</code>: {e(item.get("reason", ""))}</li>' for item in reasons)
        protocol = row.get("library_construction_protocol") or "미기재"
        detail = (f'<details><summary>{e(protocol)}</summary><p>Flags: {e(flags or "없음")}</p><ul>{details}</ul>'
                  f'<p>Library PCR: {e(row.get("library_pcr_status", ""))}</p>'
                  f'<p>WGA archive: {e(row.get("archive_wga_mda_status", ""))}; 생물학적 WGA 여부 미확인</p>'
                  f'<p>Local layout: {e(row.get("local_layout_interpretation", ""))}</p>'
                  f'{link(row.get("library_source_url", ""), "원문 metadata")}</details>')
        status = row["genome_estimation_status"]
        sample_rows.append(f'<tr data-search="{e(" ".join((row["sample_id"], row["scientific_name"], status, flags, protocol, row.get("previous_trm_classification", ""), row.get("previous_detected_motif", ""))).lower())}">'
                           f'<td>{e(row["sample_id"])}</td><td><em>{e(row["scientific_name"])}</em></td>'
                           f'<td title="{e(status)}">{e(STATUS_LABELS.get(status, status))}</td>'
                           f'<td>{e(row.get("previous_trm_classification", ""))}<br><code>{e(row.get("previous_detected_motif") or "—")}</code></td>'
                           f'<td>{number(row.get("records"))}</td><td>{number(row.get("bases"), 3, 1e9)}</td>'
                           f'<td>{number(row.get("mean_read_length"), 1)}</td><td>{e(row.get("instrument_platform", ""))}</td>'
                           f'<td>{e(row.get("library_strategy", ""))} / {e(row.get("library_source", ""))}</td>'
                           f'<td>{detail}</td></tr>')
    status_counts = " · ".join(f'{e(STATUS_LABELS.get(key, key))} {value}건' for key, value in metadata.get("status_counts", {}).items())
    overall = summary["concurrent_overall_wall_seconds"]
    overall_text = f'{number(overall, 1)}초 ({number(overall, 1, 60)}분)' if overall is not None else "pilot 전체 실행 미완료"
    resources = summary["runtime"]
    completed_label = "완료" if summary["pilot_completed"] else "진행 중인 결과 스냅샷"
    conflict_novel_ids = ", ".join(summary["protocol_conflict_previous_novel_candidate_ids"]) or "없음"
    return f'''<!doctype html>
<html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Genome size 추정 가능성 검토</title>
<style>
body{{font:15px/1.55 system-ui,sans-serif;color:#21323d;background:#f6f8fa;margin:0;padding:30px;max-width:1600px}}
h1{{font-size:28px}}h2{{font-size:21px;margin-top:32px}}a{{color:#126793}}.muted{{color:#596b77}}
.note{{padding:15px 20px;background:#fff3d8;border-left:5px solid #d49322;border-radius:4px}}
.cards{{display:flex;gap:16px;flex-wrap:wrap}}article{{background:white;border:1px solid #dce3e9;border-radius:8px;padding:18px;min-width:250px}}
.big{{font-size:24px;font-weight:650;padding:10px 0}}.scroll{{overflow:auto;max-height:70vh;background:white;border:1px solid #dce3e9}}
table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{padding:10px;border-bottom:1px solid #e1e6eb;text-align:left;vertical-align:top}}
th{{position:sticky;top:0;background:#e9eff4;white-space:nowrap}}td:nth-child(1){{white-space:nowrap}}details{{min-width:200px;max-width:520px}}
input{{padding:10px;margin:8px 0 14px;width:min(600px,85%);font-size:15px}}code{{overflow-wrap:anywhere}}.links{{line-height:2}}
</style>
<h1>현재 FASTQ의 genome size 추정 가능성</h1>
<p class="muted">{e(summary["created_utc"])} · {e(completed_label)} · 현재 TRM 결과 및 검출력 판정은 변경하지 않았습니다.</p>
<p>기존 TRM 분석은 리드 내 tandem repeat을 찾습니다. Genome size에는 별도의 전체 genome k-mer 빈도 분포가 필요합니다.
이번 pilot은 선택한 샘플의 모든 local FASTQ 리드를 KMC로 집계하고 GenomeScope2 모델을 비교합니다.</p>
<div class="note">제공된 genome size와 염색체 수는 모두 미확인입니다. 아래 수치는 가정한 ploidy·k·고빈도 cutoff에 따른 <b>탐색적 모델 결과</b>입니다.
모델 계산 완료가 genome size 정확도, 실제 ploidy·염색체 수, 균일한 sequencing coverage, 또는 telomere 미검출의 충분한 검출력을 입증하지 않습니다.</div>
<p class="links"><a href="samples.tsv">전체 샘플 metadata / 적합성 TSV</a> · <a href="pilot_estimates.tsv">모든 pilot 모델 TSV</a> ·
<a href="runtime.tsv">단계별 실행 시간 TSV</a> · <a href="report_summary.json">집계 JSON</a> ·
<a href="protocol_review.tsv">Protocol 불일치 / 기존 TRM 결과 재검토 TSV</a> ·
<a href="library_metadata_audit.json">공개 metadata audit</a> · <a href="pilot/run_configuration.json">도구·실행 설정</a> · <a href="README.md">해석 문서</a></p>
<h2>Pilot 모델 비교</h2><div class="cards">{''.join(cards)}</div>
<p>위 범위는 각 샘플의 계산 조건을 통과한 모든 모델의 최솟값–최댓값이며, <b>생물학적 신뢰구간이 아닙니다</b>.
아래 genome length min/max는 GenomeScope 출력의 parameter 불확실성 범위입니다. p=1 / p=2는 입력한 모델 가정이며 실제 ploidy 추정값이 아닙니다.
kmercov는 모델의 coverage parameter이며, 그대로 nucleotide coverage 또는 chromosome 수로 해석하지 않습니다.</p>
<div class="scroll"><table><thead><tr><th>샘플</th><th>k</th><th>p 가정</th><th>최대 빈도</th><th>상태</th><th>크기 min–max (Mb)</th>
<th>kmercov</th><th>가정상 bases / G</th><th>Fit: all-kmer / full-model</th><th>cutoff 위 occurrence</th><th>counter cap 도달</th><th>원문 / 그래프</th></tr></thead><tbody>{''.join(model_rows)}</tbody></table></div>
<p class="muted">탐색적 추정의 계산 조건: 완료된 프로세스의 exit code 0, progress 성공 표시, 유한한 양수 genome length 및 model length/kmercov.
Fit 두 값은 설치된 GenomeScope의 all-kmer / full-model 기술적 점수이며 신뢰구간 또는 정답 확률이 아닙니다.
각각 오류 구간을 제외한 전체 평가 구간 / 주 peak 모델 구간의 weighted absolute residual 점수이며, 리드가 genome에 매핑된 비율도 아닙니다.
GenomeScope의 모델 fitting과 genome length 합계 모두 cutoff 위 mass를 제외합니다. 표의 제외 비율은 KMC 전체 occurrence에서 cutoff 이하 histogram mass를 뺀 값입니다.
필터 제외가 0이고 distinct count가 일치하며 cutoff가 counter cap보다 낮을 때만 이 값을 제공합니다. Counter가 포화된 고빈도 histogram은 완전한 빈도 spectrum이 아닙니다.
고빈도 tail에는 실제 반복서열과 organelle·오염 등이 함께 들어갈 수 있습니다. 모든 모델은 생물학적 검토가 필요합니다.</p>
<p class="muted">‘가정상 bases / G’는 저장된 전체 리드의 base 수를 각 모델 genome size 범위로 나눈 값입니다.
오염·organelle·중복·overlapping mates의 base까지 포함하며, 실제 유효 coverage 또는 telomere 검출력을 의미하지 않습니다.</p>
<h2>실행 시간과 자원</h2>
<p>동시 실행 전체 wall time: <b>{overall_text}</b>. 완료된 각 단계의 wall time 합: {number(resources["sum_step_wall_seconds"], 1)}초;
user CPU 합 {number(resources["sum_user_cpu_seconds"], 1)}초, system CPU 합 {number(resources["sum_system_cpu_seconds"], 1)}초.
단일 단계 최대 RSS: {number(resources["largest_step_rss_kb"], 1, 1024)} MiB.</p>
<p class="muted">단계별 wall time은 병렬 실행으로 겹치므로 전체 경과 시간과 다릅니다. 단일 단계 RSS는 전체 batch의 peak memory가 아닙니다.
다운로드·환경 설치·metadata 조사·report 작성 시간은 pilot 시간에 포함되지 않습니다.</p>
<h2>전체 {len(samples)}개 샘플의 metadata 기반 적합성</h2>
<p>{status_counts}</p>
<p>WGS 표기와 protocol을 따로 확인합니다. 공개 library protocol에서 {audit.get("protocol_label_conflict_count", 0)}건의 WGS 표기 불일치가 발견되었습니다.
전 샘플 WGA 여부는 미확인이며, “No PCR”와 WGA 미기재는 WGA 부재의 증거가 아닙니다. Metadata는 ENA / INSDC 공유 archive 자료와 별도 NCBI 표본 확인을 사용했습니다.</p>
<p>‘spectrum 분석 후보’는 데이터 양·coverage가 충분하다는 판정이 아닙니다. 현재 표의 상태는 전체 230개 샘플에 대한 pilot 실행 전 적합성 검토입니다.
Pilot 수치는 위 표에서만 제공하며 기존 genome size unknown 및 TRM adequacy 판정에 자동 반영하지 않습니다.</p>
<p>기존 TRM 관측값도 함께 표시합니다. Protocol이 RNA PolyA / Hi-C와 충돌한 샘플은 해당 반복서열 관측을 novel telomere 근거로 해석하기 전에
실제 library 유형을 확인해야 합니다. Protocol이 충돌하는 기존 novel 후보: {e(conflict_novel_ids)}.</p>
<input id="search" aria-label="샘플 및 조건 검색" placeholder="샘플 ID, 종, protocol, flag 검색"><span id="visible">{len(samples)} / {len(samples)}</span>
<div class="scroll"><table id="samples"><thead><tr><th>샘플</th><th>종</th><th>적합성</th><th>기존 TRM 관측</th><th>리드 수</th><th>Gb</th><th>평균 bp</th><th>Platform</th><th>archive 표기</th><th>Protocol / 근거</th></tr></thead><tbody>{''.join(sample_rows)}</tbody></table></div>
<script>
const rows=Array.from(document.querySelectorAll('#samples tbody tr'));
document.getElementById('search').addEventListener('input',event=>{{const query=event.target.value.toLowerCase().trim();let count=0;
for(const row of rows){{row.hidden=!row.dataset.search.includes(query);if(!row.hidden)count++;}}
document.getElementById('visible').textContent=count+' / '+rows.length;}});
</script></html>'''


def write_tsv(path, rows, fields):
    content = io.StringIO(newline="")
    writer = csv.DictWriter(content, fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
    writer.writeheader()
    writer.writerows(rows)
    atomic_write(path, content.getvalue())


def atomic_write(path, text):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text)
    temporary.replace(path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=Path("results/genome_feasibility_20260921"))
    args = parser.parse_args(argv)
    root = args.results_dir.resolve()
    try:
        if root.name.startswith("wgs_assessed") or root.name.startswith("comparison_"):
            raise ValueError("Use a separate feasibility output directory")
        with (root / "samples.tsv").open(newline="") as handle:
            samples = list(csv.DictReader(handle, delimiter="\t"))
        metadata = read_json(root / "metadata_summary.json")
        audit = read_json(root / "library_metadata_audit.json")
        if len(samples) != metadata["sample_count"] or len({r["sample_id"] for r in samples}) != len(samples):
            raise ValueError("Metadata sample count or sample IDs disagree")
        if hashlib.sha256((root / "samples.tsv").read_bytes()).hexdigest() != metadata["samples_tsv_sha256"]:
            raise ValueError("samples.tsv does not match metadata provenance")
        estimates, runtimes, config, final, artifacts = collect_pilot(root, samples)
        previous_artifacts, protocol_review = join_previous_observations(samples)
        artifacts.update(previous_artifacts)
        stage_stats = {stage: {"completed_steps": len(rows), "wall_seconds_sum": sum(r["wall_seconds"] or 0 for r in rows)}
                       for stage in sorted({r["stage"] for r in runtimes})
                       for rows in [[r for r in runtimes if r["stage"] == stage]]}
        summary = {
            "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
            "sample_count": len(samples), "pilot_completed": final is not None,
            "pilot_model_count": len(estimates), "pilot_model_status_counts": dict(Counter(r["fit_status"] for r in estimates)),
            "pilot_sample_model_ranges": sample_model_ranges(estimates),
            "concurrent_overall_wall_seconds": final.get("elapsed_seconds") if final else None,
            "runtime": {"completed_steps": len(runtimes), "by_stage": stage_stats,
                        "sum_step_wall_seconds": sum(r["wall_seconds"] or 0 for r in runtimes),
                        "sum_user_cpu_seconds": sum(r["user_cpu_seconds"] or 0 for r in runtimes),
                        "sum_system_cpu_seconds": sum(r["system_cpu_seconds"] or 0 for r in runtimes),
                        "largest_step_rss_kb": max((r["max_rss_kb"] or 0 for r in runtimes), default=None),
                        "largest_step_rss_is_batch_peak": False},
            "ploidy_inferred": False, "karyotype_inferred": False, "trm_adequacy_modified": False,
            "trusted_genome_size_integrated": False,
            "pilot_model_range_is_biological_confidence_interval": False,
            "model_fit_guard": "Completed process exit 0; progress final token pass or done; finite positive summary genome length and model length/kmercov.",
            "configuration": config.get("arguments", {}),
            "metadata_status_counts": metadata.get("status_counts", {}),
            "protocol_conflict_count": len(protocol_review),
            "protocol_conflict_previous_classifications": dict(Counter(r["previous_trm_classification"] for r in protocol_review)),
            "protocol_conflict_previous_novel_candidate_ids": [r["sample_id"] for r in protocol_review if r["previous_trm_classification"] == "novel_trm_candidate"],
        }
        artifacts.update(root / name for name in ("samples.tsv", "metadata_summary.json", "library_metadata_audit.json"))
        artifacts.add(Path(__file__).resolve())
        summary["inputs"] = [{"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                             for path in sorted(artifacts)]
        write_tsv(root / "pilot_estimates.tsv", estimates, ESTIMATE_FIELDS)
        write_tsv(root / "runtime.tsv", runtimes, RUNTIME_FIELDS)
        write_tsv(root / "protocol_review.tsv", protocol_review, PROTOCOL_FIELDS)
        atomic_write(root / "report_summary.json", json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
        atomic_write(root / "index.html", render_html(samples, metadata, audit, estimates, runtimes, summary))
        print(json.dumps({key: summary[key] for key in ("sample_count", "pilot_completed", "pilot_model_count", "pilot_model_status_counts")}, ensure_ascii=False))
        return 0
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    raise SystemExit(main())
