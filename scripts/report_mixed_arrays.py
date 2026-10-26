#!/usr/bin/env python3
"""Report explicit mixed-unit structural observations without changing TRM calls."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import html
import json
from pathlib import Path
import re


REPORT_VERSION = "0.7.0"
GROUP_COUNTS = (
    "records", "mixed_records", "pure_only_records", "nonqualifying_records", "total_mixed_arrays",
    "full_mixed_records", "entire_mixed_records", "repeat_bases", "max_array_bp",
    "max_any_array_bp", "long_array_records", "g_right_boundary_records",
    "c_left_boundary_records", "g_left_boundary_records", "c_right_boundary_records",
    "internal_array_records", "clean_two_flank_records", "left_g_end_records",
    "left_c_end_records", "right_g_end_records", "right_c_end_records",
    "long_expected_end_records", "segmentation_ambiguous_records",
    "ambiguous_orientation_records", "ambiguous_boundary_records",
)
MIXED_RECORD_COUNTS = (
    "long_array_records", "g_right_boundary_records", "c_left_boundary_records",
    "g_left_boundary_records", "c_right_boundary_records", "internal_array_records",
    "clean_two_flank_records", "left_g_end_records", "left_c_end_records",
    "right_g_end_records", "right_c_end_records", "long_expected_end_records",
    "ambiguous_boundary_records",
)
GROUP_FIELDS = (
    "sample_id", "group", "units_g_oriented", "observation", "input_records",
    "input_bases", *GROUP_COUNTS, "expected_boundary_records_sum",
    "opposite_boundary_records_sum", "examples_saved",
)
COMPONENT_FIELDS = (
    "sample_id", "group", "unit", "canonical_family", "previous_observation",
    "previous_reads", "previous_full_repeat_reads", "previous_repeat_bases",
    "previous_max_tract_bp", "previous_result_sha256",
)
EXAMPLE_FIELDS = (
    "sample_id", "group", "record_id", "record_length", "input_index", "start",
    "end", "span_bp", "strand", "touches_record_start", "touches_record_end",
    "biological_length_censored", "unit_counts", "native_example",
)
LIMITATIONS = [
    "명시한 반복 단위 집합의 targeted 구조 검사이며 de novo motif 발견이 아닙니다.",
    "기존 v0.6 TRM 판정과 데이터 충분성 평가는 변경하지 않습니다. 이 보고서는 canonical·novel·TRM 부재·ALT를 판정하지 않습니다.",
    "Exact 단위만 연결하며 오류·삽입/결실·불규칙 배열·제공하지 않은 단위·호환되지 않는 단위 경계는 놓칠 수 있습니다. 미관측으로 기술적 검출 실패를 배제할 수 없습니다.",
    "Component family 계수와 group 계수는 서로 겹칠 수 있습니다. 이들을 합해 고유 reads 또는 분자 수로 해석하지 않습니다.",
    "최대 span은 하나의 record에서 관측된 구간입니다. 여러 짧은 reads를 합해 500 bp 연속 배열로 해석하지 않습니다.",
    "Read 끝은 라이브러리 fragment 끝입니다. Assembly 끝의 방향 일치도 독립적인 물리적 말단 검증은 아닙니다.",
    "저장된 examples는 제한된 좌표 표본입니다. 전체 지지 records나 전체 검출 감도를 나타내지 않습니다.",
    "여러 단위 분해가 가능한 서열에서는 unit_counts와 unit_counts_max가 보수적 하한·주변 상한일 수 있습니다. 예시의 segmentation_ambiguous·orientation_ambiguous와 원시 JSON을 함께 확인합니다.",
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(unit):
    unit = unit.upper()
    if not re.fullmatch("[ACGT]+", unit):
        raise ValueError(f"Invalid DNA unit: {unit!r}")
    for size in range(1, len(unit) + 1):
        if len(unit) % size == 0 and unit[:size] * (len(unit) // size) == unit:
            unit = unit[:size]
            break
    rc = unit.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return min(s[i:] + s[:i] for s in (unit, rc) for i in range(len(unit)))


def nonnegative_integer(value, label):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{label}: required nonnegative integer measurement")
    return value


def validate_source(document):
    if document.get("analysis_kind") != "targeted_mixed_arrays":
        raise ValueError("--source must be a native targeted_mixed_arrays result")
    summary = document.get("summary", {})
    for key in ("records", "bases"):
        nonnegative_integer(summary.get(key), f"summary.{key}")
    block = document.get("mixed_arrays")
    if not isinstance(block, dict) or not isinstance(block.get("settings"), dict):
        raise ValueError("mixed_arrays.settings is required")
    fraction = block["settings"].get("full_fraction")
    if isinstance(fraction, bool) or not isinstance(fraction, (float, int)) or not 0 < fraction <= 1:
        raise ValueError("mixed_arrays.settings.full_fraction must be in (0,1]")
    groups = block.get("groups")
    if not isinstance(groups, list) or not groups:
        raise ValueError("mixed_arrays.groups must contain the supplied groups, including zeros")
    names = set()
    for group in groups:
        name = group.get("name")
        if not isinstance(name, str) or not name or name in names:
            raise ValueError("Mixed group names must be nonempty and unique")
        names.add(name)
        units = group.get("units")
        if not isinstance(units, list) or len(units) < 2:
            raise ValueError(f"{name}: at least two supplied units required")
        families = [canonical(unit) for unit in units]
        if len(set(families)) != len(families):
            raise ValueError(f"{name}: supplied units must be distinct primitive families")
        for key in GROUP_COUNTS:
            nonnegative_integer(group.get(key), f"{name}.{key}")
        if not (group["entire_mixed_records"] <= group["full_mixed_records"] <=
                group["mixed_records"] <= group["records"] <= summary["records"]):
            raise ValueError(f"{name}: incompatible record subset counts")
        if sum(group[key] for key in ("mixed_records", "pure_only_records", "nonqualifying_records")) != group["records"]:
            raise ValueError(f"{name}: mixed, pure-observed and nonqualifying records must partition records")
        if group["total_mixed_arrays"] < group["mixed_records"]:
            raise ValueError(f"{name}: fewer arrays than supporting records")
        if group["long_array_records"] > group["mixed_records"]:
            raise ValueError(f"{name}: long-array count exceeds mixed-record count")
        if group["max_array_bp"] > group["max_any_array_bp"]:
            raise ValueError(f"{name}: mixed maximum exceeds any-array maximum")
        if group["repeat_bases"] > summary["bases"]:
            raise ValueError(f"{name}: mixed union bases exceed all input bases")
        for key in MIXED_RECORD_COUNTS:
            if group[key] > group["mixed_records"]:
                raise ValueError(f"{name}.{key}: count exceeds mixed-record count")
        for key in ("segmentation_ambiguous_records", "ambiguous_orientation_records"):
            if group[key] > group["records"]:
                raise ValueError(f"{name}.{key}: count exceeds any-array record count")
        if group["long_expected_end_records"] > group["long_array_records"]:
            raise ValueError(f"{name}: long expected-end count exceeds long-array record count")
        if not isinstance(group.get("examples"), list):
            raise ValueError(f"{name}: examples list required (may be empty)")
        if not group["mixed_records"] and any(group[key] for key in (
                "total_mixed_arrays", "repeat_bases", "max_array_bp")):
            raise ValueError(f"{name}: mixed-array measurements exist without mixed records")
    return document


def check_previous_inputs(current, previous):
    """Check recorded provenance; this deliberately does not rehash FASTQ inputs."""
    def signature(document):
        inputs = document.get("inputs")
        if not isinstance(inputs, list) or not inputs:
            raise ValueError("Component comparison requires recorded input provenance")
        result = []
        for entry in inputs:
            if not isinstance(entry, dict) or not entry.get("path") or entry["path"] == "-":
                raise ValueError("Component comparison requires named input paths")
            if entry.get("bytes") is None or entry.get("modified_unix_seconds") is None:
                raise ValueError("Component comparison requires input size and modification time")
            result.append((str(Path(entry["path"]).resolve()), entry["bytes"],
                           entry["modified_unix_seconds"]))
        return result
    if signature(current) != signature(previous):
        raise ValueError("Previous and mixed results have different recorded inputs")
    for key in ("records", "bases"):
        if previous.get("summary", {}).get(key) != current["summary"][key]:
            raise ValueError(f"Previous and mixed results have different analyzed {key}")
    return {"status": "recorded_provenance_and_totals_match", "sequence_rehashed": False,
            "meaning": "Ordered input paths, sizes, modification times, analyzed records and bases match; sequence identity was not independently rehashed by this reporter."}


def component_rows(groups, previous, sample_id, previous_sha):
    if not any(isinstance(previous.get(key), list) for key in ("candidates", "rejected_candidates")):
        raise ValueError("Previous result must contain homogeneous native candidate arrays")
    index = {}
    for candidate in previous.get("candidates", []) + previous.get("rejected_candidates", []):
        key = canonical(candidate["motif"])
        if key in index:
            raise ValueError(f"Duplicate previous component family: {key}")
        index[key] = candidate
    rows = []
    for group in groups:
        for unit in group["units"]:
            key = canonical(unit)
            candidate = index.get(key)
            reads = None if candidate is None else nonnegative_integer(candidate.get("reads"), f"{key}.reads")
            rows.append({
                "sample_id": sample_id, "group": group["name"], "unit": unit,
                "canonical_family": key, "previous_observation": "not_reported" if candidate is None else
                ("reported_positive" if reads else "explicit_zero"), "previous_reads": reads,
                "previous_full_repeat_reads": None if candidate is None else candidate.get("full_repeat_reads"),
                "previous_repeat_bases": None if candidate is None else candidate.get("repeat_bases"),
                "previous_max_tract_bp": None if candidate is None else candidate.get("max_tract"),
                "previous_result_sha256": previous_sha,
            })
    return rows


def example_rows(groups, sample_id):
    rows = []
    for group in groups:
        for example in group["examples"]:
            start, end, length = example.get("start"), example.get("end"), example.get("record_length")
            for key, value in (("start", start), ("end", end), ("record_length", length)):
                nonnegative_integer(value, f"{group['name']}.example.{key}")
            if not 0 <= start < end <= length:
                raise ValueError(f"{group['name']}: invalid 0-based half-open example coordinates")
            rows.append({
                "sample_id": sample_id, "group": group["name"], "record_id": example.get("record_id"),
                "record_length": length, "input_index": example.get("input_index"),
                "start": start, "end": end, "span_bp": end - start,
                "strand": "G" if example.get("strand") is True else "C" if example.get("strand") is False
                else "ambiguous" if example.get("strand") is None else example["strand"], "touches_record_start": start == 0,
                "touches_record_end": end == length,
                "biological_length_censored": start == 0 or end == length,
                "unit_counts": example.get("unit_counts"), "native_example": example,
            })
    return rows


def build_report(source_path, previous_path=None, sample_id=None):
    source_path = Path(source_path).resolve()
    native = validate_source(json.loads(source_path.read_text()))
    sample_id = sample_id or source_path.parent.name
    block = native["mixed_arrays"]
    groups = []
    for group in block["groups"]:
        groups.append({
            "sample_id": sample_id, "group": group["name"],
            "units_g_oriented": group["units"],
            "observation": "exact_mixed_array_observed" if group["mixed_records"] else "supplied_exact_mixed_pattern_not_observed",
            "input_records": native["summary"]["records"], "input_bases": native["summary"]["bases"],
            **{key: group[key] for key in GROUP_COUNTS},
            "expected_boundary_records_sum": group["g_right_boundary_records"] + group["c_left_boundary_records"],
            "opposite_boundary_records_sum": group["g_left_boundary_records"] + group["c_right_boundary_records"],
            "examples_saved": len(group["examples"]),
        })
    report = {
        "schema_version": 1, "report_version": REPORT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(), "sample_id": sample_id,
        "analysis_kind": "targeted_mixed_array_observation_report",
        "source": {"path": str(source_path), "sha256": digest(source_path), "version": native.get("version")},
        "mode": block.get("mode", native.get("config", {}).get("mode", "unspecified")),
        "threads": native.get("threads"),
        "summary": native["summary"], "inputs": native.get("inputs", []),
        "settings": block["settings"], "groups": groups,
        "native_method": block.get("method", {}),
        "definitions_source": native.get("definitions_source"),
        "examples": example_rows(block["groups"], sample_id), "previous_components": [],
        "classification_changed": False, "data_adequacy_changed": False,
        "biological_call": None, "limitations": LIMITATIONS,
        "coordinate_convention": "0-based half-open [start,end)",
        "boundary_sum_interpretation": "The sums add record counts from two orientation categories. One record can occur in more than one category; these are not unique molecule counts.",
        "span_interpretation": "Maximum observed single-record exact mixed span; read-edge observations can be censored. Not an inferred complete chromosome-telomere length.",
        "native_warnings": native.get("warnings", []),
    }
    if previous_path:
        previous_path = Path(previous_path).resolve()
        previous = json.loads(previous_path.read_text())
        report["previous"] = {"path": str(previous_path), "sha256": digest(previous_path),
                              "version": previous.get("version"),
                              "summary": previous.get("summary", {}),
                              "threads": previous.get("threads"),
                              "input_check": check_previous_inputs(native, previous)}
        report["previous_components"] = component_rows(block["groups"], previous, sample_id, report["previous"]["sha256"])
    return report


def text_value(value):
    if value is None:
        return ""
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    if isinstance(value, bool):
        return str(value).lower()
    return str(value)


def write_tsv(path, fields, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({key: text_value(row.get(key)) for key in fields} for row in rows)


def table(headers, rows):
    return "<table><thead><tr>" + "".join(f"<th>{html.escape(str(x))}</th>" for x in headers) + "</tr></thead><tbody>" + "".join(
        "<tr>" + "".join(f"<td>{html.escape(text_value(x))}</td>" for x in row) + "</tr>" for row in rows) + "</tbody></table>"


def render_html(report):
    esc = html.escape
    groups = report["groups"]
    full_percent = f"{report['settings']['full_fraction']:.0%}"
    primary = table(
        ["Group", "명시한 G 방향 단위", "혼합 records", f"≥{full_percent} 덮인 records", "100% 덮인 records", "최대 관측 span (bp)", "긴 배열 records", "혼합 기준 충족 없음·순수 단위 관측", "그 외 기준 미달 records"],
        [(g["group"], " + ".join(g["units_g_oriented"]), g["mixed_records"], g["full_mixed_records"],
          g["entire_mixed_records"], g["max_array_bp"], g["long_array_records"], g["pure_only_records"],
          g["nonqualifying_records"]) for g in groups])
    geometry = table(
        ["Group", "G → 오른쪽 경계", "C → 왼쪽 경계", "G → 왼쪽 (반대)", "C → 오른쪽 (반대)", "양끝에서 떨어진 내부 배열 records", "양쪽 clean flank records"],
        [(g["group"], g["g_right_boundary_records"], g["c_left_boundary_records"], g["g_left_boundary_records"],
          g["c_right_boundary_records"], g["internal_array_records"], g["clean_two_flank_records"]) for g in groups])
    endpoints = table(
        ["Group", "C-left 끝 근처", "G-right 끝 근처", "G-left 끝 근처", "C-right 끝 근처", "긴 배열 + 예상 말단 방향 records"],
        [(g["group"], g["left_c_end_records"], g["right_g_end_records"], g["left_g_end_records"],
          g["right_c_end_records"], g["long_expected_end_records"]) for g in groups])
    ambiguity = table(
        ["Group", "단위 분해 모호 records (전체 관측)", "방향 모호 records (전체 관측)", "방향 모호 경계 records (혼합)"],
        [(g["group"], g["segmentation_ambiguous_records"], g["ambiguous_orientation_records"],
          g["ambiguous_boundary_records"]) for g in groups])
    components = ""
    if report["previous_components"]:
        components = "<h2>이전 단일 주기 탐색과 비교</h2><p>이전 각 family 계수는 서로 겹칩니다. 혼합 배열의 고유 records 수와 합산하지 않습니다. 미보고는 0개로 확인한 검사가 아닙니다.</p>" + table(
            ["Group", "단위", "이전 관측 상태", "이전 지지 records", "이전 full-repeat records", "이전 최대 단일 주기 tract (bp)"],
            [(r["group"], r["unit"], r["previous_observation"], r["previous_reads"], r["previous_full_repeat_reads"], r["previous_max_tract_bp"]) for r in report["previous_components"]])
    settings = esc(json.dumps(report["settings"], ensure_ascii=False, indent=2))
    limits = "".join(f"<li>{esc(x)}</li>" for x in report["limitations"])
    timing = ""
    elapsed = report["summary"].get("elapsed_seconds")
    if elapsed is not None:
        timing = f"<p>Native scanner가 기록한 분석 시간: <strong>{esc(text_value(elapsed))}초</strong> · threads: {esc(text_value(report['threads']))}. 다운로드와 보고서 생성 시간은 포함하지 않습니다.</p>"
    return f"""<!doctype html>
<html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(report['sample_id'])} — 혼합 반복 구조 관측</title>
<style>body{{font:16px/1.6 system-ui,sans-serif;max-width:1500px;margin:auto;padding:24px;color:#182b38}}table{{border-collapse:collapse;display:block;overflow:auto;margin:16px 0}}th,td{{padding:9px;border:1px solid #ccd6df;text-align:left;vertical-align:top}}th{{background:#eaf0f5}}pre{{white-space:pre-wrap;background:#f4f6f8;padding:16px}}.notice{{padding:16px;background:#fff5d9;border-left:4px solid #b88010}}a{{color:#075da0}}</style>
<h1>{esc(report['sample_id'])}: 혼합 반복 구조 관측</h1>
<p class="notice">명시한 단위 집합에 대한 구조 검사입니다. 기존 v0.6 TRM 판정과 데이터 충분성 평가는 변경하지 않습니다.</p>
<p>전체 입력 <strong>{report['summary']['records']:,} records / {report['summary']['bases']:,} bp</strong> · mode: {esc(report['mode'])} · native: {esc(str(report['source']['version']))}</p>
{timing}
<p><a href="mixed_arrays.tsv">Group별 관측 TSV</a> · <a href="previous_components.tsv">이전 component 비교 TSV</a> · <a href="examples.tsv">제한된 좌표 예시 TSV</a> · <a href="report.json">전체 JSON·출처</a></p>
<h2>관측된 혼합 배열</h2>{primary}
<p>Full mixed는 하나의 혼합 배열이 record의 {full_percent} 이상을 덮는 경우이며, 100% coverage는 별도 열입니다. 최대 span은 단일 record의 관측 구간입니다. Read 끝에 닿는 배열은 길이가 잘릴 수 있으며, 150 bp read 하나로 500 bp telomere를 확인할 수 없습니다. 긴 배열의 bp 기준은 아래 실행 설정에 기록됩니다. 미관측 결과는 이 exact 혼합 규칙의 미관측입니다.</p>
<p>“그 외 기준 미달”은 순수 배열도, 기준을 충족한 혼합 배열도 보고되지 않은 records로, component copy 수 부족 또는 모호한 분해를 포함합니다. “순수 단위 관측” records에도 기준 미달 혼합 구간이 함께 있을 수 있습니다.</p>
<h2>방향성과 내부 배열</h2>{geometry}
<p>G는 제공한 G 방향 단위, C는 그 reverse complement입니다. Clean flank와 record 경계를 만족한 G-right / C-left를 예상 방향으로, 나머지를 반대 방향으로 따로 표시합니다. Full mixed records는 경계 집계에서 제외됩니다. 방향별 열은 서로 겹칠 수 있습니다. Record 끝과 물리적 염색체 말단은 구분해야 합니다.</p>
<h2>Record 끝 근처의 방향</h2>{endpoints}
<p>이 표는 end_window_bp 안에 있는 구간을 세며 clean flank를 요구하지 않습니다. 위 경계 표와 조건이 다릅니다. Assembly mode에서도 count 단위는 contig/chromosome record이며 독립 말단의 수가 아닙니다.</p>
<h2>모호성 진단</h2>{ambiguity}
{components}<h2>검사 설정</h2><pre>{settings}</pre>
<h2>해석 범위</h2><ul>{limits}</ul>
<p>Source SHA-256: <code>{esc(report['source']['sha256'])}</code> · Report {REPORT_VERSION}</p>
</html>"""


def write_report(report, output):
    output = Path(output)
    output_files = {(output / name).resolve() for name in (
        "report.json", "mixed_arrays.tsv", "previous_components.tsv", "examples.tsv", "index.html")}
    source_files = {Path(report["source"]["path"]).resolve()}
    if report.get("previous"):
        source_files.add(Path(report["previous"]["path"]).resolve())
    if source_files & output_files:
        raise ValueError("Output would overwrite a source result")
    output.mkdir(parents=True, exist_ok=True)
    (output / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    write_tsv(output / "mixed_arrays.tsv", GROUP_FIELDS, report["groups"])
    write_tsv(output / "previous_components.tsv", COMPONENT_FIELDS, report["previous_components"])
    write_tsv(output / "examples.tsv", EXAMPLE_FIELDS, report["examples"])
    (output / "index.html").write_text(render_html(report))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="Native targeted_mixed_arrays JSON")
    parser.add_argument("--previous", type=Path, help="Previous homogeneous native JSON from the same ordered inputs")
    parser.add_argument("--sample-id", help="Display label; defaults to the source parent directory name")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        report = build_report(args.source, args.previous, args.sample_id)
        write_report(report, args.output)
    except (ValueError, KeyError, TypeError, OSError) as exc:
        parser.error(str(exc))
    print(json.dumps({"sample_id": report["sample_id"], "groups": len(report["groups"]),
                      "output": str(args.output.resolve()), "classification_changed": False}, ensure_ascii=False))


if __name__ == "__main__":
    main()
