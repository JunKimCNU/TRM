#!/usr/bin/env python3
"""Report optional native flank-period annotations without reading sequence inputs."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import csv
import hashlib
import html
import json
from pathlib import Path
import re
import sys

VERSION = '0.1.0'
CLASSES = ('insufficient_window', 'neither', 'lag_only', 'family_only', 'both')
QUADRANTS = ('left_forward', 'left_reverse', 'right_forward', 'right_reverse', 'ambiguous')
LABELS = {
    'insufficient_window': '검사 구간 부족',
    'neither': '두 기준 모두 미충족',
    'lag_only': '주기 기준만 충족',
    'family_only': 'motif 기준만 충족',
    'both': '두 기준 모두 충족',
}
STATES = {'unmeasured': '미측정', 'partial': '일부 사건만 측정',
          'measured_zero': '측정됨 · 대상 사건 0', 'measured': '대상 사건 모두 측정'}
BASE_FIELDS = ['source_index', 'sample_label', 'source_path', 'source_sha256',
               'analysis_kind', 'candidate_container', 'motif', 'period_bases',
               'native_candidate_status', 'native_reads', 'native_junction_reads',
               'measurement_status', 'measured_all_events', 'required_window_bases',
               'threshold_percent', 'observed_events', 'evaluated_events']
FIELDS = BASE_FIELDS + list(CLASSES) + [f'{name}_{q}' for name in CLASSES for q in QUADRANTS]


def integer(value, name, minimum=0, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise ValueError(f'{name}: expected integer in permitted range')
    return value


def canonical(motif):
    if not isinstance(motif, str) or not 1 <= len(motif) <= 64 or set(motif) - set('ACGT'):
        raise ValueError('motif: expected 1..64 uppercase A/C/G/T bases')
    for period in range(1, len(motif)):
        if len(motif) % period == 0 and motif == motif[:period] * (len(motif) // period):
            raise ValueError('motif: expected primitive repeat unit')
    reverse = motif.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    return min(s[offset:] + s[:offset] for s in (motif, reverse) for offset in range(len(s)))


def object_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('duplicate JSON key: ' + key)
        result[key] = value
    return result


def reject_constant(value):
    raise ValueError('non-finite JSON number: ' + value)


def snapshot(path):
    s = path.stat()
    return s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns


def load_source(path):
    path = Path(path).resolve(strict=True)
    if not path.is_file():
        raise ValueError('input must be a regular JSON file: ' + str(path))
    before = snapshot(path)
    raw = path.read_bytes()
    if snapshot(path) != before:
        raise ValueError('input changed while reading: ' + str(path))
    value = json.loads(raw, object_pairs_hook=object_pairs, parse_constant=reject_constant)
    if not isinstance(value, dict):
        raise ValueError('native result must be a JSON object: ' + str(path))
    return value, {'path': str(path), 'size_bytes': len(raw),
                   'sha256': hashlib.sha256(raw).hexdigest()}


def validate_annotation(annotation, motif, junction_reads):
    """Validate stored arithmetic, not the original sequence-level calculation."""
    if not isinstance(annotation, dict):
        raise ValueError('flank_period_annotation must be an object or null')
    if integer(annotation.get('schema_version'), 'annotation schema', 1) != 1:
        raise ValueError('unsupported flank annotation schema')
    if annotation.get('candidate_motif') != motif:
        raise ValueError('annotation candidate_motif differs from native family')
    period = integer(annotation.get('primitive_period_bases'), 'annotation period', 1, 64)
    if period != len(motif):
        raise ValueError('annotation period differs from primitive motif length')
    window = integer(annotation.get('required_window_bases'), 'annotation window', 20, 192)
    if window != max(20, 3 * period):
        raise ValueError('annotation window must equal max(20,3*k)')
    if integer(annotation.get('threshold_percent'), 'annotation threshold') != 85:
        raise ValueError('annotation threshold must be 85')
    measured = annotation.get('measured_all_events')
    if type(measured) is not bool:
        raise ValueError('measured_all_events must be a boolean')
    observed = integer(annotation.get('observed_events'), 'observed_events')
    classes = annotation.get('classes')
    if not isinstance(classes, dict) or set(classes) != set(CLASSES):
        raise ValueError('annotation requires exactly the five supported classes')
    total = 0
    for name in CLASSES:
        category = classes[name]
        if not isinstance(category, dict):
            raise ValueError('annotation class must be an object: ' + name)
        events = integer(category.get('events'), name + ' events')
        counts = category.get('quadrant_counts')
        if not isinstance(counts, dict) or set(counts) != set(QUADRANTS):
            raise ValueError('annotation class requires five raw quadrant counts: ' + name)
        if sum(integer(counts[q], name + '.' + q) for q in QUADRANTS) != events:
            raise ValueError('quadrant sum differs from class events: ' + name)
        total += events
    if total != observed:
        raise ValueError('class sum differs from observed_events')
    if observed and not junction_reads:
        raise ValueError('annotated events require at least one native junction read')
    if measured and observed < junction_reads:
        raise ValueError('complete event annotation cannot contain fewer events than junction reads')
    return annotation


def sample_label(document, source):
    """A display hint, never a claim of independent biological sample identity."""
    for key in ('sample_id', 'sample'):
        if isinstance(document.get(key), str) and document[key]:
            return document[key]
    labels = set()
    inputs = document.get('inputs', [])
    if isinstance(inputs, list):
        for item in inputs:
            path = item.get('path', '') if isinstance(item, dict) else item
            if isinstance(path, str):
                labels.update(re.findall(r'(?<![A-Za-z0-9])(?:SRR|ERR|DRR)\d+(?!\d)', path))
    return ','.join(sorted(labels)) if labels else Path(source['path']).parent.name + '/' + Path(source['path']).name


def result_rows(document, source, source_index):
    kinds = [name for name in ('candidates', 'candidate_evidence') if name in document]
    if len(kinds) != 1:
        raise ValueError('expected exactly one native candidates or candidate_evidence array')
    containers = kinds + (['rejected_candidates'] if 'rejected_candidates' in document else [])
    seen, rows = set(), []
    for container in containers:
        candidates = document[container]
        if not isinstance(candidates, list):
            raise ValueError(container + ' must be an array')
        for candidate in candidates:
            if not isinstance(candidate, dict):
                raise ValueError('candidate must be an object')
            motif = candidate.get('motif')
            family = canonical(motif)
            if motif != family:
                raise ValueError('native candidate motif must be canonical')
            if family in seen:
                raise ValueError('duplicate family across native candidate containers: ' + family)
            seen.add(family)
            if 'period' in candidate and integer(candidate['period'], 'candidate period', 1, 64) != len(motif):
                raise ValueError('native candidate period differs from motif')
            reads = integer(candidate.get('reads'), 'native candidate reads')
            junction = candidate.get('junction_evidence')
            if not isinstance(junction, dict):
                raise ValueError('native candidate junction_evidence object required')
            junction_reads = integer(junction.get('junction_reads'), 'native junction_reads')
            if junction_reads > reads:
                raise ValueError('junction_reads exceeds native candidate reads')
            row = dict.fromkeys(FIELDS)
            row.update(source_index=source_index, sample_label=sample_label(document, source),
                       source_path=source['path'], source_sha256=source['sha256'],
                       analysis_kind=document.get('analysis_kind', 'native_scan'),
                       candidate_container=container, motif=motif, period_bases=len(motif),
                       native_candidate_status=candidate.get('status', ''),
                       native_reads=reads, native_junction_reads=junction_reads,
                       measurement_status='unmeasured')
            annotation = junction.get('flank_period_annotation')
            if annotation is not None:
                annotation = validate_annotation(annotation, motif, junction_reads)
                observed = annotation['observed_events']
                complete = annotation['measured_all_events']
                row.update(measurement_status='partial' if not complete else ('measured' if observed else 'measured_zero'),
                           measured_all_events=complete, observed_events=observed,
                           evaluated_events=observed - annotation['classes']['insufficient_window']['events'],
                           required_window_bases=annotation['required_window_bases'], threshold_percent=85)
                for name in CLASSES:
                    row[name] = annotation['classes'][name]['events']
                    for q in QUADRANTS:
                        row[name + '_' + q] = annotation['classes'][name]['quadrant_counts'][q]
            rows.append(row)
    return rows


def render(rows, sources):
    esc = lambda value: html.escape(str(value), quote=True)
    number = lambda value: '—' if value is None else f'{value:,}'
    body = []
    for row in rows:
        classes = ''.join('<td>' + number(row[name]) + '</td>' for name in CLASSES)
        detail = ['<table><tr><th>분류</th>' + ''.join('<th>' + esc(q) + '</th>' for q in QUADRANTS) + '</tr>']
        for name in CLASSES:
            detail.append('<tr><th>' + LABELS[name] + '</th>' + ''.join('<td>' + number(row[name + '_' + q]) + '</td>' for q in QUADRANTS) + '</tr>')
        detail.append('</table>')
        title = esc(row['sample_label']) + ' · ' + esc(row['motif'])
        body.append('<tr><td>' + title + '<br><small>입력 #' + str(row['source_index']) + ' · ' + esc(row['candidate_container']) + '</small></td>'
                    + '<td>' + STATES[row['measurement_status']] + '</td>'
                    + '<td>' + number(row['native_junction_reads']) + '</td><td>' + number(row['observed_events']) + '</td><td>' + number(row['required_window_bases']) + '</td>' + classes + '</tr>'
                    + '<tr><td colspan="10"><details><summary>' + title + ' — raw canonical 방향별 사건 수</summary>' + ''.join(detail) + '</details></td></tr>')
    source_items = ''.join('<li><strong>입력 #' + str(index) + '</strong> ' + esc(s['path'])
                           + '<br><code>' + s['sha256'] + '</code></li>' for index, s in enumerate(sources, 1))
    return '''<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>TRM 경계의 반복 주기 주석</title><style>body{font:16px/1.6 system-ui,sans-serif;max-width:1500px;margin:40px auto;padding:0 24px;color:#203047}h1{font-size:28px}table{border-collapse:collapse;width:100%;margin:18px 0;font-size:14px}td,th{padding:9px;border:1px solid #d9e0e7;text-align:right;vertical-align:top}td:first-child,th:first-child{text-align:left}th{background:#eef3f7}small{color:#536579}.table{overflow-x:auto}.note{padding:14px 18px;background:#f0f5f8;border-left:4px solid #3a6b8d}code{overflow-wrap:anywhere}li{margin:12px 0}summary{cursor:pointer;text-align:left}a{color:#216895}</style>
<h1>TRM 경계의 반복 주기 주석</h1><p>기존 native scan/evidence 결과의 선택적 설명 자료입니다. 원래 후보·read 계수·TRM 판정은 바꾸지 않습니다.</p>
<div class="note"><strong>‘두 기준 모두 미충족’은 비반복 subtelomere 또는 고유 locus의 증명이 아닙니다.</strong> 구간 부족은 음성으로 분류하지 않습니다. 누락된 주석은 미측정이며, ‘대상 사건 0’과 구분합니다. 일부 측정 결과는 전체 사건의 분포로 해석하지 마세요.</div>
<p>한 read에서 여러 경계 사건이 생길 수 있습니다. native junction reads와 annotation events는 서로 다른 단위입니다. 사건은 저장된 bottom-128 flank 표본에서 추정하지 않고 native의 전체 적격 경계에서 집계합니다. 같은 입력의 여러 분석 결과는 독립 표본으로 합산하지 않습니다.</p>
<div class="table"><table><thead><tr><th>표시 이름 · family</th><th>측정 범위</th><th>기존 junction reads</th><th>주석 events</th><th>검사 bp</th>''' + ''.join('<th>' + LABELS[name] + '</th>' for name in CLASSES) + '</tr></thead><tbody>' + ''.join(body) + '''</tbody></table></div>
<p>검사 길이 W=max(20,3k), k=primitive motif 길이. Lag는 k bp 떨어진 염기의 일치율, family는 motif의 모든 rotation과 reverse complement 중 최대 ungapped 일치율이며 각각 기준은 85%입니다. 후자는 위상을 자유롭게 선택하므로 tract와의 같은 위상 연속성을 입증하지 않습니다. 방향은 canonical family 기준이며 생물학적 G/C 방향으로 변환하지 않았습니다. ambiguous는 방향을 정할 수 없는 사건입니다.</p>
<h2>입력 근거</h2><p>표시 이름은 metadata 또는 accession 파일명에서 얻은 힌트입니다. 생물학적 독립성은 검증하지 않습니다. 입력 SHA256은 JSON 파일에 대한 값이며 WGS 파일을 다시 검사한 값이 아닙니다.</p><ol>''' + source_items + '''</ol><p><a href="samples.tsv">전체 TSV</a> · <a href="summary.json">설정·입력 해시</a> · <a href="README.md">설명</a></p></html>'''


def build_report(paths, output):
    output = Path(output).resolve()
    if output.exists():
        raise ValueError('output directory already exists; choose a new path')
    rows, sources, seen = [], [], set()
    for path in paths:
        document, source = load_source(path)
        if source['path'] in seen:
            raise ValueError('duplicate input path: ' + source['path'])
        seen.add(source['path'])
        sources.append(source)
        rows.extend(result_rows(document, source, len(sources)))
    # Validate every source before committing any report output.
    page = render(rows, sources)
    summary = {'schema_version': 1, 'reporter': 'report_flank_period', 'version': VERSION,
               'created_utc': datetime.now(timezone.utc).isoformat(), 'status': 'complete',
               'input_result_files': len(sources), 'candidate_result_rows': len(rows),
               'measurement_status_counts': dict(sorted(Counter(r['measurement_status'] for r in rows).items())),
               'source_files': sources, 'changes_native_classification': False,
               'reads_sequence_inputs': False, 'validates_stored_arithmetic_only': True,
               'source_hash_scope': 'Native result JSON bytes; original WGS is not reopened.',
               'independent_biological_sample_count': None,
               'threshold_percent': 85, 'required_window': 'max(20,3*primitive_period_bases)',
               'classes': LABELS, 'quadrants': list(QUADRANTS),
               'notes': ['Missing annotation is unmeasured, not zero.',
                         'Partial annotation cannot represent all eligible events.',
                         'Neither threshold passed does not establish a nonrepetitive or unique flank.',
                         'Family compatibility allows all phases and both strands, not same-phase continuation.',
                         'Events and junction reads are different units; no cross-result biological pooling.']}
    output.mkdir(parents=True, exist_ok=False)
    with (output / 'samples.tsv').open('x', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    (output / 'index.html').write_text(page, encoding='utf-8')
    (output / 'README.md').write_text(
        '# TRM 경계의 반복 주기 주석\n\n'
        '[읽기 쉬운 결과](index.html) · [전체 TSV](samples.tsv) · [입력 해시와 범위](summary.json)\n\n'
        f'Native 결과 {len(sources)}개 파일의 family 결과 {len(rows)}행을 정리했다. '
        '동일 입력의 여러 분석을 독립 표본으로 합산하지 않았다.\n\n'
        '기존 후보·read 수·판정은 유지한다. Annotation events는 적격 tract-side 경계의 수이며 '
        'junction reads와 단위가 다르다. 저장된 flank 표본으로 전체를 추정하지 않는다.\n\n'
        '검사 길이는 `W=max(20,3*k)`, k는 primitive motif 길이다. read 끝·N 등 비ACGT·다른 '
        'same-family tract에 막혀 W를 확보하지 못하면 `insufficient_window`다. '
        'Lag-k 일치율과 모든 phase/RC 중 최대 motif 일치율을 각각 85% 기준으로 비교한다. '
        '`lag_only`, `family_only`, `both`, `neither`는 그 두 조건의 조합이다. '
        '`neither`는 비반복 subtelomere·고유 locus의 증거가 아니며, family 일치는 같은 위상 연속성을 입증하지 않는다.\n\n'
        '미측정 값은 TSV 빈칸/HTML 대시이며 실제 0과 구분한다. `measured_all_events=false`는 일부 사건만 '
        '측정했다는 뜻이다. 방향은 canonical family 기준이며 ambiguous는 방향 불명이다. '
        '이 보고서는 JSON 내부 합계와 형식만 검사한다. 원시 서열, WGS 충분성, 물리적 염색체 말단은 검증하지 않는다.\n',
        encoding='utf-8')
    summary['outputs'] = []
    for name in ('samples.tsv', 'index.html', 'README.md'):
        raw = (output / name).read_bytes()
        summary['outputs'].append({'path': str(output / name), 'size_bytes': len(raw),
                                   'sha256': hashlib.sha256(raw).hexdigest()})
    (output / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path, nargs='+', help='native scan/evidence JSON result(s)')
    parser.add_argument('--output', type=Path, required=True, help='new report directory (must not exist)')
    args = parser.parse_args(argv)
    try:
        summary = build_report(args.results, args.output)
    except (OSError, ValueError, TypeError) as error:
        parser.exit(2, 'error: ' + str(error) + '\n')
    print(json.dumps({'status': summary['status'], 'candidate_result_rows': summary['candidate_result_rows'],
                      'output': str(args.output.resolve())}, ensure_ascii=False))


if __name__ == '__main__':
    main()
