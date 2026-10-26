#!/usr/bin/env python3
"""Report exact mixed-block compatibility of saved native flank representatives."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import html
import json
from pathlib import Path
import time

from evidence_io import EvidenceIO
from report_flank_period import canonical, integer, result_rows, sample_label, QUADRANTS
from saved_flank_grammar import analyze, checked_units, VERSION as ALGORITHM_VERSION

VERSION = '0.1.0'
IO = EvidenceIO(ValueError, 'saved flank report: ')
MAX_SOURCE_BYTES = 64 * 1024**2
MAX_SOURCES = 128
MAX_QUERY_ROWS = 32768
SUMMARY_FIELDS = ['source_index', 'sample_label', 'source_path', 'source_sha256', 'native_version',
    'source_sampled_prefix', 'family', 'candidate_container', 'native_reads', 'native_junction_reads',
    'native_candidate_status', 'measurement_status', 'source_validation_measured', 'source_sketch_saturated',
    'source_observed_anchor_events', 'source_distinct_anchor_lower_bound', 'retained_anchors',
    'representatives_evaluated', 'representatives_unavailable', 'full_grammar_compatible_queries',
    'mixed_only_queries', 'queries_with_ge20_outside_target_grammar']
QUERY_FIELDS = ['source_index', 'sample_label', 'source_path', 'source_sha256', 'family', 'source_json_pointer',
    'anchor_signature', 'anchor20', 'source_anchor_events', 'source_raw_quadrants', 'source_normalized_quadrants',
    'measurement_status', 'sequence', 'sequence_sha256', 'sequence_bases', 'full_grammar_compatible',
    'full_single_unit_compatible', 'full_compatible_units', 'full_grammar_strands', 'compatibility_class',
    'longest_grammar_span_bp', 'minimum_union_span_bp', 'grammar_union_bp', 'prefix_outside_bp',
    'suffix_outside_bp', 'max_contiguous_outside_bp', 'intervals']
LIMITS = [
    'Only retained representative sequences are evaluated; not all boundary events or raw WGS.',
    'Anchor event counts are not support counts for the longest representative, independent molecules or loci.',
    'Outside this exact target grammar does not establish a nonrepeat, unique locus or subtelomere.',
    'Partial query-edge blocks and either strand are allowed; original tract phase continuity is not established.',
    'No read ID, mate, per-base quality, physical chromosome end or biological independence is inferred.',
    'Saturated sketches omit contexts; observed proportions are not error rates, sensitivity or biological prevalence.',
    'The declared literal units are a supplied model, not automatically inferred canonical TRMs.',
    'No source count, ranking, software classification or adequacy is changed.',
]


def _bool(value, name):
    if type(value) is not bool:
        raise ValueError(name + ' must be boolean')
    return value


def _quadrants(value, events):
    if not isinstance(value, dict) or set(value) != set(QUADRANTS):
        raise ValueError('five quadrant counts required')
    if sum(integer(value[q], q) for q in QUADRANTS) != events:
        raise ValueError('quadrant sum differs from anchor events')
    return value


def _entries(sketch):
    if not isinstance(sketch, dict):
        raise ValueError('flank_contexts must be an object')
    if integer(sketch.get('capacity'), 'capacity') != 128 or integer(sketch.get('flank_bases'), 'flank_bases') != 20:
        raise ValueError('only native128-entry/20bp identity schema is supported')
    observed = integer(sketch.get('observed_events'), 'observed_events')
    lower = integer(sketch.get('distinct_lower_bound'), 'distinct_lower_bound')
    saturated = _bool(sketch.get('saturated'), 'saturated')
    entries = sketch.get('entries')
    if not isinstance(entries, list) or len(entries) > 128:
        raise ValueError('entries must be a bounded native sketch list')
    if lower != len(entries) + int(saturated) or (saturated and len(entries) != 128):
        raise ValueError('sketch lower bound/capacity is inconsistent')
    seen, keys, total = set(), [], 0
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError('entry must be an object')
        anchor = entry.get('sequence')
        if not isinstance(anchor, str) or len(anchor) != 20 or set(anchor) - set('ACGT'):
            raise ValueError('anchor must be20 uppercase ACGT bases')
        packed, rank = 0, 0x6c62272e07bb014262b821756295c58d
        for base in anchor:
            packed = (packed << 2) | 'ACGT'.index(base)
            rank = ((rank ^ ord(base)) * 0x0000000001000000000000000000013b) & ((1 << 128)-1)
        if entry.get('signature') != f'{packed:010x}' or entry.get('rank') != f'{rank:032x}':
            raise ValueError('anchor signature/rank differs from sequence')
        if anchor in seen:
            raise ValueError('duplicate retained anchor')
        seen.add(anchor); keys.append((entry['rank'], entry['signature']))
        events = integer(entry.get('events'), 'anchor events', 1); total += events
        raw = _quadrants(entry.get('quadrant_counts'), events)
        normalized = _quadrants(entry.get('normalized_quadrant_counts'), events)
        expected = dict(left_forward=raw['left_forward']+raw['right_reverse'],
            left_reverse=raw['left_reverse']+raw['right_forward'], right_forward=0, right_reverse=0,
            ambiguous=raw['ambiguous'])
        if normalized != expected:
            raise ValueError('normalized quadrant geometry differs from raw counts')
        seq = entry.get('representative_flank')
        if seq is not None and (not isinstance(seq, str) or not 20 <= len(seq) <= 80
                or set(seq)-set('ACGT') or not seq.startswith(anchor)):
            raise ValueError('representative must be20..80 ACGT bases starting with anchor')
    if keys != sorted(keys) or (not saturated and total != observed) or (saturated and total >= observed):
        raise ValueError('sketch order/observed event accounting is inconsistent')
    return entries


def extract_rows(document, source, index, units, families, min_span):
    if (document.get('program') != 'trm' or type(document.get('schema_version')) is not int
            or document['schema_version'] != 4 or document.get('version') not in ('0.8.2','0.8.3','0.8.4','0.8.5')):
        raise ValueError('supported native source is trm0.8.2..0.8.5 schema4')
    source_summary = document.get('summary')
    if source_summary is not None and not isinstance(source_summary, dict):
        raise ValueError('native summary must be an object')
    prefix = (source_summary or {}).get('sampled_prefix')
    if prefix is not None:
        _bool(prefix, 'sampled_prefix')
    native = {r['motif']: r for r in result_rows(document, source, index)}
    candidates = {c['motif']: (container, i, c) for container in ('candidates', 'candidate_evidence', 'rejected_candidates')
                  for i, c in enumerate(document.get(container, []))}
    summaries, queries = [], []
    for family in families:
        row = dict.fromkeys(SUMMARY_FIELDS)
        row.update(source_index=index, sample_label=sample_label(document, source), source_path=source['path'],
            source_sha256=source['sha256'], native_version=document['version'], source_sampled_prefix=prefix,
            family=family, measurement_status='candidate_not_reported')
        summaries.append(row)
        if family not in native:
            continue
        n = native[family]; container, ci, candidate = candidates[family]
        for key in ('candidate_container','native_reads','native_junction_reads','native_candidate_status'):
            row[key] = n[key]
        row['measurement_status'] = 'unmeasured'
        validation = candidate['junction_evidence'].get('validation')
        if validation is None:
            continue
        if not isinstance(validation, dict) or integer(validation.get('schema_version'), 'validation schema') != 1:
            raise ValueError('only native validation schema1 is supported')
        measured = validation.get('measured', False)
        row['source_validation_measured'] = _bool(measured, 'validation measured')
        sketch = validation.get('flank_contexts')
        if sketch is None:
            row['measurement_status'] = 'sketch_unavailable' if measured else 'unmeasured'
            continue
        entries = _entries(sketch)
        row.update(source_sketch_saturated=sketch['saturated'], source_observed_anchor_events=sketch['observed_events'],
            source_distinct_anchor_lower_bound=sketch['distinct_lower_bound'], retained_anchors=len(entries))
        if not measured:
            continue  # A merge can contain some saved contexts without complete measurement.
        row.update(representatives_evaluated=0, representatives_unavailable=0,
            full_grammar_compatible_queries=0, mixed_only_queries=0, queries_with_ge20_outside_target_grammar=0)
        for ei, entry in enumerate(entries):
            q = dict.fromkeys(QUERY_FIELDS)
            q.update(source_index=index, sample_label=row['sample_label'], source_path=source['path'], source_sha256=source['sha256'],
                family=family, source_json_pointer=f'/{container}/{ci}/junction_evidence/validation/flank_contexts/entries/{ei}',
                anchor_signature=entry['signature'], anchor20=entry['sequence'], source_anchor_events=entry['events'],
                source_raw_quadrants=entry['quadrant_counts'], source_normalized_quadrants=entry['normalized_quadrant_counts'],
                measurement_status='representative_unavailable')
            seq = entry.get('representative_flank')
            if seq is None:
                row['representatives_unavailable'] += 1
            else:
                q.update(analyze(seq, units, min_span))
                q.update(sequence=seq, sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(), measurement_status='evaluated_saved_representative')
                row['representatives_evaluated'] += 1
                row['full_grammar_compatible_queries'] += int(q['full_grammar_compatible'])
                row['mixed_only_queries'] += int(q['compatibility_class'] == 'mixed_only')
                row['queries_with_ge20_outside_target_grammar'] += int(q['max_contiguous_outside_bp'] >= 20)
            queries.append(q)
        row['measurement_status'] = ('measured_zero_retained' if not entries else
            'representatives_incomplete' if row['representatives_unavailable'] else 'evaluated_retained_subset')
        if entries and not row['representatives_evaluated']:
            for key in ('full_grammar_compatible_queries', 'mixed_only_queries', 'queries_with_ge20_outside_target_grammar'):
                row[key] = None
    return summaries, queries


def build_report(paths, output, units, families, min_span=20):
    start, cpu = time.monotonic(), time.process_time()
    units = list(checked_units(units))
    integer(min_span, 'minimum union span', 1, 80)
    if not families or len(families) > 32:
        raise ValueError('select1..32 candidate families explicitly')
    families = [canonical(f) for f in families]
    if len(set(families)) != len(families):
        raise ValueError('duplicate canonical selected family')
    if not paths or len(paths) > MAX_SOURCES:
        raise ValueError('select1..128 input JSON files')
    output = Path(output).resolve()
    if output.exists():
        raise ValueError('output exists; choose a new directory')
    sources, summaries, queries, seen = [], [], [], set()
    for path in paths:
        path = Path(path).resolve(strict=True)
        if path in seen:
            raise ValueError('duplicate source path')
        seen.add(path)
        raw = IO.read_bytes(path, max_bytes=MAX_SOURCE_BYTES)
        document, sha = IO.decode_json(raw, path, object_pairs_hook=IO.unique_object, parse_constant=IO.reject_constant)
        source = dict(path=str(path), size_bytes=len(raw), sha256=sha)
        sources.append(source)
        rows, q = extract_rows(document, source, len(sources), units, families, min_span)
        summaries.extend(rows); queries.extend(q)
        if len(queries) > MAX_QUERY_ROWS:
            raise ValueError('query output exceeds32768; split sources/families into separate reports')
    esc = lambda x: '—' if x is None else html.escape(str(x), quote=True)
    display = ['sample_label','family','measurement_status','retained_anchors','representatives_evaluated',
               'full_grammar_compatible_queries','mixed_only_queries','queries_with_ge20_outside_target_grammar']
    labels = ['입력 표시','family','평가 상태','저장 anchor','평가 query','전체 문법 호환','혼합만 호환','문법 밖 연속≥20bp query']
    page = ('<!doctype html><html lang="ko"><meta charset="utf-8"><title>저장 flank 혼합 반복 주석</title>'
        '<style>body{font:16px/1.6 system-ui;margin:2rem}td,th{border:1px solid #ccc;padding:.5rem}table{border-collapse:collapse}</style>'
        '<h1>저장 flank의 혼합 반복 호환성</h1><p>공급한 단위: '+esc(', '.join(units))+
        '. 기존 검출 계수·판정은 유지합니다. 아래 수치는 저장 대표 query 수이며 전체 경계 사건·독립 분자 수가 아닙니다.</p>'
        '<p>표적 문법 밖은 비반복·고유 locus의 증명이 아닙니다. 원 read·mate·quality·assembly 위치는 여기서 확인하지 않습니다.</p>'
        '<table><tr>'+''.join('<th>'+x+'</th>' for x in labels)+'</tr>'+
        ''.join('<tr>'+''.join('<td>'+esc(r[k])+'</td>' for k in display)+'</tr>' for r in summaries)+
        '</table><p>candidate_not_reported: 입력에 후보 없음 · unmeasured/sketch_unavailable: 평가 자료 없음 · '
        'measured_zero_retained: 측정된 저장 항목0 · representatives_incomplete: 일부 대표 서열 없음 · '
        'evaluated_retained_subset: 저장 대표만 평가. 미평가 수치는 대시로 표시합니다.</p>'
        '<p><a href="families.tsv">family별 원 계수·범위</a> · <a href="queries.tsv">query별 상세</a> · '
        '<a href="summary.json">방법·입력 해시·한계</a></p></html>')
    summary = dict(schema_version=1, reporter='report_saved_flank', version=VERSION, algorithm_version=ALGORITHM_VERSION,
        status='complete', completed_utc=datetime.now(timezone.utc).isoformat(), grammar_units=units,
        selected_families=families, minimum_union_span_bp=min_span, full_query_compatibility_uses_cutoff=False,
        family_rows=len(summaries), query_rows=len(queries), evaluated_queries=sum(r['representatives_evaluated'] or 0 for r in summaries),
        full_grammar_compatible_queries=sum(r['full_grammar_compatible_queries'] or 0 for r in summaries),
        mixed_only_queries=sum(r['mixed_only_queries'] or 0 for r in summaries),
        measurement_status_counts=dict(Counter(r['measurement_status'] for r in summaries)), source_files=sources,
        changes_native_counts_or_classification=False, raw_sequence_inputs_read=False, limitations=LIMITS,
        source_hash_scope='Saved JSON bytes; no raw WGS verification', independent_sample_count=None)
    if not summary['evaluated_queries']:
        summary['full_grammar_compatible_queries'] = None
        summary['mixed_only_queries'] = None
    output.mkdir(parents=True, exist_ok=False)
    for name, rows, fields in [('families.tsv',summaries,SUMMARY_FIELDS),('queries.tsv',queries,QUERY_FIELDS)]:
        with (output/name).open('x', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader()
            for row in rows:
                writer.writerow({k: json.dumps(v, sort_keys=True, separators=(',', ':')) if isinstance(v,(dict,list)) else v for k,v in row.items()})
    (output/'index.html').write_text(page, encoding='utf-8')
    summary['outputs'] = [dict(path=str(output/name), bytes=(output/name).stat().st_size,
        sha256=hashlib.sha256((output/name).read_bytes()).hexdigest()) for name in ('families.tsv','queries.tsv','index.html')]
    summary['resources'] = dict(wall_seconds=time.monotonic()-start, process_cpu_seconds=time.process_time()-cpu,
        input_json_bytes=sum(s['size_bytes'] for s in sources), scope='Saved JSON annotation and report only')
    (output/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', nargs='+', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--grammar-unit', action='append', required=True, help='repeatable literal phase-preserved block')
    parser.add_argument('--family', action='append', required=True, help='repeatable candidate family (rotation/RC normalized)')
    parser.add_argument('--min-span', type=int, default=20, help='minimum exact factor span for union annotation; default20')
    args = parser.parse_args(argv)
    try:
        result = build_report(args.results, args.output, args.grammar_unit, args.family, args.min_span)
    except (OSError, ValueError, TypeError) as error:
        parser.exit(2, 'error: '+str(error)+'\n')
    print(json.dumps({k:result[k] for k in ('status','family_rows','query_rows','full_grammar_compatible_queries')}))


if __name__ == '__main__':
    main()
