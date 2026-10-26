#!/usr/bin/env python3
"""Compare exact, phase-preserved block grammars in registered assembly windows."""
import argparse
from collections import defaultdict
from contextlib import ExitStack
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import time

import audit_irregular_arrays as audit
from anchor_store import apply_memory_limit
import telomere_grammar as grammar

PROGRAM = 'trm-exact-grammar-evidence'


def validate_plan(plan, assembly, max_window_bases, max_windows):
    if not isinstance(plan, dict) or type(plan.get('schema_version')) is not int or plan['schema_version'] != 1:
        raise ValueError('plan schema 1 required')
    expected = plan['expected_source']
    if not isinstance(expected, dict):
        raise ValueError('source metadata object required')
    if Path(expected['path']).resolve() != assembly.resolve() or not audit.SHA.fullmatch(expected['sha256']):
        raise ValueError('registered source path and SHA256 required')
    for key in ('records', 'bases', 'size_bytes'):
        grammar.natural(expected[key], 'source '+key, 1)
    if plan.get('data_role') not in ('previously_seen', 'prospective_unseen', 'unspecified'):
        raise ValueError('explicit supported source data_role required')
    profiles = plan['profiles']
    if not isinstance(profiles, list) or not 1 <= len(profiles) <= 8:
        raise ValueError('1..8 profiles required')
    ids = set()
    for profile in profiles:
        if not isinstance(profile, dict):
            raise ValueError('profile object required')
        if not isinstance(profile.get('id'), str) or not audit.SAFE_ID.fullmatch(profile['id']) or profile['id'] in ids:
            raise ValueError('unique safe profile ID required')
        ids.add(profile['id']); grammar.validate_units(profile['units'])
    windows = plan['windows']
    if not isinstance(windows, list) or not 1 <= len(windows) <= max_windows:
        raise ValueError('window count outside computation limit')
    ids, coordinates = set(), set()
    for w in windows:
        if not isinstance(w, dict):
            raise ValueError('window object required')
        if not isinstance(w.get('id'), str) or not audit.SAFE_ID.fullmatch(w['id']) or w['id'] in ids:
            raise ValueError('unique safe window ID required')
        ids.add(w['id'])
        if not isinstance(w.get('record_id'), str) or not w['record_id']:
            raise ValueError('record ID required')
        for key in ('start0', 'end0', 'record_length'):
            grammar.natural(w[key], key)
        a, b, length = (w[k] for k in ('start0', 'end0', 'record_length'))
        if not 0 <= a < b <= length or b-a > max_window_bases or w['strand'] not in ('+', '-'):
            raise ValueError('invalid window geometry')
        geometry = (w['record_id'], a, b, w['strand'])
        if geometry in coordinates:
            raise ValueError('duplicate oriented window')
        coordinates.add(geometry)
        if w['role'] == 'deposited_terminal_window':
            if (a == 0) == (b == length) or w['strand'] != ('-' if a == 0 else '+'):
                raise ValueError('terminal window must face outward at one original end')
            if w.get('side') not in (None, 'left' if a == 0 else 'right'):
                raise ValueError('terminal side disagrees with geometry')
        elif w['role'] == 'coordinate_internal_unlabelled':
            if a == 0 or b == length:
                raise ValueError('internal window touches source end')
        else:
            raise ValueError('explicit supported window role required')
    config = plan['configuration']
    if not isinstance(config, dict):
        raise ValueError('configuration object required')
    for key in ('min_span', 'min_blocks'):
        grammar.natural(config[key], key, 1)
    grammar.natural(config['end_distance'], 'end_distance')
    lengths = plan.get('record_lengths', {})
    if not isinstance(lengths, dict):
        raise ValueError('record_lengths object required')
    for name, length in lengths.items():
        if not isinstance(name, str) or not name:
            raise ValueError('nonempty record ID required')
        grammar.natural(length, 'registered record length', 1)
    references = plan.get('references', [])
    if not isinstance(references, list):
        raise ValueError('reference list required')
    for ref in references:
        if not isinstance(ref, dict) or not isinstance(ref.get('path'), str):
            raise ValueError('reference path object required')
        if not isinstance(ref.get('sha256'), str) or not audit.SHA.fullmatch(ref['sha256']):
            raise ValueError('reference SHA256 required')
        grammar.natural(ref.get('size_bytes'), 'reference size')
    return plan


def table(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v for k, v in row.items()})


def render(path, views):
    data = json.dumps(views, ensure_ascii=False, allow_nan=False).replace('<', '\\u003c')
    page = '''<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>불규칙 반복의 정확 문법 비교</title><style>body{font:16px system-ui;max-width:1180px;margin:30px auto;padding:0 20px;line-height:1.65;color:#182634}select{padding:8px;max-width:100%}svg{width:100%;height:auto}table{border-collapse:collapse;width:100%}td,th{border:1px solid #cad4dd;padding:6px;text-align:left}.note{color:#526475}a{color:#135e9e}</style>
<h1>같은 서열 창에서 비교하는 반복 규칙</h1><p>완전한 블록의 정확 일치만 연결합니다. 편집·빈 구간 연결·블록마다의 위상 변경은 허용하지 않습니다. 이 화면은 등록된 창의 설명이며 텔로미어 확정, 신규성 또는 부재 판정이 아닙니다.</p>
<p><a href="windows.tsv">창별 수치</a> · <a href="paths.tsv">연속 경로</a> · <a href="fragments.tsv">부분 블록</a> · <a href="details.jsonl">전체 근거</a> · <a href="summary.json">실행 기록</a></p>
<label>관찰 창 <select id="window"></select></label><p id="context"></p><svg id="plot" viewBox="0 0 1100 280" role="img" aria-label="규칙별 정확 coverage와 최장 연속 경로"></svg>
<p class="note">각 규칙의 위쪽 파랑 막대: 정확 블록이 덮는 염기. 아래쪽 초록 막대: 하나의 최장 연속 경로. 겹치는 블록·경로 길이를 더하지 않습니다. 0…W는 요청 방향의 창 내부 좌표입니다. 부분 블록은 경로를 연장하지 않습니다.</p>
<table><thead><tr><th>규칙</th><th>덮는 bp</th><th>최장 경로 bp</th><th>설정된 원래 끝 거리 내 최장 bp</th><th>지역 조건과 끝 거리 동시 충족</th></tr></thead><tbody id="values"></tbody></table>
<p>사용한 규칙과 거리·길이 조건은 실행 기록의 source_plan에 명시합니다. 서로 다른 규칙의 길이 차이는 민감도 차이가 아닙니다. 넓은 조성 규칙은 일반적인 microsatellite도 설명할 수 있습니다. 창과 양방향은 독립 표본이 아니며 assembly 끝의 물리적 telomere 여부도 미확정입니다.</p>
<script>const views=__DATA__;const selector=document.getElementById('window');const groups=new Map();for(const v of views){if(!groups.has(v.id))groups.set(v.id,[]);groups.get(v.id).push(v);}for(const id of groups.keys()){const o=document.createElement('option');o.value=id;o.textContent=id;selector.appendChild(o);}const svg=document.getElementById('plot');const ns='http://www.w3.org/2000/svg';function item(tag,a,text){const e=document.createElementNS(ns,tag);for(const[k,v]of Object.entries(a))e.setAttribute(k,v);if(text!==undefined)e.textContent=text;svg.appendChild(e);}
function show(){const rows=groups.get(selector.value),d=rows[0];document.getElementById('context').textContent=d.record+' ['+d.start+','+d.end+') strand '+d.strand+' · '+d.role+' · 원본 0-based half-open';svg.replaceChildren();svg.setAttribute('viewBox','0 0 1100 '+(rows.length*75+55));const body=document.getElementById('values');body.replaceChildren();const x=p=>180+880*p/d.n;rows.forEach((r,i)=>{const y=20+i*75;item('text',{x:0,y:y+16},r.profile);item('rect',{x:180,y,width:880,height:18,fill:'#d5dce1'});for(const s of r.covered)item('rect',{x:x(s[0]),y,width:Math.max(.5,880*(s[1]-s[0])/d.n),height:18,fill:'#277eb1'});if(r.longest)item('rect',{x:x(r.longest[0]),y:y+25,width:880*(r.longest[1]-r.longest[0])/d.n,height:16,fill:'#238b63'});const tr=document.createElement('tr');for(const v of[r.profile,r.summary.literal_union_bp,r.summary.max_path_bp,r.summary.near_end_max_path_bp,r.summary.has_local_gate_near_end]){const td=document.createElement('td');td.textContent=v===null?'해당 없음':String(v);tr.appendChild(td);}body.appendChild(tr);});item('text',{x:180,y:rows.length*75+35},'0');item('text',{x:1020,y:rows.length*75+35},d.n+' bp');}selector.addEventListener('change',show);show();</script></html>'''
    path.write_text(page.replace('__DATA__', data))


def run(args):
    started = time.monotonic()
    if args.output_dir.exists():
        raise ValueError('output directory exists; previous results are preserved')
    resources = apply_memory_limit(args.memory_limit_mib)
    for name in ('max_window_bases', 'max_windows', 'max_edges'):
        grammar.natural(getattr(args, name), name, 1)
    before = audit.snapshot(args.plan); plan_hash = audit.digest(args.plan)
    plan = validate_plan(audit.parse_json(args.plan.read_text()), args.assembly, args.max_window_bases, args.max_windows)
    code = [Path(__file__), Path(grammar.__file__), Path(audit.__file__), Path(__file__).with_name('anchor_store.py')]
    code_hashes = {str(p.resolve()): audit.digest(p) for p in code}
    references = {}
    for ref in plan.get('references', []):
        p = Path(ref['path']); snapshot = audit.snapshot(p)
        if snapshot['size_bytes'] != ref['size_bytes'] or audit.digest(p) != ref['sha256']:
            raise ValueError('reference hash/size mismatch')
        references[p] = snapshot
    by_record = defaultdict(list)
    for w in plan['windows']:
        by_record[w['record_id']].append(w)
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.grammar-', dir=args.output_dir.parent) as temp:
        output = Path(temp)/'completed'; output.mkdir()
        rows, paths, fragments, views, seen = [], [], [], [], set()
        records = bases = 0
        with ExitStack() as stack:
            detail = stack.enter_context((output/'details.jsonl').open('w'))
            stream = audit.VerifiedStream(stack, args.assembly, plan['expected_source'])
            while stream.stream.peek(1)[:1] in (b'\n', b'\r'):
                stream.stream.readline()
            if stream.stream.peek(1)[:1] != b'>':
                raise ValueError('assembly FASTA required')
            for name, sequence in audit.records(stream.stream):
                if name in seen:
                    raise ValueError('duplicate source record ID')
                seen.add(name); records += 1; bases += len(sequence)
                if plan.get('record_lengths') and plan['record_lengths'].get(name) != len(sequence):
                    raise ValueError('source record identity/length mismatch')
                for window in by_record.get(name, []):
                    if len(sequence) != window['record_length']:
                        raise ValueError('window record length mismatch')
                    original = sequence[window['start0']:window['end0']]
                    oriented = original[::-1].translate(grammar.COMPLEMENT) if window['strand'] == '-' else original
                    for profile in plan['profiles']:
                        config = plan['configuration']
                        result = grammar.analyze_window(oriented, profile['units'], window,
                            config['min_span'], config['min_blocks'], config['end_distance'], args.max_edges)
                        result['profile_id'] = profile['id']
                        result['original_slice_sha256'] = hashlib.sha256(original).hexdigest()
                        detail.write(json.dumps(result, ensure_ascii=False, allow_nan=False)+'\n')
                        identity = {'window_id': window['id'], 'profile_id': profile['id']}
                        rows.append({**identity, 'record_id': name, 'start0': window['start0'], 'end0': window['end0'],
                            'strand': window['strand'], 'role': window['role'], **result['summary']})
                        paths.extend({**identity, **p} for p in result['paths'])
                        fragments.extend({**identity, **p} for p in result['boundary_fragments'])
                        longest = max(result['paths'], key=lambda p: p['span_bp'], default=None)
                        views.append({'id': window['id'], 'profile': profile['id'], 'record': name,
                            'start': window['start0'], 'end': window['end0'], 'strand': window['strand'],
                            'role': window['role'], 'n': len(oriented), 'summary': result['summary'],
                            'covered': [[s['start0'], s['end0']] for s in result['coverage_segments'] if s['class'] == 'literal'],
                            'longest': [longest['start0'], longest['end0']] if longest else None})
            receipt = stream.finish()
        if records != plan['expected_source']['records'] or bases != plan['expected_source']['bases']:
            raise ValueError('full source record/base totals disagree')
        if set(by_record)-seen or len(rows) != len(plan['windows'])*len(plan['profiles']):
            raise ValueError('missing registered window/profile')
        if plan.get('record_lengths') and set(plan['record_lengths']) != seen:
            raise ValueError('source record ID set differs')
        table(output/'windows.tsv', rows)
        # Empty detail tables still retain headers by deriving the schema from
        # an explicit empty-state template rather than omitting the file.
        for filename, values, fields in [
            ('paths.tsv', paths, ['window_id','profile_id','start0','end0','span_bp','minimum_block_count','maximum_block_count','path_count_capped2','representative_tokens','original_start0','original_end0','original_end_distance_bp','meets_local_gate','near_original_end','span_at_least_500','next_base_covered','termination']),
            ('fragments.tsv', fragments, ['window_id','profile_id','start0','end0','original_start0','original_end0','unit_index','oriented_side','original_side','boundary_kind','adjacent_complete_tokens','observed_source_slice','same_interval_complete_tokens','extends_path'])]:
            if values:
                table(output/filename, values)
            else:
                (output/filename).write_text('\t'.join(fields)+'\n')
        render(output/'index.html', views)
        if audit.snapshot(args.plan) != before or audit.digest(args.plan) != plan_hash:
            raise ValueError('plan changed during execution')
        if any(audit.snapshot(p) != snap for p, snap in references.items()):
            raise ValueError('reference changed during execution')
        if audit.snapshot(args.assembly) != receipt['stat_before_after']:
            raise ValueError('source changed after EOF')
        if any(audit.digest(Path(p)) != digest for p, digest in code_hashes.items()):
            raise ValueError('source code changed during execution')
        summary = {'schema_version': 1, 'program': PROGRAM, 'version': grammar.VERSION, 'status': 'complete',
            'created_utc': datetime.now(timezone.utc).isoformat(), 'source': receipt,
            'source_records': records, 'source_bases': bases, 'source_data_role': plan['data_role'],
            'source_data_role_is_user_declared': True, 'window_rows': len(plan['windows']),
            'profile_window_rows': len(rows), 'profile_count': len(plan['profiles']),
            'plan_sha256': plan_hash, 'source_plan': plan, 'source_code_sha256': code_hashes,
            'all_source_integrity_checked': True, 'whole_source_biologically_analyzed': False,
            'changes_detector_arrays': False, 'changes_biological_classification': False,
            'physical_telomere_confirmed': False, 'wall_seconds': time.monotonic()-started,
            'execution_resources': {**resources, 'max_window_bases': args.max_window_bases,
                'max_windows': args.max_windows, 'max_edges_per_profile_window': args.max_edges},
            'limitations': ['Exact bounded block languages; no edits, gaps, phase normalization, or learned frequencies.',
                'Complete-path length differs from coverage, partial fragments and biological telomere length.',
                'Block counts depend on dictionary; parse counts are capped at two and are not probabilities.',
                'Same-window model comparison is not independently validated recall or specificity.',
                'Nested windows and opposite strands reuse bases. Internal coordinates are unlabelled.',
                '500bp is an annotation; no absence, ALT, novelty or sufficient-power inference.'],
            'outputs': {p.name: {'size_bytes': p.stat().st_size, 'sha256': audit.digest(p)} for p in sorted(output.iterdir())}}
        (output/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
        if args.output_dir.exists():
            raise ValueError('output directory appeared during execution')
        output.rename(args.output_dir)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('assembly', 'plan', 'output-dir'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--memory-limit-mib', type=int)
    for name, default in [('max-window-bases', 10000), ('max-windows', 1000), ('max-edges', 200000)]:
        parser.add_argument('--'+name, type=int, default=default)
    args = parser.parse_args()
    reserve = bytearray(1024*1024)
    try:
        result = run(args)
    except (ValueError, OSError, MemoryError, KeyError, TypeError) as exc:
        del reserve
        failure = {'program': PROGRAM, 'status': 'failed', 'complete': False,
                   'error_type': type(exc).__name__, 'error': str(exc), 'biological_interpretation_available': False}
        try:
            args.output_dir.mkdir(parents=True, exist_ok=False)
        except OSError:
            pass
        else:
            try:
                (args.output_dir/'failure.json').write_text(json.dumps(failure)+'\n')
            except OSError:
                pass
        print(json.dumps(failure), file=sys.stderr); return 1
    print(json.dumps({'status': 'complete', 'profile_window_rows': result['profile_window_rows'], 'wall_seconds': result['wall_seconds']}))
    return 0


if __name__ == '__main__':
    sys.exit(main())
