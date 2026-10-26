#!/usr/bin/env python3
"""Explain fixed actual assembly windows; preserve model gaps and ambiguity.

This optional tool checks the whole source file, analyzes the registered windows
only, and does not discover/classify telomeres or extend detector arrays.
"""
import argparse
from collections import Counter, defaultdict
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
import repeat_structure as structure

PROGRAM = 'trm-repeat-structure-evidence'
EVIDENCE_VERSION = '0.1.1'
SOURCE_ROLES = {
    'unspecified', 'prospective_unseen', 'previously_seen',
    'previously_seen_assembly_reused_explanatory_regression_not_holdout',
}


def source_role(plan):
    role = plan.get('data_role', 'unspecified')
    if not isinstance(role, str) or role not in SOURCE_ROLES:
        raise ValueError('unsupported source data_role')
    return role


def read_json(path):
    return audit.parse_json(Path(path).read_text())


def verify_reference(item):
    path = Path(item['path'])
    before = audit.snapshot(path)
    if (before['size_bytes'] != audit.integer(item['size_bytes'], 'reference size') or
            not audit.SHA.fullmatch(str(item['sha256'])) or audit.digest(path) != item['sha256'] or
            audit.snapshot(path) != before):
        raise ValueError('registered reference hash/size/stability mismatch: '+str(path))
    return path, before


def reference_items(value):
    if isinstance(value, dict):
        if {'path', 'sha256', 'size_bytes'} <= value.keys():
            yield value
        else:
            for item in value.values():
                yield from reference_items(item)
    elif isinstance(value, list):
        for item in value:
            yield from reference_items(item)


def validate_plan(plan, assembly, max_window_bases=10000, max_windows=1000):
    if not isinstance(plan, dict) or type(plan.get('schema_version')) is not int or plan['schema_version'] != 1:
        raise ValueError('supported plan schema1 required')
    source_role(plan)
    expected = plan.get('expected_source')
    if not isinstance(expected, dict) or Path(expected['path']).resolve() != assembly.resolve():
        raise ValueError('assembly and registered source path differ')
    for key in ('records', 'bases', 'size_bytes'):
        audit.integer(expected.get(key), 'source '+key, 1)
    if not audit.SHA.fullmatch(str(expected.get('sha256', ''))):
        raise ValueError('expected source SHA256 required')
    structure.validated_units(plan.get('primary_units'))
    structure.validated_units(plan.get('explanatory_units', []), allow_empty=True)
    windows = plan.get('windows')
    if not isinstance(windows, list) or not 1 <= len(windows) <= max_windows:
        raise ValueError('nonempty window list within explicit maximum required')
    limit = min(max_window_bases, audit.integer(plan.get('max_window_bases', max_window_bases), 'plan window limit', 1))
    ids, intervals = set(), set()
    for w in windows:
        if not isinstance(w, dict) or not isinstance(w.get('id'), str) or not audit.SAFE_ID.fullmatch(w['id']) or w['id'] in ids:
            raise ValueError('unique safe window id required')
        ids.add(w['id'])
        if not isinstance(w.get('record_id'), str) or not w['record_id']:
            raise ValueError('record id required')
        for key in ('record_length', 'start0', 'end0'):
            audit.integer(w.get(key), 'window '+key)
        a,b,L = w['start0'],w['end0'],w['record_length']
        if not 0 <= a < b <= L or b-a > limit or w.get('strand') not in ('+', '-'):
            raise ValueError('invalid window interval, strand or computational size')
        key = (w['record_id'], a, b, w['strand'])
        if key in intervals:
            raise ValueError('duplicate oriented window')
        intervals.add(key)
        terminal = a == 0 or b == L
        if w.get('role') == 'deposited_terminal_window':
            if not terminal or (a == 0 and b == L) or (a == 0 and w['strand'] != '-') or (b == L and w['strand'] != '+'):
                raise ValueError('terminal window must use outward orientation at one original end')
            if w.get('side') not in (None, 'left' if a == 0 else 'right'):
                raise ValueError('window side contradicts original endpoint')
        elif w.get('role') == 'coordinate_internal_unlabelled':
            if terminal:
                raise ValueError('internal window cannot touch original endpoint')
        else:
            raise ValueError('explicit supported window role required')
        for key in ('prior_original_slice_sha256', 'prior_oriented_slice_sha256'):
            if key in w and not audit.SHA.fullmatch(str(w[key])):
                raise ValueError('invalid prior slice SHA256')
    config = plan.get('frozen_model_configuration', {})
    for key, value in [('min_span',100), ('min_exact_anchors',4)]:
        audit.integer(config.get(key,value), key, 1)
    return plan


def write_row(writer, row):
    writer.writerow({k:json.dumps(v, ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in row.items()})


def render_page(path, views, case_id):
    data = json.dumps(views, ensure_ascii=False, allow_nan=False).replace('<', '\\u003c')
    text = '''<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>반복 구조 근거</title><style>body{font:16px system-ui,sans-serif;max-width:1150px;margin:28px auto;padding:0 20px;line-height:1.7;color:#182634}select{padding:9px;max-width:100%}table{border-collapse:collapse;width:100%}th,td{padding:7px;border:1px solid #cbd4dc;text-align:left}svg{width:100%;height:auto;background:#fafbfc}a{color:#135e9e}.muted{color:#576775}</style>
<h1>원본 좌표에서 보는 반복 구조 근거</h1><p>입력 파일 전체의 무결성을 확인하고 등록한 창만 분석했습니다. 아래 구간은 실제 창에 대한 설명이며, 텔로미어 확정·새로운 검출기 성능 검증이 아닙니다.</p>
<p><a href="windows.tsv">창별 요약</a> · <a href="paths.tsv">연속 경로</a> · <a href="fragments.tsv">경계 부분 단위</a> · <a href="alternatives.tsv">다른 단위 길이 해석</a> · <a href="details.jsonl">전체 단위·편집 근거</a> · <a href="summary.json">출처·제한</a></p>
<label for="window">관찰 창: </label><select id="window"></select><p id="context"></p>
<p class="muted">도표의 0…W는 요청한 방향으로 정규화한 창 내부 offset입니다. strand가 −이면 원본 좌표의 증가 방향과 반대입니다. 원본 좌표 범위는 위 설명과 TSV에 별도로 보존합니다.</p>
<svg id="plot" viewBox="0 0 1050 155" role="img" aria-label="반복 단위의 염기별 근거와 최장 연속 경로"></svg>
<p class="muted">파랑: 정확한 완전 단위가 덮음 · 주황: 편집된 완전 단위만 덮음 · 회색: 허용된 완전 단위가 덮지 않음. 초록: 하나의 최장 연속 경로. 위쪽 coverage와 경로 길이는 다르며, 회색 구간을 연결하지 않습니다. 보라 선은 별도 literal 가설이며 기본 경로에 추가하지 않습니다.</p>
<table><thead><tr><th>관측</th><th>값</th></tr></thead><tbody id="values"></tbody></table>
<p>부분 단위는 원래 record 끝과 인공 crop 경계를 구분합니다. 가까운 정확 단위가 지지해도 기존 배열을 연장하지 않습니다. 내부 좌표는 검증된 생물학적 음성으로 취급하지 않습니다. 경로 수는 단위·편집 위치 해석까지 포함해 1 또는 2 이상으로 구분하며 생물학적 확률이 아닙니다.</p>
<script>const data=__DATA__; const select=document.getElementById('window');
const ns='http://www.w3.org/2000/svg'; const plot=document.getElementById('plot');
function element(tag,attrs,text){const node=document.createElementNS(ns,tag);for(const[k,v]of Object.entries(attrs))node.setAttribute(k,v);if(text!==undefined)node.textContent=text;plot.appendChild(node);return node;}
for(let i=0;i<data.length;i++){const o=document.createElement('option');o.value=i;o.textContent=data[i].id;select.appendChild(o);}
function show(){const d=data[+select.value];plot.replaceChildren();document.getElementById('context').textContent=d.record+' ['+d.start+','+d.end+') · strand '+d.strand+' · '+d.role+' · 원본 좌표, 0-based half-open';
const x=p=>125+900*p/d.n;const color={literal:'#277eb1',edited_only:'#e78f32',uncovered:'#d3d9df'};
element('text',{x:0,y:40},'완전 단위 coverage');for(const s of d.segments)element('rect',{x:x(s[0]),y:23,width:Math.max(.5,900*(s[1]-s[0])/d.n),height:24,fill:color[s[2]]});
element('text',{x:0,y:88},'최장 연속 경로');element('line',{x1:125,x2:1025,y1:78,y2:78,stroke:'#d3d9df','stroke-width':2});if(d.longest)element('rect',{x:x(d.longest[0]),y:67,width:900*(d.longest[1]-d.longest[0])/d.n,height:22,fill:'#238b63'});
for(const h of d.hypotheses)element('line',{x1:x(h[0]),x2:x(h[0]),y1:18,y2:51,stroke:'#8b48b4','stroke-width':1});element('text',{x:125,y:125},'0');element('text',{x:990,y:125},d.n+' bp');
const body=document.getElementById('values');body.replaceChildren();for(const[k,v]of Object.entries(d.values)){const tr=document.createElement('tr');for(const value of[k,v]){const td=document.createElement('td');td.textContent=value;tr.appendChild(td);}body.appendChild(tr);}}
select.addEventListener('change',show);show();</script></html>'''
    path.write_text(text.replace('__DATA__',data))


def run(args):
    started = time.monotonic()
    if args.output_dir.exists():
        raise ValueError('output directory exists; previous results are preserved')
    resources = apply_memory_limit(args.memory_limit_mib)
    for key in ('max_window_bases','max_windows','max_edges','max_explanations'):
        audit.integer(getattr(args,key), key, 1)
    plan_before = audit.snapshot(args.plan)
    plan_hash = audit.digest(args.plan)
    source_files = [Path(__file__),Path(structure.__file__),Path(audit.__file__),Path(__file__).with_name('anchor_store.py')]
    code_hashes = {str(p.resolve()):audit.digest(p) for p in source_files}
    plan = validate_plan(read_json(args.plan), args.assembly, args.max_window_bases, args.max_windows)
    refs = {p:s for p,s in (verify_reference(r) for r in reference_items(plan.get('frozen_references', {})))}
    frozen = plan.get('frozen_references', {}).get('frozen_irregular_summary')
    if frozen:
        old = read_json(frozen['path'])
        if (old.get('status') != 'complete' or old.get('complete_input') is not True or
                old.get('validated_eof') is not True or old.get('sampled_prefix') is not False or
                old.get('resource_limit_hit') or old.get('export_truncated') or
                old.get('model_version') != '0.1.0' or old.get('program') != 'scan_irregular_arrays'):
            raise ValueError('complete frozen model source required')
        expected = plan['expected_source']
        if (len(old.get('inputs',[])) != 1 or
                Path(old['inputs'][0]['path']).resolve() != args.assembly.resolve() or
                any(old['inputs'][0][key] != expected[key] for key in ('sha256','size_bytes','records','bases'))):
            raise ValueError('frozen model and planned source provenance disagree')
        if (len(old['definitions']) != 1 or old['definitions'][0]['units'] != sorted(plan['primary_units']) or
                old['config']['min_span'] != plan.get('frozen_model_configuration',{}).get('min_span',100) or
                old['config']['min_exact_anchors'] != plan.get('frozen_model_configuration',{}).get('min_exact_anchors',4)):
            raise ValueError('frozen primary model/config and evidence plan disagree')
    by_record = defaultdict(list)
    for window in plan['windows']:
        by_record[window['record_id']].append(window)
    args.output_dir.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.structure-', dir=args.output_dir.parent) as temporary:
        output = Path(temporary)/'completed'; output.mkdir()
        views, summary_rows = [], []
        counts = Counter(); seen = set()
        source_records = source_bases = 0
        with ExitStack() as stack:
            detail = stack.enter_context((output/'details.jsonl').open('w'))
            writers = {}
            for name, fields in {
                'paths.tsv':['window_id','start0','end0','original_start0','original_end0','span_bp','minimum_edit_count','optimal_path_count_capped2','possible_path_count_capped2','literal_gate_hits','meets_existing_local_span_and_literal_gate','successor_status','next_token_start0','next_base_has_token_coverage','left_crop_context_may_be_missing','right_crop_context_may_be_missing','intervals_may_overlap_do_not_sum_spans'],
                'fragments.tsv':['window_id','start0','end0','original_start0','original_end0','span_bp','oriented_edge','original_edge_side','boundary_kind','vocabulary','unit_index','unit','observed_source_slice','adjacent_literal_anchors','anchored_original_record_fragment','complete_primary_token_explanations_also_possible','extends_detector_array','physical_telomere_confirmed'],
                'alternatives.tsv':['window_id','start0','end0','original_start0','original_end0','span_bp','hypothesis_unit_index','hypothesis_literal','hypothesis_role','primary_same_interval_explanations','local_unit_length_nonidentifiable'],
            }.items():
                handle = stack.enter_context((output/name).open('w',newline=''))
                writers[name] = csv.DictWriter(handle,fields,delimiter='\t',lineterminator='\n');writers[name].writeheader()
            stream = audit.VerifiedStream(stack,args.assembly,plan['expected_source'])
            while stream.stream.peek(1)[:1] in (b'\n',b'\r'):
                stream.stream.readline()
            if stream.stream.peek(1)[:1] != b'>':
                raise ValueError('assembly FASTA required; FASTQ record edges need a different plan')
            for name, sequence in audit.records(stream.stream):
                if name in seen:
                    raise ValueError('duplicate assembly record id')
                seen.add(name); source_records += 1; source_bases += len(sequence)
                expected_length = plan.get('record_lengths',{}).get(name)
                if expected_length is not None and len(sequence) != expected_length:
                    raise ValueError('source record length differs from plan')
                for window in by_record.get(name, []):
                    if len(sequence) != window['record_length']:
                        raise ValueError('window record length mismatch')
                    original = sequence[window['start0']:window['end0']]
                    oriented = original[::-1].translate(structure.COMPLEMENT) if window['strand'] == '-' else original
                    original_hash = hashlib.sha256(original).hexdigest()
                    oriented_hash = hashlib.sha256(oriented).hexdigest()
                    for key,value in [('prior_original_slice_sha256',original_hash),('prior_oriented_slice_sha256',oriented_hash)]:
                        if key in window and window[key] != value:
                            raise ValueError('registered existing window bytes differ: '+window['id'])
                    config = plan.get('frozen_model_configuration',{})
                    result = structure.explain_window(oriented,plan['primary_units'],plan.get('explanatory_units',[]),window,
                        config.get('min_span',100),config.get('min_exact_anchors',4),args.max_edges,args.max_explanations)
                    result['original_slice_sha256'] = original_hash
                    detail.write(json.dumps(result,ensure_ascii=False,allow_nan=False)+'\n')
                    for filename,key in [('paths.tsv','maximal_complete_paths'),('fragments.tsv','boundary_fragments'),('alternatives.tsv','literal_hypotheses')]:
                        for row in result[key]:write_row(writers[filename],{'window_id':window['id'],**row})
                    summary_rows.append({'window_id':window['id'],'record_id':name,'start0':window['start0'],'end0':window['end0'],
                        'strand':window['strand'],'role':window['role'],**result['summary']})
                    for key,value in result['summary'].items():
                        if key not in ('max_complete_path_bp',): counts[key] += value
                    longest = max(result['maximal_complete_paths'],key=lambda x:x['span_bp'],default=None)
                    views.append({'id':window['id'],'record':name,'start':window['start0'],'end':window['end0'],
                        'strand':window['strand'],'role':window['role'],'n':len(oriented),
                        'segments':[[s['start0'],s['end0'],s['class']] for s in result['coverage_segments']],
                        'longest':[longest['start0'],longest['end0']] if longest else None,
                        'hypotheses':[[h['start0'],h['end0']] for h in result['literal_hypotheses']],
                        'values':{'정확 단위 coverage bp':result['summary']['literal_union_bp'],
                          '허용 완전 단위 coverage bp':result['summary']['allowed_token_union_bp'],
                          '최장 연속 경로 bp':result['summary']['max_complete_path_bp'],
                          '최소 편집 해석이 여럿인 경로':result['summary']['optimal_parse_ambiguous_paths'],
                          '허용 해석이 여럿인 경로':result['summary']['any_parse_ambiguous_paths'],
                          '원래 끝의 정확 단위로 지지된 부분 조각':result['summary']['anchored_original_record_fragments'],
                          '인공 crop의 정확 단위로 지지된 부분 조각':result['summary']['anchored_artificial_crop_fragments']}})
            receipt = stream.finish()
        if set(by_record)-seen or len(summary_rows) != len(plan['windows']):
            raise ValueError('planned source windows missing')
        if source_records != plan['expected_source']['records'] or source_bases != plan['expected_source']['bases']:
            raise ValueError('source EOF record/base totals disagree with plan')
        if plan.get('record_lengths') and set(plan['record_lengths']) != seen:
            raise ValueError('assembly record identity set differs from plan')
        with (output/'windows.tsv').open('w',newline='') as handle:
            writer = csv.DictWriter(handle,fieldnames=list(summary_rows[0]),delimiter='\t',lineterminator='\n')
            writer.writeheader();writer.writerows(summary_rows)
        render_page(output/'index.html',views,plan.get('case_id',''))
        if audit.snapshot(args.plan) != plan_before or audit.digest(args.plan) != plan_hash:
            raise ValueError('plan changed during execution')
        for path,before in refs.items():
            if audit.snapshot(path) != before:
                raise ValueError('frozen reference changed during execution')
        if audit.snapshot(args.assembly) != receipt['stat_before_after']:
            raise ValueError('source changed after EOF')
        if any(audit.digest(Path(path)) != value for path,value in code_hashes.items()):
            raise ValueError('program source changed during execution')
        summary = {'schema_version':1,'program':PROGRAM,'evidence_version':EVIDENCE_VERSION,
            'structure_algorithm_version':structure.VERSION,'status':'complete',
            'created_utc':datetime.now(timezone.utc).isoformat(),'plan_sha256':plan_hash,'case_id':plan.get('case_id'),
            'analysis_kind':'registered_windows_not_detector_accuracy',
            'source_data_role':source_role(plan),'source_data_role_is_user_declared':True,
            'source':receipt,'source_records':source_records,'source_bases':source_bases,'window_rows':len(summary_rows),
            'all_source_integrity_checked':True,'whole_source_biologically_analyzed':False,
            'counts_with_nested_and_opposite_strand_window_reuse':dict(counts),
            'primary_units':plan['primary_units'],'explanatory_units':plan.get('explanatory_units',[]),
            'report_only_model_configuration':plan.get('frozen_model_configuration',{}),
            'execution_resources':{**resources,'max_window_bases':args.max_window_bases,'max_windows':args.max_windows,
                'max_edges_per_window':args.max_edges,'max_explanations_per_window':args.max_explanations},
            'changes_detector_arrays':False,'changes_biological_classification':False,
            'physical_telomere_confirmed':False,'source_plan':plan,
            'limitations':['Complete source integrity is different from analysis of every sequence position.',
                'Only complete primary tokens wholly inside each window enter the graph; crop boundaries can censor paths.',
                'Coverage union and overlapping maximal paths are not a tiling or total telomere length.',
                'Capped2 counts distinguish one assignment from at least two; assignment includes unit, operation and edit position.',
                'Explanatory literals and partial fragments do not add primary graph edges or extend detector arrays.',
                'Requested orientation is not joint forward/reverse strand disambiguation.',
                'Original record edges do not establish physical telomeres; internal coordinates are not known negatives.',
                'Repeated positions across nested windows and orientations are not independent observations.'],
            'source_code_sha256':code_hashes,
            'wall_seconds':time.monotonic()-started,
            'outputs':{p.name:{'size_bytes':p.stat().st_size,'sha256':audit.digest(p)} for p in sorted(output.iterdir())}}
        (output/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
        if args.output_dir.exists():raise ValueError('output directory appeared during execution')
        output.rename(args.output_dir)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('assembly','plan','output-dir'):parser.add_argument('--'+key,type=Path,required=True)
    parser.add_argument('--memory-limit-mib',type=int)
    for key,value in [('max-window-bases',10000),('max-windows',1000),('max-edges',200000),('max-explanations',1000000)]:
        parser.add_argument('--'+key,type=int,default=value)
    args = parser.parse_args()
    reserve = bytearray(1024*1024)
    try:
        result = run(args)
    except (ValueError, OSError, MemoryError, KeyError, TypeError) as exc:
        del reserve
        failure = {'program':PROGRAM,'status':'failed','error_type':type(exc).__name__,'error':str(exc),
            'complete':False,'biological_interpretation_available':False}
        try:
            args.output_dir.mkdir(parents=True,exist_ok=False)
        except OSError:pass
        else:
            try:(args.output_dir/'failure.json').write_text(json.dumps(failure)+'\n')
            except OSError:pass
        print(json.dumps(failure),file=sys.stderr);return 1
    print(json.dumps({'status':'complete','window_rows':result['window_rows'],'wall_seconds':result['wall_seconds']}))
    return 0

if __name__=='__main__':sys.exit(main())
