#!/usr/bin/env python3
"""Compare fixed-phase and bounded-edit hypotheses for registered tract pairs.

This optional evidence report preserves native observations and the original
tract-context model. A sequence path never establishes a biological telomere.
"""
import argparse
from collections import Counter,defaultdict
from contextlib import ExitStack
from datetime import datetime,timezone
import csv
import html
import json
import os
from pathlib import Path
import resource
import sys
import tempfile
import time

import audit_irregular_arrays as audit
from anchor_store import apply_memory_limit
import report_tract_context as io
import tract_context as baseline
import phase_bridge as phase

VERSION='0.1.0'
FIELDS=['id','record_id','motif','left_start0','left_end0','right_start0','right_end0',
        'gap_bp','previous_status','same_range_previous_status','current_status','previous_reasons','current_reasons',
        'edit_budget','minimum_edit_cost','optimal_path_count_capped2','alignment_status',
        'biological_tracts_merged','physical_chromosome_end_confirmed']


def flatten(spec,old,same_range,new,motif):
    bridge=new.get('bridge') or {}
    return dict(id=spec['id'],record_id=spec['record_id'],motif=motif,
        left_start0=old['left']['start0'],left_end0=old['left']['end0'],
        right_start0=old['right']['start0'],right_end0=old['right']['end0'],
        gap_bp=old['observed_gap_bp'],previous_status=old['status'],same_range_previous_status=same_range['status'],current_status=new['status'],
        previous_reasons=old.get('reason_codes',[]),current_reasons=new.get('reason_codes',[]),
        edit_budget=new.get('effective_edit_budget'),minimum_edit_cost=bridge.get('min_edits'),
        optimal_path_count_capped2=bridge.get('optimal_path_count_capped2'),
        alignment_status=old['alignment_evidence']['status'],
        biological_tracts_merged=False,physical_chromosome_end_confirmed=False)


def render(path,plan,rows,details):
    esc=lambda value:html.escape(str(value))
    statuses={'compatible_under_bounded_model':'기존 모델 조건 충족',
        'rejected_by_bounded_model':'기존 모델 조건 불충족',
        phase.SUCCESS:'편집 경로 있음(연결 가설)',phase.REJECTED:'새 모델 조건 불충족',
        phase.COMPUTATION_LIMIT:'계산 상한으로 미평가'}
    reasons={'left_anchor_not_exact':'왼쪽 경계가 정확 반복 아님',
        'right_anchor_not_exact':'오른쪽 경계가 정확 반복 아님',
        'left_anchor_too_short':'왼쪽 경계 근거 길이 부족','right_anchor_too_short':'오른쪽 경계 근거 길이 부족',
        'left_anchor_phase_ambiguous':'왼쪽 경계 위상 모호','right_anchor_phase_ambiguous':'오른쪽 경계 위상 모호',
        'no_path_within_edit_budget_and_indel_run':'편집 수·indel 길이 조건 안의 경로 없음',
        'ambiguous_gap_bases':'gap에 미확정 염기 있음','dp_state_cap_reached':'계산 상태 수 상한 도달',
        'gap_exceeds_max_bp':'gap 길이 상한 초과','span_exceeds_max_bp':'전체 비교 범위 상한 초과',
        'anchors_overlap':'조각이 겹침','anchors_out_of_order':'좌표 순서 불일치',
        'anchor_orientation_conflict':'조각 방향 불일치','expected_strand_conflict':'예상 방향 불일치',
        'anchor_template_strand_ambiguous':'motif 방향 구분 불가'}
    table=''.join('<tr>'+''.join('<td>'+esc(value if value is not None else '미평가')+'</td>' for value in (
        row['id'],row['gap_bp'],statuses.get(row['previous_status'],row['previous_status']),
        statuses.get(row['same_range_previous_status'],row['same_range_previous_status']),
        statuses.get(row['current_status'],row['current_status']),row['minimum_edit_cost'],
        row['edit_budget'],row['optimal_path_count_capped2'],
        ', '.join(reasons.get(r,r) for r in row['current_reasons']) or '해당 없음'))+'</tr>' for row in rows)
    # Whole evidence is available in JSON; the plot only shows the selected path.
    views=[]
    for record in details:
        current=record['current']; bridge=current.get('bridge') or {}
        views.append(dict(id=record['id'],record=record['record_id'],
            gap_start0=record['previous']['left']['end0'],gap_end0=record['previous']['right']['start0'],
            status=current['status'],bridge=bridge,boundary=current.get('anchor_phase_evidence'),
            alignment=record['previous']['alignment_evidence']['status']))
    payload=json.dumps(views,ensure_ascii=False,allow_nan=False).replace('<','\\u003c')
    page='''<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>위상 전환과 gap 연결 가설</title>
<style>body{font:16px system-ui,sans-serif;max-width:1300px;margin:30px auto;padding:0 20px;line-height:1.7;color:#203645}table{border-collapse:collapse;width:100%;font-size:14px}td,th{border:1px solid #c7d3dd;padding:7px;text-align:left;overflow-wrap:anywhere}select{padding:8px;max-width:100%}svg{width:100%;background:#f5f8fa}a{color:#156294}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:13px}</style>
<h1>위상 전환과 gap 연결 가설</h1><p>__CASE__ · motif <b>__MOTIF__</b> · 자료 역할 __ROLE__</p>
<p>양쪽 조각의 정확한 경계 반복으로 위상을 고정하고, gap의 모든 염기를 제한된 편집 경로로 설명하는지 검사합니다. 이전 같은 위상 모델도 같은 입력에서 실행해 나란히 보존했습니다. 추가로 이전 모델의 gap·전체 범위 상한만 새 모델과 맞춘 대조를 제공합니다. 이 대조는 나머지 exact coverage·위상 조건을 유지합니다.</p>
<p><a href="comparison.tsv">이전/현재 비교 TSV</a> · <a href="details.jsonl">양쪽 모델의 전체 근거</a> · <a href="summary.json">설정·출처·자원</a> · <a href="input_plan.json">입력 계획</a></p>
<p><b>compatible은 모델과 양립하는 서열 경로가 있다는 뜻입니다.</b> 텔로미어 위치·생물학적 변이·sequencing error·novelty를 확정하지 않습니다. rejected는 모델 조건 불충족, not_evaluated는 계산 등 미평가 상태입니다. 조각을 병합하거나 연쇄 연결하지 않습니다.</p>
<table><thead><tr><th>쌍</th><th>Gap bp</th><th>이전 모델</th><th>이전 모델: 범위 상한 통일</th><th>새 모델</th><th>최소 편집 수</th><th>허용 편집 수</th><th>최적 경로 수(2=둘 이상)</th><th>조건 불충족·미평가 사유</th></tr></thead><tbody>__ROWS__</tbody></table>
<h2>선택한 쌍의 편집 경로</h2><label>쌍 <select id="pair"></select></label><p id="context"></p>
<svg id="plot" viewBox="0 0 1100 150" role="img" aria-label="gap의 모든 염기를 소비하는 선택 편집 경로"></svg>
<p>파랑 M: 일치 · 주황 X: 치환 · 보라 I: 관측 서열에만 있는 염기 · 빨강 선 D: 반복 template에만 있는 염기. I/D 명칭은 이 정렬 모델 기준이며 실제 변이의 발생 방향을 확정하지 않습니다. 빈 경로나 기각·미평가는 상세 근거를 확인하세요.</p>
<details><summary>선택 경로·경계 위상·제한 상세</summary><pre id="detail"></pre></details>
<p>표시된 경로는 동률 중 결정적으로 고른 한 경로입니다. 최적 경로가 여럿이면 특정 편집 위치를 확정할 수 없습니다. Gap 길이 기반의 허용 편집 수에는 짧거나 빈 gap을 위한 명시적 최소 허용량이 포함됩니다. 이 수치는 개발 가정이며 생물학적 오류 확률이 아닙니다.</p>
<p>기존 전체 read 정렬 맥락은 별도 기록합니다. 분할·복수 정렬이 있어도 서열 계산을 숨기지 않으며 말단 위치를 자동 확인하지 않습니다. 입력 파일 전체의 무결성·record 수·염기 수를 검사합니다.</p>
<script>const data=__DATA__;const select=document.getElementById('pair'),plot=document.getElementById('plot');const ns='http://www.w3.org/2000/svg';
function el(tag,attrs,text){const node=document.createElementNS(ns,tag);for(const[k,v]of Object.entries(attrs))node.setAttribute(k,v);if(text!==undefined)node.textContent=text;plot.appendChild(node);return node;}
for(let i=0;i<data.length;i++){const o=document.createElement('option');o.value=i;o.textContent=data[i].id;select.appendChild(o);}
function show(){const d=data[+select.value];plot.replaceChildren();document.getElementById('context').textContent=d.record+' gap ['+d.gap_start0+','+d.gap_end0+') · '+d.status+' · '+d.alignment;document.getElementById('detail').textContent=JSON.stringify(d,null,2);
const n=Math.max(1,d.gap_end0-d.gap_start0),x=p=>80+940*p/n;el('rect',{x:80,y:45,width:940,height:35,fill:'#dbe2e7'});const colors={M:'#267fb0',X:'#d48223',I:'#8c59ad',D:'#c33942'};
for(const op of (d.bridge.path_ops||[])){const a=op.query_start0-d.gap_start0,b=op.query_end0-d.gap_start0;if(!Number.isFinite(a)||!Number.isFinite(b))continue;if(op.operation==='D')el('line',{x1:x(a),x2:x(a),y1:35,y2:90,stroke:colors.D,'stroke-width':2});else el('rect',{x:x(a),y:45,width:Math.max(.5,940*(b-a)/n),height:35,fill:colors[op.operation]||'#777'});}el('text',{x:80,y:120},'gap offset 0');el('text',{x:940,y:120},Math.max(0,d.gap_end0-d.gap_start0));}
select.onchange=show;if(data.length)show();</script></html>'''
    for key,value in {'CASE':esc(plan['case_id']),'MOTIF':esc(plan['motif']),'ROLE':esc(plan['data_role']),
                      'ROWS':table,'DATA':payload}.items(): page=page.replace('__'+key+'__',value)
    path.write_text(page)


def run(args):
    started=time.monotonic(); memory=apply_memory_limit(args.memory_limit_mib)
    for key in ('max_regions','max_record_bases','max_region_bases','max_total_region_bases'):
        audit.integer(getattr(args,key),key,1)
    source=Path(args.input).resolve(strict=True);plan_path=Path(args.plan).resolve(strict=True)
    output=Path(args.output).resolve()
    if output.exists(): raise ValueError('output already exists; choose a new directory')
    if plan_path.stat().st_size>16*1024*1024: raise ValueError('plan exceeds 16 MiB')
    code=[io.binding(Path(m.__file__)) for m in (phase,baseline,io,audit)]
    code.extend([io.binding(Path(__file__)),io.binding(Path(__file__).with_name('anchor_store.py'))])
    before=io.binding(plan_path);plan=audit.parse_json(plan_path.read_text())
    lengths,tracts=io.validate_plan(plan,args)
    if not plan.get('pairs'): raise ValueError('nonempty registered pairs required')
    config=phase.validate_config(plan.get('phase_config'))
    same_range_config=dict(baseline.validate_pair_config(plan.get('pair_config')),
        max_gap_bp=config['max_gap_bp'],max_span_bp=config['max_span_bp'])
    refs=[io.verify_reference(r,plan_path.parent) for r in plan.get('references',[])]
    pair_groups,alignments=defaultdict(list),defaultdict(list)
    for item in plan['pairs']:pair_groups[item['record_id']].append(item)
    for item in plan.get('alignments',[]):alignments[item['record_id']].append(item)
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.'+output.name+'-',dir=output.parent) as tmp:
        stage=Path(tmp)/'report';stage.mkdir()
        seen=set();record_count=bases=0;rows=[];details=[]
        with ExitStack() as stack:
            verified=audit.VerifiedStream(stack,source,plan['expected_source'])
            detail_file=stack.enter_context((stage/'details.jsonl').open('w'))
            for record,seq,quality in io.records_with_quality(verified.stream,args.max_record_bases):
                record_count+=1;bases+=len(seq)
                if record not in lengths:continue
                if record in seen:raise ValueError('duplicate selected record ID: '+record)
                seen.add(record)
                if len(seq)!=lengths[record]:raise ValueError('registered record length mismatch')
                for spec in pair_groups[record]:
                    left,right=[tracts[record][spec[k]] for k in ('left_tract_id','right_tract_id')]
                    options=dict(alignments=alignments.get(record),expected_strand=spec.get('expected_strand'))
                    old=baseline.evaluate_pair(seq,left,right,plan['motif'],config=plan.get('pair_config'),**options)
                    same_range=baseline.evaluate_pair(seq,left,right,plan['motif'],config=same_range_config,**options)
                    new=phase.evaluate_pair(seq,left,right,plan['motif'],config=config,**options)
                    item=dict(id=spec['id'],record_id=record,motif=plan['motif'],previous=old,same_range_previous=same_range,current=new)
                    detail_file.write(json.dumps(item,ensure_ascii=False,allow_nan=False)+'\n')
                    details.append(item);rows.append(flatten(spec,old,same_range,new,plan['motif']))
            input_binding=verified.finish()
        if seen!=set(lengths):raise ValueError('registered record missing')
        if record_count!=plan['expected_source']['records'] or bases!=plan['expected_source']['bases']:
            raise ValueError('whole source record/base count mismatch')
        if before!=io.binding(plan_path):raise ValueError('plan changed during execution')
        for item in refs+code:io.verify_reference(item,plan_path.parent)
        with (stage/'comparison.tsv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=FIELDS,delimiter='\t');writer.writeheader()
            for row in rows:io.tsv_row(writer,row)
        io.dump(stage/'input_plan.json',plan)
        render(stage/'index.html',plan,rows,details)
        usage=resource.getrusage(resource.RUSAGE_SELF)
        result=dict(status='complete',program='trm-phase-bridge-report',version=VERSION,
            phase_model_version=phase.VERSION,previous_model_version=baseline.VERSION,
            created_utc=datetime.now(timezone.utc).isoformat(),case_id=plan['case_id'],data_role=plan['data_role'],
            motif=plan['motif'],phase_config=config,previous_config=baseline.validate_pair_config(plan.get('pair_config')),
            same_range_previous_config=same_range_config,
            source=input_binding,plan=before,references=refs,code=code,
            input_records=record_count,input_bases=bases,selected_records=len(seen),pairs=len(rows),
            previous_status_counts=dict(Counter(r['previous_status'] for r in rows)),
            same_range_previous_status_counts=dict(Counter(r['same_range_previous_status'] for r in rows)),
            current_status_counts=dict(Counter(r['current_status'] for r in rows)),
            status_transitions=dict(Counter(r['previous_status']+' -> '+r['current_status'] for r in rows)),
            runtime=dict(wall_seconds=time.monotonic()-started,user_cpu_seconds=usage.ru_utime,
                system_cpu_seconds=usage.ru_stime,max_rss_kib=usage.ru_maxrss,**memory),
            changes_native_classification=False,biological_tracts_merged=False,physical_chromosome_end_confirmed=False,
            limitations=['Selected registered pairs only; no full-WGS motif discovery.',
                'Bounded edit paths are descriptive hypotheses, not errors, variants or physical terminal proof.',
                'Same settings are used across cases; outcome rates are not biological sensitivity/specificity.',
                'Both models are run on the same actual input; native observations and earlier files remain unchanged.'])
        result['outputs']=[io.binding(p) for p in sorted(stage.iterdir())]
        for item in result['outputs']:item['path']=str(output/Path(item['path']).name)
        io.dump(stage/'summary.json',result)
        os.rename(stage,output)
    return dict(status='complete',output=str(output),pairs=len(rows),current_status_counts=result['current_status_counts'],
                wall_seconds=result['runtime']['wall_seconds'])


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='trm-phase-bridge-report '+VERSION)
    for key in ('input','plan','output'):parser.add_argument('--'+key,required=True)
    parser.add_argument('--memory-limit-mib',type=int)
    parser.add_argument('--max-record-bases',type=int,default=1000000)
    parser.add_argument('--max-region-bases',type=int,default=100000)
    parser.add_argument('--max-total-region-bases',type=int,default=2000000)
    parser.add_argument('--max-regions',type=int,default=1000)
    args=parser.parse_args(argv)
    try:print(json.dumps(run(args),ensure_ascii=False))
    except (ValueError,OSError,KeyError,TypeError,MemoryError) as error:parser.exit(1,'phase-bridge: '+str(error)+'\n')


if __name__=='__main__':main()
