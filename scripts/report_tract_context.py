#!/usr/bin/env python3
"""Report registered repeat intervals and bounded pair hypotheses on actual reads.

This optional postprocessor streams a hash-bound FASTA/FASTQ once, preserves
original coordinates, and never changes native counts or biological calls.
"""
import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
from datetime import datetime, timezone
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
import resource
import sys
import tempfile
import time

import audit_irregular_arrays as audit
from anchor_store import apply_memory_limit
import tract_context as context

VERSION = '0.1.0'
ROLES = {'development_previously_seen', 'independent_validation', 'unspecified'}
REGION_FIELDS = ['id', 'record_id', 'role', 'start0', 'end0', 'length_bp',
                 'record_length_bp', 'left_record_distance_bp', 'right_record_distance_bp',
                 'expected_strand', 'longest_exact_segment_bp', 'native_tract_union_bp',
                 'best_fixed_phase_identity', 'best_fixed_phase_tie_count']
PAIR_FIELDS = ['id', 'record_id', 'left_tract_id', 'right_tract_id', 'status',
               'reason_codes', 'gap_bp', 'biological_tracts_merged',
               'physical_chromosome_end_confirmed', 'alignment_evidence']


def dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def binding(path):
    path = Path(path).resolve()
    before = audit.snapshot(path)
    sha = audit.digest(path)
    if before != audit.snapshot(path):
        raise ValueError('file changed during hashing: '+str(path))
    return {'path': str(path), 'size_bytes': before['size_bytes'], 'sha256': sha}


def verify_reference(item, base):
    if not isinstance(item, dict):
        raise ValueError('invalid reference binding')
    path = Path(item['path'])
    path = path if path.is_absolute() else base/path
    value = binding(path)
    if (value['size_bytes'] != audit.integer(item.get('size_bytes'), 'reference size')
            or value['sha256'] != item.get('sha256')):
        raise ValueError('reference size/hash mismatch: '+str(path))
    return value


def records_with_quality(stream, max_record_bases):
    """Wrapped records with retained Phred+33; bounded record memory."""
    pending, kind = None, None
    def line():
        raw = stream.readline(max_record_bases+4097)
        if len(raw) > max_record_bases+4096:
            raise ValueError('line exceeds record computation limit')
        return raw
    while True:
        raw = pending if pending is not None else line()
        pending = None
        if not raw:
            return
        header = audit.without_eol(raw)
        if kind is None:
            if not header:
                continue
            kind = header[:1]
            if kind not in (b'>', b'@'):
                raise ValueError('expected FASTA/FASTQ header')
        if not header.startswith(kind):
            raise ValueError('invalid or mixed source record header')
        identifier = audit.record_id(header[1:])
        pieces, size = [], 0
        while True:
            raw = line()
            if not raw:
                if kind == b'@':
                    raise ValueError('truncated FASTQ sequence/separator')
                break
            value = audit.without_eol(raw)
            if kind == b'>' and value.startswith(b'>'):
                pending = raw
                break
            if kind == b'@' and value.startswith(b'+'):
                if value[1:] and audit.record_id(value[1:]) != identifier:
                    raise ValueError('FASTQ separator ID mismatch')
                break
            if kind == b'@' and not value:
                raise ValueError('empty FASTQ sequence line')
            value = audit.normalize_sequence(value, kind == b'>')
            if not value:
                continue
            size += len(value)
            if size > max_record_bases:
                raise ValueError('record exceeds max-record-bases')
            pieces.append(value)
        if not size:
            raise ValueError('empty sequence')
        sequence = b''.join(pieces).decode('ascii')
        quality = None
        if kind == b'@':
            pieces, count = [], 0
            while count < size:
                raw = line()
                if not raw:
                    raise ValueError('truncated FASTQ quality')
                value = audit.without_eol(raw)
                if not audit.QUALITY.fullmatch(value):
                    raise ValueError('invalid Phred+33 quality')
                pieces.append(value)
                count += len(value)
            if count != size:
                raise ValueError('FASTQ sequence/quality length mismatch')
            quality = b''.join(pieces).decode('ascii')
        yield identifier, sequence, quality


def validate_plan(plan, args):
    if not isinstance(plan, dict) or type(plan.get('schema_version')) is not int or plan['schema_version'] != 1:
        raise ValueError('plan schema_version 1 required')
    if plan.get('data_role') not in ROLES:
        raise ValueError('explicit supported data_role required')
    if not isinstance(plan.get('case_id'), str) or not plan['case_id'].strip():
        raise ValueError('case_id required')
    motif = plan.get('motif')
    # Validate the selected family and cutoff contract before reading inputs.
    context.analyze_interval('N', 0, 1, motif, cutoffs=plan.get('cutoffs'))
    context.validate_pair_config(plan.get('pair_config'))
    expected = plan.get('expected_source')
    if not isinstance(expected, dict) or not audit.SHA.fullmatch(str(expected.get('sha256', ''))):
        raise ValueError('expected source hash required')
    for key in ('size_bytes', 'records', 'bases'):
        audit.integer(expected.get(key), 'expected source '+key, 1)
    regions = plan.get('regions')
    if not isinstance(regions, list) or not 1 <= len(regions) <= args.max_regions:
        raise ValueError('nonempty regions within max-regions required')
    ids, record_lengths, tracts = set(), {}, defaultdict(dict)
    total = 0
    for row in regions:
        if not isinstance(row, dict) or not isinstance(row.get('id'), str) or not row['id'] or row['id'] in ids:
            raise ValueError('unique nonempty region ID required')
        ids.add(row['id'])
        record = row.get('record_id')
        if not isinstance(record, str) or not record or any(c.isspace() for c in record):
            raise ValueError('valid record_id required')
        length = audit.integer(row.get('record_length'), 'record length', 1)
        start, end = [audit.integer(row.get(k), k) for k in ('start0', 'end0')]
        if not 0 <= start <= end <= length or end-start > args.max_region_bases or length > args.max_record_bases:
            raise ValueError('region coordinates/computation limit exceeded')
        if record in record_lengths and record_lengths[record] != length:
            raise ValueError('conflicting record lengths')
        record_lengths[record] = length
        if row.get('expected_strand') not in (None, '+', '-'):
            raise ValueError('invalid expected strand')
        total += end-start
        known = row.get('known_tracts', [])
        if not isinstance(known, list):
            raise ValueError('known_tracts must be a list')
        local_ids = set()
        for tract in known:
            if not isinstance(tract, dict):
                raise ValueError('known tract must be an object')
            tid = tract.get('tract_id')
            if not isinstance(tid, str) or not tid or tid in local_ids:
                raise ValueError('unique tract IDs required within region')
            local_ids.add(tid)
            a,b = [audit.integer(tract.get(k), 'tract '+k) for k in ('start0','end0')]
            if not 0 <= a < b <= length or tract.get('strand') not in ('+', '-'):
                raise ValueError('invalid known tract coordinates/strand')
            identity=tract.get('identity')
            if identity is not None and (isinstance(identity,bool) or not isinstance(identity,(int,float))
                    or not math.isfinite(identity) or not 0 <= identity <= 1):
                raise ValueError('known tract identity must be finite in [0,1]')
            if tid in tracts[record] and tracts[record][tid] != tract:
                raise ValueError('conflicting known tract metadata')
            tracts[record][tid] = tract
    if total > args.max_total_region_bases:
        raise ValueError('total region computation limit exceeded')
    pairs = plan.get('pairs', [])
    if not isinstance(pairs, list) or len(pairs) > args.max_regions:
        raise ValueError('pair count exceeds max-regions')
    pair_ids = set()
    for pair in pairs:
        if not isinstance(pair, dict) or not isinstance(pair.get('id'), str) or not pair['id'] or pair['id'] in pair_ids:
            raise ValueError('unique nonempty pair ID required')
        pair_ids.add(pair['id'])
        record = pair.get('record_id')
        if record not in record_lengths:
            raise ValueError('pair requires registered record')
        left,right = [tracts[record].get(pair.get(k)) for k in ('left_tract_id','right_tract_id')]
        if left is None or right is None or pair['left_tract_id'] == pair['right_tract_id']:
            raise ValueError('pair references two distinct known tracts')
        if pair.get('expected_strand') not in (None, '+', '-'):
            raise ValueError('invalid pair expected strand')
    alignments = plan.get('alignments', [])
    if not isinstance(alignments, list) or len(alignments) > 10000:
        raise ValueError('alignment metadata limit exceeded')
    alignment_ids=set()
    for row in alignments:
        if not isinstance(row, dict) or row.get('record_id') not in record_lengths:
            raise ValueError('alignment requires registered record')
        a,b = [audit.integer(row.get(k), k) for k in ('query_start0','query_end0')]
        if not 0 <= a < b <= record_lengths[row['record_id']] or row.get('strand') not in ('+','-'):
            raise ValueError('invalid alignment query coordinates/strand')
        if not isinstance(row.get('target'),str) or not row['target']:
            raise ValueError('alignment target required')
        if audit.integer(row.get('mapq'), 'MAPQ') > 255 or row.get('tp') not in (None,'P','S'):
            raise ValueError('invalid mapping metadata')
        if 'alignment_id' in row:
            aid=row['alignment_id']
            if not isinstance(aid,str) or not aid or aid in alignment_ids:
                raise ValueError('supplied alignment IDs must be unique and nonempty')
            alignment_ids.add(aid)
        if 'query_length' in row and row['query_length'] != record_lengths[row['record_id']]:
            raise ValueError('alignment query length mismatch')
        target_keys={'target_start0','target_end0','target_length'}
        if target_keys.intersection(row):
            if not target_keys <= row.keys():
                raise ValueError('all target coordinates and length required together')
            a,b,length=[audit.integer(row[k],k) for k in ('target_start0','target_end0','target_length')]
            if not 0 <= a < b <= length:
                raise ValueError('invalid alignment target coordinates')
        for key in ('NM','cs_matches','cs_mismatches','cs_insertions','cs_deletions','cs_reference_skips'):
            if key in row: audit.integer(row[key],key)
    return record_lengths, tracts


def tsv_row(writer, row):
    writer.writerow({key: json.dumps(value, ensure_ascii=False, allow_nan=False)
                     if isinstance(value,(dict,list)) else value for key,value in row.items()})


def render(path, plan, rows, pairs, views, runtime):
    esc = lambda x: html.escape(str(x))
    status_labels={'compatible_under_bounded_model':'서열 조건 충족(연결 가설)',
                   'rejected_by_bounded_model':'연결 조건 불충족'}
    reason_labels={'gap_exceeds_max_bp':'사이 구간이 길이 상한 초과',
        'span_exceeds_max_bp':'비교 범위 상한 초과','anchors_overlap':'조각이 겹침',
        'anchors_out_of_order':'조각 좌표 순서 불일치','anchor_orientation_conflict':'조각 방향 불일치',
        'expected_strand_conflict':'예상 방향 불일치','anchor_template_strand_ambiguous':'방향 구분 불가',
        'anchor_phase_ambiguous':'조각의 위상 모호','anchor_declared_strand_not_best':'지정 방향이 최적 설명 아님',
        'anchor_identity_below_threshold':'조각의 단일 위상 일치율 부족',
        'anchor_absolute_phase_mismatch':'두 조각의 위상 불일치',
        'gap_exact_coverage_below_threshold':'사이 구간의 exact coverage 부족',
        'gap_uncovered_run_exceeds_max_bp':'미포함 구간 길이 초과',
        'gap_fixed_phase_identity_below_threshold':'사이 구간의 단일 위상 일치율 부족'}
    mapping_labels={'not_provided':'정렬 자료 없음','no_primary_mapping':'주정렬 없음',
        'split_or_competing_primary_mapping':'분할 또는 경쟁 주정렬 있음',
        'single_primary_covers_pair':'한 주정렬의 query 범위에 포함',
        'single_primary_query_span_only':'한 주정렬의 query 범위만 관측(참조 위치 미제공)',
        'primary_mapping_does_not_cover_pair':'주정렬이 비교 구간 전체를 덮지 않음'}
    headers = ['구간','길이(bp)','최장 exact 조각','Native tract 합집합','최고 단일 위상 일치율','왼쪽/오른쪽 read 끝 거리']
    table = ''.join('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in (
        r['id'],r['length_bp'],r['longest_exact_segment_bp'],r['native_tract_union_bp'],
        '—' if r['best_fixed_phase_identity'] is None else f"{r['best_fixed_phase_identity']:.3f}",
        f"{r['left_record_distance_bp']} / {r['right_record_distance_bp']}"))+'</tr>' for r in rows)
    ptable = ''.join('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in (
        p['id'],status_labels.get(p['status'],p['status']),
        ', '.join(reason_labels.get(r,r) for r in p.get('reason_codes',[])) or '해당 없음',
        mapping_labels.get(p.get('alignment_evidence',{}).get('status'),p.get('alignment_evidence',{}).get('status','미평가'))))
        +'</tr>' for p in pairs)
    data = json.dumps(views,ensure_ascii=False,allow_nan=False).replace('<','\\u003c')
    page = '''<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>반복 조각의 근거</title><style>body{font:16px system-ui,sans-serif;max-width:1250px;margin:30px auto;padding:0 18px;line-height:1.65;color:#20303d}table{border-collapse:collapse;width:100%;font-size:14px}td,th{border:1px solid #ccd6df;padding:7px;text-align:left;overflow-wrap:anywhere}a{color:#126296}select{padding:8px;max-width:100%}svg{width:100%;background:#f7f9fb}details{margin:18px 0}.muted{color:#526474}</style>
<h1>반복 조각의 근거</h1><p>__CASE__ · motif <b>__MOTIF__</b> · 자료 역할: __ROLE__</p>
<p>지정한 원시 서열 구간에서 방향·위상·겹침 제거 coverage·미포함 구간을 비교합니다. 조각 합계는 연속 배열 길이가 아니며, read 끝은 염색체의 물리적 끝을 뜻하지 않습니다.</p>
<p><a href="regions.tsv">구간 요약</a> · <a href="details.jsonl">전체 구간 근거</a> · <a href="pairs.tsv">연결 가설 요약</a> · <a href="pairs.json">연결 가설 상세</a> · <a href="alignments.tsv">기존 전체 read 정렬</a> · <a href="summary.json">설정·출처·시간</a></p>
<h2>구간별 관측</h2><label>구간 선택 <select id="region"></select></label><p id="meta"></p>
<svg id="plot" viewBox="0 0 1050 210" role="img" aria-label="기존 tract와 각 방향의 정확 반복 coverage"></svg>
<p class="muted">회색: 관찰 구간 전체 · 초록: 기존 native tract 합집합 · 파랑/주황: +/− 방향 exact 조각의 합집합. 각 줄의 기준 길이는 선택 구간 안에서 측정합니다. 합집합은 중복 계수를 제거한 면적이며 조각 병합이 아닙니다. 위상 번호는 원본 record 좌표 기준입니다.</p>
<details><summary>Exact 조각의 원본 좌표·방향·위상</summary><table><thead><tr><th>시작</th><th>끝</th><th>길이</th><th>방향</th><th>절대 위상</th></tr></thead><tbody id="segments"></tbody></table></details>
<details><summary>전체 구간 수치</summary><table><thead><tr>__HEADERS__</tr></thead><tbody>__REGIONS__</tbody></table></details>
<h2>제한된 연결 가설</h2><p>두 조각만 비교하는 보수적인 단일 위상 모델입니다. compatible은 설정 범위 안에서의 서열상 양립 가능성을 뜻하며, 텔로미어 확정이나 배열 병합을 수행하지 않습니다. 단계별 통과를 연쇄 연결로 확장하지 않습니다. rejected는 이 제한 모델의 조건 불충족이며, 실제 불규칙 텔로미어의 부재 판정이 아닙니다.</p>
<table><thead><tr><th>쌍</th><th>서열 모델 결과</th><th>조건 불충족 사유</th><th>정렬 맥락</th></tr></thead><tbody>__PAIRS__</tbody></table>
<p>정렬은 입력한 기존 mapping metadata를 재사용합니다. 전체 read 정렬의 mismatch/indel 수와 지정 gap의 서열 진단은 범위가 다릅니다. 분할·복수 정렬은 별도로 표시하고 생물학적 원인을 확정하지 않습니다.</p>
<p>이 실행에서 기본 검출 계수·정책·미검출 판단은 바뀌지 않았습니다. 경과 시간: __TIME__초. 상세 자원 기록은 summary.json을 참조하세요.</p>
<script>const data=__DATA__;const select=document.getElementById('region'),plot=document.getElementById('plot');const ns='http://www.w3.org/2000/svg';
function el(tag,attrs,text){const n=document.createElementNS(ns,tag);for(const[k,v]of Object.entries(attrs))n.setAttribute(k,v);if(text!==undefined)n.textContent=text;plot.appendChild(n);return n;}
for(let i=0;i<data.length;i++){const o=document.createElement('option');o.value=i;o.textContent=data[i].id;select.appendChild(o);}
function show(){const d=data[+select.value];plot.replaceChildren();document.getElementById('meta').textContent=d.record+' ['+d.start+','+d.end+') · 예상 방향 '+(d.expected||'미지정')+' · 0-based half-open';const n=Math.max(1,d.end-d.start),x=p=>190+830*(p-d.start)/n;
const lines=[['기존 native',d.native,'#248868'],['exact ≥'+d.cutoff+' bp (+)',d.plus,'#267aac'],['exact ≥'+d.cutoff+' bp (−)',d.minus,'#d47e26']];lines.forEach((r,i)=>{const y=30+i*57;el('text',{x:0,y:y+16},r[0]);el('rect',{x:190,y,width:830,height:23,fill:'#d8dfe4'});for(const[a,b]of r[1])el('rect',{x:x(a),y,width:Math.max(.5,830*(b-a)/n),height:23,fill:r[2]});});el('text',{x:190,y:200},d.start);el('text',{x:950,y:200},d.end);
const body=document.getElementById('segments');body.replaceChildren();for(const s of d.segments){const tr=document.createElement('tr');for(const v of s){const td=document.createElement('td');td.textContent=v;tr.appendChild(td);}body.appendChild(tr);}}
select.onchange=show;if(data.length)show();</script></html>'''
    for key,value in {'CASE':esc(plan['case_id']),'MOTIF':esc(plan['motif']),'ROLE':esc(plan['data_role']),
                      'HEADERS':''.join('<th>'+h+'</th>' for h in headers),'REGIONS':table,'PAIRS':ptable,
                      'TIME':f'{runtime:.3f}','DATA':data}.items():
        page=page.replace('__'+key+'__',value)
    path.write_text(page)


def run(args):
    started=time.monotonic()
    memory=apply_memory_limit(args.memory_limit_mib)
    code_before=[binding(Path(__file__)),binding(Path(context.__file__)),binding(Path(audit.__file__)),
                 binding(Path(__file__).with_name('anchor_store.py'))]
    source=Path(args.input).resolve(strict=True)
    plan_path=Path(args.plan).resolve(strict=True)
    output=Path(args.output).resolve()
    if output.exists():
        raise ValueError('output already exists; use a new directory')
    for key in ('max_regions','max_record_bases','max_region_bases','max_total_region_bases'):
        audit.integer(getattr(args,key),key,1)
    if plan_path.stat().st_size > 16*1024*1024:
        raise ValueError('plan exceeds 16 MiB')
    plan_before=binding(plan_path)
    plan=audit.parse_json(plan_path.read_text())
    lengths,tracts=validate_plan(plan,args)
    refs=[verify_reference(r,plan_path.parent) for r in plan.get('references',[])]
    by_record,by_pair,by_alignment=defaultdict(list),defaultdict(list),defaultdict(list)
    for row in plan['regions']: by_record[row['record_id']].append(row)
    for row in plan.get('pairs',[]): by_pair[row['record_id']].append(row)
    for row in plan.get('alignments',[]): by_alignment[row['record_id']].append(row)
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.'+output.name+'-',dir=output.parent) as tmp:
        stage=Path(tmp)/'report'; stage.mkdir()
        region_rows,pair_rows,views=[],[],[]
        resolved_cutoffs=None
        seen=set(); records=bases=0
        with ExitStack() as stack:
            verified=audit.VerifiedStream(stack,source,plan['expected_source'])
            details=stack.enter_context((stage/'details.jsonl').open('w'))
            for record,seq,quality in records_with_quality(verified.stream,args.max_record_bases):
                records+=1; bases+=len(seq)
                if record not in lengths: continue
                if record in seen: raise ValueError('duplicate selected record ID')
                seen.add(record)
                if len(seq)!=lengths[record]: raise ValueError('record length mismatch: '+record)
                for spec in by_record[record]:
                    result=context.analyze_interval(seq,spec['start0'],spec['end0'],plan['motif'],
                        quality=quality,known_tracts=spec.get('known_tracts',[]),
                        expected_strand=spec.get('expected_strand'),cutoffs=plan.get('cutoffs'))
                    resolved_cutoffs=sorted(map(int,result['exact_coverage_views']))
                    details.write(json.dumps({'id':spec['id'],'record_id':record,'role':spec.get('role'),
                        'analysis':result},ensure_ascii=False,allow_nan=False)+'\n')
                    best=result['fixed_phase']['all_best_ties']
                    row={key:spec.get(key) for key in ('id','record_id','role','start0','end0','expected_strand')}
                    row.update(length_bp=result['length_bp'],record_length_bp=len(seq),
                        left_record_distance_bp=spec['start0'],right_record_distance_bp=len(seq)-spec['end0'],
                        longest_exact_segment_bp=max((s['length_bp'] for s in result['exact_segments_ge_period']),default=0),
                        native_tract_union_bp=result['known_tract_coverage']['covered_bp'],
                        best_fixed_phase_identity=best[0]['identity'] if best else None,best_fixed_phase_tie_count=len(best))
                    for cutoff,view in result['exact_coverage_views'].items():
                        row['exact'+cutoff+'_union_bp']=view['all_strands']['covered_bp']
                        row['exact'+cutoff+'_uncovered_bp']=view['all_strands']['uncovered_bp']
                    region_rows.append(row)
                    cutoff=str(max(map(int,result['exact_coverage_views'])))
                    view=result['exact_coverage_views'][cutoff]
                    views.append(dict(id=spec['id'],record=record,start=spec['start0'],end=spec['end0'],
                        expected=spec.get('expected_strand'),cutoff=int(cutoff),
                        native=result['known_tract_coverage']['coverage_union_intervals'],
                        plus=view['by_strand']['+']['coverage_union_intervals'],minus=view['by_strand']['-']['coverage_union_intervals'],
                        segments=[[s['start0'],s['end0'],s['length_bp'],s['strand'],s['absolute_template_phase']]
                                  for s in result['exact_segments_ge_period'] if s['length_bp']>=int(cutoff)]))
                for spec in by_pair[record]:
                    result=context.evaluate_pair(seq,tracts[record][spec['left_tract_id']],tracts[record][spec['right_tract_id']],
                        plan['motif'],config=plan.get('pair_config'),alignments=by_alignment.get(record),
                        expected_strand=spec.get('expected_strand'))
                    pair_rows.append(dict(spec,**result))
            input_binding=verified.finish()
        if seen != set(lengths): raise ValueError('selected records missing from source')
        if records!=plan['expected_source']['records'] or bases!=plan['expected_source']['bases']:
            raise ValueError('source record/base counts differ from plan')
        if plan_before!=binding(plan_path): raise ValueError('plan changed during execution')
        for ref in refs: verify_reference(ref,plan_path.parent)
        for ref in code_before: verify_reference(ref,plan_path.parent)
        cutoffs=sorted({k for row in region_rows for k in row if k not in REGION_FIELDS})
        with (stage/'regions.tsv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=REGION_FIELDS+cutoffs,delimiter='\t');writer.writeheader()
            for row in region_rows: tsv_row(writer,row)
        dump(stage/'pairs.json',pair_rows)
        with (stage/'pairs.tsv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=PAIR_FIELDS,delimiter='\t');writer.writeheader()
            for row in pair_rows:
                value={k:row.get(k) for k in PAIR_FIELDS}
                left=tracts[row['record_id']][row['left_tract_id']]
                right=tracts[row['record_id']][row['right_tract_id']]
                value['gap_bp']=right['start0']-left['end0']
                tsv_row(writer,value)
        alignment_rows=plan.get('alignments',[])
        afields=sorted({k for row in alignment_rows for k in row}) or ['record_id']
        with (stage/'alignments.tsv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=afields,delimiter='\t');writer.writeheader()
            for row in alignment_rows: tsv_row(writer,row)
        render(stage/'index.html',plan,region_rows,pair_rows,views,time.monotonic()-started)
        usage=resource.getrusage(resource.RUSAGE_SELF)
        summary=dict(status='complete',program='trm-tract-context',version=VERSION,
            context_version=context.VERSION,schema_version=1,created_utc=datetime.now(timezone.utc).isoformat(),
            case_id=plan['case_id'],data_role=plan['data_role'],motif=plan['motif'],
            pair_config=context.validate_pair_config(plan.get('pair_config')),
            cutoffs=resolved_cutoffs,
            source=input_binding,plan=plan_before,references=refs,
            input_records=records,input_bases=bases,selected_records=len(seen),regions=len(region_rows),pairs=len(pair_rows),
            pair_status_counts=dict(Counter(r['status'] for r in pair_rows)),
            coordinate_system='original record, zero-based half-open',
            changes_native_classification=False,biological_tracts_merged=False,physical_chromosome_end_confirmed=False,
            runtime=dict(wall_seconds=time.monotonic()-started,max_rss_kib=usage.ru_maxrss,
                user_cpu_seconds=usage.ru_utime,system_cpu_seconds=usage.ru_stime,**memory),
            limits={key:getattr(args,key) for key in ('max_regions','max_record_bases','max_region_bases','max_total_region_bases')},
            code=code_before,
            limitations=['Selected original-coordinate intervals only; no full-WGS motif discovery.',
                'Hypothesis thresholds are development settings, not calibrated sensitivity or specificity.',
                'Pairwise compatibility does not create merged arrays or transitive chains.',
                'Read edges are not verified physical chromosome ends.',
                'Whole-read alignment metadata is supplied; no new mapping or sequence-level alignment revalidation.',
                'No biological absence, novelty or non-detection power conclusion.'])
        summary['outputs']=[binding(p) for p in sorted(stage.iterdir()) if p.is_file()]
        for row in summary['outputs']: row['path']=str(output/Path(row['path']).name)
        dump(stage/'summary.json',summary)
        os.rename(stage,output)
    return {'status':'complete','output':str(output),'regions':len(region_rows),'pairs':len(pair_rows),
            'wall_seconds':summary['runtime']['wall_seconds']}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version='trm-tract-context '+VERSION)
    parser.add_argument('--input',required=True,help='Hash-bound FASTA/FASTQ, optionally gzip')
    parser.add_argument('--plan',required=True,help='Schema 1 registered regions, tracts and optional pairs')
    parser.add_argument('--output',required=True,help='New output directory')
    parser.add_argument('--memory-limit-mib',type=int,default=None)
    parser.add_argument('--max-record-bases',type=int,default=1000000)
    parser.add_argument('--max-region-bases',type=int,default=100000)
    parser.add_argument('--max-total-region-bases',type=int,default=2000000)
    parser.add_argument('--max-regions',type=int,default=1000)
    args=parser.parse_args(argv)
    try:
        print(json.dumps(run(args),ensure_ascii=False))
    except (ValueError,OSError,KeyError,TypeError,MemoryError) as error:
        parser.exit(1,'tract-context: '+str(error)+'\n')


if __name__=='__main__':
    main()
