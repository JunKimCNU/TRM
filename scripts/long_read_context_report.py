"""Shared HTML/TSV presentation for optional long-read evidence and standalone reports."""
import html
import json
from urllib.parse import urlsplit
import long_read_context_evidence as evidence_module

TABLES=('long_read_context_summary.tsv','long_read_context_reads.tsv','long_read_context_contexts.tsv')
LABELS={
    'selected_reads':'선정 read 수','contexts':'선정 reference 문맥 수',
    'full_signature_anywhere':'전체 기술 서열 일치 · tail 안 어느 위치든',
    'full_signature_at_read_end':'전체 기술 서열 일치 · 정확히 read 끝',
    'whole_tail_signature':'꼬리 전체가 기술 서열과 일치',
    'partial_signature_ge20_at_read_end':'20 bp 이상 부분 일치 · 정확히 read 끝',
    'boundary_overlap_reads':'기술 서열이 native 반복 경계와 겹치는 read',
    'structure_reviewed_reads':'지역 구조를 검토한 read',
    'other_tract_support_reads':'추가 native tract가 바깥 반복 일부에 정렬된 read',
    'distal_boundary_support_reads':'바깥 반복 시작 이후까지 정렬',
    'record_end_support_reads':'조립 record 마지막 염기 정렬',
    'exact_match_reads_ge12':'기술 서열 12 bp 이상 정확 일치',
    'exact_match_reads_ge16':'기술 서열 16 bp 이상 정확 일치',
    'exact_match_reads_ge20':'기술 서열 20 bp 이상 정확 일치',
}


def escaped(value):
    return html.escape('미평가' if value is None else str(value),quote=True)


def table_fields():
    def unique(values):return tuple(dict.fromkeys(values))
    return (
        unique(('sample_id','scientific_name','source_run','status','scope',*evidence_module.SUMMARY_FIELDS)),
        unique(('sample_id','scientific_name','source_run',*evidence_module.READ_FIELDS)),
        unique(('sample_id','scientific_name','source_run',*evidence_module.CONTEXT_FIELDS)),
    )


def tables(evidence,sample,scientific_name):
    """Return rows without inventing zero for unprovided optional measurements."""
    groups=([dict(status=evidence['status'],scope=evidence.get('scope'),**evidence['summary'])],
            evidence.get('reads',[]),evidence.get('contexts',[]))
    for name,fields,rows in zip(TABLES,table_fields(),groups):
        result=[]
        for row in rows:
            values=dict(row,sample_id=sample,scientific_name=scientific_name,source_run=evidence.get('source_run'))
            result.append({k:json.dumps(v,ensure_ascii=False,sort_keys=True,allow_nan=False)
                           if isinstance(v,(dict,list)) else v for k,v in values.items()})
        yield name,fields,result


def render_html(evidence,table_prefix='../'):
    if not evidence:return []
    out=['<h2>긴 read의 기술 서열과 바깥 반복 구조</h2>']
    if evidence.get('status') in ('unassessed','not_provided'):
        out.append('<p>선택적 긴 read 근거 미제공. 추가 관측은 미평가이며 0건이나 TRM 부재를 뜻하지 않습니다.</p>')
    else:
        out.append('<p>출처 run: <strong>'+escaped(evidence.get('source_run'))+'</strong>. '
            '미리 선정한 read의 추가 주석입니다. 원 검출 계수·판정·검출력은 그대로 유지합니다.</p>')
        signature=evidence.get('signature') or {}
        sequence=signature.get('sequence')
        if sequence:
            out.append('<p>대조한 알려진 기술 서열: <code>'+escaped(sequence)+'</code> · '+
                escaped(signature.get('length_bp',len(sequence)))+' bp.</p>')
        url=signature.get('source_url')
        if isinstance(url,str) and urlsplit(url).scheme in ('https','http'):
            out.append('<p><a href="'+escaped(url)+'">기술 서열의 원 출처</a></p>')
        out.append('<table><tr><th>관측</th><th>선정 범위의 수치</th></tr>')
        for field in evidence_module.SUMMARY_FIELDS:
            out.append('<tr><th>'+escaped(LABELS.get(field,field))+'</th><td>'+escaped(evidence['summary'][field])+'</td></tr>')
        out.append('</table><p>12/16/20 bp 기준은 서로 포함됩니다. 전체 서열 일치, read 끝의 일치, 꼬리 전체 일치를 구분합니다. '
            '검토하지 않은 위치는 미평가로 남깁니다.</p>')
        columns=(('context_id','문맥'),('target_record','Reference'),('linked_read_count','기존 연결 read'),
            ('selected_boundary_same_hsp_count','선정 경계 same-HSP'),
            ('selected_boundary_within100bp_count','선정 경계 100 bp 이내'),
            ('structure_reviewed_read_count','지역 검토 read'),
            ('distal_boundary_supported_read_count','바깥 반복 도달'),
            ('record_end_supported_read_count','조립 끝 염기 정렬'),
            ('other_native_tract_supported_read_count','추가 tract의 바깥 반복 일부 정렬'))
        out.append('<h3>대표 반복의 경계와 같은 read의 추가 구조</h3><table><tr>'+''.join('<th>'+escaped(label)+'</th>' for _,label in columns)+'</tr>')
        for row in evidence.get('contexts',[]):
            out.append('<tr>'+''.join('<td>'+escaped(row.get(key))+'</td>' for key,_ in columns)+'</tr>')
        out.append('</table><p>바깥 반복 도달은 선정 tract의 경계 지표를 대체하지 않습니다. '
            '추가 tract 지지는 바깥 배열 일부와의 정렬이며 양쪽 경계를 가로질렀다는 뜻이 아닙니다. '
            '조립 record 끝과 adapter 인접성만으로 물리적 텔로미어를 확정하지 않습니다.</p>')
    out.append('<p>'+' · '.join('<a href="'+escaped(table_prefix+name)+'">'+label+'</a>' for name,label in zip(TABLES,('요약 TSV','read별 좌표·기술 서열 TSV','문맥별 구조 TSV')))+'</p>')
    if evidence.get('limitations'):
        out.append('<details><summary>해석 범위와 남은 한계</summary><ul>'+''.join('<li>'+escaped(x)+'</li>' for x in evidence['limitations'])+'</ul></details>')
    if evidence.get('notes'):
        out.append('<details><summary>이 자료의 검토 주석</summary><ul>'+''.join('<li>'+escaped(x)+'</li>' for x in evidence['notes'])+'</ul></details>')
    return out


def page(title,body):
    return ('<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
            '<title>'+escaped(title)+'</title><style>body{font:16px/1.65 sans-serif;max-width:1200px;margin:40px auto;padding:0 22px;color:#20323c}'
            'table{border-collapse:collapse;width:100%;margin:18px 0}td,th{padding:8px;border-bottom:1px solid #dce3e7;text-align:left}'
            'th{background:#eef3f5}code{overflow-wrap:anywhere}a{color:#176c80}details{margin:24px 0}</style><body><h1>'+escaped(title)+'</h1>'+''.join(body)+'</body></html>')
