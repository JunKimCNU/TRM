#!/usr/bin/env python3
"""Render sealed saved long-read annotations without fabricating native TRM discovery."""
import argparse
import csv
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import tempfile
import time

import long_read_context_evidence as evidence_module
from long_read_context_report import TABLES,tables,render_html,page,escaped,table_fields
from observation_reporting import REPORTING_VERSION
from report_code_dependencies import local_dependency_paths


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence',type=Path,required=True,help='Sealed long-read evidence manifest')
    parser.add_argument('--output',type=Path,required=True,help='New output directory; existing reports are preserved')
    parser.add_argument('--sample',action='append',help='Optional exact sample IDs; default all manifest samples')
    args=parser.parse_args();start=time.monotonic()
    loaded=evidence_module.load_evidence(args.evidence)
    code_pins={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
               local_dependency_paths(Path(__file__).parent,entries=('report_long_read_context.py',))}
    inputs=[Path(p).resolve() for p in evidence_module.source_dependencies(args.evidence)]
    pins={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    samples=sorted(loaded) if args.sample is None else sorted(set(args.sample))
    if not samples or any(s not in loaded for s in samples):parser.error('Requested sample absent from evidence')
    if any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,119}',s) for s in samples):parser.error('Unsafe sample identifier for output filenames')
    output=args.output.absolute()
    if output.exists() or output.is_symlink():parser.error('Output already exists; choose a new directory to preserve previous results')
    if any(output.resolve()==p or output.resolve() in p.parents for p in inputs):parser.error('Output must not contain evidence input files')
    documents=[]
    for sample in samples:
        evidence=evidence_module.build_evidence(sample,None,loaded)
        documents.append({'schema_version':1,'program':'trm-long-read-context-report',
            'reporting_version':REPORTING_VERSION,'sample_id':sample,
            'scientific_name':evidence.get('identity',{}).get('scientific_name'),
            'report_scope':'saved_selected_long_read_evidence_only',
            'native_discovery_performed':False,'native_discovery_status':'not_provided',
            'classification':None,'detected_motif':None,'data_adequacy':None,
            'policy_assessment':None,'long_read_context_evidence':evidence})
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.long-read-report-',dir=output.parent) as work:
        dest=Path(work)/'report';(dest/'samples').mkdir(parents=True)
        grouped={name:[] for name in TABLES};links=[]
        for doc in documents:
            sid=doc['sample_id'];evidence=doc['long_read_context_evidence']
            (dest/'samples'/f'{sid}.json').write_text(json.dumps(doc,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
            notice='<p>이 페이지는 저장된 표적 read 근거를 보여줍니다. 일반 native 탐색·TRM 정책 판정·미검출 검출력은 미제공이며, 수치 0 또는 분석 실패로 표시하지 않습니다.</p>'
            (dest/'samples'/f'{sid}.html').write_text(page(sid+' — 긴 read 추가 근거',[notice,*render_html(evidence)]))
            links.append('<li><a href="samples/'+escaped(sid)+'.html">'+escaped(sid)+' · '+escaped(doc['scientific_name'])+'</a></li>')
            for name,_,rows in tables(evidence,sid,doc['scientific_name']):grouped[name].extend(rows)
        for name,fields in zip(TABLES,table_fields()):
            with (dest/name).open('w',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',extrasaction='ignore');writer.writeheader();writer.writerows(grouped[name])
        (dest/'index.html').write_text(page('저장된 긴 read 근거 보고',[
            '<p>기술 서열과 같은 read의 반복·reference 구조를 검토한 별도 보고입니다. '
            '새로운 native 탐색 결과나 종의 TRM 확정으로 합산하지 않습니다.</p>','<ul>'+''.join(links)+'</ul>']))
        # Validation happens before publication; do not publish a result if a source changed meanwhile.
        evidence_module.source_dependencies(args.evidence)
        for path,digest in pins.items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:raise ValueError('Evidence changed while rendering: '+path)
        for path,digest in code_pins.items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:raise ValueError('Code changed while rendering: '+path)
        artifacts=[dict(path=str(p.relative_to(dest)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(dest.rglob('*')) if p.is_file()]
        ru=resource.getrusage(resource.RUSAGE_SELF)
        completion={'status':'complete_evidence_only_report','reporting_version':REPORTING_VERSION,'sample_ids':samples,
            'completed_utc':datetime.now(timezone.utc).isoformat(),'classification_provided':False,'native_scan_runs':0,'mapping_runs':0,
            'input_sha256':pins,'code_sha256':code_pins,'artifacts':artifacts,'wall_seconds':time.monotonic()-start,
            'self_CPU_seconds':ru.ru_utime+ru.ru_stime,'max_RSS_KiB':ru.ru_maxrss}
        (dest/'completion.json').write_text(json.dumps(completion,indent=2,ensure_ascii=False)+'\n')
        if output.exists() or output.is_symlink():raise ValueError('Output appeared during report generation')
        os.rename(dest,output)
    print(json.dumps({'status':completion['status'],'samples':len(samples),'output':str(output)}))


if __name__=='__main__':main()
