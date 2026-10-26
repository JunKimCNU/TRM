"""Attach reviewed, sample-bound followup diagnostics without changing decisions.

Optional SAMPLE/{status,evidence}.json sidecars refer to small saved artifacts.
The original sequencing data and compact sequencing subsets are never reread.
Counts describe registered pairs/windows, not whole-library detection power.
"""
from __future__ import annotations

from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path
import re

SCHEMA_VERSION = 1
REPORTING_VERSION = '1.0.0'
MAX_ARTIFACT_BYTES = 64 * 1024 * 1024
SHA = re.compile(r'[0-9a-f]{64}')
ROUTES = {
    'geometry_or_orientation_blocked': '겹침·방향 등 선행 조건 미충족',
    'exact_boundary_blocked': '고정 경계의 정확 반복 미충족',
    'sequence_alphabet_blocked': 'gap의 미확정 염기로 미평가',
    'computation_limit_unevaluated': '계산 상한으로 미평가',
    'no_path_phase_length_budget_incompatible': '위상·길이·편집 예산의 필요조건 미충족',
    'no_path_despite_feasible_template_length': '가능한 길이는 있으나 현재 조건에서 경로 없음',
    'compatible_empty_gap': '빈 gap의 연결 가설 성립',
    'compatible_positive_gap': '양의 길이 gap의 연결 가설 성립',
}
SCOPE = '등록된 부분 자료의 고정 pair/window 후속 진단; 전체 샘플의 검출률·부재·충분성 평가가 아님'


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _count(value, name):
    _require(type(value) is int and value >= 0, 'invalid nonnegative count: ' + name)
    return value


def _read(path):
    path = Path(path).resolve()
    _require(path.is_file(), 'diagnostic artifact does not exist: ' + str(path))
    before = path.stat()
    _require(path.is_file() and before.st_size <= MAX_ARTIFACT_BYTES, 'missing or oversized diagnostic artifact: ' + str(path))
    content = path.read_bytes()
    after = path.stat()
    _require((before.st_size, before.st_mtime_ns, before.st_ino) ==
             (after.st_size, after.st_mtime_ns, after.st_ino), 'diagnostic artifact changed while reading')
    return content, {'path': str(path), 'size_bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest()}


def _normalize_binding(ref):
    _require(isinstance(ref, dict), 'artifact binding must be an object')
    path = Path(ref.get('path', ''))
    _require(path.is_absolute(), 'diagnostic artifact paths must be absolute')
    size = _count(ref.get('size_bytes', ref.get('bytes')), 'artifact bytes')
    digest = ref.get('sha256')
    _require(isinstance(digest, str) and SHA.fullmatch(digest) is not None, 'artifact SHA256 is invalid')
    return {'path': str(path.resolve()), 'size_bytes': size, 'sha256': digest}


def _artifact(ref, table=False):
    expected = _normalize_binding(ref)
    content, actual = _read(expected['path'])
    _require(actual == expected, 'diagnostic artifact hash or size mismatch: ' + expected['path'])
    text = content.decode('utf-8')
    value = list(csv.DictReader(io.StringIO(text), delimiter='\t')) if table else json.loads(text)
    return value, actual


def _contains_binding(value, wanted):
    if isinstance(value, dict):
        if {'path', 'sha256'} <= value.keys() and ('bytes' in value or 'size_bytes' in value):
            if _normalize_binding(value) == wanted:
                return True
        return any(_contains_binding(x, wanted) for x in value.values())
    if isinstance(value, list):
        return any(_contains_binding(x, wanted) for x in value)
    return False


def _review(ref, artifacts):
    doc, binding = _artifact(ref)
    _require(doc.get('status') == 'passed', 'diagnostic evidence requires a passed review')
    for item in artifacts:
        _require(_contains_binding(doc, item), 'review does not bind diagnostic artifact: ' + item['path'])
    return binding


def _identity(status, sample, meta, raw, source_hash):
    _require(status.get('schema_version') == 1 and status.get('sample_id') == sample,
             'followup diagnostic schema/sample mismatch')
    identity = status.get('identity', {})
    _require(isinstance(source_hash, str) and SHA.fullmatch(source_hash) is not None and
             identity.get('native_result_sha256') == source_hash, 'followup native result hash differs')
    for key in ('scientific_name', 'tax_id'):
        expected = str(meta.get(key, '')).strip()
        _require(bool(expected) and str(identity.get(key, '')).strip() == expected,
                 'followup taxon identity differs: ' + key)
    _require(bool(raw.get('version')) and identity.get('native_version') == raw['version'],
             'followup native version differs')
    _require(bool(raw.get('inputs')) and identity.get('inputs') == raw['inputs'],
             'followup native input identity differs')
    for key in ('records', 'bases'):
        value = _count(raw.get('summary', {}).get(key), 'native ' + key)
        _require(_count(identity.get(key), 'identity ' + key) == value,
                 'followup native totals differ: ' + key)
    return identity


def _relation(spec, summary, raw, kind):
    receipt, receipt_binding = _artifact(spec['receipt'])
    historical = _normalize_binding(spec['original_input_binding'])
    review_binding = _review(spec['review'], [receipt_binding,historical])
    inputs = raw['inputs']
    _require(len(inputs) == 1 and isinstance(inputs[0], dict),
             'this derivation receipt requires exactly one native source input')
    original = inputs[0]
    original_path = str(Path(original['path']).resolve())
    _count(original.get('bytes'), 'native input bytes')
    _require(historical['path'] == original_path and historical['size_bytes'] == original['bytes'],
             'historical source binding differs from native input')
    if spec.get('kind') == 'registered_full_record_extraction':
        _require(kind == 'fixed_pair_bridge', 'extraction relation does not match diagnostic kind')
        parent = receipt.get('source_identity', {})
        _require(str(Path(parent.get('path', '')).resolve()) == original_path and
                 parent.get('bytes', parent.get('size_bytes')) == original['bytes'],
                 'derived diagnostic parent input differs from native input')
        _require(receipt.get('complete_input_eof') is True, 'derived diagnostic parent EOF was not verified')
        for key in ('records', 'bases'):
            _require(_count(receipt.get('input_' + key), 'parent ' + key) == raw['summary'][key],
                     'derivation parent totals differ from native input')
        derived = receipt['output_binding']
        measured = summary['input_binding']
        _require(_normalize_binding(derived) == _normalize_binding(measured),
                 'derivation output differs from diagnostic input')
        for key in ('records', 'bases'):
            _require(_count(receipt.get(key), 'derived ' + key) == derived.get(key) == summary.get('source_' + key),
                     'derived diagnostic totals mismatch')
        caveat = receipt.get('source_hash_scope', 'Archived extraction and analysis receipts; no source checksum recomputed during report join.')
        selected_records = summary['selected_records']
    elif spec.get('kind') == 'sequence_unique_full_record_selection':
        _require(kind == 'registered_mixed_windows', 'sequence selection relation does not match diagnostic kind')
        _require(str(Path(receipt.get('source_fastq', '')).resolve()) == original_path,
                 'selected diagnostic parent input differs from native input')
        _require(receipt.get('gzip_eof_verified') is True, 'selected diagnostic parent EOF was not verified')
        for key in ('records', 'bases'):
            _require(_count(receipt.get('total_full_' + key), 'parent ' + key) == raw['summary'][key],
                     'selection parent totals differ from native input')
        measured = summary['source']
        _require(str(Path(receipt.get('extracted_fastq', '')).resolve()) == measured.get('path') and
                 receipt.get('extracted_fastq_sha256') == measured.get('sha256'),
                 'selected full-record input differs from diagnostic input')
        _require(receipt.get('distinct_sequence_unique_extracted_archive_reads') == summary.get('source_records'),
                 'selected record denominator differs')
        _require(SHA.fullmatch(str(receipt.get('source_fastq_sha256', ''))) is not None,
                 'selection receipt lacks historical source checksum')
        _require(receipt['source_fastq_sha256'] == historical['sha256'],
                 'selection historical source checksum differs')
        caveat = receipt.get('caveat', '')
        selected_records = summary['selected_records']
    else:
        raise ValueError('unsupported followup source relation')
    _require(measured.get('validated_eof') is True, 'diagnostic subset input EOF was not verified')
    records = _count(summary.get('source_records'), 'diagnostic source records')
    bases = _count(summary.get('source_bases'), 'diagnostic source bases')
    selected_records = _count(selected_records, 'diagnostic selected records')
    _require(selected_records <= records <= raw['summary']['records'] and bases <= raw['summary']['bases'],
             'diagnostic subset exceeds original input totals')
    return dict(kind=spec['kind'], original_input=original, original_records=raw['summary']['records'],
        original_bases=raw['summary']['bases'], compact_input=_normalize_binding(measured),
        compact_records=records, compact_bases=bases, diagnostic_selected_records=selected_records,
        receipt_binding=receipt_binding, review_binding=review_binding, historical_original_input_binding=historical,
        source_hash_scope='archived_provenance_chain_no_raw_or_compact_sequence_reread',
        caveat=caveat, original_molecule_identity_confirmed=receipt.get('original_name_identity_confirmed'))


def _row(metric, label, value, unit, scope):
    return dict(metric=metric, label=label, value=value, unit=unit, scope=scope)


def _phase(summary, table, source):
    _require(summary.get('source_id') == source, 'diagnostic source ID differs')
    _require(summary.get('new_model_evaluations') is False and summary.get('changes_native_classification') is False,
             'followup summary changes native decisions or reevaluates models')
    total = _count(summary.get('pairs'), 'registered pair count')
    _require(len(table) == total and len({r['pair_id'] for r in table}) == total,
             'diagnostic pair table count or unique identity differs')
    counts = Counter()
    for r in table:
        _require(r['source_id'] == source and r['diagnostic_stage'] in ROUTES,
                 'diagnostic pair table source or route differs')
        _require(r['native_status_changed'] == 'False' and r['model_recomputed'] == 'False',
                 'diagnostic pair claims a changed native decision')
        counts[r['diagnostic_stage']] += 1
    saved = summary.get('frozen_stage_counts', {})
    _require(all(k in ROUTES for k in saved) and all(_count(v, k) == counts.get(k, 0) for k, v in saved.items())
             and sum(saved.values()) == total and dict(counts) == {k:v for k,v in saved.items() if v},
             'diagnostic stage counts differ from the bound pair table')
    measures = dict(kind='fixed_pair_bridge', registered_pairs=total,
        unique_tracts=_count(summary.get('unique_tracts'), 'unique tracts'),
        unique_intervals=_count(summary.get('unique_intervals'), 'unique intervals'),
        core_unavailable=_count(summary.get('empty_cores'), 'empty cores'), stage_counts=dict(counts))
    rows = [_row('registered_pairs', '진단한 고정 pair', total, 'pairs', 'registered_pairs')]
    rows.extend(_row(k, label, counts.get(k, 0), 'pairs', 'registered_pairs') for k,label in ROUTES.items())
    rows.extend([_row('unique_tracts', '대상 고유 tract', measures['unique_tracts'], 'tracts', 'referenced_native_tracts'),
                 _row('core_unavailable', '양끝 2k 제외 후 빈 core', measures['core_unavailable'], 'tracts', 'referenced_native_tracts')])
    recommendations = []
    if counts['exact_boundary_blocked']:
        recommendations.append('고정 경계 창의 염기 차이와 품질을 검토하세요. 경계를 이동해 자동 수락하거나 이를 sequencing error로 확정하지 않습니다.')
    if counts['no_path_phase_length_budget_incompatible']:
        recommendations.append('경계 위상·gap 길이·고정 편집 예산의 불일치를 확인하세요. 이 필요조건은 최소 편집 거리나 생물학적 단절의 판정이 아닙니다.')
    if counts['no_path_despite_feasible_template_length']:
        recommendations.append('가능한 template 길이가 있어도 현재 서열·연속 indel 조건의 경로는 없었습니다. 알려진 반복 구조와 독립 long-read 자료를 검토하세요.')
    recommendations.append('선택된 pair의 진단 수를 전체 샘플의 검출 민감도·부재·새 TRM 근거로 환산하지 마세요.')
    text = f'선택된 고정 pair {total:,}쌍의 기존 판정을 보존하고 경계·gap 조건을 설명했습니다.'
    return measures, rows, recommendations, text


def _mixed(summary, table):
    _require(summary.get('changes_native_classification') is False and
             summary.get('whole_source_biologically_analyzed') is False and
             summary.get('all_orientations_and_vocabularies_retained') is True,
             'mixed diagnostic scope or decision contract differs')
    total = _count(summary.get('analyses'), 'window analyses')
    pairs = _count(summary.get('historical_pairs'), 'historical pairs')
    vocabularies = summary.get('vocabularies', {})
    _require(bool(vocabularies) and all(isinstance(v,list) and v for v in vocabularies.values()), 'missing registered vocabularies')
    _require(len(table) == total and len({r['pair_id'] for r in table}) == pairs,
             'mixed diagnostic pair/analysis count differs')
    keys = {(r['pair_id'],r['vocabulary'],r['requested_strand']) for r in table}
    _require(len(keys) == total == pairs*len(vocabularies)*2, 'mixed diagnostic combinations incomplete or duplicated')
    selected = {r['record_id'] for r in table}
    _require(len(selected) == summary['selected_records'], 'mixed selected-record count differs')
    groups = []; rows = [_row('historical_pairs','기존 고정 pair',pairs,'pairs','registered_pairs'),
        _row('window_analyses','방향·어휘별 평가',total,'analyses','registered_windows')]
    for vocabulary,units in vocabularies.items():
        for strand in ('+','-'):
            subset=[r for r in table if r['vocabulary']==vocabulary and r['requested_strand']==strand]
            _require(len(subset)==pairs, 'mixed vocabulary/orientation denominator differs')
            for r in subset:
                _require(json.loads(r['primary_units']) == units and r['changes_native_classification']=='False'
                         and r['physical_telomere_confirmed']=='False', 'mixed vocabulary or inference changed')
            group=dict(vocabulary=vocabulary,requested_strand=strand,windows=len(subset))
            for target,source in [('exact_anchor_span_paths','exact_covers_original_phase_anchor_span'),
                                  ('edited_anchor_span_paths','edited_covers_original_phase_anchor_span'),
                                  ('exact_full_window_paths','exact_full_window_path'),
                                  ('edited_full_window_paths','edited_full_window_path')]:
                _require(all(r[source] in ('True','False') for r in subset), 'invalid mixed path flag')
                group[target]=sum(r[source]=='True' for r in subset)
                label={'exact_anchor_span_paths':'정확 단위의 기존 경계 span 경로',
                       'edited_anchor_span_paths':'편집 허용 단위의 같은 span 경로',
                       'exact_full_window_paths':'정확 단위의 전체 window 경로',
                       'edited_full_window_paths':'편집 허용 단위의 전체 window 경로'}[target]
                rows.append(_row(vocabulary+':'+strand+':'+target,vocabulary+' / '+strand+' / '+label,
                    group[target],'windows',f'{vocabulary}/{strand}: {pairs} registered windows'))
            groups.append(group)
    recommendations=['등록된 어휘와 양방향 결과를 함께 검토하세요. 두 단위를 허용한 경로가 실제로 두 단위를 모두 쓰는 qualifying mixed array라는 뜻은 아닙니다.',
        '단위당 편집 모델과 고정 gap 예산 모델은 가정이 다릅니다. 기존 phase 판정을 바꾸지 않으며 구조를 확정하려면 독립 자료가 필요합니다.']
    return (dict(kind='registered_mixed_windows',registered_pairs=pairs,analyses=total,groups=groups),rows,recommendations,
        f'기존 {pairs:,}쌍을 등록된 어휘·양방향으로 {total:,}회 설명했습니다. 전체 샘플의 검출 평가는 아닙니다.')


def _result(status, display, summary):
    return dict(schema_version=SCHEMA_VERSION, reporting_version=REPORTING_VERSION,status=status,
        display=display,scope=SCOPE,summary=summary,recommendations=[],rows=[],measurements=None,provenance=None,
        changes_native_decision=False,whole_sample_inference_permitted=False)


def load_followup_diagnostics(sample, meta, raw, source_hash, directory, analysis_complete=True):
    """Return supplemental observations; never modify the supplied native result."""
    if directory is None:
        return _result('not_requested','후속 진단 미요청','후속 진단 파일을 지정하지 않았습니다. 미평가는 0건 관측이 아닙니다.')
    root=Path(directory)
    _require(root.is_dir(), 'followup diagnostic directory does not exist')
    _require(isinstance(sample,str) and sample not in ('','.','..') and Path(sample).name==sample,
             'invalid followup sample identifier')
    folder=root/sample
    if not folder.exists():
        return _result('not_available','후속 진단 자료 없음','이 샘플에 연결할 후속 진단이 없습니다. 미평가는 0건 관측이 아닙니다.')
    _require(folder.is_dir() and (folder/'status.json').is_file(), 'followup diagnostic sidecar is partial')
    content,status_binding=_read(folder/'status.json');status=json.loads(content)
    _identity(status,sample,meta,raw,source_hash)
    state=status.get('status')
    if state=='failed':
        _require(type(status.get('exit_code')) is int and status['exit_code']!=0 and
                 isinstance(status.get('failure_reason'),str) and status['failure_reason'].strip(),
                 'failed diagnostic requires an explicit failure reason and nonzero exit code')
        failure,binding=_artifact(status['failure_artifact'])
        _require(failure.get('status') in ('failed','operational_failure','incomplete_operational_failure'),
                 'declared failure artifact does not describe a failure')
        result=_result('failed','후속 진단 실행 실패','후속 진단의 처리 실패입니다. 생물학적 미검출로 해석하지 않습니다.')
        result['recommendations']=['기록된 후속 진단 실패를 해결한 뒤 해당 단계만 다시 실행하세요. 기존 native 결과와 충분성 판정은 유지됩니다.']
        result['provenance']=dict(status_binding=status_binding,failure_binding=binding,
            failure_reason=status['failure_reason'],native_analysis_complete=bool(analysis_complete))
        return result
    _require(state=='complete' and type(status.get('exit_code')) is int and status['exit_code']==0,
             'followup diagnostic status is not complete or failed')
    _require((folder/'evidence.json').is_file(),'followup diagnostic evidence is missing')
    content,evidence_binding=_read(folder/'evidence.json')
    _require(status.get('evidence_sha256')==evidence_binding['sha256'],'followup evidence hash differs')
    evidence=json.loads(content)
    _require(evidence.get('schema_version')==1 and evidence.get('kind') in ('fixed_pair_bridge','registered_mixed_windows'),
             'unsupported followup evidence schema or kind')
    summary,summary_binding=_artifact(evidence['artifacts']['summary'])
    table,table_binding=_artifact(evidence['artifacts']['table'],table=True)
    review_binding=_review(evidence['artifacts']['review'],[summary_binding,table_binding])
    _require(summary.get('status')=='complete','bound diagnostic summary is not complete')
    relation=_relation(evidence['source_relation'],summary,raw,evidence['kind'])
    if evidence['kind']=='fixed_pair_bridge':
        _require(_contains_binding(summary.get('outputs', []),table_binding),
                 'diagnostic summary does not bind its pair table')
        measurements,rows,recommendations,text=_phase(summary,table,evidence.get('diagnostic_source_id'))
    else:
        measurements,rows,recommendations,text=_mixed(summary,table)
    rows[:0]=[_row('original_records','원래 전체 입력 records',relation['original_records'],'records','original_whole_input'),
              _row('compact_records','후속 진단 입력에 포함된 full records',relation['compact_records'],'records','selected_full_records'),
              _row('diagnostic_selected_records','고정 pair/window가 있는 records',relation['diagnostic_selected_records'],'records','diagnostic_selected_records')]
    result=_result('complete','등록된 부분 자료의 후속 진단 완료',text)
    result.update(rows=rows,measurements=measurements,recommendations=recommendations,
        provenance=dict(status_binding=status_binding,evidence_binding=evidence_binding,
            diagnostic_summary_binding=summary_binding,diagnostic_table_binding=table_binding,
            diagnostic_review_binding=review_binding,source_relation=relation,
            native_result_sha256=source_hash,native_analysis_complete=bool(analysis_complete),
            historical_diagnostic_runtime=summary.get('runtime'),
            data_role=evidence.get('data_role','unspecified')))
    if not analysis_complete:
        result['recommendations'].append('Native 분석은 완료되지 않았습니다. 이 부분 자료 진단으로 입력 실패나 충분성 제한을 해소하지 않습니다.')
    return result
