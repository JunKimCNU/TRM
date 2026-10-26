"""Join saved whole-input quality/mixed observations; never classify a TRM.

All sources are small native JSON files. No FASTQ is opened by this module.
The manifest seal and source hashes establish integrity, not biological truth.
"""
from copy import deepcopy
from pathlib import Path
import re

from evidence_io import EvidenceIO
from quality_observation import describe_quality
from report_mixed_arrays import GROUP_COUNTS, check_previous_inputs, validate_source

MAX_BYTES = 16 * 1024 * 1024
MAX_SAMPLES = 1000
MAX_TOTAL_BYTES = 256 * 1024 * 1024
SUPPORTED_NATIVE_VERSION = '0.8.5'
SHA = re.compile(r'^[0-9a-f]{64}$')
FIELDS = ('saved_structure_status', 'saved_quality_display', 'saved_mixed_display')
LIMITATIONS = [
    'Whole-input observations in these saved sources are separate from previously selected interval models.',
    'A single constant quality value cannot support discriminating-base confidence; its origin is unresolved.',
    'Multiple quality bins do not establish authentic qualities, calibrated error rates or reliable variant calls.',
    'Mixed and component family counts overlap and must not be added as independent reads or molecules.',
    'Exact supplied-unit grammars can miss errors, indels or other structures; zero is model-specific.',
    'Boundary orientation and read-end length are not physical chromosome-terminal evidence.',
    'No source or displayed summary changes native counts, candidate ranking, classification or adequacy.',
]


class EvidenceError(ValueError):
    pass


IO = EvidenceIO(EvidenceError, 'saved repeat structure: ')


def payload_sha256(value):
    return IO.payload_sha256(value)


def _keys(value, keys, where):
    if not isinstance(value, dict) or set(value) != set(keys):
        IO.fail(where + ' has missing or unsupported fields')


def _sha(value):
    if not isinstance(value, str) or not SHA.fullmatch(value):
        IO.fail('invalid SHA256')


def _read(path):
    data = IO.read_bytes(path, max_bytes=MAX_BYTES)
    value, digest = IO.decode_json(data, path, object_pairs_hook=IO.unique_object,
                                   parse_constant=IO.reject_constant)
    return value, digest, len(data)


def load_evidence(path):
    if path is None:
        return None
    path = Path(path).absolute()
    manifest, digest, total = _read(path)
    _keys(manifest, ('schema_version', 'samples', 'payload_sha256'), 'manifest')
    if type(manifest['schema_version']) is not int or manifest['schema_version'] != 1:
        IO.fail('unsupported schema')
    if payload_sha256(manifest) != manifest['payload_sha256']:
        IO.fail('manifest payload hash differs')
    rows = manifest['samples']
    if not isinstance(rows, list) or not 1 <= len(rows) <= MAX_SAMPLES:
        IO.fail('sample list outside supported bounds')
    result = {'manifest_sha256': digest, 'paths': [str(path)], 'samples': {}}
    for row in rows:
        _keys(row, ('sample_id', 'identity', 'original_native_sha256',
                    'quality_result', 'mixed_result'), 'sample')
        sample = row['sample_id']
        if not isinstance(sample, str) or not sample.strip() or sample in result['samples']:
            IO.fail('empty or duplicate sample ID')
        _keys(row['identity'], ('scientific_name', 'tax_id', 'biosample'), 'identity')
        if any(not isinstance(v, str) or not v.strip() for v in row['identity'].values()):
            IO.fail('identity values must be nonempty strings')
        _sha(row['original_native_sha256'])
        value = deepcopy(row)
        if row['quality_result'] is None and row['mixed_result'] is None:
            IO.fail('registered sample has no measurements')
        for kind in ('quality_result', 'mixed_result'):
            pin = row[kind]
            if pin is None:
                continue
            _keys(pin, ('path', 'sha256'), kind)
            _sha(pin['sha256'])
            if not isinstance(pin['path'], str) or not pin['path'].strip():
                IO.fail('source path must be nonempty')
            source = Path(pin['path'])
            if not source.is_absolute():
                source = path.parent / source
            document, source_hash, size = _read(source)
            total += size
            if total > MAX_TOTAL_BYTES:
                IO.fail('total sources exceed size bound')
            if source_hash != pin['sha256']:
                IO.fail(kind + ' source hash differs')
            value[kind] = dict(path=str(source.absolute()), sha256=source_hash, document=document)
            result['paths'].append(str(source.absolute()))
        result['samples'][sample] = value
    result['paths'] = sorted(set(result['paths']))
    return result


def source_dependencies(path):
    return load_evidence(path)['paths']


def _same_input(source, raw):
    if any(report.get('version') != SUPPORTED_NATIVE_VERSION for report in (source, raw)):
        IO.fail('unsupported native version; reviewed 0.8.5 sources required')
    if source.get('program') != 'trm' or raw.get('program') != 'trm':
        IO.fail('unsupported native program')
    for report in (source, raw):
        inputs = report.get('inputs')
        if not isinstance(inputs, list) or not inputs:
            IO.fail('nonempty input provenance required')
        for entry in inputs:
            if (not isinstance(entry, dict) or not isinstance(entry.get('path'), str)
                    or not entry['path'].strip()):
                IO.fail('invalid input path provenance')
            for key in ('bytes', 'modified_unix_seconds'):
                if type(entry.get(key)) is not int or entry[key] < 0:
                    IO.fail('invalid input size or timestamp provenance')
        summary = report.get('summary') or {}
        if summary.get('sampled_prefix') is not False:
            IO.fail('full-input EOF summary required')
        for key in ('records', 'bases'):
            if type(summary.get(key)) is not int or summary[key] < 0:
                IO.fail('invalid input count')
    if any(source['summary'][key] != raw['summary'][key] for key in ('records', 'bases')):
        IO.fail('input counts differ')
    check_previous_inputs(source, raw)


def _quality(source, digest):
    if source.get('analysis_kind') not in {'targeted_evidence', 'de_novo_discovery'}:
        IO.fail('quality source is not a supported native measurement')
    quality = describe_quality(source, digest)
    value = quality['diagnostics']
    if value is not None and quality['input_scope'] != 'complete_inputs':
        IO.fail('quality input scope is not complete')
    histogram = value['input_records']['whole_read'] if value else None
    bins = ({str(i): count for i, count in enumerate(histogram['phred33_counts']) if count}
            if histogram else {})
    # Constant refers only to quality-bearing analyzed bases. Partial quality
    # must never be labelled as constant across every input base.
    distribution = ('unmeasured' if not bins else
                    'single_value' if len(bins) == 1 else 'multiple_values')
    complete = quality['status'] == 'measured'
    label = ('Quality 미측정' if not bins else
             f"Quality 관측값 {len(bins)}개 · 실제 오류율 보정 미확인" if len(bins)>1 else
             f"{'전체 입력' if complete else 'quality 보유 부분'}의 quality가 Q{next(iter(bins))} 하나로 고정 · 변이 신뢰도 근거로 사용 불가")
    return {'status': quality['status'], 'distribution': distribution,
            'conditional_phred33_counts': bins,
            'quality_bases': histogram['bases'] if histogram else None,
            'all_input_quality_available': complete,
            'constant_across_complete_input': len(bins) == 1 and complete,
            'input_scope': quality['input_scope'], 'source_sha256': digest,
            'encoding': quality['encoding'], 'quality_origin': 'unresolved',
            'quality_calibration_verified': False, 'display': label}


def build_evidence(sample, document, raw, loaded):
    base = {'schema_version': 1, 'status': 'not_provided', 'quality': None, 'mixed': None,
            'changes_classification_or_adequacy': False, 'limitations': list(LIMITATIONS),
            'display': {'saved_structure_status': '저장 추가 관측 미제공',
                'saved_quality_display': '별도 whole-input quality 미평가',
                'saved_mixed_display': '별도 혼합 배열 미평가'}}
    if loaded is None or sample not in loaded['samples']:
        return base
    row = loaded['samples'][sample]
    if not document.get('analysis_complete') or document.get('sample_id') != sample:
        IO.fail('source assessment is incomplete or sample differs')
    if row['original_native_sha256'] != document.get('provenance', {}).get('original_discovery_sha256'):
        IO.fail('original native result hash differs')
    meta = document.get('metadata') or {}
    identity = {'scientific_name': document.get('scientific_name'),
                'tax_id': str(document.get('tax_id', '')),
                'biosample': meta.get('sample_accession')}
    if row['identity'] != identity:
        IO.fail('sample/taxon/BioSample identity differs')
    base.update(status='sources_validated', manifest_sha256=loaded['manifest_sha256'],
                original_native_sha256=row['original_native_sha256'], source_artifacts={})
    base['display']['saved_structure_status'] = '저장 자료 검증됨 · 항목별 측정 범위 참조'
    for kind in ('quality_result', 'mixed_result'):
        source = row[kind]
        if source is None:
            continue
        measured = source['document']
        _same_input(measured, raw)
        base['source_artifacts'][kind] = {key: source[key] for key in ('path', 'sha256')}
        if kind == 'quality_result':
            base['quality'] = _quality(measured, source['sha256'])
            base['display']['saved_quality_display'] = base['quality']['display']
        else:
            validate_source(measured)
            block = measured['mixed_arrays']
            if block.get('mode') != 'reads' or block['settings'].get('exact_only') is not True:
                IO.fail('only saved exact mixed reads models are supported')
            if len(block['groups']) > 64:
                IO.fail('too many mixed groups')
            groups = [{key: deepcopy(g[key]) for key in ('name', 'units', *GROUP_COUNTS)} for g in block['groups']]
            base['mixed'] = {'scope': 'complete_inputs_explicit_exact_unit_grammar',
                'input_records': measured['summary']['records'], 'input_bases': measured['summary']['bases'],
                'settings': deepcopy(block['settings']), 'groups': groups}
            base['display']['saved_mixed_display'] = '; '.join(
                f"{'+'.join(g['units'])}: mixed {g['mixed_records']}, full {g['full_mixed_records']}, "
                f"max {g['max_array_bp']} bp, 양쪽 flank {g['clean_two_flank_records']}" for g in groups)
    return base
