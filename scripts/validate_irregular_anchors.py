#!/usr/bin/env python3
"""Localize exported irregular arrays without promoting TRM or absence calls.

The detector owns the edit model. This adapter checks its geometry/accounting,
preserves the model evidence, and reuses the existing adjacent-anchor mapping
rules. A long modeled span is never labelled a long exact repeat tract.
"""
import argparse
from collections import Counter
from contextlib import ExitStack
from datetime import datetime, timezone
import csv
import html
import json
import math
from pathlib import Path
import re
import sqlite3
import sys
import tempfile
import time

from anchor_store import AnchorStore, apply_memory_limit

from validate_terminal_anchors import (
    DEFAULTS, anchor_geometry, assembly_compartments, digest, fasta_lengths,
    parse_paf, select_alignment, validate_receipt,
)

PROGRAM = 'trm-irregular-anchor-context'
SCHEMA_VERSION = 1
SAFE_ID = re.compile(r'[A-Za-z0-9_.:-]+\Z')
SHA256 = re.compile(r'[0-9a-f]{64}\Z')
MODEL_COUNTS = ('matches', 'substitutions', 'insertions', 'deletions',
                'model_bases', 'exact_unit_bases')


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('duplicate JSON key: ' + key)
        result[key] = value
    return result


def read_json(path):
    return json.loads(Path(path).read_text(), object_pairs_hook=unique_object)


def integer(value, label, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f'{label} must be an integer >= {minimum}')
    return value


def identifier(value, label):
    if not isinstance(value, str) or not SAFE_ID.fullmatch(value):
        raise ValueError('invalid safe identifier: ' + label)
    return value


def verify_file(path, record, size_key):
    path = Path(path)
    if not isinstance(record, dict) or not SHA256.fullmatch(str(record.get('sha256', ''))):
        raise ValueError('missing/invalid SHA256: ' + str(path))
    expected_size = integer(record.get(size_key), str(path) + ' size')
    before = path.stat()
    if not path.is_file() or before.st_size != expected_size or digest(path) != record['sha256']:
        raise ValueError('file size/hash mismatch: ' + str(path))
    after = path.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
            after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns):
        raise ValueError('input changed during verification: ' + str(path))


def validate_array(row, sources, group_ids, anchor_lengths):
    """Validate metadata/geometry only; the scanner's trace is not reinterpreted."""
    if not isinstance(row, dict):
        raise ValueError('array row must be an object')
    array_id = identifier(row.get('array_id'), 'array_id')
    source = integer(row.get('source_index'), 'source_index', 1)
    if source > len(sources):
        raise ValueError('source_index outside input list')
    ordinal = integer(row.get('read_ordinal'), 'read_ordinal', 1)
    if ordinal > sources[source - 1]['records']:
        raise ValueError('read_ordinal outside source record count')
    if not isinstance(row.get('read_id'), str) or not row['read_id']:
        raise ValueError('read_id must be nonempty')
    n = integer(row.get('read_length'), 'read_length', 1)
    if n > sources[source - 1]['bases']:
        raise ValueError('read_length exceeds all bases in its source')
    if row.get('group_id') not in group_ids:
        raise ValueError('unknown definition group_id')
    start = integer(row.get('start0'), 'start0')
    end = integer(row.get('end0'), 'end0', 1)
    length = integer(row.get('span_bp'), 'span_bp', 1)
    if not 0 <= start < end <= n or end - start != length:
        raise ValueError('array span disagrees with read coordinates')
    if (integer(row.get('left_end_distance'), 'left_end_distance') != start or
            integer(row.get('right_end_distance'), 'right_end_distance') != n - end):
        raise ValueError('array end distances disagree with coordinates')
    if row.get('strand') not in {'+', '-', 'ambiguous'}:
        raise ValueError('invalid array strand')
    model = row.get('model')
    if not isinstance(model, dict):
        raise ValueError('model evidence object required')
    for key in MODEL_COUNTS:
        integer(model.get(key), 'model.' + key)
    if model['matches'] + model['substitutions'] + model['insertions'] != length:
        raise ValueError('model read-consuming operations disagree with span')
    if model['matches'] + model['substitutions'] + model['deletions'] != model['model_bases']:
        raise ValueError('model reference-consuming operations disagree with model_bases')
    if model['model_bases'] <= 0 or model['exact_unit_bases'] > model['matches']:
        raise ValueError('invalid exact unit/model base counts')
    identity = model.get('identity')
    if (isinstance(identity, bool) or not isinstance(identity, (int, float)) or
            not math.isfinite(identity) or not 0 <= identity <= 1):
        raise ValueError('invalid model identity')
    if (type(model.get('optimal_parse_ambiguous')) is not bool or
            not isinstance(model.get('trace'), list) or not model['trace']):
        raise ValueError('model parse ambiguity and trace list required')
    names = []
    for side, suffix in [('left', 'L'), ('right', 'R')]:
        anchor = row.get(side + '_anchor')
        if not isinstance(anchor, dict):
            raise ValueError('missing adjacent anchor descriptor: ' + side)
        a = integer(anchor.get('start0'), side + ' anchor start')
        b = integer(anchor.get('end0'), side + ' anchor end')
        size = integer(anchor.get('length'), side + ' anchor length')
        if not 0 <= a <= b <= n or b - a != size or size > 1000:
            raise ValueError('invalid adjacent anchor interval')
        if side == 'left' and b != start or side == 'right' and a != end:
            raise ValueError('anchor is not directly adjacent to array')
        if not isinstance(anchor.get('status'), str) or not isinstance(anchor.get('stop_reason'), str):
            raise ValueError('anchor status/stop reason required')
        name = anchor.get('id')
        if name:
            if (name != array_id + '_' + suffix or anchor['status'] != 'exported' or
                    not 200 <= size <= 1000 or anchor_lengths.get(name) != size):
                raise ValueError('anchor ID/status/FASTA length mismatch')
            names.append(name)
        elif name != '' or anchor['status'] == 'exported':
            raise ValueError('nonexported anchors need an empty ID and explicit status')
    return names


def load_export(directory, store):
    directory = Path(directory)
    summary = read_json(directory / 'summary.json')
    if (type(summary.get('schema_version')) is not int or summary['schema_version'] != 1 or
            summary.get('program') != 'scan_irregular_arrays' or
            summary.get('model_version') != '0.1.0' or summary.get('status') != 'complete' or
            summary.get('complete_input') is not True or summary.get('validated_eof') is not True or
            summary.get('sampled_prefix') is not False):
        raise ValueError('completed supported full-input irregular export required')
    if summary.get('export_truncated') is True or summary.get('resource_limit_hit') is True:
        raise ValueError('capped/incomplete irregular export cannot be localized as complete')
    if not isinstance(summary.get('config'), dict):
        raise ValueError('scanner config required')
    definitions = summary.get('definitions')
    if not isinstance(definitions, list) or not definitions:
        raise ValueError('explicit oriented definitions required')
    groups = set()
    for definition in definitions:
        if not isinstance(definition, dict):
            raise ValueError('invalid definition')
        name = definition.get('id')
        if not isinstance(name, str) or not name or name in groups:
            raise ValueError('empty/duplicate definition group')
        units = definition.get('units')
        if (not isinstance(units, list) or not units or
                any(not isinstance(x, str) or not re.fullmatch('[ACGT]+', x) for x in units)):
            raise ValueError('definition requires explicit oriented DNA units')
        groups.add(name)
    sources = summary.get('inputs')
    if not isinstance(sources, list) or not sources:
        raise ValueError('input provenance required')
    for source in sources:
        if not isinstance(source, dict) or not isinstance(source.get('path'), str):
            raise ValueError('invalid input source')
        integer(source.get('records'), 'input records')
        integer(source.get('bases'), 'input bases')
        verify_file(Path(source['path']), source, 'size_bytes')
    for key in ['records', 'bases']:
        if integer(summary.get(key), 'summary.' + key) != sum(x[key] for x in sources):
            raise ValueError('source counts disagree with summary: ' + key)
    if summary['records'] == 0:
        raise ValueError('empty source is not a completed biological observation')
    outputs = summary.get('outputs', {})
    for name in ['arrays.jsonl', 'anchors.fasta']:
        verify_file(directory / name, outputs.get(name), 'size_bytes')
    store.load_fasta(directory / 'anchors.fasta')
    return summary, groups


def iter_arrays(directory, summary, groups, store):
    """Validate one full trace-bearing row at a time; identities live on disk."""
    directory = Path(directory)
    with (directory / 'arrays.jsonl').open() as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line, object_pairs_hook=unique_object)
            names = validate_array(row, summary['inputs'], groups, store.lengths)
            store.add_array(row)
            for name in names:
                store.mark_anchor(name)
            yield row
    store.finish_export(integer(summary.get('arrays'), 'summary.arrays'),
                        integer(summary.get('anchors'), 'summary.anchors'))


def evaluate_array(array, alignments, settings, compartments=None):
    """Use original-read span for geometry only; expose no legacy exact-length flag."""
    geometry = {'strand': array['strand'], 'tract_length': array['span_bp'],
                'left_end_distance': array['left_end_distance'],
                'right_end_distance': array['right_end_distance']}
    anchors = []
    for side in ['left', 'right']:
        descriptor = array[side + '_anchor']
        name = descriptor['id']
        if not name:
            continue
        chosen, status = select_alignment(alignments.get(name, []), settings)
        row = {'anchor_id': name, 'anchor_side': side, 'mapping_status': status,
               'reported_alignments': len(alignments.get(name, []))}
        if chosen:
            mapped = anchor_geometry(geometry, side, chosen, settings)
            mapped['model_span_ge_long_threshold'] = mapped.pop('tract_ge_long_threshold')
            row.update(mapped)
            row['target_compartment'] = (compartments or {}).get(chosen['target'], 'unannotated')
            if row['target_compartment'] == 'organelle':
                row['assembly_end_anchor_concordant'] = False
                row['read_end_and_assembly_end_concordant'] = False
        row['unique_nonrepeat_subtelomere_confirmed'] = False
        row['anchor_repeat_contamination_assessed'] = False
        anchors.append(row)
    qualified = [x for x in anchors if x['mapping_status'] == 'qualified_within_reported_alignments']
    concordant = [x for x in qualified if x['assembly_end_anchor_concordant']]
    joint = [x for x in qualified if x['read_end_and_assembly_end_concordant']]
    if array['strand'] == 'ambiguous':
        status = 'ambiguous_array_direction'
    elif len(qualified) == 2:
        status = 'two_sided_anchor_context'
    elif joint:
        status = 'read_end_and_assembly_end_concordant'
    elif concordant:
        status = 'assembly_end_anchor_only'
    elif qualified:
        status = 'qualified_anchor_without_terminal_concordance'
    elif anchors:
        status = 'anchor_mapping_unresolved'
    else:
        status = 'no_exportable_adjacent_anchor'
    joint_flag = len(qualified) == 1 and bool(joint) and array['strand'] != 'ambiguous'
    return {**array, 'localization_status': status, 'anchor_alignments': anchors,
            'qualified_anchor_count': len(qualified), 'two_sided_qualified_anchors': len(qualified) == 2,
            'joint_terminal_context': joint_flag,
            'joint_long_model_span_context': joint_flag and array['span_bp'] >= settings['long_tract_bp'],
            'joint_nuclear_terminal_context': joint_flag and joint[0]['target_compartment'] == 'nuclear_chromosome',
            'length_evidence_kind': 'accepted_irregular_model_span_in_original_read_coordinates',
            'left_record_edge_censored': array['start0'] == 0,
            'right_record_edge_censored': array['end0'] == array['read_length'],
            'physical_chromosome_end_confirmed': False, 'changes_native_classification': False}


def file_stamp(path):
    stat = Path(path).stat()
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def flatten_array(row):
    return {key: value for key, value in row.items() if key != 'anchor_alignments'}


def flatten_anchors(row):
    for anchor in row['anchor_alignments']:
        yield {'array_id': row['array_id'], 'source_index': row['source_index'],
               'read_ordinal': row['read_ordinal'], **anchor}


def tsv_row(row):
    return {key: json.dumps(value, ensure_ascii=False) if isinstance(value, (dict, list)) else value
            for key, value in row.items()}


def write_outputs(spool, output_dir, array_fields, anchor_fields):
    """Stream ordered rows after collecting the union of even late TSV fields."""
    with ExitStack() as stack:
        details = stack.enter_context((output_dir / 'details.json').open('w'))
        arrays = stack.enter_context((output_dir / 'arrays.tsv').open('w', newline=''))
        anchors = stack.enter_context((output_dir / 'anchors.tsv').open('w', newline=''))
        page = stack.enter_context((output_dir / 'index.html').open('w'))
        array_writer = csv.DictWriter(arrays, array_fields, delimiter='\t')
        anchor_writer = csv.DictWriter(anchors, anchor_fields, delimiter='\t')
        if array_fields:
            array_writer.writeheader()
        if anchor_fields:
            anchor_writer.writeheader()
        page.write(HTML_HEADER)
        first = True
        with spool.open() as handle:
            for line in handle:
                row = json.loads(line)
                details.write('[\n' if first else ',\n')
                details.write('  ' + json.dumps(row, ensure_ascii=False, indent=2).replace('\n', '\n  '))
                first = False
                array_writer.writerow(tsv_row(flatten_array(row)))
                for anchor in flatten_anchors(row):
                    anchor_writer.writerow(tsv_row(anchor))
                page.write('<tr>' + ''.join('<td>' + html.escape(str(row[key])) + '</td>' for key in
                    ['read_id', 'span_bp', 'strand', 'localization_status', 'joint_long_model_span_context']) + '</tr>')
        details.write('[]\n' if first else '\n]\n')
        page.write('</table></html>')


HTML_HEADER = ('<!doctype html><html lang="ko"><meta charset="utf-8">'
        '<title>Irregular array anchor context</title><h1>불규칙 배열과 같은 read의 인접 서열 위치</h1>'
        '<p>길이는 원래 read에서 모델로 설명한 span입니다. 500 bp 모델 span은 500 bp 정확 반복이 아닙니다. '
        'Assembly 끝·방향 일치로 물리적 telomere, 비반복 subtelomere, 신규성 또는 미검출 검출력을 확정하지 않습니다.</p>'
        '<p><a href="summary.json">설정·검증·한계</a> · <a href="details.json">모델·편집·정렬 전체 근거</a> · '
        '<a href="arrays.tsv">배열 TSV</a> · <a href="anchors.tsv">정렬 TSV</a></p>'
        '<table><tr><th>Read</th><th>Model span bp</th><th>Strand</th><th>Context</th><th>Long model context</th></tr>')


def localize(export_dir, assembly, assembly_report, paf, alignment_receipt, output_dir, settings,
             memory_limit_mib=None, index_cache_mib=8, scratch_dir=None):
    started = time.monotonic()
    export_dir, assembly, output_dir = map(Path, [export_dir, assembly, output_dir])
    if output_dir.exists():
        raise ValueError('output directory already exists; preserve previous results')
    if type(index_cache_mib) is not int or index_cache_mib <= 0:
        raise ValueError('index cache must be a positive integer MiB')
    resources = apply_memory_limit(memory_limit_mib)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    scratch_parent = Path(scratch_dir) if scratch_dir is not None else output_dir.parent
    # Publication staging is on the output filesystem, even with separate scratch.
    with tempfile.TemporaryDirectory(prefix='.irregular-output-', dir=output_dir.parent) as staging, \
         tempfile.TemporaryDirectory(prefix='.irregular-index-', dir=scratch_parent) as scratch:
        store = AnchorStore(Path(scratch) / 'index.sqlite', index_cache_mib)
        try:
            stage = Path(staging) / 'completed'
            stage.mkdir()
            resources.update({'index_cache_mib': index_cache_mib, 'index_backend': 'sqlite_disk',
                              'scratch_parent': str(scratch_parent.resolve()),
                              'single_array_in_memory': True,
                              'whole_input_constant_memory_guaranteed': False})
            report = _localize(export_dir, assembly, assembly_report, paf, alignment_receipt,
                               stage, settings, store, Path(scratch), resources, started)
            store.close()
            if output_dir.exists():
                raise ValueError('output directory appeared during run; refusing to replace it')
            stage.rename(output_dir)
            return report
        finally:
            store.close()


def _localize(export_dir, assembly, assembly_report, paf, alignment_receipt,
              output_dir, settings, store, scratch, resources, started):
    source_paths = [export_dir / 'summary.json', export_dir / 'arrays.jsonl', export_dir / 'anchors.fasta', assembly]
    if assembly_report is not None:
        source_paths.append(Path(assembly_report))
    if (paf is None) != (alignment_receipt is None):
        raise ValueError('provide both PAF and alignment receipt, or neither for zero anchors')
    if paf is not None:
        paf, alignment_receipt = Path(paf), Path(alignment_receipt)
        source_paths += [paf, alignment_receipt]
    stamps = {path: file_stamp(path) for path in source_paths}
    # Include original reads in the before/after stability check as well as the
    # export's existing SHA256 verification, without rereading sequence records.
    for source in read_json(export_dir / 'summary.json').get('inputs', []):
        if isinstance(source, dict) and isinstance(source.get('path'), str):
            path = Path(source['path'])
            stamps[path] = file_stamp(path)
    summary, groups = load_export(export_dir, store)
    assembly_lengths = fasta_lengths(assembly)
    if not assembly_lengths:
        raise ValueError('empty assembly')
    compartments = assembly_compartments(assembly_report, assembly_lengths)
    if len(store.lengths) and paf is None:
        raise ValueError('exported anchors require PAF and completed mapping receipt')
    if paf is not None:
        validate_receipt(read_json(alignment_receipt), {
            'reference': assembly, 'queries': export_dir / 'anchors.fasta', 'paf': paf})
        with paf.open() as handle:
            for line in handle:
                if not line.strip():
                    continue
                row = parse_paf(line)
                if (store.length(row['query']) != row['query_length'] or
                        assembly_lengths.get(row['target']) != row['target_length']):
                    raise ValueError('PAF identifier/length disagrees with actual FASTA')
                store.add_alignment(row)
        mapping_status = 'complete_verified_alignment'
    else:
        mapping_status = 'not_run_no_exported_anchors'
    spool = scratch / 'outcomes.jsonl'
    array_fields, anchor_fields, status_counts = {}, {}, Counter()
    array_count = 0
    with spool.open('w') as handle:
        for row in iter_arrays(export_dir, summary, groups, store):
            outcome = evaluate_array(row, store, settings, compartments)
            store.mark_outcome(outcome)
            handle.write(json.dumps(outcome, ensure_ascii=False) + '\n')
            array_fields.update(dict.fromkeys(flatten_array(outcome)))
            for anchor in flatten_anchors(outcome):
                anchor_fields.update(dict.fromkeys(anchor))
            status_counts[outcome['localization_status']] += 1
            array_count += 1
    record_counts = store.record_counts()
    store.db.commit()
    source_hashes = {str(path.resolve()): digest(path) for path in source_paths}
    for path, stamp in stamps.items():
        if file_stamp(path) != stamp:
            raise ValueError('input changed during localization: ' + str(path))
    report = {'schema_version': SCHEMA_VERSION, 'program': PROGRAM, 'status': 'complete',
              'created_utc': datetime.now(timezone.utc).isoformat(), 'settings': settings,
              'model_version': summary['model_version'], 'source_export_summary': summary,
              'arrays': array_count, 'exported_anchors': len(store.lengths),
              'array_supporting_records': record_counts[0],
              'localization_status_counts': dict(status_counts),
              'joint_terminal_context_records': record_counts[1],
              'joint_long_model_span_context_records': record_counts[2],
              'joint_nuclear_terminal_context_records': record_counts[3],
              'length_evidence_kind': 'model_span_not_exact_repeat_length',
              'mapping_status': mapping_status,
              'source_sha256': source_hashes,
              'physical_chromosome_end_confirmed': False, 'changes_native_classification': False,
              'changes_adequacy': False, 'changes_catalogue': False,
              'limitations': [
                  'The detector defines and accepts the irregular model; this adapter does not refit it or replay its edit trace.',
                  'A >=500 bp model span is not >=500 exact repeat bases. Model edits, exact-unit bases and ambiguity remain separate.',
                  'One array must satisfy its own anchor geometry; separate arrays and records are never added into a long tract.',
                  'Read ends are fragment ends and assembly ends are coordinates, not independently confirmed physical chromosome ends.',
                  'Matched reads and their assembly are dependent evidence, not independent biological replication.',
                  'Exported anchors may contain repeats. Bounded secondary alignments and MAPQ do not prove global uniqueness.',
                  'No repeat mask is inferred for this multi-unit group; nonrepeat subtelomeric identity remains unconfirmed.',
                  'Ambiguous strand, missing anchor, competing alignment and two-sided context do not receive terminal confirmation.',
                  'No catalogue, main-TRM, novelty, ALT, absence, karyotype or data-adequacy decision is changed.'
              ]}
    write_outputs(spool, output_dir, list(array_fields), list(anchor_fields))
    # Output formatting may be long; require stable sources through completion.
    for path, stamp in stamps.items():
        if file_stamp(path) != stamp:
            raise ValueError('input changed during output writing: ' + str(path))
    resources['scratch_bytes_before_cleanup'] = sum(path.stat().st_size for path in scratch.iterdir() if path.is_file())
    report['execution_resources'] = resources
    report['wall_seconds'] = time.monotonic() - started
    report['outputs_sha256'] = {name: digest(output_dir / name) for name in ['details.json', 'arrays.tsv', 'anchors.tsv', 'index.html']}
    (output_dir / 'summary.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['export-dir', 'assembly', 'output-dir']:
        parser.add_argument('--' + name, type=Path, required=True)
    for name in ['assembly-report', 'paf', 'alignment-receipt']:
        parser.add_argument('--' + name, type=Path)
    parser.add_argument('--memory-limit-mib', type=int,
                        help='optional process virtual-address-space cap (RLIMIT_AS), not an RSS target')
    parser.add_argument('--index-cache-mib', type=int, default=8,
                        help='SQLite page cache budget in MiB (default: 8); not a total memory limit')
    parser.add_argument('--scratch-dir', type=Path,
                        help='existing directory for temporary disk index/spool; default: output parent')
    for name, value in DEFAULTS.items():
        parser.add_argument('--' + name.replace('_', '-'), type=type(value), default=value)
    args = parser.parse_args()
    settings = {name: getattr(args, name) for name in DEFAULTS}
    if (any(not math.isfinite(v) or v < 0 for v in settings.values()) or
            any(not 0 < settings[k] <= 1 for k in ['min_identity', 'min_query_coverage', 'competing_score_fraction'])):
        parser.error('invalid geometry thresholds')
    if (args.memory_limit_mib is not None and args.memory_limit_mib <= 0) or args.index_cache_mib <= 0:
        parser.error('memory and index cache limits must be positive integers')
    # Release this reserve before reporting allocation failures under RLIMIT_AS.
    reserve = bytearray(1024 * 1024)
    try:
        result = localize(args.export_dir, args.assembly, args.assembly_report, args.paf,
                          args.alignment_receipt, args.output_dir, settings,
                          args.memory_limit_mib, args.index_cache_mib, args.scratch_dir)
    except (MemoryError, ValueError, OSError, sqlite3.Error) as exc:
        del reserve
        failure = {'program': PROGRAM, 'status': 'failed', 'complete': False,
                   'error_type': type(exc).__name__, 'error': str(exc),
                   'requested_memory_limit_mib': args.memory_limit_mib,
                   'biological_interpretation_available': False}
        # Never overwrite an existing result, even when that caused the failure.
        try:
            args.output_dir.mkdir(parents=True, exist_ok=False)
        except OSError:
            pass
        else:
            try:
                (args.output_dir / 'failure.json').write_text(json.dumps(failure, indent=2) + '\n')
            except OSError:
                pass
        print(json.dumps(failure), file=sys.stderr)
        return 1
    print(json.dumps({key: result[key] for key in ['arrays', 'exported_anchors', 'joint_terminal_context_records',
                                                 'joint_long_model_span_context_records', 'wall_seconds']}))
    return 0


if __name__ == '__main__':
    sys.exit(main())
