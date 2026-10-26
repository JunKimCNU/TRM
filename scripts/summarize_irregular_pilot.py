#!/usr/bin/env python3
"""Aggregate frozen irregular-pilot outputs; never read or scan source sequences.

The exact-mixed, homogeneous and edit-model measurements deliberately remain
separate. This report does not estimate sensitivity, specificity or speedup.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BASE = ROOT / 'results/irregular_pilot_20260926'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


class Artifacts:
    def __init__(self):
        self.sources = {}

    def bind(self, path, expected=None):
        path = Path(path).resolve()
        value = digest(path)
        if expected is not None and expected != value:
            raise ValueError('artifact SHA256 mismatch: ' + str(path))
        self.sources[str(path)] = {'sha256': value, 'size_bytes': path.stat().st_size}
        return path

    def read(self, path, expected=None):
        return json.loads(self.bind(path, expected).read_text())


def union_bp(intervals):
    end, total = None, 0
    for start, stop in sorted(intervals):
        if start >= stop:
            raise ValueError('invalid interval')
        if end is None or start > end:
            total += stop - start
        elif stop > end:
            total += stop - end
        end = stop if end is None else max(end, stop)
    return total


def overlap(a, b, start, end):
    return max(0, min(b, end) - max(a, start))


def canonical(unit):
    rc = unit.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    return min(s[i:] + s[:i] for s in [unit, rc] for i in range(len(s)))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def write_table(path, rows):
    with Path(path).open('w', newline='') as handle:
        if not rows:
            return
        fields = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(handle, fields, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow({key: ('NA' if value is None else json.dumps(value, ensure_ascii=False)
                                   if isinstance(value, (list, dict)) else value)
                             for key, value in row.items()})


def check_plan(base, artifacts):
    protocol = artifacts.read(base / 'protocol.json')
    if len(protocol['datasets']) != 6:
        raise ValueError('this report expects the registered six-dataset pilot')
    for ds in protocol['datasets']:
        artifacts.read(ROOT / ds['definition']['path'], ds['definition']['sha256'])
        for source in ds['comparators']:
            prior = artifacts.read(ROOT / source['path'], source['sha256'])
            if (prior['summary']['records'], prior['summary']['bases']) != (ds['expected_records'], ds['expected_bases']):
                raise ValueError('registered comparator denominator mismatch')
    review_spec = next(x for x in protocol['registration_evidence'] if x['path'].endswith('/review_plan/prospective_review.json'))
    review = artifacts.read(ROOT / review_spec['path'], review_spec['sha256'])
    windows_path = artifacts.bind(base / 'review_plan/proposed_internal_windows.tsv', review['coordinate_file_sha256'])
    with windows_path.open() as handle:
        windows = list(csv.DictReader(handle, delimiter='\t'))
    if len(windows) != 78:
        raise ValueError('expected all 78 proposed windows, without selection by results')
    for row in windows:
        for key in ['record_length', 'start', 'end', 'quarter']:
            row[key] = int(row[key])
        if not 0 <= row['start'] < row['end'] <= row['record_length'] or row['end'] - row['start'] != 5000:
            raise ValueError('invalid registered internal window')
    return protocol, windows


def prior_rows(ds, artifacts):
    rows = []
    definitions = artifacts.read(ROOT / ds['definition']['path'], ds['definition']['sha256'])
    units = sorted({unit for d in definitions for unit in d['units']})
    for spec in ds['comparators']:
        path = ROOT / spec['path']
        result = artifacts.read(path, spec['sha256'])
        metrics_path = path.with_suffix('.metrics.json') if path.name != 'result.json' else path.parent / 'metrics.json'
        metrics = artifacts.read(metrics_path) if metrics_path.exists() else {}
        common = {'dataset_id': ds['dataset_id'], 'species': ds['species'], 'role': ds['role'],
                  'method': spec['scope'], 'version': result['version'], 'source_path': spec['path'],
                  'records': result['summary']['records'], 'bases': result['summary']['bases'],
                  'scanner_elapsed_seconds': result['summary'].get('elapsed_seconds'),
                  'invocation_wall_seconds': metrics.get('wall_seconds'),
                  'max_rss_mib': metrics.get('max_rss_mib', metrics.get('max_rss_kib', 0) / 1024 if 'max_rss_kib' in metrics else None),
                  'threads': result.get('threads'), 'minimum_identity': result.get('config', {}).get('min_identity'),
                  'memory_limit_mib': result.get('memory_policy', {}).get('requested_limit_mib'),
                  'settings': spec.get('settings'), 'measurement_is_comparable_to_irregular_model': False}
        if 'mixed_arrays' in result:
            for group in result['mixed_arrays']['groups']:
                rows.append({**common, 'unit_or_group': group['name'], 'status': 'measured_exact_mixed',
                    'positive_records': group['mixed_records'], 'arrays': group['total_mixed_arrays'],
                    'max_span_bp': group['max_array_bp'], 'long_records': group['long_array_records'],
                    'full_records': group['full_mixed_records'], 'repeat_bases_reported': group['repeat_bases'],
                    'array_interval_union_bp': None, 'exact_mixed_group_counts': {k: v for k, v in group.items() if k != 'examples'},
                    'interpretation': 'Exact grammar, >=24 bp and guaranteed two copies of each of >=2 units; bounded examples are not a whole-assembly coordinate census.'})
        else:
            available = {}
            for key in ['candidates', 'candidate_evidence', 'rejected_candidates']:
                for row in result.get(key, []):
                    motif = row['motif']
                    if motif in available:
                        raise ValueError('duplicate canonical native family across result sections')
                    available[motif] = (key, row)
            for unit in units:
                found = available.get(canonical(unit))
                section, row = found if found else (None, {})
                rows.append({**common, 'unit_or_group': unit, 'canonical_family': canonical(unit),
                    'status': 'measured_' + section if found else 'not_reported_not_zero',
                    'positive_records': row.get('reads'), 'arrays': None, 'max_span_bp': row.get('max_tract'),
                    'long_records': row.get('long_tract_reads'), 'full_records': row.get('full_repeat_reads'),
                    'repeat_bases_reported': row.get('repeat_bases'), 'array_interval_union_bp': None,
                    'prior_class': row.get('telomere_prior', {}).get('class'),
                    'interpretation': 'Homogeneous extension can allow substitutions (identity >=0.85); rejected section is catalogue placement, not biological novelty/absence. Unit rows overlap and must not be added.'})
    return rows


def longest_edited_span(trace):
    best = current = 0
    for token in trace:
        current = 0 if token['edit'] == '=' else current + token['end0'] - token['start0']
        best = max(best, current)
    return best


def model_rows(ds, base, artifacts):
    folder = base / 'analysis' / ds['dataset_id']
    summary = artifacts.read(folder / 'output/summary.json')
    runtime = artifacts.read(folder / 'runtime.json')
    verification = artifacts.read(folder / 'verification.json')
    if (summary.get('status') != 'complete' or summary.get('program') != 'scan_irregular_arrays' or
            summary.get('complete_input') is not True or summary.get('validated_eof') is not True or
            summary.get('sampled_prefix') is not False or summary.get('resource_limit_hit') is not False or
            summary.get('export_truncated') is not False or runtime.get('exit_code') != 0):
        raise ValueError('incomplete or failed model output: ' + ds['dataset_id'])
    if (summary['records'], summary['bases']) != (ds['expected_records'], ds['expected_bases']):
        raise ValueError('model denominator mismatch')
    if verification['summary_sha256'] != digest(folder / 'output/summary.json'):
        raise ValueError('primary verification does not bind this summary')
    path = artifacts.bind(folder / 'output/arrays.jsonl', summary['outputs']['arrays.jsonl']['sha256'])
    before = path.stat()
    counts = Counter()
    maxima = Counter()
    identity_min = identity_max = None
    minimum_gate = None
    previous_key, record_rows = None, []
    assembly = []

    def finish_record(rows):
        if not rows:
            return
        counts['positive_records'] += 1
        counts['interval_union_modeled_bp'] += union_bp((x['start0'], x['end0']) for x in rows)
        counts['long_modeled_records'] += any(x['span_bp'] >= 500 for x in rows)
        counts['full_modeled_records'] += any(10 * x['span_bp'] >= 9 * x['read_length'] for x in rows)

    with path.open() as handle:
        for line in handle:
            row = json.loads(line)
            key = row['source_index'], row['read_ordinal']
            if previous_key is not None and key < previous_key:
                raise ValueError('array file not in deterministic source/read order')
            if key != previous_key:
                finish_record(record_rows)
                record_rows = []
            previous_key = key
            if not 0 <= row['start0'] < row['end0'] <= row['read_length'] or row['span_bp'] != row['end0'] - row['start0']:
                raise ValueError('invalid model coordinates')
            model = row['model']
            counts['arrays'] += 1
            counts['modeled_span_bp_sum'] += row['span_bp']
            counts['long_modeled_arrays'] += row['span_bp'] >= 500
            counts['full_modeled_arrays'] += 10 * row['span_bp'] >= 9 * row['read_length']
            counts['orientation_ambiguous_arrays'] += row['strand'] == 'ambiguous'
            counts['parse_ambiguous_arrays'] += model['optimal_parse_ambiguous']
            counts['anchors'] += bool(row['left_anchor']['id']) + bool(row['right_anchor']['id'])
            for field in ['matches', 'substitutions', 'insertions', 'deletions', 'model_bases', 'exact_unit_bases']:
                counts[field + '_sum'] += model[field]
            maxima['max_model_span_bp'] = max(maxima['max_model_span_bp'], row['span_bp'])
            maxima['max_consecutive_edited_units'] = max(maxima['max_consecutive_edited_units'], model['max_consecutive_edited_units'])
            maxima['max_edited_only_model_span_bp'] = max(maxima['max_edited_only_model_span_bp'], longest_edited_span(model['trace']))
            maxima['max_unmatched_query_run'] = max(maxima['max_unmatched_query_run'], model['max_unmatched_query_run'])
            identity_min = model['identity'] if identity_min is None else min(identity_min, model['identity'])
            identity_max = model['identity'] if identity_max is None else max(identity_max, model['identity'])
            minimum_gate = model['literal_gate_hits'] if minimum_gate is None else min(minimum_gate, model['literal_gate_hits'])
            slim = {k: v for k, v in row.items() if k not in ['model', 'left_anchor', 'right_anchor']}
            slim['model'] = {k: v for k, v in model.items() if k != 'trace'}
            slim['max_edited_only_model_span_bp'] = longest_edited_span(model['trace'])
            record_rows.append(slim)
            if ds['modality'] == 'assembly':
                assembly.append(slim)
    finish_record(record_rows)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
        raise ValueError('model output changed during report aggregation')
    for key, expected in [('arrays', summary['arrays']), ('anchors', summary['anchors']),
                          ('positive_records', summary['counts']['records_with_arrays']),
                          ('long_modeled_arrays', summary['counts']['arrays_model_span_ge_500']),
                          ('full_modeled_arrays', summary['counts']['arrays_covering_ge_90_percent_of_record'])]:
        if counts[key] != expected:
            raise ValueError('recount disagrees with model summary: ' + key)
    denominator = sum(counts[x + '_sum'] for x in ['matches', 'substitutions', 'insertions', 'deletions'])
    numeric_fields = ['positive_records', 'interval_union_modeled_bp', 'long_modeled_records', 'full_modeled_records',
                      'arrays', 'modeled_span_bp_sum', 'long_modeled_arrays', 'full_modeled_arrays',
                      'orientation_ambiguous_arrays', 'parse_ambiguous_arrays', 'anchors']
    values = {key: counts[key] for key in numeric_fields}
    values.update({key + '_sum': counts[key + '_sum'] for key in ['matches', 'substitutions', 'insertions', 'deletions', 'model_bases', 'exact_unit_bases']})
    values.update({key: maxima[key] for key in ['max_model_span_bp', 'max_consecutive_edited_units', 'max_edited_only_model_span_bp', 'max_unmatched_query_run']})
    row = {'dataset_id': ds['dataset_id'], 'species': ds['species'], 'role': ds['role'], 'modality': ds['modality'],
           'status': 'complete_full_input', 'records': summary['records'], 'bases': summary['bases'], **values,
           'identity_min': identity_min, 'identity_max': identity_max,
           'identity_from_summed_model_operations': counts['matches_sum'] / denominator if denominator else None,
           'minimum_array_literal_gate_hits': minimum_gate,
           'orientations_considered': summary['counts']['orientations_considered'],
           'orientations_unscored_literal_gate': summary['counts']['orientations_unscored_literal_gate'],
           'invocation_wall_seconds': runtime['wall_seconds'], 'max_rss_mib': runtime['max_rss_mib'],
           'preflight_seconds': summary['timing']['preflight_seconds'],
           'scan_and_emit_seconds': summary['timing']['analysis_seconds'],
           'postflight_verification_seconds': summary['timing']['postflight_verification_seconds'],
           'megabases_per_scan_and_emit_second': summary['bases'] / 1e6 / summary['timing']['analysis_seconds'],
           'definitions': summary['definitions'], 'model_version': summary['model_version'],
           'group_counts': summary['counts']['groups'], 'truth_limitations': ds['truth_limitations']}
    return row, assembly


def assembly_end_rows(ds, arrays):
    records = {x['id']: x['length'] for x in ds['nuclear_records']}
    grouped = defaultdict(list)
    for row in arrays:
        if records.get(row['read_id']) != row['read_length']:
            raise ValueError('array assembly record ID/length not in registered nuclear records')
        grouped[row['read_id']].append(row)
    ends, locations = [], []
    for row in arrays:
        n, start, stop = row['read_length'], row['start0'], row['end0']
        locations.append({'dataset_id': ds['dataset_id'], **row,
            'coordinate_inset_beyond_5000bp': min(start, n - stop) > 5000,
            'crosses_left_5000_boundary': start < 5000 < stop,
            'crosses_right_5000_boundary': start < n - 5000 < stop,
            'physical_interstitial_origin_confirmed': False})
    for name, n in records.items():
        for side, strand in [('left', '-'), ('right', '+')]:
            distance = lambda row: row['start0'] if side == 'left' else n - row['end0']
            expected = [row for row in grouped[name] if row['strand'] == strand]
            near50 = [row for row in expected if distance(row) <= 50]
            near5000 = [row for row in expected if distance(row) <= 5000]
            nearest = min(expected, key=lambda row: (distance(row), -row['span_bp'], row['array_id']), default=None)
            resolution = 'not_established_by_this_protocol'
            if ds['dataset_id'] == 'GCA_032883015.1':
                resolution = ('known_unresolved_rDNA_end' if name == 'CP133361.1' and side == 'left'
                              else 'literature_resolved_end')
            ends.append({'dataset_id': ds['dataset_id'], 'role': ds['role'], 'record_id': name,
                'record_length': n, 'end': side, 'expected_strand': strand,
                'literature_endpoint_resolution':resolution, 'biological_negative_truth':False,
                'expected_oriented_within50': bool(near50), 'expected_oriented_within5000': bool(near5000),
                'expected_arrays_within50': len(near50), 'expected_arrays_within5000': len(near5000),
                'long_model_expected_arrays_within50': sum(x['span_bp'] >= 500 for x in near50),
                'long_model_expected_arrays_within5000': sum(x['span_bp'] >= 500 for x in near5000),
                'max_expected_span_within50': max((x['span_bp'] for x in near50), default=0),
                'max_expected_span_within5000': max((x['span_bp'] for x in near5000), default=0),
                'ambiguous_arrays_within50': sum(x['strand'] == 'ambiguous' and distance(x) <= 50 for x in grouped[name]),
                'opposite_arrays_within50': sum(x['strand'] not in {strand, 'ambiguous'} and distance(x) <= 50 for x in grouped[name]),
                'nearest_expected_array_distance': distance(nearest) if nearest else None,
                'nearest_expected_array': nearest, 'original_full_record_coordinates': True,
                'physical_endpoint_confirmed': False})
    return ends, locations


def window_rows(windows, protocol, assemblies):
    by_path = {str((ROOT / ds['inputs'][0]['path']).resolve()): ds for ds in protocol['datasets'] if ds['modality'] == 'assembly'}
    rows = []
    for window in windows:
        ds = by_path.get(str((ROOT / window['assembly_path']).resolve()))
        common = {**window, 'dataset_id': ds['dataset_id'] if ds else None,
                  'scope': 'coordinate-only overlap of whole-assembly outputs', 'crop_edges_are_endpoints': False,
                  'negative_biological_truth': False, 'physical_endpoint': False}
        if ds is None:
            rows.append({**common, 'status': 'unmeasured_assembly_outside_registered_six',
                         'overlapping_arrays': None, 'union_modeled_bp': None, 'union_fraction': None})
            continue
        nuclear = {record['id']:record['length'] for record in ds['nuclear_records']}
        if nuclear.get(window['record_id']) != window['record_length']:
            raise ValueError('registered window ID/length disagrees with nuclear record manifest')
        hits = [a for a in assemblies[ds['dataset_id']] if a['read_id'] == window['record_id'] and
                overlap(a['start0'], a['end0'], window['start'], window['end'])]
        for a in hits:
            if a['read_length'] != window['record_length']:
                raise ValueError('window and whole-assembly lengths disagree')
        clips = [(max(a['start0'], window['start']), min(a['end0'], window['end'])) for a in hits]
        covered = union_bp(clips)
        rows.append({**common, 'status': 'measured_original_coordinates', 'overlapping_arrays': len(hits),
            'union_modeled_bp': covered, 'union_fraction': covered / (window['end'] - window['start']),
            'arrays_crossing_window_boundary': sum(a['start0'] < window['start'] or a['end0'] > window['end'] for a in hits),
            'max_full_array_span_bp': max((a['span_bp'] for a in hits), default=0),
            'max_window_clipped_overlap_bp': max((b-a for a,b in clips), default=0),
            'plus_arrays': sum(a['strand'] == '+' for a in hits), 'minus_arrays': sum(a['strand'] == '-' for a in hits),
            'ambiguous_arrays': sum(a['strand'] == 'ambiguous' for a in hits),
            'parse_ambiguous_arrays': sum(a['model']['optimal_parse_ambiguous'] for a in hits),
            'array_ids': [a['array_id'] for a in hits]})
    return rows


def anchor_rows(protocol, base, artifacts, require):
    state_path = base / 'anchor_run_summary.json'
    state = artifacts.read(state_path) if state_path.exists() else {'status': 'not_run'}
    if require and state.get('status') != 'complete':
        raise ValueError('requested complete anchor stage is not complete')
    rows = []
    for sid, aid in protocol['optional_anchor_validation']['assembly_pairs'].items():
        path = base / 'anchor_context' / sid / 'localization/summary.json'
        common = {'dataset_id': sid, 'assembly_id': aid, 'independent_biological_replicate': False,
                  'physical_chromosome_end_confirmed': False, 'model_span_is_exact_tract': False}
        if not path.exists():
            rows.append({**common, 'status': 'unmeasured_anchor_stage_' + state.get('status', 'unknown'),
                         'joint_terminal_context_records': None, 'joint_long_model_span_context_records': None})
            continue
        runtime = next((x for x in state.get('datasets', []) if x['dataset_id'] == sid), None)
        if runtime is None:
            if require:
                raise ValueError('anchor stage has no completed dataset receipt: ' + sid)
            rows.append({**common,'status':'unmeasured_anchor_dataset_receipt_pending',
                         'joint_terminal_context_records':None,'joint_long_model_span_context_records':None})
            continue
        if (ROOT / runtime['summary_path']).resolve() != path.resolve():
            raise ValueError('anchor stage summary path disagrees with dataset')
        summary = artifacts.read(path, runtime['summary_sha256'])
        if summary.get('status') != 'complete' or summary.get('changes_native_classification') is not False:
            raise ValueError('invalid anchor-context summary')
        # Bind observed output/PAF/receipt provenance without rereading the raw WGS
        # or assembly sequence. Raw assembly identity is compared to registration.
        export_summary = base / 'analysis' / sid / 'output/summary.json'
        check_paths = [export_summary]
        if summary['mapping_status'] != 'not_run_no_exported_anchors':
            check_paths += [path.parent.parent / 'anchors.paf',path.parent.parent / 'alignment_receipt.json']
        for source in check_paths:
            expected = summary.get('source_sha256', {}).get(str(source.resolve()))
            if expected is None:
                raise ValueError('anchor summary omits required source provenance')
            artifacts.bind(source, expected)
        assembly_ds = next(ds for ds in protocol['datasets'] if ds['dataset_id'] == aid)
        for source in assembly_ds['inputs']:
            if summary.get('source_sha256', {}).get(str((ROOT / source['path']).resolve())) != source['sha256']:
                raise ValueError('anchor reference identity differs from registered assembly')
        for name, expected in summary['outputs_sha256'].items():
            artifacts.bind(path.parent / name, expected)
        details = artifacts.read(path.parent / 'details.json', summary['outputs_sha256']['details.json'])
        joint = [x for x in details if x['joint_terminal_context']]
        if len({(x['source_index'],x['read_ordinal']) for x in joint}) != summary['joint_terminal_context_records']:
            raise ValueError('anchor joint source-record count differs from summary')
        boundary_groups = defaultdict(list)
        for array in joint:
            for anchor in array['anchor_alignments']:
                if anchor.get('read_end_and_assembly_end_concordant'):
                    key = (anchor['target'],anchor['outward_target_end'],anchor['mapped_junction_boundary0'],
                           anchor['anchor_boundary_to_assembly_end_bp'])
                    boundary_groups[key].append({'array_id':array['array_id'],'source_index':array['source_index'],
                        'read_ordinal':array['read_ordinal'],'read_id':array['read_id'],'model_span_bp':array['span_bp'],
                        'read_array_strand':array['strand'],'anchor_mapping_strand':anchor['mapping_strand'],
                        'assembly_repeat_strand':anchor['assembly_repeat_strand'],
                        'repeat_outer_read_gap_bp':anchor['repeat_outer_read_gap_bp']})
        boundaries = [{'target':k[0],'outward_target_end':k[1],'mapped_junction_boundary0':k[2],
            'anchor_boundary_to_assembly_end_bp':k[3],
            'source_record_count':len({(x['source_index'],x['read_ordinal']) for x in values}),
            'model_span_min_bp':min(x['model_span_bp'] for x in values),
            'model_span_max_bp':max(x['model_span_bp'] for x in values), 'observations':values}
            for k,values in sorted(boundary_groups.items())]
        rows.append({**common, 'status': summary['mapping_status'],
            **{key: summary[key] for key in ['arrays', 'exported_anchors', 'array_supporting_records',
               'joint_terminal_context_records', 'joint_long_model_span_context_records',
               'joint_nuclear_terminal_context_records', 'localization_status_counts', 'settings']},
            'mapping_wall_seconds': (runtime.get('mapping') or {}).get('wall_seconds'),
            'adapter_wall_seconds': (runtime.get('adapter') or {}).get('wall_seconds'),
            'joint_boundary_groups':boundaries,
            'counts_are_source_records_not_ends_or_independent_molecules':True,
            'details_path':str((path.parent / 'details.json').relative_to(ROOT)),
            'summary_path': str(path.relative_to(ROOT))})
    return rows


def source_audit_rows(protocol, base, artifacts, require):
    rows = []
    for ds in protocol['datasets']:
        path = base / 'source_audits' / ds['dataset_id'] / 'audit.json'
        if not path.exists():
            if require:
                raise ValueError('required independent source audit not available: ' + ds['dataset_id'])
            rows.append({'dataset_id':ds['dataset_id'], 'status':'not_available', 'all_exported_traces_replayed':None})
            continue
        audit = artifacts.read(path)
        passed = (audit.get('status') == 'complete' and audit.get('complete_input') is True and
                  audit.get('validated_eof') is True and audit.get('all_exported_traces_replayed') is True)
        if require and not passed:
            raise ValueError('required independent source audit incomplete: ' + ds['dataset_id'])
        if passed:
            export = base / 'analysis' / ds['dataset_id'] / 'output/summary.json'
            counts = audit.get('counts', {})
            if (audit.get('export_summary', {}).get('sha256') != digest(export) or
                    counts.get('records') != ds['expected_records'] or counts.get('bases') != ds['expected_bases'] or
                    counts.get('arrays_not_replayed') != 0):
                raise ValueError('source audit does not bind the current complete export')
        rows.append({'dataset_id':ds['dataset_id'], 'status':audit.get('status'),
                     'all_exported_traces_replayed':audit.get('all_exported_traces_replayed'),
                     'counts':audit.get('counts'), 'audit_path':str(path.relative_to(ROOT)),
                     'complete_source_replay':passed})
    summary = base / 'source_audits/summary.json'
    if summary.exists():
        artifacts.bind(summary)
    return rows


def assembly_coordinate_summary(protocol, locations, ends, previous):
    rows = []
    for ds in protocol['datasets']:
        if ds['modality'] != 'assembly':
            continue
        arrays = [x for x in locations if x['dataset_id'] == ds['dataset_id']]
        endpoint = [x for x in ends if x['dataset_id'] == ds['dataset_id']]
        resolved = [x for x in endpoint if x['literature_endpoint_resolution'] == 'literature_resolved_end']
        old = [x for x in previous if x['dataset_id'] == ds['dataset_id'] and x['status'] == 'measured_exact_mixed']
        longest = max(arrays, key=lambda x:x['span_bp'], default=None)
        rows.append({'dataset_id':ds['dataset_id'], 'role':ds['role'], 'all_nuclear_ends':len(endpoint),
            'oriented_model_ends_within50':sum(x['expected_oriented_within50'] for x in endpoint),
            'oriented_model_ends_within5000':sum(x['expected_oriented_within5000'] for x in endpoint),
            'literature_resolved_end_denominator':len(resolved) if resolved else None,
            'literature_resolved_oriented_model_ends_within50':sum(x['expected_oriented_within50'] for x in resolved) if resolved else None,
            'oriented_long_model_ends_within50':sum(x['long_model_expected_arrays_within50'] > 0 for x in endpoint),
            'oriented_long_model_ends_within5000':sum(x['long_model_expected_arrays_within5000'] > 0 for x in endpoint),
            'coordinate_inset_arrays_beyond5000':sum(x['coordinate_inset_beyond_5000bp'] for x in arrays),
            'coordinate_inset_long_arrays_beyond5000':sum(x['coordinate_inset_beyond_5000bp'] and x['span_bp'] >= 500 for x in arrays),
            'arrays_crossing_either5000_boundary':sum(x['crosses_left_5000_boundary'] or x['crosses_right_5000_boundary'] for x in arrays),
            'longest_model_array':longest,
            'prior_exact_mixed_end_aggregate':[{ 'group':x['unit_or_group'],
                'left_c_end_records':x['exact_mixed_group_counts']['left_c_end_records'],
                'right_g_end_records':x['exact_mixed_group_counts']['right_g_end_records'],
                'max_exact_mixed_span_bp':x['max_span_bp']} for x in old],
            'measurement_equivalence':False,'physical_telomere_architecture_confirmed':False})
    return rows


def display(value):
    if value is None: return '미측정'
    if isinstance(value, float): return f'{value:,.3f}'
    if isinstance(value, int): return f'{value:,}'
    return str(value)


def markdown_table(rows, fields):
    lines = ['| ' + ' | '.join(label for _, label in fields) + ' |',
             '| ' + ' | '.join('---' for _ in fields) + ' |']
    lines += ['| ' + ' | '.join(('미보고' if row.get(key) is None and row.get('status') == 'not_reported_not_zero'
                                else display(row.get(key))).replace('|', '\\|') for key, _ in fields) + ' |' for row in rows]
    return '\n'.join(lines)


def render(report, output):
    models, prior, ends, windows, anchors = (report[k] for k in ['datasets', 'previous_measurements', 'assembly_ends', 'internal_windows', 'hifi_anchor_context'])
    measured = [x for x in windows if x['status'] == 'measured_original_coordinates']
    summary_ends = report['assembly_coordinate_summary']
    fields = [('dataset_id', '자료'), ('positive_records', '모델 관측 records'), ('arrays', '모델 구간'),
              ('max_model_span_bp', '최장 model span(bp)'), ('long_modeled_records', '≥500 bp 모델 records'),
              ('full_modeled_records', '≥90% 모델 records'), ('invocation_wall_seconds', '검증 포함 초'), ('max_rss_mib', 'RSS MiB')]
    text = '# 선택적 불규칙 배열 모델: 기존 결과와 별도 비교\n\n'
    if report['status'] != 'complete':
        text += '**중간 스냅샷입니다.** 아래 source audit 상태는 생성 시점의 상태입니다. 진행 중인 audit 파일은 이후 바뀔 수 있으므로 이 중간 보고서의 provenance 해시가 최신 파일과 계속 일치한다고 주장하지 않습니다. 최종 검증 보고서는 별도 comparison/에 생성합니다.\n\n'
    text += '등록한 6개 전체 입력의 결과를 집계했습니다. 기존 native 판정은 그대로이며, 여기의 길이는 **제공한 단위와 허용 편집으로 설명한 model span**입니다. 정확 반복 길이 또는 생물학적 telomere 길이로 바꾸어 읽으면 안 됩니다.\n\n'
    text += 'S. japonicus 세 자료는 개발 자료, S. osmophilus 두 자료와 Col-0은 이미 사용한 회귀 자료입니다. 새 독립 holdout은 없습니다. 동일 BioSample의 reads와 derived assembly는 독립 생물학적 반복이 아닙니다.\n\n'
    text += '이 선택적 helper는 short WGS의 R1/R2를 각각 source record로 처리하며 paired-end 방향을 다시 판정하지 않습니다. 기존 native의 paired-end 근거와 분류는 보존됩니다.\n\n'
    for row in summary_ends:
        text += f'**{row["dataset_id"]}**: 원래 핵 record 끝 {row["all_nuclear_ends"]}개 중 방향 일치 모델이 50 bp 이내인 끝은 {row["oriented_model_ends_within50"]}개, 5 kb 이내인 끝은 {row["oriented_model_ends_within5000"]}개입니다. 양 끝에서 모두 5 kb 넘게 떨어진 ≥500 bp 모델 배열은 {row["coordinate_inset_long_arrays_beyond5000"]}개입니다. 긴 모델 구간의 존재만으로 말단 구조가 복원됐다고 판단하지 않습니다.\n\n'
        if row['literature_resolved_end_denominator']:
            text += f'문헌에서 해소된 끝만 별도로 보면 ≤50 bp 방향 일치는 {row["literature_resolved_oriented_model_ends_within50"]}/{row["literature_resolved_end_denominator"]}입니다. 알려진 미해결 chr3-left rDNA 쪽 record 끝은 위의 전체 6개 좌표 분모에 보존하되 생물학적 음성으로 세지 않습니다.\n\n'
    text += markdown_table(models, fields) + '\n\n'
    text += '검증 포함 시간은 실제 invocation의 pre/post 입력 해시·읽기·scan·출력을 포함합니다. 내부 scan_and_emit_seconds는 flush/sync 전 scan 및 출력 작성 시간이며 별도 TSV에 있습니다. 이전 실행은 버전·설정·threads·메모리·검증 범위가 달라 직접적인 speedup을 계산하지 않습니다.\n\n'
    text += '## 이전 측정과의 차이\n\n'
    text += '기존 exact-mixed는 ≥24 bp, ≥2종 단위마다 ≥2개를 보장하는 정확 문법이며 record 경계에서 문법과 일치하는 partial unit을 허용합니다. 새 모델은 ≥100 bp, 방향별 literal 단위 4개를 gate로 사용하고 각 단위에 허용된 편집을 적용하며 완성된 edited unit만 사용합니다. Record 경계의 partial-unit 확장은 하지 않습니다. **새 모델은 한 component만으로도 통과할 수 있으므로 mixed 검출 증가나 민감도 개선으로 해석하지 않습니다.** Native homogeneous extension은 identity≥0.85로 치환을 허용하므로 exact-only 결과도 아닙니다.\n\n'
    text += markdown_table(prior, [('dataset_id','자료'),('method','이전 측정'),('unit_or_group','단위/그룹'),('status','보고 위치'),('positive_records','관측 records'),('max_span_bp','이전 측정 최장 bp'),('long_records','이전 ≥500 records')]) + '\n\n'
    text += '`rejected_candidates`는 global catalogue 배치이며 해당 motif의 생물학적 허위·novelty를 뜻하지 않습니다. 이전 후보/제외 목록에 모두 없으면 스캔했으나 보고 기준 미충족인 **미보고**로 남기며 숫자 0으로 바꾸지 않습니다. 수행하지 않은 분석은 미측정입니다. 서로 다른 단위의 records를 더하지 않습니다. [기존 측정 전체](previous_measurements.tsv)\n\n'
    text += '## 원래 assembly 말단 좌표\n\n'
    text += markdown_table(summary_ends, [('dataset_id','Assembly'),('oriented_model_ends_within50','방향 일치 ≤50 bp 끝'),('oriented_model_ends_within5000','방향 일치 ≤5 kb 끝'),('oriented_long_model_ends_within50','≥500 bp 모델·방향 ≤50 bp 끝'),('all_nuclear_ends','모든 핵 record 끝')]) + '\n\n'
    text += '왼쪽 C(-), 오른쪽 G(+) 방향으로 판정했습니다. 모든 원래 핵 record의 양 끝을 유지하며 미관측 끝도 분모에서 빼지 않습니다. Ambiguous 방향은 일치로 세지 않습니다. 5 kb 관측은 50 bp 관측과 별도이며 둘 다 물리적 말단 확인은 아닙니다. S. japonicus chr3 왼쪽과 S. osmophilus rDNA 관련 불완전 말단의 한계도 남습니다. [끝별 결과](assembly_ends.tsv) · [전체 원래 배열 좌표](assembly_arrays.tsv)\n\n'
    text += f'## 내부 좌표 진단\n\n사전 제안된 78개 창 중 현재 등록 assembly에 속하는 **{len(measured)}개**만 전체 assembly 출력과 겹침을 계산했습니다. 나머지 **{78-len(measured)}개는 미측정**이며 0으로 채우지 않았습니다. 측정 창에서 모델 구간과 겹친 창 {sum(x["overlapping_arrays"]>0 for x in measured)}개, union-covered bp 합 {sum(x["union_modeled_bp"] for x in measured):,} bp입니다.\n\n'
    text += '이는 염색체 내부 좌표의 모델 부담 진단이며 biological false-positive rate가 아닙니다. 원래 좌표를 사용하고 crop 경계를 말단으로 취급하지 않았습니다. 창 경계를 넘는 배열은 별도 열로 표시하며 clipped overlap 길이를 완전한 배열 길이로 해석하지 않습니다. [78개 전체 목록](internal_windows.tsv)\n\n'
    text += '## HiFi의 같은 read anchor 근거\n\n'
    text += markdown_table(anchors, [('dataset_id','HiFi'),('status','상태'),('joint_terminal_context_records','말단·방향 동시 맥락 records'),('joint_long_model_span_context_records','같은 ≥500 bp 모델 동시 맥락 records')]) + '\n\n'
    text += '**위 숫자는 서로 다른 source records 수이며, chromosome ends 수 또는 독립 molecules 수가 아닙니다.** 이 anchor 판정의 assembly 거리는 5 kb 기준이고, assembly 배열 자체의 50 bp 기준과 다릅니다.\n\n'
    for row in anchors:
        for boundary in row.get('joint_boundary_groups', []):
            detail_link = os.path.relpath(ROOT / row['details_path'], output)
            text += f'- {row["dataset_id"]}: {boundary["source_record_count"]}개 source records의 anchor 반복 쪽 정렬 경계가 {boundary["target"]}의 {boundary["outward_target_end"]} 방향, 원래 끝에서 {boundary["anchor_boundary_to_assembly_end_bp"]:,} bp에 위치합니다. 해당 model span은 {boundary["model_span_min_bp"]}–{boundary["model_span_max_bp"]} bp입니다. [read·strand·정렬 상세]({detail_link})\n'
    text += '\n'
    text += 'Read 끝·assembly 끝·방향이 같은 배열에서 함께 맞는지 별도로 평가합니다. 서로 떨어진 배열 또는 다른 reads의 길이를 합하지 않습니다. 양쪽 anchor가 정렬되는 경우와 모호한 방향은 별도 상태입니다. Anchor의 유일한 비반복 subtelomere 성격, 완전한 생물학적 telomere 길이 및 독립 복제는 아직 확인하지 않습니다. [정렬 맥락 요약](hifi_anchor_context.tsv)\n\n'
    text += '## 남는 한계와 해석\n\n'
    text += 'Model identity, 정확 단위 bp 합, substitutions/insertions/deletions, parse ambiguity, 최대 연속 edited-unit 수 및 edited-only span을 [현재 모델 결과](dataset_summary.tsv)에 분리했습니다. 정확 단위 bp 합은 연속 정확 반복 길이가 아니며, 겹치는 모델 배열의 편집 합계는 중복될 수 있습니다. Interval-union bp는 같은 원래 record 안의 겹침을 제거한 값입니다.\n\n'
    text += markdown_table(models,[('dataset_id','자료'),('identity_min','최소 model identity'),('identity_from_summed_model_operations','작업 수 합산 identity'),('exact_unit_bases_sum','정확 단위 bp 합'),('parse_ambiguous_arrays','분해 모호 배열'),('orientation_ambiguous_arrays','방향 모호 배열'),('max_edited_only_model_span_bp','최장 edited-only span bp')]) + '\n\n'
    text += '전역 literal 4개는 배열 내부 모든 구간에서 높은 정확 지지를 보장하지 않습니다. Gate로 미계산한 방향은 음성 판정이 아닙니다. 데이터 충실도·민감도·부재·ALT·신규성·karyotype를 새로 판정하지 않으며, short reads를 특정 내부 또는 말단 locus에 배정하지 않습니다. 다음 단계는 별도 실제 자료와 기준을 고정한 검증이며 이번 결과에 맞추어 thresholds를 조정하지 않았습니다.\n\n'
    text += '**과학적 권고:** 제공한 두 단위와 단위당 한 편집 모델은 S. japonicus 말단 구조를 충분히 설명하지 못합니다. 이 단계는 선택적 구조 주석으로 유지하고, 실제 말단 서열에서 누락된 단위·위상·국소 지지의 한계를 먼저 정리해야 합니다. 이후 수정할 모델은 별도 설정으로 동결한 뒤 독립 실제 자료에서 검증하는 것이 우선입니다. 현재 모델 미관측을 calibrated novel/absence 판정이나 새로운 말단 기작의 근거로 사용하지 않습니다.\n\n'
    text += '독립 source replay 상태는 [검증 목록](source_audits.tsv)에 별도로 연결했습니다. 집계 성공은 별도 sequence replay 또는 생물학적 검증을 대신하지 않습니다.\n\n[전체 근거 JSON](comparison.json) · [파일 해시 및 집계 receipt](summary.json)\n'
    (output / 'README.md').write_text(text)
    links = ' · '.join('<a href="'+name+'">'+html.escape(name)+'</a>' for name in
                       ['README.md','dataset_summary.tsv','previous_measurements.tsv','assembly_ends.tsv','assembly_arrays.tsv','internal_windows.tsv','hifi_anchor_context.tsv','source_audits.tsv','comparison.json','summary.json'])
    def html_table(rows, fields):
        header = '<tr>'+''.join('<th>'+html.escape(label)+'</th>' for _,label in fields)+'</tr>'
        body = ''.join('<tr>'+''.join('<td>'+html.escape(display(row.get(key)))+'</td>' for key,_ in fields)+'</tr>' for row in rows)
        return '<div class="scroll"><table>'+header+body+'</table></div>'
    body = '<h1>선택적 불규칙 배열 모델 비교</h1><p><strong>Model span은 정확 반복 길이가 아닙니다.</strong> 동일 측정값의 개선·민감도·speedup을 계산하지 않았으며 기존 TRM 판정은 변경하지 않았습니다.</p>'
    if report['status'] != 'complete':
        body += '<p><strong>중간 스냅샷:</strong> 진행 중인 audit 파일의 해시는 생성 시점 관측이며 이후 최신 파일과의 일치를 보장하지 않습니다. 최종 보고서는 별도 comparison/입니다.</p>'
    body += '<p>6개 전체 자료는 개발 또는 기존 회귀 자료입니다. 새로운 독립 holdout은 없습니다.</p><h2>전체 입력 결과</h2>'+html_table(models,fields)
    body += '<h2>원래 assembly 좌표의 말단 근거</h2>'+html_table(summary_ends,[('dataset_id','Assembly'),('oriented_model_ends_within50','방향 일치 ≤50 bp'),('oriented_model_ends_within5000','방향 일치 ≤5 kb'),('all_nuclear_ends','모든 핵 record 끝'),('coordinate_inset_long_arrays_beyond5000','좌표상 내부 ≥500 bp 배열')])
    body += '<p>긴 배열의 최대값과 물리적 말단의 확인은 다릅니다. 상세 위치·편집 근거 및 미관측 끝을 그대로 보존했습니다.</p>'
    body += '<h2>HiFi 같은 read의 anchor 맥락</h2>'+html_table(anchors,[('dataset_id','HiFi'),('status','상태'),('joint_terminal_context_records','말단·방향 records'),('joint_long_model_span_context_records','같은 ≥500 bp 모델 records')])
    body += '<p>수치는 source records 수이며 chromosome ends 또는 독립 molecules 수가 아닙니다. Anchor 맥락은 5 kb 기준으로, assembly의 50 bp 기준과 별도입니다.</p>'
    boundary_rows = [{'dataset_id':row['dataset_id'],**b} for row in anchors for b in row.get('joint_boundary_groups',[])]
    body += html_table(boundary_rows,[('dataset_id','자료'),('target','Reference record'),('outward_target_end','방향'),('anchor_boundary_to_assembly_end_bp','원래 끝까지 bp'),('source_record_count','source records'),('model_span_min_bp','최소 model span'),('model_span_max_bp','최대 model span')])
    body += f'<h2>내부 창 진단</h2><p>78개 제안 창 중 {len(measured)}개 측정, {78-len(measured)}개 미측정입니다. Union coverage는 원래 전체 assembly 좌표와 겹친 구간이며 biological false-positive rate가 아닙니다. Crop 경계는 말단이 아닙니다.</p>'
    body += '<p>기존 exact-mixed/native 측정, 모든 단위별 비교와 한계는 <a href="README.md">전체 설명</a> 및 아래 파일에 있습니다. 부재·ALT·신규성·데이터 충분성 판정은 하지 않습니다.</p><nav>'+links+'</nav>'
    (output / 'index.html').write_text('<!doctype html><html lang="ko"><meta charset="utf-8"><title>Irregular pilot comparison</title><style>body{max-width:1300px;margin:2em auto;padding:0 1em;font:16px/1.6 sans-serif;color:#172b3a}table{border-collapse:collapse;width:100%;font-size:14px}td,th{border:1px solid #ccd7df;padding:.4em;text-align:left}th{background:#eaf1f6}.scroll{overflow:auto}nav{line-height:2}h2{margin-top:2em}</style>'+body+'</html>')


def self_test():
    assert union_bp([]) == 0
    assert union_bp([(10, 20), (15, 30), (30, 35), (50, 60)]) == 35
    assert overlap(40, 70, 50, 60) == 10
    assert overlap(0, 50, 50, 60) == 0
    ds = {'dataset_id': 'numeric', 'role': 'numeric_test', 'nuclear_records': [{'id':'coordinate', 'length':10000}]}
    arrays = [{'array_id':'a','read_id':'coordinate','read_length':10000,'start0':50,'end0':550,'span_bp':500,'strand':'-', 'model':{}},
              {'array_id':'b','read_id':'coordinate','read_length':10000,'start0':9400,'end0':9949,'span_bp':549,'strand':'+','model':{}}]
    ends, _ = assembly_end_rows(ds, arrays)
    assert ends[0]['expected_oriented_within50'] is True
    assert ends[1]['expected_oriented_within50'] is False
    assert ends[1]['expected_oriented_within5000'] is True
    assert longest_edited_span([{'edit':'S','start0':0,'end0':6},{'edit':'I','start0':6,'end0':13},{'edit':'=','start0':13,'end0':19}]) == 13
    print(json.dumps({'status':'pass','test_kind':'numeric_coordinates_only','source_sequences_read':False}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=DEFAULT_BASE)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--require-anchors', action='store_true')
    parser.add_argument('--require-source-audits', action='store_true')
    parser.add_argument('--check', action='store_true', help='metadata/registered-comparator check only; no model outputs read')
    parser.add_argument('--self-test', action='store_true', help='numeric coordinate checks; no sequence files')
    args = parser.parse_args()
    if args.self_test:
        self_test(); return
    started = time.monotonic()
    artifacts = Artifacts()
    protocol, windows = check_plan(args.base, artifacts)
    if args.check:
        print(json.dumps({'status':'prepared','datasets':6,'proposed_windows':len(windows),'source_sequences_read':False})); return
    output = args.output_dir or args.base / 'comparison'
    if output.exists():
        raise ValueError('new report output directory required; preserve previous reports')
    state = artifacts.read(args.base / 'run_summary.json')
    if state.get('status') != 'complete' or len(state.get('datasets', [])) != 6:
        raise ValueError('all six primary full-input runs must complete before reporting')
    if state.get('protocol_sha256') != digest(args.base / 'protocol.json'):
        raise ValueError('primary run does not bind current protocol')
    models, previous, ends, locations, assemblies = [], [], [], [], {}
    for ds in protocol['datasets']:
        row, arrays = model_rows(ds, args.base, artifacts)
        models.append(row)
        previous.extend(prior_rows(ds, artifacts))
        if ds['modality'] == 'assembly':
            assemblies[ds['dataset_id']] = arrays
            endpoint, positioned = assembly_end_rows(ds, arrays)
            ends.extend(endpoint); locations.extend(positioned)
    window_results = window_rows(windows, protocol, assemblies)
    anchors = anchor_rows(protocol, args.base, artifacts, args.require_anchors)
    audits = source_audit_rows(protocol, args.base, artifacts, args.require_source_audits)
    coordinate_summary = assembly_coordinate_summary(protocol, locations, ends, previous)
    final_audits = all(x.get('complete_source_replay') is True for x in audits)
    final_anchors = all(x['status'] in {'complete_verified_alignment','not_run_no_exported_anchors'} for x in anchors)
    report_status = 'complete' if final_audits and final_anchors else 'provisional_pending_validation'
    report = {'schema_version':1, 'program':'summarize_irregular_pilot', 'status':report_status,
        'created_utc':datetime.now(timezone.utc).isoformat(), 'protocol_version':protocol['protocol_version'],
        'datasets':models,'previous_measurements':previous,'assembly_ends':ends,'assembly_arrays':locations,
        'internal_windows':window_results,'hifi_anchor_context':anchors,'inference_limits':protocol['inference_limits'],
        'assembly_coordinate_summary':coordinate_summary,'source_audits':audits,
        'measurement_equivalence':False,'sensitivity_or_speedup_calculated':False,
        'biological_false_positive_rate_calculated':False,'source_sequences_read':False,
        'independent_source_audit_status':'linked per-dataset receipts; aggregation is not an independent source replay',
        'changes_native_classification':False,'changes_catalogue':False,'changes_adequacy':False}
    output.mkdir(parents=True)
    write_json(output / 'comparison.json', report)
    for name, rows in [('dataset_summary',models),('previous_measurements',previous),('assembly_ends',ends),
                       ('assembly_arrays',locations),('internal_windows',window_results),('hifi_anchor_context',anchors)]:
        write_table(output / (name + '.tsv'), rows)
    write_table(output / 'assembly_coordinate_summary.tsv', coordinate_summary)
    write_table(output / 'source_audits.tsv', audits)
    render(report, output)
    artifacts.bind(__file__)
    receipt = {'schema_version':1,'status':report_status,'program':'summarize_irregular_pilot',
        'created_utc':report['created_utc'],'wall_seconds':time.monotonic()-started,
        'source_artifacts':artifacts.sources,'source_sequences_read':False,
        'source_artifacts_are_observation_time_snapshots':report_status != 'complete',
        'measured_internal_windows':sum(x['status']=='measured_original_coordinates' for x in window_results),
        'unmeasured_internal_windows':sum(x['status']!='measured_original_coordinates' for x in window_results),
        'outputs':{p.name:{'sha256':digest(p),'size_bytes':p.stat().st_size} for p in sorted(output.iterdir()) if p.is_file()}}
    write_json(output / 'summary.json', receipt)
    print(json.dumps({'status':report_status,'output_dir':str(output),'datasets':6,'wall_seconds':receipt['wall_seconds']}))


if __name__ == '__main__':
    main()
