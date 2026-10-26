#!/usr/bin/env python3
"""Describe actual deposited assembly windows under a frozen unit vocabulary.

This is an explanatory audit, not a new detector, fitted grammar, or biological
negative classifier. Full-record detector outputs remain authoritative. All
sequence strings in its outputs are unchanged source slices or their declared
reverse complements; no constructed sequence control is used.
"""
import argparse
from collections import Counter, defaultdict
import contextlib
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time

import audit_irregular_arrays as audit

COMPLEMENT = bytes.maketrans(b"ACGTRYSWKMBDHVN", b"TGCAYRSWMKVHDBN")


def union(intervals):
    result = []
    for start, end in sorted(intervals):
        if not 0 <= start < end:
            raise ValueError("invalid interval")
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def original_coordinates(start, end, window_start, window_end, reverse):
    if not 0 <= start < end <= window_end - window_start:
        raise ValueError("interval outside window")
    return (window_end - end, window_end - start) if reverse else (window_start + start, window_start + end)


def paths(n, edges):
    """Full endpoint table for an integer DAG; edges are (start,end,cost).

    Keep earliest start, then minimum cost, independently at every endpoint.
    Return maximal intervals as well as all endpoint states for diagnosis.
    No path-history feasibility constraint or free gap edge is added.
    """
    by_end = defaultdict(list)
    for start, end, cost in edges:
        if not 0 <= start < end <= n or cost < 0:
            raise ValueError("invalid numeric edge")
        by_end[end].append((start, cost))
    states = [None] * (n + 1)
    for end in range(1, n + 1):
        choices = []
        for start, cost in by_end[end]:
            previous = states[start]
            choices.append((previous[0] if previous else start,
                            cost + (previous[1] if previous else 0)))
        if choices:
            states[end] = min(choices)
    endpoints = [(state[0], end, state[1]) for end, state in enumerate(states) if state]
    maximal, furthest = [], -1
    for start, end, cost in sorted(endpoints, key=lambda x: (x[0], -x[1], x[2])):
        if end > furthest:
            maximal.append((start, end, cost))
            furthest = end
    return maximal, endpoints


def tokenize(sequence, units):
    lengths = sorted({len(u) + d for u in units for d in ((-1, 0, 1) if len(u) >= 5 else (0,))})
    edges, literals, cache = [], [], {}
    for end in range(1, len(sequence) + 1):
        for length in lengths:
            start = end - length
            if start < 0:
                continue
            observed = sequence[start:end]
            if observed not in cache:
                cache[observed] = audit.explanations(observed, units) if audit.ACGT.fullmatch(observed) else []
            explanations = cache[observed]
            if explanations:
                cost = min(x[0] for x in explanations)
                edges.append((start, end, cost))
                if cost == 0:
                    literals.append((start, end))
    return edges, literals


def nonoverlapping_literals(intervals):
    chosen, last = [], 0
    for start, end in sorted(set(intervals), key=lambda x: (x[1], x[0])):
        if start >= last:
            chosen.append((start, end))
            last = end
    return chosen


def describe(sequence, units, window_start, window_end, reverse, terminal, config):
    edges, literals = tokenize(sequence, units)
    maximal, endpoints = paths(len(sequence), edges)
    exact_paths, _ = paths(len(sequence), [e for e in edges if e[2] == 0])
    chosen = nonoverlapping_literals(literals)
    literal_union = union(literals)
    token_union = union((a, b) for a, b, _ in edges)
    gaps = [(a[1], b[0]) for a, b in zip(chosen, chosen[1:]) if b[0] > a[1]]
    gap_strings = Counter(sequence[a:b].decode() for a, b in gaps)

    def row(path):
        start, end, cost = path
        a, b = original_coordinates(start, end, window_start, window_end, reverse)
        hits = audit.literal_hits(sequence[start:end], units)
        return {'oriented_start0': start, 'oriented_end0': end, 'original_start0': a,
                'original_end0': b, 'span_bp': end-start, 'minimum_edit_count': cost,
                'literal_gate_hits': hits, 'passes_min_span': end-start >= config['min_span'],
                'passes_literal_gate': hits >= config['min_exact_anchors'],
                'distance_to_deposited_end': len(sequence)-end if terminal else None,
                'near_artificial_inner_edge': start < max(map(len, units))+1,
                'local_window_diagnostic_only': True}

    path_rows = [row(p) for p in maximal]
    near = [p for p in path_rows if terminal and p['distance_to_deposited_end'] <= 50]
    qualifying = [p for p in path_rows if p['passes_min_span'] and p['passes_literal_gate']]
    qualifying_near = [p for p in near if p['passes_min_span'] and p['passes_literal_gate']]
    literal_bp = sum(b-a for a, b in literal_union)
    token_bp = sum(b-a for a, b in token_union)
    covered = bytearray(len(sequence))
    for a, b in token_union:
        covered[a:b] = b'\x01' * (b-a)
    uncovered, start = [], None
    for pos, flag in enumerate([*covered, 1]):
        if not flag and start is None:
            start = pos
        elif flag and start is not None:
            a, b = original_coordinates(start, pos, window_start, window_end, reverse)
            uncovered.append({'oriented_start0':start, 'oriented_end0':pos,
                              'original_start0':a, 'original_end0':b,
                              'source_string_outward_orientation':sequence[start:pos].decode()})
            start = None
    return {'literal_occurrences':len(literals), 'nonoverlapping_literal_hits':len(chosen),
            'literal_union_bp':literal_bp, 'literal_union_fraction':literal_bp/len(sequence),
            'allowed_token_union_bp':token_bp, 'allowed_token_union_fraction':token_bp/len(sequence),
            'positions_not_in_any_allowed_token':len(sequence)-token_bp,
            'max_exact_complete_unit_path_bp':max((b-a for a,b,_ in exact_paths),default=0),
            'max_edit_model_path_bp_before_filters':max((p['span_bp'] for p in path_rows),default=0),
            'max_edit_model_path_near_end_bp_before_filters':max((p['span_bp'] for p in near),default=0) if terminal else None,
            'max_near_end_path_literal_hits_before_filters':max((p['literal_gate_hits'] for p in near),default=0) if terminal else None,
            'qualifying_local_paths':len(qualifying), 'qualifying_local_paths_near_end':len(qualifying_near) if terminal else None,
            'nearest_qualifying_local_path_distance_bp':min((p['distance_to_deposited_end'] for p in qualifying), default=None) if terminal else None,
            'near_end_paths_passing_span_only':sum(p['passes_min_span'] for p in near) if terminal else None,
            'near_end_paths_passing_gate_only':sum(p['passes_literal_gate'] for p in near) if terminal else None,
            'literal_gap_count':len(gaps),
            'literal_gap_length_histogram':dict(sorted(Counter(b-a for a,b in gaps).items())),
            'observed_literal_gap_strings':[{'source_string':s,'count':c,'length':len(s)} for s,c in gap_strings.most_common()],
            'literal_intervals_oriented':literals, 'chosen_nonoverlapping_literals_oriented':chosen,
            'allowed_token_union_intervals_oriented':token_union,
            'uncovered_source_runs':uncovered, 'maximal_local_paths':path_rows,
            'endpoint_states_count':len(endpoints), 'allowed_token_edges':len(edges),
            'sequence_outward_orientation':sequence.decode()}


def write_json(path, value):
    with path.open('x') as f:
        json.dump(value, f, indent=2, ensure_ascii=False, allow_nan=False)
        f.write('\n')


def write_tsv(path, rows):
    keys = list(rows[0]) if rows else []
    with path.open('x') as f:
        writer = csv.DictWriter(f, keys, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def run(args):
    started = time.monotonic()
    if args.output_dir.exists():
        raise ValueError('new output directory required')
    summary_path = args.model_export/'summary.json'
    arrays_path = args.model_export/'arrays.jsonl'
    frozen_inputs = {str(p):audit.digest(p) for p in (summary_path, arrays_path, args.interior_windows)}
    summary = json.loads(summary_path.read_text())
    if summary['status'] != 'complete' or summary['complete_input'] is not True or summary['validated_eof'] is not True:
        raise ValueError('complete full-record model export required')
    if summary['sampled_prefix'] or summary['resource_limit_hit'] or summary['export_truncated']:
        raise ValueError('censored model export')
    if len(summary['inputs']) != 1 or Path(summary['inputs'][0]['path']).resolve() != args.assembly.resolve():
        raise ValueError('full-record model source mismatch')
    if len(summary['definitions']) != 1:
        raise ValueError('this diagnostic accepts one supplied group at a time')
    definitions = audit.definitions_from_summary(summary)
    units = next(iter(definitions.values()))['+']
    config = summary['config']
    expected = summary['inputs'][0]
    source_stat = audit.snapshot(args.assembly)
    if audit.digest(args.assembly) != expected['sha256'] or source_stat['size_bytes'] != expected['size_bytes']:
        raise ValueError('source does not match frozen model export')
    if frozen_inputs[str(arrays_path)] != summary['outputs']['arrays.jsonl']['sha256']:
        raise ValueError('frozen arrays hash mismatch')
    arrays = [json.loads(line) for line in arrays_path.open()]
    if len(arrays) != summary['arrays']:
        raise ValueError('array count mismatch')
    with contextlib.ExitStack() as stack:
        stream = audit.VerifiedStream(stack, args.assembly, expected)
        sequences = dict(audit.records(stream.stream))
        source_receipt = stream.finish()
    if len(sequences) != expected['records'] or sum(map(len,sequences.values())) != expected['bases']:
        raise ValueError('source record/base count mismatch')
    widths = [500, 1000, 5000]
    planned = []
    for name, sequence in sequences.items():
        if len(sequence) < max(widths)*2:
            raise ValueError('nonoverlapping complete fixed terminal windows require records >=10000bp')
        for side in ['left','right']:
            for width in widths:
                start = 0 if side == 'left' else len(sequence)-width
                planned.append({'record_id':name,'record_length':len(sequence),'role':'deposited_terminal_window',
                                'side':side,'window_start0':start,'window_end0':start+width,'window_bp':width})
    with args.interior_windows.open() as f:
        for window in csv.DictReader(f, delimiter='\t'):
            if Path(window['assembly_path']).resolve() != args.assembly.resolve():
                continue
            name = window['record_id']; start,end = int(window['start']),int(window['end'])
            if name not in sequences or len(sequences[name]) != int(window['record_length']) or not 0 <= start < end <= len(sequences[name]):
                raise ValueError('invalid registered interior window')
            planned.append({'record_id':name,'record_length':len(sequences[name]),'role':'coordinate_internal_unlabelled',
                            'side':'internal','window_start0':start,'window_end0':end,'window_bp':end-start})
    details, table = [], []
    for window in planned:
        terminal = window['side'] != 'internal'
        orientations = ['-'] if window['side']=='left' else ['+'] if terminal else ['+','-']
        for strand in orientations:
            reverse = strand == '-'
            source_slice = sequences[window['record_id']][window['window_start0']:window['window_end0']]
            sequence = source_slice[::-1].translate(COMPLEMENT) if reverse else source_slice
            data = describe(sequence,units,window['window_start0'],window['window_end0'],reverse,terminal,config)
            observed = [a for a in arrays if a['read_id']==window['record_id'] and a['strand']==strand]
            nearest = min(observed,key=lambda a:a['start0'] if reverse else window['record_length']-a['end0'],default=None) if terminal else None
            distance = (nearest['start0'] if reverse else window['record_length']-nearest['end0']) if nearest else None
            accepted = bool(distance is not None and distance <= 50)
            if not terminal:
                reason = 'internal_coordinate_diagnostic_not_a_negative_control'
            elif accepted:
                reason = 'frozen_full_record_model_meets_distance_and_orientation'
            elif not data['near_end_paths_passing_span_only']:
                reason = 'no_local_allowed_path_ge_min_span_near_end'
            elif not data['qualifying_local_paths_near_end']:
                reason = 'local_paths_fail_literal_gate'
            else:
                reason = 'local_vs_full_record_discrepancy_requires_review'
            data.update(window, strand=strand, original_slice_sha256=hashlib.sha256(source_slice).hexdigest(),
                        oriented_slice_sha256=hashlib.sha256(sequence).hexdigest(),
                        frozen_nearest_expected_distance_bp=distance,
                        frozen_nearest_expected_span_bp=nearest['span_bp'] if nearest else None,
                        frozen_expected_within50=accepted if terminal else None,
                        explanatory_status=reason,
                        frozen_accepted_path_excluded_by_distance=bool(distance is not None and distance > 50) if terminal else None,
                        crop_edges_are_physical_endpoints=False,
                        physical_telomere_confirmed=False)
            details.append(data)
            row = {k:v for k,v in data.items() if not isinstance(v,(list,dict)) and k not in ['sequence_outward_orientation']}
            table.append(row)
    if audit.snapshot(args.assembly) != source_stat:
        raise ValueError('source changed during diagnosis')
    if any(audit.digest(Path(p)) != digest for p,digest in frozen_inputs.items()):
        raise ValueError('frozen model or registered window input changed during diagnosis')
    args.output_dir.mkdir(parents=True,exist_ok=False)
    write_json(args.output_dir/'windows.json',details)
    write_tsv(args.output_dir/'windows.tsv',table)
    with (args.output_dir/'actual_windows_outward.fasta').open('x') as f:
        for i, row in enumerate(details,1):
            f.write(f">window{i:03} {row['record_id']}:{row['window_start0']}-{row['window_end0']} strand={row['strand']} side={row['side']}\n{row['sequence_outward_orientation']}\n")
    receipt = {'status':'complete','analysis_kind':'data_known_explanatory_diagnosis_not_new_detector',
               'created_utc':datetime.now(timezone.utc).isoformat(),'wall_seconds':time.monotonic()-started,
               'assembly':source_receipt,'model_summary_sha256':frozen_inputs[str(summary_path)],
               'frozen_arrays_sha256':frozen_inputs[str(arrays_path)], 'definitions':summary['definitions'],'config':config,
               'frozen_inputs_verified_before_and_after':frozen_inputs,
               'record_orientations_considered':summary['counts']['orientations_considered'],
               'record_orientations_unscored_literal_gate':summary['counts']['orientations_unscored_literal_gate'],
               'source_records':len(sequences),'source_bases':sum(map(len,sequences.values())),
               'terminal_windows':sum(w['role']=='deposited_terminal_window' for w in planned),
               'interior_coordinate_windows':sum(w['role']=='coordinate_internal_unlabelled' for w in planned),
               'oriented_window_rows':len(details),
               'script_sha256':audit.digest(Path(__file__)), 'audit_helper_sha256':audit.digest(Path(audit.__file__)),
               'interior_windows_sha256':frozen_inputs[str(args.interior_windows)],
               'changes_detector_or_thresholds':False,'proposes_new_units':False,'biological_classification_changes':False,
               'limits':['Known-development-data diagnosis; no blind or prospective accuracy claim.',
                         'Token union and literal coverage are not a contiguous repeat path.',
                         'Strings between literal hits are observed gaps in this vocabulary, not inferred biological spacers or new repeat units.',
                         'Cropped-window paths can differ at artificial boundaries; full-record outputs remain authoritative.',
                         'No alternative phase, unit vocabulary or edit threshold is fitted.',
                         'Per-unit abundance, phase histograms and partial-boundary token recovery are not measured here; boundary proximity only flags possible clipping.',
                         'Physical terminal sequence, mutation mechanism, sensitivity, biological absence and short-read loci remain unconfirmed.'],
               'outputs':{p.name:{'size_bytes':p.stat().st_size,'sha256':audit.digest(p)} for p in sorted(args.output_dir.iterdir()) if p.is_file()}}
    write_json(args.output_dir/'summary.json',receipt)
    print(json.dumps({'status':'complete','windows':len(details),'wall_seconds':receipt['wall_seconds']}))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assembly',required=True,type=Path)
    p.add_argument('--model-export',required=True,type=Path)
    p.add_argument('--interior-windows',required=True,type=Path)
    p.add_argument('--output-dir',required=True,type=Path)
    run(p.parse_args())


if __name__=='__main__':
    try:
        main()
    except (ValueError,OSError,KeyError) as exc:
        print(str(exc),file=sys.stderr)
        raise SystemExit(1)
