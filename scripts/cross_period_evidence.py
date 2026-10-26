"""Sealed, optional saved cross-period diagnostics; never changes TRM calls.

This module validates observations and joins identities. It does not run an edit
model, aligner, FASTQ reader, biological classifier or quality probability model.
"""
from collections import defaultdict
from collections.abc import Mapping
from copy import deepcopy
import hashlib
import json
# Retain this module binding for callers instrumenting descriptor reads.
import os
from pathlib import Path
import re
import stat

from evidence_io import EvidenceIO

SCHEMA_VERSION = 1
SCOPE = 'conditional_saved_cross_period_junction_subset'
MAX_ARTIFACT_BYTES = 64 * 1024 * 1024
MAX_SOURCE_ARTIFACTS = 64
MAX_SOURCE_BYTES = 512 * 1024 * 1024
SHA = re.compile(r'^[0-9a-f]{64}$')
DNA = re.compile(r'^[ACGT]+$')
IDENTITY_FIELDS = ('scientific_name', 'tax_id', 'biosample')
SUMMARY_FIELDS = ('selected_event_count', 'selected_pair_count', 'selected_read_count',
    'selected_interval_count', 'strict_contrast_position_count', 'qualified_context_count',
    'boundary_supported_event_count', 'mate_geometry_evaluable_event_count',
    'inward_mate_event_count', 'joint_repeat_inward_event_count',
    'joint_continuous_matching_repeat_inward_event_count',
    'candidate_lower_cost_interval_count', 'comparator_lower_cost_interval_count',
    'equal_cost_interval_count')
FAMILY_FIELDS = ('family', 'declared_motif', 'role', *SUMMARY_FIELDS,
    'max_adjacent_reference_run_bp', 'max_joint_adjacent_reference_run_bp')
INTERVAL_FIELDS = ('interval_id', 'pair_query_id', 'mate', 'start0', 'end0', 'interval_bp',
    'candidate_motif', 'comparator_motif', 'candidate_cost', 'comparator_cost',
    'comparator_minus_candidate_cost', 'lower_cost_model', 'candidate_operations',
    'comparator_operations', 'candidate_ambiguous_operation_positions',
    'comparator_ambiguous_operation_positions', 'strict_contrast_position_count', 'deletion_quality')
EVENT_FIELDS = ('event_id', 'interval_id', 'family', 'role', 'declared_motif', 'pair_query_id',
    'mate', 'qualified_context_count', 'mate_geometry_evaluable_context_count',
    'inward_mate_context_count', 'joint_repeat_inward_context_count',
    'joint_continuous_matching_repeat_inward_context_count',
    'strict_positions_in_full_hsp_count', 'max_adjacent_reference_run_bp',
    'max_joint_adjacent_reference_run_bp')
CONTEXT_FIELDS = ('joint_id', 'event_id', 'interval_id', 'family', 'role', 'reference_accession',
    'mode', 'target', 'flank_hsp_id', 'full_read_hsp_id', 'full_read_retained_hit_index0',
    'full_read_strand', 'full_read_query_start0', 'full_read_query_end0',
    'full_read_reference_start0', 'full_read_reference_end0', 'boundary_read0',
    'boundary_reference0', 'selected_adjacent_reference_run_bp',
    'full_hsp_observed_interval_bases', 'strict_positions_in_full_hsp',
    'mate_geometry_evaluable', 'same_hsp_inward_witness', 'joint_repeat_inward',
    'joint_continuous_matching_repeat_inward', 'boundary_reference_continuous',
    'boundary_adjacent_acgt_matches', 'technical_boundary_warning', 'search_capped',
    'closest_mate', 'compartment', 'flank_fcs', 'full_read_fcs', 'mate_fcs', 'aligned_technical_overlap')
POSITION_FIELDS = ('interval_id', 'pair_query_id', 'mate', 'absolute_position0', 'query_base',
    'contrast', 'quality_character', 'quality_ascii', 'conditional_phred33',
    'candidate_operations', 'comparator_operations', 'candidate_template_bases',
    'comparator_template_bases')
_CONTRAST_FIELDS = set(POSITION_FIELDS) - {'interval_id', 'pair_query_id', 'mate'}
_INTERVAL_FIELDS = {'interval_id', 'pair_query_id', 'mate', 'start0', 'end0', 'full_read_length',
    'sequence', 'quality', 'candidate_cost', 'comparator_cost', 'candidate_operations',
    'comparator_operations', 'candidate_ambiguous_operation_positions',
    'comparator_ambiguous_operation_positions', 'deletion_quality', 'strict_contrasts'}
_EVENT_FIELDS = {'event_id', 'interval_id', 'family', 'role', 'declared_motif'}
_CONTEXT_FIELDS = set(CONTEXT_FIELDS) - {'interval_id', 'family', 'role',
    'selected_adjacent_reference_run_bp', 'joint_repeat_inward',
    'joint_continuous_matching_repeat_inward'} | {'exact_flank_hit', 'selected_expected_run'}
_CONFIRMATIONS = ('biological_origin_confirmed', 'physical_terminal_confirmed',
    'main_trm_confirmed', 'novel_trm_confirmed', 'library_insert_compatibility')
LIMITATIONS = [
    'Only previously selected saved intervals and junctions are assessed; no full-library frequency or independent biological replication is estimated.',
    'A lower periodic edit cost is conditioned on earlier motif-dependent tract selection; it is not an independent motif accuracy test, likelihood or biological confirmation.',
    'Operation ranges cover complete optimal paths; marginal minima and maxima are not necessarily jointly attainable.',
    'Strict contrasts require an observed base to match on every optimum of one model and be edited on every optimum of the other model. Other positions are not asserted to be nondiscriminating.',
    'Raw quality characters are preserved; Phred+33 is conditional and encoding/authenticity are unverified. Quality is not a calibrated variant probability.',
    'Model deletions have no observed query base and no quality value; mandatory deletion counts cannot supply base-quality support.',
    'Joint evidence requires one retained full-read alignment and its own flank context; independent maxima or mate placements from another HSP are not combined.',
    'Adjacent reference run length and the number of observed interval bases covered by that HSP are separate quantities.',
    'Geometry is null when unassessed. A zero joint-positive count means no observed positive association, not an evaluated biological absence.',
    'Reference and method copies are deduplicated by event, read, interval and pair in summaries; mates, PCR copies and accessions are not independent biological replicates.',
    'Strict contrast counts use interval-plus-query-position identities. Overlapping intervals may describe the same physical base more than once; counts are not independent base observations.',
    'Capped or heuristic alignments, different reference specimens, compartment assignments and FCS overlaps limit interpretation; any-query or alternate-mate overlap is not junction contamination.',
    'Reference ends, inward mates and exact repeat length do not establish physical termini, insert compatibility, main TRM, novel TRM or read origin. No universal 500 bp requirement is applied.',
    'Payload and source hashes establish integrity relative to supplied artifacts, not authorship or biological truth; native counts, classification and adequacy are unchanged.',
]


class EvidenceError(ValueError):
    """Explicit evidence is invalid, stale or bound to another result."""


_IO = EvidenceIO(EvidenceError, 'cross-period evidence: ')


def _fail(message):
    raise EvidenceError('cross-period evidence: ' + message)


def _keys(value, keys, where):
    if not isinstance(value, dict) or set(value) != set(keys):
        _fail(where + ' has missing or unsupported fields')


def _text(value, where):
    if not isinstance(value, str) or not value.strip():_fail(where + ' must be a nonempty string')
    return value


def _int(value, where, minimum=0, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        _fail(where + ' must be a bounded integer')
    return value


def _bool(value, where, nullable=False):
    if type(value) is not bool and not (nullable and value is None):_fail(where + ' must be boolean' + (' or null' if nullable else ''))


def _list(value, where):
    if not isinstance(value, list):_fail(where + ' must be a list')
    return value


def _sha(value, where):
    if not isinstance(value, str) or not SHA.fullmatch(value):_fail(where + ' must be lowercase SHA256')
    return value


def _canonical(motif):
    if not isinstance(motif, str) or not 1 <= len(motif) <= 64 or not DNA.fullmatch(motif):
        _fail('motif must contain 1..64 uppercase ACGT bases')
    for n in range(1, len(motif) + 1):
        if len(motif) % n == 0 and motif == motif[:n] * (len(motif) // n):
            motif = motif[:n];break
    rc = motif.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    return min(s[i:] + s[:i] for s in (motif, rc) for i in range(len(s)))


def _unique_object(pairs):
    return _IO.unique_object(pairs)


def _constant(value):
    return _IO.reject_constant(value)


def evidence_payload_sha256(document):
    return _IO.payload_sha256(document)


def _read_bytes(path):
    return _IO.read_bytes(path, max_bytes=MAX_ARTIFACT_BYTES)


def _read_json(path):
    return _IO.decode_json(_read_bytes(path), path,
        object_pairs_hook=_unique_object, parse_constant=_constant)


def _source_pins(doc, base):
    pins = _list(doc.get('source_artifacts'), 'source_artifacts')
    if not pins or len(pins) > MAX_SOURCE_ARTIFACTS:_fail('source artifact count outside supported bounds')
    result = []; seen = set()
    for pin in pins:
        _keys(pin, {'path', 'sha256'}, 'source artifact')
        path = Path(_text(pin['path'], 'source path'))
        path = (base / path).absolute() if not path.is_absolute() else path
        canonical_path = str(path.resolve())
        if canonical_path in seen:_fail('duplicate source artifact path')
        seen.add(canonical_path)
        result.append(dict(path=str(path), sha256=_sha(pin['sha256'], 'source SHA256')))
    return result


def source_dependencies(path):
    """Cache dependencies only. Call load_evidence for validation."""
    path = Path(path).absolute(); doc, _ = _read_json(path)
    pins = _source_pins(doc, path.parent)
    total = 0
    for pin in pins:
        try:
            stamp = Path(pin['path']).stat()
        except OSError as error:
            raise EvidenceError('cross-period evidence: cannot inspect source ' + pin['path']) from error
        if not stat.S_ISREG(stamp.st_mode) or stamp.st_size > MAX_ARTIFACT_BYTES:
            _fail('source must be a small regular file: ' + pin['path'])
        total += stamp.st_size
        if total > MAX_SOURCE_BYTES:_fail('total source artifact bytes exceed supported bounds')
    return [str(path), *(p['path'] for p in pins)]


def _operations(ranges, cost, n):
    _keys(ranges, {'match', 'substitution', 'insertion', 'deletion'}, 'operation ranges')
    for operation, bounds in ranges.items():
        _keys(bounds, {'min', 'max'}, 'operation range')
        _int(bounds['min'], 'operation min', maximum=n)
        _int(bounds['max'], 'operation max', bounds['min'], n)
        if operation != 'match' and bounds['max'] > cost:_fail('operation maximum exceeds edit cost')
    # Necessary aggregate bounds only: marginal ranges must not be combined as one path.
    if not sum(ranges[k]['min'] for k in ('match','substitution','insertion')) <= n <= sum(ranges[k]['max'] for k in ('match','substitution','insertion')):
        _fail('operation ranges cannot consume the query')
    if not sum(ranges[k]['min'] for k in ('substitution','insertion','deletion')) <= cost <= sum(ranges[k]['max'] for k in ('substitution','insertion','deletion')):
        _fail('operation ranges inconsistent with cost')


def _validate_interval(row):
    _keys(row, _INTERVAL_FIELDS, 'interval')
    for key in ('interval_id', 'pair_query_id'):_text(row[key], 'interval.' + key)
    _int(row['mate'], 'mate', 1, 2)
    start = _int(row['start0'], 'start0'); end = _int(row['end0'], 'end0', start + 1)
    _int(row['full_read_length'], 'full_read_length', end)
    n = end - start
    if n > 1000:_fail('interval exceeds validated diagnostic query bound')
    if not isinstance(row['sequence'], str) or not DNA.fullmatch(row['sequence']) or len(row['sequence']) != n:_fail('interval sequence/coordinate mismatch')
    q = row['quality']
    if not isinstance(q, str) or len(q) != n or any(not 33 <= ord(c) <= 126 for c in q):_fail('interval quality must be printable ASCII of query length')
    if row['deletion_quality'] is not None:_fail('model deletion quality must be null')
    for role in ('candidate', 'comparator'):
        cost = _int(row[role + '_cost'], role + ' cost', maximum=n)
        _operations(row[role + '_operations'], cost, n)
        _int(row[role + '_ambiguous_operation_positions'], 'ambiguous operation positions', maximum=n)
    seen = set()
    for pos in _list(row['strict_contrasts'], 'strict contrasts'):
        _keys(pos, _CONTRAST_FIELDS, 'strict contrast')
        at = _int(pos['absolute_position0'], 'contrast coordinate', start, end - 1)
        if at in seen:_fail('duplicate strict contrast coordinate')
        seen.add(at)
        if pos['query_base'] != row['sequence'][at - start]:_fail('contrast observed base mismatch')
        if pos['quality_character'] != q[at - start]:_fail('contrast quality character mismatch')
        if _int(pos['quality_ascii'], 'quality ASCII', 33, 126) != ord(q[at - start]):_fail('contrast quality ASCII mismatch')
        if _int(pos['conditional_phred33'], 'conditional Phred33', maximum=93) != ord(q[at - start]) - 33:_fail('conditional Phred33 mismatch')
        for role in ('candidate', 'comparator'):
            ops = _list(pos[role + '_operations'], 'contrast operations')
            if not ops or len(set(ops)) != len(ops) or any(x not in ('match', 'substitution', 'insertion') for x in ops):_fail('invalid observed-base operation set')
            if any(row[role + '_operations'][op]['max'] == 0 for op in ops):_fail('position operation absent from every complete optimal path')
            bases = _list(pos[role + '_template_bases'], 'contrast template bases')
            if len(set(bases)) != len(bases) or any(x not in ('A','C','G','T') for x in bases):_fail('invalid template base set')
            if ('match' in ops or 'substitution' in ops) and not bases:_fail('template-consuming operation requires a template base')
            if 'match' in ops and pos['query_base'] not in bases:_fail('match operation lacks observed template base')
            if 'match' not in ops and pos['query_base'] in bases:_fail('nonmatching operation contains observed template base')
            if ops == ['insertion'] and bases:_fail('insertion-only position cannot have template bases')
            if 'substitution' in ops and not any(b != pos['query_base'] for b in bases):_fail('substitution lacks a different template base')
            if 'substitution' not in ops and any(b != pos['query_base'] for b in bases):_fail('template base requires an unreported substitution')
        ca, co = pos['candidate_operations'], pos['comparator_operations']
        expected = ('candidate_always_match_comparator_always_observed_edit' if ca == ['match'] and 'match' not in co else
                    'comparator_always_match_candidate_always_observed_edit' if co == ['match'] and 'match' not in ca else None)
        if expected is None or pos['contrast'] != expected:_fail('contrast is not strict across all optimal operations')


def _fcs(value):
    _keys(value, {'assessed','span_overlap','aligned_base_overlap_bp','published_actions'}, 'FCS')
    _bool(value['assessed'], 'FCS assessed'); _bool(value['span_overlap'], 'FCS span overlap', True)
    if value['aligned_base_overlap_bp'] is not None:_int(value['aligned_base_overlap_bp'], 'FCS overlap')
    if value['published_actions'] is not None:
        for x in _list(value['published_actions'], 'FCS actions'):_text(x, 'FCS action')
    if not value['assessed'] and any(value[k] is not None for k in ('span_overlap','aligned_base_overlap_bp','published_actions')):_fail('unassessed FCS must be null')


def _validate_context(row, events, intervals):
    _keys(row, _CONTEXT_FIELDS, 'context')
    for key in ('joint_id','event_id','reference_accession','target'):_text(row[key], 'context.' + key)
    if row['event_id'] not in events:_fail('context refers to unknown event')
    event = events[row['event_id']]; interval = intervals[event['interval_id']]
    if row['mode'] not in ('exact','gapped'):_fail('unsupported context mode')
    if row['full_read_strand'] not in ('+','-'):_fail('invalid alignment strand')
    for key in ('flank_hsp_id','full_read_hsp_id'):
        if row['mode'] == 'gapped':_text(row[key], key)
        elif row[key] is not None:_fail('exact HSP IDs must be null')
    _int(row['full_read_retained_hit_index0'], 'retained hit index')
    for prefix in ('full_read_query_', 'full_read_reference_'):
        a = _int(row[prefix+'start0'], prefix+'start0'); b = _int(row[prefix+'end0'], prefix+'end0', a+1)
        if prefix == 'full_read_query_' and b > interval['full_read_length']:_fail('full-read HSP outside query')
    br = _int(row['boundary_read0'], 'boundary read coordinate', maximum=interval['full_read_length'])
    _int(row['boundary_reference0'], 'boundary reference coordinate')
    if br not in (interval['start0'], interval['end0']):_fail('boundary must touch saved native tract')
    repeat_left_in_query = br == interval['end0']
    repeat_left_in_reference = repeat_left_in_query == (row['full_read_strand'] == '+')
    if not row['full_read_query_start0'] < br < row['full_read_query_end0']:_fail('full-read HSP must cross the boundary')
    if not row['full_read_reference_start0'] <= row['boundary_reference0'] <= row['full_read_reference_end0']:_fail('boundary outside reference HSP')
    if row['mode'] == 'exact':
        x=row['exact_flank_hit'];_keys(x, {'start0','end0','strand','target'}, 'exact flank hit')
        _int(x['start0'], 'exact flank start');_int(x['end0'], 'exact flank end', x['start0']+1)
        if x['target'] != row['target'] or x['strand'] not in ('+','-'):_fail('exact flank/full-read target mismatch')
        # Exact flank sequence can use the outward normalized frame, opposite to
        # the original full read. Do not equate their strand labels.
        flank_boundary = x['start0'] if repeat_left_in_reference else x['end0']
        if row['boundary_reference0'] != flank_boundary:_fail('exact flank must touch the nonrepeat side of this boundary')
    elif row['exact_flank_hit'] is not None:_fail('gapped context cannot have exact flank hit')
    for key in ('mate_geometry_evaluable','boundary_reference_continuous','boundary_adjacent_acgt_matches','technical_boundary_warning','search_capped'):_bool(row[key], key)
    _bool(row['same_hsp_inward_witness'], 'inward witness', True)
    if row['mate_geometry_evaluable'] != (row['same_hsp_inward_witness'] is not None):_fail('geometry evaluability/null mismatch')
    mate=row['closest_mate']
    if row['same_hsp_inward_witness'] is not True:
        if mate is not None:_fail('closest inward mate supplied without inward witness')
    else:
        _keys(mate, {'five_prime_boundary0','five_prime_separation_bp','hsp_id','library_insert_compatibility','retained_hit_index0','strand','support_five_prime_boundary0','target'}, 'inward mate')
        for k in ('five_prime_boundary0','five_prime_separation_bp','retained_hit_index0','support_five_prime_boundary0'):_int(mate[k], 'mate.'+k)
        if mate['hsp_id'] is not None:_text(mate['hsp_id'], 'mate HSP')
        if mate['library_insert_compatibility'] is not None:_fail('insert compatibility is unassessed')
        support=row['full_read_reference_start0'] if row['full_read_strand']=='+' else row['full_read_reference_end0']
        if mate['target'] != row['target'] or mate['support_five_prime_boundary0'] != support:_fail('mate witness belongs to another full-read HSP')
        expected='-' if row['full_read_strand']=='+' else '+'
        inward=(support <= mate['five_prime_boundary0']) if row['full_read_strand']=='+' else (mate['five_prime_boundary0'] <= support)
        if mate['strand'] != expected or not inward or mate['five_prime_separation_bp'] != abs(support-mate['five_prime_boundary0']):_fail('invalid same-HSP inward geometry')
    n=interval['end0']-interval['start0'];covered=_int(row['full_hsp_observed_interval_bases'], 'observed interval HSP bases', maximum=n)
    overlap=max(0,min(row['full_read_query_end0'],interval['end0'])-max(row['full_read_query_start0'],interval['start0']))
    if covered > overlap:_fail('observed interval HSP bases exceed coordinate overlap')
    seen=set();valid={p['absolute_position0'] for p in interval['strict_contrasts']}
    for x in _list(row['strict_positions_in_full_hsp'], 'HSP strict positions'):
        _int(x, 'HSP strict position')
        if x in seen or x not in valid or not row['full_read_query_start0'] <= x < row['full_read_query_end0']:_fail('invalid or duplicate strict HSP position')
        seen.add(x)
    if len(seen)>covered:_fail('strict positions exceed observed interval HSP coverage')
    comp=row['compartment']
    _keys(comp, {'assembly_unit','assigned_molecule','assigned_molecule_location_type','reference_sequence_bp','sequence_role'}, 'compartment')
    _int(comp['reference_sequence_bp'], 'reference sequence bp', row['full_read_reference_end0'])
    for key in set(comp)-{'reference_sequence_bp'}:_text(comp[key], 'compartment.'+key)
    run=row['selected_expected_run']
    if run is not None:
        _keys(run, {'complete_periods','contiguous_exact_bp','contiguous_exact_repeat_ge_500','end0','family_sequence_in_reference_forward_frame','matching_phase_indices_at_boundary_base','run_reaches_context_limit','run_reaches_sequence_edge','sequence_beyond_run_unresolved','start0','stop_position0','stop_reason'}, 'adjacent run')
        a=_int(run['start0'],'run start');b=_int(run['end0'],'run end',a+1,comp['reference_sequence_bp'])
        length=_int(run['contiguous_exact_bp'],'run bp',1)
        if length != b-a or row['boundary_reference0'] not in (a,b):_fail('selected run must touch this boundary')
        repeat_boundary = b if repeat_left_in_reference else a
        if row['boundary_reference0'] != repeat_boundary:_fail('selected run is on the wrong side of this boundary')
        if _canonical(run['family_sequence_in_reference_forward_frame']) != event['family']:_fail('adjacent run family mismatch')
        if _int(run['complete_periods'],'complete periods') != length//len(event['family']):_fail('adjacent run periods mismatch')
        for key in ('contiguous_exact_repeat_ge_500','run_reaches_context_limit','run_reaches_sequence_edge','sequence_beyond_run_unresolved'):_bool(run[key], key)
        if run['contiguous_exact_repeat_ge_500'] != (length>=500):_fail('run length flag mismatch')
        phases=_list(run['matching_phase_indices_at_boundary_base'],'boundary phases')
        if not phases or len(set(phases))!=len(phases):_fail('invalid boundary phases')
        for p in phases:_int(p,'phase',maximum=len(event['family'])-1)
        if run['stop_position0'] is not None:_int(run['stop_position0'],'run stop',-1,comp['reference_sequence_bp'])
        _text(run['stop_reason'],'run stop reason')
    _keys(row['aligned_technical_overlap'], {'12','16','20'}, 'technical overlap')
    for x in row['aligned_technical_overlap'].values():_int(x, 'technical overlap', maximum=interval['full_read_length'])
    for key in ('flank_fcs','full_read_fcs'):_fcs(row[key])
    if row['mate_fcs'] is not None:_fcs(row['mate_fcs'])


def _validate_sample(sample):
    _keys(sample, {'sample_id','source_result_sha256',*IDENTITY_FIELDS,'motifs','review_notes','intervals','events','contexts'}, 'sample')
    _text(sample['sample_id'],'sample_id');_sha(sample['source_result_sha256'],'native SHA256')
    for field in IDENTITY_FIELDS:_text(sample[field],field)
    if not sample['tax_id'].isdigit():_fail('tax_id must be a digit string')
    _keys(sample['motifs'], {'candidate','comparator'}, 'motifs')
    families={k:_canonical(v) for k,v in sample['motifs'].items()}
    if len(set(families.values()))!=2:_fail('candidate/comparator must represent distinct motif families')
    if len(families['candidate'])==len(families['comparator']):_fail('cross-period evidence requires unequal primitive periods')
    for note in _list(sample['review_notes'], 'review notes'):_text(note,'review note')
    intervals={};physical=set();read_lengths={};read_bases={}
    for row in _list(sample['intervals'],'intervals'):
        _validate_interval(row)
        key=(row['pair_query_id'],row['mate'],row['start0'],row['end0'])
        if row['interval_id'] in intervals or key in physical:_fail('duplicate interval or physical interval identity')
        read=(row['pair_query_id'],row['mate'])
        if read in read_lengths and read_lengths[read]!=row['full_read_length']:_fail('inconsistent full read length')
        read_lengths[read]=row['full_read_length']
        for offset,(base,quality) in enumerate(zip(row['sequence'],row['quality'])):
            coordinate=(*read,row['start0']+offset)
            if coordinate in read_bases and read_bases[coordinate]!=(base,quality):_fail('overlapping intervals disagree on observed base or quality')
            read_bases[coordinate]=(base,quality)
        physical.add(key);intervals[row['interval_id']]=row
    events={};event_identity=set()
    for row in _list(sample['events'],'events'):
        _keys(row,_EVENT_FIELDS,'event')
        for key in _EVENT_FIELDS:_text(row[key],'event.'+key)
        if row['event_id'] in events:_fail('duplicate event_id')
        if row['interval_id'] not in intervals:_fail('event refers to unknown interval')
        if row['role'] not in families or _canonical(row['declared_motif'])!=row['family'] or row['family']!=families[row['role']]:_fail('event role/motif/family mismatch')
        key=(row['interval_id'],row['family'],row['role'])
        if key in event_identity:_fail('duplicate physical event identity')
        event_identity.add(key);events[row['event_id']]=row
    if set(intervals)!={e['interval_id'] for e in events.values()}:_fail('unreferenced intervals')
    joints=set();identities=set()
    for row in _list(sample['contexts'],'contexts'):
        _validate_context(row,events,intervals)
        identity=(row['event_id'],row['reference_accession'],row['mode'],row['target'],row['flank_hsp_id'],row['full_read_hsp_id'],row['full_read_retained_hit_index0'],row['full_read_strand'],row['full_read_reference_start0'],row['full_read_reference_end0'],json.dumps(row['exact_flank_hit'],sort_keys=True))
        if row['joint_id'] in joints or identity in identities:_fail('duplicate context identity')
        joints.add(row['joint_id']);identities.add(identity)


class ValidatedEvidence(Mapping):
    """Validated sample mapping; every exposed mutable object is an isolated copy."""
    def __init__(self,samples,*,path,sha256,pins,payload_sha256):
        self._samples=deepcopy(samples);self._pins=deepcopy(pins)
        self.path=str(path);self.sha256=sha256;self.payload_sha256=payload_sha256
        self.source_paths=tuple(p['path'] for p in pins)
    def __getitem__(self,key):return deepcopy(self._samples[key])
    def __iter__(self):return iter(self._samples)
    def __len__(self):return len(self._samples)
    @property
    def provenance(self):
        return dict(sidecar_path=self.path,sidecar_sha256=self.sha256,payload_sha256=self.payload_sha256,source_artifacts=deepcopy(self._pins))


def load_evidence(path):
    """Validate seal, compact observations and bounded stable source hashes."""
    path=Path(path).absolute();doc,sha=_read_json(path)
    _keys(doc, {'schema_version','scope','source_artifacts','payload_sha256','samples'}, 'document')
    if type(doc['schema_version']) is not int or doc['schema_version']!=SCHEMA_VERSION:_fail('unsupported schema_version')
    if doc['scope']!=SCOPE:_fail('unexpected scope')
    if _sha(doc['payload_sha256'],'payload seal') != evidence_payload_sha256(doc):_fail('payload seal mismatch')
    samples={}
    for sample in _list(doc['samples'],'samples'):
        try:
            _validate_sample(sample)
        except (TypeError, KeyError, OverflowError) as error:
            raise EvidenceError('cross-period evidence: malformed observation field shape') from error
        if sample['sample_id'] in samples:_fail('duplicate sample_id')
        samples[sample['sample_id']]=sample
    if not samples:_fail('samples cannot be empty')
    pins=_source_pins(doc,path.parent)
    source_bytes=0
    for pin in pins:
        data=_read_bytes(pin['path']);source_bytes+=len(data)
        if source_bytes>MAX_SOURCE_BYTES:_fail('total source artifact bytes exceed supported bounds')
        if hashlib.sha256(data).hexdigest()!=pin['sha256']:_fail('source artifact SHA256 mismatch: '+pin['path'])
    return ValidatedEvidence(samples,path=path,sha256=sha,pins=pins,payload_sha256=doc['payload_sha256'])


def _positive(row):
    return row['selected_expected_run'] is not None and row['same_hsp_inward_witness'] is True


def _continuous(row):
    return _positive(row) and row['boundary_reference_continuous'] and row['boundary_adjacent_acgt_matches']


def _maxrun(rows, joint=False):
    return max((r['selected_expected_run']['contiguous_exact_bp'] for r in rows if r['selected_expected_run'] is not None and (not joint or _positive(r))),default=None)


def _summary(events,intervals,contexts):
    selected={e['interval_id']:intervals[e['interval_id']] for e in events}
    def count(predicate):return len({r['event_id'] for r in contexts if predicate(r)})
    evaluable=count(lambda r:r['mate_geometry_evaluable'])
    return dict(selected_event_count=len(events),selected_pair_count=len({r['pair_query_id'] for r in selected.values()}),
        selected_read_count=len({(r['pair_query_id'],r['mate']) for r in selected.values()}),selected_interval_count=len(selected),
        strict_contrast_position_count=sum(len(r['strict_contrasts']) for r in selected.values()),qualified_context_count=len(contexts),
        boundary_supported_event_count=count(lambda r:True),mate_geometry_evaluable_event_count=evaluable,
        inward_mate_event_count=count(lambda r:r['same_hsp_inward_witness'] is True) if evaluable else None,
        joint_repeat_inward_event_count=count(_positive),joint_continuous_matching_repeat_inward_event_count=count(_continuous),
        candidate_lower_cost_interval_count=sum(r['candidate_cost']<r['comparator_cost'] for r in selected.values()),
        comparator_lower_cost_interval_count=sum(r['comparator_cost']<r['candidate_cost'] for r in selected.values()),
        equal_cost_interval_count=sum(r['candidate_cost']==r['comparator_cost'] for r in selected.values()))


def build_cross_period_evidence(sample_id, source_report, evidence=None):
    """Return independent reporting fields bound to native hash and specimen."""
    _text(sample_id,'requested sample_id')
    if evidence is not None and not isinstance(evidence,ValidatedEvidence):_fail('call load_evidence before building diagnostics')
    result=dict(schema_version=SCHEMA_VERSION,status='unassessed',reason='evidence_not_supplied',scope=SCOPE,
        sample_id=sample_id,display_label='주기 간·동일 위치 진단 미평가',
        summary='검증된 저장 구간의 주기 간 비교·동일 위치 결합 근거가 없습니다.',
        recommendation='저장된 원래 구간·품질·동일 HSP 위치 근거를 검증한 선택 진단 파일을 제공하세요.',
        **dict.fromkeys(SUMMARY_FIELDS),family_summaries=[],intervals=[],events=[],contexts=[],positions=[],
        review_notes=[],limitations=list(LIMITATIONS),quality_encoding_verified=False,
        quality_authenticity_verified=False,deletion_quality=None,
        provenance=evidence.provenance if evidence is not None else None,
        native_counts_changed=False,classification_changed=False,adequacy_changed=False,**dict.fromkeys(_CONFIRMATIONS))
    if evidence is None:return result
    if sample_id not in evidence:
        result['reason']='sample_not_in_evidence';return result
    sample=evidence[sample_id]
    if not isinstance(source_report,dict) or source_report.get('sample_id')!=sample_id:_fail('current report sample_id mismatch')
    provenance=source_report.get('provenance',{})
    if not isinstance(provenance,dict):_fail('current report provenance must be an object')
    sha=provenance.get('original_discovery_sha256',source_report.get('source_result_sha256'))
    if _sha(sha,'native result SHA256')!=sample['source_result_sha256']:_fail('native result SHA256 mismatch')
    for field in IDENTITY_FIELDS:
        metadata=source_report.get('metadata',{})
        value=metadata.get('sample_accession') if field=='biosample' and isinstance(metadata,dict) else source_report.get(field)
        if value is None or str(value)!=sample[field]:_fail('current report '+field+' mismatch')
    if 'candidates' in source_report:
        present=set()
        for candidate in _list(source_report['candidates'],'current report candidates'):
            if not isinstance(candidate,dict):_fail('current candidate must be an object')
            if isinstance(candidate.get('canonical_family'),str):present.add(candidate['canonical_family'])
        if not {e['family'] for e in sample['events']}<=present:_fail('evidence family absent from current report')
    result['provenance'].update(source_result_sha256=sha,sample_identity={k:sample[k] for k in ('sample_id',*IDENTITY_FIELDS)})
    result['review_notes']=sample['review_notes']
    intervals={r['interval_id']:r for r in sample['intervals']};events={r['event_id']:r for r in sample['events']}
    contexts=sample['contexts'];byevent=defaultdict(list)
    for c in contexts:byevent[c['event_id']].append(c)
    result.update(_summary(list(events.values()),intervals,contexts))
    for row in sample['intervals']:
        out={k:row[k] for k in INTERVAL_FIELDS if k in row}
        delta=row['comparator_cost']-row['candidate_cost']
        out.update(interval_bp=row['end0']-row['start0'],candidate_motif=sample['motifs']['candidate'],comparator_motif=sample['motifs']['comparator'],
            comparator_minus_candidate_cost=delta,lower_cost_model='candidate' if delta>0 else 'comparator' if delta<0 else 'tie',strict_contrast_position_count=len(row['strict_contrasts']))
        result['intervals'].append(out)
        result['positions'].extend(dict(interval_id=row['interval_id'],pair_query_id=row['pair_query_id'],mate=row['mate'],**p) for p in row['strict_contrasts'])
    for event in sample['events']:
        rows=byevent[event['event_id']];interval=intervals[event['interval_id']]
        result['events'].append(dict(**event,pair_query_id=interval['pair_query_id'],mate=interval['mate'],qualified_context_count=len(rows),
            mate_geometry_evaluable_context_count=sum(r['mate_geometry_evaluable'] for r in rows),
            inward_mate_context_count=sum(r['same_hsp_inward_witness'] is True for r in rows) if any(r['mate_geometry_evaluable'] for r in rows) else None,
            joint_repeat_inward_context_count=sum(_positive(r) for r in rows),joint_continuous_matching_repeat_inward_context_count=sum(_continuous(r) for r in rows),
            strict_positions_in_full_hsp_count=len({p for r in rows for p in r['strict_positions_in_full_hsp']}),
            max_adjacent_reference_run_bp=_maxrun(rows),max_joint_adjacent_reference_run_bp=_maxrun(rows,True)))
    for row in contexts:
        out=deepcopy(row);event=events[row['event_id']]
        out.update(interval_id=event['interval_id'],family=event['family'],role=event['role'],selected_adjacent_reference_run_bp=row['selected_expected_run']['contiguous_exact_bp'] if row['selected_expected_run'] else None,
            joint_repeat_inward=_positive(row),joint_continuous_matching_repeat_inward=_continuous(row))
        result['contexts'].append(out)
    for family in sorted({e['family'] for e in events.values()}):
        selected=[e for e in events.values() if e['family']==family];ids={e['event_id'] for e in selected};rows=[c for c in contexts if c['event_id'] in ids]
        result['family_summaries'].append(dict(family=family,declared_motif=selected[0]['declared_motif'],role=selected[0]['role'],**_summary(selected,intervals,rows),
            max_adjacent_reference_run_bp=_maxrun(rows),max_joint_adjacent_reference_run_bp=_maxrun(rows,True)))
    result.update(status='reviewed_saved_subset',reason=None,display_label='저장 구간의 주기 간·동일 위치 진단 완료',
        summary=f"선정 {len(events)}사건·{len(intervals)}구간에서 구간별 관측 염기의 엄격한 모델 대비 {result['strict_contrast_position_count']}개를 기록했습니다. 같은 위치의 반복·inward mate 양성 사건은 {result['joint_repeat_inward_event_count']}건이며 생물학적 확정은 아닙니다.",
        recommendation='구간 선택 편향과 품질 미검증을 고려해 비용·실제 염기·동일 위치 근거를 별도로 검토하고, 적합한 긴 서열 또는 독립 표본으로 검증하세요.')
    if not events:result.update(status='unassessed',reason='no_selected_junction_events',display_label='선정 구간 없음')
    return deepcopy(result)
