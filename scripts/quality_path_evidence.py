"""Optional sealed, same-read quality/path annotations; never changes TRM calls.

This module reads only bounded completed JSON/TSV artifacts. It does not align,
reoptimize models, infer likelihoods, calibrate qualities or rescan sequences.
"""
from collections import defaultdict
from copy import deepcopy
import csv
import hashlib
import io
import json
from pathlib import Path
import re

from evidence_io import EvidenceIO

SCHEMA_VERSION = 1
SCOPE = 'conditional_saved_same_read_normalized_quality_and_all_optimal_paths'
MAX_BYTES = 8 * 1024 * 1024
MAX_TOTAL_BYTES = 64 * 1024 * 1024
MAX_SAMPLES = 100
MAX_ROWS = 10000
SHA = re.compile(r'^[0-9a-f]{64}$')
FIELDS = ('quality_path_status', 'quality_path_read_display', 'quality_path_scope_display',
          'quality_path_model_display', 'quality_path_limit_display')
SUMMARY_FIELDS = ('normalized_same_read_comparisons', 'modeled_reads', 'modeled_intervals',
    'ordinary_reads_without_models', 'focused_reads', 'focused_mixed_reads',
    'focused_mate_only_reads', 'focused_pairs', 'focused_mixed_intervals',
    'focused_mixed_union_bases', 'strict_fixed_position_contrast_intervals',
    'quality_origin_status', 'quality_calibration_status', 'physical_terminal_status',
    'main_trm_status', 'novel_trm_status', 'independent_validation_status')
MODEL_FIELDS = ('model', 'motif', 'focused_mixed_intervals', 'reads', 'pairs',
    'intervals_requiring_observed_edit', 'intervals_requiring_Q_ge30_observed_edit',
    'intervals_permitting_zero_observed_edit', 'intervals_with_exclusively_gap_edits',
    'intervals_requiring_deletion', 'deletion_required_without_fixed_mandatory_boundary',
    'strict_position_contrast_intervals', 'gap_base_quality', 'quality_used_to_optimize')
READ_FIELDS = ('read_query_id', 'spot_id', 'mate', 'bases', 'normalized_quality_bins',
    'normalized_variable', 'quality_differing_positions', 'model_quality_evaluated',
    'focused', 'focused_mixed_bearing', 'independent_biological_replication',
    'original_submitted_quality_confirmed')
SOURCE_ROLES = ('assessment_completion', 'prior_assessment', 'comparison_completion',
    'read_comparison', 'path_completion', 'read_coverage', 'model_summary', 'interval_pair_summary')
JSON_ROLES = {'assessment_completion', 'prior_assessment', 'comparison_completion', 'path_completion'}
MODELS = {'period5': 'TTAGG', 'period6': 'TTAGGG', 'period11': 'AGGGTTAGGTT'}
IDENTITY_FIELDS = ('scientific_name', 'tax_id', 'biosample')
PAYLOAD_FIELDS = ('schema_version', 'scope', 'sample_id', 'identity', 'original_native_sha256',
    'representative_motif', 'summary', 'models', 'read_comparisons', 'review_notes',
    'source_artifacts', 'payload_sha256')
LIMITATIONS = [
    'Selected motif-conditioned reads are a conditional subset, not a whole-library sensitivity or absence denominator.',
    'Normalized qualities come from the same reads and specimens; mates, models and repeated intervals are not independent biological replication.',
    'The saved constant-Q30 values differ from normalized same-read values; original conversion provenance and empirical quality calibration remain unresolved.',
    'A positive minimum Q>=30 edited-base count across optimal paths does not identify a fixed discriminating base position.',
    'The fixed-position contrast count applies to focused mixed intervals only; other saved intervals can have such contrasts.',
    'Gap/deletion boundaries have no observed base quality: their quality is null, never zero or Q30.',
    'All-optimal edit-count ranges describe equally minimal unweighted edit paths, not likelihoods, posterior probabilities or calibrated sequencing errors.',
    'Per-operation and quality-bin extrema need not be simultaneously attainable; qualities were not used to reoptimize paths.',
    'The period-11 model is one specific periodic sequence, not a model of every TTAGG/TTAGGG mixture.',
    'Allowing zero observed-base edits on one optimal path differs from requiring exclusively gap edits on every optimal path.',
    'This optional annotation confirms neither physical telomeres, host/compartment origin, main TRM nor novel TRM, and changes no native counts, ranking, classification or adequacy.',
]
ABSENT_LIMITATIONS = [
    'No optional quality/path evidence was provided for this sample; all additional measurements remain unassessed.',
    'Missing optional evidence does not indicate zero observations, absence, failed detection or a biological conclusion.',
    'This optional annotation changes no native counts, ranking, classification or adequacy.',
]


class EvidenceError(ValueError):
    pass


IO = EvidenceIO(EvidenceError, 'quality path evidence: ')
payload_sha256 = IO.payload_sha256


def _keys(value, keys, where):
    if not isinstance(value, dict) or set(value) != set(keys):
        IO.fail(where + ' has missing or unsupported fields')


def _text(value, where, maximum=4096):
    if not isinstance(value, str) or not value.strip() or len(value) > maximum:
        IO.fail(where + ' must be bounded nonempty text')
    return value


def _int(value, where, maximum=1000000000):
    if type(value) is not int or not 0 <= value <= maximum:
        IO.fail(where + ' must be a bounded nonnegative integer')
    return value


def _bool(value, where, expected=None):
    if type(value) is not bool or (expected is not None and value is not expected):
        IO.fail(where + ' must be the required boolean')
    return value


def _sha(value):
    if not isinstance(value, str) or not SHA.fullmatch(value): IO.fail('invalid SHA256')
    return value


def _tsv_int(value):
    if not isinstance(value, str) or not re.fullmatch(r'0|[1-9][0-9]*', value): IO.fail('invalid source TSV integer')
    return _int(int(value), 'source TSV integer')


def _tsv_bool(value):
    if value not in ('True', 'False'): IO.fail('invalid source TSV boolean')
    return value == 'True'


def _stamp(path):
    try:
        s = path.stat()
        return (str(path.resolve(strict=True)), s.st_dev, s.st_ino, s.st_size,
                s.st_mtime_ns, s.st_ctime_ns, s.st_mode)
    except (OSError, RuntimeError) as error:
        raise EvidenceError('quality path evidence: source unavailable: ' + str(path)) from error


def _json(data, path):
    return IO.decode_json(data, path, object_pairs_hook=IO.unique_object,
                          parse_constant=IO.reject_constant)[0]


def _nested(value):
    try:
        return json.loads(value, object_pairs_hook=IO.unique_object, parse_constant=IO.reject_constant)
    except (ValueError, TypeError) as error:
        raise EvidenceError('quality path evidence: invalid nested source JSON') from error


def _table(data, path):
    try:
        reader = csv.DictReader(io.StringIO(data.decode('utf-8')), delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames): IO.fail('duplicate/missing TSV columns')
        rows = list(reader)
        if len(rows) > MAX_ROWS or any(None in row or any(v is None for v in row.values()) for row in rows):
            IO.fail('source TSV shape or count invalid')
        return rows
    except (UnicodeError, csv.Error) as error:
        raise EvidenceError('quality path evidence: invalid TSV: ' + str(path)) from error


def _identity(document):
    if not isinstance(document, dict) or not isinstance(document.get('metadata') or {}, dict): IO.fail('assessment/metadata must be objects')
    return dict(scientific_name=document.get('scientific_name'), tax_id=str(document.get('tax_id', '')),
                biosample=(document.get('metadata') or {}).get('sample_accession'))


def _assessment(document, sample, identity, digest):
    if not isinstance(document, dict) or not isinstance(document.get('provenance') or {}, dict): IO.fail('assessment/provenance must be objects')
    if document.get('analysis_complete') is not True or document.get('sample_id') != sample: IO.fail('source assessment incomplete or sample differs')
    if document.get('native_scan_version') != '0.8.5': IO.fail('only reviewed native 0.8.5 assessments supported')
    if _identity(document) != identity: IO.fail('sample/taxon/BioSample identity differs')
    if (document.get('provenance') or {}).get('original_discovery_sha256') != digest: IO.fail('original native SHA differs')


def _index(rows, fields):
    result = {}
    for row in rows:
        key = tuple(_text(row[field], field) for field in fields)
        if key in result: IO.fail('duplicate source row identity')
        result[key] = row
    return result


def _range(value, fields, maximum):
    _keys(value, fields, 'range collection')
    for item in value.values():
        _keys(item, ('min', 'max'), 'range')
        low = _int(item['min'], 'range minimum', maximum)
        if not low <= _int(item['max'], 'range maximum', maximum): IO.fail('range minimum exceeds maximum')
    return value


def _derive(sample, sources):
    try:
        return _derive_checked(sample, sources)
    except (KeyError, TypeError, AttributeError) as error:
        raise EvidenceError('quality path evidence: detailed source shape is invalid') from error


def _derive_checked(sample, sources):
    prior = sources['prior_assessment']; identity = _identity(prior)
    native = (prior.get('provenance') or {}).get('original_discovery_sha256'); _sha(native)
    _assessment(prior, sample, identity, native)
    comparisons = _index(sources['read_comparison'], ('sample_id', 'query_id'))
    coverage = _index(sources['read_coverage'], ('sample_id', 'read_query_id'))
    all_models = _index(sources['model_summary'], ('sample_id', 'interval_id', 'model'))
    all_pairs = _index(sources['interval_pair_summary'], ('sample_id', 'interval_id', 'left_model', 'right_model'))
    comp = {key[1]: row for key, row in comparisons.items() if key[0] == sample}
    cov = {key[1]: row for key, row in coverage.items() if key[0] == sample}
    if not comp or set(comp) != set(cov): IO.fail('comparison/coverage read identities differ')
    if len({(r['spot_id'], r['mate']) for r in comp.values()}) != len(comp): IO.fail('duplicate same-read spot/mate identity')
    intervals = {}; selected_models = {}; focused = set(); mixed_reads = set()
    for (sid, iid, model), row in all_models.items():
        if sid != sample: continue
        if model not in MODELS or row['motif'] != MODELS[model]: IO.fail('unsupported reviewed model/motif')
        read = row['read_query_id']
        if read not in comp: IO.fail('modeled read lacks same-read quality comparison')
        c = comp[read]; mate = _tsv_int(row['mate']); spot = _tsv_int(row['spot_id'])
        if mate not in (1, 2) or str(mate) != c['mate'] or str(spot) != c['spot_id']: IO.fail('model spot/mate binding differs')
        expected_pair = sample + '__pair' + str(spot)
        if row['pair_query_id'] != expected_pair or read != expected_pair + '__R' + str(mate): IO.fail('model pair/read binding differs')
        a, b = _tsv_int(row['start0']), _tsv_int(row['end0'])
        if not 0 <= a < b <= _tsv_int(c['bases']): IO.fail('model interval outside read')
        if iid != f'{expected_pair}__mate{mate}__{a}_{b}': IO.fail('model interval identity differs')
        focus, mixed = _tsv_bool(row['focused']), _tsv_bool(row['mixed'])
        roles = _nested(cov[read]['roles'])
        if not isinstance(roles, list) or any(not isinstance(x, str) for x in roles): IO.fail('invalid read coverage roles')
        if focus != ('focused_oriented_mixed_boundary_pair' in roles): IO.fail('focused read binding differs')
        if _tsv_bool(row['quality_used_to_optimize']) or not _tsv_bool(row['count_extrema_jointly_attainable_not_asserted']): IO.fail('source optimization scope differs')
        state = _nested(row['states'])
        if not isinstance(state, dict) or not set(state) <= {'M', 'E', 'A'}: IO.fail('invalid observed-position states')
        if sum(_int(v, 'state count', b-a) for v in state.values()) != b-a: IO.fail('state count length differs')
        cost = _tsv_int(row['minimum_edit_cost'])
        op = _range(_nested(row['operation_count_ranges']), ('match', 'substitution', 'insertion', 'deletion'), 2*(b-a)+32)
        q = _range(_nested(row['path_quality_count_ranges']), ('observed_edits', 'observed_edits_Q_lt20', 'observed_edits_Q20_29', 'observed_edits_Q_ge30'), b-a)
        if cost > 2*(b-a)+32 or any(op[k]['min'] > cost or op[k]['max'] > cost for k in ('substitution', 'insertion', 'deletion')): IO.fail('operation count exceeds edit cost')
        if any(q[k]['min'] > q['observed_edits']['min'] or q[k]['max'] > q['observed_edits']['max'] for k in q if k != 'observed_edits'): IO.fail('quality-bin range exceeds observed edits')
        observed = q['observed_edits']
        if (observed['min'] != cost-op['deletion']['max'] or observed['max'] != cost-op['deletion']['min']
                or op['match']['min'] != b-a-observed['max'] or op['match']['max'] != b-a-observed['min']):
            IO.fail('fixed-cost observed-edit/deletion ranges differ')
        if state.get('E',0)>observed['min'] or state.get('M',0)>b-a-observed['max']:
            IO.fail('fixed observed-position state count exceeds path range')
        gap_only = _tsv_bool(row['exclusively_gap_edits'])
        if gap_only != (cost > 0 and q['observed_edits']['max'] == 0): IO.fail('gap-only claim differs from all-path range')
        always = _tsv_bool(row['all_observed_positions_always_match'])
        if always != (state.get('M', 0) == b-a): IO.fail('all-match state differs')
        mandatory = _tsv_int(row['mandatory_deletion_boundaries'])
        possible = _tsv_int(row['possible_deletion_boundaries'])
        if not mandatory <= possible <= b-a+1: IO.fail('deletion boundary bounds differ')
        no_fixed = _tsv_bool(row['deletion_required_but_no_mandatory_boundary'])
        if no_fixed != (op['deletion']['min'] > 0 and mandatory == 0): IO.fail('required deletion scope differs')
        binding = (read, expected_pair, a, b, focus, mixed)
        if iid in intervals and intervals[iid] != binding: IO.fail('model interval binding differs')
        intervals[iid] = binding
        selected_models[(iid, model)] = dict(row=row, op=op, q=q, cost=cost, gap_only=gap_only, no_fixed=no_fixed, state=state)
        if focus: focused.add(read)
        if focus and mixed: mixed_reads.add(read)
    if not intervals or set(selected_models) != {(iid, model) for iid in intervals for model in MODELS}: IO.fail('incomplete three-model interval coverage')
    strict = {model:set() for model in MODELS}; focus_strict = set()
    required_pairs = {('period5', 'period6'), ('period5', 'period11'), ('period6', 'period11')}
    selected_pairs = {(iid, left, right):row for (sid,iid,left,right),row in all_pairs.items() if sid == sample}
    if set(selected_pairs) != {(iid,left,right) for iid in intervals for left,right in required_pairs}: IO.fail('incomplete model-pair interval coverage')
    for (iid, left, right), row in selected_pairs.items():
        read, pair, a, b, focus, mixed = intervals[iid]
        if (row['read_query_id'], row['pair_query_id'], _tsv_bool(row['focused']), _tsv_bool(row['mixed'])) != (read,pair,focus,mixed): IO.fail('model-pair source binding differs')
        if _tsv_int(row['bases']) != b-a: IO.fail('model-pair bases differ')
        for side, model in (('left',left),('right',right)):
            saved = selected_models[(iid,model)]
            if _tsv_int(row[side+'_cost']) != saved['cost'] or _tsv_bool(row[side+'_only_gap_edits']) != saved['gap_only']: IO.fail('model-pair cost/gap binding differs')
        hist = _nested(row['state_pair_counts'])
        if not isinstance(hist,dict) or not set(hist) <= {a+b for a in 'MEA' for b in 'MEA'}: IO.fail('invalid paired states')
        if sum(_int(v,'paired-state count',b-a) for v in hist.values()) != b-a: IO.fail('paired-state length differs')
        for index,model in enumerate((left,right)):
            marginal={state:sum(count for pair,count in hist.items() if pair[index]==state) for state in 'MEA'}
            if any(marginal[state]!=selected_models[(iid,model)]['state'].get(state,0) for state in 'MEA'):
                IO.fail('paired-state marginals differ from model states')
        if focus and mixed and (hist.get('ME',0) or hist.get('EM',0)):
            focus_strict.add(iid); strict[left].add(iid); strict[right].add(iid)
    read_rows=[]; modeled={v[0] for v in intervals.values()}
    for read,c in sorted(comp.items()):
        bases = _tsv_int(c['bases']); diff = _tsv_int(c['quality_differing_positions'])
        mate = _tsv_int(c['mate']); spot = _tsv_int(c['spot_id']); bins=[_tsv_int(x) for x in c['normalized_quality_bins'].split(',')]
        if bases < 1 or mate not in (1,2) or diff > bases or not bins or len(set(bins)) != len(bins) or any(q>93 for q in bins): IO.fail('quality comparison counts/bins invalid')
        if not _tsv_bool(c['sequence_and_instrument_name_match']) or c['orientation']!='direct' or c['saved_quality']!='constant_Q30': IO.fail('same-read identity/quality scope differs')
        if _tsv_bool(c['independent_biological_replication']) or _tsv_bool(c['original_submitted_quality_confirmed']): IO.fail('unsupported quality/independence confirmation')
        variable = _tsv_bool(c['normalized_variable'])
        if variable != (len(bins)>1): IO.fail('variable quality claim differs')
        cv = cov[read]; assessed = _tsv_bool(cv['model_quality_evaluated'])
        count = sum(v[0]==read for v in intervals.values())
        if assessed != (read in modeled) or _tsv_int(cv['saved_intervals']) != count: IO.fail('model/read coverage count differs')
        roles = _nested(cv['roles'])
        if not isinstance(roles,list) or any(not isinstance(x,str) for x in roles): IO.fail('invalid read roles')
        if not assessed and 'ordinary_unselected_prefix_control_not_biological_negative' not in roles: IO.fail('unsupported unmodeled control scope')
        if ('focused_oriented_mixed_boundary_pair' in roles) != (read in focused): IO.fail('focused coverage/model mismatch')
        read_rows.append(dict(read_query_id=read,spot_id=spot,mate=mate,bases=bases,normalized_quality_bins=bins,
            normalized_variable=variable,quality_differing_positions=diff,model_quality_evaluated=assessed,
            focused=read in focused,focused_mixed_bearing=read in mixed_reads,
            independent_biological_replication=False,original_submitted_quality_confirmed=False))
    fm = {iid:v for iid,v in intervals.items() if v[4] and v[5]}
    union = defaultdict(list)
    for read,pair,a,b,focus,mixed in fm.values(): union[read].append((a,b))
    union_bases = 0
    for segments in union.values():
        previous_end = -1
        for a,b in sorted(segments):
            union_bases += max(0, b-max(a,previous_end))
            previous_end = max(previous_end,b)
    model_rows=[]
    for model,motif in MODELS.items():
        rows=[selected_models[(iid,model)] for iid in fm]
        model_rows.append(dict(model=model,motif=motif,focused_mixed_intervals=len(fm),reads=len(mixed_reads),pairs=len({v[1] for v in fm.values()}),
            intervals_requiring_observed_edit=sum(r['q']['observed_edits']['min']>0 for r in rows),
            intervals_requiring_Q_ge30_observed_edit=sum(r['q']['observed_edits_Q_ge30']['min']>0 for r in rows),
            intervals_permitting_zero_observed_edit=sum(r['q']['observed_edits']['min']==0 for r in rows),
            intervals_with_exclusively_gap_edits=sum(r['gap_only'] for r in rows),
            intervals_requiring_deletion=sum(r['op']['deletion']['min']>0 for r in rows),
            deletion_required_without_fixed_mandatory_boundary=sum(r['no_fixed'] for r in rows),
            strict_position_contrast_intervals=len(strict[model]),gap_base_quality=None,quality_used_to_optimize=False))
    summary=dict(normalized_same_read_comparisons=len(comp),modeled_reads=len(modeled),modeled_intervals=len(intervals),
        ordinary_reads_without_models=len(comp)-len(modeled),focused_reads=len(focused),focused_mixed_reads=len(mixed_reads),
        focused_mate_only_reads=len(focused-mixed_reads),focused_pairs=len({v[1] for v in intervals.values() if v[4]}),
        focused_mixed_intervals=len(fm),focused_mixed_union_bases=union_bases,
        strict_fixed_position_contrast_intervals=len(focus_strict),quality_origin_status='normalized_same_read_recovered_original_conversion_unresolved',
        quality_calibration_status='not_calibrated',physical_terminal_status='not_confirmed',main_trm_status='not_confirmed',
        novel_trm_status='not_confirmed',independent_validation_status='not_confirmed')
    return dict(identity=identity,original_native_sha256=native,representative_motif=prior.get('detected_motif'),
                summary=summary,models=model_rows,read_comparisons=read_rows)


def _validate_payload(doc):
    _keys(doc,PAYLOAD_FIELDS,'payload')
    if type(doc['schema_version']) is not int or doc['schema_version']!=SCHEMA_VERSION or doc['scope']!=SCOPE: IO.fail('unsupported payload schema/scope')
    if payload_sha256(doc)!=_sha(doc['payload_sha256']): IO.fail('payload seal differs')
    _text(doc['sample_id'],'sample ID'); _sha(doc['original_native_sha256'])
    _keys(doc['identity'],IDENTITY_FIELDS,'identity')
    for key,value in doc['identity'].items(): _text(value,key)
    if not re.fullmatch('[ACGT]{1,64}',_text(doc['representative_motif'],'representative motif')): IO.fail('invalid representative motif')
    _keys(doc['summary'],SUMMARY_FIELDS,'summary')
    for key,value in doc['summary'].items():
        if key.endswith('_status'): _text(value,key)
        else: _int(value,key)
    if not isinstance(doc['models'],list) or len(doc['models'])!=3: IO.fail('three model rows required')
    if {r.get('model') for r in doc['models'] if isinstance(r,dict)}!=set(MODELS): IO.fail('duplicate/missing model ID')
    for row in doc['models']:
        _keys(row,MODEL_FIELDS,'model row')
        for key,value in row.items():
            if key in ('model','motif'): _text(value,key)
            elif key=='gap_base_quality':
                if value is not None: IO.fail('gap base quality must be null')
            elif key=='quality_used_to_optimize': _bool(value,key,False)
            else: _int(value,key)
    if not isinstance(doc['read_comparisons'],list) or not 1<=len(doc['read_comparisons'])<=MAX_ROWS: IO.fail('read count outside bound')
    seen=set()
    for row in doc['read_comparisons']:
        _keys(row,READ_FIELDS,'read comparison')
        read=_text(row['read_query_id'],'read ID')
        if read in seen: IO.fail('duplicate read ID')
        seen.add(read)
        for key in ('spot_id','mate','bases','quality_differing_positions'): _int(row[key],key)
        for key in ('normalized_variable','model_quality_evaluated','focused','focused_mixed_bearing'): _bool(row[key],key)
        for key in ('independent_biological_replication','original_submitted_quality_confirmed'): _bool(row[key],key,False)
        if not isinstance(row['normalized_quality_bins'],list) or not 1<=len(row['normalized_quality_bins'])<=94: IO.fail('quality bins outside bound')
        for value in row['normalized_quality_bins']: _int(value,'quality bin',93)
    if not isinstance(doc['review_notes'],list) or len(doc['review_notes'])>16: IO.fail('review notes outside bound')
    for note in doc['review_notes']: _text(note,'review note',2000)


def _bind_sources(sources,paths,hashes):
    groups={'assessment_completion':('prior_assessment',),'comparison_completion':('read_comparison',),
            'path_completion':('read_coverage','model_summary','interval_pair_summary')}
    for receipt_role,roles in groups.items():
        receipt=sources[receipt_role]
        if not isinstance(receipt,dict) or receipt.get('status')!='complete' or not isinstance(receipt.get('artifacts_sha256'),dict): IO.fail('source completion receipt invalid')
        for role in roles:
            try: key=str(paths[role].resolve().relative_to(paths[receipt_role].parent.resolve()))
            except ValueError: IO.fail('source outside completion receipt root')
            if receipt['artifacts_sha256'].get(key)!=hashes[role]: IO.fail('source differs from completed artifact seal')


def load_evidence(path):
    if path is None: return None
    path=Path(path).absolute(); cache={}; stamps={}; total=0
    def read(source):
        nonlocal total
        source=Path(source).absolute(); key=str(source)
        if key not in cache:
            before=_stamp(source); data=IO.read_bytes(source,max_bytes=MAX_BYTES); after=_stamp(source)
            if before!=after: IO.fail('source changed around bounded read')
            total+=len(data)
            if total>MAX_TOTAL_BYTES: IO.fail('total evidence bytes exceed bound')
            cache[key]=(data,hashlib.sha256(data).hexdigest()); stamps[key]=after
        return cache[key]
    raw,digest=read(path); manifest=_json(raw,path)
    _keys(manifest,('schema_version','scope','samples','payload_sha256'),'manifest')
    if type(manifest['schema_version']) is not int or manifest['schema_version']!=SCHEMA_VERSION or manifest['scope']!=SCOPE: IO.fail('manifest schema/scope unsupported')
    if payload_sha256(manifest)!=_sha(manifest['payload_sha256']): IO.fail('manifest seal differs')
    if not isinstance(manifest['samples'],list) or not 1<=len(manifest['samples'])<=MAX_SAMPLES: IO.fail('manifest sample count outside bound')
    result=dict(manifest_sha256=digest,samples={},paths=[],stamps=stamps)
    for row in manifest['samples']:
        _keys(row,('sample_id','identity','original_native_sha256','evidence'),'manifest sample')
        sample=_text(row['sample_id'],'sample ID')
        if sample in result['samples']: IO.fail('duplicate sample ID')
        _keys(row['identity'],IDENTITY_FIELDS,'manifest identity')
        for value in row['identity'].values(): _text(value,'manifest identity')
        _sha(row['original_native_sha256']); _keys(row['evidence'],('path','sha256'),'evidence pin')
        pp=Path(_text(row['evidence']['path'],'payload path')); pp=pp if pp.is_absolute() else path.parent/pp
        data,actual=read(pp)
        if actual!=_sha(row['evidence']['sha256']): IO.fail('sample evidence hash differs')
        doc=_json(data,pp); _validate_payload(doc)
        for key in ('sample_id','identity','original_native_sha256'):
            if row[key]!=doc[key]: IO.fail('manifest/payload identity differs')
        pins=doc['source_artifacts']
        if not isinstance(pins,list) or len(pins)!=len(SOURCE_ROLES): IO.fail('source artifact count differs')
        sources={}; paths={}; hashes={}; seen=set()
        for pin in pins:
            _keys(pin,('role','path','sha256'),'source pin'); role=pin['role']
            if not isinstance(role,str) or role not in SOURCE_ROLES or role in sources: IO.fail('missing/duplicate/unsupported source role')
            source=Path(_text(pin['path'],'source path')); source=source if source.is_absolute() else pp.parent/source
            if str(source.resolve()) in seen: IO.fail('duplicate source path')
            seen.add(str(source.resolve()))
            if source.suffix!=('.json' if role in JSON_ROLES else '.tsv'): IO.fail('source format not supported')
            data,actual=read(source)
            if actual!=_sha(pin['sha256']): IO.fail('detailed source hash differs')
            sources[role]=_json(data,source) if role in JSON_ROLES else _table(data,source)
            paths[role]=source.absolute(); hashes[role]=actual
        _bind_sources(sources,paths,hashes)
        rebuilt=_derive(sample,sources)
        for key,expected in rebuilt.items():
            if doc[key]!=expected: IO.fail('exported '+key+' differs from detailed sealed sources')
        result['samples'][sample]=dict(document=doc,validated_payload_sha256=doc['payload_sha256'],path=str(pp.absolute()),sha256=row['evidence']['sha256'],
            source_artifacts=[dict(role=role,path=str(paths[role]),sha256=hashes[role]) for role in SOURCE_ROLES])
    result['paths']=sorted(stamps)
    return result


def source_dependencies(path):
    return [] if path is None else load_evidence(path)['paths']


def build_evidence(sample,document,loaded):
    result=dict(status='not_provided',display=dict(zip(FIELDS,('품질·경로 추가 근거 미제공','같은 read 품질 비교 미평가',
        '선정 mixed 구간 경로 미평가','모델별 품질 제약 미평가','별도 검토 자료가 제공되지 않음'))),
        summary={key:None for key in SUMMARY_FIELDS},models=[],read_comparisons=[],scope=SCOPE,
        limitations=list(ABSENT_LIMITATIONS),changes_classification_or_adequacy=False,source_artifacts=[])
    if loaded is None or sample not in loaded['samples']: return result
    for path,stamp in loaded['stamps'].items():
        if _stamp(Path(path))!=stamp: IO.fail('loaded evidence is stale')
    row=loaded['samples'][sample]; doc=row['document']; _validate_payload(doc)
    # Recheck the canonical seal captured on load, not just a caller-resealable object.
    if doc['payload_sha256']!=row['validated_payload_sha256']: IO.fail('loaded payload changed since source validation')
    _assessment(document,sample,doc['identity'],doc['original_native_sha256'])
    if document.get('detected_motif')!=doc['representative_motif']: IO.fail('representative motif differs; overlay cannot replace it')
    summary=doc['summary']; models=doc['models']
    display=dict(quality_path_status='파일·표본 연결 무결성 검증됨 · 생물학적 확정 아님',
        quality_path_read_display=f"정규화된 같은 read 품질 {summary['normalized_same_read_comparisons']}건 비교 · 기존 모델 포함 {summary['modeled_reads']} reads · 원본 품질 변환 경위·오류율 보정 미해결",
        quality_path_scope_display=f"선정 mixed {summary['focused_mixed_intervals']}구간/{summary['focused_mixed_reads']} reads/{summary['focused_pairs']} pairs · 고정 위치 대조 {summary['strict_fixed_position_contrast_intervals']}구간 · mate-only {summary['focused_mate_only_reads']} reads 별도",
        quality_path_model_display='; '.join(f"{r['motif']}: 모든 최적 경로에서 Q≥30 관측 염기 편집 필요 {r['intervals_requiring_Q_ge30_observed_edit']}/{r['focused_mixed_intervals']}구간" for r in models),
        quality_path_limit_display='고정 위치 대조와 경로별 필수 편집 수는 다름 · gap 품질 없음 · 우도/오류율/물리적 말단/main/novel 확정 아님' + ('; '+'; '.join(doc['review_notes']) if doc['review_notes'] else ''))
    result.update(status='sources_validated',display=display,summary=deepcopy(summary),models=deepcopy(models),
        read_comparisons=deepcopy(doc['read_comparisons']),limitations=list(LIMITATIONS),
        source_artifacts=deepcopy(row['source_artifacts']),manifest_sha256=loaded['manifest_sha256'],evidence_sha256=row['sha256'],
        original_native_sha256=doc['original_native_sha256'],representative_motif=doc['representative_motif'],review_notes=list(doc['review_notes']))
    return result
