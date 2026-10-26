"""Validate sealed saved-junction reference evidence as a reporting-only overlay.

No FASTQ, assembly, aligner, native count or classification is touched.  The
canonical payload seal checks integrity, not authorship or biological truth.
"""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
from copy import deepcopy
import hashlib
# Retain this module binding for callers instrumenting descriptor reads.
import os
from pathlib import Path
import re

from evidence_io import EvidenceIO

SCHEMA_VERSION = 1
SCOPE = 'conditional_saved_junction_subset'
MAX_ARTIFACT_BYTES = 64 * 1024 * 1024
MAX_MOTIF_LABEL_BP = 256
_SHA = re.compile(r'^[0-9a-f]{64}$')
_DNA = re.compile(r'^[ACGT]+$')
_TOP = {'schema_version', 'scope', 'source_artifacts', 'payload_sha256', 'references', 'samples'}
_REFERENCE = {'reference_accession', 'sample_id', 'repository_taxon', 'role',
              'sample_linkage', 'verification', 'linkage_limit'}
_SAMPLE_REQUIRED = {'sample_id', 'source_result_sha256', 'selected_event_count',
                    'selected_pair_count', 'review_notes', 'events'}
_IDENTITY = ('scientific_name', 'tax_id', 'biosample')
_COUNTS = ('retained_flank_hit_count', 'retained_context_rows', 'assessed_context_rows',
           'full_read_coordinate_supported_assessed_context_rows',
           'raw_ungated_bridge_associations', 'qualified_bridge_associations',
           'qualified_bridge_adjacent_reference_continuous_associations',
           'qualified_bridge_adjacent_reference_acgt_matching_associations')
_BOOLEANS = ('assessed_flank_context_present', 'same_flank_locus_mate_geometry_evaluable',
             'flank_search_cap_flag', 'any_event_query_search_cap_flag', 'technical_boundary_warning')
_NULLABLE_BOOLEANS = ('event_has_same_flank_locus_inward_mate_witness',
    'pair_has_any_inward_witness_any_retained_locus', 'adjacent_expected_exact_run_ge500',
    'selected_expected_run_reaches_reference_sequence_edge', 'any_ambiguous_query_not_exact_assessable',
    'fcs_overlay_assessed', 'flank_fcs_overlap', 'any_event_query_fcs_overlap')
_RUN_LENGTHS = ('selected_expected_run_max_bp',
               'selected_expected_run_max_bp_with_full_read_coordinate_support',
               'expected_orientation_measured_run_max_bp')
_CONFIRMATIONS = ('biological_origin_confirmed', 'physical_terminal_confirmed',
                  'main_trm_confirmed', 'novel_trm_confirmed')
EVENT_FIELDS = ('reference_accession', 'sample_id', 'mode', 'family', 'role', 'declared_motif',
    'event_id', 'pair_query_id', 'state', *_COUNTS[:3], 'assessed_flank_context_present',
    *_COUNTS[3:], 'same_flank_locus_mate_geometry_evaluable',
    'event_has_same_flank_locus_inward_mate_witness',
    'smallest_retained_same_flank_locus_inward_five_prime_separation_bp',
    'pair_has_any_inward_witness_any_retained_locus', *_RUN_LENGTHS,
    'adjacent_expected_exact_run_ge500', 'selected_expected_run_reaches_reference_sequence_edge',
    'flank_search_cap_flag', 'any_event_query_search_cap_flag',
    'any_ambiguous_query_not_exact_assessable', 'fcs_overlay_assessed', 'flank_fcs_overlap',
    'any_event_query_fcs_overlap', 'any_event_query_published_fcs_actions',
    'technical_boundary_warning', 'query_selection_scope', *_CONFIRMATIONS)
_SUMMARY_COUNTS = ('selected_event_count', 'selected_pair_count', 'assessed_event_count',
    'coordinate_supported_event_count', 'boundary_supported_event_count',
    'mate_geometry_evaluable_event_count', 'inward_mate_event_count', 'capped_event_count')
_EXTRA_SUMMARY = ('max_selected_expected_exact_run_bp',
    'max_selected_expected_exact_run_bp_with_coordinate_support',
    'max_expected_orientation_measured_run_bp', 'exact_ge500_event_count',
    'reference_edge_event_count', 'flank_fcs_overlap_event_count',
    'any_query_fcs_overlap_event_count', 'published_fcs_actions')
FAMILY_SUMMARY_FIELDS = ('family', 'declared_motif', 'role', *_SUMMARY_COUNTS, *_EXTRA_SUMMARY)
REFERENCE_SUMMARY_FIELDS = ('reference_accession', 'mode', *FAMILY_SUMMARY_FIELDS)
LIMITATIONS = [
    'Only a conditional saved-junction subset was reviewed; counts are not library-wide frequencies or independent biological replicates.',
    'Reference/method copies of the same event are deduplicated; mates and shared BioSamples are not independent biological replication.',
    'Retained alignments can be capped and gapped search is heuristic; absence of retained support does not establish telomere absence.',
    'Database size, the fixed E-value threshold, short-flank length and sequence composition affect reporting sensitivity across references.',
    'Full-read coordinate support, repeat–flank boundary crossing and mate orientation are separate evidence fields.',
    'An inward mate witness does not establish library insert compatibility; the reported nearest separation uses only retained placements.',
    'Inward-mate event counts are null if every geometry value is unassessed; zero requires an evaluated geometry denominator.',
    'A reference edge is not a physical chromosome end, and an exact run shorter than 500 bp is not a universal failure criterion.',
    'Reference taxonomy, primary-assembly labels and FCS overlaps do not confirm read origin, a clean junction, or candidate contamination.',
    'Any-query FCS overlap can be a distant alternative mate placement; it is distinct from flank overlap and is not a filter.',
    'Main TRM, novel TRM, physical termini and biological origin remain unconfirmed; native counts, classification and adequacy never change.',
]


class EvidenceError(ValueError):
    """Explicit evidence is malformed, stale, or bound to a different result."""


_IO = EvidenceIO(EvidenceError, 'reference localization evidence: ')


def _fail(message):
    raise EvidenceError('reference localization evidence: ' + message)


def _object(value, where):
    if not isinstance(value, dict):_fail(where + ' must be an object')
    return value


def _list(value, where):
    if not isinstance(value, list):_fail(where + ' must be a list')
    return value


def _text(value, where):
    if not isinstance(value, str) or not value.strip():_fail(where + ' must be a nonempty string')
    return value


def _integer(value, where, minimum=0):
    if type(value) is not int or value < minimum:_fail(where + ' must be an integer >= ' + str(minimum))
    return value


def _boolean(value, where, nullable=False):
    if not (type(value) is bool or (nullable and value is None)):_fail(where + ' must be boolean' + (' or null' if nullable else ''))
    return value


def _keys(value, required, where, optional=()):
    _object(value, where)
    if not required <= set(value) or set(value) - required - set(optional):
        _fail(where + ' has missing or unsupported fields')


def _sha(value, where):
    if not isinstance(value, str) or not _SHA.fullmatch(value):_fail(where + ' must be a lowercase SHA256 digest')
    return value


def _canonical(sequence):
    """Bounded metadata equivalence only; never a motif filter or search."""
    if not 1 <= len(sequence) <= MAX_MOTIF_LABEL_BP or not _DNA.fullmatch(sequence):
        _fail('motif labels must contain 1..256 uppercase DNA bases')
    for period in range(1,len(sequence)+1):
        if len(sequence)%period==0 and sequence==sequence[:period]*(len(sequence)//period):
            unit=sequence[:period];break
    reverse=unit.translate(str.maketrans('ACGT','TGCA'))[::-1]
    return min(s[i:]+s[:i] for s in (unit,reverse) for i in range(len(s)))


def _unique_object(pairs):
    return _IO.unique_object(pairs)


def _reject_constant(value):
    return _IO.reject_constant(value)


def evidence_payload_sha256(document):
    """Canonical UTF-8 SHA256 excluding only the top-level payload seal."""
    return _IO.payload_sha256(document)


def _read_bytes(path):
    return _IO.read_bytes(path, max_bytes=MAX_ARTIFACT_BYTES)


def _read_json(path):
    return _IO.decode_json(_read_bytes(path), path,
        object_pairs_hook=_unique_object, parse_constant=_reject_constant)


def _source_pins(document, base):
    pins = _list(document.get('source_artifacts'), 'source_artifacts')
    if not pins:_fail('at least one pinned source artifact is required')
    normalized, seen = [], set()
    for pin in pins:
        _keys(pin, {'path','sha256'}, 'source_artifact')
        path = Path(_text(pin['path'], 'source_artifact.path'))
        path = (base / path).absolute() if not path.is_absolute() else path
        lexical = str(path)
        if lexical in seen:_fail('duplicate source_artifact path: ' + lexical)
        seen.add(lexical)
        normalized.append({'path':lexical,'sha256':_sha(pin['sha256'],'source_artifact.sha256')})
    return normalized


def source_dependencies(path):
    """Return lexical sidecar/source paths for report-cache guards, not validation."""
    path = Path(path).absolute()
    document, _ = _read_json(path)
    return [str(path), *(p['path'] for p in _source_pins(document, path.parent))]


def _validate_event(event, sample_id, references):
    _keys(event, set(EVENT_FIELDS), 'event')
    for name in ('reference_accession','sample_id','family','declared_motif','event_id','pair_query_id'):
        _text(event[name], 'event.' + name)
    if event['sample_id'] != sample_id:_fail('event sample mismatch')
    if (sample_id,event['reference_accession']) not in references:_fail('event reference/sample linkage mismatch')
    if event['mode'] not in ('exact','gapped'):_fail('unsupported event mode')
    if event['role'] not in ('candidate','comparator'):_fail('unsupported event role')
    if event['state'] != 'completed_independently_audited':_fail('event is not completed and independently audited')
    if event['query_selection_scope'] != 'conditional_orientation_stratified_saved_events':_fail('unexpected event query selection scope')
    if _canonical(event['family']) != event['family']:_fail('family must be a primitive canonical DNA repeat')
    if _canonical(event['declared_motif']) != event['family']:_fail('declared_motif is not equivalent to its family')
    for name in _COUNTS:_integer(event[name], 'event.' + name)
    for name in _BOOLEANS:_boolean(event[name], 'event.' + name)
    for name in _NULLABLE_BOOLEANS:_boolean(event[name], 'event.' + name, nullable=True)
    for name in _RUN_LENGTHS:
        if event[name] is not None:_integer(event[name], 'event.' + name)
    span = event['smallest_retained_same_flank_locus_inward_five_prime_separation_bp']
    if span is not None:_integer(span, 'event inward separation', 1)
    for name in _CONFIRMATIONS:
        if event[name] is not None:_fail('biological confirmation fields must be null: ' + name)
    if not (event['retained_flank_hit_count'] >= event['retained_context_rows'] >= event['assessed_context_rows']
            >= event['full_read_coordinate_supported_assessed_context_rows']):_fail('inconsistent context count hierarchy')
    if event['assessed_flank_context_present'] != bool(event['assessed_context_rows']):_fail('assessed context flag/count mismatch')
    if not (event['raw_ungated_bridge_associations'] >= event['qualified_bridge_associations']
            >= event['qualified_bridge_adjacent_reference_continuous_associations']
            >= event['qualified_bridge_adjacent_reference_acgt_matching_associations']):_fail('inconsistent bridge count hierarchy')
    if event['qualified_bridge_associations'] and not event['full_read_coordinate_supported_assessed_context_rows']:
        _fail('qualified boundary bridge lacks assessed coordinate support')
    if event['flank_search_cap_flag'] and not event['any_event_query_search_cap_flag']:_fail('flank cap missing from any-query cap')
    if not event['assessed_context_rows'] and any(event[k] is not None for k in (*_RUN_LENGTHS,
            'adjacent_expected_exact_run_ge500','selected_expected_run_reaches_reference_sequence_edge')):
        _fail('unassessed repeat context requires null run/edge fields')
    measured = event['expected_orientation_measured_run_max_bp']
    if event['adjacent_expected_exact_run_ge500'] != (None if measured is None else measured >= 500):
        _fail('ge500 flag must agree with measured expected run or remain null')
    if measured is None and event['selected_expected_run_reaches_reference_sequence_edge'] is not None:
        _fail('unmeasured expected run requires null edge field')
    supported = event['selected_expected_run_max_bp_with_full_read_coordinate_support']
    selected = event['selected_expected_run_max_bp']
    if supported is not None and (not event['full_read_coordinate_supported_assessed_context_rows']
            or selected is None or supported > selected):_fail('coordinate-supported run is inconsistent')
    if selected is not None and (measured is None or selected > measured):_fail('selected run exceeds measured expected run')
    evaluable = event['same_flank_locus_mate_geometry_evaluable']
    witness = event['event_has_same_flank_locus_inward_mate_witness']
    if (not evaluable and witness is not None) or (evaluable and witness is None):_fail('mate geometry null/evaluable mismatch')
    if (witness is True) != (span is not None):_fail('inward witness/separation mismatch')
    # A broad pair witness is intentionally not required to support a boundary.
    if event['mode']=='gapped' and event['any_ambiguous_query_not_exact_assessable'] is not None:
        _fail('exact-query ambiguity field must be null for gapped mode')
    if event['mode']=='exact' and event['any_ambiguous_query_not_exact_assessable'] is None:
        _fail('exact-query ambiguity field must be evaluated for exact mode')
    actions = event['any_event_query_published_fcs_actions']
    if event['fcs_overlay_assessed'] is not True:
        if any(event[k] is not None for k in ('flank_fcs_overlap','any_event_query_fcs_overlap','any_event_query_published_fcs_actions')):
            _fail('unassessed FCS requires null overlaps/actions')
    else:
        _boolean(event['flank_fcs_overlap'],'flank_fcs_overlap')
        _boolean(event['any_event_query_fcs_overlap'],'any_event_query_fcs_overlap')
        _list(actions,'FCS actions')
        for action in actions:_text(action,'FCS action')
        if len(set(actions)) != len(actions):_fail('duplicate FCS actions')
        if bool(actions) != event['any_event_query_fcs_overlap']:_fail('FCS action/overlap mismatch')
        if event['flank_fcs_overlap'] and not event['any_event_query_fcs_overlap']:_fail('flank FCS overlap missing any-query overlap')


def _validate_sample(sample, references):
    _keys(sample, _SAMPLE_REQUIRED, 'sample', _IDENTITY)
    sid = _text(sample['sample_id'], 'sample_id')
    _sha(sample['source_result_sha256'], 'source_result_sha256')
    for name in ('scientific_name','biosample'):
        if name in sample:_text(sample[name], name)
    if 'tax_id' in sample:
        value=sample['tax_id']
        if not ((type(value) is int and value > 0) or (type(value) is str and value.isascii() and value.isdigit() and int(value)>0)):
            _fail('tax_id must be a positive numeric taxon ID')
    for name in ('selected_event_count','selected_pair_count'):_integer(sample[name],name)
    for note in _list(sample['review_notes'],'review_notes'):_text(note,'review note')
    seen, signatures, pairs = set(), {}, set()
    family_labels = {}
    for event in _list(sample['events'],'events'):
        _validate_event(event,sid,references)
        key=(event['reference_accession'],event['mode'],event['event_id'])
        if key in seen:_fail('duplicate reference/mode/event row')
        seen.add(key)
        signature=tuple(event[k] for k in ('family','pair_query_id','role','declared_motif'))
        eid=event['event_id']
        if eid in signatures and signatures[eid]!=signature:_fail('inconsistent event identity across references/methods')
        signatures[eid]=signature
        pairs.add(event['pair_query_id'])
        label=(event['declared_motif'],event['role'])
        if event['family'] in family_labels and family_labels[event['family']]!=label:_fail('inconsistent family motif/role labels')
        family_labels[event['family']]=label
    if len(signatures)!=sample['selected_event_count']:_fail('selected_event_count does not match unique saved events')
    if len(pairs)!=sample['selected_pair_count']:_fail('selected_pair_count does not match unique saved pairs')


class ValidatedEvidence(Mapping):
    """Validated sample map; returned values are copies, never the stored payload."""
    def __init__(self,samples,*,path,sha256,pins,references,payload_sha256):
        self._samples=deepcopy(samples)
        self.path=str(path);self.sha256=sha256;self.payload_sha256=payload_sha256
        self.source_paths=tuple(p['path'] for p in pins)
        self._pins=deepcopy(pins);self._references=deepcopy(references)
    def __getitem__(self,key):return deepcopy(self._samples[key])
    def __iter__(self):return iter(self._samples)
    def __len__(self):return len(self._samples)
    @property
    def provenance(self):
        return dict(sidecar_path=self.path,sidecar_sha256=self.sha256,
                    payload_sha256=self.payload_sha256,source_artifacts=deepcopy(self._pins))
    @property
    def references(self):return deepcopy(self._references)


def load_evidence(path):
    """Fail closed on explicit malformed evidence, seal or small-source mismatch."""
    path=Path(path).absolute()
    document,sidecar_sha=_read_json(path)
    _keys(document,_TOP,'document')
    if type(document['schema_version']) is not int or document['schema_version']!=SCHEMA_VERSION:_fail('unsupported schema_version')
    if document['scope']!=SCOPE:_fail('unexpected scope')
    seal=_sha(document['payload_sha256'],'payload_sha256')
    if seal!=evidence_payload_sha256(document):_fail('sidecar payload seal mismatch')
    refs=_list(document['references'],'references')
    if not refs:_fail('references cannot be empty')
    indexed={}
    for row in refs:
        _keys(row,_REFERENCE,'reference')
        for k in _REFERENCE:_text(row[k],'reference.'+k)
        if row['verification']!='fully_verified':_fail('reference must be fully_verified')
        key=(row['sample_id'],row['reference_accession'])
        if key in indexed:_fail('duplicate sample/reference row')
        indexed[key]=row
    samples={}
    for sample in _list(document['samples'],'samples'):
        _validate_sample(sample,indexed)
        if sample['sample_id'] in samples:_fail('duplicate sample_id')
        samples[sample['sample_id']]=sample
    if not samples:_fail('samples cannot be empty')
    if any(sid not in samples for sid,_ in indexed):_fail('reference links an undeclared sample')
    pins=_source_pins(document,path.parent)
    for pin in pins:
        if hashlib.sha256(_read_bytes(pin['path'])).hexdigest()!=pin['sha256']:
            _fail('source artifact SHA256 mismatch: '+pin['path'])
    return ValidatedEvidence(samples,path=path,sha256=sidecar_sha,pins=pins,references=refs,payload_sha256=seal)


def _summarize(rows):
    def count(field):return len({r['event_id'] for r in rows if r[field]})
    def optional_count(field):return count(field) if any(r[field] is not None for r in rows) else None
    def maximum(field):return max((r[field] for r in rows if r[field] is not None),default=None)
    return dict(selected_event_count=len({r['event_id'] for r in rows}),
        selected_pair_count=len({r['pair_query_id'] for r in rows}),
        assessed_event_count=count('assessed_context_rows'),
        coordinate_supported_event_count=count('full_read_coordinate_supported_assessed_context_rows'),
        boundary_supported_event_count=count('qualified_bridge_associations'),
        mate_geometry_evaluable_event_count=count('same_flank_locus_mate_geometry_evaluable'),
        inward_mate_event_count=optional_count('event_has_same_flank_locus_inward_mate_witness'),
        capped_event_count=count('any_event_query_search_cap_flag'),
        max_selected_expected_exact_run_bp=maximum('selected_expected_run_max_bp'),
        max_selected_expected_exact_run_bp_with_coordinate_support=maximum('selected_expected_run_max_bp_with_full_read_coordinate_support'),
        max_expected_orientation_measured_run_bp=maximum('expected_orientation_measured_run_max_bp'),
        exact_ge500_event_count=optional_count('adjacent_expected_exact_run_ge500'),
        reference_edge_event_count=optional_count('selected_expected_run_reaches_reference_sequence_edge'),
        flank_fcs_overlap_event_count=optional_count('flank_fcs_overlap'),
        any_query_fcs_overlap_event_count=optional_count('any_event_query_fcs_overlap'),
        published_fcs_actions=sorted({a for r in rows for a in r['any_event_query_published_fcs_actions'] or []})
            if any(r['any_event_query_published_fcs_actions'] is not None for r in rows) else None)


def build_reference_localization_evidence(sample_id,source_report,evidence=None):
    """Build a diagnostic bound to the native result, without mutating its report.

    Accept an assessed report with provenance.original_discovery_sha256, or a
    minimal sample_id/source_result_sha256 binding. Candidate membership is
    checked whenever candidates are present; optional identities bind exactly.
    """
    _text(sample_id,'requested sample_id')
    if evidence is not None and not isinstance(evidence,ValidatedEvidence):_fail('call load_evidence before building a diagnostic')
    result=dict(schema_version=SCHEMA_VERSION,status='unassessed',reason='evidence_not_supplied',
        display_label='참조 위치 근거 미평가',scope=SCOPE,
        summary='선정된 반복–flank 경계의 검증된 참조 위치 진단이 없습니다.',
        recommendation='저장된 경계·mate 자료와 적합한 참조의 검증된 위치 근거를 준비하세요.',
        sample_id=sample_id,**{k:None for k in (*_SUMMARY_COUNTS,'reference_count')},
        family_summaries=[],reference_summaries=[],events=[],references=[],review_notes=[],
        limitations=list(LIMITATIONS),provenance=evidence.provenance if evidence is not None else None,
        native_counts_changed=False,classification_changed=False,adequacy_changed=False,
        **dict.fromkeys((*_CONFIRMATIONS,'library_insert_compatibility')))
    if evidence is None:return result
    if sample_id not in evidence:
        result['reason']='sample_not_in_evidence';return result
    sample=evidence[sample_id]
    if not isinstance(source_report,dict) or source_report.get('sample_id')!=sample_id:_fail('current report sample_id mismatch')
    provenance=source_report.get('provenance',{})
    _object(provenance,'current report provenance')
    native_sha=provenance.get('original_discovery_sha256',source_report.get('source_result_sha256'))
    _sha(native_sha,'current report native result SHA256')
    if native_sha!=sample['source_result_sha256']:_fail('native result SHA256 mismatch for '+sample_id)
    families={e['family'] for e in sample['events']}
    if 'candidates' in source_report:
        candidates=_list(source_report['candidates'],'current report candidates')
        for c in candidates:_object(c,'current report candidate')
        if not families <= {c.get('canonical_family') for c in candidates if isinstance(c.get('canonical_family'),str)}:
            _fail('evidence family is absent from current report candidates for '+sample_id)
    for field in _IDENTITY:
        if field not in sample:continue
        metadata=source_report.get('metadata',{})
        current=(metadata.get('sample_accession') if isinstance(metadata,dict) else None) if field=='biosample' else source_report.get(field)
        if current is None or str(current)!=str(sample[field]):_fail('current report '+field+' mismatch for '+sample_id)
    result['provenance'].update(source_result_sha256=native_sha,
        sample_identity={k:sample[k] for k in ('sample_id',*_IDENTITY) if k in sample})
    if not sample['selected_event_count']:
        result['reason']='no_selected_junction_events';return result
    rows=sample['events'];groups=defaultdict(list);refgroups=defaultdict(list)
    for row in rows:
        groups[row['family']].append(row)
        refgroups[(row['reference_accession'],row['mode'],row['family'])].append(row)
    def labelled(rs):return dict(family=rs[0]['family'],declared_motif=rs[0]['declared_motif'],role=rs[0]['role'],**_summarize(rs))
    result['family_summaries']=[labelled(groups[k]) for k in sorted(groups)]
    result['reference_summaries']=[dict(reference_accession=k[0],mode=k[1],**labelled(refgroups[k])) for k in sorted(refgroups)]
    result['events']=rows
    result['references']=[r for r in evidence.references if r['sample_id']==sample_id]
    result['reference_count']=len(result['references'])
    result['review_notes']=sample['review_notes']
    all_counts=_summarize(rows)
    result.update({k:all_counts[k] for k in _SUMMARY_COUNTS})
    supported=result['boundary_supported_event_count']
    result.update(reason=None,
        status='reviewed_subset_boundary_support' if supported else 'reviewed_subset_no_boundary_support',
        display_label='선정 경계의 참조 연결 근거 있음' if supported else '선정 경계에서 참조 연결 근거 미확보',
        summary=f"선정 경계 {result['selected_event_count']}건 중 {supported}건에서 평가된 flank와 연결되는 반복–flank 경계 정렬 근거를 확인했습니다. 물리적 말단·생물학적 기원은 미확정입니다.",
        recommendation='참조·구획·경계·mate 방향 및 FCS 경고를 구분해 검토하고, 적합한 긴 서열이나 독립 표본으로 검증하세요. 미지지 결과를 텔로미어 부재로 해석하지 마세요.')
    return result
