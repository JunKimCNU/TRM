"""Read sealed, bounded technical-boundary evidence without opening FASTQ files.

This is a reporting overlay, never a motif filter or a scanner. Schema 1 stores
selected junction events, three nested exact-match thresholds (12/16/20 bp), a
small reference panel, and SHA256 bindings to both source artifacts and the
native result used by the current report. ``evidence_payload_sha256`` supplies
the sidecar's canonical content seal; it is an integrity check, not a signature.

Unavailable full reads or mates must have null counts, tiers and events. A
complete sample has one event row per selected junction and two saved full
reads per distinct original pair. Counts never describe the whole library.
"""
from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re

SCHEMA_VERSION = 1
THRESHOLDS = (12, 16, 20)
SCOPE = 'selected_junction_subset'
MAX_ARTIFACT_BYTES = 64 * 1024 * 1024
_SHA = re.compile(r'^[0-9a-f]{64}$')
_FAMILY = re.compile(r'^[ACGT]+$')
_COUNT_FIELDS = ('selected_event_count', 'selected_pair_count', 'full_read_count')
_TIER_FIELDS = ('exact_match_associations', 'full_reads_with_matches',
                'pairs_with_matches', 'events_with_flank_overlap')
LIMITATIONS = [
    'Only the selected saved junction subset was assessed; counts are not whole-library fractions.',
    'No exact hit to this limited panel does not establish adapter-free reads or genomic origin.',
    'Reference aliases, reverse complements and overlapping matches are associations, not independent molecules.',
    'Thresholds 12, 16 and 20 bp are nested, not additive; a 19-bp match remains visible at 12 and 16 bp.',
    'Read-through-compatible pair geometry does not establish technical origin or physical telomere location.',
    'Remaining flank bases are not established as genomic, unique or terminal.',
    'This evidence never changes native motif counts, classification, adequacy or sequence.',
]


class EvidenceError(ValueError):
    """An explicitly supplied evidence input is invalid, stale or mismatched."""


def _fail(message):
    raise EvidenceError('technical boundary evidence: ' + message)


def _object(value, where):
    if not isinstance(value, dict):
        _fail(where + ' must be an object')
    return value


def _list(value, where):
    if not isinstance(value, list):
        _fail(where + ' must be a list')
    return value


def _text(value, where):
    if not isinstance(value, str) or not value.strip():
        _fail(where + ' must be a nonempty string')
    return value


def _integer(value, where, minimum=0):
    # bool is an int subclass, but must not become a count.
    if type(value) is not int or value < minimum:
        _fail(where + ' must be an integer >= ' + str(minimum))
    return value


def _boolean(value, where):
    if type(value) is not bool:
        _fail(where + ' must be boolean')
    return value


def _sha(value, where):
    if not isinstance(value, str) or not _SHA.fullmatch(value):
        _fail(where + ' must be a lowercase SHA256 digest')
    return value


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            _fail('duplicate JSON key ' + repr(key))
        result[key] = value
    return result


def _reject_constant(value):
    _fail('non-finite JSON number ' + value)


def evidence_payload_sha256(document):
    """Return the canonical seal excluding only the top-level seal itself."""
    _object(document, 'document')
    value = {key: item for key, item in document.items() if key != 'payload_sha256'}
    try:
        data = json.dumps(value, sort_keys=True, ensure_ascii=False,
                          separators=(',', ':'), allow_nan=False).encode('utf-8')
    except (TypeError, ValueError) as error:
        raise EvidenceError('technical boundary evidence: non-JSON payload') from error
    return hashlib.sha256(data).hexdigest()


def _read_bytes(path):
    """Bounded stable read, including detection of a symlink target change."""
    path = Path(path).absolute()
    try:
        resolved = path.resolve(strict=True)
        before = path.stat()
        if not path.is_file() or before.st_size > MAX_ARTIFACT_BYTES:
            _fail('source must be a small regular file: ' + str(path))
        with path.open('rb') as handle:
            data = handle.read(MAX_ARTIFACT_BYTES + 1)
        after = path.stat()
    except OSError as error:
        raise EvidenceError('technical boundary evidence: cannot read ' + str(path)) from error
    identity = lambda stat: (stat.st_dev, stat.st_ino, stat.st_size,
                              stat.st_mtime_ns, stat.st_ctime_ns, stat.st_mode)
    if (len(data) > MAX_ARTIFACT_BYTES or identity(before) != identity(after)
            or resolved != path.resolve(strict=True)):
        _fail('source changed while reading: ' + str(path))
    return data


def _read_json(path):
    data = _read_bytes(path)
    try:
        document = json.loads(data, object_pairs_hook=_unique_object,
                              parse_constant=_reject_constant)
    except (ValueError, UnicodeError) as error:
        raise EvidenceError('technical boundary evidence: malformed JSON in ' + str(path)) from error
    return _object(document, 'document'), hashlib.sha256(data).hexdigest()


def _source_pins(document, base):
    pins = _list(document.get('source_artifacts'), 'source_artifacts')
    if not pins:
        _fail('at least one pinned source artifact is required')
    normalized, paths = [], set()
    for pin in pins:
        _object(pin, 'source_artifact')
        path = Path(_text(pin.get('path'), 'source_artifact.path'))
        path = (base / path).absolute() if not path.is_absolute() else path
        # Keep the lexical path for the caller's symlink/race guard.
        lexical = str(path)
        if lexical in paths:
            _fail('duplicate source_artifact path: ' + lexical)
        paths.add(lexical)
        normalized.append({'path': lexical, 'sha256': _sha(pin.get('sha256'), 'source_artifact.sha256')})
    return normalized


def source_dependencies(path):
    """List sidecar and pinned source paths for a report cache's own guards.

    This only parses the sidecar. ``load_evidence`` performs full validation;
    callers must not consider this dependency list validation of the evidence.
    """
    path = Path(path).absolute()
    document, _ = _read_json(path)
    return [str(path), *(pin['path'] for pin in _source_pins(document, path.parent))]


def _threshold_map(value, where, maximum):
    _object(value, where)
    if set(value) != {str(t) for t in THRESHOLDS}:
        _fail(where + ' must contain exactly 12, 16 and 20')
    result = {}
    for threshold in THRESHOLDS:
        count = _integer(value[str(threshold)], where + '.' + str(threshold))
        if count > maximum:
            _fail(where + ' exceeds flank length')
        result[threshold] = count
    return result


def _validate_sample(sample, reference_lengths):
    _object(sample, 'sample')
    sample_id = _text(sample.get('sample_id'), 'sample_id')
    _sha(sample.get('source_result_sha256'), sample_id + '.source_result_sha256')
    for identity_field in ('scientific_name', 'biosample'):
        if identity_field in sample:
            _text(sample[identity_field], sample_id + '.' + identity_field)
    if 'tax_id' in sample and (isinstance(sample['tax_id'], bool)
                              or not str(sample['tax_id']).isdigit()):
        _fail(sample_id + '.tax_id must be a numeric taxon ID')
    families = _list(sample.get('families'), sample_id + '.families')
    if not families or any(not isinstance(f, str) or not _FAMILY.fullmatch(f) for f in families):
        _fail(sample_id + '.families must contain DNA repeat families')
    if len(set(families)) != len(families):
        _fail(sample_id + '.families contains duplicates')
    availability = _object(sample.get('availability'), sample_id + '.availability')
    full = _boolean(availability.get('full_reads'), 'availability.full_reads')
    mates = _boolean(availability.get('mates'), 'availability.mates')
    if not full or not mates:
        for field in (*_COUNT_FIELDS, 'tiers', 'events'):
            if field not in sample or sample[field] is not None:
                _fail(sample_id + ': unavailable evidence requires null ' + field)
        return
    counts = {field: _integer(sample.get(field), sample_id + '.' + field) for field in _COUNT_FIELDS}
    if counts['full_read_count'] != counts['selected_pair_count'] * 2:
        _fail(sample_id + ': full_read_count must be twice selected_pair_count')
    events = _list(sample.get('events'), sample_id + '.events')
    if len(events) != counts['selected_event_count']:
        _fail(sample_id + ': selected_event_count does not match events')
    identifiers, pair_ids = set(), set()
    overlap_event_ids = {threshold: set() for threshold in THRESHOLDS}
    overlap_pair_ids = {threshold: set() for threshold in THRESHOLDS}
    overlap_read_ids = {threshold: set() for threshold in THRESHOLDS}
    for event in events:
        _object(event, 'event')
        event_id = _text(event.get('event_id'), 'event_id')
        if event_id in identifiers:
            _fail(sample_id + ': duplicate event_id ' + event_id)
        identifiers.add(event_id)
        pair = _integer(event.get('original_pair_index_1based'), 'original_pair_index_1based', 1)
        mate = _integer(event.get('mate'), 'mate', 1)
        if mate not in (1, 2):
            _fail('mate must be 1 or 2')
        pair_ids.add(pair)
        if event.get('family') not in families:
            _fail(sample_id + ': event family is absent from declared families')
        flank = _integer(event.get('flank_bp'), 'flank_bp', 1)
        maximum = _integer(event.get('maximum_exact_match_bp_overlapping_flank'),
                           'maximum_exact_match_bp_overlapping_flank')
        if maximum and (maximum < min(THRESHOLDS) or maximum > max(reference_lengths)):
            _fail(sample_id + ': maximum match is outside panel length/minimum threshold')
        union = _threshold_map(event.get('union_flank_overlap_by_threshold'),
                               'union_flank_overlap_by_threshold', flank)
        adjacent = _threshold_map(event.get('adjacent_unmatched_flank_by_threshold'),
                                  'adjacent_unmatched_flank_by_threshold', flank)
        paired = _boolean(event.get('paired_readthrough_compatible'), 'paired_readthrough_compatible')
        for threshold in THRESHOLDS:
            overlap = union[threshold]
            if bool(overlap) != (maximum >= threshold):
                _fail(sample_id + ': maximum match and threshold flank overlap disagree')
            if adjacent[threshold] > flank - overlap:
                _fail(sample_id + ': adjacent unmatched bases exceed remaining flank')
            if not overlap and adjacent[threshold] != flank:
                _fail(sample_id + ': no-overlap event must retain its complete flank')
            if overlap:
                overlap_event_ids[threshold].add(event_id)
                overlap_pair_ids[threshold].add(pair)
                overlap_read_ids[threshold].add((pair, mate))
        if not (union[12] >= union[16] >= union[20] and adjacent[12] <= adjacent[16] <= adjacent[20]):
            _fail(sample_id + ': threshold flank quantities are not nested')
        if paired and not union[12]:
            _fail(sample_id + ': paired read-through flag requires a matched flank')
    if len(pair_ids) != counts['selected_pair_count']:
        _fail(sample_id + ': selected_pair_count does not match distinct event pairs')
    tiers = _list(sample.get('tiers'), sample_id + '.tiers')
    if len(tiers) != len(THRESHOLDS):
        _fail(sample_id + ': exactly three threshold tiers are required')
    by_threshold = {}
    for tier in tiers:
        _object(tier, 'tier')
        threshold = _integer(tier.get('threshold_bp'), 'threshold_bp')
        if threshold not in THRESHOLDS or threshold in by_threshold:
            _fail(sample_id + ': duplicate or unsupported threshold tier')
        by_threshold[threshold] = tier
        values = {field: _integer(tier.get(field), field) for field in _TIER_FIELDS}
        associations, reads, pairs, event_count = (values[field] for field in _TIER_FIELDS)
        if (pairs > counts['selected_pair_count'] or reads > counts['full_read_count']
                or not pairs <= reads <= 2 * pairs or associations < reads
                or bool(associations) != bool(reads)):
            _fail(sample_id + ': inconsistent tier association/read/pair counts')
        if event_count != len(overlap_event_ids[threshold]):
            _fail(sample_id + ': tier flank event count differs from event rows')
        if reads < len(overlap_read_ids[threshold]) or pairs < len(overlap_pair_ids[threshold]):
            _fail(sample_id + ': tier matches omit a flank-overlap read/pair')
        if associations and not any(length >= threshold for length in reference_lengths):
            _fail(sample_id + ': no reference is eligible for a populated threshold')
    for field in _TIER_FIELDS:
        if not by_threshold[12][field] >= by_threshold[16][field] >= by_threshold[20][field]:
            _fail(sample_id + ': threshold counts are not nested: ' + field)


class ValidatedEvidence(Mapping):
    """Validated sample map; access returns copies so caller edits cannot persist."""
    def __init__(self, samples, *, path, sha256, pins, panel, payload_sha256):
        self._samples = deepcopy(samples)
        self.path = str(path)
        self.sha256 = sha256
        self.payload_sha256 = payload_sha256
        self.source_paths = tuple(pin['path'] for pin in pins)
        self._pins = deepcopy(pins)
        self._panel = deepcopy(panel)

    def __getitem__(self, key):
        return deepcopy(self._samples[key])

    def __iter__(self):
        return iter(self._samples)

    def __len__(self):
        return len(self._samples)

    @property
    def provenance(self):
        return {'sidecar_path': self.path, 'sidecar_sha256': self.sha256,
                'payload_sha256': self.payload_sha256,
                'source_artifacts': deepcopy(self._pins)}

    @property
    def reference_panel(self):
        return deepcopy(self._panel)


def load_evidence(path):
    """Validate a sidecar and all small source pins, or raise EvidenceError."""
    path = Path(path).absolute()
    document, sidecar_sha = _read_json(path)
    if type(document.get('schema_version')) is not int or document['schema_version'] != SCHEMA_VERSION:
        _fail('unsupported schema_version')
    if document.get('scope') != SCOPE:
        _fail('scope must be selected_junction_subset')
    if document.get('thresholds_bp') != list(THRESHOLDS) or any(type(t) is not int for t in document['thresholds_bp']):
        _fail('thresholds_bp must be [12, 16, 20]')
    seal = _sha(document.get('payload_sha256'), 'payload_sha256')
    if seal != evidence_payload_sha256(document):
        _fail('sidecar payload seal mismatch')
    panel = _list(document.get('reference_panel'), 'reference_panel')
    if not panel:
        _fail('reference_panel cannot be empty')
    reference_ids, lengths = set(), []
    for reference in panel:
        _object(reference, 'reference')
        reference_id = _text(reference.get('reference_id'), 'reference_id')
        if reference_id in reference_ids:
            _fail('duplicate reference_id ' + reference_id)
        reference_ids.add(reference_id)
        lengths.append(_integer(reference.get('length_bp'), 'length_bp', 1))
    samples = {}
    for sample in _list(document.get('samples'), 'samples'):
        _validate_sample(sample, lengths)
        sample_id = sample['sample_id']
        if sample_id in samples:
            _fail('duplicate sample_id ' + sample_id)
        samples[sample_id] = sample
    pins = _source_pins(document, path.parent)
    for pin in pins:
        if hashlib.sha256(_read_bytes(pin['path'])).hexdigest() != pin['sha256']:
            _fail('source artifact SHA256 mismatch: ' + pin['path'])
    return ValidatedEvidence(samples, path=path, sha256=sidecar_sha, pins=pins,
                             panel=panel, payload_sha256=seal)


def build_technical_boundary_evidence(sample_id, source_report, evidence=None):
    """Return a conservative diagnostic bound to the current native result SHA.

    ``source_report`` is an assessed sample document, with ``sample_id`` and
    ``provenance.original_discovery_sha256``. A minimal explicit binding with
    ``sample_id`` and ``source_result_sha256`` is also accepted. It is not mutated.
    A supplied raw dict is rejected: evidence must first pass ``load_evidence``.
    """
    _text(sample_id, 'requested sample_id')
    if evidence is not None and not isinstance(evidence, ValidatedEvidence):
        _fail('call load_evidence before building a diagnostic')
    result = {
        'schema_version': SCHEMA_VERSION,
        'status': 'unassessed', 'reason': 'evidence_not_supplied',
        'display_label': '기술서열 경계 미평가', 'scope': SCOPE,
        'summary': '저장된 full-read·mate 기반 기술서열 경계 진단이 없습니다.',
        'recommendation': '선정 경계의 full read와 mate를 확보해 제한된 기술서열 패널과 비교하세요.',
        'sample_id': sample_id, 'families': None,
        'selected_event_count': None, 'selected_pair_count': None,
        'full_read_count': None, 'flagged_event_count': None,
        'paired_readthrough_compatible_event_count': None,
        'events': [], 'tiers': [], 'thresholds_bp': list(THRESHOLDS),
        'reference_panel': None,
        'provenance': evidence.provenance if evidence is not None else None,
        'native_counts_changed': False, 'classification_changed': False,
        'adequacy_changed': False, 'sequence_trimmed': False,
        'technical_origin_confirmed': None, 'physical_terminal_confirmed': None,
        'novel_trm_confirmed': None, 'limitations': list(LIMITATIONS),
    }
    if evidence is None:
        return result
    if sample_id not in evidence:
        result['reason'] = 'sample_not_in_evidence'
        return result
    sample = evidence[sample_id]
    if not isinstance(source_report, dict) or source_report.get('sample_id') != sample_id:
        _fail('current report sample_id does not match ' + sample_id)
    provenance = source_report.get('provenance', {})
    if not isinstance(provenance, dict):
        _fail('current report provenance must be an object')
    source_sha = provenance.get('original_discovery_sha256', source_report.get('source_result_sha256'))
    _sha(source_sha, 'current report native result SHA256')
    if source_sha != sample['source_result_sha256']:
        _fail('native result SHA256 mismatch for ' + sample_id)
    if 'candidates' in source_report:
        candidates = source_report['candidates']
        if not isinstance(candidates, list) or any(not isinstance(row, dict) for row in candidates):
            _fail('current report candidates must be a list of objects')
        current_families = {row.get('canonical_family') for row in candidates}
        if not set(sample['families']) <= current_families:
            _fail('evidence family is absent from current report candidates for ' + sample_id)
    for field in ('scientific_name', 'tax_id', 'biosample'):
        if field not in sample:
            continue
        if field == 'biosample':
            metadata = source_report.get('metadata', {})
            current_value = metadata.get('sample_accession') if isinstance(metadata, dict) else None
        else:
            current_value = source_report.get(field)
        if current_value is None or str(current_value) != str(sample[field]):
            _fail('current report ' + field + ' mismatch for ' + sample_id)
    result['provenance']['source_result_sha256'] = source_sha
    result['provenance']['sample_identity'] = {
        key: sample[key] for key in ('sample_id', 'scientific_name', 'tax_id', 'biosample') if key in sample}
    result['families'] = sample['families']
    result['reference_panel'] = [
        {**reference, 'eligible_thresholds_bp': [t for t in THRESHOLDS if t <= reference['length_bp']]}
        for reference in evidence.reference_panel]
    availability = sample['availability']
    result['availability'] = deepcopy(availability)
    if not availability['full_reads'] or not availability['mates']:
        result['reason'] = ('full_reads_unavailable' if not availability['full_reads'] else 'mates_unavailable')
        return result
    if not sample['selected_event_count']:
        result['reason'] = 'no_selected_junction_events'
        return result
    for field in _COUNT_FIELDS:
        result[field] = sample[field]
    result['events'] = sample['events']
    result['tiers'] = sorted(sample['tiers'], key=lambda tier: tier['threshold_bp'])
    result['flagged_event_count'] = sum(bool(event['union_flank_overlap_by_threshold']['12'])
                                        for event in sample['events'])
    result['paired_readthrough_compatible_event_count'] = sum(event['paired_readthrough_compatible']
                                                               for event in sample['events'])
    result['reason'] = None
    n = result['selected_event_count']
    flagged = result['flagged_event_count']
    if flagged:
        result.update(
            status='technical_boundary_warning', display_label='선정 경계의 기술서열 일치 경고',
            summary=f'선정 경계 {n}건 중 {flagged}건에서 기술서열과 12 bp 이상 정확히 일치하는 구간이 비반복 flank와 겹칩니다.',
            recommendation='해당 flank 전체를 독립적인 genomic junction 근거로 사용하지 말고, 기술서열 구간과 남은 경계를 분리해 위치·독립 자료를 검증하세요. Motif 자체는 제외하지 않습니다.')
    elif result['tiers'][0]['exact_match_associations']:
        result.update(
            status='technical_reference_match_without_flank_overlap', display_label='Full read 기술서열 일치·flank 겹침 없음',
            summary=f'선정 경계 {n}건의 full read 또는 mate에 기술서열 일치가 있으나 flank와 겹치지 않습니다.',
            recommendation='일치 위치를 검토하고 genomic 위치 검증을 계속하세요. 경계의 생물학적 기원은 아직 확정되지 않았습니다.')
    else:
        result.update(
            status='no_exact_flank_overlap_in_selected_subset', display_label='선정 경계·제한 패널에서 정확 일치 없음',
            summary=f'선정 경계 {n}건의 저장된 full read·mate에서 제한 패널과 12 bp 이상 정확한 일치를 찾지 못했습니다.',
            recommendation='다른 기술서열·불일치·전처리 영향은 배제되지 않습니다. Assembly 위치 및 독립 자료 검증을 계속하세요.')
    return result
