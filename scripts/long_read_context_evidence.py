"""Bounded saved long-read annotations; no discovery, trimming or classification.

Schema 1 is deliberately separate from paired short-read technical evidence.
All coordinates are zero-based, half-open. Counts describe a selected read
subset and are derived from validated rows, never a whole-library denominator.
"""
from __future__ import annotations

from copy import deepcopy
import csv
import hashlib
import io
import json
from pathlib import Path
import re
from urllib.parse import urlparse

from evidence_io import EvidenceIO

SCHEMA_VERSION = 1
MAX_ARTIFACT_BYTES = 64 * 1024 * 1024
MAX_TOTAL_SOURCE_BYTES = 256 * 1024 * 1024
MAX_MANIFEST_BYTES = 4 * 1024 * 1024
MAX_TABLE_BYTES = 16 * 1024 * 1024
THRESHOLDS = (12, 16, 20)
SUMMARY_FIELDS = (
    'selected_reads', 'contexts', 'exact_match_reads_ge12', 'exact_match_reads_ge16',
    'exact_match_reads_ge20', 'full_signature_anywhere', 'full_signature_at_read_end',
    'whole_tail_signature', 'partial_signature_ge20_at_read_end',
    'boundary_overlap_reads', 'structure_reviewed_reads', 'other_tract_support_reads',
    'distal_boundary_support_reads', 'record_end_support_reads',
)
FIELDS = tuple('long_read_context_' + key for key in ('status', 'source_run', *SUMMARY_FIELDS))
SAMPLE_FIELDS = FIELDS
READ_FIELDS = (
    'sample_id', 'source_run', 'read_key', 'read_id', 'context_id', 'tract_id',
    'read_length', 'tract_start0', 'tract_end0', 'query_outward_direction',
    'outward_tail_bp', 'native_tract_bp', 'HSP_outward_clip_bp',
    'outward_tail_maximum_exact_match_bp', 'outward_tail_matched_bp_ge12',
    'outward_tail_matched_bp_ge16', 'outward_tail_matched_bp_ge20',
    'outward_tail_full26_anywhere', 'outward_tail_full26_at_read_end',
    'terminal20_match_bp', 'boundary_crossing_exact12_hits',
    'boundary_crossing_max_native_overlap_bp',
)
CONTEXT_FIELDS = (
    'sample_id', 'source_run', 'context_id', 'target_record', 'reference_side',
    'reference_length', 'linked_read_count', 'selected_boundary_start0',
    'selected_boundary_end0', 'selected_boundary_same_hsp_count',
    'selected_boundary_within100bp_count', 'local_structure_reviewed',
    'structure_reviewed_read_count', 'distal_boundary_supported_read_count',
    'record_end_supported_read_count', 'other_native_tract_supported_read_count',
)
STRUCTURE_FIELDS = (
    'sample_id', 'source_run', 'read_key', 'context_id', 'alignment_id',
    'mapping_strand', 'query_length', 'target_start0', 'target_end0',
    'target_record_end', 'unmapped_reference_suffix_bp',
    'selected_native_tract_id', 'selected_native_start0', 'selected_native_end0',
    'distal_boundary_reached', 'outermost_array_reached',
    'reference_record_last_base_aligned', 'other_native_tract_count',
    'other_native_tract_distal_boundary_support_count',
    'selected_native_tract_distal_boundary_support',
)
LIMITATIONS = [
    'Saved selected long-read subset only; not whole-library prevalence or discovery sensitivity.',
    'The technical signature is externally documented and protocol-compatible; exact kit use is not established.',
    'Exact signature matches do not establish RNA origin, complete adapter removal, or original molecule ends.',
    'The 12/16/20-bp exact-match thresholds are nested and must not be added.',
    'Partial terminal signature counts cover matches of at least 20 bp; shorter terminal matches are not assessed.',
    'Selected reference repeat envelopes may contain interruptions; they are not necessarily continuous exact arrays.',
    'Additional tracts on the same read are support, not independent molecules or biological replicates.',
    'Reference record-end alignment and distal-array support do not confirm a physical telomere, main TRM or novel TRM.',
    'Unreviewed local structure counts are null, not zero; only the saved reviewed subset is assessed.',
    'Native counts, ranks, classification, adequacy and sequences are unchanged.',
]
_SHA = re.compile(r'^[0-9a-f]{64}$')
_DNA = re.compile(r'^[ACGT]+$')


class EvidenceError(ValueError):
    """Invalid, stale, inconsistent or mismatched saved evidence."""


_IO = EvidenceIO(EvidenceError, 'long-read context evidence: ')
_fail = _IO.fail
payload_sha256 = _IO.payload_sha256


def _object(value, label):
    if not isinstance(value, dict):
        _fail(label + ' must be an object')
    return value


def _text(value, label):
    if not isinstance(value, str) or not value.strip():
        _fail(label + ' must be nonempty text')
    return value


def _integer(value, label, minimum=0):
    if type(value) is not int or value < minimum:
        _fail(label + ' must be an integer >= ' + str(minimum))
    return value


def _sha(value, label):
    if not isinstance(value, str) or not _SHA.fullmatch(value):
        _fail(label + ' must be a lowercase SHA256 digest')
    return value


def _read_json(path):
    data = _IO.read_bytes(path, max_bytes=MAX_MANIFEST_BYTES)
    return _IO.decode_json(data, path, object_pairs_hook=_IO.unique_object,
                           parse_constant=_IO.reject_constant)


def _resolve(path, base):
    value = Path(_text(path, 'path'))
    return str((base / value).absolute() if not value.is_absolute() else value)


def _table(path, sources, fields):
    if path not in sources:
        _fail('table has no source pin: ' + path)
    data = sources[path]
    if len(data) > MAX_TABLE_BYTES:
        _fail('table exceeds size bound: ' + path)
    try:
        reader = csv.DictReader(io.StringIO(data.decode('utf-8-sig')), delimiter='\t')
        if (not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames))
                or not set(fields) <= set(reader.fieldnames)):
            _fail('table missing required or duplicate columns: ' + path)
        rows = list(reader)
    except (UnicodeError, csv.Error) as error:
        raise EvidenceError('long-read context evidence: malformed TSV ' + path) from error
    if len(rows) > 100000:
        _fail('table row bound exceeded')
    for row in rows:
        if None in row or any(value is None for value in row.values()):
            _fail('malformed TSV row: ' + path)
    return rows


def _count(row, field, nullable=False):
    value = row.get(field)
    if nullable and value in ('', 'null', None):
        row[field] = None
        return None
    if not isinstance(value, str) or not re.fullmatch(r'0|[1-9][0-9]*', value):
        _fail(field + ' must be a nonnegative TSV integer')
    row[field] = int(value)
    return row[field]


def _boolean(row, field):
    value = row.get(field)
    if value not in ('True', 'False', 'true', 'false'):
        _fail(field + ' must be a TSV boolean')
    row[field] = value.lower() == 'true'
    return row[field]


def _row_identity(row, sample):
    if row.get('sample_id') != sample['sample_id'] or row.get('source_run') != sample['source_run']:
        _fail('table sample/run identity mismatch')


def _validate_sample(sample, base, sources, populate_summary=False):
    _object(sample, 'sample')
    for key in ('sample_id', 'source_run', 'scientific_name', 'biosample'):
        _text(sample.get(key), key)
    tax = sample.get('tax_id')
    if isinstance(tax, bool) or not str(tax).isdigit() or int(tax) < 1:
        _fail('tax_id must be positive numeric taxonomy ID')
    if sample.get('original_native_sha256') is not None:
        _sha(sample['original_native_sha256'], 'original_native_sha256')
    signature = _object(sample.get('signature'), 'signature')
    _text(signature.get('signature_id'), 'signature_id')
    sequence = signature.get('sequence')
    if not isinstance(sequence, str) or not _DNA.fullmatch(sequence):
        _fail('signature sequence must be uppercase DNA')
    length = _integer(signature.get('length_bp'), 'signature length', 20)
    # Schema 1 field names retain the frozen 26-bp study column convention.
    if length != len(sequence) or length != 26:
        _fail('schema 1 requires the documented 26-bp signature length')
    parsed_url = urlparse(_text(signature.get('source_url'), 'signature source_url'))
    if parsed_url.scheme != 'https' or not parsed_url.netloc:
        _fail('signature requires an external HTTPS source URL')
    applicable = signature.get('applicable_runs')
    if (not isinstance(applicable, list) or not applicable
            or any(not isinstance(run, str) or not run for run in applicable)
            or len(applicable) != len(set(applicable)) or sample['source_run'] not in applicable):
        _fail('signature is not applicable to this source run')
    if signature.get('sequence_derived_from_reads', False) is not False:
        _fail('technical signature must not be derived from selected reads')
    if signature.get('known_signature', True) is not True:
        _fail('signature must be externally known')
    signature_path = _resolve(sample.get('signature_source'), base)
    if signature_path not in sources:
        _fail('signature source requires a bounded source pin')
    panel, _ = _IO.decode_json(sources[signature_path], signature_path,
                              object_pairs_hook=_IO.unique_object,
                              parse_constant=_IO.reject_constant)
    entries = panel.get('signatures')
    if not isinstance(entries, list):
        _fail('signature source panel must contain signatures')
    selected = [entry for entry in entries if isinstance(entry, dict)
                and entry.get('signature_id', entry.get('id')) == signature['signature_id']]
    if len(selected) != 1:
        _fail('signature identifier must match one frozen panel entry')
    frozen = selected[0]
    for field in ('sequence', 'length_bp', 'source_url'):
        if frozen.get(field) != signature[field]:
            _fail('signature differs from pinned source panel: ' + field)
    if frozen.get('known_signature') is not True or frozen.get('sequence_derived_from_reads') is not False:
        _fail('pinned signature panel requires externally known, non-read-derived sequence')
    reverse = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    for record in (signature, frozen):
        if 'reverse_complement' in record and record['reverse_complement'] != reverse:
            _fail('signature reverse complement is inconsistent')
    tables = {}
    for key, fields in (('reads_source', READ_FIELDS), ('contexts_source', CONTEXT_FIELDS),
                        ('structures_source', STRUCTURE_FIELDS)):
        if key == 'structures_source' and sample.get(key) is None:
            tables[key] = []
        else:
            path = _resolve(sample.get(key), base)
            tables[key] = _table(path, sources, fields)
    reads, contexts, structures = (tables[k] for k in ('reads_source', 'contexts_source', 'structures_source'))
    if not reads or not contexts:
        _fail('annotated sample requires nonempty reads and contexts')
    context_map = {}
    for row in contexts:
        _row_identity(row, sample)
        key = _text(row.get('context_id'), 'context_id')
        if key in context_map:
            _fail('duplicate context_id')
        context_map[key] = row
        _text(row.get('target_record'), 'target_record')
        if row.get('reference_side') not in ('left', 'right'):
            _fail('reference_side must be left or right')
        for field in ('reference_length', 'linked_read_count', 'selected_boundary_start0',
                      'selected_boundary_end0', 'selected_boundary_same_hsp_count',
                      'selected_boundary_within100bp_count'):
            _count(row, field)
        if row['linked_read_count'] < 1:
            _fail('context must contain at least one selected read')
        if not 0 <= row['selected_boundary_start0'] < row['selected_boundary_end0'] <= row['reference_length']:
            _fail('selected boundary coordinates outside reference')
        for field in ('selected_boundary_same_hsp_count', 'selected_boundary_within100bp_count'):
            if row[field] > row['linked_read_count']:
                _fail('selected-boundary count exceeds linked reads')
        reviewed = _boolean(row, 'local_structure_reviewed')
        for field in ('structure_reviewed_read_count', 'distal_boundary_supported_read_count',
                      'record_end_supported_read_count', 'other_native_tract_supported_read_count'):
            value = _count(row, field, nullable=not reviewed)
            if not reviewed and value is not None:
                _fail('unreviewed structure counts must be null')
            if reviewed and value > row['linked_read_count']:
                _fail('structure count exceeds linked reads')
        if reviewed and row['structure_reviewed_read_count'] < 1:
            _fail('reviewed local structure requires at least one reviewed read')
    read_map, read_ids, tract_ids = {}, set(), set()
    text_fields = {'sample_id', 'source_run', 'read_key', 'read_id', 'context_id', 'tract_id',
                   'query_outward_direction'}
    bool_fields = {'outward_tail_full26_anywhere', 'outward_tail_full26_at_read_end'}
    for row in reads:
        _row_identity(row, sample)
        for field in text_fields:
            _text(row.get(field), field)
        if row['context_id'] not in context_map:
            _fail('read context is undeclared')
        if (row['read_key'] in read_map or row['read_id'] in read_ids or row['tract_id'] in tract_ids):
            _fail('duplicate read or selected tract identity')
        if not row['read_id'].startswith(sample['source_run'] + '.'):
            _fail('read_id does not belong to source_run')
        read_map[row['read_key']] = row
        read_ids.add(row['read_id'])
        tract_ids.add(row['tract_id'])
        for field in READ_FIELDS:
            if field in bool_fields:
                _boolean(row, field)
            elif field not in text_fields:
                _count(row, field)
        start, end, read_length = row['tract_start0'], row['tract_end0'], row['read_length']
        if not 0 <= start < end <= read_length or end - start != row['native_tract_bp']:
            _fail('native tract coordinates/length disagree')
        direction = row['query_outward_direction']
        if direction not in ('increasing', 'decreasing'):
            _fail('query_outward_direction must be increasing or decreasing')
        tail = read_length - end if direction == 'increasing' else start
        if tail != row['outward_tail_bp'] or row['HSP_outward_clip_bp'] > read_length:
            _fail('outward tail or clipping coordinates disagree')
        maximum = row['outward_tail_maximum_exact_match_bp']
        if maximum > length or maximum > tail or (0 < maximum < min(THRESHOLDS)):
            _fail('maximum exact signature match outside valid lengths')
        unions = [row['outward_tail_matched_bp_ge' + str(t)] for t in THRESHOLDS]
        if not tail >= unions[0] >= unions[1] >= unions[2] >= 0:
            _fail('signature threshold unions are not nested or exceed tail')
        for threshold, count in zip(THRESHOLDS, unions):
            if bool(count) != (maximum >= threshold) or (count and count < maximum):
                _fail('signature threshold union disagrees with maximum exact match')
        if row['outward_tail_full26_anywhere'] != (maximum == length):
            _fail('full signature flag disagrees with exact match length')
        terminal = row['terminal20_match_bp']
        if terminal != 0 and not 20 <= terminal <= maximum:
            _fail('terminal signature match outside evaluated threshold/maximum')
        if row['outward_tail_full26_at_read_end'] != (terminal == length):
            _fail('full terminal signature flag disagrees with terminal match')
        crossings, overlap = row['boundary_crossing_exact12_hits'], row['boundary_crossing_max_native_overlap_bp']
        if bool(crossings) != bool(overlap) or overlap > min(length - 1, end - start):
            _fail('boundary crossing count/overlap disagree')
    structure_keys = set()
    for row in structures:
        _row_identity(row, sample)
        key = row.get('read_key')
        if key in structure_keys or key not in read_map:
            _fail('duplicate or unknown structure read')
        structure_keys.add(key)
        read = read_map[key]
        if row.get('context_id') != read['context_id']:
            _fail('structure/read context mismatch')
        context = context_map[row['context_id']]
        if not context['local_structure_reviewed']:
            _fail('structure row in unreviewed context')
        _text(row.get('alignment_id'), 'alignment_id')
        if row.get('mapping_strand') not in ('+', '-'):
            _fail('mapping_strand must be + or -')
        expected_direction = ('increasing' if (context['reference_side'] == 'right') == (row['mapping_strand'] == '+')
                              else 'decreasing')
        if read['query_outward_direction'] != expected_direction:
            _fail('query/target reverse-complement direction mismatch')
        for field in ('query_length', 'target_start0', 'target_end0', 'target_record_end',
                      'unmapped_reference_suffix_bp', 'selected_native_start0', 'selected_native_end0',
                      'other_native_tract_count', 'other_native_tract_distal_boundary_support_count'):
            _count(row, field)
        for field in ('distal_boundary_reached', 'outermost_array_reached',
                      'reference_record_last_base_aligned', 'selected_native_tract_distal_boundary_support'):
            _boolean(row, field)
        if (row['query_length'] != read['read_length']
                or row.get('selected_native_tract_id') != read['tract_id']
                or (row['selected_native_start0'], row['selected_native_end0']) != (read['tract_start0'], read['tract_end0'])):
            _fail('structure/read selected tract mismatch')
        if (row['target_record_end'] != context['reference_length']
                or not 0 <= row['target_start0'] < row['target_end0'] <= row['target_record_end']
                or row['unmapped_reference_suffix_bp'] != row['target_record_end'] - row['target_end0']):
            _fail('structure target coordinates disagree')
        if (row['other_native_tract_distal_boundary_support_count'] > row['other_native_tract_count']
                or (row['reference_record_last_base_aligned'] and row['target_end0'] != row['target_record_end'])):
            _fail('structure support is inconsistent with counts/coordinates')
        if row['outermost_array_reached'] and not row['distal_boundary_reached']:
            _fail('outermost-array support requires distal-array support')
        if row['reference_record_last_base_aligned'] and not row['outermost_array_reached']:
            _fail('last reference base requires outermost array support')
        if ((row['other_native_tract_distal_boundary_support_count']
             or row['selected_native_tract_distal_boundary_support']) and not row['distal_boundary_reached']):
            _fail('tract distal support requires full-read distal support')
    for context in contexts:
        group = [row for row in reads if row['context_id'] == context['context_id']]
        if context['linked_read_count'] != len(group):
            _fail('context linked_read_count does not match read rows')
        if context['local_structure_reviewed']:
            group = [row for row in structures if row['context_id'] == context['context_id']]
            expected = {
                'structure_reviewed_read_count': len(group),
                'distal_boundary_supported_read_count': sum(row['distal_boundary_reached'] for row in group),
                'record_end_supported_read_count': sum(row['reference_record_last_base_aligned'] for row in group),
                'other_native_tract_supported_read_count': sum(row['other_native_tract_distal_boundary_support_count'] > 0 for row in group),
            }
            for field, count in expected.items():
                if context[field] != count:
                    _fail('context ' + field + ' does not match structure rows')
    summary = {
        'selected_reads': len(reads), 'contexts': len(contexts),
        **{'exact_match_reads_ge' + str(t): sum(row['outward_tail_maximum_exact_match_bp'] >= t for row in reads)
           for t in THRESHOLDS},
        'full_signature_anywhere': sum(row['outward_tail_full26_anywhere'] for row in reads),
        'full_signature_at_read_end': sum(row['outward_tail_full26_at_read_end'] for row in reads),
        'whole_tail_signature': sum(row['outward_tail_full26_at_read_end'] and row['outward_tail_bp'] == length for row in reads),
        'partial_signature_ge20_at_read_end': sum(20 <= row['terminal20_match_bp'] < length for row in reads),
        'boundary_overlap_reads': sum(row['boundary_crossing_exact12_hits'] > 0 for row in reads),
        'structure_reviewed_reads': len(structures) if structures else None,
        'other_tract_support_reads': sum(row['other_native_tract_distal_boundary_support_count'] > 0 for row in structures) if structures else None,
        'distal_boundary_support_reads': sum(row['distal_boundary_reached'] for row in structures) if structures else None,
        'record_end_support_reads': sum(row['reference_record_last_base_aligned'] for row in structures) if structures else None,
    }
    if populate_summary:
        sample['summary'] = summary
    else:
        declared = _object(sample.get('summary'), 'summary')
        if set(declared) != set(SUMMARY_FIELDS):
            _fail('summary fields are incomplete or unknown')
        for key, value in summary.items():
            if value is not None:
                _integer(declared.get(key), 'summary.' + key)
            if declared.get(key) != value:
                _fail('declared summary disagrees with source rows: ' + key)
    return {'sample': deepcopy(sample), 'summary': summary, 'reads': reads,
            'contexts': contexts, 'structures': structures}


def _validate_manifest(document, path, manifest_sha, *, populate_summary=False):
    if (document.get('schema_version') != SCHEMA_VERSION
            or type(document.get('schema_version')) is not int
            or document.get('evidence_type') != 'long_read_context_evidence'):
        _fail('unsupported schema or evidence_type')
    pins = document.get('source_artifacts')
    if not isinstance(pins, list) or not 1 <= len(pins) <= 64:
        _fail('source_artifacts must contain 1..64 pinned small sources')
    base = path.parent
    sources, normalized, total = {}, [], 0
    for pin in pins:
        _object(pin, 'source_artifact')
        source_path = _resolve(pin.get('path'), base)
        if source_path in sources:
            _fail('duplicate source path')
        digest = _sha(pin.get('sha256'), 'source sha256')
        data = _IO.read_bytes(source_path, max_bytes=MAX_ARTIFACT_BYTES)
        total += len(data)
        if total > MAX_TOTAL_SOURCE_BYTES:
            _fail('total bounded source size exceeded')
        if hashlib.sha256(data).hexdigest() != digest:
            _fail('source SHA256 mismatch: ' + source_path)
        sources[source_path] = data
        normalized.append({'path': source_path, 'sha256': digest})
    samples = document.get('samples')
    if not isinstance(samples, list) or not samples:
        _fail('samples must be a nonempty list')
    result = {}
    for sample in samples:
        validated = _validate_sample(sample, base, sources, populate_summary=populate_summary)
        key = sample['sample_id']
        if key in result:
            _fail('duplicate sample_id')
        validated['provenance'] = {
            'manifest_path': str(path), 'manifest_sha256': manifest_sha,
            'payload_sha256': document.get('payload_sha256'), 'source_artifacts': deepcopy(normalized),
            'original_native_sha256': sample.get('original_native_sha256'),
            'new_raw_passes': 0, 'new_mapping_runs': 0,
        }
        result[key] = validated
    return result


class LoadedEvidence(dict):
    """Dictionary-compatible validated snapshot with a captured integrity seal."""

    def __init__(self, values, manifest_path):
        super().__init__(values)
        self._manifest_path = manifest_path
        self._validated_sha256 = _IO.payload_sha256({'validated_samples': dict(self)})

    def checked_snapshot(self):
        if _IO.payload_sha256({'validated_samples': dict(self)}) != self._validated_sha256:
            _fail('loaded evidence changed after validation')
        current = load_evidence(self._manifest_path)
        if current._validated_sha256 != self._validated_sha256:
            _fail('evidence sources changed after loading')
        return current


def load_evidence(path=None):
    """Fully validate current source bytes and return samples keyed by sample ID."""
    if path is None:
        return {}
    path = Path(path).absolute()
    document, digest = _read_json(path)
    if document.get('payload_sha256') != payload_sha256(document):
        _fail('manifest canonical payload SHA256 mismatch')
    return LoadedEvidence(_validate_manifest(document, path, digest), str(path))


def source_dependencies(path=None):
    """Revalidate sources even on a warm report cache before listing dependencies."""
    if path is None:
        return []
    loaded = load_evidence(path)
    first = next(iter(loaded.values()))
    return [first['provenance']['manifest_path'],
            *(pin['path'] for pin in first['provenance']['source_artifacts'])]


def create_manifest(path, samples, source_paths):
    """Write a sealed manifest from pinned tables, deriving all summary counts.

    ``samples`` contains identity, signature and the three source-table paths;
    ``source_paths`` also includes the frozen source artifacts used by a producer.
    The helper validates first and never invents native-result identity bindings.
    """
    path = Path(path).absolute()
    pins = []
    for item in source_paths:
        source = Path(item).absolute()
        data = _IO.read_bytes(source, max_bytes=MAX_ARTIFACT_BYTES)
        pins.append({'path': str(source), 'sha256': hashlib.sha256(data).hexdigest()})
    document = {'schema_version': SCHEMA_VERSION, 'evidence_type': 'long_read_context_evidence',
                'source_artifacts': pins, 'samples': deepcopy(samples)}
    _validate_manifest(document, path, None, populate_summary=True)
    document['payload_sha256'] = payload_sha256(document)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    return document


def build_evidence(sample, document=None, loaded=None):
    """Build a read-only reporting overlay, or a standalone saved-evidence report.

    A native document requires its own complete result, original hash and exact
    sample/taxon/BioSample identity. An unbound study sidecar is standalone only.
    """
    result = {
        'status': 'unassessed', 'reason': 'not_supplied' if not loaded else 'sample_not_in_evidence',
        'scope': 'selected_saved_long_read_contexts', 'source_run': None, 'identity': None,
        'signature': None, 'summary': dict.fromkeys(SUMMARY_FIELDS), 'reads': [], 'contexts': [],
        'structures': [], 'notes': [], 'limitations': list(LIMITATIONS), 'provenance': {},
        'changes_classification_or_adequacy': False, 'native_counts_changed': False,
        'sequence_trimmed': False, 'physical_terminal_confirmed': None,
        'main_trm_confirmed': None, 'novel_trm_confirmed': None,
    }
    if loaded:
        if not isinstance(loaded, LoadedEvidence):
            _fail('loaded evidence must originate from load_evidence')
        loaded = loaded.checked_snapshot()
    if loaded and sample in loaded:
        entry = loaded[sample]
        identity = entry['sample']
        if document is not None:
            if identity['sample_id'] != identity['source_run']:
                _fail('native sample_id must equal source_run; aliases cannot bind this evidence')
            _object(document, 'native report document')
            if document.get('analysis_complete') is not True:
                _fail('native report analysis_complete must be true')
            native = identity.get('original_native_sha256')
            if native is None:
                _fail('unbound evidence is standalone only; native original hash required')
            if _object(document.get('provenance'), 'native provenance').get('original_discovery_sha256') != native:
                _fail('native original discovery SHA256 mismatch')
            meta = document.get('metadata', {})
            if not isinstance(meta, dict):
                _fail('native metadata must be an object')
            current = {
                'sample_id': document.get('sample_id'),
                'scientific_name': document.get('scientific_name'),
                'tax_id': document.get('tax_id'),
                'biosample': meta.get('sample_accession'),
            }
            for key, value in current.items():
                if value is None or str(value) != str(identity[key]):
                    _fail('native ' + key + ' identity mismatch')
        result.update(
            status='annotated_saved_subset', reason='saved_rows_validated', source_run=identity['source_run'],
            identity={key: identity[key] for key in ('sample_id', 'source_run', 'scientific_name', 'tax_id', 'biosample')},
            signature=deepcopy(identity['signature']), summary=deepcopy(entry['summary']),
            reads=deepcopy(entry['reads']), contexts=deepcopy(entry['contexts']),
            structures=deepcopy(entry['structures']), provenance=deepcopy(entry['provenance']),
            notes=deepcopy(identity.get('notes', [])),
        )
    result['display'] = {'long_read_context_status': result['status'],
                         'long_read_context_source_run': result['source_run'],
                         **{'long_read_context_' + key: value for key, value in result['summary'].items()}}
    return result
