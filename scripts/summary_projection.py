"""Small report-produced projections, bound to complete assessment bytes.

Consumers still stream-hash the full source. They avoid loading/decoding its
candidate array. Independent recalculation remains available via --full-json.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path

SCHEMA = 1
MANIFEST = 'summary_projections.tsv'
FIELDS = ['sample_id', 'source_sha256', 'source_bytes', 'projection_sha256', 'projection_bytes']


def file_stamp(path):
    value = Path(path).stat()
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)


def stream_digest(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_projection(document, saved_document, public_document, output, *, summarize_sample, version):
    before = file_stamp(saved_document)
    digest = stream_digest(saved_document)
    if file_stamp(saved_document) != before:
        raise ValueError('Assessment changed while preparing summary projection')
    row, expected, queue = summarize_sample(document, public_document, digest)
    projection = {'schema_version': SCHEMA, 'program': 'trm-summary-projection',
                  'summarizer_version': version, 'sample_id': document['sample_id'],
                  'source_sha256': digest, 'source_bytes': before[2],
                  'row': row, 'expected': expected, 'queue': queue,
                  'scope': 'Producer-computed projection; source bytes require verification. Not an independent recalculation.'}
    Path(output).write_text(json.dumps(projection, ensure_ascii=False, allow_nan=False,
                                     separators=(',', ':')) + '\n', encoding='utf-8')


def write_manifest(receipts, output):
    output = Path(output)
    with (output / MANIFEST).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for receipt in receipts:
            sample = receipt['sample_id']
            source = receipt['artifacts'][str((output / 'samples' / (sample + '.json')).resolve())]
            projection = receipt['artifacts'][str((output / 'summaries' / (sample + '.json')).resolve())]
            writer.writerow({'sample_id': sample, 'source_sha256': source['sha256'],
                             'source_bytes': source['size_bytes'], 'projection_sha256': projection['sha256'],
                             'projection_bytes': projection['size_bytes']})


def manifest_entries(controls, paths):
    payload = controls.get(MANIFEST)
    metrics = json.loads(controls['metrics.json']) if controls.get('metrics.json') else {}
    expected_hash = (metrics.get('files') or {}).get(MANIFEST)
    if payload is None and expected_hash is None:
        return {}
    if payload is None or expected_hash != hashlib.sha256(payload).hexdigest():
        raise ValueError('Summary projection manifest checksum differs from completed metrics')
    reader = csv.DictReader(io.StringIO(payload.decode('utf-8')), delimiter='\t')
    if reader.fieldnames != FIELDS:
        raise ValueError('Invalid summary projection manifest columns')
    entries = {}
    for row in reader:
        sample = row['sample_id']
        if sample in entries:
            raise ValueError('Duplicate sample in summary projection manifest')
        for key in ('source_bytes', 'projection_bytes'):
            row[key] = int(row[key])
            if row[key] < 0:
                raise ValueError('Negative summary projection size')
        entries[sample] = row
    if set(entries) != {path.stem for path in paths}:
        raise ValueError('Summary projection cohort differs from samples.tsv')
    return entries


def read_projection(source, entry, *, version, row_fields, expected_fields, queue_fields):
    source = Path(source)
    projection_path = source.parent.parent / 'summaries' / source.name
    before = file_stamp(source)
    projection_before = file_stamp(projection_path)
    payload = projection_path.read_bytes()
    if (len(payload) != entry['projection_bytes']
            or hashlib.sha256(payload).hexdigest() != entry['projection_sha256']
            or file_stamp(projection_path) != projection_before):
        raise ValueError('Summary projection bytes differ from manifest: ' + str(projection_path))
    value = json.loads(payload)
    if (value.get('program') != 'trm-summary-projection' or value.get('schema_version') != SCHEMA
            or value.get('summarizer_version') != version or value.get('sample_id') != source.stem):
        raise ValueError('Unsupported or mismatched summary projection; use --full-json for independent recalculation')
    if (value.get('source_sha256') != entry['source_sha256']
            or value.get('source_bytes') != entry['source_bytes'] or before[2] != entry['source_bytes']):
        raise ValueError('Summary projection source binding differs from manifest')
    digest = stream_digest(source)
    if digest != entry['source_sha256'] or before != file_stamp(source):
        raise ValueError('Assessment bytes changed since summary projection was generated: ' + str(source))
    row, expected, queue = value['row'], value['expected'], value['queue']
    if not isinstance(expected, list) or not isinstance(queue, list):
        raise ValueError('Summary projection tables must be arrays')
    for rows, fields in (([row], row_fields), (expected, expected_fields), (queue, queue_fields)):
        for item in rows:
            if (not isinstance(item, dict) or set(item) != set(fields)
                    or item.get('sample_id') != source.stem or item.get('source_sha256') != digest):
                raise ValueError('Summary projection row identity or schema differs')
            # Report directories can be copied; retain the verified bytes and
            # link to their current location just as the full-JSON path does.
            item['source_json'] = str(source)
    return row, expected, queue, {'path': str(source), 'sha256': digest, 'bytes': before[2]}
