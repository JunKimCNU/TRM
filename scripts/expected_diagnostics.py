"""Join separately measured expected families without replacing discovery evidence.

An optional diagnostic directory contains SAMPLE/evidence.json and status.json.
Only the expectation table is supplemented. Candidate ranking, geometry, paired
annotations, observed counts and adequacy retain their existing evidence source.
"""
from __future__ import annotations

import json
from pathlib import Path

from evidence_batch import canonical, sha256


DETECTION_FIELDS = ('mode', 'seed_bases', 'min_copies', 'min_identity', 'min_tract', 'end_window')
COUNT_FIELDS = ('reads', 'full_repeat_reads', 'repeat_bases', 'max_tract')


def supplement_expected_observations(sample, observations, directory, raw, source_hash,
                                     analysis_complete=True):
    """Return a new expectation table and provenance, or fail on a stale/partial join."""
    if directory is None:
        return observations, None
    folder = Path(directory) / sample
    status_path, evidence_path = folder / 'status.json', folder / 'evidence.json'
    if not folder.exists():
        return observations, None
    if not status_path.is_file() or not evidence_path.is_file():
        raise ValueError(f'{sample}: expected-family diagnostic is incomplete')
    if not analysis_complete:
        raise ValueError(f'{sample}: expected-family diagnostic cannot repair an incomplete source analysis')
    status, evidence = (json.loads(p.read_text()) for p in (status_path, evidence_path))
    identity = status.get('identity', {})
    if status.get('sample_id') != sample or status.get('status') != 'complete' or status.get('exit_code', 0) != 0:
        raise ValueError(f'{sample}: expected-family diagnostic is not a completed matching sample')
    if identity.get('source_sha256') != source_hash:
        raise ValueError(f'{sample}: expected-family diagnostic source hash differs')
    evidence_hash = sha256(evidence_path)
    if status.get('result_sha256') != evidence_hash:
        raise ValueError(f'{sample}: expected-family diagnostic result hash differs')
    if evidence.get('analysis_kind') != 'targeted_evidence':
        raise ValueError(f'{sample}: expected-family diagnostic must be targeted evidence')
    if not raw.get('version') or evidence.get('version') != raw['version']:
        raise ValueError(f'{sample}: expected-family diagnostic detector version differs')
    if not raw.get('inputs') or evidence.get('inputs') != raw['inputs']:
        raise ValueError(f'{sample}: expected-family diagnostic input identity differs')
    if evidence.get('summary', {}).get('sampled_prefix') is not False:
        raise ValueError(f'{sample}: expected-family diagnostic must cover the full input')
    for field in ('records', 'bases'):
        if evidence.get('summary', {}).get(field) != raw.get('summary', {}).get(field):
            raise ValueError(f'{sample}: expected-family diagnostic {field} differs')
    for field in DETECTION_FIELDS:
        if field not in raw.get('config', {}) or evidence.get('config', {}).get(field) != raw['config'][field]:
            raise ValueError(f'{sample}: expected-family diagnostic detection setting {field} differs')
    if evidence.get('junction_options') != raw.get('junction_options'):
        raise ValueError(f'{sample}: expected-family diagnostic junction settings differ')
    targets = identity.get('motifs', identity.get('targets', {}).get('motifs'))
    if not isinstance(targets, list) or not targets:
        raise ValueError(f'{sample}: expected-family diagnostic lacks registered targets')
    target_families = {canonical(m) for m in targets}
    measured = {}
    for row in evidence.get('candidate_evidence', []):
        family = canonical(row['motif'])
        if family in measured:
            raise ValueError(f'{sample}: duplicate expected-family diagnostic family')
        for key in COUNT_FIELDS:
            value = row.get(key)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f'{sample}: expected-family diagnostic has invalid {key}')
        if not 0 <= row['full_repeat_reads'] <= row['reads'] <= evidence['summary']['records']:
            raise ValueError(f'{sample}: expected-family diagnostic count ordering is invalid')
        if row['reads'] == 0 and any(row[key] for key in COUNT_FIELDS):
            raise ValueError(f'{sample}: zero diagnostic reads have nonzero supporting counts')
        measured[family] = row
    if set(measured) != target_families:
        raise ValueError(f'{sample}: expected-family diagnostic target set differs')
    native = {canonical(c['motif']): c for c in raw.get('candidates', []) + raw.get('rejected_candidates', [])}
    updated, used = [], []
    for observation in observations:
        family = canonical(observation['motif'])
        if family not in measured:
            updated.append(observation)
            continue
        row = measured[family]
        previous = native.get(family)
        disagreement = previous is not None and any(previous.get(k) != row[k] for k in COUNT_FIELDS)
        recovery = ('measurement_disagreement' if disagreement else
                    'reported_positive' if row['reads'] else 'explicit_zero')
        updated.append({**observation, 'recovery': recovery, 'reads': row['reads'],
                        'full_repeat_reads': row['full_repeat_reads'],
                        'repeat_bases': row['repeat_bases'], 'max_tract_bp': row['max_tract'],
                        'source': 'expected_family_diagnostic',
                        'original_discovery_reads': previous.get('reads') if previous else None,
                        'original_discovery_reported': previous is not None,
                        'measurement_disagreement': disagreement,
                        'interpretation': 'Full-input targeted measurement under the recorded detector settings. '
                        'Zero is not biological absence. Discovery candidates, ranking and adequacy are unchanged.'})
        used.append(family)
    provenance = {'status': 'complete', 'result_file': str(evidence_path.resolve()),
                  'result_sha256': evidence_hash, 'status_file': str(status_path.resolve()),
                  'status_sha256': sha256(status_path), 'source_sha256': source_hash,
                  'detector_version': evidence['version'], 'registered_families': sorted(target_families),
                  'used_expected_families': sorted(used),
                  'unused_registered_families': sorted(target_families - set(used)),
                  'wall_seconds': status.get('wall_seconds'), 'max_rss_kb': status.get('max_rss_kb'),
                  'scope': 'Expected-family table only; no change to discovery candidates, ranking, geometry or adequacy.'}
    return updated, provenance
