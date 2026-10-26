#!/usr/bin/env python3
"""Summarize saved assessments on separate observation and confirmation axes.

Use verified producer projections when available, or read one sample JSON at
a time (--full-json forces independent recalculation). Never rescan reads, reclassify candidates, infer
taxonomy, or turn a clade prior into species validation. Every candidate stays
in its immutable source; every stored provisional/geometry nomination enters
the separate review queue, including nominations outside the representative.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re


VERSION = '1.1.0'
SAMPLE_FIELDS = ['sample_id', 'scientific_name', 'tax_id', 'classification',
    'representative_motif', 'analysis_complete', 'candidate_count', 'positive_current_families',
    'known_positive_families', 'known_positive_motifs', 'known_observed_families',
    'known_observed_motifs', 'known_pattern_families',
    'expected_family_status', 'expected_positive_families', 'expected_unreported_families',
    'expected_explicit_zero_families', 'expected_disagreement_families',
    'reference_scope', 'reference_status', 'reference_strength',
    'provisional_candidate_families', 'geometry_only_families', 'secondary_provisional_families',
    'secondary_geometry_only_families', 'main_trm_status', 'externally_confirmed_recorded',
    'adequacy_status', 'missing_power_assumptions', 'missing_power_provenance',
    'library_limited', 'library_reasons', 'uniform_wgs_model_blockers',
    'library_selection', 'protocol_available', 'ordinary_pcr_selection',
    'whole_genome_amplification_metadata', 'explicit_wga_blocker',
    'source_json', 'source_sha256']
SAMPLE_FIELDS += ['clade_variant_status', 'clade_variant_observed_families',
                 'clade_canonical_observed_families', 'published_clade_exception_families',
                 'clade_variant_review_required']
EXPECTED_FIELDS = ['sample_id', 'scientific_name', 'canonical_family', 'expected_motifs',
    'reference_scope', 'reference_status', 'reference_strength', 'family_status',
    'exact_measurement_status', 'current_status', 'current_reads', 'current_read_counts',
    'original_status', 'original_reads', 'explicit_zero', 'measurement_disagreement',
    'nearby_family_count', 'source_json', 'source_sha256']
REVIEW_FIELDS = ['priority', 'sample_id', 'scientific_name', 'reason', 'canonical_family',
    'motif', 'current_reads', 'is_representative', 'reference_strength', 'action',
    'source_json', 'source_sha256']
LIMITATIONS = [
    'Known-family observations, provisional repeats and physical terminal confirmation are separate axes.',
    'Catalogue-external nominations do not establish literature-wide novelty or a replacement telomeric motif.',
    'Clade/genus expectations and reused species references are not held-out accuracy measurements.',
    'Missing expected-family reports are not measured zeros; an explicit zero does not establish biological absence.',
    'Ordinary PCR selection does not establish whole-genome amplification; unknown preparation remains unknown.',
    'All candidates remain in the linked source JSON. No source candidate is edited, removed or promoted.',
    'Queue priority is a deterministic evidence-review workflow, not a calibrated probability or biological ranking.',
]


def reference_strength(scope, status):
    if status == 'supported':
        return {'species': 'species_reference', 'isolate': 'isolate_reference',
                'genus': 'genus_expectation', 'clade': 'weak_clade_expectation'}.get(scope, 'unspecified_scope')
    return {'candidate': 'published_candidate_nonbinding', 'out_of_scope': 'out_of_scope'}.get(status, 'unknown_reference')


def positive(value):
    return type(value) is int and value > 0


def joined(values):
    return ';'.join(sorted(set(str(value) for value in values)))


def summarize_sample(document, source, digest):
    sample = document.get('sample_id')
    if (document.get('program') != 'trm-assessment' or not isinstance(sample, str)
            or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample)
            or not isinstance(document.get('candidates'), list)):
        raise ValueError(f'{source}: not a valid saved assessment')
    candidates = document['candidates']
    families = [candidate.get('canonical_family') for candidate in candidates]
    if any(not isinstance(family, str) or not family for family in families) or len(set(families)) != len(families):
        raise ValueError(f'{sample}: missing or duplicate candidate family')
    assignment = document.get('reference_assignment') or {}
    agreement = document.get('reference_agreement') or {}
    adequacy = document.get('data_adequacy') or {}
    library = document.get('library_context') or {}
    clade = document.get('clade_variant_observation') or {}
    metadata = document.get('metadata') or {}
    scope, status = assignment.get('reference_scope', ''), assignment.get('reference_status', '')
    strength = reference_strength(scope, status)
    representative = (document.get('top_candidate') or {}).get('canonical_family')
    base = {'sample_id': sample, 'scientific_name': document.get('scientific_name'),
            'source_json': str(source), 'source_sha256': digest}
    queue = []

    def review(priority, reason, action, candidate=None, family=None):
        candidate = candidate or {}
        key = family or candidate.get('canonical_family', '')
        queue.append({**base, 'priority': priority, 'reason': reason, 'canonical_family': key,
            'motif': candidate.get('motif', ''), 'current_reads': candidate.get('reads'),
            'is_representative': bool(key and key == representative),
            'reference_strength': strength, 'action': action})

    expected, expected_counts = [], Counter()
    family_view = document.get('family_observation') or {}
    for family in family_view.get('families', []):
        current, original = family.get('exact_current') or {}, family.get('exact_original') or {}
        row_strength = reference_strength(family.get('reference_scope'), family.get('reference_status'))
        explicit_zero = current.get('explicit_zero_observed') is True
        disagreement = family.get('measurement_disagreement') is True
        unreported = current.get('status') == 'not_reported'
        expected_counts.update(positive=int(current.get('positive_observed') is True),
                               unreported=int(unreported), explicit_zero=int(explicit_zero),
                               disagreement=int(disagreement))
        expected.append({**base, 'canonical_family': family.get('canonical_family'),
            'expected_motifs': joined(family.get('expected_motifs', [])),
            'reference_scope': family.get('reference_scope'), 'reference_status': family.get('reference_status'),
            'reference_strength': row_strength, 'family_status': family.get('status'),
            'exact_measurement_status': family.get('exact_measurement_status'),
            'current_status': current.get('status'), 'current_reads': current.get('reads'),
            'current_read_counts': json.dumps(current.get('read_counts', []), separators=(',', ':')),
            'original_status': original.get('status'), 'original_reads': original.get('reads'),
            'explicit_zero': explicit_zero, 'measurement_disagreement': disagreement,
            'nearby_family_count': family.get('nearby_family_count')})
        key = family.get('canonical_family')
        if disagreement:
            review(20, 'expected_measurement_disagreement', 'Compare source/settings; retain both observations.', family=key)
        elif unreported:
            review(30, 'expected_family_not_reported', 'Review reference applicability and consider registered full-input measurement; missing is not zero.', family=key)
        elif explicit_zero:
            review(35, 'expected_family_explicit_zero', 'Review detector scope and conditional power before biological interpretation.', family=key)

    provisional, geometry, secondary_provisional, secondary_geometry = 0, 0, 0, 0
    for candidate in candidates:
        if candidate.get('novel_candidate') is True:
            provisional += 1
            secondary_provisional += candidate['canonical_family'] != representative
            review(50, 'provisional_repeat_requires_external_validation',
                   'Check chromosome-terminal location/polarity, independent biological library and literature novelty.', candidate)
        elif candidate.get('geometry_nomination') is True:
            geometry += 1
            secondary_geometry += candidate['canonical_family'] != representative
            review(60, 'geometry_nomination_requires_evidence',
                   'Review stored sequence/flank limitations, internal-array evidence and independent terminal location.', candidate)
    if document.get('analysis_complete') is not True or adequacy.get('status') == 'observed_qc_failure':
        review(10, 'analysis_or_qc_incomplete', 'Resolve input/analysis failure before biological attribution.')
    if document.get('classification') == 'taxonomic_scope_unresolved':
        review(25, 'taxonomic_attribution_unresolved', 'Resolve official metadata identity/lineage; retain observed repeats.')
    if library.get('limited') is True:
        review(40, 'library_context_limited', 'Review construction and genomic representation; ordinary PCR is not WGA.')
    if not expected:
        review(70, 'no_applicable_expected_family', 'Review species/strain primary literature without forcing a broader default.')
    if adequacy.get('missing_assumptions') or adequacy.get('missing_assumption_provenance'):
        review(80, 'conditional_power_inputs_missing', 'Document sourced genome/library assumptions; do not infer absence from read totals.')
    by_family = {candidate['canonical_family']: candidate for candidate in candidates}
    for item in clade.get('review_items', []):
        review(45, 'published_clade_exception_requires_sample_validation' if item.get('published_scoped_exception')
               else 'clade_variant_requires_scoped_validation',
               item['next_action'], by_family.get(item['canonical_family']), item['canonical_family'])
    row = {**base, 'tax_id': document.get('tax_id'), 'classification': document.get('classification'),
        'representative_motif': document.get('detected_motif'),
        'analysis_complete': document.get('analysis_complete') is True,
        'candidate_count': len(candidates),
        'positive_current_families': sum(positive(c.get('reads')) for c in candidates),
        'known_positive_families': sum(c.get('reference_status') == 'established' and positive(c.get('reads')) for c in candidates),
        'known_positive_motifs': joined(c.get('motif', c['canonical_family']) for c in candidates
                                       if c.get('reference_status') == 'established' and positive(c.get('reads'))),
        'known_observed_families': sum(c.get('known_observed') is True for c in candidates),
        'known_observed_motifs': joined(c.get('motif', c['canonical_family']) for c in candidates
                                       if c.get('known_observed') is True),
        'known_pattern_families': sum(c.get('candidate_status') == 'known_trm_pattern_supported' for c in candidates),
        'expected_family_status': family_view.get('status', 'unavailable'),
        'expected_positive_families': expected_counts['positive'],
        'expected_unreported_families': expected_counts['unreported'],
        'expected_explicit_zero_families': expected_counts['explicit_zero'],
        'expected_disagreement_families': expected_counts['disagreement'],
        'reference_scope': scope, 'reference_status': status, 'reference_strength': strength,
        'provisional_candidate_families': provisional, 'geometry_only_families': geometry,
        'secondary_provisional_families': secondary_provisional, 'secondary_geometry_only_families': secondary_geometry,
        'main_trm_status': agreement.get('main_trm_status', 'unavailable'),
        'externally_confirmed_recorded': document.get('externally_confirmed') is True,
        'adequacy_status': adequacy.get('status', 'unavailable'),
        'missing_power_assumptions': joined(adequacy.get('missing_assumptions', [])),
        'missing_power_provenance': joined(adequacy.get('missing_assumption_provenance', [])),
        'library_limited': library.get('limited') is True, 'library_reasons': joined(library.get('reasons', [])),
        'uniform_wgs_model_blockers': joined((adequacy.get('library') or {}).get('uniform_wgs_model_blockers', [])),
        'library_selection': metadata.get('library_selection', ''),
        'protocol_available': bool(str(metadata.get('library_construction_protocol') or '').strip()),
        'ordinary_pcr_selection': str(metadata.get('library_selection', '')).strip().upper() in {'PCR', 'RANDOM PCR'},
        'whole_genome_amplification_metadata': metadata.get('whole_genome_amplification', 'unknown'),
        'explicit_wga_blocker': 'explicit_whole_genome_amplification' in (adequacy.get('library') or {}).get('uniform_wgs_model_blockers', []),
        'clade_variant_status': clade.get('status', 'not_assessed'),
        'clade_variant_observed_families': clade.get('variant_observed_family_count'),
        'clade_canonical_observed_families': clade.get('canonical_observed_family_count'),
        'published_clade_exception_families': clade.get('published_exception_family_count'),
        'clade_variant_review_required': clade.get('review_required')}
    return row, expected, queue


def write_table(path, rows, fields):
    with path.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def assessment_sources(assessment, allow_json_only=False):
    """Choose the current cohort from its public manifest, never stale JSONs."""
    controls, pins = {}, []
    for name in ('interpretation_status.json', 'metrics.json', 'samples.tsv', 'summary_projections.tsv'):
        path = assessment / name
        if path.is_file():
            data = path.read_bytes()
            controls[name] = data
            pins.append({'path': str(path), 'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)})
        else:
            controls[name] = None
    status = json.loads(controls['interpretation_status.json']) if controls['interpretation_status.json'] is not None else None
    metrics = json.loads(controls['metrics.json']) if controls['metrics.json'] is not None else None
    if status is not None and status.get('status') != 'complete':
        raise ValueError('Assessment interpretation status is not complete')
    manifest = controls['samples.tsv']
    if manifest is None:
        if status is not None or metrics is not None or not allow_json_only:
            raise ValueError('Assessment requires authoritative samples.tsv; JSON-only fixtures need explicit --allow-json-only')
        paths = sorted((assessment / 'samples').glob('*.json'))
        authority = 'explicit_json_only_fallback_no_cohort_manifest'
    else:
        table = csv.DictReader(io.StringIO(manifest.decode('utf-8')), delimiter='\t')
        if 'sample_id' not in (table.fieldnames or []):
            raise ValueError('samples.tsv has no sample_id column')
        ids = [row['sample_id'] for row in table]
        if any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value or '') for value in ids) or len(set(ids)) != len(ids):
            raise ValueError('samples.tsv has invalid or duplicate sample IDs')
        if metrics is None:
            raise ValueError('Manifested assessment requires completed metrics.json')
        if type(metrics.get('samples')) is not int or metrics['samples'] != len(ids):
            raise ValueError('metrics.samples differs from samples.tsv')
        manifest_hash = (metrics.get('files') or {}).get('samples.tsv')
        if manifest_hash is not None and manifest_hash != hashlib.sha256(manifest).hexdigest():
            raise ValueError('samples.tsv checksum differs from metrics.json')
        for label, control in [('metrics', metrics), ('interpretation status', status)]:
            if control is None:
                continue
            if 'sample_ids' in control:
                values = control['sample_ids']
                if (not isinstance(values, list) or any(not isinstance(value, str) for value in values)
                        or len(set(values)) != len(values) or sorted(values) != sorted(ids)):
                    raise ValueError(f'{label} sample_ids differ from samples.tsv')
            if 'samples' in control and (type(control['samples']) is not int or control['samples'] != len(ids)):
                raise ValueError(f'{label} sample count differs from samples.tsv')
        if status is not None:
            if ('completed' in status and status['completed'] != len(ids)) or status.get('failed', 0) != 0:
                raise ValueError('Complete interpretation status has incomplete/failed samples')
            if metrics.get('execution_id') and metrics['execution_id'] != status.get('execution_id'):
                raise ValueError('Metrics and interpretation status execution IDs differ')
        paths = [assessment / 'samples' / f'{sample}.json' for sample in sorted(ids)]
        authority = 'samples.tsv'
    if not paths:
        raise ValueError('Assessment cohort has no samples')
    return paths, controls, pins, authority


def verify_controls_unchanged(assessment, controls):
    for name, before in controls.items():
        path = assessment / name
        after = path.read_bytes() if path.is_file() else None
        if before != after:
            raise ValueError(f'Assessment cohort control changed during summary: {name}')


def summarize(assessment, output, allow_json_only=False, *, use_projections=True):
    assessment, output = Path(assessment).resolve(), Path(output).resolve()
    paths, controls, control_pins, authority = assessment_sources(assessment, allow_json_only)
    from summary_projection import manifest_entries, read_projection
    projections = manifest_entries(controls, paths) if use_projections else {}
    if output == assessment or assessment in output.parents:
        raise ValueError('Summary output must be separate from the input assessment')
    output.mkdir(parents=True, exist_ok=False)
    rows, expected, queue, pins = [], [], [], []
    seen = set()
    projected_count = 0
    for path in paths:
        if path.stem in projections:
            row, sample_expected, sample_queue, pin = read_projection(path, projections[path.stem],
                version=VERSION, row_fields=SAMPLE_FIELDS, expected_fields=EXPECTED_FIELDS, queue_fields=REVIEW_FIELDS)
            if row['sample_id'] in seen:
                raise ValueError('Duplicate sample identity in summary projection')
            seen.add(row['sample_id'])
            rows.append(row); expected.extend(sample_expected); queue.extend(sample_queue); pins.append(pin)
            projected_count += 1
            continue
        before = path.stat()
        payload = path.read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        document = json.loads(payload)
        del payload
        after = path.stat()
        if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
                after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns):
            raise ValueError(f'Source changed while reading: {path}')
        row, sample_expected, sample_queue = summarize_sample(document, path, digest)
        if row['sample_id'] != path.stem or row['sample_id'] in seen:
            raise ValueError('Sample filename/identity differs or duplicates another sample')
        seen.add(row['sample_id'])
        rows.append(row); expected.extend(sample_expected); queue.extend(sample_queue)
        pins.append({'path': str(path), 'sha256': digest, 'bytes': after.st_size})
        del document
    verify_controls_unchanged(assessment, controls)
    queue.sort(key=lambda row: (row['priority'], row['sample_id'], row['reason'], row['canonical_family'] or ''))
    expected.sort(key=lambda row: (row['sample_id'], row['canonical_family'] or ''))
    write_table(output / 'samples.tsv', rows, SAMPLE_FIELDS)
    write_table(output / 'expected_families.tsv', expected, EXPECTED_FIELDS)
    write_table(output / 'review_queue.tsv', queue, REVIEW_FIELDS)
    counts = lambda field: dict(sorted(Counter(row[field] for row in rows).items()))
    summary = {'schema_version': 1, 'program': 'trm-interpretation-summary', 'version': VERSION,
        'created_utc': datetime.now(timezone.utc).isoformat(), 'status': 'complete',
        'assessment': str(assessment), 'samples': len(rows),
        'cohort_authority': authority, 'cohort_control_pins': control_pins,
        'producer_projection_samples': projected_count,
        'independently_recalculated_samples': len(paths) - projected_count,
        'source_verification': 'full_file_sha256_including_projection_mode',
        'classification_counts': counts('classification'), 'main_trm_status_counts': counts('main_trm_status'),
        'family_observation_counts': counts('expected_family_status'),
        'reference_strength_counts': counts('reference_strength'), 'adequacy_status_counts': counts('adequacy_status'),
        'candidate_rows_examined': sum(row['candidate_count'] for row in rows),
        'candidate_count_scope': 'All source candidates represented; projection-mode rows were examined by the report producer, not redecoded here.',
        'clade_variant_status_counts': counts('clade_variant_status'),
        'clade_variant_review_samples': sum(row['clade_variant_review_required'] is True for row in rows),
        'samples_with_known_positive_family': sum(row['known_positive_families'] > 0 for row in rows),
        'samples_with_provisional_candidates': sum(row['provisional_candidate_families'] > 0 for row in rows),
        'provisional_candidate_families_including_secondary': sum(row['provisional_candidate_families'] for row in rows),
        'geometry_only_families_including_secondary': sum(row['geometry_only_families'] for row in rows),
        'samples_with_known_and_provisional_families': sum(row['known_positive_families'] > 0 and row['provisional_candidate_families'] > 0 for row in rows),
        'expected_positive_by_reference_strength': dict(sorted(Counter(row['reference_strength'] for row in expected if row['current_status'] == 'reported_positive').items())),
        'library_limited_samples': sum(row['library_limited'] for row in rows),
        'ordinary_pcr_selection_samples': sum(row['ordinary_pcr_selection'] for row in rows),
        'explicit_wga_blocker_samples': sum(row['explicit_wga_blocker'] for row in rows),
        'missing_library_protocol_samples': sum(not row['protocol_available'] for row in rows),
        'missing_power_assumption_samples': sum(bool(row['missing_power_assumptions']) for row in rows),
        'missing_power_assumption_counts': dict(sorted(Counter(key for row in rows
            for key in row['missing_power_assumptions'].split(';') if key).items())),
        'missing_power_provenance_counts': dict(sorted(Counter(key for row in rows
            for key in row['missing_power_provenance'].split(';') if key).items())),
        'uniform_wgs_model_blocker_counts': dict(sorted(Counter(key for row in rows
            for key in row['uniform_wgs_model_blockers'].split(';') if key).items())),
        'review_queue_rows': len(queue), 'review_reason_counts': dict(sorted(Counter(row['reason'] for row in queue).items())),
        'source_pins': pins, 'limitations': LIMITATIONS,
        'classification_changed': False, 'candidate_filtering_performed': False,
        'production_raw_inputs_opened': False, 'downloads_started': 0}
    with (output / 'summary.json').open('x') as handle:
        json.dump(summary, handle, indent=2, allow_nan=False); handle.write('\n')
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assessment', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--allow-json-only', action='store_true',
                        help='Explicit fixture/standalone fallback only when no cohort controls exist; cannot bypass report status/metrics')
    parser.add_argument('--full-json', action='store_true',
                        help='Independently recalculate from complete candidate JSON instead of using verified producer projections')
    args = parser.parse_args(argv)
    try:
        result = summarize(args.assessment, args.output, args.allow_json_only, use_projections=not args.full_json)
    except (OSError, ValueError, TypeError, KeyError) as error:
        parser.error(str(error))
    print(json.dumps({key: result[key] for key in ('status', 'samples', 'candidate_rows_examined', 'review_queue_rows')}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
