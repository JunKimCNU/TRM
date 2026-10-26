"""Validate and display optional native Phred+33 observations without filtering."""
from __future__ import annotations

from copy import deepcopy


def _count(value, label):
    if type(value) is not int or value < 0:
        raise ValueError('Invalid quality count: ' + label)
    return value


def _histogram(value):
    if not isinstance(value, dict):
        raise ValueError('Quality histogram must be an object')
    counts = value.get('phred33_counts')
    if not isinstance(counts, list) or len(counts) != 94:
        raise ValueError('Quality histogram requires 94 Phred+33 bins')
    for item in counts:
        _count(item, 'histogram bin')
    for key in ('bases', 'below_q20_bases', 'below_q30_bases'):
        _count(value.get(key), key)
    if (sum(counts) != value['bases'] or sum(counts[:20]) != value['below_q20_bases']
            or sum(counts[:30]) != value['below_q30_bases']):
        raise ValueError('Quality histogram totals disagree')
    return value['bases']


def _records(value, expected=None):
    if not isinstance(value, dict):
        raise ValueError('Quality record counts must be an object')
    for key in ('observed_records', 'quality_records', 'fasta_records_without_quality',
                'fastq_records_above_retention_limit'):
        _count(value.get(key), key)
    if value['observed_records'] != sum(value[key] for key in (
            'quality_records', 'fasta_records_without_quality', 'fastq_records_above_retention_limit')):
        raise ValueError('Quality record availability counts disagree')
    if expected is not None and value['observed_records'] != expected:
        raise ValueError('Quality record count differs from sequence observations')
    bases = _histogram(value.get('whole_read'))
    if not value['quality_records'] and bases:
        raise ValueError('Quality bases are present without quality-bearing records')
    # Native input records are nonempty, and retained quality covers the whole read.
    if bases < value['quality_records']:
        raise ValueError('Quality bases are fewer than quality-bearing records')
    return bases


def describe_quality(report, source_sha256):
    value = report.get('quality_diagnostics')
    result = {'status': 'not_measured', 'encoding': None, 'source_sha256': source_sha256,
              'input_scope': None, 'diagnostics': None,
              'filtering_applied': False, 'changes_classification_or_adequacy': False,
              'variant_discriminating_bases_assessed': False,
              'quality_authenticity_or_library_retention_verified': False,
              'scope': 'Supplied FASTQ qualities under explicitly declared encoding; no inference of variant-base accuracy, independent molecules, or telomere retention.'}
    if value is None:
        return result
    if (not isinstance(value, dict) or value.get('schema_version') != 1 or value.get('encoding') != 'phred33'
            or value.get('encoding_provenance') != 'explicit_user_declaration_not_inferred'
            or value.get('input_scope') not in {'complete_inputs', 'analyzed_prefix'}
            or value.get('filtering_applied') is not False
            or value.get('variant_discriminating_bases_assessed') is not False):
        raise ValueError('Unsupported native quality diagnostic contract')
    summary = report.get('summary')
    if not isinstance(summary, dict):
        raise ValueError('Quality diagnostics require sequence summary counts')
    sequence_records = _count(summary.get('records'), 'summary.records')
    sequence_bases = _count(summary.get('bases'), 'summary.bases')
    records = value.get('input_records')
    bases = _records(records, sequence_records)
    if bases > sequence_bases:
        raise ValueError('Quality bases exceed analyzed sequence bases')
    missing = records['fasta_records_without_quality'] + records['fastq_records_above_retention_limit']
    if not missing and bases != sequence_bases:
        raise ValueError('Complete quality bases differ from analyzed sequence bases')
    status = ('unavailable_in_analyzed_records' if not records['quality_records'] else
              'partially_measured' if missing else 'measured')
    result.update(status=status, encoding='phred33', input_scope=value['input_scope'], diagnostics=deepcopy(value))
    return result


def describe_candidate_quality(candidate, context):
    """Do not fill a targeted measurement's missing quality from old discovery."""
    value = candidate.get('quality_diagnostics')
    if value is None or context['diagnostics'] is None:
        return None
    if not isinstance(value, dict):
        raise ValueError('Candidate quality diagnostics must be an object')
    _records(value.get('supporting_records'), candidate.get('reads'))
    measured = _histogram(value.get('union_tract'))
    union = _count(value.get('union_tract_bases'), 'union_tract_bases')
    missing = _count(value.get('union_tract_bases_without_quality'), 'union_tract_bases_without_quality')
    if measured + missing != union or union != candidate.get('repeat_bases'):
        raise ValueError('Candidate quality tract counts differ from repeat observations')
    return {'status': 'measured' if measured and not missing else 'partially_measured' if measured else 'unavailable',
            'encoding': context['encoding'], 'source_sha256': context['source_sha256'],
            'input_scope': context['input_scope'], 'filtering_applied': False,
            'variant_discriminating_bases_assessed': False, 'diagnostics': deepcopy(value)}
