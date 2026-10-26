#!/usr/bin/env python3
"""Independently check v0.4/v0.5/v0.6 report artifacts without opening sequence inputs.

Checks consistency and preserved evidence, not the biological truth of a TRM.
Default expectations describe the current 230-sample batch; explicit count/path
options allow use on smaller fixtures. No report implementation is imported.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import json
import math
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit


REQUIRED_ASSUMPTIONS = ('genome_size_bp', 'haploid_chromosome_count', 'random_wgs_assumed',
                        'telomere_retention_fraction', 'effective_independent_read_fraction',
                        'detector_sensitivity')
KNOWN_CLASSIFICATIONS = {'known_trm_detected', 'known_repeat_detected_boundary_unresolved'}
OBSERVATION_CLASSIFICATIONS = {'known_repeat_observed_pattern_unresolved',
                              'published_candidate_observed_pattern_unresolved',
                              'published_repeat_trace_observed'}
ESTABLISHED_EVIDENCE = {'experimental', 'experimental_terminal', 'assembly_terminal',
                        'published_terminal_long_read', 'long_read_terminal', 'clade_consensus'}
ALLOWED_CLASSIFICATIONS = KNOWN_CLASSIFICATIONS | {
    'analysis_failed', 'detection_inconclusive_technical', 'detection_inconclusive_measurement',
    'taxonomic_scope_unresolved',
    'novel_trm_candidate', 'novel_geometry_nomination', 'novel_provisional_candidate',
    'known_alternative_architecture', 'library_context_limited',
    'published_candidate_supported', 'detection_inconclusive_depth',
    'canonical_not_detected_conditional_power', 'detection_inconclusive'} | OBSERVATION_CLASSIFICATIONS


def read_json(path):
    return json.loads(Path(path).read_text())


def read_tsv(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def canonical(motif):
    """Independent primitive/rotation/strand normalization of a DNA repeat."""
    assert motif and not set(motif) - set('ACGT'), f'invalid motif {motif!r}'
    unit = next(motif[:k] for k in range(1, len(motif) + 1)
                if len(motif) % k == 0 and motif[:k] * (len(motif) // k) == motif)
    reverse = unit.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    return min(sequence[k:] + sequence[:k] for sequence in (unit, reverse)
               for k in range(len(unit)))


def index_rows(rows, key='sample_id'):
    result = {}
    for row in rows:
        identity = row[key]
        assert identity and identity not in result, f'duplicate/empty {key}: {identity!r}'
        if key == 'sample_id':
            assert re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}', identity), 'unsafe sample ID'
        result[identity] = row
    return result


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = []

    def handle_starttag(self, tag, attrs):
        self.hrefs.extend(value for key, value in attrs if tag == 'a' and key == 'href')


def check_links(path, required=()):
    links = Links()
    links.feed(path.read_text())
    for href in links.hrefs:
        parsed = urlsplit(href)
        assert parsed.scheme not in {'javascript', 'data'}, f'unsafe HTML link {href!r}'
        if parsed.scheme or parsed.netloc or not parsed.path:
            continue
        assert (path.parent / unquote(parsed.path)).is_file(), f'broken HTML link {href!r}'
    assert set(required) <= set(links.hrefs), 'required sample/navigation HTML links missing'
    return links.hrefs


def equal(left, right, label):
    assert left == right, f'{label}: {left!r} != {right!r}'


def digest_matches(path, expected, label):
    assert expected and len(expected) == 64, f'{label}: missing SHA256'
    equal(sha256(path), expected, label)


class Audit:
    def __init__(self):
        self.checks = []

    def run(self, name, action, sample=None):
        row = {'check': name, 'status': 'passed'}
        if sample:
            row['sample_id'] = sample
        try:
            detail = action()
            if detail is not None:
                row['detail'] = detail
        except (AssertionError, OSError, ValueError, KeyError, TypeError, IndexError) as error:
            row.update(status='failed', detail=f'{type(error).__name__}: {error}')
        self.checks.append(row)

    def finish(self, args):
        counts = Counter(row['status'] for row in self.checks)
        return {'schema_version': 1, 'created_utc': datetime.now(timezone.utc).isoformat(),
                'scope': 'Artifact consistency/provenance audit only; no FASTQ opened and no biological absence inference.',
                'assessment': str(args.assessment.resolve()), 'comparison': str(args.comparison.resolve()),
                'passed': not counts['failed'], 'summary': dict(counts), 'checks': self.checks}


def check_sample_tables(document, summary, comparison, candidate_rows, adequacy_row, raw, path):
    sample = document['sample_id']
    for field in ('sample_id', 'scientific_name', 'result', 'reference_context'):
        equal(summary[field], str(document[field]), f'summary {field}')
    equal(summary['detected_motif'], document['detected_motif'] or '', 'summary motif')
    equal(summary['expected_motifs'], ','.join(document['expected_motifs']), 'expected motifs')
    equal(adequacy_row['status'], document['data_adequacy']['status'], 'adequacy table')
    equal(summary['details'], f'samples/{sample}.html', 'sample HTML path')
    for field in ('records', 'bases'):
        equal(summary[field], str(raw.get('summary', {}).get(field, '')), f'summary {field}')
        equal(document['input_summary'].get(field), raw.get('summary', {}).get(field), f'document {field}')
    for field in ('reference_status', 'repeat_architecture'):
        equal(summary[field], document['reference_assignment'].get(field, ''), f'reference {field}')
    top = document['top_candidate']
    assert not top or top in document['candidates'], 'top candidate absent from candidates'
    equal(document['detected_motif'], top['motif'] if top else None, 'top/detected motif')
    previous = raw.get('candidates', [])
    previous = previous[0] if previous else {}
    equal(comparison['previous_canonical_key'], previous.get('motif', ''), 'previous family')
    equal(comparison['previous_reference_anchor'], previous.get('telomere_prior', {}).get('matched_reference', ''), 'previous anchor')
    if comparison['previous_motif_display']:
        equal(canonical(comparison['previous_motif_display']), canonical(previous['motif']), 'previous display family')
    equal(comparison['updated_result'], document['result'], 'comparison result')
    equal(comparison['updated_motif_display'], document['detected_motif'] or '', 'comparison motif')
    equal(comparison['updated_canonical_key'], top['canonical_family'] if top else '', 'comparison family')
    equal(Path(comparison['updated_result_file']).resolve(), path.resolve(), 'comparison file')
    equal(Path(comparison['previous_result']).resolve(), Path(document['provenance']['original_discovery_result']).resolve(), 'comparison previous file')
    if document['analysis_complete']:
        changed = previous.get('motif', '') != (top['canonical_family'] if top else '')
        equal(comparison['motif_family_changed'], str(changed), 'family change flag')
    by_family = index_rows(candidate_rows, 'canonical_family')
    equal(set(by_family), {c['canonical_family'] for c in document['candidates']}, 'candidate family set')
    for candidate in document['candidates']:
        family = candidate['canonical_family']
        equal(canonical(candidate['motif']), family, 'display/canonical family')
        row = by_family[family]
        for field in ('motif', 'reference_status', 'candidate_status', 'reads', 'full_repeat_reads', 'repeat_bases'):
            equal(row[field], str(candidate[field]), f'candidate {family} {field}')
        for field in ('left_consistent', 'right_consistent', 'opposite_junctions'):
            table_field = field + '_junctions' if field != 'opposite_junctions' else field
            value = candidate['directionality'][field]
            equal(row[table_field], '' if value is None else str(value), f'direction {family} {field}')
        assert not (candidate['reference_status'] == 'established' and candidate['novel_candidate']), 'established family called novel'


def observation_schema(document):
    version = str(document.get('version', '')).split('.')
    version_pair = tuple(int(piece) for piece in version[:2]) if all(piece.isdigit() for piece in version[:2]) else ()
    return document.get('schema_version', 0) >= 6 or version_pair >= (0, 6)


def _nonnegative_integer(value, label):
    assert isinstance(value, int) and not isinstance(value, bool) and value >= 0, f'invalid {label}'
    return value


def _saved_families(document):
    rows = document.get('candidate_evidence', document.get('candidates', []) + document.get('rejected_candidates', []))
    return index_rows([dict(row, canonical_family=canonical(row['motif'])) for row in rows], 'canonical_family')


def check_observation_contract(document, raw=None, evidence=None):
    """v0.6 evidence contracts; may be called directly for native-only reports.

    With raw/evidence JSON arguments, verify counts against saved scanner
    records. No production decision or canonicalization helper is imported.
    Older schemas retain their existing audit and do not require new fields.
    """
    if not observation_schema(document):
        return 'legacy schema: new observation fields not required'
    policy = document['policy']
    equal(policy['minimum_known_reads'], 10, 'unchanged representative count threshold')
    equal(policy['minimum_full_repeat_reads'], 5, 'unchanged full-repeat pattern threshold')
    candidates = index_rows(document['candidates'], 'canonical_family')
    originals = _saved_families(raw) if raw is not None else None
    measured = _saved_families(evidence) if evidence is not None else None
    expected_families = {canonical(motif) for motif in document['expected_motifs']}
    for family, candidate in candidates.items():
        equal(canonical(candidate['motif']), family, 'observed motif identity')
        reads = _nonnegative_integer(candidate['reads'], family + ' reads')
        full = _nonnegative_integer(candidate['full_repeat_reads'], family + ' full-repeat reads')
        original = _nonnegative_integer(candidate['original_discovery_reads'], family + ' original reads')
        assert full <= reads, 'full-repeat records exceed current reads'
        equal(candidate.get('expected_family_match'), family in expected_families if expected_families else None,
              'expected family flag')
        if originals is not None and measured is not None:
            initial = originals.get(family, {})
            target = measured.get(family, {})
            complete = (all(key in target for key in ('reads', 'repeat_bases', 'max_tract'))
                        and isinstance(target.get('junction_evidence'), dict)
                        and 'full_repeat_records' in target['junction_evidence'])
            selected = target if complete else initial
            for output, source in (('reads', 'reads'), ('repeat_bases', 'repeat_bases'), ('max_tract_bp', 'max_tract')):
                equal(candidate[output], selected.get(source, 0), 'saved selected ' + output)
            equal(full, target['junction_evidence']['full_repeat_records'] if complete else initial.get('full_repeat_reads', 0),
                  'saved selected full-repeat count')
            equal(original, initial.get('reads', 0), 'saved original reads')
            count_source = ('targeted_only' if complete and original == 0 and reads > 0 else
                            'measured_evidence' if complete else 'original_discovery')
            equal(candidate['evidence_count_source'], count_source, 'complete count-source selection')
            equal(candidate['measurement_disagreement'], complete and original > 0 and reads == 0,
                  'saved measurement disagreement')
        disagreement = candidate['measurement_disagreement']
        eligible = reads >= 10 and not disagreement
        status = ('measurement_disagreement' if disagreement else 'sequence_observed' if eligible else
                  'trace_sequence_observed' if reads else 'no_sequence_observed')
        observation = candidate['sequence_observation']
        for field, value in (('status', status), ('observed', bool(reads or original)), ('current_reads', reads),
                             ('original_discovery_reads', original), ('count_source', candidate['evidence_count_source']),
                             ('minimum_representative_reads', 10), ('representative_eligible', eligible),
                             ('measurement_disagreement', disagreement)):
            equal(observation[field], value, 'sequence observation ' + field)
        established = candidate['reference_status'] == 'established'
        referenced = candidate['reference_status'] in {'established', 'published_candidate'}
        equal(candidate['known_observed'], established and eligible, 'known observation independent of pattern')
        equal(candidate['reference_observed'], referenced and eligible, 'published observation independent of pattern')
        direction = candidate['directionality']
        known_pattern_counts = established and reads >= 10 and full >= 5 and direction['status'] != 'orientation_conflict'
        equal(candidate['known_detected'], known_pattern_counts, 'known pattern unchanged gates')
        if direction['status'] == 'consistent_both_ends':
            left, right, opposite = (direction[field] for field in ('left_consistent', 'right_consistent', 'opposite_junctions'))
            total = left + right + opposite
            unique = direction['unique_junction_reads']
            assert unique is not None and unique >= 10 and min(left, right) >= 3, 'insufficient both-end junction counts'
            fraction = (left + right) / total
            z2 = 1.96 ** 2
            lower = (fraction + z2 / (2 * total) - 1.96 * math.sqrt(fraction * (1 - fraction) / total + z2 / (4 * total ** 2))) / (1 + z2 / total)
            assert fraction >= 0.8 and lower >= 0.6, 'insufficient directional pattern support'
        full_and_boundary = full >= 5 and direction['status'] == 'consistent_both_ends'
        if established:
            expected_status = ('known_trm_pattern_supported' if known_pattern_counts and full_and_boundary else
                               'known_repeat_detected_boundary_unresolved' if known_pattern_counts else 'weak_or_conflicting_known_repeat')
            equal(candidate['candidate_status'], expected_status, 'known candidate status from unchanged pattern gates')
        elif referenced:
            equal(candidate['candidate_status'], 'published_candidate_pattern_supported' if full_and_boundary else 'published_candidate_unresolved',
                  'published candidate evidence remains unconfirmed')
            equal(candidate['known_detected'], False, 'candidate cannot become known')
        assert not (referenced and candidate['novel_candidate']), 'published family called uncatalogued novel'

    top = document['top_candidate']
    assert not top or top == candidates.get(top['canonical_family']), 'representative differs from candidate evidence'
    equal(document['detected_motif'], top['motif'] if top else None, 'representative displayed identity')
    if top:
        equal(document['sequence_observation'], top['sequence_observation'], 'representative sequence observation')
    else:
        equal(document['sequence_observation']['observed'], False, 'no representative observation')
        equal(document['sequence_observation']['current_reads'], None, 'no representative is not explicit zero')
    code = document['classification']
    ranked = sorted(candidates.values(), key=lambda candidate: (-candidate['ranking_score'], candidate['canonical_family']))
    qualifying_known = [candidate for candidate in ranked if candidate['known_detected']]
    if code in KNOWN_CLASSIFICATIONS:
        assert top and top['known_detected'], 'known sample call without qualifying known pattern'
        equal(top, qualifying_known[0], 'representative strongest qualifying known family')
        equal(top['candidate_status'] == 'known_trm_pattern_supported', code == 'known_trm_detected', 'known sample boundary state')
    if code in OBSERVATION_CLASSIFICATIONS:
        assert top and not top['known_detected'], 'observation-only state promoted to known pattern'
        assert not qualifying_known, 'observation-only state ignores qualifying known pattern'
        assert document['analysis_complete'], 'observation-only biological call on incomplete analysis'
        assert document['data_adequacy']['status'] != 'observed_qc_failure', 'observation call bypassed QC block'
        assert not document.get('library_context', {}).get('limited'), 'observation call bypassed library block'
        if code == 'published_repeat_trace_observed':
            assert 0 < top['reads'] < 10 and top['expected_family_match'], 'trace state lacks a small expected-family observation'
            assert document['reference_assignment'].get('reference_scope') in {'species', 'isolate'}, 'trace state lacks exact taxon applicability'
        else:
            equal(top['sequence_observation']['representative_eligible'], True, 'observation count gate')
            equal(top['reference_status'], 'established' if code == 'known_repeat_observed_pattern_unresolved' else 'published_candidate',
                  'observation reference certainty')
            assert top['candidate_status'] != 'published_candidate_pattern_supported', 'supported candidate mislabeled unresolved observation'
    exact_scope = document['reference_assignment'].get('reference_scope') in {'species', 'isolate'}
    expected_observed = [candidate for candidate in ranked if exact_scope and candidate['expected_family_match']
                         and candidate['reference_status'] in {'established', 'published_candidate'}
                         and candidate['sequence_observation']['representative_eligible']]
    limited_states = {'analysis_failed', 'detection_inconclusive_technical', 'taxonomic_scope_unresolved',
                      'library_context_limited', 'known_alternative_architecture'}
    if expected_observed and not qualifying_known and code not in limited_states:
        equal(top, expected_observed[0], 'applicable observed reference precedes unrelated novelty')
        assert code in OBSERVATION_CLASSIFICATIONS | {'published_candidate_supported'}, 'observed expected family erased by sample conclusion'

    observations = document['expected_family_observations']
    equal([row['motif'] for row in observations], document['expected_motifs'], 'expected observation motif list')
    for row in observations:
        family = canonical(row['motif'])
        equal(row['canonical_family'], family, 'expected observation family')
        candidate = candidates.get(family)
        if candidate is None:
            equal(row['recovery'], 'not_reported', 'unsaved expected family recovery')
            equal(row['reads'], None, 'not-reported expected count')
            equal(row['full_repeat_reads'], None, 'not-reported full-repeat count')
            equal(row['candidate_status'], None, 'unsaved family candidate status')
            equal(row['reference_status'], None, 'unsaved family reference status')
            if originals is not None and measured is not None:
                assert family not in originals and family not in measured, 'saved expected family lost from candidates'
            continue
        equal(row['candidate_status'], candidate['candidate_status'], 'expected observation candidate state')
        equal(row['reference_status'], candidate['reference_status'], 'expected observation reference certainty')
        recovery = row['recovery']
        assert recovery in {'reported_positive', 'measurement_disagreement', 'explicit_zero', 'not_reported'}, 'invalid expected recovery state'
        if recovery == 'not_reported':
            equal(row['reads'], None, 'not-reported expected count')
            equal(row['full_repeat_reads'], None, 'not-reported full-repeat count')
            assert candidate['reads'] == 0 and not candidate['measurement_disagreement'], 'positive evidence hidden as not reported'
        else:
            equal(row['reads'], candidate['reads'], 'expected observation read count')
            equal(row['full_repeat_reads'], candidate['full_repeat_reads'], 'expected observation full-repeat count')
            if recovery == 'reported_positive':
                assert row['reads'] > 0, 'reported-positive expected family has no reads'
            elif recovery == 'measurement_disagreement':
                assert candidate['measurement_disagreement'], 'disagreement recovery without conflicting observations'
            else:
                equal(row['reads'], 0, 'explicit-zero expected count')
                equal(row['full_repeat_reads'], 0, 'explicit-zero full-repeat count')
                assert not candidate['measurement_disagreement'], 'explicit zero hides a positive original observation'
                assert candidate['evidence_count_source'] in {'measured_evidence', 'targeted_only'}, 'explicit zero without measured complete counts'
        if originals is not None and measured is not None:
            saved = family in originals or family in measured
            if not saved:
                equal(recovery, 'not_reported', 'unsaved family is not an explicit zero')
            elif candidate['reads'] > 0:
                equal(recovery, 'reported_positive', 'saved positive family recovery')
            elif candidate['measurement_disagreement']:
                equal(recovery, 'measurement_disagreement', 'conflicting family recovery')
            elif candidate['evidence_count_source'] in {'measured_evidence', 'targeted_only'}:
                equal(recovery, 'explicit_zero', 'saved measured zero family recovery')
    return {'candidate_observations': len(candidates), 'expected_families': len(observations)}


def check_interpretation(document, adequacy_row, raw=None, evidence=None):
    adequacy = document['data_adequacy']
    assert document['classification'] in ALLOWED_CLASSIFICATIONS, 'unsupported biological conclusion/classification'
    for field in ('true_absence_inferred', 'alternative_maintenance_inferred'):
        equal(adequacy[field], False, field)
    for field in ('absence_inferred', 'alt_inferred'):
        equal(adequacy_row[field], 'False', f'adequacy {field}')
    if adequacy['status'] == 'conditional_power_supported':
        assumptions = adequacy['assumptions']
        assert all(field in assumptions for field in REQUIRED_ASSUMPTIONS), 'power supported with missing assumptions'
        assert all(str(assumptions.get('provenance', {}).get(field, '')).strip() for field in REQUIRED_ASSUMPTIONS), 'power supported without sourced assumptions'
        assert assumptions['random_wgs_assumed'] is True, 'uniform WGS assumption not explicit'
        assert not adequacy.get('missing_assumptions'), 'power supported despite missing assumptions'
        assert not adequacy.get('missing_assumption_provenance'), 'power supported despite missing provenance'
        assert not adequacy.get('library', {}).get('uniform_wgs_model_blockers'), 'power supported despite library blockers'
        assert adequacy.get('power'), 'power supported without calculation'
    if document['classification'] == 'analysis_failed':
        equal(document['analysis_complete'], False, 'failed completion')
        equal(document['top_candidate'], None, 'failed top candidate')
        equal(adequacy['status'], 'observed_qc_failure', 'failed adequacy')
    check_observation_contract(document, raw, evidence)


def check_identity(document, metadata, assignment, taxon_rows, args):
    equal(document['reference_assignment'], assignment, 'complete assignment')
    for field in ('scientific_name', 'tax_id', 'lineage'):
        equal(str(document['metadata'].get(field, '')).strip(), str(metadata.get(field, '')).strip(), f'metadata {field}')
        equal(str(assignment.get(field, '')).strip(), str(metadata.get(field, '')).strip(), f'assignment {field}')
    equal(document['metadata']['status'], metadata['status'], 'taxonomy status preserved')
    provenance = document['provenance']
    for field, file in (('catalogue_sha256', args.catalog), ('assignment_sha256', args.assignments),
                        ('taxon_reference_sha256', args.taxon_references)):
        digest_matches(file, provenance[field], field)
    if assignment.get('reference_status') == 'supported':
        ids = [key for key in assignment.get('reference_ids', '').split(',') if key]
        assert ids, 'supported assignment has no reference IDs'
        motifs = set()
        for key in ids:
            row = taxon_rows[key]
            equal(row['reference_status'], 'supported', f'reference status {key}')
            assert row.get('evidence') and row.get('source'), f'undocumented supported reference {key}'
            motifs.update(canonical(m) for m in row.get('related_motifs', '').split(',') if m)
        expected = {canonical(m) for m in assignment.get('related_motifs', '').split(',') if m}
        assert expected <= motifs, 'assignment families absent from referenced table rows'


def check_reference_evidence(document, catalogue_rows, taxon_rows):
    """Rebuild exact-family reference certainty from cited input tables."""
    if not observation_schema(document):
        return 'legacy schema: retained original reference audit'
    assignment = document['reference_assignment']
    by_family = {}
    for row in catalogue_rows:
        by_family.setdefault(canonical(row['motif']), []).append(row)
    if assignment.get('reference_status') in {'supported', 'candidate'}:
        supported_ids = set(filter(None, assignment.get('reference_ids', '').split(',')))
        candidate_ids = set(filter(None, assignment.get('review_reference_ids', '').split(',')))
        supported_families = {canonical(m) for m in assignment.get('related_motifs', '').split(',') if m}
        candidate_families = {canonical(m) for m in assignment.get('review_motifs', '').split(',') if m}
        for key in supported_ids | candidate_ids:
            assert key in taxon_rows, 'assignment references an absent reference ID'
            row = taxon_rows[key]
            supported = (key in supported_ids and assignment.get('reference_status') == 'supported'
                         and row.get('reference_status') == 'supported' and row.get('evidence') in ESTABLISHED_EVIDENCE)
            reviewed = (key in candidate_ids and row.get('reference_status') == 'candidate'
                        and row.get('evidence') and row.get('source', '').startswith(('https://', 'http://')))
            if not (supported or reviewed) or row.get('repeat_architecture') not in {'exact_short_unit', 'mixed_short_units'}:
                continue
            for motif in row.get('related_motifs', '').split(','):
                if not 3 <= len(motif) <= 32 or set(motif) - set('ACGT'):
                    continue
                family = canonical(motif)
                if family in (supported_families if supported else candidate_families):
                    by_family.setdefault(family, []).append(dict(row, motif=motif))
    for candidate in document['candidates']:
        rows = by_family.get(candidate['canonical_family'], [])
        established = any(row.get('evidence') in ESTABLISHED_EVIDENCE and row.get('reference_status') != 'candidate' for row in rows)
        expected = 'established' if established else 'published_candidate' if rows else 'not_in_catalogue'
        equal(candidate['reference_status'], expected, 'reference certainty from exact cited family')
        equal(candidate['reference_sources'], sorted({row['source'] for row in rows}), 'exact reference sources')
        equal(candidate['reference_evidence'], sorted({row['evidence'] for row in rows}), 'exact reference evidence')
        equal(candidate['published_motifs'], sorted({row['motif'] for row in rows}), 'exact published motif displays')
    return {'exact_reference_families': len(by_family)}


def check_provenance(document, source_path, raw, status, evidence_path, args, decision):
    provenance = document['provenance']
    equal(Path(provenance['original_discovery_result']).resolve(), source_path.resolve(), 'discovery path')
    if not source_path.is_file():
        equal(provenance['original_discovery_sha256'], None, 'missing source hash')
        equal(status['status'], 'source_failed', 'missing source status')
        equal(document['classification'], 'analysis_failed', 'source failure is not biological negative')
        return 'source_failed: no biological inference'
    digest_matches(source_path, provenance['original_discovery_sha256'], 'discovery hash')
    digest_matches(source_path, decision['sha256'], 'frozen decision-audit discovery hash')
    equal(Path(decision['path']).resolve(), source_path.resolve(), 'decision audit source path')
    source_state = read_json(source_path.parent / 'state.json')
    digest_matches(source_path, source_state['result_sha256'], 'discovery state hash')
    prior = read_json(source_path.parent / 'provenance.json')
    digest_matches(source_path, prior['result_sha256'], 'prior result hash')
    raw_path = Path(prior['source_result_json'])
    digest_matches(raw_path, prior['original_result_sha256'], 'original raw discovery hash')
    digest_matches(raw_path, prior['source_state']['result_sha256'], 'embedded original raw state hash')
    digest_matches(raw_path, read_json(raw_path.parent / 'state.json')['result_sha256'], 'original raw state hash')
    original = read_json(raw_path)
    for field in ('records', 'bases'):
        equal(raw['summary'][field], original['summary'][field], f'raw/prior {field}')
    raw_families = {canonical(c['motif']) for c in raw.get('candidates', []) + raw.get('rejected_candidates', [])}
    original_families = {canonical(c['motif']) for c in original.get('candidates', []) + original.get('rejected_candidates', [])}
    equal(raw_families, original_families, 'all raw families preserved')
    equal(provenance['all_original_repeat_families_retained_in_source'], len(raw_families), 'retained-family count')
    equal(status['status'], 'complete', 'final evidence completed')
    equal(document['analysis_complete'], True, 'completed analysis')
    digest_matches(source_path, status['identity']['source_sha256'], 'evidence source hash')
    # A view may select saved evidence from several collector versions. Bind
    # each row to its original executable; never substitute the view's binary.
    evidence_binary = evidence_path.resolve().parent.parent / 'trm-evidence'
    digest_matches(evidence_binary, status['identity']['binary_sha256'], 'frozen evidence executable hash')
    digest_matches(evidence_path, status['result_sha256'], 'completed evidence hash')
    digest_matches(evidence_path, provenance['targeted_evidence_sha256'], 'document evidence hash')
    equal(Path(provenance['targeted_evidence_result']).resolve(), evidence_path.resolve(), 'evidence path')
    evidence = read_json(evidence_path)
    equal(status['identity']['max_records'], 0, 'no prefix cap')
    equal(evidence['summary']['sampled_prefix'], False, 'full evidence input')
    equal(raw['summary']['sampled_prefix'], False, 'full original input')
    for field in ('records', 'bases'):
        equal(evidence['summary'][field], raw['summary'][field], f'full-scan {field}')
        equal(status[field], raw['summary'][field], f'completed status {field}')
        equal(evidence['evidence_qc'][field], raw['summary'][field], f'evidence QC {field}')
    targets = status['identity']['targets']
    equal(read_json(evidence_path.parent / 'targets.json'), targets, 'saved target selection')
    equal(set(targets['motifs']), set(evidence['target_motifs']), 'requested evidence motifs')
    equal(set(targets['motifs']), {c['motif'] for c in evidence['candidate_evidence']}, 'output evidence families')


def run_audit(args):
    audit = Audit()
    try:
        summaries = index_rows(read_tsv(args.assessment / 'samples.tsv'))
        comparisons = index_rows(read_tsv(args.comparison / 'comparison.tsv'))
        candidates = read_tsv(args.assessment / 'candidates.tsv')
        adequacy = index_rows(read_tsv(args.assessment / 'adequacy.tsv'))
        metadata = index_rows(read_tsv(args.metadata))
        assignments = index_rows(read_tsv(args.assignments))
        taxon_rows = index_rows(read_tsv(args.taxon_references), 'reference_id')
        catalogue_rows = read_tsv(args.catalog)
        metrics = read_json(args.assessment / 'metrics.json')
        batch_metrics = read_json(args.evidence / 'metrics.json')
        decision = index_rows(read_json(args.decision_audit)['input_result_sha256'])
        pilot = index_rows(read_json(args.pilot_check)['samples'])
        documents = {path.stem: read_json(path) for path in sorted((args.assessment / 'samples').glob('*.json'))}
    except (OSError, ValueError, KeyError, AssertionError) as error:
        audit.checks.append({'check': 'required_artifacts', 'status': 'failed', 'detail': str(error)})
        return audit.finish(args)

    def sample_sets():
        for label, rows in (('comparison', comparisons), ('adequacy', adequacy), ('metadata', metadata),
                            ('assignments', assignments), ('sample JSON', documents)):
            equal(set(rows), set(summaries), label + ' sample set')
        equal(len(summaries), args.expected_samples, 'expected sample count')
        assert {r['sample_id'] for r in candidates} <= summaries.keys(), 'candidate sample absent from summary'
        equal(set(decision), {s for s in summaries if (args.source / s / 'result.json').is_file()}, 'frozen discovery hash sample set')
    audit.run('sample_sets_and_counts', sample_sets)
    grouped = {sample: [] for sample in summaries}
    for row in candidates:
        grouped.setdefault(row['sample_id'], []).append(row)
    failed_sources = []
    for sample, document in documents.items():
        if sample not in summaries or sample not in comparisons or sample not in metadata or sample not in assignments or sample not in adequacy:
            continue
        source_path = args.source / sample / 'result.json'
        try:
            raw = read_json(source_path) if source_path.is_file() else {'summary': {}, 'candidates': []}
            status = read_json(args.evidence / sample / 'status.json')
            evidence_file = args.evidence / sample / 'evidence.json'
            evidence = read_json(evidence_file) if evidence_file.is_file() else None
        except (OSError, ValueError) as error:
            audit.checks.append({'check': 'source_artifacts', 'sample_id': sample, 'status': 'failed', 'detail': str(error)})
            continue
        if status.get('status') == 'source_failed':
            failed_sources.append(sample)
        audit.run('sample_table_json_consistency', lambda: check_sample_tables(document, summaries[sample], comparisons[sample], grouped[sample], adequacy[sample], raw, args.assessment / 'samples' / f'{sample}.json'), sample)
        audit.run('metadata_assignment_reference_integrity', lambda: check_identity(document, metadata[sample], assignments[sample], taxon_rows, args), sample)
        audit.run('reference_evidence_certainty', lambda: check_reference_evidence(document, catalogue_rows, taxon_rows), sample)
        audit.run('original_and_evidence_provenance', lambda: check_provenance(document, source_path, raw, status, args.evidence / sample / 'evidence.json', args, decision.get(sample)), sample)
        audit.run('no_absence_or_unjustified_power_upgrade', lambda: check_interpretation(document, adequacy[sample], raw, evidence), sample)
        audit.run('sample_html_links', lambda: check_links(args.assessment / 'samples' / f'{sample}.html',
                                                         (f'{sample}.json', '../index.html')), sample)

    def aggregate_metrics():
        equal(metrics['samples'], len(documents), 'metric samples')
        equal(metrics['candidate_rows'], len(candidates), 'metric candidate rows')
        equal(metrics['classifications'], dict(Counter(d['classification'] for d in documents.values())), 'classification counts')
        equal(metrics['data_adequacy_counts'], dict(Counter(d['data_adequacy']['status'] for d in documents.values())), 'adequacy counts')
        equal(metrics['top_family_changes'], sum(r['motif_family_changed'] == 'True' for r in comparisons.values()), 'family changes')
        equal(metrics['novel_candidate_rows_including_secondary'], sum(c['novel_candidate'] for d in documents.values() for c in d['candidates']), 'novel candidate count')
        for name, path in (('samples.tsv', args.assessment / 'samples.tsv'), ('candidates.tsv', args.assessment / 'candidates.tsv'), ('comparison.tsv', args.comparison / 'comparison.tsv')):
            digest_matches(path, metrics['files'][name], 'metric ' + name)
        equal(len(failed_sources), args.expected_source_failures, 'source failure count')
        equal(batch_metrics['sample_count'], len(documents), 'evidence batch samples')
        equal(batch_metrics['source_failed'], len(failed_sources), 'evidence source failures')
        completed = [d for d in documents.values() if d['analysis_complete']]
        equal(batch_metrics['completed'], len(completed), 'evidence completed samples')
        equal(batch_metrics['failed'], 0, 'no incomplete/failed evidence scans')
        for field in ('records', 'bases'):
            equal(batch_metrics[field], sum(d['input_summary'][field] for d in completed), 'evidence aggregate ' + field)
        return {'samples': len(documents), 'candidates': len(candidates), 'source_failed_samples': failed_sources}
    audit.run('aggregate_metrics_and_source_failures', aggregate_metrics)
    audit.run('summary_html_links', lambda: check_links(args.assessment / 'index.html',
                                                       [row['details'] for row in summaries.values()]))
    audit.run('comparison_html_links', lambda: check_links(args.comparison / 'index.html'))

    def biological_regressions():
        equal(len(pilot), args.expected_pilot_samples, 'pilot sample count')
        for sample, expected in pilot.items():
            motif = expected['detected_motif']
            assert motif in ('TTAGAC', 'TTAGGT') or args.expected_pilot_samples != 6, 'unexpected pilot reference'
            document = documents[sample]
            equal(document['detected_motif'], motif, f'pilot {sample} observed motif')
            assert document['classification'] in KNOWN_CLASSIFICATIONS, f'{sample} known family lost'
            equal(document['top_candidate']['reference_status'], 'established', 'pilot reference status')
            assert not document['top_candidate']['novel_candidate'], 'pilot family called novel'
        hapla = [d for d in documents.values() if d['scientific_name'].strip() == 'Meloidogyne hapla']
        equal(len(hapla), args.expected_mhapla_samples, 'M. hapla sample count')
        family = canonical('CCCAAGGTTTAAAAGG')
        for document in hapla:
            candidate = next((c for c in document['candidates'] if c['canonical_family'] == family), None)
            assert candidate, 'M. hapla published 16-bp family absent from candidate report'
            equal(candidate['reference_status'], 'established', 'M. hapla reference status')
            equal(candidate['novel_candidate'], False, 'M. hapla cannot be novel')
        return {'pilot_samples': sorted(pilot), 'M_hapla_samples': [d['sample_id'] for d in hapla]}
    audit.run('known_family_reporting_regressions', biological_regressions)
    return audit.finish(args)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for flag, default in (
        ('assessment', 'results/wgs_assessed_v040'), ('comparison', 'results/comparison_v030_v040'),
        ('source', 'results/wgs_taxonomy_v030'), ('evidence', 'results/evidence_v040_20260921'),
        ('metadata', 'data/sample_metadata_20260920/sample_metadata.tsv'),
        ('assignments', 'results/priors_v040_20260921/assignments.tsv'),
        ('catalog', 'data/telomere_motifs.tsv'), ('taxon-references', 'data/taxon_motif_references.tsv'),
        ('decision-audit', 'results/motif_review_20260920/decision_audit.json'),
        ('pilot-check', 'results/evidence_v040_pilot/interpretation_check.json')):
        parser.add_argument('--' + flag, type=Path, default=Path(default))
    for flag, default in (('expected-samples', 230), ('expected-source-failures', 2),
                          ('expected-pilot-samples', 6), ('expected-mhapla-samples', 2)):
        parser.add_argument('--' + flag, type=int, default=default)
    parser.add_argument('--output', type=Path, help='Write audit JSON here as well as stdout')
    return parser.parse_args(argv)


def main():
    args = parse_args()
    result = run_audit(args)
    text = json.dumps(result, indent=2, ensure_ascii=False) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end='')
    raise SystemExit(0 if result['passed'] else 1)


if __name__ == '__main__':
    main()
