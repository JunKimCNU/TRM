#!/usr/bin/env python3
"""Explicitly executed, pinned two-family full-input Volvariella diagnostic.

Default mode validates metadata only. No downloads, input rewriting, alternate
input discovery, prefix sampling, or modification of original assessments.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
sys.path.insert(0, str(PROJECT / 'scripts'))
from evidence_batch import canonical, resources
from expected_diagnostics import supplement_expected_observations, DETECTION_FIELDS

SAMPLE = 'SRR1283432'
TARGETS = ['TTAGGG', 'ACACCCACCTCCACC']


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''): digest.update(block)
    return digest.hexdigest()


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha(path)}


def checked(item):
    path = Path(item['path'])
    if path.stat().st_size != item['bytes'] or sha(path) != item['sha256']:
        raise ValueError('Pinned file changed: ' + str(path))
    return path


def fingerprint(path):
    stat = Path(path).stat()
    return {'path': str(Path(path).resolve()), 'device': stat.st_dev, 'inode': stat.st_ino,
            'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns, 'ctime_ns': stat.st_ctime_ns}


def atomic(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    os.replace(temporary, path)


def utc(): return datetime.now(timezone.utc).isoformat()


def preflight(config_path, config_sha256):
    if sha(config_path) != config_sha256:
        raise ValueError('Configuration SHA256 differs')
    config = json.loads(Path(config_path).read_text())
    if config['sample_id'] != SAMPLE or config['motifs'] != TARGETS or config['max_records'] != 0:
        raise ValueError('Only the two registered full-input Volvariella targets are allowed')
    if not 1 <= config['threads'] <= 8 or not 512 <= config['memory_limit_mib'] <= 65536:
        raise ValueError('Diagnostic resource budget exceeds registration')
    if config.get('downloads_enabled') is not False:
        raise ValueError('Downloads are not allowed')
    checked(config['runner'])
    for source in config['dependencies']: checked(source)
    source = checked(config['source']); assessment_path = checked(config['assessment']); binary = checked(config['binary'])
    raw = json.loads(source.read_text())
    assessment = json.loads(assessment_path.read_text())
    if (assessment.get('sample_id') != SAMPLE or assessment.get('detected_motif') != TARGETS[1]
            or assessment.get('provenance', {}).get('original_discovery_sha256') != config['source']['sha256']):
        raise ValueError('Assessment sample, representative or source identity differs')
    if raw['version'] != config['detector_version']:
        raise ValueError('Native detector version differs')
    version = subprocess.run([str(binary), '--version'], capture_output=True, text=True, check=True).stdout.strip()
    if version != 'trm ' + raw['version']:
        raise ValueError('Executable detector version differs')
    if raw['summary'].get('sampled_prefix') is not False or any(
            raw['summary'].get(key) != config['expected_summary'][key] for key in ('records', 'bases')):
        raise ValueError('Source is not the registered full input')
    if raw['config'] != config['native_config'] or raw['junction_options'] != config['junction_options']:
        raise ValueError('Original detector settings changed')
    if raw['inputs'] != config['native_inputs'] or len(raw['inputs']) != 2 or len(config['inputs']) != 2:
        raise ValueError('Paired input identity differs')
    if (raw.get('paired_end') or {}).get('layout') != 'fr':
        raise ValueError('Source is not the registered paired FR analysis')
    for stored, identity in zip(config['inputs'], raw['inputs']):
        current = fingerprint(stored['path'])
        if current != stored or current['path'] != identity['path'] or current['bytes'] != identity['bytes'] or current['mtime_ns'] // 1000000000 != identity['modified_unix_seconds']:
            raise ValueError('Full-input archive fingerprint changed')
    return config, raw


def command(config, directory):
    argv = [config['binary']['path'], 'evidence', *[row['path'] for row in config['native_inputs']],
            '--targets', str(directory / 'targets.json'), '--output', str(directory / 'evidence'),
            '--threads', str(config['threads']), '--memory-limit-mib', str(config['memory_limit_mib']),
            '--max-records', '0', '--paired-layout', 'fr']
    # Target families replace the discovery min/max-period search range.
    for field in ('mode', 'seed_bases', 'min_copies', 'min_identity', 'min_tract', 'min_reads', 'end_window', 'max_candidates'):
        argv += ['--' + field.replace('_', '-'), str(config['native_config'][field])]
    argv += ['--junction-flank-bases', str(config['junction_options']['flank_bases']),
             '--junction-end-slack', str(config['junction_options']['end_slack']),
             '--flank-period-annotation', 'on' if config['junction_options']['flank_period_annotation'] else 'off']
    return argv


def validate_join(config, raw, directory, status):
    assessment = json.loads(Path(config['assessment']['path']).read_text())
    # Validate through the production join before publishing a completed status.
    with tempfile.TemporaryDirectory(prefix='join-validation-', dir=directory) as temporary:
        folder = Path(temporary) / SAMPLE; folder.mkdir()
        os.link(directory / 'evidence.json', folder / 'evidence.json')
        (folder / 'status.json').write_text(json.dumps({**status, 'status': 'complete'}))
        observed, provenance = supplement_expected_observations(SAMPLE,
            assessment['expected_family_observations'], temporary, raw,
            config['source']['sha256'], analysis_complete=True)
    return {'join_accepted': True, 'registered_families': provenance['registered_families'],
            'used_expected_families': provenance['used_expected_families'],
            'observations': observed}


def execute(config_path, config_sha256, *, run=subprocess.run):
    config, raw = preflight(config_path, config_sha256)
    directory = Path(config['output_root']) / SAMPLE
    directory.mkdir(parents=True, exist_ok=False)
    status = {'sample_id': SAMPLE, 'status': 'running', 'started_utc': utc(),
              'operation': 'registered_full_input_expected_family_diagnostic',
              'config_sha256': config_sha256, 'downloads_started': 0,
              'identity': {'source_sha256': config['source']['sha256'], 'motifs': [canonical(m) for m in TARGETS],
                           'inputs': config['native_inputs'], 'binary_sha256': config['binary']['sha256'],
                           'detector_version': config['detector_version'], 'max_records': 0,
                           'threads': config['threads'], 'memory_limit_mib': config['memory_limit_mib'],
                           'native_config': config['native_config'], 'junction_options': config['junction_options']}}
    atomic(directory / 'status.json', status)
    started = time.monotonic()
    try:
        atomic(directory / 'targets.json', {'motifs': [canonical(m) for m in TARGETS]})
        argv = command(config, directory)
        status['command'] = argv
        atomic(directory / 'status.json', status)
        timed = ['/usr/bin/time', '-f', '%e\t%U\t%S\t%M\t%x', '-o', str(directory / 'resource.tsv'), *argv]
        with (directory / 'stdout.log').open('wb') as stdout, (directory / 'stderr.log').open('wb') as stderr:
            completed = run(timed, stdout=stdout, stderr=stderr, check=False)
        status.update(exit_code=completed.returncode, wall_seconds=time.monotonic() - started)
        status.update(resources(directory / 'resource.tsv'))
        if completed.returncode or status['resource_exit_code'] != completed.returncode:
            raise RuntimeError('Diagnostic process failed; inspect stderr.log')
        # Check both source identities and all input fingerprints after the full pass.
        preflight(config_path, config_sha256)
        evidence = json.loads((directory / 'evidence.json').read_text())
        status['result_sha256'] = sha(directory / 'evidence.json')
        status['result_json'] = str(directory / 'evidence.json')
        validation = validate_join(config, raw, directory, status)
        if evidence.get('paired_end', {}).get('layout') != 'fr' or evidence.get('paired_end', {}).get('scope') != 'complete_inputs':
            raise ValueError('Paired evidence does not cover complete inputs')
        representative = canonical(TARGETS[1])
        native = next(row for row in raw['candidates'] + raw['rejected_candidates'] if canonical(row['motif']) == representative)
        current = next(row for row in evidence['candidate_evidence'] if canonical(row['motif']) == representative)
        validation['representative_control'] = {key: {'original': native[key], 'targeted': current[key],
            'equal': native[key] == current[key]} for key in ('reads', 'full_repeat_reads', 'repeat_bases', 'max_tract')}
        validation['source_inputs_unchanged'] = True
        atomic(directory / 'validation.json', validation)
        status.update(status='complete', finished_utc=utc(), wall_seconds=time.monotonic() - started,
                      full_input_records=evidence['summary']['records'], full_input_bases=evidence['summary']['bases'],
                      expected_diagnostics_join_accepted=True, validation=pin(directory / 'validation.json'))
        atomic(directory / 'status.json', status)
        return status
    except BaseException as error:
        status.update(status='failed', error=type(error).__name__ + ': ' + str(error),
                      finished_utc=utc(), wall_seconds=time.monotonic() - started)
        atomic(directory / 'status.json', status)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=HERE / 'volvariella_config.json')
    parser.add_argument('--config-sha256', required=True)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    if args.execute:
        value = execute(args.config, args.config_sha256)
    else:
        config, raw = preflight(args.config, args.config_sha256)
        value = {'status': 'preflight_passed_no_scan', 'sample_id': SAMPLE,
                 'records': raw['summary']['records'], 'bases': raw['summary']['bases'],
                 'command': command(config, Path(config['output_root']) / SAMPLE)}
    print(json.dumps(value, indent=2))


if __name__ == '__main__': main()
