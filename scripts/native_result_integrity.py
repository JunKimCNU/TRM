"""Bind native evidence to an available batch/reclassification completion record.

Standalone native JSON has no batch state; retain that supported input format
without claiming that its execution receipt has been verified. A present state
is authoritative and may never be ignored in favor of an older result file.
"""
from __future__ import annotations

import json
from pathlib import Path


def validate_native_result(source, sample, result_sha256):
    source = Path(source)
    state_path = source.with_name('state.json')
    base = {'state_path': str(state_path.resolve()),
            'result_sha256': result_sha256,
            'execution_receipt_verified': False}
    if not state_path.exists() and not state_path.is_symlink():
        return {**base, 'status': 'standalone_state_unavailable' if source.is_file() else 'missing_result',
                'analysis_eligible': source.is_file(),
                'reason': 'No batch state supplied; only native document completeness can be evaluated.'}
    if not state_path.is_file():
        raise ValueError(f'{sample}: native state is not a readable file: {state_path}')
    state = json.loads(state_path.read_text(encoding='utf-8'))
    if not isinstance(state, dict) or not isinstance(state.get('status'), str):
        raise ValueError(f'{sample}: native state lacks a valid status')
    if state.get('sample_id') is not None and state['sample_id'] != sample:
        raise ValueError(f'{sample}: native state sample identity differs')
    if state['status'] != 'complete':
        return {**base, 'status': 'source_not_complete', 'analysis_eligible': False,
                'recorded_status': state['status'],
                'reason': 'Native execution is ' + state['status'] + '; retained previous results are not current successful evidence.'}
    expected = state.get('result_sha256')
    if not source.is_file() or not isinstance(expected, str) or expected != result_sha256:
        raise ValueError(f'{sample}: native result hash differs from recorded completed scan')
    return {**base, 'status': 'completed_state_verified', 'analysis_eligible': True,
            'execution_receipt_verified': True, 'recorded_status': 'complete',
            'reason': 'Native result bytes match the available completed execution record.'}
