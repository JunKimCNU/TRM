"""Resolve the current managed scanner without changing historical executables.

Precedence: explicit argument, TRM_BINARY, then <project>/bin/trm. Only the
managed default is tied to Cargo.toml and its installation receipt. Explicit
overrides remain available for older reproducibility runs and test doubles.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def project_version(root: Path = ROOT) -> str:
    text = (Path(root) / 'Cargo.toml').read_text(encoding='utf-8')
    package = re.search(r'(?ms)^\[package\]\s*\n(.*?)(?=^\[|\Z)', text)
    version = re.search(r'(?m)^version\s*=\s*"([^"]+)"\s*$', package[1] if package else '')
    if not version:
        raise ValueError('Cannot read the current scanner version from Cargo.toml')
    return version[1]


def binary_version(path: Path) -> str:
    try:
        result = subprocess.run([str(path), '--version'], check=True, text=True,
                                capture_output=True, timeout=10)
    except (OSError, subprocess.SubprocessError) as error:
        raise ValueError(f'Cannot identify scanner version: {path}: {error}') from error
    match = re.fullmatch(r'trm ([^\s]+)', result.stdout.strip())
    if not match:
        raise ValueError(f'Unexpected scanner version response: {path}: {result.stdout.strip()!r}')
    return match[1]


def installed_identity(root: Path = ROOT) -> dict:
    root = Path(root).resolve()
    binary = root / 'bin' / 'trm'
    receipt = root / 'bin' / 'runtime.json'
    install_hint = f'Install the current runtime with: python3 {root / "scripts/install_runtime.py"}'
    try:
        if binary.parent.is_symlink() or binary.is_symlink() or receipt.is_symlink():
            raise ValueError('Managed runtime files must not be symlinks')
        if not binary.is_file() or not os.access(binary, os.X_OK):
            raise ValueError(f'Current managed scanner is absent or not executable: {binary}')
        recorded = json.loads(receipt.read_text(encoding='utf-8'))
        if (not isinstance(recorded, dict) or recorded.get('schema_version') != 1 or recorded.get('binary') != 'trm'
                or not re.fullmatch('[0-9a-f]{64}', str(recorded.get('sha256', '')))):
            raise ValueError(f'Invalid runtime receipt: {receipt}')
        expected_version = project_version(root)
        if recorded.get('version') != expected_version:
            raise ValueError(f'Managed scanner version {recorded.get("version")} differs from project {expected_version}')
        before = binary.stat()
        actual_hash = sha256(binary)
        if actual_hash != recorded['sha256']:
            raise ValueError(f'Managed scanner SHA256 mismatch: {binary}')
        actual_version = binary_version(binary)
        after = binary.stat()
        identity = lambda info: (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns)
        if identity(before) != identity(after):
            raise ValueError(f'Managed scanner changed during verification: {binary}')
        if actual_version != expected_version:
            raise ValueError(f'Managed scanner reports {actual_version}; expected {expected_version}')
        return {'path': str(binary), 'version': actual_version, 'sha256': actual_hash,
                'receipt': str(receipt), 'selection': 'managed_runtime'}
    except (OSError, ValueError, TypeError) as error:
        raise ValueError(f'{error}. {install_hint}. Historical target/release/trm is never selected automatically.') from error


def resolve_binary(explicit: str | Path | None = None, *, root: Path = ROOT,
                   environ: dict | None = None) -> Path:
    """Resolve an executable; never silently fall back after an invalid choice."""
    env = os.environ if environ is None else environ
    selected = explicit if explicit is not None else env.get('TRM_BINARY')
    if selected is None:
        return Path(installed_identity(root)['path'])
    value = os.fspath(selected)
    if not value.strip():
        raise ValueError('The explicit scanner/TRM_BINARY must not be empty')
    expanded = os.path.expanduser(value)
    executable = shutil.which(expanded, path=env.get('PATH'))
    if executable is None:
        raise ValueError(f'Scanner is not an executable file or PATH command: {value}')
    path = Path(executable).resolve()
    managed = Path(root).resolve() / 'bin' / 'trm'
    if path == managed or Path(expanded).absolute() == managed:
        return Path(installed_identity(root)['path'])
    return path


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', help='Explicit scanner override; older versions are permitted')
    args = parser.parse_args()
    try:
        binary = resolve_binary(args.binary)
        print(json.dumps({'path': str(binary), 'version': binary_version(binary), 'sha256': sha256(binary)}, indent=2))
    except ValueError as error:
        parser.error(str(error))
