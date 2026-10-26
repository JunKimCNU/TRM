#!/usr/bin/env python3
"""Install the current scanner at bin/trm without touching target/release/trm.

With no --from-binary, build from Cargo.lock in a separate temporary target.
An existing installation with different bytes is refused, never overwritten.
Historical provenance can be embedded using --provenance; it is informational,
and the installed runtime does not depend on that original path remaining.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

from runtime_binary import ROOT, binary_version, installed_identity, project_version, sha256


def install(source: Path, *, root: Path = ROOT, expected_sha256: str | None = None,
            provenance: dict | None = None) -> dict:
    root, source = Path(root).resolve(), Path(source).resolve(strict=True)
    version = project_version(root)
    digest = sha256(source)
    if expected_sha256 is not None and digest != expected_sha256:
        raise ValueError('Source scanner SHA256 differs from --expected-sha256')
    if not source.is_file() or not os.access(source, os.X_OK) or binary_version(source) != version:
        raise ValueError(f'Source scanner must be executable and report trm {version}')
    bindir = root / 'bin'
    if bindir.is_symlink():
        raise ValueError('Managed bin directory must not be a symlink')
    bindir.mkdir(parents=True, exist_ok=True)
    binary, manifest = bindir / 'trm', bindir / 'runtime.json'
    lockpath = bindir / '.install.lock'
    if lockpath.is_symlink():
        raise ValueError('Runtime installation lock must not be a symlink')
    with lockpath.open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if binary.exists() or binary.is_symlink() or manifest.exists() or manifest.is_symlink():
            current = installed_identity(root)
            if current['sha256'] != digest:
                raise ValueError('A different managed runtime already exists; preserve it before explicitly installing another release')
            return {**current, 'status': 'already_installed'}
        record = {'schema_version': 1, 'binary': 'trm', 'version': version, 'sha256': digest,
                  'bytes': source.stat().st_size, 'installed_utc': datetime.now(timezone.utc).isoformat(),
                  'source_binary': str(source), 'source_binary_sha256': digest,
                  'build_provenance': provenance,
                  'provenance_scope': 'Source paths are historical annotations, not runtime dependencies'}
        temporary = []
        try:
            fd, name = tempfile.mkstemp(prefix='.trm-', dir=bindir)
            staged = Path(name); temporary.append(staged)
            with os.fdopen(fd, 'wb') as out, source.open('rb') as original:
                shutil.copyfileobj(original, out, 1024 * 1024)
                out.flush(); os.fsync(out.fileno())
            staged.chmod(0o755)
            if sha256(staged) != digest or sha256(source) != digest or binary_version(staged) != version:
                raise ValueError('Source scanner changed during installation')
            fd, name = tempfile.mkstemp(prefix='.runtime-', dir=bindir)
            staged_manifest = Path(name); temporary.append(staged_manifest)
            with os.fdopen(fd, 'w', encoding='utf-8') as out:
                json.dump(record, out, indent=2, sort_keys=True); out.write('\n')
                out.flush(); os.fsync(out.fileno())
            # Any crash between the two replacements fails default validation
            # closed; it cannot cause a fallback to the preserved old scanner.
            os.replace(staged, binary)
            os.replace(staged_manifest, manifest)
            directory_fd = os.open(bindir, os.O_RDONLY)
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
        finally:
            for path in temporary:
                if path.exists():
                    path.unlink()
        return {**installed_identity(root), 'status': 'installed'}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--from-binary', type=Path, help='Install an already built scanner instead of building')
    parser.add_argument('--expected-sha256', help='Required source digest when supplied')
    parser.add_argument('--provenance', type=Path, help='Embed the existing build receipt without changing it')
    parser.add_argument('--offline', action='store_true', help='Use only already cached Cargo dependencies')
    parser.add_argument('--jobs', type=int, default=min(8, os.cpu_count() or 1), help='Cargo build workers (default: up to 8)')
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    try:
        if args.from_binary:
            provenance = json.loads(args.provenance.read_text()) if args.provenance else None
            result = install(args.from_binary, root=ROOT, expected_sha256=args.expected_sha256, provenance=provenance)
        else:
            if args.provenance:
                parser.error('--provenance is only for --from-binary; builds record their own provenance')
            with tempfile.TemporaryDirectory(prefix='trm-current-build-') as builddir:
                command = ['cargo', 'build', '--release', '--locked', '--target-dir', builddir, '--jobs', str(args.jobs)]
                if args.offline:
                    command.append('--offline')
                source_paths = [ROOT/'Cargo.toml', ROOT/'Cargo.lock', *sorted((ROOT/'src').glob('*.rs'))]
                inputs = {str(p.relative_to(ROOT)): sha256(p) for p in source_paths}
                subprocess.run(command, cwd=ROOT, check=True)
                if inputs != {str(p.relative_to(ROOT)): sha256(p) for p in source_paths}:
                    raise ValueError('Build sources changed while compiling; runtime was not installed')
                provenance = {'command': command, 'source_sha256': inputs,
                              'rustc': subprocess.check_output(['rustc', '--version'], text=True).strip(),
                              'cargo': subprocess.check_output(['cargo', '--version'], text=True).strip()}
                result = install(Path(builddir)/'release/trm', root=ROOT, expected_sha256=args.expected_sha256, provenance=provenance)
        print(json.dumps(result, indent=2))
        return 0
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        parser.error(str(error))


if __name__ == '__main__':
    raise SystemExit(main())
