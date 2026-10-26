#!/usr/bin/env python3
"""Run Python regressions with one explicitly verified scanner selection.

Use --binary for a development or historical build. The normal invocation uses
the current managed runtime, never the preserved target/release/trm artifact.
This runner does not build or rewrite any native binary.
"""
import argparse
import os
import sys
import unittest

from runtime_binary import ROOT, binary_version, resolve_binary, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary')
    parser.add_argument('--pattern', default='test_*.py')
    parser.add_argument('--quiet', action='store_true')
    parser.add_argument('--include-windows-tests', action='store_true',
                        help='Explicitly include the separately maintained Windows downloader tests')
    args = parser.parse_args()
    try:
        binary = resolve_binary(args.binary)
        version = binary_version(binary)
    except ValueError as error:
        parser.error(str(error))
    print(f'Python regressions: trm {version}; {binary}; SHA256 {sha256(binary)}', flush=True)
    # Import only selected modules: filtering an already-discovered suite would
    # still import the separate Windows project's tests and their dependencies.
    os.environ['TRM_BINARY'] = str(binary)
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    os.chdir(ROOT)
    sys.path.insert(0, str(ROOT / 'tests'))
    selected = sorted((ROOT / 'tests').glob(args.pattern))
    selected = [p for p in selected if p.is_file() and p.suffix == '.py'
                and (args.include_windows_tests or not p.name.startswith('test_windows'))]
    if not selected:
        parser.error('No test modules match the selected pattern')
    loader = unittest.TestLoader()
    suite = unittest.TestSuite(loader.loadTestsFromName(p.stem) for p in selected)
    result = unittest.TextTestRunner(verbosity=1 if args.quiet else 2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    raise SystemExit(main())
