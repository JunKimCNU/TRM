"""Shared native selection for unittest discovery and individual test files."""
import importlib.util
from pathlib import Path

_path = Path(__file__).resolve().parents[1] / 'scripts/runtime_binary.py'
_spec = importlib.util.spec_from_file_location('_trm_runtime_for_tests', _path)
_runtime = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_runtime)


def runtime_binary():
    return _runtime.resolve_binary()
