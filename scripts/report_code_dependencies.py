"""Conservative transitive local import closure for report cache fingerprints.

Walk all imports, including function-local branches. This deliberately tracks
whole modules rather than guessing which function bodies can affect policy.
Unrelated acquisition/install scripts are not report dependencies.
"""
from __future__ import annotations

import ast
from pathlib import Path


def local_dependency_paths(directory, entries=('report_assessment.py', 'report_runtime.py')):
    directory = Path(directory).resolve()
    pending = [directory / name for name in entries]
    seen = set()
    while pending:
        path = pending.pop().resolve()
        if path in seen:
            continue
        if path.parent != directory or not path.is_file():
            raise ValueError('Report dependency is missing or outside the scripts directory: ' + str(path))
        seen.add(path)
        tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
        names = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names.update(alias.name.split('.')[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                names.add(node.module.split('.')[0])
            elif isinstance(node, ast.Call):
                dynamic = ((isinstance(node.func, ast.Name) and node.func.id == '__import__')
                           or (isinstance(node.func, ast.Attribute) and node.func.attr == 'import_module'))
                if dynamic:
                    if not node.args or not isinstance(node.args[0], ast.Constant) or not isinstance(node.args[0].value, str):
                        raise ValueError('Dynamic report imports require an explicit literal dependency: ' + str(path))
                    names.add(node.args[0].value.split('.')[0])
        for name in names:
            candidate = directory / (name + '.py')
            if candidate.is_file():
                pending.append(candidate)
    return sorted(seen)
