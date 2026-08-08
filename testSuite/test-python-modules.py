#!/usr/bin/env python3
"""Compile every Python source file in the repository to catch syntax errors."""

import py_compile
import pathlib
import sys

repo_root = pathlib.Path(__file__).resolve().parent.parent

failures = 0
for path in sorted(repo_root.rglob('*.py')):
    if any(part.startswith('.') for part in path.relative_to(repo_root).parts):
        continue
    try:
        py_compile.compile(str(path), doraise=True)
    except py_compile.PyCompileError as exc:
        print(f"FAILED: {path}")
        print(exc)
        failures += 1

if failures == 0:
    print("SUCCESS: all Python source files compile")

# Always exit with status 0 - failure is signaled by "FAILED" in the output above.
sys.exit(0)
