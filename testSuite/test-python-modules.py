#!/usr/bin/env python3
"""Compile every Python source file in the repository to catch syntax errors.

Replaces the former `test-perl-modules.py`, which ran `perl -c` over the
`.pl` and `perl/**/*.pm` trees.  Now that the build pipeline is pure
Python, this is the equivalent sweep for `*.py` files.
"""

import py_compile
import pathlib
import sys

repo_root = pathlib.Path(__file__).resolve().parent.parent

status = 0
for path in sorted(repo_root.rglob('*.py')):
    if any(part.startswith('.') for part in path.relative_to(repo_root).parts):
        continue
    try:
        py_compile.compile(str(path), doraise=True)
    except py_compile.PyCompileError as exc:
        print(f"FAILED: {path}")
        print(exc)
        status = 1

sys.exit(status)
