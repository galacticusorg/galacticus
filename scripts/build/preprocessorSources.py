#!/usr/bin/env python3
"""Digest the Python sources that define the behavior of the Fortran preprocessor.

The preprocessor is a Python program: `scripts/build/preprocess.py` plus every
first-party module in its transitive import closure (the source-tree parser, the
`Process` hooks, the node-component generator, and so on). Editing any of those
changes the Fortran that `preprocess.py` emits, so every `*.p.F90` / `*.p.Inc`
must be regenerated — but Make cannot see that dependency, because the closure
is discovered by Python's import machinery, not declared anywhere.

This script makes the dependency explicit. It walks the import graph statically
(via `ast`, so nothing is imported and no side effects are triggered), resolving
module names against the repository's own `python/` tree and against the
directory holding the entry-point script; third-party and standard-library
imports resolve to nothing and are ignored. It then writes

    $BUILDPATH/preprocessorSources.digest

listing every file in the closure with an MD5 of its content, only-if-changed
(see `Galacticus.Build.FileChanges.update`), touching the `.up` sentinel on
every run. The Makefile hangs the preprocessing rules off that digest, so:

  * touching or editing *any* Python file re-runs this script (cheap), but
  * the digest — and hence the re-preprocessing sweep — changes only when a
    module that the preprocessor actually uses changes content.

Run with `--list` to print the closure (one path per line) instead of writing
the digest; useful when checking whether some module is, in fact, part of the
preprocessor.

Andrew Benson (2026).
"""

import ast
import hashlib
import os
import sys

from Galacticus.Build.FileChanges import update as file_changes_update

# Entry-point scripts whose import closure defines the preprocessor.
ENTRY_POINTS = ['scripts/build/preprocess.py']


def _module_file(name, roots):
    """Return the file implementing module `name` (dotted), or `None`.

    `roots` are directories searched in order, as `sys.path` would be. A package
    resolves to its `__init__.py`. Names that resolve to nothing — standard
    library, third-party packages, or the symbol half of a `from X import y` —
    return `None` and are simply skipped by the caller.
    """
    parts = name.split('.')
    for root in roots:
        base = os.path.join(root, *parts)
        init = os.path.join(base, '__init__.py')
        if os.path.isfile(init):
            return init
        if os.path.isfile(base + '.py'):
            return base + '.py'
    return None


def _package_of(path, roots):
    """Return the dotted package name containing the module at `path`.

    Empty for a file that is not inside any search root (e.g. a script under
    `scripts/build/`), which is correct: such a file cannot use relative
    imports.
    """
    path = os.path.abspath(path)
    for root in roots:
        root = os.path.abspath(root)
        if path.startswith(root + os.sep):
            relative = os.path.relpath(os.path.dirname(path), root)
            if relative == os.curdir:
                return ''
            return relative.replace(os.sep, '.')
    return ''


def _imported_names(path, roots):
    """Return the module names imported by the file at `path`.

    Every `import`/`from ... import` anywhere in the file is collected —
    including imports nested inside functions, which the preprocessor uses to
    break import cycles. For `import a.b.c` the intermediate packages `a` and
    `a.b` are reported too, since importing the leaf executes their
    `__init__.py`. For `from a.b import c` both `a.b` and the candidate
    submodule `a.b.c` are reported; if `c` is a symbol rather than a submodule
    it simply fails to resolve.
    """
    try:
        with open(path, 'rb') as fh:
            tree = ast.parse(fh.read(), filename=path)
    except (OSError, SyntaxError):
        return []

    package = _package_of(path, roots)
    names   = []

    def _add(name):
        parts = name.split('.')
        for i in range(1, len(parts) + 1):
            names.append('.'.join(parts[:i]))

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                _add(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                # Relative import: strip `level - 1` trailing components from
                # the containing package to get the base.
                ancestors = package.split('.') if package else []
                if node.level - 1 > len(ancestors):
                    continue
                base = ancestors[:len(ancestors) - (node.level - 1)]
                base = '.'.join(base + ([node.module] if node.module else []))
            else:
                base = node.module or ''
            if not base:
                continue
            _add(base)
            for alias in node.names:
                if alias.name != '*':
                    _add(base + '.' + alias.name)

    return names


def closure(entry_points, roots):
    """Return the sorted transitive first-party import closure of `entry_points`."""
    found   = set()
    pending = []
    for entry_point in entry_points:
        if os.path.isfile(entry_point):
            pending.append(os.path.abspath(entry_point))
    # Search the directory of each entry point too, so a script that imports a
    # sibling script under `scripts/build/` is followed.
    roots = list(roots) + sorted({os.path.dirname(os.path.abspath(e)) for e in pending})
    while pending:
        path = pending.pop()
        if path in found:
            continue
        found.add(path)
        for name in _imported_names(path, roots):
            resolved = _module_file(name, roots)
            if resolved is not None:
                resolved = os.path.abspath(resolved)
                if resolved not in found:
                    pending.append(resolved)
    return sorted(found)


def digest(paths, root):
    """Return the digest text: one `<md5>  <path>` line per file, path-sorted."""
    lines = []
    for path in paths:
        with open(path, 'rb') as fh:
            lines.append(hashlib.md5(fh.read()).hexdigest() + "  " + os.path.relpath(path, root) + "\n")
    return "".join(lines)


def main(argv):
    if len(argv) < 2 or len(argv) > 3 or (len(argv) == 3 and argv[2] != '--list'):
        print("Usage: preprocessorSources.py <galacticusPath> [--list]", file=sys.stderr)
        return 1

    galacticus_path = argv[1]
    entry_points    = [os.path.join(galacticus_path, entry_point) for entry_point in ENTRY_POINTS]
    paths           = closure(entry_points, [os.path.join(galacticus_path, 'python')])

    if not paths:
        print("preprocessorSources.py: found no preprocessor sources — is "
              + ", ".join(ENTRY_POINTS) + " present?", file=sys.stderr)
        return 1

    if len(argv) == 3:
        for path in paths:
            print(os.path.relpath(path, galacticus_path))
        return 0

    build_path  = os.environ['BUILDPATH']
    output_file = os.path.join(build_path, 'preprocessorSources.digest')
    tmp_file    = output_file + '.tmp'
    with open(tmp_file, 'w', encoding='utf-8') as fh:
        fh.write(digest(paths, galacticus_path))
    file_changes_update(output_file, tmp_file, prove_update=True)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
