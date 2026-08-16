#!/usr/bin/env python3
# Format variable declaration blocks of Fortran source files.
# Andrew Benson (2026)

import os
import re
import sys
import tempfile
import argparse


from Galacticus.Build.SourceTree import (
    parse_code_for_format,
    serialize_for_format,
    walk_tree,
)
from Galacticus.Build.SourceTree.Parse.Declarations import format_declarations_tree

parser = argparse.ArgumentParser(
    description='Format variable declaration blocks of Fortran source files')
parser.add_argument('infile', help='Input Fortran source file')
parser.add_argument('--suffix', default='~', help='Suffix for backup file (default: ~)')
parser.add_argument('--check', action='store_true',
                    help='Report whether the file is correctly formatted without modifying it; '
                         'exit status 1 if it would change')
parser.add_argument('--no-backup', action='store_true',
                    help='Overwrite the input file without leaving a backup')
parser.add_argument('--include-external', action='store_true',
                    help='Also format files under source/external/ (vendored third-party '
                         'code, skipped by default)')
args = parser.parse_args()


def is_vendored(path):
    """True for third-party code vendored under `source/external/`.

    Reformatting vendored sources makes them impossible to diff against
    upstream, so they are left alone by default — the same exclusion the build
    applies when discovering executables (`scripts/build/findExecutables.py`).
    """
    parts = os.path.normpath(os.path.abspath(path)).split(os.sep)
    for index in range(len(parts) - 1):
        if parts[index] == 'source' and parts[index + 1] == 'external':
            return True
    return False


def declaration_summary(tree):
    """Reduce a tree to the declaration content that formatting must preserve.

    For each declaration: the intrinsic type, the type-spec, the *set* of
    attributes (the formatter reorders them, so order is deliberately dropped)
    and the variables with their initializers.  Whitespace and case are
    normalized, since Fortran is insensitive to both and the emitter adjusts
    both.  Comparing this before and after catches a dropped variable, a lost
    attribute or a mangled initializer.
    """
    summary = []
    for node in walk_tree(tree):
        if node.get('type') != 'declaration':
            continue
        for declaration in node.get('declarations', []):
            attributes = frozenset(
                re.sub(r'\s+', '', a).lower()
                for a in (declaration.get('attributes') or []))
            type_spec = declaration.get('type')
            summary.append((
                declaration['intrinsic'],
                re.sub(r'\s+', '', type_spec).lower() if type_spec else '',
                attributes,
                tuple(re.sub(r'\s+', '', v).lower()
                      for v in declaration.get('variables', [])),
                bool(declaration.get('openMP')),
            ))
    return sorted(summary)


if is_vendored(args.infile) and not args.include_external:
    sys.exit(0)

with open(args.infile, 'r', errors='replace') as fh:
    original = fh.read()

tree = parse_code_for_format(original, source=args.infile)

# Safety check 1 — the parse/serialize round-trip must be exact *before* any
# formatting is applied.  If it is not, the parser has failed to represent
# something in this file and any output would silently corrupt it.
unmodified = serialize_for_format(tree['firstChild'])
if unmodified != original:
    sys.stderr.write(
        f'{args.infile}: parse/serialize round-trip is not exact — refusing to '
        f'reformat.\n')
    sys.exit(2)

before = declaration_summary(tree)

format_declarations_tree(tree)
content = serialize_for_format(tree['firstChild'])

# Safety check 2 — reformatting must not change what is declared.  This
# re-parses the formatted text rather than trusting the tree it came from, so a
# bug in the emitter is caught too.
after = declaration_summary(parse_code_for_format(content, source=args.infile))
if after != before:
    sys.stderr.write(
        f'{args.infile}: reformatting would change the declarations — refusing '
        f'to write.\n')
    sys.exit(2)

if content == original:
    sys.exit(0)

if args.check:
    print(f'{args.infile}: would be reformatted')
    sys.exit(1)

dirpath  = os.path.dirname(os.path.abspath(args.infile))
tmp_path = None
try:
    fd, tmp_path = tempfile.mkstemp(dir=dirpath)
    with os.fdopen(fd, 'w') as fh:
        fh.write(content)
        fh.flush()
        os.fsync(fh.fileno())
    if args.no_backup:
        os.replace(tmp_path, args.infile)
    else:
        os.replace(args.infile, args.infile + args.suffix)
        os.replace(tmp_path, args.infile)
    tmp_path = None
except Exception:
    if tmp_path is not None:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
    raise
