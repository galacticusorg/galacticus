#!/usr/bin/env python3
# Format module use sections of Fortran source files.
# Andrew Benson (2026)

import os
import sys
import tempfile
import argparse


from Galacticus.Build.SourceTree import (
    parse_code_for_format,
    serialize_for_format,
    preprocessor_guard_balance,
    walk_tree,
)
from Galacticus.Build.SourceTree.Parse.ModuleUses import update_uses

parser = argparse.ArgumentParser(description='Format module use sections of Fortran source files')
parser.add_argument('infile', help='Input Fortran source file')
parser.add_argument('--suffix', default='~', help='Suffix for backup file (default: ~)')
parser.add_argument('--check', action='store_true',
                    help='Report whether the file is correctly formatted without modifying it; '
                         'exit status 1 if it would change')
parser.add_argument('--no-backup', action='store_true',
                    help='Overwrite the input file without leaving a backup')
args = parser.parse_args()



def module_use_summary(tree):
    """Reduce a tree to the module-use content that formatting must preserve.

    Maps each module name to the set of `(condition set, symbol set)` pairs
    declared for it, ignoring order and layout — exactly the axes the
    formatter is allowed to change.  Comparing this before and after catches
    a dropped, duplicated, or mangled symbol even though the formatter sorts
    symbols and merges repeated `use` statements for the same module.
    """
    summary = {}
    for node in walk_tree(tree):
        if node.get('type') != 'moduleUse':
            continue
        for name, entries in node.get('moduleUse', {}).items():
            if not isinstance(entries, list):
                entries = [entries]
            for entry in entries:
                conditions = tuple(sorted(
                    (c['name'], c['invert']) for c in entry.get('conditions', [])))
                symbols = frozenset(entry.get('only', {}).keys())
                summary.setdefault(name.lower(), set()).add(
                    (conditions, symbols, bool(entry.get('all')),
                     bool(entry.get('openMP')), bool(entry.get('intrinsic'))))
    return summary


with open(args.infile, 'r', errors='replace') as fh:
    original = fh.read()

tree = parse_code_for_format(original, source=args.infile)

# Safety check 1 — the parse/serialize round-trip must be exact *before* any
# formatting is applied.  If it is not, the parser has failed to represent
# something in this file and any output would silently corrupt it, so refuse
# to write rather than guess.
unmodified = serialize_for_format(tree['firstChild'])
if unmodified != original:
    sys.stderr.write(
        f'{args.infile}: parse/serialize round-trip is not exact — refusing to '
        f'reformat.\n')
    sys.exit(2)

before = module_use_summary(tree)

for node in walk_tree(tree):
    if node.get('type') == 'moduleUse':
        update_uses(node)

content = serialize_for_format(tree['firstChild'])

# Safety check 2 — reformatting must not change which symbols are imported
# from which modules under which preprocessor conditions.  This re-parses the
# formatted text rather than trusting the tree it came from, so a bug in the
# emitter is caught too.
after = module_use_summary(parse_code_for_format(content, source=args.infile))
if preprocessor_guard_balance(content) != preprocessor_guard_balance(original):
    sys.stderr.write(
        f'{args.infile}: reformatting would leave the preprocessor conditionals '
        f'unbalanced — refusing to write.\n')
    sys.exit(2)

if after != before:
    sys.stderr.write(
        f'{args.infile}: reformatting would change the module uses — refusing to '
        f'write.\n')
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
