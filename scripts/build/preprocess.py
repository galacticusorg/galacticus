#!/usr/bin/env python3
"""Preprocess a Galacticus Fortran source file.

Parses the input source into a Galacticus::Build::SourceTree-compatible AST,
runs every registered process hook, optionally analyzes the tree, and writes
the resulting source — together with a `.lmap` line-number mapping file — to
the requested output path.

Setting `GALACTICUS_PREPROCESSOR_ANALYZE=yes` in the environment runs the
tree analyzer after processing. This is a manual debugging hook only — no
Makefile rule or CI workflow sets it.

Mirrors scripts/build/preprocess.pl.
Andrew Benson (ported to Python 2026).
"""

import os
import sys


from Galacticus.Build.FileChanges                       import update as file_changes_update
from Galacticus.Build.SourceTree              import parse_file, serialize, analyze_tree
from Galacticus.Build.SourceTree.Process      import process_tree
from Galacticus._logging                      import configure_default as _configure_default

# Show INFO-level diagnostic output from the library modules (mirrors the
# verbose `print()`-driven output of the Perl-era driver).
_configure_default()

# Register every source-tree process hook (the full set — see Process/all.py).
import Galacticus.Build.SourceTree.Process.all  # noqa: F401


def main(argv):
    if len(argv) != 3:
        print("Usage: preprocess.py <infile> <outfile>", file=sys.stderr)
        return 1

    input_file  = argv[1]
    output_file = argv[2]

    tree = parse_file(input_file)
    process_tree(tree)

    if os.environ.get('GALACTICUS_PREPROCESSOR_ANALYZE') == 'yes':
        analyze_tree(tree)

    source, mappings = serialize(tree, annotate=True, strip_mappings=True)

    tmp_file = output_file + '.tmp'
    with open(tmp_file, 'wb') as fh:
        fh.write(source.encode('utf-8'))
    file_changes_update(output_file, tmp_file, prove_update=True)

    with open(output_file + '.lmap', 'w', encoding='utf-8') as fh:
        fh.write(mappings)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
