#!/usr/bin/env python3
"""Preprocess a Galacticus Fortran source file.

Parses the input source into a Galacticus::Build::SourceTree-compatible AST,
runs every registered process hook, optionally analyzes the tree, and writes
the resulting source — together with a `.lmap` line-number mapping file — to
the requested output path.

Mirrors scripts/build/preprocess.pl.
Andrew Benson (ported to Python 2026).
"""

import os
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from build.file_changes                       import update as file_changes_update
from Galacticus.Build.SourceTree              import parse_file, serialize, analyze_tree
from Galacticus.Build.SourceTree.Process      import process_tree

# Importing each Process submodule registers its hook.  Mirrors the `use`
# statements at the top of perl/Galacticus/Build/SourceTree.pm.
import Galacticus.Build.SourceTree.Process.AddMetaProperty          # noqa: F401
import Galacticus.Build.SourceTree.Process.Allocate                 # noqa: F401
import Galacticus.Build.SourceTree.Process.ClassDocumentation       # noqa: F401
import Galacticus.Build.SourceTree.Process.ConditionalCall          # noqa: F401
import Galacticus.Build.SourceTree.Process.Constants                # noqa: F401
import Galacticus.Build.SourceTree.Process.Constructors             # noqa: F401
import Galacticus.Build.SourceTree.Process.DebugHDF5                # noqa: F401
import Galacticus.Build.SourceTree.Process.DebugMPI                 # noqa: F401
import Galacticus.Build.SourceTree.Process.DeepCopyActions          # noqa: F401
import Galacticus.Build.SourceTree.Process.DeepCopyFinalize         # noqa: F401
import Galacticus.Build.SourceTree.Process.DeepCopyReset            # noqa: F401
import Galacticus.Build.SourceTree.Process.Dependencies             # noqa: F401
import Galacticus.Build.SourceTree.Process.Enumeration              # noqa: F401
import Galacticus.Build.SourceTree.Process.EventHooks               # noqa: F401
import Galacticus.Build.SourceTree.Process.EventHooksStatic         # noqa: F401
import Galacticus.Build.SourceTree.Process.ForEach                  # noqa: F401
import Galacticus.Build.SourceTree.Process.FunctionClass            # noqa: F401
import Galacticus.Build.SourceTree.Process.FunctionsGlobal          # noqa: F401
import Galacticus.Build.SourceTree.Process.Generics                 # noqa: F401
import Galacticus.Build.SourceTree.Process.HDF5FCInterop            # noqa: F401
import Galacticus.Build.SourceTree.Process.InputParameter           # noqa: F401
import Galacticus.Build.SourceTree.Process.InputParametersValidate  # noqa: F401
import Galacticus.Build.SourceTree.Process.MetaPropertyDatabase     # noqa: F401
import Galacticus.Build.SourceTree.Process.NonProcessed             # noqa: F401
import Galacticus.Build.SourceTree.Process.ObjectBuilder            # noqa: F401
import Galacticus.Build.SourceTree.Process.OptionalArgument         # noqa: F401
import Galacticus.Build.SourceTree.Process.ParameterMigration       # noqa: F401
import Galacticus.Build.SourceTree.Process.ProfileOpenMP            # noqa: F401
import Galacticus.Build.SourceTree.Process.SourceDigest             # noqa: F401
import Galacticus.Build.SourceTree.Process.SourceIntrospection      # noqa: F401
import Galacticus.Build.SourceTree.Process.StateStorable            # noqa: F401
import Galacticus.Build.SourceTree.Process.StateStore               # noqa: F401
import Galacticus.Build.SourceTree.Process.ThreadSafeIO             # noqa: F401


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
