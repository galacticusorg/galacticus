# If the `-DDEBUGMPI` flag is in $GALACTICUS_FCFLAGS, prefixes every
# `mpiSelf%...` call and `call mpiBarrier()` with a debug write.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/DebugMPI.pm

import io
import re
import os


from Galacticus.Build.FortranUtils                              import get_fortran_line
from Galacticus.Build.SourceTree                      import walk_tree
from Galacticus.Build.SourceTree.Process              import register_process
from Galacticus.Build.SourceTree.Process.SourceIntrospection import location


def _debug_enabled():
    flags = os.environ.get('GALACTICUS_FCFLAGS', '')
    return '-DDEBUGMPI' in flags.split()


def process_debug_mpi(tree, options):
    """Mirrors Process_DebugMPI() from DebugMPI.pm."""
    if not _debug_enabled():
        return

    for node in walk_tree(tree):
        if node.get('type') != 'code':
            continue
        new_content = []
        line_number = node.get('line', 0)
        fh = io.StringIO(node.get('content', ''))
        while True:
            raw_line, processed_line, _ = get_fortran_line(fh)
            if not raw_line and not processed_line:
                break
            line_number += 1

            m = re.search(r'mpiSelf%([a-zA-Z0-9_]+)', processed_line)
            if m:
                method = m.group(1)
                raw_line = (
                    "write (0,*) 'MPI DEBUG ['//char(mpiSelf%rankLabel())//"
                    f"': mpiSelf call to method \"{method}\"'//"
                    + location(node, line_number)
                    + "\n" + raw_line
                )
            if re.match(r'^\s*call\s+mpiBarrier\s*\(\s*\)', processed_line):
                raw_line = (
                    "write (0,*) 'MPI DEBUG ['//char(mpiSelf%rankLabel())//"
                    "']: mpiBarrier'//"
                    + location(node, line_number)
                    + "\n" + raw_line
                )
            new_content.append(raw_line)

        node['content'] = ''.join(new_content)


register_process('debugMPI', process_debug_mpi)
