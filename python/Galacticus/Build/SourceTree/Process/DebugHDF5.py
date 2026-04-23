# If the `-DDEBUGHDF5` flag is in $GALACTICUS_FCFLAGS, replaces calls to
# `hdf5Access%set()` / `hdf5Access%unset()` with explicit `IO_HDF5_Start_Locked`
# / `IO_HDF5_End_Locked` calls, and injects the corresponding `use IO_HDF5`
# in the enclosing subprogram.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/DebugHDF5.pm

import re
import os
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from Galacticus.Build.SourceTree                  import walk_tree
from Galacticus.Build.SourceTree.Process          import register_process
from Galacticus.Build.SourceTree.Parse.ModuleUses import add_uses


def _debug_enabled():
    return '-DDEBUGHDF5' in os.environ.get('GALACTICUS_FCFLAGS', '').split()


def _enclosing_subprogram(node):
    """Return the enclosing function/subroutine/program ancestor of `node`.

    Mirrors the loop at DebugHDF5.pm:56-58.
    """
    cur = node
    while cur is not None and cur.get('type') not in ('function', 'subroutine', 'program'):
        cur = cur.get('parent')
    return cur


def _enclosing_module_name(node):
    """Return the enclosing `module` node's name (or '' if not in a module)."""
    cur = node
    while cur is not None:
        if cur.get('type') == 'module':
            return cur.get('name', '')
        cur = cur.get('parent')
    return ''


def _add_hdf5_use(subprogram, module_name):
    """Inject `use IO_HDF5, only : IO_HDF5_Start_Locked, IO_HDF5_End_Locked`
    into `subprogram` unless its enclosing module is one we must not add
    a use-dependency to (IO_HDF5 or Error).
    """
    if module_name in ('IO_HDF5', 'Error'):
        return
    add_uses(subprogram, {
        'moduleUse': {
            'IO_HDF5': {
                'openMP':    False,
                'intrinsic': False,
                'only': {
                    'IO_HDF5_Start_Locked': True,
                    'IO_HDF5_End_Locked':   True,
                },
            },
        },
        'moduleOrder': ['IO_HDF5'],
    })


def process_debug_hdf5(tree, options):
    """Mirrors Process_DebugHDF5() from DebugHDF5.pm."""
    if not _debug_enabled():
        return

    set_re   = re.compile(r'call\s+hdf5Access\s*%\s*set\s*\(\s*\)')
    unset_re = re.compile(r'call\s+hdf5Access\s*%\s*unset\s*\(\s*\)')

    for node in walk_tree(tree):
        if node.get('type') != 'code':
            continue
        content = node.get('content', '')
        if 'hdf5Access' not in content:
            continue

        subprogram = _enclosing_subprogram(node)
        module_name = _enclosing_module_name(node)
        # Perl emits the call substitution only when the enclosing module is
        # neither IO_HDF5 nor Error (those define the locked routines); the
        # module-use statement is suppressed in those modules too.
        add_call = module_name not in ('IO_HDF5', 'Error')

        new_lines = []
        use_added = False
        for line in content.splitlines(keepends=True):
            if set_re.search(line) or unset_re.search(line):
                if not use_added and subprogram is not None:
                    _add_hdf5_use(subprogram, module_name)
                    use_added = True
                if add_call:
                    line = set_re.sub('call IO_HDF5_Start_Locked()', line)
                    line = unset_re.sub('call IO_HDF5_End_Locked()', line)
            new_lines.append(line)
        node['content'] = ''.join(new_lines)


register_process('debugHDF5', process_debug_hdf5)
