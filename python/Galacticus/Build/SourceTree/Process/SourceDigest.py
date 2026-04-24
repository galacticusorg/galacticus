# Processes `sourceDigest` directives: emits a C-interop character array
# declaration bound to the per-build-target MD5 hash symbol (defined by
# the Perl build system's `Find_Hash` chain), plus the matching
# `use :: ISO_C_Binding, only : C_Char` import.
# Andrew Benson (ported to Python 2026)
#
# Mirrors the `Process_SourceDigests` hook from
# perl/Galacticus/Build/SourceTree/Process/SourceDigest.pm.  The Perl
# module also ships `Find_Hash` / `Hash_Data_Files` / `modificationTime`
# helpers used by the Makefile to actually compute the MD5 hashes; those
# are not part of the Process pipeline and are deferred to a follow-up
# port.

import os
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from Galacticus.Build.SourceTree                    import walk_tree
from Galacticus.Build.SourceTree.Process            import register_process
from Galacticus.Build.SourceTree.Parse.Declarations import add_declarations
from Galacticus.Build.SourceTree.Parse.ModuleUses   import add_uses


def process_source_digests(tree, options):
    """Mirrors Process_SourceDigests() from SourceDigest.pm:24-58."""
    for node in walk_tree(tree):
        if node.get('type') != 'sourceDigest':
            continue
        directive = node.get('directive') or {}
        if directive.get('processed'):
            continue
        directive['processed'] = True

        name = directive['name']
        add_declarations(node['parent'], [{
            'intrinsic':     'character',
            'type':          'C_Char',
            'openMP':        False,
            'attributes':    [
                'dimension(23)',
                f'bind(C, name="{name}MD5")',
            ],
            'variables':     [name],
            'variableNames': [name],
        }])
        add_uses(node['parent'], {
            'moduleUse': {
                'ISO_C_Binding': {
                    'intrinsic': True,
                    'only':      {'C_Char': True},
                },
            },
            'moduleOrder': ['ISO_C_Binding'],
        })


register_process('sourceDigests', process_source_digests)
