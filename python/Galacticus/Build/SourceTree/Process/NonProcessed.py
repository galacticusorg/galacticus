# Marks a fixed list of directive types as processed so they do not trigger
# the `post_process_directives` check — these directives are consumed by
# other parts of the build (e.g. direct XML extractors), not by the
# SourceTree Process pipeline.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/NonProcessed.pm

import os
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from Galacticus.Build.SourceTree         import walk_tree
from Galacticus.Build.SourceTree.Process import register_process

# Directive types that we simply mark as processed.  Mirrors the
# @nonProcessedDirectives list at NonProcessed.pm:17.
_NON_PROCESSED_DIRECTIVES = frozenset((
    'methods',
    'workaround',
    'include',
    'functionGlobal',
    'component',
    'radiusSolverPlausibility',
    'interTreePositionInsert',
    'expiry',
    'scoping',
))


def is_non_processed_type(node_type):
    """Return True if `node_type` is a directive that does not require any
    Process hook to handle it (mirrors the NonProcessed.pm exemption list,
    plus any `*Task` type).

    Used both by `process_non_processed` (to mark such directives as
    processed early in the pipeline) and by `post_process_directives` (to
    forgive any later-injected directive of the same type — e.g. a
    `<methods>` block emitted by a code-generating Process hook into the
    tree after `nonProcessed` has already run).
    """
    return bool(node_type) and (
        node_type.endswith('Task') or node_type in _NON_PROCESSED_DIRECTIVES
    )


def process_non_processed(tree, options):
    """Mirrors Process_NonProcessed() from NonProcessed.pm."""
    for node in walk_tree(tree):
        if is_non_processed_type(node.get('type', '')):
            directive = node.get('directive')
            if directive is not None:
                directive['processed'] = True


register_process('nonProcessed', process_non_processed, before=['generics'])
