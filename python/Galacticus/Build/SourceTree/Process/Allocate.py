# Processes `allocate` directives: emits an `allocate(var(…))` statement
# whose rank is inferred from the variable's declaration and whose shape is
# drawn from either a `shape=` or `size=` argument on the directive.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/Allocate.pm

import os
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from Galacticus.Build.SourceTree                    import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process            import register_process
from Galacticus.Build.SourceTree.Parse.Declarations import get_declaration


def _declaration_rank(parent, variable):
    """Return the rank (number of dimensions) of `variable` declared under
    `parent`, or 0 if the variable has no `dimension(...)` attribute.
    """
    declaration = get_declaration(parent, variable)
    for attr in declaration.get('attributes') or []:
        if attr.startswith('dimension'):
            # Count commas inside dimension(...).  Perl: ($attribute =~ tr/,//)+1
            return attr.count(',') + 1
    return 0


def process_allocate(tree, options):
    """Mirrors Process_Allocate() from Allocate.pm."""
    for node in walk_tree(tree):
        if node.get('type') != 'allocate':
            continue
        directive = node.get('directive') or {}
        if directive.get('processed'):
            continue
        directive['processed'] = True

        variable = directive['variable']
        if 'rank' in directive:
            rank = int(directive['rank'])
        else:
            rank = _declaration_rank(node['parent'], variable)

        if 'shape' in directive:
            source = directive['shape']
            indices = "" if rank == 0 else "(" + ",".join(
                f"{source}({i})" for i in range(1, rank + 1)) + ")"
        elif 'size' in directive:
            source = directive['size']
            indices = "" if rank == 0 else "(" + ",".join(
                f"size({source},dim={i})" for i in range(1, rank + 1)) + ")"
        else:
            raise RuntimeError("process_allocate: no source given for allocation")

        code  = "! Auto-generated allocation\n"
        code += f"allocate({variable}{indices})\n"
        code += "! End auto-generated allocation\n"

        insert_after_node(node, [{
            'type':       'code',
            'content':    code,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':     node.get('source', 'unknown'),
            'line':       node.get('line', 0),
        }])


register_process('allocate', process_allocate, before=['generics'])
