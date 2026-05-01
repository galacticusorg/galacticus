# Processes `deepCopyReset` directives: emits a call to `%deepCopyReset()` for
# each named object and marks the directive as processed.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/DeepCopyReset.pm



from Galacticus.Build.SourceTree         import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process import register_process


def process_deep_copy_reset(tree, options):
    """Mirrors Process_DeepCopyReset() from DeepCopyReset.pm."""
    for node in walk_tree(tree):
        if node.get('type') != 'deepCopyReset':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        objects = (directive.get('variables') or '').split()
        code = ''.join(f"call {obj}%deepCopyReset()\n" for obj in objects)
        insert_after_node(node, [{
            'type':       'code',
            'content':    code,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':     node.get('source', 'unknown'),
            'line':       node.get('line', 0),
        }])
        directive['processed'] = True


register_process('deepCopyReset', process_deep_copy_reset)
