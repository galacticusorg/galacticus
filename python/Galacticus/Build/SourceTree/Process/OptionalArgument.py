"""Processes `optionalArgument` directives: clones the named argument's
declaration with an underscore-suffixed name, strips `optional`/`intent(...)`
attributes, and emits a setter that defaults the `_`-variable and copies in
the caller-supplied value when `present()`.

Andrew Benson (ported to Python 2026)
"""

import re


from Galacticus.Build.SourceTree                      import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process              import register_process
from Galacticus.Build.SourceTree.Parse.Declarations   import get_declaration, add_declarations


def _is_argument_attribute(attr):
    """True if attr is `optional` or any `intent(...)` form (including `intent()`)."""
    return attr == 'optional' or re.match(r'^intent\((in)?(out)?\)', attr) is not None


def process_optional_arguments(tree, options):
    """Process `optionalArgument` directives in the tree."""
    for node in walk_tree(tree):
        if node.get('type') != 'optionalArgument':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        directive['processed'] = True

        name = directive['name']
        parent = node['parent']

        declaration = get_declaration(parent, name)
        declaration['variables']  = [name + '_']
        declaration['attributes'] = [
            a for a in (declaration.get('attributes') or [])
            if not _is_argument_attribute(a)
        ]
        if 'if' in directive:
            declaration['preprocessor'] = directive['if']

        add_declarations(parent, [declaration])

        setter  = "   ! Auto-generated optional argument setter\n"
        setter += f"   {name}_={directive['defaultsTo']}\n"
        setter += f"   if (present({name})) {name}_={name}\n"
        setter += "   ! End auto-generated optional argument setter\n"

        insert_after_node(node, [{
            'type':       'code',
            'content':    setter,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':     node.get('source', 'unknown'),
            'line':       node.get('line', 0),
        }])


register_process('optionalArguments', process_optional_arguments)
