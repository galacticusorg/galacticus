# Processes `addMetaProperty` directives: emits a call to the
# `add<Type>Rank<N>MetaProperty` method of the enclosing component
# object, with a `<component>:<name>` prefix and the isCreator /
# isEvolvable flags set per the directive.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/AddMetaProperty.pm



from Galacticus.Build.SourceTree                    import (
    walk_tree, insert_after_node,
)
from Galacticus.Build.SourceTree.Process            import register_process
from Galacticus.Build.SourceTree.Parse.ModuleUses   import add_uses


_BOOLEAN = {'no': '.false.', 'yes': '.true.'}


def process_add_meta_property(tree, options):
    """Mirrors Process_AddMetaProperty() from AddMetaProperty.pm."""
    for node in walk_tree(tree):
        if node.get('type') != 'addMetaProperty':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        directive['processed'] = True

        directive.setdefault('type',        'float')
        directive.setdefault('rank',        0)
        directive.setdefault('isEvolvable', 'no')
        directive.setdefault('isCreator',   'no')

        # Perl stores directive values as strings; normalise `rank` for the
        # numeric comparisons below.
        try:
            rank = int(directive['rank'])
        except (TypeError, ValueError):
            rank = directive['rank']

        if directive['isEvolvable'] == 'yes' and directive['type'] != 'float':
            raise RuntimeError(
                "non-float meta-properties can not be evolvable")
        if directive['isEvolvable'] == 'yes' and rank > 0:
            raise RuntimeError(
                "rank > 0 meta-properties can not be evolvable")
        if rank > 1:
            raise RuntimeError(
                "rank > 1 meta-properties are not supported")

        component = (
            'default' + directive['component'][:1].upper()
            + directive['component'][1:] + 'Component')
        type_prefix = (
            directive['type'][:1].upper() + directive['type'][1:])
        type_suffix = f"Rank{rank}"

        name_attr = directive['name']
        if name_attr.startswith("'"):
            # Fortran character literal — pass through var_str(...).
            name         = f"var_str({name_attr})"
            prefixed_name = (
                f"'{directive['component']}:'//{name_attr}")
        else:
            name         = f"var_str('{name_attr}')"
            prefixed_name = (
                f"'{directive['component']}:{name_attr}'")

        evolvable_arg = ''
        if directive['type'] == 'float' and rank == 0:
            evolvable_arg = (
                f",isEvolvable={_BOOLEAN[directive['isEvolvable']]}")

        code = (
            f"{directive['id']}={component}%add{type_prefix}{type_suffix}"
            f"MetaProperty({name},{prefixed_name},isCreator="
            f"{_BOOLEAN[directive['isCreator']]}{evolvable_arg})\n"
        )

        add_uses(node['parent'], {
            'moduleUse': {
                'ISO_Varying_String': {
                    'intrinsic': False,
                    'only':      {'var_str': True},
                },
                'Galacticus_Nodes':   {
                    'intrinsic': False,
                    'only':      {component: True},
                },
            },
            'moduleOrder': ['ISO_Varying_String', 'Galacticus_Nodes'],
        })

        insert_after_node(node, [{
            'type':       'code',
            'content':    code,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':
                'Galacticus.Build.SourceTree.Process.AddMetaProperty'
                '.process_add_meta_property()',
            'line':       1,
        }])


register_process('addMetaProperty', process_add_meta_property)
