# Generate `<class><Property>AttributeMatch` accessors per property
# of a component class — return the list of implementations matching a
# requested attribute requirement.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Properties/Attributes.pm.



from Galacticus.Build.Components.Utils                import (
    register,
    _component_properties,
)
from Galacticus.Build.Components.Properties.Utils     import attribute_adjective


def Attributes_Match(build, class_dict):
    """Mirrors `Attributes_Match`."""
    cap_class = _ucfirst(class_dict['name'])
    type_name = 'nodeComponent' + cap_class

    # Build a {prop_name → {member_name → {get/set/rate flags}}} map.
    properties = {}
    for member in class_dict.get('members') or []:
        # Walk member + ancestors so an inherited property is recorded
        # against this member too.
        cursor = member
        while isinstance(cursor, dict):
            for prop in _component_properties(cursor):
                slot = (
                    properties
                    .setdefault(prop['name'], {'name': prop['name'], 'members': {}})
                    ['members']
                    .setdefault(member['name'], {})
                )
                attrs = prop.get('attributes') or {}
                for verb in ('set', 'get', 'rate'):
                    slot[verb] = bool(attrs.get(attribute_adjective[verb]))
            ext = cursor.get('extends') or {}
            cursor = ext.get('implementation') if isinstance(ext, dict) else None

    for prop_name in sorted(properties.keys()):
        prop_entry = properties[prop_name]
        cap_prop = _ucfirst(prop_name)

        function = {
            'type':        'type(varying_string), allocatable, dimension(:) => matches',
            'name':        f"{class_dict['name']}{cap_prop}AttributeMatch",
            'description': (
                f"Return a text list of component implementations in the "
                f"\\mono{{{class_dict['name']}}} class that have the desired "
                f"attributes for the \\mono{{{prop_name}}} property"
            ),
            'modules':     ['ISO_Varying_String'],
            'variables':   [
                {
                    'intrinsic':  'logical',
                    'attributes': ['intent(in   )', 'optional'],
                    'variables':  [
                        'requireSettable', 'requireGettable', 'requireEvolvable',
                    ],
                },
                {
                    'intrinsic':  'logical',
                    'variables':  [
                        'requireSettableActual', 'requireGettableActual',
                        'requireEvolvableActual',
                    ],
                },
                {
                    'intrinsic':  'type',
                    'type':       'varying_string',
                    'attributes': ['allocatable', 'dimension(:)'],
                    'variables':  ['temporaryList'],
                },
            ],
        }
        content = (
            "requireSettableActual =.false.\n"
            "requireGettableActual =.false.\n"
            "requireEvolvableActual=.false.\n"
            "if (present(requireSettable )) requireSettableActual =requireSettable\n"
            "if (present(requireGettable )) requireGettableActual =requireGettable\n"
            "if (present(requireEvolvable)) requireEvolvableActual=requireEvolvable\n"
        )

        # Iterate over members in sorted order (Perl `hashList(..., keyAs)`
        # iterates the underlying hash without sorting; we sort for
        # determinism).
        for member_name in sorted(prop_entry['members'].keys()):
            flags = prop_entry['members'][member_name]
            logic = []
            for verb in ('get', 'set', 'rate'):
                if not flags.get(verb):
                    suffix = (
                        'Evolvable' if verb == 'rate'
                        else _ucfirst(verb) + 'table'
                    )
                    logic.append(f".not.require{suffix}Actual")
            if logic:
                content += f"if ({'.and.'.join(logic)}) then\n"
            content += (
                "if (allocated(matches)) then\n"
                "   call Move_Alloc(matches,temporaryList)\n"
                "   allocate(matches(size(temporaryList)+1))\n"
                "   matches(1:size(temporaryList))=temporaryList\n"
                "   deallocate(temporaryList)\n"
                "else\n"
                "   allocate(matches(1))\n"
                "end if\n"
                f"matches(size(matches))='{member_name}'\n"
            )
            if logic:
                content += "end if\n"

        function['content'] = content

        build.setdefault('types', {}).setdefault(type_name, {}) \
                                      .setdefault('boundFunctions', []) \
                                      .append({
            'type':       'procedure',
            'descriptor': function,
            'pass':       'nopass',
            'name':       prop_name + 'AttributeMatch',
        })


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


register('propertiesAttributes', 'classIteratedFunctions', Attributes_Match)
