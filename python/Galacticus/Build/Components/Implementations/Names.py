# Per-implementation type-name accessor.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Implementations/Names.pm.
from __future__ import annotations

from Galacticus.Build.Components.Utils import register


def Implementation_Type(build: dict, class_dict: dict, member: dict) -> None:
    """Generate `nodeComponent<Class><Member>Type` returning a static
    `'nodeComponent:<class>:<member>'` label.

    Mirrors `Implementation_Type`.
    """
    cap_class  = _ucfirst(class_dict['name'])
    cap_member = _ucfirst(member['name'])
    type_name  = 'nodeComponent' + cap_class + cap_member
    function = {
        'type':        'type(varying_string) => name',
        'name':        type_name + 'Type',
        'description': (
            f"Returns the type name for the {member['name']} implementation "
            f"of the {class_dict['name']} component class."
        ),
        'modules':     ['ISO_Varying_String'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       type_name,
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
        ],
        'content': (
            "!$GLC attributes unused :: self\n"
            f"name='nodeComponent:{class_dict['name']}:{member['name']}'\n"
        ),
    }
    build.setdefault('types', {}).setdefault(type_name, {}) \
                                  .setdefault('boundFunctions', []) \
                                  .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       'type',
    })


def _ucfirst(text: str) -> str:
    return text[:1].upper() + text[1:] if text else text


register('implementationNames', 'implementationIteratedFunctions', Implementation_Type)
