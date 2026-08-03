"""Per-class `<class>Type` name accessor.

Andrew Benson (ported to Python 2026)
"""
from __future__ import annotations

from Galacticus.Build.Components.Utils import register


def Class_Type(build: dict, class_dict: dict) -> None:
    """Generate `nodeComponent<Class>Type` returning the class label.

    Body just returns `'nodeComponent:<class>'`.
    """
    name      = class_dict['name']
    type_name = 'nodeComponent' + _ucfirst(name)
    function  = {
        'type':        'type(varying_string) => name',
        'name':        type_name + 'Type',
        'description': f"Returns the type name for the {name} component class.",
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
            f"name='nodeComponent:{name}'\n"
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


register('classNames', 'classIteratedFunctions', Class_Type)
