# Per-class `<class>Type` name accessor.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Classes/Names.pm.



from Galacticus.Build.Components.Utils import register


def Class_Type(build, class_dict):
    """Generate `nodeComponent<Class>Type` returning the class label.

    Mirrors `Class_Type`.  Body just returns `'nodeComponent:<class>'`.
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


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


register('classNames', 'classIteratedFunctions', Class_Type)
