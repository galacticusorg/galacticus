# Per-class generic ASCII-serialization stub.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Classes/Serialization.pm.



from Galacticus.Build.Components       import Utils as _Utils
from Galacticus.Build.Components.Utils import register


def Class_Serialize_ASCII(build, class_dict):
    """Generate `nodeComponent<Class>SerializeASCII` — the abstract
    parent-class stub that just emits a single `'<class>: generic'`
    line.  Concrete implementations override this in the
    Implementations stage.

    Mirrors `Class_Serialize_ASCII`.
    """
    name      = class_dict['name']
    cap       = _ucfirst(name)
    type_name = 'nodeComponent' + cap

    padding_len = max(
        (_Utils.fully_qualified_name_length_max or 0) - len(name), 0
    )
    padding = ' ' * padding_len

    function = {
        'type':        'void',
        'name':        type_name + 'SerializeASCII',
        'description': (
            f"Serialize the content of a \\mono{{{name}}} component to ASCII."
        ),
        'modules':     ['Display', 'ISO_Varying_String'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       type_name,
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'type',
                'type':       'enumerationVerbosityLevelType',
                'variables':  ['verbosityLevel'],
                'attributes': ['intent(in   )'],
            },
        ],
        'content': (
            "!$GLC attributes unused :: self\n"
            f"call displayIndent('{name}: {padding}generic',verbosityLevel)\n"
            "call displayUnindent('done',verbosityLevel)\n"
        ),
    }
    build.setdefault('types', {}).setdefault(type_name, {}) \
                                  .setdefault('boundFunctions', []) \
                                  .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       'serializeASCII',
    })


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


register('classesSerialization', 'classIteratedFunctions', Class_Serialize_ASCII)
