"""Per-property `<prop>Set` setters bound at the component level + the
class-level `<prop>IsSettable` boolean stub.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Components/Properties/Set.pm.
"""
from __future__ import annotations

from Galacticus.Build.Components.Utils      import (
    register,
    _component_properties,
)
from Galacticus.Build.Components.DataTypes  import data_object_definition


def _is_deferred(prop: dict, verb: str) -> bool:
    deferred = (prop.get('attributes') or {}).get('isDeferred')
    if not deferred:
        return False
    return verb in str(deferred).split(':')


def Build_Class_Setters(build: dict) -> None:
    """Insert one `<prop>IsSettable` Boolean_False stub per property
    into the matching `nodeComponent<Class>` type.

    Mirrors `Build_Class_Setters`.
    """
    for component in (build.get('components') or {}).values():
        cap_class = _ucfirst(component['class'])
        type_name = 'nodeComponent' + cap_class
        bound = build.setdefault('types', {}).setdefault(type_name, {}) \
                                              .setdefault('boundFunctions', [])
        for prop in _component_properties(component):
            fn_name = prop['name'] + 'IsSettable'
            if any(b.get('name') == fn_name for b in bound):
                continue
            bound.append({
                'type':        'procedure',
                'pass':        'nopass',
                'name':        fn_name,
                'function':    'Boolean_False',
                'returnType':  r"\logicalzero",
                'arguments':   "",
                'description': (
                    f"Specify whether the \\mono{{{prop['name']}}} "
                    f"property of the \\mono{{{component['class']}}} "
                    "component is settable."
                ),
            })


def Bind_Set_Functions(build: dict, class_dict: dict, member: dict, prop: dict) -> None:
    """Bind a compile-time custom set function to the component
    implementation when the user supplied one and the property's `set`
    is not deferred.

    Mirrors `Bind_Set_Functions`.
    """
    attrs = prop.get('attributes') or {}
    set_function = prop.get('setFunction') or {}
    if (
        attrs.get('isSettable')
        and not set_function.get('build')
        and not _is_deferred(prop, 'set')
    ):
        impl_type = (
            'nodeComponent'
            + _ucfirst(class_dict['name'])
            + _ucfirst(member['name'])
        )
        build.setdefault('types', {}).setdefault(impl_type, {}) \
                                      .setdefault('boundFunctions', []) \
                                      .append({
            'type':     'procedure',
            'name':     prop['name'] + 'Set',
            'function': set_function.get('content'),
        })


def Build_Set_Functions(build: dict, class_dict: dict, member: dict, prop: dict) -> None:
    """Build the auto-generated set function for a non-deferred,
    non-virtual settable property.  Mirrors `Build_Set_Functions`.
    """
    attrs = prop.get('attributes') or {}
    set_function = prop.get('setFunction') or {}
    if (
        attrs.get('isVirtual')
        or not attrs.get('isSettable')
        or not set_function.get('build')
    ):
        return

    suffix = 'Value' if _is_deferred(prop, 'set') else ''
    data = prop.get('data') or {}
    rank = int(data.get('rank') or 0)

    type_descriptor, _label = data_object_definition(data)
    type_descriptor['variables']  = ['setValue']
    type_descriptor['attributes'] = ['intent(in   )']
    if rank > 0:
        dims = ','.join([':'] * rank)
        type_descriptor['attributes'].append(f"dimension({dims})")

    impl_type = (
        'nodeComponent'
        + _ucfirst(class_dict['name'])
        + _ucfirst(member['name'])
    )
    fn_name = (
        class_dict['name']
        + _ucfirst(member['name'])
        + _ucfirst(prop['name'])
        + 'Set' + suffix
    )

    function = {
        'type':        'void',
        'name':        fn_name,
        'description': (
            f"Set the \\mono{{{prop['name']}}} property of an "
            f"\\mono{{{member['name']}}} implementation of the "
            f"\\mono{{{class_dict['name']}}} component class."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            type_descriptor,
        ],
    }

    if rank == 0:
        function['content'] = f"self%{prop['name']}Data=setValue\n"
    elif rank == 1:
        function['content'] = (
            f"if (.not.allocated(self%{prop['name']}Data)) then\n"
            "      !![\n"
            f"      <allocate variable=\"self%{prop['name']}Data\" "
            f"size=\"setValue\" rank=\"{rank}\"/>\n"
            "      !!]\n"
            "else\n"
            f"   if (size(self%{prop['name']}Data) /= size(setValue)) then\n"
            f"      deallocate(self%{prop['name']}Data)\n"
            "      !![\n"
            f"      <allocate variable=\"self%{prop['name']}Data\" "
            f"size=\"setValue\" rank=\"{rank}\"/>\n"
            "      !!]\n"
            "   end if\n"
            "end if\n"
            f"self%{prop['name']}Data=setValue\n"
        )

    build.setdefault('types', {}).setdefault(impl_type, {}) \
                                  .setdefault('boundFunctions', []) \
                                  .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       prop['name'] + 'Set' + suffix,
    })


def _ucfirst(text: str) -> str:
    return text[:1].upper() + text[1:] if text else text


register('propertiesSet', 'functions', Build_Class_Setters)
register('propertiesSet', 'propertyIteratedFunctions', Bind_Set_Functions)
register('propertiesSet', 'propertyIteratedFunctions', Build_Set_Functions)
