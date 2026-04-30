# Per-property `<prop>` getters bound at the component level.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Properties/Get.pm.

import os
import sys

sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))

from Galacticus.Build.Components.Utils      import register
from Galacticus.Build.Components.DataTypes  import data_object_definition


def _is_deferred(prop, verb):
    """Return True if `verb` ('get' / 'set') appears in the property's
    `attributes.isDeferred` (a colon-separated list).
    """
    deferred = (prop.get('attributes') or {}).get('isDeferred')
    if not deferred:
        return False
    return verb in str(deferred).split(':')


def Bind_Get_Functions(build, class_dict, member, prop):
    """Bind a compile-time custom get function to the component
    implementation when the user supplied one (`getFunction.build` is
    False) and the property's `get` is not deferred.

    Mirrors `Bind_Get_Functions`.
    """
    attrs = prop.get('attributes') or {}
    get_function = prop.get('getFunction') or {}
    if (
        attrs.get('isGettable')
        and not get_function.get('build')
        and not _is_deferred(prop, 'get')
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
            'name':     prop['name'],
            'function': get_function.get('content'),
        })


def Build_Get_Functions(build, class_dict, member, prop):
    """Build the auto-generated get function for a non-deferred,
    non-virtual gettable property.  Mirrors `Build_Get_Functions`.
    """
    attrs = prop.get('attributes') or {}
    get_function = prop.get('getFunction') or {}
    if (
        attrs.get('isVirtual')
        or not attrs.get('isGettable')
        or not get_function.get('build')
    ):
        return

    suffix = 'Value' if _is_deferred(prop, 'get') else ''

    type_descriptor, _label = data_object_definition(prop.get('data') or {})
    function_type = type_descriptor['intrinsic']
    if 'type' in type_descriptor:
        function_type += f"({type_descriptor['type']})"
    if type_descriptor.get('attributes'):
        function_type += ', ' + ', '.join(type_descriptor['attributes'])

    impl_type = (
        'nodeComponent'
        + _ucfirst(class_dict['name'])
        + _ucfirst(member['name'])
    )
    fn_name = (
        class_dict['name']
        + _ucfirst(member['name'])
        + _ucfirst(prop['name'])
        + 'Get' + suffix
    )

    function = {
        'type':        function_type + ' => propertyValue',
        'name':        fn_name,
        'description': (
            f"Get the \\mono{{{prop['name']}}} property of an "
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
        ],
        'content': f"propertyValue=self%{prop['name']}Data\n",
    }
    build.setdefault('types', {}).setdefault(impl_type, {}) \
                                  .setdefault('boundFunctions', []) \
                                  .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       prop['name'] + suffix,
    })


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


register('propertiesGet', 'propertyIteratedFunctions', Bind_Get_Functions)
register('propertiesGet', 'propertyIteratedFunctions', Build_Get_Functions)
