# Per-property deferred-attribute hooks.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Properties/Deferred.pm.
# Four propertyIteratedFunctions hooks emit module-scope procedure
# pointers + IsAttached flags + wrapper functions + attacher methods
# for each `isDeferred` attribute (`get` / `set` / `rate`) on every
# component property.

import os
import sys

sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))

from Galacticus.Build.Components.Utils         import register
from Galacticus.Build.Components.DataTypes     import data_object_definition
from Galacticus.Build.Components.NullFunctions import create_null_function
from Galacticus.Build.Components.Properties.Utils import attribute_adjective


# Module-level pointer-creation cache.  Mirrors `%createdPointers` at
# Properties/Deferred.pm:33 — across all calls to
# `Properties_Deferred_Pointers` we emit each unique pointer once.
_created_pointers = set()


def _deferred_attributes(prop):
    """Return the list of deferred attribute names declared on `prop`."""
    raw = (prop.get('attributes') or {}).get('isDeferred')
    if not raw:
        return []
    return [a for a in str(raw).split(':') if a]


def Properties_Deferred_Pointers(build, class_dict, member, prop):
    """Mirrors `Properties_Deferred_Pointers`."""
    if not _deferred_attributes(prop):
        return
    attach_to = class_dict['name'] + _ucfirst(member['name'])
    for attribute in _deferred_attributes(prop):
        function_label = (
            _lcfirst(attach_to)
            + _ucfirst(prop['name'])
            + _ucfirst(attribute)
        )
        adjective = attribute_adjective.get(attribute)
        if (
            adjective is None
            or not (prop.get('attributes') or {}).get(adjective)
            or function_label in _created_pointers
        ):
            continue
        if attribute == 'get':
            template = (
                class_dict['name']
                + _ucfirst(member['name'])
                + _ucfirst(prop['name'])
                + _ucfirst(attribute)
            )
        else:
            template = create_null_function(build, {
                'selfType':  class_dict['name'],
                'attribute': attribute,
                'property':  prop,
                'intent':    'inout',
            })
        build.setdefault('variables', []).extend([
            {
                'intrinsic':  'procedure',
                'type':       template,
                'attributes': ['pointer'],
                'variables':  [function_label + 'Deferred'],
            },
            {
                'intrinsic':  'logical',
                'variables':  [function_label + 'IsAttchdVl=.false.'],
            },
        ])
        _created_pointers.add(function_label)


def reset_pointer_cache():
    """Clear `_created_pointers` — useful between test runs."""
    _created_pointers.clear()


def Properties_Deferred_Get_Functions(build, class_dict, member, prop):
    """Mirrors `Properties_Deferred_Get_Functions`."""
    attrs        = prop.get('attributes')  or {}
    get_function = prop.get('getFunction') or {}
    if not (
        'get' in _deferred_attributes(prop)
        and attrs.get('isGettable')
        and get_function.get('build')
    ):
        return

    type_descriptor, _ = data_object_definition(prop)
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
    function = {
        'type':        function_type + ' => propertyValue',
        'name':        (
            class_dict['name']
            + _ucfirst(member['name'])
            + _ucfirst(prop['name'])
            + 'Get'
        ),
        'description': (
            f"Get the value of the \\mono{{{prop['name']}}} property of "
            f"the \\mono{{{member['name']}}} implementation of the "
            f"\\mono{{{class_dict['name']}}} component using a deferred "
            "function."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
        ],
        'content': (
            f"propertyValue={class_dict['name']}{_ucfirst(member['name'])}"
            f"{_ucfirst(prop['name'])}GetDeferred(self)\n"
        ),
    }
    build.setdefault('types', {}).setdefault(impl_type, {}) \
                                  .setdefault('boundFunctions', []) \
                                  .append({
        'type':       'procedure',
        'name':       prop['name'],
        'descriptor': function,
    })
    _generate_deferred_attacher(member, prop, build, 'get')


def Properties_Deferred_Set_Functions(build, class_dict, member, prop):
    """Mirrors `Properties_Deferred_Set_Functions`."""
    attrs        = prop.get('attributes')  or {}
    set_function = prop.get('setFunction') or {}
    if not (
        'set' in _deferred_attributes(prop)
        and attrs.get('isSettable')
        and set_function.get('build')
    ):
        return

    type_descriptor, _ = data_object_definition(prop, match_only=True)
    type_descriptor.setdefault('variables',  []).append('setValue')
    type_descriptor.setdefault('attributes', []).append('intent(in   )')

    impl_type = (
        'nodeComponent'
        + _ucfirst(class_dict['name'])
        + _ucfirst(member['name'])
    )
    function = {
        'type':        'void',
        'name':        (
            class_dict['name']
            + _ucfirst(member['name'])
            + _ucfirst(prop['name'])
            + 'Set'
        ),
        'description': (
            f"Set the value of the \\mono{{{prop['name']}}} property of "
            f"the \\mono{{{member['name']}}} implementation of the "
            f"\\mono{{{class_dict['name']}}} component using a deferred "
            "function."
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
        'content': (
            f"call {class_dict['name']}{_ucfirst(member['name'])}"
            f"{_ucfirst(prop['name'])}SetDeferred(self,setValue)\n"
        ),
    }
    build.setdefault('types', {}).setdefault(impl_type, {}) \
                                  .setdefault('boundFunctions', []) \
                                  .append({
        'type':       'procedure',
        'name':       prop['name'] + 'Set',
        'descriptor': function,
    })
    _generate_deferred_attacher(member, prop, build, 'set')


def Properties_Deferred_Rate_Functions(build, class_dict, member, prop):
    """Mirrors `Properties_Deferred_Rate_Functions`."""
    attrs = prop.get('attributes') or {}
    if not (
        'rate' in _deferred_attributes(prop)
        and attrs.get('isEvolvable')
    ):
        return

    type_descriptor, _ = data_object_definition(prop, match_only=True)
    type_descriptor.setdefault('variables',  []).append('setValue')
    type_descriptor.setdefault('attributes', []).append('intent(in   )')

    cls_member  = class_dict['name'] + _ucfirst(member['name'])
    impl_type   = 'nodeComponent' + _ucfirst(class_dict['name']) + _ucfirst(member['name'])

    # Skip if a `<prop>Rate` already bound (e.g. from Set or another
    # path).  Mirrors Perl's `grep {...} boundFunctions` check.
    bound = build.setdefault('types', {}).setdefault(impl_type, {}) \
                                          .setdefault('boundFunctions', [])
    if any(b.get('name') == prop['name'] + 'Rate' for b in bound):
        return

    function = {
        'type':        'void',
        'name':        (
            class_dict['name']
            + _ucfirst(member['name'])
            + _ucfirst(prop['name'])
            + 'Rate'
        ),
        'description': (
            f"Accumulate the rate of change of the \\mono{{{prop['name']}}} "
            f"property of the \\mono{{{member['name']}}} implementation of "
            f"the \\mono{{{class_dict['name']}}} component using a deferred "
            "function."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            type_descriptor,
            {
                'intrinsic':  'logical',
                'attributes': ['intent(inout)', 'optional'],
                'variables':  ['interrupt'],
            },
            {
                'intrinsic':  'procedure',
                'type':       'interruptTask',
                'attributes': ['pointer', 'optional', 'intent(inout)'],
                'variables':  ['interruptProcedure'],
            },
        ],
        'content': (
            f"call {cls_member}{_ucfirst(prop['name'])}RateDeferred"
            "(self,setValue,interrupt,interruptProcedure)\n"
        ),
    }
    bound.append({
        'type':       'procedure',
        'name':       prop['name'] + 'Rate',
        'descriptor': function,
    })
    _generate_deferred_attacher(member, prop, build, 'rate')


def _generate_deferred_attacher(component, prop, build, method):
    """Emit `<...>Function` attacher and `<...>IsAttached` query
    methods on the component's nodeComponent type.

    Mirrors `Generate_Deferred_Function_Attacher`.
    """
    method_suffix = '' if method == 'get' else _ucfirst(method)
    component_class_name = component['class']
    component_name       = component['fullyQualifiedName']
    property_name        = prop['name']

    function_label = (
        _lcfirst(component_name)
        + _ucfirst(property_name)
        + _ucfirst(method)
    )

    type_name = 'nodeComponent' + _ucfirst(component_name)
    bound = build.setdefault('types', {}).setdefault(type_name, {}) \
                                          .setdefault('boundFunctions', [])
    target_name = property_name + method_suffix + 'Function'
    if any(b.get('name') == target_name for b in bound):
        return

    if method == 'get':
        attach_template = (
            component_name
            + _ucfirst(property_name)
            + _ucfirst(method)
        )
    else:
        attach_template = create_null_function(build, {
            'selfType':  component_class_name,
            'attribute': method,
            'property':  prop,
            'intent':    'inout',
        })

    attach_function = {
        'type':        'void',
        'name':        function_label + 'Function',
        'description': (
            f"Set the function to be used for the \\mono{{{method}}} "
            f"method of the \\mono{{{property_name}}} property of the "
            f"\\mono{{{component_name}}} component."
        ),
        'variables':   [
            {
                'intrinsic':  'procedure',
                'type':       attach_template,
                'variables':  ['deferredFunction'],
                'isArgument': True,
            },
        ],
        'content': (
            f"{function_label}Deferred        => deferredFunction\n"
            f"{function_label}IsAttchdVl =  .true.\n"
        ),
    }

    attach_status_function = {
        'type':        'logical',
        'name':        function_label + 'IsAttached',
        'description': (
            f"Return true if the deferred function used to {method} the "
            f"\\mono{{{property_name}}} property of the "
            f"\\mono{{{component_name}}} component class has been attached."
        ),
        'content': (
            f"{function_label}IsAttached={function_label}IsAttchdVl\n"
        ),
    }

    bound.append({
        'type':       'procedure',
        'pass':       'nopass',
        'name':       property_name + method_suffix + 'Function',
        'descriptor': attach_function,
    })
    bound.append({
        'type':       'procedure',
        'pass':       'nopass',
        'name':       property_name + method_suffix + 'IsAttached',
        'descriptor': attach_status_function,
    })


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


def _lcfirst(text):
    return text[:1].lower() + text[1:] if text else text


# ---------------------------------------------------------------------------
# Hook registration.  Order matches Perl Properties/Deferred.pm:24-28.
# ---------------------------------------------------------------------------

register('propertiesDeferred', 'propertyIteratedFunctions',
         Properties_Deferred_Pointers)
register('propertiesDeferred', 'propertyIteratedFunctions',
         Properties_Deferred_Get_Functions)
register('propertiesDeferred', 'propertyIteratedFunctions',
         Properties_Deferred_Set_Functions)
register('propertiesDeferred', 'propertyIteratedFunctions',
         Properties_Deferred_Rate_Functions)
