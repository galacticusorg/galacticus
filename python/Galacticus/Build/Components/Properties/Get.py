"""Per-property `<prop>` getters bound at the component level.

Andrew Benson (ported to Python 2026)
"""



from Galacticus.Build.Components.Utils      import register, is_intrinsic, offset_name
from Galacticus.Build.Components.DataTypes  import data_object_definition
from Galacticus.Build.Components.Properties.Evolve import _is_debugging, _evolvable


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
    non-virtual gettable property.
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
    # In debugging builds, prepend a guard to evolvable-property getters that aborts if the property's value is read
    # while active-property derivatives are being evaluated - during that phase an inactive property's stored value
    # is stale, so any such read is a bug (see issue #128). The guard is emitted only under -DDEBUGGING because value
    # getters are the hottest functions in the code and it must be free in production builds.
    if _is_debugging() and _evolvable(prop):
        guard, needs_count = _inactive_read_guard(class_dict, member, prop)
        function['content'] = guard + function['content']
        function['modules'] = ['Error']
        if needs_count:
            function['variables'].append({
                'intrinsic': 'integer',
                'variables': ['count'],
            })
    build.setdefault('types', {}).setdefault(impl_type, {}) \
                                  .setdefault('boundFunctions', []) \
                                  .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       prop['name'] + suffix,
    })


def _inactive_read_guard(class_dict, member, prop):
    """Build the Fortran guard block (as a string) that aborts if this
    evolvable property's value is read while active-property derivatives are
    being evaluated, together with a flag indicating whether an integer
    `count` local variable is required.

    The rank-1 / object-typed patterns mirror those in
    `Build_Analytic_Functions` (`Properties/Evolve.py`).
    """
    data  = prop.get('data') or {}
    ptype = data.get('type')
    rank  = int(data.get('rank') or 0)

    cls_member = class_dict['name'] + _ucfirst(member['name'])
    offset     = offset_name('all', cls_member, prop['name'])

    message = (
        f"'value of inactive property \"{prop['name']}\" of the "
        f"\"{member['name']}\" {class_dict['name']} component was read during "
        "active property evolution'//{introspection:location}"
    )

    needs_count = False
    if is_intrinsic(ptype):
        if rank == 0:
            condition = f"nodeInactives({offset})"
        else:
            condition = (
                f"any(nodeInactives({offset}:{offset}"
                f"+size(self%{prop['name']}Data)-1))"
            )
        guard = (
            "if (evaluationActiveRHS) then\n"
            f" if ({condition}) call Error_Report({message})\n"
            "end if\n"
        )
    else:
        needs_count = True
        guard = (
            "if (evaluationActiveRHS) then\n"
            f" count=self%{prop['name']}Data%serializeCount()\n"
            f" if (count > 0) then\n"
            f"  if (any(nodeInactives({offset}:{offset}+count-1))) "
            f"call Error_Report({message})\n"
            " end if\n"
            "end if\n"
        )
    return guard, needs_count


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


register('propertiesGet', 'propertyIteratedFunctions', Bind_Get_Functions)
register('propertiesGet', 'propertyIteratedFunctions', Build_Get_Functions)
