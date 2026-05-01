# Component-to-component assignment operator.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Hierarchy/Utils.pm.



from Galacticus.Build.Components.Utils      import register, is_intrinsic


def Component_Assign(build):
    """Generate the `assign` method (overloads `assignment(=)`) on
    `nodeComponent`.

    Mirrors `Component_Assign`.  Walks every (class, implementation,
    property) triple to emit the right per-property copy code.
    """
    function = {
        'type':        'void',
        'name':        'nodeComponentAssign',
        'methodName':  'assignment(=)',
        'description': "Assign a node component to another node component.",
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'nodeComponent',
                'attributes': ['intent(  out)'],
                'variables':  ['to'],
            },
            {
                'intrinsic':  'class',
                'type':       'nodeComponent',
                'attributes': ['intent(in   )'],
                'variables':  ['from'],
            },
        ],
    }

    content = (
        "to%hostNode => from%hostNode\n"
        "select type (to)\n"
    )
    for class_dict in (build.get('componentClasses') or {}).values():
        for member in class_dict.get('members') or []:
            impl_type = (
                'nodeComponent'
                + _ucfirst(class_dict['name'])
                + _ucfirst(member['name'])
            )
            content += (
                f"type is ({impl_type})\n"
                f"   select type (from)\n"
                f"   type is ({impl_type})\n"
            )
            for prop in _properties(member):
                if prop.get('attributes', {}).get('isVirtual'):
                    continue
                ptype = (prop.get('data') or {}).get('type')
                rank  = int((prop.get('data') or {}).get('rank') or 0)
                if is_intrinsic(ptype):
                    if ptype == 'double' and rank > 0:
                        content += (
                            f"      if (allocated(to%{prop['name']}Data)) "
                            f"deallocate(to%{prop['name']}Data)\n"
                        )
                else:
                    content += f"      call to%{prop['name']}Data%destroy()\n"
                content += (
                    f"      to%{prop['name']}Data=from%{prop['name']}Data\n"
                )
            content += "   end select\n"
    content += "end select\n"

    function['content'] = content

    bound = build.setdefault('types', {}).setdefault('nodeComponent', {}) \
                                          .setdefault('boundFunctions', [])
    bound.append({
        'type':       'procedure',
        'descriptor': function,
        'name':       'assign',
    })
    bound.append({
        'type':     'generic',
        'name':     'assignment(=)',
        'function': 'assign',
    })


def _properties(member):
    """Yield every property dict declared on `member`."""
    props = (member.get('properties') or {}).get('property')
    if props is None:
        return []
    if isinstance(props, list):
        return props
    if isinstance(props, dict):
        if all(isinstance(v, dict) for v in props.values()):
            return list(props.values())
        return [props]
    return []


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


register('hierarchyUtils', 'functions', Component_Assign)
