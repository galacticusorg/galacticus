# Per-implementation state variable + sizeOf accessor.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Implementations/State.pm.

import os
import sys


from Galacticus.Build.Components.Utils import (
    register,
    is_intrinsic,
    _component_properties,
)


def Implementation_State(build, class_dict, member):
    """Declare `nodeComponent<Class><Member>IsActiveValue` module variable.

    Mirrors `Implementation_State`.
    """
    cap_class  = _ucfirst(class_dict['name'])
    cap_member = _ucfirst(member['name'])
    name       = f'nodeComponent{cap_class}{cap_member}IsActiveValue=.false.'
    build.setdefault('variables', []).append({
        'intrinsic':  'logical',
        'variables':  [name],
    })


def Implementation_Size_Of(build, class_dict, member):
    """Generate `nodeComponent<Class><Member>SizeOf` reporting the
    in-memory size.  Mirrors `Implementation_Size_Of`.
    """
    cap_class  = _ucfirst(class_dict['name'])
    cap_member = _ucfirst(member['name'])
    type_name  = f'nodeComponent{cap_class}{cap_member}'

    function = {
        'type':        'integer(c_size_t)',
        'name':        type_name + 'SizeOf',
        'description': (
            f"Return the size in bytes of a {member['name']} implementation "
            f"of the {class_dict['name']} component."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       type_name,
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
        ],
        'content': '',
    }
    content = f"{type_name}SizeOf=sizeof(self)\n"

    loop_iterator_required = False
    for prop in _component_properties(member):
        attrs = prop.get('attributes') or {}
        if attrs.get('isVirtual'):
            continue
        data = prop.get('data') or {}
        ptype = data.get('type')
        rank  = int(data.get('rank') or 0)
        if is_intrinsic(ptype):
            if rank > 0:
                content += (
                    f"{type_name}SizeOf={type_name}SizeOf"
                    f"+sizeof(self%{prop['name']}Data)\n"
                )
        else:
            if rank == 0:
                content += (
                    f"{type_name}SizeOf={type_name}SizeOf"
                    f"+self%{prop['name']}Data%nonStaticSizeOf()\n"
                )
            else:
                loop_iterator_required = True
                content += (
                    f"do i=1,size(self%{prop['name']}Data)\n"
                    f"   {type_name}SizeOf={type_name}SizeOf"
                    f"+self%{prop['name']}Data(i)%nonStaticSizeOf()\n"
                    "end do\n"
                )

    if loop_iterator_required:
        function['variables'].append({
            'intrinsic':  'integer',
            'type':       'c_size_t',
            'attributes': [],
            'variables':  ['i'],
        })

    function['content'] = content
    build.setdefault('types', {}).setdefault(type_name, {}) \
                                  .setdefault('boundFunctions', []) \
                                  .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       'sizeOf',
    })


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


register('implementationsState', 'implementationIteratedFunctions', Implementation_State)
register('implementationsState', 'implementationIteratedFunctions', Implementation_Size_Of)
