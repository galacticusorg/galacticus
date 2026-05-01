# Per-implementation deferred-binding hooks.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Implementations/Deferred.pm.
#
# Four `implementationIteratedFunctions` hooks per (class, member):
#
#   Implementation_Deferred_Binding_Pointers      — module-scope
#       procedure pointer + IsSet flag for each deferred binding.
#   Implementation_Deferred_Binding_Attachers     — `<...>DfrrdFnctnSet`
#       method to attach a function to a deferred binding.
#   Implementation_Deferred_Binding_Attach_Status — `<...>DfrrdFnctnIsSet`
#       method to query attachment status.
#   Implementation_Deferred_Binding_Wrappers      — wrapper method that
#       calls through the deferred function pointer when set, falls
#       back to a parent class's binding (if any), or emits an
#       `Error_Report` if no fallback exists.

import os
import sys


from Fortran.Utils                                  import unformat_variables
from Galacticus.Build.Components.Utils              import (
    register,
    is_intrinsic,
    intrinsic_types,
    intrinsic_nulls,
    argument_list,
)


def _bindings(member):
    """Return the list of bindings on `member`, normalised to a list."""
    bindings = (member.get('bindings') or {}).get('binding') or []
    if not isinstance(bindings, list):
        bindings = [bindings]
    return [b for b in bindings if isinstance(b, dict)]


def _walk_to_template_member(member, method_name):
    """Walk the `extends.implementation` chain upward as long as each
    parent has a binding for `method_name`; return the topmost member.

    Mirrors the `while ( exists($baseMember->{'extends'}) && grep ... )`
    loop in Deferred.pm.  The Perl code uses `$baseMember = $baseMember->{'extends'}`
    on the assignment side (a quirk in the original — `$baseMember` ends
    up as the `extends` *dict*, not the parent member dict, so subsequent
    `$baseMember->{'name'}` reads the parent's class entry's `name`
    rather than the parent member's `name`).  Mirrored here verbatim.
    """
    base_member = member
    while True:
        ext = base_member.get('extends')
        if not isinstance(ext, dict):
            break
        parent = ext.get('implementation')
        if not isinstance(parent, dict):
            break
        if not any(b.get('method') == method_name for b in _bindings(parent)):
            break
        base_member = ext  # Faithful Perl-ism — see docstring.
    return base_member


def Implementation_Deferred_Binding_Pointers(build, class_dict, member):
    """Mirrors `Implementation_Deferred_Binding_Pointers`."""
    name      = class_dict['name']
    cap_member = _ucfirst(member['name'])
    for binding in _bindings(member):
        if not binding.get('isDeferred'):
            continue
        method_name = binding['method']
        base_member = _walk_to_template_member(member, method_name)
        component_function_name = (
            name + cap_member + _ucfirst(method_name)
        )
        template_type = (
            name + _ucfirst(base_member.get('name', '')) + _ucfirst(method_name)
        )
        build.setdefault('variables', []).extend([
            {
                'intrinsic':  'procedure',
                'type':       template_type,
                'attributes': ['pointer'],
                'variables':  [component_function_name + 'Deferred'],
            },
            {
                'intrinsic':  'logical',
                'variables':  [component_function_name + 'IsSetValue=.false.'],
            },
        ])


def Implementation_Deferred_Binding_Attachers(build, class_dict, member):
    """Mirrors `Implementation_Deferred_Binding_Attachers`."""
    name        = class_dict['name']
    cap_member  = _ucfirst(member['name'])
    impl_type   = 'nodeComponent' + _ucfirst(name) + cap_member
    for binding in _bindings(member):
        if not binding.get('isDeferred'):
            continue
        method_name = binding['method']
        base_member = _walk_to_template_member(member, method_name)
        member_function_name = name + cap_member + _ucfirst(method_name)
        template_type = (
            name + _ucfirst(base_member.get('name', '')) + _ucfirst(method_name)
        )

        function = {
            'type':        'void',
            'name':        member_function_name + 'DfrrdFnctnSet',
            'description': (
                f"Set the function to be used for the \\mono{{{method_name}}} "
                f"method of the \\mono{{{member['name']}}} implementation of "
                f"the \\mono{{{name}}} component class."
            ),
            'variables':   [
                {
                    'intrinsic':  'procedure',
                    'type':       template_type,
                    'isArgument': True,
                    'variables':  ['deferredFunction'],
                },
            ],
            'content': (
                f"{member_function_name}Deferred   => deferredFunction\n"
                f"{member_function_name}IsSetValue =  .true.\n"
            ),
        }
        build.setdefault('types', {}).setdefault(impl_type, {}) \
                                      .setdefault('boundFunctions', []) \
                                      .append({
            'type':       'procedure',
            'pass':       'nopass',
            'descriptor': function,
            'name':       method_name + 'Function',
        })


def Implementation_Deferred_Binding_Attach_Status(build, class_dict, member):
    """Mirrors `Implementation_Deferred_Binding_Attach_Status`."""
    name        = class_dict['name']
    cap_member  = _ucfirst(member['name'])
    impl_type   = 'nodeComponent' + _ucfirst(name) + cap_member
    for binding in _bindings(member):
        if not binding.get('isDeferred'):
            continue
        method_name = binding['method']
        member_function_name = name + cap_member + _ucfirst(method_name)
        function = {
            'type':        'logical',
            'name':        member_function_name + 'DfrrdFnctnIsSet',
            'description': (
                f"Return true if the deferred function for the "
                f"\\mono{{{method_name}}} method of the "
                f"\\mono{{{member['name']}}} implementation of the "
                f"\\mono{{{name}}} component class has been set."
            ),
            'content': (
                f"{member_function_name}DfrrdFnctnIsSet="
                f"{member_function_name}IsSetValue\n"
            ),
        }
        build.setdefault('types', {}).setdefault(impl_type, {}) \
                                      .setdefault('boundFunctions', []) \
                                      .append({
            'type':       'procedure',
            'pass':       'nopass',
            'descriptor': function,
            'name':       method_name + 'FunctionIsSet',
        })


def Implementation_Deferred_Binding_Wrappers(build, class_dict, member):
    """Mirrors `Implementation_Deferred_Binding_Wrappers`."""
    name        = class_dict['name']
    cap         = _ucfirst(name)
    cap_member  = _ucfirst(member['name'])
    impl_type   = 'nodeComponent' + cap + cap_member

    for binding in _bindings(member):
        if not binding.get('isDeferred'):
            continue
        method_name              = binding['method']
        member_function_name     = name + cap_member + _ucfirst(method_name)
        return_name              = member_function_name
        interface                = binding.get('interface') or {}
        interface_type           = interface.get('type', 'void')
        specific_type            = (
            interface_type != 'void'
            and interface_type not in intrinsic_types
        )
        is_pointer               = (
            specific_type and ', pointer' in interface_type
        )
        if specific_type:
            return_name += '_'

        if interface_type == 'void':
            function_type = 'void'
        elif specific_type:
            function_type = f"{interface_type} => {return_name}"
        else:
            function_type = intrinsic_types[interface_type]

        function = {
            'type':        function_type,
            'name':        member_function_name,
            'description': (
                f"Call the deferred function for the \\mono{{{method_name}}} "
                f"method of the \\mono{{{name}}} component class if it has "
                "been set."
            ),
            'modules':     ['Error'],
            'variables':   [],
        }
        # Optional extra modules from the binding's interface.
        for module in _as_list(interface.get('module')):
            function['modules'].append(module)

        # Handle rank/shape on the return type.
        rank  = interface.get('rank')
        shape = interface.get('shape')
        if rank is not None and int(rank) > 0:
            if shape is not None:
                raise RuntimeError('can not specify both "rank" and "shape"')
            return_name += '_'
            dims = ','.join([':'] * int(rank))
            function['type'] += (
                f", allocatable, dimension({dims}) => {return_name}"
            )
        elif shape is not None:
            if rank is not None:
                raise RuntimeError('can not specify both "rank" and "shape"')
            return_name += '_'
            function['type'] += f", dimension({shape}) => {return_name}"

        # `self` argument if interface.self.pass is set.
        if (interface.get('self') or {}).get('pass'):
            self_intent = (interface['self'].get('intent') or 'in')
            function['variables'].append({
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': [f'intent({self_intent})'],
                'variables':  ['self'],
            })

        # Other arguments.
        for arg_string in _as_list(interface.get('argument')):
            decl = unformat_variables(arg_string)
            if decl is not None:
                function['variables'].append(decl)

        # Argument-list extracted via the standard helper.
        arguments = argument_list(*function['variables'])

        # Body.
        content = (
            "select type (self)\n"
            f"class is ({impl_type})\n"
            f"   if (self%{method_name}FunctionIsSet()) then\n"
        )
        if interface_type == 'void':
            content += (
                f"      call {member_function_name}Deferred("
                + ",".join(arguments) + ")\n"
            )
        else:
            assigner = ' => ' if is_pointer else '='
            content += (
                f"      {return_name}{assigner}{member_function_name}Deferred("
                + ",".join(arguments) + ")\n"
            )
        content += "   else\n"

        # Parent-class fallback.
        parent_type = None
        ext = member.get('extends')
        if isinstance(ext, dict):
            parent_impl = ext.get('implementation') or {}
            if any(b.get('method') == method_name for b in _bindings(parent_impl)):
                parent_type = (
                    'nodeComponent'
                    + _ucfirst(ext.get('class', ''))
                    + _ucfirst(ext.get('name', ''))
                )

        if parent_type is not None:
            args_no_self = ",".join(a for a in arguments if a != 'self')
            if interface_type == 'void':
                content += (
                    f"      call self%{parent_type}%{method_name}({args_no_self})\n"
                )
            else:
                assigner = ' => ' if is_pointer else '='
                content += (
                    f"      {return_name}{assigner}self%{parent_type}"
                    f"%{method_name}({args_no_self})\n"
                )
        else:
            if (
                is_intrinsic(interface_type)
                and interface_type != 'void'
            ):
                null_value = intrinsic_nulls[interface_type]
                content += f"      {return_name}={null_value}\n"
            content += (
                "      call Error_Report('deferred function has not been "
                "assigned'//{introspection:location})\n"
            )
        content += "   end if\nclass default\n"

        if (
            is_intrinsic(interface_type)
            and interface_type != 'void'
        ):
            null_value = intrinsic_nulls[interface_type]
            content += f"      {return_name}={null_value}\n"
        content += (
            "   call Error_Report('incorrect class - this should not happen'"
            "//{introspection:location})\n"
            "end select\n"
        )

        function['content'] = content

        build.setdefault('types', {}).setdefault(impl_type, {}) \
                                      .setdefault('boundFunctions', []) \
                                      .append({
            'type':       'procedure',
            'descriptor': function,
            'name':       method_name,
        })


def _as_list(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


# ---------------------------------------------------------------------------
# Hook registration.  Order matches Perl Implementations/Deferred.pm:23-27.
# ---------------------------------------------------------------------------

register('implementationsDeferred', 'implementationIteratedFunctions',
         Implementation_Deferred_Binding_Pointers)
register('implementationsDeferred', 'implementationIteratedFunctions',
         Implementation_Deferred_Binding_Attachers)
register('implementationsDeferred', 'implementationIteratedFunctions',
         Implementation_Deferred_Binding_Attach_Status)
register('implementationsDeferred', 'implementationIteratedFunctions',
         Implementation_Deferred_Binding_Wrappers)
