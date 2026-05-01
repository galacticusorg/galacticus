# Driver for the `classIteratedFunctions` sub-iteration phase.
# Andrew Benson (ported to Python 2026)
#
# Mirrors the small portion of perl/Galacticus/Build/Components/Classes/Utils.pm
# that runs as part of the components-build pipeline:
# `Class_Function_Iterator`.  This is itself a `functions`-phase hook that,
# when invoked, walks the `component_utils` registry looking for
# `classIteratedFunctions` entries and dispatches each function once per
# component class.
#
# The remaining functions in the Perl Classes::Utils module
# (`Class_Move`, `Class_Remove`, the main `functions`-phase group) belong
# to a later port stage.

import os
import sys


from Galacticus.Build.Components.Utils import register, component_utils, verbosity_level


def Class_Function_Iterator(build):
    """Run every `classIteratedFunctions` hook once per component class.

    Mirrors `Class_Function_Iterator`.  Owners may register a list of
    functions under the `classIteratedFunctions` key; each function is
    called as `fn(build, class_dict)` for every entry in
    `build['componentClasses']`.
    """
    if verbosity_level >= 1:
        # Match Perl's leading-six-space report indent.
        pass

    for owner_name in sorted(component_utils.keys()):
        owner = component_utils[owner_name]
        functions = owner.get('classIteratedFunctions')
        if not functions:
            continue
        if not isinstance(functions, list):
            functions = [functions]
        for fn in functions:
            marker = (
                f" {{{getattr(fn, '__name__', '<fn>')}}}"
                if len(functions) > 1 else ''
            )
            print(f"         --> {owner_name}{marker}")
            for class_dict in (build.get('componentClasses') or {}).values():
                fn(build, class_dict)


def Class_Move(build):
    """Generate `nodeComponent<Class>Move` per class — moves all
    instances of a component class from one node to another, optionally
    overwriting any pre-existing target instances.

    Mirrors `Class_Move`.  Bound to `treeNode`, not the class type
    itself, since the method takes a target `treeNode` argument.
    """
    active = set(build.get('componentClassListActive') or [])
    for class_dict in (build.get('componentClasses') or {}).values():
        name      = class_dict['name']
        cap       = _ucfirst(name)
        type_name = 'nodeComponent' + cap

        function = {
            'type':        'void',
            'name':        type_name + 'Move',
            'description': (
                f"Move instances of the \\mono{{{name}}} component, "
                "from one node to another."
            ),
            'modules':     ['Error'],
            'variables':   [
                {
                    'intrinsic':  'class',
                    'type':       'treeNode',
                    'attributes': ['intent(inout)'],
                    'variables':  ['self'],
                },
                {
                    'intrinsic':  'type',
                    'type':       'treeNode',
                    'attributes': ['intent(inout)', 'target'],
                    'variables':  ['targetNode'],
                },
                {
                    'intrinsic':  'logical',
                    'attributes': ['intent(in   )', 'optional'],
                    'variables':  ['overwrite'],
                },
                {
                    'intrinsic':  'integer',
                    'variables':  ['instanceCount', 'targetCount', 'i'],
                },
                {
                    'intrinsic':  'class',
                    'type':       type_name,
                    'attributes': ['allocatable, dimension(:)'],
                    'variables':  ['instancesTemporary'],
                },
                {
                    'intrinsic':  'logical',
                    'variables':  ['overwrite_'],
                },
            ],
        }

        if name in active:
            content  = (
                "overwrite_=.false.\n"
                "if (present(overwrite)) overwrite_=overwrite\n"
                f"instanceCount=self      %{name}count()\n"
                f"targetCount  =targetNode%{name}count()\n"
                "if (overwrite_ .and. targetCount > 0) then\n"
                f"  do i=1,targetCount\n"
                f"    call targetNode%component{cap}(i)%destroy()\n"
                "  end do \n"
                "  targetCount=0\n"
                f"  deallocate(targetNode%component{cap})\n"
                f"  allocate(targetNode%component{cap}(1))\n"
                "end if\t\n"
                "if (instanceCount == 0) return\n"
                "if (targetCount == 0) then\n"
                f"  deallocate(targetNode%component{cap})\n"
                f"  call Move_Alloc(self%component{cap},targetNode%component{cap})\n"
                "else\n"
                "  ! Multiple instances, so remove the specified instance.\n"
                f"  allocate(instancesTemporary(instanceCount+targetCount),"
                f"source=self%component{cap}(1))\n"
            )
            for member in class_dict.get('members') or []:
                impl_type = type_name + _ucfirst(member['name'])
                content += (
                    f"  select type (from => targetNode%component{cap})\n"
                    f"  type is ({impl_type})\n"
                    f"    select type (to => instancesTemporary)\n"
                    f"    type is ({impl_type})\n"
                    f"      to(1:targetCount)=from\n"
                    f"    end select\n"
                    f"  end select\n"
                )
            for member in class_dict.get('members') or []:
                impl_type = type_name + _ucfirst(member['name'])
                content += (
                    f"   select type (from => self%component{cap})\n"
                    f"   type is ({impl_type})\n"
                    f"     select type (to => instancesTemporary)\n"
                    f"     type is ({impl_type})\n"
                    f"       to(targetCount+1:targetCount+instanceCount)=from\n"
                    f"     end select\n"
                    f"   end select\n"
                )
            content += (
                f"  call targetNode%{name}Destroy()\n"
                f"  call self      %{name}Destroy()\n"
                f"  call Move_Alloc(instancesTemporary,targetNode%component{cap})\n"
                f"  allocate(self%component{cap}(1))\n"
                "end if\n"
                f"do i=1,size(targetNode%component{cap})\n"
                f"   targetNode%component{cap}(i)%hostNode => targetNode\n"
                "end do\n"
            )
        else:
            content = (
                "!$GLC attributes unused :: self, targetNode, overwrite, "
                "instanceCount, targetCount, i, instancesTemporary, overwrite_\n"
                "call Error_Report('Galacticus was not compiled with support "
                "for this class'//{introspection:location})\n"
            )
        function['content'] = content

        build.setdefault('types', {}).setdefault('treeNode', {}) \
                                     .setdefault('boundFunctions', []) \
                                     .append({
            'type':       'procedure',
            'descriptor': function,
            'name':       name + 'Move',
            'returnType': r"\void",
            'arguments':  (
                r"\textcolor{red}{\textless type(treeNode)\textgreater} "
                r"targetNode\arginout"
            ),
        })


def Class_Remove(build):
    """Generate `nodeComponent<Class>Remove` per class — removes one
    indexed instance of a component class from a node.

    Mirrors `Class_Remove`.
    """
    active = set(build.get('componentClassListActive') or [])
    for class_dict in (build.get('componentClasses') or {}).values():
        name      = class_dict['name']
        cap       = _ucfirst(name)
        type_name = 'nodeComponent' + cap

        function = {
            'type':        'void',
            'name':        type_name + 'Remove',
            'description': (
                f"Remove an instance of the \\mono{{{name}}} component "
                "from a node."
            ),
            'modules':     ['Error'],
            'variables':   [
                {
                    'intrinsic':  'class',
                    'type':       'treeNode',
                    'attributes': ['intent(inout)'],
                    'variables':  ['self'],
                },
                {
                    'intrinsic':  'integer',
                    'attributes': ['intent(in   )'],
                    'variables':  ['instance'],
                },
                {
                    'intrinsic':  'integer',
                    'variables':  ['instanceCount'],
                },
                {
                    'intrinsic':  'class',
                    'type':       type_name,
                    'attributes': ['allocatable, dimension(:)'],
                    'variables':  ['instancesTemporary'],
                },
            ],
        }

        if name in active:
            content = (
                f"instanceCount=self%{name}count()\n"
                "if (instance < 1 .or. instance > instanceCount) "
                "call Error_Report('instance out of range'"
                "//{introspection:location})\n"
                f"call self%component{cap}(instance)%destroy()\n"
                "if (instanceCount == 1) then\n"
                "  ! Only one instance of this component. Deallocate it and "
                "reallocate with generic type.\n"
                f"  deallocate(self%component{cap})\n"
                f"  allocate(self%component{cap}(1))\n"
                "else\n"
                "  ! Multiple instances, so remove the specified instance.\n"
                f"  allocate(instancesTemporary(instanceCount-1),"
                f"source=self%component{cap}(1))\n"
            )
            for member in class_dict.get('members') or []:
                impl_type = type_name + _ucfirst(member['name'])
                content += (
                    f"  select type (from => self%component{cap})\n"
                    f"  type is ({impl_type})\n"
                    f"    select type (to => instancesTemporary)\n"
                    f"    type is ({impl_type})\n"
                    "      if (instance >             1) "
                    "to(       1:instance     -1)=from(         1:instance     -1)\n"
                    "      if (instance < instanceCount) "
                    "to(instance:instanceCount-1)=from(instance+1:instanceCount  )\n"
                    f"    end select\n"
                    f"  end select\n"
                )
            content += (
                f"  deallocate(self%component{cap})\n"
                f"  call Move_Alloc(instancesTemporary,self%component{cap})\n"
                "end if\n"
            )
        else:
            content = (
                "!$GLC attributes unused :: self, instance, instanceCount, "
                "instancesTemporary\n"
                "call Error_Report('Galacticus was not compiled with support "
                "for this class'//{introspection:location})\n"
            )
        function['content'] = content

        build.setdefault('types', {}).setdefault('treeNode', {}) \
                                     .setdefault('boundFunctions', []) \
                                     .append({
            'type':       'procedure',
            'descriptor': function,
            'name':       name + 'Remove',
        })


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


# Hook registration order matches Perl Classes/Utils.pm:23-27.
register('classUtils', 'functions', Class_Move)
register('classUtils', 'functions', Class_Remove)
register('classUtils', 'functions', Class_Function_Iterator)
