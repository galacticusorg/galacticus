# Per-implementation lifecycle methods: initialize / destroy / builder /
# createFunctionSet.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Implementations/CreateDestroy.pm.

import re


from Galacticus.Build.Components.Utils                  import (
    register,
    is_intrinsic,
    intrinsic_nulls,
    _component_properties,
)
from Galacticus.Build.Components.Classes.MetaProperties import meta_property_types


_SELF_REF_RE = re.compile(r'self([a-zA-Z]+)\s*%')


def Implementation_Creation(build, class_dict, member):
    """Generate `nodeComponent<Class><Member>Initialize`.

    Mirrors `Implementation_Creation`.  Allocates parent state, walks
    `selfXxx%` references in classDefault.code to add cross-component
    pointer locals, and either initialises every non-virtual property
    from its `classDefault.code` or sets it to the matching null value.
    Active classes also allocate the meta-property storage arrays.
    """
    name        = class_dict['name']
    cap_class   = _ucfirst(name)
    cap_member  = _ucfirst(member['name'])
    impl_type   = 'nodeComponent' + cap_class + cap_member

    function = {
        'type':        'void',
        'name':        impl_type + 'Initialize',
        'description': (
            f"Initialize a \\mono{{{member['name']}}} member of the "
            f"\\mono{{{name}}} component."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
        ],
    }

    component_initialization_lines = []
    required_components = set()
    for prop in _component_properties(member):
        attrs = prop.get('attributes') or {}
        class_default = prop.get('classDefault') or {}
        if attrs.get('isVirtual') or 'code' not in class_default:
            continue
        code = class_default['code']
        # Strip and discover every `selfXxx%` reference.  Match order
        # matches Perl's iterative `s///` which short-circuits on
        # already-seen names via `%requiredComponents`.
        cursor = code
        while True:
            m = _SELF_REF_RE.search(cursor)
            if m is None:
                break
            component_name = m.group(1)
            cursor = cursor[:m.start()] + cursor[m.end():]
            if component_name.lower() in required_components:
                continue
            required_components.add(component_name.lower())
            component_initialization_lines.append(
                f"self{component_name} => self%hostNode%"
                f"{component_name.lower()}()\n"
            )
            function['variables'].append({
                'intrinsic':  'class',
                'type':       'nodeComponent' + _ucfirst(component_name),
                'attributes': ['pointer'],
                'variables':  ['self' + _ucfirst(component_name)],
            })

    # Aggregate every classDefault.modules entry across this member's
    # properties.
    modules = []
    for prop in _component_properties(member):
        cd = prop.get('classDefault') or {}
        for m in cd.get('modules') or []:
            modules.append(m)
    if modules:
        function.setdefault('modules', []).extend(modules)

    has_extends    = isinstance(member.get('extends'), dict)
    has_real_props = any(
        not (p.get('attributes') or {}).get('isVirtual')
        for p in _component_properties(member)
    )

    content = ''
    if (
        not component_initialization_lines
        and not has_extends
        and not has_real_props
    ):
        content += "!$GLC attributes unused :: self\n"

    if has_extends:
        ext = member['extends']
        parent_type = (
            'nodeComponent'
            + _ucfirst(ext.get('class', ''))
            + _ucfirst(ext.get('name', ''))
        )
        content += f"call self%{parent_type}%initialize()\n"

    content += ''.join(component_initialization_lines)

    for prop in _component_properties(member):
        attrs = prop.get('attributes') or {}
        if attrs.get('isVirtual'):
            continue
        cd    = prop.get('classDefault') or {}
        data  = prop.get('data')         or {}
        ptype = data.get('type')
        rank  = int(data.get('rank') or 0)
        if 'code' in cd:
            if 'count' in cd:
                content += (
                    f"allocate(self%{prop['name']}Data({cd['count']}))\n"
                )
            content += f"self%{prop['name']}Data={cd['code']}\n"
        else:
            if (
                ptype in ('double', 'longInteger')
                and rank == 1
            ):
                # Mirror the Perl `join(",", "0" x $rank)`: when rank==1
                # this evaluates to `"0"` (a single zero) — so the Fortran
                # emits `allocate(...Data(0))`.  Faithfully reproduced.
                zeros = ','.join(['0'] * rank)
                content += (
                    f"allocate(self%{prop['name']}Data({zeros}))\n"
                )
            if is_intrinsic(ptype):
                content += (
                    f"self%{prop['name']}Data={intrinsic_nulls[ptype]}\n"
                )
            else:
                content += f"call self%{prop['name']}Data%reset()\n"

    if name in (build.get('componentClassListActive') or []):
        for mpt in meta_property_types:
            cap_label = _ucfirst(mpt['label'])
            rank      = mpt['rank']
            prefix    = f"{cap_label}Rank{rank}"
            var_pref  = f"{name}{prefix}"
            if rank == 0:
                if mpt['intrinsic'] == 'double precision':
                    base = '0.0'
                elif mpt['intrinsic'] == 'integer':
                    base = '0'
                else:
                    raise RuntimeError(
                        "Implementation_Creation: unknown meta-property type"
                    )
                if 'type' in mpt:
                    base += f"_{mpt['type']}"
                elif mpt['intrinsic'] == 'double precision':
                    base += "d0"
                initializer = f" self%{prefix}MetaProperties={base}"
            else:
                initializer = (
                    f" do i=1,size({var_pref}MetaPropertyNames)\n"
                    f"  allocate(self%{prefix}MetaProperties(i)%values(0))\n"
                    " end do\n"
                )
            content += (
                f"if (allocated({var_pref}MetaPropertyNames)"
                f".and..not.allocated(self%{prefix}MetaProperties)) then\n"
                f" allocate(self%{prefix}MetaProperties(size("
                f"{var_pref}MetaPropertyNames)))\n"
                f"{initializer}\n"
                "end if\n"
            )

    function['content'] = content
    _bind(build, impl_type, function, 'initialize')


def Implementation_Finalization(build, class_dict, member):
    """Generate `nodeComponent<Class><Member>Finalize`.

    Mirrors `Implementation_Finalization`.
    """
    cap_class   = _ucfirst(class_dict['name'])
    cap_member  = _ucfirst(member['name'])
    impl_type   = 'nodeComponent' + cap_class + cap_member

    function = {
        'type':        'void',
        'name':        impl_type + 'Finalize',
        'description': (
            f"Finalize a \\mono{{{member['name']}}} implementation of the "
            f"\\mono{{{class_dict['name']}}} component."
        ),
        'content':     '',
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
        ],
    }

    needs_real_finalize = any(
        not (p.get('attributes') or {}).get('isVirtual')
        and (
            int((p.get('data') or {}).get('rank') or 0) > 0
            or not is_intrinsic((p.get('data') or {}).get('type'))
        )
        for p in _component_properties(member)
    )
    if not needs_real_finalize:
        function['content'] += "!$GLC attributes unused :: self\n"

    if isinstance(member.get('extends'), dict):
        ext = member['extends']
        parent_type = (
            'nodeComponent'
            + _ucfirst(ext.get('class', ''))
            + _ucfirst(ext.get('name', ''))
        )
        function['content'] += f"call self%{parent_type}%destroy()\n"

    for prop in _component_properties(member):
        attrs = prop.get('attributes') or {}
        if attrs.get('isVirtual'):
            continue
        data  = prop.get('data') or {}
        ptype = data.get('type')
        rank  = int(data.get('rank') or 0)
        if not is_intrinsic(ptype):
            function['content'] += f"call self%{prop['name']}Data%destroy()\n"
        if rank > 0:
            function['content'] += (
                f"if (allocated(self%{prop['name']}Data)) "
                f"deallocate(self%{prop['name']}Data)\n"
            )

    _bind(build, impl_type, function, 'destroy')


def Implementation_Builder(build, class_dict, member):
    """Generate `nodeComponent<Class><Member>Builder`.

    Mirrors `Implementation_Builder`.  Walks the FoX_DOM XML
    sub-elements named after each non-virtual property and emits the
    appropriate `extractDataContent` / `%builder()` calls.  Active
    classes also walk the matching meta-property name list.
    """
    name        = class_dict['name']
    cap_class   = _ucfirst(name)
    cap_member  = _ucfirst(member['name'])
    impl_type   = 'nodeComponent' + cap_class + cap_member

    function = {
        'type':        'void',
        'name':        impl_type + 'Builder',
        'description': (
            f"Build a \\mono{{{member['name']}}} implementation of the "
            f"\\mono{{{class_dict['name']}}} component from a supplied "
            "XML definition."
        ),
        'modules':     [
            'Error, only : Error_Report',
            ('FoX_DOM, only : node, nodeList, getLength, '
             'extractDataContent, getElementsByTagName, item'),
            'ISO_Varying_String',
        ],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'type',
                'type':       'node',
                'attributes': ['intent(in   )', 'pointer'],
                'variables':  ['componentDefinition'],
            },
            {
                'intrinsic':  'type',
                'type':       'node',
                'attributes': ['pointer'],
                'variables':  ['property'],
            },
            {
                'intrinsic':  'type',
                'type':       'nodeList',
                'attributes': ['pointer'],
                'variables':  ['propertyList'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i', 'j', 'propertyListLength'],
            },
        ],
    }

    if member['name'] == 'null':
        content = "!$GLC attributes unused :: componentDefinition\n"
    else:
        content = "call self%initialize()\n"
        if isinstance(member.get('extends'), dict):
            ext = member['extends']
            parent_type = (
                'nodeComponent'
                + _ucfirst(ext.get('class', ''))
                + _ucfirst(ext.get('name', ''))
            )
            content += (
                f"call self%{parent_type}%builder(componentDefinition)\n"
            )

        for prop in _component_properties(member):
            attrs = prop.get('attributes') or {}
            if attrs.get('isVirtual'):
                continue
            data  = prop.get('data') or {}
            ptype = data.get('type')
            rank  = int(data.get('rank') or 0)

            content += (
                "!$omp critical (FoX_DOM_Access)\n"
                f"propertyList => getElementsByTagName(componentDefinition,"
                f"'{prop['name']}')\n"
                "!$omp end critical (FoX_DOM_Access)\n"
            )

            if rank == 0:
                content += (
                    "!$omp critical (FoX_DOM_Access)\n"
                    "propertyListLength=getLength(propertyList)\n"
                    "!$omp end critical (FoX_DOM_Access)\n"
                    "if (propertyListLength > 1) call Error_Report("
                    "'scalar property must have precisely one value'"
                    "//{introspection:location})\n"
                    "if (propertyListLength == 1) then\n"
                    "  !$omp critical (FoX_DOM_Access)\n"
                    "  property => item(propertyList,0)\n"
                    "  !$omp end critical (FoX_DOM_Access)\n"
                )
                if ptype in ('double', 'integer', 'logical'):
                    content += (
                        "  !$omp critical (FoX_DOM_Access)\n"
                        f"  call extractDataContent(property,self%{prop['name']}Data)\n"
                        "  !$omp end critical (FoX_DOM_Access)\n"
                    )
                elif ptype == 'longInteger':
                    content += (
                        "  call Error_Report('building of long integer "
                        "properties currently not supported'"
                        "//{introspection:location})\n"
                    )
                else:
                    content += (
                        f"  call self%{prop['name']}Data%builder(property)\n"
                    )
                content += "end if\n"
            elif rank == 1:
                content += (
                    "!$omp critical (FoX_DOM_Access)\n"
                    "propertyListLength=getLength(propertyList)\n"
                    "!$omp end critical (FoX_DOM_Access)\n"
                    "if (propertyListLength >= 1) then\n"
                )
                if ptype in ('double', 'integer', 'logical'):
                    content += (
                        f"  if (allocated(self%{prop['name']}Data)) "
                        f"deallocate(self%{prop['name']}Data)\n"
                        f"  allocate(self%{prop['name']}Data(propertyListLength))\n"
                        "  do i=1,propertyListLength\n"
                        "    !$omp critical (FoX_DOM_Access)\n"
                        "    property => item(propertyList,i-1)\n"
                        f"    call extractDataContent(property,self%{prop['name']}Data(i))\n"
                        "    !$omp end critical (FoX_DOM_Access)\n"
                        "  end do\n"
                    )
                elif ptype == 'longInteger':
                    content += (
                        "  call Error_Report('building of long integer "
                        "properties currently not supported'"
                        "//{introspection:location})\n"
                    )
                else:
                    content += (
                        f"  if (allocated(self%{prop['name']}Data)) "
                        f"deallocate(self%{prop['name']}Data)\n"
                        f"  allocate(self%{prop['name']}Data(propertyListLength))\n"
                        "  do i=1,propertyListLength\n"
                        "    !$omp critical (FoX_DOM_Access)\n"
                        "    property => item(propertyList,i-1)\n"
                        "    !$omp end critical (FoX_DOM_Access)\n"
                        f"    call self%{prop['name']}Data(i)%builder(property)\n"
                        "  end do\n"
                    )
                content += "end if\n"

    if name in (build.get('componentClassListActive') or []):
        for mpt in meta_property_types:
            cap_label = _ucfirst(mpt['label'])
            rank      = mpt['rank']
            prefix    = f"{cap_label}Rank{rank}"
            var_pref  = f"{name}{prefix}"
            content += (
                f"if (allocated({var_pref}MetaPropertyNames)) then\n"
                f" do i=1,size(({var_pref}MetaPropertyNames))\n"
                "  !$omp critical (FoX_DOM_Access)\n"
                f"  propertyList => getElementsByTagName(componentDefinition,"
                f"char({var_pref}MetaPropertyNames(i)))\n"
                "  propertyListLength=getLength(propertyList)\n"
                "  !$omp end critical (FoX_DOM_Access)\n"
            )
            is_long_integer = (
                mpt['intrinsic'] == 'integer'
                and mpt.get('type') == 'kind_int8'
            )
            if rank == 0:
                content += (
                    "  if (propertyListLength > 1) call Error_Report("
                    "'meta-property must have precisely one value'"
                    "//{introspection:location})\n"
                )
                if is_long_integer:
                    content += (
                        "  call Error_Report('building of long integer "
                        "properties currently not supported'"
                        "//{introspection:location})\n"
                    )
                else:
                    content += (
                        "  if (propertyListLength == 1) then\n"
                        "    !$omp critical (FoX_DOM_Access)\n"
                        "    property => item(propertyList,0)\n"
                        f"    call extractDataContent(property,self%{prefix}MetaProperties(i))\n"
                        "    !$omp end critical (FoX_DOM_Access)\n"
                        "  end if\n"
                    )
            elif rank == 1:
                if is_long_integer:
                    content += (
                        "  call Error_Report('building of long integer "
                        "properties currently not supported'"
                        "//{introspection:location})\n"
                    )
                else:
                    content += (
                        "  if (propertyListLength >= 1) then\n"
                        f"    if (allocated(self%{prefix}MetaProperties(i)%values)) "
                        f"deallocate(self%{prefix}MetaProperties(i)%values)\n"
                        f"    allocate(self%{prefix}MetaProperties(i)%values(propertyListLength))\n"
                        "    do j=1,propertyListLength\n"
                        "     !$omp critical (FoX_DOM_Access)\n"
                        "     property => item(propertyList,j-1)\n"
                        f"     call extractDataContent(property,self%{prefix}MetaProperties(i)%values(j))\n"
                        "     !$omp end critical (FoX_DOM_Access)\n"
                        "    end do\n"
                        "  end if\n"
                    )
            else:
                raise RuntimeError(
                    "Implementation_Builder: unsupported meta-property rank"
                )
            content += " end do\nend if\n"

    function['content'] = content
    _bind(build, impl_type, function, 'builder')


def Implementation_Deferred_Create_Set(build, class_dict, member):
    """Generate `<class><Member>CreateFunctionSet` when the member's
    `createFunction` is deferred.  Mirrors `Implementation_Deferred_Create_Set`.
    """
    if not (member.get('createFunction') or {}).get('isDeferred'):
        return

    name       = class_dict['name']
    cap_class  = _ucfirst(name)
    cap_member = _ucfirst(member['name'])

    function = {
        'type':        'void',
        'name':        f"{name}{cap_member}CreateFunctionSet",
        'description': (
            f"Set the create function for the \\mono{{{member['name']}}} "
            f"implementation of the \\mono{{{name}}} component class."
        ),
        'variables':   [
            {
                'intrinsic': 'external',
                'variables': ['createFunction'],
            },
        ],
        'content': f"{name}{cap_member}CreateFunction => createFunction\n",
    }

    impl_type = 'nodeComponent' + cap_class + cap_member
    build.setdefault('types', {}).setdefault(impl_type, {}) \
                                  .setdefault('boundFunctions', []) \
                                  .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       'createFunctionSet',
        'pass':       'nopass',
    })

    # Module-scope function pointer that the setter assigns.
    build.setdefault('variables', []).append({
        'intrinsic':  'procedure',
        'type':       '',
        'attributes': ['pointer'],
        'variables':  [f"{name}{cap_member}CreateFunction"],
    })


def _bind(build, type_name, function, method_name):
    build.setdefault('types', {}).setdefault(type_name, {}) \
                                  .setdefault('boundFunctions', []) \
                                  .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       method_name,
    })


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


# ---------------------------------------------------------------------------
# Hook registration.  Order matches Perl Implementations/CreateDestroy.pm:24-28.
# ---------------------------------------------------------------------------

register('implementationsCreateDestroy', 'implementationIteratedFunctions',
         Implementation_Creation)
register('implementationsCreateDestroy', 'implementationIteratedFunctions',
         Implementation_Builder)
register('implementationsCreateDestroy', 'implementationIteratedFunctions',
         Implementation_Finalization)
register('implementationsCreateDestroy', 'implementationIteratedFunctions',
         Implementation_Deferred_Create_Set)
