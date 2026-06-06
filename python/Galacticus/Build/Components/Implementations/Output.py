"""Per-implementation output methods: outputCount / outputNames / postOutput.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Components/Implementations/Output.pm.
"""

import re


from Galacticus.Build.Components.Utils import (
    register,
    is_output_intrinsic,
    output_type_map,
    _component_properties,
)


def _output_common_tasks(build, class_dict, member, function, type_suffixes):
    """Append common bookkeeping (modules, derived-type locals, unused-args
    line) to `function`.  Mirrors `Implementation_Output_Common_Tasks`.
    """
    # Determine which arguments are read by the body.  `extends` adds a
    # parent-call line that uses every argument, so all are flagged used.
    has_extends = isinstance(member.get('extends'), dict)
    arguments_used = {
        k: has_extends
        for k in ('self', 'time', 'instance', 'integer', 'double')
    }
    if (member.get('output') or {}).get('instances') == 'first':
        arguments_used['instance'] = True

    modules_required = []
    types_required   = []
    for prop in _component_properties(member):
        if 'output' not in prop or prop.get('definedInParent'):
            continue
        output = prop['output'] if isinstance(prop['output'], dict) else {}
        if 'modules' in output:
            module_list = re.sub(r'\s', '', str(output['modules']))
            for m in module_list.split(','):
                modules_required.append(m)
        ptype = (prop.get('data') or {}).get('type')
        if not is_output_intrinsic(ptype):
            types_required.append(ptype)
            arguments_used['self'] = True
            arguments_used['time'] = True
            arguments_used['integer'] = True
            arguments_used['double']  = True
        else:
            arguments_used[output_type_map[ptype]] = True

    # Append modules (sorted, uniqued).
    function.setdefault('modules', []).extend(sorted(set(modules_required)))

    # Append derived-type locals (sorted, uniqued).
    for tname in sorted(set(types_required)):
        function.setdefault('variables', []).append({
            'intrinsic': 'type',
            'type':      tname,
            'variables': [f'output{_ucfirst(tname)}'],
        })

    # Build the unused-arguments line.  For `integer`/`double`, expand to
    # one entry per `type_suffixes` element; for other names, just the
    # name itself.  Mirrors the inner nestedmap at Output.pm:378.
    unused_keys = [k for k in sorted(arguments_used.keys())
                   if not arguments_used[k]]
    unused_args = []
    for key in unused_keys:
        if key in ('integer', 'double'):
            for suffix in type_suffixes:
                unused_args.append(key + suffix)
        else:
            unused_args.append(key)
    if unused_args:
        function['content'] += (
            "!$GLC attributes unused :: " + ",".join(unused_args) + "\n"
        )


def Implementation_Output_Count(build, class_dict, member):
    """Generate `<class><Member>OutputCount`.

    Mirrors `Implementation_Output_Count`.
    """
    cap_class  = _ucfirst(class_dict['name'])
    cap_member = _ucfirst(member['name'])
    impl_type  = 'nodeComponent' + cap_class + cap_member

    function = {
        'type':        'void',
        'name':        class_dict['name'] + cap_member + 'OutputCount',
        'description': (
            f"Increment the count of properties to output for a "
            f"\\mono{{{member['name']}}} implementation of the "
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
            {
                'intrinsic':  'integer',
                'attributes': ['intent(inout)'],
                'variables':  ['integerPropertyCount', 'doublePropertyCount'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['intent(in   )'],
                'variables':  ['time'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['instance'],
            },
        ],
    }

    _output_common_tasks(build, class_dict, member, function, ['PropertyCount'])

    ext = member.get('extends')
    if isinstance(ext, dict):
        parent_type = (
            'nodeComponent'
            + _ucfirst(ext.get('class', ''))
            + _ucfirst(ext.get('name', ''))
        )
        function['content'] += (
            f"call self%{parent_type}%outputCount("
            "integerPropertyCount,doublePropertyCount,time,instance)\n"
        )

    if (member.get('output') or {}).get('instances') == 'first':
        function['content'] += "if (instance > 1) return\n"

    fixed_size_count = {'integer': 0, 'double': 0}

    for prop in _component_properties(member):
        if 'output' not in prop or prop.get('definedInParent'):
            continue
        data = prop.get('data') or {}
        ptype = data.get('type')
        rank  = int(data.get('rank') or 0)
        if is_output_intrinsic(ptype):
            output_type = output_type_map[ptype]
            if rank == 0:
                fixed_size_count[output_type] += 1
            elif rank == 1:
                output = prop['output'] if isinstance(prop['output'], dict) else {}
                labels = output.get('labels')
                m = re.match(r'^\[(.*)\]$', labels or '')
                if m:
                    count = m.group(1).count(',') + 1
                else:
                    count = output.get('count')
                # If count starts with a digit, it's a literal — fold it
                # into the fixed-size accumulator; otherwise it's a
                # variable expression and emit a runtime add.
                if isinstance(count, int) or (
                    isinstance(count, str) and count and count[0].isdigit()
                ):
                    fixed_size_count[output_type] += int(count)
                else:
                    function['content'] += (
                        f"{output_type}PropertyCount={output_type}PropertyCount"
                        f"+{count}\n"
                    )
        else:
            cap_type = _ucfirst(ptype)
            function['content'] += (
                f"output{cap_type}=self%{prop['name']}()\n"
                f"call output{cap_type}%outputCount(integerPropertyCount,"
                "doublePropertyCount,time)\n"
            )

    for tname, count in fixed_size_count.items():
        if count > 0:
            function['content'] += (
                f"{tname}PropertyCount={tname}PropertyCount+{count}\n"
            )

    _bind(build, impl_type, function, 'outputCount')


def Implementation_Output_Names(build, class_dict, member):
    """Generate `<class><Member>OutputNames`.

    Mirrors `Implementation_Output_Names`.
    """
    cap_class  = _ucfirst(class_dict['name'])
    cap_member = _ucfirst(member['name'])
    impl_type  = 'nodeComponent' + cap_class + cap_member

    function = {
        'type':        'void',
        'name':        class_dict['name'] + cap_member + 'OutputNames',
        'description': (
            f"Return the names of properties to output for a "
            f"\\mono{{{member['name']}}} implementation of the "
            f"\\mono{{{class_dict['name']}}} component."
        ),
        'content':     '',
        'modules':     [
            'Merger_Tree_Outputter_Buffer_Types',
            'Units_MetaData',
        ],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(inout)'],
                'variables':  ['integerProperty'],
            },
            {
                'intrinsic':  'type',
                'type':       'outputPropertyInteger',
                'attributes': ['intent(inout)', 'dimension(:)'],
                'variables':  ['integerProperties'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(inout)'],
                'variables':  ['doubleProperty'],
            },
            {
                'intrinsic':  'type',
                'type':       'outputPropertyDouble',
                'attributes': ['intent(inout)', 'dimension(:)'],
                'variables':  ['doubleProperties'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['intent(in   )'],
                'variables':  ['time'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['instance'],
            },
        ],
    }

    _output_common_tasks(build, class_dict, member, function, ['Property', 'Properties'])

    ext = member.get('extends')
    if isinstance(ext, dict):
        parent_type = (
            'nodeComponent'
            + _ucfirst(ext.get('class', ''))
            + _ucfirst(ext.get('name', ''))
        )
        function['content'] += (
            f"call self%{parent_type}%outputNames("
            "integerProperty,integerProperties,doubleProperty,"
            "doubleProperties,time,instance)\n"
        )

    if (member.get('output') or {}).get('instances') == 'first':
        function['content'] += "if (instance > 1) return\n"

    for prop in _component_properties(member):
        if 'output' not in prop or prop.get('definedInParent'):
            continue
        data = prop.get('data') or {}
        ptype = data.get('type')
        rank  = int(data.get('rank') or 0)
        output = prop['output'] if isinstance(prop['output'], dict) else {}
        units_in_si  = output.get('unitsInSI', '1.0d0')
        units_desc   = output.get('unitsDescription') or ''
        units_qty    = output.get('unitsQuantity')    or ''
        is_comoving  = output.get('isComoving')       or '.false.'
        comment      = output.get('comment', '')

        if is_output_intrinsic(ptype):
            ot = output_type_map[ptype]
            cap_prop = _ucfirst(prop['name'])
            if rank == 0:
                function['content'] += (
                    f"{ot}Property                                   = {ot}Property+1\n"
                    f"{ot}Properties({ot}Property)%name     ='{class_dict['name']}{cap_prop}'\n"
                    f"{ot}Properties({ot}Property)%comment  ='{comment}'\n"
                    f"{ot}Properties({ot}Property)%units    ="
                    f"unitType(unitsInSI={units_in_si},"
                    f"description='{units_desc}',quantity='{units_qty}',"
                    f"isComoving={is_comoving})\n"
                )
            elif rank == 1:
                labels = output.get('labels', '')
                m = re.match(r'^\[(.*)\]$', labels)
                if m:
                    label_list = re.sub(r'\s', '', m.group(1)).split(',')
                    for label in label_list:
                        function['content'] += (
                            f"{ot}Property                                   ={ot}Property+1\n"
                            f"{ot}Properties({ot}Property)%name     ='{class_dict['name']}{cap_prop}{label}'\n"
                            f"{ot}Properties({ot}Property)%comment  ='{comment} [{label}]'\n"
                            f"{ot}Properties({ot}Property)%units    ="
                            f"unitType(unitsInSI={units_in_si},"
                            f"description='{units_desc}',quantity='{units_qty}',"
                            f"isComoving={is_comoving})\n"
                        )
                elif 'count' in output:
                    label = re.sub(r'\{i\}', 'i', output.get('labels', ''))
                    function['content'] += (
                        f"do i=1,{output['count']}\n"
                        f"   {ot}Property                                   ={ot}Property+1\n"
                        f"   {ot}Properties({ot}Property)%name     ='{class_dict['name']}{cap_prop}'//{label}\n"
                        f"   {ot}Properties({ot}Property)%comment  ='{comment} [' //{label}//']'\n"
                        f"   {ot}Properties({ot}Property)%units    ="
                        f"unitType(unitsInSI={units_in_si},"
                        f"description='{units_desc}',quantity='{units_qty}',"
                        f"isComoving={is_comoving})\n"
                        "end do\n"
                    )
        else:
            units_in_si_d  = output.get('unitsInSI',         '1.0d0')
            units_desc_d   = output.get('unitsDescription', "''")
            units_qty_d    = output.get('unitsQuantity',    "''")
            cap_type       = _ucfirst(ptype)
            cap_prop       = _ucfirst(prop['name'])
            function['content'] += (
                f"output{cap_type}=self%{prop['name']}()\t\t\t   \n"
                f"call output{cap_type}%outputNames(integerProperty,"
                "integerProperties,doubleProperty,doubleProperties,time,"
                f"'{class_dict['name']}{cap_prop}','{comment}',"
                f"{units_in_si_d},'{units_desc_d}','{units_qty_d}')\n"
            )

    _bind(build, impl_type, function, 'outputNames')


def Implementation_Post_Output(build, class_dict, member):
    """Generate `<class><Member>PostOutput`.

    Mirrors `Implementation_Post_Output`.  Skipped entirely when no
    body would be emitted (no parent extends, no derived-type
    properties).
    """
    cap_class  = _ucfirst(class_dict['name'])
    cap_member = _ucfirst(member['name'])
    impl_type  = 'nodeComponent' + cap_class + cap_member

    function = {
        'type':        'void',
        'name':        class_dict['name'] + cap_member + 'PostOutput',
        'description': (
            f"Perform post-output processing for a \\mono{{{member['name']}}} "
            f"implementation of the \\mono{{{class_dict['name']}}} component."
        ),
        'content':     '',
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['intent(in   )'],
                'variables':  ['time'],
            },
        ],
    }

    ext = member.get('extends')
    if isinstance(ext, dict):
        parent_type = (
            'nodeComponent'
            + _ucfirst(ext.get('class', ''))
            + _ucfirst(ext.get('name', ''))
        )
        function['content'] += f"call self%{parent_type}%postOutput(time)\n"

    for prop in _component_properties(member):
        if 'output' not in prop or prop.get('definedInParent'):
            continue
        data = prop.get('data') or {}
        if not is_output_intrinsic(data.get('type')):
            function['content'] += (
                f"call self%{prop['name']}Data%postOutput(time)\n"
            )

    if function['content']:
        _bind(build, impl_type, function, 'postOutput')


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
# Hook registration.  Order matches Perl Implementations/Output.pm:25-28.
# ---------------------------------------------------------------------------

register('implementationsOutput', 'implementationIteratedFunctions', Implementation_Output_Count)
register('implementationsOutput', 'implementationIteratedFunctions', Implementation_Output_Names)
register('implementationsOutput', 'implementationIteratedFunctions', Implementation_Post_Output)
