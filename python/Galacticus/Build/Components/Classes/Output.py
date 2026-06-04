"""Per-class output methods: dumpASCII, outputCount, outputNames,
postOutput, output.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Components/Classes/Output.pm.
"""

import re


from Galacticus.Build.Components       import Utils as _Utils
from Galacticus.Build.Components.Utils import (
    register,
    is_output_intrinsic,
    _component_properties,
)


# Output.pm:278-283 — map from property type to (label, intrinsic, type).
_INTRINSIC_TYPE_MAP = {
    'double':       {'label': 'double',  'intrinsic': 'double precision'                       },
    'integer':      {'label': 'integer', 'intrinsic': 'integer',          'type': 'kind_int8'  },
    'longInteger':  {'label': 'integer', 'intrinsic': 'integer',          'type': 'kind_int8'  },
}


def Class_Dump_ASCII(build, class_dict):
    """Generate `nodeComponent<Class>DumpASCII`.  Mirrors `Class_Dump_ASCII`."""
    name = class_dict['name']
    if name not in (build.get('componentClassListActive') or []):
        return
    type_name = 'nodeComponent' + _ucfirst(name)
    padding   = ' ' * max(
        (_Utils.fully_qualified_name_length_max or 0) - len(name), 0
    )
    function = {
        'type':        'void',
        'name':        type_name + 'DumpASCII',
        'description': f"Dump the content of a \\mono{{{name}}} component.",
        'modules':     ['Display', 'ISO_Varying_String'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       type_name,
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
        ],
        'content': (
            "!$GLC attributes unused :: self\n"
            f"call displayIndent('{name}: {padding}generic')\n"
            "call displayUnindent('done')\n"
        ),
    }
    _bind(build, type_name, function, 'dumpASCII')


def Class_Output_Count(build, class_dict):
    """Generate `<class>OutputCount` — delegates to the default
    component instance.  Mirrors `Class_Output_Count`.
    """
    name = class_dict['name']
    if name not in (build.get('componentClassListActive') or []):
        return
    cap       = _ucfirst(name)
    type_name = 'nodeComponent' + cap
    function = {
        'type':        'void',
        'name':        name + 'OutputCount',
        'description': (
            f"Increment the count of properties to output for a generic "
            f"\\mono{{{name}}} component."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       type_name,
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
            {
                'intrinsic':  'class',
                'type':       type_name,
                'attributes': ['allocatable'],
                'variables':  ['selfDefault'],
            },
        ],
        'content': (
            f"allocate(selfDefault,source=default{cap}Component)\n"
            "selfDefault%hostNode => self%hostNode\n"
            "call selfDefault%outputCount(integerPropertyCount,doublePropertyCount,"
            "time,instance)\n"
        ),
    }
    _bind(build, type_name, function, 'outputCount')


def Class_Output_Names(build, class_dict):
    """Generate `<class>OutputNames`.  Mirrors `Class_Output_Names`."""
    name = class_dict['name']
    if name not in (build.get('componentClassListActive') or []):
        return
    cap       = _ucfirst(name)
    type_name = 'nodeComponent' + cap
    function = {
        'type':        'void',
        'name':        name + 'OutputNames',
        'description': (
            f"Establish the names of properties to output for a generic "
            f"\\mono{{{name}}} component."
        ),
        'modules':     ['Merger_Tree_Outputter_Buffer_Types'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       type_name,
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
            {
                'intrinsic':  'class',
                'type':       type_name,
                'attributes': ['allocatable'],
                'variables':  ['selfDefault'],
            },
        ],
        'content': (
            f"allocate(selfDefault,source=default{cap}Component)\n"
            "selfDefault%hostNode => self%hostNode\n"
            "call selfDefault%outputNames(integerProperty,integerProperties,"
            "doubleProperty,doubleProperties,time,instance)\n"
        ),
    }
    _bind(build, type_name, function, 'outputNames')


def Class_Output(build, class_dict):
    """Generate `<class>Output` — populate output buffers from the
    properties of every active member.  Mirrors `Class_Output`.
    """
    name = class_dict['name']
    if name not in (build.get('componentClassListActive') or []):
        return
    cap       = _ucfirst(name)
    type_name = 'nodeComponent' + cap

    base_variables = [
        {
            'intrinsic':  'class',
            'type':       type_name,
            'attributes': ['intent(inout)'],
            'variables':  ['self'],
        },
        {
            'intrinsic':  'integer',
            'attributes': ['intent(inout)'],
            'variables':  ['integerProperty', 'integerBufferCount'],
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
            'variables':  ['doubleProperty', 'doubleBufferCount'],
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
            'intrinsic':  'type',
            'type':       'multiCounter',
            'attributes': ['intent(in   )'],
            'variables':  ['outputInstance'],
        },
        {
            'intrinsic':  'integer',
            'attributes': ['intent(in   )'],
            'variables':  ['instance'],
        },
    ]

    function = {
        'type':        'void',
        'name':        name + 'Output',
        'description': (
            f"Populate output buffers with properties to output for a "
            f"\\mono{{{name}}} component."
        ),
        'content':     '',
        'modules':     [
            'Multi_Counters',
            'Merger_Tree_Outputter_Buffer_Types',
        ],
        'variables':   list(base_variables),
    }

    # Discover output-derived types, modules, and arg-usage.
    argument_usage = {k: False for k in
                      ('self', 'time', 'instance', 'integer', 'double')}
    output_derived_types = set()
    modules_required     = set()

    for member in class_dict.get('members') or []:
        if (
            isinstance(member.get('output'), dict)
            and member['output'].get('instances') == 'first'
        ):
            argument_usage['instance'] = True
        for prop in _component_properties(member):
            if 'output' not in prop:
                continue
            ptype = (prop.get('data') or {}).get('type')
            argument_usage['self'] = True
            if not is_output_intrinsic(ptype):
                argument_usage['time'] = True
                argument_usage['integer'] = True
                argument_usage['double']  = True
                output_derived_types.add(ptype)
            else:
                argument_usage[
                    _INTRINSIC_TYPE_MAP[ptype]['label']
                ] = True
            output_block = prop['output'] if isinstance(prop['output'], dict) else {}
            if 'modules' in output_block:
                for m in str(output_block['modules']).split(','):
                    modules_required.add(m)

    # Per-derived-type output variable.
    for type_label in sorted(output_derived_types):
        function['variables'].append({
            'intrinsic': 'type',
            'type':      type_label,
            'variables': [f'output{_ucfirst(type_label)}'],
        })

    # Required modules — appended sorted.
    for m in sorted(modules_required):
        function['modules'].append(m)

    # Unused argument list — mirrors the nestedmap at Output.pm:344-359.
    arguments_unused = []
    for key in sorted(argument_usage.keys()):
        if argument_usage[key]:
            continue
        if key in ('integer', 'double'):
            for suffix in ('Property', 'BufferCount', 'Properties'):
                arguments_unused.append(key + suffix)
        else:
            arguments_unused.append(key)
    arguments_unused.append('outputInstance')

    content = ''
    if arguments_unused:
        content += (
            "!$GLC attributes unused :: " + ", ".join(arguments_unused) + "\n"
        )

    # Per-member body.
    tmps_added = {'double': False, 'integer': False, 'index': False}

    for member in class_dict.get('members') or []:
        # Skip members with no outputs.
        if not any('output' in p for p in _component_properties(member)):
            continue
        guard_extra = ''
        m_output = member.get('output') if isinstance(member.get('output'), dict) else {}
        if m_output.get('instances') == 'first':
            guard_extra = " .and. instance == 1"
        content += (
            f"if (default{cap}Component%{member['name']}IsActive() {guard_extra}) then\n"
        )

        for prop in _component_properties(member):
            if 'output' not in prop or prop.get('definedInParent'):
                continue
            data = prop.get('data') or {}
            ptype = data.get('type')
            rank  = int(data.get('rank') or 0)

            if rank == 0:
                count = 1
            elif rank == 1:
                output = prop['output'] if isinstance(prop['output'], dict) else {}
                labels = output.get('labels')
                m = re.match(r'^\[(.*)\]$', labels or '')
                if m:
                    count = m.group(1).count(',') + 1
                else:
                    count = output.get('count')
            else:
                count = None

            buffer_type = _INTRINSIC_TYPE_MAP[ptype]['label'] \
                if ptype in _INTRINSIC_TYPE_MAP else None

            if is_output_intrinsic(ptype):
                if rank == 0:
                    content += (
                        f"{buffer_type}Property={buffer_type}Property+1\n"
                        f"{buffer_type}Properties({buffer_type}Property)"
                        f"%scalar({buffer_type}BufferCount)=self%{prop['name']}()\n"
                    )
                else:
                    if buffer_type == 'integer' and not tmps_added['integer']:
                        tmps_added['integer'] = True
                        function['variables'].append({
                            'intrinsic':  'integer',
                            'type':       'kind_int8',
                            'attributes': ['allocatable', 'dimension(:)'],
                            'variables':  ['integerOutputTmp'],
                        })
                    if buffer_type == 'double' and not tmps_added['double']:
                        tmps_added['double'] = True
                        function['variables'].append({
                            'intrinsic':  'double precision',
                            'attributes': ['allocatable', 'dimension(:)'],
                            'variables':  ['doubleOutputTmp'],
                        })
                    if not tmps_added['index']:
                        tmps_added['index'] = True
                        function['variables'].append({
                            'intrinsic': 'integer',
                            'variables': ['i'],
                        })
                    content += (
                        f"{buffer_type}OutputTmp=reshape(self%{prop['name']}(),"
                        f"[{count}])\n"
                        f"do i=1,{count}\n"
                        f"  {buffer_type}Properties({buffer_type}Property+i)"
                        f"%scalar({buffer_type}BufferCount)={buffer_type}OutputTmp(i)\n"
                        f"end do\n"
                        f"deallocate({buffer_type}OutputTmp)\n"
                        f"{buffer_type}Property={buffer_type}Property+{count}\n"
                    )
            else:
                cap_type = _ucfirst(ptype)
                content += (
                    f"output{cap_type}=self%{prop['name']}()\n"
                    f"call output{cap_type}%output(integerProperty,"
                    "integerBufferCount,integerProperties,doubleProperty,"
                    "doubleBufferCount,doubleProperties,time,outputInstance)\n"
                    f"if (.not.same_type_as(self,{name}Class)) "
                    f"call self%{prop['name']}Set(output{cap_type})\n"
                )

        content += "end if\n"

    function['content'] = content
    _bind(build, type_name, function, 'output')


def Class_Post_Output(build, class_dict):
    """Generate `<class>PostOutput` — no-op stub.
    Mirrors `Class_Post_Output`.

    Note this hook does NOT skip inactive classes (matches Perl).
    """
    name      = class_dict['name']
    cap       = _ucfirst(name)
    type_name = 'nodeComponent' + cap
    function = {
        'type':        'void',
        'name':        name + 'PostOutput',
        'description': (
            f"Perform post-output processing of a \\mono{{{name}}} component."
        ),
        'content':     "!$GLC attributes unused :: self, time\n",
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       type_name,
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
    _bind(build, type_name, function, 'postOutput')


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
# Hook registration.  Order matches Perl Classes/Output.pm:23-28.
# ---------------------------------------------------------------------------

register('classesOutput', 'classIteratedFunctions', Class_Dump_ASCII)
register('classesOutput', 'classIteratedFunctions', Class_Output_Count)
register('classesOutput', 'classIteratedFunctions', Class_Output_Names)
register('classesOutput', 'classIteratedFunctions', Class_Post_Output)
register('classesOutput', 'classIteratedFunctions', Class_Output)
