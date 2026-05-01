# Per-property evolution methods: count / rateGet / Rate / scale /
# inactive / analytic + auto-create rate functions on the class type.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Properties/Evolve.pm.

import os
import sys


from Galacticus.Build.Components.Utils      import (
    register,
    is_intrinsic,
    offset_name,
)
from Galacticus.Build.Components.DataTypes  import data_object_definition


def _is_debugging():
    flags = os.environ.get('GALACTICUS_FCFLAGS', '')
    return any(tok == '-DDEBUGGING' for tok in flags.split())


def _evolvable(prop):
    """Return True if this property is non-virtual and evolvable."""
    attrs = prop.get('attributes') or {}
    return (
        not attrs.get('isVirtual')
        and bool(attrs.get('isEvolvable'))
    )


def _impl_type(class_dict, member):
    return (
        'nodeComponent'
        + _ucfirst(class_dict['name'])
        + _ucfirst(member['name'])
    )


def Build_Count_Functions(build, class_dict, member, prop):
    """Generate `<class><Member><Prop>Count`.  Mirrors `Build_Count_Functions`."""
    if not _evolvable(prop):
        return
    impl_type = _impl_type(class_dict, member)
    type_prefix = class_dict['name'] + _ucfirst(member['name'])
    cap_prop  = _ucfirst(prop['name'])
    fn_name   = type_prefix + cap_prop + 'Count'

    function = {
        'type':        'integer',
        'name':        fn_name,
        'description': (
            f"Return a count of the number of scalar properties in the "
            f"\\mono{{{prop['name']}}} property of an "
            f"\\mono{{{member['name']}}} implementation of the "
            f"\\mono{{{class_dict['name']}}} component class."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
        ],
    }
    rank = int((prop.get('data') or {}).get('rank') or 0)
    if rank == 0:
        function['content'] = (
            "!$GLC attributes unused :: self\n"
            f"{type_prefix}{cap_prop}Count=1\n"
        )
    elif rank == 1:
        function['content'] = (
            f"if (allocated(self%{prop['name']}Data)) then\n"
            f"   {type_prefix}{cap_prop}Count=size(self%{prop['name']}Data)\n"
            "else\n"
            f"   {type_prefix}{cap_prop}Count=0\n"
            "end if\n"
        )
    else:
        return  # rank > 1 not supported

    _bind(build, impl_type, function, prop['name'] + 'Count')


def Build_Rate_Get_Functions(build, class_dict, member, prop):
    """Generate `<class><Member><Prop>RateGet`.  Mirrors `Build_Rate_Get_Functions`."""
    if not _evolvable(prop):
        return

    type_descriptor, _ = data_object_definition(prop.get('data') or {})
    function_type = type_descriptor['intrinsic']
    if 'type' in type_descriptor:
        function_type += f"({type_descriptor['type']})"
    if type_descriptor.get('attributes'):
        function_type += ', ' + ', '.join(type_descriptor['attributes'])

    impl_type = _impl_type(class_dict, member)
    data  = prop.get('data') or {}
    ptype = data.get('type')
    rank  = int(data.get('rank') or 0)

    function = {
        'type':        function_type + ' => propertyRate',
        'name':        (
            class_dict['name']
            + _ucfirst(member['name'])
            + _ucfirst(prop['name'])
            + 'RateGet'
        ),
        'description': (
            f"Get the rate of change of the \\mono{{{prop['name']}}} "
            f"property of an \\mono{{{member['name']}}} implementation "
            f"of the \\mono{{{class_dict['name']}}} component class."
        ),
        'modules':     ['Error'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['offset'],
            },
        ],
    }
    if (is_intrinsic(ptype) and rank > 0) or not is_intrinsic(ptype):
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': ['count'],
        })

    unused = []
    if is_intrinsic(ptype):
        unused.append('self')

    cls_member = class_dict['name'] + _ucfirst(member['name'])
    offset_all      = offset_name('all',      cls_member, prop['name'])
    offset_active   = offset_name('active',   cls_member, prop['name'])

    offset_block = (
        f"if (nodeAnalytics({offset_all})) call Error_Report("
        "'rates are gettable only for numerically-solved variables'"
        "//{introspection:location})\n"
        "if (rateComputeState == propertyTypeInactive) then\n"
        f" if (nodeInactives({offset_all})) call Error_Report("
        "'rates are gettable only for active variables'"
        "//{introspection:location})\n"
        f" offset={offset_active}\n"
        "else\n"
        " offset=0\n"
        " call Error_Report("
        "'rates are gettable only during inactive variable integration'"
        "//{introspection:location})\n"
        "end if\n"
    )

    content = ''
    if unused:
        content += "!$GLC attributes unused :: " + ",".join(unused) + "\n"

    if is_intrinsic(ptype):
        if rank == 0:
            content += offset_block
            content += "propertyRate=nodeRatesActives(offset)\n"
        else:
            content += offset_block
            content += (
                f"count=size(self%{prop['name']}Data)\n"
                "propertyRate=nodeRatesActives(offset:offset+count-1)\n"
            )
    else:
        content += (
            f"count=self%{prop['name']}Data%serializeCount()\n"
            "if (count > 0) then\n"
            f"{offset_block}"
            f"   propertyRate=self%{prop['name']}Data\n"
            "   call propertyRate%deserialize"
            "(nodeRatesActives(offset:offset+count-1))\n"
            "end if\n"
        )

    function['content'] = content
    _bind(build, impl_type, function, prop['name'] + 'RateGet')


def Build_Rate_Functions(build, class_dict, member, prop):
    """Generate `<class><Member><Prop>Rate(Intrinsic)?`.  Mirrors `Build_Rate_Functions`."""
    if not _evolvable(prop):
        return

    debugging = _is_debugging()
    deferred  = (prop.get('attributes') or {}).get('isDeferred')
    intrinsic_rate = bool(deferred) and 'rate' in str(deferred).split(':')
    suffix = 'Intrinsic' if intrinsic_rate else ''

    type_descriptor, _ = data_object_definition(prop.get('data') or {}, match_only=True)
    type_descriptor.setdefault('variables',  []).append('setValue')
    type_descriptor.setdefault('attributes', []).append('intent(in   )')

    impl_type = _impl_type(class_dict, member)
    data  = prop.get('data') or {}
    ptype = data.get('type')
    rank  = int(data.get('rank') or 0)

    function = {
        'type':        'void',
        'name':        (
            class_dict['name']
            + _ucfirst(member['name'])
            + _ucfirst(prop['name'])
            + 'Rate'
            + suffix
        ),
        'description': (
            "Accumulate"
            + (" directly (i.e. circumventing any deferred function binding)" if intrinsic_rate else "")
            + f" to the rate of change of the \\mono{{{prop['name']}}} "
            f"property of an \\mono{{{member['name']}}} implementation of "
            f"the \\mono{{{class_dict['name']}}} component class."
        ),
        'modules':     ['Error'],
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
                'attributes': ['optional', 'intent(inout)'],
                'variables':  ['interrupt'],
            },
            {
                'intrinsic':  'procedure',
                'type':       'interruptTask',
                'attributes': ['optional', 'intent(inout)', 'pointer'],
                'variables':  ['interruptProcedure'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['offset'],
            },
        ],
    }
    if debugging:
        function['modules'].extend(['Debugging', 'ISO_Varying_String'])
        function['variables'].extend([
            {
                'intrinsic': 'type',
                'type':      'varying_string',
                'variables': ['message'],
            },
            {
                'intrinsic': 'character',
                'type':      'len=32',
                'variables': ['label'],
            },
        ])

    debug_iterator_required = False
    if not is_intrinsic(ptype):
        current_descriptor, _ = data_object_definition(prop.get('data') or {}, match_only=True)
        current_descriptor.setdefault('variables', []).append('current')
        function['variables'].append(current_descriptor)
    if (is_intrinsic(ptype) and rank > 0) or not is_intrinsic(ptype):
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': ['count'],
        })

    unused = ['interrupt', 'interruptProcedure']
    if is_intrinsic(ptype):
        unused.append('self')

    cls_member = class_dict['name'] + _ucfirst(member['name'])
    offset_all      = offset_name('all',      cls_member, prop['name'])
    offset_active   = offset_name('active',   cls_member, prop['name'])
    offset_inactive = offset_name('inactive', cls_member, prop['name'])

    offset_block = (
        f"if (nodeAnalytics({offset_all})) call Error_Report("
        "'rates are only settable for numerically-solved variables'"
        "//{introspection:location})\n"
        "if (rateComputeState == propertyTypeAll          ) then\n"
        f" offset={offset_all}\n"
        "else if (rateComputeState == propertyTypeActive  ) then\n"
        f" if (     nodeInactives({offset_all})) return\n"
        f" offset={offset_active}\n"
        "else if (rateComputeState == propertyTypeInactive) then\n"
        f" if (.not.nodeInactives({offset_all})) return\n"
        f" offset={offset_inactive}\n"
        "else if (rateComputeState == propertyTypeNumerics) then\n"
        f" offset={offset_active}\n"
        "else\n"
        " return\n"
        "end if\n"
    )

    content = "!$GLC attributes unused :: " + ",".join(unused) + "\n"

    if is_intrinsic(ptype):
        if rank == 0:
            content += offset_block
            content += "nodeRates(offset)=nodeRates(offset)+setValue\n"
            if debugging:
                content += (
                    "if (isDebugging()) then\n"
                    " write (label,'(e12.6)') setValue\n"
                    f" message=\"   rate: (\"//getCaller()//\") "
                    f"{class_dict['name']}:{member['name']}:{prop['name']} "
                    "\"//trim(adjustl(label))\n"
                    " call debugLog(message)\n"
                    "end if\n"
                )
        else:
            content += offset_block
            content += (
                "count=size(setValue)\n"
                "nodeRates(offset:offset+count-1)="
                "nodeRates(offset:offset+count-1)+setValue\n"
            )
            if debugging:
                debug_iterator_required = True
                content += (
                    "if (isDebugging()) then\n"
                    " do i=1,count\n"
                    "  write (label,'(a1,i4.4,a2,e12.6)') \"[\",i,\"] \",setValue(i)\n"
                    f"  message=\"   rate: (\"//getCaller()//\") "
                    f"{class_dict['name']}:{member['name']}:{prop['name']}\""
                    "//trim(adjustl(label))\n"
                    "  call debugLog(message)\n"
                    " end do\n"
                    "end if\n"
                )
    else:
        content += (
            f"count=self%{prop['name']}Data%serializeCount()\n"
            "if (count > 0) then\n"
            f"{offset_block}"
            f"   current=self%{prop['name']}Data\n"
            "   call current%deserialize(nodeRates(offset:offset+count-1))\n"
            "   call current%increment(setValue)\n"
            "   call current%serialize(nodeRates(offset:offset+count-1))\n"
            "end if\n"
        )
        if debugging:
            function['variables'].append({
                'intrinsic':  'double precision',
                'attributes': ['allocatable', 'dimension(:)'],
                'variables':  ['rates'],
            })
            debug_iterator_required = True
            content += (
                "if (isDebugging() .and. count > 0) then\n"
                " allocate(rates(count))\n"
                " call setValue%serialize(rates)\n"
                " do i=1,count\n"
                "  write (label,'(a1,i4.4,a2,e12.6)') \"[\",i,\"] \",rates(i)\n"
                f"  message=\"   rate: (\"//getCaller()//\") "
                f"{class_dict['name']}:{member['name']}:{prop['name']}\""
                "//trim(adjustl(label))\n"
                "  call debugLog(message)\n"
                " end do\n"
                "end if\n"
            )

    if debug_iterator_required:
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': ['i'],
        })

    function['content'] = content
    _bind(build, impl_type, function, prop['name'] + 'Rate' + suffix)


def Build_Auto_Create_Rate_Functions(build, class_dict, member, prop):
    """Mirrors `Build_Auto_Create_Rate_Functions`.  Bound to the
    component class type, not the implementation.
    """
    attrs = prop.get('attributes') or {}
    if (
        attrs.get('isVirtual')
        or not attrs.get('isEvolvable')
        or not attrs.get('createIfNeeded')
    ):
        return

    cap_class = _ucfirst(class_dict['name'])
    type_name = 'nodeComponent' + cap_class
    method_name = prop['name'] + 'Rate'
    bound = build.setdefault('types', {}).setdefault(type_name, {}) \
                                          .setdefault('boundFunctions', [])
    if any(b.get('name') == method_name for b in bound):
        return

    type_descriptor, _ = data_object_definition(prop.get('data') or {}, match_only=True)
    type_descriptor.setdefault('variables',  []).append('setValue')
    type_descriptor.setdefault('attributes', []).append('intent(in   )')

    function = {
        'type':        'void',
        'name':        class_dict['name'] + _ucfirst(prop['name']) + 'Rate',
        'description': (
            f"Accept a rate set for the \\mono{{{prop['name']}}} property "
            f"of the \\mono{{{class_dict['name']}}} component class. "
            "Trigger an interrupt to create the component."
        ),
        'modules':     ['Error'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       type_name,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            type_descriptor,
            {
                'intrinsic':  'logical',
                'attributes': ['optional', 'intent(inout)'],
                'variables':  ['interrupt'],
            },
            {
                'intrinsic':  'procedure',
                'type':       'interruptTask',
                'attributes': ['optional', 'intent(inout)', 'pointer'],
                'variables':  ['interruptProcedure'],
            },
        ],
    }

    data  = prop.get('data') or {}
    ptype = data.get('type')
    rank  = int(data.get('rank') or 0)

    content = (
        "!$GLC attributes unused :: self\n"
        "! No specific component exists, so we must interrupt and create "
        "one unless the rate is zero.\n"
    )
    if rank == 0:
        if is_intrinsic(ptype):
            content += "if (setValue == 0.0d0) return\n"
        else:
            content += "if (setValue%isZero()) return\n"
    else:
        if is_intrinsic(ptype):
            content += "if (all(setValue == 0.0d0)) return\n"
    content += (
        "if (.not.(present(interrupt).and.present(interruptProcedure))) "
        "call Error_Report('interrupt required, but optional arguments "
        "missing'//{introspection:location})\n"
        "interrupt=.true.\n"
        f"interruptProcedure => {class_dict['name']}CreateByInterrupt\n"
    )
    function['content'] = content
    bound.append({
        'type':       'procedure',
        'descriptor': function,
        'name':       method_name,
    })


def Build_Scale_Functions(build, class_dict, member, prop):
    """Mirrors `Build_Scale_Functions`."""
    if not _evolvable(prop):
        return
    type_descriptor, _ = data_object_definition(prop.get('data') or {}, match_only=True)
    type_descriptor.setdefault('variables',  []).append('setValue')
    type_descriptor.setdefault('attributes', []).append('intent(in   )')

    impl_type = _impl_type(class_dict, member)
    data  = prop.get('data') or {}
    ptype = data.get('type')
    rank  = int(data.get('rank') or 0)

    function = {
        'type':        'void',
        'name':        (
            class_dict['name']
            + _ucfirst(member['name'])
            + _ucfirst(prop['name'])
            + 'Scale'
        ),
        'description': (
            f"Set the absolute scale of the \\mono{{{prop['name']}}} "
            f"property of an \\mono{{{member['name']}}} implementation of "
            f"the \\mono{{{class_dict['name']}}} component class."
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
    }
    if not is_intrinsic(ptype):
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': ['count'],
        })

    cls_member = class_dict['name'] + _ucfirst(member['name'])
    offset = offset_name('all', cls_member, prop['name'])
    content = "!$GLC attributes unused :: self\n"
    if is_intrinsic(ptype):
        if rank == 0:
            content += f"nodeScales({offset})=setValue\n"
        else:
            content += f"nodeScales({offset}:{offset}+size(setValue))=setValue\n"
    else:
        content += (
            "count=setValue%serializeCount()\n"
            f"if (count > 0) call setValue%serialize(nodeScales({offset}:{offset}+count-1))\n"
        )
    function['content'] = content
    _bind(build, impl_type, function, prop['name'] + 'Scale')


def Build_Inactive_Functions(build, class_dict, member, prop):
    """Mirrors `Build_Inactive_Functions`."""
    if not _evolvable(prop):
        return
    impl_type = _impl_type(class_dict, member)
    data  = prop.get('data') or {}
    ptype = data.get('type')
    rank  = int(data.get('rank') or 0)

    function = {
        'type':        'void',
        'name':        (
            class_dict['name']
            + _ucfirst(member['name'])
            + _ucfirst(prop['name'])
            + 'JcbnZr'
        ),
        'description': (
            f"Indicate that the \\mono{{{prop['name']}}} property of an "
            f"\\mono{{{member['name']}}} implementation of the "
            f"\\mono{{{class_dict['name']}}} component class is inactive "
            "for differential equation solving."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
        ],
    }
    if not is_intrinsic(ptype):
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': ['count'],
        })

    cls_member = class_dict['name'] + _ucfirst(member['name'])
    offset = offset_name('all', cls_member, prop['name'])
    content = "!$GLC attributes unused :: self\n"
    if is_intrinsic(ptype):
        if rank == 0:
            content += f"nodeInactives({offset})=.true.\n"
        else:
            content += (
                f"nodeInactives({offset}:{offset}+size(self%{prop['name']}Data)-1)=.true.\n"
            )
    else:
        content += (
            f"count=self%{prop['name']}Data%serializeCount()\n"
            f"if (count > 0) nodeInactives({offset}:{offset}+count-1)=.true.\n"
        )
    function['content'] = content
    _bind(build, impl_type, function, prop['name'] + 'Inactive')


def Build_Analytic_Functions(build, class_dict, member, prop):
    """Mirrors `Build_Analytic_Functions`."""
    if not _evolvable(prop):
        return
    impl_type = _impl_type(class_dict, member)
    data  = prop.get('data') or {}
    ptype = data.get('type')
    rank  = int(data.get('rank') or 0)

    function = {
        'type':        'void',
        'name':        (
            class_dict['name']
            + _ucfirst(member['name'])
            + _ucfirst(prop['name'])
            + 'Alytc'
        ),
        'description': (
            f"Indicate that the \\mono{{{prop['name']}}} property of an "
            f"\\mono{{{member['name']}}} implementation of the "
            f"\\mono{{{class_dict['name']}}} component class is to be "
            "solved analytically during differential evolution."
        ),
        'modules':     ['Error'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
        ],
    }
    if not is_intrinsic(ptype):
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': ['count'],
        })

    cls_member = class_dict['name'] + _ucfirst(member['name'])
    offset = offset_name('all', cls_member, prop['name'])
    content = "!$GLC attributes unused :: self\n"
    if is_intrinsic(ptype):
        if rank == 0:
            content += (
                f"if (nodeAnalytics({offset})) call Error_Report("
                "'property is already marked analytically-solvable'"
                "//{introspection:location})\n"
                f"nodeAnalytics({offset})=.true.\n"
            )
        else:
            content += (
                f"if (any(nodeAnalytics({offset}:{offset}+size(self%{prop['name']}Data)-1))) "
                "call Error_Report('property is already marked analytically-solvable'"
                "//{introspection:location})\n"
                f"nodeAnalytics({offset}:{offset}+size(self%{prop['name']}Data)-1)=.true.\n"
            )
    else:
        content += (
            f"count=self%{prop['name']}Data%serializeCount()\n"
            "if (count > 0) then\n"
            f" if (any(nodeAnalytics({offset}:{offset}+count-1))) "
            "call Error_Report('property is already marked analytically-solvable'"
            "//{introspection:location})\n"
            f" nodeAnalytics({offset}:{offset}+count-1)=.true.\n"
            "end if\n"
        )
    function['content'] = content
    _bind(build, impl_type, function, prop['name'] + 'Analytic')


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
# Hook registration.  Order matches Perl Properties/Evolve.pm:21-29.
# ---------------------------------------------------------------------------

register('propertiesEvolve', 'propertyIteratedFunctions', Build_Count_Functions)
register('propertiesEvolve', 'propertyIteratedFunctions', Build_Rate_Get_Functions)
register('propertiesEvolve', 'propertyIteratedFunctions', Build_Rate_Functions)
register('propertiesEvolve', 'propertyIteratedFunctions', Build_Auto_Create_Rate_Functions)
register('propertiesEvolve', 'propertyIteratedFunctions', Build_Scale_Functions)
register('propertiesEvolve', 'propertyIteratedFunctions', Build_Inactive_Functions)
register('propertiesEvolve', 'propertyIteratedFunctions', Build_Analytic_Functions)
