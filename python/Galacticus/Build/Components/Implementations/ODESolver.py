"""Per-implementation ODE-solver hooks.

Andrew Benson (ported to Python 2026)

Eight hooks total — seven implementationIteratedFunctions plus one
top-level `functions` hook (Implementation_ODE_Rate_Variables).
"""



from Galacticus.Build.Components.Utils                  import (
    register,
    offset_name,
    _component_properties,
)
from Galacticus.Build.Components.Implementations.Utils  import (
    has_real_evolvers,
    has_real_non_trivial_evolvers,
    list_real_evolvers,
)


def _impl_type(class_dict, member):
    return (
        'nodeComponent'
        + _ucfirst(class_dict['name'])
        + _ucfirst(member['name'])
    )


def _evolver_count_expr(prop):
    """Return the Fortran expression for the serialization count of a
    single property — `size(self%<X>Data)` for rank>0, `1` for rank-0
    doubles, `self%<X>Data%serializeCount()` otherwise.

    Mirrors the count-builder ternary repeated across most of these
    hooks.
    """
    data  = prop.get('data') or {}
    rank  = int(data.get('rank') or 0)
    if rank > 0:
        return f"size(self%{prop['name']}Data)"
    if data.get('type') == 'double':
        return "1"
    return f"self%{prop['name']}Data%serializeCount()"


def Implementation_ODE_Name_From_Index(build, class_dict, member):
    """Generate `<impl>NameFromIndex`.
    """
    impl_type = _impl_type(class_dict, member)
    has_extends = isinstance(member.get('extends'), dict)

    function = {
        'type':        'type(varying_string) => name',
        'name':        impl_type + 'NameFromIndex',
        'description': (
            f"Return the name of the property of given index for a "
            f"\\mono{{{member['name']}}} implementation of the "
            f"\\mono{{{class_dict['name']}}} component class."
        ),
        'modules':     ['ISO_Varying_String'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(inout)'],
                'variables':  ['count'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['propertyType'],
            },
        ],
    }

    unused = []
    if not (has_extends or has_real_non_trivial_evolvers(member)):
        unused.append('self')
    # A second append to `unused` is deliberately never emitted here: the
    # historical guard for it was always true, so the append was never
    # reached.  That behavior is preserved (by simply omitting the append)
    # to keep the generated code identical.

    if not has_extends:
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': ['i'],
        })

    content = ''
    if unused:
        content += "!$GLC attributes unused :: " + ",".join(unused) + "\n"

    if has_extends:
        ext = member['extends']
        parent_type = (
            'nodeComponent'
            + _ucfirst(ext.get('class', ''))
            + _ucfirst(ext.get('name', ''))
        )
        content += (
            f"name=self%{parent_type}%nameFromIndex(count,propertyType)\n"
            "if (count <= 0) return\n"
        )
    elif class_dict['name'] in (build.get('componentClassListActive') or []):
        offset = offset_name('all', class_dict['name'], 'floatRank0MetaProperties')
        content += (
            f"if (allocated({class_dict['name']}FloatRank0MetaPropertyNames)) then\n"
            f" do i=1,size({class_dict['name']}FloatRank0MetaPropertyNames)\n"
            "   if (                                                                                    &\n"
            f"    &  {class_dict['name']}FloatRank0MetaPropertyEvolvable(i)                                &\n"
            "    & .and.                                                                                &\n"
            f"    &                                                 .not.nodeAnalytics({offset}(i)) &\n"
            "    & .and.                                                                                &\n"
            "    &  (                                                                                   &\n"
            "    &     propertyType == propertyTypeAll                                                  &\n"
            "    &   .or.                                                                               &\n"
            f"    &    (propertyType == propertyTypeActive   .and. .not.nodeInactives({offset}(i))) &\n"
            "    &   .or.                                                                               &\n"
            f"    &    (propertyType == propertyTypeInactive .and.      nodeInactives({offset}(i))) &\n"
            "    &   .or.                                                                               &\n"
            "    &     propertyType == propertyTypeNumerics                                             &\n"
            "    &  )                                                                                   &\n"
            "    & ) count=count-1\n"
            "  if (count <= 0) then\n"
            f"   name={class_dict['name']}FloatRank0MetaPropertyNames(i)\n"
            "   return\n"
            "  end if\n"
            " end do\n"
            "end if\n"
        )

    for prop in list_real_evolvers(member):
        data  = prop.get('data') or {}
        rank  = int(data.get('rank') or 0)
        count_expr = _evolver_count_expr(prop)
        if rank > 0:
            cond_open  = f"if (allocated(self%{prop['name']}Data)) then"
            cond_close = "end if"
        else:
            cond_open  = ""
            cond_close = ""
        offset = offset_name('all', class_dict, member, prop)
        content += (
            f"{cond_open}\n"
            " if (                                                                                                               &\n"
            f"  &                                                all(.not.nodeAnalytics({offset}:{offset}+{count_expr}-1))  &\n"
            "  &  .and.                                                                                                          &\n"
            "  &  (                                                                                                              &\n"
            "  &     propertyType == propertyTypeAll                                                                             &\n"
            "  &   .or.                                                                                                          &\n"
            f"  &    (propertyType == propertyTypeActive   .and. all(.not.nodeInactives({offset}:{offset}+{count_expr}-1))) &\n"
            "  &   .or.                                                                                                          &\n"
            f"  &    (propertyType == propertyTypeInactive .and. all(     nodeInactives({offset}:{offset}+{count_expr}-1))) &\n"
            "  &   .or.                                                                                                          &\n"
            "  &     propertyType == propertyTypeNumerics                                                                        &\n"
            "  &  )                                                                                                              &\n"
            f"  & ) count=count-{count_expr}\n"
            f"{cond_close}\n"
            "if (count <= 0) then\n"
            f"  name='{class_dict['name']}:{member['name']}:{prop['name']}'\n"
            "  return\n"
            "end if\n"
        )

    content += "name='?'\n"
    function['content'] = content

    _bind(build, impl_type, function, 'nameFromIndex')


def Implementation_ODE_Serialize_Count(build, class_dict, member):
    """Generate `<impl>SerializeCount`.
    """
    impl_type   = _impl_type(class_dict, member)
    has_extends = isinstance(member.get('extends'), dict)
    has_evolve  = has_real_evolvers(member)
    has_nontriv = has_real_non_trivial_evolvers(member)

    function = {
        'type':        'integer',
        'name':        impl_type + 'SerializeCount',
        'description': (
            f"Return a count of the serialization of a {member['name']} "
            f"implementation of the {class_dict['name']} component."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['propertyType'],
            },
        ],
    }
    if has_evolve or not has_extends:
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': ['count'],
        })
    if not has_extends:
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': ['i'],
        })

    unused = []
    if not (has_extends or has_nontriv):
        unused.append('self')
    if not has_evolve:
        unused.append('propertyType')

    content = ''
    if unused:
        content += "!$GLC attributes unused :: " + ",".join(unused) + "\n"

    if has_extends:
        ext = member['extends']
        parent_type = (
            'nodeComponent'
            + _ucfirst(ext.get('class', ''))
            + _ucfirst(ext.get('name', ''))
        )
        content += (
            f"{impl_type}SerializeCount=self%{parent_type}%serializeCount(propertyType)\n"
        )
    elif class_dict['name'] in (build.get('componentClassListActive') or []):
        cn = class_dict['name']
        content += (
            f"{impl_type}SerializeCount=0\n"
            f"if (allocated({cn}FloatRank0MetaPropertyNames)) then\n"
            " count=0\n"
            f" do i=1,size({cn}FloatRank0MetaPropertyNames)\n"
            f"  if (.not.{cn}FloatRank0MetaPropertyEvolvable(i)) cycle\n"
            "  if      (propertyType == propertyTypeAll     ) then\n"
            f"                                                                                                                                                                  {impl_type}SerializeCount={impl_type}SerializeCount+1\n"
            "  else if (propertyType == propertyTypeInactive) then\n"
            f"   if (     nodeInactives(offsetAll{cn}FloatRank0MetaProperties(i)).and..not.nodeAnalytics(offsetAll{cn}FloatRank0MetaProperties(i))) {impl_type}SerializeCount={impl_type}SerializeCount+1\n"
            "  else if (propertyType == propertyTypeActive  ) then\n"
            f"   if (.not.nodeInactives(offsetAll{cn}FloatRank0MetaProperties(i)).and..not.nodeAnalytics(offsetAll{cn}FloatRank0MetaProperties(i))) {impl_type}SerializeCount={impl_type}SerializeCount+1\n"
            "  else if (propertyType == propertyTypeNumerics) then\n"
            f"   if (                                                                               .not.nodeAnalytics(offsetAll{cn}FloatRank0MetaProperties(i))) {impl_type}SerializeCount={impl_type}SerializeCount+1\n"
            "  end if\n"
            "  count=count+1\n"
            " end do\n"
            "end if\n"
        )

    for prop in list_real_evolvers(member):
        data  = prop.get('data') or {}
        rank  = int(data.get('rank') or 0)
        count_expr = _evolver_count_expr(prop)
        offset = offset_name('all', class_dict, member, prop)
        content += f"count={count_expr}\n"
        if rank > 0:
            cond_open  = f"if (allocated(self%{prop['name']}Data)) then"
            cond_close = "end if"
        else:
            cond_open  = ""
            cond_close = ""
        content += (
            f"{cond_open}\n"
            "if (propertyType == propertyTypeAll) then\n"
            f" {impl_type}SerializeCount={impl_type}SerializeCount+count\n"
            f"else if (((propertyType == propertyTypeActive .and. all(.not.nodeInactives({offset}:{offset}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({offset}:{offset}+count-1))) .or. propertyType == propertyTypeNumerics) .and. all(.not.nodeAnalytics({offset}:{offset}+count-1))) then\n"
            f" {impl_type}SerializeCount={impl_type}SerializeCount+count\n"
            "end if\n"
            f"{cond_close}\n"
        )

    function['content'] = content
    _bind(build, impl_type, function, 'serializeCount')


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


def Implementation_ODE_Serialize_Values(build, class_dict, member):
    """Generate `<impl>SerializeValues`.
    """
    impl_type   = _impl_type(class_dict, member)
    has_extends = isinstance(member.get('extends'), dict)
    has_evolve  = has_real_evolvers(member)
    has_nontriv = has_real_non_trivial_evolvers(member)

    function = {
        'type':        'void',
        'name':        impl_type + 'SerializeValues',
        'description': (
            f"Serialize evolvable properties of a "
            f"\\mono{{{member['name']}}} implementation of the "
            f"\\mono{{{class_dict['name']}}} component to array."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['intent(  out)', 'dimension(:)'],
                'variables':  ['array'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['propertyType'],
            },
        ],
    }

    required = ['offset']
    if has_extends or has_nontriv:
        required.append('count')
    if not has_extends:
        required.append('i')
    if required:
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': list(required),
        })

    unused = [] if (has_extends or has_evolve) else ['array', 'propertyType']

    content = ''
    if unused:
        content += "!$GLC attributes unused :: " + ",".join(unused) + "\n"
    if 'offset' in required:
        content += "offset=1\n"

    if has_extends:
        ext = member['extends']
        parent_type = (
            'nodeComponent'
            + _ucfirst(ext.get('class', ''))
            + _ucfirst(ext.get('name', ''))
        )
        content += (
            f"count=self%{parent_type}%serializeCount (      propertyType)\n"
            "if (count > 0) then\n"
            f" call self%{parent_type}%serializeValues(array,propertyType)\n"
            " offset=offset+count\n"
            "end if\n"
        )
    elif class_dict['name'] in (build.get('componentClassListActive') or []):
        offset = offset_name('all', class_dict['name'], 'floatRank0MetaProperties')
        cn = class_dict['name']
        content += (
            f"if (allocated({cn}FloatRank0MetaPropertyNames)) then\n"
            f" do i=1,size({cn}FloatRank0MetaPropertyNames)\n"
            f"  if ({cn}FloatRank0MetaPropertyEvolvable(i).and..not.nodeAnalytics({offset}(i)).and.(propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. .not.nodeInactives({offset}(i))) .or. (propertyType == propertyTypeInactive .and. nodeInactives({offset}(i))) .or. propertyType == propertyTypeNumerics)) then\n"
            "   array(offset)=self%floatRank0MetaProperties(i)\n"
            "   offset=offset+1\n"
            "  end if\n"
            " end do\n"
            "end if\n"
        )

    for prop in list_real_evolvers(member):
        offset = offset_name('all', class_dict, member, prop)
        data  = prop.get('data') or {}
        rank  = int(data.get('rank') or 0)
        if rank == 0:
            if data.get('type') == 'double':
                content += (
                    f"if (.not.nodeAnalytics({offset}) .and. (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. .not.nodeInactives({offset})) .or. (propertyType == propertyTypeInactive .and. nodeInactives({offset})) .or. propertyType == propertyTypeNumerics)) then\n"
                    f" array(offset)=self%{prop['name']}Data\n"
                    " offset=offset+1\n"
                    "end if\n"
                )
            else:
                content += (
                    f"count=self%{prop['name']}Data%serializeCount()\n"
                    f"if (all(.not.nodeAnalytics({offset}:{offset}+count-1)) .and. (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. all(.not.nodeInactives({offset}:{offset}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({offset}:{offset}+count-1))) .or. propertyType == propertyTypeNumerics)) then\n"
                    f" if (count > 0) call self%{prop['name']}Data%serialize(array(offset:offset+count-1))\n"
                    " offset=offset+count\n"
                    "end if\n"
                )
        else:
            content += (
                f"if (allocated(self%{prop['name']}Data)) then\n"
                f"   count=size(self%{prop['name']}Data)\n"
                f"   if (all(.not.nodeAnalytics({offset}:{offset}+count-1)) .and. (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. all(.not.nodeInactives({offset}:{offset}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({offset}:{offset}+count-1))) .or. propertyType == propertyTypeNumerics)) then\n"
                f"    array(offset:offset+count-1)=reshape(self%{prop['name']}Data,[count])\n"
                "    offset=offset+count\n"
                "   end if\n"
                "end if\n"
            )

    function['content'] = content
    _bind(build, impl_type, function, 'serializeValues')


def Implementation_ODE_Deserialize_Values(build, class_dict, member):
    """Generate `<impl>DeserializeValues`.
    """
    impl_type   = _impl_type(class_dict, member)
    has_extends = isinstance(member.get('extends'), dict)
    has_evolve  = has_real_evolvers(member)
    has_nontriv = has_real_non_trivial_evolvers(member)

    function = {
        'type':        'void',
        'name':        impl_type + 'DeserializeValues',
        'description': (
            f"Deserialize evolvable properties of a {member['name']} "
            f"implementation of the {class_dict['name']} component from "
            "array."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['intent(in   )', 'dimension(:)'],
                'variables':  ['array'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['propertyType'],
            },
        ],
    }

    required = ['offset']
    if has_extends or has_nontriv:
        required.append('count')
    if required:
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': list(required),
        })
    if not has_extends:
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': ['i'],
        })

    unused = [] if (has_extends or has_evolve) else ['array', 'propertyType']

    content = ''
    if unused:
        content += "!$GLC attributes unused :: " + ",".join(unused) + "\n"
    if 'offset' in required:
        content += "offset=1\n"

    if has_extends:
        ext = member['extends']
        parent_type = (
            'nodeComponent'
            + _ucfirst(ext.get('class', ''))
            + _ucfirst(ext.get('name', ''))
        )
        content += (
            f"count=self%{parent_type}%serializeCount   (      propertyType)\n"
            "if (count > 0) then\n"
            f" call self%{parent_type}%deserializeValues(array,propertyType)\n"
            " offset=offset+count\n"
            "end if\n"
        )
    elif class_dict['name'] in (build.get('componentClassListActive') or []):
        offset = offset_name('all', class_dict['name'], 'floatRank0MetaProperties')
        cn = class_dict['name']
        content += (
            f"if (allocated({cn}FloatRank0MetaPropertyNames)) then\n"
            f" do i=1,size({cn}FloatRank0MetaPropertyNames)\n"
            f"  if ({cn}FloatRank0MetaPropertyEvolvable(i).and..not.nodeAnalytics({offset}(i)).and.(propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. .not.nodeInactives({offset}(i))) .or. (propertyType == propertyTypeInactive .and. nodeInactives({offset}(i))) .or. propertyType == propertyTypeNumerics)) then\n"
            "   self%floatRank0MetaProperties(i)=array(offset)\n"
            "   offset=offset+1\n"
            "  end if\n"
            " end do\n"
            "end if\n"
        )

    for prop in list_real_evolvers(member):
        offset = offset_name('all', class_dict, member, prop)
        data  = prop.get('data') or {}
        rank  = int(data.get('rank') or 0)
        if rank == 0:
            if data.get('type') == 'double':
                content += (
                    f"if (.not.nodeAnalytics({offset}) .and. (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. .not.nodeInactives({offset})) .or. (propertyType == propertyTypeInactive .and. nodeInactives({offset})) .or. propertyType == propertyTypeNumerics)) then\n"
                    f" self%{prop['name']}Data=array(offset)\n"
                    " offset=offset+1\n"
                    "end if\n"
                )
            else:
                content += (
                    f"count=self%{prop['name']}Data%serializeCount()\n"
                    f"if (all(.not.nodeAnalytics({offset}:{offset}+count-1)) .and. (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. all(.not.nodeInactives({offset}:{offset}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({offset}:{offset}+count-1))) .or. propertyType == propertyTypeNumerics)) then\n"
                    f" if (count > 0) call self%{prop['name']}Data%deserialize(array(offset:offset+count-1))\n"
                    " offset=offset+count\n"
                    "end if\n"
                )
        else:
            content += (
                f"if (allocated(self%{prop['name']}Data)) then\n"
                f"   count=size(self%{prop['name']}Data)\n"
                f"   if (all(.not.nodeAnalytics({offset}:{offset}+count-1)) .and. (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. all(.not.nodeInactives({offset}:{offset}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({offset}:{offset}+count-1))) .or. propertyType == propertyTypeNumerics)) then\n"
                f"    self%{prop['name']}Data=reshape(array(offset:offset+count-1),shape(self%{prop['name']}Data))\n"
                "    offset=offset+count\n"
                "   end if\n"
                "end if\n"
            )

    function['content'] = content
    _bind(build, impl_type, function, 'deserializeValues')


def Implementation_ODE_Serialize_NonNegative(build, class_dict, member):
    """Generate `<impl>SerializeNonNegative`.
    """
    impl_type   = _impl_type(class_dict, member)
    has_extends = isinstance(member.get('extends'), dict)
    has_evolve  = has_real_evolvers(member)
    has_nontriv = has_real_non_trivial_evolvers(member)

    function = {
        'type':        'void',
        'name':        impl_type + 'SerializeNonNegative',
        'description': (
            f"Serialize non-negative status of evolvable properties of a "
            f"\\mono{{{member['name']}}} implementation of the "
            f"\\mono{{{class_dict['name']}}} component to array."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'logical',
                'attributes': ['intent(  out)', 'dimension(:)'],
                'variables':  ['array'],
            },
        ],
    }

    required = ['offset']
    if has_extends or has_nontriv:
        required.append('count')
    if not has_extends:
        required.append('i')
    if required:
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': list(required),
        })

    unused = [] if (has_extends or has_evolve) else ['array']

    content = ''
    if unused:
        content += "!$GLC attributes unused :: " + ",".join(unused) + "\n"
    if 'offset' in required:
        content += "offset=1\n"

    if has_extends:
        ext = member['extends']
        parent_type = (
            'nodeComponent'
            + _ucfirst(ext.get('class', ''))
            + _ucfirst(ext.get('name', ''))
        )
        content += (
            f"count=self%{parent_type}%serializeCount(propertyTypeActive)\n"
            "if (count > 0) then\n"
            f" call self%{parent_type}%serializeNonNegative(array)\n"
            " offset=offset+count\n"
            "end if\n"
        )
    elif class_dict['name'] in (build.get('componentClassListActive') or []):
        offset = offset_name('all', class_dict['name'], 'floatRank0MetaProperties')
        cn = class_dict['name']
        content += (
            f"if (allocated({cn}FloatRank0MetaPropertyNames)) then\n"
            f" do i=1,size({cn}FloatRank0MetaPropertyNames)\n"
            f"  if ({cn}FloatRank0MetaPropertyEvolvable(i).and..not.nodeAnalytics({offset}(i)) .and. .not.nodeInactives({offset}(i))) then\n"
            "   array(offset)=.false. ! Currently no support for non-negative meta-properties.\n"
            "   offset=offset+1\n"
            "  end if\n"
            " end do\n"
            "end if\n"
        )

    for prop in list_real_evolvers(member):
        offset = offset_name('all', class_dict, member, prop)
        is_non_negative = (
            ".true." if (prop.get('attributes') or {}).get('isNonNegative') else ".false."
        )
        data  = prop.get('data') or {}
        rank  = int(data.get('rank') or 0)
        if rank == 0:
            if data.get('type') == 'double':
                content += (
                    f"if (.not.nodeAnalytics({offset}) .and. .not.nodeInactives({offset})) then\n"
                    f" array(offset)={is_non_negative}\n"
                    " offset=offset+1\n"
                    "end if\n"
                )
            else:
                content += (
                    f"count=self%{prop['name']}Data%serializeCount()\n"
                    f"if (all(.not.nodeAnalytics({offset}:{offset}+count-1)) .and. all(.not.nodeInactives({offset}:{offset}+count-1))) then\n"
                    f" if (count > 0) array(offset:offset+count-1)={is_non_negative}\n"
                    " offset=offset+count\n"
                    "end if\n"
                )
        else:
            content += (
                f"if (allocated(self%{prop['name']}Data)) then\n"
                f"   count=size(self%{prop['name']}Data)\n"
                f"   if (all(.not.nodeAnalytics({offset}:{offset}+count-1)) .and. all(.not.nodeInactives({offset}:{offset}+count-1))) then\n"
                f"    array(offset:offset+count-1)={is_non_negative}\n"
                "    offset=offset+count\n"
                "   end if\n"
                "end if\n"
            )

    function['content'] = content
    _bind(build, impl_type, function, 'serializeNonNegative')


def Implementation_ODE_Offsets(build, class_dict, member):
    """Generate `<impl>SerializeOffsets`.
    """
    impl_type   = _impl_type(class_dict, member)
    has_extends = isinstance(member.get('extends'), dict)
    has_evolve  = has_real_evolvers(member)
    has_nontriv = has_real_non_trivial_evolvers(member)

    function = {
        'type':        'void',
        'name':        impl_type + 'SerializeOffsets',
        'description': (
            f"Compute offsets into serialization arrays for a "
            f"\\mono{{{member['name']}}} implementation of the "
            f"\\mono{{{class_dict['name']}}} component."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(inout)'],
                'variables':  ['count', 'countSubset'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['propertyType'],
            },
        ],
    }
    if not has_extends:
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': ['i'],
        })

    unused = []
    if not (has_extends or has_nontriv):
        unused.append('self')
    if not (has_extends or has_evolve):
        unused.extend(['count', 'countSubset', 'propertyType'])

    if has_evolve:
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': ['propertySize'],
        })

    content = ''
    if unused:
        content += "!$GLC attributes unused :: " + ",".join(unused) + "\n"

    if has_extends:
        ext = member['extends']
        parent_type = (
            'nodeComponent'
            + _ucfirst(ext.get('class', ''))
            + _ucfirst(ext.get('name', ''))
        )
        content += (
            f"call self%{parent_type}%serializationOffsets(count,countSubset,propertyType)\n"
        )
    elif class_dict['name'] in (build.get('componentClassListActive') or []):
        cn = class_dict['name']
        offset_all      = offset_name('all',      cn, 'floatRank0MetaProperties')
        offset_active   = offset_name('active',   cn, 'floatRank0MetaProperties')
        offset_inactive = offset_name('inactive', cn, 'floatRank0MetaProperties')
        content += f"if (allocated({cn}FloatRank0MetaPropertyNames)) then\n"
        for status in ('all', 'active', 'inactive'):
            offset = offset_name(status, cn, 'floatRank0MetaProperties')
            content += (
                f" if (.not.allocated({offset})) then\n"
                f"  allocate({offset}({cn}FloatRank0MetaPropertyCount))\n"
                f" else if (size({offset}) /= {cn}FloatRank0MetaPropertyCount) then\n"
                f"  deallocate({offset})\n"
                f"  allocate({offset}({cn}FloatRank0MetaPropertyCount))\n"
                " end if\n"
            )
        content += (
            f" do i=1,size({cn}FloatRank0MetaPropertyNames)\n"
            f"  if (.not.{cn}FloatRank0MetaPropertyEvolvable(i)) cycle\n"
            "  if      (propertyType == propertyTypeAll     ) then\n"
            f"                                    {offset_all}     (i)=count      +1\n"
            "  else if (propertyType == propertyTypeInactive) then\n"
            f"   if (     nodeInactives(count+1).and..not.nodeAnalytics(count+1)) {offset_inactive}(i)=countSubset+1\n"
            "  else if (propertyType == propertyTypeActive  ) then\n"
            f"   if (.not.nodeInactives(count+1).and..not.nodeAnalytics(count+1)) {offset_active}  (i)=countSubset+1\n"
            "  else if (propertyType == propertyTypeNumerics) then\n"
            f"   if (                                .not.nodeAnalytics(count+1)) {offset_active}  (i)=countSubset+1\n"
            "  end if\n"
            "  if (propertyType /= propertyTypeAll) then\n"
            "   if (((nodeInactives(count+1) .and. propertyType == propertyTypeInactive) .or. (.not.nodeInactives(count+1) .and. propertyType == propertyTypeActive) .or. propertyType == propertyTypeNumerics) .and. .not.nodeAnalytics(count+1)) countSubset=countSubset+1\n"
            "  end if\n"
            "  count=count+1\n"
            " end do\n"
            "end if\n"
        )

    for prop in list_real_evolvers(member):
        offset_all_p      = offset_name('all',      class_dict, member, prop)
        offset_active_p   = offset_name('active',   class_dict, member, prop)
        offset_inactive_p = offset_name('inactive', class_dict, member, prop)
        data  = prop.get('data') or {}
        rank  = int(data.get('rank') or 0)
        if rank == 0:
            if data.get('type') == 'double':
                content += "propertySize=1\n"
            else:
                content += (
                    f"propertySize=self%{prop['name']}Data%serializeCount()\n"
                )
        else:
            content += f"propertySize=size(self%{prop['name']}Data)\n"
        content += (
            "if (propertySize > 0) then\n"
            " if      (propertyType == propertyTypeAll     ) then\n"
            f"                                   {offset_all_p}     =count      +1\n"
            " else if (propertyType == propertyTypeInactive) then\n"
            f"  if (     nodeInactives(count+1).and..not.nodeAnalytics(count+1)) {offset_inactive_p}=countSubset+1\n"
            " else if (propertyType == propertyTypeActive  ) then\n"
            f"  if (.not.nodeInactives(count+1).and..not.nodeAnalytics(count+1)) {offset_active_p}  =countSubset+1\n"
            " else if (propertyType == propertyTypeNumerics) then\n"
            f"  if (                                .not.nodeAnalytics(count+1)) {offset_active_p}  =countSubset+1\n"
            " end if\n"
            " if (propertyType /= propertyTypeAll) then\n"
            "  if (((nodeInactives(count+1) .and. propertyType == propertyTypeInactive) .or. (.not.nodeInactives(count+1) .and. propertyType == propertyTypeActive) .or. propertyType == propertyTypeNumerics) .and. .not.nodeAnalytics(count+1)) countSubset=countSubset+propertySize\n"
            " end if\n"
            " count=count+propertySize\n"
            "end if\n"
        )

    function['content'] = content
    _bind(build, impl_type, function, 'serializationOffsets')


def Implementation_ODE_Offset_Variables(build, class_dict, member):
    """Declare the `offsetAll<...>` / `offsetActive<...>` / `offsetInactive<...>`
    module-scope variables.  Meta-property offsets are emitted only for
    the `null` member of each class so that we get exactly one copy.
    """
    name = class_dict['name']
    if name in (build.get('componentClassListActive') or []):
        if member.get('name') == 'null':
            for status in ('all', 'active', 'inactive'):
                offset = offset_name(status, name, 'floatRank0MetaProperties')
                build.setdefault('variables', []).append({
                    'intrinsic':  'integer',
                    'ompPrivate': True,
                    'attributes': ['allocatable', 'dimension(:)'],
                    'variables':  [offset],
                })

    for prop in list_real_evolvers(member):
        for status in ('all', 'active', 'inactive'):
            # The concatenated class+member is intentionally passed as the
            # "componentName" arg of the 3-arg form — this fixes the
            # generated offset variable names.
            offset = offset_name(status, name + member['name'], prop['name'])
            build.setdefault('variables', []).append({
                'intrinsic':  'integer',
                'ompPrivate': True,
                'variables':  [offset],
            })


def Implementation_ODE_Rate_Variables(build):
    """Declare the module-scope state arrays consumed by the
    per-implementation ODE methods (rates, scales, inactive flags,
    analytic flags).
    """
    build.setdefault('variables', []).extend([
        {
            'intrinsic':  'integer',
            'ompPrivate': True,
            'variables':  [
                'nodeSerializationCount',
                'nodeSerializationCountActive',
                'nodeSerializationCountInactive',
            ],
        },
        {
            'intrinsic':  'double precision',
            'attributes': ['allocatable', 'dimension(:)'],
            'ompPrivate': True,
            'variables':  ['nodeScales', 'nodeRates', 'nodeRatesActives'],
        },
        {
            'intrinsic':  'logical',
            'attributes': ['allocatable', 'dimension(:)'],
            'ompPrivate': True,
            'variables':  ['nodeInactives', 'nodeAnalytics'],
        },
    ])


# Hook registrations land at the end of this file once every hook is
# defined.
register('implementationODESolver', 'implementationIteratedFunctions',
         Implementation_ODE_Serialize_Count)
register('implementationODESolver', 'implementationIteratedFunctions',
         Implementation_ODE_Serialize_Values)
register('implementationODESolver', 'implementationIteratedFunctions',
         Implementation_ODE_Deserialize_Values)
register('implementationODESolver', 'implementationIteratedFunctions',
         Implementation_ODE_Serialize_NonNegative)
register('implementationODESolver', 'implementationIteratedFunctions',
         Implementation_ODE_Name_From_Index)
register('implementationODESolver', 'implementationIteratedFunctions',
         Implementation_ODE_Offsets)
register('implementationODESolver', 'implementationIteratedFunctions',
         Implementation_ODE_Offset_Variables)
register('implementationODESolver', 'functions',
         Implementation_ODE_Rate_Variables)
