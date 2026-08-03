"""Per-implementation (de)serialization methods (ASCII / XML / raw / deserialize raw).

Andrew Benson (ported to Python 2026)

Note: three of the four hooks deliberately emit no `!$GLC attributes
unused` lines — the unused-variable list was never populated for them
historically, so those branches are omitted here.
"""



from Galacticus.Build.Components                       import Utils as _Utils
from Galacticus.Build.Components.Utils                 import (
    register,
    is_intrinsic,
    _component_properties,
)
from Galacticus.Build.Components.Classes.MetaProperties import meta_property_types


_FORMAT_LABEL_QUOTED = {
    'double':      "'(e22.16)'",
    'integer':     "'(i8)'",
    'longInteger': "'(i16)'",
    'logical':     "'(l1)'",
}

_FORMAT_LABEL_BARE = {
    'double':      "(e12.6)",
    'integer':     "(i8)",
    'longInteger': "(i16)",
    'logical':     "(l1)",
}


def _meta_format(mpt, format_table):
    """Return the format spec for a meta-property type, given a format
    table.  Mirrors the four-way conditional repeated in each hook.
    """
    if mpt['label'] == 'float':
        return format_table['double']
    if mpt['label'] == 'integer':
        return format_table['integer']
    if mpt['label'] == 'longInteger':
        return format_table['longInteger']
    raise RuntimeError(
        "Implementations/Serialization: unknown meta-property type"
    )


def Implementation_Serialize_ASCII(build, class_dict, member):
    """Generate `nodeComponent<Class><Member>SerializeASCII`."""
    cap_class   = _ucfirst(class_dict['name'])
    cap_member  = _ucfirst(member['name'])
    impl_type   = 'nodeComponent' + cap_class + cap_member

    function = {
        'type':        'void',
        'name':        impl_type + 'SerializeASCII',
        'description': (
            f"Serialize the contents of a {member['name']} implementation of "
            f"the {class_dict['name']} component to ASCII."
        ),
        'modules':     ['Display', 'ISO_Varying_String', 'String_Handling'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'type',
                'type':       'enumerationVerbosityLevelType',
                'variables':  ['verbosityLevel'],
                'attributes': ['intent(in   )'],
            },
            {
                'intrinsic':  'type',
                'type':       'varying_string',
                'variables':  ['message'],
            },
            {
                'intrinsic':  'character',
                'type':       'len=22',
                'variables':  ['label'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i', 'j'],
            },
        ],
        'content':     '',
    }

    content = ''

    if member['name'] != 'null':
        if isinstance(member.get('extends'), dict):
            ext = member['extends']
            parent_type = (
                'nodeComponent'
                + _ucfirst(ext.get('class', ''))
                + _ucfirst(ext.get('name', ''))
            )
            content += (
                f"call self%{parent_type}%serializeASCII(verbosityLevel)\n"
            )
        padding_len = max(
            (_Utils.fully_qualified_name_length_max or 0)
            - len(class_dict['name']), 0
        )
        padding = ' ' * padding_len
        content += (
            f"call displayIndent('{class_dict['name']}: "
            f"{padding}{member['name']}',verbosityLevel)\n"
        )
        for prop in _component_properties(member):
            attrs = prop.get('attributes') or {}
            if attrs.get('isVirtual'):
                continue
            data  = prop.get('data') or {}
            ptype = data.get('type')
            rank  = int(data.get('rank') or 0)
            name_length = len(prop['name'])
            if rank == 0:
                if is_intrinsic(ptype):
                    fmt = _FORMAT_LABEL_QUOTED[ptype]
                    content += (
                        f"write (label,{fmt}) self%{prop['name']}Data\n"
                        f"message='{prop['name']}: '"
                        f"//repeat(' ',propertyNameLengthMax-{name_length})"
                        f"//label\n"
                        "call displayMessage(message,verbosityLevel)\n"
                    )
                else:
                    content += (
                        f"message='{prop['name']}:'\n"
                        "call displayIndent(message,verbosityLevel)\n"
                        f"call self%{prop['name']}Data%dump(verbosityLevel)\n"
                        "call displayUnindent('end',verbosityLevel)\n"
                    )
            elif rank == 1:
                if is_intrinsic(ptype):
                    fmt = _FORMAT_LABEL_QUOTED[ptype]
                    content += (
                        f"do i=1,size(self%{prop['name']}Data)\n"
                        "   write (label,'(i3)') i\n"
                        f"   message='{prop['name']}: '"
                        f"//repeat(' ',propertyNameLengthMax-{name_length})"
                        f"//trim(label)\n"
                        f"   write (label,{fmt}) self%{prop['name']}Data(i)\n"
                        "   message=message//': '//label\n"
                        "   call displayMessage(message,verbosityLevel)\n"
                        "end do\n"
                    )
                else:
                    content += (
                        f"do i=1,size(self%{prop['name']}Data)\n"
                        "   write (label,'(i3)') i\n"
                        f"   message='{prop['name']}: '"
                        f"//repeat(' ',propertyNameLengthMax-{name_length})"
                        f"//trim(label)\n"
                        "   call displayIndent(message,verbosityLevel)\n"
                        f"   call self%{prop['name']}Data(i)%dump(verbosityLevel)\n"
                        "   call displayUnindent('end',verbosityLevel)\n"
                        "end do\n"
                    )

    if class_dict['name'] in (build.get('componentClassListActive') or []):
        for mpt in meta_property_types:
            cap_label = _ucfirst(mpt['label'])
            rank      = mpt['rank']
            prefix    = f"{cap_label}Rank{rank}"
            var_pref  = f"{class_dict['name']}{prefix}"
            fmt       = _meta_format(mpt, _FORMAT_LABEL_QUOTED)
            if rank == 0:
                content += (
                    f"if (allocated({var_pref}MetaPropertyNames)) then\n"
                    f" do i=1,size({var_pref}MetaPropertyNames)\n"
                    f"  write (label,{fmt}) self%{prefix}MetaProperties(i)\n"
                    f"  message=trim({var_pref}MetaPropertyNames(i))//': '"
                    f"//repeat(' ',propertyNameLengthMax-len_trim("
                    f"{var_pref}MetaPropertyNames(i)))//label\n"
                    "  call displayMessage(message,verbosityLevel)\n"
                    " end do\n"
                    "end if\n"
                )
            elif rank == 1:
                content += (
                    f"if (allocated({var_pref}MetaPropertyNames)) then\n"
                    f" do i=1,size({var_pref}MetaPropertyNames)\n"
                    f"  do j=1,size( self%{prefix}MetaProperties(i)%values)\n"
                    "   write (label,'(i3)') j\n"
                    f"   message=trim({var_pref}MetaPropertyNames(i))//': '"
                    f"//repeat(' ',propertyNameLengthMax-len_trim("
                    f"{var_pref}MetaPropertyNames(i)))//trim(label)\n"
                    f"   write (label,{fmt}) self%{prefix}MetaProperties(i)%values(j)\n"
                    "   message=message//': '//label\n"
                    "   call displayMessage(message,verbosityLevel)\n"
                    "  end do\n"
                    " end do\n"
                    "end if\n"
                )
            else:
                raise RuntimeError(
                    "Serialize_ASCII: unsupported meta-property rank"
                )

    content += "call displayUnindent('done',verbosityLevel)\n"
    function['content'] = content
    _bind(build, impl_type, function, 'serializeASCII')


def Implementation_Serialize_XML(build, class_dict, member):
    """Generate `nodeComponent<Class><Member>SerializeXML`."""
    cap_class  = _ucfirst(class_dict['name'])
    cap_member = _ucfirst(member['name'])
    impl_type  = 'nodeComponent' + cap_class + cap_member

    function = {
        'type':        'void',
        'name':        impl_type + 'SerializeXML',
        'description': (
            f"Serialize the contents of a {member['name']} implementation of "
            f"the {class_dict['name']} component to XML."
        ),
        'modules':     ['ISO_Varying_String'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['fileHandle'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i', 'j'],
            },
        ],
        'content':     '',
    }

    content = ''

    if isinstance(member.get('extends'), dict):
        ext = member['extends']
        parent_type = (
            'nodeComponent'
            + _ucfirst(ext.get('class', ''))
            + _ucfirst(ext.get('name', ''))
        )
        content += (
            f"call self%{parent_type}%serializeXML(fileHandle)\n"
        )
    content += (
        f"write (fileHandle,'(a)') '  <{class_dict['name']} "
        f"type=\"{member['name']}\">'\n"
    )

    for prop in _component_properties(member):
        attrs = prop.get('attributes') or {}
        data  = prop.get('data') or {}
        ptype = data.get('type')
        rank  = int(data.get('rank') or 0)
        # Keep iff `not isVirtual or data.rank == 0`.
        if attrs.get('isVirtual') and rank != 0:
            continue
        if not attrs.get('isVirtual'):
            if rank == 0:
                if is_intrinsic(ptype):
                    fmt = _FORMAT_LABEL_BARE[ptype]
                    content += (
                        f"write (fileHandle,'(a,{fmt},a)') "
                        f"'   <{prop['name']}>',self%{prop['name']}Data,"
                        f"'</{prop['name']}>'\n"
                    )
                else:
                    content += (
                        f"write (fileHandle,'(a)') '   <{prop['name']}>'\n"
                        f"write (fileHandle,'(a)') '   </{prop['name']}>'\n"
                    )
            elif rank == 1:
                if is_intrinsic(ptype):
                    fmt = _FORMAT_LABEL_BARE[ptype]
                    content += (
                        f"do i=1,size(self%{prop['name']}Data)\n"
                        f"   write (fileHandle,'(a,{fmt},a)') "
                        f"'   <{prop['name']}>',self%{prop['name']}Data(i),"
                        f"'</{prop['name']}>'\n"
                        "end do\n"
                    )
                else:
                    content += (
                        f"do i=1,size(self%{prop['name']}Data)\n"
                        f"   write (fileHandle,'(a)') '   <{prop['name']}>'\n"
                        f"   write (fileHandle,'(a)') '   </{prop['name']}>'\n"
                        "end do\n"
                    )
        else:
            # Virtual + rank-0 + isGettable + intrinsic: call the getter.
            if (
                attrs.get('isGettable')
                and rank == 0
                and is_intrinsic(ptype)
            ):
                fmt = _FORMAT_LABEL_BARE[ptype]
                content += (
                    f"write (fileHandle,'(a,{fmt},a)') "
                    f"'   <{prop['name']}>',self%{prop['name']}(),"
                    f"'</{prop['name']}>'\n"
                )

    if class_dict['name'] in (build.get('componentClassListActive') or []):
        for mpt in meta_property_types:
            cap_label = _ucfirst(mpt['label'])
            rank      = mpt['rank']
            prefix    = f"{cap_label}Rank{rank}"
            var_pref  = f"{class_dict['name']}{prefix}"
            fmt       = _meta_format(mpt, _FORMAT_LABEL_BARE)
            if rank == 0:
                content += (
                    f"if (allocated({var_pref}MetaPropertyNames)) then\n"
                    f" do i=1,size(({var_pref}MetaPropertyNames))\n"
                    f"  write (fileHandle,'(a,a,a,{fmt},a,a,a)') "
                    f"'   <'//char({var_pref}MetaPropertyNames(i))//'>',"
                    f"self%{prefix}MetaProperties(i),"
                    f"'</'//char({var_pref}MetaPropertyNames(i))//'>'\n"
                    " end do\n"
                    "end if\n"
                )
            elif rank == 1:
                content += (
                    f"if (allocated({var_pref}MetaPropertyNames)) then\n"
                    f" do i=1,size(({var_pref}MetaPropertyNames))\n"
                    f"  do j=1,size(self%{prefix}MetaProperties(i)%values)\n"
                    f"   write (fileHandle,'(a,a,a,{fmt},a,a,a)') "
                    f"'   <'//char({var_pref}MetaPropertyNames(i))//'>',"
                    f"self%{prefix}MetaProperties(i)%values(j),"
                    f"'</'//char({var_pref}MetaPropertyNames(i))//'>'\n"
                    "  end do\n"
                    " end do\n"
                    "end if\n"
                )
            else:
                raise RuntimeError(
                    "Serialize_XML: unsupported meta-property rank"
                )

    content += (
        f"write (fileHandle,'(a)') '  </{class_dict['name']}>'\n"
    )
    function['content'] = content
    _bind(build, impl_type, function, 'serializeXML')


def Implementation_Serialize_Raw(build, class_dict, member):
    """Generate `nodeComponent<Class><Member>SerializeRaw`."""
    cap_class  = _ucfirst(class_dict['name'])
    cap_member = _ucfirst(member['name'])
    impl_type  = 'nodeComponent' + cap_class + cap_member

    function = {
        'type':        'void',
        'name':        impl_type + 'SerializeRaw',
        'description': (
            f"Serialize the contents of a {member['name']} implementation "
            f"of the {class_dict['name']} component to raw (binary) file."
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
                'attributes': ['intent(in   )'],
                'variables':  ['fileHandle'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
        ],
        'content':     '',
    }

    content = ''
    if isinstance(member.get('extends'), dict):
        ext = member['extends']
        parent_type = (
            'nodeComponent'
            + _ucfirst(ext.get('class', ''))
            + _ucfirst(ext.get('name', ''))
        )
        content += (
            f"call self%{parent_type}%serializeRaw(fileHandle)\n"
        )

    for prop in _component_properties(member):
        attrs = prop.get('attributes') or {}
        if attrs.get('isVirtual'):
            continue
        data  = prop.get('data') or {}
        ptype = data.get('type')
        rank  = int(data.get('rank') or 0)
        if rank == 0:
            if is_intrinsic(ptype):
                content += f"write (fileHandle) self%{prop['name']}Data\n"
            else:
                content += (
                    f"call self%{prop['name']}Data%dumpRaw(fileHandle)\n"
                )
        elif rank == 1:
            content += (
                f"write (fileHandle) allocated(self%{prop['name']}Data)\n"
                f"if (allocated(self%{prop['name']}Data)) then\n"
                f"   write (fileHandle) size(self%{prop['name']}Data)\n"
            )
            if is_intrinsic(ptype):
                content += f"write (fileHandle) self%{prop['name']}Data\n"
            else:
                content += (
                    f"   do i=1,size(self%{prop['name']}Data)\n"
                    f"      call self%{prop['name']}Data(i)%dumpRaw(fileHandle)\n"
                    "   end do\n"
                )
            content += "end if\n"

    if class_dict['name'] in (build.get('componentClassListActive') or []):
        # Meta-properties are written in a self-describing form: for each
        # meta-property type we write a count, and then for each meta-property
        # its name (length followed by characters) alongside its value(s). This
        # makes the on-disk format independent of the set of meta-properties
        # registered by the (task-dependent) run performing the store, so that a
        # run performing the restore can map values back by name regardless of
        # which meta-properties it has itself registered. See the matching
        # deserialization in `Implementation_Deserialize_Raw`.
        for mpt in meta_property_types:
            cap_label = _ucfirst(mpt['label'])
            rank      = mpt['rank']
            prefix    = f"{cap_label}Rank{rank}"
            var_pref  = f"{class_dict['name']}{prefix}"
            if rank == 0:
                content += (
                    f"if (allocated({var_pref}MetaPropertyNames)) then\n"
                    f" write (fileHandle) size({var_pref}MetaPropertyNames)\n"
                    f" do i=1,size({var_pref}MetaPropertyNames)\n"
                    f"  write (fileHandle) len({var_pref}MetaPropertyNames(i))\n"
                    f"  write (fileHandle) char({var_pref}MetaPropertyNames(i))\n"
                    f"  write (fileHandle) self%{prefix}MetaProperties(i)\n"
                    " end do\n"
                    "else\n"
                    " write (fileHandle) 0\n"
                    "end if\n"
                )
            elif rank == 1:
                content += (
                    f"if (allocated({var_pref}MetaPropertyNames)) then\n"
                    f" write (fileHandle) size({var_pref}MetaPropertyNames)\n"
                    f" do i=1,size({var_pref}MetaPropertyNames)\n"
                    f"  write (fileHandle) len({var_pref}MetaPropertyNames(i))\n"
                    f"  write (fileHandle) char({var_pref}MetaPropertyNames(i))\n"
                    f"  write (fileHandle) size(self%{prefix}MetaProperties(i)%values)\n"
                    f"  write (fileHandle) self%{prefix}MetaProperties(i)%values\n"
                    " end do\n"
                    "else\n"
                    " write (fileHandle) 0\n"
                    "end if\n"
                )
            else:
                raise RuntimeError(
                    "Serialize_Raw: unsupported meta-property rank"
                )

    function['content'] = content
    _bind(build, impl_type, function, 'serializeRaw')


def Implementation_Deserialize_Raw(build, class_dict, member):
    """Generate `nodeComponent<Class><Member>DeserializeRaw`."""
    cap_class  = _ucfirst(class_dict['name'])
    cap_member = _ucfirst(member['name'])
    impl_type  = 'nodeComponent' + cap_class + cap_member

    function = {
        'type':        'void',
        'name':        impl_type + 'DeserializeRaw',
        'description': (
            f"Deserialize the contents of a {member['name']} implementation "
            f"of the {class_dict['name']} component from raw (binary) file."
        ),
        'modules':     ['ISO_Varying_String', 'Kind_Numbers'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       impl_type,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['fileHandle'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i', 'metaPropertySize'],
            },
        ],
        'content':     '',
    }

    active = class_dict['name'] in (build.get('componentClassListActive') or [])
    if active:
        # Scratch variables used to read the self-describing meta-property
        # block and map stored values back onto the currently-registered
        # meta-properties by name (see `Implementation_Serialize_Raw`).
        function['variables'].append({
            'intrinsic': 'integer',
            'variables': ['j', 'metaPropertyCount', 'metaPropertyNameLength',
                          'metaPropertyMatch'],
        })
        function['variables'].append({
            'intrinsic':  'character',
            'type':       'len=:',
            'attributes': ['allocatable'],
            'variables':  ['metaPropertyName'],
        })
        for mpt in meta_property_types:
            if mpt['rank'] != 0:
                continue
            cap_label = _ucfirst(mpt['label'])
            var = {
                'intrinsic': mpt['intrinsic'],
                'variables': [f"metaPropertyScalar{cap_label}"],
            }
            if mpt.get('type'):
                var['type'] = mpt['type']
            function['variables'].append(var)
        for mpt in meta_property_types:
            if mpt['rank'] != 1:
                continue
            cap_label = _ucfirst(mpt['label'])
            var = {
                'intrinsic':  mpt['intrinsic'],
                'attributes': ['allocatable', 'dimension(:)'],
                'variables':  [f"metaPropertyBuffer{cap_label}"],
            }
            if mpt.get('type'):
                var['type'] = mpt['type']
            function['variables'].append(var)

    has_rank1_real = any(
        not (p.get('attributes') or {}).get('isVirtual')
        and int((p.get('data') or {}).get('rank') or 0) == 1
        for p in _component_properties(member)
    )
    if has_rank1_real:
        function['variables'].extend([
            {
                'intrinsic': 'integer',
                'variables': ['arraySize'],
            },
            {
                'intrinsic': 'logical',
                'variables': ['isAllocated'],
            },
        ])

    content = ''
    if isinstance(member.get('extends'), dict):
        ext = member['extends']
        parent_type = (
            'nodeComponent'
            + _ucfirst(ext.get('class', ''))
            + _ucfirst(ext.get('name', ''))
        )
        content += (
            f"call self%{parent_type}%deserializeRaw(fileHandle)\n"
        )

    for prop in _component_properties(member):
        attrs = prop.get('attributes') or {}
        if attrs.get('isVirtual'):
            continue
        data  = prop.get('data') or {}
        ptype = data.get('type')
        rank  = int(data.get('rank') or 0)
        if rank == 0:
            if is_intrinsic(ptype):
                content += f"read (fileHandle) self%{prop['name']}Data\n"
            else:
                content += (
                    f"call self%{prop['name']}Data%readRaw(fileHandle)\n"
                )
        elif rank == 1:
            content += (
                "read (fileHandle) isAllocated\n"
                "if (isAllocated) then\n"
                "   read (fileHandle) arraySize\n"
            )
            if is_intrinsic(ptype):
                content += (
                    f"   allocate(self%{prop['name']}Data(arraySize))\n"
                    f"   read (fileHandle) self%{prop['name']}Data\n"
                )
            else:
                content += (
                    f"   allocate(self%{prop['name']}Data(arraySize))\n"
                    "   do i=1,arraySize\n"
                    f"      call self%{prop['name']}Data(i)%readRaw(fileHandle)\n"
                    "   end do\n"
                )
            content += "end if\n"

    if active:
        # Read the self-describing meta-property block written by
        # `Implementation_Serialize_Raw`. For each meta-property type we read a
        # count and then, for each stored meta-property, its name and value(s).
        # Every stored entry is consumed from the stream regardless of whether
        # this run has the meta-property registered (so the stream never
        # desyncs), and values are assigned by matching the stored name against
        # the currently-registered names. Registered meta-properties absent from
        # the file retain their default (zero / unallocated) value.
        zero_literal = {
            'float':       '0.0d0',
            'longInteger': '0_kind_int8',
            'integer':     '0',
        }
        for mpt in meta_property_types:
            cap_label = _ucfirst(mpt['label'])
            rank      = mpt['rank']
            prefix    = f"{cap_label}Rank{rank}"
            var_pref  = f"{class_dict['name']}{prefix}"
            match_head = (
                "read (fileHandle) metaPropertyNameLength\n"
                "if (allocated(metaPropertyName)) deallocate(metaPropertyName)\n"
                "allocate(character(len=metaPropertyNameLength) :: "
                "metaPropertyName)\n"
                "read (fileHandle) metaPropertyName\n"
                "metaPropertyMatch=0\n"
                f"if (allocated({var_pref}MetaPropertyNames)) then\n"
                f" do j=1,size({var_pref}MetaPropertyNames)\n"
                f"  if ({var_pref}MetaPropertyNames(j) == metaPropertyName) then\n"
                "   metaPropertyMatch=j\n"
                "   exit\n"
                "  end if\n"
                " end do\n"
                "end if\n"
            )
            if rank == 0:
                content += (
                    f"if (allocated({var_pref}MetaPropertyNames)) then\n"
                    f" if (allocated(self%{prefix}MetaProperties)) "
                    f"deallocate(self%{prefix}MetaProperties)\n"
                    f" allocate(self%{prefix}MetaProperties(size("
                    f"{var_pref}MetaPropertyNames)))\n"
                    f" self%{prefix}MetaProperties={zero_literal[mpt['label']]}\n"
                    "end if\n"
                    "read (fileHandle) metaPropertyCount\n"
                    "do i=1,metaPropertyCount\n"
                    + match_head +
                    " if (metaPropertyMatch > 0) then\n"
                    f"  read (fileHandle) self%{prefix}MetaProperties(metaPropertyMatch)\n"
                    " else\n"
                    f"  read (fileHandle) metaPropertyScalar{cap_label}\n"
                    " end if\n"
                    "end do\n"
                )
            elif rank == 1:
                content += (
                    f"if (allocated({var_pref}MetaPropertyNames)) then\n"
                    f" if (allocated(self%{prefix}MetaProperties)) "
                    f"deallocate(self%{prefix}MetaProperties)\n"
                    f" allocate(self%{prefix}MetaProperties(size("
                    f"{var_pref}MetaPropertyNames)))\n"
                    "end if\n"
                    "read (fileHandle) metaPropertyCount\n"
                    "do i=1,metaPropertyCount\n"
                    + match_head +
                    " read (fileHandle) metaPropertySize\n"
                    " if (metaPropertyMatch > 0) then\n"
                    f"  if (allocated(self%{prefix}MetaProperties(metaPropertyMatch)%values)) "
                    f"deallocate(self%{prefix}MetaProperties(metaPropertyMatch)%values)\n"
                    f"  allocate(self%{prefix}MetaProperties(metaPropertyMatch)%values(metaPropertySize))\n"
                    f"  read (fileHandle) self%{prefix}MetaProperties(metaPropertyMatch)%values\n"
                    " else\n"
                    f"  if (allocated(metaPropertyBuffer{cap_label})) "
                    f"deallocate(metaPropertyBuffer{cap_label})\n"
                    f"  allocate(metaPropertyBuffer{cap_label}(metaPropertySize))\n"
                    f"  read (fileHandle) metaPropertyBuffer{cap_label}\n"
                    " end if\n"
                    "end do\n"
                )
            else:
                raise RuntimeError(
                    "Deserialize_Raw: unsupported meta-property rank"
                )

    function['content'] = content
    _bind(build, impl_type, function, 'deserializeRaw')


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
# Hook registration.  Registration order determines the order of generated
# code — do not reorder.
# ---------------------------------------------------------------------------

register('implementationsSerialization', 'implementationIteratedFunctions',
         Implementation_Serialize_ASCII)
register('implementationsSerialization', 'implementationIteratedFunctions',
         Implementation_Serialize_XML)
register('implementationsSerialization', 'implementationIteratedFunctions',
         Implementation_Serialize_Raw)
register('implementationsSerialization', 'implementationIteratedFunctions',
         Implementation_Deserialize_Raw)
