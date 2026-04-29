# Per-implementation (de)serialization methods (ASCII / XML / raw / deserialize raw).
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Implementations/Serialization.pm.
# Note: the Perl source references `@code::unused` in three of the four
# hooks, but never populates it — so the `!$GLC attributes unused` lines
# in those functions are effectively dead code in the original.  We omit
# those branches here.

import os
import sys

sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))

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
    """Mirrors `Implementation_Serialize_ASCII`."""
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
    """Mirrors `Implementation_Serialize_XML`."""
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
        # Perl filter: keep iff `! isVirtual || data.rank == 0`.
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
            # WART preserved: Perl checks `$property->{'rank'} == 0` (the
            # raw-input rank field) rather than `data.rank` here.
            raw_rank = int(prop.get('rank') or 0)
            if (
                attrs.get('isGettable')
                and raw_rank == 0
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
    """Mirrors `Implementation_Serialize_Raw`."""
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
        for mpt in meta_property_types:
            cap_label = _ucfirst(mpt['label'])
            rank      = mpt['rank']
            prefix    = f"{cap_label}Rank{rank}"
            var_pref  = f"{class_dict['name']}{prefix}"
            if rank == 0:
                content += (
                    f"if (allocated({var_pref}MetaPropertyNames)) "
                    f"write (fileHandle) self%{prefix}MetaProperties\n"
                )
            elif rank == 1:
                content += (
                    f"if (allocated({var_pref}MetaPropertyNames)) then\n"
                    f" do i=1,size({var_pref}MetaPropertyNames)\n"
                    f"  write (fileHandle) size(self%{prefix}MetaProperties(i)%values)\n"
                    f"  write (fileHandle) self%{prefix}MetaProperties(i)%values\n"
                    " end do\n"
                    "end if\n"
                )
            else:
                raise RuntimeError(
                    "Serialize_Raw: unsupported meta-property rank"
                )

    function['content'] = content
    _bind(build, impl_type, function, 'serializeRaw')


def Implementation_Deserialize_Raw(build, class_dict, member):
    """Mirrors `Implementation_Deserialize_Raw`.

    WART preserved: Perl Serialization.pm:571 has `do i=1,arraySize)` —
    a stray closing paren that would produce syntactically invalid
    Fortran for any rank-1 derived-type property in this hook.  This
    branch appears to be never exercised in the current build (no
    rank-1 derived-type properties in `objects.nodes.components.*`),
    so the bug never surfaces.  Mirrored verbatim.
    """
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
                    "   do i=1,arraySize)\n"
                    f"      call self%{prop['name']}Data(i)%readRaw(fileHandle)\n"
                    "   end do\n"
                )
            content += "end if\n"

    if class_dict['name'] in (build.get('componentClassListActive') or []):
        for mpt in meta_property_types:
            cap_label = _ucfirst(mpt['label'])
            rank      = mpt['rank']
            prefix    = f"{cap_label}Rank{rank}"
            var_pref  = f"{class_dict['name']}{prefix}"
            if rank == 0:
                content += (
                    f"if (allocated({var_pref}MetaPropertyNames)) then\n"
                    f" allocate(self%{prefix}MetaProperties(size("
                    f"{var_pref}MetaPropertyNames)))\n"
                    f" read (fileHandle) self%{prefix}MetaProperties\n"
                    "end if\n"
                )
            elif rank == 1:
                content += (
                    f"if (allocated({var_pref}MetaPropertyNames  )) then\n"
                    f" allocate(self%{prefix}MetaProperties(size("
                    f"{var_pref}MetaPropertyNames)))\n"
                    f" do i=1,size({var_pref}MetaPropertyNames)\n"
                    "  read (fileHandle) metaPropertySize\n"
                    f"  allocate(self%{prefix}MetaProperties(i)%values(metaPropertySize))\n"
                    f"  read (fileHandle) self%{prefix}MetaProperties(i)%values\n"
                    " end do\n"
                    "end if\n"
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
# Hook registration.  Order matches Perl Implementations/Serialization.pm:23-27.
# ---------------------------------------------------------------------------

register('implementationsSerialization', 'implementationIteratedFunctions',
         Implementation_Serialize_ASCII)
register('implementationsSerialization', 'implementationIteratedFunctions',
         Implementation_Serialize_XML)
register('implementationsSerialization', 'implementationIteratedFunctions',
         Implementation_Serialize_Raw)
register('implementationsSerialization', 'implementationIteratedFunctions',
         Implementation_Deserialize_Raw)
