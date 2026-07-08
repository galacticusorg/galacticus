"""`treeNode` (de)serialization methods (ASCII / XML / raw binary).

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Components/TreeNodes/Serialization.pm.
Four hooks on the `functions` phase: `serializeASCII`, `serializeXML`,
`serializeRaw`, and `deserializeRaw`.
"""



from Galacticus.Build.Components.Utils import register


_TREE_NODE_POINTERS = (
    'parent', 'firstChild', 'sibling', 'firstSatellite',
    'mergeTarget', 'firstMergee', 'siblingMergee', 'formationNode',
)
_TREE_NODE_STATES = ('isPhysicallyPlausible', 'isSolvable')

# Magic marker written to the raw stream after each component class, and checked
# on restore to detect desynchronization (see Tree_Node_Serialize_Raw /
# Tree_Node_Deserialize_Raw).
_RAW_SENTINEL = 987654321


def Tree_Node_Serialize_ASCII(build):
    """Generate `treeNodeSerializeASCII`.  Mirrors `Tree_Node_Serialize_ASCII`."""
    function = {
        'type':        'void',
        'name':        'treeNodeSerializeASCII',
        'description': "Serialize node content to ASCII.",
        'modules':     [
            'ISO_Varying_String',
            'Display',
            'String_Handling',
        ],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'type',
                'type':       'enumerationVerbosityLevelType',
                'variables':  ['verbosityLevel'],
                'attributes': ['intent(in   )', 'optional'],
            },
            {
                'intrinsic':  'type',
                'type':       'enumerationVerbosityLevelType',
                'variables':  ['verbosityLevel_'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
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
        ],
    }

    content = (
        "verbosityLevel_=verbosityLevelStandard\n"
        "if (present(verbosityLevel)) verbosityLevel_=verbosityLevel\n"
        "message='Dumping node '\n"
        "message=message//self%index()\n"
        "call displayIndent(message,verbosityLevel_)\n"
        "message='host tree: '\n"
        "if (associated(self%hostTree)) then\n"
        " message=message//self%hostTree%index\n"
        "else\n"
        " message=message//'unhosted'\n"
        "end if\n"
        "call displayMessage(message,verbosityLevel_)\n"
        "call displayIndent('pointers',verbosityLevel_)\n"
    )
    for pointer in _TREE_NODE_POINTERS:
        pad = ' ' * (14 - len(pointer))
        content += (
            f"message='{pad}{pointer}: '\n"
            f"message=message//self%{pointer}%index()\n"
            "call displayMessage(message,verbosityLevel_)\n"
        )
    content += (
        "call displayUnindent('done',verbosityLevel_)\n"
        "call displayIndent('state',verbosityLevel_)\n"
    )
    for state in _TREE_NODE_STATES:
        pad = ' ' * (22 - len(state))
        content += (
            f"message='{pad}{state}: '\n"
            f"write (label,'(l1)') self%{state}\n"
            "message=message//trim(adjustl(label))\n"
        )
    content += (
        "call displayMessage(message,verbosityLevel_)\n"
        "call displayUnindent('done',verbosityLevel_)\n"
    )
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"if (allocated(self%component{cap})) then\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    call self%component{cap}(i)%serializeASCII(verbosityLevel_)\n"
            f"  end do\n"
            f"end if\n"
        )
    content += "call displayUnindent('done',verbosityLevel_)\n"

    function['content'] = content
    _bind(build, 'serializeASCII', function)


def Tree_Node_Serialize_XML(build):
    """Generate `treeNodeSerializeXML`.  Mirrors `Tree_Node_Serialize_XML`."""
    function = {
        'type':        'void',
        'name':        'treeNodeSerializeXML',
        'description': "Serialize tree node content as XML.",
        'modules':     [
            'ISO_Varying_String',
            'Display',
            'String_Handling',
        ],
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
                'variables':  ['fileHandle'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
            {
                'intrinsic':  'character',
                'type':       'len=20',
                'variables':  ['idLabel', 'treeLabel'],
            },
        ],
    }

    content = (
        "!$omp critical(Node_XML_Dump)\n"
        "write (  idLabel,'(i20)') self         %index()\n"
        "write (treeLabel,'(i20)') self%hostTree%index\n"
        "write (fileHandle,'(a,a,a,a,a)') ' <node tree=\"',"
        "trim(adjustl(treeLabel)),'\" id=\"',trim(adjustl(idLabel)),'\" >'\n"
        "write (fileHandle,'(a)') '  <pointer>'\n"
    )
    for pointer in _TREE_NODE_POINTERS:
        content += (
            f"write (idLabel,'(i20)') self%{pointer}%index()\n"
            f"write (fileHandle,'(a,a,a)') '   <{pointer}>',"
            f"trim(adjustl(idLabel)),'</{pointer}>'\n"
        )
    content += "write (fileHandle,'(a)') '  </pointer>'\n"
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"if (allocated(self%component{cap})) then\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    call self%component{cap}(i)%serializeXML(fileHandle)\n"
            f"  end do\n"
            f"end if\n"
        )
    content += (
        "write (fileHandle,*) ' </node>'\n"
        "!$omp end critical(Node_XML_Dump)\n"
    )

    function['content'] = content
    _bind(build, 'serializeXML', function)


def Tree_Node_Serialize_Raw(build):
    """Generate `treeNodeSerializeRaw`.  Mirrors `Tree_Node_Serialize_Raw`."""
    function = {
        'type':        'void',
        'name':        'treeNodeSerializeRaw',
        'description': (
            "Serialize all content of a tree node to a raw (binary) file."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
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
    }

    content = "write (fileHandle) self%isPhysicallyPlausible\n"
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"write (fileHandle) allocated(self%component{cap})\n"
            f"if (allocated(self%component{cap})) then\n"
            f"  select type (component => self%component{cap}(1))\n"
            f"  type is (nodeComponent{cap})\n"
            f"    write (fileHandle) .false.\n"
            f"  class is (nodeComponent{cap})\n"
            f"    write (fileHandle) .true.\n"
            f"  end select\n"
            f"  write (fileHandle) size(self%component{cap})\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    call self%component{cap}(i)%serializeRaw(fileHandle)\n"
            f"  end do\n"
            f"end if\n"
            # Sentinel written after each component class. On restore this is
            # checked to detect (and report, rather than silently corrupting a
            # later read) any desynchronization of the raw stream — e.g. from a
            # component property whose stored size depends on globally-registered
            # state that differs between the storing and restoring runs.
            f"write (fileHandle) {_RAW_SENTINEL}\n"
        )

    function['content'] = content
    _bind(build, 'serializeRaw', function)


def Tree_Node_Deserialize_Raw(build):
    """Generate `treeNodeDeserializeRaw`.  Mirrors `Tree_Node_Deserialize_Raw`."""
    function = {
        'type':        'void',
        'name':        'treeNodeDeserializeRaw',
        'description': "Deserialize a tree node object from a raw (binary) file.",
        'modules':     ['Error'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)', 'target'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['fileHandle'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i', 'componentCount', 'nodeSentinel'],
            },
            {
                'intrinsic':  'logical',
                'variables':  ['isAllocated'],
            },
        ],
    }

    content = "read (fileHandle) self%isPhysicallyPlausible\n"
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"read (fileHandle) isAllocated\n"
            f"if (isAllocated) then\n"
            f"  read (fileHandle) isAllocated\n"
            f"  read (fileHandle) componentCount\n"
            f"  if (allocated(self%component{cap})) "
            f"deallocate(self%component{cap})\n"
            # gfortran bug workaround (array-descriptor `dtype`/`rank` not set):
            # `allocate(class_array(n), mold=/source=<non-polymorphic type
            # scalar>)` sets the descriptor's data/bounds/span but does NOT write
            # its `dtype` (which holds the rank). For these `class(...),
            # allocatable, dimension(:)` component arrays the descriptor then
            # keeps whatever garbage was in the (heap-allocated) node's memory,
            # and when the array is finalized during tree destruction the
            # compiler-generated finalizer reads a bogus (possibly negative)
            # rank, under-allocates an internal temporary (malloc of 1 byte) and
            # overflows it — corrupting the heap and hanging. Allocation forms
            # that DO set `dtype`: a polymorphic `mold=`, a type-spec allocate,
            # and a bare allocate. So: use the polymorphic `mold=default...` for
            # the default-implementation branch, and a bare allocate (which
            # allocates the declared base type) for the base-type/null branches.
            f"  if (isAllocated) then\n"
            f"    allocate(self%component{cap}(componentCount),"
            f"mold=default{cap}Component)\n"
            f"  else\n"
            f"    allocate(self%component{cap}(componentCount))\n"
            f"  end if\n"
            f"  do i=1,componentCount\n"
            f"    self%component{cap}(i)%hostNode => self\n"
            f"  end do\n"
            f"  do i=1,componentCount\n"
            f"    call self%component{cap}(i)%deserializeRaw(fileHandle)\n"
            f"  end do\n"
            f"else\n"
            f"   if (allocated(self%component{cap})) "
            f"deallocate(self%component{cap})\n"
            # Bare allocate of the size-1 "null" placeholder: allocates the
            # declared base type AND sets the descriptor `dtype`/rank (see the
            # gfortran-bug note above; `mold=<base type scalar>` would leave the
            # rank as garbage and hang on finalization during tree destruction).
            f"   allocate(self%component{cap}(1))\n"
            f"end if\n"
            # Verify the sentinel written after this component class. A mismatch
            # means the raw stream desynchronized while reading this component,
            # so fail loudly (with the offending component named) rather than
            # letting a corrupted array descriptor cause an eventual hang in
            # tree destruction.
            f"read (fileHandle) nodeSentinel\n"
            f"if (nodeSentinel /= {_RAW_SENTINEL}) call Error_Report("
            f"'raw tree-node deserialization desynchronized while reading the "
            f"\"{class_dict['name']}\" component (stored file is incompatible "
            f"with the currently-registered component/property configuration)'"
            f"//char(10)//'   file:objects.nodes.components.Inc')\n"
        )

    function['content'] = content
    _bind(build, 'deserializeRaw', function)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _active_classes(build):
    active = set(build.get('componentClassListActive') or [])
    for class_dict in (build.get('componentClasses') or {}).values():
        if class_dict['name'] in active:
            yield class_dict


def _bind(build, method_name, function):
    build.setdefault('types', {}).setdefault('treeNode', {}) \
                                 .setdefault('boundFunctions', []) \
                                 .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       method_name,
    })


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


# ---------------------------------------------------------------------------
# Hook registration.  Order matches Perl Serialization.pm:21-26.
# ---------------------------------------------------------------------------

register('treeNodeSerialization', 'functions', Tree_Node_Serialize_ASCII)
register('treeNodeSerialization', 'functions', Tree_Node_Serialize_XML)
register('treeNodeSerialization', 'functions', Tree_Node_Serialize_Raw)
register('treeNodeSerialization', 'functions', Tree_Node_Deserialize_Raw)
