"""`treeNode` lifecycle methods: initialize, build-from-XML, destroy, and
per-class create/destroy.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Components/TreeNodes/CreateDestroy.pm.
Three `functions`-phase hooks emit the global lifecycle methods, plus
two `classIteratedFunctions` hooks emit per-class create/destroy
methods.
"""



from Galacticus.Build.Components.Utils import register


def Tree_Node_Creation(build):
    """Generate `treeNodeInitialize`.

    Mirrors `Tree_Node_Creation`.  Initialises pointers, allocates one
    instance per active component class, then sets index, unique ID,
    timestep, subsampling weight, and physical-state flags.
    """
    pointers = ('parent', 'firstChild', 'sibling', 'firstSatellite',
                'mergeTarget', 'firstMergee', 'siblingMergee',
                'formationNode', 'event')
    nullify_block = ''.join(f"nullify (self%{p})\n" for p in pointers)

    active = build.get('componentClassListActive') or []
    allocate_block = ''.join(
        f"allocate(self%component{c}(1))\n" for c in active
    )
    host_block = ''.join(
        f"   self%component{c}(1)%hostNode => self\n" for c in active
    )

    content = (
        "! Ensure pointers are nullified.\n"
        f"{nullify_block}"
        f"{allocate_block}"
        "select type (self)\n"
        "type is (treeNode)\n"
        f"{host_block}"
        "end select\n"
        "! Assign a host tree if supplied.\n"
        "if (present(hostTree)) then\n"
        "   self%hostTree => hostTree\n"
        "else\n"
        "   self%hostTree => null()\n"
        "end if\n"
        "! Assign index if supplied.\n"
        "if (present(index)) then\n"
        " call self%indexSet(index)\n"
        "else\n"
        " call self%indexSet(-1_kind_int8)\n"
        "end if\n"
        "! Assign a unique ID.\n"
        "call self%uniqueIDSet()\n"
        "! Assign a timestep and subsampling weight.\n"
        "self%timeStepValue         =-1.0d0\n"
        "self%subsamplingWeightValue= 1.0d0\n"
        "! Initialize physical state.\n"
        "self%isSolvable            =.true.\n"
        "self%isPhysicallyPlausible =.true.\n"
    )

    function = {
        'type':        'void',
        'name':        'treeNodeInitialize',
        'description': r"Initialize a \mono{treeNode} object.",
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
                'type':       'kind=kind_int8',
                'attributes': ['intent(in   )', 'optional'],
                'variables':  ['index'],
            },
            {
                'intrinsic':  'type',
                'type':       'mergerTree',
                'attributes': ['intent(in   )', 'optional', 'target'],
                'variables':  ['hostTree'],
            },
        ],
        'content':     content,
    }
    _bind(build, 'initialize', function)


def Tree_Node_Builder(build):
    """Generate `treeNodeComponentBuilder` — populate components from a
    FoX_DOM XML node.

    Mirrors `Tree_Node_Builder`.
    """
    active = set(build.get('componentClassListActive') or [])
    component_classes = list(
        (build.get('componentClasses') or {}).values()
    )

    head = (
        "select type (self)\n"
        "type is (treeNode)\n"
        "   call componentIndex%initialize()\n"
        "   !$omp critical (FoX_DOM_Access)\n"
    )

    counts = ""
    for component in component_classes:
        if component['name'] not in active:
            continue
        cap   = _ucfirst(component['name'])
        cname = component['name']
        counts += (
            f"    componentList => getChildNodes(nodeDefinition)\n"
            f"    componentCount=0\n"
            f"    do i=0,getLength(componentList)-1\n"
            f"      componentDefinition => item(componentList,i)\n"
            f"      if (getNodeName(componentDefinition) == '{cname}') "
            f"componentCount=componentCount+1\n"
            f"    end do\n"
            f"    if (componentCount > 0) then\n"
            f"      if (allocated(self%component{cap})) "
            f"deallocate(self%component{cap})\n"
            f"      allocate(self%component{cap}(componentCount),"
            f"source=default{cap}Component)\n"
            f"      call componentIndex%set('{cname}',0)\n"
            f"    end if\n"
        )

    middle = (
        "   componentCount=getLength(componentList)\n"
        "   !$omp end critical (FoX_DOM_Access)\n"
        "   do i=0,componentCount-1\n"
    )

    body = ""
    for component in component_classes:
        if component['name'] not in active:
            continue
        cap   = _ucfirst(component['name'])
        cname = component['name']
        body += (
            f"     !$omp critical (FoX_DOM_Access)\n"
            f"     componentDefinition => item(componentList,i)\n"
            f"     nodeName=getNodeName(componentDefinition)\n"
            f"     !$omp end critical (FoX_DOM_Access)\n"
            f"     if (trim(nodeName) == '{cname}') then\n"
            f"       j=componentIndex%value('{cname}')\n"
            f"       j=j+1\n"
            f"       self%component{cap}(j)%hostNode => self\n"
            f"       call self%component{cap}(j)%builder(componentDefinition)\n"
            f"       call componentIndex%set('{cname}',j)\n"
            f"     end if\n"
        )

    tail = (
        "  end do\n"
        "  call componentIndex%destroy()\n"
        "end select\n"
    )

    function = {
        'type':        'void',
        'name':        'treeNodeComponentBuilder',
        'description': (
            r"Build components in a \mono{treeNode} object given an XML "
            r"definition."
        ),
        'modules':     [
            ("FoX_DOM, only : node, nodeList, getChildNodes, getLength, "
             "getNodeName, item"),
            'Hashes',
        ],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'type',
                'type':       'node',
                'attributes': ['intent(in   )', 'pointer'],
                'variables':  ['nodeDefinition'],
            },
            {
                'intrinsic':  'type',
                'type':       'node',
                'attributes': ['pointer'],
                'variables':  ['componentDefinition'],
            },
            {
                'intrinsic':  'type',
                'type':       'nodeList',
                'attributes': ['pointer'],
                'variables':  ['componentList'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i', 'j', 'componentCount'],
            },
            {
                'intrinsic':  'type',
                'type':       'integerHash',
                'variables':  ['componentIndex'],
            },
            {
                'intrinsic':  'character',
                'type':       'len=128',
                'variables':  ['nodeName'],
            },
        ],
        'content':     head + counts + middle + body + tail,
    }
    _bind(build, 'componentBuilder', function)


def Tree_Node_Finalization(build):
    """Generate `treeNodeDestroy`.

    Mirrors `Tree_Node_Finalization`.  Destroys all per-class component
    arrays, the formation node, and walks the linked list of attached
    events to free each (and its paired event on a partner node).
    """
    active = build.get('componentClassListActive') or []
    destroy_block = ''.join(
        f"call self%{c}Destroy()\n" for c in active
    )
    if destroy_block:
        # Match Perl `join(" ", ...)` — a single space between entries
        # rather than newlines.  The `\n` per line is already there.
        destroy_block = " ".join(destroy_block.splitlines(keepends=True))

    content = (
        "! Destroy all components.\n"
        f"{destroy_block}"
        "! Destroy any formation node.\n"
        "if (associated(self%formationNode)) then\n"
        "   call self%formationNode%destroy()\n"
        "   deallocate(self%formationNode)\n"
        "   nullify(self%formationNode)\n"
        "end if\n"
        "! Remove any events attached to the node, along with their paired "
        "event in other nodes.\n"
        "thisEvent => self%event\n"
        "do while (associated(thisEvent))\n"
        "   ! If a paired node is given, remove any paired event from it.\n"
        "   if (associated(thisEvent%node)) then\n"
        "      ! Locate the paired event and remove it.\n"
        "      pairEvent => thisEvent%node%event\n"
        "      lastEvent => thisEvent%node%event\n"
        "      ! Iterate over all events.\n"
        "      pairMatched=.false.\n"
        "      do while (associated(pairEvent).and..not.pairMatched)\n"
        "         ! Match the paired event ID with the current event ID.\n"
        "         if (pairEvent%ID == thisEvent%ID) then\n"
        "            pairMatched=.true.\n"
        "            if (associated(pairEvent,thisEvent%node%event)) then\n"
        "               thisEvent%node  %event => pairEvent%next\n"
        "               lastEvent       => thisEvent%node %event\n"
        "            else\n"
        "               lastEvent%next  => pairEvent%next\n"
        "            end if\n"
        "            nextEvent => pairEvent%next\n"
        "            deallocate(pairEvent)\n"
        "            pairEvent => nextEvent\n"
        "         else\n"
        "            lastEvent => pairEvent\n"
        "            pairEvent => pairEvent%next\n"
        "         end if\n"
        "      end do \n"
        "      if (.not.pairMatched) call Error_Report('unable to find paired"
        " event'//{introspection:location})\n"
        "   end if\n"
        "   nextEvent => thisEvent%next\n"
        "   deallocate(thisEvent)\n"
        "   thisEvent => nextEvent\n"
        "end do\n"
    )

    function = {
        'type':        'void',
        'name':        'treeNodeDestroy',
        'description': r"Destroy a \mono{treeNode} object.",
        'modules':     ['Error'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'class',
                'type':       'nodeEvent',
                'attributes': ['pointer'],
                'variables':  ['thisEvent', 'pairEvent', 'lastEvent', 'nextEvent'],
            },
            {
                'intrinsic':  'logical',
                'variables':  ['pairMatched'],
            },
        ],
        'content':     content,
    }
    _bind(build, 'destroy', function)


def Tree_Node_Class_Creation(build, class_dict):
    """Generate `nodeComponent<Class>Create` — adds one instance of the
    class on `self`.  Skipped entirely for inactive classes (matches Perl).

    Mirrors `Tree_Node_Class_Creation`.
    """
    name = class_dict['name']
    if name not in (build.get('componentClassListActive') or []):
        return
    cap = _ucfirst(name)

    non_null_components = [
        m['name'] for m in class_dict.get('members') or []
        if m['name'] != 'null'
    ]
    non_null_components.sort()
    if non_null_components:
        non_null_block = (
            "char(10)//"
            + "//char(10)//".join(f"'   {n}'" for n in non_null_components)
        )
    else:
        non_null_block = "char(10)"

    content = (
        "if (displayVerbosity() >= verbosityLevelInfo) then\n"
        "  block\n"
        "    type(varying_string) :: message\n"
        f"    message='Creating {name} in node '\n"
        "    message=message//self%index()\n"
        "    call displayMessage(message,verbosityLevelInfo)\n"
        "  end block\n"
        "end if\n"
        "if (present(template)) then\n"
        f"   allocate(self%component{cap}(1),source=template)\n"
        "else\n"
        f"   select type (default{cap}Component)\n"
        f"   type is (nodeComponent{cap}Null)\n"
        "      block\n"
        "         type(varying_string) :: message\n"
        f"         message=         'creation of the {name} component "
        f"requested, but that component is null'//char(10)\n"
        f"         message=message//'please select a non-null {name} "
        f"component - available options are:'\n"
        f"         message=message//{non_null_block}\n"
        "         call displayMessage(message,verbosityLevelSilent)\n"
        "         call Error_Report('refusing to create null instance'"
        "//{introspection:location})\n"
        "      end block\n"
        "   class default\n"
        f"      allocate(self%component{cap}(1),source=default{cap}Component)\n"
        "   end select\n"
        "end if\n"
        "select type (self)\n"
        "type is (treeNode)\n"
        f"  do i=1,size(self%component{cap})\n"
        f"    self%component{cap}(i)%hostNode => self\n"
        f"    call self%component{cap}(i)%initialize()\n"
        "  end do\n"
        "end select\n"
    )

    function = {
        'type':        'void',
        'name':        f"nodeComponent{cap}Create",
        'description': (
            f"Create the \\mono{{{name}}} component of \\mono{{self}}."
        ),
        'modules':     [
            'ISO_Varying_String',
            'Display',
            'Error',
            'String_Handling',
        ],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['target', 'intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'class',
                'type':       f"nodeComponent{cap}",
                'attributes': ['intent(in   )', 'optional'],
                'variables':  ['template'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
        ],
        'content':     content,
    }
    _bind(build, f"{name}Create", function)


def Tree_Node_Class_Destruction(build, class_dict):
    """Generate `nodeComponent<Class>Destroy` — frees the per-class array.

    Mirrors `Tree_Node_Class_Destruction`.  Inactive classes get a stub
    that raises an `Error_Report`.
    """
    name = class_dict['name']
    cap  = _ucfirst(name)

    function = {
        'type':        'void',
        'name':        f"nodeComponent{cap}Destroy",
        'description': (
            f"Destroy the \\mono{{{name}}} component of \\mono{{self}}"
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
        ],
    }

    if name in (build.get('componentClassListActive') or []):
        function['content'] = (
            f"if (allocated(self%component{cap})) then\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    call self%component{cap}(i)%destroy()\n"
            "  end do\n"
            f"  deallocate (self%component{cap})\n"
            "end if\n"
        )
    else:
        function['modules'] = ['Error']
        function['content'] = (
            f"call Error_Report(\"component '{name}' is not active\""
            "//{introspection:location})\n"
        )

    _bind(build, f"{name}Destroy", function)


def _bind(build, method_name, function):
    """Append a `procedure :: <method_name> => <function>` binding to
    `treeNode`.
    """
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
# Hook registration
# ---------------------------------------------------------------------------

register('treeNodesCreateDestroy', 'functions', Tree_Node_Creation)
register('treeNodesCreateDestroy', 'functions', Tree_Node_Builder)
register('treeNodesCreateDestroy', 'functions', Tree_Node_Finalization)
register('treeNodesCreateDestroy', 'classIteratedFunctions', Tree_Node_Class_Creation)
register('treeNodesCreateDestroy', 'classIteratedFunctions', Tree_Node_Class_Destruction)
