"""Per-class `<class>Count` and `<class>Get` methods on `treeNode`.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Components/TreeNodes/Classes.pm.  Two
`classIteratedFunctions` hooks — each runs once per component class
during the `Class_Function_Iterator` sub-iteration of the `functions`
phase.
"""



from Galacticus.Build.Components.Utils import register


def Tree_Node_Class_Count(build, class_dict):
    """Generate `treeNode.<class>Count()` returning the number of
    instances of a given component class on a node.

    Mirrors `Tree_Node_Class_Count`.  Components.pm:50 — for inactive
    classes the body is just an `Error_Report` call.
    """
    name = class_dict['name']
    cap  = _ucfirst(name)

    function = {
        'type':        'integer',
        'name':        f"treeNode{cap}Count",
        'description': (
            f"Returns the number of \\mono{{{name}}} components in the node."
        ),
        'modules':     ['Error'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)', 'target'],
                'variables':  ['self'],
            },
        ],
    }

    if name in (build.get('componentClassListActive') or []):
        function['content'] = (
            "select type (self)\n"
            "class is (treeNode)\n"
            f" if (allocated(self%component{cap})) then\n"
            f"   select type (component => self%component{cap}(1))\n"
            f"   type is (nodeComponent{cap})\n"
            f"     treeNode{cap}Count=0\n"
            f"   class default\n"
            f"     treeNode{cap}Count=size(self%component{cap})\n"
            "   end select\n"
            " else\n"
            f"    treeNode{cap}Count=0\n"
            " end if\n"
            "class default\n"
            f" treeNode{cap}Count=0\n"
            " call Error_Report('treeNode of unknown class'//{introspection:location})\n"
            "end select\n"
        )
    else:
        function['content'] = (
            "!$GLC attributes unused :: self\n"
            f"treeNode{cap}Count=0\n"
            "call Error_Report('Galacticus was not compiled with support for "
            "this class'//{introspection:location})\n"
        )

    build.setdefault('types', {}).setdefault('treeNode', {}) \
                                 .setdefault('boundFunctions', []) \
                                 .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       f"{name}Count",
    })


def Tree_Node_Class_Get(build, class_dict):
    """Generate `treeNode.<class>()` returning the component instance,
    optionally auto-creating one if none exists.

    Mirrors `Tree_Node_Class_Get`.  Recursive — the function body calls
    `self%<class>Create()` when `autoCreate` is requested.
    """
    name = class_dict['name']
    cap  = _ucfirst(name)

    function = {
        'type':        f"class(nodeComponent{cap}), pointer => component",
        'name':        f"treeNode{cap}Get",
        'description': (
            f"Return a \\mono{{{name}}} component member of the node. "
            "If no \\mono{instance} is specified, return the first instance. "
            "If \\mono{autoCreate} is \\mono{true} then create a single "
            "instance of the component if none exists in the node."
        ),
        'recursive':   True,
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
                'attributes': ['intent(in   )', 'optional'],
                'variables':  ['instance'],
            },
            {
                'intrinsic':  'logical',
                'attributes': ['intent(in   )', 'optional'],
                'variables':  ['autoCreate'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['instance_'],
            },
            {
                'intrinsic':  'logical',
                'variables':  ['autoCreate_'],
            },
        ],
    }

    if name in (build.get('componentClassListActive') or []):
        function['content'] = (
            f"if (.not.present(autoCreate).and..not.present(instance).and.allocated(self%component{cap})) then\n"
            "    ! Rapid path for the most common call pattern - component exists, no auto-creation is\n"
            "    ! allowed, and no instance is specified (so the first instance will be used).\n"
            f"   component => self%component{cap}(1)\n"
            "else\n"
            "   ! Set default instance and auto-creation options.\n"
            "   instance_=1\n"
            "   if (present(instance)) instance_=instance\n"
            "   autoCreate_=.false.\n"
            "   if (present(autoCreate)) autoCreate_=autoCreate\n"
            "   ! If we are allowed to auto-create the component and it has generic type then deallocate it\n"
            "   ! to force it to be created later.\n"
            f"   if (autoCreate_.and.allocated(self%component{cap}).and.same_type_as(self%component{cap}(1),{cap}Class)) &\n"
            f"        & deallocate(self%component{cap})\n"
            "   ! If the component is not yet allocated, allocate it now.\n"
            f"   if (.not.allocated(self%component{cap})) then\n"
            "      if (autoCreate_) then\n"
            "         ! Auto-creation is allowed - create the component.\n"
            f"         call self%{name}Create()\n"
            "      else\n"
            "         ! Auto-creation is not allowed - create a generic type.\n"
            f"         allocate(self%component{cap}(1))\n"
            f"         self%component{cap}(1)%hostNode => self\n"
            "      end if\n"
            "   end if\n"
            "   ! Return the requested component.\n"
            f"   component => self%component{cap}(instance_)\n"
            "end if\n"
        )
    else:
        function['content'] = (
            "!$GLC attributes unused :: self, instance, instance_\n"
            "autoCreate_=.false.\n"
            "if (present(autoCreate)) autoCreate_=autoCreate\n"
            "if (autoCreate_) then\n"
            " ! Support for this component was not compiled, so we can not "
            "create it.\n"
            " component => null()\n"
            " call Error_Report('Galacticus was compiled without support for "
            "this class'//{introspection:location})\n"
            "else\n"
            " ! Support for this component was not compiled, return the "
            "default of the class - and trust that the user knows what they "
            "are doing.\n"
            f" component => default{cap}Component\n"
            "end if\n"
        )

    build.setdefault('types', {}).setdefault('treeNode', {}) \
                                 .setdefault('boundFunctions', []) \
                                 .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       name,
    })


register('treeNodesClasses', 'classIteratedFunctions', Tree_Node_Class_Count)
register('treeNodesClasses', 'classIteratedFunctions', Tree_Node_Class_Get)


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text
