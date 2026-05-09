"""Per-class create/destroy/builder + meta-property add/count/name hooks.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Components/Classes/CreateDestroy.pm.
Seven `classIteratedFunctions` hooks: Initialization, Builder,
Finalization (small stubs); Class_Create_By_Interrupt (top-level
helper for createIfNeeded); plus Add/Count/Name_Meta_Property (one
function per meta-property type per class).
"""



from Galacticus.Build.Components.Utils                  import (
    register,
    _component_properties,
)
from Galacticus.Build.Components.Classes.MetaProperties import meta_property_types


def Class_Initialization(build, class_dict):
    """Generate `nodeComponent<Class>Initialize` — abstract stub that
    errors when called.  Mirrors `Class_Initialization`.
    """
    name      = class_dict['name']
    type_name = 'nodeComponent' + _ucfirst(name)
    function  = {
        'type':        'void',
        'name':        type_name + 'Initialize',
        'description': f"Initialize a generic \\mono{{{name}}} component.",
        'modules':     ['Error'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       type_name,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
        ],
        'content': (
            "!$GLC attributes unused :: self\n"
            "\n"
            "call Error_Report('can not initialize a generic component'"
            "//{introspection:location})\n"
        ),
    }
    _bind(build, type_name, function, 'initialize')


def Class_Finalization(build, class_dict):
    """Generate `nodeComponent<Class>Finalize` — no-op stub.
    Mirrors `Class_Finalization`.
    """
    name      = class_dict['name']
    type_name = 'nodeComponent' + _ucfirst(name)
    function  = {
        'type':        'void',
        'name':        type_name + 'Finalize',
        'description': f"Finalize a generic \\mono{{{name}}} component.",
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       type_name,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
        ],
        'content': (
            "!$GLC attributes unused :: self\n"
            "\n"
            "! Nothing to do.\n"
        ),
    }
    _bind(build, type_name, function, 'destroy')


def Class_Builder(build, class_dict):
    """Generate `nodeComponent<Class>Builder` — abstract stub.
    Mirrors `Class_Builder`.
    """
    name      = class_dict['name']
    type_name = 'nodeComponent' + _ucfirst(name)
    function  = {
        'type':        'void',
        'name':        type_name + 'Builder',
        'description': (
            f"Build a generic \\mono{{{name}}} component from a supplied "
            "XML definition."
        ),
        'modules':     ['Error', 'FoX_DOM, only : node'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       type_name,
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'type',
                'type':       'node',
                'attributes': ['intent(in   )', 'pointer'],
                'variables':  ['componentDefinition'],
            },
        ],
        'content': (
            "!$GLC attributes unused :: self, componentDefinition\n"
            "\n"
            "call Error_Report('can not build a generic component'"
            "//{introspection:location})\n"
        ),
    }
    _bind(build, type_name, function, 'builder')


def Class_Create_By_Interrupt(build, class_dict):
    """Generate `<class>CreateByInterrupt` (a free function, not bound)
    when at least one member has at least one property with the
    `createIfNeeded` attribute.  Mirrors `Class_Create_By_Interrupt`.
    """
    name = class_dict['name']

    # Skip unless any property requests createIfNeeded.
    function_required = False
    for member in class_dict.get('members') or []:
        for prop in _component_properties(member):
            if (prop.get('attributes') or {}).get('createIfNeeded'):
                function_required = True
                break
        if function_required:
            break
    if not function_required:
        return

    cap       = _ucfirst(name)
    type_name = 'nodeComponent' + cap

    content = f"{name} => self%{name}(autoCreate=.true.)\n"

    # Iterate through members carrying explicit createFunction overrides.
    members_with_create = [
        m for m in (class_dict.get('members') or [])
        if 'createFunction' in m
    ]
    if members_with_create:
        content += f"select type ({name})\n"
        for member in members_with_create:
            create_function = _resolve_create_function(name, member)
            impl_type = type_name + _ucfirst(member['name'])
            content += (
                f"type is ({impl_type})\n"
                f"   call {create_function}({name},timeEnd)\n"
            )
        content += "end select\n"

    function = {
        'type':        'void',
        'name':        name + 'CreateByInterrupt',
        'description': (
            f"Create the \\mono{{{name}}} component of \\mono{{self}} via "
            "an interrupt."
        ),
        'variables':   [
            {
                'intrinsic':  'type',
                'type':       'treeNode',
                'attributes': ['target', 'intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['intent(in   )', 'optional'],
                'variables':  ['timeEnd'],
            },
            {
                'intrinsic':  'class',
                'type':       type_name,
                'attributes': ['pointer'],
                'variables':  [name],
            },
        ],
        'content':     content,
    }
    build.setdefault('functions', []).append(function)


def _resolve_create_function(class_name, member):
    """Return the Fortran symbol to call from a `createIfNeeded`
    interrupt branch for `member`.  Mirrors the conditional at
    Classes/CreateDestroy.pm:201-213.
    """
    cf = member['createFunction']
    if isinstance(cf, dict):
        if cf.get('isDeferred'):
            return (
                class_name
                + _ucfirst(member['name'])
                + 'CreateFunction'
            )
        return cf.get('content', cf)
    return cf


def Class_Add_Meta_Property(build, class_dict):
    """Generate `nodeComponent<Class>Add<Label>Rank<N>MetaProperty` per
    meta-property type.  Mirrors `Class_Add_Meta_Property`.

    Body uses one `!$omp critical` to serialize meta-property updates
    and grows the class's bookkeeping arrays in lock-step.  The
    `float`+rank-0 form gets an extra evolvable-tracking branch.
    """
    name      = class_dict['name']
    cap       = _ucfirst(name)
    type_name = 'nodeComponent' + cap
    active    = name in (build.get('componentClassListActive') or [])

    for mpt in meta_property_types:
        cap_label = _ucfirst(mpt['label'])
        rank      = mpt['rank']
        is_evolvable = (mpt['label'] == 'float' and rank == 0)
        prefix    = f"{cap_label}Rank{rank}"
        var_pref  = f"{name}{prefix}"
        fn_id     = f"{type_name}Add{prefix}MetaProperty"

        options = ['isCreator']
        tmps    = ['creatorTmp']
        if is_evolvable:
            options.insert(0, 'isEvolvable')
            tmps.insert(0, 'evolvableTmp')

        function = {
            'type':        'integer',
            'name':        fn_id,
            'description': (
                f"Add a rank-{rank} {mpt['label']}meta-property to the "
                f"generic \\mono{{{name}}} component."
            ),
            'modules':     ['ISO_Varying_String', 'Error'],
            'variables':   [
                {
                    'intrinsic':  'class',
                    'type':       type_name,
                    'attributes': ['intent(inout)'],
                    'variables':  ['self'],
                },
                {
                    'intrinsic':  'type',
                    'type':       'varying_string',
                    'attributes': ['intent(in   )'],
                    'variables':  ['label'],
                },
                {
                    'intrinsic':  'character',
                    'type':       'len=*',
                    'attributes': ['intent(in   )'],
                    'variables':  ['name'],
                },
                {
                    'intrinsic':  'logical',
                    'attributes': ['intent(in   )', 'optional'],
                    'variables':  options,
                },
                {
                    'intrinsic':  'type',
                    'type':       'varying_string',
                    'attributes': ['allocatable', 'dimension(:)'],
                    'variables':  ['labelsTmp', 'namesTmp'],
                },
                {
                    'intrinsic':  'logical',
                    'attributes': ['allocatable', 'dimension(:)'],
                    'variables':  tmps,
                },
                {
                    'intrinsic':  'logical',
                    'variables':  ['found'],
                },
            ],
        }

        if active:
            content  = (
                "!$GLC attributes unused :: self\n"
                "\n"
                f"!$omp critical ({var_pref}MetaPropertyUpdate)\n"
                "found=.false.\n"
                f"if (allocated({var_pref}MetaPropertyLabels)) then\n"
                f" do {fn_id}=1,size({var_pref}MetaPropertyLabels)\n"
                f"  if ({var_pref}MetaPropertyLabels({fn_id}) == label) then\n"
                f"   found=.true.\n"
                f"   exit\n"
                f"  end if\n"
                f" end do\n"
                f" if (.not.found) then\n"
                f"  call move_alloc({var_pref}MetaPropertyLabels , labelsTmp)\n"
                f"  call move_alloc({var_pref}MetaPropertyNames  ,  namesTmp)\n"
                f"  call move_alloc({var_pref}MetaPropertyCreator,creatorTmp)\n"
                f"  allocate({var_pref}MetaPropertyLabels (size( labelsTmp)+1))\n"
                f"  allocate({var_pref}MetaPropertyNames  (size(  namesTmp)+1))\n"
                f"  allocate({var_pref}MetaPropertyCreator(size(creatorTmp)+1))\n"
                f"  {var_pref}MetaPropertyLabels (1:size( labelsTmp))= labelsTmp\n"
                f"  {var_pref}MetaPropertyNames  (1:size(  namesTmp))=  namesTmp\n"
                f"  {var_pref}MetaPropertyCreator(1:size(creatorTmp))=creatorTmp\n"
                f"  deallocate( labelsTmp)\n"
                f"  deallocate(  namesTmp)\n"
                f"  deallocate(creatorTmp)\n"
            )
            if is_evolvable:
                content += (
                    f"  call move_alloc({var_pref}MetaPropertyEvolvable,evolvableTmp)\n"
                    f"  allocate({var_pref}MetaPropertyEvolvable(size(evolvableTmp)+1))\n"
                    f"  {var_pref}MetaPropertyEvolvable(1:size(evolvableTmp))=evolvableTmp\n"
                    f"  deallocate(evolvableTmp)\n"
                )
            content += (
                " end if\n"
                "else\n"
                f" allocate({var_pref}MetaPropertyLabels (1))\n"
                f" allocate({var_pref}MetaPropertyNames  (1))\n"
                f" allocate({var_pref}MetaPropertyCreator(1))\n"
            )
            if is_evolvable:
                content += f" allocate({var_pref}MetaPropertyEvolvable(1))\n"
            content += "end if\n"
            content += (
                "if (.not.found) then\n"
                f" {fn_id}=size({var_pref}MetaPropertyLabels)\n"
                f" {var_pref}MetaPropertyLabels   ({fn_id})=label\n"
                f" {var_pref}MetaPropertyNames    ({fn_id})=name\n"
            )
            if is_evolvable:
                content += (
                    " if (present(isEvolvable)) then\n"
                    f"  {var_pref}MetaPropertyEvolvable({fn_id})=isEvolvable\n"
                    " else\n"
                    f"  {var_pref}MetaPropertyEvolvable({fn_id})=.false.\n"
                    " end if\n"
                )
            content += (
                " if (present(isCreator  )) then\n"
                f"  {var_pref}MetaPropertyCreator  ({fn_id})=isCreator\n"
                " else\n"
                f"  {var_pref}MetaPropertyCreator  ({fn_id})=.false.\n"
                " end if\n"
                f" {var_pref}MetaPropertyCount={var_pref}MetaPropertyCount+1\n"
            )
            if is_evolvable:
                content += (
                    f" if ({var_pref}MetaPropertyEvolvable({fn_id})) "
                    f"{var_pref}MetaPropertyEvolvableCount="
                    f"{var_pref}MetaPropertyEvolvableCount+1\n"
                )
            content += (
                " propertyNameLengthMax=max(len(name),propertyNameLengthMax) \n"
                "else\n"
                " if (present(isCreator)) then\n"
                f"  if (isCreator) {var_pref}MetaPropertyCreator({fn_id})=.true.\n"
                " end if\n"
            )
            if is_evolvable:
                content += (
                    " if (present(isEvolvable)) then\n"
                    f"  if ({var_pref}MetaPropertyEvolvable({fn_id}) .neqv. isEvolvable) "
                    "call Error_Report('inconsistent evolvability for meta-property'"
                    "//{introspection:location})\n"
                    " else\n"
                    f"  if ({var_pref}MetaPropertyEvolvable({fn_id})                   ) "
                    "call Error_Report('inconsistent evolvability for meta-property'"
                    "//{introspection:location})\n"
                    " end if\n"
                )
            content += (
                "end if\n"
                f"!$omp end critical ({var_pref}MetaPropertyUpdate)\n"
            )
            function['content'] = content
        # Inactive class → no body; the function returns whatever the
        # implicit value of the result is.  Match Perl, which simply
        # leaves $function->{content} unset.

        _bind(build, type_name, function, f"add{prefix}MetaProperty")


def Class_Count_Meta_Property(build, class_dict):
    """Generate `component<Class>Count<Label>Rank<N>MetaProperties`.
    Mirrors `Class_Count_Meta_Property`.
    """
    name      = class_dict['name']
    cap       = _ucfirst(name)
    type_name = 'nodeComponent' + cap
    active    = name in (build.get('componentClassListActive') or [])

    for mpt in meta_property_types:
        cap_label = _ucfirst(mpt['label'])
        rank      = mpt['rank']
        prefix    = f"{cap_label}Rank{rank}"
        var_pref  = f"{name}{prefix}"

        function = {
            'type':        'integer => countMetaProperties',
            'name':        f"component{cap}Count{prefix}MetaProperties",
            'description': (
                f"Return the number of rank-{rank} {mpt['label']}meta-properties "
                f"associated with the generic \\mono{{{name}}} component."
            ),
            'variables':   [
                {
                    'intrinsic':  'class',
                    'type':       type_name,
                    'attributes': ['intent(inout)'],
                    'variables':  ['self'],
                },
            ],
        }
        if active:
            function['content'] = (
                "!$GLC attributes unused :: self\n"
                "\n"
                f"!$omp critical ({var_pref}MetaPropertyUpdate)\n"
                f"if (allocated({var_pref}MetaPropertyNames)) then\n"
                f" countMetaProperties=size({var_pref}MetaPropertyNames)\n"
                "else\n"
                " countMetaProperties=0\n"
                "end if\n"
                f"!$omp end critical ({var_pref}MetaPropertyUpdate)\n"
            )
        _bind(build, type_name, function, f"count{prefix}MetaProperties")


def Class_Name_Meta_Property(build, class_dict):
    """Generate `component<Class>Name<Label>Rank<N>MetaProperty`.
    Mirrors `Class_Name_Meta_Property`.
    """
    name      = class_dict['name']
    cap       = _ucfirst(name)
    type_name = 'nodeComponent' + cap
    active    = name in (build.get('componentClassListActive') or [])

    for mpt in meta_property_types:
        cap_label = _ucfirst(mpt['label'])
        rank      = mpt['rank']
        prefix    = f"{cap_label}Rank{rank}"
        var_pref  = f"{name}{prefix}"

        function = {
            'type':        'type(varying_string) => nameMetaProperty',
            'name':        f"component{cap}Name{prefix}MetaProperty",
            'description': (
                f"Return the name of the indexed of rank-{rank} {mpt['label']} "
                f"meta-property associated with the generic \\mono{{{name}}} "
                "component."
            ),
            'modules':     ['ISO_Varying_String'],
            'variables':   [
                {
                    'intrinsic':  'class',
                    'type':       type_name,
                    'attributes': ['intent(inout)'],
                    'variables':  ['self'],
                },
                {
                    'intrinsic':  'integer',
                    'attributes': ['intent(in   )'],
                    'variables':  ['index'],
                },
            ],
        }
        if active:
            function['modules'].append('Error')
            function['content'] = (
                "!$GLC attributes unused :: self\n"
                "\n"
                f"!$omp critical ({var_pref}MetaPropertyUpdate)\n"
                f"if (index > 0 .and. index <= size({var_pref}MetaPropertyNames)) then\n"
                f" nameMetaProperty={var_pref}MetaPropertyNames(index)\n"
                "else\n"
                " nameMetaProperty=var_str('')\n"
                " call Error_Report('meta-property index is out of range'"
                "//{introspection:location})\n"
                "end if\n"
                f"!$omp end critical ({var_pref}MetaPropertyUpdate)\n"
            )
        _bind(build, type_name, function, f"name{prefix}MetaProperty")


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
# Hook registration.  Order matches Perl Classes/CreateDestroy.pm:23-30.
# ---------------------------------------------------------------------------

register('classesCreateDestroy', 'classIteratedFunctions', Class_Initialization)
register('classesCreateDestroy', 'classIteratedFunctions', Class_Builder)
register('classesCreateDestroy', 'classIteratedFunctions', Class_Finalization)
register('classesCreateDestroy', 'classIteratedFunctions', Class_Create_By_Interrupt)
register('classesCreateDestroy', 'classIteratedFunctions', Class_Add_Meta_Property)
register('classesCreateDestroy', 'classIteratedFunctions', Class_Count_Meta_Property)
register('classesCreateDestroy', 'classIteratedFunctions', Class_Name_Meta_Property)
