# `treeNode` ODE-solver methods: serialize/deserialize count, values,
# rates, scales, inactive/non-negative flags, plus the offsets builder.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/TreeNodes/ODESolver.pm.  Eleven
# `functions`-phase hooks (one of which — `Tree_Node_ODE_Step_Initialize`
# — itself emits four bound methods, one per quantity).



from Galacticus.Build.Components.Utils import register


# ---------------------------------------------------------------------------
# Tree_Node_ODE_Step_Initialize — emits four `odeStep<X>sInitialize` methods
# ---------------------------------------------------------------------------

_STEP_INIT_QUANTITIES = (
    {'name': 'rate',     'value': '0.0d0'},
    {'name': 'scale',    'value': '1.0d0'},
    {'name': 'inactive', 'value': '.false.'},
    {'name': 'analytic', 'value': '.false.'},
)


def Tree_Node_ODE_Step_Initialize(build):
    """Emit four `odeStep<X>sInitialize` methods (rates / scales /
    inactives / analytics).  Mirrors `Tree_Node_ODE_Step_Initialize`.
    """
    for q in _STEP_INIT_QUANTITIES:
        cap = _ucfirst(q['name'])
        function = {
            'type':        'void',
            'name':        f"treeNodeODEStep{cap}sInitialize",
            'description': (
                f"Initialize the {q['name']}s in components of tree node "
                r"\mono{self} in preparation for an ODE solver step."
            ),
            'variables':   [
                {
                    'intrinsic':  'class',
                    'type':       'treeNode',
                    'attributes': ['intent(in   )'],
                    'variables':  ['self'],
                },
            ],
            'content':     (
                "!$GLC attributes unused :: self\n"
                f"node{cap}s={q['value']}\n"
            ),
        }
        _bind(build, f"odeStep{cap}sInitialize", function)


def Tree_Node_ODE_Serialize_Count(build):
    """Generate `treeNodeSerializeCount`.  Mirrors `Tree_Node_ODE_Serialize_Count`."""
    function = {
        'type':        'integer',
        'name':        'treeNodeSerializeCount',
        'description': "Return a count of the size of the node when serialized to an array.",
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
                'variables':  ['propertyType'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
        ],
    }
    content = "treeNodeSerializeCount=0\n"
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"if (allocated(self%component{cap})) then\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    treeNodeSerializeCount=treeNodeSerializeCount"
            f"+self%component{cap}(i)%serializeCount(propertyType)\n"
            f"  end do\n"
            f"end if\n"
        )
    function['content'] = content
    _bind(build, 'serializeCount', function)


def Tree_Node_ODE_Serialize_Values(build):
    """Generate `treeNodeSerializeValuesToArray`.  Mirrors `Tree_Node_ODE_Serialize_Values`."""
    function = {
        'type':        'void',
        'name':        'treeNodeSerializeValuesToArray',
        'description': "Serialize evolvable properties of a node into an array.",
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['dimension(:)', 'intent(  out)'],
                'variables':  ['array'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['propertyType'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['count', 'offset', 'i'],
            },
        ],
    }
    content = "offset=1\n"
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"if (allocated(self%component{cap})) then\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    count=self%component{cap}(i)%serializeCount(propertyType)\n"
            f"    if (count > 0) call self%component{cap}(i)%serializeValues"
            f"(array(offset:),propertyType)\n"
            f"    offset=offset+count\n"
            f"  end do\n"
            f"end if\n"
        )
    function['content'] = content
    _bind(build, 'serializeValues', function)


def Tree_Node_ODE_Deserialize_Values(build):
    """Generate `treeNodeDeserializeValuesFromArray`.  Mirrors `Tree_Node_ODE_Deserialize_Values`."""
    function = {
        'type':        'void',
        'name':        'treeNodeDeserializeValuesFromArray',
        'description': "Deserialize evolvable properties of a node from an array.",
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['dimension(:)', 'intent(in   )'],
                'variables':  ['array'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['propertyType'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['count', 'offset', 'i'],
            },
        ],
    }
    content = "offset=1\n"
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"if (allocated(self%component{cap})) then\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    count=self%component{cap}(i)%serializeCount(propertyType)\n"
            f"    if (count > 0) then\n"
            f"       call self%component{cap}(i)%deserializeValues"
            f"(array(offset:),propertyType)\n"
            f"       offset=offset+count\n"
            f"    end if\n"
            f"  end do\n"
            f"end if\n"
        )
    function['content'] = content
    _bind(build, 'deserializeValues', function)


def Tree_Node_ODE_Serialize_NonNegative(build):
    """Generate `treeNodeSerializeNonNegativeToArray`.  Mirrors
    `Tree_Node_ODE_Serialize_NonNegative`.
    """
    function = {
        'type':        'void',
        'name':        'treeNodeSerializeNonNegativeToArray',
        'description': (
            "Serialize non-negative status of evolvable properties of a "
            "node into an array."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'logical',
                'attributes': ['dimension(:)', 'intent(  out)'],
                'variables':  ['array'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['count', 'offset', 'i'],
            },
        ],
    }
    content = "offset=1\n"
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"if (allocated(self%component{cap})) then\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    count=self%component{cap}(i)%serializeCount(propertyTypeActive)\n"
            f"    if (count > 0) call self%component{cap}(i)"
            f"%serializeNonNegative(array(offset:))\n"
            f"    offset=offset+count\n"
            f"  end do\n"
            f"end if\n"
        )
    function['content'] = content
    _bind(build, 'serializeNonNegative', function)


def Tree_Node_ODE_Serialize_Scales(build):
    """Generate `treeNodeSerializeScaleToArray`.  Mirrors
    `Tree_Node_ODE_Serialize_Scales`.  The body is fully static — it
    just packs from module-level `nodeScales` based on `propertyType`.
    """
    function = {
        'type':        'void',
        'name':        'treeNodeSerializeScaleToArray',
        'description': "Serialize scales of evolvable properties of a node into an array.",
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['dimension(:)', 'intent(  out)'],
                'variables':  ['array'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['propertyType'],
            },
        ],
        'content': (
            "!$GLC attributes unused :: self\n"
            "select case (propertyType)\n"
            "case (propertyTypeAll     )\n"
            " array(1:nodeSerializationCount        )=pack(nodeScales(1:nodeSerializationCount),"
            "                                                   .not.nodeAnalytics(1:nodeSerializationCount))\n"
            "case (propertyTypeActive  )\n"
            " array(1:nodeSerializationCountActive  )=pack(nodeScales(1:nodeSerializationCount),"
            ".not.nodeInactives(1:nodeSerializationCount) .and. .not.nodeAnalytics(1:nodeSerializationCount))\n"
            "case (propertyTypeInactive)\n"
            " array(1:nodeSerializationCountInactive)=pack(nodeScales(1:nodeSerializationCount),"
            "     nodeInactives(1:nodeSerializationCount) .and. .not.nodeAnalytics(1:nodeSerializationCount))\n"
            "case (propertyTypeNumerics)\n"
            " array(1:nodeSerializationCountActive  )=pack(nodeScales(1:nodeSerializationCount),"
            "                                                   .not.nodeAnalytics(1:nodeSerializationCount))\n"
            "end select\n"
        ),
    }
    _bind(build, 'serializeScales', function)


def Tree_Node_ODE_Serialize_Rates(build):
    """Generate `treeNodeSerializeRatesToArray`.  Mirrors `Tree_Node_ODE_Serialize_Rates`."""
    function = {
        'type':        'void',
        'name':        'treeNodeSerializeRatesToArray',
        'description': "Serialize rates of evolvable properties of a node into an array.",
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['dimension(:)', 'intent(  out)'],
                'variables':  ['array'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['propertyType'],
            },
        ],
        'content': (
            "!$GLC attributes unused :: self\n"
            "select case (propertyType)\n"
            "case (propertyTypeAll                          )\n"
            " array(1:nodeSerializationCount        )=nodeRates(1:nodeSerializationCount        )\n"
            "case (propertyTypeActive  ,propertyTypeNumerics)\n"
            " array(1:nodeSerializationCountActive  )=nodeRates(1:nodeSerializationCountActive  )\n"
            "case (propertyTypeInactive                     )\n"
            " array(1:nodeSerializationCountInactive)=nodeRates(1:nodeSerializationCountInactive)\n"
            "end select\n"
        ),
    }
    _bind(build, 'serializeRates', function)


def Tree_Node_ODE_Deserialize_Rates(build):
    """Generate `treeNodeDeserializeRatesToArray`.  Mirrors `Tree_Node_ODE_Deserialize_Rates`."""
    function = {
        'type':        'void',
        'name':        'treeNodeDeserializeRatesToArray',
        'description': "Deserialize rates of evolvable properties of a node from an array.",
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['dimension(:)', 'intent(in   )'],
                'variables':  ['array'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['propertyType'],
            },
        ],
        'content': (
            "!$GLC attributes unused :: self\n"
            "select case (propertyType)\n"
            "case (propertyTypeAll                          )\n"
            " nodeRatesActives(1:nodeSerializationCount        )=array(1:nodeSerializationCount        )\n"
            "case (propertyTypeActive  ,propertyTypeNumerics)\n"
            " nodeRatesActives(1:nodeSerializationCountActive  )=array(1:nodeSerializationCountActive  )\n"
            "case (propertyTypeInactive                     )\n"
            " nodeRatesActives(1:nodeSerializationCountInactive)=array(1:nodeSerializationCountInactive)\n"
            "end select\n"
        ),
    }
    _bind(build, 'deserializeRates', function)


def Tree_Node_ODE_Serialize_Inactive(build):
    """Generate `treeNodeSerializeInactiveToArray`.  Mirrors `Tree_Node_ODE_Serialize_Inactive`."""
    function = {
        'type':        'void',
        'name':        'treeNodeSerializeInactiveToArray',
        'description': "Serialize inactive statuses of evolvable properties of a node into an array.",
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'logical',
                'attributes': ['dimension(:)', 'intent(  out)'],
                'variables':  ['array'],
            },
        ],
        'content': (
            "!$GLC attributes unused :: self\n"
            "array(1:nodeSerializationCount)=nodeInactives(1:nodeSerializationCount)\n"
        ),
    }
    _bind(build, 'serializeInactives', function)


def Tree_Node_ODE_Name_From_Index(build):
    """Generate `treeNodePropertyNameFromIndex`.  Mirrors `Tree_Node_ODE_Name_From_Index`."""
    function = {
        'type':        'type(varying_string) => name',
        'name':        'treeNodePropertyNameFromIndex',
        'description': "Return the name of a property given its index within an array.",
        'modules':     ['ISO_Varying_String', 'Error'],
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
                'variables':  ['index', 'propertyType'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['count', 'i'],
            },
        ],
    }
    content = "name='unknown'\ncount =index\n"
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"if (allocated(self%component{cap})) then\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    name=self%component{cap}(i)%nameFromIndex(count,propertyType)\n"
            f"    if (count <= 0) return\n"
            f"  end do\n"
            f"end if\n"
        )
    content += (
        "if (name == 'unknown') call Error_Report('property index out of range'"
        "//{introspection:location})\n"
    )
    function['content'] = content
    _bind(build, 'nameFromIndex', function)


def Tree_Node_ODE_Offsets(build):
    """Generate `treeNodeSerializeOffsets`.  Mirrors `Tree_Node_ODE_Offsets`.

    Walks each active class to delegate `serializationOffsets` per
    component, then ensures the module-level `nodeScales`, `nodeRates`,
    `nodeRatesActives`, `nodeInactives`, `nodeAnalytics` arrays are
    sized to fit.
    """
    function = {
        'type':        'void',
        'name':        'treeNodeSerializeOffsets',
        'description': (
            r"Compute offsets into serialization arrays for \mono{treeNode} object."
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
                'variables':  ['propertyType'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i', 'count', 'countSubset'],
            },
        ],
    }
    content = "count      =0\ncountSubset=0\n"
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"if (allocated(self%component{cap})) then\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    call self%component{cap}(i)%serializationOffsets"
            f"(count,countSubset,propertyType)\n"
            f"  end do\n"
            f"end if\n"
        )
    content += (
        "if (.not.allocated(nodeScales)) then\n"
        "   allocate  (nodeScales      (count))\n"
        "   allocate  (nodeRates       (count))\n"
        "   allocate  (nodeRatesActives(count))\n"
        "   allocate  (nodeInactives   (count))\n"
        "   allocate  (nodeAnalytics   (count))\n"
        "else if (size(nodeScales) < count) then\n"
        "   deallocate(nodeScales             )\n"
        "   deallocate(nodeRates              )\n"
        "   deallocate(nodeRatesActives       )\n"
        "   deallocate(nodeInactives          )\n"
        "   deallocate(nodeAnalytics          )\n"
        "   allocate  (nodeScales      (count))\n"
        "   allocate  (nodeRates       (count))\n"
        "   allocate  (nodeRatesActives(count))\n"
        "   allocate  (nodeInactives   (count))\n"
        "   allocate  (nodeAnalytics   (count))\n"
        "end if\n"
        "nodeSerializationCount         =count\n"
        "select case (propertyType)\n"
        "case (propertyTypeInactive)\n"
        " nodeSerializationCountInactive=countSubset\n"
        "case (propertyTypeActive  )\n"
        " nodeSerializationCountActive  =countSubset\n"
        "case (propertyTypeNumerics)\n"
        " nodeSerializationCountActive  =countSubset\n"
        "end select\n"
    )
    function['content'] = content
    _bind(build, 'serializationOffsets', function)


# ---------------------------------------------------------------------------
# Helpers (shared with sister modules)
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
# Hook registration.  Order matches Perl ODESolver.pm:21-32.
# ---------------------------------------------------------------------------

register('treeNodeODESolver', 'functions', Tree_Node_ODE_Step_Initialize)
register('treeNodeODESolver', 'functions', Tree_Node_ODE_Serialize_Count)
register('treeNodeODESolver', 'functions', Tree_Node_ODE_Serialize_Values)
register('treeNodeODESolver', 'functions', Tree_Node_ODE_Deserialize_Values)
register('treeNodeODESolver', 'functions', Tree_Node_ODE_Serialize_Rates)
register('treeNodeODESolver', 'functions', Tree_Node_ODE_Deserialize_Rates)
register('treeNodeODESolver', 'functions', Tree_Node_ODE_Serialize_Scales)
register('treeNodeODESolver', 'functions', Tree_Node_ODE_Serialize_Inactive)
register('treeNodeODESolver', 'functions', Tree_Node_ODE_Serialize_NonNegative)
register('treeNodeODESolver', 'functions', Tree_Node_ODE_Name_From_Index)
register('treeNodeODESolver', 'functions', Tree_Node_ODE_Offsets)
