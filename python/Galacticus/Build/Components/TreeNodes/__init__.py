"""Build the `treeNode` class at the start of the components-build pipeline.

Andrew Benson (ported to Python 2026)

Two hooks:

  types     → Build_Tree_Node_Class    — declares the `treeNode` Fortran
                                         type with all of its
                                         type-bound methods and data
                                         content (one allocatable array
                                         per active component class).
  functions → Insert_Interrupt_Interface — emits the `interruptTask`
                                           abstract interface used by
                                           differential-evolution
                                           interrupt handlers.
"""



from Galacticus.Build.Components.Utils import register


def Build_Tree_Node_Class(build):
    """Define `treeNode` on `build['types']`.

    Almost entirely declarative; only
    the per-class data-content section depends on
    `componentClassListActive` (populated by Classes/Gather_Classes).
    """
    type_bound_functions = list(_BASE_BOUND_FUNCTIONS)

    data_content = list(_BASE_DATA_CONTENT)
    for class_name in build.get('componentClassListActive') or []:
        data_content.append({
            'intrinsic':  'class',
            'type':       f"nodeComponent{_ucfirst(class_name)}",
            'attributes': ['allocatable', 'dimension(:)'],
            'variables':  [f"component{_ucfirst(class_name)}"],
            'comment':    f"A generic {class_name} object.",
        })

    build.setdefault('types', {})['treeNode'] = {
        'name':           'treeNode',
        'comment':        r"A class for \glspl{node} in merger trees.",
        'isPublic':       True,
        'boundFunctions': type_bound_functions,
        'dataContent':    data_content,
    }


def Insert_Interrupt_Interface(build):
    """Emit the `interruptTask` abstract interface.
    """
    build.setdefault('interfaces', {})['interruptTask'] = {
        'name':      'interruptTask',
        'comment':   "Interface for differential evolution interrupt tasks.",
        'intrinsic': 'void',
        'data':      [
            {
                'intrinsic':  'type',
                'type':       'treeNode',
                'attributes': ['target', 'intent(inout)'],
                'variables':  ['node'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['intent(in   )', 'optional'],
                'variables':  ['timeEnd'],
            },
        ],
    }


# ---------------------------------------------------------------------------
# Standard data content for the `treeNode` type.
# ---------------------------------------------------------------------------

_BASE_DATA_CONTENT = [
    {
        'intrinsic':  'integer',
        'type':       'kind=kind_int8',
        'variables':  ['indexValue', 'uniqueIdValue'],
    },
    {
        'intrinsic':  'double precision',
        'variables':  ['timeStepValue'],
    },
    {
        'intrinsic':  'double precision',
        'variables':  ['subsamplingWeightValue'],
    },
    {
        'intrinsic':  'type',
        'type':       'treeNode',
        'attributes': ['pointer', 'public'],
        'variables':  ['parent', 'firstChild', 'sibling', 'firstSatellite',
                       'mergeTarget', 'firstMergee', 'siblingMergee',
                       'formationNode'],
    },
    {
        'intrinsic':  'logical',
        'attributes': ['public'],
        'variables':  ['isPhysicallyPlausible', 'isSolvable'],
    },
    {
        'intrinsic':  'class',
        'type':       'nodeEvent',
        'attributes': ['public', 'pointer'],
        'variables':  ['event'],
    },
    {
        'intrinsic':  'type',
        'type':       'mergerTree',
        'attributes': ['public', 'pointer'],
        'variables':  ['hostTree'],
    },
]

# ---------------------------------------------------------------------------
# Standard type-bound function table.
# ---------------------------------------------------------------------------

_BASE_BOUND_FUNCTIONS = [
    {
        'type':        'procedure',
        'name':        'type',
        'function':    'Tree_Node_Type',
        'description': "Return the type of this node.",
        'returnType':  r"\textcolor{red}{\textless type(varying\_string)\textgreater}",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'index',
        'function':    'Tree_Node_Index',
        'description': "Return the index of this node.",
        'returnType':  r"\textcolor{red}{\textless integer(kind\_int8)\textgreater}",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'indexSet',
        'function':    'Tree_Node_Index_Set',
        'description': "Set the index of this node.",
        'returnType':  r"\void",
        'arguments':   r"\textcolor{red}{\textless integer(kind\_int8)\textgreater} index\argin",
    },
    {
        'type':        'procedure',
        'name':        'timeStep',
        'function':    'Tree_Node_Time_Step',
        'description': "Return the time-step last used by this node.",
        'returnType':  r"\doublezero",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'timeStepSet',
        'function':    'Tree_Node_Time_Step_Set',
        'description': "Set the time-step used by this node.",
        'returnType':  r"\void",
        'arguments':   r"\doublezero\ timeStep\argin",
    },
    {
        'type':        'procedure',
        'name':        'subsamplingWeight',
        'function':    'Tree_Node_Subsampling_Weight',
        'description': "Return the subsampling weight of this node.",
        'returnType':  r"\doublezero",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'subsamplingWeightSet',
        'function':    'Tree_Node_Subsampling_Weight_Set',
        'description': "Set the subsampling weight of this node.",
        'returnType':  r"\void",
        'arguments':   r"\doublezero\ subsamplingWeight\argin",
    },
    {
        'type':        'procedure',
        'name':        'uniqueID',
        'function':    'Tree_Node_Unique_ID',
        'description': "Return the unique identifier for this node.",
        'returnType':  r"\textcolor{red}{\textless integer(kind\_int8)\textgreater}",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'uniqueIDSet',
        'function':    'Tree_Node_Unique_ID_Set',
        'description': "Set the unique identifier for this node.",
        'returnType':  r"\void",
        'arguments':   r"\textcolor{red}{\textless integer(kind\_int8)\textgreater} uniqueID\argin",
    },
    {
        'type':        'procedure',
        'name':        'removeFromHost',
        'function':    'Tree_Node_Remove_From_Host',
        'description': "Remove this node from the satellite population of its host halo.",
        'returnType':  r"\void",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'removeFromMergee',
        'function':    'Tree_Node_Remove_From_Mergee',
        'description': "Remove this node from the list of mergees associated with its merge target.",
        'returnType':  r"\void",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'isPrimaryProgenitor',
        'function':    'treeNodeIsPrimaryProgenitor',
        'description': "Return true if this node is the primary progenitor of its descendant, false otherwise.",
        'returnType':  r"\logicalzero",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'function':    'Tree_Node_Is_Primary_Progenitor_Of_Index',
    },
    {
        'type':        'procedure',
        'function':    'Tree_Node_Is_Primary_Progenitor_Of_Node',
    },
    {
        'type':        'generic',
        'name':        'isPrimaryProgenitorOf',
        'function':    ['Tree_Node_Is_Primary_Progenitor_Of_Index',
                        'Tree_Node_Is_Primary_Progenitor_Of_Node'],
        'description': "Return true is this node is the primary progenitor of the specified (by index or pointer) node, false otherwise.",
        'returnType':  r"\logicalzero",
        'arguments':   r"\textcolor{red}{\textless integer(kind\_int8)\textgreater} targetNodeIndex\argin|\textcolor{red}{\textless *type(treeNode)\textgreater} targetNode\argin",
    },
    {
        'type':        'procedure',
        'function':    'Tree_Node_Is_Progenitor_Of_Index',
    },
    {
        'type':        'procedure',
        'function':    'Tree_Node_Is_Progenitor_Of_Node',
    },
    {
        'type':        'generic',
        'name':        'isProgenitorOf',
        'function':    ['Tree_Node_Is_Progenitor_Of_Index',
                        'Tree_Node_Is_Progenitor_Of_Node'],
        'description': "Return true is this node is a progenitor of the specified (by index or pointer) node, false otherwise.",
        'returnType':  r"\logicalzero",
        'arguments':   r"\textcolor{red}{\textless integer(kind\_int8)\textgreater} targetNodeIndex\argin|\textcolor{red}{\textless *type(treeNode)\textgreater} targetNode\argin",
    },
    {
        'type':        'procedure',
        'name':        'isOnMainBranch',
        'function':    'Tree_Node_Is_On_Main_Branch',
        'description': "Return true if this node is on the main branch of its tree, false otherwise.",
        'returnType':  r"\logicalzero",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'isSatellite',
        'function':    'Tree_Node_Is_Satellite',
        'description': "Return true if this node is a satellite, false otherwise.",
        'returnType':  r"\logicalzero",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'isolatedParent',
        'function':    'Tree_Node_Get_Isolated_Parent',
        'description': "Return a pointer to the isolated parent node of this node.",
        'returnType':  r"\textcolor{red}{\textless *type(treeNode)\textgreater}",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'lastSatellite',
        'function':    'Tree_Node_Get_Last_Satellite',
        'description': "Return a pointer to the last satellite in the list of satellites beloning to this node.",
        'returnType':  r"\textcolor{red}{\textless *type(treeNode)\textgreater}",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'earliestProgenitor',
        'function':    'Tree_Node_Get_Earliest_Progenitor',
        'description': "Return a pointer to the earliest progenitor (along the main branch) of this node.",
        'returnType':  r"\textcolor{red}{\textless *type(treeNode)\textgreater}",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'mergesWith',
        'function':    'Tree_Node_Merges_With_Node',
        'description': "Return a pointer to the node with which this node will merge.",
        'returnType':  r"\textcolor{red}{\textless *type(treeNode)\textgreater}",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'walkBranchWithSatellites',
        'function':    'treeNodeWalkBranchWithSatellites',
        'description': "Return a pointer to the next node when performing a walk of a single branch of the tree, including satellites.",
        'returnType':  r"\void",
        'arguments':   r"\textcolor{red}{\textless *type(treeNode)\textgreater} startNode\arginout",
    },
    {
        'type':        'procedure',
        'name':        'walkTreeWithSatellites',
        'function':    'treeNodeWalkTreeWithSatellites',
        'description': "Return a pointer to the next node when performing a walk of the entire tree, including satellites.",
        'returnType':  r"\void",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'destroyBranch',
        'function':    'treeNodeDestroyBranch',
        'description': "Destroy a branch of a merger tree rooted at this node.",
        'returnType':  r"\void",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'attachEvent',
        'function':    'Tree_Node_Attach_Event',
        'description': r"Attach a \mono{nodeEvent} object to this node.",
        'returnType':  r"\void",
        'arguments':   r"\textcolor{red}{\textless *class(nodeEvent)\textgreater} newEvent\arginout",
    },
    {
        'type':        'procedure',
        'name':        'removePairedEvent',
        'function':    'Tree_Node_Remove_Paired_Event',
        'description': r"Remove a paired \mono{nodeEvent} from this node.",
        'returnType':  r"\void",
        'arguments':   r"\textcolor{red}{\textless class(nodeEvent)\textgreater} event\argin",
    },
]


# ---------------------------------------------------------------------------
# Hook registration
# ---------------------------------------------------------------------------

register('treeNodes', 'types',     Build_Tree_Node_Class)
register('treeNodes', 'functions', Insert_Interrupt_Interface)


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text
