"""Build the `nodeComponent` base type at the start of the components-build
pipeline.

Andrew Benson (ported to Python 2026)

A single hook on
the `types` phase that defines `nodeComponent`, the abstract base type
every per-class component derives from.  The type carries the standard
set of ~25 type-bound methods (lifecycle, ODE solver, output, mass
distribution, …) plus six meta-property `add…` methods, and stores a
pointer back to its owning `treeNode` along with one allocatable array
per meta-property kind.
"""



from Galacticus.Build.Components.Utils import register
from Galacticus.Build.Components.Classes.MetaProperties import meta_property_types


def Build_Node_Component_Class(build):
    """Define `nodeComponent` on `build['types']`.

    The body is almost entirely
    declarative — the data structure here matches the dict passed to
    Components/__init__.py's `derived_types_serialize`.
    """
    type_bound_functions = list(_BASE_BOUND_FUNCTIONS)

    # Add one `addFloatRank0MetaProperty` / `addLongIntegerRank1MetaProperty`
    # / … entry per meta-property type.
    for mpt in meta_property_types:
        label = _ucfirst(mpt['label'])
        rank  = mpt['rank']
        type_bound_functions.append({
            'type':        'procedure',
            'name':        f"add{label}Rank{rank}MetaProperty",
            'function':    (
                f"Node_Component_Generic_Add_{label}"
                f"_Rank{rank}_Meta_Property"
            ),
            'description': (
                f"Add a rank-{rank} {mpt['label']} meta-property to this class."
            ),
            'returnType':  r"``integer``",
            'arguments':   _meta_property_arg_doc(mpt),
        })

    # Data content: a back-pointer to the host treeNode plus one
    # allocatable array per meta-property kind.
    data_content = [
        {
            'intrinsic':  'type',
            'type':       'treeNode',
            'attributes': ['pointer', 'public'],
            'variables':  ['hostNode => null()'],
        },
    ]
    for mpt in meta_property_types:
        decl = {}
        if mpt['rank'] == 0:
            decl['intrinsic'] = mpt['intrinsic']
            if 'type' in mpt:
                decl['type'] = mpt['type']
        else:
            decl['intrinsic'] = 'type'
            decl['type'     ] = (
                f"{mpt['label']}Rank{mpt['rank']}MetaProperty"
            )
        decl['attributes'] = ['allocatable', 'dimension(:)']
        decl['variables' ] = [
            f"{mpt['label']}Rank{mpt['rank']}MetaProperties"
        ]
        data_content.append(decl)

    build.setdefault('types', {})['nodeComponent'] = {
        'name':           'nodeComponent',
        'comment':        r"A class for components in :term:`node`\ s.",
        'isPublic':       True,
        'boundFunctions': type_bound_functions,
        'dataContent':    data_content,
    }


def _meta_property_arg_doc(mpt):
    """Return the documentation-only argument string for an `add…MetaProperty`
    method.
    """
    is_evolvable = r", ``logical`` [isEvolvable]" \
        if (mpt['label'] == 'float' and mpt['rank'] == 0) else ""
    return (
        r"``type(varying_string)`` label, "
        r"``character(len=*)`` name"
        + is_evolvable
        + r", ``logical`` [isCreator]"
    )


# ---------------------------------------------------------------------------
# Standard type-bound function table.  Each entry is a dict with fields
# `type`, `name`, `function`, `description`, `returnType`, `arguments`,
# and optional `mappable`.
# ---------------------------------------------------------------------------

_BASE_BOUND_FUNCTIONS = [
    {
        'type':        'procedure',
        'name':        'type',
        'function':    'Node_Component_Generic_Type',
        'description': "Return the type of this object.",
        'returnType':  r"``type(varying_string)``",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'host',
        'function':    'Node_Component_Host_Node',
        'description': r"Return a pointer to the host ``treeNode`` object.",
        'returnType':  r"``*type(treeNode)``",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'destroy',
        'function':    'Node_Component_Generic_Destroy',
        'description': "Destroy the object.",
        'returnType':  r"``void``",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'serializeCount',
        'function':    'Node_Component_Serialize_Count_Zero',
        'description': "Return a count of the number of evolvable quantities to be evolved.",
        'returnType':  r"``integer``",
        'arguments':   r"``integer`` propertyType [in]",
    },
    {
        'type':        'procedure',
        'name':        'serializationOffsets',
        'function':    'Node_Component_Serialization_Offsets',
        'description': "Set offsets into serialization arrays.",
        'returnType':  r"``void``",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'serializeValues',
        'function':    'Node_Component_Serialize_Null',
        'description': "Serialize the evolvable quantities to an array.",
        'returnType':  r"``void``",
        'arguments':   r"``double(:)`` array [out], ``integer`` propertyType [in]",
    },
    {
        'type':        'procedure',
        'name':        'serializeNonNegative',
        'function':    'Node_Component_Serialize_NonNegative_Null',
        'description': "Serialize the non-negative status of evolvable quantities to an array.",
        'returnType':  r"``void``",
        'arguments':   r"``logical(:)`` array [out]",
    },
    {
        'type':        'procedure',
        'name':        'deserializeRaw',
        'function':    'Node_Component_Read_Raw_Null',
        'description': "Read properties from raw file.",
        'returnType':  r"``void``",
        'arguments':   r"``integer`` fileHandle [in]",
    },
    {
        'type':        'procedure',
        'name':        'deserializeValues',
        'function':    'Node_Component_Deserialize_Null',
        'description': "Deserialize the evolvable quantities from an array.",
        'returnType':  r"``void``",
        'arguments':   r"``double(:)`` array [in], ``integer`` propertyType [in]",
    },
    {
        'type':        'procedure',
        'name':        'odeStepRatesInitialize',
        'function':    'Node_Component_ODE_Step_Initialize_Null',
        'description': "Initialize rates for evolvable properties.",
        'returnType':  r"``void``",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'odeStepScalesInitialize',
        'function':    'Node_Component_ODE_Step_Initialize_Null',
        'description': "Initialize scales for evolvable properties.",
        'returnType':  r"``void``",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'serializeASCII',
        'function':    'Node_Component_Dump_Null',
        'description': "Generate an ASCII dump of all properties.",
        'returnType':  r"``void``",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'serializeXML',
        'function':    'Node_Component_Dump_XML_Null',
        'description': "Generate an XML dump of all properties.",
        'returnType':  r"``void``",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'serializeRaw',
        'function':    'Node_Component_Dump_Raw_Null',
        'description': "Generate a binary dump of all properties.",
        'returnType':  r"``void``",
        'arguments':   r"``integer`` fileHandle [in]",
    },
    {
        'type':        'procedure',
        'name':        'outputCount',
        'function':    'Node_Component_Output_Count_Null',
        'description': "Compute a count of outputtable properties.",
        'returnType':  r"``void``",
        'arguments':   (
            r"``integer`` integerPropertyCount [inout], "
            r"``integer`` doublePropertyCount [inout], "
            r"``double`` time [in], "
            r"``integer`` instance [in]"
        ),
    },
    {
        'type':        'procedure',
        'name':        'outputNames',
        'function':    'Node_Component_Output_Names_Null',
        'description': "Generate names of outputtable properties.",
        'returnType':  r"``void``",
        'arguments':   (
            r"``integer`` integerProperty [inout], "
            r"``type(outputPropertyInteger)(:)`` integerProperties [inout], "
            r"``integer`` doubleProperty [inout], "
            r"``type(otuputPropertyDouble)(:)`` doubleProperties [inout], "
            r"``double`` time [in], "
            r"``integer`` instance [in]"
        ),
    },
    {
        'type':        'procedure',
        'name':        'output',
        'function':    'Node_Component_Output_Null',
        'description': "Generate values of outputtable properties.",
        'returnType':  r"``void``",
        'arguments':   (
            r"``integer`` integerProperty [inout], "
            r"``integer`` integerBufferCount [inout], "
            r"``type(outputPropertyInteger)(:)`` integerProperties [inout], "
            r"``integer`` doubleProperty [inout], "
            r"``integer`` doubleBufferCount [inout], "
            r"``type(outputPropertyDouble)(:)`` doubleProperties [inout], "
            r"``double`` time [in], "
            r"``integer`` instance [in]"
        ),
    },
    {
        'type':        'procedure',
        'name':        'massDistribution',
        'function':    'Node_Component_Mass_Distribution_Null',
        'description': "Return the mass distribution for this component.",
        'returnType':  r"``class(massDistribution)``",
        'arguments':   (
            r"``type(enumerationComponentTypeType)`` [componentType] [in], "
            r"``type(enumeratioMassTypeType)`` [massType] [in], "
            r"``type(enumeratioWeightByType)`` [weightBy] [in], "
            r"``integer`` [weightIndex] [in]"
        ),
    },
    {
        'type':        'procedure',
        'name':        'massBaryonic',
        'function':    'Node_Component_Mass_Baryonic_Null',
        'description': "Return the total baryonic mass for this component.",
        'returnType':  r"``double``",
        'arguments':   "",
    },
    {
        'type':        'procedure',
        'name':        'density',
        'function':    'Node_Component_Density_Null',
        'description': "Compute the density.",
        'mappable':    "summation",
        'returnType':  r"``double``",
        'arguments':   (
            r"``double(3)`` positionSpherical [in], "
            r"``enumerationComponentTypeType`` [componentType] [in], "
            r"``enumerationMassTypeType`` [massType] [in], "
            r"``enumerationWeightByType`` [weightBy] [in], "
            r"``integer`` [weightIndex] [in]"
        ),
    },
    {
        'type':        'procedure',
        'name':        'densitySphericalAverage',
        'function':    'Node_Component_Density_Spherical_Average_Null',
        'description': "Compute the spherically-averaged density.",
        'mappable':    "summation",
        'returnType':  r"``double``",
        'arguments':   (
            r"``double`` radius [in], "
            r"``enumerationComponentTypeType`` [componentType] [in], "
            r"``enumerationMassTypeType`` [massType] [in], "
            r"``enumerationWeightByType`` [weightBy] [in], "
            r"``integer`` [weightIndex] [in]"
        ),
    },
    {
        'type':        'procedure',
        'name':        'surfaceDensity',
        'function':    'Node_Component_Surface_Density_Null',
        'description': "Compute the surface density.",
        'mappable':    "summation",
        'returnType':  r"``double``",
        'arguments':   (
            r"``double(3)`` positionCylindrical [in], "
            r"``enumerationComponentTypeType`` [componentType] [in], "
            r"``enumerationMassTypeType`` [massType] [in], "
            r"``enumerationWeightByType`` [weightBy] [in], "
            r"``integer`` [weightIndex] [in]"
        ),
    },
]


# ---------------------------------------------------------------------------
# Hook registration
# ---------------------------------------------------------------------------

register('baseTypes', 'types', Build_Node_Component_Class)


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text
