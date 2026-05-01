# Hierarchy initialization and finalization.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Hierarchy.pm.  Two
# `functions`-phase hooks emitting `nodeClassHierarchyInitialize` and
# `nodeClassHierarchyFinalize`.

import re


from Galacticus.Build.Components.Utils import register


_OUTPUT_CONDITION_RE = re.compile(r'\[\[([^\]]+)\]\]')


def Hierarchy_Initialization(build):
    """Generate `nodeClassHierarchyInitialize`.

    Mirrors `Hierarchy_Initialization`.  Reads run-time parameters that
    select which implementation each component class uses, then
    allocates `default<Class>Component` from the matching template and
    flips on `nodeComponent<Full>IsActiveValue` for the chosen
    implementation and every ancestor in its `extends` chain.
    """
    component_id_list = build.get('componentIdList') or []
    components        = build.get('components')      or {}
    component_classes = build.get('componentClasses') or {}

    variables = [
        {
            'intrinsic':  'type',
            'type':       'inputParameters',
            'attributes': ['intent(inout)'],
            'variables':  ['parameters_'],
        },
        {
            'intrinsic':  'type',
            'type':       'varying_string',
            'variables':  ['componentSelection', 'message'],
        },
    ]
    for impl_id in component_id_list:
        variables.append({
            'intrinsic':  'type',
            'type':       f'nodeComponent{impl_id}',
            'variables':  [f'default{impl_id}Component'],
        })

    function = {
        'type':        'void',
        'name':        'nodeClassHierarchyInitialize',
        'description': r"Initialize the \glc\ node/component class hierarchy.",
        'modules':     [
            'Input_Parameters',
            'ISO_Varying_String',
            'Error',
        ],
        'variables':   variables,
    }

    content = (
        "if (hierarchyInitialized == 0) then\n"
        "  !$omp critical (Node_Class_Hierarchy_Initialize)\n"
        "  if (hierarchyInitialized == 0) then\n"
        "    ! Parameters controlling output.\n"
    )

    output_conditions = set()

    for class_name in sorted(component_classes.keys()):
        class_dict = component_classes[class_name]
        cap_class  = _ucfirst(class_name)
        default_impl = class_dict.get('defaultImplementation', '')
        content += (
            "    ! Insert a function call to get the parameter controlling "
            "the choice of implementation for this class.\n"
            f"    if (.not.parameters_%isPresent('component{cap_class}')) "
            f"call parameters_%addParameter('component{cap_class}',"
            f"'{default_impl}')\n"
            "    !![\n"
            f"    <inputParameter>\n"
            f"      <name>component{cap_class}</name>\n"
            f"      <variable>componentSelection</variable>\n"
            f"      <source>parameters_</source>\n"
            f"      <defaultValue>var_str('{default_impl}')</defaultValue>\n"
            f"      <description>Specifies the implementation to be used "
            f"for the {class_name} component of nodes.</description>\n"
            f"      <type>string</type>\n"
            f"      <cardinality>1</cardinality>\n"
            f"    </inputParameter>\n"
            f"    !!]\n"
        )
        for member in class_dict.get('members') or []:
            full_name = cap_class + _ucfirst(member['name'])
            content += (
                f"    if (componentSelection == '{member['name']}') then\n"
                f"       allocate(default{cap_class}Component,"
                f"source=default{full_name}Component)\n"
                f"        nodeComponent{full_name}IsActiveValue=.true.\n"
            )
            # Walk up the `extends` chain, flipping on every ancestor.
            cursor = full_name
            while cursor:
                comp = components.get(cursor)
                if comp and 'extends' in comp:
                    ext = comp['extends']
                    cursor = _ucfirst(ext['class']) + _ucfirst(ext['name'])
                    content += (
                        f"\tnodeComponent{cursor}IsActiveValue=.true.\n"
                    )
                else:
                    cursor = ""
            content += "    end if\n"

            # Emit input parameters for any output-condition references.
            for prop in _properties(member):
                output = prop.get('output')
                if not isinstance(output, dict):
                    continue
                condition = output.get('condition')
                if not condition:
                    continue
                m = _OUTPUT_CONDITION_RE.search(condition)
                if not m:
                    continue
                parameter_name = m.group(1)
                if parameter_name in output_conditions:
                    continue
                output_conditions.add(parameter_name)
                content += (
                    "     !![\n"
                    f"     <inputParameter>\n"
                    f"      <name>{parameter_name}</name>\n"
                    f"      <variable>{parameter_name}</variable>\n"
                    f"      <source>parameters_</source>\n"
                    f"      <defaultValue>.false.</defaultValue>\n"
                    f"      <attachedTo>module</attachedTo>\n"
                    f"      <description>Specifies whether the "
                    f"\\mono{{{prop['name']}}} method of the "
                    f"\\mono{{{member['name']}}} implementation of the "
                    f"\\mono{{{class_name}}} component class should be output."
                    f"</description>\n"
                    f"      <type>string</type>\n"
                    f"      <cardinality>1</cardinality>\n"
                    f"     </inputParameter>\n"
                    "     !!]\n"
                )
                # Add a module-scope variable to store the parameter.
                build.setdefault('variables', []).append({
                    'intrinsic': 'logical',
                    'variables': [parameter_name],
                })

        # Validation block: error out if no implementation matched.
        members_list = '\n      '.join(
            f"message=message//char(10)//'    {m['name']}'"
            for m in class_dict.get('members') or []
        )
        content += (
            f"    if (.not.allocated(default{cap_class}Component)) then\n"
            f"      message='unrecognized class \"'//componentSelection"
            f"//'\" for \"{class_name}\" component'\n"
            f"      message=message//char(10)//'  available methods are:'\n"
            f"      {members_list}\n"
            f"      call Error_Report(message//{{introspection:location}})\n"
            f"    end if\n"
        )

    content += (
        "    ! Increment the count of hierarchy initializations.\n"
        "    hierarchyInitialized=hierarchyInitialized+1\n"
        "  end if\n"
        "  !$omp end critical (Node_Class_Hierarchy_Initialize)\n"
        "end if\n"
    )

    function['content'] = content
    build.setdefault('functions', []).append(function)


def Hierarchy_Finalization(build):
    """Generate `nodeClassHierarchyFinalize`.

    Mirrors `Hierarchy_Finalization`.  Decrements the init counter and,
    if zero, deallocates every `default<Class>Component` and tears down
    the mass-distributions cache.
    """
    component_classes = build.get('componentClasses') or {}
    deallocations = ' '.join(
        f"deallocate(default{c['name']}Component)\n"
        for c in component_classes.values()
    )

    function = {
        'type':        'void',
        'name':        'nodeClassHierarchyFinalize',
        'description': r"Finalize the \glc\ node/component class hierarchy.",
        'content':     (
            "!$omp critical (Node_Class_Hierarchy_Initialize)\n"
            "hierarchyInitialized=hierarchyInitialized-1\n"
            "if (hierarchyInitialized == 0) then\n"
            f" {deallocations}\n"
            " call massDistributionsDestroy()\n"
            "end if\n"
            "!$omp end critical (Node_Class_Hierarchy_Initialize)\n"
        ),
    }
    build.setdefault('functions', []).append(function)


def _properties(member):
    """Yield every property dict declared on `member`."""
    props = (member.get('properties') or {}).get('property')
    if props is None:
        return []
    if isinstance(props, list):
        return props
    if isinstance(props, dict):
        if all(isinstance(v, dict) for v in props.values()):
            return list(props.values())
        return [props]
    return []


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


# ---------------------------------------------------------------------------
# Hook registration.  Order matches Perl Hierarchy.pm:23-26.
# ---------------------------------------------------------------------------

register('hierarchy', 'functions', Hierarchy_Initialization)
register('hierarchy', 'functions', Hierarchy_Finalization)
