# Processes the functionClass lifecycle directives — objectBuilder,
# objectDestructor, referenceCountIncrement, referenceAcquire,
# referenceConstruct, deepCopy — by synthesising the appropriate
# reference-counted get / release / build / copy code and wiring in the
# required moduleUse imports and local declarations.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/ObjectBuilder.pm

import os
import re
import xml.etree.ElementTree as ET


from XML.Utils                                      import xml_to_dict
from Galacticus.Build.StateStorables                import function_class_entries
from Galacticus.Build.SourceTree                    import (
    walk_tree, insert_after_node,
)
from Galacticus.Build.SourceTree.Process            import register_process
from Galacticus.Build.SourceTree.Parse.Declarations import (
    add_declarations, add_attributes, declaration_exists, get_declaration,
)
from Galacticus.Build.SourceTree.Parse.ModuleUses   import add_uses
from Galacticus.Build.SourceTree.Process.SourceIntrospection import location


_STATE_STORABLES = None


def _load_state_storables():
    global _STATE_STORABLES
    if _STATE_STORABLES is not None:
        return _STATE_STORABLES
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        raise RuntimeError("process_object_builder: BUILDPATH is not set")
    path = os.path.join(build_path, 'stateStorables.xml')
    _STATE_STORABLES = xml_to_dict(ET.parse(path).getroot())
    return _STATE_STORABLES


def _function_classes_by_name(state_storables):
    """Return `{<name>Class: {module: …, …}}` keyed by class name."""
    out = {}
    for entry in function_class_entries(state_storables):
        name = entry.get('name')
        if name:
            out[name] = {k: v for k, v in entry.items() if k != 'name'}
    return out


def _result_name(parent_function):
    opener = parent_function.get('opener') or ''
    m = re.search(r'result\s*\(\s*([a-zA-Z0-9_]+)\s*\)\s*$', opener)
    if m:
        return m.group(1)
    return parent_function.get('name', '')


def _code_node(content, source, line=1):
    return {
        'type':       'code',
        'content':    content,
        'parent':     None,
        'firstChild': None,
        'sibling':    None,
        'source':     source,
        'line':       line,
    }


def _dict_to_parameters_xml(data):
    """Serialize a parsed `<default>` directive block back to a compact
    `<parameters …/>` XML string.

    Mirrors the Perl idiom
        my $xml = XML::Simple->XMLout($dict, RootName => 'parameters');
        $xml =~ s/\\s*\\n\\s*//g; $xml =~ s/\\s{2,}/ /g;
    with XML::Simple's default attribute-promotion semantics (scalar leaves
    become attributes on the enclosing element).
    """
    root = ET.Element('parameters')
    _fill_element(root, data)
    text = ET.tostring(root, encoding='unicode')
    # Compact whitespace to match Perl's normalised single-line output.
    text = re.sub(r'\s*\n\s*', '', text)
    text = re.sub(r'\s{2,}', ' ', text)
    return text


def _fill_element(elem, data):
    if not isinstance(data, dict):
        if data is not None:
            elem.text = str(data)
        return
    for key, value in data.items():
        if isinstance(value, dict):
            child = ET.SubElement(elem, key)
            _fill_element(child, value)
        elif isinstance(value, list):
            for item in value:
                child = ET.SubElement(elem, key)
                _fill_element(child, item)
        else:
            elem.set(key, str(value))


# ---------------------------------------------------------------------------
# Per-directive handlers
# ---------------------------------------------------------------------------

_SOURCE_TAG = ('Galacticus.Build.SourceTree.Process.ObjectBuilder'
               '.process_object_builder()')


def _handle_object_builder(node, state_storables, function_classes):
    directive = node['directive']

    if 'default' in directive and 'parameterName' not in directive:
        raise RuntimeError(
            "process_object_builder: objects with defaults must have explicit "
            "parameter names")

    parent = node['parent']
    parameter_name = directive.get('parameterName') or directive['class']
    default_name   = parameter_name == directive['class']
    loc_expr       = location(node, node.get('line', 0))
    parameters_default_required = False

    lines  = "   ! Determine where to build+store or point to the required object....\n"
    lines += f"   parametersCurrent => {directive['source']}\n"
    if 'parameterName' in directive:
        if 'default' in directive:
            parameters_default_required = True
            default_xml = _dict_to_parameters_xml(directive['default'])
            lines += (
                f"   if (.not.parametersCurrent%isPresent('{parameter_name}')) then\n"
                "    block\n"
                "     type(varying_string), dimension(1) :: allowedParameterName__\n"
                "     type(varying_string)               :: defaultXML__\n"
                f"     defaultXML__=var_str('{default_xml}')\n"
                f"     allowedParameterName__(1)=var_str('{parameter_name}')\n"
                "     parametersDefault=inputParameters(defaultXML__,"
                "allowedParameterNames=allowedParameterName__,noOutput=.true.)\n"
                "    end block\n"
                "    parametersCurrent => parametersDefault\n"
                "    parametersDefaultCreated=.true.\n"
                "  else\n"
                "    parametersDefaultCreated=.false.\n"
                "  end if\n"
            )
        else:
            lines += (
                f"   do while (.not.parametersCurrent%isPresent('{parameter_name}')"
                ".and.associated(parametersCurrent%parent))\n"
                "      parametersCurrent => parametersCurrent%parent\n"
                "   end do\n"
                f"   if (.not.parametersCurrent%isPresent('{parameter_name}')) "
                f"call Error_Report('[{parameter_name}] object is undefined'//{loc_expr})\n"
            )
    else:
        lines += (
            f"   do while (.not.parametersCurrent%isPresent('{parameter_name}')"
            ".and.associated(parametersCurrent%parent))\n"
            "      parametersCurrent => parametersCurrent%parent\n"
            "   end do\n"
        )

    if default_name:
        lines += f"   if (parametersCurrent%isPresent('{parameter_name}')) then\n"

    # Handle copy=<token=range> or copy=<scalar>.
    copy_instance  = ""
    copy_loop_open = ""
    copy_loop_close = ""
    if 'copy' in directive:
        copy_raw = directive['copy']
        m = re.match(r'^([a-zA-Z0-9_]+)=', copy_raw)
        if m:
            copy_instance   = f",copyInstance={m.group(1)}"
            copy_loop_open  = f"      do {copy_raw}\n"
            copy_loop_close = "      end do\n"
        else:
            copy_instance = f",copyInstance={copy_raw}"

    # Recursive-build guard when `self` is in scope.
    if declaration_exists(parent, 'self'):
        lines += (
            "      ! Test for a recursive build.\n"
            "      genericObject => self\n"
            "      select type (genericObject)\n"
            f"      class is ({directive['class']}Class)\n"
            f"         if (associated(parametersCurrent,{directive['source']}%parent)) "
            f"call Error_Report('recursive build of [{directive['class']}] class detected'//{loc_expr})\n"
            "      end select\n"
        )

    lines += copy_loop_open
    lines += "      ! Object should belong to the parameter node. Get the node and test whether the object has already been created in it.\n"
    lines += (
        f"      parameterNode => parametersCurrent%node('{parameter_name}'"
        f"{copy_instance})\n"
    )
    lines += "      if (parameterNode%objectCreated()) then\n"
    lines += "         ! Object already exists - simply get a pointer to it. Increment the reference counter as this is a new reference to an existing object.\n"
    lines += "         genericObject => parameterNode%objectGet()\n"
    lines += "         select type (genericObject)\n"
    lines += f"         class is ({directive['class']}Class)\n"
    lines += f"            {directive['name']} => genericObject\n"
    lines += f"            call {directive['name']}%referenceCountIncrement()\n"
    lines += "         class default\n"
    lines += (
        "            call Error_Report('parameter-stored object is not of "
        f"[{directive['class']}] class'//{loc_expr})\n"
    )
    lines += "         end select\n"
    lines += "      else\n"
    lines += "         ! Object does not yet exist - build it and store in the parameter node. Increment reference counter here as this is a newly constructed object.\n"
    name_arg = (f",parameterName='{parameter_name}'"
                if 'parameterName' in directive else '')
    lines += (
        f"         {directive['name']} => {directive['class']}"
        f"(parametersCurrent{copy_instance}{name_arg})\n"
    )
    lines += f"            call {directive['name']}%referenceCountIncrement()\n"
    lines += f"         call parameterNode%objectSet({directive['name']})\n"
    lines += f"         call {directive['name']}%autoHook()\n"
    lines += "      end if\n"
    lines += copy_loop_close

    if default_name:
        # Find a free warnObjectBuilderN__ name.
        i = 0
        warn_status = "warnObjectBuilder0__"
        while declaration_exists(parent, warn_status):
            i += 1
            warn_status = f"warnObjectBuilder{i}__"
        add_declarations(parent, [{
            'intrinsic':     'logical',
            'type':          None,
            'openMP':        False,
            'attributes':    ['save'],
            'variables':     [f"{warn_status}=.false."],
            'variableNames': [warn_status],
            'threadprivate': True,
        }])
        lines += "   else\n"
        lines += "      ! Object is not explicitly defined. Cause a default object of the class to be added to the parameters. Increment the reference count here as this is a new object.\n"
        lines += copy_loop_open
        lines += (
            f"      {directive['name']} => "
            f"{directive['class']}(parametersCurrent)\n"
        )
        lines += f"      call {directive['name']}%referenceCountIncrement()\n"
        lines += f"      call {directive['name']}%autoHook()\n"
        lines += copy_loop_close
        lines += (
            f"      if (mpiSelf%isMaster() .and. .not.{warn_status}) then\n"
            "         block\n"
            "            type(varying_string) :: parametersPath\n"
            "            parametersPath=parametersCurrent%path()\n"
            "            call Warn('Using default class for parameter "
            f"''['//char(parametersPath)//'{parameter_name}]''')\n"
            f"            {warn_status}=.true.\n"
            "         end block\n"
            "      end if\n"
            "   end if\n"
        )

    if 'parameterName' in directive and 'default' in directive:
        lines += "   if (parametersDefaultCreated) call parametersDefault%destroy()\n"

    insert_after_node(node, [_code_node(lines, _SOURCE_TAG)])

    # MPI + varying_string module uses when using a default-named object.
    if default_name:
        add_uses(parent, {
            'moduleUse': {
                'MPI_Utilities':      {'intrinsic': False,
                                       'only': {'mpiSelf': True}},
                'ISO_Varying_String': {'intrinsic': False,
                                       'only': {'varying_string': True}},
            },
            'moduleOrder': ['MPI_Utilities', 'ISO_Varying_String'],
        })

    # Main module-use block.  The class-providing module is skipped when the
    # enclosing node is part of the same module that declares the class
    # (`isSelf` in Perl).
    class_key = directive['class'] + 'Class'
    class_entry = function_classes.get(class_key) or {}
    module_name = class_entry.get('module')
    if module_name:
        is_self = False
        for n in walk_tree(node):   # walks from node DOWN, not ideal but OK
            pass
        # Use a tree-wide search from the file root.
        root = node
        while root.get('parent') is not None:
            root = root['parent']
        for n in walk_tree(root):
            if 'directive' in n:
                match_key = (n.get('type') or '') + 'Class'
                if match_key in function_classes and \
                        function_classes[match_key].get('module') == module_name:
                    is_self = True
                    break

        module_uses = {
            'Input_Parameters':   {'intrinsic': False,
                                   'only': {'inputParameter': True}},
            'Error':              {'intrinsic': False,
                                   'only': {'Warn': True}},
            'ISO_Varying_String': {'intrinsic': False,
                                   'only': {'char': True}},
        }
        if parameters_default_required:
            module_uses['ISO_Varying_String']['only']['var_str'] = True
        if not is_self:
            module_uses[module_name] = {
                'intrinsic': False,
                'only':      {directive['class']:               True,
                              directive['class'] + 'Class':     True},
            }
        add_uses(parent, {
            'moduleUse':   module_uses,
            'moduleOrder': sorted(module_uses.keys()),
        })

    # Add required local declarations once per parent.
    if not parent.get('objectBuilderDeclarations'):
        add_declarations(parent, [
            {
                'intrinsic':     'type',
                'type':          'inputParameters',
                'openMP':        False,
                'attributes':    ['pointer'],
                'variables':     ['parametersCurrent'],
                'variableNames': ['parametersCurrent'],
            },
            {
                'intrinsic':     'type',
                'type':          'inputParameter',
                'openMP':        False,
                'attributes':    ['pointer'],
                'variables':     ['parameterNode'],
                'variableNames': ['parameterNode'],
            },
            {
                'intrinsic':     'class',
                'type':          '*',
                'openMP':        False,
                'attributes':    ['pointer'],
                'variables':     ['genericObject'],
                'variableNames': ['genericObject'],
            },
        ])
        add_uses(parent, {
            'moduleUse':   {'Error': {'intrinsic': False, 'all': True}},
            'moduleOrder': ['Error'],
        })
        parent['objectBuilderDeclarations'] = True

    if parameters_default_required and not parent.get(
            'objectBuilderDefaultDeclarations'):
        add_declarations(parent, [
            {
                'intrinsic':     'type',
                'type':          'inputParameters',
                'openMP':        False,
                'attributes':    ['target'],
                'variables':     ['parametersDefault'],
                'variableNames': ['parametersDefault'],
            },
            {
                'intrinsic':     'logical',
                'type':          None,
                'openMP':        False,
                'attributes':    [],
                'variables':     ['parametersDefaultCreated'],
                'variableNames': ['parametersDefaultCreated'],
            },
        ])
        add_uses(parent, {
            'moduleUse': {
                'ISO_Varying_String': {
                    'intrinsic': False,
                    'only': {'varying_string': True, 'assignment(=)': True},
                },
            },
            'moduleOrder': ['ISO_Varying_String'],
        })
        parent['objectBuilderDefaultDeclarations'] = True

    # Make `source` and `self` have a `target` attribute if they don't
    # already have `target`/`pointer`, so that `=>` / `isAssociated(...)`
    # work on them.  Perl tracks done-work via
    # `$parent->{objectBuilderAttributes}->{<source>}`.
    attrs_map = parent.setdefault('objectBuilderAttributes', {})
    src_name = directive['source']
    if src_name not in attrs_map:
        if declaration_exists(parent, src_name):
            decl = get_declaration(parent, src_name)
            if not any(a in ('target', 'pointer')
                       for a in decl.get('attributes') or []):
                add_attributes(parent, src_name, ['target'])
        if declaration_exists(parent, 'self'):
            decl = get_declaration(parent, 'self')
            if not any(a in ('target', 'pointer')
                       for a in decl.get('attributes') or []):
                add_attributes(parent, 'self', ['target'])
        attrs_map[src_name] = True


def _handle_object_destructor(node):
    directive = node['directive']
    loc_expr = location(node, node.get('line', 0))
    name = directive['name']
    lines  = f"if (associated({name})) then\n"
    lines += "   ! Decrement the reference count, and decide if this object can be destroyed.\n"
    lines += f"   referenceCount_={name}%referenceCountDecrement()\n"
    lines += "   if (referenceCount_ == 0) then\n"
    lines += "      ! Deallocate the pointer.\n"
    lines += f"      deallocate({name})\n"
    lines += "   else if (referenceCount_ < 0) then\n"
    lines += "      ! Negative counter - should not happen.\n"
    lines += (
        f"      call Error_Report('negative reference counter in object "
        f"\"{name}\"'//{loc_expr})\n"
    )
    lines += "   else\n"
    lines += "      ! Nullify the pointer.\n"
    if directive.get('nullify') != 'no':
        lines += f"      nullify({name})\n"
    lines += "   end if\n"
    lines += "end if\n"

    insert_after_node(node, [_code_node(lines, _SOURCE_TAG)])

    parent = node['parent']
    if not parent.get('objectDestructorDeclarations'):
        add_declarations(parent, [{
            'intrinsic':     'integer',
            'type':          None,
            'openMP':        False,
            'attributes':    ['save'],
            'variables':     ['referenceCount_'],
            'variableNames': ['referenceCount_'],
            'threadprivate': True,
        }])
        parent['objectDestructorDeclarations'] = True

    add_uses(parent, {
        'moduleUse':   {'Error': {'intrinsic': False, 'all': True}},
        'moduleOrder': ['Error'],
    })


def _handle_reference_count_increment(node):
    directive = node['directive']
    owner = directive.get('owner')
    prefix = f"{owner}%" if owner else ''
    code = f"call {prefix}{directive['object']}%referenceCountIncrement()\n"
    insert_after_node(node, [_code_node(code, _SOURCE_TAG)])


def _handle_reference_acquire(node):
    directive = node['directive']
    owner = directive.get('owner')
    prefix = f"{owner}%" if owner else ''
    code = (
        f"{prefix}{directive['target']} => {directive['source']}\n"
        f"call {prefix}{directive['target']}%referenceCountIncrement()\n"
    )
    insert_after_node(node, [_code_node(code, _SOURCE_TAG)])


def _handle_reference_construct(node):
    directive = node['directive']
    constructor = (directive.get('constructor') or '').strip()
    name_associated = directive.get('nameAssociated')
    owner = directive.get('owner')
    if name_associated:
        object_name = name_associated
    else:
        prefix = f"{owner}%" if owner else ''
        object_name = f"{prefix}{directive['object']}"
    code = (
        f"{object_name}={constructor}\n"
        f"call {object_name}%referenceCountIncrement()\n"
        f"call {object_name}%autoHook()\n"
    )
    insert_after_node(node, [_code_node(code, _SOURCE_TAG)])


def _handle_deep_copy(node):
    directive = node['directive']
    code = (
        f"call {directive['source']}%deepCopy({directive['destination']})\n"
        f"{directive['source']}%copiedSelf => {directive['destination']}\n"
        f"call {directive['destination']}%autoHook()\n"
    )
    insert_after_node(node, [_code_node(code, _SOURCE_TAG)])


# ---------------------------------------------------------------------------
# Main dispatcher
# ---------------------------------------------------------------------------

def process_object_builder(tree, options):
    """Mirrors Process_ObjectBuilder() from ObjectBuilder.pm."""
    state_storables  = _load_state_storables()
    function_classes = _function_classes_by_name(state_storables)

    handled_types = (
        'objectBuilder', 'objectDestructor', 'referenceCountIncrement',
        'referenceAcquire', 'referenceConstruct', 'deepCopy',
    )
    for node in list(walk_tree(tree)):
        ntype = node.get('type')
        if ntype not in handled_types:
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue

        if ntype == 'objectBuilder':
            _handle_object_builder(node, state_storables, function_classes)
        elif ntype == 'objectDestructor':
            _handle_object_destructor(node)
        elif ntype == 'referenceCountIncrement':
            _handle_reference_count_increment(node)
        elif ntype == 'referenceAcquire':
            _handle_reference_acquire(node)
        elif ntype == 'referenceConstruct':
            _handle_reference_construct(node)
        elif ntype == 'deepCopy':
            _handle_deep_copy(node)

        directive['processed'] = True


register_process('objectBuilder', process_object_builder)
