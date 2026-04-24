# Processes `inputParametersValidate` directives: emits a runtime check
# against the `allowedParameters` list of the enclosing function-class
# object, optionally with multi-parameter and extra-allowed-name lists.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/InputParametersValidate.pm

import os
import re
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from XML.Utils                                      import xml_to_dict
from Galacticus.Build.StateStorables                import (
    function_class_names    as _shared_function_class_names,
)
from Galacticus.Build.SourceTree                    import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process            import register_process
from Galacticus.Build.SourceTree.Parse.Declarations import (
    add_declarations, declaration_exists,
)
from Galacticus.Build.SourceTree.Parse.ModuleUses   import add_uses


def _load_state_storables():
    """Load `$BUILDPATH/stateStorables.xml`.  Required for this hook."""
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        raise RuntimeError(
            "process_input_parameters_validate: BUILDPATH is not set")
    path = os.path.join(build_path, 'stateStorables.xml')
    if not os.path.exists(path):
        raise RuntimeError(
            f"process_input_parameters_validate: stateStorables.xml not found at {path}")
    root = ET.parse(path).getroot()
    return xml_to_dict(root)


def _function_class_result_name(parent_function):
    """Return the `result(...)` identifier for a function opener, or the
    function's own name if no explicit `result(...)` is present.

    Mirrors the regex at InputParametersValidate.pm:105-109.
    """
    opener = parent_function.get('opener') or ''
    m = re.search(r'result\s*\(\s*([a-zA-Z0-9_]+)\s*\)\s*$', opener)
    if m:
        return m.group(1)
    return parent_function.get('name', '')


def process_input_parameters_validate(tree, options):
    """Mirrors Process_InputParametersValidate() from InputParametersValidate.pm."""
    state_storables = _load_state_storables()
    class_names     = _shared_function_class_names(state_storables)

    function_class_name = None

    for node in walk_tree(tree):
        ntype = node.get('type', '')
        if ntype + 'Class' in class_names:
            function_class_name = ntype

        if ntype != 'inputParametersValidate':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        directive['processed'] = True

        if 'source' not in directive:
            raise RuntimeError(
                "process_input_parameters_validate: no `source` given")
        source = directive['source']
        parent = node['parent']

        variable_name = 'allowedParameterNames_'
        if not declaration_exists(parent, variable_name):
            add_declarations(parent, [{
                'intrinsic':     'type',
                'type':          'varying_string',
                'openMP':        False,
                'attributes':    ['dimension(:)', 'allocatable'],
                'variables':     [variable_name],
                'variableNames': [variable_name],
            }])

        code = ''

        multi_names = 'allowedMultiParameterNames_'
        if 'multiParameters' in directive:
            multi_parameter_names = [
                m.strip()
                for m in re.split(r'\s*,\s*', directive['multiParameters'])
                if m.strip()
            ]
            if not declaration_exists(parent, multi_names):
                add_declarations(parent, [{
                    'intrinsic':     'type',
                    'type':          'varying_string',
                    'openMP':        False,
                    'attributes':    [f'dimension({len(multi_parameter_names)})'],
                    'variables':     [multi_names],
                    'variableNames': [multi_names],
                }])
            for i, name in enumerate(multi_parameter_names, start=1):
                code += f"{multi_names}({i})='{name}'\n"

        add_uses(parent, {
            'moduleUse': {
                'ISO_Varying_String': {
                    'intrinsic': False,
                    'all':       True,
                },
            },
            'moduleOrder': ['ISO_Varying_String'],
        })

        if 'extraAllowedNames' in directive:
            extra_names = directive['extraAllowedNames'].split()
            code += f"allocate({variable_name}({len(extra_names)}))\n"
            for i, name in enumerate(extra_names, start=1):
                code += f"{variable_name}({i})='{name}'\n"

        if parent.get('type') != 'function':
            raise RuntimeError(
                "process_input_parameters_validate: parent is not a function")
        result = _function_class_result_name(parent)

        code += f"   call {result}%allowedParameters({variable_name},'{source}',.false.)\n"
        extra_multi = (f",allowedMultiParameterNames={multi_names}"
                       if 'multiParameters' in directive else '')
        code += (
            f"   if ({function_class_name}DsblVldtn == 0) call {source}"
            f"%checkParameters(allowedParameterNames={variable_name}{extra_multi})\n"
        )
        code += f"   if (allocated({variable_name})) deallocate({variable_name})\n"

        insert_after_node(node, [{
            'type':       'code',
            'content':    code,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':
                'Galacticus.Build.SourceTree.Process.InputParametersValidate.process_input_parameters_validate()',
            'line':       1,
        }])


register_process('inputParametersValidate', process_input_parameters_validate)
