# Processes `inputParameter` directives: emits the
# `source%value('name', var, defaultValue=…, writeOutput=…)` call that binds
# a Fortran variable to a named run-time parameter, and ensures the enclosing
# subprogram imports `Input_Parameters`.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/InputParameter.pm

import os
import re
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from XML.Utils                                    import xml_to_dict
from Galacticus.Build.SourceTree                  import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process          import register_process
from Galacticus.Build.SourceTree.Parse.ModuleUses import add_uses


def _read_executables(build_path):
    """Return the list of non-test executables listed on `all_exes = …` of
    `$BUILDPATH/Makefile_All_Execs`.  Mirrors InputParameter.pm:25-34.
    """
    path = os.path.join(build_path, 'Makefile_All_Execs')
    executables = []
    with open(path, 'r', errors='replace') as fh:
        for line in fh:
            m = re.match(r'^all_exes\s+=\s+(.*)', line)
            if m:
                executables = [
                    token for token in m.group(1).split()
                    if not token.startswith('tests.')
                ]
                break
    return executables


def _read_dependencies(build_path, executables):
    """For each executable, parse `<exe>.d` and collect its object-file stems.

    Mirrors InputParameter.pm:36-48.  Missing `.d` files are silently skipped
    (Perl's `-e` guard).  The returned shape — `{exe: {stem: True}}` — matches
    Perl's `$dependencies->{$exe}->{$stem} = 1` idiom.
    """
    dependencies = {}
    for executable in executables:
        dep_name = re.sub(r'\.exe$', '.d', executable)
        dep_path = os.path.join(build_path, dep_name)
        if not os.path.exists(dep_path):
            continue
        dependencies[executable] = {}
        with open(dep_path, 'r', errors='replace') as fh:
            for line in fh:
                m = re.search(r'([^/]+)\.o$', line.rstrip('\n'))
                if m:
                    dependencies[executable][m.group(1)] = True
    return dependencies


def process_input_parameters(tree, options):
    """Mirrors Process_InputParameters() from InputParameter.pm."""
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        raise RuntimeError(
            "process_input_parameters: BUILDPATH is not set")

    # The executables/dependencies/directiveLocations reads are performed even
    # when no inputParameter directives exist — matches Perl's eager load.
    executables  = _read_executables(build_path)
    _dependencies = _read_dependencies(build_path, executables)
    _directive_locations_path = os.path.join(build_path, 'directiveLocations.xml')
    if os.path.exists(_directive_locations_path):
        xml_to_dict(ET.parse(_directive_locations_path).getroot())

    for node in walk_tree(tree):
        if node.get('type') != 'inputParameter':
            continue
        directive = node.get('directive') or {}
        if directive.get('processed'):
            continue
        directive['processed'] = True

        code  = "  ! Auto-generated input parameter\n"
        if 'name' in directive:
            source        = directive['source']
            parameter     = directive['name']
            # Add delimiters unless the name is actually a function call.
            parameter_arg = parameter if '(' in parameter else f"'{parameter}'"
            variable      = directive.get('variable', directive['name'])
            code += f"  call {source}%value({parameter_arg},{variable}"
            if 'defaultValue' in directive:
                code += f",defaultValue={directive['defaultValue']}"
            if 'writeOutput' in directive:
                write_flag = ('.false.' if directive['writeOutput'] == 'no'
                              else '.true.')
                code += f",writeOutput={write_flag}"
            code += ")\n"
        code += "  ! End auto-generated input parameter\n\n"

        insert_after_node(node, [{
            'type':       'code',
            'content':    code,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':
                'Galacticus.Build.SourceTree.Process.InputParameter.process_input_parameters()',
            'line':       1,
        }])

        add_uses(node['parent'], {
            'moduleUse': {
                'Input_Parameters': {
                    'intrinsic': False,
                    'all':       True,
                },
            },
            'moduleOrder': ['Input_Parameters'],
        })


register_process('inputParameters', process_input_parameters)
