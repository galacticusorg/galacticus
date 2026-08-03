"""Processes `functionsGlobal` directives: either synthesizes the
procedure-pointer declarations + `_Null` default implementations
(`type="pointers"`) or emits the runtime pointer assignments
(`type="establish"`) that wire named global functions into the pointer
table.

Andrew Benson (ported to Python 2026)
"""

import os
import re
import xml.etree.ElementTree as ET


from List.ExtraUtils                                 import as_array
from XML.Utils                                       import xml_to_dict
from Fortran.Utils                                   import extract_variables
from Galacticus.Build.StateStorables                 import (
    function_class_module_map as _shared_function_class_module_map,
)
from Galacticus.Build.Directives                     import extract_directives
from Galacticus.Build.SourceTree                     import (
    walk_tree, parse_file, parse_code,
    insert_after_node, insert_pre_contains, insert_post_contains,
)
from Galacticus.Build.SourceTree.Process             import register_process, process_tree
from Galacticus.Build.SourceTree.Parse.ModuleUses    import add_uses


_DIRECTIVE_LOCATIONS = None
_STATE_STORABLES     = None


def _load_xml_once():
    """Load directiveLocations.xml + stateStorables.xml once, caching at
    module scope.
    """
    global _DIRECTIVE_LOCATIONS, _STATE_STORABLES
    if _DIRECTIVE_LOCATIONS is not None and _STATE_STORABLES is not None:
        return _DIRECTIVE_LOCATIONS, _STATE_STORABLES
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        raise RuntimeError(
            "process_functions_global: BUILDPATH is not set")
    if _DIRECTIVE_LOCATIONS is None:
        path = os.path.join(build_path, 'directiveLocations.xml')
        _DIRECTIVE_LOCATIONS = xml_to_dict(ET.parse(path).getroot())
    if _STATE_STORABLES is None:
        path = os.path.join(build_path, 'stateStorables.xml')
        _STATE_STORABLES = xml_to_dict(ET.parse(path).getroot())
    return _DIRECTIVE_LOCATIONS, _STATE_STORABLES


def _collect_function_globals(directive_locations):
    """Walk every file listed under directiveLocations/functionGlobal/file,
    extracting every `<functionGlobal>` directive and tagging each with the
    file it came from.
    """
    block = directive_locations.get('functionGlobal') or {}
    files = list(as_array(block.get('file')))
    collected = []
    for f in files:
        for d in extract_directives(f, 'functionGlobal'):
            d['file'] = f
            collected.append(d)
    return collected


def _class_module_map(state_storables):
    """Return `{<root>+'Class': <module>}` for every functionClass entry so the
    establish path can look up which module provides the named class.
    """
    return _shared_function_class_module_map(state_storables)


def _argument_names(arguments):
    """Extract the declared variable names from every `<arguments>::` entry."""
    names = []
    for arg in as_array(arguments):
        m = re.search(r'::(.*)$', arg)
        if m:
            names.extend(extract_variables(m.group(1)))
    return names


def _render_null_definitions(function_globals):
    """Render the `_Null` subprogram bodies for every functionGlobal."""
    out = ""
    for fg in function_globals:
        void = fg.get('type') == 'void'
        opener = 'subroutine' if void else 'function'
        names  = _argument_names(fg.get('arguments'))
        out += (
            f" {opener} {fg['unitName']}_Null(" + ",".join(names) + ")\n"
            "  use :: Error, only : Error_Report\n"
        )
        for module_entry in as_array(fg.get('module')):
            module_name = re.sub(r'\s*,.*', '', module_entry)
            intrinsic_tag = ', intrinsic' if module_name == 'ISO_C_Binding' else ''
            out += f"use{intrinsic_tag} :: {module_entry}\n"
        if not void:
            out += f"  {fg['type']} :: {fg['unitName']}_Null\n"
        for arg in as_array(fg.get('arguments')):
            out += f"  {arg}\n"
        if names:
            out += "  !$GLC attributes unused :: " + ",".join(names) + "\n"
        if fg.get('type') == 'double precision':
            out += f"  {fg['unitName']}_Null=0.0d0\n"
        elif fg.get('type') and re.search(r',\s*pointer', fg['type']):
            out += f"  {fg['unitName']}_Null => null()\n"
        out += (
            "  call Error_Report('global functions have not been "
            "initialized'//{introspection:location})\n"
        )
        out += f" end {opener} {fg['unitName']}_Null\n"
    return out


def _collect_module_uses_for_establish(function_globals, class_module_map):
    """Build the `moduleUse` dict required by the `establish` path.

    For each functionGlobal, parse the
    containing file, and for every module at file-scope (or every
    functionClass whose class module is recorded in stateStorables.xml)
    add a `use <module>, only : <unitName>` entry.
    """
    module_uses = {}
    for fg in function_globals:
        tree = parse_file(fg['file'])
        child = tree.get('firstChild')
        while child is not None:
            if child.get('type') == 'module':
                name = child.get('name')
                if name:
                    entry = module_uses.setdefault(
                        name, {'intrinsic': False})
                    entry.setdefault('only', {})[fg['unitName']] = True
            # functionClass instances declared as sibling types pull in the
            # module that provides the class itself.
            class_key = (child.get('type') or '') + 'Class'
            if class_key in class_module_map:
                provider = class_module_map[class_key]
                if provider:
                    entry = module_uses.setdefault(
                        provider, {'intrinsic': False})
                    entry.setdefault('only', {})[fg['unitName']] = True
            child = child.get('sibling')
    return module_uses


def process_functions_global(tree, options):
    """Process `functionsGlobal` directives in the tree."""
    directive_locations, state_storables = _load_xml_once()
    class_module_map = _class_module_map(state_storables)

    for node in walk_tree(tree):
        if node.get('type') != 'functionsGlobal':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        directive['processed'] = True

        function_globals = _collect_function_globals(directive_locations)
        kind = directive.get('type')

        if kind == 'pointers':
            # Procedure-pointer declarations, inserted before `contains`.
            pointer_lines = "\n".join(
                f"    procedure({fg['unitName']}_Null), pointer :: "
                f"{fg['unitName']}_ => {fg['unitName']}_Null"
                for fg in function_globals
            ) + "\n"
            insert_pre_contains(node['parent'], [{
                'type':       'code',
                'content':    pointer_lines,
                'parent':     None,
                'firstChild': None,
                'sibling':    None,
                'source':
                    'Galacticus.Build.SourceTree.Process.FunctionsGlobal'
                    '.process_functions_global()',
                'line':       1,
            }])

            # Null-default subprograms — synthesized, re-parsed so nested
            # directives (e.g. {introspection:location}) are picked up by
            # downstream passes, then inserted after `contains`.
            definitions = _render_null_definitions(function_globals)
            sub_tree = parse_code(
                definitions,
                name='Galacticus.Build.SourceTree.Process.FunctionsGlobal'
                     '.process_functions_global()',
            )
            process_tree(sub_tree)
            insert_post_contains(node['parent'], [sub_tree])

        elif kind == 'establish':
            # Pointer assignments at the directive site, plus module-use
            # injection so every callee in the assignments is visible.
            assignments = "\n".join(
                f"    {fg['unitName']}_ => {fg['unitName']}"
                for fg in function_globals
            ) + "\n"
            module_uses = _collect_module_uses_for_establish(
                function_globals, class_module_map)
            if module_uses:
                add_uses(node['parent'], {
                    'moduleUse':   module_uses,
                    'moduleOrder': sorted(module_uses.keys()),
                })
            insert_after_node(node, [{
                'type':       'code',
                'content':    assignments,
                'parent':     None,
                'firstChild': None,
                'sibling':    None,
                'source':
                    'Galacticus.Build.SourceTree.Process.FunctionsGlobal'
                    '.process_functions_global()',
                'line':       1,
            }])
        else:
            raise RuntimeError(
                f"process_functions_global: unknown type '{kind}'")


register_process('functionsGlobal', process_functions_global)
