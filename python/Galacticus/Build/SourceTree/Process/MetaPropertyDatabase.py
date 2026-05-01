# Processes `metaPropertyDatabase` directives: walks every source file that
# declares an `addMetaProperty` directive and synthesizes a
# `subroutine metaPropertyNoCreator(…)` whose if/else-if ladder identifies
# which functionClass implementation would create a given meta-property.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/MetaPropertyDatabase.pm

import os
import re
import xml.etree.ElementTree as ET


from List.ExtraUtils                                         import as_array
from XML.Utils                                               import xml_to_dict
from Galacticus.Build.StateStorables                         import (
    function_class_names as _shared_function_class_names,
)
from Galacticus.Build.Directives                             import extract_directives
from Galacticus.Build.SourceTree                             import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process                     import register_process
from Galacticus.Build.SourceTree.Process.SourceIntrospection import location


def _load_xml(build_path, name):
    path = os.path.join(build_path, name)
    if not os.path.exists(path):
        raise RuntimeError(
            f"process_meta_property_database: {name} not found at {path}")
    return xml_to_dict(ET.parse(path).getroot())


def _function_class_names(state_storables):
    """Return the set of `<name>Class` keys advertised in stateStorables.xml."""
    return _shared_function_class_names(state_storables)


def _locate_function_class(file_name, class_names):
    """Return `(functionClassName, implementationName)` for the `*` directive
    in `file_name` whose root tag + 'Class' matches one of `class_names`.

    Mirrors MetaPropertyDatabase.pm:45-58.  Raises if the implementation's
    `name` attribute does not begin with the function class name (the Perl
    code `die`s in the same spot).
    """
    for directive in extract_directives(file_name, '*', set_root_element_type=True):
        root_type = directive.get('rootElementType')
        if (root_type or '') + 'Class' not in class_names:
            continue
        implementation_name = directive.get('name', '')
        if not implementation_name.startswith(root_type):
            raise RuntimeError(
                "process_meta_property_database: functionClass implementation "
                "name has incorrect prefix")
        # Perl: lcfirst(stripped)
        stripped = implementation_name[len(root_type):]
        implementation_name = stripped[:1].lower() + stripped[1:] if stripped else stripped
        return root_type, implementation_name
    return None, None


def process_meta_property_database(tree, options):
    """Mirrors Process_MetaPropertyDatabase() from MetaPropertyDatabase.pm."""
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        raise RuntimeError(
            "process_meta_property_database: BUILDPATH is not set")

    directive_locations = _load_xml(build_path, 'directiveLocations.xml')
    state_storables     = _load_xml(build_path, 'stateStorables.xml')
    class_names         = _function_class_names(state_storables)

    for node in walk_tree(tree):
        if node.get('type') != 'metaPropertyDatabase':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        directive['processed'] = True

        # Collect every <addMetaProperty isCreator="yes"/> from every file
        # declared under directiveLocations/addMetaProperty/file.
        amp_block = (directive_locations.get('addMetaProperty') or {})
        amp_files = list(as_array(amp_block.get('file')))
        creators = []
        for file_name in amp_files:
            file_creators = extract_directives(
                file_name, 'addMetaProperty',
                conditions={'isCreator': 'yes'},
            )
            if not file_creators:
                continue
            func_class, impl_name = _locate_function_class(file_name, class_names)
            for creator in file_creators:
                creator['functionClass']      = func_class
                creator['implementationName'] = impl_name
                creator.setdefault('type', 'float')
                creator.setdefault('rank', 0)
                creators.append(creator)

        loc_expr = location(node, node.get('line', 0))

        code  = "subroutine metaPropertyNoCreator(component_,name_,type_,rank_)\n"
        code += " use :: Error, only : Error_Report\n"
        code += " implicit none\n"
        code += " character(len=* ), intent(in   ) :: component_, name_             , type_\n"
        code += " integer          , intent(in   ) :: rank_\n"
        code += " character(len=64)                :: className , implementationName, rankLabel\n"
        code += "\n"
        code += " write (rankLabel,'(i1)') rank_\n"

        join = ""
        for creator in creators:
            name = creator.get('name', '')
            if name.startswith("'"):
                continue
            code += (
                join
                + f"if (component_ == '{creator.get('component', '')}'"
                  f" .and. name_ == '{name}'"
                  f" .and. type_ == '{creator.get('type', '')}'"
                  f" .and. rank_ == {creator.get('rank', 0)}) then\n"
            )
            code += f"          className='{creator.get('functionClass') or 'unknown'}'\n"
            code += f" implementationName='{creator.get('implementationName') or 'unknown'}'\n"
            join = "else "

        code += " else\n"
        code += "  className         =\"\"\n"
        code += "  implementationName=\"\"\n"
        code += (
            '  call Error_Report("no class creates the rank-"'
            '//trim(rankLabel)//" \'"//trim(type_)//"\' type meta-property \'"'
            '//trim(name_)//"\' in component \'"//trim(component_)//"\'"//'
            + loc_expr + ")\n"
        )
        code += " end if\n"
        code += (
            ' call Error_Report("the rank-"//trim(rankLabel)//" \'"//trim(type_)'
            '//"\' type meta-property \'"//trim(name_)//"\' in component \'"'
            '//trim(component_)//"\' is required"//char(10)'
            '//"it is created by the \'"//trim(implementationName)//"\' '
            'implementation of the \'"//trim(className)//"\' class"//char(10)'
            '//"to create this meta-property include the following in your '
            'parameter file:"//char(10)//" <"//trim(className)//" value="""'
            '//trim(implementationName)//"""/>"//' + loc_expr + ")\n"
        )
        code += " return\n"
        code += "end subroutine metaPropertyNoCreator\n"

        insert_after_node(node, [{
            'type':       'code',
            'content':    code,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':
                'Galacticus.Build.SourceTree.Process.MetaPropertyDatabase'
                '.process_meta_property_database()',
            'line':       1,
        }])


register_process('metaPropertyDatabase', process_meta_property_database)
