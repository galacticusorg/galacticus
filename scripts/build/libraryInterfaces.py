#!/usr/bin/env python3
# Generates C/Fortran/Python interfaces to Galacticus library classes.
# Andrew Benson (ported to Python 2026)

import sys
import os
import re
from pathlib import Path

# Set up path for imports from python/
sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))

import xml.etree.ElementTree as ET
from Galacticus.Build import SourceTree
from Galacticus.Build.SourceTree.Parse import Declarations
from List.ExtraUtils import as_array, hash_list, sorted_keys
from Sort.Topo import sort as topo_sort

# Placeholder: full implementation follows the Perl script structure.
# This is skeleton code to be completed.

def main():
    """Main entry point — mirrors libraryInterfaces.pl."""

    # Initialize code and Python interface structures
    code = {'main': []}
    python = {'c_lib': []}

    # Load XML configuration files
    build_path = os.environ.get('BUILDPATH', './work/build')
    exec_path = os.environ['GALACTICUS_EXEC_PATH']

    directive_locations = _load_xml(os.path.join(build_path, 'directiveLocations.xml'))
    state_storables = _load_xml(os.path.join(build_path, 'stateStorables.xml'))
    library_classes = _load_xml(os.path.join(exec_path, 'source', 'libraryClasses.xml'))

    # Process function classes
    lib_function_classes = library_classes.get('classes', {})
    for class_name in sorted_keys(lib_function_classes):
        func_class = lib_function_classes[class_name]
        # Augment with module information
        class_key = class_name + 'Class'
        if 'functionClasses' in state_storables and class_key in state_storables['functionClasses']:
            func_class['module'] = state_storables['functionClasses'][class_key].get('module')

    # Process each function class
    if 'functionClass' in directive_locations and 'file' in directive_locations['functionClass']:
        for file_name in as_array(directive_locations['functionClass']['file']):
            _process_function_class_file(
                file_name, code, python, lib_function_classes,
                directive_locations, state_storables
            )

    # Write output files
    _write_fortran_code(code, build_path)
    _write_python_interface(python)


def _load_xml(path):
    """Load and parse XML file into nested dict structure."""
    if not os.path.exists(path):
        return {}
    try:
        root = ET.parse(path).getroot()
        return _xml_elem_to_dict(root)
    except ET.ParseError:
        return {}


def _xml_elem_to_dict(elem):
    """Convert XML element to nested dict (simple version)."""
    result = dict(elem.attrib)
    children_by_tag = {}

    for child in elem:
        tag = child.tag
        val = _xml_elem_to_dict(child)
        if tag in children_by_tag:
            if not isinstance(children_by_tag[tag], list):
                children_by_tag[tag] = [children_by_tag[tag]]
            children_by_tag[tag].append(val)
        else:
            children_by_tag[tag] = val

    result.update(children_by_tag)
    return result


def _process_function_class_file(file_name, code, python, lib_function_classes,
                                  directive_locations, state_storables):
    """Process a single file containing function class definitions."""
    # TODO: Port the main loop from libraryInterfaces.pl lines 43-198
    # This involves:
    # 1. Parsing the file with SourceTree.parse_file()
    # 2. Walking tree for functionClass directives and implementations
    # 3. Calling interfaces_*() generator functions
    pass


def _write_fortran_code(code, build_path):
    """Write generated Fortran code to files."""
    out_dir = os.path.join(build_path, 'libgalacticus')
    os.makedirs(out_dir, exist_ok=True)

    for class_name in sorted(code.keys()):
        if class_name == 'main':
            out_file = os.path.join(build_path, 'libgalacticus.Inc')
        else:
            out_file = os.path.join(out_dir, f'{class_name}.F90')

        with open(out_file, 'w') as fh:
            fh.write('\n'.join(code[class_name]) + '\n')


def _write_python_interface(python):
    """Write generated Python code to galacticus.py."""
    # TODO: Port topological sort and Python file generation
    # from libraryInterfaces.pl lines 250-273
    pass


if __name__ == '__main__':
    main()
