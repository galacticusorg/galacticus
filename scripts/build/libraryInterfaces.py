#!/usr/bin/env python3
# Generates C/Fortran/Python interfaces to Galacticus library classes.
# Andrew Benson (ported to Python with assistance from Claude 2026)

import sys
import os
import re
from pathlib import Path

# Set up path for imports from python/

import xml.etree.ElementTree as ET
from Galacticus.Build import SourceTree
from Galacticus.Build.SourceTree.Parse import Declarations
from List.ExtraUtils import as_array, hash_list, sorted_keys
from Sort.Topo import sort as topo_sort
from XML.Utils import xml_to_dict

from LibraryInterfaces.ArgSpec import ArgSpec
from LibraryInterfaces.Pipeline import (
    assign_c_types,
    assign_c_attributes,
    build_python_reassignments,
    build_fortran_reassignments,
)
from LibraryInterfaces.Emitters import (
    ctypes_arg_types,
    fortran_arg_list,
    fortran_declarations,
    fortran_reassignments,
    fortran_module_uses,
    fortran_call_code,
    iso_c_binding_import,
    python_arg_list,
    python_reassignments,
    python_call_code,
)
def main():
    """Main entry point — mirrors libraryInterfaces.pl."""

    # Initialize code and Python interface structures
    code = {'main': []}
    python = {'c_lib': [], 'units': {}}

    # Load XML configuration files
    build_path = os.environ.get('BUILDPATH', './work/build')
    exec_path = os.environ['GALACTICUS_EXEC_PATH']

    directive_locations = _load_xml(os.path.join(build_path, 'directiveLocations.xml'), required=True)
    state_storables = _load_xml(os.path.join(build_path, 'stateStorables.xml'), required=True)
    library_classes = _load_xml(os.path.join(exec_path, 'source', 'libraryClasses.xml'), required=True)

    # stateStorables.xml stores functionClasses as a flat list of elements each
    # carrying name= and module= attributes.  XML::Simple re-keys these by name
    # automatically; we replicate that here.
    fc_storables = {
        fc['name']: fc
        for fc in as_array(state_storables.get('functionClasses', []))
        if isinstance(fc, dict) and 'name' in fc
    }

    # Process function classes
    # libraryClasses.xml uses self-closing tags as class names, so _xml_elem_to_dict
    # returns '' for empty entries.  Normalize to dicts in-place so downstream code
    # can always call .get() on each entry and modifications persist.
    raw_classes = library_classes.get('classes', {})
    lib_function_classes = {}
    for name, val in raw_classes.items():
        fc = val if isinstance(val, dict) else {}
        fc['name'] = name
        lib_function_classes[name] = fc
        
    # Augment with module information from stateStorables.
    for name, fc in lib_function_classes.items():
        class_key = name + 'Class'
        if fc_storables.get(class_key):
            fc['module'] = fc_storables[class_key].get('module')

    # Process each function class
    directive_fc = directive_locations.get('functionClass', {})
    if directive_fc.get('file'):
        for file_name in as_array(directive_fc['file']):
            _process_function_class_file(
                file_name, code, python, lib_function_classes,
                directive_locations, state_storables
            )

    # Append initialization code
    _append_init_code(code)

    # Write output files
    _write_fortran_code(code, build_path)
    _write_python_interface(python)


def _load_xml(path, required=False):
    """Load and parse XML file into nested dict structure.

    If *required* is True, exits with a descriptive message when the file is
    missing or cannot be parsed rather than silently returning an empty dict.
    """
    if not os.path.exists(path):
        if required:
            sys.exit(f"libraryInterfaces.py: required XML file not found: {path}")
        return {}
    try:
        root = ET.parse(path).getroot()
        return xml_to_dict(root)
    except ET.ParseError as exc:
        if required:
            sys.exit(
                f"libraryInterfaces.py: failed to parse required XML file"
                f" '{path}': {exc}"
            )
        return {}




def _process_function_class_file(file_name, code, python, lib_function_classes,
                                  directive_locations, state_storables):
    """Process a single file containing function class definitions."""
    tree = SourceTree.parse_file(file_name)

    module_uses = []
    for node in SourceTree.walk_tree(tree):
        # Collect module uses
        if node['type'] == 'moduleUse':
            module_uses.append(node.get('moduleUse', {}))

        # Process functionClass directives
        if node['type'] != 'functionClass':
            continue

        directive = node.get('directive', {})
        class_name = directive.get('name')
        if not class_name or class_name not in lib_function_classes:
            continue

        func_class = lib_function_classes[class_name]

        # Augment methods
        methods = directive.get('method', {})

        if isinstance(methods, dict) and 'name' in methods:
            # Single method
            func_class['methods'] = {methods['name']: methods}
        elif isinstance(methods, list):
            # Multiple methods
            func_class['methods'] = {methods[i]['name']: method for i, method in enumerate(methods)}
        else:
            raise ValueError(
                f"incomprehensible 'method' field in '{file_name}' "
                f"for functionClass '{class_name}': expected a dict with a "
                f"'name' key or a list of such dicts, got "
                f"{type(methods).__name__}: {methods!r}"
            )

        func_class['moduleUses'] = module_uses

        # Parse implementations
        _process_implementations(
            func_class, directive_locations, state_storables,
            code, python, lib_function_classes
        )

def _unsupported_constructor_arg(args):
    """If any constructor argument has a type the pipeline can't translate,
    return ``(name, reason)``; otherwise ``None``.

    Today this catches:

    * ``complex`` / ``double complex`` — ctypes has no built-in c_complex.
    * Any source-level ``dimension(...)`` attribute.  Arrays of basic types
      need a size-passing convention we don't have yet, and arrays of
      ``character`` / ``type(varying_string)`` would also collide with the
      ``dimension(*)`` that ``assign_c_attributes`` auto-appends for
      c_char_p, producing invalid two-``dimension``-clause declarations.

    Allocatable scalars (no ``dimension``) are still permitted.
    """
    for arg in args:
        intrinsic = arg.get('intrinsic')
        if intrinsic in ('complex', 'double complex'):
            return arg['name'], f"{intrinsic}({arg.get('type','')})"
        for attr in arg.get('attributes', []):
            if attr.startswith('dimension'):
                return arg['name'], 'dimensioned argument'
    return None


def _process_implementations(func_class, directive_locations, state_storables,
                             code, python, lib_function_classes):
    """Process all implementations of a function class."""
    class_name = func_class['name']
    impls = directive_locations.get(class_name, {}).get('file', [])

    extensions = {}
    module_uses_impls = {}

    # First pass: collect extensions and module uses per implementation file.
    for impl_file in as_array(impls):
        tree = SourceTree.parse_file(impl_file)
        impl_name = None
        local_module_uses = []

        for node in SourceTree.walk_tree(tree):
            if node['type'] == class_name:
                impl_name = node.get('directive', {}).get('name')
            elif (impl_name and node['type'] == 'type'
                  and node.get('name') == impl_name):
                opener = node.get('opener', '')
                m = re.search(r',\s*extends\s*\(\s*([a-zA-Z0-9_]+)\s*\)', opener)
                if m:
                    extensions[node['name']] = m.group(1)
            elif node['type'] == 'moduleUse':
                module_use = node.get('moduleUse', {})
                if module_use:
                    local_module_uses.append(module_use)

        if impl_name:
            impl_conf = func_class.get(impl_name)
            is_excluded = (isinstance(impl_conf, dict)
                           and impl_conf.get('exclude') == 'yes')
            if not is_excluded:
                module_uses_impls[impl_name] = local_module_uses

    # Second pass: find constructors and build the implementations list.
    class_id = 0
    impls_list = []
    for impl_file in as_array(impls):
        tree = SourceTree.parse_file(impl_file)
        impl_name = None
        is_abstract = False
        name_constructor = None
        args_constructor = []

        for node in SourceTree.walk_tree(tree):
            if node['type'] == class_name:
                impl_name = node.get('directive', {}).get('name')
                is_abstract = (node.get('directive', {}).get('abstract', 'no') == 'yes')

            elif (impl_name
                  and node['type'] == 'interface'
                  and node.get('name', '').lower() == impl_name.lower()):
                # Scan moduleProcedure children for an Internal-suffixed constructor.
                child = node.get('firstChild')
                while child:
                    if child['type'] == 'moduleProcedure':
                        internal = [n for n in child.get('names', [])
                                    if n.lower().endswith('internal')]
                        if len(internal) == 1:
                            name_constructor = internal[0]
                    child = child.get('sibling')

            elif (name_constructor
                  and node['type'] == 'function'
                  and node.get('name', '').lower() == name_constructor.lower()):
                # Extract argument names from the function opener.
                opener = node.get('opener', '')
                m = re.search(
                    r'function\s+' + re.escape(name_constructor) + r'\s*\(([^)]+)\)',
                    opener, re.IGNORECASE)
                if m:
                    args_constructor = [{'name': a.strip()}
                                        for a in m.group(1).split(',')]
                # Enrich each argument with its declared type from child declaration nodes.
                child = node.get('firstChild')
                while child:
                    if child['type'] == 'declaration':
                        for decl in child.get('declarations', []):
                            for var_name in decl.get('variableNames', []):
                                for arg in args_constructor:
                                    if arg['name'].lower() == var_name.lower():
                                        arg['intrinsic'] = decl.get('intrinsic')
                                        arg['type'] = decl.get('type')
                                        arg['attributes'] = decl.get('attributes', [])
                    child = child.get('sibling')

        if impl_name is None:
            raise ValueError(
                f"Unable to find implementation of '{class_name}' in '{impl_file}'")

        # Fall back to impl_name as constructor if no Internal constructor was found.
        if not name_constructor:
            name_constructor = impl_name

        # classID is assigned to every file, even abstract/excluded ones (mirrors Perl).
        class_id += 1

        impl_conf = func_class.get(impl_name)
        is_excluded = (isinstance(impl_conf, dict)
                       and impl_conf.get('exclude') == 'yes')
        if not is_abstract and not is_excluded:
            unsupported = _unsupported_constructor_arg(args_constructor)
            if unsupported:
                # Skip implementations whose constructor takes an argument the
                # pipeline can't translate (today: source-level dimensioned
                # args of any intrinsic, and complex / double complex).
                # The class itself is still exposed; only the offending
                # implementation is omitted from impls_list, so it won't
                # appear in the Python class hierarchy or GetPtr dispatcher.
                sys.stderr.write(
                    f"libraryInterfaces.py: caution: implementation"
                    f" '{impl_name}' of class '{class_name}' has constructor"
                    f" argument '{unsupported[0]}' of unsupported kind"
                    f" ({unsupported[1]}) — skipping implementation\n"
                )
            else:
                impls_list.append({
                    'name':      impl_name,
                    'classID':   class_id,
                    'fileName':  impl_file,
                    'moduleUses': module_uses_impls.get(impl_name, []),
                    'arguments': args_constructor,
                })

    func_class['implementations'] = impls_list

    # Generate interfaces.
    interfaces_python_classes(python, func_class)
    interfaces_pointer_get(code, func_class)
    interfaces_constructors(code, python, func_class, lib_function_classes,
                           extensions, module_uses_impls)
    interfaces_methods(code, python, func_class, extensions, module_uses_impls,
                       lib_function_classes)
    interfaces_destructor(code, python, func_class)


def interfaces_pointer_get(code, func_class):
    """Generate pointer getter function for a function class."""
    class_name = func_class['name']
    impls = func_class.get('implementations', [])

    symbols = [class_name + 'Class'] + [impl['name'] for impl in impls]

    code.setdefault(class_name, []).append(f'''function {class_name}GetPtr({class_name}_,classID)
  use, intrinsic :: ISO_C_Binding, only : c_ptr, c_int, c_f_pointer
  use :: Error, only : Error_Report
  use :: {func_class.get('module', 'Unknown')}, only : {', '.join(symbols)}
  implicit none
  class({class_name}Class), pointer :: {class_name}GetPtr
  type(c_ptr), intent(in) :: {class_name}_
  integer(c_int), intent(in) :: classID
{chr(10).join(f'  type({impl["name"]}), pointer :: {impl["name"]}_' for impl in impls)}

  select case (classID)
{chr(10).join(f'  case ({impl["classID"]}){chr(10)}     call c_f_pointer({class_name}_, {impl["name"]}_){chr(10)}     {class_name}GetPtr => {impl["name"]}_' for impl in impls)}
  case default
     {class_name}GetPtr => null()
     call Error_Report('unknown classID'//{{introspection:location}})
  end select
  return
end function {class_name}GetPtr
''')


def interfaces_python_classes(python, func_class):
    """Generate Python class hierarchy for a function class."""
    class_name = func_class['name']

    # Parent class
    parent_code = f'''class {class_name}:

    # Constructor
    def __init__(self):
        # Assign class ID to negative (not a concrete class)
        self._classID = -1
'''
    python['units'][class_name] = {
        'content': parent_code,
        'indent': 0,
        'dependencies': ['init'],
    }

    # Child classes (implementations)
    for impl in func_class.get('implementations', []):
        child_code = f"class {impl['name']}({class_name}):"
        python['units'][impl['name']] = {
            'content': child_code,
            'indent': 0,
            'dependencies': [class_name],
        }


def interfaces_constructors(code, python, func_class, lib_function_classes,
                           extensions, module_uses_impls):
    """Generate constructor wrappers."""
    class_name = func_class['name']

    for impl in func_class.get('implementations', []):
        # Process argument list
        arg_list = impl.get('arguments', [])
        arg_list = assign_c_types(arg_list, lib_function_classes)
        arg_list = assign_c_attributes(arg_list)
        arg_list = build_python_reassignments(arg_list)
        arg_list = build_fortran_reassignments(
            arg_list, func_class, impl, extensions, module_uses_impls,
            lib_function_classes)

	# Construct pre- and post-arguments content for the call from Fortran to Galacticus.
        preArguments = f'''  !![
  <referenceConstruct object="self">
   <constructor>
    {impl.get('name', [])}( &amp;
'''
        postArguments = '''     &amp;                     )
   </constructor>
  </referenceConstruct>
  !!]
'''

        # Generate Fortran constructor
        iso_imports = iso_c_binding_import(arg_list, 'c_ptr', 'c_loc')
        fort_args = fortran_arg_list(arg_list)
        declarations = fortran_declarations(arg_list)
        reassignments = fortran_reassignments(arg_list)
        module_uses = fortran_module_uses(arg_list)
        call_code = fortran_call_code(arg_list, preArguments, postArguments, '&amp;')

        fort_constructor = f'''function {impl["name"]}L({','.join(fort_args)}) bind(c,name='{impl["name"]}L')
  use :: {func_class.get('module')}, only : {impl["name"]}
{iso_imports}
{module_uses}
  implicit none
  type(c_ptr) :: {impl["name"]}L
  type({impl["name"]}), pointer :: self
{declarations}

{reassignments}  allocate(self)
{call_code}  {impl["name"]}L=c_loc(self)
  return
end function {impl["name"]}L
'''
        code.setdefault(class_name, []).append(fort_constructor)

        # Add c_lib interface
        arg_types = ctypes_arg_types(arg_list)
        python['c_lib'].append({
            'name': impl['name'] + 'L',
            'restype': 'c_void_p',
            'argtypes': arg_types,
        })

        # Generate Python constructor
        py_args = python_arg_list(arg_list)
        py_reassignments = python_reassignments(arg_list)
        py_call = python_call_code(arg_list, 'self._glcObj = c_lib.' + impl['name'] + 'L')

        py_constructor = f'''# Constructor
def __init__({','.join(py_args)}):
    self._classID = {impl["classID"]}
{py_reassignments}{py_call}
'''
        python['units'][impl['name']].setdefault('subUnits', []).append({
            'content': py_constructor,
        })


def interfaces_methods(code, python, func_class, extensions, module_uses_impls,
                       lib_function_classes=None):
    """Generate method wrappers."""
    class_name = func_class['name']

    methods_to_delete = []
    for method_name, method_spec in func_class.get('methods', {}).items():
        method_type = method_spec.get('type', 'void')

        # Build argument list for method
        arg_list = [
            {
                'intrinsic': 'class',
                'type': class_name + 'Class',
                'attributes': ['intent(inout)'],
                'name': 'self',
            }
        ]

        # Generate method name (used below by some result-type conversions).
        method_name_c = class_name + method_name[0].upper() + method_name[1:] + 'L'

        # Determine any ISO_C_Binding imports needed.
        isoImports                = {}
        result_conversion_open    = ""
        result_conversion_close   = ""
        result_extra_module_uses  = ""
        result_extra_declarations = ""
        result_post_call_code     = ""
        result_call_target        = method_name_c
        result_python_decode      = False
        if method_type == "double precision":
            method_type_c                        = "real(c_double)"
            clib_res_type                        = "c_double";
            isoImports['c_double'] = 1
        elif method_type == "integer":
            method_type_c                        = "integer(c_int)"
            clib_res_type                        = "c_int";
            isoImports['c_int'] = 1
        elif method_type == "integer(c_long)":
            method_type_c                        = "integer(c_long)"
            clib_res_type                        = "c_long";
            isoImports['c_long'] = 1
        elif method_type == "integer(c_size_t)":
            method_type_c                        = "integer(c_size_t)"
            clib_res_type                        = "c_size_t";
            isoImports['c_size_t'] = 1
        elif method_type == "logical":
            method_type_c                        = "logical(c_bool)"
            clib_res_type                        = "c_bool";
            result_conversion_open               = "logical(";
            result_conversion_close              = ",kind=c_bool)";
            isoImports['c_bool'] = 1
        elif method_type == "type(varying_string)":
            # Returned varying_string is copied into a per-function static C
            # buffer (deallocated on each call so only the most recent result
            # persists), then a c_ptr to that buffer is returned.  Python's
            # ctypes c_char_p restype copies the bytes; we then decode to str.
            method_type_c                        = "type(c_ptr)"
            clib_res_type                        = "c_char_p"
            isoImports['c_ptr']       = 1
            isoImports['c_loc']       = 1
            isoImports['c_char']      = 1
            isoImports['c_null_char'] = 1
            result_extra_module_uses = (
                f'  use :: ISO_Varying_String, only : varying_string, char\n'
            )
            result_call_target = f'{method_name_c}_result_'
            result_extra_declarations = (
                f'  type     (varying_string)                                       :: {method_name_c}_result_\n'
                f'  character(kind=c_char   ), dimension(:), allocatable, save, target :: {method_name_c}_buffer_\n'
                f'  character(len=:         ), allocatable                          :: {method_name_c}_chars_\n'
                f'  integer                                                         :: {method_name_c}_i_\n'
            )
            result_post_call_code = (
                f'  {method_name_c}_chars_ = char({method_name_c}_result_)\n'
                f'  if (allocated({method_name_c}_buffer_)) deallocate({method_name_c}_buffer_)\n'
                f'  allocate({method_name_c}_buffer_(len({method_name_c}_chars_)+1))\n'
                f'  do {method_name_c}_i_ = 1, len({method_name_c}_chars_)\n'
                f'     {method_name_c}_buffer_({method_name_c}_i_) = {method_name_c}_chars_({method_name_c}_i_:{method_name_c}_i_)\n'
                f'  end do\n'
                f'  {method_name_c}_buffer_(len({method_name_c}_chars_)+1) = c_null_char\n'
                f'  {method_name_c} = c_loc({method_name_c}_buffer_)\n'
            )
            result_python_decode = True
        elif method_type == "void":
            pass
        else:
            sys.stderr.write(
                f"libraryInterfaces.py: caution: unsupported method return type"
                f" '{method_type}' in class '{class_name}', method"
                f" '{method_name}' — skipping\n"
            )
            methods_to_delete.append(method_name)
            continue

        # Add method arguments
        for arg_spec in as_array(method_spec.get('argument', [])):
            decl = Declarations.parse_declaration(arg_spec)
            if decl:
                for var_name in decl.get('variableNames', []):
                    arg_list.append({
                        'intrinsic': decl['intrinsic'],
                        'type': decl['type'],
                        'attributes': decl['attributes'],
                        'name': var_name,
                    })

        # Skip the method if any argument has a Fortran intrinsic the pipeline
        # can't translate.  Currently this is `complex` / `double complex`:
        # ctypes has no built-in c_complex, and the only known callers (e.g.
        # surveyGeometry::windowFunctions) combine complex with multi-dimensional
        # runtime-sized arrays — together a substantial design effort that
        # isn't justified for a single method.  Methods are dropped with a
        # warning rather than emitted as broken Fortran/Python.
        unsupported_arg = next(
            (arg for arg in arg_list[1:]
             if arg['intrinsic'] in ('complex', 'double complex')),
            None,
        )
        if unsupported_arg:
            sys.stderr.write(
                f"libraryInterfaces.py: caution: method '{method_name}' in"
                f" class '{class_name}' has argument '{unsupported_arg['name']}'"
                f" of unsupported type '{unsupported_arg['intrinsic']}"
                f"({unsupported_arg.get('type','')})' — skipping method\n"
            )
            methods_to_delete.append(method_name)
            continue

        # Process arguments
        arg_list = assign_c_types(arg_list, lib_function_classes or {})
        arg_list = assign_c_attributes(arg_list)
        arg_list = build_python_reassignments(arg_list)
        arg_list = build_fortran_reassignments(arg_list, func_class, None,
                                              extensions, module_uses_impls,
                                              lib_function_classes)

        # Generate Fortran method
        procedure = 'subroutine' if method_type == 'void' else 'function'
        func_decl = '' if method_type == 'void' else f'{method_type_c} :: {method_name_c}\n'

        iso_imports = iso_c_binding_import(arg_list, *isoImports.keys())
        fort_args = fortran_arg_list(arg_list)
        declarations = fortran_declarations(arg_list)
        reassignments = fortran_reassignments(arg_list)
        module_uses = fortran_module_uses(arg_list)
        call_lhs = "call" if method_type == "void" else f'{result_call_target} ='
        call_code = fortran_call_code(arg_list,
                                     f'{call_lhs} {result_conversion_open} self_%{method_name}( &\n',
                                     f'&){result_conversion_close}\n', '&')
        call_code += result_post_call_code

        fort_method = f'''{procedure} {method_name_c}({','.join(fort_args)}) bind(c,name='{method_name_c}')
  use :: {func_class.get('module')}, only : {class_name}Class
{module_uses}{result_extra_module_uses}{iso_imports}
  implicit none
{func_decl}{declarations}{result_extra_declarations}
{reassignments}{call_code}  return
end {procedure} {method_name_c}
'''
        code.setdefault(class_name, []).append(fort_method)

        # Add c_lib interface
        arg_types = ctypes_arg_types(arg_list)
        restype = None if method_type == 'void' else clib_res_type
        python['c_lib'].append({
            'name': method_name_c,
            'restype': restype,
            'argtypes': arg_types,
        })

        # Generate Python method
        py_args = python_arg_list(arg_list)
        py_call = python_call_code(arg_list, f'return c_lib.{method_name_c}')
        if result_python_decode:
            # Append .decode("utf-8") to each call line so the bytes returned by
            # ctypes c_char_p are converted to a Python str.
            py_call = re.sub(
                r'(c_lib\.' + re.escape(method_name_c) + r'\([^\n]*\))(\n)',
                r'\1.decode("utf-8")\2',
                py_call,
            )

        py_method = f'''def {method_name}({','.join(py_args)}):
{py_call}
'''
        python['units'][class_name].setdefault('subUnits', []).append({
            'content': py_method,
        })

    # Delete any unsupported methods.
    for key in methods_to_delete:
        del func_class['methods'][key]

def interfaces_destructor(code, python, func_class):
    """Generate destructor wrapper for a function class."""
    class_name = func_class['name']

    destructor_code = f'''subroutine {class_name}DestructorL(self,classID) bind(c,name='{class_name}DestructorL')
  use, intrinsic :: ISO_C_Binding, only : c_ptr, c_int
  use :: {func_class.get('module', 'Unknown')}, only : {class_name}Class
  implicit none
  type(c_ptr), value, intent(in) :: self
  integer(c_int), value, intent(in) :: classID
  class({class_name}Class), pointer :: self_, {class_name}GetPtr

  self_ => {class_name}GetPtr(self,classID)
  !![
  <objectDestructor name="self_"/>
  !!]
  return
end subroutine {class_name}DestructorL
'''

    code.setdefault(class_name, []).append(destructor_code)

    # Add c_lib interface
    python['c_lib'].append({
        'name': class_name + 'DestructorL',
        'restype': None,
        'argtypes': ['c_void_p', 'c_int'],
    })

    # Add Python destructor
    py_destructor = f'''# Destructor
def __del__(self):
    c_lib.{class_name}DestructorL(self._glcObj,self._classID)
'''
    python['units'].setdefault(class_name, {}).setdefault('subUnits', []).append({
        'content': py_destructor,
    })


def _append_init_code(code):
    """Append initialization code to the main code block."""
    init_code = '''subroutine libGalacticusInitL() bind(c,name='libGalacticusInitL')
  use:: Events_Hooks, only : eventsHooksInitialize
  use :: IO_HDF5, only : ioHDF5AccessInitialize

  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Initialize HDF5 library access lock.
  call ioHDF5AccessInitialize()
end subroutine libGalacticusInitL

program libGalacticusInit
end program libGalacticusInit
'''
    code['main'].append(init_code)


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
    # Initialize the init unit
    init_content = '''from ctypes import *
# Load the shared library into ctypes.
import os
cwd = os.getcwd()
libname = os.path.join(cwd, "galacticus/lib/libgalacticus.so")
c_lib = CDLL(libname)
c_lib.libGalacticusInitL()
'''

    # Add c_lib restype/argtypes
    for func in python['c_lib']:
        if func.get('restype'):
            init_content += f"c_lib.{func['name']}.restype = {func['restype']}\n"
        if func.get('argtypes'):
            argtypes_str = '[' + ', '.join(func['argtypes']) + ']'
            init_content += f"c_lib.{func['name']}.argtypes = {argtypes_str}\n"

    python['units']['init'] = {
        'content': init_content,
        'indent': 0,
    }

    # Topologically sort units
    dependencies = {}
    for unit_name, unit_data in python['units'].items():
        if 'dependencies' in unit_data:
            dependencies[unit_name] = unit_data['dependencies']

    unit_names = sorted(python['units'].keys())
    try:
        sorted_names = topo_sort(unit_names, dependencies)
    except RuntimeError:
        sorted_names = unit_names

    # Write galacticus.py
    with open('galacticus.py', 'w') as fh:
        stack = [python['units'][name] for name in sorted_names]
        while stack:
            unit = stack.pop(0)
            if 'subUnits' in unit:
                for subUnit in unit['subUnits']:
                    subUnit['indent'] = unit['indent']+1
                stack = unit['subUnits'] + stack
            indent = '    ' * unit.get('indent', 0)
            content = unit['content']
            for line in content.splitlines():
                fh.write(indent + line + '\n')
            fh.write('\n')


if __name__ == '__main__':
    main()

