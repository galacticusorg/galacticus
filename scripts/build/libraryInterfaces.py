#!/usr/bin/env python3
# Generates C/Fortran/Python interfaces to Galacticus library classes.
# Andrew Benson (ported to Python with assistance from Claude 2026)

import sys
import os
import re
from pathlib import Path
from itertools import chain, combinations

# Set up path for imports from python/
sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))

import xml.etree.ElementTree as ET
from Galacticus.Build import SourceTree
from Galacticus.Build.SourceTree.Parse import Declarations
from List.ExtraUtils import as_array, hash_list, sorted_keys
from Sort.Topo import sort as topo_sort


def main():
    """Main entry point — mirrors libraryInterfaces.pl."""

    # Initialize code and Python interface structures
    code = {'main': []}
    python = {'c_lib': [], 'units': {}}

    # Load XML configuration files
    build_path = os.environ.get('BUILDPATH', './work/build')
    exec_path = os.environ['GALACTICUS_EXEC_PATH']

    directive_locations = _load_xml(os.path.join(build_path, 'directiveLocations.xml'))
    state_storables = _load_xml(os.path.join(build_path, 'stateStorables.xml'))
    library_classes = _load_xml(os.path.join(exec_path, 'source', 'libraryClasses.xml'))

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

def _powerset(iterable):
    """powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"""
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

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
    """Convert XML element to nested dict (mirrors XML::Simple behaviour).

    Text-only elements (no attributes, no children) return the text string
    directly, matching XML::Simple's default behaviour.
    """
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
    # If element has no attributes and no children, return text content directly.
    # This mirrors XML::Simple which returns the text string for such elements.
    if not result:
        return (elem.text or '').strip()
    return result


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
            raise ValueError("incomprehensible methods")

        func_class['moduleUses'] = module_uses

        # Parse implementations
        _process_implementations(
            func_class, directive_locations, state_storables,
            code, python, lib_function_classes
        )

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

        # Determine any ISO_C_Binding imports needed.
        isoImports = {}
        result_conversion_open = ""
        result_conversion_close = ""
        if method_type == "double precision":
            method_type_c                        = "real(c_double)"
            clib_res_type                        = "c_double";
            isoImports['c_double'] = 1
        elif method_type == "logical":
            method_type_c                        = "logical(c_bool)"
            clib_res_type                        = "c_bool";
            result_conversion_open               = "logical(";
            result_conversion_close              = ",kind=c_bool)";
            isoImports['c_bool'] = 1
        elif method_type == "void":
            pass
        else:
            print(f"unsupported type '{method_type}'")
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

        # Process arguments
        arg_list = assign_c_types(arg_list, lib_function_classes or {})
        arg_list = assign_c_attributes(arg_list)
        arg_list = build_python_reassignments(arg_list)
        arg_list = build_fortran_reassignments(arg_list, func_class, None,
                                              extensions, module_uses_impls,
                                              lib_function_classes)

        # Generate method name
        method_name_c = class_name + method_name[0].upper() + method_name[1:] + 'L'

        # Generate Fortran method
        procedure = 'subroutine' if method_type == 'void' else 'function'
        func_decl = '' if method_type == 'void' else f'{method_type_c} :: {method_name_c}\n'

        iso_imports = iso_c_binding_import(arg_list, *isoImports.keys())
        fort_args = fortran_arg_list(arg_list)
        declarations = fortran_declarations(arg_list)
        reassignments = fortran_reassignments(arg_list)
        module_uses = fortran_module_uses(arg_list)
        call_code = fortran_call_code(arg_list,
                                     f'{"call" if method_type == "void" else method_name_c + "="} {result_conversion_open} self_%{method_name}( &\n',
                                     f'&){result_conversion_close}\n', '&')

        fort_method = f'''{procedure} {method_name_c}({','.join(fort_args)}) bind(c,name='{method_name_c}')
  use :: {func_class.get('module')}, only : {class_name}Class
{module_uses}
{iso_imports}
  implicit none
{func_decl}{declarations}
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


def assign_c_types(argument_list, lib_function_classes):
    """Assign appropriate C types for each argument.

    Mirrors Perl assignCTypes().  Processes the list in reverse and builds a
    new one so that the _ID companion argument for functionClass parameters
    can be inserted immediately after its parent without disturbing the rest
    of the list.
    """
    new_list = []
    for arg in reversed(argument_list):
        arg.setdefault('ctypes', {})
        arg.setdefault('fortran', {})
        arg.setdefault('python', {})
        arg.setdefault('galacticus', {})

        arg['isOptional'] = bool('optional' in arg.get('attributes', []))
        arg['fortran']['isPresent'] = 1
        arg['python']['isPresent'] = 1
        arg['galacticus']['isPresent'] = 1
        arg['isFunctionClass'] = 0

        intrinsic  = arg.get('intrinsic', '')
        type_spec  = arg.get('type') or ''

        if intrinsic == 'double precision':
            arg['ctypes']['type'] = 'c_double'
            arg['fortran']['type'] = 'real(c_double)'
        elif intrinsic == 'integer':
            arg['ctypes']['type'] = 'c_int'
            arg['fortran']['type'] = 'integer(c_int)'
        elif intrinsic == 'logical':
            arg['ctypes']['type'] = 'c_bool'
            arg['fortran']['type'] = 'logical(c_bool)'
        elif intrinsic == 'character':
            arg['ctypes']['type'] = 'c_char_p'
            arg['fortran']['type'] = 'character(c_char)'
        elif intrinsic == 'type':
            if type_spec == 'varying_string':
                arg['ctypes']['type'] = 'c_char_p'
                arg['fortran']['type'] = 'character(c_char)'
            elif re.match(r'^enumeration[a-z0-9_]+type$', type_spec, re.IGNORECASE):
                # Enumeration types map to C int.
                arg['ctypes']['type'] = 'c_int'
                arg['fortran']['type'] = 'integer(c_int)'
            else:
                arg['ctypes']['type'] = 'c_void_p'
                arg['fortran']['type'] = 'type(c_ptr)'
        elif intrinsic == 'class':
            arg['ctypes']['type'] = 'c_void_p'
            arg['fortran']['type'] = 'type(c_ptr)'
            # Check whether this is a functionClass argument.
            if type_spec.endswith('Class'):
                class_key = type_spec[:-5]  # strip trailing 'Class'
                if class_key in lib_function_classes:
                    arg['isFunctionClass'] = 1
                    # 'self' is dispatched via the method binding, not passed directly.
                    if arg['name'] == 'self':
                        arg['galacticus']['isPresent'] = 0
                    # Insert a companion _ID argument (carries the concrete class ID).
                    arg_id = {
                        'name':       arg['name'] + '_ID',
                        'intrinsic':  'integer',
                        'type':       None,
                        'attributes': ['intent(in)'],
                        'ctypes':     {'type': 'c_int'},
                        'fortran':    {'type': 'integer(c_int)', 'isPresent': 1},
                        'python':     {'isPresent': 0},
                        'galacticus': {'isPresent': 0},
                        'isOptional': False,
                        'isFunctionClass': False,
                    }
                    if arg['isOptional']:
                        arg_id['attributes'].append('optional')
                        arg_id['isOptional'] = True
                        arg_id['python']['present'] = arg['name']
                    # Insert _ID before the current front, then arg before that.
                    new_list.insert(0, arg_id)

        new_list.insert(0, arg)

    return new_list


def assign_c_attributes(argument_list):
    """Assign C attributes to arguments."""
    for arg in argument_list:
        arg.setdefault('fortran', {})
        arg['fortran']['attributes'] = []

        if arg.get('isOptional'):
            arg['fortran']['attributes'].append('optional')

        attr_filters = [a for a in arg.get('attributes', [])
                       if a.startswith('dimension') or a == 'allocatable']
        arg['fortran']['attributes'].extend(attr_filters)

        if arg.get('ctypes', {}).get('type') == 'c_char_p':
            arg['fortran']['attributes'].append('dimension(*)')

        # Determine pass-by method
        is_ptr_type = arg.get('ctypes', {}).get('type', '').endswith('_p')
        is_intent_in = any('intent(in)' in a for a in arg.get('attributes', []))
        is_non_scalar = any(a.startswith('dimension') for a in arg.get('fortran', {}).get('attributes', []))

        if (is_ptr_type or is_intent_in) and not arg.get('isOptional') and not is_non_scalar:
            arg['passBy'] = 'value'
        else:
            arg['passBy'] = 'reference'

        if arg['passBy'] == 'value':
            arg['fortran']['attributes'].append('value')

        arg['ctypes']['pointer'] = (arg['passBy'] == 'reference' and
                                   arg.get('ctypes', {}).get('type') != 'c_char_p')

    return argument_list


def build_python_reassignments(argument_list):
    """Set python['passAs'] and python['reassignment'] for functionClass args.

    Mirrors Perl buildPythonReassignments().  Processes in reverse so that when
    a functionClass arg is encountered, its _ID companion is already sitting at
    the front of new_list (it was the immediately preceding arg in forward
    order, so the last one pushed in reverse order).

    Non-optional:  passAs = 'name._glcObj' / 'name._classID'
    Optional:      passAs = 'name_glcObj'  / 'name_classID' plus an
                   if/else reassignment block that extracts the values or
                   sets them to None when the argument is absent.
    """
    new_list = []
    for arg in reversed(argument_list):
        if arg.get('isFunctionClass'):
            arg_id = new_list.pop(0)          # shift _ID off front of new list
            name = arg['name']
            if arg.get('isOptional'):
                arg['python']['passAs']        = name + '_glcObj'
                arg_id['python']['passAs']     = name + '_classID'
                arg['python']['reassignment']  = (
                    f'    if {name}:\n'
                    f'        {name}_glcObj ={name}._glcObj\n'
                    f'        {name}_classID={name}._classID\n'
                    f'    else:\n'
                    f'        {name}_glcObj =None\n'
                    f'        {name}_classID=None\n'
                )
            else:
                arg['python']['passAs']    = name + '._glcObj'
                arg_id['python']['passAs'] = name + '._classID'
            new_list.insert(0, arg_id)        # unshift _ID back
        new_list.insert(0, arg)               # unshift current arg
    return new_list


def build_fortran_reassignments(argument_list, func_class, implementation,
                               extensions, module_uses_impls,
                               lib_function_classes=None):
    """Generate Fortran reassignments for cross-language type conversions.

    Mirrors Perl buildFortranReassignments().  Processes in reverse (same
    pop/unshift skeleton as the other builders); no new args are inserted so
    the order is unchanged.

    Cases handled:
      logical       — c_bool  → logical via logical() cast
      character     — c_char_p → varying_string via char(String_C_to_Fortran())
      varying_string — c_char_p → varying_string via String_C_to_Fortran()
      enumeration   — c_int   → type(enumXxx) via %ID assignment; module located
                                by walking implementation/extension/functionClass uses
      treeNode      — c_ptr   → type(treeNode) via c_f_pointer
      other types   — c_ptr   → type(X) via c_f_pointer
      isFunctionClass — c_ptr → class(XClass) via XGetPtr(ptr, ID)
      class(*)      — c_ptr   → concrete type via c_f_pointer (type from libraryClasses config)
    """
    if lib_function_classes is None:
        lib_function_classes = {}

    new_list = []
    for arg in reversed(argument_list):
        arg.setdefault('fortran', {})
        intrinsic  = arg.get('intrinsic', '')
        type_spec  = arg.get('type') or ''
        name       = arg['name']
        is_optional = arg.get('isOptional', False)
        opt_prefix  = f'if (present({name})) ' if is_optional else ''

        if intrinsic == 'logical':
            # c_bool must be recast to a plain Fortran logical.
            arg['fortran']['reassignment'] = f'{opt_prefix}{name}_=logical({name})\n'
            arg['fortran']['declarations'] = f'logical :: {name}_\n'
            arg['fortran']['passAs']       = name + '_'

        elif intrinsic == 'character':
            # c_char_p → character via String_C_to_Fortran then char()
            arg['fortran'].setdefault('modules', {})
            arg['fortran']['modules'].setdefault('String_Handling',    {})['String_C_to_Fortran'] = 1
            arg['fortran']['modules'].setdefault('ISO_Varying_String', {})['char']                = 1
            arg['fortran']['passAs'] = f'char(String_C_to_Fortran({name}))'

        elif intrinsic == 'type':
            if type_spec == 'varying_string':
                # c_char_p → varying_string via String_C_to_Fortran
                arg['fortran'].setdefault('modules', {})
                arg['fortran']['modules'].setdefault('String_Handling', {})['String_C_to_Fortran'] = 1
                arg['fortran']['passAs'] = f'String_C_to_Fortran({name})'

            elif re.match(r'^enumeration[a-z0-9_]+type$', type_spec, re.IGNORECASE):
                # c_int → type(enumXxx) via %ID assignment.
                arg['fortran']['declarations'] = f'type({type_spec}) :: {name}_\n'
                arg['fortran']['passAs']       = name + '_'
                arg['fortran']['reassignment'] = f'{opt_prefix}{name}_%ID={name}\n'
                # Locate the module that imports this enumeration type:
                # 1. walk implementation's module uses, following the extends chain.
                import_module = None
                if implementation:
                    cls = implementation['name']
                    while cls and not import_module:
                        for use_block in module_uses_impls.get(cls, []):
                            for mod_name, mod_data in use_block.items():
                                if (isinstance(mod_data, dict)
                                        and type_spec in mod_data.get('only', {})):
                                    import_module = mod_name
                                    break
                            if import_module:
                                break
                        cls = extensions.get(cls)
                # 2. fall back to the functionClass file's own module uses.
                if not import_module:
                    for use_block in func_class.get('moduleUses', []):
                        for mod_name, mod_data in use_block.items():
                            if (isinstance(mod_data, dict)
                                    and type_spec in mod_data.get('only', {})):
                                import_module = mod_name
                                break
                        if import_module:
                            break
                # 3. last resort: the functionClass's own module.
                if not import_module:
                    import_module = func_class.get('module')
                if import_module:
                    arg['fortran'].setdefault('modules', {})
                    arg['fortran']['modules'].setdefault(import_module, {})[type_spec] = 1

            elif type_spec == 'treeNode':
                # c_ptr → type(treeNode) via c_f_pointer.
                arg['fortran']['declarations'] = f'type({type_spec}), pointer :: {name}_\n'
                arg['fortran']['passAs']       = name + '_'
                reassign = f'call c_f_pointer({name},{name}_)\n'
                if is_optional:
                    reassign = (f'if (present({name})) then\n '
                                f'{reassign}else\n {name}_ => null()\nend if\n')
                arg['fortran']['reassignment']       = reassign
                arg['fortran']['isoCBindingSymbols'] = ['c_f_pointer']
                arg['fortran'].setdefault('modules', {})
                arg['fortran']['modules'].setdefault('Galacticus_Nodes', {})['treeNode'] = 1

            else:
                # Other derived types: c_ptr → type(X) via c_f_pointer.
                arg['fortran']['declarations'] = f'type({type_spec}), pointer :: {name}_\n'
                arg['fortran']['passAs']       = name + '_'
                arg['fortran']['reassignment'] = f'{opt_prefix}call c_f_pointer({name},{name}_)\n'
                arg['fortran']['isoCBindingSymbols'] = ['c_f_pointer']
                arg['fortran'].setdefault('modules', {})
                if type_spec == 'inputParameters':
                    arg['fortran']['modules'].setdefault('Input_Parameters', {})['inputParameters'] = 1
                else:
                    mod = func_class.get('module', '')
                    if mod:
                        arg['fortran']['modules'].setdefault(mod, {})[type_spec] = 1

        elif arg.get('isFunctionClass'):
            # class(FooClass) → class(FooClass) pointer via FooGetPtr(ptr, ID).
            class_key = type_spec[:-5] if type_spec.endswith('Class') else type_spec
            mod_name  = lib_function_classes.get(class_key, {}).get('module', '')
            arg['fortran']['declarations'] = f'class({type_spec}), pointer :: {name}_\n'
            arg['fortran']['passAs']       = name + '_'
            arg['fortran']['reassignment'] = (
                f'{opt_prefix}{name}_ => {class_key}GetPtr({name},{name}_ID)\n'
            )
            arg['fortran']['functionClass'] = class_key
            if mod_name:
                arg['fortran'].setdefault('modules', {})
                arg['fortran']['modules'].setdefault(mod_name, {})[type_spec] = 1

        elif intrinsic == 'class' and type_spec == '*':
            # Unlimited polymorphic: look up the concrete type in libraryClasses config.
            impl_name = implementation['name'] if implementation else None
            fc_name   = func_class.get('name', '')
            impl_info = (lib_function_classes.get(fc_name, {}).get(impl_name, {})
                         if impl_name else {})
            ctor_args = as_array(impl_info.get('constructor', {}).get('argument', []))
            concrete  = [a for a in ctor_args
                         if isinstance(a, dict) and a.get('name') == name]
            if concrete:
                ct      = concrete[0]
                ct_type = ct.get('type', '')
                ct_mod  = ct.get('module', '')
                arg['fortran']['declarations']       = f'type({ct_type}), pointer :: {name}_\n'
                arg['fortran']['passAs']             = name + '_'
                arg['fortran']['reassignment']       = f'{opt_prefix}call c_f_pointer({name},{name}_)\n'
                arg['fortran']['isoCBindingSymbols'] = ['c_f_pointer']
                if ct_mod and ct_type:
                    arg['fortran'].setdefault('modules', {})
                    arg['fortran']['modules'].setdefault(ct_mod, {})[ct_type] = 1

        new_list.insert(0, arg)

    return new_list


def ctypes_arg_types(argument_list):
    """Generate ctypes argument type list."""
    types = []
    for arg in argument_list:
        ctype = arg.get('ctypes', {}).get('type', 'c_int')
        if arg.get('ctypes', {}).get('pointer'):
            ctype = f'POINTER({ctype})'
        types.append(ctype)
    return types


def fortran_arg_list(argument_list):
    """Generate Fortran argument list."""
    return [arg['name'] for arg in argument_list
            if arg.get('fortran', {}).get('isPresent')]


def fortran_declarations(argument_list):
    """Generate Fortran declarations.

    Mirrors Perl fortranDeclarations(): emits one declaration per argument
    plus any extra declarations stored in fortran['declarations'].  Also
    emits an interface block for every distinct functionClass GetPtr function
    that is referenced by an isFunctionClass argument.
    """
    code = ''
    function_classes = {}   # {className: True}  — deduplicates interface blocks

    for arg in argument_list:
        attrs    = arg.get('fortran', {}).get('attributes', [])
        attr_str = (', ' + ', '.join(attrs)) if attrs else ''
        fort_type = arg.get('fortran', {}).get('type', 'integer(c_int)')
        code += f'  {fort_type}{attr_str} :: {arg["name"]}\n'
        if arg.get('fortran', {}).get('declarations'):
            code += arg['fortran']['declarations']
        if arg.get('fortran', {}).get('functionClass'):
            function_classes[arg['fortran']['functionClass']] = True

    # Emit interface blocks so the compiler knows the GetPtr signatures.
    for fc in sorted(function_classes):
        code += (
            f'interface\n'
            f' function {fc}GetPtr(ptr_,classID)\n'
            f'  import c_int, c_ptr, {fc}Class\n'
            f'  class({fc}Class), pointer :: {fc}GetPtr\n'
            f'  type   (c_ptr), intent(in   ) :: ptr_\n'
            f'  integer(c_int), intent(in   ) :: classID\n'
            f' end function {fc}GetPtr\n'
            f'end interface\n'
        )
    return code


def fortran_reassignments(argument_list):
    """Generate Fortran reassignments."""
    return ''.join(arg.get('fortran', {}).get('reassignment', '')
                  for arg in argument_list)


def fortran_module_uses(argument_list):
    """Generate Fortran module use statements.

    Mirrors Perl fortranModuleUses(): accumulates {module: {symbol: 1}} dicts
    from every arg's fortran['modules'] field, then emits one 'use' line per
    module with a sorted 'only' list.
    """
    modules = {}
    for arg in argument_list:
        for mod_name, symbols in arg.get('fortran', {}).get('modules', {}).items():
            modules.setdefault(mod_name, {}).update(symbols)

    code = ''
    for mod_name in sorted(modules):
        syms = ', '.join(sorted(modules[mod_name]))
        code += f'  use :: {mod_name}, only : {syms}\n'
    return code


def fortran_call_code(argument_list, pre_arguments, post_arguments, continuation):
    """Generate Fortran call code, with optional-argument branching.

    Mirrors Perl fortranCallCode().  When N optional args are present, emits
    2^N if/else-if branches using present() intrinsic, plus an unconditional
    else fallback to suppress compiler warnings about unset function results.
    """
    def make_call(present_set):
        # present_set=None → all galacticus-present args included
        # otherwise only include optionals whose name is in present_set
        args = []
        for arg in argument_list:
            if not arg.get('galacticus', {}).get('isPresent'):
                continue
            if (arg.get('isOptional') and present_set is not None
                    and arg['name'] not in present_set):
                continue
            pass_name = arg.get('fortran', {}).get('passAs', arg['name'])
            args.append(f"{arg['name']}={pass_name}")
        return (f"{pre_arguments}{continuation} {(',' + continuation + chr(10) + '  ' + continuation + ' ').join(args)} {continuation}\n"
                f"{post_arguments}")

    optional_names = sorted(
        arg['name'] for arg in argument_list
        if arg.get('isOptional') and arg.get('galacticus', {}).get('isPresent')
    )

    if not optional_names:
        return make_call(None)

    code = ''
    first = True
    for subset in _powerset(optional_names):
        present_set = set(subset)
        conditions = [
            (f'present({n})' if n in present_set else f'.not.present({n})')
            for n in optional_names
        ]
        prefix = '' if first else 'else '
        code += f'{prefix}if ({" .and. ".join(conditions)}) then\n'
        code += make_call(present_set)
        first = False
    # Unconditional else avoids compiler warnings about unset function results.
    code += 'else\n'
    code += make_call(None)
    code += 'end if\n'
    return code


def iso_c_binding_import(argument_list, *extra_symbols):
    """Generate ISO_C_Binding import statement.

    Mirrors Perl isoCBindingImport(): collects the kind symbol from each
    argument's Fortran type (e.g. 'c_double' from 'real(c_double)') plus
    any extra symbols stored in fortran['isoCBindingSymbols'] (e.g.
    'c_f_pointer' added for pointer-dereference reassignments).
    """
    symbols = set(extra_symbols)
    for arg in argument_list:
        fort_type = arg.get('fortran', {}).get('type', '')
        m = re.search(r'\(([a-z_]+)\)', fort_type)
        if m:
            symbols.add(m.group(1))
        for sym in arg.get('fortran', {}).get('isoCBindingSymbols', []):
            symbols.add(sym)
    return f"  use, intrinsic :: ISO_C_Binding, only : {', '.join(sorted(symbols))}\n"


def python_arg_list(argument_list):
    """Generate Python argument list.

    Mirrors Perl pythonArgList(): if the first argument is not named 'self'
    (i.e. for constructors) prepend 'self' explicitly.  For methods the first
    argument already is 'self' (python.isPresent=1) so the loop adds it and
    no explicit prepend is needed.
    """
    first_name = argument_list[0]['name'] if argument_list else None
    args = [] if first_name == 'self' else ['self']

    first_optional = False
    for arg in argument_list:
        if not arg.get('python', {}).get('isPresent'):
            continue
        name = arg['name']
        if arg.get('isOptional') or first_optional:
            name += '=None'
            first_optional = True
        args.append(name)
    return args


def python_reassignments(argument_list):
    """Generate Python argument reassignments."""
    lines = []
    for arg in argument_list:
        if arg.get('python', {}).get('reassignment'):
            lines.append(arg['python']['reassignment'])
    return ''.join(lines)


def python_call_code(argument_list, call):
    """Generate Python call code, with optional-argument branching.

    Mirrors Perl pythonCallCode().  When N optional args are present, emits
    2^N if/elif branches.  Once the first optional arg is encountered, ALL
    subsequent args need explicit ctype wrapping (ctypes cannot infer types
    when some args may be None).  Absent optional args are passed as None.
    """
    def make_py_call(present_set, indent='        '):
        # present_set=None → all present; else set of presence-variable names
        args = []
        first_optional = False
        for arg in argument_list:
            if not arg.get('fortran', {}).get('isPresent'):
                continue
            is_opt = arg.get('isOptional', False)
            pv     = arg.get('python', {}).get('present', arg['name'])
            pa     = arg.get('python', {}).get('passAs', arg['name'])
            if is_opt:
                first_optional = True
            if first_optional:
                if is_opt and present_set is not None and pv not in present_set:
                    args.append('None')
                else:
                    ctype = arg.get('ctypes', {}).get('type', 'c_void_p')
                    args.append(f'{ctype}({pa})')
            else:
                args.append(pa)
        return f"{indent}{call}({','.join(args)})\n"

    # Collect the unique presence-variable names for optional fortran-present args.
    pv_seen = {}
    for arg in argument_list:
        if not arg.get('isOptional'):
            continue
        if not arg.get('fortran', {}).get('isPresent'):
            continue
        pv = arg.get('python', {}).get('present', arg['name'])
        pv_seen[pv] = True
    optional_pvs = sorted(pv_seen.keys())

    if not optional_pvs:
        return make_py_call(None, indent='    ')

    code = ''
    first = True
    for subset in _powerset(optional_pvs):
        present_set = set(subset)
        conditions = [
            (f'{pv} is not None' if pv in present_set else f'{pv} is None')
            for pv in optional_pvs
        ]
        prefix = '' if first else 'el'
        code += f"    {prefix}if {' and '.join(conditions)}:\n"
        code += make_py_call(present_set)
        first = False
    return code


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
