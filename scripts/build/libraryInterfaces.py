#!/usr/bin/env python3
# Generates C/Fortran/Python interfaces to Galacticus library classes.
# Andrew Benson (ported to Python with assistance from Claude 2026)

import sys
import os
import re
import keyword
import hashlib

import xml.etree.ElementTree as ET
from Galacticus.Build import SourceTree
from Galacticus.Build.FileChanges import update as file_changes_update
from Galacticus.Build.SourceTree.Parse import Declarations
from List.ExtraUtils import as_array
from Sort.Topo import sort as topo_sort
from XML.Utils import xml_to_dict

from LibraryInterfaces.Classification import (
    ENUM_RETURN_RX             as _ENUM_RETURN_RX,
    CLASS_RETURN_RX            as _CLASS_RETURN_RX,
    ARRAY_RETURN_RX            as _ARRAY_RETURN_RX,
    DYNAMIC_ARRAY_RETURN_RX    as _DYNAMIC_ARRAY_RETURN_RX,
    DYNAMIC_ARRAY_RETURN_2D_RX as _DYNAMIC_ARRAY_RETURN_2D_RX,
    TYPE_POINTER_RETURN_RX     as _TYPE_POINTER_RETURN_RX,
    normalize_method_return_type as _normalize_method_return_type,
    unsupported_arg            as _unsupported_arg,
    unsupported_output_array_method as _unsupported_output_array_method,
)
from LibraryInterfaces.Pipeline import (
    assign_c_types,
    assign_c_attributes,
    build_python_reassignments,
    build_fortran_reassignments,
    _SHARED_TYPE_MODULES,
    _CALLBACK_PROCEDURE_INTERFACES,
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
    python_safe_name,
)
_CLASS_HIERARCHY = {}


def main():
    """Main entry point — mirrors libraryInterfaces.pl."""

    # Initialize code and Python interface structures
    code = {'main': []}
    python = {'c_lib': [], 'units': {}}

    # Load XML configuration files
    build_path = os.environ.get('BUILDPATH', './work/build')
    exec_path = os.environ['GALACTICUS_EXEC_PATH']

    # Scan source/ for the derived-type hierarchy so the pipeline can
    # recognise `class(<intermediate>)` constructor args whose parent
    # chain reaches a registered functionClass.  Built once, read by
    # _unsupported_arg and passed into assign_c_types via the global.
    global _CLASS_HIERARCHY
    from LibraryInterfaces.Hierarchy import build_type_hierarchy
    _CLASS_HIERARCHY = build_type_hierarchy(
        os.path.join(exec_path, 'source'))

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


def _load_xml(path):
    """Load and parse an XML file into a nested dict structure, exiting with
    a descriptive message when the file is missing or cannot be parsed.
    """
    if not os.path.exists(path):
        sys.exit(f"libraryInterfaces.py: required XML file not found: {path}")
    try:
        root = ET.parse(path).getroot()
        return xml_to_dict(root)
    except ET.ParseError as exc:
        sys.exit(
            f"libraryInterfaces.py: failed to parse required XML file"
            f" '{path}': {exc}"
        )




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

def _find_enum_module(enum_type, func_class):
    """Locate the Fortran module that exports *enum_type* (an
    `enumerationXxxType` derived-type name).

    Order of resolution:
      0. explicit override in :data:`_SHARED_TYPE_MODULES` (wins outright;
         used for enums whose home module isn't imported by name in the
         functionClass's own ``use``-blocks),
      1. the functionClass file's own ``use``-block imports
         (collected during parsing as ``func_class['moduleUses']``),
      2. fall back to the functionClass's own module.
    Mirrors the scan ``build_fortran_reassignments`` already does for
    enumeration *arguments*.
    """
    explicit = _SHARED_TYPE_MODULES.get(enum_type)
    if explicit:
        return explicit
    for use_block in func_class.get('moduleUses', []):
        for mod_name, mod_data in use_block.items():
            # A parsed moduleUse node maps each module to a *list* of entries
            # (one per preprocessor condition set); tolerate a bare dict too.
            entries = mod_data if isinstance(mod_data, list) else [mod_data]
            for entry in entries:
                if (isinstance(entry, dict)
                        and enum_type in entry.get('only', {})):
                    return mod_name
    return func_class.get('module')


def _unsupported_constructor_arg(args, lib_function_classes,
                                 constructor_overrides=(),
                                 class_hierarchy=None):
    """If any constructor argument is unsupported, return ``(name, reason)``;
    otherwise ``None``.  See :func:`_unsupported_arg` for the predicate."""
    if class_hierarchy is None:
        class_hierarchy = _CLASS_HIERARCHY
    for arg in args:
        reason = _unsupported_arg(
            arg, lib_function_classes,
            constructor_overrides=constructor_overrides,
            class_hierarchy=class_hierarchy,
        )
        if reason:
            return arg['name'], reason
    return None


def _unsupported_method_arg(args, lib_function_classes, class_hierarchy=None):
    """If any method argument is unsupported, return ``(name, reason)``;
    otherwise ``None``.  Methods don't have constructor-style overrides for
    ``class(*)``, so the override list is empty and any ``class(*)`` arg is
    rejected.  `allow_pointer_writeback` is method-only: the method wrapper
    has a post-call hook for the pointer write-back protocol, the
    constructor wrapper does not."""
    if class_hierarchy is None:
        class_hierarchy = _CLASS_HIERARCHY
    for arg in args:
        reason = _unsupported_arg(arg, lib_function_classes,
                                  class_hierarchy=class_hierarchy,
                                  allow_pointer_writeback=True)
        if reason:
            return arg['name'], reason
    return None


def _process_implementations(func_class, directive_locations, state_storables,
                             code, python, lib_function_classes):
    """Process all implementations of a function class."""
    class_name = func_class['name']
    impls = directive_locations.get(class_name, {}).get('file', [])

    extensions = {}
    module_uses_impls = {}
    class_id = 0
    impls_list = []

    # Each implementation file is parsed once and its tree walked twice —
    # first for extensions/module-uses, then for constructor discovery.
    # (Previously two sequential passes each re-parsed every file; parsing
    # dominates this script's cost.)
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

        # Constructor discovery (second walk of the same tree).
        impl_name = None
        is_abstract = False
        name_constructor = None
        ambiguous_internals = None
        args_constructor = []
        impl_conf = None

        for node in SourceTree.walk_tree(tree):
            if node['type'] == class_name:
                impl_name = node.get('directive', {}).get('name')
                is_abstract = (node.get('directive', {}).get('abstract', 'no') == 'yes')
                # Fetch the libraryClasses.xml override up front so the
                # interface-block walk below can consult its
                # `<constructor internal="..."/>` hint to disambiguate
                # multiple Internal-suffixed module procedures.
                impl_conf = func_class.get(impl_name)

            elif (impl_name
                  and node['type'] == 'interface'
                  and node.get('name', '').lower() == impl_name.lower()):
                # Collect every Internal-marked module procedure across this
                # interface's children.  The Galacticus convention names
                # constructors `<short>Constructor<Variant>` where Variant
                # is `Parameters` (XML-driven) or `Internal[Suffix]`; we
                # want only the Internal flavour.  Some classes (the merger
                # tree walkers, for example) use the shorter `<short>Internal`
                # form without a `Constructor` infix.  Accepting either
                # `endswith('internal')` or `'constructorinternal' in name`
                # covers both conventions plus the rarer
                # `ConstructorInternalType` / `ConstructorInternalDefined`
                # disambiguation suffixes; `<short>ConstructorParameters`
                # satisfies neither rule and is correctly rejected.
                candidates = []
                child = node.get('firstChild')
                while child:
                    if child['type'] == 'moduleProcedure':
                        candidates.extend(
                            n for n in child.get('names', [])
                            if (n.lower().endswith('internal')
                                or 'constructorinternal' in n.lower())
                        )
                    child = child.get('sibling')
                if len(candidates) == 1:
                    name_constructor = candidates[0]
                elif len(candidates) > 1:
                    # Multiple Internal-suffixed constructors (e.g.
                    # darkMatterProfileConcentrationDuttonMaccio2014's
                    # InternalType vs InternalDefined).  libraryClasses.xml
                    # can disambiguate via
                    # `<constructor internal="<chosen-name>"/>`; if a
                    # valid hint is present, use it.  Otherwise stash
                    # for the warning-and-skip below; leave
                    # name_constructor unset so the impl is dropped
                    # from impls_list.
                    chosen = None
                    if isinstance(impl_conf, dict):
                        chosen = (impl_conf.get('constructor', {}) or {}) \
                                 .get('internal')
                    if chosen and chosen in candidates:
                        name_constructor = chosen
                    else:
                        ambiguous_internals = candidates

            elif (name_constructor
                  and node['type'] == 'function'
                  and node.get('name', '').lower() == name_constructor.lower()):
                # Extract argument names from the function opener.
                opener = node.get('opener', '')
                m = re.search(
                    r'function\s+' + re.escape(name_constructor) + r'\s*\(([^)]+)\)',
                    opener, re.IGNORECASE)
                if m:
                    # Strip Fortran line-continuation characters (`&`) and any
                    # whitespace (including embedded newlines from multi-line
                    # function openers) from each captured argument name —
                    # otherwise tokens like "&\n   &  delta_0" leak into
                    # downstream emitters (e.g. as text in <referenceConstruct>)
                    # and break XML parsing because the literal `&` isn't escaped.
                    # Drop entries that strip to empty: an argument-less opener
                    # written with whitespace between the parentheses (`( )`, as
                    # opposed to `()`) captures a lone space here, which would
                    # otherwise become a phantom empty-named argument and emit a
                    # broken `integer(c_int) ::` declaration and `& = &` call.
                    args_constructor = [
                        {'name': name}
                        for name in (re.sub(r'[\s&]+', '', a)
                                     for a in m.group(1).split(','))
                        if name
                    ]
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

        # impl_conf was already fetched up front so the interface-block
        # walk could consult the optional `<constructor internal=…/>`
        # disambiguation hint.
        is_excluded = (isinstance(impl_conf, dict)
                       and impl_conf.get('exclude') == 'yes')
        if not is_abstract and not is_excluded:
            if ambiguous_internals:
                # Drop the impl: we can't pick one of N Internal constructors
                # without a hint, and emitting the previous fall-back
                # (constructor call with no args) crashes gfortran with
                # "No initializer for component 'X' given …" whenever the
                # type has any required components.
                sys.stderr.write(
                    f"libraryInterfaces.py: caution: implementation"
                    f" '{impl_name}' of class '{class_name}' has multiple"
                    f" Internal-suffixed constructors"
                    f" ({', '.join(ambiguous_internals)}) — ambiguous,"
                    f" skipping implementation\n"
                )
            else:
                constructor_overrides = ()
                if isinstance(impl_conf, dict):
                    constructor_overrides = as_array(
                        impl_conf.get('constructor', {}).get('argument', []))
                unsupported = _unsupported_constructor_arg(
                    args_constructor, lib_function_classes, constructor_overrides)
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


def _shared_bucket(code, class_name):
    """Return the per-class 'shared' code list, creating the bucket if needed.

    ``code[class_name]`` is a dict ``{'shared': […], 'per_impl': {…}}`` —
    'shared' holds the per-class pieces (GetPtr, GetIdAndPtr, methods,
    destructor) that all live together in ``<class>.F90``, and
    'per_impl[impl_name]' holds the constructor wrapper for one impl,
    written to its own ``<class>__<impl>.F90`` so that gfortran doesn't
    have to compile every impl's bind(c) wrapper in a single
    compilation unit.
    """
    return code.setdefault(class_name, {'shared': [], 'per_impl': {}})['shared']


def _impl_bucket(code, class_name, impl_name):
    """Return the per-impl code list, creating the bucket if needed.

    See :func:`_shared_bucket` for the structure.
    """
    return (code.setdefault(class_name, {'shared': [], 'per_impl': {}})
                ['per_impl'].setdefault(impl_name, []))


def interfaces_pointer_get(code, func_class):
    """Generate pointer getters for a function class.

    Emits two Fortran helpers:

    * ``<class>GetPtr(c_ptr, classID) -> class(<class>Class), pointer`` —
      the inverse, used by every wrapper that takes a ``class(<class>Class)``
      argument from the Python side: c_ptr + classID → typed pointer.

    * ``<class>GetIdAndPtr(class_obj, classID) -> c_ptr`` — the *forward*
      direction, used by every wrapper whose method *returns* a
      ``class(<class>Class)``.  Walks the class's concrete impls with a
      `select type` block, sets ``classID`` to the matching impl's id,
      and returns ``c_loc(obj)`` for the concrete typed pointer (which
      ``c_loc`` requires — it can't take a polymorphic target directly).
    """
    class_name = func_class['name']
    impls = func_class.get('implementations', [])

    symbols = [class_name + 'Class'] + [impl['name'] for impl in impls]

    _shared_bucket(code, class_name).append(f'''function {class_name}GetPtr({class_name}_,classID)
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

    # Forward helper used by class(...) return-type method wrappers.
    select_branches = chr(10).join(
        f'  type is ({impl["name"]})\n'
        f'     ptr     = c_loc(obj)\n'
        f'     classID = {impl["classID"]}'
        for impl in impls
    )
    _shared_bucket(code, class_name).append(f'''function {class_name}GetIdAndPtr(obj,classID) result(ptr)
  use, intrinsic :: ISO_C_Binding, only : c_ptr, c_int, c_loc, c_null_ptr
  use :: {func_class.get('module', 'Unknown')}, only : {', '.join(symbols)}
  implicit none
  type(c_ptr) :: ptr
  class({class_name}Class), pointer, intent(in) :: obj
  integer(c_int), intent(out) :: classID

  if (.not.associated(obj)) then
     ptr     = c_null_ptr
     classID = -1
     return
  end if
  select type (obj)
{select_branches}
  class default
     ptr     = c_null_ptr
     classID = -1
  end select
  return
end function {class_name}GetIdAndPtr
''')


def interfaces_python_classes(python, func_class):
    """Generate Python class hierarchy for a function class."""
    class_name = func_class['name']

    # Parent class.  _from_classID() is the entry point used by methods
    # that return class(...)-typed objects: it walks the direct subclasses,
    # picks the one whose _classIDStatic matches the runtime classID, and
    # constructs an instance via __new__ (skipping __init__ so we don't
    # build a *new* Galacticus object — we just wrap the existing pointer).
    # The returned wrapper is marked _owned=False so the destructor
    # doesn't try to free Galacticus's object.
    parent_code = f'''class {class_name}:

    # Constructor
    def __init__(self):
        # Assign class ID to negative (not a concrete class)
        self._classID = -1

    @classmethod
    def _from_classID(cls, classID, ptr):
        for subcls in cls.__subclasses__():
            if getattr(subcls, '_classIDStatic', None) == classID:
                obj = subcls.__new__(subcls)
                obj._glcObj  = ptr
                obj._classID = classID
                obj._owned   = False
                return obj
        raise ValueError(f"Unknown classID {{classID}} for {{cls.__name__}}")
'''
    python['units'][class_name] = {
        'content': parent_code,
        'indent': 0,
        'dependencies': ['init'],
    }

    # Child classes.  _classIDStatic is the class-level companion to the
    # per-instance _classID set by the constructor, used by the parent
    # class's _from_classID to pick the right subclass without having
    # to instantiate one.
    for impl in func_class.get('implementations', []):
        child_code = (
            f"class {impl['name']}({class_name}):\n"
            f"    _classIDStatic = {impl['classID']}"
        )
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
        # Pull the impl's libraryClasses.xml overrides (if any) into
        # assign_c_types so it can recognise `value="null"` directives
        # and drop the matching args from both wrappers before the
        # Python/Fortran emitters see them.
        impl_conf = func_class.get(impl['name'])
        constructor_overrides = ()
        if isinstance(impl_conf, dict):
            constructor_overrides = as_array(
                impl_conf.get('constructor', {}).get('argument', []))
        arg_list = assign_c_types(arg_list, lib_function_classes,
                                  class_hierarchy=_CLASS_HIERARCHY,
                                  constructor_overrides=constructor_overrides)
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
        _impl_bucket(code, class_name, impl['name']).append(fort_constructor)

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
        method_type = _normalize_method_return_type(
            method_spec.get('type', 'void'))

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
        result_pre_call_code      = ""    # emitted before the inner call
        result_post_call_code     = ""    # emitted after the inner call
        result_call_target        = method_name_c
        result_assign_op          = '='   # `=>` for pointer-returning methods (class(...))
        result_python_decode      = False
        result_extra_fort_args    = []    # extra args appended to bind(c) signature
        result_extra_clib_argtypes = []   # matching ctypes wrappers
        result_python_class_wrap  = None  # (parent_class, out_id_var) or None
        result_python_array_wrap  = None  # (size, elem_ctype, elem_dtype) or None
        result_python_dyn_array_wrap = None  # (elem_ctype, elem_dtype) for save-buffer 1D returns
        result_python_dyn_array_2d_wrap = None  # (elem_ctype, elem_dtype) for save-buffer 2D returns
        result_is_subroutine      = False # force `subroutine` even when method_type != "void"
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
            #
            # Local-variable names are short, generic, and scoped to the
            # bind(c) function — the Galacticus convention is a `glc` prefix
            # plus a trailing underscore so they don't collide with user-named
            # method arguments.  This also avoids the Fortran 63-character
            # identifier limit, which long method names would otherwise trip
            # if we used `<method>_result_` etc.
            method_type_c                        = "type(c_ptr)"
            clib_res_type                        = "c_char_p"
            isoImports['c_ptr']       = 1
            isoImports['c_loc']       = 1
            isoImports['c_char']      = 1
            isoImports['c_null_char'] = 1
            result_extra_module_uses = (
                f'  use :: ISO_Varying_String, only : varying_string, char\n'
            )
            result_call_target = 'glcResult_'
            result_extra_declarations = (
                f'  type     (varying_string)                                          :: glcResult_\n'
                f'  character(kind=c_char   ), dimension(:), allocatable, save, target :: glcBuffer_\n'
                f'  character(len=:         ), allocatable                             :: glcChars_\n'
                f'  integer                                                            :: glcI_\n'
            )
            result_post_call_code = (
                f'  glcChars_ = char(glcResult_)\n'
                f'  if (allocated(glcBuffer_)) deallocate(glcBuffer_)\n'
                f'  allocate(glcBuffer_(len(glcChars_)+1))\n'
                f'  do glcI_ = 1, len(glcChars_)\n'
                f'     glcBuffer_(glcI_) = glcChars_(glcI_:glcI_)\n'
                f'  end do\n'
                f'  glcBuffer_(len(glcChars_)+1) = c_null_char\n'
                f'  {method_name_c} = c_loc(glcBuffer_)\n'
            )
            result_python_decode = True
        elif _ENUM_RETURN_RX.match(method_type):
            # Returned type(enumerationXxxType): the inner method gives us a
            # derived type whose %ID component holds the c_int value we want
            # to surface to Python.  Same scaffolding shape as varying_string
            # — call into a temporary, then pull the c_int out — minus the
            # buffer/decode dance.
            enum_type = _ENUM_RETURN_RX.match(method_type).group(1)
            enum_module = _find_enum_module(enum_type, func_class)
            method_type_c                        = "integer(c_int)"
            clib_res_type                        = "c_int"
            isoImports['c_int'] = 1
            if enum_module:
                result_extra_module_uses = (
                    f'  use :: {enum_module}, only : {enum_type}\n'
                )
            result_call_target = 'glcResult_'
            result_extra_declarations = (
                f'  type({enum_type}) :: glcResult_\n'
            )
            result_post_call_code = (
                f'  {method_name_c} = glcResult_%ID\n'
            )
        elif _CLASS_RETURN_RX.match(method_type):
            # Returned class(FooClass): the inner method gives us a
            # polymorphic pointer.  We can't c_loc() it directly (Fortran
            # forbids that on a polymorphic target), so we hand it to
            # FooGetIdAndPtr — generated by interfaces_pointer_get for
            # every registered class — which dispatches via select type
            # to extract a typed c_loc + the matching classID.  Python
            # gets back (c_void_p, c_int via byref) and dispatches into
            # the right Foo subclass via Foo._from_classID.
            return_class_type = _CLASS_RETURN_RX.match(method_type).group(1)
            return_stem = (return_class_type[:-5]
                           if return_class_type.endswith('Class')
                           else return_class_type)
            if return_stem not in (lib_function_classes or {}):
                sys.stderr.write(
                    f"libraryInterfaces.py: caution: method '{method_name}' in"
                    f" class '{class_name}' returns class({return_class_type}) —"
                    f" '{return_stem}' is not a registered functionClass in"
                    f" libraryClasses.xml — skipping method\n"
                )
                methods_to_delete.append(method_name)
                continue
            return_module = (lib_function_classes or {}) \
                            .get(return_stem, {}).get('module')
            # Short, scoped Fortran identifiers — see the varying_string
            # branch above for the same rationale (Fortran 63-char limit
            # plus collision-avoidance with user-named method arguments).
            out_classID_name = 'glcCidOut_'
            method_type_c   = "type(c_ptr)"
            clib_res_type   = "c_void_p"
            isoImports['c_ptr'] = 1
            isoImports['c_int'] = 1
            if return_module:
                result_extra_module_uses = (
                    f'  use :: {return_module}, only : {return_class_type}\n'
                )
            result_call_target = 'glcResult_'
            result_assign_op   = '=>'   # method returns a polymorphic pointer
            result_extra_declarations = (
                f'  class({return_class_type}), pointer :: glcResult_\n'
                f'  integer(c_int), intent(out) :: {out_classID_name}\n'
                f'  interface\n'
                f'    function {return_stem}GetIdAndPtr(obj,classID) result(ptr)\n'
                f'      import :: c_ptr, c_int, {return_class_type}\n'
                f'      type(c_ptr) :: ptr\n'
                f'      class({return_class_type}), pointer, intent(in) :: obj\n'
                f'      integer(c_int), intent(out) :: classID\n'
                f'    end function {return_stem}GetIdAndPtr\n'
                f'  end interface\n'
            )
            result_post_call_code = (
                f'  {method_name_c} = {return_stem}GetIdAndPtr('
                f'glcResult_,{out_classID_name})\n'
            )
            result_extra_fort_args     = [out_classID_name]
            result_extra_clib_argtypes = ['POINTER(c_int)']
            result_python_class_wrap   = (return_stem, out_classID_name)
        elif _ARRAY_RETURN_RX.match(method_type):
            # Returned 1D fixed-size numeric array (e.g.
            # `double precision, dimension(3)`).  bind(c) functions can't
            # return arrays directly, so we lower to a subroutine with an
            # extra `intent(out), dimension(N)` arg and have the inner
            # method's result assigned into it.  Python pre-allocates a
            # numpy array of the right shape/dtype, passes it via
            # data_as(POINTER(...)), and returns it.
            m = _ARRAY_RETURN_RX.match(method_type)
            arr_intrinsic_raw = m.group(1).lower()
            arr_intrinsic     = re.sub(r'\s+', ' ', arr_intrinsic_raw)
            arr_size = int(m.group(2))
            if arr_intrinsic == 'double precision':
                elem_ctype, elem_fort, elem_dtype = ('c_double',
                                                    'real(c_double)',
                                                    'float64')
            else:                                # integer (default kind)
                elem_ctype, elem_fort, elem_dtype = ('c_int',
                                                    'integer(c_int)',
                                                    'int32')
            isoImports[elem_ctype] = 1
            method_type_c              = ''      # subroutine, no return type
            clib_res_type              = None
            result_call_target         = 'glcResult_'
            result_extra_declarations  = (
                f'  {elem_fort}, dimension({arr_size}), intent(out) ::'
                f' glcResult_\n'
            )
            result_extra_fort_args     = ['glcResult_']
            result_extra_clib_argtypes = [f'POINTER({elem_ctype})']
            result_python_array_wrap   = (arr_size, elem_ctype, elem_dtype)
            result_is_subroutine       = True
        elif _DYNAMIC_ARRAY_RETURN_RX.match(method_type):
            # Returned 1D allocatable or dynamic-size numeric array.
            # Covers both source-level shapes:
            #
            #   <type>double precision, allocatable, dimension(:)</type>
            #   <type>double precision, dimension(self%parameterCount)</type>
            #
            # bind(c) functions can't return arrays, so we lower to a
            # subroutine with two intent(out) companions — a c_ptr that
            # the wrapper sets to c_loc(glcResult_) and a c_size_t that
            # the wrapper sets to size(glcResult_).  glcResult_ is a
            # function-local `save, target` allocatable; Fortran 2003
            # auto-realloc on assignment resizes it to the inner method's
            # result shape on every call (both for genuinely allocatable
            # results and for fixed-shape-but-runtime-extent results).
            #
            # On the Python side we wrap the (ptr, size) pair into a
            # numpy view and *copy* before returning — the save buffer
            # is overwritten on the next call to the same method, so
            # holding a reference to it past the next call would alias
            # to fresh data.  Mirrors the `varying_string` return
            # convention (same buffer-lifetime contract; different
            # element type).
            m = _DYNAMIC_ARRAY_RETURN_RX.match(method_type)
            arr_intrinsic_raw = m.group(1).lower()
            arr_intrinsic     = re.sub(r'\s+', ' ', arr_intrinsic_raw)
            if arr_intrinsic == 'double precision':
                elem_ctype, elem_fort, elem_dtype = ('c_double',
                                                    'real(c_double)',
                                                    'float64')
            else:
                elem_ctype, elem_fort, elem_dtype = ('c_int',
                                                    'integer(c_int)',
                                                    'int32')
            isoImports[elem_ctype]    = 1
            isoImports['c_ptr']       = 1
            isoImports['c_size_t']    = 1
            isoImports['c_loc']       = 1
            method_type_c              = ''      # subroutine, no return type
            clib_res_type              = None
            result_call_target         = 'glcResult_'
            # The save-target buffer is a function-local allocatable; on
            # entry we deallocate any prior contents (the previous call's
            # data is no longer needed — Python copied it) and then let
            # the assignment auto-allocate to the new shape.
            result_extra_declarations  = (
                f'  {elem_fort}, dimension(:), allocatable, save, target'
                f' :: glcResult_\n'
                f'  type(c_ptr),       intent(out) :: glcDataPtr_\n'
                f'  integer(c_size_t), intent(out) :: glcSize_\n'
            )
            result_post_call_code      = (
                f'  glcSize_    = size(glcResult_, kind=c_size_t)\n'
                f'  glcDataPtr_ = c_loc(glcResult_)\n'
            )
            # Buffer hygiene: a previous call may have left glcResult_
            # allocated to a different shape; deallocate so the next
            # assignment auto-allocates cleanly (Fortran 2003 auto-realloc
            # would also handle a same-shape reuse, but we want to free
            # any unused-memory excess too).
            result_pre_call_code       = (
                f'  if (allocated(glcResult_)) deallocate(glcResult_)\n'
            )
            result_extra_fort_args     = ['glcDataPtr_', 'glcSize_']
            result_extra_clib_argtypes = ['POINTER(c_void_p)',
                                          'POINTER(c_size_t)']
            result_python_dyn_array_wrap = (elem_ctype, elem_dtype)
            result_is_subroutine       = True
        elif _DYNAMIC_ARRAY_RETURN_2D_RX.match(method_type):
            # 2D variant of the allocatable / dynamic-size array return
            # — same save-target buffer trick as the 1D case above, with
            # one extra `c_size_t` size companion so Python can reshape
            # the flat byte buffer to the right (size1, size2).  Fortran
            # stores 2D arrays in column-major order; the Python wrapper
            # uses `reshape((size1, size2), order='F')` so the indexing
            # convention matches the inner method.
            m = _DYNAMIC_ARRAY_RETURN_2D_RX.match(method_type)
            arr_intrinsic_raw = m.group(1).lower()
            arr_intrinsic     = re.sub(r'\s+', ' ', arr_intrinsic_raw)
            if arr_intrinsic == 'double precision':
                elem_ctype, elem_fort, elem_dtype = ('c_double',
                                                    'real(c_double)',
                                                    'float64')
            else:
                elem_ctype, elem_fort, elem_dtype = ('c_int',
                                                    'integer(c_int)',
                                                    'int32')
            isoImports[elem_ctype]    = 1
            isoImports['c_ptr']       = 1
            isoImports['c_size_t']    = 1
            isoImports['c_loc']       = 1
            method_type_c              = ''      # subroutine, no return type
            clib_res_type              = None
            result_call_target         = 'glcResult_'
            result_extra_declarations  = (
                f'  {elem_fort}, dimension(:,:), allocatable, save, target'
                f' :: glcResult_\n'
                f'  type(c_ptr),       intent(out) :: glcDataPtr_\n'
                f'  integer(c_size_t), intent(out) :: glcSize1_\n'
                f'  integer(c_size_t), intent(out) :: glcSize2_\n'
            )
            result_post_call_code      = (
                f'  glcSize1_   = size(glcResult_, dim=1, kind=c_size_t)\n'
                f'  glcSize2_   = size(glcResult_, dim=2, kind=c_size_t)\n'
                f'  glcDataPtr_ = c_loc(glcResult_)\n'
            )
            result_pre_call_code       = (
                f'  if (allocated(glcResult_)) deallocate(glcResult_)\n'
            )
            result_extra_fort_args     = ['glcDataPtr_', 'glcSize1_', 'glcSize2_']
            result_extra_clib_argtypes = ['POINTER(c_void_p)',
                                          'POINTER(c_size_t)',
                                          'POINTER(c_size_t)']
            result_python_dyn_array_2d_wrap = (elem_ctype, elem_dtype)
            result_is_subroutine       = True
        elif (_TYPE_POINTER_RETURN_RX.match(method_type)
              and _TYPE_POINTER_RETURN_RX.match(method_type).group(1)
                  in _SHARED_TYPE_MODULES):
            # `type(<X>), pointer` return for a shared type: returned to
            # Python as an OPAQUE HANDLE. The wrapper captures the result in
            # a local pointer and returns c_loc of its target (c_null_ptr
            # when disassociated); Python receives the address as an integer
            # handle that can be passed back into any `type(<X>)` argument
            # (e.g. a mergerTree handle into the tree-walker constructors).
            # Ownership stays with Fortran: the handle does not free the
            # object.
            handle_type = _TYPE_POINTER_RETURN_RX.match(method_type).group(1)
            method_type_c              = 'type(c_ptr)'
            clib_res_type              = 'c_void_p'
            result_call_target         = 'glcResult_'
            result_assign_op           = '=>'
            result_extra_declarations  = (
                f'  type({handle_type}), pointer :: glcResult_\n')
            result_extra_module_uses   = (
                f'  use :: {_SHARED_TYPE_MODULES[handle_type]}, only :'
                f' {handle_type}\n')
            result_post_call_code      = (
                f'  if (associated(glcResult_)) then\n'
                f'    {method_name_c} = c_loc(glcResult_)\n'
                f'  else\n'
                f'    {method_name_c} = c_null_ptr\n'
                f'  end if\n'
            )
            isoImports['c_ptr']      = 1
            isoImports['c_loc']      = 1
            isoImports['c_null_ptr'] = 1
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

        # Skip the method if any argument has a type the pipeline can't
        # translate (complex/double complex, class(non-registered), class(*)).
        # See :func:`_unsupported_arg` for the predicate.  Self (arg_list[0])
        # is class(<the current class>Class) — registered by definition — so
        # the slice [1:] is just to avoid noise in the iteration.
        unsupported_arg = _unsupported_method_arg(
            arg_list[1:], lib_function_classes or {})
        if unsupported_arg:
            sys.stderr.write(
                f"libraryInterfaces.py: caution: method '{method_name}' in"
                f" class '{class_name}' has argument '{unsupported_arg[0]}'"
                f" of unsupported kind ({unsupported_arg[1]}) — skipping"
                f" method\n"
            )
            methods_to_delete.append(method_name)
            continue

        # Whole-method gate for output-array args (`intent(out),
        # allocatable, dimension(:)`).  A per-arg check accepts each output
        # array, but the current codegen only handles them on void-returning
        # methods with no optionals whose other args are supported inputs;
        # anything else (scalar intent(out) companion, inout allocatable,
        # non-void return) would have an output the wrapper silently drops.
        # Shared with the audit so their verdicts can't drift.
        output_array_block = _unsupported_output_array_method(
            [{'name': a['name'], 'intrinsic': a['intrinsic'],
              'type': a['type'], 'attributes': a['attributes']}
             for a in arg_list[1:]],
            method_type)
        if output_array_block:
            sys.stderr.write(
                f"libraryInterfaces.py: caution: method '{method_name}' in"
                f" class '{class_name}': {output_array_block} — skipping"
                f" method\n"
            )
            methods_to_delete.append(method_name)
            continue

        # Process arguments
        arg_list = assign_c_types(arg_list, lib_function_classes or {},
                                  class_hierarchy=_CLASS_HIERARCHY)
        arg_list = assign_c_attributes(arg_list)
        arg_list = build_python_reassignments(arg_list)
        arg_list = build_fortran_reassignments(arg_list, func_class, None,
                                              extensions, module_uses_impls,
                                              lib_function_classes)

        # Output-array args (`intent(out), allocatable, dimension(:)`): each
        # is a wrapper-local `save, target` allocatable (declared + pre-call
        # deallocated by assign_c_types) that the inner call fills.  Append a
        # `(c_ptr, c_size_t)` intent(out) companion pair per output array to
        # the bind(c) signature and set them, after the call, to
        # `c_loc(buffer)` and `size(buffer)`.  The Python wrapper (bespoke
        # branch below) wraps each pair into a numpy array.  Companions are
        # appended at the end of the signature; the Python call passes them
        # in this same order, so ctypes/bind(c) stay aligned.  Mirrors the
        # single-slot dynamic-array *return* path (the _DYNAMIC_ARRAY_RETURN
        # branch), generalised to N per-argument outputs.
        # Pointer write-back args: after the inner call, write the local
        # Fortran pointer's (possibly repointed) target address back
        # through the by-reference c_ptr handle — c_null_ptr when the
        # method left it disassociated (end-of-iteration for a tree
        # walker).  The pre-call guarded c_f_pointer lives in the arg's
        # fort_reassignment (see assign_c_types).
        for a in (x for x in arg_list if x.is_pointer_writeback):
            local = f'{a.name}_'
            result_post_call_code += (
                f'  if (associated({local})) then\n'
                f'    {a.name} = c_loc({local})\n'
                f'  else\n'
                f'    {a.name} = c_null_ptr\n'
                f'  end if\n'
            )

        # Per-argument copy-back code (e.g. kind-narrowing an
        # intent(out) logical's default-kind local back into the c_bool
        # bind(c) dummy) — set by build_fortran_reassignments.
        for a in arg_list:
            if a.fort_postcall:
                result_post_call_code += f'  {a.fort_postcall}'

        output_arrays = [a for a in arg_list if a.is_output_array]
        sized_outputs = [a for a in arg_list if a.is_output_sized]
        for a in output_arrays:
            local     = f'glcOut_{a.name}_'
            ptr_name  = f'{a.name}DataPtr_'
            result_extra_declarations += (
                f'  type(c_ptr),       intent(out) :: {ptr_name}\n'
            )
            # Logical output arrays: the inner method filled the
            # default-kind `glcOutInner_` local (see assign_c_types);
            # kind-narrow it into the c_bool export buffer before taking
            # its address (auto-realloc sizes the LHS — both buffers were
            # deallocated pre-call).
            if a.intrinsic == 'logical':
                inner = f'glcOutInner_{a.name}_'
                result_post_call_code += (
                    f'  if (allocated({inner})) then\n'
                    f'    {local} = logical({inner}, c_bool)\n'
                    f'  end if\n'
                )
            # Guard c_loc: it is non-conforming on a zero-size or unallocated
            # array (F2008 requires a contiguous, nonzero-size target), and
            # an intent(out) allocatable may be left unallocated (or
            # allocated to length 0 — e.g. transferFunction's default
            # `allocate(wavenumbers(0))`) by the inner method.  In those
            # cases report size 0 and a null pointer; the Python wrapper
            # returns an empty array without dereferencing it.
            if a.array_rank == 2:
                # 2D: one size companion per axis; the Python wrapper
                # reshapes the flat buffer column-major (order='F').
                s1 = f'{a.name}Size1_'
                s2 = f'{a.name}Size2_'
                result_extra_declarations += (
                    f'  integer(c_size_t), intent(out) :: {s1}\n'
                    f'  integer(c_size_t), intent(out) :: {s2}\n'
                )
                result_post_call_code += (
                    f'  if (allocated({local})) then\n'
                    f'    {s1} = size({local}, dim=1, kind=c_size_t)\n'
                    f'    {s2} = size({local}, dim=2, kind=c_size_t)\n'
                    f'    if ({s1}*{s2} > 0_c_size_t) then\n'
                    f'      {ptr_name} = c_loc({local})\n'
                    f'    else\n'
                    f'      {ptr_name} = c_null_ptr\n'
                    f'    end if\n'
                    f'  else\n'
                    f'    {s1} = 0_c_size_t\n'
                    f'    {s2} = 0_c_size_t\n'
                    f'    {ptr_name} = c_null_ptr\n'
                    f'  end if\n'
                )
                result_extra_fort_args     += [ptr_name, s1, s2]
                result_extra_clib_argtypes += ['POINTER(c_void_p)',
                                               'POINTER(c_size_t)',
                                               'POINTER(c_size_t)']
            else:
                size_name = f'{a.name}Size_'
                result_extra_declarations += (
                    f'  integer(c_size_t), intent(out) :: {size_name}\n'
                )
                result_post_call_code += (
                    f'  if (allocated({local})) then\n'
                    f'    {size_name} = size({local}, kind=c_size_t)\n'
                    f'    if ({size_name} > 0_c_size_t) then\n'
                    f'      {ptr_name} = c_loc({local})\n'
                    f'    else\n'
                    f'      {ptr_name} = c_null_ptr\n'
                    f'    end if\n'
                    f'  else\n'
                    f'    {size_name} = 0_c_size_t\n'
                    f'    {ptr_name} = c_null_ptr\n'
                    f'  end if\n'
                )
                result_extra_fort_args     += [ptr_name, size_name]
                result_extra_clib_argtypes += ['POINTER(c_void_p)',
                                               'POINTER(c_size_t)']
            isoImports['c_ptr']            = 1
            isoImports['c_size_t']         = 1
            isoImports['c_loc']            = 1
            isoImports['c_null_ptr']       = 1
            isoImports[a.output_elem_ctype] = 1   # kind for the local buffer

        # Inbound callback args (`procedure(<iface>)` with <iface> registered
        # in _CALLBACK_PROCEDURE_INTERFACES; marked is_callback by
        # assign_c_types).  Fortran has no closures, so per callback-taking
        # method we emit a small module holding a `type(c_funptr)` slot and a
        # shim function matching the Galacticus-side interface; the bind(c)
        # wrapper stores the incoming funptr in the slot (fort_reassignment)
        # and passes the shim to the inner method (fort_pass_as).  The slot
        # is process-global: the callback is valid for the duration of the
        # wrapped call only, matching the save-buffer conventions.  On the
        # Python side the user's callable is adapted (py_adapter) and wrapped
        # in the registry's CFUNCTYPE; the local reference keeps it alive for
        # the duration of the call.
        callback_args = [a for a in arg_list if a.is_callback]
        if callback_args:
            module_name = f'glcCB_{class_name}_{method_name}'
            if len(module_name) > 63:
                # Fortran caps identifiers at 63 chars; hash-truncate the
                # tail deterministically (same scheme the interface names
                # would need — nothing hits this today).
                digest      = hashlib.md5(module_name.encode()).hexdigest()[:8]
                module_name = module_name[:55] + digest
            iface_blocks = ''
            var_decls    = ''
            shim_blocks  = ''
            shim_uses    = {}
            for a in callback_args:
                spec  = _CALLBACK_PROCEDURE_INTERFACES[a.type_spec]
                var   = f'glcCBPtr_{a.name}_'
                shim  = f'glcCBShim_{a.name}_'
                iface = f'glcCBIface_{a.name}_'
                iface_blocks += spec['c_iface'].format(iface=iface)
                var_decls    += f'  type(c_funptr) :: {var} = c_null_funptr\n'
                shim_blocks  += spec['shim'].format(shim=shim, iface=iface,
                                                    var=var)
                for mod, symbols in spec['shim_uses'].items():
                    shim_uses.setdefault(mod, set()).update(symbols)
                if 'pointer' in a.attributes:
                    # `procedure(...), pointer` dummy: the actual must
                    # itself be a procedure pointer, so the module also
                    # carries a slot declared with the registry's
                    # Galacticus-side abstract interface (a module
                    # procedure from the same contains part can't serve as
                    # the interface-name in the specification part).  The
                    # wrapper re-aims it at the shim before every call, so
                    # a callee that repoints it does no lasting harm.
                    giface = f'glcCBGIface_{a.name}_'
                    pp     = f'glcCBPP_{a.name}_'
                    iface_blocks += spec['g_iface'].format(giface=giface)
                    var_decls    += (f'  procedure({giface}), pointer ::'
                                     f' {pp} => null()\n')
                    a.fort_reassignment = (f'  {var} = {a.name}\n'
                                           f'  {pp} => {shim}\n')
                    a.fort_pass_as      = pp
                    a.fort_modules      = {module_name:
                                           {var: 1, shim: 1, pp: 1}}
                else:
                    # Plain procedure dummy: pass the shim directly (module
                    # procedures have explicit interfaces).
                    a.fort_reassignment = f'  {var} = {a.name}\n'
                    a.fort_pass_as      = shim
                    a.fort_modules      = {module_name: {var: 1, shim: 1}}
                # Python-side wiring: adapt the user callable, wrap in the
                # CFUNCTYPE, pass as a void pointer.
                pv        = python_safe_name(a.name)
                adapted   = spec['py_adapter'].format(fn=pv)
                a.py_reassignment = (
                    f'    _glcCBType_{pv}_ = {spec["cfunctype"]}\n'
                    f'    _glcCB_{pv}_ = _glcCBType_{pv}_({adapted})\n'
                )
                a.py_pass_as = f'cast(_glcCB_{pv}_, c_void_p)'
            use_lines = ''.join(
                f'  use :: {mod}, only : {", ".join(sorted(syms))}\n'
                for mod, syms in sorted(shim_uses.items()))
            callback_module = (
                f'module {module_name}\n'
                f'  use, intrinsic :: ISO_C_Binding, only : c_double,'
                f' c_f_procpointer, c_funptr, c_null_funptr\n'
                f'{use_lines}'
                f'  implicit none\n'
                f'  public\n'
                f'{iface_blocks}'
                f'{var_decls}'
                f'contains\n'
                f'{shim_blocks}'
                f'end module {module_name}\n'
            )
            # The module must precede the wrapper subroutine (appended
            # below) in the class's compilation unit.
            _shared_bucket(code, class_name).append(callback_module)

        # Generate Fortran method.  An array-return method is lowered to a
        # subroutine with an extra intent(out) array arg (see the
        # _ARRAY_RETURN_RX branch above), so result_is_subroutine forces
        # `subroutine` even though method_type is non-void.
        is_subroutine = method_type == 'void' or result_is_subroutine
        procedure = 'subroutine' if is_subroutine else 'function'
        func_decl = '' if is_subroutine else f'{method_type_c} :: {method_name_c}\n'

        iso_imports = iso_c_binding_import(arg_list, *isoImports.keys())
        fort_args = list(fortran_arg_list(arg_list)) + list(result_extra_fort_args)
        declarations = fortran_declarations(arg_list)
        reassignments = fortran_reassignments(arg_list)
        module_uses = fortran_module_uses(arg_list)
        call_lhs = ("call" if method_type == "void"
                    else f'{result_call_target} {result_assign_op}')
        call_code = fortran_call_code(arg_list,
                                     f'{call_lhs} {result_conversion_open} self_%{method_name}( &\n',
                                     f'&){result_conversion_close}\n', '&')
        # `result_pre_call_code` is for setup that must happen between the
        # arg-reassignments block and the inner call (currently: the
        # save-buffer deallocate for the dynamic-array-return path).  Goes
        # *before* call_code; `result_post_call_code` goes after as usual.
        call_code = result_pre_call_code + call_code + result_post_call_code

        fort_method = f'''{procedure} {method_name_c}({','.join(fort_args)}) bind(c,name='{method_name_c}')
  use :: {func_class.get('module')}, only : {class_name}Class
{module_uses}{result_extra_module_uses}{iso_imports}
  implicit none
{func_decl}{declarations}{result_extra_declarations}
{reassignments}{call_code}  return
end {procedure} {method_name_c}
'''
        _shared_bucket(code, class_name).append(fort_method)

        # Add c_lib interface
        arg_types = (list(ctypes_arg_types(arg_list))
                     + list(result_extra_clib_argtypes))
        restype = None if method_type == 'void' else clib_res_type
        python['c_lib'].append({
            'name': method_name_c,
            'restype': restype,
            'argtypes': arg_types,
        })

        # Generate Python method
        py_args = python_arg_list(arg_list)
        if result_python_class_wrap:
            # class(...) return: bypass python_call_code's "return c_lib(args)"
            # template — we need a setup statement (the c_int that ctypes
            # will fill), then the call, then a wrap-into-Python-class step.
            # Optional-arg branching is intentionally not supported here;
            # no current method needs it.
            parent_class, out_id_var = result_python_class_wrap
            py_call_args = []
            for a in arg_list:
                if not a.fort_is_present:
                    continue
                py_call_args.append(a.py_pass_as if a.py_pass_as else python_safe_name(a.name))
            py_call_args.append(f'byref({out_id_var})')
            reassignments_block = ''.join(a.py_reassignment for a in arg_list)
            py_call = (
                reassignments_block
                + f'    {out_id_var} = c_int(-1)\n'
                + f'    _ptr_ = c_lib.{method_name_c}({",".join(py_call_args)})\n'
                + f'    return {parent_class}._from_classID({out_id_var}.value, _ptr_)\n'
            )
        elif result_python_array_wrap:
            # Fixed-size array return: pre-allocate a numpy array of the
            # right shape/dtype, hand its data pointer as the synthetic
            # intent(out) arg, then return the (now-filled) numpy array.
            # As with class(...) returns, optional-arg branching isn't
            # supported here — none of today's array-return methods need it.
            arr_size, elem_ctype, elem_dtype = result_python_array_wrap
            py_call_args = []
            for a in arg_list:
                if not a.fort_is_present:
                    continue
                py_call_args.append(a.py_pass_as if a.py_pass_as else python_safe_name(a.name))
            py_call_args.append(
                f'_glcArr_.ctypes.data_as(POINTER({elem_ctype}))'
            )
            reassignments_block = ''.join(a.py_reassignment for a in arg_list)
            py_call = (
                reassignments_block
                + f'    _glcArr_ = np.zeros({arr_size}, dtype=np.{elem_dtype})\n'
                + f'    c_lib.{method_name_c}({",".join(py_call_args)})\n'
                + f'    return _glcArr_\n'
            )
        elif result_python_dyn_array_2d_wrap and not output_arrays:
            # 2D allocatable / dynamic-size array return.  Same shape as
            # the 1D path below, with one extra size companion and a
            # `reshape((s1, s2), order='F')` on the Python side so the
            # column-major Fortran storage maps to the right numpy
            # indexing.  (When output arrays are also present, the
            # output-array branch below handles the return's companions
            # too — they precede the per-argument ones in the signature.)
            elem_ctype, elem_dtype = result_python_dyn_array_2d_wrap
            py_call_args = []
            for a in arg_list:
                if not a.fort_is_present:
                    continue
                py_call_args.append(a.py_pass_as if a.py_pass_as else python_safe_name(a.name))
            py_call_args.append('byref(_glcDataPtr_)')
            py_call_args.append('byref(_glcSize1_)')
            py_call_args.append('byref(_glcSize2_)')
            reassignments_block = ''.join(a.py_reassignment for a in arg_list)
            py_call = (
                reassignments_block
                + f'    _glcDataPtr_ = c_void_p()\n'
                + f'    _glcSize1_   = c_size_t()\n'
                + f'    _glcSize2_   = c_size_t()\n'
                + f'    c_lib.{method_name_c}({",".join(py_call_args)})\n'
                # `from_address` over the save buffer; reshape in
                # column-major order to match the Fortran-side
                # (size1, size2) layout; `.copy()` materialises a fresh
                # numpy-owned array so the caller's reference survives
                # subsequent calls to the same method.
                + f'    return np.ctypeslib.as_array(\n'
                + f'        ({elem_ctype} * (_glcSize1_.value * _glcSize2_.value))'
                + f'.from_address(_glcDataPtr_.value)\n'
                + f'    ).reshape((_glcSize1_.value, _glcSize2_.value), order="F").copy()\n'
            )
        elif result_python_dyn_array_wrap and not output_arrays:
            # Allocatable / dynamic-size 1D array return.  The Fortran
            # wrapper writes the buffer's c_loc into a c_void_p out-param
            # and the element count into a c_size_t out-param; we wrap
            # them into a numpy view and *copy* before returning, because
            # the save-target buffer on the Fortran side is overwritten
            # on the next call to this method (mirrors the
            # `varying_string`-return convention: bytes are valid until
            # the next call).
            #
            # Optional-arg branching isn't supported on this path either;
            # the current set of dynamic-array-return methods have no
            # optional args.
            elem_ctype, elem_dtype = result_python_dyn_array_wrap
            py_call_args = []
            for a in arg_list:
                if not a.fort_is_present:
                    continue
                py_call_args.append(a.py_pass_as if a.py_pass_as else python_safe_name(a.name))
            py_call_args.append('byref(_glcDataPtr_)')
            py_call_args.append('byref(_glcSize_)')
            reassignments_block = ''.join(a.py_reassignment for a in arg_list)
            py_call = (
                reassignments_block
                + f'    _glcDataPtr_ = c_void_p()\n'
                + f'    _glcSize_    = c_size_t()\n'
                + f'    c_lib.{method_name_c}({",".join(py_call_args)})\n'
                # `from_address` builds a numpy array view over the
                # save buffer's bytes; `.copy()` materialises a fresh
                # numpy-owned array so the caller's reference survives
                # subsequent calls to the same method.
                + f'    return np.ctypeslib.as_array(\n'
                + f'        ({elem_ctype} * _glcSize_.value).from_address(_glcDataPtr_.value)\n'
                + f'    ).copy()\n'
            )
        elif output_arrays or sized_outputs:
            # Output-array method: one or more `intent(out), allocatable,
            # dimension(:)` args, each filled by the inner call into a
            # save-target buffer and conveyed back as a (c_ptr, c_size_t)
            # companion pair appended to the bind(c) signature above.  The
            # method may also carry scalar `intent(out)` numeric/logical
            # companions (is_output_scalar) — ordinary by-reference dummies
            # in the signature whose filled values are returned too.
            #
            # ctypes call order matches the bind(c) signature: fort_is_present
            # args first (regular inputs pass their value/pointer expr; scalar
            # outputs pass byref of a fresh ctype local), then the array
            # (c_ptr, c_size_t) companion pairs — mirroring
            # result_extra_fort_args.  After the call each array pair is
            # wrapped into a fresh numpy array (`.copy()`, since the save
            # buffer is overwritten on the next call; a null pointer / zero
            # size from the Fortran guard yields an empty array without
            # dereferencing).  The return interleaves scalar (.value) and
            # array outputs in declaration order — a single value bare, more
            # than one as a tuple.  The gate guarantees a void return and no
            # optional args, so no optional-arg branching is needed.
            setup_lines   = ''
            collect_lines = ''
            py_call_args  = []
            for a in arg_list:
                if not a.fort_is_present:
                    continue
                if a.is_output_scalar:
                    pv    = python_safe_name(a.name)
                    local = f'_{pv}_'
                    setup_lines += f'    {local} = {a.ctype or "c_int"}()\n'
                    py_call_args.append(f'byref({local})')
                elif a.is_output_sized:
                    # Sized output buffer: pre-allocate a flat numpy array
                    # of the product of the extent args (they're this
                    # function's own parameters) and pass its data pointer;
                    # the reshape-to-declared-shape happens in the collect
                    # block below.  No copy is needed — Python owns the
                    # buffer.  Extents were lowercased by the declaration
                    # parser; resolve them to the actual (case-preserved)
                    # parameter names case-insensitively.
                    name_map  = {x.name.lower(): python_safe_name(x.name)
                                 for x in arg_list}
                    extents   = [name_map.get(e.lower(), e)
                                 for e in a.output_extents]
                    pv        = python_safe_name(a.name)
                    arr_local = f'_{pv}Arr_'
                    size_expr = '*'.join(extents)
                    setup_lines += (
                        f'    {arr_local} = np.zeros(({size_expr},),'
                        f' dtype=np.{a.output_elem_dtype})\n')
                    py_call_args.append(
                        f'{arr_local}.ctypes.data_as('
                        f'POINTER({a.output_elem_ctype}))')
                    if len(extents) > 1:
                        shape = ', '.join(extents)
                        collect_lines += (
                            f'    {arr_local} = {arr_local}.reshape('
                            f'({shape}), order="F")\n')
                else:
                    py_call_args.append(a.py_pass_as if a.py_pass_as
                                        else python_safe_name(a.name))
            # A dynamic-array RETURN's companions were appended to
            # result_extra_fort_args by the return-type switch, BEFORE the
            # per-argument companions below — mirror that order here.
            ret_token = None
            if result_python_dyn_array_2d_wrap:
                elem_ctype, elem_dtype = result_python_dyn_array_2d_wrap
                py_call_args += ['byref(_glcDataPtr_)', 'byref(_glcSize1_)',
                                 'byref(_glcSize2_)']
                setup_lines += (
                    '    _glcDataPtr_ = c_void_p()\n'
                    '    _glcSize1_   = c_size_t()\n'
                    '    _glcSize2_   = c_size_t()\n'
                )
                collect_lines += (
                    f'    if _glcSize1_.value and _glcSize2_.value:\n'
                    f'        _glcRetArr_ = np.ctypeslib.as_array(\n'
                    f'            ({elem_ctype} * (_glcSize1_.value *'
                    f' _glcSize2_.value)).from_address(_glcDataPtr_.value)\n'
                    f'        ).reshape((_glcSize1_.value, _glcSize2_.value),'
                    f' order="F").copy()\n'
                    f'    else:\n'
                    f'        _glcRetArr_ = np.empty((_glcSize1_.value,'
                    f' _glcSize2_.value), dtype=np.{elem_dtype})\n'
                )
                ret_token = '_glcRetArr_'
            elif result_python_dyn_array_wrap:
                elem_ctype, elem_dtype = result_python_dyn_array_wrap
                py_call_args += ['byref(_glcDataPtr_)', 'byref(_glcSize_)']
                setup_lines += (
                    '    _glcDataPtr_ = c_void_p()\n'
                    '    _glcSize_    = c_size_t()\n'
                )
                collect_lines += (
                    f'    if _glcSize_.value:\n'
                    f'        _glcRetArr_ = np.ctypeslib.as_array(\n'
                    f'            ({elem_ctype} * _glcSize_.value)'
                    f'.from_address(_glcDataPtr_.value)\n'
                    f'        ).copy()\n'
                    f'    else:\n'
                    f'        _glcRetArr_ = np.empty(0,'
                    f' dtype=np.{elem_dtype})\n'
                )
                ret_token = '_glcRetArr_'
            elif method_type != 'void' and not result_is_subroutine:
                # Direct-scalar function result; ctypes converts the
                # restype to a plain Python scalar, so no .value.
                ret_token = '_glcRet_'
            for a in output_arrays:
                pv        = python_safe_name(a.name)
                ptr_local = f'_{pv}Ptr_'
                arr_local = f'_{pv}Arr_'
                setup_lines += f'    {ptr_local}  = c_void_p()\n'
                py_call_args.append(f'byref({ptr_local})')
                if a.array_rank == 2:
                    s1 = f'_{pv}Size1_'
                    s2 = f'_{pv}Size2_'
                    setup_lines += (f'    {s1} = c_size_t()\n'
                                    f'    {s2} = c_size_t()\n')
                    py_call_args += [f'byref({s1})', f'byref({s2})']
                    collect_lines += (
                        f'    if {s1}.value and {s2}.value:\n'
                        f'        {arr_local} = np.ctypeslib.as_array(\n'
                        f'            ({a.output_elem_ctype} * ({s1}.value *'
                        f' {s2}.value)).from_address({ptr_local}.value)\n'
                        f'        ).reshape(({s1}.value, {s2}.value),'
                        f' order="F").copy()\n'
                        f'    else:\n'
                        f'        {arr_local} = np.empty(({s1}.value,'
                        f' {s2}.value), dtype=np.{a.output_elem_dtype})\n'
                    )
                else:
                    size_local = f'_{pv}Size_'
                    setup_lines += f'    {size_local} = c_size_t()\n'
                    py_call_args.append(f'byref({size_local})')
                    collect_lines += (
                        f'    if {size_local}.value:\n'
                        f'        {arr_local} = np.ctypeslib.as_array(\n'
                        f'            ({a.output_elem_ctype} *'
                        f' {size_local}.value)'
                        f'.from_address({ptr_local}.value)\n'
                        f'        ).copy()\n'
                        f'    else:\n'
                        f'        {arr_local} = np.empty(0,'
                        f' dtype=np.{a.output_elem_dtype})\n'
                    )
            # Return tokens in declaration (arg_list) order so the tuple
            # matches the method's signature left-to-right; the function
            # result (scalar or dynamic-array) leads the tuple.
            result_tokens = [ret_token] if ret_token else []
            for a in arg_list:
                pv = python_safe_name(a.name)
                if a.is_output_scalar:
                    result_tokens.append(f'_{pv}_.value')
                elif a.is_output_array or a.is_output_sized:
                    result_tokens.append(f'_{pv}Arr_')
            reassignments_block = ''.join(a.py_reassignment for a in arg_list)
            call_lhs_py = '_glcRet_ = ' if ret_token == '_glcRet_' else ''
            if len(result_tokens) == 1:
                return_stmt = f'    return {result_tokens[0]}\n'
            else:
                return_stmt = f'    return ({", ".join(result_tokens)})\n'
            py_call = (
                reassignments_block
                + setup_lines
                + f'    {call_lhs_py}c_lib.{method_name_c}'
                + f'({",".join(py_call_args)})\n'
                + collect_lines
                + return_stmt
            )
        else:
            # Reassignments (numpy conversion, optional-arg unpacking, …)
            # belong before the call, mirroring the constructor template.
            # Without this, methods with array args would see an unconverted
            # input passed straight to data_as(), and methods with optional
            # functionClass args would skip their presence-check block.
            py_call = (python_reassignments(arg_list)
                       + python_call_code(arg_list,
                                          f'return c_lib.{method_name_c}'))
            if result_python_decode:
                # Append .decode("utf-8") to each call line so the bytes returned by
                # ctypes c_char_p are converted to a Python str.
                py_call = re.sub(
                    r'(c_lib\.' + re.escape(method_name_c) + r'\([^\n]*\))(\n)',
                    r'\1.decode("utf-8")\2',
                    py_call,
                )

        # Python-side identifier for the method.  Galacticus method names
        # (`yield`, `class`, ...) can collide with Python reserved words,
        # which would emit a `def yield(...)` and break import.  Follow
        # PEP 8's trailing-underscore convention to escape; the Fortran
        # call (`self_%{method_name}`) and the bind(c) symbol
        # (`{method_name_c}`) are unaffected.
        py_method_name = (method_name + '_'
                          if keyword.iskeyword(method_name)
                          else method_name)
        py_method = f'''def {py_method_name}({','.join(py_args)}):
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

    _shared_bucket(code, class_name).append(destructor_code)

    # Add c_lib interface
    python['c_lib'].append({
        'name': class_name + 'DestructorL',
        'restype': None,
        'argtypes': ['c_void_p', 'c_int'],
    })

    # Add Python destructor.  The _owned guard lets class(...)-returned
    # objects (which Galacticus owns and frees as part of their parent's
    # lifecycle) skip the destructor — see _from_classID in the parent
    # class body emitted by interfaces_python_classes.  Default-True via
    # getattr keeps every constructor-built object on the destroy path
    # without each constructor having to set the flag.
    #
    # The `_glcObj` guard handles partially-constructed objects: if
    # `__init__` raised (e.g. a size-validator ValueError before the
    # bind(c) constructor ran), Python still GCs the half-built
    # instance and calls __del__; without the guard we'd dereference
    # the missing attribute and surface an unrelated AttributeError
    # at interpreter shutdown.
    py_destructor = f'''# Destructor
def __del__(self):
    if getattr(self, '_owned', True) and hasattr(self, '_glcObj'):
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
  use :: Functions_Global_Utilities, only : Functions_Global_Set

  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Initialize HDF5 library access lock.
  call ioHDF5AccessInitialize()
  ! Connect global function pointers (mirrors Galacticus.exe startup) -
  ! required by any code path that crosses the deferred-binding seams
  ! between modules (e.g. merger-tree construction's state store).
  call Functions_Global_Set()
end subroutine libGalacticusInitL

subroutine libGalacticusNodesInitL() bind(c,name='libGalacticusNodesInitL')
  ! Initialize the node-component class system with default component
  ! selections - required before any tree nodes or merger trees can be
  ! created (mirrors the startup sequence of Galacticus.exe and of the
  ! unit tests). Kept separate from libGalacticusInitL because most
  ! library use never touches the node system, and this initialization
  ! selects component implementations. Idempotent.
  use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize
  use :: Node_Components , only : Node_Components_Initialize, Node_Components_Thread_Initialize
  use :: Input_Parameters, only : inputParameters
  type   (inputParameters), save :: parameters
  logical                 , save :: initialized=.false.

  if (initialized) return
  parameters=inputParameters()
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  initialized=.true.
end subroutine libGalacticusNodesInitL

subroutine libGalacticusVerbositySetL(verbosity) bind(c,name='libGalacticusVerbositySetL')
  ! Set the display verbosity level from a C string naming the level
  ! ("silent", "standard", "working", "warn", "info", "debug").  Exposed so
  ! library users can quiet Galacticus's progress bars and messages (which
  ! write directly to the process's stdout, bypassing Python-level stream
  ! redirection).
  use :: Display        , only : displayVerbositySet   , enumerationVerbosityLevelEncode
  use :: String_Handling, only : String_C_to_Fortran
  use :: ISO_Varying_String, only : char
  use, intrinsic :: ISO_C_Binding, only : c_char
  character(c_char), dimension(*), intent(in   ) :: verbosity

  call displayVerbositySet(enumerationVerbosityLevelEncode(char(String_C_to_Fortran(verbosity)),includesPrefix=.false.))
end subroutine libGalacticusVerbositySetL

program libGalacticusInit
end program libGalacticusInit
'''
    code['main'].append(init_code)


def _write_fortran_code(code, build_path):
    """Write generated Fortran code to files.

    Per class, the per-impl constructor wrappers are split into separate
    .F90 files (``<class>__<impl>.F90``) so that gfortran has a small
    compilation unit per impl.  This works around a memory-blow-up
    pathology that triggers OOM kills on classes with many impls
    (galacticFilter and nodePropertyExtractor were the original
    offenders, ~12 GB for a single .F90).  The class's shared pieces —
    GetPtr / GetIdAndPtr / methods / destructor — go to the original
    ``<class>.F90`` filename so existing dependency rules continue to
    apply to it; the per-impl files are picked up by
    libraryInterfacesDependencies.py via a directory listing.

    Stale per-impl files from a previous run (whose impl was since
    excluded or renamed) are cleaned up so they don't get linked into
    the .so.
    """
    out_dir = os.path.join(build_path, 'libgalacticus')
    os.makedirs(out_dir, exist_ok=True)

    written = set()

    if 'main' in code:
        # Written unconditionally, NOT only-if-changed: libgalacticus.Inc is
        # the target of the Makefile rule that runs this script, so its mtime
        # must advance past the rule's prerequisites or make would re-run the
        # (slow) generator on every invocation. The downstream re-preprocess
        # of this single file is cheap.
        main_file = os.path.join(build_path, 'libgalacticus.Inc')
        with open(main_file, 'w') as fh:
            fh.write('\n'.join(code['main']) + '\n')

    # The per-class/per-impl wrapper units are written only-if-changed
    # (mtime preserved when identical): their generated make rules run this
    # script only when the file is missing, so a stable mtime is safe here —
    # and it stops an unrelated catalog change from cascading into a
    # re-preprocess and recompile of every one of the ~1500 wrapper units.
    def _write_unit(path, lines):
        tmp = path + '.tmp'
        with open(tmp, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')
        file_changes_update(path, tmp)

    for class_name in sorted(k for k in code if k != 'main'):
        bucket = code[class_name]
        # Shared pieces (GetPtr, GetIdAndPtr, methods, destructor) go to
        # <class>.F90.
        _write_unit(os.path.join(out_dir, f'{class_name}.F90'),
                    bucket['shared'])
        written.add(f'{class_name}.F90')
        # One file per concrete impl's constructor wrapper.
        for impl_name, blocks in sorted(bucket['per_impl'].items()):
            _write_unit(os.path.join(out_dir, f'{class_name}__{impl_name}.F90'),
                        blocks)
            written.add(f'{class_name}__{impl_name}.F90')

    # Remove any stale .F90 (and matching .p.F90, .o, .d, .m) files that
    # weren't regenerated this run — e.g. an impl that picked up an
    # exclude="yes" or whose constructor newly fails the predicate would
    # otherwise leave a dangling object the linker still pulls in.
    for fname in os.listdir(out_dir):
        if (fname.endswith('.F90') and not fname.endswith('.p.F90')
                and fname not in written):
            stem = fname[:-len('.F90')]
            for ext in ('.F90', '.p.F90', '.p.F90.up', '.o', '.d', '.m'):
                stale = os.path.join(out_dir, stem + ext)
                if os.path.exists(stale):
                    os.remove(stale)



def _write_python_interface(python):
    """Write generated Python code to galacticus.py."""
    # Initialize the init unit
    init_content = '''from ctypes import *
import numpy as np
# Load the shared library into ctypes.
import os
cwd = os.getcwd()
libname = os.path.join(cwd, "galacticus/lib/libgalacticus.so")
c_lib = CDLL(libname)
c_lib.libGalacticusInitL()
c_lib.libGalacticusVerbositySetL.argtypes = [c_char_p]

def nodesInitialize():
    """Initialize the node-component class system (default component
    selections).  Required once before creating tree nodes or merger trees
    (e.g. before mergerTreeConstructor.construct); idempotent."""
    c_lib.libGalacticusNodesInitL()

def verbositySet(level):
    """Set Galacticus's display verbosity: one of 'silent', 'standard',
    'working', 'warn', 'info', or 'debug'.  Useful to quiet progress bars
    and messages, which Galacticus writes directly to the process's stdout
    (bypassing Python-level stream redirection)."""
    c_lib.libGalacticusVerbositySetL(level.encode('ascii'))
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

