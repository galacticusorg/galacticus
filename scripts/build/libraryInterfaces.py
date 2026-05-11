#!/usr/bin/env python3
# Generates C/Fortran/Python interfaces to Galacticus library classes.
# Andrew Benson (ported to Python with assistance from Claude 2026)

import sys
import os
import re
import keyword
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
    _SHARED_TYPE_MODULES,
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

# Match `type(enumerationFooType)` — the Galacticus convention for an
# enumeration kind named "foo".  Used by interfaces_methods to detect
# enumeration return types and lift the inner %ID component out as c_int.
_ENUM_RETURN_RX = re.compile(
    r'^type\s*\(\s*(enumeration[a-z0-9_]+type)\s*\)$',
    re.IGNORECASE,
)

# Match `class(fooClass)` — the Galacticus convention for a polymorphic
# pointer to a registered functionClass.  Used by interfaces_methods to
# detect class(...) return types and route them through the per-class
# XGetIdAndPtr helper.
_CLASS_RETURN_RX = re.compile(
    r'^class\s*\(\s*([a-z][a-zA-Z0-9_]*Class)\s*\)$',
    re.IGNORECASE,
)

# Match a fixed-size 1D dimension(N) attribute on a constructor/method arg.
# Mirrors Pipeline._DIM_FIXED_RX (kept locally so libraryInterfaces.py
# stays usable standalone) — accepts both `dimension(N)` and the
# explicit-lower-bound form `dimension(L:U)`.  Only the *match* is
# checked here; this file doesn't compute the element count.
# doesn't need to import a private name).
_DIM_FIXED_RX_INTERFACES = re.compile(
    r'^dimension\s*\(\s*(?:\d+\s*:\s*)?\d+\s*\)$'
)

# Match `len=N` literal in a character type-spec.  Same regex as
# Pipeline.py's `_CHAR_LEN_RX`; duplicated here to keep this module
# self-contained for the validation pass that runs before the pipeline.
_CHAR_LEN_RX_INTERFACES = re.compile(r'^len\s*=\s*(\d+)$')

# Match a 1D fixed-size array RETURN type, e.g. `double precision, dimension(3)`
# or `integer, dimension(3)`.  Captures (intrinsic, size).
_ARRAY_RETURN_RX = re.compile(
    r'^(double\s+precision|integer)\s*,\s*dimension\s*\(\s*(\d+)\s*\)\s*$',
    re.IGNORECASE,
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
            if (isinstance(mod_data, dict)
                    and enum_type in mod_data.get('only', {})):
                return mod_name
    return func_class.get('module')


def _unsupported_arg(arg, lib_function_classes, *,
                     constructor_overrides=(), class_hierarchy=None):
    """Return a human-readable reason if ``arg`` has a type the pipeline
    can't translate, otherwise ``None``.

    Shared core used by both :func:`_unsupported_constructor_arg` and
    :func:`_unsupported_method_arg`.

    Always rejects:

    * ``complex`` / ``double complex`` — ctypes has no built-in c_complex.
    * ``class(FooClass)`` whose stem ``Foo`` is not listed in
      ``libraryClasses.xml``.  Without that registration, ``assign_c_types``
      can't set ``is_function_class=True``, so no ``_ID`` companion or
      ``FooGetPtr``-based reassignment is generated and the bind(c) wrapper
      passes a raw c_ptr into a Fortran routine that wants ``class(...)``.
    * ``class(SomethingElse)`` that isn't a registered functionClass — same
      mismatch (e.g. ``class(nodeComponentBlackHole)``).
    * ``class(*)`` unless ``constructor_overrides`` contains a matching
      ``<argument name=... type=... module=.../>`` hint from
      libraryClasses.xml.
    * Any source-level ``dimension(...)`` other than 1D deferred-shape
      ``dimension(:)`` of double precision / integer (the array path
      plumbed through assign_c_types et al. — see Pipeline.py).  And even
      among those, ``allocatable`` and ``intent(out)`` / ``intent(inout)``
      array args are rejected: those are output arrays whose shape is
      determined by the inner Galacticus call, which needs a different
      bind(c) shape (allocatable can't be assumed-size, the count has to
      come back, etc.) than the input-array path does.
    """
    intrinsic = arg.get('intrinsic')
    # A `<argument name="..." value="null"/>` override in libraryClasses.xml
    # tells the wrapper to drop the arg and pass a local null pointer to
    # the inner constructor instead — sufficient for the
    # procedure-pointer + class(*) callback-injection escape hatch the
    # parameter-driven paths of some impls already null out (see
    # assign_c_types / build_fortran_reassignments).  Accept it for
    # those two intrinsics regardless of the rest of this predicate.
    if any(isinstance(o, dict)
           and o.get('name') == arg.get('name')
           and o.get('value') == 'null'
           for o in constructor_overrides):
        type_spec = (arg.get('type') or '').strip()
        if intrinsic == 'procedure' \
                or (intrinsic == 'class' and type_spec == '*'):
            return None
        return (f"value='null' override on {intrinsic}({type_spec}) — "
                f"only procedure and class(*) args are supported")
    if intrinsic in ('complex', 'double complex'):
        return f"{intrinsic}({arg.get('type','')})"
    if intrinsic == 'procedure':
        # Procedure-pointer args (e.g. `procedure(integrand) :: f`) have
        # no ctypes counterpart in this pipeline — assign_c_types has no
        # branch for them, so without an early reject the arg falls
        # through with empty ctype/fort_type and the bind(c) wrapper
        # emits a broken declaration that mismatches the inner method
        # signature.  Skip the surrounding constructor or method instead.
        return f"procedure({arg.get('type','')}) — procedure-pointer args are not supported"
    if intrinsic == 'integer' \
            and (arg.get('type') or '').strip() == 'omp_lock_kind' \
            and 'optional' not in arg.get('attributes', []):
        # `integer(omp_lock_kind)` has a platform-dependent size
        # (INTEGER(4) on Linux GCC vs INTEGER(8) on macOS GCC), so the
        # wrapper can't pick a single C-interop kind that's correct on
        # both.  When the arg is `optional`, assign_c_types drops it
        # silently and the inner method's optional default takes over.
        # When it's non-optional we have to reject the surrounding
        # method — silently dropping a required arg would mis-call the
        # inner Galacticus routine.
        return ("integer(omp_lock_kind) non-optional argument — kind is "
                "platform-dependent (INTEGER(4) on Linux vs INTEGER(8) on "
                "macOS), no portable C-interop kind exists")
    if intrinsic == 'class':
        type_spec = (arg.get('type') or '').strip()
        if type_spec == '*':
            has_override = any(
                isinstance(o, dict) and o.get('name') == arg.get('name')
                for o in constructor_overrides
            )
            if not has_override:
                return ("class(*) without a libraryClasses.xml override "
                        "(<constructor><argument name=… type=… module=…/>)")
        elif type_spec.endswith('Class'):
            stem = type_spec[:-5]
            if stem not in lib_function_classes:
                return (f"class({type_spec}) — '{stem}' is not a registered "
                        f"functionClass in libraryClasses.xml")
        else:
            # Abstract intermediate: a class whose extends-chain reaches
            # a registered <base>Class is acceptable — the pipeline routes
            # through <base>GetPtr and narrows to <intermediate> at
            # runtime (see build_fortran_reassignments).
            base = None
            if class_hierarchy:
                from LibraryInterfaces.Hierarchy import resolve_function_class_base
                base, _ = resolve_function_class_base(
                    type_spec, class_hierarchy, set(lib_function_classes.keys()))
            if base is None:
                return (f"class({type_spec}) — only registered functionClasses "
                        f"and class(*) (with override) are supported")
    # Dimension check — same predicate for constructor and method args now
    # that both use the array-arg pipeline.  (Previously constructor-only,
    # which left method args producing broken-but-not-fatal Fortran when
    # they had array shapes.  After array support landed, *some* of those
    # method args now compile, but the cases that don't fit our
    # input-only pipeline — allocatable / intent(out) — produce hard
    # compile errors instead of silently broken code, so we now reject
    # them cleanly here for both contexts.)
    attrs = arg.get('attributes', [])
    for attr in attrs:
        if attr.startswith('dimension'):
            # 1D numeric arrays are plumbed through with a numpy ndarray
            # conversion at the Python boundary; deferred-shape gets a
            # count companion (see assign_c_types /
            # build_python_reassignments in Pipeline.py), fixed-size
            # uses the literal length from the dimension spec.  Output /
            # allocatable arrays still need the bind(c) function to
            # communicate the size BACK to Python plus a non-assumed-size
            # shape — that's a separate implementation we don't have yet.
            is_supported_dim = (
                intrinsic in ('double precision', 'integer')
                and (attr == 'dimension(:)'
                     or attr == 'dimension(:,:)'
                     or _DIM_FIXED_RX_INTERFACES.match(attr))
            )
            # 1D deferred-shape logical: same plumbing as numeric 1D
            # arrays, except build_fortran_reassignments emits a
            # kind-narrowing copy from `logical(c_bool)` to the default
            # `logical` kind that the inner method's dummy declares.
            is_supported_logical = (
                intrinsic == 'logical' and attr == 'dimension(:)'
            )
            # Fixed-length character arrays (`character(len=N),
            # dimension(:)`) are supported at deferred shape: the bind(c)
            # boundary receives a contiguous count*N byte buffer plus a
            # count companion, and the wrapper repacks into
            # `character(len=N), dimension(:)` for the inner call.
            is_supported_char_array = (
                intrinsic == 'character'
                and attr == 'dimension(:)'
                and _CHAR_LEN_RX_INTERFACES.match((arg.get('type') or '').strip())
            )
            # `type(varying_string), dimension(:)`: shipped via a flat
            # contiguous byte buffer plus count + per-element-length
            # companions; the wrapper builds a `type(varying_string),
            # dimension(:)` for the inner call (Pipeline.py's
            # build_python_reassignments / build_fortran_reassignments).
            is_supported_vstring_array = (
                intrinsic == 'type'
                and (arg.get('type') or '').strip() == 'varying_string'
                and attr == 'dimension(:)'
            )
            # `type(<class>List), dimension(:)`: Galacticus's idiom for
            # an array of polymorphic pointers — each element is a
            # wrapper struct whose `<class>_` component holds a
            # `class(<class>Class), pointer`.  We support these when
            # `<class>` is a registered functionClass; the wrapper ships
            # parallel arrays of c_ptr + classID and rebuilds the list
            # via `<class>GetPtr` element-by-element (same machinery as
            # the scalar `class(...)` arg path).
            type_spec = (arg.get('type') or '').strip()
            # The shape is only recognised when the wrapper type is
            # registered in `_SHARED_TYPE_MODULES` (Pipeline.py's
            # gatekeeper); see the corresponding branch in
            # `assign_c_types` for the rationale.  Other locally-defined
            # `<class>List` structs that happen to share the naming
            # convention fall through to the generic-rejection path.
            is_supported_list_array = (
                intrinsic == 'type'
                and attr == 'dimension(:)'
                and type_spec.endswith('List')
                and type_spec[:-4] in lib_function_classes
                and type_spec in _SHARED_TYPE_MODULES
            )
            if (is_supported_dim or is_supported_logical
                    or is_supported_char_array
                    or is_supported_vstring_array or is_supported_list_array):
                if 'allocatable' in attrs:
                    return ('1D allocatable array argument'
                            ' (output arrays not yet supported)')
                if 'intent(in)' not in attrs:
                    return ('non-intent(in) array argument'
                            ' (output arrays not yet supported)')
                continue
            return ('dimensioned argument'
                    ' (only 1D deferred-shape or fixed-size numeric input,'
                    ' 2D deferred-shape numeric input, or 1D deferred-shape'
                    ' fixed-length character arrays, are supported)')
    return None


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
    rejected."""
    if class_hierarchy is None:
        class_hierarchy = _CLASS_HIERARCHY
    for arg in args:
        reason = _unsupported_arg(arg, lib_function_classes,
                                  class_hierarchy=class_hierarchy)
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
        ambiguous_internals = None
        args_constructor = []

        for node in SourceTree.walk_tree(tree):
            if node['type'] == class_name:
                impl_name = node.get('directive', {}).get('name')
                is_abstract = (node.get('directive', {}).get('abstract', 'no') == 'yes')

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
                    # InternalType vs InternalDefined).  Stash for the
                    # warning-and-skip below; leave name_constructor unset
                    # so the impl is dropped from impls_list.
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
                    args_constructor = [{'name': re.sub(r'[\s&]+', '', a)}
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
        result_assign_op          = '='   # `=>` for pointer-returning methods (class(...))
        result_python_decode      = False
        result_extra_fort_args    = []    # extra args appended to bind(c) signature
        result_extra_clib_argtypes = []   # matching ctypes wrappers
        result_python_class_wrap  = None  # (parent_class, out_id_var) or None
        result_python_array_wrap  = None  # (size, elem_ctype, elem_dtype) or None
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

        # Process arguments
        arg_list = assign_c_types(arg_list, lib_function_classes or {},
                                  class_hierarchy=_CLASS_HIERARCHY)
        arg_list = assign_c_attributes(arg_list)
        arg_list = build_python_reassignments(arg_list)
        arg_list = build_fortran_reassignments(arg_list, func_class, None,
                                              extensions, module_uses_impls,
                                              lib_function_classes)

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
        call_code += result_post_call_code

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
        main_file = os.path.join(build_path, 'libgalacticus.Inc')
        with open(main_file, 'w') as fh:
            fh.write('\n'.join(code['main']) + '\n')

    for class_name in sorted(k for k in code if k != 'main'):
        bucket = code[class_name]
        # Shared pieces (GetPtr, GetIdAndPtr, methods, destructor) go to
        # <class>.F90.
        shared_file = os.path.join(out_dir, f'{class_name}.F90')
        with open(shared_file, 'w') as fh:
            fh.write('\n'.join(bucket['shared']) + '\n')
        written.add(f'{class_name}.F90')
        # One file per concrete impl's constructor wrapper.
        for impl_name, blocks in sorted(bucket['per_impl'].items()):
            impl_file = os.path.join(out_dir, f'{class_name}__{impl_name}.F90')
            with open(impl_file, 'w') as fh:
                fh.write('\n'.join(blocks) + '\n')
            written.add(f'{class_name}__{impl_name}.F90')

    # Remove any stale .F90 (and matching .p.F90, .o, .d, .m) files that
    # weren't regenerated this run — e.g. an impl that picked up an
    # exclude="yes" or whose constructor newly fails the predicate would
    # otherwise leave a dangling object the linker still pulls in.
    for fname in os.listdir(out_dir):
        if fname.endswith('.F90') and fname not in written:
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

