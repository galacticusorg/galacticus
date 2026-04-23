#!/usr/bin/env python3
# Golden-output regression tests for the LibraryInterfaces pipeline and emitters.
# Tests the following modules:
#   - LibraryInterfaces.ArgSpec     (ArgSpec dataclass / from_raw)
#   - LibraryInterfaces.Pipeline    (assign_c_types, assign_c_attributes,
#                                    build_python_reassignments,
#                                    build_fortran_reassignments)
#   - LibraryInterfaces.Emitters    (ctypes_arg_types, fortran_arg_list,
#                                    fortran_declarations, fortran_reassignments,
#                                    fortran_module_uses, fortran_call_code,
#                                    iso_c_binding_import, python_arg_list,
#                                    python_reassignments, python_call_code)

import sys
import os

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', '.'), 'python'))

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

# ============================================================================
# Test framework (mirrors test-python-utils.py)
# ============================================================================

PASS_COUNT = 0
FAIL_COUNT = 0


def assert_equal(actual, expected, msg):
    global PASS_COUNT, FAIL_COUNT
    if actual == expected:
        PASS_COUNT += 1
        print(f"✓ {msg}")
    else:
        FAIL_COUNT += 1
        print(f"✗ {msg}")
        print(f"  Expected: {expected!r}")
        print(f"  Got:      {actual!r}")


def assert_in(needle, haystack, msg):
    global PASS_COUNT, FAIL_COUNT
    if needle in haystack:
        PASS_COUNT += 1
        print(f"✓ {msg}")
    else:
        FAIL_COUNT += 1
        print(f"✗ {msg}")
        print(f"  Expected {needle!r} to be in {haystack!r}")


def assert_not_in(needle, haystack, msg):
    global PASS_COUNT, FAIL_COUNT
    if needle not in haystack:
        PASS_COUNT += 1
        print(f"✓ {msg}")
    else:
        FAIL_COUNT += 1
        print(f"✗ {msg}")
        print(f"  Expected {needle!r} NOT to be in {haystack!r}")


# ============================================================================
# Helpers
# ============================================================================

def run_pipeline(raw_args, lib_fc=None, func_class=None,
                 implementation=None, extensions=None, module_uses_impls=None):
    """Run all four pipeline stages and return the enriched ArgSpec list."""
    if lib_fc is None:
        lib_fc = {}
    if func_class is None:
        func_class = {}
    if extensions is None:
        extensions = {}
    if module_uses_impls is None:
        module_uses_impls = {}
    al = assign_c_types(raw_args, lib_fc)
    al = assign_c_attributes(al)
    al = build_python_reassignments(al)
    al = build_fortran_reassignments(al, func_class, implementation,
                                     extensions, module_uses_impls, lib_fc)
    return al


# ============================================================================
# ArgSpec tests
# ============================================================================

def test_argspec_from_raw():
    print("\n=== Testing ArgSpec.from_raw() ===")

    d = {'name': 'omegaMatter', 'intrinsic': 'double precision',
         'type': None, 'attributes': ['intent(in)']}
    arg = ArgSpec.from_raw(d)
    assert_equal(arg.name, 'omegaMatter', "from_raw: name")
    assert_equal(arg.intrinsic, 'double precision', "from_raw: intrinsic")
    assert_equal(arg.type_spec, '', "from_raw: type_spec for None type")
    assert_equal(arg.attributes, ['intent(in)'], "from_raw: attributes")

    # Defaults are set correctly
    assert_equal(arg.ctype, '', "from_raw: ctype empty initially")
    assert_equal(arg.fort_is_present, True, "from_raw: fort_is_present default True")
    assert_equal(arg.py_is_present, True, "from_raw: py_is_present default True")
    assert_equal(arg.galacticus_is_present, True, "from_raw: galacticus_is_present default True")
    assert_equal(arg.fort_modules, {}, "from_raw: fort_modules empty dict")
    assert_equal(arg.fort_iso_c_symbols, [], "from_raw: fort_iso_c_symbols empty list")

    # type field is mapped to type_spec
    d2 = {'name': 'node', 'intrinsic': 'type', 'type': 'treeNode', 'attributes': []}
    arg2 = ArgSpec.from_raw(d2)
    assert_equal(arg2.type_spec, 'treeNode', "from_raw: type_spec from 'type' key")


# ============================================================================
# assign_c_types tests
# ============================================================================

def test_assign_c_types_scalar_intrinsics():
    print("\n=== Testing assign_c_types: scalar intrinsics ===")

    cases = [
        ({'name': 'x', 'intrinsic': 'double precision', 'type': None,
          'attributes': ['intent(in)']},
         'c_double', 'real(c_double)', "double precision"),
        ({'name': 'n', 'intrinsic': 'integer', 'type': None,
          'attributes': ['intent(in)']},
         'c_int', 'integer(c_int)', "integer"),
        ({'name': 'flag', 'intrinsic': 'logical', 'type': None,
          'attributes': ['intent(in)']},
         'c_bool', 'logical(c_bool)', "logical"),
        ({'name': 'label', 'intrinsic': 'character', 'type': None,
          'attributes': ['intent(in)']},
         'c_char_p', 'character(c_char)', "character"),
    ]
    for raw, exp_ctype, exp_ftype, label in cases:
        al = assign_c_types([raw], {})
        assert_equal(len(al), 1, f"{label}: single arg produced")
        assert_equal(al[0].ctype, exp_ctype, f"{label}: ctype")
        assert_equal(al[0].fort_type, exp_ftype, f"{label}: fort_type")
        assert_equal(al[0].is_optional, False, f"{label}: not optional")
        assert_equal(al[0].is_function_class, False, f"{label}: not functionClass")


def test_assign_c_types_type_variants():
    print("\n=== Testing assign_c_types: derived types ===")

    # varying_string
    raw = {'name': 'desc', 'intrinsic': 'type', 'type': 'varying_string',
           'attributes': ['intent(in)']}
    al = assign_c_types([raw], {})
    assert_equal(al[0].ctype, 'c_char_p', "varying_string: ctype")
    assert_equal(al[0].fort_type, 'character(c_char)', "varying_string: fort_type")

    # enumeration type
    raw = {'name': 'method', 'intrinsic': 'type',
           'type': 'enumerationIntegrationMethodType',
           'attributes': ['intent(in)']}
    al = assign_c_types([raw], {})
    assert_equal(al[0].ctype, 'c_int', "enumeration: ctype")
    assert_equal(al[0].fort_type, 'integer(c_int)', "enumeration: fort_type")

    # other derived type (treeNode)
    raw = {'name': 'node', 'intrinsic': 'type', 'type': 'treeNode',
           'attributes': ['intent(inout)']}
    al = assign_c_types([raw], {})
    assert_equal(al[0].ctype, 'c_void_p', "treeNode: ctype")
    assert_equal(al[0].fort_type, 'type(c_ptr)', "treeNode: fort_type")


def test_assign_c_types_optional():
    print("\n=== Testing assign_c_types: optional args ===")

    raw = {'name': 'time', 'intrinsic': 'double precision', 'type': None,
           'attributes': ['intent(in)', 'optional']}
    al = assign_c_types([raw], {})
    assert_equal(al[0].is_optional, True, "optional: is_optional flag set")


def test_assign_c_types_function_class():
    print("\n=== Testing assign_c_types: functionClass arg ===")

    lib_fc = {'cosmologyParameters': {'module': 'Cosmology_Parameters'}}
    raw = {'name': 'cosmology', 'intrinsic': 'class',
           'type': 'cosmologyParametersClass', 'attributes': ['intent(in)']}
    al = assign_c_types([raw], lib_fc)

    # Two args: the original + the _ID companion
    assert_equal(len(al), 2, "functionClass: two args produced")
    main_arg = al[0]
    id_arg   = al[1]

    assert_equal(main_arg.name, 'cosmology', "functionClass: main arg name")
    assert_equal(main_arg.is_function_class, True, "functionClass: is_function_class=True")
    assert_equal(main_arg.ctype, 'c_void_p', "functionClass: main ctype")
    assert_equal(main_arg.fort_type, 'type(c_ptr)', "functionClass: main fort_type")

    assert_equal(id_arg.name, 'cosmology_ID', "functionClass: _ID arg name")
    assert_equal(id_arg.ctype, 'c_int', "functionClass: _ID ctype")
    assert_equal(id_arg.fort_type, 'integer(c_int)', "functionClass: _ID fort_type")
    assert_equal(id_arg.py_is_present, False, "functionClass: _ID py_is_present=False")
    assert_equal(id_arg.galacticus_is_present, False, "functionClass: _ID galacticus_is_present=False")


def test_assign_c_types_self_exclusion():
    print("\n=== Testing assign_c_types: 'self' arg galacticus exclusion ===")

    lib_fc = {'cosmologyParameters': {'module': 'Cosmology_Parameters'}}
    raw = {'name': 'self', 'intrinsic': 'class',
           'type': 'cosmologyParametersClass', 'attributes': ['intent(inout)']}
    al = assign_c_types([raw], lib_fc)
    self_arg = al[0]
    assert_equal(self_arg.galacticus_is_present, False,
                 "self functionClass: galacticus_is_present=False")


def test_assign_c_types_optional_function_class():
    print("\n=== Testing assign_c_types: optional functionClass arg ===")

    lib_fc = {'cosmologyParameters': {'module': 'Cosmology_Parameters'}}
    raw = {'name': 'cosmo', 'intrinsic': 'class',
           'type': 'cosmologyParametersClass', 'attributes': ['intent(in)', 'optional']}
    al = assign_c_types([raw], lib_fc)
    assert_equal(len(al), 2, "optional functionClass: two args produced")
    id_arg = al[1]
    assert_equal(id_arg.is_optional, True, "optional functionClass: _ID is_optional=True")
    assert_equal(id_arg.py_present, 'cosmo', "optional functionClass: _ID py_present set to parent name")


# ============================================================================
# assign_c_attributes tests
# ============================================================================

def test_assign_c_attributes():
    print("\n=== Testing assign_c_attributes ===")

    # intent(in) double → pass-by-value, no pointer wrapper
    raw = {'name': 'x', 'intrinsic': 'double precision', 'type': None,
           'attributes': ['intent(in)']}
    al = assign_c_types([raw], {})
    al = assign_c_attributes(al)
    assert_equal(al[0].pass_by, 'value', "intent(in) double: pass_by=value")
    assert_equal('value' in al[0].fort_attributes, True, "intent(in) double: 'value' in fort_attributes")
    assert_equal(al[0].ctype_pointer, False, "intent(in) double: ctype_pointer=False")

    # intent(out) double → pass-by-reference, pointer wrapper
    raw = {'name': 'y', 'intrinsic': 'double precision', 'type': None,
           'attributes': ['intent(out)']}
    al = assign_c_types([raw], {})
    al = assign_c_attributes(al)
    assert_equal(al[0].pass_by, 'reference', "intent(out) double: pass_by=reference")
    assert_equal(al[0].ctype_pointer, True, "intent(out) double: ctype_pointer=True")
    assert_equal('value' in al[0].fort_attributes, False, "intent(out) double: no 'value' attribute")

    # optional arg → pass-by-reference
    raw = {'name': 'z', 'intrinsic': 'double precision', 'type': None,
           'attributes': ['intent(in)', 'optional']}
    al = assign_c_types([raw], {})
    al = assign_c_attributes(al)
    assert_equal(al[0].pass_by, 'reference', "optional double: pass_by=reference")
    assert_equal('optional' in al[0].fort_attributes, True, "optional double: 'optional' in fort_attributes")

    # character → dimension(*) appended, pass-by-value (c_char_p)
    raw = {'name': 'label', 'intrinsic': 'character', 'type': None,
           'attributes': ['intent(in)']}
    al = assign_c_types([raw], {})
    al = assign_c_attributes(al)
    assert_equal('dimension(*)' in al[0].fort_attributes, True, "character: dimension(*) added")
    assert_equal(al[0].ctype_pointer, False, "character: ctype_pointer=False (c_char_p excluded)")


# ============================================================================
# build_python_reassignments tests
# ============================================================================

def test_build_python_reassignments_non_optional():
    print("\n=== Testing build_python_reassignments: non-optional functionClass ===")

    lib_fc = {'cosmologyParameters': {'module': 'Cosmology_Parameters'}}
    raw = {'name': 'cosmology', 'intrinsic': 'class',
           'type': 'cosmologyParametersClass', 'attributes': ['intent(in)']}
    al = assign_c_types([raw], lib_fc)
    al = assign_c_attributes(al)
    al = build_python_reassignments(al)

    main_arg = al[0]
    id_arg   = al[1]
    assert_equal(main_arg.py_pass_as, 'cosmology._glcObj',  "non-opt FC: py_pass_as of main arg")
    assert_equal(id_arg.py_pass_as,   'cosmology._classID', "non-opt FC: py_pass_as of _ID arg")
    assert_equal(main_arg.py_reassignment, '', "non-opt FC: no py_reassignment for main")


def test_build_python_reassignments_optional():
    print("\n=== Testing build_python_reassignments: optional functionClass ===")

    lib_fc = {'cosmologyParameters': {'module': 'Cosmology_Parameters'}}
    raw = {'name': 'cosmo', 'intrinsic': 'class',
           'type': 'cosmologyParametersClass', 'attributes': ['intent(in)', 'optional']}
    al = assign_c_types([raw], lib_fc)
    al = assign_c_attributes(al)
    al = build_python_reassignments(al)

    main_arg = al[0]
    id_arg   = al[1]
    assert_equal(main_arg.py_pass_as, 'cosmo_glcObj',  "opt FC: py_pass_as of main arg")
    assert_equal(id_arg.py_pass_as,   'cosmo_classID', "opt FC: py_pass_as of _ID arg")
    assert_in('if cosmo:', main_arg.py_reassignment, "opt FC: if-branch in py_reassignment")
    assert_in('cosmo_glcObj', main_arg.py_reassignment, "opt FC: glcObj var in py_reassignment")
    assert_in('else:', main_arg.py_reassignment, "opt FC: else-branch in py_reassignment")
    assert_in('None', main_arg.py_reassignment, "opt FC: None in absent branch")


# ============================================================================
# build_fortran_reassignments tests
# ============================================================================

def test_fort_reassign_logical():
    print("\n=== Testing build_fortran_reassignments: logical ===")

    fc = {'name': 'foo', 'module': 'Foo_Module'}
    raw = {'name': 'flag', 'intrinsic': 'logical', 'type': None,
           'attributes': ['intent(in)']}
    al = run_pipeline([raw], func_class=fc)

    assert_equal(al[0].fort_pass_as,      'flag_',                  "logical: fort_pass_as")
    assert_equal(al[0].fort_reassignment, 'flag_=logical(flag)\n',  "logical: fort_reassignment")
    assert_equal(al[0].fort_declarations, 'logical :: flag_\n',     "logical: fort_declarations")


def test_fort_reassign_character():
    print("\n=== Testing build_fortran_reassignments: character ===")

    fc = {'name': 'foo', 'module': 'Foo_Module'}
    raw = {'name': 'label', 'intrinsic': 'character', 'type': None,
           'attributes': ['intent(in)']}
    al = run_pipeline([raw], func_class=fc)

    assert_equal(al[0].fort_pass_as, 'char(String_C_to_Fortran(label))', "character: fort_pass_as")
    assert_in('String_Handling',    al[0].fort_modules, "character: String_Handling module")
    assert_in('ISO_Varying_String', al[0].fort_modules, "character: ISO_Varying_String module")


def test_fort_reassign_varying_string():
    print("\n=== Testing build_fortran_reassignments: varying_string ===")

    fc = {'name': 'foo', 'module': 'Foo_Module'}
    raw = {'name': 'desc', 'intrinsic': 'type', 'type': 'varying_string',
           'attributes': ['intent(in)']}
    al = run_pipeline([raw], func_class=fc)

    assert_equal(al[0].fort_pass_as, 'String_C_to_Fortran(desc)', "varying_string: fort_pass_as")
    assert_in('String_Handling', al[0].fort_modules, "varying_string: String_Handling module")


def test_fort_reassign_enumeration():
    print("\n=== Testing build_fortran_reassignments: enumeration ===")

    fc = {
        'name': 'integratorODE',
        'module': 'Integrator_ODE',
        'moduleUses': [
            {'Integrator_ODE_Enumerations': {
                'only': {'enumerationIntegrationMethodType': 1}
            }}
        ],
    }
    raw = {'name': 'method', 'intrinsic': 'type',
           'type': 'enumerationIntegrationMethodType',
           'attributes': ['intent(in)']}
    al = run_pipeline([raw], func_class=fc)

    assert_equal(al[0].fort_pass_as, 'method_', "enumeration: fort_pass_as")
    assert_in('method_%ID=method\n', al[0].fort_reassignment, "enumeration: %ID reassignment")
    assert_in('type(enumerationIntegrationMethodType) :: method_\n',
              al[0].fort_declarations, "enumeration: fort_declarations")
    assert_in('Integrator_ODE_Enumerations', al[0].fort_modules, "enumeration: module located via moduleUses")


def test_fort_reassign_tree_node():
    print("\n=== Testing build_fortran_reassignments: treeNode ===")

    fc = {'name': 'foo', 'module': 'Foo_Module'}
    raw = {'name': 'node', 'intrinsic': 'type', 'type': 'treeNode',
           'attributes': ['intent(inout)']}
    al = run_pipeline([raw], func_class=fc)

    assert_in('type(treeNode), pointer :: node_\n', al[0].fort_declarations, "treeNode: declaration")
    assert_equal(al[0].fort_pass_as, 'node_', "treeNode: fort_pass_as")
    assert_in('c_f_pointer', al[0].fort_reassignment, "treeNode: c_f_pointer in reassignment")
    assert_in('c_f_pointer', al[0].fort_iso_c_symbols, "treeNode: c_f_pointer in iso_c_symbols")
    assert_in('Galacticus_Nodes', al[0].fort_modules, "treeNode: Galacticus_Nodes module")
    # Not optional: no if/else block
    assert_not_in('null()', al[0].fort_reassignment, "treeNode non-optional: no null() fallback")


def test_fort_reassign_tree_node_optional():
    print("\n=== Testing build_fortran_reassignments: optional treeNode ===")

    fc = {'name': 'foo', 'module': 'Foo_Module'}
    raw = {'name': 'node', 'intrinsic': 'type', 'type': 'treeNode',
           'attributes': ['intent(inout)', 'optional']}
    al = run_pipeline([raw], func_class=fc)

    assert_in('null()', al[0].fort_reassignment, "optional treeNode: null() fallback")
    assert_in('if (present(node)) then', al[0].fort_reassignment, "optional treeNode: if-present guard")


def test_fort_reassign_function_class():
    print("\n=== Testing build_fortran_reassignments: functionClass ===")

    lib_fc = {'cosmologyParameters': {'module': 'Cosmology_Parameters'}}
    fc     = {'name': 'cosmologyParameters', 'module': 'Cosmology_Parameters'}
    raw = {'name': 'cosmology', 'intrinsic': 'class',
           'type': 'cosmologyParametersClass', 'attributes': ['intent(in)']}
    al = run_pipeline([raw], lib_fc=lib_fc, func_class=fc)

    fc_arg = al[0]
    assert_equal(fc_arg.fort_pass_as, 'cosmology_', "FC: fort_pass_as")
    assert_in('class(cosmologyParametersClass), pointer :: cosmology_\n',
              fc_arg.fort_declarations, "FC: class pointer declaration")
    assert_in('cosmologyParametersGetPtr(cosmology,cosmology_ID)',
              fc_arg.fort_reassignment, "FC: GetPtr call in reassignment")
    assert_equal(fc_arg.fort_function_class, 'cosmologyParameters', "FC: fort_function_class set")


# ============================================================================
# Emitter tests
# ============================================================================

def test_ctypes_arg_types():
    print("\n=== Testing ctypes_arg_types() ===")

    # Mixed types: value double, reference double, char_p
    raw_args = [
        {'name': 'x', 'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)']},
        {'name': 'y', 'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(out)']},
        {'name': 'label', 'intrinsic': 'character', 'type': None, 'attributes': ['intent(in)']},
    ]
    al = assign_c_types(raw_args, {})
    al = assign_c_attributes(al)
    result = ctypes_arg_types(al)
    assert_equal(result[0], 'c_double',             "ctypes_arg_types: intent(in) double → c_double")
    assert_equal(result[1], 'POINTER(c_double)',     "ctypes_arg_types: intent(out) double → POINTER")
    assert_equal(result[2], 'c_char_p',              "ctypes_arg_types: char → c_char_p (no pointer wrap)")


def test_fortran_arg_list():
    print("\n=== Testing fortran_arg_list() ===")

    lib_fc = {'cosmologyParameters': {'module': 'Cosmology_Parameters'}}
    raw_args = [
        {'name': 'x',         'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)']},
        {'name': 'cosmology', 'intrinsic': 'class', 'type': 'cosmologyParametersClass', 'attributes': ['intent(in)']},
    ]
    al = assign_c_types(raw_args, lib_fc)
    # x, cosmology, cosmology_ID — all fort_is_present=True
    result = fortran_arg_list(al)
    assert_equal(result, ['x', 'cosmology', 'cosmology_ID'], "fortran_arg_list: includes _ID companion")


def test_fortran_declarations():
    print("\n=== Testing fortran_declarations() ===")

    lib_fc = {'cosmologyParameters': {'module': 'Cosmology_Parameters'}}
    fc = {'name': 'cosmologyParameters', 'module': 'Cosmology_Parameters'}
    raw_args = [
        {'name': 'x',         'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)']},
        {'name': 'cosmology', 'intrinsic': 'class', 'type': 'cosmologyParametersClass', 'attributes': ['intent(in)']},
    ]
    al = run_pipeline(raw_args, lib_fc=lib_fc, func_class=fc)
    decls = fortran_declarations(al)

    assert_in('real(c_double), value :: x\n',         decls, "declarations: double precision line")
    assert_in('type(c_ptr), value :: cosmology\n',    decls, "declarations: functionClass c_ptr line")
    assert_in('integer(c_int), value :: cosmology_ID\n', decls, "declarations: _ID line")
    # interface block
    assert_in('function cosmologyParametersGetPtr(ptr_,classID)', decls, "declarations: GetPtr interface block")
    assert_in(f'import c_int, c_ptr, cosmologyParametersClass',   decls, "declarations: interface import")


def test_fortran_reassignments():
    print("\n=== Testing fortran_reassignments() ===")

    fc = {'name': 'foo', 'module': 'Foo_Module'}
    raw_args = [
        {'name': 'flag',  'intrinsic': 'logical',          'type': None,    'attributes': ['intent(in)']},
        {'name': 'value', 'intrinsic': 'double precision',  'type': None,    'attributes': ['intent(in)']},
    ]
    al = run_pipeline(raw_args, func_class=fc)
    reassign = fortran_reassignments(al)

    assert_in('flag_=logical(flag)\n', reassign, "fortran_reassignments: logical cast present")
    assert_not_in('value', reassign,              "fortran_reassignments: no reassignment for double")


def test_fortran_module_uses():
    print("\n=== Testing fortran_module_uses() ===")

    fc = {'name': 'foo', 'module': 'Foo_Module'}
    raw_args = [
        {'name': 'node',  'intrinsic': 'type', 'type': 'treeNode',       'attributes': ['intent(inout)']},
        {'name': 'desc',  'intrinsic': 'type', 'type': 'varying_string', 'attributes': ['intent(in)']},
    ]
    al = run_pipeline(raw_args, func_class=fc)
    uses = fortran_module_uses(al)

    assert_in('use :: Galacticus_Nodes, only : treeNode', uses, "module_uses: Galacticus_Nodes")
    assert_in('use :: String_Handling',                   uses, "module_uses: String_Handling")


def test_iso_c_binding_import():
    print("\n=== Testing iso_c_binding_import() ===")

    fc = {'name': 'foo', 'module': 'Foo_Module'}
    raw_args = [
        {'name': 'x',    'intrinsic': 'double precision', 'type': None,       'attributes': ['intent(in)']},
        {'name': 'node', 'intrinsic': 'type',             'type': 'treeNode', 'attributes': ['intent(inout)']},
    ]
    al = run_pipeline(raw_args, func_class=fc)
    stmt = iso_c_binding_import(al, 'c_ptr', 'c_loc')

    assert_in('c_double',    stmt, "iso_c_binding: c_double from double precision")
    assert_in('c_f_pointer', stmt, "iso_c_binding: c_f_pointer from treeNode")
    assert_in('c_loc',       stmt, "iso_c_binding: extra symbol passed in")
    assert_in('c_ptr',       stmt, "iso_c_binding: extra c_ptr passed in")


def test_fortran_call_code_no_optional():
    print("\n=== Testing fortran_call_code(): no optional args ===")

    fc = {'name': 'foo', 'module': 'Foo_Module'}
    raw_args = [
        {'name': 'x', 'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)']},
        {'name': 'y', 'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)']},
    ]
    al = run_pipeline(raw_args, func_class=fc)
    call = fortran_call_code(al, 'result=myFunc(', ')\n', '&')

    assert_in('x=x',     call, "no-optional call: x argument")
    assert_in('y=y',     call, "no-optional call: y argument")
    assert_not_in('if (', call, "no-optional call: no branching")


def test_fortran_call_code_with_optional():
    print("\n=== Testing fortran_call_code(): with optional args ===")

    fc = {'name': 'foo', 'module': 'Foo_Module'}
    raw_args = [
        {'name': 'x',    'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)']},
        {'name': 'time', 'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)', 'optional']},
        {'name': 'z',    'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)', 'optional']},
    ]
    al = run_pipeline(raw_args, func_class=fc)
    call = fortran_call_code(al, 'result=myFunc(', ')\n', '&')

    # 2 optional args → 4 branches (2^2 subsets from powerset).
    # Each pv appears once per branch condition (both as present(x) and .not.present(x));
    # since .not.present(x) contains present(x) as a substring, the total count is 4.
    assert_equal(call.count('present(time)'), 4, "optional call: present(time) in each of 4 branch conditions")
    assert_equal(call.count('present(z)'),    4, "optional call: present(z) in each of 4 branch conditions")
    assert_in('.not.present(time)',   call, "optional call: negated present check present")
    assert_in('else\n',               call, "optional call: unconditional else fallback")
    assert_in('end if\n',             call, "optional call: end if terminator")


def test_python_arg_list():
    print("\n=== Testing python_arg_list() ===")

    fc = {'name': 'foo', 'module': 'Foo_Module'}

    # Constructor-style: first arg is not 'self' → 'self' prepended
    raw_args = [
        {'name': 'x', 'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)']},
        {'name': 'y', 'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)']},
    ]
    al = run_pipeline(raw_args, func_class=fc)
    args = python_arg_list(al)
    assert_equal(args[0], 'self',   "constructor: 'self' prepended")
    assert_equal(args[1], 'x',     "constructor: x follows")
    assert_equal(args[2], 'y',     "constructor: y follows")

    # With optional arg: gets '=None' suffix
    raw_args = [
        {'name': 'x',    'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)']},
        {'name': 'time', 'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)', 'optional']},
    ]
    al = run_pipeline(raw_args, func_class=fc)
    args = python_arg_list(al)
    assert_equal(args[-1], 'time=None', "optional arg gets '=None' suffix")
    assert_not_in('x=None', args,        "non-optional arg has no '=None' suffix")

    # functionClass _ID arg: py_is_present=False → excluded from Python list
    lib_fc = {'cosmologyParameters': {'module': 'Cosmology_Parameters'}}
    raw_args = [
        {'name': 'cosmology', 'intrinsic': 'class',
         'type': 'cosmologyParametersClass', 'attributes': ['intent(in)']},
    ]
    al = run_pipeline(raw_args, lib_fc=lib_fc, func_class=fc)
    args = python_arg_list(al)
    names = [a.split('=')[0] for a in args]
    assert_not_in('cosmology_ID', names, "functionClass _ID hidden from Python arg list")


def test_python_reassignments():
    print("\n=== Testing python_reassignments() ===")

    lib_fc = {'cosmologyParameters': {'module': 'Cosmology_Parameters'}}
    fc = {'name': 'cosmologyParameters', 'module': 'Cosmology_Parameters'}

    # Non-optional FC: no reassignment block
    raw_args = [
        {'name': 'cosmology', 'intrinsic': 'class',
         'type': 'cosmologyParametersClass', 'attributes': ['intent(in)']},
    ]
    al = run_pipeline(raw_args, lib_fc=lib_fc, func_class=fc)
    assert_equal(python_reassignments(al), '', "non-optional FC: empty python_reassignments")

    # Optional FC: reassignment block present
    raw_args = [
        {'name': 'cosmo', 'intrinsic': 'class',
         'type': 'cosmologyParametersClass', 'attributes': ['intent(in)', 'optional']},
    ]
    al = run_pipeline(raw_args, lib_fc=lib_fc, func_class=fc)
    block = python_reassignments(al)
    assert_in('if cosmo:', block, "optional FC: if-branch in python_reassignments")
    assert_in('cosmo_glcObj', block, "optional FC: glcObj assignment")
    assert_in('None',         block, "optional FC: None in absent branch")


def test_python_call_code_no_optional():
    print("\n=== Testing python_call_code(): no optional args ===")

    fc = {'name': 'foo', 'module': 'Foo_Module'}
    raw_args = [
        {'name': 'x', 'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)']},
        {'name': 'y', 'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)']},
    ]
    al = run_pipeline(raw_args, func_class=fc)
    call = python_call_code(al, 'c_lib.fooL')

    assert_in('c_lib.fooL(x,y)', call, "no-optional python call: correct invocation")
    assert_not_in('if ', call,         "no-optional python call: no branching")


def test_python_call_code_with_optional():
    print("\n=== Testing python_call_code(): with optional args ===")

    fc = {'name': 'foo', 'module': 'Foo_Module'}
    raw_args = [
        {'name': 'x',    'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)']},
        {'name': 'time', 'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)', 'optional']},
        {'name': 'z',    'intrinsic': 'double precision', 'type': None, 'attributes': ['intent(in)', 'optional']},
    ]
    al = run_pipeline(raw_args, func_class=fc)
    call = python_call_code(al, 'c_lib.fooL')

    # 2 optional pvs → 2^2 = 4 if/elif branches.
    # Each pv appears in every branch condition (positive or negative),
    # so 'is not None' appears: 0+1+1+2 = 4 times across all conditions,
    # and 'is None' appears: 2+1+1+0 = 4 times.
    assert_equal(call.count('is not None'), 4, "python call: 'is not None' appears across 4 branches")
    assert_equal(call.count('is None'),     4, "python call: 'is None' appears across 4 branches")
    assert_in('None',        call,             "python call: absent optionals passed as None")
    assert_in('c_double(z)', call,             "python call: explicit ctype cast after first optional")


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 70)
    print("LibraryInterfaces Pipeline + Emitters Regression Tests")
    print("=" * 70)

    test_argspec_from_raw()
    test_assign_c_types_scalar_intrinsics()
    test_assign_c_types_type_variants()
    test_assign_c_types_optional()
    test_assign_c_types_function_class()
    test_assign_c_types_self_exclusion()
    test_assign_c_types_optional_function_class()
    test_assign_c_attributes()
    test_build_python_reassignments_non_optional()
    test_build_python_reassignments_optional()
    test_fort_reassign_logical()
    test_fort_reassign_character()
    test_fort_reassign_varying_string()
    test_fort_reassign_enumeration()
    test_fort_reassign_tree_node()
    test_fort_reassign_tree_node_optional()
    test_fort_reassign_function_class()
    test_ctypes_arg_types()
    test_fortran_arg_list()
    test_fortran_declarations()
    test_fortran_reassignments()
    test_fortran_module_uses()
    test_iso_c_binding_import()
    test_fortran_call_code_no_optional()
    test_fortran_call_code_with_optional()
    test_python_arg_list()
    test_python_reassignments()
    test_python_call_code_no_optional()
    test_python_call_code_with_optional()

    print("\n" + "=" * 70)
    print(f"Results: {PASS_COUNT} passed, {FAIL_COUNT} failed")
    print("=" * 70)

    if FAIL_COUNT > 0:
        print("FAIL: Some tests failed")
        sys.exit(1)
    else:
        print("PASS: All tests passed")
        sys.exit(0)


if __name__ == '__main__':
    main()
