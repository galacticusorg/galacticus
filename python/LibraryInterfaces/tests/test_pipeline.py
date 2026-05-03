# Tests for `LibraryInterfaces.Pipeline`.
#
# Four pipeline stages enrich `ArgSpec` lists for the cross-language Galacticus
# library wrapper.  Each test exercises one type-mapping branch end-to-end
# rather than every internal helper -- the type mappings are the contract that
# the rest of the build relies on.

import pytest

from LibraryInterfaces.ArgSpec  import ArgSpec
from LibraryInterfaces.Pipeline import (
    assign_c_types,
    assign_c_attributes,
    build_python_reassignments,
    build_fortran_reassignments,
)


# ---------------------------------------------------------------------------
# assign_c_types — Fortran intrinsic / class type → ctype mapping
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("intrinsic,expected_ctype,expected_fort", [
    ('double precision', 'c_double', 'real(c_double)'),
    ('integer',          'c_int',    'integer(c_int)'),
    ('logical',          'c_bool',   'logical(c_bool)'),
    ('character',        'c_char_p', 'character(c_char)'),
])
def test_assign_c_types_intrinsic_scalars(intrinsic, expected_ctype, expected_fort):
    raw = [{'name': 'x', 'intrinsic': intrinsic, 'attributes': ['intent(in)']}]
    out = assign_c_types(raw, lib_function_classes={})
    assert len(out) == 1
    assert out[0].ctype     == expected_ctype
    assert out[0].fort_type == expected_fort


@pytest.mark.parametrize("kind,expected_ctype,expected_fort", [
    ('c_long',   'c_long',   'integer(c_long)'),
    ('c_size_t', 'c_size_t', 'integer(c_size_t)'),
])
def test_assign_c_types_integer_kinds(kind, expected_ctype, expected_fort):
    """integer(c_long) and integer(c_size_t) pass through with matching
    ctypes wrappers rather than being silently demoted to c_int."""
    raw = [{'name': 'n', 'intrinsic': 'integer', 'type': kind,
            'attributes': ['intent(in)']}]
    out = assign_c_types(raw, lib_function_classes={})
    assert out[0].ctype     == expected_ctype
    assert out[0].fort_type == expected_fort


def test_assign_c_types_varying_string_maps_to_c_char_p():
    """type(varying_string) is treated like character — passes as c_char_p."""
    raw = [{'name': 'name', 'intrinsic': 'type', 'type': 'varying_string',
            'attributes': ['intent(in)']}]
    out = assign_c_types(raw, lib_function_classes={})
    assert out[0].ctype     == 'c_char_p'
    assert out[0].fort_type == 'character(c_char)'


def test_assign_c_types_enumeration_maps_to_c_int():
    """type(enumerationXxxType) maps to c_int (case-insensitive)."""
    raw = [{'name': 'kind', 'intrinsic': 'type',
            'type': 'enumerationStarFormationModeType', 'attributes': []}]
    out = assign_c_types(raw, lib_function_classes={})
    assert out[0].ctype     == 'c_int'
    assert out[0].fort_type == 'integer(c_int)'


def test_assign_c_types_other_derived_type_maps_to_c_void_p():
    raw = [{'name': 'n', 'intrinsic': 'type', 'type': 'treeNode',
            'attributes': ['intent(in)']}]
    out = assign_c_types(raw, lib_function_classes={})
    assert out[0].ctype     == 'c_void_p'
    assert out[0].fort_type == 'type(c_ptr)'


def test_assign_c_types_function_class_inserts_id_companion():
    """class(fooClass) where 'foo' is a registered functionClass must
    insert a `name_ID` companion arg of integer type immediately after."""
    raw = [{'name': 'engine', 'intrinsic': 'class', 'type': 'fooClass',
            'attributes': ['intent(in)']}]
    out = assign_c_types(raw, lib_function_classes={'foo': {'module': 'Foo'}})

    # Should produce TWO ArgSpecs: the original arg, then the _ID companion.
    assert [a.name for a in out] == ['engine', 'engine_ID']
    assert out[0].is_function_class is True
    assert out[1].name              == 'engine_ID'
    assert out[1].ctype             == 'c_int'
    assert out[1].py_is_present     is False
    assert out[1].galacticus_is_present is False


def test_assign_c_types_self_function_class_skips_galacticus_pass_through():
    """A `self` argument of class(XClass) is dispatched via method binding,
    not passed to the Galacticus constructor — galacticus_is_present should
    be False."""
    raw = [{'name': 'self', 'intrinsic': 'class', 'type': 'fooClass',
            'attributes': ['intent(in)']}]
    out = assign_c_types(raw, lib_function_classes={'foo': {'module': 'Foo'}})
    assert out[0].name                  == 'self'
    assert out[0].galacticus_is_present is False


def test_assign_c_types_optional_function_class_marks_id_companion_optional():
    """When the parent functionClass arg is optional, the inserted _ID
    companion is also marked optional and gets the parent's name as its
    py_present trigger."""
    raw = [{'name': 'engine', 'intrinsic': 'class', 'type': 'fooClass',
            'attributes': ['intent(in)', 'optional']}]
    out = assign_c_types(raw, lib_function_classes={'foo': {'module': 'Foo'}})
    parent, companion = out
    assert parent.is_optional      is True
    assert companion.is_optional   is True
    assert 'optional' in companion.attributes
    assert companion.py_present    == 'engine'


def test_assign_c_types_class_not_in_lib_classes_is_plain_pointer():
    """class(X) where X isn't a registered functionClass is just c_void_p,
    no _ID companion."""
    raw = [{'name': 'engine', 'intrinsic': 'class', 'type': 'fooClass',
            'attributes': ['intent(in)']}]
    out = assign_c_types(raw, lib_function_classes={})  # foo not registered
    assert len(out) == 1
    assert out[0].is_function_class is False
    assert out[0].ctype             == 'c_void_p'


def test_assign_c_types_input_order_preserved():
    """Multiple args without _ID expansion stay in input order."""
    raw = [
        {'name': 'a', 'intrinsic': 'integer',          'attributes': ['intent(in)']},
        {'name': 'b', 'intrinsic': 'double precision', 'attributes': ['intent(in)']},
        {'name': 'c', 'intrinsic': 'logical',          'attributes': ['intent(in)']},
    ]
    out = assign_c_types(raw, lib_function_classes={})
    assert [a.name for a in out] == ['a', 'b', 'c']


# ---------------------------------------------------------------------------
# assign_c_attributes — pass-by determination
# ---------------------------------------------------------------------------

def test_assign_c_attributes_intent_in_scalar_is_pass_by_value():
    """A non-optional scalar with intent(in) passes by value."""
    args = [ArgSpec(name='x', intrinsic='integer', ctype='c_int',
                    attributes=['intent(in)'])]
    assign_c_attributes(args)
    assert args[0].pass_by         == 'value'
    assert 'value' in args[0].fort_attributes
    assert args[0].ctype_pointer   is False


def test_assign_c_attributes_intent_out_passes_by_reference():
    """intent(out) cannot be value-passed; must be POINTER on Python side."""
    args = [ArgSpec(name='x', intrinsic='integer', ctype='c_int',
                    attributes=['intent(out)'])]
    assign_c_attributes(args)
    assert args[0].pass_by       == 'reference'
    assert args[0].ctype_pointer is True


def test_assign_c_attributes_optional_uses_reference_passing():
    """Optional args can't be value-passed (Fortran limitation)."""
    args = [ArgSpec(name='x', intrinsic='integer', ctype='c_int',
                    attributes=['intent(in)', 'optional'], is_optional=True)]
    assign_c_attributes(args)
    assert args[0].pass_by         == 'reference'
    assert 'optional' in args[0].fort_attributes


def test_assign_c_attributes_dimension_attr_propagates_and_forces_reference():
    """A `dimension(...)` attribute is preserved on fort_attributes and
    forces reference passing (arrays can't be value-passed)."""
    args = [ArgSpec(name='vec', intrinsic='double precision', ctype='c_double',
                    attributes=['intent(in)', 'dimension(:)'])]
    assign_c_attributes(args)
    assert 'dimension(:)' in args[0].fort_attributes
    assert args[0].pass_by == 'reference'


def test_assign_c_attributes_c_char_p_gets_dimension_star():
    """c_char_p strings get an implicit `dimension(*)` attribute on the
    Fortran side and ctype_pointer is False (no POINTER wrapping)."""
    args = [ArgSpec(name='s', intrinsic='character', ctype='c_char_p',
                    attributes=['intent(in)'])]
    assign_c_attributes(args)
    assert 'dimension(*)' in args[0].fort_attributes
    assert args[0].ctype_pointer is False


# ---------------------------------------------------------------------------
# build_python_reassignments
# ---------------------------------------------------------------------------

def test_python_reassignments_non_optional_function_class():
    """Non-optional functionClass args use attribute access on the Python
    object: `arg._glcObj` / `arg._classID`."""
    parent = ArgSpec(name='engine', is_function_class=True, type_spec='fooClass')
    companion = ArgSpec(name='engine_ID')
    out = build_python_reassignments([parent, companion])
    assert out[0].py_pass_as == 'engine._glcObj'
    assert out[1].py_pass_as == 'engine._classID'
    # Non-optional case: no py_reassignment block.
    assert out[0].py_reassignment == ''


def test_python_reassignments_optional_function_class_synthesises_block():
    """Optional functionClass args get a Python if/else block that
    extracts _glcObj/_classID or sets them to None."""
    parent = ArgSpec(name='engine', is_function_class=True, is_optional=True,
                     type_spec='fooClass')
    companion = ArgSpec(name='engine_ID')
    out = build_python_reassignments([parent, companion])
    assert out[0].py_pass_as == 'engine_glcObj'
    assert out[1].py_pass_as == 'engine_classID'
    assert 'if engine:' in out[0].py_reassignment
    assert 'engine._glcObj' in out[0].py_reassignment
    assert 'else:' in out[0].py_reassignment


# ---------------------------------------------------------------------------
# build_fortran_reassignments — the heart of cross-language type conversion
# ---------------------------------------------------------------------------

def test_fortran_reassignments_logical_recasts_via_logical_call():
    """c_bool -> logical via `logical()` cast."""
    args = [ArgSpec(name='flag', intrinsic='logical')]
    out = build_fortran_reassignments(
        args, func_class={}, implementation=None,
        extensions={}, module_uses_impls={},
    )
    assert 'flag_=logical(flag)' in out[0].fort_reassignment
    assert out[0].fort_pass_as == 'flag_'


def test_fortran_reassignments_character_routes_through_String_C_to_Fortran():
    args = [ArgSpec(name='s', intrinsic='character')]
    out = build_fortran_reassignments(
        args, func_class={}, implementation=None,
        extensions={}, module_uses_impls={},
    )
    assert 'String_C_to_Fortran(s)' in out[0].fort_pass_as
    # Module declarations are recorded for the Fortran-emitter to use.
    assert 'String_C_to_Fortran' in out[0].fort_modules['String_Handling']
    assert 'char' in out[0].fort_modules['ISO_Varying_String']


def test_fortran_reassignments_varying_string_routes_through_String_C_to_Fortran():
    args = [ArgSpec(name='label', intrinsic='type', type_spec='varying_string')]
    out = build_fortran_reassignments(
        args, func_class={}, implementation=None,
        extensions={}, module_uses_impls={},
    )
    assert out[0].fort_pass_as == 'String_C_to_Fortran(label)'
    assert 'String_C_to_Fortran' in out[0].fort_modules['String_Handling']


def test_fortran_reassignments_treeNode_uses_c_f_pointer():
    args = [ArgSpec(name='n', intrinsic='type', type_spec='treeNode')]
    out = build_fortran_reassignments(
        args, func_class={}, implementation=None,
        extensions={}, module_uses_impls={},
    )
    assert 'c_f_pointer(n,n_)' in out[0].fort_reassignment
    assert out[0].fort_pass_as == 'n_'
    assert 'c_f_pointer' in out[0].fort_iso_c_symbols
    assert 'treeNode' in out[0].fort_modules['Galacticus_Nodes']


def test_fortran_reassignments_mergerTree_imported_from_Galacticus_Nodes():
    """type(mergerTree) must come from Galacticus_Nodes — without the
    explicit special case the fall-back branch would import it from the
    functionClass's own module (e.g. Merger_Trees_Build_Mass_Resolution),
    where it isn't defined."""
    args = [ArgSpec(name='tree', intrinsic='type', type_spec='mergerTree')]
    out = build_fortran_reassignments(
        args,
        func_class={'module': 'Merger_Trees_Build_Mass_Resolution'},
        implementation=None, extensions={}, module_uses_impls={},
    )
    assert 'mergerTree' in out[0].fort_modules['Galacticus_Nodes']
    assert 'Merger_Trees_Build_Mass_Resolution' not in out[0].fort_modules
    assert 'c_f_pointer(tree,tree_)' in out[0].fort_reassignment


def test_fortran_reassignments_multiCounter_imported_from_Multi_Counters():
    """type(multiCounter) must come from the Multi_Counters module rather
    than the functionClass's own module."""
    args = [ArgSpec(name='counter', intrinsic='type', type_spec='multiCounter')]
    out = build_fortran_reassignments(
        args,
        func_class={'module': 'Some_Other_Module'},
        implementation=None, extensions={}, module_uses_impls={},
    )
    assert 'multiCounter' in out[0].fort_modules['Multi_Counters']
    assert 'Some_Other_Module' not in out[0].fort_modules
    assert 'c_f_pointer(counter,counter_)' in out[0].fort_reassignment


def test_fortran_reassignments_optional_treeNode_emits_present_branch():
    """An optional treeNode argument gets an explicit `if (present(...))`
    block in fort_reassignment so the absent case sets the pointer to null."""
    args = [ArgSpec(name='n', intrinsic='type', type_spec='treeNode',
                    is_optional=True)]
    out = build_fortran_reassignments(
        args, func_class={}, implementation=None,
        extensions={}, module_uses_impls={},
    )
    assert 'if (present(n))' in out[0].fort_reassignment
    assert 'else' in out[0].fort_reassignment
    assert '=> null()' in out[0].fort_reassignment


def test_fortran_reassignments_function_class_uses_GetPtr_function():
    """isFunctionClass args use FooGetPtr(ptr, ID) to recover a typed pointer."""
    parent = ArgSpec(name='engine', is_function_class=True, type_spec='fooClass')
    out = build_fortran_reassignments(
        [parent],
        func_class={}, implementation=None,
        extensions={}, module_uses_impls={},
        lib_function_classes={'foo': {'module': 'Foo'}},
    )
    assert 'fooGetPtr(engine,engine_ID)' in out[0].fort_reassignment
    assert out[0].fort_pass_as == 'engine_'
    assert out[0].fort_function_class == 'foo'
    assert 'fooClass' in out[0].fort_modules['Foo']


def test_fortran_reassignments_arbitrary_derived_type_uses_c_f_pointer():
    """Unknown derived types get c_f_pointer + a module declaration from
    func_class.module."""
    args = [ArgSpec(name='cosmo', intrinsic='type', type_spec='cosmologyParameters')]
    out = build_fortran_reassignments(
        args,
        func_class={'module': 'Cosmology_Parameters'},
        implementation=None, extensions={}, module_uses_impls={},
    )
    assert 'call c_f_pointer(cosmo,cosmo_)' in out[0].fort_reassignment
    assert 'cosmologyParameters' in out[0].fort_modules['Cosmology_Parameters']
