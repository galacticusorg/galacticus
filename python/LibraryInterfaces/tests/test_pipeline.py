# Tests for `LibraryInterfaces.Pipeline`.
#
# Four pipeline stages enrich `ArgSpec` lists for the cross-language Galacticus
# library wrapper.  Each test exercises one type-mapping branch end-to-end
# rather than every internal helper -- the type mappings are the contract that
# the rest of the build relies on.

import pytest

from LibraryInterfaces.ArgSpec  import ArgSpec
from LibraryInterfaces.Emitters import python_call_code
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


def test_assign_c_types_optional_omp_lock_kind_is_dropped():
    """`integer(omp_lock_kind)` has a platform-dependent size — INTEGER(4)
    on Linux GCC, INTEGER(8) on macOS GCC — so the wrapper can't pick a
    single C-interop kind that matches the inner method on both
    platforms.  Optional ones are silently dropped from the wrapper
    (the inner method's optional default kicks in); the surrounding
    method's other args still flow through unchanged."""
    raw = [
        {'name': 'tree',  'intrinsic': 'type', 'type': 'mergerTree',
         'attributes': ['intent(inout)', 'target']},
        {'name': 'initializationLock', 'intrinsic': 'integer',
         'type': 'omp_lock_kind',
         'attributes': ['intent(inout)', 'optional']},
    ]
    out = assign_c_types(raw, lib_function_classes={})
    names = [a.name for a in out]
    assert 'initializationLock' not in names
    assert 'tree' in names


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


def test_assign_c_types_fixed_array_no_count_companion():
    """A double precision, dimension(3) argument is is_array=True with
    array_size=3, but does NOT get a count companion (the size is in the
    dimension spec)."""
    raw = [{'name': 'point', 'intrinsic': 'double precision', 'type': None,
            'attributes': ['intent(in)', 'dimension(3)']}]
    out = assign_c_types(raw, lib_function_classes={})
    assert [a.name for a in out] == ['point']
    assert out[0].is_array   is True
    assert out[0].array_size == 3
    assert out[0].ctype      == 'c_double'


def test_python_reassignments_fixed_array_validates_size():
    """A fixed-size is_array arg's py_reassignment converts via numpy AND
    raises ValueError if the input's size doesn't match the declared
    length — a wrong-size input would otherwise silently corrupt
    Fortran-side memory."""
    arr = ArgSpec(name='point', intrinsic='double precision',
                  ctype='c_double', is_array=True, array_size=3)
    out = build_python_reassignments([arr])
    assert 'np.ascontiguousarray(point' in out[0].py_reassignment
    assert 'point.size != 3'            in out[0].py_reassignment
    assert 'raise ValueError'           in out[0].py_reassignment
    assert out[0].py_pass_as == 'point.ctypes.data_as(POINTER(c_double))'


def test_fortran_reassignments_fixed_array_passed_directly():
    """A fixed-size is_array arg's fort_pass_as is unset (defaults to the
    arg name) — no slicing is needed because the bind(c) declaration is
    already explicit-shape `dimension(N)` and matches the inner method."""
    arr = ArgSpec(name='point', intrinsic='double precision', ctype='c_double',
                  fort_type='real(c_double)', is_array=True, array_size=3)
    out = build_fortran_reassignments(
        [arr], func_class={}, implementation=None,
        extensions={}, module_uses_impls={},
    )
    assert out[0].fort_pass_as == ''  # no slicing — name passes through


def test_assign_c_types_deferred_array_inserts_count_companion():
    """A double precision, dimension(:) argument gets a hidden c_size_t
    count companion immediately after it; both are flagged so the rest of
    the pipeline can recognise them."""
    raw = [{'name': 'times', 'intrinsic': 'double precision', 'type': None,
            'attributes': ['intent(in)', 'dimension(:)']}]
    out = assign_c_types(raw, lib_function_classes={})
    assert [a.name for a in out] == ['times', 'times_count']
    assert out[0].is_array is True
    assert out[0].ctype     == 'c_double'
    assert out[1].name      == 'times_count'
    assert out[1].ctype     == 'c_size_t'
    assert out[1].fort_type == 'integer(c_size_t)'
    # Hidden from the user-facing Python signature; never passed to the
    # inner Galacticus call.
    assert out[1].py_is_present         is False
    assert out[1].galacticus_is_present is False


@pytest.mark.parametrize("intrinsic,kind,expected_ctype", [
    ('double precision', None      , 'c_double'),
    ('integer'         , None      , 'c_int'   ),
    ('integer'         , 'c_long'  , 'c_long'  ),
    ('integer'         , 'c_size_t', 'c_size_t'),
])
def test_assign_c_types_deferred_array_ctype_matches_scalar_kind(
        intrinsic, kind, expected_ctype):
    """The array's element ctype follows the same kind-mapping as a scalar
    of the same intrinsic."""
    raw = [{'name': 'a', 'intrinsic': intrinsic, 'type': kind,
            'attributes': ['intent(in)', 'dimension(:)']}]
    out = assign_c_types(raw, lib_function_classes={})
    assert out[0].is_array is True
    assert out[0].ctype     == expected_ctype


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
    """A `dimension(:)` attribute (assumed-shape) is rewritten to
    `dimension(*)` (assumed-size) for bind(c) compatibility, propagated to
    fort_attributes, and forces reference passing (arrays can't be
    value-passed)."""
    args = [ArgSpec(name='vec', intrinsic='double precision', ctype='c_double',
                    attributes=['intent(in)', 'dimension(:)'])]
    assign_c_attributes(args)
    assert 'dimension(*)' in args[0].fort_attributes
    assert 'dimension(:)' not in args[0].fort_attributes
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


def test_fortran_reassignments_array_arg_passed_as_slice():
    """An is_array arg's fort_pass_as is `name(1:name_count)` so the inner
    Galacticus method receives a proper rank-1 array section bounded by
    the count companion's runtime value."""
    arr = ArgSpec(name='times', intrinsic='double precision', ctype='c_double',
                  fort_type='real(c_double)', is_array=True)
    out = build_fortran_reassignments(
        [arr], func_class={}, implementation=None,
        extensions={}, module_uses_impls={},
    )
    assert out[0].fort_pass_as == 'times(1:times_count)'


def test_python_reassignments_array_arg_converts_via_numpy():
    """An is_array arg's py_reassignment converts the user's input through
    np.ascontiguousarray of the right dtype, and py_pass_as hands ctypes
    the array's data pointer.  The count companion's py_pass_as reads
    .size off the (now numpy-array) input."""
    arr = ArgSpec(name='times', intrinsic='double precision',
                  ctype='c_double', is_array=True)
    cnt = ArgSpec(name='times_count', intrinsic='integer',
                  type_spec='c_size_t', ctype='c_size_t',
                  fort_is_present=True, py_is_present=False,
                  galacticus_is_present=False)
    out = build_python_reassignments([arr, cnt])
    assert 'np.ascontiguousarray(times' in out[0].py_reassignment
    assert 'dtype=np.float64'           in out[0].py_reassignment
    assert out[0].py_pass_as == 'times.ctypes.data_as(POINTER(c_double))'
    assert out[1].py_pass_as == 'c_size_t(times.size)'


def test_python_reassignments_optional_array_gates_conversion_on_None():
    """An optional 1D array gates the np.ascontiguousarray conversion on
    `is not None`.  Without the gate, a default value of None propagates
    into `np.ascontiguousarray(None, dtype=...)` which yields a 0-D
    scalar in modern numpy and breaks downstream `.size` / `.ctypes`
    access.  The pass expression and count are likewise None-aware so
    the absent-arg branch of python_call_code passes a NULL pointer
    and a zero count cleanly."""
    arr = ArgSpec(name='valueTarget', intrinsic='double precision',
                  ctype='c_double', is_array=True, is_optional=True)
    cnt = ArgSpec(name='valueTarget_count', intrinsic='integer',
                  type_spec='c_size_t', ctype='c_size_t',
                  fort_is_present=True, py_is_present=False,
                  galacticus_is_present=False)
    out = build_python_reassignments([arr, cnt])
    assert 'if valueTarget is not None:' in out[0].py_reassignment
    assert 'np.ascontiguousarray(valueTarget'  in out[0].py_reassignment
    assert out[0].py_pass_as == \
        'valueTarget.ctypes.data_as(POINTER(c_double)) if valueTarget is not None else None'
    # Count companion: NOT wrapped in c_size_t(...) — python_call_code
    # adds that once the arg crosses the first-optional boundary, and a
    # second wrap (c_size_t(c_size_t(...))) is rejected by ctypes.
    assert out[1].py_pass_as == \
        'valueTarget.size if valueTarget is not None else 0'


def test_python_reassignments_optional_2d_array_gates_conversion_on_None():
    """Same gating for the 2D case.  Without it,
    `np.asfortranarray(np.asarray(None, dtype=...))` produces a 1-D
    array (asarray gives 0-D, asfortranarray promotes) and trips the
    ndim==2 check below the conversion."""
    arr = ArgSpec(name='covarianceTarget', intrinsic='double precision',
                  ctype='c_double', is_array=True, array_rank=2,
                  is_optional=True)
    c1  = ArgSpec(name='covarianceTarget_count_1', intrinsic='integer',
                  type_spec='c_size_t', ctype='c_size_t',
                  fort_is_present=True, py_is_present=False,
                  galacticus_is_present=False)
    c2  = ArgSpec(name='covarianceTarget_count_2', intrinsic='integer',
                  type_spec='c_size_t', ctype='c_size_t',
                  fort_is_present=True, py_is_present=False,
                  galacticus_is_present=False)
    out = build_python_reassignments([arr, c1, c2])
    assert 'if covarianceTarget is not None:'        in out[0].py_reassignment
    assert 'covarianceTarget is not None else None'  in out[0].py_pass_as
    # Same c_size_t double-wrap concern as the 1D test above — leave
    # the wrap to python_call_code; emit a plain Python int / 0 here.
    assert 'covarianceTarget is not None else 0'     in out[1].py_pass_as
    assert 'covarianceTarget is not None else 0'     in out[2].py_pass_as
    assert 'c_size_t' not in out[1].py_pass_as
    assert 'c_size_t' not in out[2].py_pass_as


def test_fortran_reassignments_optional_2d_array_gated_on_present():
    """Optional 2D array's fort_reassignment is wrapped in
    `if (present(...)) then ... end if` so the bind(c) wrapper doesn't
    slice a NULL pointer when the user passed None on the Python side."""
    arr = ArgSpec(name='arr', intrinsic='double precision',
                  ctype='c_double', fort_type='real(c_double)',
                  is_array=True, array_rank=2, is_optional=True)
    out = build_fortran_reassignments(
        [arr], func_class={}, implementation=None,
        extensions={}, module_uses_impls={},
    )
    assert 'if (present(arr)) then'   in out[0].fort_reassignment
    assert 'allocate(arr_F_'          in out[0].fort_reassignment
    assert 'reshape(arr(1:'           in out[0].fort_reassignment
    assert 'end if'                   in out[0].fort_reassignment


def test_python_reassignments_c_char_p_encodes_str_to_bytes():
    """ctypes' c_char_p expects bytes (it implements `char *` as a
    NUL-terminated byte string).  Python users naturally pass `str`,
    and ctypes won't auto-encode — it raises 'bytes or integer address
    expected instead of str instance'.  build_python_reassignments
    must emit an isinstance-gated `.encode()` so str inputs work and
    None / bytes pass through unchanged."""
    s = ArgSpec(name='xAxisLabel', intrinsic='type', type_spec='varying_string',
                ctype='c_char_p', is_optional=True)
    out = build_python_reassignments([s])
    assert 'isinstance(xAxisLabel, str)'        in out[0].py_reassignment
    assert 'xAxisLabel = xAxisLabel.encode('    in out[0].py_reassignment


def test_python_reassignments_c_char_p_uses_safe_name_for_keyword_args():
    """A `class`/`lambda`/... arg name has been escaped via python_safe_name
    in the Python signature, so the encode reassignment must reference the
    escaped identifier — otherwise the generated code raises NameError."""
    s = ArgSpec(name='class', intrinsic='character', ctype='c_char_p')
    out = build_python_reassignments([s])
    assert 'isinstance(class_, str)' in out[0].py_reassignment
    assert 'class_ = class_.encode(' in out[0].py_reassignment


def test_python_call_code_array_arg_not_wrapped_in_ctype():
    """Array args are passed as POINTER(...) expressions; wrapping them
    in `{ctype}(pa)` (e.g. `c_double(<pointer>)`) — which is what
    python_call_code does for scalar optional args — would crash at
    runtime.  python_call_code must recognise is_array and pass the
    expression directly."""
    arr = ArgSpec(name='valueTarget', intrinsic='double precision',
                  ctype='c_double', is_array=True, is_optional=True,
                  py_pass_as='valueTarget.ctypes.data_as(POINTER(c_double)) if valueTarget is not None else None',
                  fort_is_present=True)
    out = python_call_code([arr], 'self._glcObj = c_lib.fooL')
    assert 'c_double(valueTarget' not in out
    assert 'valueTarget.ctypes.data_as(POINTER(c_double))' in out


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


# ---------------------------------------------------------------------------
# 2D deferred-shape numeric arrays — `dimension(:,:)`
# ---------------------------------------------------------------------------

def test_assign_c_types_2d_array_inserts_two_count_companions():
    """A `double precision, dimension(:,:)` argument expands into the
    arg followed by two c_size_t count companions (one per axis).
    array_rank=2 marks it for the 2D-specific Python and Fortran paths."""
    raw = [{'name': 'covariance', 'intrinsic': 'double precision', 'type': None,
            'attributes': ['intent(in)', 'dimension(:,:)']}]
    out = assign_c_types(raw, lib_function_classes={})
    assert [a.name for a in out] == [
        'covariance', 'covariance_count_1', 'covariance_count_2',
    ]
    assert out[0].is_array   is True
    assert out[0].array_rank == 2
    assert out[0].ctype      == 'c_double'
    assert all(c.ctype == 'c_size_t' for c in out[1:])
    assert all(c.py_is_present is False and c.galacticus_is_present is False
               for c in out[1:])


def test_assign_c_attributes_2d_array_collapsed_to_dimension_star():
    """Both `dimension(:)` and `dimension(:,:)` collapse to `dimension(*)`
    on the bind(c) signature — assumed-shape isn't interoperable, and the
    wrapper rebuilds the rank-2 view inside via `reshape`."""
    raw = [{'name': 'm', 'intrinsic': 'double precision', 'type': None,
            'attributes': ['intent(in)', 'dimension(:,:)']}]
    out = assign_c_types(raw, lib_function_classes={})
    out = assign_c_attributes(out)
    assert 'dimension(*)' in out[0].fort_attributes
    assert 'dimension(:,:)' not in out[0].fort_attributes


def test_python_reassignments_2d_array_uses_asfortranarray_with_shape_counts():
    """The Python wrapper converts inputs through np.asfortranarray (so
    the underlying buffer matches Fortran's column-major layout) and
    fills both count companions from `.shape[0]` / `.shape[1]`."""
    arr = ArgSpec(name='m', intrinsic='double precision',
                  ctype='c_double', is_array=True, array_rank=2)
    c1  = ArgSpec(name='m_count_1', intrinsic='integer',
                  type_spec='c_size_t', ctype='c_size_t',
                  fort_is_present=True, py_is_present=False,
                  galacticus_is_present=False)
    c2  = ArgSpec(name='m_count_2', intrinsic='integer',
                  type_spec='c_size_t', ctype='c_size_t',
                  fort_is_present=True, py_is_present=False,
                  galacticus_is_present=False)
    out = build_python_reassignments([arr, c1, c2])
    assert 'np.asfortranarray' in out[0].py_reassignment
    assert 'm.ndim != 2'       in out[0].py_reassignment
    assert out[1].py_pass_as == 'c_size_t(m.shape[0])'
    assert out[2].py_pass_as == 'c_size_t(m.shape[1])'
    assert out[0].py_pass_as == 'm.ctypes.data_as(POINTER(c_double))'


def test_fortran_reassignments_2d_array_reshapes_into_local():
    """The wrapper allocates a `dimension(:,:)` local and reshapes the
    flat bind(c) buffer into it before passing to the inner Galacticus
    method; the inner call receives the rank-2 local, not the flat
    rank-1 dummy."""
    arr = ArgSpec(name='m', intrinsic='double precision',
                  ctype='c_double', fort_type='real(c_double)',
                  is_array=True, array_rank=2)
    out = build_fortran_reassignments(
        [arr], func_class={}, implementation=None,
        extensions={}, module_uses_impls={},
    )
    decls    = out[0].fort_declarations
    reassign = out[0].fort_reassignment
    assert 'real(c_double), dimension(:,:), allocatable :: m_F_' in decls
    assert 'allocate(m_F_(m_count_1, m_count_2))' in reassign
    assert 'reshape(m(1:m_count_1*m_count_2), [m_count_1, m_count_2])' in reassign
    assert out[0].fort_pass_as == 'm_F_'


# ---------------------------------------------------------------------------
# 1D fixed-length character arrays — `character(len=N), dimension(:)`
# ---------------------------------------------------------------------------

def test_assign_c_types_char_array_inserts_count_companion():
    """`character(len=2), dimension(:)` is treated as a 1D array: a hidden
    c_size_t count companion is inserted, the element ctype is c_char,
    and char_len is captured for the downstream emitters."""
    raw = [{'name': 'elements', 'intrinsic': 'character', 'type': 'len=2',
            'attributes': ['intent(in)', 'dimension(:)']}]
    out = assign_c_types(raw, lib_function_classes={})
    assert [a.name for a in out] == ['elements', 'elements_count']
    assert out[0].is_array  is True
    assert out[0].ctype     == 'c_char'
    assert out[0].fort_type == 'character(c_char)'
    assert out[0].char_len  == 2
    assert out[1].ctype     == 'c_size_t'


def test_assign_c_types_scalar_character_unaffected():
    """Plain `character(len=*)` / `character(c_char)` scalars still go
    through the null-terminated-string path (ctype=c_char_p), not the new
    array branch."""
    raw = [{'name': 'name', 'intrinsic': 'character', 'type': 'len=*',
            'attributes': ['intent(in)']}]
    out = assign_c_types(raw, lib_function_classes={})
    assert len(out) == 1
    assert out[0].is_array is False
    assert out[0].char_len == 0
    assert out[0].ctype    == 'c_char_p'


def test_python_reassignments_char_array_packs_to_fixed_byte_buffer():
    """The Python wrapper space-pads + clips each input string to N
    characters, encodes ASCII, and lays them out contiguously as `S{N}`
    so .ctypes.data_as gives a valid count*N byte buffer to bind(c)."""
    arr = ArgSpec(name='elements', intrinsic='character',
                  ctype='c_char', char_len=2, is_array=True)
    cnt = ArgSpec(name='elements_count', intrinsic='integer',
                  type_spec='c_size_t', ctype='c_size_t',
                  fort_is_present=True, py_is_present=False,
                  galacticus_is_present=False)
    out = build_python_reassignments([arr, cnt])
    assert "ljust(2)" in out[0].py_reassignment
    assert "[:2]"     in out[0].py_reassignment
    assert "dtype='S2'" in out[0].py_reassignment
    assert out[0].py_pass_as == "elements.ctypes.data_as(POINTER(c_char))"
    # Count companion reads .size off the now-numpy buffer.
    assert out[1].py_pass_as == 'c_size_t(elements.size)'


def test_fortran_reassignments_char_array_repacks_into_local():
    """The wrapper allocates a `character(len=N), dimension(:)` local,
    copies the flat bind(c) byte buffer in element by element, and
    passes that local to the inner Galacticus call."""
    arr = ArgSpec(name='elements', intrinsic='character',
                  ctype='c_char', fort_type='character(c_char)',
                  char_len=2, is_array=True)
    out = build_fortran_reassignments(
        [arr], func_class={}, implementation=None,
        extensions={}, module_uses_impls={},
    )
    decls = out[0].fort_declarations
    reassign = out[0].fort_reassignment
    assert 'character(len=2), dimension(:), allocatable :: elements_F_' in decls
    assert 'allocate(elements_F_(elements_count))' in reassign
    assert 'elements_F_(elements_glcI_)(elements_glcJ_:elements_glcJ_)' in reassign
    assert 'elements((elements_glcI_-1)*2 + elements_glcJ_)' in reassign
    assert out[0].fort_pass_as == 'elements_F_'


# ---------------------------------------------------------------------------
# 1D varying_string arrays — `type(varying_string), dimension(:)`
# ---------------------------------------------------------------------------

def test_assign_c_types_vstring_array_inserts_count_and_charlen_companions():
    """`type(varying_string), dimension(:)` ships as a flat byte buffer
    plus two c_size_t companions: count (number of elements) and
    charLen (per-element ASCII byte length, runtime-decided)."""
    raw = [{'name': 'names', 'intrinsic': 'type', 'type': 'varying_string',
            'attributes': ['intent(in)', 'dimension(:)']}]
    out = assign_c_types(raw, lib_function_classes={})
    assert [a.name for a in out] == ['names', 'names_count', 'names_charLen']
    assert out[0].is_array             is True
    assert out[0].varying_string_array is True
    assert out[0].ctype                == 'c_char'
    assert out[0].fort_type            == 'character(c_char)'
    assert out[1].ctype                == 'c_size_t'
    assert out[2].ctype                == 'c_size_t'


def test_python_reassignments_vstring_array_picks_runtime_max_length():
    """The Python wrapper ASCII-encodes each input element, picks the
    max length at runtime, and lays the result out as a contiguous
    `S{N}` numpy buffer.  The count companion reads `.size`, the
    char-length companion reads `.itemsize` so Fortran sees the
    runtime-chosen N."""
    arr = ArgSpec(name='names', intrinsic='type', type_spec='varying_string',
                  ctype='c_char', is_array=True, varying_string_array=True)
    cnt = ArgSpec(name='names_count', intrinsic='integer',
                  type_spec='c_size_t', ctype='c_size_t',
                  fort_is_present=True, py_is_present=False,
                  galacticus_is_present=False)
    cl  = ArgSpec(name='names_charLen', intrinsic='integer',
                  type_spec='c_size_t', ctype='c_size_t',
                  fort_is_present=True, py_is_present=False,
                  galacticus_is_present=False)
    out = build_python_reassignments([arr, cnt, cl])
    assert "names_glcStrs_" in out[0].py_reassignment
    assert "names_glcLen_"  in out[0].py_reassignment
    assert "dtype=f'S{names_glcLen_}'" in out[0].py_reassignment
    assert out[0].py_pass_as == "names.ctypes.data_as(POINTER(c_char))"
    assert out[1].py_pass_as == 'c_size_t(names.size)'
    assert out[2].py_pass_as == 'c_size_t(names.itemsize)'


def test_fortran_reassignments_vstring_array_repacks_into_varying_string_local():
    """The wrapper allocates a `type(varying_string), dimension(:)`
    local, builds a per-element scratch character buffer of the
    runtime-decided length, copies bytes in, then trims and assigns to
    each varying_string element.  ISO_Varying_String must supply both
    the type and the assignment(=) overload."""
    arr = ArgSpec(name='names', intrinsic='type', type_spec='varying_string',
                  ctype='c_char', fort_type='character(c_char)',
                  is_array=True, varying_string_array=True)
    out = build_fortran_reassignments(
        [arr], func_class={}, implementation=None,
        extensions={}, module_uses_impls={},
    )
    decls    = out[0].fort_declarations
    reassign = out[0].fort_reassignment
    assert 'type(varying_string), dimension(:), allocatable :: names_F_' in decls
    assert 'character(len=:), allocatable :: names_buf_'                  in decls
    assert 'allocate(names_F_(names_count))'                              in reassign
    assert 'allocate(character(len=int(names_charLen)) :: names_buf_)'    in reassign
    assert 'names_F_(names_glcI_) = trim(names_buf_)'                     in reassign
    assert out[0].fort_pass_as                                            == 'names_F_'
    assert 'varying_string'  in out[0].fort_modules.get('ISO_Varying_String', {})
    assert 'assignment(=)'   in out[0].fort_modules.get('ISO_Varying_String', {})


# ---------------------------------------------------------------------------
# 1D polymorphic-list arrays — `type(<class>List), dimension(:)`
# ---------------------------------------------------------------------------

def test_assign_c_types_list_array_inserts_count_and_id_companions():
    """`type(modelParameterList), dimension(:)` (the Galacticus idiom for
    "array of class(modelParameterClass)") expands into a c_void_p
    pointer buffer + a count companion + a parallel c_int IDs buffer."""
    raw = [{'name': 'parms', 'intrinsic': 'type',
            'type': 'modelParameterList',
            'attributes': ['intent(in)', 'dimension(:)']}]
    out = assign_c_types(
        raw, lib_function_classes={'modelParameter': {'module': 'Model_Parameters'}})
    assert [a.name for a in out] == ['parms', 'parms_count', 'parms_IDs']
    assert out[0].is_array                   is True
    assert out[0].polymorphic_list_array     is True
    assert out[0].polymorphic_list_class     == 'modelParameter'
    assert out[0].polymorphic_list_type      == 'modelParameterList'
    assert out[0].polymorphic_list_component == 'modelParameter_'
    assert out[0].ctype                      == 'c_void_p'
    assert out[0].fort_type                  == 'type(c_ptr)'
    assert out[1].ctype                      == 'c_size_t'   # count
    assert out[2].ctype                      == 'c_int'      # IDs (array)
    # IDs companion is a deferred-shape array so the bind(c) emitter
    # rewrites it to dimension(*) (see assign_c_attributes).
    assert 'dimension(:)' in out[2].attributes


def test_assign_c_types_list_array_unregistered_stem_falls_back_to_void_p():
    """`type(fooList), dimension(:)` where `foo` is not a registered
    functionClass takes the generic-derived-type path (ctype=c_void_p)
    rather than the polymorphic-list-array path — _unsupported_arg
    will reject the surrounding constructor in libraryInterfaces.py."""
    raw = [{'name': 'xs', 'intrinsic': 'type', 'type': 'fooList',
            'attributes': ['intent(in)', 'dimension(:)']}]
    out = assign_c_types(raw, lib_function_classes={})
    assert out[0].polymorphic_list_array is False


def test_assign_c_types_list_array_unregistered_wrapper_falls_back_to_void_p():
    """`type(<class>List), dimension(:)` where `<class>` IS a registered
    functionClass but the wrapper type isn't in `_SHARED_TYPE_MODULES`
    must NOT take the polymorphic-list-array path — locally-defined
    wrappers can carry extra members (e.g.
    `virialDensityContrastList` in tasks.halo_mass_function.F90 has a
    `label` field alongside the polymorphic pointer) which our
    rebuild loop would silently drop.  Falling through to the
    generic-derived-type path means `_unsupported_arg` rejects the
    surrounding constructor cleanly."""
    raw = [{'name': 'xs', 'intrinsic': 'type',
            'type': 'virialDensityContrastList',
            'attributes': ['intent(in)', 'dimension(:)']}]
    out = assign_c_types(
        raw, lib_function_classes={'virialDensityContrast': {'module': 'X'}})
    assert out[0].polymorphic_list_array is False
    # Single-element output: no companion args were inserted.
    assert len(out) == 1


def test_python_reassignments_list_array_builds_ptr_and_id_buffers():
    """The Python wrapper builds parallel ctypes arrays of object
    pointers + class IDs, one per element, off the `_glcObj` /
    `_classID` fields the wrapper attaches to every functionClass
    instance."""
    arr = ArgSpec(name='parms', intrinsic='type',
                  type_spec='modelParameterList',
                  ctype='c_void_p', is_array=True,
                  polymorphic_list_array=True,
                  polymorphic_list_class='modelParameter',
                  polymorphic_list_type='modelParameterList',
                  polymorphic_list_component='modelParameter_')
    cnt = ArgSpec(name='parms_count', intrinsic='integer',
                  type_spec='c_size_t', ctype='c_size_t',
                  fort_is_present=True, py_is_present=False,
                  galacticus_is_present=False)
    ids = ArgSpec(name='parms_IDs', intrinsic='integer',
                  ctype='c_int', fort_is_present=True,
                  py_is_present=False, galacticus_is_present=False,
                  attributes=['intent(in)', 'dimension(:)'])
    out = build_python_reassignments([arr, cnt, ids])
    reassign = out[0].py_reassignment
    assert 'parms_glcN_'    in reassign
    assert '_glcObj'        in reassign
    assert '_classID'       in reassign
    assert '(c_void_p * parms_glcN_)' in reassign
    assert '(c_int    * parms_glcN_)' in reassign
    assert out[0].py_pass_as == 'parms_glcPtrs_'
    assert out[1].py_pass_as == 'c_size_t(parms_glcN_)'
    assert out[2].py_pass_as == 'parms_glcIDs_'


def test_fortran_reassignments_list_array_loops_via_get_ptr():
    """The Fortran wrapper allocates a `type(<class>List), dimension(:)`
    local of the right size, then loops over each element and uses the
    registered class's `<class>GetPtr` helper to recover the polymorphic
    pointer that the wrapper struct's `<class>_` component should hold.
    The `_SHARED_TYPE_MODULES` table supplies the wrapper struct's home
    module for the `use ::` line."""
    arr = ArgSpec(name='parms', intrinsic='type',
                  type_spec='modelParameterList',
                  ctype='c_void_p', fort_type='type(c_ptr)',
                  is_array=True, polymorphic_list_array=True,
                  polymorphic_list_class='modelParameter',
                  polymorphic_list_type='modelParameterList',
                  polymorphic_list_component='modelParameter_')
    out = build_fortran_reassignments(
        [arr], func_class={}, implementation=None,
        extensions={}, module_uses_impls={},
        lib_function_classes={'modelParameter': {'module': 'Model_Parameters'}},
    )
    decls    = out[0].fort_declarations
    reassign = out[0].fort_reassignment
    assert 'type(modelParameterList), dimension(:), allocatable :: parms_F_' in decls
    assert 'integer(c_size_t) :: parms_glcI_'                                in decls
    assert 'allocate(parms_F_(parms_count))'                                 in reassign
    assert ('parms_F_(parms_glcI_)%modelParameter_ =>'
            ' modelParameterGetPtr(parms(parms_glcI_), parms_IDs(parms_glcI_))') in reassign
    assert out[0].fort_pass_as        == 'parms_F_'
    assert out[0].fort_function_class == 'modelParameter'
    # `_SHARED_TYPE_MODULES` supplies modelParameterList's home module.
    # The Class type ALSO has to be imported into the bind(c) host
    # scope — otherwise the GetPtr interface block emitted by
    # `fortran_declarations` (from `fort_function_class`) hits
    # "Cannot IMPORT 'modelparameterclass' from host scoping unit"
    # because the host only `use`s the List type, not the Class.
    assert 'modelParameterList'  in out[0].fort_modules.get('Model_Parameters', {})
    assert 'modelParameterClass' in out[0].fort_modules.get('Model_Parameters', {})
