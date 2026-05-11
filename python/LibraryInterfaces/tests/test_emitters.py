"""Tests for `LibraryInterfaces.Emitters`.

Ten emitters that turn enriched ArgSpec lists into the Fortran/Python code
strings consumed by the build of `libgalacticus.so` and its `galacticus.py`
wrapper.  These tests focus on the cross-cutting behaviours -- ctypes
wrapping, optional-arg branching, ISO_C_Binding symbol collection, and
functionClass interface-block emission -- rather than every byte of
generated text.
"""

from LibraryInterfaces.ArgSpec  import ArgSpec
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


# ---------------------------------------------------------------------------
# ctypes_arg_types
# ---------------------------------------------------------------------------

def test_ctypes_arg_types_wraps_pointer_args_in_POINTER():
    args = [
        ArgSpec(name='x', ctype='c_double', ctype_pointer=False),
        ArgSpec(name='y', ctype='c_int',    ctype_pointer=True),
    ]
    assert ctypes_arg_types(args) == ['c_double', 'POINTER(c_int)']


def test_ctypes_arg_types_falls_back_to_c_int_when_ctype_unset():
    """An ArgSpec with no ctype gets `c_int` — defensive default for any
    misconfigured arg slipping through the pipeline."""
    args = [ArgSpec(name='x')]
    assert ctypes_arg_types(args) == ['c_int']


# ---------------------------------------------------------------------------
# fortran_arg_list
# ---------------------------------------------------------------------------

def test_fortran_arg_list_filters_by_fort_is_present():
    """Only args with fort_is_present=True appear in the bind(C) signature."""
    args = [
        ArgSpec(name='a',                                    fort_is_present=True),
        ArgSpec(name='hidden_internal_param',                fort_is_present=False),
        ArgSpec(name='b',                                    fort_is_present=True),
    ]
    assert fortran_arg_list(args) == ['a', 'b']


# ---------------------------------------------------------------------------
# fortran_declarations
# ---------------------------------------------------------------------------

def test_fortran_declarations_emits_per_arg_decl_lines():
    args = [
        ArgSpec(name='x', fort_type='real(c_double)',
                fort_attributes=['intent(in)', 'value']),
    ]
    out = fortran_declarations(args)
    assert '  real(c_double), intent(in), value :: x' in out


def test_fortran_declarations_appends_extra_fort_declarations():
    """Extra per-arg fort_declarations (e.g. local pointer vars introduced by
    Pipeline reassignments) are emitted right after the main declaration."""
    args = [
        ArgSpec(name='n', fort_type='type(c_ptr)',
                fort_attributes=['intent(in)', 'value'],
                fort_declarations='type(treeNode), pointer :: n_\n'),
    ]
    out = fortran_declarations(args)
    assert 'type(c_ptr), intent(in), value :: n' in out
    assert 'type(treeNode), pointer :: n_' in out


def test_fortran_declarations_skips_signature_decl_for_absent_args():
    """Args with fort_is_present=False (e.g. null-filled overrides) must
    not get the default `integer(c_int)` signature declaration — only
    their `fort_declarations` block, which carries the local
    null-pointer the inner constructor receives.  Emitting both would
    declare the same name twice and trip gfortran with `VALUE attribute
    conflicts with POINTER attribute`."""
    args = [
        ArgSpec(name='initializationFunction',
                fort_is_present=False,
                fort_declarations='procedure(someInitializor), pointer :: '
                                  'initializationFunction => null()\n'),
    ]
    out = fortran_declarations(args)
    assert 'integer(c_int)' not in out
    assert ':: initializationFunction\n' not in out
    assert ('procedure(someInitializor), pointer :: '
            'initializationFunction => null()') in out


def test_fortran_declarations_emits_GetPtr_interface_block():
    """Each distinct fort_function_class produces ONE interface block
    declaring the GetPtr function used to recover a typed pointer."""
    args = [
        ArgSpec(name='engine', fort_type='type(c_ptr)',
                fort_attributes=['intent(in)', 'value'],
                fort_function_class='cooling'),
    ]
    out = fortran_declarations(args)
    assert 'function coolingGetPtr(ptr_,classID)' in out
    assert 'class(coolingClass), pointer :: coolingGetPtr' in out


def test_fortran_declarations_dedups_GetPtr_blocks():
    """Two args citing the same fort_function_class produce ONE
    interface block, not two.  (`coolingGetPtr` appears twice per block:
    once in `function coolingGetPtr(...)` and once in `end function
    coolingGetPtr` — so a single block has count==2.)"""
    args = [
        ArgSpec(name='a', fort_function_class='cooling'),
        ArgSpec(name='b', fort_function_class='cooling'),
    ]
    out = fortran_declarations(args)
    # `function coolingGetPtr(ptr_,classID)` appears exactly once per block.
    assert out.count('function coolingGetPtr(ptr_,classID)') == 1


# ---------------------------------------------------------------------------
# fortran_reassignments  /  python_reassignments — concatenation
# ---------------------------------------------------------------------------

def test_fortran_reassignments_concatenates_in_order():
    args = [
        ArgSpec(name='a', fort_reassignment='line_a\n'),
        ArgSpec(name='b'),  # empty fort_reassignment
        ArgSpec(name='c', fort_reassignment='line_c\n'),
    ]
    assert fortran_reassignments(args) == 'line_a\nline_c\n'


def test_python_reassignments_concatenates_in_order():
    args = [
        ArgSpec(name='a', py_reassignment='py_a\n'),
        ArgSpec(name='b', py_reassignment='py_b\n'),
    ]
    assert python_reassignments(args) == 'py_a\npy_b\n'


# ---------------------------------------------------------------------------
# fortran_module_uses — accumulate + dedupe + sort
# ---------------------------------------------------------------------------

def test_fortran_module_uses_emits_one_use_per_module_with_sorted_only_list():
    args = [
        ArgSpec(name='x', fort_modules={'String_Handling': {'String_C_to_Fortran': 1}}),
        ArgSpec(name='y', fort_modules={'String_Handling': {'char': 1}}),
        ArgSpec(name='z', fort_modules={'Galacticus_Nodes': {'treeNode': 1}}),
    ]
    out = fortran_module_uses(args)
    # Modules emitted in sorted order; symbols within each `only:` sorted.
    assert out == (
        '  use :: Galacticus_Nodes, only : treeNode\n'
        '  use :: String_Handling, only : String_C_to_Fortran, char\n'
    )


def test_fortran_module_uses_no_args_returns_empty():
    assert fortran_module_uses([]) == ''


# ---------------------------------------------------------------------------
# iso_c_binding_import
# ---------------------------------------------------------------------------

def test_iso_c_binding_import_extracts_kind_symbols_from_fort_types():
    args = [
        ArgSpec(name='a', fort_type='real(c_double)'),
        ArgSpec(name='b', fort_type='integer(c_int)'),
    ]
    out = iso_c_binding_import(args)
    # Both symbols present, sorted alphabetically.
    assert out == '  use, intrinsic :: ISO_C_Binding, only : c_double, c_int\n'


def test_iso_c_binding_import_includes_extra_symbols():
    """`iso_c_binding_import(args, 'c_loc')` adds c_loc to the only-list."""
    args = [ArgSpec(name='a', fort_type='real(c_double)')]
    out = iso_c_binding_import(args, 'c_loc')
    assert 'c_double' in out
    assert 'c_loc' in out


def test_iso_c_binding_import_includes_per_arg_iso_c_symbols():
    """Symbols pushed into fort_iso_c_symbols (e.g. `c_f_pointer` from a
    treeNode reassignment) are picked up automatically."""
    args = [ArgSpec(name='n', fort_type='type(c_ptr)',
                    fort_iso_c_symbols=['c_f_pointer'])]
    out = iso_c_binding_import(args)
    assert 'c_f_pointer' in out
    assert 'c_ptr'       in out


# ---------------------------------------------------------------------------
# python_arg_list
# ---------------------------------------------------------------------------

def test_python_arg_list_constructor_prepends_self():
    """A constructor (first arg not named 'self') gets an explicit 'self'
    prepended to the Python argument list."""
    args = [
        ArgSpec(name='count', py_is_present=True),
        ArgSpec(name='value', py_is_present=True),
    ]
    assert python_arg_list(args) == ['self', 'count', 'value']


def test_python_arg_list_method_does_not_double_prepend_self():
    """A method (first arg already named 'self') doesn't get a duplicate
    'self' prepended."""
    args = [
        ArgSpec(name='self', py_is_present=True),
        ArgSpec(name='x',    py_is_present=True),
    ]
    assert python_arg_list(args) == ['self', 'x']


def test_python_arg_list_marks_optional_with_default_None():
    """Once an optional arg appears, IT and every subsequent included arg
    get `=None`."""
    args = [
        ArgSpec(name='self', py_is_present=True),
        ArgSpec(name='req',  py_is_present=True),
        ArgSpec(name='opt',  py_is_present=True, is_optional=True),
        ArgSpec(name='also', py_is_present=True),
    ]
    out = python_arg_list(args)
    assert out == ['self', 'req', 'opt=None', 'also=None']


def test_python_arg_list_filters_by_py_is_present():
    """Args with py_is_present=False (e.g. _ID companions) are excluded."""
    args = [
        ArgSpec(name='self', py_is_present=True),
        ArgSpec(name='hidden_id', py_is_present=False),
        ArgSpec(name='x', py_is_present=True),
    ]
    assert python_arg_list(args) == ['self', 'x']


# ---------------------------------------------------------------------------
# fortran_call_code — optional-argument branching
# ---------------------------------------------------------------------------

def test_fortran_call_code_no_optionals_emits_unconditional_call():
    args = [
        ArgSpec(name='a', galacticus_is_present=True),
        ArgSpec(name='b', galacticus_is_present=True),
    ]
    out = fortran_call_code(args, pre_arguments='self%method(',
                            post_arguments=')\n', continuation='&')
    assert 'self%method(' in out
    assert 'a=a' in out
    assert 'b=b' in out
    # No present()-conditional in the no-optionals branch.
    assert 'present(' not in out


def test_fortran_call_code_one_optional_emits_2_branches_plus_else():
    """N optional args → 2^N if/else-if branches plus an unconditional
    else fallback."""
    args = [
        ArgSpec(name='req', galacticus_is_present=True),
        ArgSpec(name='opt', galacticus_is_present=True, is_optional=True),
    ]
    out = fortran_call_code(args, pre_arguments='self%method(',
                            post_arguments=')\n', continuation='&')
    assert out.count('if (') >= 2  # one 'if', possibly 'else if' or 'else'
    assert '.not.present(opt)' in out
    assert 'present(opt)' in out
    assert 'else' in out
    assert 'end if' in out


def test_fortran_call_code_skips_galacticus_is_present_false():
    """Args with galacticus_is_present=False (e.g. self via dispatch) are
    omitted from the call."""
    args = [
        ArgSpec(name='self', galacticus_is_present=False),
        ArgSpec(name='x',    galacticus_is_present=True),
    ]
    out = fortran_call_code(args, pre_arguments='self%method(',
                            post_arguments=')\n', continuation='&')
    assert 'x=x' in out
    assert 'self=self' not in out


# ---------------------------------------------------------------------------
# python_call_code — optional-argument branching
# ---------------------------------------------------------------------------

def test_python_call_code_no_optionals_unconditional():
    args = [
        ArgSpec(name='a', fort_is_present=True),
        ArgSpec(name='b', fort_is_present=True),
    ]
    out = python_call_code(args, call='c_lib.method')
    # Single line with no if-branching.
    assert 'c_lib.method(a,b)' in out
    assert 'is None' not in out


def test_python_call_code_optional_emits_branches_with_ctype_wrapping():
    """Once the first optional is hit, every subsequent included arg gets
    explicit ctype wrapping (ctypes can't infer types when some args may be
    None).  Absent optionals are passed as `None`."""
    args = [
        ArgSpec(name='req', fort_is_present=True),
        ArgSpec(name='opt', fort_is_present=True, is_optional=True,
                ctype='c_int'),
    ]
    out = python_call_code(args, call='c_lib.method')
    # Branches over `opt is not None` / `opt is None`.
    assert 'opt is not None' in out
    assert 'opt is None'     in out
    # When opt is present: ctype-wrapped.
    assert 'c_int(opt)' in out
    # When opt is absent: passed as None.
    assert ',None' in out
