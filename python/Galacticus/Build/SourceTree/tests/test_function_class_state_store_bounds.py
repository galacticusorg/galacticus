"""Regression test for the bounds of allocatable arrays in the functionClass
state store (`Galacticus.Build.SourceTree.Process.FunctionClass.StateStore`).

Bug (#1786): only the shape of an allocatable member array was written to the
state file, and on restore the array was reallocated with that shape - so with
lower bounds of 1. The tabulations of the Farahi excursion set solvers are
declared with zero lower bounds (`dimension(0:n,0:m,:)`), so a restored solver
read every rate from the neighboring tabulated point, and merger trees built
after a state retrieve differed from those of the original run.
"""

from Galacticus.Build.SourceTree.Process.FunctionClass.StateStore import (
    bounds_expression, bounds_specification,
    generate_allocatable_state_store_code,
)


def test_bounds_record_and_specification_agree():
    assert bounds_expression('x', 2) == (
        'lbound(self%x,kind=c_size_t),ubound(self%x,kind=c_size_t)')
    # Lower bounds come first in the record, then upper bounds.
    assert bounds_specification(3) == (
        'storedBounds(1):storedBounds(4),'
        'storedBounds(2):storedBounds(5),'
        'storedBounds(3):storedBounds(6)')


def test_allocatable_member_is_restored_with_its_bounds():
    state_store, state_stores = {}, {}
    generate_allocatable_state_store_code(
        state_store, state_stores, 'rate', 3, 'self%rate', 'self%rate')
    assert 'shape(self%rate' not in state_store['outputCode']
    assert (' write (stateFile) lbound(self%rate,kind=c_size_t),'
            'ubound(self%rate,kind=c_size_t)') in state_store['outputCode']
    assert 'allocate(storedBounds(6))' in state_store['inputCode']
    assert ('allocate(self%rate(storedBounds(1):storedBounds(4),'
            'storedBounds(2):storedBounds(5),'
            'storedBounds(3):storedBounds(6)))') in state_store['inputCode']
