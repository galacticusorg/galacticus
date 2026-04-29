# Regression test for `_process_allocatable_intrinsic` in
# `Galacticus.Build.SourceTree.Process.StateStorable`.
#
# Bug: `_process_allocatable_intrinsic` set `fragments['rank_seen']` to the
# allocatable's rank, but its emitted body uses Fortran whole-array I/O
# (`write (stateFile) self%name`) and never opens a `do iN=...` loop.
# The caller bumped `rank_maximum` from `rank_seen` and therefore declared
# `integer(c_size_t) :: i1, i2, i3` — `i2`/`i3` then triggered
# `[-Wunused-variable]` warnings whenever the type's only multi-rank
# member was an allocatable intrinsic.

from Galacticus.Build.SourceTree.Process.StateStorable import (
    _process_allocatable_intrinsic,
)


def _decl(intrinsic, attributes, variables):
    return {
        'intrinsic':  intrinsic,
        'type':       None,
        'attributes': list(attributes),
        'variables':  list(variables),
    }


def test_allocatable_intrinsic_rank3_does_not_request_indices():
    """A rank-3 allocatable intrinsic should report `rank_seen == 0` so
    the caller doesn't declare unused index variables."""
    decl = _decl('real',
                 ['allocatable', 'dimension(:,:,:)'],
                 ['matrix'])
    frag = _process_allocatable_intrinsic(decl, exclude=[])
    assert frag['rank_seen'] == 0, frag
    # The actual body still works for rank-3 — it uses whole-array I/O —
    # so `stored_shape_required` and `was_allocated_required` are set and
    # the body references the variable.
    assert frag['stored_shape_required']  is True
    assert frag['was_allocated_required'] is True
    assert 'self%matrix' in frag['output']
    assert 'self%matrix' in frag['input']
    # And critically, NO `do iN=…` loop got emitted.
    assert 'do i1' not in frag['output']
    assert 'do i1' not in frag['input']


def test_allocatable_intrinsic_rank1_does_not_request_indices():
    decl = _decl('real', ['allocatable', 'dimension(:)'], ['vec'])
    frag = _process_allocatable_intrinsic(decl, exclude=[])
    assert frag['rank_seen'] == 0


def test_allocatable_scalar_intrinsic_unchanged():
    """A rank-0 allocatable scalar already had `rank_seen == 0` — make sure
    the fix didn't perturb that path."""
    decl = _decl('real', ['allocatable'], ['x'])
    frag = _process_allocatable_intrinsic(decl, exclude=[])
    assert frag['rank_seen'] == 0
    assert frag['stored_shape_required'] is False
