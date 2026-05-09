"""Regression test for `_parse_interface_arguments` in
`Galacticus.Build.SourceTree.Process.EventHooks`.

Bug: the original implementation split each declaration line on the FIRST
`::`.  But Fortran attribute lists may contain `dimension(:)` (or
`dimension(:,:)`), and the colon inside the parens looks just like the
`::` separator to a naive split — so the parser locked onto the inner
colon and failed to find the real `::`, dropping the argument list.

Concretely, the abstract-interface body for `events.hooks.p.F90` contains

    real(c_double), dimension(:), intent(in   ) :: time
    class(treeNode), intent(inout) :: node

The old splitter dropped `time`, leaving the rendered abstract interface
missing one of its arguments and yielding a Fortran syntax error.  The
fix splits on the LAST `::` instead.
"""

from Galacticus.Build.SourceTree.Process.EventHooks import _parse_interface_arguments


def test_parse_interface_arguments_with_dimension_colon():
    text = (
        "real(c_double), dimension(:), intent(in   ) :: time\n"
        "class(treeNode), intent(inout)              :: node\n"
    )
    args = _parse_interface_arguments(text)
    assert args == ['time', 'node']


def test_parse_interface_arguments_handles_multi_dim():
    text = "real, dimension(:,:), intent(in) :: matrix\n"
    args = _parse_interface_arguments(text)
    assert args == ['matrix']


def test_parse_interface_arguments_multiple_per_line():
    """A single declaration with multiple variables exposes them all."""
    text = "integer, intent(in) :: i, j, k\n"
    args = _parse_interface_arguments(text)
    assert args == ['i', 'j', 'k']


def test_parse_interface_arguments_skips_non_declarations():
    text = (
        "! a comment line\n"
        "use, intrinsic :: iso_c_binding\n"
        "integer, intent(in) :: i\n"
    )
    args = _parse_interface_arguments(text)
    assert args == ['i']
