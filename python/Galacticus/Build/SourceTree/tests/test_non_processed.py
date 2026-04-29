# Regression test for `is_non_processed_type`.
#
# Bug: directives like `<methods>` or `<optionalArgument>Task` are emitted by
# code-generating Process hooks AFTER `nonProcessed` has run.  Without the
# `is_non_processed_type` exemption they tripped `post_process_directives`'s
# "directive 'foo' was not processed" check at the end of the pipeline.
#
# The exemption list is fixed (Methods, Workaround, Include, FunctionGlobal,
# Component, RadiusSolverPlausibility, InterTreePositionInsert, Expiry,
# Scoping) plus any directive whose type ends in `Task`.

from Galacticus.Build.SourceTree.Process.NonProcessed import is_non_processed_type


def test_known_non_processed_types():
    for t in ('methods', 'workaround', 'include', 'functionGlobal',
              'component', 'radiusSolverPlausibility',
              'interTreePositionInsert', 'expiry', 'scoping'):
        assert is_non_processed_type(t), t


def test_task_suffix_is_non_processed():
    """Any `*Task` directive is exempt — covers `optionalArgumentTask`,
    `nodeOperatorTask`, etc., which would otherwise need to be enumerated."""
    assert is_non_processed_type('optionalArgumentTask')
    assert is_non_processed_type('nodeOperatorTask')


def test_unknown_type_is_processed():
    """Random directives are NOT exempt — the post-process check catches
    them."""
    assert not is_non_processed_type('functionClass')
    assert not is_non_processed_type('eventHook')
    assert not is_non_processed_type('')
    assert not is_non_processed_type(None)
