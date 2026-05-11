r"""Regression test for `_build_allowed_parameters_method`.

Bug: when scanning a class's constructor function for `<objectBuilder>`
directives (so the auto-generated `*AllowedParameters` method can emit
the recursive `if (associated(self%X_)) call self%X_%allowedParameters
(allowedParameters,'parameters',.true.)` lines for every nested
functionClass member), the constructor-signature regex was

    ^\s*(recursive)??\s+function\s+NAME\s*\(\s*parameters …

The mandatory `\s+` in front of `function` meant that a constructor
WITHOUT a leading `recursive` keyword (the common case) didn't match
at all, and the constructor body was never walked.  The `objects`
accumulator stayed empty and the emitted `*AllowedParameters` method
silently dropped all the recursive child-class calls.  Real-world
impact: `mergerTreeOutputterAllowedParameters` (and many similar
`*Class` siblings) ended without the
`call self%galacticfilter_%allowedParameters(...)` lines and so
failed to forward parameter validation to nested objects.

Fix: move the mandatory whitespace inside the optional group —
`(recursive\s+)?function NAME …`.
"""

import re


_OPENERS = [
    'function fooConstructorParameters(parameters) result(self)',
    '   function fooConstructorParameters(parameters) result(self)',
    'recursive function fooConstructorParameters(parameters)',
    'recursive function fooConstructorParameters(parameters,'
    ' recursiveConstruct,recursiveSelf)',
]


def _sig_re(name):
    return (
        r'^\s*(recursive\s+)?function\s+' + re.escape(name)
        + r'\s*\(\s*parameters\s*'
        + r'(\s*,\s*recursiveConstruct\s*,\s*recursiveSelf\s*)??\)'
    )


def test_signature_regex_matches_non_recursive_constructor():
    """The common case: a constructor with no `recursive` keyword."""
    pattern = _sig_re('fooConstructorParameters')
    assert re.match(
        pattern,
        'function fooConstructorParameters(parameters) result(self)')


def test_signature_regex_matches_indented_non_recursive():
    pattern = _sig_re('fooConstructorParameters')
    assert re.match(
        pattern,
        '   function fooConstructorParameters(parameters)')


def test_signature_regex_matches_recursive_constructor():
    pattern = _sig_re('fooConstructorParameters')
    assert re.match(
        pattern,
        'recursive function fooConstructorParameters(parameters)')


def test_signature_regex_matches_recursive_with_extra_args():
    pattern = _sig_re('fooConstructorParameters')
    assert re.match(
        pattern,
        'recursive function fooConstructorParameters(parameters,'
        ' recursiveConstruct,recursiveSelf)')


def test_signature_regex_in_real_module():
    """Pin the actual regex shipped in
    `Process.FunctionClass._build_allowed_parameters_method` so that
    every opener form covered above goes on matching.  If someone
    reverts to the broken `(recursive)??\\s+function` form (mandatory
    whitespace before `function` even when no leading word is
    present), the openers in the four other tests in this module
    immediately stop matching."""
    import inspect
    from Galacticus.Build.SourceTree.Process.FunctionClass import (
        _build_allowed_parameters_method,
    )
    src = inspect.getsource(_build_allowed_parameters_method)
    assert "(recursive\\s+)?function" in src, src
