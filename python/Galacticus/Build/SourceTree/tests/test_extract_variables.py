"""Regression test for `Fortran.Utils.extract_variables`.

Bug: `extract_variables` lowercased the entire post-`::` text up front,
which mangled string-literal initializers — a `character(len=56) ::
fileNameCoolingFunction='cooling/cooling_function_Atomic_CIE_Cloudy.hdf5'`
became
`character(len=56) :: filenamecoolingfunction='cooling/cooling_function_atomic_cie_cloudy.hdf5'`
in the emitted `.p.F90`.  Fortran is case-insensitive for identifiers but
the bytes inside a string literal are not, and on a case-sensitive
filesystem the lowercased filename does not exist, so the build runs but
the model fails to open its data file at runtime.

The fix stashes string literals as opaque placeholders before any
lowercase/whitespace pass, and restores them at the end — which also
fixes the latent bug where a comma inside a string literal would split
the variable list incorrectly.
"""

from Fortran.Utils import extract_variables


def test_string_literal_initializer_keeps_case():
    out = extract_variables(
        "fileNameCoolingFunction="
        "'cooling/cooling_function_Atomic_CIE_Cloudy.hdf5'",
        keep_qualifiers=True, lower_case=True)
    assert out == [
        "filenamecoolingfunction="
        "'cooling/cooling_function_Atomic_CIE_Cloudy.hdf5'"
    ]


def test_multiple_string_literal_initializers_keep_case():
    out = extract_variables(
        "fileA='Foo.hdf5', fileB='Bar.hdf5', n=42",
        keep_qualifiers=True, lower_case=True)
    assert out == ["filea='Foo.hdf5'", "fileb='Bar.hdf5'", "n=42"]


def test_comma_inside_string_does_not_split():
    """Pre-fix, the comma inside `'Hello, World'` would split the variable
    list into two bogus entries `label='Hello` and `World'`."""
    out = extract_variables(
        "label='Hello, World'", keep_qualifiers=True, lower_case=True)
    assert out == ["label='Hello, World'"]


def test_doubled_quote_inside_single_quoted_string():
    """Fortran escapes a single quote inside a `'…'` string by doubling it
    (`''`).  The literal extractor must consume both quotes as part of the
    same string."""
    out = extract_variables(
        "name='Don''t'", keep_qualifiers=True, lower_case=True)
    assert out == ["name='Don''t'"]


def test_double_quoted_string_keeps_case():
    out = extract_variables(
        'msg="Hello World"', keep_qualifiers=True, lower_case=True)
    assert out == ['msg="Hello World"']


def test_no_strings_unchanged():
    """The fix is purely additive — non-string inputs go through the same
    pipeline as before."""
    out = extract_variables(
        'foo, bar(3,4), baz', keep_qualifiers=True, lower_case=True)
    assert out == ['foo', 'bar(3,4)', 'baz']


def test_drop_qualifiers_strips_initializer():
    """The `keep_qualifiers=False` path must continue to strip the
    `=initializer` text entirely (string literal or not)."""
    out = extract_variables(
        "fileNameA='Foo.hdf5', n=42",
        keep_qualifiers=False, lower_case=True)
    assert out == ['filenamea', 'n']
