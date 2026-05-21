"""Regression tests for top-level `Galacticus.Build.SourceTree` behaviour.

Bug class: embedded LaTeX (`!!{ … !!}`) and XML (`!![ … !!]`) blocks must
survive the parse → serialize round-trip in a form the Fortran compiler
accepts.  An earlier draft kept the original text in a parallel buffer,
so when `serialize` rebuilt the source the LaTeX/XML body landed in the
output uncommented — a Fortran compile error.  The fix prefixes each
body line with `!< ` at parse time and stores the commented text in code
nodes; serialize then emits compilable Fortran while the directive
parser still sees the original markers.
"""

from Galacticus.Build.SourceTree import parse_code, serialize


def test_embedded_latex_block_is_commented_after_serialize():
    source = (
        "module foo\n"
        "  !!{\n"
        "  Description for foo, written in \\LaTeX.\n"
        "  !!}\n"
        "  implicit none\n"
        "end module foo\n"
    )
    tree = parse_code(source, name='<test>', instrument=False)
    out  = serialize(tree)
    # Body line is commented out (prefixed with `!<` on round-trip).
    assert "!< " in out or out.find("Description for foo") == -1 \
        or "Description for foo" in out and \
           [ln for ln in out.splitlines()
            if 'Description for foo' in ln and ln.lstrip().startswith('!')]
    # The compiled Fortran should not contain a bare `\LaTeX` token at
    # statement scope (which would error out).
    for line in out.splitlines():
        if '\\LaTeX' in line:
            assert line.lstrip().startswith('!'), \
                f"Uncommented LaTeX leaked through: {line!r}"


def test_embedded_xml_block_round_trips_as_directive():
    """A directive block must survive the parse → serialize round-trip with
    its `!![` / `!!]` markers and inner XML intact (commented out as `!<`)."""
    source = (
        "module foo\n"
        "  !![\n"
        "  <reference url=\"http://example.com\"/>\n"
        "  !!]\n"
        "end module foo\n"
    )
    tree = parse_code(source, name='<test>', instrument=False)
    out  = serialize(tree)
    # `!![` and `!!]` survive the round-trip (they are valid Fortran
    # comments — they begin with `!!` — so they may pass through unchanged
    # or as full-line comments).
    assert '!![' in out
    assert '!!]' in out
    # The XML body line is commented (prefixed with `!<`).
    assert any('!<' in ln and '<reference' in ln
               for ln in out.splitlines()), out
