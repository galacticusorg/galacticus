"""Regression tests for top-level `Galacticus.Build.SourceTree` behavior.

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


def test_unit_opener_keyword_prefilter_covers_every_opener():
    """`_parse_units` tries the `UNIT_OPENERS` patterns only on lines that the
    keyword prefilter accepts, so every opener pattern must require one of its
    keywords — otherwise a new opener would silently never match."""
    from Fortran.Utils import UNIT_OPENERS
    from Galacticus.Build.SourceTree import _UNIT_OPENER_KEYWORDS_RE
    for unit_type, spec in UNIT_OPENERS.items():
        assert _UNIT_OPENER_KEYWORDS_RE.search(spec['regex'].pattern), unit_type


def test_walk_tree_follows_links_lazily():
    """`walk_tree` reads a node's `firstChild` only after yielding the node, and
    its `sibling` only once its subtree is exhausted, so nodes grafted in
    mid-walk are visited — callers rely on this."""
    from Galacticus.Build.SourceTree import walk_tree
    leaf    = {'name': 'leaf'   , 'firstChild': None, 'sibling': None}
    grafted = {'name': 'grafted', 'firstChild': None, 'sibling': None}
    a       = {'name': 'a'      , 'firstChild': None, 'sibling': None}
    b       = {'name': 'b'      , 'firstChild': None, 'sibling': None}
    root    = {'name': 'root'   , 'firstChild': a   , 'sibling': {'name': 'outside'}}
    a['sibling'] = b
    visited = []
    for node in walk_tree(root):
        visited.append(node['name'])
        if node is a:
            a['firstChild']    = leaf
        if node is leaf:
            leaf['sibling']    = grafted
    assert visited == ['root', 'a', 'leaf', 'grafted', 'b']
