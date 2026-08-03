"""Tests for `Galacticus.Build.Directives` — file-based directive scanner.

This module is the on-disk counterpart to the SourceTree-walking
`parse_directives`; it scans a Fortran source file for `!![...!!]` XML
directive blocks using a small state machine (see Directives.py:80-110).

The canonical format in Galacticus source is multi-line:

    !![
    <directive .../>
    !!]

i.e. `!![` and `!!]` are block markers on their own lines, not an inline
wrapping syntax.  Several edge cases of that state machine — !-->
instrumentation lines, leading `!<` strip, &nbsp; substitution, multiple
directives per file, wildcard matching, conditions filter,
set_root_element_type — were previously untested.
"""

import pytest

from Galacticus.Build.Directives import (
    extract_directives,
    extract_directive,
)


def _write_source(tmp_path, body, name='source.F90'):
    """Write `body` to `tmp_path/name` and return the absolute path."""
    p = tmp_path / name
    p.write_text(body)
    return str(p)


# ---------------------------------------------------------------------------
# Empty / missing file
# ---------------------------------------------------------------------------

def test_missing_file_returns_empty_list(tmp_path):
    """A path that doesn't exist must NOT raise, just return []."""
    assert extract_directives(str(tmp_path / 'nope.F90'), 'foo') == []


def test_file_with_no_directives_returns_empty(tmp_path):
    src = _write_source(tmp_path, "module foo\nend module foo\n")
    assert extract_directives(src, 'anything') == []


# ---------------------------------------------------------------------------
# Basic directive extraction
# ---------------------------------------------------------------------------

def test_extracts_self_closing_directive(tmp_path):
    src = _write_source(tmp_path,
        "subroutine foo()\n"
        "  !![\n"
        "  <verbatim name=\"x\" value=\"42\"/>\n"
        "  !!]\n"
        "end subroutine foo\n"
    )
    out = extract_directives(src, 'verbatim')
    assert out == [{'name': 'x', 'value': '42'}]


def test_extracts_directive_with_child_elements(tmp_path):
    src = _write_source(tmp_path,
        "  !![\n"
        "  <foo>\n"
        "    <bar>inner</bar>\n"
        "  </foo>\n"
        "  !!]\n"
    )
    out = extract_directives(src, 'foo')
    assert out == [{'bar': 'inner'}]


def test_extracts_multiple_directives(tmp_path):
    src = _write_source(tmp_path,
        "  !![\n"
        "  <a id=\"1\"/>\n"
        "  !!]\n"
        "code\n"
        "  !![\n"
        "  <a id=\"2\"/>\n"
        "  !!]\n"
        "more code\n"
        "  !![\n"
        "  <a id=\"3\"/>\n"
        "  !!]\n"
    )
    out = extract_directives(src, 'a')
    assert [d['id'] for d in out] == ['1', '2', '3']


def test_strips_bang_lt_prefix_inside_xml_block(tmp_path):
    """Lines inside a directive that start with `!<` (at column 0) have
    that prefix stripped before XML parsing.  Legacy convention from
    the F77-fixed-column-comment era; harmless on free-form sources."""
    src = _write_source(tmp_path,
        "  !![\n"
        "!< <wrap><inner attr=\"v\"/></wrap>\n"
        "  !!]\n"
    )
    out = extract_directives(src, 'wrap')
    assert out == [{'inner': {'attr': 'v'}}]


def test_replaces_nbsp_entity_with_space(tmp_path):
    """`&nbsp;` inside a directive is replaced with a literal space.  Some
    legacy content uses it to preserve formatting; ParseError would
    be the silent failure mode otherwise (XML doesn't define `&nbsp;`)."""
    src = _write_source(tmp_path,
        "  !![\n"
        "  <note>hello&nbsp;world</note>\n"
        "  !!]\n"
    )
    out = extract_directives(src, 'note')
    # `xml_to_dict` turns a text-only element into the bare string; the
    # scanner wraps that in {'content': ...} so callers can attach
    # rootElementType / conditions uniformly.
    assert out == [{'content': 'hello world'}]


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def test_filters_by_directive_name(tmp_path):
    """Directives whose root element does NOT match `directive_name` are
    excluded."""
    src = _write_source(tmp_path,
        "  !![\n  <wanted id=\"1\"/>\n  !!]\n"
        "  !![\n  <other  id=\"2\"/>\n  !!]\n"
        "  !![\n  <wanted id=\"3\"/>\n  !!]\n"
    )
    out = extract_directives(src, 'wanted')
    assert [d['id'] for d in out] == ['1', '3']


def test_wildcard_directive_name_matches_all(tmp_path):
    """`'*'` matches every root element."""
    src = _write_source(tmp_path,
        "  !![\n  <a id=\"1\"/>\n  !!]\n"
        "  !![\n  <b id=\"2\"/>\n  !!]\n"
    )
    out = extract_directives(src, '*')
    assert len(out) == 2


def test_conditions_filter_includes_only_matching(tmp_path):
    src = _write_source(tmp_path,
        "  !![\n  <thing kind=\"A\" id=\"1\"/>\n  !!]\n"
        "  !![\n  <thing kind=\"B\" id=\"2\"/>\n  !!]\n"
        "  !![\n  <thing kind=\"A\" id=\"3\"/>\n  !!]\n"
    )
    out = extract_directives(src, 'thing', conditions={'kind': 'A'})
    assert [d['id'] for d in out] == ['1', '3']


def test_conditions_with_no_match_returns_empty(tmp_path):
    src = _write_source(tmp_path,
        "  !![\n  <thing kind=\"A\"/>\n  !!]\n"
    )
    assert extract_directives(src, 'thing', conditions={'kind': 'Z'}) == []


def test_set_root_element_type_stamps_root_tag(tmp_path):
    """When set_root_element_type=True, each returned dict gains a
    `rootElementType` key carrying the actual XML root tag — needed by
    callers that pass `directive_name='*'` and need to know which element
    they got back."""
    src = _write_source(tmp_path,
        "  !![\n  <alpha v=\"1\"/>\n  !!]\n"
        "  !![\n  <beta  v=\"2\"/>\n  !!]\n"
    )
    out = extract_directives(src, '*', set_root_element_type=True)
    by_type = {d['rootElementType']: d for d in out}
    assert set(by_type) == {'alpha', 'beta'}
    assert by_type['alpha']['v'] == '1'
    assert by_type['beta']['v']  == '2'


# ---------------------------------------------------------------------------
# Instrumentation lines
# ---------------------------------------------------------------------------

def test_skips_instrument_marker_lines(tmp_path):
    """Lines starting with `!-->` are SourceIntrospection markers and must
    NOT be included in the accumulated XML — otherwise they'd corrupt
    parsing of the surrounding directive."""
    src = _write_source(tmp_path,
        "  !![\n"
        "  <foo>\n"
        "!--> some_introspection_marker\n"
        "    <bar>x</bar>\n"
        "  </foo>\n"
        "  !!]\n"
    )
    out = extract_directives(src, 'foo')
    assert out == [{'bar': 'x'}]


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

def test_malformed_xml_raises_RuntimeError(tmp_path):
    """Parse failure is fatal — the caller has no way to
    recover from a malformed embedded XML, so it raises with a message
    naming the file."""
    src = _write_source(tmp_path,
        "  !![\n  <foo>\n  </bar>\n  !!]\n",
        name='broken.F90',
    )
    with pytest.raises(RuntimeError, match=r"failed parsing XML.*broken\.F90"):
        extract_directives(src, 'foo')


# ---------------------------------------------------------------------------
# extract_directive (singular convenience wrapper)
# ---------------------------------------------------------------------------

def test_extract_directive_returns_first_match(tmp_path):
    src = _write_source(tmp_path,
        "  !![\n  <x id=\"1\"/>\n  !!]\n"
        "  !![\n  <x id=\"2\"/>\n  !!]\n"
    )
    out = extract_directive(src, 'x')
    assert out == {'id': '1'}


def test_extract_directive_returns_None_when_no_match(tmp_path):
    src = _write_source(tmp_path,
        "  !![\n  <x id=\"1\"/>\n  !!]\n"
    )
    assert extract_directive(src, 'absent') is None
