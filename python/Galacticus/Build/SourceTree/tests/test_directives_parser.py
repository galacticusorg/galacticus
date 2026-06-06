r"""Regression tests for `Galacticus.Build.SourceTree.Parse.Directives`.

Three bug classes:

  1. A single `!![ … !!]` block with multiple sibling tags must yield ONE
     directive node per tag.  An earlier draft passed the whole block to
     the XML parser, which auto-wrapped multi-sibling content in a
     synthetic `<root>` element — leaving a bogus `root` directive in the
     tree that `post_process_directives` then complained about.

  2. A self-closing tag whose attributes contain `/` (typically URLs in
     `<reference>`) was missed by the original `<tag\s[^>]*\/>` regex,
     because `[^>]` was inadvertently overrestrictive on text containing
     `/`.  The fix uses `.` (any char incl. `/`) and matches Perl's
     regex.

  3. `post_process_directives` must forgive directives whose type is in
     the NonProcessed exemption list (e.g. `<methods>` injected late by
     a code-generating hook) but must still flag any other unprocessed
     directive — that error is the build's safety net against typos in
     Process hook names.
"""

import pytest

from Galacticus.Build.SourceTree                   import (
    parse_code, walk_tree, serialize,
)
from Galacticus.Build.SourceTree.Parse.Directives  import post_process_directives


def _directive_types(tree):
    return [n.get('type') for n in walk_tree(tree) if n.get('directive')]


def test_multi_sibling_directive_block_emits_one_node_per_tag():
    """Two `<constant/>` siblings inside one `!![ … !!]` block produce two
    directive nodes — not a synthetic `root` wrapper directive."""
    source = (
        "module foo\n"
        "  !![\n"
        "  <constant variable=\"a\" value=\"1\"/>\n"
        "  <constant variable=\"b\" value=\"2\"/>\n"
        "  !!]\n"
        "end module foo\n"
    )
    tree = parse_code(source, name='<test>', instrument=False)
    types = _directive_types(tree)
    assert types == ['constant', 'constant'], types
    assert 'root' not in types


def test_self_closing_with_slash_in_attribute_value():
    """A self-closing tag with `/` in its attribute value (e.g. a URL) closes
    correctly."""
    source = (
        "module foo\n"
        "  !![\n"
        "  <reference url=\"http://example.com/papers/foo\"/>\n"
        "  !!]\n"
        "end module foo\n"
    )
    tree = parse_code(source, name='<test>', instrument=False)
    types = _directive_types(tree)
    assert types == ['reference'], types


def test_post_process_forgives_non_processed_methods_directive():
    """A `<methods>` block injected after `nonProcessed` ran (and so
    carrying no `processed` flag) is forgiven by the post-process check."""
    source = (
        "module foo\n"
        "  !![\n"
        "  <methods>\n"
        "    <method name=\"foo\" description=\"do foo\"/>\n"
        "  </methods>\n"
        "  !!]\n"
        "end module foo\n"
    )
    tree = parse_code(source, name='<test>', instrument=False)
    # Force the directive into the unprocessed state.
    for n in walk_tree(tree):
        if n.get('type') == 'methods' and n.get('directive') is not None:
            n['directive'].pop('processed', None)
    # Must NOT raise.
    post_process_directives(tree)


def test_post_process_flags_unhandled_directive():
    source = (
        "module foo\n"
        "  !![\n"
        "  <bogusDirective name=\"x\"/>\n"
        "  !!]\n"
        "end module foo\n"
    )
    tree = parse_code(source, name='<test>', instrument=False)
    with pytest.raises(RuntimeError, match="bogusDirective.*was not processed"):
        post_process_directives(tree)
