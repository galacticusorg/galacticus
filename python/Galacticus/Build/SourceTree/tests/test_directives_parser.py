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
     `/`.  The fix uses `.` (any char incl. `/`).

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
    directive nodes — not a synthetic `root` wrapper directive. (Directives are
    validated against their schemas, so the fixtures carry the attributes
    `constant.xsd` requires.)"""
    source = (
        "module foo\n"
        "  !![\n"
        "  <constant variable=\"a\" value=\"1\" description=\"a\" reference=\"none\"/>\n"
        "  <constant variable=\"b\" value=\"2\" description=\"b\" reference=\"none\"/>\n"
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
        "    <method method=\"foo\" description=\"do foo\"/>\n"
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


def test_file_based_directive_schemas_are_applied():
    """Directive schemas such as `objectBuilder.xsd` live in the repository's
    `schema/` directory. The directory was once built as
    `os.path.join(EXEC_PATH, "/schemas")`, which discards `EXEC_PATH` and names
    a directory that does not exist, so every file-based schema was silently
    skipped. A directive violating its schema must be rejected."""
    import os
    import Galacticus.Build.SourceTree.Parse.Directives as directives
    if not directives._HAS_LXML:
        pytest.skip("lxml is not installed")
    assert directives.SCHEMAS_DIR is not None
    assert os.path.isfile(os.path.join(directives.SCHEMAS_DIR, 'objectBuilder.xsd'))
    context = {'type': 'file', 'name': 'test.F90'}
    directives._validate_directive(
        'objectBuilder',
        '<objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>',
        context, 1)
    with pytest.raises(RuntimeError, match="source"):
        directives._validate_directive(
            'objectBuilder',
            '<objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_"/>',
            context, 1)
