# Regression test for `_strip_processed_directive_content` in
# `Galacticus.Build.SourceTree.Process.Generics`.
#
# Bug: when `_expand_subtree` clones a generic-templated subtree, it
# blanks the firstChild content of every `processed=True` directive so
# the reparse step doesn't recreate a fresh directive node and have
# the inner `process_tree` re-emit its generated code (e.g. an
# `<optionalArgument>` setter).  But that strip was over-eager — it
# also blanked `<methods>` directives, which are tagged `processed=True`
# by `nonProcessed.process_non_processed` purely to silence
# `post_process_directives`'s safety check, NOT because they emit
# code.  `<methods>` blocks carry the per-class method descriptions
# that the inner `process_tree`'s `classDocumentation` pass relies on,
# so stripping them caused every generic-expanded class
# (e.g. `polyRankdouble`) to arrive at the doc consumer with empty
# `<descriptions>` — producing spurious "missing method descriptions"
# warnings for every type-bound procedure.
#
# Fix: short-circuit the strip when the directive's type is in
# `is_non_processed_type`'s exemption list.

from Galacticus.Build.SourceTree.Process.Generics import (
    _strip_processed_directive_content,
)


def _directive_node(node_type, content, processed=True):
    first_child = {
        'type':       'code',
        'content':    content,
        'firstChild': None,
        'sibling':    None,
        'parent':     None,
    }
    node = {
        'type':       node_type,
        'directive':  {'processed': processed},
        'firstChild': first_child,
        'sibling':    None,
        'parent':     None,
    }
    first_child['parent'] = node
    return node


def _root_with_children(*children):
    root = {
        'type':       'unit',
        'firstChild': None,
        'sibling':    None,
        'parent':     None,
    }
    for i, child in enumerate(children):
        child['parent'] = root
        child['sibling'] = children[i + 1] if i + 1 < len(children) else None
    root['firstChild'] = children[0] if children else None
    return root


def test_methods_directive_content_preserved():
    """`<methods>` is in the `nonProcessed` exemption list — its
    firstChild content must survive the strip step so the inner
    `process_tree` can populate `<descriptions>` from it."""
    methods_node = _directive_node(
        'methods',
        '!![\n<methods><method method="rank"/></methods>\n!!]\n')
    root = _root_with_children(methods_node)
    _strip_processed_directive_content(root)
    assert methods_node['firstChild']['content'] != ''
    assert '<methods>' in methods_node['firstChild']['content']


def test_code_generating_directive_content_blanked():
    """A code-generating directive whose hook has already emitted its
    output (e.g. `<optionalArgument>`) must have its firstChild content
    blanked, so the reparse below doesn't recreate the directive and
    the inner `process_tree` re-emits the same code."""
    optional_arg = _directive_node(
        'optionalArgument',
        '!![\n<optionalArgument name="foo" defaultsTo="0"/>\n!!]\n')
    root = _root_with_children(optional_arg)
    _strip_processed_directive_content(root)
    assert optional_arg['firstChild']['content'] == ''


def test_unprocessed_directive_content_preserved():
    """A directive that hasn't been processed yet (`processed != True`)
    is left alone — `_strip_processed_directive_content` only strips the
    ones whose hook has already run."""
    pending = _directive_node(
        'allocate',
        '!![\n<allocate variable="x" size="y"/>\n!!]\n',
        processed=False)
    root = _root_with_children(pending)
    _strip_processed_directive_content(root)
    assert pending['firstChild']['content'] != ''


def test_task_suffix_directives_also_preserved():
    """`is_non_processed_type` also forgives `*Task` directives.  Verify
    a `nodeOperatorTask` survives the strip even though its `processed`
    flag is set."""
    task_node = _directive_node(
        'nodeOperatorTask',
        '!![\n<nodeOperatorTask name="foo"/>\n!!]\n')
    root = _root_with_children(task_node)
    _strip_processed_directive_content(root)
    assert '<nodeOperatorTask' in task_node['firstChild']['content']


def test_mixed_tree_strips_only_code_generating_directives():
    """Sanity-check the integrated behaviour: a root with one
    `<methods>` and one `<optionalArgument>` directive sees only the
    latter blanked."""
    methods_node = _directive_node(
        'methods', '!![ <methods><method method="rank"/></methods> !!]')
    optional_arg = _directive_node(
        'optionalArgument', '!![ <optionalArgument name="x"/> !!]')
    root = _root_with_children(methods_node, optional_arg)
    _strip_processed_directive_content(root)
    assert methods_node['firstChild']['content'] != ''
    assert optional_arg['firstChild']['content'] == ''
