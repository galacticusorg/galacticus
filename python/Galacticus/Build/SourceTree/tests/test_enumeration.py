"""Regression test for `_insert_code_tree` in
`Galacticus.Build.SourceTree.Process.Enumeration`.

Bug: `_insert_code_tree` parsed its synthetic source text with a
fixed `name='Enumeration'`, so every node inserted by the enumeration
emitter (the `enumeration<N>Type` derived type, the
`enumeration<N>IsEqual` function, validators, encoders, decoders, …)
carried the literal string `'Enumeration'` in its `source` field
instead of the real `.F90` filename.

Fix: parse the synthetic enumeration source text with the OUTER
source's name, so its nodes inherit the outer file's `source` field.
"""

from Galacticus.Build.SourceTree.Process.Enumeration import _insert_code_tree


def _outer_module_node(source):
    """Build a minimal `module` parent whose `source` attribute matches
    `source` — the shape `_insert_code_tree` is called against in the
    real pipeline.  We use a self-closing `contains` marker so the
    `insert_post_contains` inserter has somewhere to anchor."""
    contains = {
        'type':       'contains',
        'firstChild': None,
        'sibling':    None,
        'parent':     None,
        'source':     source,
    }
    parent = {
        'type':       'module',
        'name':       'foo',
        'firstChild': contains,
        'sibling':    None,
        'parent':     None,
        'source':     source,
    }
    contains['parent'] = parent
    return parent


def test_insert_code_tree_inherits_parents_source():
    """Inserted nodes carry the OUTER source name, not 'Enumeration'."""
    parent = _outer_module_node('utility.foo.F90')
    _insert_code_tree(
        parent,
        '\n  type, public :: barType\n  end type barType\n',
        inserter=lambda p, kids: setattr(_inserted_holder, 'kids', kids),
    )
    inserted = _inserted_holder.kids
    assert inserted, "inserter callback didn't capture children"
    assert all(node.get('source') == 'utility.foo.F90' for node in inserted), \
        [node.get('source') for node in inserted]


def test_insert_code_tree_falls_back_when_parent_has_no_source():
    """If the parent has no `source` attribute (defensive), keep the
    legacy `'Enumeration'` placeholder rather than `None`."""
    parent = {
        'type':       'module',
        'name':       'foo',
        'firstChild': None,
        'sibling':    None,
        'parent':     None,
    }
    _insert_code_tree(
        parent,
        '\n  type :: barType\n  end type barType\n',
        inserter=lambda p, kids: setattr(_inserted_holder, 'kids', kids),
    )
    assert all(node.get('source') == 'Enumeration'
               for node in _inserted_holder.kids)


class _inserted_holder:
    """Module-level container so the test's `inserter` callback can
    stash its `kids` argument for the assert step."""
    kids = ()
