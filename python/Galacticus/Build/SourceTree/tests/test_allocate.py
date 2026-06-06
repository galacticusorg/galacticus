"""Regression test for `process_allocate` in
`Galacticus.Build.SourceTree.Process.Allocate`.

Bug: `process_allocate` is registered with `before=['generics']` so it
runs before generic-template expansion.  For a generic-templated
declaration like

    logical, allocatable {Type¦rank} :: passed
    !![ <allocate variable="passed" size="value1"/> !!]

the placeholder `{Type¦rank}` does NOT match `dimension(...)` at this
stage, so `_declaration_rank` returned 0 and emitted `allocate(passed)`
— no parentheses, no size.  Generics then cloned the entire subroutine
subtree (and the emitted code with it) per instance, including
rank-1+ instances where the cloned `, dimension(:) :: passed`
declaration leaves `allocate(passed)` invalid:

    Error: Array specification required in ALLOCATE statement at (1)

Fix: when allocate sees an unresolved `¦`-bearing placeholder in the
declaration's attributes / type / variables, DEFER (don't mark the
directive as processed and don't emit code).  Generics' clone-then-
reparse step preserves the `!![…!!]` markers of un-processed
directives, so the inner `process_tree` call on each cloned subtree
re-invokes `process_allocate` with a fully-resolved declaration.
"""

import pytest


def _wrap_in_parent(declaration_attrs, allocate_directive):
    """Build the minimal AST shape `process_allocate` walks: a parent unit
    with a declaration child and an `<allocate>` directive child.
    """
    parent = {
        'type':       'unit',
        'firstChild': None,
        'sibling':    None,
        'parent':     None,
    }
    decl_node = {
        'type':         'declaration',
        'declarations': [{
            'intrinsic':     'logical',
            'type':          None,
            'attributes':    list(declaration_attrs),
            'variables':     ['passed'],
            'variableNames': ['passed'],
        }],
        'firstChild': None,
        'sibling':    None,
        'parent':     parent,
    }
    allocate_node = {
        'type':       'allocate',
        'directive':  dict(allocate_directive),
        'firstChild': None,
        'sibling':    None,
        'parent':     parent,
        'source':     '<test>',
        'line':       1,
    }
    decl_node['sibling']  = allocate_node
    parent['firstChild']  = decl_node
    return parent, allocate_node


def _emitted_allocate_lines(parent_node):
    """Return the auto-generated `allocate(...)` lines inserted by
    `process_allocate` after the directive, or `[]` if none were inserted.
    """
    out = []
    child = parent_node['firstChild']
    while child is not None:
        if child.get('type') == 'code' and 'allocate(' in (child.get('content') or ''):
            for line in child['content'].splitlines():
                if 'allocate(' in line:
                    out.append(line.strip())
        child = child.get('sibling')
    return out


def test_generic_placeholder_defers_allocate():
    """`{Type¦rank}` in attributes ⇒ no emission, directive stays unprocessed."""
    from Galacticus.Build.SourceTree.Process.Allocate import process_allocate

    parent, alloc = _wrap_in_parent(
        declaration_attrs=['allocatable', '{Type¦rank}'],
        allocate_directive={'variable': 'passed', 'size': 'value1'})
    process_allocate(parent, options={})

    assert alloc['directive'].get('processed') is not True, alloc['directive']
    assert _emitted_allocate_lines(parent) == [], _emitted_allocate_lines(parent)


def test_resolved_scalar_emits_no_size():
    """After expansion to a scalar (`{Type¦rank}` → empty), allocate emits
    the bare `allocate(passed)` form."""
    from Galacticus.Build.SourceTree.Process.Allocate import process_allocate

    parent, alloc = _wrap_in_parent(
        declaration_attrs=['allocatable'],
        allocate_directive={'variable': 'passed', 'size': 'value1'})
    process_allocate(parent, options={})

    assert alloc['directive'].get('processed') is True
    assert _emitted_allocate_lines(parent) == ['allocate(passed)']


def test_resolved_rank1_emits_size_form():
    """After expansion to rank-1 (`{Type¦rank}` → `dimension(:)`), allocate
    emits `allocate(passed(size(value1,dim=1)))`."""
    from Galacticus.Build.SourceTree.Process.Allocate import process_allocate

    parent, alloc = _wrap_in_parent(
        declaration_attrs=['allocatable', 'dimension(:)'],
        allocate_directive={'variable': 'passed', 'size': 'value1'})
    process_allocate(parent, options={})

    assert alloc['directive'].get('processed') is True
    assert _emitted_allocate_lines(parent) == [
        'allocate(passed(size(value1,dim=1)))']


def test_explicit_rank_overrides_declaration_check():
    """If the directive itself carries `rank=N`, the declaration is not
    consulted, and emission proceeds even with an unresolved placeholder."""
    from Galacticus.Build.SourceTree.Process.Allocate import process_allocate

    parent, alloc = _wrap_in_parent(
        declaration_attrs=['allocatable', '{Type¦rank}'],
        allocate_directive={'variable': 'passed', 'size': 'value1', 'rank': 2})
    process_allocate(parent, options={})

    assert alloc['directive'].get('processed') is True
    assert _emitted_allocate_lines(parent) == [
        'allocate(passed(size(value1,dim=1),size(value1,dim=2)))']
