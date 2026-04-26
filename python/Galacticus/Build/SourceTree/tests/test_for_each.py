# Regression test for `process_for_each` in
# `Galacticus.Build.SourceTree.Process.ForEach`.
#
# Bug: same shape as the `<allocate>` bug.  `<forEach variable="passed">`
# is registered to run before `process_generics`, and uses the variable's
# declaration to determine its rank.  For a generic-templated declaration
#
#     logical, allocatable {Type¦rank} :: passed
#
# the placeholder `{Type¦rank}` does not match `dimension(...)` at first-
# pass time, so `_rank_from_declaration` returned 0 and the emitted
# expansion of the directive's body had:
#
#   * no `do foreach__1=…` index loops
#   * `{index}` replaced with an empty string
#
# Generics then cloned the subroutine subtree (and the emitted code with
# it) for every instance.  The rank-1+ clones got
# `result%note=result%note//value1` — `value1` is the rank-1 dummy
# argument and `//` (string concat) requires scalar operands, so gfortran
# rejected it.
#
# Fix: defer when the declaration still carries an unresolved
# `¦`-bearing placeholder, mirroring the Allocate fix.

import pytest


def _wrap_in_parent(declaration_attrs, forEach_directive_kwargs):
    """Build the minimal AST shape `process_for_each` walks: a parent unit
    with a declaration child and a `<forEach>` directive child."""
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
    forEach_node = {
        'type':       'forEach',
        'directive':  dict(forEach_directive_kwargs),
        'firstChild': None,
        'sibling':    None,
        'parent':     parent,
        'source':     '<test>',
        'line':       1,
    }
    decl_node['sibling']  = forEach_node
    parent['firstChild']  = decl_node
    return parent, forEach_node


def _emitted_iteration(parent_node):
    """Return the auto-generated iteration block inserted after the
    forEach directive, or None if nothing was inserted."""
    child = parent_node['firstChild']
    while child is not None:
        if (child.get('type') == 'code'
                and 'Auto-generated iteration' in (child.get('content') or '')):
            return child['content']
        child = child.get('sibling')
    return None


def test_generic_placeholder_defers_for_each():
    """`{Type¦rank}` in attributes ⇒ no emission, directive stays unprocessed."""
    from Galacticus.Build.SourceTree.Process.ForEach import process_for_each

    parent, fe = _wrap_in_parent(
        declaration_attrs=['allocatable', '{Type¦rank}'],
        forEach_directive_kwargs={
            'variable': 'passed',
            'content':  ' if (.not.passed{index}) result=.false.\n',
        })
    process_for_each(parent, options={})

    assert fe['directive'].get('processed') is not True, fe['directive']
    assert _emitted_iteration(parent) is None


def test_resolved_scalar_emits_no_loop():
    """For a scalar declaration, no `do` loop is emitted and `{index}`
    becomes the empty string."""
    from Galacticus.Build.SourceTree.Process.ForEach import process_for_each

    parent, fe = _wrap_in_parent(
        declaration_attrs=['allocatable'],
        forEach_directive_kwargs={
            'variable': 'passed',
            'content':  ' x=passed{index}\n',
        })
    process_for_each(parent, options={})

    assert fe['directive'].get('processed') is True
    body = _emitted_iteration(parent)
    assert body is not None
    assert 'do foreach__' not in body
    assert ' x=passed\n' in body, body


def test_resolved_rank1_emits_loop_and_index():
    """For a rank-1 declaration, a single `do` loop is emitted and
    `{index}` becomes `(foreach__1)`."""
    from Galacticus.Build.SourceTree.Process.ForEach import process_for_each

    parent, fe = _wrap_in_parent(
        declaration_attrs=['allocatable', 'dimension(:)'],
        forEach_directive_kwargs={
            'variable': 'passed',
            'content':  ' x=passed{index}\n',
        })
    process_for_each(parent, options={})

    assert fe['directive'].get('processed') is True
    body = _emitted_iteration(parent)
    assert body is not None
    assert 'do foreach__1=1,size(passed,dim=1)' in body
    assert ' x=passed(foreach__1)\n' in body, body
    assert 'end do' in body


def test_format_meta_write_paren_count():
    """The auto-generated `write (indexesFormatMeta_,…) v1, v2` line must
    have balanced parentheses.  An earlier draft had an extra `)` that
    gfortran rejected with `Syntax error in WRITE statement`."""
    from Galacticus.Build.SourceTree.Process.ForEach import process_for_each

    parent, fe = _wrap_in_parent(
        declaration_attrs=['allocatable', 'dimension(:)'],
        forEach_directive_kwargs={
            'variable': 'passed',
            'content':  ' x=passed%index%\n',   # `%index%` triggers the format setup
        })
    process_for_each(parent, options={})

    body = _emitted_iteration(parent)
    assert body is not None
    for line in body.splitlines():
        if line.startswith('write (indexesFormatMeta_,'):
            opens  = sum(1 for ch in line if ch == '(')
            closes = sum(1 for ch in line if ch == ')')
            assert opens == closes, (
                f"unbalanced parens ({opens} ( vs {closes} )): {line!r}")
            break
    else:
        raise AssertionError("expected a write(indexesFormatMeta_,…) line, got:\n"
                             + body)


def test_trailing_end_if_in_body_is_separated_from_end_do():
    """If the directive's body ends with `end if\\n`, the emitted iterator
    must put `end do` on its OWN line — not concatenate `end ifend do`."""
    from Galacticus.Build.SourceTree.Process.ForEach import process_for_each

    parent, fe = _wrap_in_parent(
        declaration_attrs=['allocatable', 'dimension(:)'],
        forEach_directive_kwargs={
            'variable': 'passed',
            'content':  (
                ' if (.not.passed{index}) then\n'
                ' x=1\n'
                ' end if\n'
            ),
        })
    process_for_each(parent, options={})

    body = _emitted_iteration(parent)
    assert body is not None
    assert 'end ifend do' not in body, body
    assert 'end if\nend do' in body, body
