"""Regression tests for `Galacticus.Build.SourceTree.Parse.Declarations`.

Covers two bug classes in the declarations parser:

  1. `parse_declaration` mishandled type-specs whose body itself contained
     `(`/`)` — e.g. `character(len=len(tagName))` was clipped at the first
     `)` to leave `len=len(tagName` as the type.  The fix moved to
     `extract_bracketed` so the OUTER paren group is matched correctly.

  2. `declaration_exists` missed declarations stored as `name=initializer`
     tokens (e.g. `warnObjectBuilder0__=.false.`).  Without stripping the
     `=value` tail the lookup of `warnObjectBuilder0__` failed and the
     caller emitted a duplicate declaration in the same scope.
"""

import pytest

from Galacticus.Build.SourceTree.Parse.Declarations import (
    parse_declaration, declaration_exists,
)


def test_parse_declaration_balanced_parens_in_type_spec():
    """`character(len=len(tagName))` keeps its inner parens — the parser must
    consume the OUTER balanced paren group, not stop at the first `)`."""
    decl = parse_declaration(
        "character(len=len(tagName)) :: foo")
    assert decl is not None
    assert decl['intrinsic']  == 'character'
    assert decl['type']       == 'len=len(tagName)'
    assert decl['attributes'] == []
    assert decl['variables']  == ['foo']


def test_parse_declaration_dimension_attribute_with_colons():
    """`dimension(:,:)` carries colons inside parens; the type-spec/attribute
    split must not get confused by them."""
    decl = parse_declaration(
        "real(kind=8), dimension(:,:), allocatable :: matrix")
    assert decl is not None
    assert decl['intrinsic']  == 'real'
    assert decl['type']       == 'kind=8'
    assert 'dimension(:,:)' in decl['attributes']
    assert 'allocatable'    in decl['attributes']
    assert decl['variables'] == ['matrix']


def test_parse_declaration_procedure_pointer():
    """`procedure(template), pointer, nopass :: foo` — used by every
    procedure-pointer data member.  The `pointer` attribute is what tells
    `_build_assignment_method` to use `=>` rather than `=`."""
    decl = parse_declaration(
        "procedure(myInterface), pointer, nopass :: callback")
    assert decl is not None
    assert decl['intrinsic']            == 'procedure'
    assert decl['type']                 == 'myInterface'
    assert 'pointer' in decl['attributes']
    assert 'nopass'  in decl['attributes']
    assert decl['variables']            == ['callback']


def test_parse_declaration_returns_none_on_non_declaration():
    assert parse_declaration("call doSomething()")          is None
    assert parse_declaration("if (x .gt. 0) y = 1")          is None
    assert parse_declaration("integer x")                    is None  # no `::`


def _decl_node(declarations):
    """Wrap a list of declaration dicts in a minimal AST node."""
    return {
        'firstChild': {
            'type':         'declaration',
            'declarations': declarations,
            'sibling':      None,
            'firstChild':   None,
        },
    }


def test_declaration_exists_strips_initializer():
    """`warnObjectBuilder0__=.false.` should be findable by the bare name."""
    parent = _decl_node([{
        'intrinsic':     'logical',
        'type':          None,
        'attributes':    [],
        'variables':     ['warnObjectBuilder0__=.false.'],
        'variableNames': [],
    }])
    assert declaration_exists(parent, 'warnObjectBuilder0__')
    assert declaration_exists(parent, 'WarnObjectBuilder0__')   # case-insensitive
    assert not declaration_exists(parent, 'somethingElse')


def test_declaration_exists_uses_variableNames():
    """When the parser populated `variableNames` we should fall back to it."""
    parent = _decl_node([{
        'intrinsic':     'integer',
        'type':          None,
        'attributes':    [],
        'variables':     ['n=0'],
        'variableNames': ['n'],
    }])
    assert declaration_exists(parent, 'n')
