"""Regression tests for derived-type definition opener recognition.

Bug class: the opener pattern omitted `bind(c)` from its attribute
alternation and admitted no whitespace inside `extends(...)`, so
`type, bind(c) :: unitType` and `type, extends(hdf5Group             ) ::
hdf5File` both failed to match.

A missed opener is not merely an unlabelled node.  The unit parser never opens
a `type` node, so the type's components are absorbed into the enclosing
scope's declaration node as though they were ordinary variables, and the
`end type` is stranded as an orphan code node — silently, with nothing
raised.  Four types in the tree were affected.

The pattern had been copied into eight modules which then drifted apart, so
these tests also pin every consumer to the shared definition in
`Fortran.Utils`: fixing one copy previously just moved the failure to the
next one (`DeepCopyActions` raised "unable to parse type definition opener"
once the type nodes it had never seen started existing).
"""

import re

import pytest

from Fortran.Utils import (
    TYPE_ATTRIBUTES_TAGGED, UNIT_OPENERS, type_opener_regex,
)
from Galacticus.Build.SourceTree import parse_code, children, walk_tree


TYPE_OPENER = UNIT_OPENERS['type']['regex']

# The four openers in the tree that the old pattern missed, and why.
MISSED_OPENERS = [
    ("  type, bind(C) :: unitType",                          'unitType'),
    ("  type, public, bind(c) :: gsl_sf_result",             'gsl_sf_result'),
    ("  type, bind(C) :: hdf5VlenC",                         'hdf5VlenC'),
    ("  type, extends(hdf5Group             ) :: hdf5File",  'hdf5File'),
]


@pytest.mark.parametrize("opener,name", MISSED_OPENERS)
def test_previously_missed_openers_now_match(opener, name):
    m = TYPE_OPENER.match(opener)
    assert m is not None, f"opener not recognized: {opener!r}"
    assert m.group('name') == name


@pytest.mark.parametrize("opener,name", [
    ("  type :: plain",                             'plain'),
    ("  type plain",                                'plain'),
    ("  type, public :: pub",                       'pub'),
    ("  type, abstract, extends(foo) :: bar",       'bar'),
    ("  type, extends(foo), abstract :: bar",       'bar'),
    ("  type, bind( c ) :: spaced",                 'spaced'),
    ("  TYPE, BIND(C) :: shouty",                   'shouty'),
])
def test_opener_forms_match(opener, name):
    m = TYPE_OPENER.match(opener)
    assert m is not None, f"opener not recognized: {opener!r}"
    assert m.group('name') == name


@pytest.mark.parametrize("line", [
    "  type(unitType) :: x",                # a variable, not a definition
    "  type(c_ptr), pointer :: p",
    "  type(varying_string), dimension(:) :: names",
])
def test_variable_declarations_are_not_openers(line):
    assert TYPE_OPENER.match(line) is None, \
        f"variable declaration wrongly claimed as a type opener: {line!r}"


def test_bind_c_type_gets_its_own_node_with_its_components():
    """The whole point of the fix: components stay inside the type.

    Before, `numerator`/`denominator` were recorded as *module-scope*
    declarations, `unitType` had no node at all, and `end type unitType`
    survived as a bare code node.
    """
    source = (
        "module demo\n"
        "  implicit none\n"
        "  type, bind(c) :: unitType\n"
        "     integer :: numerator\n"
        "     integer :: denominator\n"
        "  end type unitType\n"
        "  integer :: moduleVariable\n"
        "end module demo\n"
    )
    tree = parse_code(source, name='<test>', instrument=False)

    def declared_names(node):
        """Every variable declared directly in `node`'s own scope.

        The parser lower-cases variable names (Fortran is case-insensitive),
        so compare against lower-case expectations.
        """
        names = []
        for child in children(node):
            if child.get('type') == 'declaration':
                for declaration in child.get('declarations') or []:
                    names.extend(declaration.get('variables') or [])
        return [name.lower() for name in names]

    types = [n for n in walk_tree(tree) if n.get('type') == 'type']
    assert [n.get('name') for n in types] == ['unitType']

    # The components belong to the type ...
    assert declared_names(types[0]) == ['numerator', 'denominator']

    # ... and are *not* leaked into the enclosing module, which keeps only its
    # own variable.  This is the assertion that fails on the unfixed parser.
    module = next(n for n in walk_tree(tree) if n.get('type') == 'module')
    assert declared_names(module) == ['modulevariable']

    # `end type` was consumed as the closer, not stranded as orphan code.
    for node in walk_tree(tree):
        if node.get('type') == 'code':
            assert 'end type' not in node.get('content', '').lower(), \
                "`end type` left as an orphan code node"


def test_padded_extends_records_the_parent():
    """`extends(hdf5Group             )` must yield the parent, not None.

    The padding is what hid `hdf5File` from the parser; a pattern that
    matched the opener but dropped the parent would silently break the
    inheritance map instead.
    """
    rx = type_opener_regex(attributes=TYPE_ATTRIBUTES_TAGGED)
    m  = rx.match("  type, extends(hdf5Group             ) :: hdf5File")
    assert m is not None
    assert m.group('name')     == 'hdf5File'
    assert m.group('extends')  == 'hdf5Group'
    assert m.group('abstract') is None


def test_tagged_attributes_capture_abstract_in_any_order():
    rx = type_opener_regex(attributes=TYPE_ATTRIBUTES_TAGGED)
    for opener in ("type, abstract, extends(foo) :: bar",
                   "type, extends(foo), abstract :: bar"):
        m = rx.match(opener)
        assert m is not None,               opener
        assert m.group('abstract')is not None, opener
        assert m.group('extends') == 'foo', opener
        assert m.group('name')    == 'bar', opener


def test_abstract_is_matched_case_insensitively():
    """Fortran is case-insensitive; an `ABSTRACT` type is still abstract.

    The hand-rolled attribute splitting this replaced compared with
    `attr == 'abstract'`, so it silently mis-classified an upper-case
    `ABSTRACT` type as concrete.
    """
    rx = type_opener_regex(attributes=TYPE_ATTRIBUTES_TAGGED)
    m  = rx.match("TYPE, ABSTRACT, EXTENDS(Foo) :: Bar")
    assert m is not None
    assert m.group('abstract') is not None


def test_narrowed_name_pattern_captures_the_whole_name():
    """Callers scanning one functionClass family pass a prefixed name."""
    rx = type_opener_regex(name=re.escape('massDistribution') + r'[a-z0-9_]+',
                           attributes=TYPE_ATTRIBUTES_TAGGED)
    m = rx.match("  type, extends(massDistributionClass) :: massDistributionGaussian")
    assert m is not None
    assert m.group('name')    == 'massDistributionGaussian'
    assert m.group('extends') == 'massDistributionClass'
    assert rx.match("  type, extends(foo) :: someOtherThing") is None


@pytest.mark.parametrize("opener,name", MISSED_OPENERS)
def test_every_consumer_parses_the_previously_missed_openers(opener, name):
    """Each site that carried its own copy of the pattern must cope.

    Fixing only the canonical regex broke the build: `DeepCopyActions` and
    `StateStorable` raise on an opener they cannot parse, and they only began
    seeing these types once the canonical fix made their nodes exist.
    """
    from Galacticus.Build.TypeScan import TYPE_OPENER_RE
    from Galacticus.Build.SourceTree.Process.DeepCopyActions import (
        _parse_type_opener as deep_copy_parse,
    )
    from Galacticus.Build.SourceTree.Process.StateStorable import (
        _parse_type_opener as state_store_parse,
    )

    assert TYPE_OPENER_RE.match(opener) is not None

    # Neither generator may raise; both must agree on the name.
    assert deep_copy_parse(opener)[0]   == name
    assert state_store_parse(opener)[0] == name


def test_bind_c_types_have_no_parent_and_so_reach_no_generator():
    """`bind(c)` types are inert for the hierarchy-walking generators.

    A C-interoperable derived type can neither extend another type nor be
    extended, so it can never join a functionClass / stateStorable /
    deepCopyActions chain.  That is why making these types visible needs no
    explicit exclusion — the generators emit code by walking parents, and
    these have none.
    """
    from Galacticus.Build.SourceTree.Process.DeepCopyActions import (
        _parse_type_opener,
    )
    for opener, name in MISSED_OPENERS:
        parsed = _parse_type_opener(opener)
        if 'bind' in opener.lower():
            assert parsed == (name, None, False)
