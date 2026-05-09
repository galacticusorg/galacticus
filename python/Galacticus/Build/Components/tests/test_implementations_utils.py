"""Tests for `Galacticus.Build.Components.Implementations.Utils`.

Two pipeline hooks (`Implementation_Is_Active`,
`Implementation_Function_Iterator`) plus three free helpers shared with
sister modules (`has_real_evolvers`, `has_real_non_trivial_evolvers`,
`list_real_evolvers`).  The free helpers are pure predicates and the
main test target -- they're called from many emit sites and any
regression here cascades through the build.
"""

from Galacticus.Build.Components.Implementations.Utils import (
    Implementation_Is_Active,
    has_real_evolvers,
    has_real_non_trivial_evolvers,
    list_real_evolvers,
)


def _member(props):
    """Wrap a list-of-property-dicts as a `member` (component implementation)
    in the shape `_component_properties` expects."""
    return {'name': 'standard', 'properties': {'property': props}}


# ---------------------------------------------------------------------------
# has_real_evolvers
# ---------------------------------------------------------------------------

def test_has_real_evolvers_true_for_evolvable_non_virtual_property():
    member = _member([{
        'name':       'mass',
        'data':       {'isEvolvable': True, 'type': 'double', 'rank': 0},
        'attributes': {},
    }])
    assert has_real_evolvers(member) is True


def test_has_real_evolvers_false_for_virtual_evolvable_property():
    """`isVirtual=True` excludes the property from "real" evolvers (its
    derivative comes from the get function, not state)."""
    member = _member([{
        'name':       'mass',
        'data':       {'isEvolvable': True, 'type': 'double', 'rank': 0},
        'attributes': {'isVirtual': True},
    }])
    assert has_real_evolvers(member) is False


def test_has_real_evolvers_false_for_non_evolvable_property():
    member = _member([{
        'name':       'mass',
        'data':       {'isEvolvable': False, 'type': 'double', 'rank': 0},
        'attributes': {},
    }])
    assert has_real_evolvers(member) is False


def test_has_real_evolvers_false_for_member_with_no_properties():
    assert has_real_evolvers({'name': 'X'}) is False
    assert has_real_evolvers({'name': 'X', 'properties': {}}) is False


def test_has_real_evolvers_true_when_any_property_qualifies():
    """Mixed list: virtual + non-evolvable + real evolver → True."""
    member = _member([
        {'name': 'mass', 'data': {'isEvolvable': True}, 'attributes': {'isVirtual': True}},
        {'name': 'spin', 'data': {'isEvolvable': False}, 'attributes': {}},
        {'name': 'energy', 'data': {'isEvolvable': True, 'type': 'double', 'rank': 0}, 'attributes': {}},
    ])
    assert has_real_evolvers(member) is True


# ---------------------------------------------------------------------------
# has_real_non_trivial_evolvers — rank > 0 OR non-double
# ---------------------------------------------------------------------------

def test_non_trivial_evolvers_false_for_double_rank0():
    """Plain rank-0 double is the trivial case."""
    member = _member([{
        'name':       'mass',
        'data':       {'isEvolvable': True, 'type': 'double', 'rank': 0},
        'attributes': {},
    }])
    assert has_real_non_trivial_evolvers(member) is False


def test_non_trivial_evolvers_true_for_rank1_array():
    member = _member([{
        'name':       'positions',
        'data':       {'isEvolvable': True, 'type': 'double', 'rank': 1},
        'attributes': {},
    }])
    assert has_real_non_trivial_evolvers(member) is True


def test_non_trivial_evolvers_true_for_non_double_type():
    member = _member([{
        'name':       'kind',
        'data':       {'isEvolvable': True, 'type': 'integer', 'rank': 0},
        'attributes': {},
    }])
    assert has_real_non_trivial_evolvers(member) is True


def test_non_trivial_evolvers_skips_virtual_and_non_evolvable():
    """Virtual or non-evolvable properties don't trigger the predicate."""
    member = _member([
        {'name': 'p1', 'data': {'isEvolvable': True, 'type': 'integer'}, 'attributes': {'isVirtual': True}},
        {'name': 'p2', 'data': {'isEvolvable': False, 'type': 'integer'}, 'attributes': {}},
    ])
    assert has_real_non_trivial_evolvers(member) is False


# ---------------------------------------------------------------------------
# list_real_evolvers
# ---------------------------------------------------------------------------

def test_list_real_evolvers_returns_only_qualifying_properties():
    p_real    = {'name': 'energy', 'data': {'isEvolvable': True}, 'attributes': {}}
    p_virtual = {'name': 'mass',   'data': {'isEvolvable': True}, 'attributes': {'isVirtual': True}}
    p_static  = {'name': 'spin',   'data': {'isEvolvable': False}, 'attributes': {}}
    member = _member([p_real, p_virtual, p_static])
    out = list_real_evolvers(member)
    assert [p['name'] for p in out] == ['energy']


def test_list_real_evolvers_empty_for_no_properties():
    assert list_real_evolvers({'name': 'X'}) == []


# ---------------------------------------------------------------------------
# Implementation_Is_Active — emits one IsActive accessor per (class, member)
# ---------------------------------------------------------------------------

def test_is_active_emits_per_member_accessor():
    build = {
        'componentClasses': {
            'darkMatter': {
                'name':    'darkMatter',
                'members': [
                    {'name': 'standard'},
                    {'name': 'alternative'},
                ],
            },
        },
    }
    Implementation_Is_Active(build)
    bound = build['types']['nodeComponentDarkMatter']['boundFunctions']
    names = [b['name'] for b in bound]
    assert 'standardIsActive'    in names
    assert 'alternativeIsActive' in names


def test_is_active_function_body_reads_module_variable():
    """The generated function returns the matching `<Type>IsActiveValue`
    module variable (set elsewhere in the build)."""
    build = {
        'componentClasses': {
            'darkMatter': {
                'name':    'darkMatter',
                'members': [{'name': 'standard'}],
            },
        },
    }
    Implementation_Is_Active(build)
    desc = build['types']['nodeComponentDarkMatter']['boundFunctions'][0]['descriptor']
    assert 'nodeComponentDarkMatterStandardIsActive=nodeComponentDarkMatterStandardIsActiveValue' in desc['content']
    assert desc['type'] == 'logical'


def test_is_active_no_componentClasses_is_noop():
    build = {}
    Implementation_Is_Active(build)
    assert 'types' not in build


def test_is_active_member_with_no_members_list_is_noop():
    """A class with no `members` list contributes no IsActive functions
    but does NOT raise."""
    build = {
        'componentClasses': {
            'X': {'name': 'X'},   # no 'members' key
        },
    }
    Implementation_Is_Active(build)
    assert 'types' not in build
