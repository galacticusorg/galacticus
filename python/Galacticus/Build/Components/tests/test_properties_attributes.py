# Tests for `Galacticus.Build.Components.Properties.Attributes`.
#
# Single hook `Attributes_Match(build, class_dict)` that emits one
# `<class><Property>AttributeMatch` accessor per property of the class.
# The accessor returns the list of implementations matching a requested
# get/set/rate attribute requirement, used at runtime for filtering.

from Galacticus.Build.Components.Properties.Attributes import Attributes_Match


# ---------------------------------------------------------------------------
# Attributes_Match — emits per-property accessor functions
# ---------------------------------------------------------------------------

def test_attributes_match_emits_one_accessor_per_property():
    class_dict = {
        'name':    'darkMatter',
        'members': [
            {
                'name': 'standard',
                'properties': {'property': [
                    {'name': 'mass', 'attributes': {'isGettable': True}},
                    {'name': 'spin', 'attributes': {'isGettable': True, 'isSettable': True}},
                ]},
            },
        ],
    }
    build = {}
    Attributes_Match(build, class_dict)

    bound = build['types']['nodeComponentDarkMatter']['boundFunctions']
    names = [b['name'] for b in bound]
    assert 'massAttributeMatch' in names
    assert 'spinAttributeMatch' in names


def test_attributes_match_descriptor_has_expected_signature():
    class_dict = {
        'name':    'darkMatter',
        'members': [
            {
                'name': 'standard',
                'properties': {'property': [
                    {'name': 'mass', 'attributes': {'isGettable': True}},
                ]},
            },
        ],
    }
    build = {}
    Attributes_Match(build, class_dict)

    desc = build['types']['nodeComponentDarkMatter']['boundFunctions'][0]['descriptor']
    assert desc['name']    == 'darkMatterMassAttributeMatch'
    assert desc['type']    == 'type(varying_string), allocatable, dimension(:) => matches'
    assert desc['modules'] == ['ISO_Varying_String']
    # Signature accepts three optional logical filters.
    var_blocks = desc['variables']
    optional_block = next(v for v in var_blocks if 'optional' in (v.get('attributes') or []))
    assert set(optional_block['variables']) == {
        'requireSettable', 'requireGettable', 'requireEvolvable',
    }


def test_attributes_match_aggregates_property_across_members():
    """A property defined on TWO members of the same class produces ONE
    accessor whose body covers both members."""
    class_dict = {
        'name':    'darkMatter',
        'members': [
            {
                'name': 'standard',
                'properties': {'property': [
                    {'name': 'mass', 'attributes': {'isGettable': True}},
                ]},
            },
            {
                'name': 'alternative',
                'properties': {'property': [
                    {'name': 'mass', 'attributes': {'isGettable': True}},
                ]},
            },
        ],
    }
    build = {}
    Attributes_Match(build, class_dict)

    bound = build['types']['nodeComponentDarkMatter']['boundFunctions']
    # Just one accessor for `mass`, not two.
    mass_accessors = [b for b in bound if b['name'] == 'massAttributeMatch']
    assert len(mass_accessors) == 1
    body = mass_accessors[0]['descriptor']['content']
    # Body lists BOTH member names.
    assert "matches(size(matches))='standard'"    in body
    assert "matches(size(matches))='alternative'" in body


def test_attributes_match_walks_inherited_properties_via_extends():
    """An implementation that `extends` another inherits its parent's
    properties (via the `extends.implementation` chain).  The accessor
    should record the property against the child member's name."""
    parent_member = {
        'name': 'parent',
        'properties': {'property': [
            {'name': 'mass', 'attributes': {'isGettable': True}},
        ]},
    }
    child_member = {
        'name': 'child',
        'extends': {'implementation': parent_member},
        'properties': {'property': []},
    }
    class_dict = {
        'name':    'darkMatter',
        'members': [child_member],
    }
    build = {}
    Attributes_Match(build, class_dict)

    bound = build['types']['nodeComponentDarkMatter']['boundFunctions']
    names = [b['name'] for b in bound]
    assert 'massAttributeMatch' in names
    body = bound[0]['descriptor']['content']
    # The inherited property gets recorded under the CHILD member's name.
    assert "matches(size(matches))='child'" in body


def test_attributes_match_emits_filter_logic_for_requested_attributes():
    """A property that's gettable but not settable produces an
    `if (.not.requireSettableActual.and..not.requireEvolvableActual)` guard
    around the matches-append (because the property doesn't satisfy
    requestable settability/evolvability)."""
    class_dict = {
        'name':    'darkMatter',
        'members': [
            {
                'name': 'standard',
                'properties': {'property': [
                    {'name': 'mass', 'attributes': {'isGettable': True}},  # NOT settable, NOT evolvable
                ]},
            },
        ],
    }
    build = {}
    Attributes_Match(build, class_dict)
    body = build['types']['nodeComponentDarkMatter']['boundFunctions'][0]['descriptor']['content']
    # The body guards the matches-append on the missing capabilities.
    assert '.not.requireSettableActual' in body
    assert '.not.requireEvolvableActual' in body
    # A property with all three attributes wouldn't generate any guards.


def test_attributes_match_no_guards_for_fully_capable_property():
    """A property that has get/set/rate ALL = True produces an unguarded
    matches-append (every `require...Actual` requirement is satisfied)."""
    class_dict = {
        'name':    'darkMatter',
        'members': [
            {
                'name': 'standard',
                'properties': {'property': [{
                    'name':       'mass',
                    'attributes': {
                        'isGettable':  True,
                        'isSettable':  True,
                        'isEvolvable': True,
                    },
                }]},
            },
        ],
    }
    build = {}
    Attributes_Match(build, class_dict)
    body = build['types']['nodeComponentDarkMatter']['boundFunctions'][0]['descriptor']['content']
    # No `.not.require*Actual` guards in the per-member section (we still
    # have the prologue lines initialising the requireXxxActual variables).
    assert 'if (' not in body or '.not.require' not in body
    # The matches-append is unguarded, i.e. it appears at the top level
    # of the per-member section.
    assert "matches(size(matches))='standard'" in body
