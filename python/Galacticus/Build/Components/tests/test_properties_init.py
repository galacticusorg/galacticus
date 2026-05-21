"""Tests for `Galacticus.Build.Components.Properties` (the package's
`__init__.py` validators and default-population hooks).

Three of the seven hooks here are pure validators that the audit
specifically called out as needing regression coverage -- one of them
(`Property_Output_Validate`) was the subject of fix commit `83da4762`
("Repair Property_Output_Validate rank>1 check").
"""

import pytest

from Galacticus.Build.Components import Utils as _Utils
from Galacticus.Build.Components.Properties import (
    Property_Defaults,
    Data_Validate,
    Property_Output_Validate,
    Construct_Data,
)


def _build(props, *, class_name='darkMatter', impl_name='standard'):
    """Build with `props` in list-of-dicts shape (the form
    `_component_properties` handles).  Used by Data_Validate /
    Property_Output_Validate tests."""
    return {
        'components': {
            f'{class_name}{impl_name}': {
                'class': class_name,
                'name':  impl_name,
                'properties': {'property': props},
            },
        },
    }


def _build_keyed(props_keyed, *, class_name='darkMatter', impl_name='standard'):
    """Build with `props_keyed` in KeyAttr-keyed shape `{name: dict}`.
    `Property_Defaults` walks this form via `apply_defaults`."""
    return {
        'components': {
            f'{class_name}{impl_name}': {
                'class': class_name,
                'name':  impl_name,
                'properties': {'property': props_keyed},
            },
        },
    }


# ---------------------------------------------------------------------------
# Property_Defaults
# ---------------------------------------------------------------------------

def test_property_defaults_fills_in_all_four_boolean_attributes():
    """Property_Defaults expects the KeyAttr-style XML shape
    `properties.property = {name: dict}` with each property carrying an
    existing `attributes` dict; missing booleans then default to False."""
    build = _build_keyed({'mass': {'attributes': {}}})
    Property_Defaults(build)
    attrs = build['components']['darkMatterstandard']['properties']['property']['mass']['attributes']
    assert attrs['isVirtual']      is False
    assert attrs['createIfNeeded'] is False
    assert attrs['isNonNegative']  is False
    # `isDeferred` is a colon-separated method list, not a bool.
    assert attrs['isDeferred']     == 'false'


def test_property_defaults_does_not_override_existing_attributes():
    """Already-set attributes are coerced to bool but NOT overwritten."""
    build = _build_keyed({
        'mass': {'attributes': {'isVirtual': 'true', 'isNonNegative': 'true'}},
    })
    Property_Defaults(build)
    attrs = build['components']['darkMatterstandard']['properties']['property']['mass']['attributes']
    # apply_defaults' booleanFalse default coerces existing 'true' to True.
    assert attrs['isVirtual']     is True
    assert attrs['isNonNegative'] is True
    # Missing attributes still get the default.
    assert attrs['createIfNeeded'] is False


def test_property_defaults_no_components_is_noop():
    """No components → silently does nothing."""
    build = {'components': {}}
    Property_Defaults(build)


def test_property_defaults_handles_multiple_properties():
    build = _build_keyed({
        'mass': {'attributes': {}},
        'spin': {'attributes': {}},
    })
    Property_Defaults(build)
    props_keyed = build['components']['darkMatterstandard']['properties']['property']
    for prop in props_keyed.values():
        assert prop['attributes']['isVirtual'] is False


# ---------------------------------------------------------------------------
# Data_Validate
# ---------------------------------------------------------------------------

def test_data_validate_passes_when_type_and_rank_present():
    build = _build([{'name': 'mass', 'type': 'double', 'rank': 0}])
    Data_Validate(build)


def test_data_validate_raises_on_missing_type():
    build = _build([{'name': 'mass', 'rank': 0}])
    with pytest.raises(ValueError, match=r"no type.*'mass'.*standard.*darkMatter"):
        Data_Validate(build)


def test_data_validate_raises_on_missing_rank():
    build = _build([{'name': 'mass', 'type': 'double'}])
    with pytest.raises(ValueError, match=r"no rank.*'mass'"):
        Data_Validate(build)


def test_data_validate_forbids_createIfNeeded_on_rank_gt_0():
    """Rank > 0 + createIfNeeded is the explicit unsupported case."""
    build = _build([{
        'name': 'positions', 'type': 'double', 'rank': 1,
        'attributes': {'createIfNeeded': True},
    }])
    with pytest.raises(ValueError, match=r"auto-creation of rank > 0"):
        Data_Validate(build)


def test_data_validate_allows_createIfNeeded_on_rank_0():
    """Rank-0 properties can use createIfNeeded freely."""
    build = _build([{
        'name': 'mass', 'type': 'double', 'rank': 0,
        'attributes': {'createIfNeeded': True},
    }])
    Data_Validate(build)  # no exception


# ---------------------------------------------------------------------------
# Property_Output_Validate (recent bug-fix area: commit 83da4762)
# ---------------------------------------------------------------------------

def test_output_validate_skips_when_no_output_block():
    """Properties without an `output` dict are silently skipped."""
    build = _build([{'name': 'mass', 'attributes': {}}])
    Property_Output_Validate(build)


def test_output_validate_raises_for_virtual_non_gettable():
    """A virtual property without isGettable cannot be output."""
    build = _build([{
        'name': 'mass',
        'attributes': {'isVirtual': True, 'isGettable': False},
        'output': {},
    }])
    with pytest.raises(ValueError, match=r"non-gettable, virtual.*can not be output"):
        Property_Output_Validate(build)


def test_output_validate_raises_for_rank_gt_1_arrays():
    """The recent fix (commit 83da4762): `data.rank > 1` is the unsupported
    case; `rank == 1` must still be allowed."""
    build = _build([{
        'name':       'tensor',
        'attributes': {},
        'data':       {'rank': 2},
        'output':     {'labels': '[a,b]'},
    }])
    with pytest.raises(ValueError, match=r"rank>1 arrays.*is not supported"):
        Property_Output_Validate(build)


def test_output_validate_allows_rank_1_with_labels():
    """rank == 1 with parseable labels (`[a,b]`) is the standard case."""
    build = _build([{
        'name':       'positions',
        'attributes': {},
        'data':       {'rank': 1},
        'output':     {'labels': '[x,y,z]'},
    }])
    Property_Output_Validate(build)  # no exception


def test_output_validate_raises_when_rank_gt_0_lacks_labels():
    """rank > 0 without `labels` in the output block is forbidden."""
    build = _build([{
        'name':       'positions',
        'attributes': {},
        'data':       {'rank': 1},
        'output':     {},
    }])
    with pytest.raises(ValueError, match=r"requires a labels attribute"):
        Property_Output_Validate(build)


def test_output_validate_strips_whitespace_from_labels_and_modules():
    """After validation, whitespace inside `labels` and `modules` is stripped."""
    build = _build([{
        'name':       'positions',
        'attributes': {},
        'data':       {'rank': 1},
        'output':     {
            'labels':  '[ x , y , z ]',
            'modules': 'mod_a , mod_b',
        },
    }])
    Property_Output_Validate(build)
    out = build['components']['darkMatterstandard']['properties']['property'][0]['output']
    assert out['labels']  == '[x,y,z]'
    assert out['modules'] == 'mod_a,mod_b'


def test_output_validate_rank_0_allowed_without_labels():
    """Rank-0 (scalar) outputs don't need a labels attribute."""
    build = _build([{
        'name':       'mass',
        'attributes': {},
        'data':       {'rank': 0},
        'output':     {},
    }])
    Property_Output_Validate(build)


# ---------------------------------------------------------------------------
# Construct_Data
#
# Construct_Data is large and validates parent-implementation chains, but
# it ALSO has a small but mistake-prone responsibility: maintaining the
# `Components.Utils.linked_data_name_length_max` global as it visits each
# property.  That global is referenced as `_Utils.linked_data_name_length_max`
# -- and an earlier refactor accidentally dropped the `from ... import Utils
# as _Utils` import while leaving the references in place.  These tests
# exercise that codepath end-to-end so any future drop of the import (or any
# rename of the global) breaks here rather than at runtime in the build.
# ---------------------------------------------------------------------------

def test_construct_data_seeds_property_data_dict():
    """Construct_Data populates each property's `data` sub-dict from its
    type / rank / isEvolvable attribute."""
    build = _build([{
        'name': 'mass', 'type': 'double', 'rank': 0,
        'attributes': {'isEvolvable': True},
    }])
    Construct_Data(build)
    prop = build['components']['darkMatterstandard']['properties']['property'][0]
    assert prop['data'] == {'type': 'double', 'rank': 0, 'isEvolvable': True}
    assert prop['definedInParent'] is False


def test_construct_data_creates_linked_data_entry_in_component_content():
    """Each non-virtual, non-inherited property's `<name>Data` becomes a
    key under `component['content']['data']`."""
    build = _build([{
        'name': 'mass', 'type': 'double', 'rank': 0,
        'attributes': {'isEvolvable': True, 'isVirtual': False},
    }])
    Construct_Data(build)
    component = build['components']['darkMatterstandard']
    assert 'massData' in component['content']['data']


def test_construct_data_updates_linked_data_name_length_max():
    """The linked-data name's length is tracked on
    `Components.Utils.linked_data_name_length_max` so downstream column-
    aligned output uses a consistent width.  Regression test for the
    `_Utils.linked_data_name_length_max` reference in Construct_Data
    (its import was accidentally dropped during the logging migration)."""
    saved = _Utils.linked_data_name_length_max
    try:
        _Utils.linked_data_name_length_max = 0
        build = _build([{
            'name': 'angularMomentum', 'type': 'double', 'rank': 0,
            'attributes': {'isEvolvable': True, 'isVirtual': False},
        }])
        Construct_Data(build)
        # `angularMomentumData` is 19 characters.
        assert _Utils.linked_data_name_length_max == len('angularMomentumData')
    finally:
        _Utils.linked_data_name_length_max = saved


def test_construct_data_skips_linked_data_for_virtual_properties():
    """Virtual properties don't get a linked-data entry (their values come
    from the get function, not state)."""
    build = _build([{
        'name': 'massVirtual', 'type': 'double', 'rank': 0,
        'attributes': {'isEvolvable': False, 'isVirtual': True},
    }])
    Construct_Data(build)
    component = build['components']['darkMatterstandard']
    # `massVirtualData` should NOT be in the linked-data registry.
    assert 'massVirtualData' not in component.get('content', {}).get('data', {})
