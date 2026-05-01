# Tests for `Galacticus.Build.Components.Attributes`.
#
# Four pipeline-hook validators that the components-build runs in
# preValidate / default / postValidate phases:
#
#   Validate_Deferreds_Functionless  — forbid `xxxFunction` when xxx is deferred
#   Default_Functions                — fill in default get/set/rate function dicts
#   Validate_Boolean                 — coerce true/false strings to bools
#   Validate_Evolvable_Intrinsics    — forbid evolvable non-double intrinsics
#
# This module was the subject of commit `f83dc402` which converted its
# `sys.exit(...)` calls to `raise ValueError(...)`.  These tests pin both the
# happy paths and every documented validation error.

import pytest

from Galacticus.Build.Components.Attributes import (
    Validate_Deferreds_Functionless,
    Default_Functions,
    Validate_Boolean,
    Validate_Evolvable_Intrinsics,
)


def _build(class_name, impl_name, props):
    """Build a minimal `build` dict with one component carrying `props`."""
    return {
        'components': {
            f'{class_name}{impl_name}': {
                'class': class_name,
                'name':  impl_name,
                'properties': {'property': props},
            },
        },
    }


# ---------------------------------------------------------------------------
# Validate_Deferreds_Functionless
# ---------------------------------------------------------------------------

def test_deferreds_functionless_passes_when_no_deferred():
    """A property without `isDeferred` is unconstrained."""
    build = _build('darkMatter', 'standard', [
        {'name': 'mass', 'attributes': {}},
    ])
    Validate_Deferreds_Functionless(build)  # no exception expected


def test_deferreds_functionless_passes_when_deferred_but_no_function_supplied():
    """Deferred property with no compile-time function is the legitimate
    deferred case."""
    build = _build('darkMatter', 'standard', [
        {'name': 'mass', 'attributes': {'isDeferred': 'set'}},
    ])
    Validate_Deferreds_Functionless(build)


def test_deferreds_functionless_raises_when_function_supplied_for_deferred():
    """`isDeferred=set` AND a `setFunction` together is a contradiction."""
    build = _build('darkMatter', 'standard', [
        {'name': 'mass', 'attributes': {'isDeferred': 'set'}, 'setFunction': 'mySet'},
    ])
    with pytest.raises(ValueError, match=r"setFunction.*set.*deferred.*'mass'"):
        Validate_Deferreds_Functionless(build)


def test_deferreds_functionless_handles_colon_separated_methods():
    """`isDeferred='get:set'` checks each method independently."""
    build = _build('darkMatter', 'standard', [
        {'name': 'mass', 'attributes': {'isDeferred': 'get:set'}, 'getFunction': 'g'},
    ])
    with pytest.raises(ValueError, match=r"getFunction.*get.*deferred"):
        Validate_Deferreds_Functionless(build)


def test_deferreds_functionless_skips_components_without_properties():
    build = {'components': {'X': {'class': 'X', 'name': 'Y'}}}
    Validate_Deferreds_Functionless(build)


# ---------------------------------------------------------------------------
# Default_Functions
# ---------------------------------------------------------------------------

def test_default_functions_fills_in_rateFunction():
    """Missing `rateFunction` defaults to `<class><Impl><Prop>Rate`."""
    build = _build('darkMatter', 'standard', [
        {'name': 'mass'},
    ])
    Default_Functions(build)
    prop = build['components']['darkMatterstandard']['properties']['property'][0]
    assert prop['rateFunction'] == 'DarkMatterStandardMassRate'


def test_default_functions_synthesises_get_and_set():
    """Missing `getFunction` and `setFunction` get a synthesised content + build flag."""
    build = _build('darkMatter', 'standard', [
        {'name': 'mass'},
    ])
    Default_Functions(build)
    prop = build['components']['darkMatterstandard']['properties']['property'][0]
    assert prop['getFunction'] == {'content': 'darkMatterStandardMassGet', 'build': True}
    assert prop['setFunction'] == {'content': 'darkMatterStandardMassSet', 'build': True}


def test_default_functions_promotes_string_to_dict_with_build_false():
    """A user-supplied string `getFunction` becomes a dict with `build=False`
    (signalling the build system not to emit a body)."""
    build = _build('darkMatter', 'standard', [
        {'name': 'mass', 'getFunction': 'myCustomGet'},
    ])
    Default_Functions(build)
    prop = build['components']['darkMatterstandard']['properties']['property'][0]
    assert prop['getFunction'] == {'content': 'myCustomGet', 'build': False}


def test_default_functions_preserves_user_dict_and_sets_build_false():
    """A user-supplied dict gets a `build=False` flag added; existing fields
    survive."""
    build = _build('darkMatter', 'standard', [
        {'name': 'mass', 'setFunction': {'content': 'X', 'extra': 'keep'}},
    ])
    Default_Functions(build)
    prop = build['components']['darkMatterstandard']['properties']['property'][0]
    assert prop['setFunction'] == {'content': 'X', 'extra': 'keep', 'build': False}


# ---------------------------------------------------------------------------
# Validate_Boolean
# ---------------------------------------------------------------------------

def test_validate_boolean_coerces_true_string():
    build = _build('cls', 'impl', [
        {'name': 'p', 'attributes': {'isSettable': 'true'}},
    ])
    Validate_Boolean(build)
    prop = build['components']['clsimpl']['properties']['property'][0]
    assert prop['attributes']['isSettable'] is True


def test_validate_boolean_coerces_false_string():
    build = _build('cls', 'impl', [
        {'name': 'p', 'attributes': {'isGettable': 'false'}},
    ])
    Validate_Boolean(build)
    prop = build['components']['clsimpl']['properties']['property'][0]
    assert prop['attributes']['isGettable'] is False


def test_validate_boolean_passes_existing_bool_through():
    """If the value is already a Python bool (e.g. Default_Functions ran twice),
    coercion is idempotent."""
    build = _build('cls', 'impl', [
        {'name': 'p', 'attributes': {'isEvolvable': True}},
    ])
    Validate_Boolean(build)
    assert build['components']['clsimpl']['properties']['property'][0]['attributes']['isEvolvable'] is True


def test_validate_boolean_raises_on_typo():
    """`isSettable='ture'` (typo) must raise — the bug this guard catches."""
    build = _build('darkMatter', 'standard', [
        {'name': 'mass', 'attributes': {'isSettable': 'ture'}},
    ])
    with pytest.raises(ValueError, match=r"isSettable.*'mass'.*'true' or 'false'"):
        Validate_Boolean(build)


def test_validate_boolean_skips_when_attribute_missing():
    """An attribute that's simply absent is left alone (no default applied)."""
    build = _build('cls', 'impl', [
        {'name': 'p', 'attributes': {}},
    ])
    Validate_Boolean(build)


def test_validate_boolean_skips_when_attributes_not_dict():
    """A property with non-dict `attributes` is silently skipped."""
    build = _build('cls', 'impl', [
        {'name': 'p', 'attributes': None},
    ])
    Validate_Boolean(build)  # no exception


# ---------------------------------------------------------------------------
# Validate_Evolvable_Intrinsics
# ---------------------------------------------------------------------------

def test_evolvable_double_is_allowed():
    build = _build('darkMatter', 'standard', [
        {'name': 'mass', 'type': 'double', 'attributes': {'isEvolvable': True}},
    ])
    Validate_Evolvable_Intrinsics(build)  # no exception


def test_evolvable_non_intrinsic_object_is_allowed():
    """A non-intrinsic object type is unconstrained (it has its own evolution mechanism)."""
    build = _build('darkMatter', 'standard', [
        {'name': 'profile', 'type': 'darkMatterProfile', 'attributes': {'isEvolvable': True}},
    ])
    Validate_Evolvable_Intrinsics(build)


def test_evolvable_integer_intrinsic_raises():
    """An evolvable integer is forbidden — only `double` is supported."""
    build = _build('darkMatter', 'standard', [
        {'name': 'count', 'type': 'integer', 'attributes': {'isEvolvable': True}},
    ])
    with pytest.raises(ValueError, match=r"non-real intrinsic.*'count'.*can not be evolvable"):
        Validate_Evolvable_Intrinsics(build)


def test_evolvable_logical_intrinsic_raises():
    build = _build('cls', 'impl', [
        {'name': 'flag', 'type': 'logical', 'attributes': {'isEvolvable': True}},
    ])
    with pytest.raises(ValueError, match=r"non-real intrinsic"):
        Validate_Evolvable_Intrinsics(build)


def test_non_evolvable_non_double_intrinsic_passes():
    """A non-evolvable intrinsic of any type is fine."""
    build = _build('cls', 'impl', [
        {'name': 'count', 'type': 'integer', 'attributes': {'isEvolvable': False}},
        {'name': 'flag',  'type': 'logical', 'attributes': {}},  # no isEvolvable at all
    ])
    Validate_Evolvable_Intrinsics(build)
