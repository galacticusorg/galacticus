"""Tests for `Galacticus.Build.Components.Properties.Set`.

Three pipeline-hook functions emit Fortran type-bound procedures for
component property setters:

  Build_Class_Setters  — `<prop>IsSettable` Boolean_False stub on the class type
  Bind_Set_Functions   — bind a user-supplied set function to the impl type
  Build_Set_Functions  — auto-generate a set function (rank 0 / rank 1 forms)

Plus the helper `_is_deferred(prop, verb)` that all three rely on for
parsing the colon-separated `isDeferred` attribute string.

This file was the subject of commit `9f98deae` ("Repair Group B WARTs in
Properties/Set"); the audit specifically called it out as needing
regression coverage.
"""

from Galacticus.Build.Components.Properties.Set import (
    Build_Class_Setters,
    Bind_Set_Functions,
    _is_deferred,
)


# ---------------------------------------------------------------------------
# _is_deferred helper
# ---------------------------------------------------------------------------

def test_is_deferred_no_attributes_returns_false():
    assert _is_deferred({}, 'set') is False
    assert _is_deferred({'attributes': None}, 'set') is False
    assert _is_deferred({'attributes': {}}, 'set') is False


def test_is_deferred_no_isDeferred_field_returns_false():
    assert _is_deferred({'attributes': {'isSettable': True}}, 'set') is False


def test_is_deferred_single_method_match():
    assert _is_deferred({'attributes': {'isDeferred': 'set'}}, 'set') is True
    assert _is_deferred({'attributes': {'isDeferred': 'set'}}, 'get') is False


def test_is_deferred_colon_separated_methods():
    """`isDeferred='get:set:rate'` means all three are deferred."""
    prop = {'attributes': {'isDeferred': 'get:set:rate'}}
    assert _is_deferred(prop, 'get')  is True
    assert _is_deferred(prop, 'set')  is True
    assert _is_deferred(prop, 'rate') is True
    assert _is_deferred(prop, 'sub')  is False  # substring of 'set' must not match


def test_is_deferred_empty_string_returns_false():
    """An explicit empty `isDeferred=''` is treated as not-deferred."""
    assert _is_deferred({'attributes': {'isDeferred': ''}}, 'set') is False


# ---------------------------------------------------------------------------
# Build_Class_Setters
# ---------------------------------------------------------------------------

def _build(class_name, impl_name, props):
    return {
        'components': {
            f'{class_name}{impl_name}': {
                'class': class_name,
                'name':  impl_name,
                'properties': {'property': props},
            },
        },
    }


def test_build_class_setters_emits_one_stub_per_property():
    """One `<prop>IsSettable` Boolean_False stub per declared property."""
    build = _build('darkMatter', 'standard', [
        {'name': 'mass'},
        {'name': 'spin'},
    ])
    Build_Class_Setters(build)
    bound = build['types']['nodeComponentDarkMatter']['boundFunctions']
    names = [b['name'] for b in bound]
    assert 'massIsSettable' in names
    assert 'spinIsSettable' in names
    # Both stubs go to the SAME class-level type (not impl-level).
    assert 'nodeComponentDarkMatterStandard' not in build.get('types', {})


def test_build_class_setters_stub_uses_Boolean_False_function():
    build = _build('darkMatter', 'standard', [{'name': 'mass'}])
    Build_Class_Setters(build)
    stub = build['types']['nodeComponentDarkMatter']['boundFunctions'][0]
    assert stub['function']   == 'Boolean_False'
    assert stub['returnType'] == r"\logicalzero"
    assert stub['arguments']  == ""
    assert stub['type']       == 'procedure'
    assert stub['pass']       == 'nopass'


def test_build_class_setters_is_idempotent_per_property():
    """Running Build_Class_Setters twice must NOT emit duplicate stubs.
    Mirrors the dedup check at the top of the loop."""
    build = _build('darkMatter', 'standard', [{'name': 'mass'}])
    Build_Class_Setters(build)
    Build_Class_Setters(build)
    bound = build['types']['nodeComponentDarkMatter']['boundFunctions']
    names = [b['name'] for b in bound]
    assert names.count('massIsSettable') == 1


def test_build_class_setters_aggregates_across_implementations():
    """Two impls of the same class share the class-level type, so a property
    common to both produces a single stub."""
    build = {
        'components': {
            'darkMatterStandard': {
                'class': 'darkMatter', 'name': 'standard',
                'properties': {'property': [{'name': 'mass'}]},
            },
            'darkMatterAlt': {
                'class': 'darkMatter', 'name': 'alt',
                'properties': {'property': [{'name': 'mass'}]},
            },
        },
    }
    Build_Class_Setters(build)
    bound = build['types']['nodeComponentDarkMatter']['boundFunctions']
    names = [b['name'] for b in bound]
    assert names.count('massIsSettable') == 1


def test_build_class_setters_handles_component_without_properties():
    """A component with `properties` absent contributes no stubs and does not raise."""
    build = {'components': {'X': {'class': 'spin', 'name': 'simple'}}}
    Build_Class_Setters(build)


# ---------------------------------------------------------------------------
# Bind_Set_Functions
# ---------------------------------------------------------------------------

def test_bind_set_binds_user_function_to_impl_type():
    """isSettable + user-supplied (build=False) setFunction → bind to impl type."""
    build = {}
    class_dict = {'name': 'darkMatter'}
    member     = {'name': 'standard'}
    prop = {
        'name': 'mass',
        'attributes':  {'isSettable': True},
        'setFunction': {'content': 'myCustomSet', 'build': False},
    }
    Bind_Set_Functions(build, class_dict, member, prop)
    bound = build['types']['nodeComponentDarkMatterStandard']['boundFunctions']
    assert bound == [{
        'type':     'procedure',
        'name':     'massSet',
        'function': 'myCustomSet',
    }]


def test_bind_set_skips_when_not_settable():
    build = {}
    Bind_Set_Functions(
        build,
        {'name': 'darkMatter'}, {'name': 'standard'},
        {
            'name': 'mass',
            'attributes':  {'isSettable': False},
            'setFunction': {'content': 'X', 'build': False},
        },
    )
    assert build == {}


def test_bind_set_skips_when_setFunction_is_to_be_built():
    """`build=True` means the build system will emit the function body --
    don't bind a user-supplied function."""
    build = {}
    Bind_Set_Functions(
        build,
        {'name': 'darkMatter'}, {'name': 'standard'},
        {
            'name': 'mass',
            'attributes':  {'isSettable': True},
            'setFunction': {'content': 'autoSet', 'build': True},
        },
    )
    assert build == {}


def test_bind_set_skips_when_set_is_deferred():
    """Deferred-set properties handle binding via a separate deferred path,
    not this one."""
    build = {}
    Bind_Set_Functions(
        build,
        {'name': 'darkMatter'}, {'name': 'standard'},
        {
            'name': 'mass',
            'attributes': {
                'isSettable':  True,
                'isDeferred':  'set',
            },
            'setFunction': {'content': 'X', 'build': False},
        },
    )
    assert build == {}


def test_bind_set_skips_when_isDeferred_includes_set_in_chain():
    """`isDeferred='get:set'` should also skip set-binding."""
    build = {}
    Bind_Set_Functions(
        build,
        {'name': 'darkMatter'}, {'name': 'standard'},
        {
            'name': 'mass',
            'attributes': {
                'isSettable':  True,
                'isDeferred':  'get:set',
            },
            'setFunction': {'content': 'X', 'build': False},
        },
    )
    assert build == {}


def test_bind_set_does_not_skip_when_only_get_deferred():
    """`isDeferred='get'` does NOT affect the set binding."""
    build = {}
    Bind_Set_Functions(
        build,
        {'name': 'darkMatter'}, {'name': 'standard'},
        {
            'name': 'mass',
            'attributes': {
                'isSettable':  True,
                'isDeferred':  'get',
            },
            'setFunction': {'content': 'myCustomSet', 'build': False},
        },
    )
    assert 'nodeComponentDarkMatterStandard' in build['types']
