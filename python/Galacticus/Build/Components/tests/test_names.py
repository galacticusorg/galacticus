# Tests for `Galacticus.Build.Components.Classes.Names` and
# `Galacticus.Build.Components.Implementations.Names`.
#
# Each module owns one `<X>_Type` function that produces a Fortran type-name
# accessor:
#
#   Class_Type           — `nodeComponent<Class>Type`
#   Implementation_Type  — `nodeComponent<Class><Member>Type`
#
# Both register on the `classNames` / `implementationNames` hook owners.
# The functions mutate `build['types']` to append a `procedure :: type`
# binding whose body returns a static string label.
#
# Pure mutation; no I/O; small enough to exercise exhaustively.

from Galacticus.Build.Components.Classes.Names import Class_Type
from Galacticus.Build.Components.Implementations.Names import Implementation_Type


# ---------------------------------------------------------------------------
# Class_Type
# ---------------------------------------------------------------------------

def test_class_type_creates_type_in_build_types_dict():
    """`Class_Type({}, {'name': 'darkMatter'})` adds an entry under
    `build['types']['nodeComponentDarkMatter']` with one bound function."""
    build = {}
    Class_Type(build, {'name': 'darkMatter'})
    assert 'types' in build
    assert 'nodeComponentDarkMatter' in build['types']
    bound = build['types']['nodeComponentDarkMatter']['boundFunctions']
    assert len(bound) == 1


def test_class_type_bound_function_has_expected_shape():
    build = {}
    Class_Type(build, {'name': 'darkMatter'})
    entry = build['types']['nodeComponentDarkMatter']['boundFunctions'][0]
    assert entry['type'] == 'procedure'
    assert entry['name'] == 'type'
    desc = entry['descriptor']
    assert desc['name'] == 'nodeComponentDarkMatterType'
    assert desc['description'] == "Returns the type name for the darkMatter component class."


def test_class_type_body_returns_static_label():
    """The Fortran body should set `name='nodeComponent:<class>'` (static
    string, no dispatch needed)."""
    build = {}
    Class_Type(build, {'name': 'spin'})
    desc = build['types']['nodeComponentSpin']['boundFunctions'][0]['descriptor']
    assert "name='nodeComponent:spin'" in desc['content']


def test_class_type_lowercase_name_is_capitalised_in_type_name():
    """The class name appears in two cases: original-case in the body label,
    capitalised in the Fortran type identifier."""
    build = {}
    Class_Type(build, {'name': 'darkMatter'})
    desc = build['types']['nodeComponentDarkMatter']['boundFunctions'][0]['descriptor']
    assert desc['name'] == 'nodeComponentDarkMatterType'  # capitalised D
    assert "name='nodeComponent:darkMatter'" in desc['content']  # original case


def test_class_type_appends_to_existing_bound_functions():
    """If `boundFunctions` already has an entry, the new one is appended
    rather than overwriting."""
    build = {
        'types': {
            'nodeComponentDarkMatter': {
                'boundFunctions': [
                    {'type': 'procedure', 'name': 'preExisting'},
                ],
            },
        },
    }
    Class_Type(build, {'name': 'darkMatter'})
    bound = build['types']['nodeComponentDarkMatter']['boundFunctions']
    assert [b['name'] for b in bound] == ['preExisting', 'type']


def test_class_type_uses_modern_string_module():
    """Description binding declares `ISO_Varying_String` import."""
    build = {}
    Class_Type(build, {'name': 'X'})
    desc = build['types']['nodeComponentX']['boundFunctions'][0]['descriptor']
    assert desc['modules'] == ['ISO_Varying_String']


# ---------------------------------------------------------------------------
# Implementation_Type
# ---------------------------------------------------------------------------

def test_implementation_type_creates_per_implementation_type():
    """Implementation type names concatenate `<Class><Member>` with both
    capitalised."""
    build = {}
    Implementation_Type(
        build,
        {'name': 'darkMatter'},
        {'name': 'standard'},
    )
    assert 'nodeComponentDarkMatterStandard' in build['types']


def test_implementation_type_body_includes_class_and_member_in_label():
    build = {}
    Implementation_Type(
        build,
        {'name': 'darkMatter'},
        {'name': 'standard'},
    )
    desc = build['types']['nodeComponentDarkMatterStandard']['boundFunctions'][0]['descriptor']
    assert "name='nodeComponent:darkMatter:standard'" in desc['content']


def test_implementation_type_description_mentions_both_levels():
    build = {}
    Implementation_Type(
        build,
        {'name': 'spin'},
        {'name': 'simple'},
    )
    desc = build['types']['nodeComponentSpinSimple']['boundFunctions'][0]['descriptor']
    assert "Returns the type name for the simple implementation of the spin component class." == desc['description']


def test_implementation_type_appends_to_existing_bound_functions():
    build = {
        'types': {
            'nodeComponentSpinSimple': {
                'boundFunctions': [
                    {'type': 'procedure', 'name': 'somethingElse'},
                ],
            },
        },
    }
    Implementation_Type(build, {'name': 'spin'}, {'name': 'simple'})
    bound = build['types']['nodeComponentSpinSimple']['boundFunctions']
    assert [b['name'] for b in bound] == ['somethingElse', 'type']
