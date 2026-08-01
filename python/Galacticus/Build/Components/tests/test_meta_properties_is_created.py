"""Tests for `Class_Meta_Property_Is_Created`.

Registering a meta-property with `isCreator="no"` always yields a valid ID, but
the meta-property has storage only if some class registered it with
`isCreator="yes"`. The generated `...MetaPropertyCrtd` functions report whether
that is so, letting callers treat an uncreated meta-property as absent rather
than triggering the fatal `metaPropertyNoCreator` path taken by the accessors.
"""

from Galacticus.Build.Components.Classes.MetaProperties import (
    Class_Meta_Property_Is_Created,
    meta_property_types,
)


def _bound(build, type_name='nodeComponentDarkMatterProfile'):
    return build['types'][type_name]['boundFunctions']


def test_one_function_generated_per_meta_property_type():
    build = {}
    Class_Meta_Property_Is_Created(build, {'name': 'darkMatterProfile'})
    assert len(_bound(build)) == len(meta_property_types)


def test_bound_names_are_spelled_out_in_full():
    """The caller-facing type-bound name is `...MetaPropertyIsCreated`."""
    build = {}
    Class_Meta_Property_Is_Created(build, {'name': 'darkMatterProfile'})
    names = {entry['name'] for entry in _bound(build)}
    assert 'floatRank0MetaPropertyIsCreated'       in names
    assert 'longIntegerRank1MetaPropertyIsCreated' in names
    for entry in _bound(build):
        assert entry['type'] == 'procedure'


def test_procedure_names_stay_within_the_fortran_identifier_limit():
    """Fortran caps identifiers at 63 characters. The longest combination is
    `darkMatterProfile` + `longInteger`, which reaches 62 with the abbreviated
    `...MetaPropertyCrtd` suffix - spelling the suffix out in full would
    overflow, which is why the internal name is abbreviated."""
    build = {}
    Class_Meta_Property_Is_Created(build, {'name': 'darkMatterProfile'})
    for entry in _bound(build):
        assert len(entry['descriptor']['name']) <= 63, entry['descriptor']['name']


def test_active_class_reads_the_creator_flag_guarded_by_allocation():
    """For an active class the body consults `...MetaPropertyCreator`, but only
    after checking it is allocated - it is unallocated until at least one
    meta-property of that kind has been registered."""
    build = {'componentClassListActive': ['darkMatterProfile']}
    Class_Meta_Property_Is_Created(build, {'name': 'darkMatterProfile'})
    contents = {
        entry['descriptor']['name']: entry['descriptor']['content']
        for entry in _bound(build)
    }
    body = contents['nodeComponentDarkMatterProfileFloatRank0MetaPropertyCrtd']
    assert 'allocated(darkMatterProfileFloatRank0MetaPropertyCreator)' in body
    assert 'isCreated_=darkMatterProfileFloatRank0MetaPropertyCreator(metaPropertyID)' \
        in body
    assert 'isCreated_=.false.' in body
    # Must never take the fatal path used by the accessors.
    assert 'metaPropertyNoCreator' not in body


def test_inactive_class_always_reports_not_created():
    build = {}                       # no class listed as active
    Class_Meta_Property_Is_Created(build, {'name': 'darkMatterProfile'})
    for entry in _bound(build):
        body = entry['descriptor']['content']
        assert 'isCreated_=.false.' in body
        assert 'MetaPropertyCreator' not in body


def test_return_type_is_a_logical_result():
    build = {'componentClassListActive': ['darkMatterProfile']}
    Class_Meta_Property_Is_Created(build, {'name': 'darkMatterProfile'})
    for entry in _bound(build):
        assert entry['descriptor']['type'] == 'logical => isCreated_'


def test_arguments_are_self_and_the_meta_property_id():
    build = {'componentClassListActive': ['darkMatterProfile']}
    Class_Meta_Property_Is_Created(build, {'name': 'darkMatterProfile'})
    descriptor = _bound(build)[0]['descriptor']
    variables = [v for decl in descriptor['variables'] for v in decl['variables']]
    assert variables == ['self', 'metaPropertyID']
    # `self` is not referenced by the body, so it must be marked unused.
    assert '!$GLC attributes unused :: self' in descriptor['content']
