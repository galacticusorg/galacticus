# Tests for `Galacticus.Build.StateStorables` shape-bridging helpers.
#
# `stateStorables.xml` is read back through `xml_to_dict()` and arrives in
# any of five different shapes depending on which writer produced it
# (XML::Simple ForceArray vs KeyAttr vs current Python writer; single vs
# multiple entries).  These functions normalise all five shapes to a
# consistent list-of-dicts representation; pinning that down with tests
# is high-leverage because the function-class generation loop iterates
# over their output.

from Galacticus.Build.StateStorables import (
    function_class_entries,
    function_class_names,
    function_class_module_map,
    function_class_entry,
    function_class_instances,
    event_hook_static_names,
)


# ---------------------------------------------------------------------------
# function_class_entries — five accepted XML shapes
# ---------------------------------------------------------------------------

def test_entries_none_returns_empty():
    assert function_class_entries(None) == []
    assert function_class_entries({}) == []


def test_entries_missing_functionClasses_key_returns_empty():
    """A `state_storables` dict with no `functionClasses` key is valid input
    (e.g. a build that registers no function classes)."""
    assert function_class_entries({'somethingElse': 1}) == []


def test_entries_list_of_dicts_passthrough():
    """Current Python writer with multiple entries: list of dicts."""
    state = {'functionClasses': [
        {'name': 'massDistributionClass', 'module': 'Mass_Distributions'},
        {'name': 'cosmologyParameters',   'module': 'Cosmology_Parameters'},
    ]}
    out = function_class_entries(state)
    assert len(out) == 2
    assert out[0]['name'] == 'massDistributionClass'
    assert out[1]['module'] == 'Cosmology_Parameters'


def test_entries_list_filters_non_dict_elements():
    """Defensive: stray non-dict items in the list are silently skipped
    rather than causing a TypeError downstream."""
    state = {'functionClasses': [
        {'name': 'X', 'module': 'M_X'},
        'invalid_string_entry',
        {'name': 'Y', 'module': 'M_Y'},
    ]}
    out = function_class_entries(state)
    assert [e['name'] for e in out] == ['X', 'Y']


def test_entries_single_dict_with_name_is_singleton():
    """Current writer with exactly one entry: xml_to_dict collapses
    `<functionClasses>` to a single dict instead of a list."""
    state = {'functionClasses': {'name': 'oneClass', 'module': 'OneMod'}}
    out = function_class_entries(state)
    assert out == [{'name': 'oneClass', 'module': 'OneMod'}]


def test_entries_force_array_wrapper_with_list():
    """XML::Simple ForceArray=>['functionClass'] writer:
    `<functionClasses><functionClass>...</functionClass>...</functionClasses>`."""
    state = {'functionClasses': {'functionClass': [
        {'name': 'A', 'module': 'M_A'},
        {'name': 'B', 'module': 'M_B'},
    ]}}
    out = function_class_entries(state)
    assert [e['name'] for e in out] == ['A', 'B']


def test_entries_force_array_wrapper_with_single_dict():
    """Force-array variant whose single entry was collapsed to a bare dict
    by xml_to_dict — must still be returned as a singleton list."""
    state = {'functionClasses': {'functionClass': {'name': 'solo', 'module': 'M'}}}
    out = function_class_entries(state)
    assert out == [{'name': 'solo', 'module': 'M'}]


def test_entries_keyAttr_mapping_synthesises_name_from_key():
    """XML::Simple KeyAttr=>['name'] writer: `{'X': {'module': 'M_X'}}`.
    The class name lives in the *key*, not in a `name` field; the helper
    must synthesise `name` from the key."""
    state = {'functionClasses': {
        'classA': {'module': 'M_A'},
        'classB': {'module': 'M_B'},
    }}
    out = function_class_entries(state)
    by_name = {e['name']: e for e in out}
    assert set(by_name) == {'classA', 'classB'}
    assert by_name['classA']['module'] == 'M_A'


def test_entries_keyAttr_does_not_overwrite_explicit_name():
    """If the inner dict already has a `name` field, it must take precedence
    over the synthesised key (the explicit name is what callers wrote)."""
    state = {'functionClasses': {
        'keyName': {'name': 'explicitName', 'module': 'M'},
    }}
    out = function_class_entries(state)
    assert out == [{'name': 'explicitName', 'module': 'M'}]


def test_entries_unknown_shape_returns_empty():
    """Anything that isn't list/dict/None (e.g. a stray string) is returned
    as empty rather than crashing."""
    state = {'functionClasses': 'unexpected_string'}
    assert function_class_entries(state) == []


# ---------------------------------------------------------------------------
# function_class_names / module_map / entry
# ---------------------------------------------------------------------------

def test_names_returns_set_of_class_names():
    state = {'functionClasses': [
        {'name': 'A', 'module': 'M_A'},
        {'name': 'B', 'module': 'M_B'},
    ]}
    assert function_class_names(state) == {'A', 'B'}


def test_names_skips_entries_missing_a_name():
    """An entry without `name` is silently skipped — better than raising,
    since legacy XMLs occasionally contain placeholder entries."""
    state = {'functionClasses': [
        {'name': 'A', 'module': 'M_A'},
        {'module': 'orphan'},
    ]}
    assert function_class_names(state) == {'A'}


def test_module_map_returns_name_to_module_dict():
    state = {'functionClasses': [
        {'name': 'A', 'module': 'M_A'},
        {'name': 'B', 'module': 'M_B'},
    ]}
    assert function_class_module_map(state) == {'A': 'M_A', 'B': 'M_B'}


def test_module_map_handles_missing_module_field():
    """An entry with no `module` field maps to None — callers can detect
    and skip rather than crashing on KeyError."""
    state = {'functionClasses': [{'name': 'A'}]}
    assert function_class_module_map(state) == {'A': None}


def test_entry_returns_match_when_present():
    state = {'functionClasses': [
        {'name': 'A', 'module': 'M_A'},
        {'name': 'B', 'module': 'M_B'},
    ]}
    assert function_class_entry(state, 'B') == {'name': 'B', 'module': 'M_B'}


def test_entry_returns_None_when_absent():
    state = {'functionClasses': [{'name': 'A', 'module': 'M_A'}]}
    assert function_class_entry(state, 'missing') is None


# ---------------------------------------------------------------------------
# function_class_instances / event_hook_static_names
# ---------------------------------------------------------------------------

def test_instances_none_returns_empty():
    assert function_class_instances(None) == []
    assert function_class_instances({}) == []


def test_instances_string_form_returns_singleton():
    """A single instance written as bare text in the XML."""
    assert function_class_instances({'functionClassInstances': 'soloInstance'}) == ['soloInstance']


def test_instances_empty_string_returns_empty():
    assert function_class_instances({'functionClassInstances': ''}) == []


def test_instances_list_of_strings():
    state = {'functionClassInstances': ['inst1', 'inst2', 'inst3']}
    assert function_class_instances(state) == ['inst1', 'inst2', 'inst3']


def test_instances_list_with_dict_entries():
    """Mixed list: instance written as dict gets its `name` (or `content`)
    field extracted; bare strings pass through."""
    state = {'functionClassInstances': [
        'plainInstance',
        {'name': 'fromName'},
        {'content': 'fromContent'},
    ]}
    assert function_class_instances(state) == ['plainInstance', 'fromName', 'fromContent']


def test_instances_single_dict_form():
    """A single instance written as a dict (xml_to_dict collapse)."""
    assert function_class_instances({'functionClassInstances': {'name': 'solo'}}) == ['solo']
    assert function_class_instances({'functionClassInstances': {'content': 'tx'}}) == ['tx']


def test_instances_dict_with_neither_name_nor_content_returns_empty():
    """A dict instance with neither `name` nor `content` is silently
    skipped (reflects writer-side defects)."""
    assert function_class_instances({'functionClassInstances': {'other': 'x'}}) == []


def test_event_hook_static_names_uses_same_logic():
    """`event_hook_static_names` mirrors `function_class_instances` for the
    `eventHookStatics` key — pin both shapes work."""
    state = {'eventHookStatics': ['hookA', {'name': 'hookB'}]}
    assert event_hook_static_names(state) == ['hookA', 'hookB']


def test_event_hook_static_names_none_returns_empty():
    assert event_hook_static_names(None) == []
    assert event_hook_static_names({}) == []
