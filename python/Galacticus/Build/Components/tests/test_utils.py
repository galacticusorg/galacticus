# Tests for `Galacticus.Build.Components.Utils`.
#
# This module is the central registry and helper-set for the components-build
# pipeline.  It owns:
#   - the `component_utils` hook registry (mutated at import time by every
#     sister module via `register(...)`),
#   - intrinsic / output-intrinsic type recognisers,
#   - the `offset_name` Fortran-symbol builder (two arities),
#   - the `pad_*` column-aligners (read module-level length-max globals),
#   - `argument_list` for extracting Fortran argument names,
#   - `apply_defaults` for recursive default-attribute population, and
#   - `Label_Lengths` for computing those length-max globals.
#
# It carried zero coverage despite being on the critical path of every
# components-build phase.  These tests pin the documented behaviour against
# the Perl reference so any accidental change to a frequently-called primitive
# is caught.

import pytest

from Galacticus.Build.Components import Utils
from Galacticus.Build.Components.Utils import (
    register,
    is_intrinsic,
    is_output_intrinsic,
    offset_name,
    pad_class,
    pad_implementation,
    pad_property_name,
    apply_defaults,
    argument_list,
    Label_Lengths,
    _as_list,
    _component_properties,
)


# ---------------------------------------------------------------------------
# is_intrinsic / is_output_intrinsic
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("type_name,expected", [
    ('integer',     True),
    ('longInteger', True),
    ('logical',     True),
    ('double',      True),
    ('void',        True),
    ('treeNode',    False),
    ('',            False),
    (None,          False),
])
def test_is_intrinsic(type_name, expected):
    assert is_intrinsic(type_name) is expected


@pytest.mark.parametrize("type_name,expected", [
    ('double',      True),
    ('integer',     True),
    ('longInteger', True),
    # `logical` is intrinsic but NOT directly outputtable.
    ('logical',     False),
    ('void',        False),
    ('treeNode',    False),
])
def test_is_output_intrinsic(type_name, expected):
    assert is_output_intrinsic(type_name) is expected


# ---------------------------------------------------------------------------
# register / component_utils registry
# ---------------------------------------------------------------------------

def test_register_appends_to_owner_phase_list():
    """Multiple registrations under the same owner+phase preserve registration
    order (Perl array semantics)."""
    fn1 = lambda build: None
    fn2 = lambda build: None
    register('test_owner_unique_xyz', 'test_phase', fn1)
    register('test_owner_unique_xyz', 'test_phase', fn2)
    hooks = Utils.component_utils['test_owner_unique_xyz']['test_phase']
    assert hooks[-2:] == [fn1, fn2]


# ---------------------------------------------------------------------------
# offset_name
# ---------------------------------------------------------------------------

def test_offset_name_three_arg_form():
    """offset_name('all', 'darkMatter', 'mass') -> 'offsetAllDarkMatterMass'."""
    assert offset_name('all', 'darkMatter', 'mass') == 'offsetAllDarkMatterMass'
    assert offset_name('active', 'spin', 'angularMomentum') == 'offsetAtvSpinAngularMomentum'
    assert offset_name('inactive', 'a', 'b') == 'offsetItvAB'


def test_offset_name_four_arg_form():
    """offset_name(status, class_dict, member_dict, prop_dict)."""
    out = offset_name(
        'all',
        {'name': 'darkMatter'},
        {'name': 'standard'},
        {'name': 'mass'},
    )
    assert out == 'offsetAllDarkMatterStandardMass'


def test_offset_name_invalid_status_raises():
    with pytest.raises(ValueError, match="unrecognized status"):
        offset_name('typo', 'cls', 'prop')


def test_offset_name_wrong_arity_raises():
    """Mirrors Perl's `die` on an unsupported signature."""
    with pytest.raises(TypeError, match="incorrect number of arguments"):
        offset_name('all')
    with pytest.raises(TypeError, match="incorrect number of arguments"):
        offset_name('all', 'a', 'b', 'c', 'd')


# ---------------------------------------------------------------------------
# argument_list
# ---------------------------------------------------------------------------

def test_argument_list_collects_intent_variables():
    desc = {
        'intrinsic':  'integer',
        'attributes': ['intent(in)'],
        'variables':  ['count', 'depth'],
    }
    assert argument_list(desc) == ['count', 'depth']


def test_argument_list_collects_procedure_pointer_arguments():
    """Procedure-pointer descriptor with `pointer` attribute counts as an
    argument."""
    desc = {
        'intrinsic':  'procedure',
        'attributes': ['pointer'],
        'variables':  ['callback'],
    }
    assert argument_list(desc) == ['callback']


def test_argument_list_collects_isArgument_procedure():
    """Procedure-pointer with `isArgument: True` flag also counts."""
    desc = {
        'intrinsic':  'procedure',
        'attributes': [],
        'variables':  ['fn'],
        'isArgument': True,
    }
    assert argument_list(desc) == ['fn']


def test_argument_list_collects_external():
    """External declarations are always arguments."""
    desc = {
        'intrinsic': 'external',
        'variables': ['ext_fn'],
    }
    assert argument_list(desc) == ['ext_fn']


def test_argument_list_skips_non_argument_locals():
    """A plain integer with no intent/procedure-pointer is NOT an argument."""
    desc = {
        'intrinsic':  'integer',
        'attributes': [],
        'variables':  ['local_counter'],
    }
    assert argument_list(desc) == []


def test_argument_list_handles_multiple_descriptors():
    """argument_list flattens across all argument-bearing descriptors."""
    args = argument_list(
        {'intrinsic': 'integer', 'attributes': ['intent(in   )'], 'variables': ['n']},
        {'intrinsic': 'integer', 'attributes': [],                'variables': ['local']},
        {'intrinsic': 'double',  'attributes': ['intent(out)'],   'variables': ['x', 'y']},
    )
    assert args == ['n', 'x', 'y']


def test_argument_list_skips_non_dict_inputs():
    """A None or string in the input list is silently skipped."""
    assert argument_list(None, 'invalid', 42) == []


# ---------------------------------------------------------------------------
# apply_defaults
# ---------------------------------------------------------------------------

def test_apply_defaults_scalar_setdefault():
    """`apply_defaults({'x': 1}, 'y', 7)` sets `y=7`; `apply_defaults({'x': 1},
    'x', 7)` leaves `x=1` unchanged (setdefault semantics)."""
    obj = {'x': 1}
    apply_defaults(obj, 'y', 7)
    assert obj == {'x': 1, 'y': 7}
    apply_defaults(obj, 'x', 7)
    assert obj == {'x': 1, 'y': 7}


def test_apply_defaults_booleanFalse_default_when_missing():
    obj = {}
    apply_defaults(obj, 'isVirtual', 'booleanFalse')
    assert obj == {'isVirtual': False}


def test_apply_defaults_booleanTrue_default_when_missing():
    obj = {}
    apply_defaults(obj, 'flag', 'booleanTrue')
    assert obj == {'flag': True}


def test_apply_defaults_boolean_coerces_existing_string():
    """If the attribute is already set, a `boolean*` default coerces the
    string to a bool."""
    obj_t = {'flag': 'true'}
    apply_defaults(obj_t, 'flag', 'booleanFalse')  # default ignored, but coerced
    assert obj_t == {'flag': True}

    obj_f = {'flag': 'false'}
    apply_defaults(obj_f, 'flag', 'booleanTrue')
    assert obj_f == {'flag': False}


def test_apply_defaults_nested_dict_walks_into_subobjects():
    """A dict default `{'k': v, ...}` recurses into `obj[name]`."""
    obj = {'props': {'a': {}}}
    apply_defaults(obj, 'props', {'a': {'kind': 'scalar'}})
    # The nested 'kind' default should land on obj['props']['a'].
    assert obj['props']['a']['kind'] == 'scalar'


def test_apply_defaults_ALL_pseudokey_walks_every_value():
    """`apply_defaults(obj, 'ALL', {...})` applies the defaults to every
    value in `obj`."""
    obj = {'a': {}, 'b': {}}
    apply_defaults(obj, 'ALL', {'kind': 'scalar'})
    assert obj == {'a': {'kind': 'scalar'}, 'b': {'kind': 'scalar'}}


# ---------------------------------------------------------------------------
# pad_* (depend on module-level globals — save/restore via fixture)
# ---------------------------------------------------------------------------

@pytest.fixture
def with_lengths():
    """Snapshot and restore the module-level *_length_max globals."""
    saved = (
        Utils.class_name_length_max,
        Utils.implementation_name_length_max,
        Utils.property_name_length_max,
    )
    try:
        Utils.class_name_length_max          = 10
        Utils.implementation_name_length_max = 8
        Utils.property_name_length_max       = 12
        yield
    finally:
        (Utils.class_name_length_max,
         Utils.implementation_name_length_max,
         Utils.property_name_length_max) = saved


def test_pad_class_pads_to_length_max(with_lengths):
    assert pad_class('foo') == 'foo' + ' ' * 7  # 10 - 3 = 7 spaces
    assert pad_class('darkMatter') == 'darkMatter'  # exact length, no padding


def test_pad_implementation_pads_to_length_max(with_lengths):
    assert pad_implementation('std') == 'std' + ' ' * 5


def test_pad_property_name_pads_to_length_max(with_lengths):
    assert pad_property_name('mass') == 'mass' + ' ' * 8


def test_pad_extra_pad_extends_length(with_lengths):
    """`extra_pad=(n, 0)` adds `n` to the pad length."""
    assert pad_class('foo', extra_pad=(2, 0)) == 'foo' + ' ' * 9  # 12-3 = 9


def test_pad_extra_pad_minimum_floor(with_lengths):
    """`extra_pad=(0, m)` enforces a minimum pad length of `m`."""
    assert pad_class('xy', extra_pad=(0, 20)) == 'xy' + ' ' * 18


# ---------------------------------------------------------------------------
# _as_list
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("input_,expected", [
    (None,           []),
    ([],             []),
    ([1, 2, 3],      [1, 2, 3]),
    ('hello',        ['hello']),
    (42,             [42]),
    ({'a': 1},       [{'a': 1}]),
])
def test_as_list_normalises_xmlin_idiom(input_, expected):
    assert _as_list(input_) == expected


# ---------------------------------------------------------------------------
# _component_properties
# ---------------------------------------------------------------------------

def test_component_properties_no_properties_returns_empty():
    assert _component_properties({'name': 'X'}) == []
    assert _component_properties({'name': 'X', 'properties': None}) == []
    assert _component_properties({'name': 'X', 'properties': {}}) == []


def test_component_properties_list_of_dicts_passthrough():
    component = {'properties': {'property': [
        {'name': 'mass'},
        {'name': 'spin'},
    ]}}
    out = _component_properties(component)
    assert [p['name'] for p in out] == ['mass', 'spin']


def test_component_properties_single_dict_with_name():
    """xml_to_dict collapse: a single property is a dict, not a list."""
    component = {'properties': {'property': {'name': 'mass'}}}
    assert _component_properties(component) == [{'name': 'mass'}]


def test_component_properties_keyAttr_dict_form():
    """KeyAttr-style: properties.property is `{name: subdict}`."""
    component = {'properties': {'property': {
        'mass': {'type': 'double'},
        'spin': {'type': 'double'},
    }}}
    out = _component_properties(component)
    types = sorted(p['type'] for p in out)
    assert types == ['double', 'double']
    assert len(out) == 2


# ---------------------------------------------------------------------------
# Label_Lengths
# ---------------------------------------------------------------------------

def test_Label_Lengths_computes_max_lengths_and_emits_variable():
    """Label_Lengths sets module globals and pushes a propertyNameLengthMax
    variable into build['variables']."""
    saved = (
        Utils.class_name_length_max,
        Utils.implementation_name_length_max,
        Utils.fully_qualified_name_length_max,
        Utils.property_name_length_max,
        Utils.implementation_property_name_length_max,
    )
    try:
        Utils.verbosity_level = 0  # suppress diagnostic prints
        build = {
            'components': {
                'darkMatterStandard': {
                    'class': 'darkMatter',
                    'name':  'standard',
                    'properties': {'property': [
                        {'name': 'angularMomentum'},
                        {'name': 'mass'},
                    ]},
                },
                'spinSimple': {
                    'class': 'spin',
                    'name':  'simple',
                    'properties': {'property': [{'name': 'value'}]},
                },
            },
        }
        Label_Lengths(build)
        assert Utils.class_name_length_max          == len('darkMatter')   # 10
        assert Utils.implementation_name_length_max == len('standard')     # 8
        assert Utils.fully_qualified_name_length_max == len('darkMatter') + len('standard')  # 18
        assert Utils.property_name_length_max       == len('angularMomentum')  # 15
        assert Utils.implementation_property_name_length_max == \
            len('darkMatter') + len('standard') + len('angularMomentum')  # 33

        # Build-emitted variable.
        emitted = build['variables']
        assert any(
            v.get('intrinsic') == 'integer'
            and any('propertyNameLengthMax=15' in s for s in v.get('variables', []))
            for v in emitted
        )
    finally:
        Utils.verbosity_level = 1
        (Utils.class_name_length_max,
         Utils.implementation_name_length_max,
         Utils.fully_qualified_name_length_max,
         Utils.property_name_length_max,
         Utils.implementation_property_name_length_max) = saved


def test_Label_Lengths_handles_empty_components():
    """No components → all length-maxes go to 0 and no exception."""
    saved = (
        Utils.class_name_length_max,
        Utils.implementation_name_length_max,
        Utils.fully_qualified_name_length_max,
        Utils.property_name_length_max,
    )
    try:
        Utils.verbosity_level = 0
        build = {'components': {}}
        Label_Lengths(build)
        assert Utils.class_name_length_max == 0
        assert Utils.implementation_name_length_max == 0
        assert Utils.fully_qualified_name_length_max == 0
        assert Utils.property_name_length_max == 0
    finally:
        Utils.verbosity_level = 1
        (Utils.class_name_length_max,
         Utils.implementation_name_length_max,
         Utils.fully_qualified_name_length_max,
         Utils.property_name_length_max) = saved
