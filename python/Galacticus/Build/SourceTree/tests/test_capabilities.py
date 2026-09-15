"""Tests for the capability markers: `<requires>` and `<forwards>` on an
implementation's directive, the `requires` method generated from them, and the
check an `objectBuilder` directive emits for its `withholds` attribute.

Andrew Benson (2026).
"""

import pytest

from Galacticus.Build.Capabilities import parse_requires, parse_withholds
from Galacticus.Build.SourceTree.Process.FunctionClass import _build_requires_method
from Galacticus.Build.SourceTree.Process.ObjectBuilder import _withholds_checks


METHODS = {
    'differential': {'argument': [
        'double precision          , intent(in   )           :: time, mass',
        'type            (treeNode), intent(inout), optional :: node',
    ]},
    'integrated': {'argument': [
        'double precision          , intent(in   )                   :: time, massLow, massHigh',
        'type            (treeNode), intent(inout), target, optional :: node',
        'integer                   , intent(  out)        , optional :: status',
    ]},
}


def _classes(**extra):
    classes = {
        'massFunctionPressSchechter': {
            'extends': 'massFunctionClass',
            'requires': {'method': 'differential', 'argument': 'node'},
        },
        # Extends the above, so must inherit its requirement as well as adding its own.
        'massFunctionPressSchechterBinned': {
            'extends': 'massFunctionPressSchechter',
            'requires': [{'method': 'integrated', 'argument': 'status'}],
        },
        'massFunctionMultiplier': {
            'extends': 'massFunctionClass',
            'forwards': {'object': 'massFunction_'},
        },
        'massFunctionTinker': {'extends': 'massFunctionClass'},
    }
    classes.update(extra)
    return classes


def _requires_method(classes):
    methods = {name: dict(method) for name, method in METHODS.items()}
    _build_requires_method({'name': 'massFunction'}, classes, methods, 'location__')
    return methods['requires']


def _branch(code, type_name):
    """The body of the `class is` branch for `type_name`, or None."""
    marker = f"class is ({type_name})\n"
    if marker not in code:
        return None
    body = code.split(marker, 1)[1]
    return body.split("class is (", 1)[0].split("end select\nend select", 1)[0]


def test_parse_withholds():
    assert parse_withholds('differential:node, integrated:node') == \
        [('differential', 'node'), ('integrated', 'node')]
    assert parse_withholds(' value:mass\n timeOfCollapse:mass ') == \
        [('value', 'mass'), ('timeOfCollapse', 'mass')]
    assert parse_withholds(None) == []
    assert parse_withholds('') == []


@pytest.mark.parametrize('text', ['differential', 'differential:node:extra', ':node', 'differential:'])
def test_parse_withholds_rejects_malformed(text):
    with pytest.raises(ValueError):
        parse_withholds(text)


def test_parse_requires():
    assert parse_requires({}) == ([], [])
    assert parse_requires({'requires': {'method': 'differential', 'argument': 'node'},
                           'forwards': [{'object': 'a_'}, {'object': 'b_'}]}) == \
        ([('differential', 'node')], ['a_', 'b_'])
    with pytest.raises(ValueError):
        parse_requires({'requires': {'method': 'differential'}})
    with pytest.raises(ValueError):
        parse_requires({'forwards': {}})


def test_requires_method_signature():
    method = _requires_method(_classes())
    assert method['type'] == 'logical'
    assert method['recursive'] == 'yes'   # a forwarding implementation calls it on another object
    assert method['argument'] == ['character(len=*), intent(in   ) :: method, argument']


def test_requires_method_recognizes_only_optional_arguments():
    code = _requires_method(_classes())['code']
    recognized = code.split('case default', 1)[0]
    for pair in ('differential:node', 'integrated:node', 'integrated:status'):
        assert f"'{pair}'" in recognized, code
    for pair in ('differential:mass', 'differential:time', 'integrated:massLow'):
        assert f"'{pair}'" not in recognized, code
    assert 'call Error_Report(' in code.split('case default', 1)[1].split('end select', 1)[0]


def test_requires_method_branches():
    code = _requires_method(_classes())['code']
    own = _branch(code, 'massFunctionPressSchechter')
    assert "'differential:node'" in own and "'integrated:status'" not in own, code
    extended = _branch(code, 'massFunctionPressSchechterBinned')
    assert "'differential:node'" in extended and "'integrated:status'" in extended, code
    forwarding = _branch(code, 'massFunctionMultiplier')
    assert 'associated(self%massFunction_)' in forwarding, code
    assert 'self%massFunction_%requires(method,argument)' in forwarding, code
    # An implementation which declares nothing needs no branch: it answers false.
    assert _branch(code, 'massFunctionTinker') is None, code
    assert 'massFunctionRequires=.false.' in code


def test_requires_method_without_declarations_has_no_select_type():
    code = _requires_method({'massFunctionTinker': {'extends': 'massFunctionClass'}})['code']
    assert 'select type' not in code


@pytest.mark.parametrize('pair', [
    {'method': 'differential', 'argument': 'mass'},     # not optional
    {'method': 'derivative',   'argument': 'node'},     # no such method
])
def test_requires_method_rejects_a_pair_which_can_not_be_withheld(pair):
    with pytest.raises(RuntimeError):
        _requires_method(_classes(massFunctionBad={'extends': 'massFunctionClass', 'requires': pair}))


def test_withholds_checks():
    lines = _withholds_checks(
        {'name': 'haloMassFunction_', 'withholds': 'differential:node integrated:node'},
        'haloMassFunction', 'location__')
    assert lines.count('call Error_Report(') == 2
    assert "haloMassFunction_%requires('differential','node')" in lines
    assert "haloMassFunction_%requires('integrated','node')" in lines
    assert _withholds_checks({'name': 'haloMassFunction_'}, 'haloMassFunction', 'location__') == ''


def test_withholds_checks_rejects_malformed():
    with pytest.raises(RuntimeError):
        _withholds_checks({'name': 'x_', 'withholds': 'differential'}, 'x', 'location__')
