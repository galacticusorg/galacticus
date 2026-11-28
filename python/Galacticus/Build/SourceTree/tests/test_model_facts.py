"""Tests for `Galacticus.Build.SourceTree.Process.ModelFacts`.

The directives are pure code generation, so the tests assert on the emitted
Fortran text:

  * `modelFact` emits a `modelFactProvide` call on the named parameters;
  * `modelRequire` emits a `modelFactRequire` call whose help text names the
    class and implementation providing the fact, and the XML that selects it;
  * a `modelRequire` naming a fact that nothing provides is a build-time error
    (the token vocabulary is checked by the build);
  * `modelExclusive` emits a `modelFactExclusive` call carrying its description
    and location, and a token which no other `modelExclusive` directive claims
    is a build-time error;
  * `modelConsistent` emits code passing the descriptor of an object (or a
    value) to the registry, declaring the variables it needs, and a token which
    no other `modelConsistent` directive declares is a build-time error;
  * malformed directives raise rather than emitting broken Fortran.
"""

import pytest

from Galacticus.Build.SourceTree.Process.ModelFacts import (
    help_text, process_model_consistent, process_model_exclusive, process_model_fact, process_model_require,
)
from Galacticus.Build.SourceTree import serialize


_PROVIDERS = {
    'satellite:tidalHeatingNormalized:evolved': [
        ('nodeOperator', 'satelliteTidalHeating', 'source/tidal.F90'),
    ],
    'component:shared': [
        ('nodeOperator', 'first',  'source/first.F90' ),
        (None          , None   ,  'source/disk.F90'  ),
    ],
}

_CLAIMANTS = {
    'spheroid:tidalHeating': 2,
    'component:lonely':      1,
}

_CONSISTENT = {
    'disk:galacticComponentTruncation': 4,
    'component:lonely':                 1,
}


def _build(directive_type, directive):
    """Minimal AST: a function unit holding a single directive node."""
    parent = {
        'type':       'function',
        'name':       'testConstructor',
        'firstChild': None,
        'sibling':    None,
        'parent':     None,
    }
    node = {
        'type':       directive_type,
        'directive':  dict(directive),
        'firstChild': None,
        'sibling':    None,
        'parent':     parent,
        'source':     '<test>',
        'line':       17,
    }
    parent['firstChild'] = node
    return parent


def _emitted(parent):
    """Concatenate the auto-generated code nodes hanging off `parent`."""
    text  = ""
    child = parent['firstChild']
    while child is not None:
        if child.get('type') == 'code':
            text += child.get('content') or ""
        child = child.get('sibling')
    return text


def _fact(directive):
    parent = _build('modelFact', directive)
    process_model_fact(parent, {})
    return _emitted(parent)


def _require(directive, providers=None):
    parent = _build('modelRequire', directive)
    process_model_require(parent, {'modelFactProviders': _PROVIDERS if providers is None else providers})
    return _emitted(parent)


def _exclusive(directive, claimants=None):
    parent = _build('modelExclusive', directive)
    process_model_exclusive(parent, {'modelExclusiveClaimants': _CLAIMANTS if claimants is None else claimants})
    return _emitted(parent)


def test_fact_default_source():
    code = _fact({'name': 'satellite:tidalHeatingNormalized:evolved'})
    assert code == "call parameters%modelFactProvide('satellite:tidalHeatingNormalized:evolved')\n"


def test_fact_explicit_source():
    code = _fact({'name': 'a:b', 'source': 'parametersRoot'})
    assert code.startswith("call parametersRoot%modelFactProvide('a:b')")


@pytest.mark.parametrize('name', ['', 'a b', 'a::b', ':a', 'a:', "a'b", '1a'])
def test_fact_rejects_bad_token(name):
    with pytest.raises(RuntimeError, match='not a valid token'):
        _fact({'name': name})


def test_fact_rejects_bad_source():
    with pytest.raises(RuntimeError, match='not a variable name'):
        _fact({'name': 'a:b', 'source': 'parameters; stop'})


def test_require_emits_call_with_help_and_location():
    code = _require({
        'name':   'satellite:tidalHeatingNormalized:evolved',
        'reason': "the disk's heating needs it",
    })
    assert code.startswith(
        "call parameters%modelFactRequire('satellite:tidalHeatingNormalized:evolved',")
    # Apostrophes in the reason are doubled for the Fortran literal.
    assert "'the disk''s heating needs it'" in code
    # Help names the provider and the XML that selects it (apostrophes doubled
    # inside the Fortran literal).
    assert "''satelliteTidalHeating'' implementation of the ''nodeOperator'' class" in code
    assert '<nodeOperator value="satelliteTidalHeating"/>' in code
    # The location expression names the enclosing unit.
    assert 'function:testConstructor' in code
    assert code.rstrip().endswith(')')


def test_require_unknown_fact_fails_the_build():
    with pytest.raises(RuntimeError, match="no `modelFact` directive.*'satellite:tidalHeating:evolvd'"):
        _require({'name': 'satellite:tidalHeating:evolvd', 'reason': 'x'})


def test_require_needs_reason():
    with pytest.raises(RuntimeError, match='`reason` attribute is required'):
        _require({'name': 'satellite:tidalHeatingNormalized:evolved'})


def test_help_lists_every_provider():
    text = help_text('component:shared', _PROVIDERS['component:shared'])
    assert text.startswith("the fact 'component:shared' is established by one of: ")
    assert '<nodeOperator value="first"/>' in text
    assert 'code in disk.F90' in text


def test_directive_processed_once():
    parent = _build('modelFact', {'name': 'a:b'})
    process_model_fact(parent, {})
    process_model_fact(parent, {})
    assert _emitted(parent).count('modelFactProvide') == 1


def test_exclusive_emits_call_with_description_and_location():
    code = _exclusive({
        'name':        'spheroid:tidalHeating',
        'description': "the spheroid's tidal heating",
    })
    assert code.startswith("call parameters%modelFactExclusive('spheroid:tidalHeating',")
    assert "'the spheroid''s tidal heating'" in code
    assert 'function:testConstructor' in code
    assert code.rstrip().endswith(')')


def test_exclusive_explicit_source():
    code = _exclusive({'name': 'spheroid:tidalHeating', 'description': 'x', 'source': 'subParameters'})
    assert code.startswith("call subParameters%modelFactExclusive('spheroid:tidalHeating',")


def test_exclusive_single_claimant_fails_the_build():
    with pytest.raises(RuntimeError, match="'component:lonely' is claimed by only one"):
        _exclusive({'name': 'component:lonely', 'description': 'x'})
    with pytest.raises(RuntimeError, match="'spheroid:tidalHeatng' is claimed by only one"):
        _exclusive({'name': 'spheroid:tidalHeatng', 'description': 'x'})


def test_exclusive_needs_description():
    with pytest.raises(RuntimeError, match='`description` attribute is required'):
        _exclusive({'name': 'spheroid:tidalHeating'})


@pytest.mark.parametrize('name', ['', 'a b', 'a::b'])
def test_exclusive_rejects_bad_token(name):
    with pytest.raises(RuntimeError, match='not a valid token'):
        _exclusive({'name': name, 'description': 'x'})


def _consistent(directive, claimants=None):
    parent = _build('modelConsistent', directive)
    process_model_consistent(parent, {'modelConsistentClaimants': _CONSISTENT if claimants is None else claimants})
    return parent


def test_consistent_object_passes_descriptor():
    parent = _consistent({
        'name':        'disk:galacticComponentTruncation',
        'object':      'galacticComponentTruncation_',
        'description': "the disk's truncation",
    })
    code = serialize(parent)
    # The descriptor variable and the module providing its type are declared.
    assert 'use :: Input_Parameters, only : inputParameters' in code
    assert 'type(inputParameters) :: modelConsistentDescriptor_' in code
    # The descriptor of the object is built, passed to the registry, and destroyed.
    assert 'modelConsistentDescriptor_=inputParameters()' in code
    assert 'call galacticComponentTruncation_%descriptor(modelConsistentDescriptor_,includeClass=.true.)' in code
    assert ("call parameters%modelFactConsistentDescriptor('disk:galacticComponentTruncation',"
            "modelConsistentDescriptor_,'the disk''s truncation',") in code
    assert 'function:testConstructor' in code
    assert code.index('%modelFactConsistentDescriptor(') < code.index('call modelConsistentDescriptor_%destroy()')


def test_consistent_value_assigned_before_call():
    code = serialize(_consistent({
        'name':        'disk:galacticComponentTruncation',
        'value':       'char(self%label)',
        'description': 'x',
        'source':      'subParameters',
    }))
    assert 'character(len=:), allocatable :: modelConsistentValue_' in code
    # The value is assigned to a variable, not passed as a function result.
    assert 'modelConsistentValue_=char(self%label)' in code
    assert "call subParameters%modelFactConsistent('disk:galacticComponentTruncation',modelConsistentValue_,'x'," in code


@pytest.mark.parametrize('directive, message', [
    ({'name': 'disk:galacticComponentTruncation', 'object': 'a_'},                                  '`description` attribute is required'),
    ({'name': 'disk:galacticComponentTruncation', 'description': 'x'},                              'exactly one of the `object` and `value`'),
    ({'name': 'disk:galacticComponentTruncation', 'description': 'x', 'object': 'a_', 'value': 'b'}, 'exactly one of the `object` and `value`'),
    ({'name': 'disk:galacticComponentTruncation', 'description': 'x', 'object': 'a_; stop'},        'is not a variable name'),
    ({'name': 'component:lonely',                 'description': 'x', 'object': 'a_'},              "'component:lonely' is declared by only one"),
    ({'name': 'disk:galacticComponentTrunction',  'description': 'x', 'object': 'a_'},              "'disk:galacticComponentTrunction' is declared by only one"),
])
def test_consistent_rejects_bad_directives(directive, message):
    with pytest.raises(RuntimeError, match=message):
        _consistent(directive)
