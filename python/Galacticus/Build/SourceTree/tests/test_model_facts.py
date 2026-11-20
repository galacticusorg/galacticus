"""Tests for `Galacticus.Build.SourceTree.Process.ModelFacts`.

The directives are pure code generation, so the tests assert on the emitted
Fortran text:

  * `modelFact` emits a `modelFactProvide` call on the named parameters;
  * `modelRequire` emits a `modelFactRequire` call whose help text names the
    class and implementation providing the fact, and the XML that selects it;
  * a `modelRequire` naming a fact that nothing provides is a build-time error
    (the token vocabulary is checked by the build);
  * malformed directives raise rather than emitting broken Fortran.
"""

import pytest

from Galacticus.Build.SourceTree.Process.ModelFacts import (
    help_text, process_model_fact, process_model_require,
)


_PROVIDERS = {
    'satellite:tidalHeatingNormalized:evolved': [
        ('nodeOperator', 'satelliteTidalHeating', 'source/tidal.F90'),
    ],
    'component:shared': [
        ('nodeOperator', 'first',  'source/first.F90' ),
        (None          , None   ,  'source/disk.F90'  ),
    ],
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
