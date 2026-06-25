"""Tests for `Galacticus.Parameters.query` (inheritance resolution for consumers)."""

from Galacticus.Parameters.query import (
    inheritance_chain, enumeration_values, resolved_parameters,
)

CATALOG = {
    'enumerations': {'densityKind': ['critical', 'mean']},
    'implementations': {
        'baseSimple': {
            'parent': 'baseClass',
            'parameters': [
                {'name': 'alpha', 'type': 'real'},
                {'name': 'mode', 'type': 'string', 'enumeration': 'densityKind'},
            ],
        },
        'baseChild': {            # extends baseSimple; redeclares alpha
            'parent': 'baseSimple',
            'parameters': [
                {'name': 'beta', 'type': 'integer'},
                {'name': 'alpha', 'type': 'real'},
            ],
        },
    },
}


def test_inheritance_chain():
    assert inheritance_chain(CATALOG, 'baseChild') == ['baseChild', 'baseSimple']
    assert inheritance_chain(CATALOG, 'baseSimple') == ['baseSimple']
    assert inheritance_chain(CATALOG, 'missing') == []


def test_enumeration_values():
    assert enumeration_values(CATALOG, 'densityKind') == ['critical', 'mean']
    assert enumeration_values(CATALOG, None) is None
    assert enumeration_values(CATALOG, 'nope') is None


def test_resolved_parameters_inheritance_dedup_and_enum():
    params = resolved_parameters(CATALOG, 'baseChild')
    assert [p['name'] for p in params] == ['beta', 'alpha', 'mode']  # own first, dedup
    by_name = {p['name']: p for p in params}
    assert by_name['beta']['inheritedFrom'] is None
    assert by_name['alpha']['inheritedFrom'] is None        # child's redeclaration wins
    assert by_name['mode']['inheritedFrom'] == 'baseSimple'  # inherited
    assert by_name['mode']['allowedValues'] == ['critical', 'mean']
    assert by_name['alpha']['allowedValues'] is None
