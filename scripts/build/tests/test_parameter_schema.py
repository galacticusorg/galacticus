"""Tests for scripts/build/parameterSchema.py (parameter-file XSD generator)."""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
from parameterSchema import build_schema   # noqa: E402

CATALOG = {
    'functionClasses': {
        'accretionHalo':        {'implementations': ['coldMode', 'simple']},
        'criticalOverdensity':  {'implementations': ['fixed']},   # also a parameter -> collision
    },
    'enumerations': {'densityKind': ['critical', 'mean']},
    'implementations': {
        'accretionHaloSimple': {'parameters': [
            {'name': 'alpha', 'enumeration': None},
            {'name': 'densityType', 'enumeration': 'densityKind'},
        ]},
        # makes `criticalOverdensity` also a scalar parameter name
        'haloMassFunctionFixed': {'parameters': [
            {'name': 'criticalOverdensity', 'enumeration': None},
        ]},
    },
}


def test_selector_enum_emitted():
    xsd = build_schema(CATALOG)
    assert '<xs:element name="accretionHalo">' in xsd
    assert '<xs:enumeration value="simple"/>' in xsd
    assert '<xs:enumeration value="coldMode"/>' in xsd


def test_enum_parameter_emitted():
    xsd = build_schema(CATALOG)
    assert '<xs:element name="densityType">' in xsd
    assert '<xs:enumeration value="critical"/>' in xsd
    assert '<xs:enumeration value="mean"/>' in xsd


def test_collision_name_is_unconstrained():
    """A name that is both a functionClass base and a scalar parameter must not
    have its value enumerated (its value may be a label or an ordinary value)."""
    xsd = build_schema(CATALOG)
    assert '<xs:element name="criticalOverdensity">' in xsd
    # 'fixed' is criticalOverdensity's only implementation; it must NOT appear as
    # an enumeration (the element is generic).
    assert '<xs:enumeration value="fixed"/>' not in xsd


def test_tolerated_roots_declared():
    """Roots of non-parameter files that share parameter directories (changes,
    merger trees, ...) are declared so a directory-scoped editor association does
    not flag them."""
    xsd = build_schema(CATALOG)
    for root in ('changes', 'mergerTrees', 'parameterGrid', 'tree'):
        assert f'<xs:element name="{root}">' in xsd


def test_schema_compiles_and_validates():
    etree = pytest.importorskip('lxml.etree')
    schema = etree.XMLSchema(etree.fromstring(build_schema(CATALOG).encode()))
    assert schema.validate(etree.fromstring(
        b'<parameters><accretionHalo value="simple">'
        b'<densityType value="critical"/></accretionHalo></parameters>'))
    # bad selector / bad enum value are rejected
    assert not schema.validate(etree.fromstring(
        b'<parameters><accretionHalo value="bogus"/></parameters>'))
    assert not schema.validate(etree.fromstring(
        b'<parameters><densityType value="bogus"/></parameters>'))
    # collision name accepts an ordinary (numeric) value
    assert schema.validate(etree.fromstring(
        b'<parameters><criticalOverdensity value="1.686"/></parameters>'))
    # a non-parameter file root (e.g. a `<changes>` file co-located with
    # parameter files) is tolerated rather than rejected
    assert schema.validate(etree.fromstring(
        b'<changes><change type="update"/></changes>'))
