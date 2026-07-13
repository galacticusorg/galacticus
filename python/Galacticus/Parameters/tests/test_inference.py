"""Tests for `Galacticus.Parameters.inference`.

Only ~1.5% of `<inputParameter>` directives declare an explicit `<type>`, so
type inference from the default-value literal and the variable's declaration is
the workhorse.  These tests pin the three inference sources and their
provenance, plus the Fortran-literal and declaration classifiers they rely on.
"""

import pytest

from Galacticus.Parameters.inference import (
    classify_literal, canonical_from_declaration, normalize_explicit_type,
    infer_parameter_type,
    TYPE_BOOLEAN, TYPE_INTEGER, TYPE_REAL, TYPE_STRING, TYPE_OBJECT, TYPE_UNKNOWN,
    PROV_EXPLICIT, PROV_DEFAULT, PROV_DECLARATION, PROV_UNKNOWN,
)


@pytest.mark.parametrize("literal, expected_type, expected_kind", [
    ('.true.',                  TYPE_BOOLEAN, None),
    ('.FALSE.',                 TYPE_BOOLEAN, None),
    ("var_str('none')",         TYPE_STRING,  None),
    ("'a string'",              TYPE_STRING,  None),
    ('"a string"',              TYPE_STRING,  None),
    ('42',                      TYPE_INTEGER, None),
    ('-7',                      TYPE_INTEGER, None),
    ('0_c_long',                TYPE_INTEGER, 'long'),
    ('1_kind_int8',             TYPE_INTEGER, 'long'),
    ('9.97d0',                  TYPE_REAL,    'double'),
    ('1.0d4',                   TYPE_REAL,    'double'),
    ('35.0d0',                  TYPE_REAL,    'double'),
    ('3.14',                    TYPE_REAL,    None),
    ('2.5e3',                   TYPE_REAL,    None),
])
def test_classify_literal(literal, expected_type, expected_kind):
    """Fortran literal defaults classify by syntax; `d` exponents read as double."""
    assert classify_literal(literal) == (expected_type, expected_kind)


@pytest.mark.parametrize("literal", [
    'huge(1_c_size_t)',                       # an intrinsic call, not a literal
    "inputPath(pathTypeDataDynamic)//'x'",    # a string-concatenation expression
    '',
    None,
])
def test_classify_literal_unrecognized(literal):
    """Expressions (function calls, concatenations) are not literals."""
    assert classify_literal(literal) == (None, None)


@pytest.mark.parametrize("declaration, expected_type, expected_kind", [
    ({'intrinsic': 'logical'},                                 TYPE_BOOLEAN, None),
    ({'intrinsic': 'integer'},                                 TYPE_INTEGER, None),
    ({'intrinsic': 'integer', 'type': 'c_long'},               TYPE_INTEGER, 'long'),
    ({'intrinsic': 'double precision'},                        TYPE_REAL,    'double'),
    ({'intrinsic': 'real'},                                    TYPE_REAL,    None),
    ({'intrinsic': 'character', 'type': 'len=32'},             TYPE_STRING,  None),
    ({'intrinsic': 'type', 'type': 'varying_string'},          TYPE_STRING,  None),
    ({'intrinsic': 'type', 'type': 'cosmologyParameters'},     TYPE_OBJECT,  'cosmologyparameters'),
])
def test_canonical_from_declaration(declaration, expected_type, expected_kind):
    """A scalar parameter's declared intrinsic maps to a canonical type."""
    assert canonical_from_declaration(declaration) == (expected_type, expected_kind)


def test_normalize_explicit_type():
    assert normalize_explicit_type('real')    == TYPE_REAL
    assert normalize_explicit_type('boolean') == TYPE_BOOLEAN
    assert normalize_explicit_type('integer') == TYPE_INTEGER
    assert normalize_explicit_type('string')  == TYPE_STRING
    assert normalize_explicit_type('nonsense') is None


def test_infer_explicit_type_wins():
    """An explicit `<type>` takes precedence over everything else."""
    result = infer_parameter_type({'name': 'x', 'type': 'real',
                                   'defaultValue': '.true.'})
    assert result == {'type': TYPE_REAL, 'kind': None, 'provenance': PROV_EXPLICIT}


def test_infer_from_default_literal():
    """With no explicit type, the default-value literal is used."""
    result = infer_parameter_type({'name': 'redshiftReionization',
                                   'defaultValue': '9.97d0'})
    assert result == {'type': TYPE_REAL, 'kind': 'double', 'provenance': PROV_DEFAULT}


def test_infer_from_declaration_when_no_default():
    """No explicit type and no usable default -> look up the variable's
    declaration via the supplied lookup (the path that types ~39% of params)."""
    def lookup(variable):
        assert variable == 'opticalDepthReionization'
        return {'intrinsic': 'double precision'}

    result = infer_parameter_type(
        {'name': 'opticalDepthReionization'}, declaration_lookup=lookup)
    assert result == {'type': TYPE_REAL, 'kind': 'double',
                      'provenance': PROV_DECLARATION}


def test_infer_uses_variable_over_name_for_lookup():
    """When `<variable>` differs from `<name>`, the declaration lookup keys on
    the variable actually written to."""
    seen = []

    def lookup(variable):
        seen.append(variable)
        return {'intrinsic': 'integer'}

    result = infer_parameter_type(
        {'name': 'countPerDecade', 'variable': 'countLocal'},
        declaration_lookup=lookup)
    assert seen == ['countLocal']
    assert result['type'] == TYPE_INTEGER


def test_infer_unknown_when_nothing_resolves():
    """An expression default with no declaration lookup is genuinely unknown."""
    result = infer_parameter_type({'name': 'limit', 'defaultValue': 'huge(0.0d0)'})
    assert result == {'type': TYPE_UNKNOWN, 'kind': None, 'provenance': PROV_UNKNOWN}


def test_infer_skips_declaration_lookup_for_array_element_variable():
    """A `name(i)` array-element variable is not a plain declaration key, so the
    lookup is skipped and the result is unknown (documents a known limitation)."""
    called = []
    result = infer_parameter_type(
        {'name': 'roots', 'variable': 'roots(1)'},
        declaration_lookup=lambda v: called.append(v))
    assert called == []
    assert result['provenance'] == PROV_UNKNOWN
