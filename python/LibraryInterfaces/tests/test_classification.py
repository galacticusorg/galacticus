"""Tests for LibraryInterfaces.Classification — the argument/return-type
rules shared by the library-interface generator and the audit tool. These
pin the rules that decide which functionClass implementations get exposed
through libgalacticus, including the generator/audit divergences the shared
module was created to prevent (e.g. the omp_lock_kind rejection the audit
had silently lost).
"""

import pytest

from LibraryInterfaces.Classification import (
    classify_arg,
    is_internal_constructor_name,
    normalize_method_return_type,
    unsupported_arg,
)

REGISTERED = {'cosmologyFunctions': {}, 'darkMatterHaloScale': {}}
KNOWN      = set(REGISTERED) | {'unregisteredThing'}


def _arg(intrinsic, type_spec='', attributes=(), name='x'):
    return {'name': name, 'intrinsic': intrinsic, 'type': type_spec,
            'attributes': list(attributes)}


# ---------------------------------------------------------------------------
# Scalars and outright rejections
# ---------------------------------------------------------------------------

def test_plain_scalar_supported():
    assert classify_arg(_arg('double precision'), REGISTERED) is None
    assert classify_arg(_arg('integer'), REGISTERED) is None


@pytest.mark.parametrize('intrinsic', ['complex', 'double complex'])
def test_complex_rejected(intrinsic):
    verdict = classify_arg(_arg(intrinsic), REGISTERED)
    assert verdict[0] == 'blocked'


def test_procedure_rejected():
    verdict = classify_arg(_arg('procedure', 'integrand'), REGISTERED)
    assert verdict[0] == 'blocked'
    assert 'procedure' in verdict[1]


def test_omp_lock_kind_nonoptional_rejected():
    # The rule whose absence from the audit's hand-synced copy was the
    # original drift this module exists to prevent.
    verdict = classify_arg(_arg('integer', 'omp_lock_kind'), REGISTERED)
    assert verdict[0] == 'blocked'
    assert 'omp_lock_kind' in verdict[1]


def test_omp_lock_kind_optional_accepted():
    arg = _arg('integer', 'omp_lock_kind', attributes=['optional'])
    assert classify_arg(arg, REGISTERED) is None


# ---------------------------------------------------------------------------
# class(...) arguments
# ---------------------------------------------------------------------------

def test_registered_class_supported():
    arg = _arg('class', 'cosmologyFunctionsClass')
    assert classify_arg(arg, REGISTERED) is None


def test_unregistered_class_blocked_in_generator_mode():
    arg = _arg('class', 'unregisteredThingClass')
    verdict = classify_arg(arg, REGISTERED)
    assert verdict[0] == 'blocked'
    assert unsupported_arg(arg, REGISTERED) == verdict[1]


def test_unregistered_class_is_missing_dep_in_audit_mode():
    arg = _arg('class', 'unregisteredThingClass')
    verdict = classify_arg(arg, REGISTERED, known_function_classes=KNOWN)
    assert verdict == ('missing-dep', {'unregisteredThing'})


def test_unknown_class_blocked_in_both_modes():
    arg = _arg('class', 'treeNode')
    assert classify_arg(arg, REGISTERED)[0] == 'blocked'
    assert classify_arg(arg, REGISTERED,
                        known_function_classes=KNOWN)[0] == 'blocked'


def test_class_star_needs_override():
    arg = _arg('class', '*')
    assert classify_arg(arg, REGISTERED)[0] == 'blocked'
    assert classify_arg(arg, REGISTERED,
                        constructor_overrides=[{'name': 'x'}]) is None


def test_abstract_intermediate_resolved_through_hierarchy():
    # An intermediate type whose extends-chain reaches a registered
    # <base>Class is routed through <base>GetPtr.
    hierarchy = {'cosmologyFunctionsMatterLambda':
                 {'parent': 'cosmologyFunctionsClass', 'module': 'm'}}
    arg = _arg('class', 'cosmologyFunctionsMatterLambda')
    assert classify_arg(arg, REGISTERED, class_hierarchy=hierarchy) is None
    # ...but with no hierarchy the same arg is blocked.
    assert classify_arg(arg, REGISTERED)[0] == 'blocked'


# ---------------------------------------------------------------------------
# Overrides
# ---------------------------------------------------------------------------

def test_null_override_accepts_procedure_and_class_star():
    overrides = [{'name': 'x', 'value': 'null'}]
    assert classify_arg(_arg('procedure', 'integrand'), REGISTERED,
                        constructor_overrides=overrides) is None
    assert classify_arg(_arg('class', '*'), REGISTERED,
                        constructor_overrides=overrides) is None
    verdict = classify_arg(_arg('integer'), REGISTERED,
                           constructor_overrides=overrides)
    assert verdict[0] == 'blocked'


def test_absent_override_requires_optional():
    overrides = [{'name': 'x', 'value': 'absent'}]
    assert classify_arg(_arg('integer', attributes=['optional']), REGISTERED,
                        constructor_overrides=overrides) is None
    assert classify_arg(_arg('integer'), REGISTERED,
                        constructor_overrides=overrides)[0] == 'blocked'


# ---------------------------------------------------------------------------
# Array shapes
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('intrinsic,type_spec,dim', [
    ('double precision', '',      'dimension(:)'),
    ('double precision', '',      'dimension(:,:)'),
    ('integer',          '',      'dimension(3)'),
    ('logical',          '',      'dimension(:)'),
    ('character',        'len=5', 'dimension(:)'),
    ('type',             'varying_string', 'dimension(:)'),
])
def test_supported_array_shapes(intrinsic, type_spec, dim):
    arg = _arg(intrinsic, type_spec, attributes=[dim, 'intent(in   )'])
    assert classify_arg(arg, REGISTERED) is None


def test_allocatable_array_rejected():
    arg = _arg('double precision', '',
               attributes=['dimension(:)', 'allocatable'])
    verdict = classify_arg(arg, REGISTERED)
    assert verdict[0] == 'blocked'
    assert 'allocatable' in verdict[1]


@pytest.mark.parametrize('intrinsic,type_spec,dim', [
    ('character', 'len=*', 'dimension(:)'),   # no fixed byte stride
    ('logical',   '',      'dimension(3)'),   # fixed-size logical unsupported
    ('class',     'cosmologyFunctionsClass', 'dimension(:)'),
])
def test_unsupported_array_shapes(intrinsic, type_spec, dim):
    arg = _arg(intrinsic, type_spec, attributes=[dim])
    verdict = classify_arg(arg, REGISTERED)
    assert verdict[0] == 'blocked'
    assert 'dimensioned argument' in verdict[1]


# ---------------------------------------------------------------------------
# Small shared helpers
# ---------------------------------------------------------------------------

def test_return_type_aliases():
    assert normalize_method_return_type('integer(kind=c_size_t)') \
        == 'integer(c_size_t)'
    assert normalize_method_return_type('integer(kind_int8)') \
        == 'integer(c_long)'
    assert normalize_method_return_type('double precision') \
        == 'double precision'


def test_internal_constructor_names():
    assert is_internal_constructor_name('coleConstructorInternal')
    assert is_internal_constructor_name('allAndFormationNodesInternal')
    assert is_internal_constructor_name('fooConstructorInternalType')
    assert not is_internal_constructor_name('fooConstructorParameters')
