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
    is_output_array_arg,
    is_output_scalar_arg,
    normalize_method_return_type,
    unsupported_arg,
    unsupported_output_array_method,
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
])
def test_unsupported_array_shapes(intrinsic, type_spec, dim):
    arg = _arg(intrinsic, type_spec, attributes=[dim])
    verdict = classify_arg(arg, REGISTERED)
    assert verdict[0] == 'blocked'
    assert 'dimensioned argument' in verdict[1]


def test_class_array_argument_gets_specific_reject():
    # Arrays of polymorphic objects get a distinct, deferrable message
    # (one dynamic type per array; no intrinsic assignment into
    # polymorphic elements) rather than the generic dimensioned reject.
    arg = _arg('class', 'cosmologyFunctionsClass',
               attributes=['dimension(:)'])
    verdict = classify_arg(arg, REGISTERED)
    assert verdict[0] == 'blocked'
    assert 'array argument' in verdict[1]
    assert verdict[1].startswith('class(cosmologyFunctionsClass)')


def test_derived_type_array_reject_names_the_type():
    # `type(nBodyData), dimension(:)` — the reject must name the element
    # type so the audit's internal-derived-type rule can defer it.
    arg = _arg('type', 'nBodyData',
               attributes=['intent(inout)', 'dimension(:)'])
    verdict = classify_arg(arg, REGISTERED)
    assert verdict[0] == 'blocked'
    assert 'type(nBodyData)' in verdict[1]


def test_object_pointer_dummies_rejected():
    # class(*) pointer output — unsupportable in principle (deferred).
    star = _arg('class', '*', attributes=['intent(out)', 'pointer'],
                name='taskSelf')
    verdict = classify_arg(star, REGISTERED)
    assert verdict[0] == 'blocked'
    assert 'pointer output' in verdict[1]
    # type(treeNode) pointer inout — repointing would be silently lost;
    # actionable via a write-back protocol, so the wording must NOT
    # contain "type(" (which would misfile it as a deferred derived type).
    node = _arg('type', 'treeNode', attributes=['intent(inout)', 'pointer'],
                name='node')
    verdict = classify_arg(node, REGISTERED)
    assert verdict[0] == 'blocked'
    assert 'write-back' in verdict[1]
    assert 'type(' not in verdict[1]
    # A plain (non-pointer) treeNode scalar stays supported.
    plain = _arg('type', 'treeNode', attributes=['intent(inout)'],
                 name='node')
    assert classify_arg(plain, REGISTERED) is None


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


# ---------------------------------------------------------------------------
# Output-array arguments (`intent(out), allocatable, dimension(:)`)
# ---------------------------------------------------------------------------

def test_output_array_arg_recognised_and_accepted():
    oa = _arg('double precision',
              attributes=['intent(out)', 'allocatable', 'dimension(:)'],
              name='wavenumbers')
    assert is_output_array_arg(oa)
    # A per-arg classify accepts it (the method-level gate handles the rest).
    assert classify_arg(oa, REGISTERED) is None


@pytest.mark.parametrize('intrinsic', ['double precision', 'integer'])
def test_output_array_numeric_kinds(intrinsic):
    assert is_output_array_arg(
        _arg(intrinsic,
             attributes=['intent(out)', 'allocatable', 'dimension(:)']))


def test_inout_allocatable_is_output_array():
    # intent(inout) allocatable is Galacticus's allocation-reuse idiom —
    # treated as an output array (contents never read), same as intent(out).
    io = _arg('integer', type_spec='c_size_t',
              attributes=['allocatable', 'dimension(:)', 'intent(inout)'],
              name='indices')
    assert is_output_array_arg(io)
    assert classify_arg(io, REGISTERED) is None


def test_inout_non_allocatable_array_not_output():
    # A non-allocatable intent(inout) dimension(:) array is the in-place
    # mutable-buffer input path, not an output array.
    io = _arg('double precision',
              attributes=['intent(inout)', 'dimension(:)'], name='buf')
    assert not is_output_array_arg(io)


def test_logical_and_2d_output_arrays_accepted():
    # 1D logical (kind-narrowing copy into a c_bool export buffer) and 2D
    # numeric (second size companion + column-major reshape) both qualify.
    lg = _arg('logical',
              attributes=['intent(out)', 'allocatable', 'dimension(:)'])
    assert is_output_array_arg(lg)
    assert classify_arg(lg, REGISTERED) is None
    d2 = _arg('double precision',
              attributes=['intent(out)', 'allocatable', 'dimension(:,:)'])
    assert is_output_array_arg(d2)
    assert classify_arg(d2, REGISTERED) is None
    # 2D logical and character allocatables stay out.
    lg2 = _arg('logical',
               attributes=['intent(out)', 'allocatable', 'dimension(:,:)'])
    assert not is_output_array_arg(lg2)
    ch = _arg('character', type_spec='len=8',
              attributes=['intent(out)', 'allocatable', 'dimension(:)'])
    assert not is_output_array_arg(ch)


def test_output_array_method_gate_clean_cases():
    oa = _arg('double precision',
              attributes=['intent(out)', 'allocatable', 'dimension(:)'],
              name='y')
    inp = _arg('double precision', attributes=['intent(in)'], name='time')
    scalar_out = _arg('integer', attributes=['intent(out)'], name='count')
    dbl_out = _arg('double precision', attributes=['intent(out)'], name='f')
    # Void return, only inputs alongside the output → allowed.
    assert unsupported_output_array_method([oa], 'void') is None
    assert unsupported_output_array_method([inp, oa], 'void') is None
    # No output array at all → gate is a no-op.
    assert unsupported_output_array_method([inp], 'void') is None
    # Scalar intent(out) numeric companions are allowed alongside an array.
    assert is_output_scalar_arg(scalar_out)
    assert is_output_scalar_arg(dbl_out)
    assert unsupported_output_array_method([scalar_out, oa], 'void') is None
    assert unsupported_output_array_method([dbl_out, oa], 'void') is None


def test_output_array_method_gate_blocks():
    oa = _arg('double precision',
              attributes=['intent(out)', 'allocatable', 'dimension(:)'],
              name='y')
    optional_in = _arg('double precision',
                       attributes=['intent(in)', 'optional'], name='z')
    inout_scalar = _arg('double precision',
                        attributes=['intent(inout)'], name='x')
    # Subroutine-lowered returns, optional args, and scalar intent(inout)
    # numerics (an output the wrapper would silently drop) stay blocked.
    assert unsupported_output_array_method([oa], 'type(varying_string)')
    assert unsupported_output_array_method([optional_in, oa], 'void')
    assert unsupported_output_array_method([inout_scalar, oa], 'void')


def test_output_array_method_gate_allows_scalar_returns():
    # Direct-restype scalar returns may accompany output arrays; the
    # wrapper stays a bind(c) function and Python leads the tuple with it.
    oa = _arg('double precision',
              attributes=['intent(out)', 'allocatable', 'dimension(:)'],
              name='y')
    for ret in ('double precision', 'integer', 'integer(c_size_t)',
                'logical', 'integer(kind=c_size_t)',    # alias normalized
                # Dynamic-array returns compose with the output-array
                # companions (their own companions come first).
                'double precision, allocatable, dimension(:)',
                'double precision, allocatable, dimension(:,:)'):
        assert unsupported_output_array_method([oa], ret) is None, ret
    # Conversion-buffer returns stay blocked.
    for ret in ('type(varying_string)', 'class(fooClass)'):
        assert unsupported_output_array_method([oa], ret), ret


def test_output_array_method_gate_allows_in_place_and_class_args():
    # lengthToCellBoundary's shape: class intent(inout) (pointer-passed),
    # fixed-size intent(out) array (in-place buffer) — both allowed.
    oa = _arg('integer', type_spec='c_size_t', name='indicesNeighbor',
              attributes=['allocatable', 'dimension(:)', 'intent(inout)'])
    photon = _arg('class', type_spec='radiativeTransferPhotonPacketClass',
                  attributes=['intent(inout)'], name='photonPacket')
    pos_out = _arg('double precision', name='positionBoundary',
                   attributes=['dimension(3)', 'intent(out)'])
    assert unsupported_output_array_method(
        [photon, pos_out, oa], 'double precision') is None


def test_output_array_method_gate_allows_inout_allocatable():
    # An intent(inout) allocatable array is itself an output array, so a
    # method carrying only inputs + such arrays is allowed.
    oa_out   = _arg('double precision',
                    attributes=['intent(out)', 'allocatable', 'dimension(:)'],
                    name='y')
    oa_inout = _arg('integer', type_spec='c_size_t',
                    attributes=['intent(inout)', 'allocatable',
                                'dimension(:)'], name='indices')
    inp = _arg('double precision',
               attributes=['intent(in)', 'dimension(3)'], name='position')
    assert unsupported_output_array_method([inp, oa_inout], 'void') is None
    assert unsupported_output_array_method([oa_out, oa_inout], 'void') is None


# ---------------------------------------------------------------------------
# Procedure (callback) arguments
# ---------------------------------------------------------------------------

def test_registered_callback_accepted():
    # procedure(<iface>) with <iface> in _CALLBACK_PROCEDURE_INTERFACES:
    # marshalled as a Python callable via CFUNCTYPE → c_funptr → shim.
    cb = _arg('procedure', type_spec='computationalDomainVolumeIntegrand',
              name='integrand')
    assert classify_arg(cb, REGISTERED) is None


def test_unregistered_procedure_blocked_in_scope():
    cb = _arg('procedure', type_spec='someUnregisteredTemplate',
              attributes=['pointer'], name='fn')
    verdict = classify_arg(cb, REGISTERED)
    assert verdict[0] == 'blocked'
    assert 'procedure-pointer args are not supported' in verdict[1]


def test_registered_callback_pointer_dummy_accepted():
    # A `procedure(...), pointer` dummy (no intent) is accepted for a
    # registered interface: the wrapper passes the shim module's
    # procedure-pointer slot (a pointer actual for a pointer dummy).
    cb = _arg('procedure', type_spec='crossSectionFunctionTemplate',
              attributes=['pointer'], name='crossSectionFunction')
    assert classify_arg(cb, REGISTERED) is None


def test_registered_callback_pointer_output_still_blocked():
    # intent(out|inout) pointer procedures are outputs — blocked even for
    # a registered interface.
    cb = _arg('procedure', type_spec='crossSectionFunctionTemplate',
              attributes=['pointer', 'intent(out)'], name='fn')
    verdict = classify_arg(cb, REGISTERED)
    assert verdict[0] == 'blocked'
    assert 'pointer output' in verdict[1]


@pytest.mark.parametrize('intent', ['intent(out)', 'intent(inout)'])
def test_procedure_pointer_output_blocked_distinctly(intent):
    # procedure(...), pointer, intent(out|inout): a Fortran procedure
    # pointer handed back to the caller — unsupportable in principle; the
    # message is distinct so the audit buckets it out-of-scope.
    cb = _arg('procedure', type_spec='interruptTask',
              attributes=[intent, 'pointer'], name='functionInterrupt')
    verdict = classify_arg(cb, REGISTERED)
    assert verdict[0] == 'blocked'
    assert 'pointer output' in verdict[1]


def test_output_scalar_excludes_derived_and_character():
    # Only numeric/logical scalars are output companions; derived types,
    # class, and character intent(out) scalars are not.
    assert not is_output_scalar_arg(
        _arg('type', type_spec='abundances', attributes=['intent(out)']))
    assert not is_output_scalar_arg(
        _arg('character', type_spec='len=32', attributes=['intent(out)']))
    assert not is_output_scalar_arg(
        _arg('double precision',
             attributes=['intent(out)', 'dimension(3)']))  # fixed array
