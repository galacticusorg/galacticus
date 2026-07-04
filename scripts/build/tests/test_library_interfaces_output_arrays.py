"""Integration tests for output-array method arguments
(`intent(out), allocatable, dimension(:)`) in the libgalacticus generator.

Drives ``libraryInterfaces.interfaces_methods`` on synthetic functionClasses
and pins the generated Fortran wrapper, the ctypes signature, and the Python
wrapper, plus the whole-method gate that keeps disqualified methods blocked.
The runtime behaviour of the emitted Python (byref out-params, from_address,
copy, tuple assembly, zero-size guard) is validated separately against a C
shim that mimics the Fortran ABI; these tests keep the *shape* of the
generated code from regressing without needing a compiler in CI.
"""

import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
import libraryInterfaces as G  # noqa: E402


def _generate(fc, lib_function_classes):
    """Run interfaces_methods on *fc* and return (fortran, c_lib, python)."""
    code   = {}
    python = {'c_lib': [], 'units': {fc['name']: {'subUnits': []}}}
    G.interfaces_methods(code, python, fc, {}, {}, lib_function_classes)
    fortran = '\n'.join(code.get(fc['name'], {}).get('shared', []))
    py      = '\n'.join(u['content']
                        for u in python['units'][fc['name']]['subUnits'])
    return fortran, python['c_lib'], py


def _method_fc(name, module, method_name, return_type, arguments):
    return {'name': name, 'module': module,
            'methods': {method_name: {'type': return_type,
                                      'argument': arguments}}}


# ---------------------------------------------------------------------------
# Single output array — the flagship transferFunction.wavenumbersLocalMinima
# ---------------------------------------------------------------------------

WAVENUMBERS = _method_fc(
    'transferFunction', 'Transfer_Functions', 'wavenumbersLocalMinima', 'void',
    ['double precision, intent(  out), allocatable, dimension(:) :: wavenumbers'])


def test_single_output_array_fortran():
    fort, _, _ = _generate(WAVENUMBERS, {'transferFunction': {}})
    # A void-returning subroutine (not a function).
    assert 'subroutine transferFunctionWavenumbersLocalMinimaL' in fort
    # The output array becomes a save/target local, not a bind(c) formal.
    assert ('real(c_double), dimension(:), allocatable, save, target'
            ' :: glcOut_wavenumbers_' in fort)
    # (c_ptr, c_size_t) intent(out) companions in the signature.
    assert 'type(c_ptr),       intent(out) :: wavenumbersDataPtr_' in fort
    assert 'integer(c_size_t), intent(out) :: wavenumbersSize_' in fort
    # Pre-call hygiene + inner call fills the local.
    assert 'if (allocated(glcOut_wavenumbers_)) deallocate(glcOut_wavenumbers_)' in fort
    assert 'wavenumbers=glcOut_wavenumbers_' in fort
    # Post-call guarded c_loc/size (zero-size / unallocated → null + 0).
    assert 'wavenumbersSize_ = size(glcOut_wavenumbers_, kind=c_size_t)' in fort
    assert 'wavenumbersDataPtr_ = c_loc(glcOut_wavenumbers_)' in fort
    assert 'wavenumbersDataPtr_ = c_null_ptr' in fort
    # Required ISO_C_Binding symbols imported.
    for sym in ('c_ptr', 'c_size_t', 'c_loc', 'c_null_ptr', 'c_double'):
        assert sym in fort


def test_single_output_array_clib_signature():
    _, c_lib, _ = _generate(WAVENUMBERS, {'transferFunction': {}})
    spec = next(s for s in c_lib
                if s['name'] == 'transferFunctionWavenumbersLocalMinimaL')
    assert spec['restype'] is None            # subroutine
    # self, self_ID, then the output companion pair.
    assert spec['argtypes'] == [
        'c_void_p', 'c_int', 'POINTER(c_void_p)', 'POINTER(c_size_t)']


def test_single_output_array_python():
    _, _, py = _generate(WAVENUMBERS, {'transferFunction': {}})
    # No output-array parameter in the Python signature.
    assert 'def wavenumbersLocalMinima(self):' in py
    assert '_wavenumbersPtr_  = c_void_p()' in py
    assert '_wavenumbersSize_ = c_size_t()' in py
    assert 'byref(_wavenumbersPtr_),byref(_wavenumbersSize_)' in py
    # Zero-size guard returns an empty typed array without dereferencing.
    assert 'if _wavenumbersSize_.value:' in py
    assert 'np.empty(0, dtype=np.float64)' in py
    assert '.from_address(_wavenumbersPtr_.value)' in py
    assert '.copy()' in py
    assert 'return _wavenumbersArr_' in py


# ---------------------------------------------------------------------------
# Multiple output arrays — mergerTreeBuildMasses.construct (returns a tuple)
# ---------------------------------------------------------------------------

CONSTRUCT = _method_fc(
    'mergerTreeBuildMasses', 'Merger_Tree_Build_Masses', 'construct', 'void',
    ['double precision, intent(in   ) :: time',
     'double precision, intent(  out), allocatable, dimension(:) ::'
     ' mass, massMinimum, massMaximum, weight'])


def test_multiple_output_arrays_return_tuple():
    fort, c_lib, py = _generate(CONSTRUCT, {'mergerTreeBuildMasses': {}})
    # Four independent save/target buffers.
    for nm in ('mass', 'massMinimum', 'massMaximum', 'weight'):
        assert f'glcOut_{nm}_' in fort
        assert f'{nm}DataPtr_' in fort
    # The ordinary input arg is still a normal Python parameter.
    assert 'def construct(self,time):' in py
    # Return is a 4-tuple in declaration order.
    assert ('return (_massArr_, _massMinimumArr_, _massMaximumArr_,'
            ' _weightArr_)' in py)
    # ctypes signature: self, self_ID, time, then 4 companion pairs.
    spec = next(s for s in c_lib if s['name'] == 'mergerTreeBuildMassesConstructL')
    assert spec['argtypes'].count('POINTER(c_void_p)') == 4
    assert spec['argtypes'].count('POINTER(c_size_t)') == 4
    assert spec['argtypes'][:3] == ['c_void_p', 'c_int', 'c_double']


# ---------------------------------------------------------------------------
# Whole-method gate: disqualified output-array methods stay blocked (the
# generator deletes them from the class, emitting a caution).
# ---------------------------------------------------------------------------

def _method_deleted(fc, method_name, lib_function_classes):
    code   = {}
    python = {'c_lib': [], 'units': {fc['name']: {'subUnits': []}}}
    G.interfaces_methods(code, python, fc, {}, {}, lib_function_classes)
    return method_name not in fc['methods']


def test_scalar_out_companion_supported():
    # accretionDiskSpectra.wavelengths: scalar intent(out) count + array —
    # the scalar is returned alongside the array, not a Python input.
    fc = _method_fc(
        'accretionDiskSpectra', 'M', 'wavelengths', 'void',
        ['integer, intent(  out) :: wavelengthsCount',
         'double precision, allocatable, dimension(:), intent(  out)'
         ' :: wavelengths'])
    fort, c_lib, py = _generate(fc, {'accretionDiskSpectra': {}})
    assert 'wavelengths' in fc['methods']          # not deleted
    # Scalar out is a by-reference intent(out) dummy passed straight through.
    assert 'integer(c_int) :: wavelengthsCount' in fort
    assert 'wavelengthsCount=wavelengthsCount' in fort
    # No Python input parameter for the outputs; both returned in order.
    assert 'def wavelengths(self):' in py
    assert '_wavelengthsCount_ = c_int()' in py
    assert 'byref(_wavelengthsCount_)' in py
    assert 'return (_wavelengthsCount_.value, _wavelengthsArr_)' in py
    spec = next(s for s in c_lib if s['name'] == 'accretionDiskSpectraWavelengthsL')
    # self, self_ID, POINTER(c_int) scalar-out, then the array companion pair.
    assert spec['argtypes'] == [
        'c_void_p', 'c_int', 'POINTER(c_int)',
        'POINTER(c_void_p)', 'POINTER(c_size_t)']


def test_scalar_and_array_outputs_declaration_order():
    # variogram.modelFDF: input arrays + scalar out `f` + output array `dfdC`.
    fc = _method_fc(
        'variogram', 'M', 'modelFDF', 'void',
        ['double precision, intent(in   ), dimension(:) :: C, separations,'
         ' semivariances',
         'double precision, intent(  out) :: f',
         'double precision, intent(  out), allocatable, dimension(:) :: dfdC'])
    _, _, py = _generate(fc, {'variogram': {}})
    assert 'def modelFDF(self,C,separations,semivariances):' in py
    assert '_f_ = c_double()' in py
    # Scalar output precedes the array output in declaration order.
    assert 'return (_f_.value, _dfdCArr_)' in py


def test_non_void_return_blocks_output_array_method():
    fc = _method_fc(
        'x', 'M', 'f', 'double precision',
        ['double precision, intent(  out), allocatable, dimension(:) :: y'])
    assert _method_deleted(fc, 'f', {'x': {}})


def test_optional_arg_blocks_output_array_method():
    fc = _method_fc(
        'x', 'M', 'f', 'void',
        ['double precision, intent(in   ), optional :: z',
         'double precision, intent(  out), allocatable, dimension(:) :: y'])
    assert _method_deleted(fc, 'f', {'x': {}})


def test_inout_allocatable_is_not_an_output_array():
    # intent(inout) allocatable stays blocked (per-arg), method dropped.
    fc = _method_fc(
        'x', 'M', 'f', 'void',
        ['double precision, intent(inout), allocatable, dimension(:) :: y'])
    assert _method_deleted(fc, 'f', {'x': {}})


def test_plain_method_still_generated():
    # A method with no output array is unaffected and still emitted.
    fc = _method_fc(
        'transferFunction', 'M', 'value', 'double precision',
        ['double precision, intent(in   ) :: wavenumber'])
    assert not _method_deleted(fc, 'value', {'transferFunction': {}})
