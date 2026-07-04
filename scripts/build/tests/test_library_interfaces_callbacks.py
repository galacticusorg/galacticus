"""Integration tests for inbound procedure (callback) method arguments in
the libgalacticus generator.

Drives ``libraryInterfaces.interfaces_methods`` on a synthetic
functionClass carrying a `procedure(computationalDomainVolumeIntegrand)`
argument — the registered-callback shape — and pins the emitted Fortran
(per-method storage+shim module, funptr store, shim pass-through), the
ctypes signature, and the Python wrapper (CFUNCTYPE adaptation).  The
runtime path (Python callable → CFUNCTYPE → c_funptr → module slot → shim →
coordinate%toCartesian() → Python) is validated separately against
gfortran-16-built stubs; these tests keep the generated code's shape from
regressing without needing a compiler in CI.
"""

import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
import libraryInterfaces as G  # noqa: E402


INTEGRATE = {
    'name'   : 'computationalDomainVolumeIntegrator',
    'module' : 'Computational_Domain_Volume_Integrators',
    'methods': {'integrate': {
        'type'    : 'double precision',
        'argument': ['procedure(computationalDomainVolumeIntegrand)'
                     ' :: integrand'],
    }},
}


def _generate(fc):
    code   = {}
    python = {'c_lib': [], 'units': {fc['name']: {'subUnits': []}}}
    G.interfaces_methods(code, python, fc, {}, {},
                         {fc['name']: {}})
    fortran = '\n'.join(code.get(fc['name'], {}).get('shared', []))
    py      = '\n'.join(u['content']
                        for u in python['units'][fc['name']]['subUnits'])
    return fortran, python['c_lib'], py


def test_callback_method_generated():
    fc = {k: (dict(v) if isinstance(v, dict) else v)
          for k, v in INTEGRATE.items()}
    fc['methods'] = {'integrate': dict(INTEGRATE['methods']['integrate'])}
    _generate(fc)
    assert 'integrate' in fc['methods']       # not deleted


def test_callback_fortran_module_and_wrapper():
    fort, _, _ = _generate(dict(INTEGRATE,
                                methods={'integrate': dict(
                                    INTEGRATE['methods']['integrate'])}))
    # Storage+shim module precedes the wrapper in the compilation unit.
    mod = 'glcCB_computationalDomainVolumeIntegrator_integrate'
    assert fort.index(f'module {mod}') < fort.index(
        'function computationalDomainVolumeIntegratorIntegrateL')
    # c_funptr slot, C-side abstract interface, and the coordinate shim.
    assert 'type(c_funptr) :: glcCBPtr_integrand_ = c_null_funptr' in fort
    assert 'function glcCBIface_integrand_(position) bind(c)' in fort
    assert 'position_ = coordinates%toCartesian()' in fort
    assert 'call c_f_procpointer(glcCBPtr_integrand_, fn_)' in fort
    # Wrapper: funptr received by value, stored, shim passed to the inner.
    assert 'type(c_funptr), value :: integrand' in fort
    assert 'glcCBPtr_integrand_ = integrand' in fort
    assert 'integrand=glcCBShim_integrand_' in fort
    assert f'use :: {mod}, only : glcCBPtr_integrand_, glcCBShim_integrand_' \
        in fort


def test_callback_ctypes_signature():
    _, c_lib, _ = _generate(dict(INTEGRATE,
                                 methods={'integrate': dict(
                                     INTEGRATE['methods']['integrate'])}))
    spec = next(s for s in c_lib
                if s['name'] == 'computationalDomainVolumeIntegratorIntegrateL')
    assert spec['restype'] == 'c_double'
    # self, self_ID, then the funptr by value.
    assert spec['argtypes'] == ['c_void_p', 'c_int', 'c_void_p']


def test_callback_python_wrapper():
    _, _, py = _generate(dict(INTEGRATE,
                              methods={'integrate': dict(
                                  INTEGRATE['methods']['integrate'])}))
    assert 'def integrate(self,integrand):' in py
    # The user callable is adapted to a numpy view and wrapped in CFUNCTYPE.
    assert '_glcCBType_integrand_ = CFUNCTYPE(c_double, POINTER(c_double))' \
        in py
    assert 'np.ctypeslib.as_array(_pos_, (3,))' in py
    assert 'cast(_glcCB_integrand_, c_void_p)' in py


def test_unregistered_procedure_method_still_dropped():
    fc = {
        'name'   : 'radiationField',
        'module' : 'M',
        'methods': {'integrateOverCrossSection': {
            'type'    : 'double precision',
            'argument': ['procedure(crossSectionFunctionTemplate), pointer'
                         ' :: crossSectionFunction'],
        }},
    }
    _generate(fc)
    assert 'integrateOverCrossSection' not in fc['methods']
