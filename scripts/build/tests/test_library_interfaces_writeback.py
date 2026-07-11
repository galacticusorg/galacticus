"""Integration tests for the pointer write-back protocol: `type(X),
pointer, intent(out|inout)` method arguments (X a shared type, e.g.
treeNode).  The bind(c) wrapper takes the handle by reference, converts a
c_associated handle to a local Fortran pointer (null → disassociated — the
tree-walker start-of-iteration idiom), and writes the (re)pointed target's
c_loc back after the call; Python passes a ctypes.c_void_p updated in
place, so `while walker.next(node): …` iterates.  Runtime behaviour is
validated against gfortran-16 stubs (a 5-node chain walked from Python);
these tests pin the generated code's shape.
"""

import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
import libraryInterfaces as G  # noqa: E402


def _generate(fc):
    code   = {}
    python = {'c_lib': [], 'units': {fc['name']: {'subUnits': []}}}
    G.interfaces_methods(code, python, fc, {}, {}, {fc['name']: {}})
    fortran = '\n'.join(code.get(fc['name'], {}).get('shared', []))
    py      = '\n'.join(u['content']
                        for u in python['units'][fc['name']]['subUnits'])
    return fortran, python['c_lib'], py


NEXT = {
    'name'   : 'mergerTreeWalker',
    'module' : 'Merger_Tree_Walkers',
    'methods': {'next': {
        'type'    : 'logical',
        'argument': ['type(treeNode), intent(inout), pointer :: node'],
    }},
}


def test_writeback_method_generated():
    fc = dict(NEXT, methods={'next': dict(NEXT['methods']['next'])})
    fort, c_lib, py = _generate(fc)
    assert 'next' in fc['methods']                     # no longer dropped
    # Handle crosses by reference, not by value.
    assert 'type(c_ptr), intent(inout) :: node' in fort
    assert 'type(c_ptr), value :: node' not in fort
    spec = next(s for s in c_lib if s['name'] == 'mergerTreeWalkerNextL')
    assert spec['argtypes'] == ['c_void_p', 'c_int', 'POINTER(c_void_p)']
    # Guarded conversion in, write-back out.
    assert 'if (c_associated(node)) then' in fort
    assert 'call c_f_pointer(node, node_)' in fort
    assert 'node_ => null()' in fort
    assert 'node = c_loc(node_)' in fort
    assert 'node = c_null_ptr' in fort
    # Python: in-place c_void_p contract with a type guard.
    assert 'byref(node)' in py
    assert 'isinstance(node, c_void_p)' in py


def test_writeback_is_method_only():
    # The same shape on a *constructor* stays rejected (the constructor
    # wrapper has no post-call write-back hook) — via the shared predicate.
    from LibraryInterfaces.Classification import classify_arg
    arg = {'name': 'nodeTip', 'intrinsic': 'type', 'type': 'treeNode',
           'attributes': ['intent(inout)', 'pointer']}
    # Constructor context (default): blocked.
    assert classify_arg(arg, {})[0] == 'blocked'
    # Method context: supported.
    assert classify_arg(arg, {}, allow_pointer_writeback=True) is None


def test_writeback_requires_shared_type():
    # A locally-defined derived type has no known module to import the
    # type from — stays rejected even in method context.
    from LibraryInterfaces.Classification import classify_arg
    arg = {'name': 'x', 'intrinsic': 'type', 'type': 'someLocalType',
           'attributes': ['intent(inout)', 'pointer']}
    assert classify_arg(arg, {}, allow_pointer_writeback=True)[0] == 'blocked'


def test_optional_writeback_still_rejected():
    from LibraryInterfaces.Classification import classify_arg
    arg = {'name': 'lockNode', 'intrinsic': 'type', 'type': 'treeNode',
           'attributes': ['intent(out)', 'optional', 'pointer']}
    assert classify_arg(arg, {}, allow_pointer_writeback=True)[0] == 'blocked'
