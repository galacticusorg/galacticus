"""Tests for useDependencies.py's conditional-compilation handling.

Covers the two state machines that previously mis-parsed real inputs:

* the Makefile-conditional scanner (`_collect_preprocessor_directives`),
  which desynchronized its stack on `ifneq` (no push, but its `endif`
  still popped) and ignored `else`;
* the source-side preprocessor stack, which treated `#if <expression>`
  blocks as *inactive* (silently dropping the `use` statements inside
  them), did not recognize `#elif`, and rejected lowercase `#ifndef`
  macro names.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
import useDependencies as ud  # noqa: E402


# ---------------------------------------------------------------------------
# _update_conditional_compile
# ---------------------------------------------------------------------------

def test_ifdef_defined_active():
    stack = [{'name': 'USEMPI', 'state': 1}]
    assert ud._update_conditional_compile(stack, {'USEMPI'})


def test_ifdef_undefined_inactive():
    stack = [{'name': 'USEMPI', 'state': 1}]
    assert not ud._update_conditional_compile(stack, set())


def test_ifndef_undefined_active():
    stack = [{'name': 'USEMPI', 'state': 0}]
    assert ud._update_conditional_compile(stack, set())


def test_unknown_entry_always_active():
    # `#if <expression>` blocks are unevaluable: every branch is scanned.
    stack = [{'name': None, 'state': 1, 'unknown': True}]
    assert ud._update_conditional_compile(stack, set())
    # ...even after an #else flip.
    stack[0]['state'] = 0
    assert ud._update_conditional_compile(stack, set())


def test_nested_false_branch_wins():
    stack = [
        {'name': None, 'state': 1, 'unknown': True},
        {'name': 'USEMPI', 'state': 1},
    ]
    assert not ud._update_conditional_compile(stack, set())


# ---------------------------------------------------------------------------
# Source-side stack, driven through a real file scan
# ---------------------------------------------------------------------------

def _scan_modules(tmp_path, source_text, defined=frozenset()):
    """Run _scan_source_file over `source_text` and return the module names
    recorded in modulesUsed (stripped of the work-dir prefix)."""
    src = tmp_path / 'probe.F90'
    src.write_text(source_text)
    entry = {
        'files':                [str(src)],
        'modulesUsed':          [],
        'dependenciesExplicit': [],
        'modulesProvided':      {},
        'submodules':           [],
        'submodulesProvided':   [],
        'libraryDependencies':  {},
    }
    ud._scan_source_file(
        entry, [str(src)],
        {'fileName': 'probe.F90', 'subDirectoryName': ''},
        {name: [] for name in ud._DIRECTIVE_NAMES},
        {}, str(tmp_path), 'WORK/', frozenset(defined),
    )
    return {m[len('WORK/'):-len('.mod')] for m in entry['modulesUsed']}


def test_plain_use_recorded(tmp_path):
    mods = _scan_modules(tmp_path, (
        'module probe\n'
        '  use Module_A\n'
        'end module probe\n'
    ))
    assert 'module_a' in mods


def test_ifdef_branches(tmp_path):
    text = (
        'module probe\n'
        '#ifdef USEMPI\n'
        '  use Module_Mpi_Side\n'
        '#else\n'
        '  use Module_Serial_Side\n'
        '#endif\n'
        'end module probe\n'
    )
    assert _scan_modules(tmp_path, text, {'USEMPI'}) == {'module_mpi_side'}
    assert _scan_modules(tmp_path, text) == {'module_serial_side'}


def test_hash_if_block_is_scanned(tmp_path):
    # Previously a `use` inside `#if defined(A) || defined(B)` was silently
    # dropped from the dependency graph.
    mods = _scan_modules(tmp_path, (
        'module probe\n'
        '#if defined(OFDAVAIL) || defined(OFDLOCKS)\n'
        '  use Module_Locks\n'
        '#else\n'
        '  use Module_No_Locks\n'
        '#endif\n'
        'end module probe\n'
    ))
    # Both branches of an unevaluable conditional are scanned.
    assert {'module_locks', 'module_no_locks'} <= mods


def test_elif_degrades_chain_to_all_branches(tmp_path):
    text = (
        'module probe\n'
        '#ifdef USEMPI\n'
        '  use Module_A\n'
        '#elif defined(OTHER)\n'
        '  use Module_B\n'
        '#else\n'
        '  use Module_C\n'
        '#endif\n'
        'end module probe\n'
    )
    # The leading #ifdef is evaluable and stays exact; from the #elif onward
    # the chain is unevaluable, so every later branch is scanned.
    mods = _scan_modules(tmp_path, text)
    assert 'module_a' not in mods         # USEMPI known-undefined
    assert {'module_b', 'module_c'} <= mods
    mods = _scan_modules(tmp_path, text, {'USEMPI'})
    assert {'module_a', 'module_b', 'module_c'} <= mods


def test_lowercase_ifndef(tmp_path):
    # `#ifndef __aarch64__` was previously unmatched (uppercase-only regex),
    # leaving its `#endif` to unbalance the stack.
    mods = _scan_modules(tmp_path, (
        'module probe\n'
        '#ifndef __aarch64__\n'
        '  use Module_X86\n'
        '#endif\n'
        '  use Module_Always\n'
        'end module probe\n'
    ))
    assert 'module_always' in mods
    assert 'module_x86' in mods  # __aarch64__ not in the defined set


# ---------------------------------------------------------------------------
# Makefile-conditional scanner
# ---------------------------------------------------------------------------

def _collect_from(tmp_path, monkeypatch, makefile_text, env=()):
    (tmp_path / 'Makefile').write_text(makefile_text)
    build = tmp_path / 'build'
    build.mkdir(exist_ok=True)
    monkeypatch.chdir(tmp_path)
    for name in env:
        monkeypatch.setenv(name, '1')
    monkeypatch.setenv('GALACTICUS_FCFLAGS', '')
    # Neutralize the compiler-macro merge so tests see only Makefile flags.
    monkeypatch.setenv('CCOMPILER', '/bin/true')
    return set(ud._collect_preprocessor_directives(str(build)))


def test_ifneq_no_longer_desynchronizes(tmp_path, monkeypatch):
    # Previously `ifneq` did not push but its `endif` popped, so the endif
    # here re-activated the enclosing false `ifdef` branch and -DWRONG was
    # collected.
    flags = _collect_from(tmp_path, monkeypatch, (
        'ifdef NOT_SET_ANYWHERE_XYZ\n'
        'ifneq ($(THING),)\n'
        'FCFLAGS += -DINNER\n'
        'endif\n'
        'FCFLAGS += -DWRONG\n'
        'endif\n'
        'FCFLAGS += -DALWAYS\n'
    ))
    assert 'ALWAYS' in flags
    assert 'WRONG' not in flags
    assert 'INNER' not in flags


def test_else_of_env_ifdef(tmp_path, monkeypatch):
    text = (
        'ifdef SOME_TEST_ENV_FLAG\n'
        'FCFLAGS += -DWITH\n'
        'else\n'
        'FCFLAGS += -DWITHOUT\n'
        'endif\n'
    )
    assert 'WITHOUT' in _collect_from(tmp_path, monkeypatch, text)
    assert 'WITH' not in _collect_from(tmp_path, monkeypatch, text)
    with_env = _collect_from(tmp_path, monkeypatch, text,
                             env=('SOME_TEST_ENV_FLAG',))
    assert 'WITH' in with_env
    assert 'WITHOUT' not in with_env


def test_ifeq_collects_all_branches(tmp_path, monkeypatch):
    # Unevaluable make-variable conditions: deliberately over-approximate.
    flags = _collect_from(tmp_path, monkeypatch, (
        "ifeq '$(OPT)' 'a'\n"
        'FCFLAGS += -DBRANCH_A\n'
        "else ifeq '$(OPT)' 'b'\n"
        'FCFLAGS += -DBRANCH_B\n'
        'else\n'
        'FCFLAGS += -DBRANCH_C\n'
        'endif\n'
    ))
    assert {'BRANCH_A', 'BRANCH_B', 'BRANCH_C'} <= flags


def test_dash_d_with_underscore(tmp_path, monkeypatch):
    flags = _collect_from(tmp_path, monkeypatch,
                          'FCFLAGS += -DUSE_MPI -Dmixed_Case1\n')
    assert 'USE_MPI' in flags
    assert 'mixed_Case1' in flags
