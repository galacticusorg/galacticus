"""Regression tests for `Galacticus.Build.SourceTree.Parse.ModuleUses`.

Bug classes covered:

  1. `_flush_code` used to absorb `raw_pp_buf` whenever flushing a code
     buffer.  That broke the case where a `#ifdef USEMPI` line sat
     immediately before a `use mpi` line: the `#ifdef` got copied into the
     surrounding code node, and again into the moduleUse block when the
     use was rebuilt — yielding a duplicated `#ifdef … #endif` wrapper
     (one of which was unmatched and produced "unterminated #ifdef" at
     compile time).

  2. Issue #1030 — preprocessor condition conflation.  A module was keyed
     by name alone, so when the same module was `use`d twice in one unit
     under different `#ifdef` guards (or under a guard *and* unconditionally),
     the second occurrence's condition set overwrote — or absorbed the
     symbols of — the first, moving imports into the wrong preprocessor
     block.  Entries are now keyed by name *and* condition set: a module
     maps to a list of entries, one per distinct guard, each re-emitted
     under its own `#ifdef … #endif`.

  3. Issue #1030 — the use-statement regex matched assignments to variables
     whose name starts with "use" (e.g. `useCache=lastCache`), fabricating
     a spurious module-use node mid-procedure-body.
"""

import re

from Galacticus.Build.SourceTree                  import parse_code, serialize, walk_tree
from Galacticus.Build.SourceTree.Parse.ModuleUses import parse_module_uses, add_uses


def _parse_with_module_uses(source):
    """parse_code already runs `_pass_module_uses` for us, but `parse_code`
    expects either a real file or pre-instrumented content.  Use the
    documented `instrument=False` form so the test source is taken
    verbatim."""
    return parse_code(source, name='<test>', instrument=False)


def _module_use_node(tree):
    return next(n for n in walk_tree(tree) if n.get('type') == 'moduleUse')


def test_pp_guarded_use_does_not_duplicate_when_rebuilt():
    """When `update_uses` rebuilds the moduleUse content (which `add_uses`
    triggers), the `#ifdef USEMPI` wrapper must come from exactly ONE
    place — the rebuild — not from a surrounding code node that absorbed
    the original `#ifdef` line.  An earlier draft of `_flush_code`
    extended `raw_code_buf` with `raw_pp_buf` whenever it ran, so the
    `#ifdef` ended up duplicated between the code node and the rebuild,
    yielding "unterminated #ifdef" at compile time."""
    from Galacticus.Build.SourceTree.Parse.ModuleUses import update_uses

    source = (
        "module foo\n"
        "#ifdef USEMPI\n"
        "  use mpi\n"
        "#endif\n"
        "  implicit none\n"
        "end module foo\n"
    )
    tree = _parse_with_module_uses(source)

    update_uses(_module_use_node(tree))

    out = serialize(tree)
    assert out.count('#ifdef USEMPI') == 1, out
    assert out.count('#endif')        == 1, out
    # The `use mpi` survives.
    assert re.search(r'\buse\b\s*(::)?\s*mpi\b', out), out


def test_unconditional_add_uses_coexists_with_conditional():
    """If the tree already imports `mpi` under `#ifdef USEMPI`, an
    unconditional `add_uses(...)` for `mpi` must make `mpi` available in
    *every* build — otherwise `mpiSelf` is undefined in serial builds.

    Under the issue #1030 fix this no longer means deleting the existing
    `#ifdef USEMPI` entry; instead a separate, unconditional entry is added
    alongside it.  The guarded entry is harmless (a whole-module `use mpi`
    is idempotent), and crucially the module is now imported unconditionally."""
    source = (
        "module foo\n"
        "#ifdef USEMPI\n"
        "  use mpi\n"
        "#endif\n"
        "  implicit none\n"
        "end module foo\n"
    )
    tree = _parse_with_module_uses(source)
    module_use_node = _module_use_node(tree)

    # Sanity: the existing entry is keyed under the USEMPI condition.
    entries = module_use_node['moduleUse']['mpi']
    assert isinstance(entries, list)
    assert entries[0].get('conditions') == [{'name': 'USEMPI', 'invert': False}]

    # Now add an unconditional `use mpi`.  The new entry carries no
    # `conditions`, mirroring the call sites that simply want the module
    # imported in every build.
    add_uses(module_use_node['parent'], {
        'moduleUse':   {'mpi': {'openMP': False, 'intrinsic': False, 'all': True}},
        'moduleOrder': ['mpi'],
    })

    # An unconditional entry now exists in the structure.
    assert any('conditions' not in e
               for e in module_use_node['moduleUse']['mpi'])

    out = serialize(tree)
    # With every `#ifdef USEMPI … #endif` block stripped, `mpi` is still
    # imported — i.e. it is available unconditionally.
    unguarded = re.sub(r'#ifdef USEMPI.*?#endif', '', out, flags=re.S)
    assert re.search(r'\buse\b\s*(::)?\s*mpi\b', unguarded), out


def test_pp_guarded_use_attaches_conditions_to_entry():
    """Round-trip identity check on the parsed structure."""
    source = (
        "module foo\n"
        "#ifdef USEMPI\n"
        "  use mpi\n"
        "#endif\n"
        "end module foo\n"
    )
    tree = _parse_with_module_uses(source)

    entries = _module_use_node(tree)['moduleUse']['mpi']
    assert entries[0].get('conditions') == [{'name': 'USEMPI', 'invert': False}]


def test_dual_condition_module_kept_in_separate_blocks():
    """Issue #1030, conditional-first case (the `output.version.F90` pattern):
    a module used both inside an `#else` branch *and* unconditionally must
    keep each occurrence under its own guard.  An `add_uses` (as the OpenMP
    profiler does) triggers a rebuild, which is where the conflation used to
    surface."""
    source = (
        "subroutine demo\n"
        "#ifdef GIT2AVAIL\n"
        "    use, intrinsic :: ISO_C_Binding, only : c_null_char\n"
        "#else\n"
        "    use :: Input_Paths, only : pathTypeDataDynamic\n"
        "#endif\n"
        "    use :: Input_Paths, only : inputPath, pathTypeDataStatic\n"
        "    implicit none\n"
        "end subroutine demo\n"
    )
    tree = _parse_with_module_uses(source)
    subroutine = next(n for n in walk_tree(tree) if n.get('type') == 'subroutine')
    add_uses(subroutine, {
        'moduleUse':   {'OMP_Lib': {'intrinsic': False, 'all': True}},
        'moduleOrder': ['OMP_Lib'],
    })

    out = serialize(tree)

    # `pathTypeDataStatic` and `inputPath` are imported unconditionally; only
    # `pathTypeDataDynamic` belongs under `#ifndef GIT2AVAIL` (`#else`).
    unguarded = re.sub(r'#if[n]?def \w+.*?#endif', '', out, flags=re.S)
    assert 'pathTypeDataStatic' in unguarded, out
    assert 'inputPath'          in unguarded, out
    assert 'pathTypeDataDynamic' not in unguarded, out

    # Input_Paths appears under two distinct condition sets.
    entries = _module_use_node(tree)['moduleUse']['Input_Paths']
    keys = sorted(tuple(sorted((c['name'], c['invert']) for c in e.get('conditions', [])))
                  for e in entries)
    assert keys == [(), (('GIT2AVAIL', True),)], entries


def test_unconditional_first_conditional_later_not_conflated():
    """Issue #1030, unconditional-first case (the `utility.input_parameters.F90`
    pattern): `use Error, only : Error_Report` unconditionally, then
    `use Error, only : Warn` under `#else`.  `Error_Report` must stay
    unconditional — it used to be dragged under the later `#ifndef` guard."""
    source = (
        "subroutine demo2\n"
        "    use :: Error, only : Error_Report\n"
        "#ifdef GIT2AVAIL\n"
        "    use :: Input_Paths, only : pathTypeExec\n"
        "#else\n"
        "    use :: Error, only : Warn\n"
        "#endif\n"
        "    implicit none\n"
        "end subroutine demo2\n"
    )
    tree = _parse_with_module_uses(source)
    subroutine = next(n for n in walk_tree(tree) if n.get('type') == 'subroutine')
    add_uses(subroutine, {
        'moduleUse':   {'OMP_Lib': {'intrinsic': False, 'all': True}},
        'moduleOrder': ['OMP_Lib'],
    })

    out = serialize(tree)
    unguarded = re.sub(r'#if[n]?def \w+.*?#endif', '', out, flags=re.S)
    assert 'Error_Report' in unguarded, out      # available in every build
    assert 'Warn' not in unguarded, out          # only under #ifndef GIT2AVAIL


def test_assignment_to_use_prefixed_var_not_parsed_as_use():
    """Issue #1030: an assignment to a variable named `use…` is code, not a
    `use` statement.  Triggering a rebuild via `add_uses` must not destroy
    the assignment or fabricate a `Cache` module."""
    source = (
        "subroutine demo3\n"
        "    use :: Error, only : Error_Report\n"
        "    implicit none\n"
        "    integer :: useCache, lastCache\n"
        "    useCache  =lastCache\n"
        "    return\n"
        "end subroutine demo3\n"
    )
    tree = _parse_with_module_uses(source)
    subroutine = next(n for n in walk_tree(tree) if n.get('type') == 'subroutine')
    add_uses(subroutine, {
        'moduleUse':   {'OMP_Lib': {'intrinsic': False, 'all': True}},
        'moduleOrder': ['OMP_Lib'],
    })

    out = serialize(tree)
    assert 'useCache  =lastCache' in out, out
    # No spurious module parsed out of the assignment.
    for node in walk_tree(tree):
        if node.get('type') == 'moduleUse':
            assert 'Cache' not in node['moduleUse'], node['moduleUse']


def test_use_double_colon_no_space_accepted():
    """The boundary fix must still accept the `use::module` spelling (which
    the Perl regex wrongly rejected)."""
    source = (
        "module foo\n"
        "  use::iso_fortran_env\n"
        "  implicit none\n"
        "end module foo\n"
    )
    tree = _parse_with_module_uses(source)
    assert 'iso_fortran_env' in _module_use_node(tree)['moduleUse']
