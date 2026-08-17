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

  4. Issue #1385 — a `#ifdef` was recorded as a condition on the `use`
     statements it preceded whether or not it *closed* around them.  Where
     the guard also covered the code after the `use` statements, the rebuilt
     block emitted its own `#ifdef … #endif` around just the `use` while the
     original `#endif` survived in the following code node — one `#endif` too
     many, which fails to compile.  A guard is now claimed only when the
     region it opens contains nothing but `use` statements and directives
     that can themselves be rebuilt.
"""

import re

from Galacticus.Build.SourceTree                  import parse_code, serialize, walk_tree
from Galacticus.Build.SourceTree.Parse.ModuleUses import parse_module_uses, add_uses, update_uses


def _parse_with_module_uses(source):
    """parse_code already runs `_pass_module_uses` for us, but `parse_code`
    expects either a real file or pre-instrumented content.  Use the
    documented `instrument=False` form so the test source is taken
    verbatim."""
    return parse_code(source, name='<test>', instrument=False)


def _module_use_node(tree):
    return next(n for n in walk_tree(tree) if n.get('type') == 'moduleUse')


def _module_use_nodes(tree):
    return [n for n in walk_tree(tree) if n.get('type') == 'moduleUse']


def _rebuild(tree):
    """Rebuild every moduleUse block and return the serialized result.

    This is what `scripts/aux/formatModuleUses.py` does, and what a build
    does to any unit an `add_uses` touches.
    """
    for node in _module_use_nodes(tree):
        update_uses(node)
    return serialize(tree)


def _guard_balance(text):
    """`(closing balance, deepest negative excursion)` of the conditional
    directives in `text`.  Both are zero for well-formed source; either
    going negative means an `#endif` without a matching `#if`."""
    depth   = 0
    minimum = 0
    for line in text.splitlines():
        if re.match(r'^#\s*if', line):
            depth += 1
        elif re.match(r'^#\s*endif', line):
            depth -= 1
        minimum = min(minimum, depth)
    return depth, minimum


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
    """The boundary fix must still accept the `use::module` spelling
    (historically rejected)."""
    source = (
        "module foo\n"
        "  use::iso_fortran_env\n"
        "  implicit none\n"
        "end module foo\n"
    )
    tree = _parse_with_module_uses(source)
    assert 'iso_fortran_env' in _module_use_node(tree)['moduleUse']


# ---------------------------------------------------------------------------
# Issue #1385 — guards that span more than the `use` statements they precede
# ---------------------------------------------------------------------------

# The `source/error/_module.F90` pattern: one `#ifdef` covers a `use` *and*
# the `implicit none` and declarations that follow it.
_SPANNING_GUARD = (
    "subroutine handler()\n"
    "    use, intrinsic :: ISO_Fortran_Env, only : error_unit\n"
    "    use            :: System_Output  , only : stdOutIsATTY\n"
    "#ifdef USEMPI\n"
    "    use            :: MPI_F08        , only : MPI_Comm_Rank, MPI_Comm_World\n"
    "    implicit none\n"
    "    integer            :: mpiRank , error\n"
    "    character(len=128) :: hostName\n"
    "#endif\n"
    "    call backtrace()\n"
    "end subroutine handler\n"
)


def test_guard_spanning_code_leaves_directives_balanced_when_rebuilt():
    """The bug itself: rebuilding must not add an `#ifdef … #endif` of its
    own around the guarded `use`, because the original `#endif` — which
    closes *after* the declarations — is still there in the code node."""
    tree = _parse_with_module_uses(_SPANNING_GUARD)
    out  = _rebuild(tree)

    assert out.count('#ifdef USEMPI') == 1, out
    assert out.count('#endif')        == 1, out
    assert _guard_balance(out) == (0, 0), out
    # The guarded declarations still sit inside the guard, in order.
    assert re.search(r'#ifdef USEMPI\n.*MPI_F08.*\n\s*implicit none\n.*'
                     r'hostName.*\n#endif\n', out, re.S), out


def test_guard_spanning_code_is_not_recorded_as_a_condition():
    """A guard that is left in place must not also be recorded on the `use`
    statements it covers — that is what makes the rebuild emit it twice."""
    tree  = _parse_with_module_uses(_SPANNING_GUARD)
    nodes = _module_use_nodes(tree)

    entries = next(n['moduleUse']['MPI_F08'] for n in nodes
                   if 'MPI_F08' in n['moduleUse'])
    assert 'conditions' not in entries[0], entries
    # …and the block holding it is flagged as sitting inside a guard.
    assert any(n.get('inUnclaimedGuard') for n in nodes), nodes


def test_guard_over_uses_only_is_still_claimed():
    """The narrow case the condition machinery exists for is unaffected: a
    guard wrapping nothing but `use` statements is still absorbed, and
    rebuilt from the entry's conditions."""
    source = (
        "module foo\n"
        "#ifdef USEMPI\n"
        "  use :: mpi\n"
        "  use :: mpi_helpers\n"
        "#endif\n"
        "  implicit none\n"
        "end module foo\n"
    )
    tree = _parse_with_module_uses(source)
    node = _module_use_node(tree)
    assert node['moduleUse']['mpi'][0]['conditions'] == \
        [{'name': 'USEMPI', 'invert': False}]
    assert not node.get('inUnclaimedGuard')

    out = _rebuild(tree)
    assert out.count('#ifdef USEMPI') == 2, out   # one wrapper per statement
    assert _guard_balance(out) == (0, 0), out


def test_add_uses_does_not_inject_into_a_left_in_place_guard():
    """An unconditional import added to a unit whose only other `use` sits
    inside a guard left in place must land *outside* that guard — otherwise
    it is compiled only when the guard happens to be satisfied."""
    tree       = _parse_with_module_uses(_SPANNING_GUARD)
    subroutine = next(n for n in walk_tree(tree) if n.get('type') == 'subroutine')
    add_uses(subroutine, {
        'moduleUse':   {'OMP_Lib': {'intrinsic': False, 'all': True}},
        'moduleOrder': ['OMP_Lib'],
    })

    out       = serialize(tree)
    unguarded = re.sub(r'#ifdef USEMPI.*?#endif', '', out, flags=re.S)
    assert 'OMP_Lib' in unguarded, out
    assert _guard_balance(out) == (0, 0), out


def test_guard_closing_outside_the_code_node_is_not_claimed():
    """A guard opened in one code node and closed in another (here, across a
    `contains`) cannot be closed by a rebuilt block, so it is not claimed."""
    source = (
        "module foo\n"
        "#ifdef USEMPI\n"
        "  use :: mpi\n"
        "  implicit none\n"
        "contains\n"
        "  subroutine s()\n"
        "    return\n"
        "  end subroutine s\n"
        "#endif\n"
        "end module foo\n"
    )
    tree = _parse_with_module_uses(source)
    assert 'conditions' not in _module_use_node(tree)['moduleUse']['mpi'][0]

    out = _rebuild(tree)
    assert out.count('#ifdef USEMPI') == 1, out
    assert out.count('#endif')        == 1, out


def test_directive_that_cannot_be_rebuilt_is_left_in_place():
    """`#include` (and `#define`, `#if …`, …) cannot be reconstructed from an
    entry's conditions, so it must not be absorbed into the block — a rebuild
    would silently delete it."""
    source = (
        "subroutine s()\n"
        "#include \"os.inc\"\n"
        "  use :: Error, only : Error_Report\n"
        "  implicit none\n"
        "end subroutine s\n"
    )
    tree = _parse_with_module_uses(source)
    out  = _rebuild(tree)
    assert '#include "os.inc"' in out, out


def test_if_expression_guard_is_not_claimed():
    """`#if defined(…)` has no `conditions` representation, so neither it nor
    any guard enclosing it may be claimed."""
    source = (
        "subroutine s()\n"
        "#ifdef USEMPI\n"
        "#if defined(HAVE_F08)\n"
        "  use :: mpi_f08\n"
        "#endif\n"
        "#endif\n"
        "  implicit none\n"
        "end subroutine s\n"
    )
    tree = _parse_with_module_uses(source)
    assert 'conditions' not in _module_use_node(tree)['moduleUse']['mpi_f08'][0]

    out = _rebuild(tree)
    assert out.count('#ifdef USEMPI')        == 1, out
    assert out.count('#if defined(HAVE_F08)') == 1, out
    assert out.count('#endif')                == 2, out
    assert _guard_balance(out) == (0, 0), out


def test_nested_guards_over_uses_only_are_claimed():
    """Nesting is fine as long as every region holds only `use` statements."""
    source = (
        "subroutine s()\n"
        "#ifdef USEMPI\n"
        "  use :: mpi\n"
        "#ifdef DEBUGGING\n"
        "  use :: mpi_debug\n"
        "#endif\n"
        "#endif\n"
        "  implicit none\n"
        "end subroutine s\n"
    )
    tree = _parse_with_module_uses(source)
    node = _module_use_node(tree)
    assert node['moduleUse']['mpi_debug'][0]['conditions'] == [
        {'name': 'USEMPI',    'invert': False},
        {'name': 'DEBUGGING', 'invert': False},
    ]

    out = _rebuild(tree)
    assert _guard_balance(out) == (0, 0), out
    # `mpi_debug` is still doubly guarded after the rebuild.
    assert re.search(r'#ifdef USEMPI\n#ifdef DEBUGGING\n\s*use.*mpi_debug', out), out


def test_stray_endif_is_not_swallowed_by_a_following_use():
    """An `#endif` closing a region opened in an earlier node used to be
    buffered and absorbed by the next `use`, which the rebuild then dropped."""
    source = (
        "module foo\n"
        "  implicit none\n"
        "contains\n"
        "  subroutine s()\n"
        "#endif\n"
        "    use :: Error, only : Error_Report\n"
        "    implicit none\n"
        "  end subroutine s\n"
        "end module foo\n"
    )
    tree = _parse_with_module_uses(source)
    out  = _rebuild(tree)
    assert out.count('#endif') == 1, out


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

def test_openmp_conditional_use_alignment():
    """`use` starts at the block's indent on every statement; the `!$ `
    sentinel sits to its left and the padding that realigns `::` goes after
    the module attributes, not before `use`."""
    source = (
        "module foo\n"
        "    use :: module1, only : a\n"
        "    !$ use :: module2, only : b\n"
        "    use :: module3, only : c\n"
        "    implicit none\n"
        "end module foo\n"
    )
    tree = _parse_with_module_uses(source)
    update_uses(_module_use_node(tree))
    lines = [line for line in
             _module_use_node(tree)['firstChild']['content'].splitlines()]

    assert lines == [
        "    use    :: module1, only : a",
        "    !$ use :: module2, only : b",
        "    use    :: module3, only : c",
    ], lines


def test_openmp_conditional_use_alignment_with_intrinsic():
    """The same, with an intrinsic module in the block: the `::` column is
    shared, and no statement is pushed off the block's indent."""
    source = (
        "module foo\n"
        "    use, intrinsic :: iso_c_binding, only : c_int\n"
        "    !$ use :: omp_lib, only : omp_get_thread_num\n"
        "    use :: module1, only : a\n"
        "    implicit none\n"
        "end module foo\n"
    )
    tree = _parse_with_module_uses(source)
    update_uses(_module_use_node(tree))
    lines = _module_use_node(tree)['firstChild']['content'].splitlines()

    assert all(line.startswith("    use") or line.startswith("    !$ use")
               for line in lines), lines
    columns = {line.index("::") for line in lines}
    assert len(columns) == 1, lines


def test_openmp_conditional_use_layout_is_idempotent():
    """Reformatting an already-formatted block must not shift it again."""
    source = (
        "module foo\n"
        "    use :: module1, only : a\n"
        "    !$ use :: module2, only : b\n"
        "    implicit none\n"
        "end module foo\n"
    )
    tree = _parse_with_module_uses(source)
    once = _rebuild(tree)
    twice = _rebuild(_parse_with_module_uses(once))
    assert twice == once, (once, twice)


def test_guard_whose_else_branch_holds_code_is_not_claimed():
    """The code need not follow the `use` inside the *same* branch: an
    `#else` branch holding anything but `use` statements disqualifies the
    whole region, since the rebuilt block would close the guard before the
    original `#endif`."""
    source = (
        "subroutine s()\n"
        "#ifdef USEMPI\n"
        "  use :: mpi\n"
        "#else\n"
        "  implicit none\n"
        "#endif\n"
        "  return\n"
        "end subroutine s\n"
    )
    tree = _parse_with_module_uses(source)
    assert 'conditions' not in _module_use_node(tree)['moduleUse']['mpi'][0]

    out = _rebuild(tree)
    assert out.count('#ifdef USEMPI') == 1, out
    assert out.count('#else')         == 1, out
    assert out.count('#endif')        == 1, out
    assert _guard_balance(out) == (0, 0), out
