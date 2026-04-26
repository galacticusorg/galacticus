# Regression tests for `Galacticus.Build.SourceTree.Parse.ModuleUses`.
#
# Two bug classes:
#
#   1. `_flush_code` used to absorb `raw_pp_buf` whenever flushing a code
#      buffer.  That broke the case where a `#ifdef USEMPI` line sat
#      immediately before a `use mpi` line: the `#ifdef` got copied into the
#      surrounding code node, and again into the moduleUse block when the
#      use was rebuilt — yielding a duplicated `#ifdef … #endif` wrapper
#      (one of which was unmatched and produced "unterminated #ifdef" at
#      compile time).
#
#   2. `add_uses` merged a new unconditional `use` onto an existing
#      conditional entry while preserving the existing entry's `conditions`
#      list.  This meant code that simply asked to `use mpiSelf` ended up
#      only importing it under `#ifdef USEMPI`, and was undefined in serial
#      builds.

from Galacticus.Build.SourceTree                  import parse_code, serialize
from Galacticus.Build.SourceTree.Parse.ModuleUses import parse_module_uses, add_uses


def _parse_with_module_uses(source):
    """parse_code already runs `_pass_module_uses` for us, but `parse_code`
    expects either a real file or pre-instrumented content.  Use the
    documented `instrument=False` form so the test source is taken
    verbatim."""
    return parse_code(source, name='<test>', instrument=False)


def test_pp_guarded_use_does_not_duplicate_when_rebuilt():
    """When `update_uses` rebuilds the moduleUse content (which `add_uses`
    triggers), the `#ifdef USEMPI` wrapper must come from exactly ONE
    place — the rebuild — not from a surrounding code node that absorbed
    the original `#ifdef` line.  An earlier draft of `_flush_code`
    extended `raw_code_buf` with `raw_pp_buf` whenever it ran, so the
    `#ifdef` ended up duplicated between the code node and the rebuild,
    yielding "unterminated #ifdef" at compile time."""
    from Galacticus.Build.SourceTree.Parse.ModuleUses import update_uses
    from Galacticus.Build.SourceTree                  import walk_tree

    source = (
        "module foo\n"
        "#ifdef USEMPI\n"
        "  use mpi\n"
        "#endif\n"
        "  implicit none\n"
        "end module foo\n"
    )
    tree = _parse_with_module_uses(source)

    module_use_node = next(
        n for n in walk_tree(tree) if n.get('type') == 'moduleUse')
    update_uses(module_use_node)

    out = serialize(tree)
    assert out.count('#ifdef USEMPI') == 1, out
    assert out.count('#endif')        == 1, out
    # The `use mpi` survives.
    import re as _re
    assert _re.search(r'\buse\b\s*(::)?\s*mpi\b', out), out


def test_unconditional_add_uses_drops_existing_conditions():
    """If the tree already imports `mpi` under `#ifdef USEMPI`, an unconditional
    `add_uses(...)` for `mpi` must drop the wrapper — otherwise `mpiSelf` is
    only available in MPI builds."""
    source = (
        "module foo\n"
        "#ifdef USEMPI\n"
        "  use mpi\n"
        "#endif\n"
        "  implicit none\n"
        "end module foo\n"
    )
    tree = _parse_with_module_uses(source)

    # Find the moduleUse node.
    from Galacticus.Build.SourceTree import walk_tree
    module_use_node = next(
        n for n in walk_tree(tree) if n.get('type') == 'moduleUse')

    # Sanity: the existing entry has a `conditions` list.
    assert 'conditions' in module_use_node['moduleUse']['mpi']

    # Now add an unconditional `use mpi`.  The new node carries no
    # `conditions`, mirroring the call sites that simply want the module
    # imported in every build.
    add_uses(module_use_node['parent'], {
        'moduleUse':   {'mpi': {'openMP': False, 'intrinsic': False, 'all': True}},
        'moduleOrder': ['mpi'],
    })

    out = serialize(tree)
    # The conditional wrapper is gone.
    assert '#ifdef USEMPI' not in out
    assert '#endif'        not in out
    # `use mpi` survives in some emit form (`use mpi` or `use :: mpi`).
    import re as _re
    assert _re.search(r'\buse\b\s*(::)?\s*mpi\b', out), out


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

    from Galacticus.Build.SourceTree import walk_tree
    module_use_node = next(
        n for n in walk_tree(tree) if n.get('type') == 'moduleUse')

    entry = module_use_node['moduleUse']['mpi']
    assert entry.get('conditions') == [{'name': 'USEMPI', 'invert': False}]
