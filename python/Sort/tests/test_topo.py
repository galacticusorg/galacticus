"""Tests for `Sort.Topo.sort` (port of Perl Sort::Topo).

Recent history (commits 4ca9596d, 62314942) shows the algorithm direction
was toggled twice during the port -- the canonical ordering convention is
subtle ("dependencies[X] = [Y]" reads as "X depends on Y, so Y emits
first"), and the related Perl module's docstring/algorithm mismatch
masked the issue for years.  These tests pin down the current convention
so any future regression is caught.
"""

import pytest

from Sort.Topo import sort


def _index(seq, item):
    """Return the position of `item` in `seq` (helper for ordering asserts)."""
    return seq.index(item)


# ---------------------------------------------------------------------------
# Convention: dependencies[X] = [Y, ...]  means  Y must come BEFORE X.
# ---------------------------------------------------------------------------

def test_no_dependencies_returns_all_objects():
    """With no dependencies, every input object appears in the output exactly
    once.  Order between independent nodes is not asserted -- graphlib's
    static_order is deterministic but the contract here is set-equality."""
    out = sort(['a', 'b', 'c'], {})
    assert set(out) == {'a', 'b', 'c'}
    assert len(out) == 3


def test_simple_dependency_orders_dependency_first():
    """`{B: [A]}` reads as 'B depends on A', so A must precede B."""
    out = sort(['A', 'B'], {'B': ['A']})
    assert _index(out, 'A') < _index(out, 'B')


def test_chain_dependency_orders_transitively():
    """A -> B -> C: dependencies[B]=[A], dependencies[C]=[B] yields A, B, C."""
    out = sort(['A', 'B', 'C'], {'B': ['A'], 'C': ['B']})
    assert _index(out, 'A') < _index(out, 'B') < _index(out, 'C')


def test_diamond_dependency_orders_root_first():
    """A -> {B, C} -> D.  A first, D last; B and C between in either order."""
    out = sort(
        ['A', 'B', 'C', 'D'],
        {'B': ['A'], 'C': ['A'], 'D': ['B', 'C']},
    )
    assert _index(out, 'A') < _index(out, 'B')
    assert _index(out, 'A') < _index(out, 'C')
    assert _index(out, 'B') < _index(out, 'D')
    assert _index(out, 'C') < _index(out, 'D')


def test_dependencies_outside_objects_set_are_ignored():
    """Names referenced in `dependencies` that aren't in `objects` must NOT
    be added to the output, and must NOT cause an error.  Mirrors
    Sort::Topo's tolerance of dangling references in `static_link_dependencies`-
    style configurations."""
    out = sort(['A', 'B'], {'B': ['A', 'X', 'Y']})
    assert set(out) == {'A', 'B'}
    assert _index(out, 'A') < _index(out, 'B')


def test_object_with_no_entry_in_dependencies_dict_is_handled():
    """An object that is not a key of `dependencies` is treated as having no
    dependencies."""
    out = sort(['A', 'B', 'C'], {'B': ['A']})
    assert set(out) == {'A', 'B', 'C'}
    assert _index(out, 'A') < _index(out, 'B')


def test_empty_objects_returns_empty_list():
    assert sort([], {}) == []
    assert sort([], {'B': ['A']}) == []


def test_self_dependency_raises_RuntimeError():
    """`{A: [A]}` is a degenerate cycle; the wrapper translates graphlib's
    CycleError into a `RuntimeError("circular dependency")` to match the
    Perl `die`."""
    with pytest.raises(RuntimeError, match="circular"):
        sort(['A'], {'A': ['A']})


def test_two_node_cycle_raises_RuntimeError():
    """A <-> B forms a cycle and must error rather than silently producing
    a partial order."""
    with pytest.raises(RuntimeError, match="circular"):
        sort(['A', 'B'], {'A': ['B'], 'B': ['A']})


def test_three_node_cycle_raises_RuntimeError():
    with pytest.raises(RuntimeError, match="circular"):
        sort(
            ['A', 'B', 'C'],
            {'A': ['B'], 'B': ['C'], 'C': ['A']},
        )


def test_objects_not_mutated():
    """The input list/dict must not be modified in place (callers in the
    build system reuse them across invocations)."""
    objects = ['A', 'B', 'C']
    deps    = {'B': ['A'], 'C': ['B']}
    sort(objects, deps)
    assert objects == ['A', 'B', 'C']
    assert deps    == {'B': ['A'], 'C': ['B']}
