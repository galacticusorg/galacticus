"""Topological sort of tasks described by `after` / `before` relations.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Dependencies.pm — specifically Dependency_Sort,
in its minimal subset used by SourceTree.Process.EventHooksStatic: no
`sortName` aliasing, no `re:` regex keys.
"""
from __future__ import annotations

from Sort.Topo        import sort as topo_sort
from List.ExtraUtils  import as_array

__all__ = ['dependency_sort']


def dependency_sort(tasks: dict[str, dict]) -> list[str]:
    """Sort task names into dependency order.

    Parameters
    ----------
    tasks : dict[str, dict]
        Each value may declare an `after` and/or `before` key whose value is
        a task name (str) or a list of task names.

    Returns
    -------
    list[str]
        Task names in an order consistent with what Perl `Dependency_Sort`
        actually emits (which is the opposite of what its docstring says
        — see the Notes below).

    Notes
    -----
    Reading is the natural one:
      - `X.after  = Y` →  X comes AFTER  Y in the output (Y precedes X).
      - `X.before = Y` →  X comes BEFORE Y in the output (X precedes Y).

    Implementation note (subtle).  Perl's `Sort::Topo` documents its
    `dependencies` argument as "key X → array of names X depends on",
    but its actual algorithm emits the *array element* later than the
    key (so `dependencies[X] = [Y]` makes X precede Y in the output —
    the opposite of what its docstring says).  Perl `Dependency_Sort`
    relied on that inversion: for `X.after = Y` it pushed
    `dependencies[Y] += [X]`, and the broken docstring + broken
    algorithm cancelled out to give the right answer.

    Our Python `Sort.Topo` wraps `graphlib.TopologicalSorter`, which
    matches the standard "X depends on Y → Y first" convention.  So
    we build edges the way the standard reading dictates:
      - `X.after  = Y` → `deps[X] += [Y]` ("X depends on Y", Y first).
      - `X.before = Y` → `deps[Y] += [X]` ("Y depends on X", X first).
    """
    names = list(tasks.keys())
    deps: dict[str, list[str]] = {name: [] for name in names}

    for name, info in tasks.items():
        for after in as_array(info.get('after')):
            if after in deps:
                deps[name].append(after)
        for before in as_array(info.get('before')):
            if before in deps:
                deps[before].append(name)

    return topo_sort(sorted(names), deps)
