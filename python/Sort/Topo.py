# Provides topological sort functionality.
# Andrew Benson (ported to Python 2026)
#
# Wraps Python's standard-library graphlib.TopologicalSorter (Python >= 3.9),
# exposing an interface that matches the Perl Sort::Topo module.
from __future__ import annotations

from graphlib import TopologicalSorter, CycleError

__all__ = ['sort']


def sort(objects: list, dependencies: dict) -> list:
    """Topologically sort objects by their dependencies.

    Mirrors Perl Sort::Topo::sort().

    Parameters
    ----------
    objects : list
        Names of the objects to sort.
    dependencies : dict
        Maps each object name to a list of names it depends on (i.e. names
        that must appear *before* it in the output).

    Returns
    -------
    list
        A sorted copy of objects where every dependency appears before the
        object that depends on it.

    Raises
    ------
    RuntimeError
        If a circular dependency is detected (matches the Perl 'die').
    """
    objects_set = set(objects)
    ts: TopologicalSorter = TopologicalSorter()
    for obj in objects:
        deps = [d for d in dependencies.get(obj, []) if d in objects_set]
        ts.add(obj, *deps)
    try:
        return list(ts.static_order())
    except CycleError:
        raise RuntimeError("circular dependency")
