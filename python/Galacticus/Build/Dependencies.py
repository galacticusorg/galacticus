# Topological sort of tasks described by `after` / `before` relations.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Dependencies.pm — specifically Dependency_Sort,
# in its minimal subset used by SourceTree.Process.EventHooksStatic: no
# `sortName` aliasing, no `re:` regex keys.

import os
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from Sort.Topo        import sort as topo_sort
from List.ExtraUtils  import as_array


def dependency_sort(tasks):
    """Sort task names into dependency order.

    Parameters
    ----------
    tasks : dict[str, dict]
        Each value may declare an `after` and/or `before` key whose value is
        a task name (str) or a list of task names.

    Returns
    -------
    list[str]
        Task names in an order consistent with Perl `Dependency_Sort`.

    Notes
    -----
    The Perl implementation builds edges as follows:
      - `X.after = Y` → Sort::Topo is told that **Y depends on X**, i.e. X
        precedes Y.
      - `X.before = Y` → Sort::Topo is told that **X depends on Y**, i.e. Y
        precedes X.
    Both conventions read "inverted" relative to plain English but are
    internally consistent and are what callers rely on.  We mirror them
    exactly here.
    """
    names = list(tasks.keys())
    deps  = {name: [] for name in names}

    for name, info in tasks.items():
        for after in as_array(info.get('after')):
            if after in deps:
                deps[after].append(name)
        for before in as_array(info.get('before')):
            if name in deps:
                deps[name].append(before)

    return topo_sort(sorted(names), deps)
