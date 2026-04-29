# Provides topological sort functionality.
# Andrew Benson (ported to Python 2026)
#
# Faithful port of perl/Sort/Topo.pm's home-grown topological sort, NOT
# a graphlib wrapper.  The Perl module's behaviour and its doc-comment
# disagree about the meaning of `dependencies[X] = [Y]`:
#
#   - The doc-comment claims: X depends on Y → Y comes before X.
#   - The actual algorithm:   Y depends on X → Y comes AFTER X.
#
# Galacticus relies on the actual algorithm semantics throughout
# (for example `Implementation_Dependencies` builds
# `dependencies[parent] += [child]` to make the child come after its
# extends-parent), and the specific tie-breaking order Perl produces
# affects which iteration of `Build_Component_Classes` first claims
# the `propertiesCreated` slot for a property shared between members.
# A graphlib-based port produced a different but equally-valid topo
# order, which silently changed downstream output.  This faithful
# port preserves Perl's exact ordering.


def sort(objects, dependencies):
    """Topologically sort objects by their dependencies.

    Mirrors Perl Sort::Topo::sort() — the algorithm at
    perl/Sort/Topo.pm:7-78, based on the example at
    https://rosettacode.org/wiki/Topological_sort .

    Parameters
    ----------
    objects : list
        Names of the objects to sort.
    dependencies : dict
        Maps each object name to a list of names that *depend on* it —
        i.e. names that must appear AFTER the key in the output.
        (Despite its name, this is the inverse of the relationship
        suggested by Perl's misleading doc-comment; Galacticus relies
        on the actual algorithm semantics.)

    Returns
    -------
    list
        A topologically-sorted copy of `objects` where for every
        `(key, value)` pair in `dependencies`, every entry of `value`
        appears AFTER `key`.

    Raises
    ------
    RuntimeError
        If a circular dependency is detected.
    """
    objects = list(objects)
    n_objects = len(objects)

    # Build the dependency list: each entry is [dependentIndex, objectIndex].
    # Perl iterates the dependency keys in sorted order (line 16) so the
    # deterministic tie-breaking matches.
    dep_pairs = []
    for object_name in sorted(dependencies.keys()):
        if object_name not in objects:
            continue
        object_index = objects.index(object_name)
        for dependent in dependencies[object_name]:
            if dependent not in objects:
                continue
            dependent_index = objects.index(dependent)
            dep_pairs.append([dependent_index, object_index])

    n_deps = len(dep_pairs)
    order    = list(range(n_objects))
    position = list(range(n_objects))

    # Bubble-from-the-end algorithm, ported faithfully.
    k = 0
    j = None
    while True:
        j = k
        k = n_objects
        for i in range(n_deps):
            dep_left  = dep_pairs[i][0]
            dep_right = dep_pairs[i][1]
            pos_left  = position[dep_left]
            pos_right = position[dep_right]
            if (
                dep_left == dep_right
                or pos_left  >= k
                or pos_left  <  j
                or pos_right <  j
            ):
                continue
            k -= 1
            position[order[k]] = pos_left
            position[dep_left] = k
            order[pos_left]    = order[k]
            order[k]           = dep_left
        if k <= j:
            break

    # Cycle detection: any unsorted suffix indicates a cycle.
    if j < n_objects:
        raise RuntimeError("circular dependency")

    return [objects[order[i]] for i in range(n_objects)]
