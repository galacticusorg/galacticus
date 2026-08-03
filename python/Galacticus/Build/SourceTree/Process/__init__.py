"""Process hook registry and tree-processing driver.

Andrew Benson (ported to Python 2026)

Each Process submodule registers itself at import time via
register_process() or register_postprocess().
"""



from Sort.Topo import sort as topo_sort


# Populated by Process submodules at import time.
PROCESS_HOOKS       = {}   # name -> callable(tree, options)
PROCESS_DEPENDENCIES = {}  # name -> list of hook names that should run AFTER name
POSTPROCESS_HOOKS   = {}   # name -> callable(tree, options)
ANALYZE_HOOKS       = {}   # name -> callable(tree, options)


def register_analyze(name, fn):
    """Register an analyze hook.

    Analyze hooks run in
    alphabetical order after process_tree() when analyze_tree() is invoked.
    """
    ANALYZE_HOOKS[name] = fn


def register_process(name, fn, before=()):
    """Register a process hook.

    Parameters
    ----------
    name : str
        Hook name; must be unique across all Process submodules.
    fn : callable(tree, options)
        Hook function invoked with the root node and an options dict.
    before : iterable of str
        Other hook names that MUST run after this one — i.e. each entry
        depends on `name`.
    """
    PROCESS_HOOKS[name] = fn
    if before:
        PROCESS_DEPENDENCIES[name] = list(before)


def register_postprocess(name, fn):
    """Register a post-process hook.

    Post-process hooks run
    in alphabetical order after every process hook has completed.
    """
    POSTPROCESS_HOOKS[name] = fn


def process_tree(tree, options=None):
    """Run every registered process hook on the tree, in dependency order.

    The hook-execution order is the topological sort of PROCESS_HOOKS using
    PROCESS_DEPENDENCIES as the precedence graph: if
    `PROCESS_DEPENDENCIES[X] = [Y, Z]`, X runs before Y and before Z.
    """
    if options is None:
        options = {}

    # Build the edge set that Sort.Topo.sort expects: edges[A] = list of
    # names that must precede A.
    edges = {}
    for name, dependents in PROCESS_DEPENDENCIES.items():
        for dependent in dependents:
            edges.setdefault(dependent, []).append(name)

    for name in topo_sort(sorted(PROCESS_HOOKS.keys()), edges):
        PROCESS_HOOKS[name](tree, options)

    for name in sorted(POSTPROCESS_HOOKS.keys()):
        POSTPROCESS_HOOKS[name](tree, options)

    return tree


# Wire up the one existing post-process hook (Parse.Directives).  This keeps
# callers who just invoke process_tree() from having to remember to run it
# themselves.
def _register_directives_postprocess():
    from Galacticus.Build.SourceTree.Parse.Directives import post_process_directives

    def _wrapper(tree, options):
        post_process_directives(tree)

    register_postprocess('directives', _wrapper)


_register_directives_postprocess()
