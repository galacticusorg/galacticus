# Process hook registry and tree-processing driver.
# Andrew Benson (ported to Python 2026)
#
# Mirrors the Galacticus::Build::SourceTree::ProcessTree() driver and the
# %processHooks / %processDependencies / %postprocessHooks registries from
# perl/Galacticus/Build/SourceTree.pm.  Each Process submodule registers
# itself at import time via register_process() or register_postprocess().

import os
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from Sort.Topo import sort as topo_sort


# Populated by Process submodules at import time.
PROCESS_HOOKS       = {}   # name -> callable(tree, options)
PROCESS_DEPENDENCIES = {}  # name -> list of hook names that should run AFTER name
POSTPROCESS_HOOKS   = {}   # name -> callable(tree, options)


def register_process(name, fn, before=()):
    """Register a process hook.

    Mirrors the Perl idiom
        $Hooks::processHooks{name}        = \&fn;
        $Hooks::processDependencies{name} = [before, ...];

    Parameters
    ----------
    name : str
        Hook name; must be unique across all Process submodules.
    fn : callable(tree, options)
        Hook function invoked with the root node and an options dict.
    before : iterable of str
        Other hook names that MUST run after this one — i.e. each entry
        depends on `name`.  Matches the meaning of Perl's
        `processDependencies{name} = [...]`.
    """
    PROCESS_HOOKS[name] = fn
    if before:
        PROCESS_DEPENDENCIES[name] = list(before)


def register_postprocess(name, fn):
    """Register a post-process hook.

    Mirrors `$Hooks::postprocessHooks{name} = \&fn`.  Post-process hooks run
    in alphabetical order after every process hook has completed.
    """
    POSTPROCESS_HOOKS[name] = fn


def process_tree(tree, options=None):
    """Run every registered process hook on the tree, in dependency order.

    Mirrors Perl Galacticus::Build::SourceTree::ProcessTree().

    The hook-execution order is the topological sort of PROCESS_HOOKS using
    PROCESS_DEPENDENCIES as the precedence graph: if
    `PROCESS_DEPENDENCIES[X] = [Y, Z]`, X runs before Y and before Z.
    """
    if options is None:
        options = {}

    # Build the edge set that Sort.Topo.sort expects: edges[A] = list of
    # names that must precede A.  Perl does exactly this (SourceTree.pm:190-195):
    #     push(@{$dependencies{$dependent}}, $processor);
    edges = {}
    for name, dependents in PROCESS_DEPENDENCIES.items():
        for dependent in dependents:
            edges.setdefault(dependent, []).append(name)

    for name in topo_sort(sorted(PROCESS_HOOKS.keys()), edges):
        PROCESS_HOOKS[name](tree, options)

    for name in sorted(POSTPROCESS_HOOKS.keys()):
        POSTPROCESS_HOOKS[name](tree, options)

    return tree


# Wire up the one existing post-process hook (Parse.Directives).  This matches
# the Perl registration at Parse/Directives.pm:14 and keeps callers who just
# invoke process_tree() from having to remember to run it themselves.
def _register_directives_postprocess():
    from Galacticus.Build.SourceTree.Parse.Directives import post_process_directives

    def _wrapper(tree, options):
        post_process_directives(tree)

    register_postprocess('directives', _wrapper)


_register_directives_postprocess()
