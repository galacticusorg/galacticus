# Driver for the `classIteratedFunctions` sub-iteration phase.
# Andrew Benson (ported to Python 2026)
#
# Mirrors the small portion of perl/Galacticus/Build/Components/Classes/Utils.pm
# that runs as part of the components-build pipeline:
# `Class_Function_Iterator`.  This is itself a `functions`-phase hook that,
# when invoked, walks the `component_utils` registry looking for
# `classIteratedFunctions` entries and dispatches each function once per
# component class.
#
# The remaining functions in the Perl Classes::Utils module
# (`Class_Move`, `Class_Remove`, the main `functions`-phase group) belong
# to a later port stage.

import os
import sys

sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))

from Galacticus.Build.Components.Utils import register, component_utils, verbosity_level


def Class_Function_Iterator(build):
    """Run every `classIteratedFunctions` hook once per component class.

    Mirrors `Class_Function_Iterator`.  Owners may register a list of
    functions under the `classIteratedFunctions` key; each function is
    called as `fn(build, class_dict)` for every entry in
    `build['componentClasses']`.
    """
    if verbosity_level >= 1:
        # Match Perl's leading-six-space report indent.
        pass

    for owner_name in sorted(component_utils.keys()):
        owner = component_utils[owner_name]
        functions = owner.get('classIteratedFunctions')
        if not functions:
            continue
        if not isinstance(functions, list):
            functions = [functions]
        for fn in functions:
            marker = (
                f" {{{getattr(fn, '__name__', '<fn>')}}}"
                if len(functions) > 1 else ''
            )
            print(f"         --> {owner_name}{marker}")
            for class_dict in (build.get('componentClasses') or {}).values():
                fn(build, class_dict)


register('classUtils', 'functions', Class_Function_Iterator)
