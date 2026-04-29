# Components/Classes — per-class generators.
# Andrew Benson (ported to Python 2026)
#
# Mirrors the perl/Galacticus/Build/Components/Classes/ directory.  The
# parent module Classes.pm contributes `Gather_Classes` (a `gather`-phase
# hook that groups component implementations by their class name and
# populates `componentClasses` / `componentClassList` /
# `componentClassListActive`) plus `Build_Component_Classes` (a `types`
# phase hook that emits one `nodeComponent<Class>` Fortran type per class).
# Only the small `Gather_Classes` piece is ported here for now — the type
# builder and the per-class iterated hooks come with the full Classes/
# port.

import os
import sys

sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))

from Galacticus.Build.Components.Utils import register, verbosity_level


def Gather_Classes(build):
    """Group component implementations by class name and build active-class
    lists.  Mirrors `Gather_Classes`.

    After this hook runs, `build` carries:

    * `componentClasses[<className>]` — dict with `name`, `memberNames`,
      and `members` keys.
    * `componentClassList` — sorted list of all class names.
    * `componentClassListActive` — same as `componentClassList`, unless
      the `GALACTICUS_ACTIVE_COMPONENTS` env var is set, in which case
      only listed names are kept.
    """
    component_classes = build.setdefault('componentClasses', {})

    for implementation in (build.get('components') or {}).values():
        class_name = implementation['class']
        impl_name  = implementation['name']
        # Implementation_Defaults may have already created the entry
        # carrying just `defaultImplementation`; fill in the rest
        # rather than overwrite.
        entry = component_classes.setdefault(class_name, {})
        entry.setdefault('name',        class_name)
        entry.setdefault('memberNames', [])
        entry.setdefault('members',     [])
        entry['memberNames'].append(impl_name)
        entry['members'    ].append(implementation)

    component_class_list = sorted(component_classes.keys())
    build['componentClassList'] = component_class_list

    active_env = os.environ.get('GALACTICUS_ACTIVE_COMPONENTS')
    if active_env is not None:
        active = active_env.split()
        build['componentClassListActive'] = [
            n for n in component_class_list if n in active
        ]
    else:
        build['componentClassListActive'] = list(component_class_list)

    if verbosity_level > 0:
        print("         --> Found the following component classes and implementations:")
        for class_name in component_class_list:
            print(f"            --> {class_name}")
            for impl_name in component_classes[class_name]['memberNames']:
                print(f"               --> {impl_name}")


register('classes', 'gather', Gather_Classes)
