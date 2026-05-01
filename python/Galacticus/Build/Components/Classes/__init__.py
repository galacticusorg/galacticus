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


from Galacticus.Build.Components.Utils         import (
    register,
    verbosity_level,
    _component_properties,
)
from Galacticus.Build.Components.DataTypes     import (
    data_object_primitive_name,
    data_object_doc_name,
)
from Galacticus.Build.Components.NullFunctions import create_null_function


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


def Build_Component_Classes(build):
    """Define one `nodeComponent<Class>` Fortran type per component class.

    Mirrors `Build_Component_Classes` (Classes.pm:75-217).  For every
    `(implementation, property)` pair under each class that has at
    least one of `isGettable` / `isSettable` / `isEvolvable`, emit
    null-function bindings for `Set` / `Count` / `Rate` / `Analytic`
    / `Inactive` / `Scale` as appropriate.  Each function name is
    deduplicated per-class via `properties_created`.
    """
    components        = build.get('components')        or {}
    component_classes = build.get('componentClasses')  or {}

    for class_name in sorted(component_classes.keys()):
        class_dict = component_classes[class_name]
        properties_created = set()
        type_bound_functions = []

        for impl_name in class_dict.get('memberNames') or []:
            impl_id = _ucfirst(class_name) + _ucfirst(impl_name)
            implementation = components.get(impl_id) or {}

            for prop in _component_properties(implementation):
                attrs = prop.get('attributes') or {}
                if not (attrs.get('isGettable')
                        or attrs.get('isSettable')
                        or attrs.get('isEvolvable')):
                    continue

                # Set function.
                if attrs.get('isSettable'):
                    fn_name = prop['name'] + 'Set'
                    if fn_name not in properties_created:
                        bound_to = create_null_function(build, {
                            'selfType':  class_name,
                            'attribute': 'set',
                            'property':  prop,
                            'intent':    'inout',
                        })
                        type_bound_functions.append({
                            'type':        'procedure',
                            'name':        fn_name,
                            'function':    bound_to,
                            'returnType':  r"\void",
                            'arguments':   (
                                data_object_doc_name(prop) + r"\ value"
                            ),
                            'description': (
                                f"Set the \\mono{{{prop['name']}}} property "
                                f"of the \\mono{{{class_name}}} component."
                            ),
                        })
                        properties_created.add(fn_name)

                # Evolve functions.
                if attrs.get('isEvolvable'):
                    # Count function.
                    fn_name = prop['name'] + 'Count'
                    if fn_name not in properties_created:
                        type_bound_functions.append({
                            'type':        'procedure',
                            'name':        fn_name,
                            'function':    create_null_function(build, {
                                'selfType':  class_name,
                                'attribute': 'get',
                                'property':  {'type': 'integer', 'rank': 0},
                                'intent':    'in',
                            }),
                            'returnType':  r"\intzero",
                            'arguments':   "",
                            'description': (
                                f"Compute the count of evolvable quantities in the "
                                f"\\mono{{{prop['name']}}} property of the "
                                f"\\mono{{{impl_id}}} component."
                            ),
                        })
                        properties_created.add(fn_name)

                    # Rate function (skipped for createIfNeeded).
                    fn_name = prop['name'] + 'Rate'
                    if fn_name not in properties_created:
                        if not attrs.get('createIfNeeded'):
                            type_bound_functions.append({
                                'type':        'procedure',
                                'name':        fn_name,
                                'function':    create_null_function(build, {
                                    'selfType':  class_name,
                                    'attribute': 'rate',
                                    'property':  prop,
                                    'intent':    'inout',
                                }),
                                'returnType':  r"\void",
                                'arguments':   (
                                    data_object_doc_name(prop) + r"\ value"
                                ),
                                'description': (
                                    f"Cumulate to the rate of the "
                                    f"\\mono{{{prop['name']}}} property of the "
                                    f"\\mono{{{impl_id}}} component."
                                ),
                            })

                        # Analytic / Inactive / Scale (only if not virtual).
                        if not attrs.get('isVirtual'):
                            type_bound_functions.append({
                                'type':        'procedure',
                                'name':        prop['name'] + 'Analytic',
                                'function':    create_null_function(build, {
                                    'selfType':  class_name,
                                    'attribute': 'analytic',
                                    'property':  prop,
                                    'intent':    'inout',
                                }),
                                'returnType':  r"\void",
                                'arguments':   "",
                                'description': (
                                    f"Mark the \\mono{{{prop['name']}}} property of "
                                    f"the \\mono{{{impl_id}}} component as "
                                    f"analtyically-solvable."
                                ),
                            })
                            type_bound_functions.append({
                                'type':        'procedure',
                                'name':        prop['name'] + 'Inactive',
                                'function':    create_null_function(build, {
                                    'selfType':  class_name,
                                    'attribute': 'inactive',
                                    'property':  prop,
                                    'intent':    'inout',
                                }),
                                'returnType':  r"\void",
                                'arguments':   "",
                                'description': (
                                    f"Mark the \\mono{{{prop['name']}}} property of "
                                    f"the \\mono{{{impl_id}}} component as inactive."
                                ),
                            })
                            type_bound_functions.append({
                                'type':        'procedure',
                                'name':        prop['name'] + 'Scale',
                                'function':    create_null_function(build, {
                                    'selfType':  class_name,
                                    'attribute': 'scale',
                                    'property':  prop,
                                    'intent':    'inout',
                                }),
                                'returnType':  r"\void",
                                'arguments':   (
                                    data_object_doc_name(prop) + r"\ value"
                                ),
                                'description': (
                                    f"Set the scale of the \\mono{{{prop['name']}}} "
                                    f"property of the \\mono{{{impl_id}}} component."
                                ),
                            })
                        properties_created.add(fn_name)

        type_name = 'nodeComponent' + _ucfirst(class_name)
        build.setdefault('types', {})[type_name] = {
            'name':           type_name,
            'comment':        (
                f"Type for the \\mono{{{class_name}}} component class."
            ),
            'isPublic':       True,
            'extends':        'nodeComponent',
            'boundFunctions': type_bound_functions,
        }


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


# ---------------------------------------------------------------------------
# Hook registration.  Order matches Perl Classes.pm:21-28.
# ---------------------------------------------------------------------------

register('classes', 'gather', Gather_Classes)
register('classes', 'types',  Build_Component_Classes)
