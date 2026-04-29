# Components/Implementations — per-implementation generators.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Implementations.pm.  This file
# currently ports only the small hooks needed to unblock other stages
# (`Implementation_ID_List` plus the two `default`-phase hooks that
# populate `defaultImplementation` and synthesise null components).  The
# remainder of Implementations.pm — `Implementation_Defaults` body
# extension, `Implementation_Dependencies`, `Implementation_Parents`,
# `Build_Component_Implementations`, `Default_Full_Name` — comes with
# the full Implementations/ port.

import os
import sys

sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))

from Galacticus.Build.Components.Utils import (
    register,
    apply_defaults,
    verbosity_level,
)


def Default_Full_Name(build):
    """Set `fullyQualifiedName = <Class><Name>` on every component.

    Mirrors `Default_Full_Name`.
    """
    for impl in (build.get('components') or {}).values():
        impl['fullyQualifiedName'] = (
            _ucfirst(impl['class']) + _ucfirst(impl['name'])
        )


def Implementation_ID_List(build):
    """Populate `build['componentIdList']` with a sorted list of component
    implementation IDs (e.g. `["BasicStandard", "BlackHoleSimple", …]`).

    Mirrors `Implementation_ID_List`.
    """
    components = build.get('components') or {}
    build['componentIdList'] = sorted(components.keys())


# Default-attribute table for each implementation.  Mirrors the
# `%defaults` hash inside `Implementation_Defaults` (Implementations.pm:64+).
# We start with the small subset Hierarchy.pm relies on; the full table
# (bindings / interfaces / properties etc.) will land with the rest of
# Implementations.pm.
_IMPLEMENTATION_DEFAULTS = {
    'isDefault': 'booleanFalse',
}


def Implementation_Defaults(build):
    """Apply per-implementation default attributes and record the default
    implementation name on each class.

    Mirrors the relevant top portion of `Implementation_Defaults`.
    """
    for impl in (build.get('components') or {}).values():
        for key, default in _IMPLEMENTATION_DEFAULTS.items():
            apply_defaults(impl, key, default)

    # Record `defaultImplementation` on each component class.
    component_classes = build.setdefault('componentClasses', {})
    for impl in (build.get('components') or {}).values():
        if impl.get('isDefault'):
            component_classes.setdefault(impl['class'], {})['defaultImplementation'] = impl['name']


def Null_Implementations(build):
    """Synthesise a `null` implementation for any class that lacks one.

    Mirrors `Null_Implementations`.  When a class has no `null` member,
    we add one; if no other implementation was marked `isDefault`, the
    new null member becomes the default.
    """
    components        = build.setdefault('components', {})
    component_classes = build.setdefault('componentClasses', {})

    classes = {}
    for impl in components.values():
        c = classes.setdefault(impl['class'], {'hasNull': False, 'hasDefault': False})
        if impl.get('name') == 'null':
            c['hasNull'] = True
        if impl.get('isDefault'):
            c['hasDefault'] = True

    for class_name in sorted(classes.keys()):
        info = classes[class_name]
        if info['hasNull']:
            continue
        impl_name = _ucfirst(class_name) + 'Null'
        is_default = not info['hasDefault']
        components[impl_name] = {
            'class':     class_name,
            'name':      'null',
            'isDefault': is_default,
        }
        if is_default:
            component_classes.setdefault(class_name, {})['defaultImplementation'] = 'null'
        # Append to the ID list (Implementation_ID_List ran in
        # `preValidate`, so we keep the trailing-append pattern).
        build.setdefault('componentIdList', []).append(impl_name)
        if verbosity_level >= 1:
            prefix = "         --> Adding null implementation "
            qualifier = "as default " if not info['hasDefault'] else ""
            print(f"{prefix}{qualifier}for {class_name} class")


# ---------------------------------------------------------------------------
# Hook registration
# ---------------------------------------------------------------------------

# `Implementation_ID_List` registers under `preValidate` so it runs before
# `default`-phase hooks that walk `components` and may mutate the dict.
register('implementations', 'preValidate', Implementation_ID_List)
register('implementations', 'default',     Implementation_Defaults)
register('implementations', 'default',     Null_Implementations)
register('implementations', 'default',     Default_Full_Name)


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text
