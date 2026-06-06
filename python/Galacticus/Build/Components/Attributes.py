"""Validation and default-population hooks for component property
attributes.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Components/Attributes.pm.  Registers four
pipeline hooks:

  preValidate  → Validate_Deferreds_Functionless
  default      → Default_Functions
  postValidate → Validate_Boolean, Validate_Evolvable_Intrinsics
"""
from __future__ import annotations



from Galacticus.Build.Components.Utils import (
    register,
    is_intrinsic,
    _component_properties,
)


def Validate_Deferreds_Functionless(build: dict) -> None:
    """Forbid build-time `xxxFunction` elements on a property whose `xxx`
    method is also flagged deferred.  Mirrors `Validate_Deferreds_Functionless`.
    """
    for component in (build.get('components') or {}).values():
        if 'properties' not in component:
            continue
        for prop in _component_properties(component):
            attrs = prop.get('attributes') or {}
            deferred_field = attrs.get('isDeferred')
            if not deferred_field:
                continue
            for method in deferred_field.split(':'):
                if (method + 'Function') in prop:
                    raise ValueError(
                        "Validate_Deferreds_Functionless(): cannot specify '"
                        f"{method}Function' when '{method}' method is "
                        f"deferred for property '{prop['name']}' of component "
                        f"'{component['class']}{_ucfirst(component['name'])}'"
                    )


def Default_Functions(build: dict) -> None:
    """Fill in default `rateFunction` / `getFunction` / `setFunction`
    entries on every component property.

    Mirrors `Default_Functions`.  Each `getFunction` / `setFunction` is
    normalised to a dict with `content` (the Fortran symbol name) and
    `build` (True if the build system needs to emit a body, False if the
    user supplied one in the XML).
    """
    for implementation in (build.get('components') or {}).values():
        component_id = (
            _ucfirst(implementation['class'])
            + _ucfirst(implementation['name'])
        )
        if 'properties' not in implementation:
            continue
        for prop in _component_properties(implementation):
            prop.setdefault(
                'rateFunction',
                component_id + _ucfirst(prop['name']) + "Rate",
            )
            for verb in ('get', 'set'):
                key = verb + 'Function'
                if key in prop:
                    if not isinstance(prop[key], dict):
                        prop[key] = {'content': prop[key]}
                    prop[key]['build'] = False
                else:
                    prop[key] = {
                        'content': (
                            _lcfirst(component_id)
                            + _ucfirst(prop['name'])
                            + _ucfirst(verb)
                        ),
                        'build':   True,
                    }


def Validate_Boolean(build: dict) -> None:
    """Require `isSettable` / `isGettable` / `isEvolvable` to be `"true"`
    or `"false"`, and convert each to a Python `bool`.

    Mirrors `Validate_Boolean`.
    """
    for component in (build.get('components') or {}).values():
        if 'properties' not in component:
            continue
        for prop in _component_properties(component):
            attrs = prop.get('attributes')
            if not isinstance(attrs, dict):
                continue
            for key in ('isSettable', 'isGettable', 'isEvolvable'):
                if key not in attrs:
                    continue
                value = attrs[key]
                if value == 'true':
                    attrs[key] = True
                elif value == 'false':
                    attrs[key] = False
                elif isinstance(value, bool):
                    # Already coerced (e.g. Default_Functions ran twice).
                    pass
                else:
                    raise ValueError(
                        f"Validate_Boolean: value of '{key}' attribute of '"
                        f"{prop['name']}' of '"
                        f"{component['class']}{_ucfirst(component['name'])}' "
                        "must be either 'true' or 'false'"
                    )


def Validate_Evolvable_Intrinsics(build: dict) -> None:
    """Forbid evolvable intrinsic properties whose type is not `double`.

    Mirrors `Validate_Evolvable_Intrinsics`.
    """
    for component in (build.get('components') or {}).values():
        if 'properties' not in component:
            continue
        for prop in _component_properties(component):
            attrs = prop.get('attributes') or {}
            if (is_intrinsic(prop.get('type'))
                    and prop['type'] != 'double'
                    and attrs.get('isEvolvable')):
                raise ValueError(
                    "Validate_Evolvable_Intrinsics: non-real intrinsic "
                    f"property '{prop['name']}' of '"
                    f"{component['class']}{_ucfirst(component['name'])}' "
                    "component can not be evolvable"
                )


# ---------------------------------------------------------------------------
# Hook registration.  Phase ordering matches Perl.
# ---------------------------------------------------------------------------

register('attributes', 'preValidate',  Validate_Deferreds_Functionless)
register('attributes', 'default',      Default_Functions)
register('attributes', 'postValidate', Validate_Boolean)
register('attributes', 'postValidate', Validate_Evolvable_Intrinsics)


def _ucfirst(text: str) -> str:
    return text[:1].upper() + text[1:] if text else text


def _lcfirst(text: str) -> str:
    return text[:1].lower() + text[1:] if text else text
