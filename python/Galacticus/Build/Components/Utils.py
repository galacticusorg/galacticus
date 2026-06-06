"""Shared helpers and global state for the components-build pipeline.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Components/Utils.pm.  Sub-modules register
themselves on the `component_utils` registry; the driver in
Galacticus.Build.Components walks that registry phase by phase.

Module-level globals reproduce Perl's `our $foo` declarations
(`class_name_length_max`, etc.).  `Label_Lengths` populates them during
the `gather` phase; later phases read them for column-aligned output.
"""
from __future__ import annotations

import logging

from typing import Any, Callable

logger = logging.getLogger(__name__)



# ---------------------------------------------------------------------------
# Hook registry
#
# component_utils['<owner>'] = {
#     '<phase>': [<callable(build)>, ...],
#     ...
# }
# ---------------------------------------------------------------------------

component_utils: dict[str, dict[str, list[Callable]]] = {}


def register(owner: str, phase: str, function: Callable) -> None:
    """Append `function` to the list of hooks for `(owner, phase)`.

    Mirrors the Perl idiom of mutating
    `%Galacticus::Build::Component::Utils::componentUtils` at module-load
    time.  Multiple registrations under the same owner+phase run in
    registration order, matching Perl's array iteration.
    """
    component_utils.setdefault(owner, {}).setdefault(phase, []).append(function)


# ---------------------------------------------------------------------------
# Global state.  Set by Label_Lengths during the `gather` phase.
# ---------------------------------------------------------------------------

boolean_label                            = ('false', 'true')

# Maximum lengths of various labels.  Initialised to 0 and then set by
# `Label_Lengths` during the `gather` phase; downstream code never reads
# them before that, matching Perl's `our $foo` being undef until populated
# (in Python we choose 0 over None so the pad_* helpers' arithmetic stays
# well-typed; reading any of these before Label_Lengths simply produces a
# zero-padded label rather than a TypeError on arithmetic).
class_name_length_max                    = 0
implementation_name_length_max           = 0
fully_qualified_name_length_max          = 0
property_name_length_max                 = 0
implementation_property_name_length_max  = 0
linked_data_name_length_max              = 0


# ---------------------------------------------------------------------------
# Type registries.  Plain dicts/sets — read constants only.
# ---------------------------------------------------------------------------

intrinsic_types = {
    'integer':     'integer',
    'longInteger': 'integer(kind=kind_int8)',
    'logical':     'logical',
    'double':      'double precision',
    'void':        'void',
}

intrinsic_nulls = {
    'double':      '0.0d0',
    'integer':     '0',
    'longInteger': '0_kind_int8',
    'logical':     '.false.',
}

output_type_map = {
    'double':      'double',
    'integer':     'integer',
    'longInteger': 'integer',
}


def is_intrinsic(type_name: str | None) -> bool:
    """Return True if `type_name` is one of the recognised intrinsic types.

    Mirrors `Galacticus::Build::Components::Utils::isIntrinsic`.
    """
    return type_name in intrinsic_types


def is_output_intrinsic(type_name: str | None) -> bool:
    """Return True if `type_name` is intrinsic *and* directly outputtable.

    Mirrors `Galacticus::Build::Components::Utils::isOutputIntrinsic`.
    """
    return type_name in ('double', 'integer', 'longInteger')


# ---------------------------------------------------------------------------
# Offset-name builder.  Two arities mirror the Perl module's two call
# patterns:
#   offset_name(status, component_name, property_name)
#   offset_name(status, class_dict,    member_dict, property_dict)
# ---------------------------------------------------------------------------

_OFFSET_SHORT_NAME = {
    'all':      'all',
    'active':   'atv',
    'inactive': 'itv',
}


def offset_name(*args: Any) -> str:
    """Return the variable name used to store the ODE solver offset of a
    property.  Mirrors `Components::Utils::offsetName`.
    """
    if len(args) == 3:
        status, component_name, property_name = args
        if status not in _OFFSET_SHORT_NAME:
            raise ValueError(f"unrecognized status '{status}'")
        return ('offset' + _ucfirst(_OFFSET_SHORT_NAME[status])
                + _ucfirst(component_name) + _ucfirst(property_name))
    if len(args) == 4:
        status, class_, member, prop = args
        if status not in _OFFSET_SHORT_NAME:
            raise ValueError(f"unrecognized status '{status}'")
        return ('offset' + _ucfirst(_OFFSET_SHORT_NAME[status])
                + _ucfirst(class_['name'])
                + _ucfirst(member['name'])
                + _ucfirst(prop  ['name']))
    raise TypeError(
        f"offset_name(): incorrect number of arguments ({len(args)})"
    )


# ---------------------------------------------------------------------------
# Padding.  Each `pad_*` helper rounds a label to the matching length-max.
# ---------------------------------------------------------------------------

def pad_class(text: str, extra_pad: tuple[int, int] = (0, 0)) -> str:
    return _pad_generic(class_name_length_max,                    text, extra_pad)


def pad_implementation(text: str, extra_pad: tuple[int, int] = (0, 0)) -> str:
    return _pad_generic(implementation_name_length_max,           text, extra_pad)


def pad_fully_qualified(text: str, extra_pad: tuple[int, int] = (0, 0)) -> str:
    return _pad_generic(fully_qualified_name_length_max,          text, extra_pad)


def pad_property_name(text: str, extra_pad: tuple[int, int] = (0, 0)) -> str:
    return _pad_generic(property_name_length_max,                 text, extra_pad)


def pad_implementation_property_name(text: str, extra_pad: tuple[int, int] = (0, 0)) -> str:
    return _pad_generic(implementation_property_name_length_max,  text, extra_pad)


def pad_linked_data(text: str, extra_pad: tuple[int, int] = (0, 0)) -> str:
    return _pad_generic(linked_data_name_length_max,              text, extra_pad)


def _pad_generic(length: int, text: str, extra_pad: tuple[int, int]) -> str:
    """Right-pad `text` to `length + extra_pad[0]`, but never shorter than
    `extra_pad[1]`.
    """
    pad_length = length + extra_pad[0]
    if extra_pad[1] > pad_length:
        pad_length = extra_pad[1]
    return text + ' ' * (pad_length - len(text))


# ---------------------------------------------------------------------------
# Apply-defaults: walk a nested structure, filling in default values where
# none are set.  Used by sub-modules to coerce missing XML attributes into
# a canonical form.  Booleans get translated from `"true"`/`"false"` strings
# into Python `True`/`False`.
# ---------------------------------------------------------------------------

def apply_defaults(obj: Any, name: str, default: Any) -> None:
    """Mirror Perl `applyDefaults($object, $name, $default)`.

    `default` may be either a plain scalar (or string starting with
    `"boolean"`) or a nested dict of further defaults.
    """
    if isinstance(default, dict):
        if isinstance(obj, dict) and name in obj:
            for sub_obj in _as_list(obj[name]):
                for k, v in default.items():
                    apply_defaults(sub_obj, k, v)
        elif name == 'ALL' and isinstance(obj, dict):
            for sub_obj in obj.values():
                for k, v in default.items():
                    apply_defaults(sub_obj, k, v)
        return

    for entry in _as_list(obj):
        if not isinstance(entry, dict):
            continue
        if isinstance(default, str) and default.startswith('boolean'):
            if name in entry:
                entry[name] = (entry[name] == 'true')
            else:
                entry[name] = (default == 'booleanTrue')
        else:
            entry.setdefault(name, default)


# ---------------------------------------------------------------------------
# argument_list: extract every argument-bearing variable from a list of
# variable descriptors.  Used when building Fortran function signatures.
# ---------------------------------------------------------------------------

def argument_list(*descriptors: Any) -> list[str]:
    """Return the flat list of argument names taken from `descriptors`.

    A descriptor contributes its `variables` to the result when it is:

    * a regular variable with an `intent(...)` attribute,
    * a procedure pointer (either via `pointer` attribute or `is_argument`
      flag), or
    * declared `external`.
    """
    arguments: list[str] = []
    for desc in descriptors:
        if not isinstance(desc, dict):
            continue
        attrs     = desc.get('attributes') or []
        intrinsic = desc.get('intrinsic')

        is_intent_var = any(
            (isinstance(a, str) and a.lstrip().startswith('intent('))
            for a in attrs
        )

        is_proc_ptr = (
            intrinsic == 'procedure'
            and (
                'pointer' in attrs
                or desc.get('isArgument') in (True, 1, 'true')
            )
        )

        is_external = (intrinsic == 'external')

        if is_intent_var or is_proc_ptr or is_external:
            arguments.extend(desc.get('variables') or [])

    return arguments


# ---------------------------------------------------------------------------
# Label_Lengths: compute the maximum length for each label kind, for
# column-aligned output.  Registered in the `gather` phase below.
# ---------------------------------------------------------------------------

def Label_Lengths(build: dict) -> None:
    """Set the module-level `*_length_max` globals from the components in
    `build`.  Pushes `propertyNameLengthMax` into `build['variables']` so
    the generated Fortran can reference it at run-time.
    """
    global class_name_length_max, implementation_name_length_max
    global fully_qualified_name_length_max
    global implementation_property_name_length_max
    global property_name_length_max

    components = list((build.get('components') or {}).values())

    class_name_length_max = max(
        (len(c.get('class', '')) for c in components), default=0
    )
    implementation_name_length_max = max(
        (len(c.get('name', '')) for c in components), default=0
    )
    fully_qualified_name_length_max = max(
        (len(c.get('name', '')) + len(c.get('class', '')) for c in components),
        default=0,
    )

    impl_prop_max = 0
    prop_max      = 0
    for component in components:
        for prop in _component_properties(component):
            prop_name = prop.get('name', '')
            impl_prop_max = max(
                impl_prop_max,
                len(component.get('name', '')) + len(component.get('class', ''))
                + len(prop_name),
            )
            prop_max = max(prop_max, len(prop_name))
    implementation_property_name_length_max = impl_prop_max
    property_name_length_max                = prop_max

    build.setdefault('variables', []).append({
        'intrinsic': 'integer',
        'variables': [f'propertyNameLengthMax={prop_max}'],
    })

    logger.info("         --> Maximum label lengths:")
    logger.info(f"            -->           Class: {class_name_length_max}")
    logger.info(f"            -->  Implementation: {implementation_name_length_max}")
    logger.info(f"            --> Fully-qualified: {fully_qualified_name_length_max}")
    logger.info(f"            -->        Property: {prop_max}")


# `Label_Lengths` is the only hook this module registers — the rest of the
# globals/helpers above are consumed by sister modules.
register('utils', 'gather', Label_Lengths)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _as_list(value: Any) -> list:
    """Normalise XMLin's "scalar or list" idiom into a list."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _component_properties(component: dict) -> list[dict]:
    """Yield every property declared on `component`, regardless of whether
    `properties.property` is a dict, a single dict, or a list of dicts.
    """
    props = (component.get('properties') or {}).get('property')
    if props is None:
        return []
    if isinstance(props, list):
        return props
    if isinstance(props, dict):
        # XMLin with KeyAttr left to defaults sometimes returns a dict of
        # name → property.  Force-array on `property` keeps it a list, but
        # be defensive in case a downstream caller hands us the keyed form.
        if all(isinstance(v, dict) for v in props.values()):
            return list(props.values())
        return [props]
    return []


def _ucfirst(text: str) -> str:
    return text[:1].upper() + text[1:] if text else text
