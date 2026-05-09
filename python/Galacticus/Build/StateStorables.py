"""Helpers for reading the structure of `$BUILDPATH/stateStorables.xml`.

Andrew Benson (ported to Python 2026)
"""
from __future__ import annotations
#
# `dict_to_xml_string()` writes `stateStorables.xml` in the no-attribute
# (`NoAttr => 1`) style: every value is a child element rather than an
# attribute.  A list value (e.g. the list of functionClasses) becomes a
# series of sibling elements with the same tag — for example
#
#     <storables>
#       <functionClasses>
#         <name>massDistributionClass</name>
#         <module>Mass_Distributions</module>
#       </functionClasses>
#       <functionClasses>...</functionClasses>
#       ...
#     </storables>
#
# When read back through `xml_to_dict()`, that becomes
# `{'functionClasses': [{'name': ..., 'module': ...}, ...]}` — a *list* of
# entry dicts.  The original Perl scripts produced (and consumed) the
# attribute-keyed shape `{'functionClasses': {'X': {'module': ...}, ...}}`,
# so most callers in the codebase still test `isinstance(fc, dict)` first
# and then look for a `functionClass` wrapper key.  Centralise the shape-
# bridging here so every reader handles the actual on-disk format.

__all__ = [
    'function_class_entries', 'function_class_names',
    'function_class_module_map', 'function_class_entry',
    'function_class_instances', 'event_hook_static_names',
]


def function_class_entries(state_storables: dict | None) -> list[dict]:
    """Return a list of ``{'name': ..., 'module': ..., ...}`` dicts for every
    functionClass listed in ``stateStorables.xml``.

    Accepts every XML shape that `xml_to_dict()` can produce and that the
    legacy Perl writers ever produced:

    * list of dicts                            (current writer, multiple entries)
    * single dict with a ``name`` key          (current writer, single entry)
    * ``{'functionClass': [...]}`` wrapper     (XML::Simple ForceArray style)
    * ``{'X': {'module': ...}, ...}`` mapping  (XML::Simple KeyAttr style)
    * ``None``/missing                          (returns ``[]``)
    """
    fc = (state_storables or {}).get('functionClasses')
    if fc is None:
        return []

    if isinstance(fc, list):
        return [e for e in fc if isinstance(e, dict)]

    if isinstance(fc, dict):
        # XML::Simple KeyAttr-style: ``{'X': {'module': ...}, 'Y': {...}}``.
        entries = fc.get('functionClass')
        if entries is not None:
            if isinstance(entries, dict):
                return [entries]
            return [e for e in entries if isinstance(e, dict)]
        # Single ``<functionClasses>`` element with ``<name>`` child.
        if 'name' in fc:
            return [fc]
        # KeyAttr-keyed mapping: synthesise entry dicts from {key: rest}.
        out = []
        for key, value in fc.items():
            if isinstance(value, dict):
                merged = dict(value)
                merged.setdefault('name', key)
                out.append(merged)
        return out

    return []


def function_class_names(state_storables: dict | None) -> set[str]:
    """Set of `<name>Class` strings from ``stateStorables.xml``."""
    return {e['name'] for e in function_class_entries(state_storables)
            if e.get('name')}


def function_class_module_map(state_storables: dict | None) -> dict[str, str | None]:
    """Map ``{<name>Class: providing_module}`` from ``stateStorables.xml``."""
    return {e['name']: e.get('module')
            for e in function_class_entries(state_storables)
            if e.get('name')}


def function_class_entry(state_storables: dict | None, class_name: str) -> dict | None:
    """Return the entry dict for ``class_name`` (e.g. ``'massDistributionClass'``)
    or ``None`` if no such entry exists.
    """
    for e in function_class_entries(state_storables):
        if e.get('name') == class_name:
            return e
    return None


def function_class_instances(state_storables: dict | None) -> list[str]:
    """List of every concrete ``functionClass`` instance name."""
    raw = (state_storables or {}).get('functionClassInstances')
    if raw is None:
        return []
    if isinstance(raw, str):
        return [raw] if raw else []
    if isinstance(raw, dict):
        # Single instance written as a dict (with optional 'name' key, or just
        # text content).
        name = raw.get('name') or raw.get('content')
        return [name] if name else []
    out = []
    for item in raw:
        if isinstance(item, str) and item:
            out.append(item)
        elif isinstance(item, dict):
            name = item.get('name') or item.get('content')
            if name:
                out.append(name)
    return out


def event_hook_static_names(state_storables: dict | None) -> list[str]:
    """List of every ``eventHookStatic`` name."""
    raw = (state_storables or {}).get('eventHookStatics')
    if raw is None:
        return []
    if isinstance(raw, str):
        return [raw] if raw else []
    if isinstance(raw, dict):
        name = raw.get('name') or raw.get('content')
        return [name] if name else []
    out = []
    for item in raw:
        if isinstance(item, str) and item:
            out.append(item)
        elif isinstance(item, dict):
            name = item.get('name') or item.get('content')
            if name:
                out.append(name)
    return out
