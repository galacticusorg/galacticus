"""Helpers for reading the structure of `$BUILDPATH/stateStorables.xml`.

Andrew Benson (ported to Python 2026)
"""
from __future__ import annotations

import sys
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
# entry dicts.  Legacy writers produced (and consumed) the attribute-keyed
# shape `{'functionClasses': {'X': {'module': ...}, ...}}`,
# so most callers in the codebase still test `isinstance(fc, dict)` first
# and then look for a `functionClass` wrapper key.  Centralize the shape-
# bridging here so every reader handles the actual on-disk format.

__all__ = [
    'function_class_entries', 'function_class_names',
    'function_class_module_map', 'function_class_entry',
    'function_class_instances', 'event_hook_static_names',
    'warn_unstored_member', 'is_self_storing_type', 'is_unstored_type',
]

# Types that are not registered as state-storable, but which provide their own `stateStore(stateFile)` and
# `stateRestore(stateFile)` type-bound procedures (with no GSL state file and no store identifier). The state store
# generators store members of these types by calling those procedures. Names are lower case.
SELF_STORING_TYPES = frozenset({
    'tabulationinverse',
})

# Types whose members the state store generators skip without warning, each with the reason. These are either handles to
# resources that exist only in the running process (a class holding one must make sure that, after a restore, the
# handle is re-created when needed, e.g. with a `<restoreTo>` that resets an "initialized" flag), or types whose objects
# are fully configured by their constructor and never change afterwards (so the freshly constructed object in the
# restoring run already holds the correct value). Names are lower case.
UNSTORED_TYPES = {
    'omplock'           : 'an OpenMP lock',
    'ompreadwritelock'  : 'an OpenMP lock',
    'lockdescriptor'    : 'a file lock',
    'hdf5file'          : 'a handle to an open HDF5 file',
    'hdf5group'         : 'a handle to an open HDF5 group',
    'file'              : 'a handle to an open file',
    'c_ptr'             : 'a C pointer',
    'c_funptr'          : 'a C function pointer',
    'resourcemanager'   : 'a reference-counted handle to a shared resource',
    'timer'             : 'a wall-clock timer',
    'inputparameters'   : 'a handle to a parameter set, fixed at construction',
    'rootfinder'        : 'configured at construction, holding only settings and GSL handles',
    'fastexponentiator' : 'a lookup table built by its constructor and never changed',
    # Strings are names, labels and file names, set by the constructor from parameters. Restoring them would be wrong,
    # since some may legitimately differ between the storing and restoring runs (e.g. output file names). A class which
    # changes a string after construction must arrange for it to be consistent after a restore (e.g. by excluding any
    # flag recording that the string has been changed).
    'varying_string'    : 'a name, label, or file name, set by the constructor from parameters',
}


def is_self_storing_type(type_name: str | None) -> bool:
    """True if `type_name` provides its own `stateStore(stateFile)`/`stateRestore(stateFile)` procedures."""
    return (type_name or '').replace(' ', '').lower() in SELF_STORING_TYPES


def is_unstored_type(type_name: str | None) -> bool:
    """True if members of `type_name` are deliberately not stored (see `UNSTORED_TYPES`)."""
    return (type_name or '').replace(' ', '').lower() in UNSTORED_TYPES


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
        # KeyAttr-keyed mapping: synthesize entry dicts from {key: rest}.
        out = []
        for key, value in fc.items():
            if isinstance(value, dict):
                merged = dict(value)
                merged.setdefault('name', key)
                out.append(merged)
        return out

    return []


# Lookup tables derived from the most recently queried ``stateStorables`` dict.  The preprocessor hooks query the same
# (never-mutated) dict once per tree node — tens of thousands of times per file — so rebuilding these on every call
# dominated their cost.  The dict itself is held so that the identity check cannot be fooled by `id()` reuse.
_lookup_cache = (None, None, None)


def _lookups(state_storables: dict | None) -> tuple[frozenset[str], dict[str, dict]]:
    """Return ``(names, entries_by_name)`` for ``state_storables``, cached by identity."""
    global _lookup_cache
    cached, names, by_name = _lookup_cache
    if cached is None or cached is not state_storables:
        by_name = {}
        for e in function_class_entries(state_storables):
            # The first entry of a given name wins, matching a linear scan.
            if e.get('name'):
                by_name.setdefault(e['name'], e)
        names         = frozenset(by_name)
        _lookup_cache = (state_storables, names, by_name)
    return names, by_name


def function_class_names(state_storables: dict | None) -> frozenset[str]:
    """Set of `<name>Class` strings from ``stateStorables.xml``."""
    return _lookups(state_storables)[0]


def function_class_module_map(state_storables: dict | None) -> dict[str, str | None]:
    """Map ``{<name>Class: providing_module}`` from ``stateStorables.xml``."""
    return {e['name']: e.get('module')
            for e in function_class_entries(state_storables)
            if e.get('name')}


def function_class_entry(state_storables: dict | None, class_name: str) -> dict | None:
    """Return the entry dict for ``class_name`` (e.g. ``'massDistributionClass'``)
    or ``None`` if no such entry exists.
    """
    return _lookups(state_storables)[1].get(class_name)


# The most recently queried ``stateStorables`` dict and its instance names (see `_lookup_cache` above).
_instances_cache = (None, None)


def function_class_instances(state_storables: dict | None) -> list[str]:
    """List of every concrete ``functionClass`` instance name."""
    global _instances_cache
    cached, instances = _instances_cache
    if cached is None or cached is not state_storables:
        instances        = _function_class_instances(state_storables)
        _instances_cache = (state_storables, instances)
    # Return a copy, so that a caller modifying its list cannot corrupt the cache.
    return list(instances)


def _function_class_instances(state_storables: dict | None) -> list[str]:
    """Uncached implementation of `function_class_instances()`."""
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


def warn_unstored_member(owner: str, member: str, type_name: str,
                         source: str | None = None) -> None:
    """Report (on standard error) a derived-type member that state store/restore code will not store.

    The generators can store only members whose type is registered as state-storable (a functionClass, or a type with a
    `stateStorable` directive). A member of any other type is skipped, so after a state restore it keeps whatever value
    the object had before the restore, which silently breaks reproducibility (see issue #1786). Each such member must
    therefore either have its type registered, or be listed in the owning class' `<stateStorable><exclude
    variables="..."/></stateStorable>` to record that it need not be stored.
    """
    where = f" ({source})" if source else ''
    sys.stderr.write(
        f"Warning: state store for '{owner}'{where} skips member '{member}' of type '{type_name}', which is not "
        "state-storable - register the type as state-storable, or list the member in "
        f"<stateStorable><exclude variables=\"{member}\"/></stateStorable> if it need not be stored\n"
    )
