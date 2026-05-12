"""Scan Galacticus Fortran sources to build a derived-type hierarchy.

The bind(c) wrapper generator treats ``class(<base>Class)`` constructor
arguments as polymorphic functionClass pointers and routes them through
``<base>GetPtr``.  Some constructors instead take an *abstract intermediate*
such as ``class(massDistributionSpherical)`` — an abstract type that
itself extends a registered functionClass base (``massDistributionClass``).
The wrapper can still accept these via the base's ``GetPtr``, then narrow
the polymorphic pointer to the intermediate with a ``select type`` block.

To know which intermediates trace back to which registered base, we walk
every ``type, [...] extends(<parent>) [...] :: <name>`` declaration in
``source/*.F90`` and record ``name -> parent``.  Callers walk that map
upward until they hit a ``<base>Class`` entry whose ``<base>`` is in the
registered functionClass set.
"""

import re
from pathlib import Path

# `type, [optional attrs (abstract, public, private, extends(...))]
#  [::] <name>` — same shape Fortran.Utils.UNIT_OPENERS['type'] recognises,
# but with the parent name captured.  Only matches when the `extends(...)`
# attribute is present; types without it have no parent edge to record.
_TYPE_EXTENDS_RX = re.compile(
    r'^\s*type\s*'
    r'(?:,\s*(?:abstract|public|private|extends\s*\([A-Za-z][A-Za-z0-9_]*\))\s*)*'
    r',\s*extends\s*\(\s*([A-Za-z][A-Za-z0-9_]*)\s*\)\s*'
    r'(?:,\s*(?:abstract|public|private)\s*)*'
    r'(?:::)?\s*([A-Za-z][A-Za-z0-9_]*)\s*$',
    re.IGNORECASE | re.MULTILINE,
)
_MODULE_RX = re.compile(
    r'^\s*module\s+([A-Za-z][A-Za-z0-9_]*)\s*$',
    re.IGNORECASE | re.MULTILINE,
)


def build_type_hierarchy(source_dir):
    """Return ``{type_name: {'parent': parent_name, 'module': module_name}}``.

    *source_dir* is a path-like to the ``source/`` directory (or any tree
    containing ``*.F90`` files).  Types declared without an ``extends(...)``
    clause are absent from the returned map — callers walk the map by
    chasing ``parent`` upward and stop when the lookup misses.

    Order within a file matters only for the per-file module assignment:
    each ``type`` is recorded against the most recently-seen ``module`` line
    in that file.  Submodules are ignored (they re-open a parent module
    rather than introducing new types in practice).
    """
    hierarchy = {}
    for path in sorted(Path(source_dir).glob('*.F90')):
        text = path.read_text(errors='replace')
        # Build a list of (offset, module_name) so each later `type ...`
        # match can be attributed to the enclosing module via a bisect.
        module_anchors = [(m.start(), m.group(1))
                          for m in _MODULE_RX.finditer(text)]
        for tm in _TYPE_EXTENDS_RX.finditer(text):
            parent, child = tm.group(1), tm.group(2)
            mod = ''
            for off, name in module_anchors:
                if off < tm.start():
                    mod = name
                else:
                    break
            # Earlier files declaring the same type name (unusual but
            # possible via includes) lose to the first hit; that matches
            # what the SourceTree parse would surface for the canonical
            # declaration site.
            hierarchy.setdefault(child, {'parent': parent, 'module': mod})
    return hierarchy


def resolve_function_class_base(type_name, hierarchy, registered):
    """Walk *hierarchy* upward from *type_name* to find a registered base.

    Returns ``(base_name, intermediate)`` where:
      * ``base_name``   — the stem of a ``<base>Class`` ancestor whose
                          ``<base>`` is in *registered*, or ``None`` if no
                          such ancestor exists.
      * ``intermediate`` — the original *type_name* (echoed back so the
                          caller can record it as the narrowing target).

    The walk stops at the first match; if the supplied *type_name* itself
    ends in ``Class`` and its stem is registered, *base_name* is the stem
    and *intermediate* is ``None`` (no narrowing needed — this is the
    plain ``class(<base>Class)`` case the generator already handled).

    Cycle-safe: terminates if the parent chain ever revisits a name,
    which shouldn't happen in valid Fortran but would otherwise loop.
    """
    if type_name.endswith('Class') and type_name[:-5] in registered:
        return type_name[:-5], None

    visited = {type_name}
    current = type_name
    while True:
        parent = (hierarchy.get(current) or {}).get('parent')
        if not parent or parent in visited:
            return None, None
        if parent.endswith('Class') and parent[:-5] in registered:
            return parent[:-5], type_name
        visited.add(parent)
        current = parent
