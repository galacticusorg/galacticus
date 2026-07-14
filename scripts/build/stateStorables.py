#!/usr/bin/env python3
"""Catalog classes (plus `functionClass` instances and `eventHookStatic`
objects) that support state storage/restoration.

For each file listed in `directiveLocations.xml` under a relevant directive
name, enumerates:

* `functionClasses`        - `{name: '<functionClass>Class', module: ...}`
                             for every `functionClass` directive in the tree.
* `functionClassInstances` - every concrete instance name declared via the
                             per-functionClass directive (one per class).
* `stateStorables`         - every derived type `{type, class, module}`
                             produced by walking each `stateStorable`
                             directive's inheritance tree.
* `eventHookStatics`       - names of every `eventHookStatic` directive.

Writes `$BUILDPATH/stateStorables.xml` (atomic update) and caches per-file
results in `$BUILDPATH/stateStorables.blob` (also atomic).

Mirrors scripts/build/stateStorables.pl.
Andrew Benson (ported to Python 2026).
"""

import os
import pickle
import re
import sys
import xml.etree.ElementTree as ET


from Galacticus.Build.FileChanges               import update as file_changes_update
from Galacticus.Build.FortranUtils              import get_fortran_line
from Fortran.Utils                    import UNIT_OPENERS
from Galacticus.Build.Directives      import extract_directives
from Galacticus.Build.ParallelScan    import scan as parallel_scan
from Galacticus.Build.SourceTree      import parse_file
from Galacticus.Build.TypeScan        import collect_types_and_directives, inherits_from
from List.ExtraUtils                  import as_array
from XML.Utils                        import dict_to_xml_string, xml_to_dict
from Galacticus.Build.ScanCache import (
    file_identifier as _file_identifier,
    load_cache      as _load_cache,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_MODULE_OPENER_RE = UNIT_OPENERS['module']['regex']
_MODULE_NAME_GROUP = UNIT_OPENERS['module']['unit_name']  # 0-based index

_MODULE_INLINE_RE = re.compile(r'^module\s+([a-zA-Z0-9_]+)')


def _parse_directive_locations(path):
    tree = ET.parse(path)
    return xml_to_dict(tree.getroot())


def _fresh(cache_mtime, entry_exists, file_name):
    """Perl's `$havePerFile && exists(...) && -M $file > $updateTime` idiom.

    Returns True when the cached entry is still current and we should skip
    rescanning this file.
    """
    if cache_mtime is None or not entry_exists:
        return False
    try:
        return os.stat(file_name).st_mtime < cache_mtime
    except FileNotFoundError:
        return False


def _collect_module_openers(file_name):
    """Return every module name opened in `file_name` (one per `module X`
    unit-opener line).

    Uses `get_fortran_line` + the module opener regex -- the same
    `$Fortran::Utils::unitOpeners{'module'}->{'regEx'}` the Perl uses.
    """
    modules = []
    with open(file_name, 'r', errors='replace') as fh:
        while True:
            raw, processed, _ = get_fortran_line(fh)
            if not raw and not processed:
                break
            m = _MODULE_OPENER_RE.match(processed)
            if m:
                modules.append(m.group(_MODULE_NAME_GROUP + 1))
    return modules


def _enclosing_module_name(node):
    """Walk up `node`'s parents, returning the first enclosing `module`
    name (or `""` if none).  Mirrors the loop at stateStorables.pl:119-126.
    """
    cur = node
    while cur is not None:
        if cur.get('type') == 'module':
            m = _MODULE_INLINE_RE.match(cur.get('opener') or '')
            if m:
                return m.group(1)
        cur = cur.get('parent')
    return ''


# ---------------------------------------------------------------------------
# Per-file workers (run concurrently; each reads files only, no shared state)
# ---------------------------------------------------------------------------
#
# Each pass's slow part is the per-file open/read/parse over NFS. The workers
# below do exactly that and return plain data; main() merges the results in the
# original order, reproducing the serial pop/setdefault sequence (and hence the
# cache blob) exactly. The final XML is built from fully sorted lists, so it is
# insensitive to ordering regardless.

def _scan_function_class(task):
    """Pass 1: module openers + functionClass directive names for one file."""
    file_identifier, file_name = task
    modules  = _collect_module_openers(file_name)
    fc_names = [fc['name'] for fc in extract_directives(file_name, 'functionClass')]
    return file_identifier, modules, fc_names


def _scan_instance(task):
    """Pass 2: concrete instance names for one (instance-file, functionClass)."""
    instance_identifier, instance_file, function_class_name = task
    names = [inst['name']
             for inst in extract_directives(instance_file, function_class_name)]
    return instance_identifier, names


def _scan_state_storable(task):
    """Pass 3: derived types matching each stateStorable directive's base."""
    file_identifier, file_name = task
    tree = parse_file(file_name)
    classes, directives = collect_types_and_directives(tree, 'stateStorable')
    entries = []
    for node in directives:
        directive  = node.get('directive') or {}
        base_class = directive.get('class')
        module     = _enclosing_module_name(node)
        for class_name in sorted(classes):
            if inherits_from(classes, class_name, base_class):
                entries.append({
                    'type':   class_name,
                    'class':  base_class,
                    'module': module,
                })
    return file_identifier, entries


def _scan_event_static(task):
    """Pass 4: eventHookStatic directive names for one file."""
    file_identifier, file_name = task
    names = [d['name'] for d in extract_directives(file_name, 'eventHookStatic')]
    return file_identifier, names


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv):
    if len(argv) != 2:
        print("Usage: stateStorables.py <installDirectory>", file=sys.stderr)
        sys.exit(1)

    build_path = os.environ['BUILDPATH']
    blob_path  = os.path.join(build_path, 'stateStorables.blob')

    directive_locations = _parse_directive_locations(
        os.path.join(build_path, 'directiveLocations.xml'),
    )

    storables_per_file, cache_mtime = _load_cache(blob_path)
    have_per_file = cache_mtime is not None

    # -----------------------------------------------------------------------
    # Pass 1: every functionClass file -- extract modules + directives.
    # Parsed concurrently, merged in file order (reproducing the serial
    # pop/setdefault/append sequence exactly).
    # -----------------------------------------------------------------------
    function_class_files = as_array(
        (directive_locations.get('functionClass') or {}).get('file')
    )
    p1_tasks = [
        (_file_identifier(f), f) for f in function_class_files
        if not _fresh(cache_mtime,
                      _file_identifier(f) in storables_per_file, f)
    ]
    p1_results = {
        fid: (modules, fc_names)
        for fid, modules, fc_names in parallel_scan(
            p1_tasks, _scan_function_class, 'stateStorables.py')
    }
    for file_name in function_class_files:
        file_identifier = _file_identifier(file_name)
        if file_identifier not in p1_results:
            continue                      # fresh -- keep the cached entry
        modules, fc_names = p1_results[file_identifier]
        storables_per_file.pop(file_identifier, None)
        entry = storables_per_file.setdefault(file_identifier, {})
        for name in fc_names:
            descriptor = {'name': name + 'Class'}
            if len(modules) == 1:
                descriptor['module'] = modules[0]
            entry.setdefault('functionClasses', []).append(descriptor)
            entry.setdefault('functionClassDirectives', []).append(name)

    # -----------------------------------------------------------------------
    # Pass 2: for each functionClass, enumerate its concrete instances. Tasks
    # are built in the original nested order (so the merge's pop/append, which
    # is order-sensitive when an instance file recurs, matches the serial run);
    # Pass 1 must be merged first since this reads functionClassDirectives.
    # -----------------------------------------------------------------------
    p2_tasks = []
    for file_name in function_class_files:
        file_identifier = _file_identifier(file_name)
        directive_names = (
            storables_per_file.get(file_identifier, {})
                              .get('functionClassDirectives') or []
        )
        for function_class_name in directive_names:
            instance_files = as_array(
                (directive_locations.get(function_class_name) or {}).get('file')
            )
            for instance_file in instance_files:
                instance_identifier = _file_identifier(instance_file)
                if _fresh(cache_mtime,
                          instance_identifier in storables_per_file,
                          instance_file):
                    continue
                p2_tasks.append(
                    (instance_identifier, instance_file, function_class_name)
                )
    # The cache entry for a freshly-scanned instance file is popped ONCE (to
    # discard any stale instances left from a previous run), then every task's
    # names are appended. A single instance file that implements instances of
    # more than one functionClass yields one task per class; popping only on
    # the first of them lets those tasks accumulate, instead of a later task's
    # pop discarding an earlier one's names (which silently dropped that file's
    # instances of every class but the last).
    popped = set()
    for instance_identifier, names in parallel_scan(
            p2_tasks, _scan_instance, 'stateStorables.py'):
        if instance_identifier not in popped:
            storables_per_file.pop(instance_identifier, None)
            popped.add(instance_identifier)
        instance_entry = storables_per_file.setdefault(instance_identifier, {})
        for name in names:
            instance_entry.setdefault('functionClassInstances', []).append(name)

    # -----------------------------------------------------------------------
    # Pass 3: every stateStorable file -- walk its AST and record derived
    # types that match the directive's declared base class.
    # -----------------------------------------------------------------------
    state_storable_files = as_array(
        (directive_locations.get('stateStorable') or {}).get('file')
    )
    p3_tasks = [
        (_file_identifier(f), f) for f in state_storable_files
        if not _fresh(cache_mtime,
                      _file_identifier(f) in storables_per_file, f)
    ]
    for file_identifier, entries in parallel_scan(
            p3_tasks, _scan_state_storable, 'stateStorables.py'):
        storables_per_file.pop(file_identifier, None)
        if entries:
            storables_per_file.setdefault(
                file_identifier, {},
            )['stateStorables'] = entries

    # -----------------------------------------------------------------------
    # Pass 4: eventHookStatic files -- record directive names.
    # -----------------------------------------------------------------------
    event_static_files = as_array(
        (directive_locations.get('eventHookStatic') or {}).get('file')
    )
    p4_tasks = [
        (_file_identifier(f), f) for f in event_static_files
        if not _fresh(cache_mtime,
                      _file_identifier(f) in storables_per_file, f)
    ]
    for file_identifier, names in parallel_scan(
            p4_tasks, _scan_event_static, 'stateStorables.py'):
        storables_per_file.pop(file_identifier, None)
        if names:
            storables_per_file.setdefault(
                file_identifier, {},
            )['eventHookStatics'] = names

    # -----------------------------------------------------------------------
    # Aggregate across files and sort.
    # -----------------------------------------------------------------------
    event_hook_statics       = []
    function_classes         = []
    function_class_instances = []
    state_storables          = []
    for entry in storables_per_file.values():
        event_hook_statics      .extend(entry.get('eventHookStatics')       or [])
        function_classes        .extend(entry.get('functionClasses')        or [])
        function_class_instances.extend(entry.get('functionClassInstances') or [])
        state_storables         .extend(entry.get('stateStorables')         or [])

    event_hook_statics.sort()
    function_classes.sort(key=lambda d: d['name'])
    function_class_instances.sort()
    state_storables.sort(key=lambda d: d['type'])

    # -----------------------------------------------------------------------
    # Write outputs (atomic).
    # -----------------------------------------------------------------------
    xml_tmp = os.path.join(build_path, 'stateStorables.xml.tmp')
    with open(xml_tmp, 'w') as fh:
        fh.write(dict_to_xml_string('storables', {
            'eventHookStatics':       event_hook_statics,
            'functionClasses':        function_classes,
            'functionClassInstances': function_class_instances,
            'stateStorables':         state_storables,
        }))
    file_changes_update(
        os.path.join(build_path, 'stateStorables.xml'),
        xml_tmp,
    )

    blob_tmp = blob_path + '.tmp'
    with open(blob_tmp, 'wb') as fh:
        pickle.dump(storables_per_file, fh, protocol=pickle.HIGHEST_PROTOCOL)
    file_changes_update(blob_path, blob_tmp)


if __name__ == '__main__':
    main(sys.argv)
