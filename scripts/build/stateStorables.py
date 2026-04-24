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

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from build.file_changes               import update as file_changes_update
from build.fortran_utils              import get_fortran_line
from Fortran.Utils                    import UNIT_OPENERS
from Galacticus.Build.Directives      import extract_directives
from Galacticus.Build.SourceTree      import parse_file, walk_tree
from List.ExtraUtils                  import as_array
from XML.Utils                        import dict_to_xml_string, xml_to_dict


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_TYPE_OPENER_RE = re.compile(
    r'^\s*type\s*'
    r'(?:,\s*(?:(abstract)|public|private|extends\s*\(([a-zA-Z0-9_]+)\))\s*)*'
    r'(?:::)?\s*([a-zA-Z0-9_]+)\s*$',
    re.IGNORECASE,
)

_MODULE_OPENER_RE = UNIT_OPENERS['module']['regex']
_MODULE_NAME_GROUP = UNIT_OPENERS['module']['unit_name']  # 0-based index

_MODULE_INLINE_RE = re.compile(r'^module\s+([a-zA-Z0-9_]+)')


def _file_identifier(path):
    """Perl `(my $id = $path) =~ s/\\//_/g; $id =~ s/^\\._??//;`."""
    return re.sub(r'^\._?', '', path.replace('/', '_'))


def _load_cache(blob_path):
    if not os.path.exists(blob_path):
        return {}, None
    try:
        with open(blob_path, 'rb') as fh:
            cache = pickle.load(fh)
    except (pickle.UnpicklingError, EOFError, AttributeError, ValueError,
            ImportError, ModuleNotFoundError):
        return {}, None
    if not isinstance(cache, dict):
        return {}, None
    return cache, os.stat(blob_path).st_mtime


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


def _collect_types_and_directives(tree, directive_type):
    """Walk `tree`, returning (type-dict, directive-node-list).  See the same
    helper in deepCopyActions.py for the structure.
    """
    classes    = {}
    directives = []
    for node in walk_tree(tree):
        if node.get('type') == directive_type:
            directives.append(node)
        if node.get('type') == 'type':
            m = _TYPE_OPENER_RE.match(node.get('opener') or '')
            if not m:
                raise RuntimeError(
                    "stateStorables.py: unable to parse type definition "
                    f"opener: {node.get('opener')!r}"
                )
            abstract, extends, name = m.group(1), m.group(2), m.group(3)
            classes[name] = {
                'extends':  extends,
                'abstract': abstract is not None,
            }
    return classes, directives


def _inherits_from(classes, class_name, base_class):
    current = class_name
    while current is not None:
        if current == base_class:
            return True
        current = classes.get(current, {}).get('extends')
    return False


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
    # -----------------------------------------------------------------------
    function_class_files = as_array(
        (directive_locations.get('functionClass') or {}).get('file')
    )
    for file_name in function_class_files:
        file_identifier = _file_identifier(file_name)
        if not _fresh(cache_mtime,
                      file_identifier in storables_per_file, file_name):
            storables_per_file.pop(file_identifier, None)
            modules = _collect_module_openers(file_name)
            entry   = storables_per_file.setdefault(file_identifier, {})
            for fc in extract_directives(file_name, 'functionClass'):
                descriptor = {'name': fc['name'] + 'Class'}
                if len(modules) == 1:
                    descriptor['module'] = modules[0]
                entry.setdefault('functionClasses', []).append(descriptor)
                entry.setdefault('functionClassDirectives', []).append(fc['name'])

        # Pass 2: for each functionClass declared here, enumerate its
        # concrete instances -- one per file under the directive's name.
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
                storables_per_file.pop(instance_identifier, None)
                instance_entry = storables_per_file.setdefault(
                    instance_identifier, {},
                )
                for inst in extract_directives(instance_file,
                                               function_class_name):
                    instance_entry.setdefault(
                        'functionClassInstances', [],
                    ).append(inst['name'])

    # -----------------------------------------------------------------------
    # Pass 3: every stateStorable file -- walk its AST and record derived
    # types that match the directive's declared base class.
    # -----------------------------------------------------------------------
    state_storable_files = as_array(
        (directive_locations.get('stateStorable') or {}).get('file')
    )
    for file_name in state_storable_files:
        file_identifier = _file_identifier(file_name)
        if _fresh(cache_mtime,
                  file_identifier in storables_per_file, file_name):
            continue
        storables_per_file.pop(file_identifier, None)

        tree = parse_file(file_name)
        classes, directives = _collect_types_and_directives(
            tree, 'stateStorable',
        )
        entries = []
        for node in directives:
            directive  = node.get('directive') or {}
            base_class = directive.get('class')
            module     = _enclosing_module_name(node)
            for class_name in sorted(classes):
                if _inherits_from(classes, class_name, base_class):
                    entries.append({
                        'type':   class_name,
                        'class':  base_class,
                        'module': module,
                    })
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
    for file_name in event_static_files:
        file_identifier = _file_identifier(file_name)
        if _fresh(cache_mtime,
                  file_identifier in storables_per_file, file_name):
            continue
        storables_per_file.pop(file_identifier, None)
        names = [
            d['name']
            for d in extract_directives(file_name, 'eventHookStatic')
        ]
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
