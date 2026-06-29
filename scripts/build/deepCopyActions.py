#!/usr/bin/env python3
"""Catalog every derived type that takes part in a `deepCopyActions` directive.

Walks every source file listed under `directiveLocations.deepCopyActions.file`,
parses it with `Galacticus.Build.SourceTree.parse_file`, collects every
`deepCopyActions` directive node and every `type ... :: name` definition, and
for each directive emits one `(type, class)` record per derived class in the
inheritance chain rooted at the directive's `class`.

Writes `$BUILDPATH/deepCopyActions.xml` and caches per-file results in
`$BUILDPATH/deepCopyActions.blob`.

Mirrors scripts/build/deepCopyActions.pl.
Andrew Benson (ported to Python 2026).
"""

import os
import pickle
import re
import sys
import xml.etree.ElementTree as ET


from Galacticus.Build.ParallelScan import scan as parallel_scan
from Galacticus.Build.SourceTree   import parse_file, walk_tree
from List.ExtraUtils              import as_array
from XML.Utils                    import dict_to_xml_string, xml_to_dict


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_TYPE_OPENER_RE = re.compile(
    r'^\s*type\s*'
    r'(?:,\s*(?:(abstract)|public|private|extends\s*\(([a-zA-Z0-9_]+)\))\s*)*'
    r'(?:::)?\s*([a-zA-Z0-9_]+)\s*$',
    re.IGNORECASE,
)


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
    """Read `directiveLocations.xml` and return the xml_to_dict structure."""
    tree = ET.parse(path)
    return xml_to_dict(tree.getroot())


def _collect_types_and_directives(tree, directive_type):
    """Walk `tree`, returning (type-dict, directive-node-list).

    type-dict: {type_name: {'extends': parent-or-None, 'abstract': bool}}
    directive-node-list: every AST node whose `type` attribute equals
                         `directive_type`.
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
                    "deepCopyActions.py: unable to parse type definition "
                    f"opener: {node.get('opener')!r}"
                )
            abstract, extends, name = m.group(1), m.group(2), m.group(3)
            classes[name] = {
                'extends':  extends,
                'abstract': abstract is not None,
            }
    return classes, directives


def _inherits_from(classes, class_name, base_class):
    """True if `class_name` derives (directly or transitively) from `base_class`."""
    current = class_name
    while current is not None:
        if current == base_class:
            return True
        current = classes.get(current, {}).get('extends')
    return False


def _scan_one(task):
    """Worker: parse one file and return `(file_identifier, entries)`. Mirrors
    the per-file body of the serial loop; writes no shared state. The file read
    + AST parse (the slow, NFS-blocking part) is what runs concurrently.
    """
    file_identifier, file_name = task
    tree = parse_file(file_name)
    classes, directives = _collect_types_and_directives(tree, 'deepCopyActions')
    entries = []
    for node in directives:
        directive  = node.get('directive') or {}
        base_class = directive.get('class')
        for class_name in sorted(classes):
            if _inherits_from(classes, class_name, base_class):
                entries.append({'type': class_name, 'class': base_class})
    return file_identifier, entries


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv):
    if len(argv) != 2:
        print("Usage: deepCopyActions.py <installDirectory>", file=sys.stderr)
        sys.exit(1)

    build_path = os.environ['BUILDPATH']
    blob_path  = os.path.join(build_path, 'deepCopyActions.blob')

    directive_locations = _parse_directive_locations(
        os.path.join(build_path, 'directiveLocations.xml'),
    )

    actions_per_file, cache_mtime = _load_cache(blob_path)
    have_per_file = cache_mtime is not None

    # Select the files that need a (re)scan, preserving directive-location order.
    files = as_array((directive_locations.get('deepCopyActions') or {}).get('file'))
    scan_list = []
    for file_name in files:
        file_identifier = _file_identifier(file_name)
        if (have_per_file
                and file_identifier in actions_per_file
                and os.stat(file_name).st_mtime < cache_mtime):
            continue
        scan_list.append((file_identifier, file_name))

    # Parse the files concurrently; merge in scan order, reproducing the serial
    # pop-then-maybe-set behaviour exactly.
    for file_identifier, entries in parallel_scan(
            scan_list, _scan_one, 'deepCopyActions.py'):
        actions_per_file.pop(file_identifier, None)
        if entries:
            actions_per_file[file_identifier] = {'deepCopyActions': entries}

    # Aggregate + sort by `type`.
    all_entries = []
    for entry in actions_per_file.values():
        all_entries.extend(entry.get('deepCopyActions') or [])
    all_entries.sort(key=lambda e: e['type'])

    out_path = os.path.join(build_path, 'deepCopyActions.xml')
    with open(out_path, 'w') as fh:
        fh.write(
            dict_to_xml_string('deepCopyActions',
                               {'deepCopyActions': all_entries})
        )

    with open(blob_path, 'wb') as fh:
        pickle.dump(actions_per_file, fh, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    main(sys.argv)
