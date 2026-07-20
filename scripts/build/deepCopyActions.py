#!/usr/bin/env python3
"""Catalog every derived type that takes part in a `deepCopyActions` directive.

Walks every source file listed under `directiveLocations.deepCopyActions.file`,
parses it with `Galacticus.Build.SourceTree.parse_file`, collects every
`deepCopyActions` directive node and every `type ... :: name` definition, and
for each directive emits one `(type, class)` record per derived class in the
inheritance chain rooted at the directive's `class`.

Writes `$BUILDPATH/deepCopyActions.xml` and caches per-file results in
`$BUILDPATH/deepCopyActions.blob`.

Andrew Benson (2026).
"""

import os
import pickle
import re
import sys
import xml.etree.ElementTree as ET


from Galacticus.Build.FileChanges  import update as file_changes_update
from Galacticus.Build.ParallelScan import scan as parallel_scan
from Galacticus.Build.SourceTree   import parse_file
from Galacticus.Build.TypeScan     import collect_types_and_directives, inherits_from
from List.ExtraUtils              import as_array
from XML.Utils                    import dict_to_xml_string, xml_to_dict
from Galacticus.Build.ScanCache import (
    file_identifier as _file_identifier,
    load_cache      as _load_cache,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_directive_locations(path):
    """Read `directiveLocations.xml` and return the xml_to_dict structure."""
    tree = ET.parse(path)
    return xml_to_dict(tree.getroot())


def _scan_one(task):
    """Worker: parse one file and return `(file_identifier, entries)`. Mirrors
    the per-file body of the serial loop; writes no shared state. The file read
    + AST parse (the slow, NFS-blocking part) is what runs concurrently.
    """
    file_identifier, file_name = task
    tree = parse_file(file_name)
    classes, directives = collect_types_and_directives(tree, 'deepCopyActions')
    entries = []
    for node in directives:
        directive  = node.get('directive') or {}
        base_class = directive.get('class')
        for class_name in sorted(classes):
            if inherits_from(classes, class_name, base_class):
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
        if not os.path.exists(file_name):
            # Stale directiveLocations.xml entry (the file was deleted since
            # the catalog was generated): drop any cached entry and skip it
            # rather than crashing in os.stat()/parse_file().
            actions_per_file.pop(file_identifier, None)
            continue
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

    # Write both outputs atomically and only-if-changed (mtime preserved when
    # content is identical, matching stateStorables.py). The xml's mtime
    # stability matters: it is a prerequisite of every preprocessing rule, so
    # an unconditional rewrite would re-preprocess the entire source tree on
    # every catalog regeneration.
    out_path = os.path.join(build_path, 'deepCopyActions.xml')
    xml_tmp  = out_path + '.tmp'
    with open(xml_tmp, 'w') as fh:
        fh.write(
            dict_to_xml_string('deepCopyActions',
                               {'deepCopyActions': all_entries})
        )
    file_changes_update(out_path, xml_tmp)

    blob_tmp = blob_path + '.tmp'
    with open(blob_tmp, 'wb') as fh:
        pickle.dump(actions_per_file, fh, protocol=pickle.HIGHEST_PROTOCOL)
    file_changes_update(blob_path, blob_tmp)


if __name__ == '__main__':
    main(sys.argv)
