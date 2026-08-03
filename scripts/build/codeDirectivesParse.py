#!/usr/bin/env python3
"""Scan source code for `!![...!!]` directives and emit the build-time catalog.

Writes, under `$BUILDPATH`:
    * `directiveLocations.xml` -- XML map of directive name -> list of source
      files that contain the directive.
    * `Makefile_Directives`    -- extra prerequisites for functionClass- and
      componentBuilder-preprocessed files, plus the ordering of
      `Makefile_Use_Dependencies` and `Makefile_Component_Includes` after the
      preprocessed `_class` file.
    * `codeDirectives.blob`    -- pickle cache of per-file directive info
      so subsequent runs skip unchanged sources.

Andrew Benson (2026).
"""

import os
import pickle
import re
import sys


from Galacticus.Build.FileChanges               import update as file_changes_update
from Galacticus.Build.ScanCache                 import (
    file_identifier as _file_identifier,
    load_cache      as _load_cache,
    prune           as _prune_cache,
)
from Galacticus.Build.Directives      import extract_directives
from Galacticus.Build.ParallelScan    import scan as parallel_scan
from XML.Utils                        import dict_to_xml_string


# ---------------------------------------------------------------------------
# Command-line / environment wiring
# ---------------------------------------------------------------------------

def _usage():
    print("Usage: codeDirectivesParse.py <installDirectory>", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
# Pickle cache helpers
# ---------------------------------------------------------------------------

def _save_cache(blob_path, cache):
    """Write `cache` to `blob_path` via a temp file + atomic replace (only when
    the pickled bytes have actually changed).
    """
    tmp_path = blob_path + '.tmp'
    with open(tmp_path, 'wb') as fh:
        pickle.dump(cache, fh, protocol=pickle.HIGHEST_PROTOCOL)
    file_changes_update(blob_path, tmp_path)


# ---------------------------------------------------------------------------
# Directive extraction helpers
# ---------------------------------------------------------------------------

_INCLUDE_LINE_RE  = re.compile(r"""^\s*include\s*['"]([^'"]+)['"]\s*$""")
_SOURCE_SUFFIX_RE = re.compile(r'\.(f|f90|c|cpp|h)$', re.IGNORECASE)


def _collect_included_files(file_path, source_directory):
    """Return the list of `include '<leaf>'` files referenced by `file_path`
    that exist under `source_directory` (after `.inc` -> `.Inc` fix-up).
    """
    included = []
    with open(file_path, 'r', errors='replace') as fh:
        for line in fh:
            m = _INCLUDE_LINE_RE.match(line)
            if not m:
                continue
            leaf = m.group(1)
            candidate = os.path.join(source_directory, leaf)
            candidate = re.sub(r'\.inc$', '.Inc', candidate)
            if os.path.exists(candidate):
                included.append(candidate)
    return included


def _add_implicit_directives(directive, per_file_entry, file_name, file_path):
    """For functionClass directives, inject implicit `stateful` /
    `functionClassDestroy` task dependencies into `per_file_entry`.
    """
    is_function_class = directive.get('rootElementType') == 'functionClass'
    implicit_map = {
        'stateful': {
            'always': is_function_class,
            'tasks':  ['galacticusStateRetrieveTask',
                       'galacticusStateStoreTask'],
        },
        'functionClassDestroy': {
            'always': is_function_class and 'functionClassDestroy' not in directive,
            'tasks':  ['functionClassDestroyTask'],
        },
    }
    non_include = per_file_entry.setdefault('nonIncludeDirectives', {})
    for key, spec in implicit_map.items():
        if spec['always'] or directive.get(key) == 'yes':
            for task in spec['tasks']:
                entry = non_include.setdefault(task, {'files': [], 'dependency': []})
                entry['files'].append(file_path)
                entry['dependency'].append(file_name)


# ---------------------------------------------------------------------------
# Per-file processing (parallelised across files)
# ---------------------------------------------------------------------------
#
# Each source file is scanned independently: it reads the file (and any
# `include`d files) and produces a single per-file `entry` dict that writes to
# no shared state. The per-file cost is dominated by blocking file I/O
# (open/read of the source on NFS), which on a loaded node stretches from
# microseconds to milliseconds per file; running the scans concurrently overlaps
# those waits (and the directive-parsing CPU work) instead of paying them one at
# a time.
#
# A fork-based process pool is used rather than threads so the workers are fully
# isolated -- no assumptions are made about the thread-safety of the parser.
# Workers inherit the read-only `source_directory`/`build_path` via copy-on-write
# fork (published to `_WORKER` before the scan), so nothing large is pickled per
# task; only the small per-file result dicts come back.

_WORKER = {}


def _scan_one(task):
    """Worker: scan one source file (and its include tree) and return
    `(file_identifier, entry)`. Mirrors the per-file body of the serial scan
    loop exactly, writing to no shared state.
    """
    source_file_name, file_path, file_identifier = task
    source_directory = _WORKER['source_directory']
    build_path       = _WORKER['build_path']

    entry = {'files': [file_path]}

    # Walk include files depth-first.
    pending = [file_path]
    while pending:
        current = pending.pop()
        included = _collect_included_files(current, source_directory)
        pending.extend(included)
        entry['files'].extend(included)

        # Extract every directive in the current file.
        for directive in extract_directives(
            current, '*', set_root_element_type=True,
        ):
            root_type = directive.get('rootElementType')
            # Remember the source file for every directive.  (The `<include>`
            # directive mechanism was retired when the `component` builder moved
            # into the SourceTree preprocessor — see
            # Galacticus.Build.SourceTree.Process.ComponentBuilder — so every
            # directive is now handled the same way.)
            non_include = entry.setdefault('nonIncludeDirectives', {})
            slot = non_include.setdefault(
                root_type, {'files': [], 'dependency': []},
            )
            slot['files'].append(file_path)

            if root_type == 'functionClass':
                preprocessed = re.sub(
                    r'\.F90$', '.p.F90',
                    os.path.join(build_path, source_file_name),
                )
                entry.setdefault(
                    'functionClasses', {},
                )[directive['name']] = preprocessed
                _add_implicit_directives(
                    directive, entry,
                    preprocessed, preprocessed,
                )

    return file_identifier, entry


# Literal (compile-time-interned) strings used as dict keys / fixed values when
# building a per-file `entry` in `_scan_one`. When the serial loop builds entries
# in-process these literals are a single interned object shared across every
# entry, so `pickle` writes them once and back-references thereafter. Entries
# returned from forked workers are unpickled into fresh, de-interned objects, so
# without help the saved blob would grow (lost back-references) and differ
# byte-for-byte from the serial blob. `_canonicalize` re-interns exactly these
# strings (only) so the merged structure pickles identically -- interning more
# (e.g. every identifier-like string) would wrongly share parser-derived strings
# that the serial run keeps distinct. Keep this set in sync with the literal
# keys/values produced by `_scan_one`.
_ENTRY_LITERALS = {
    'files',
    'nonIncludeDirectives', 'dependency', 'functionClasses',
    'galacticusStateRetrieveTask', 'galacticusStateStoreTask',
    'functionClassDestroyTask',
}
_ENTRY_CANON = {s: s for s in _ENTRY_LITERALS}


def _canonicalize(obj):
    """Return a copy of `obj` with the `_ENTRY_LITERALS` strings replaced by
    their canonical interned instances (preserving dict insertion order and all
    other object identities), so a worker-returned entry pickles identically to
    one the serial loop built in-process. A no-op in the serial path (those
    strings are already the interned literals), so output is unchanged there.
    """
    if isinstance(obj, dict):
        return {
            _ENTRY_CANON.get(k, k) if isinstance(k, str) else k:
                _canonicalize(v)
            for k, v in obj.items()
        }
    if isinstance(obj, list):
        return [_canonicalize(x) for x in obj]
    if isinstance(obj, str):
        return _ENTRY_CANON.get(obj, obj)
    return obj


def _scan_files(scan_list, source_directory, build_path):
    """Scan every file in `scan_list` (a list of
    `(source_file_name, file_path, file_identifier)` preserving the serial
    loop's order) and return the results in that SAME order, so the merge into
    `directives_per_file` -- and therefore the cross-file reduction and emitted
    Makefile -- is identical to the serial version. The read-only
    `source_directory`/`build_path` are published to `_WORKER` first so forked
    workers inherit them via copy-on-write instead of re-pickling per task.
    """
    _WORKER['source_directory'] = source_directory
    _WORKER['build_path']       = build_path
    return parallel_scan(scan_list, _scan_one, 'codeDirectivesParse.py')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv):
    if len(argv) != 2:
        _usage()
    install_directory = argv[1]
    build_path        = os.environ['BUILDPATH']
    blob_path         = os.path.join(build_path, 'codeDirectives.blob')

    # Load per-file cache and capture its mtime for freshness checks.
    directives_per_file, cache_mtime = _load_cache(blob_path)
    have_per_file = bool(directives_per_file)

    # List source files under `<installDir>/source/`, recursing into
    # subdirectories.  Names are kept relative to `source/`.
    source_directory   = os.path.join(install_directory, 'source')
    source_file_names  = []
    for dirpath, dirnames, filenames in os.walk(source_directory):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith('.'))
        for f in filenames:
            if _SOURCE_SUFFIX_RE.search(f) and not f.startswith('.#'):
                source_file_names.append(
                    os.path.relpath(os.path.join(dirpath, f), source_directory))
    source_file_names.sort()

    # Build the identifier list used for the new/removed-file checks and the
    # stale-entry prune below. Canonicalized via `file_identifier()` so it uses
    # the SAME key form as the cache entries; were it built unstripped it would
    # agree with the cache keys only for the absolute paths make passes, and a
    # `./`-relative invocation would make every cached key look removed — and,
    # with pruning, silently discard the whole cache.
    file_identifiers = [
        _file_identifier(source_directory + '/' + name)
        for name in source_file_names
    ]
    file_identifier_set = set(file_identifiers)

    # Force a global rescan if any file was added or removed since the cache
    # was written.
    force_rescan = False
    if any(fid not in directives_per_file for fid in file_identifiers):
        force_rescan = True
    if any(fid not in file_identifier_set for fid in directives_per_file):
        force_rescan = True

    # Decide which files need a (re)scan, preserving the deterministic sorted
    # `source_file_names` order so the merge below reproduces the serial
    # insertion order (and hence the cross-file reduction) exactly.
    scan_list = []
    for source_file_name in source_file_names:
        file_path       = source_directory + '/' + source_file_name
        file_identifier = _file_identifier(file_path)

        # Decide whether this file needs a rescan.
        rescan = True
        if have_per_file and file_identifier in directives_per_file:
            tracked = directives_per_file[file_identifier].get('files') or []
            stale   = any(
                os.path.exists(t) and os.stat(t).st_mtime > cache_mtime
                for t in tracked
            )
            rescan = bool(stale)
        if not (rescan or force_rescan):
            continue

        scan_list.append((source_file_name, file_path, file_identifier))

    # Scan the files (concurrently) and merge results in scan-list order so the
    # output is identical to a serial run: each rescanned file's stale cache
    # entry is dropped and the freshly-built entry reinserted in the same order
    # the serial loop would have produced.
    for file_identifier, entry in _scan_files(
            scan_list, source_directory, build_path):
        directives_per_file.pop(file_identifier, None)
        directives_per_file[file_identifier] = _canonicalize(entry)

    # Drop cache entries for source files that no longer exist. The rescan
    # above only ever replaces entries for surviving files, so without this a
    # deleted file's directives would be re-emitted into directiveLocations.xml
    # and Makefile_Directives forever — a phantom make dependency on a file
    # that is gone, which survives even deleting the generated fragments
    # (they are rebuilt from this same cache).
    _prune_cache(directives_per_file, file_identifier_set)

    # -----------------------------------------------------------------------
    # Reduce across files.
    # -----------------------------------------------------------------------
    non_include_directives  = {}
    function_classes        = {}
    for entry in directives_per_file.values():
        for directive, data in (entry.get('nonIncludeDirectives') or {}).items():
            merged = non_include_directives.setdefault(
                directive, {'files': [], 'dependency': []},
            )
            for key in ('files', 'dependency'):
                merged[key].extend(data.get(key) or [])
        for k, v in (entry.get('functionClasses') or {}).items():
            function_classes[k] = v

    # Uniquify (sort + dedupe) every list.
    for data in non_include_directives.values():
        for key in ('files', 'dependency'):
            if data.get(key):
                data[key] = sorted(set(data[key]))

    # -----------------------------------------------------------------------
    # directiveLocations.xml (atomic, only-if-changed).
    # -----------------------------------------------------------------------
    output = {
        name: {'file': data['files']}
        for name, data in non_include_directives.items()
    }
    locations_tmp = os.path.join(build_path, 'directiveLocations.xml.tmp')
    with open(locations_tmp, 'w') as fh:
        fh.write(dict_to_xml_string('directives', output))
    file_changes_update(
        os.path.join(build_path, 'directiveLocations.xml'),
        locations_tmp,
    )

    # -----------------------------------------------------------------------
    # Makefile_Directives.
    # -----------------------------------------------------------------------
    # Written UNCONDITIONALLY (fresh mtime on every run), NOT only-if-changed.
    # Makefile_Directives is `-include`d by the main Makefile, whose rule for
    # it runs this script; make re-executes itself (re-reading the regenerated
    # Makefile_Directives) only when the file is updated by that rule. A stable
    # mtime would suppress that restart on a clean build, so make would never
    # re-read the functionClass / componentBuilder extra-dependency lines nor
    # the ordering lines below — see the Makefile_Directives rule comment in the
    # main Makefile.
    makefile_path = os.path.join(build_path, 'Makefile_Directives')
    with open(makefile_path, 'w') as mk:
        # Extra dependencies: object files of functionClass-bearing sources
        # pick up every other source that carries the same directive.
        for directive_name in sorted(function_classes):
            preprocessed = function_classes[directive_name]
            deps = sorted(
                non_include_directives.get(directive_name, {}).get('files', [])
            )
            mk.write(f"{preprocessed}.up: {' '.join(deps)}\n\n")

        # The `componentBuilder` directive (in objects/nodes/_class.F90)
        # synthesizes the node-component class hierarchy at preprocess time (see
        # Galacticus.Build.SourceTree.Process.ComponentBuilder). Make the
        # preprocessed `_class` file depend on every `<component>` directive file
        # and on the component generator's Python sources (excluding the
        # generators' own unit tests), so editing any of them re-preprocesses
        # `_class.F90`. Also add stateStorables.xml / deepCopyActions.xml, which
        # the process pipeline reads while expanding the grafted component code
        # (the generic `%.p.F90.up` rule does not list them).
        component_builder_targets = []
        builder_files = sorted(
            non_include_directives.get('componentBuilder', {}).get('files', [])
        )
        if builder_files:
            component_files = sorted(
                non_include_directives.get('component', {}).get('files', [])
            )
            components_root = os.path.join(
                install_directory,
                'python', 'Galacticus', 'Build', 'Components',
            )
            generator_sources = []
            for root, _, files in os.walk(components_root):
                parts = root.split(os.sep)
                if '__pycache__' in parts or 'tests' in parts:
                    continue
                generator_sources.extend(
                    os.path.join(root, name)
                    for name in files if name.endswith('.py')
                )
            builder_deps = (
                component_files
                + sorted(generator_sources)
                + [
                    os.path.join(build_path, 'stateStorables.xml'),
                    os.path.join(build_path, 'deepCopyActions.xml'),
                ]
            )
            for builder_file in builder_files:
                relative     = os.path.relpath(builder_file, source_directory)
                preprocessed = os.path.join(
                    build_path, re.sub(r'\.F90$', '.p.F90', relative),
                )
                component_builder_targets.append(preprocessed)
                mk.write(f"{preprocessed}.up: {' '.join(builder_deps)}\n\n")

        component_builder_targets.sort()

        # Makefile_Use_Dependencies must be rebuilt after the component code has
        # been grafted into the preprocessed `_class` file: useDependencies scans
        # that preprocessed file to recover the generated `use` statements (of
        # the per-component `*_Data` modules) and the `include`d bound-functions
        # sources.
        if component_builder_targets:
            mk.write(
                os.path.join(build_path, 'Makefile_Use_Dependencies')
                + ": " + ' '.join(component_builder_targets) + "\n\n"
            )

        # Makefile_Component_Includes carries the dependency of `_class.o` /
        # `_class.p.F90` on the per-component bound-functions include files; it is
        # (re)written by the component generator as `_class.F90` is preprocessed,
        # so it must be built after the preprocessed `_class` file.
        mk.write(
            f"-include {os.path.join(build_path, 'Makefile_Component_Includes')}\n"
        )
        if component_builder_targets:
            mk.write(
                os.path.join(build_path, 'Makefile_Component_Includes')
                + ": " + ' '.join(component_builder_targets) + "\n\n"
            )

    # -----------------------------------------------------------------------
    # Persist per-file cache.
    # -----------------------------------------------------------------------
    _save_cache(blob_path, directives_per_file)


if __name__ == '__main__':
    main(sys.argv)
