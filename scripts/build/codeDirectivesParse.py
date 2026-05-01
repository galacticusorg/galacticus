#!/usr/bin/env python3
"""Scan source code for `!![...!!]` directives and emit the build-time catalog.

Writes, under `$BUILDPATH`:
    * `directiveLocations.xml` -- XML map of non-include directive name ->
      list of source files that contain the directive.
    * `Makefile_Directives`    -- build rules for every include directive
      (`<file>.Inc.up` / `<file>.Inc` targets) plus implicit dependencies
      for functionClass-preprocessed files and `Makefile_Use_Dependencies`.
    * `<directive>.xml`        -- one file per include directive, atomically
      updated only when its XML content actually changes.
    * `codeDirectives.blob`    -- pickle cache of per-file directive info
      so subsequent runs skip unchanged sources.

Mirrors scripts/build/codeDirectivesParse.pl.
Andrew Benson (ported to Python 2026).
"""

import os
import pickle
import re
import sys


from build.file_changes               import update as file_changes_update
from Galacticus.Build.Directives      import extract_directives
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

def _load_cache(blob_path):
    """Return the per-file cache dict, or `{}` if no readable cache exists.

    A pre-existing Perl Storable blob at `blob_path` will fail to unpickle;
    we catch that and return an empty cache, forcing a full rescan.  Matches
    the `$havePerFile` logic in the Perl script.
    """
    if not os.path.exists(blob_path):
        return {}
    try:
        with open(blob_path, 'rb') as fh:
            cache = pickle.load(fh)
    except (pickle.UnpicklingError, EOFError, AttributeError, ValueError,
            ImportError, ModuleNotFoundError):
        return {}
    return cache if isinstance(cache, dict) else {}


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
_INCLUDE_BODY_RE  = re.compile(r"""^\s*\#?\s*include\s*["'<]([^"'>]+)["'>]""",
                               re.IGNORECASE | re.MULTILINE)
_SOURCE_SUFFIX_RE = re.compile(r'\.(f|f90|c|cpp|h)$', re.IGNORECASE)


def _collect_included_files(file_path, source_directory):
    """Return the list of `include '<leaf>'` files referenced by `file_path`
    that exist under `source_directory` (after `.inc` -> `.Inc` fix-up).

    Mirrors the `map { … include … }` block at codeDirectivesParse.pl:84-92.
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

    Mirrors addImplicitDirectives() at codeDirectivesParse.pl:214-242.
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
# Main
# ---------------------------------------------------------------------------

def main(argv):
    if len(argv) != 2:
        _usage()
    install_directory = argv[1]
    build_path        = os.environ['BUILDPATH']
    blob_path         = os.path.join(build_path, 'codeDirectives.blob')

    # Load per-file cache and capture its mtime for freshness checks.
    directives_per_file = _load_cache(blob_path)
    have_per_file       = bool(directives_per_file) and os.path.exists(blob_path)
    cache_mtime         = os.stat(blob_path).st_mtime if have_per_file else None

    # List source files under `<installDir>/source/`.
    source_directory   = os.path.join(install_directory, 'source')
    source_file_names  = sorted(
        f for f in os.listdir(source_directory)
        if _SOURCE_SUFFIX_RE.search(f) and not f.startswith('.#')
    )

    # Build the UNSTRIPPED identifier list used for the new/removed-file
    # checks below.  Matches the first `@fileIdentifiers` loop in the Perl.
    file_identifiers = [
        (source_directory + '/' + name).replace('/', '_')
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

    # Iterate over source files.
    for source_file_name in source_file_names:
        file_path       = source_directory + '/' + source_file_name
        file_identifier = file_path.replace('/', '_')
        # Strip leading `.` (optionally followed by `_`) -- mirrors the
        # Perl `s/^\._??//`.  A no-op for absolute install paths.
        file_identifier = re.sub(r'^\._?', '', file_identifier)

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

        # Drop stale cache entry and reinitialise.
        directives_per_file.pop(file_identifier, None)
        entry = directives_per_file.setdefault(
            file_identifier,
            {'files': [file_path]},
        )

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
                if root_type == 'include':
                    # Include directive: record the source file and the
                    # include file name, drop `content`, and XML-serialize
                    # the remaining attributes for the per-directive file.
                    directive['source'] = current
                    content = directive.get('content', '')
                    if isinstance(content, str):
                        m = _INCLUDE_BODY_RE.search(content)
                        if m:
                            include_leaf = m.group(1)
                            include_leaf = re.sub(r'\.inc$', '.Inc',
                                                  include_leaf)
                            directive['fileName'] = os.path.join(
                                build_path, include_leaf,
                            )
                    directive.pop('content', None)

                    directive_name = (
                        directive.get('name')
                        or directive.get('directive')
                    )
                    key = f"{directive_name}.{directive.get('type')}"
                    entry.setdefault('includeDirectives', {})[key] = {
                        'source':   current,
                        'fileName': directive.get('fileName'),
                        'xml':      dict_to_xml_string(root_type, directive),
                    }
                else:
                    # Non-include directive: remember the source file.
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

    # -----------------------------------------------------------------------
    # Reduce across files.
    # -----------------------------------------------------------------------
    include_directives      = {}
    non_include_directives  = {}
    function_classes        = {}
    for entry in directives_per_file.values():
        for k, v in (entry.get('includeDirectives') or {}).items():
            include_directives[k] = v
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
    # Makefile_Directives plus per-directive XML files.
    # -----------------------------------------------------------------------
    makefile_path = os.path.join(build_path, 'Makefile_Directives')
    with open(makefile_path, 'w') as mk:
        for directive in sorted(include_directives):
            info      = include_directives[directive]
            file_name = re.sub(r'\.inc$', '.Inc', info['fileName'])

            # Extra dependencies per directive-key suffix.
            extra_deps = []
            m = re.match(r'^([a-zA-Z0-9_]+)\.function$', directive)
            if m:
                base = m.group(1)
                extra_deps.extend(
                    sorted(non_include_directives.get(base, {}).get('files', []))
                )
            m = re.match(r'^([a-zA-Z0-9_]+)\.(moduleUse|functionCall)$',
                         directive)
            if m:
                base = m.group(1)
                deps = non_include_directives.get(base, {}).get('dependency')
                if deps:
                    extra_deps.extend(deps)

            directive_xml = os.path.join(build_path, directive + '.xml')
            mk.write(
                f"{file_name}.up: {directive_xml} {' '.join(extra_deps)}"
                " $(BUILDPATH)/hdf5FCInterop.dat"
                " $(BUILDPATH)/openMPCriticalSections.xml\n"
            )
            mk.write(
                f"\t./scripts/build/buildCode.py {install_directory}"
                f" {directive_xml}\n"
            )
            mk.write(f"{file_name}: {file_name}.up\n\n")

            # Per-directive XML, atomic update.
            tmp = directive_xml + '.tmp'
            with open(tmp, 'w') as fh:
                fh.write(info['xml'])
            file_changes_update(directive_xml, tmp)

        # Extra dependencies: object files of functionClass-bearing sources
        # pick up every other source that carries the same directive.
        for directive_name in sorted(function_classes):
            preprocessed = function_classes[directive_name]
            deps = sorted(
                non_include_directives.get(directive_name, {}).get('files', [])
            )
            mk.write(f"{preprocessed}.up: {' '.join(deps)}\n\n")

        # Makefile_Use_Dependencies must be rebuilt after the include files
        # have been constructed.
        use_deps = sorted(
            re.sub(r'\.inc$', '.Inc', info['fileName'])
            for info in include_directives.values()
        )
        mk.write(
            os.path.join(build_path, 'Makefile_Use_Dependencies')
            + ": " + ' '.join(use_deps) + "\n\n"
        )

        # Makefile_Component_Includes include (must live here because it
        # depends on objects.nodes.components.Inc built via this Makefile).
        mk.write(
            f"-include {os.path.join(build_path, 'Makefile_Component_Includes')}\n"
        )
        mk.write(
            os.path.join(build_path, 'Makefile_Component_Includes')
            + ": "
            + os.path.join(build_path, 'objects.nodes.components.Inc')
            + "\n\n"
        )

    # -----------------------------------------------------------------------
    # Persist per-file cache.
    # -----------------------------------------------------------------------
    _save_cache(blob_path, directives_per_file)


if __name__ == '__main__':
    main(sys.argv)
