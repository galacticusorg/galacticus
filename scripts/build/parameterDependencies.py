#!/usr/bin/env python3
"""Emit `<target>.parameters.F90` -- a Fortran subroutine that enumerates
every input parameter a given executable can accept.

Reads `<target>.d` to discover which object files participate in the
executable, finds each object file's corresponding F90/cpp source, scans
it (plus any pre-processed `<stem>.p.F90` variant and the per-functionClass
`<stem>.p` parameter listings the source-tree processor drops into
`$BUILDPATH`) for `<inputParameter>` and `<objectBuilder>` directives, and
aggregates every parameter name into a de-duplicated, sorted list.

Outputs:
  * `$BUILDPATH/<target>.parameters.F90` -- a `knownParameterNames`
    subroutine allocating a `varying_string` array of the names.
  * `$BUILDPATH/<target>.blob` -- pickle cache of per-file parameter sets.

Andrew Benson (2026).
"""

import os
import pickle
import re
import sys


from Galacticus.Build.Directives import extract_directives
from Galacticus.Build.ScanCache import (
    file_identifier as _file_identifier,
    load_cache      as _load_cache,
)


# Include files to exclude from the parameter search (third-party code).
_EXCLUDED_INCLUDES = frozenset({'fftw3.f03'})


def _object_files_from_dep(dependency_file_name, build_path):
    """Extract every `<build>/<name>.o` path from `<target>.d` and return the
    bare `<name>.o` list.
    """
    if not os.path.exists(dependency_file_name):
        return []
    prefix_re = re.compile(r'^' + re.escape(build_path) + r'/(.+\.o)$')
    out = []
    with open(dependency_file_name, 'r') as fh:
        for line in fh:
            m = prefix_re.match(line.rstrip('\n'))
            if m:
                out.append(m.group(1))
    return out


def _slurp_lines_if_present(path):
    """Return the list of chomped lines from `path`, or `[]` if it does not
    exist.
    """
    if not os.path.exists(path):
        return []
    with open(path, 'r', errors='replace') as fh:
        return fh.read().splitlines()


_INCLUDE_LINE_RE = re.compile(r"""^\s*include\s*['"]([^'"]+)['"]""")


def _find_included_files(file_path, build_path):
    """Return every `<build>/<leaf>` referenced via `include '<leaf>'` in
    `file_path`, skipping any leaf in `_EXCLUDED_INCLUDES`.
    """
    includes = []
    for line in _slurp_lines_if_present(file_path):
        m = _INCLUDE_LINE_RE.match(line)
        if m and m.group(1) not in _EXCLUDED_INCLUDES:
            includes.append(os.path.join(build_path, m.group(1)))
    return includes


def main(argv):
    if len(argv) != 3:
        print("Usage: parameterDependencies.py <sourceDirectory> <target>",
              file=sys.stderr)
        sys.exit(1)
    source_directory_name = argv[1]
    target_name           = argv[2]
    build_path            = os.environ['BUILDPATH']

    # Work out the per-target side-file paths.
    dep_name = re.sub(r'\.(exe|o)$', '.d', target_name)
    blob_name = re.sub(r'\.(exe|o)$', '.blob', target_name)
    out_name  = re.sub(r'\.(exe|o)$', '.parameters.F90', target_name)

    dependency_file_name = os.path.join(build_path, dep_name)
    blob_path            = os.path.join(build_path, blob_name)
    output_path          = os.path.join(build_path, out_name)

    object_files = set(_object_files_from_dep(dependency_file_name,
                                              build_path))

    parameters_per_file, cache_mtime = _load_cache(blob_path)
    seen_identifiers = set()

    source_root = os.path.join(source_directory_name, 'source')
    if not os.path.isdir(source_root):
        sys.exit("parameterDependencies.py: can not open the source "
                 "directory: " + source_root)
    source_names = []
    for dirpath, dirnames, filenames in os.walk(source_root):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith('.'))
        for f in filenames:
            source_names.append(
                os.path.relpath(os.path.join(dirpath, f), source_root))
    source_names.sort()

    for file_name in source_names:
        if os.path.basename(file_name).startswith('.#'):
            continue
        if not re.search(r'\.(F90|cpp)$', file_name):
            continue
        object_file_name = re.sub(r'\.(F90|cpp)$', '.o', file_name)
        if object_file_name not in object_files:
            continue

        # Prefer the pre-processed `.p.F90` variant when present.
        root_file_name = re.sub(
            r'\.F90$', '.',
            os.path.join(build_path, file_name),
        )
        if file_name.endswith('.F90') and os.path.exists(
                root_file_name + 'p.F90'):
            file_stack = [root_file_name + 'p.F90']
        else:
            file_stack = [os.path.join(source_root, file_name)]

        file_identifier = _file_identifier(file_stack[0])
        seen_identifiers.add(file_identifier)

        rescan = True
        if file_identifier in parameters_per_file:
            tracked = parameters_per_file[file_identifier].get('files') or []
            stale = any(
                os.path.exists(t) and os.stat(t).st_mtime > cache_mtime
                for t in tracked
            )
            rescan = bool(stale)
        if not rescan:
            continue

        parameters_per_file.pop(file_identifier, None)
        entry = parameters_per_file.setdefault(
            file_identifier, {'files': list(file_stack), 'parameter': []},
        )

        # Pull any pre-emitted `.p` parameter-list file (one parameter
        # name per line, dropped by the FunctionClass source-tree
        # processor).  Only for .F90 source.  Track the `.p` file itself
        # so a change to it alone can trigger a rescan of this entry.
        if file_name.endswith('.F90'):
            p_list_path = root_file_name + 'p'
            entry['files'].append(p_list_path)
            entry['parameter'].extend(
                _slurp_lines_if_present(p_list_path)
            )

        while file_stack:
            file_to_process = file_stack.pop(0)
            includes = _find_included_files(file_to_process, build_path)
            file_stack.extend(includes)
            entry['files'].extend(includes)

            # inputParameter directives contribute their `name` attribute.
            for d in extract_directives(file_to_process, 'inputParameter'):
                if 'name' in d:
                    entry['parameter'].append(d['name'])

            # objectBuilder directives contribute `parameterName` when it
            # differs from the default `class` naming.
            for d in extract_directives(file_to_process, 'objectBuilder'):
                pname = d.get('parameterName')
                if pname is not None and pname != d.get('class'):
                    entry['parameter'].append(pname)

    # Prune cache entries for files that no longer participate in this target
    # (deleted, renamed, dropped from the `.d` list, or whose identifier
    # flipped between the source and preprocessed path). Without this, their
    # parameters would keep contributing to `knownParameterNames` forever,
    # weakening the executable's unknown-parameter detection.
    parameters_per_file = {
        k: v for k, v in parameters_per_file.items() if k in seen_identifiers
    }

    # Aggregate + deduplicate (sort + unique).
    all_parameters = []
    for entry in parameters_per_file.values():
        all_parameters.extend(entry.get('parameter') or [])
    all_parameters = sorted(set(all_parameters))

    with open(output_path, 'w') as fh:
        fh.write("subroutine knownParameterNames(names)\n")
        fh.write("  use ISO_Varying_String\n")
        fh.write("  implicit none\n")
        fh.write("  type(varying_string), dimension(:), allocatable, "
                 "intent(inout) :: names \n")
        fh.write(f"  allocate(names({len(all_parameters)}))\n")
        for i, name in enumerate(all_parameters):
            fh.write(f"  names({i + 1})='{name}'\n")
        fh.write("end subroutine knownParameterNames\n")

    with open(blob_path, 'wb') as fh:
        pickle.dump(parameters_per_file, fh,
                    protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    main(sys.argv)
