#!/usr/bin/env python3
import os
import re
import subprocess
import sys

from Galacticus.Build.ParallelScan import scan as parallel_scan

# Construct Makefile rules for source files which have dependencies on include files.
# Andrew Benson (ported to Python 2026)

# ---------------------------------------------------------------------------
# Read-only data shared with forked workers via copy-on-write. These are
# populated in `main()` BEFORE `parallel_scan` is called, so each forked worker
# inherits them without anything being pickled per task.
# ---------------------------------------------------------------------------
_WORKER = {}


def _source_files_recursive(root):
    """Yield every file path under `root`, relative to `root`, walking
    subdirectories to any depth."""
    rel = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith('.'))
        for fn in filenames:
            rel.append(os.path.relpath(os.path.join(dirpath, fn), root))
    return sorted(rel)


def _c_include_dirs(source_directory):
    """Build a list of pre-existing C include file paths: source dir, system
    GCC paths, and any paths from GALACTICUS_CFLAGS / GALACTICUS_CPPFLAGS."""
    c_include_dirs = [source_directory]
    try:
        cpp_output = subprocess.run(
            ["cpp", "-Wp,-v"],
            input="",
            capture_output=True,
            text=True
        ).stderr
        for line in cpp_output.splitlines():
            m = re.match(r'^\s*(/\S+)', line)
            if m:
                c_include_dirs.append(m.group(1))
    except Exception:
        pass

    for env_var in ("GALACTICUS_CFLAGS", "GALACTICUS_CPPFLAGS"):
        val = os.environ.get(env_var, "")
        for token in val.split():
            m = re.match(r'^-I(.*)', token)
            if m:
                c_include_dirs.append(m.group(1))
    return c_include_dirs


def _scan_one(file_name):
    """Worker: read and parse one source file, returning its contribution as a
    `(rule_text, dependency_file_names)` tuple. Writes NOTHING to shared state.
    Mirrors the per-file body of the original serial loop exactly. The read-only
    `source_directory`, `build_path` and `c_include_dirs` are inherited from the
    `_WORKER` globals via copy-on-write fork.
    """
    source_directory = _WORKER['source_directory']
    build_path       = _WORKER['build_path']
    c_include_dirs   = _WORKER['c_include_dirs']

    file_full = os.path.join(source_directory, file_name)
    included_files = []

    try:
        with open(file_full, 'r', errors='replace') as fh:
            for line in fh:
                m = re.match(
                    r'''^\s*#??include\s*[<'"](.+\.(inc|h)\d*)['">](\s*!\s*NO_USES)?''',
                    line
                )
                if not m:
                    continue
                inc_name  = m.group(1)
                ext       = m.group(2)
                no_uses   = m.group(3) is not None
                if ext == 'inc':
                    # The `automatic` flag records a `! NO_USES` marker, which
                    # opts a generated include out of the auto-generated
                    # dependency chain fed to Makefile_Use_Dependencies below.
                    included_files.append({'fileName': inc_name, 'automatic': not no_uses})
                elif ext == 'h':
                    # Only include .h files not already present in source or system dirs.
                    in_source = os.path.exists(os.path.join(source_directory, inc_name))
                    in_system = any(
                        os.path.exists(os.path.join(d, inc_name)) for d in c_include_dirs
                    )
                    if not in_source and not in_system:
                        included_files.append({'fileName': inc_name})
    except OSError:
        return "", []

    rule_text = ""
    if included_files:
        obj_name = re.sub(r'\.[^.]+$', '.o', file_name)
        deps     = ' '.join(sorted(
            os.path.join(build_path, f['fileName']) for f in included_files
        ))
        rule_text += f"{os.path.join(build_path, obj_name)}: {deps}\n\n"
        if file_name.endswith('.F90'):
            pp_name = file_name.replace('.F90', '.p.F90')
            rule_text += f"{os.path.join(build_path, pp_name)}: {deps}\n\n"

    # Collect auto-generated .inc dependencies for Makefile_Use_Dependencies.
    dependency_file_names = []
    for inc in included_files:
        if inc['fileName'].endswith('.inc') and inc.get('automatic'):
            unpp = re.sub(r'\.inc$', '.Inc', inc['fileName'])
            if os.path.exists(os.path.join(source_directory, unpp)):
                dependency_file_names.append(unpp)
            else:
                dependency_file_names.append(inc['fileName'])

    return rule_text, dependency_file_names


def main(argv):
    if len(argv) != 2:
        print("Usage: includeDependencies.py <galacticusDirectory>", file=sys.stderr)
        sys.exit(1)

    galacticus_directory = argv[1]
    source_directory     = os.path.join(galacticus_directory, "source")
    build_path           = os.environ['BUILDPATH']

    c_include_dirs = _c_include_dirs(source_directory)

    # Build the ordered list of files to process, preserving the original loop's
    # order and skip rules (temporary files, non source/include extensions).
    tasks = []
    for file_name in _source_files_recursive(source_directory):
        # Skip temporary files.
        if os.path.basename(file_name).startswith('.#'):
            continue
        # Only process Fortran, Fortran include, C, C++, and header files.
        if not re.search(r'\.(f(90)?|inc|c(pp)?|h)$', file_name, re.IGNORECASE):
            continue
        tasks.append(file_name)

    # Publish read-only data for the forked workers, then scan all files
    # concurrently. Results come back in the SAME order as `tasks`.
    _WORKER['source_directory'] = source_directory
    _WORKER['build_path']       = build_path
    _WORKER['c_include_dirs']   = c_include_dirs

    results = parallel_scan(tasks, _scan_one, "includeDependencies.py")

    # Write the output Makefile SEQUENTIALLY, consuming results in order so the
    # per-file rule lines keep their original ordering.
    makefile_path = os.path.join(build_path, "Makefile_Include_Dependencies")
    dependency_file_names = []
    with open(makefile_path, 'w') as makefile:
        for rule_text, dep_names in results:
            if rule_text:
                makefile.write(rule_text)
            dependency_file_names.extend(dep_names)

        # Dependency line for Makefile_Use_Dependencies.
        dep_targets = ' '.join(
            os.path.join(build_path, f) for f in sorted(set(dependency_file_names))
        )
        makefile.write(
            f"\n{os.path.join(build_path, 'Makefile_Use_Dependencies')}: {dep_targets}\n"
        )


if __name__ == '__main__':
    main(sys.argv)
