#!/usr/bin/env python3
import os
import re
import sys

# Locate files which contain programs and append to a list of executables.
# Andrew Benson (ported to Python 2026)

if len(sys.argv) != 2:
    print("Usage: findExecutables.py <galacticusDirectory>", file=sys.stderr)
    sys.exit(1)

galacticus_directory = sys.argv[1]
source_directory     = os.path.join(galacticus_directory, "source")
build_path           = os.environ['BUILDPATH']
work_dir             = build_path.rstrip('/') + '/'

executable_names = []

def _source_files_recursive(root):
    """Yield every file path under `root`, relative to `root`, walking
    subdirectories to any depth."""
    rel = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith('.'))
        for fn in filenames:
            rel.append(os.path.relpath(os.path.join(dirpath, fn), root))
    return sorted(rel)


with open(os.path.join(build_path, "Makefile_All_Execs"), 'w') as out:
    for file_name in _source_files_recursive(source_directory):
        # Skip temporary files.
        if os.path.basename(file_name).startswith('.#'):
            continue
        # Skip vendored third-party code: it may carry its own test `program`
        # units (e.g. Genz's `PROGRAM TSTNRM`) that are not Galacticus
        # executables.  The flat, pre-hierarchy scan never reached these
        # subdirectories; preserve that by excluding `external/`.
        if file_name.split('/', 1)[0] == 'external':
            continue
        # Only Fortran source files.
        if not re.search(r'\.[fF](90)?$', file_name):
            continue

        file_full           = os.path.join(source_directory, file_name)
        exclude_from_all    = False
        found_program       = False

        in_doc_block = False

        try:
            with open(file_full, 'r', errors='replace') as fh:
                for line in fh:
                    # Skip lines inside `!!{ ... !!}` LaTeX documentation blocks,
                    # whose unprefixed prose can otherwise look like a `program`
                    # statement.
                    if in_doc_block:
                        if '!!}' in line:
                            in_doc_block = False
                        continue
                    if re.match(r'^\s*!!\{', line) and '!!}' not in line.split('!!{', 1)[1]:
                        in_doc_block = True
                        continue
                    if re.match(r'^\s*!/\s+exclude', line):
                        exclude_from_all = True
                    if re.match(r'^\s*program\s', line, re.IGNORECASE):
                        found_program = True
                        break
        except OSError:
            continue

        if not found_program:
            continue

        # `obj_root` is the path-based stem (relative to source/, e.g.
        # `tests/nodes`) used for every build artifact, which mirrors the
        # source hierarchy.  `exe_root` is the historical flat, dot-separated
        # name (e.g. `tests.nodes`) used for the user-facing executable so that
        # `make tests.nodes.exe`, CI matrices, and test harnesses keep working
        # unchanged after the source tree was made hierarchical.
        obj_root = re.sub(r'\.[fF](90)?t?$', '', file_name)
        exe_root = obj_root.replace('/', '.')

        if not exclude_from_all:
            executable_names.append(exe_root + '.exe')

        rule = f"""\
{exe_root}.exe: {work_dir}{obj_root}.o {work_dir}{obj_root}.d $(MAKE_DEPS) $(UPDATE_DEPS)
\t./scripts/build/parameterDependencies.py `pwd` {obj_root}.exe
\t$(FCCOMPILER) -c {work_dir}{obj_root}.parameters.F90 -o {work_dir}{obj_root}.parameters.o $(FCFLAGS)
\t@if echo "$(MAKEFLAGS)" | grep -q -E -- ' -j1( |$$)'; then \\
\t useLocks=no; \\
\telif echo "$(MAKEFLAGS)" | grep -q -E -- ' -j( |$$)'; then \\
\t useLocks=$(LOCKMD5); \\
\telif echo "$(MAKEFLAGS)" | grep -q -E -- ' -j[0-9]+( |$$)'; then \\
\t useLocks=$(LOCKMD5); \\
\telse \\
\t useLocks=no; \\
\tfi; \\
\t./scripts/build/sourceDigests.py `pwd` {obj_root}.exe $$useLocks
\t$(CCOMPILER) -c {work_dir}{obj_root}.md5s.c -o {work_dir}{obj_root}.md5s.o $(CFLAGS)
\t$(FCCOMPILER) `cat {work_dir}{obj_root}.d` {work_dir}{obj_root}.parameters.o {work_dir}{obj_root}.md5s.o -o {exe_root}.exe$(SUFFIX) $(FCFLAGS) $(FCFLAGS_LINK) `./scripts/build/libraryDependencies.py {obj_root}.exe $(FCFLAGS)` 2>&1 | ./scripts/build/postprocessLinker.py

"""
        out.write(rule)

    if executable_names:
        out.write("all_exes = " + " ".join(executable_names) + "\n")
