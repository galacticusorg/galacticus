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

with open(os.path.join(build_path, "Makefile_All_Execs"), 'w') as out:
    for file_name in sorted(os.listdir(source_directory)):
        # Skip temporary files.
        if file_name.startswith('.#'):
            continue
        # Only Fortran source files.
        if not re.search(r'\.[fF](90)?$', file_name):
            continue

        file_full           = os.path.join(source_directory, file_name)
        exclude_from_all    = False
        found_program       = False

        try:
            with open(file_full, 'r', errors='replace') as fh:
                for line in fh:
                    if re.match(r'^\s*!/\s+exclude', line):
                        exclude_from_all = True
                    if re.match(r'^\s*program\s', line, re.IGNORECASE):
                        found_program = True
                        break
        except OSError:
            continue

        if not found_program:
            continue

        file_name_root = re.sub(r'\.[fF](90)?t?$', '', file_name)

        if not exclude_from_all:
            executable_names.append(file_name_root + '.exe')

        rule = f"""\
{file_name_root}.exe: {work_dir}{file_name_root}.o {work_dir}{file_name_root}.d $(MAKE_DEPS) $(UPDATE_DEPS)
\t./scripts/build/parameterDependencies.pl `pwd` {file_name_root}.exe
\t$(FCCOMPILER) -c {work_dir}{file_name_root}.parameters.F90 -o {work_dir}{file_name_root}.parameters.o $(FCFLAGS)
\t@if echo "$(MAKEFLAGS)" | grep -q -E -- ' -j1( |$$)'; then \\
\t useLocks=no; \\
\telif echo "$(MAKEFLAGS)" | grep -q -E -- ' -j( |$$)'; then \\
\t useLocks=$(LOCKMD5); \\
\telif echo "$(MAKEFLAGS)" | grep -q -E -- ' -j[0-9]+( |$$)'; then \\
\t useLocks=$(LOCKMD5); \\
\telse \\
\t useLocks=no; \\
\tfi; \\
\t./scripts/build/sourceDigests.pl `pwd` {file_name_root}.exe $$useLocks
\t$(CCOMPILER) -c {work_dir}{file_name_root}.md5s.c -o {work_dir}{file_name_root}.md5s.o $(CFLAGS)
\t$(FCCOMPILER) `cat {work_dir}{file_name_root}.d` {work_dir}{file_name_root}.parameters.o {work_dir}{file_name_root}.md5s.o -o {file_name_root}.exe$(SUFFIX) $(FCFLAGS) `./scripts/build/libraryDependencies.py {file_name_root}.exe $(FCFLAGS)` 2>&1 | ./scripts/build/postprocessLinker.py

"""
        out.write(rule)

    if executable_names:
        out.write("all_exes = " + " ".join(executable_names) + "\n")
