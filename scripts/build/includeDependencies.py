#!/usr/bin/env python3
import os
import re
import subprocess
import sys

# Construct Makefile rules for source files which have dependencies on include files.
# Andrew Benson (ported to Python 2026)

if len(sys.argv) != 2:
    print("Usage: includeDependencies.py <galacticusDirectory>", file=sys.stderr)
    sys.exit(1)

galacticus_directory = sys.argv[1]
source_directory     = os.path.join(galacticus_directory, "source")
build_path           = os.environ['BUILDPATH']

# Build a list of pre-existing C include file paths: source dir, system GCC paths,
# and any paths from GALACTICUS_CFLAGS / GALACTICUS_CPPFLAGS.
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

# Open output Makefile.
makefile_path = os.path.join(build_path, "Makefile_Include_Dependencies")
dependency_file_names = []

with open(makefile_path, 'w') as makefile:
    for file_name in sorted(os.listdir(source_directory)):
        # Skip temporary files.
        if file_name.startswith('.#'):
            continue
        # Only process Fortran, Fortran include, C, C++, and header files.
        if not re.search(r'\.(f(90)?|inc|c(pp)?|h)$', file_name, re.IGNORECASE):
            continue

        file_full = os.path.join(source_directory, file_name)
        included_files = []

        try:
            with open(file_full, 'r', errors='replace') as fh:
                for line in fh:
                    m = re.match(
                        r'''^\s*#??include\s*[<'"]((.+)\.(inc|h)\d*)['">](\s*!\s*NO_USES)?''',
                        line
                    )
                    if not m:
                        continue
                    inc_name  = m.group(1)
                    ext       = m.group(3)
                    no_uses   = m.group(4) is not None
                    if ext == 'inc':
                        included_files.append({'fileName': inc_name, 'automatic': not no_uses})
                    elif ext == 'h':
                        # Only include .h files not already present in source or system dirs.
                        in_source = os.path.exists(os.path.join(source_directory, inc_name))
                        in_system = any(
                            os.path.exists(os.path.join(d, inc_name)) for d in c_include_dirs
                        )
                        if not in_source and not in_system:
                            included_files.append({'fileName': inc_name, 'automatic': not no_uses})
        except OSError:
            continue

        if included_files:
            obj_name = re.sub(r'\.[^.]+$', '.o', file_name)
            deps     = ' '.join(sorted(
                os.path.join(build_path, f['fileName']) for f in included_files
            ))
            makefile.write(f"{os.path.join(build_path, obj_name)}: {deps}\n\n")
            if file_name.endswith('.F90'):
                pp_name = file_name.replace('.F90', '.p.F90')
                makefile.write(f"{os.path.join(build_path, pp_name)}: {deps}\n\n")

        # Collect auto-generated .inc dependencies for Makefile_Use_Dependencies.
        for inc in included_files:
            if inc['fileName'].endswith('.inc') and inc['automatic']:
                unpp = re.sub(r'\.inc$', '.Inc', inc['fileName'])
                if os.path.exists(os.path.join(source_directory, unpp)):
                    dependency_file_names.append(unpp)
                else:
                    dependency_file_names.append(inc['fileName'])

    # Dependency line for Makefile_Use_Dependencies.
    dep_targets = ' '.join(
        os.path.join(build_path, f) for f in sorted(set(dependency_file_names))
    )
    makefile.write(
        f"\n{os.path.join(build_path, 'Makefile_Use_Dependencies')}: {dep_targets}\n"
    )
