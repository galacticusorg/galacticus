#!/usr/bin/env python3
import os
import re
import shutil
import subprocess
import sys
from graphlib import TopologicalSorter

# Output linker options to link required libraries for building an executable.
# Andrew Benson (ported to Python 2026)

if len(sys.argv) < 2:
    print("Usage: libraryDependencies.py <executable> [compiler options...]", file=sys.stderr)
    sys.exit(1)

executable       = sys.argv[1]
compiler_options = sys.argv[2:]

# Library dependency graph: key depends on values.
dependencies = {
    'hdf5hl_fortran': ['hdf5_hl'],
    'hdf5_hl'       : ['hdf5'],
    'hdf5_fortran'  : ['hdf5'],
    'hdf5'          : ['z'],
    'gsl'           : ['gslcblas'],
    'FoX_dom'       : ['FoX_fsys', 'FoX_utils', 'FoX_sax'],
    'FoX_sax'       : ['FoX_common'],
    'FoX_utils'     : ['FoX_wxml'],
    'qhullcpp'      : ['qhull_r', 'stdc++'],
}

# Static-link ordering: key must appear before values.
static_link_dependencies = {
    'hdf5'          : ['z', 'dl'],
    'hdf5_hl'       : ['hdf5'],
    'hdf5_fortran'  : ['hdf5'],
    'hdf5hl_fortran': ['hdf5_hl'],
    'gsl'           : ['gslcblas'],
    'FoX_dom'       : ['FoX_fsys', 'FoX_utils', 'FoX_sax', 'FoX_wxml'],
    'FoX_sax'       : ['FoX_common'],
    'FoX_wxml'      : ['FoX_utils'],
    'FoX_common'    : ['FoX_fsys'],
    'qhullcpp'      : ['stdc++'],
}

# Detect compiler preprocessor defines.
c_compiler = os.environ.get('CCOMPILER', 'gcc')
try:
    defs_out = subprocess.run(
        [c_compiler, '-dM', '-E', '-'],
        input='', capture_output=True, text=True
    ).stdout
except FileNotFoundError:
    defs_out = ''
preprocessor_directives = [
    line.split()[1] for line in defs_out.splitlines()
    if len(line.split()) >= 2 and line.startswith('#define')
]

is_macos  = '__APPLE__' in preprocessor_directives
is_static = '-static'   in compiler_options

# For static builds, HDF5 needs an explicit libdl dependency.
if is_static:
    dependencies.setdefault('hdf5', []).append('dl')

pthread_included = '-lpthread' in compiler_options

build_path = os.environ.get('BUILDPATH')
if not build_path:
    print('libraryDependencies.py: "BUILDPATH" environment variable is not set', file=sys.stderr)
    sys.exit(1)

# Locate the main dependency file (.d) for the executable.
main_dep = re.sub(r'\.(exe|o)$', '.d', os.path.join(build_path, executable))
if not os.path.exists(main_dep):
    print(f'libraryDependencies.py: dependency file is missing: {main_dep}', file=sys.stderr)
    sys.exit(1)

dep_files = [main_dep]
if executable == 'libgalacticus.o':
    dep_files.append(os.path.join(build_path, 'libgalacticus_classes.d'))

object_files = []
for dep_file in dep_files:
    with open(dep_file) as fh:
        object_files.extend(line.rstrip('\n') for line in fh)

# Read .fl (library dependency) files for each object file.
libraries = {}
for obj in object_files:
    fl_file = re.sub(r'\.o$', '.fl', obj)
    if os.path.exists(fl_file):
        with open(fl_file) as fh:
            for lib in fh:
                lib = lib.rstrip('\n')
                if lib:
                    libraries[lib] = libraries.get(lib, 0) + 1

# Transitively expand library dependencies.
lib_count = -1
while len(libraries) != lib_count:
    lib_count = len(libraries)
    for lib in sorted(libraries):
        for dep in dependencies.get(lib, []):
            libraries[dep] = libraries.get(dep, 0) + 1

# Remove conditionally compiled libraries if the corresponding flag is absent.
conditional = {
    'fftw3'   : '-DFFTW3AVAIL',
    'ANN'     : '-DANNAVAIL',
    'qhullcpp': '-DQHULLAVAIL',
    'qhull_r' : '-DQHULLAVAIL',
    'matheval': '-DMATHEVALAVAIL',
    'git2'    : '-DGIT2AVAIL',
}
for lib, flag in conditional.items():
    if lib in libraries and flag not in compiler_options:
        del libraries[lib]

# Topological sort using graphlib (Python 3.9+).
# Build a graph where an edge A→B means A must come before B (A depends on B in link order).
# We invert static_link_dependencies: if key must appear before values, then values depend on key.
ts_graph = {lib: set() for lib in libraries}
for lib, before_these in static_link_dependencies.items():
    if lib not in libraries:
        continue
    for dep in before_these:
        if dep in libraries:
            # dep must come after lib, so dep has a dependency on lib.
            ts_graph.setdefault(dep, set()).add(lib)

ts = TopologicalSorter(ts_graph)
try:
    sorted_libraries = list(ts.static_order())
except Exception:
    sorted_libraries = sorted(libraries)
# Filter to only the libraries we actually need.
sorted_libraries = [lib for lib in sorted_libraries if lib in libraries]

# Rename qhull_r for static/macOS builds.
if is_static or is_macos:
    sorted_libraries = [
        'qhullstatic_r' if lib == 'qhull_r' else lib
        for lib in sorted_libraries
    ]

# Build static-link option.
static_opts = ' -Wl,--whole-archive -lpthread -Wl,--no-whole-archive' \
    if (is_static and not pthread_included) else ''

# macOS-specific option.
import platform
os_opts = ' -Wl,-commons,use_dylibs' if platform.system() == 'Darwin' else ''

# Determine whether to link librt (needed for older glibc).
link_librt = True
ldd = shutil.which('ldd')
if ldd:
    try:
        ldd_out = subprocess.run(['ldd', '--version'], capture_output=True, text=True).stdout
        for line in ldd_out.splitlines():
            m = re.match(r'^ldd \(GNU libc\) (\d+)\.(\d+)', line)
            if m:
                major, minor = int(m.group(1)), int(m.group(2))
                link_librt = major < 2 or (major == 2 and minor <= 16)
    except Exception:
        pass
else:
    link_librt = False
    
# If ldd not found, link_librt stays False.

lib_flags    = ' '.join(f'-l{lib}' for lib in sorted_libraries)
debug_flag   = ' -rdynamic'       if '-DDEBUGGING' in compiler_options else ''
blas_flag    = ' -fexternal-blas' if 'blas' in sorted_libraries        else ''
rt_flag      = ' -lrt'            if link_librt                         else ''

print(lib_flags + static_opts + os_opts + debug_flag + blas_flag + rt_flag)
