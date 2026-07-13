#!/usr/bin/env python3
import os
import re
import sys

from Galacticus.Build.ParallelScan import scan as parallel_scan

# Locate files which contain programs and append to a list of executables.
# Andrew Benson (ported to Python 2026)

if len(sys.argv) != 2:
    print("Usage: findExecutables.py <galacticusDirectory>", file=sys.stderr)
    sys.exit(1)

galacticus_directory = sys.argv[1]
source_directory     = os.path.join(galacticus_directory, "source")
build_path           = os.environ['BUILDPATH']
work_dir             = build_path.rstrip('/') + '/'


def _source_files_recursive(root):
    """Yield every file path under `root`, relative to `root`, walking
    subdirectories to any depth."""
    rel = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith('.'))
        for fn in filenames:
            rel.append(os.path.relpath(os.path.join(dirpath, fn), root))
    return sorted(rel)


# ---------------------------------------------------------------------------
# Per-file processing (parallelised across files)
# ---------------------------------------------------------------------------
#
# Reading each candidate source file (to detect a `program` statement) is the
# slow, serial part on NFS: each open blocks on close-to-open metadata
# revalidation, so doing them one at a time serialises thousands of millisecond
# waits. Each file is scanned independently and produces only its own result
# (the Makefile rule text plus the executable name it contributes, if any),
# writing to no shared state. The read-only `source_directory`/`work_dir` are
# published as module-level globals BEFORE the scan so forked workers inherit
# them via copy-on-write rather than re-pickling per task. Results come back in
# task order, so the emitted Makefile is byte-identical to the serial version.


def _scan_one(file_name):
    """Worker: read+parse ONE source file. Returns `None` if the file
    contributes nothing, otherwise `(rule_text, exe_name_or_None)` where
    `exe_name_or_None` is the executable name to add to `all_exes` (or `None`
    when the file is excluded from `all`). Writes to no shared state.
    """
    file_full        = os.path.join(source_directory, file_name)
    exclude_from_all = False
    found_program    = False

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
        return None

    if not found_program:
        return None

    # `obj_root` is the path-based stem (relative to source/, e.g.
    # `tests/nodes`) used for every build artifact, which mirrors the
    # source hierarchy.  `exe_root` is the historical flat, dot-separated
    # name (e.g. `tests.nodes`) used for the user-facing executable so that
    # `make tests.nodes.exe`, CI matrices, and test harnesses keep working
    # unchanged after the source tree was made hierarchical.
    obj_root = re.sub(r'\.[fF](90)?$', '', file_name)
    exe_root = obj_root.replace('/', '.')

    exe_name = None if exclude_from_all else exe_root + '.exe'

    # The link procedure itself lives in the main Makefile's LINK_EXECUTABLE canned recipe (which
    # has a single owner there); each generated rule contributes only the target's dependency line
    # and a $(call) to that recipe.
    rule = f"""\
{exe_root}.exe: {work_dir}{obj_root}.o {work_dir}{obj_root}.d $(MAKE_DEPS)
\t$(call LINK_EXECUTABLE,{exe_root},{obj_root})

"""
    return rule, exe_name


# Build the ordered list of candidate files, applying the cheap filters that
# never required reading the file (so they impose no I/O cost), then scan the
# survivors concurrently.
candidate_files = []
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
    candidate_files.append(file_name)

results = parallel_scan(candidate_files, _scan_one, "findExecutables.py")

# Write the Makefile SEQUENTIALLY, consuming results in order so the output is
# byte-identical to the serial version.
executable_names = []
with open(os.path.join(build_path, "Makefile_All_Execs"), 'w') as out:
    for result in results:
        if result is None:
            continue
        rule, exe_name = result
        if exe_name is not None:
            executable_names.append(exe_name)
        out.write(rule)

    if executable_names:
        out.write("all_exes = " + " ".join(executable_names) + "\n")
