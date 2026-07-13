#!/usr/bin/env python3
import os
import re
import sys
import xml.etree.ElementTree as ET

from Galacticus.Build.FileChanges  import update as file_changes_update
from Galacticus.Build.FortranUtils import get_fortran_line
from Galacticus.Build.ParallelScan import scan as parallel_scan

# Locate all OpenMP critical sections, and build an enumeration of them for use
# in source code instrumentation.
# Andrew Benson (ported to Python 2026)


def _scan_one(full_path):
    """Worker: scan one source file and return a dict mapping each OpenMP
    critical-section name found to its occurrence count within that file.

    Each file is read independently and writes to no shared state, so the scans
    run concurrently across a fork-based process pool; the per-file cost is
    dominated by blocking file I/O (open/read over NFS), which overlaps when run
    in parallel. The small per-file count dicts are then merged in task order.
    """
    counts = {}
    try:
        with open(full_path, 'r', errors='replace') as fh:
            while True:
                raw, processed, _ = get_fortran_line(fh)
                if not raw:
                    break
                m = re.match(r'^\s*!\$omp\s+critical\s*\(([a-z0-9_]+)\)', processed, re.IGNORECASE)
                if m:
                    name = m.group(1).lower()
                    counts[name] = counts.get(name, 0) + 1
    except OSError:
        return {}
    return counts


if len(sys.argv) != 2:
    print("Usage: enumerateOpenMPCriticalSections.py <sourceDirectory>", file=sys.stderr)
    sys.exit(1)

source_directory = sys.argv[1]
build_path       = os.environ['BUILDPATH']
src_path         = os.path.join(source_directory, "source")

critical_section_names = {}

source_file_paths = []
for dirpath, dirnames, filenames in os.walk(src_path):
    dirnames[:] = sorted(d for d in dirnames if not d.startswith('.'))
    for file_name in filenames:
        if re.search(r'\.f(90)?$', file_name, re.IGNORECASE):
            source_file_paths.append(os.path.join(dirpath, file_name))

# Scan every file (concurrently) and merge results in the original sorted order
# so the accumulated counts match a serial run exactly.
tasks = sorted(source_file_paths)
for counts in parallel_scan(tasks, _scan_one, "enumerateOpenMPCriticalSections.py"):
    for name, count in counts.items():
        critical_section_names[name] = critical_section_names.get(name, 0) + count


# --- openMPCriticalSections.xml ---
xml_tmp = os.path.join(build_path, "openMPCriticalSections.xml.tmp")
root = ET.Element("criticalSections")
for idx, name in enumerate(sorted(critical_section_names), start=1):
    child = ET.SubElement(root, "critical")
    child.set("name",      name)
    child.set("frequency", str(critical_section_names[name]))
    child.set("id",        str(idx))
tree = ET.ElementTree(root)
ET.indent(tree, space='  ')
tree.write(xml_tmp, encoding='unicode', xml_declaration=False)
with open(xml_tmp, 'a') as fh:
    fh.write('\n')
# `prove_update=True` touches the `openMPCriticalSections.xml.up` sentinel — the Makefile's rule
# target — so make records that the enumeration ran even when the .xml itself (whose mtime drives
# the re-preprocessing cascade) was left untouched because its content did not change.
file_changes_update(os.path.join(build_path, "openMPCriticalSections.xml"), xml_tmp,
                    prove_update=True)

# --- openMPCriticalSections.count.inc ---
count_tmp = os.path.join(build_path, "openMPCriticalSections.count.inc.tmp")
with open(count_tmp, 'w') as fh:
    fh.write("! Number of named OpenMP critical sections in the source.\n")
    fh.write(f"integer, public, parameter :: criticalSectionCount={len(critical_section_names)}\n")
file_changes_update(os.path.join(build_path, "openMPCriticalSections.count.inc"), count_tmp)

# --- openMPCriticalSections.enumerate.inc ---
# A fixed-length `character` array is used (rather than `varying_string`) to
# mirror the event-hook wait-time code in EventHooks.py: a local automatic
# `varying_string` array segfaults gfortran when it is finalized on return
# from OpenMP_Critical_Wait_Times (it is written to the output file and never
# otherwise needs the dynamic length).
enum_tmp = os.path.join(build_path, "openMPCriticalSections.enumerate.inc.tmp")
sorted_names = sorted(critical_section_names)
name_length_max = max((len(name) for name in sorted_names), default=1)
with open(enum_tmp, 'w') as fh:
    fh.write(f"character(len={name_length_max}), dimension(criticalSectionCount) :: criticalSectionNames\n")
    fh.write(f"criticalSectionNames=[character(len={name_length_max}) :: &\n")
    joined = "', &\n & '".join(sorted_names)
    fh.write(f" & '{joined}' &\n")
    fh.write(" & ]\n")
file_changes_update(os.path.join(build_path, "openMPCriticalSections.enumerate.inc"), enum_tmp)
