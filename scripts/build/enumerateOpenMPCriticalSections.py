#!/usr/bin/env python3
import os
import re
import sys
import xml.etree.ElementTree as ET

from Galacticus.Build.FortranUtils import get_fortran_line

# Locate all OpenMP critical sections, and build an enumeration of them for use
# in source code instrumentation.
# Andrew Benson (ported to Python 2026)

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

for full_path in sorted(source_file_paths):
    try:
        with open(full_path, 'r', errors='replace') as fh:
            while True:
                raw, processed, _ = get_fortran_line(fh)
                if not raw:
                    break
                m = re.match(r'^\s*!\$omp\s+critical\s*\(([a-z0-9_]+)\)', processed, re.IGNORECASE)
                if m:
                    name = m.group(1).lower()
                    critical_section_names[name] = critical_section_names.get(name, 0) + 1
    except OSError:
        continue


def _update_file(old_path, new_path):
    """Replace old_path with new_path only if their contents differ."""
    import shutil
    if not os.path.exists(old_path):
        shutil.move(new_path, old_path)
    else:
        with open(old_path, 'rb') as f1, open(new_path, 'rb') as f2:
            if f1.read() == f2.read():
                os.unlink(new_path)
            else:
                shutil.move(new_path, old_path)


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
_update_file(os.path.join(build_path, "openMPCriticalSections.xml"), xml_tmp)

# --- openMPCriticalSections.count.inc ---
count_tmp = os.path.join(build_path, "openMPCriticalSections.count.inc.tmp")
with open(count_tmp, 'w') as fh:
    fh.write("! Number of named OpenMP critical sections in the source.\n")
    fh.write(f"integer, public, parameter :: criticalSectionCount={len(critical_section_names)}\n")
_update_file(os.path.join(build_path, "openMPCriticalSections.count.inc"), count_tmp)

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
_update_file(os.path.join(build_path, "openMPCriticalSections.enumerate.inc"), enum_tmp)
