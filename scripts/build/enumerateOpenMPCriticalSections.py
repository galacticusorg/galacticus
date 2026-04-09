#!/usr/bin/env python3
import os
import re
import sys
import xml.etree.ElementTree as ET

galacticus_exec_path = os.environ.get('GALACTICUS_EXEC_PATH')
if not galacticus_exec_path:
    print("Error: GALACTICUS_EXEC_PATH environment variable is not set.", file=sys.stderr)
    sys.exit(1)
python_path = os.path.abspath(os.path.join(galacticus_exec_path, 'python'))
sys.path.insert(0, python_path)
from build.fortran_utils import get_fortran_line

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

for file_name in os.listdir(src_path):
    if not re.search(r'\.f(90)?$', file_name, re.IGNORECASE):
        continue
    full_path = os.path.join(src_path, file_name)
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
enum_tmp = os.path.join(build_path, "openMPCriticalSections.enumerate.inc.tmp")
sorted_names = sorted(critical_section_names)
with open(enum_tmp, 'w') as fh:
    fh.write("type(varying_string), dimension(criticalSectionCount) :: criticalSectionNames\n")
    fh.write("criticalSectionNames=[ &\n")
    joined = "'), &\n & var_str('".join(sorted_names)
    fh.write(f" & var_str('{joined}') &\n")
    fh.write(" & ]\n")
_update_file(os.path.join(build_path, "openMPCriticalSections.enumerate.inc"), enum_tmp)
