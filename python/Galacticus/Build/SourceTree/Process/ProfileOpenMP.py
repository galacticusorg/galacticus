# If the `-DOMPPROFILE` flag is in $GALACTICUS_FCFLAGS, wraps every
# `!$omp critical(name)` directive listed in `$BUILDPATH/openMPCriticalSections.xml`
# with timing calls that accumulate wait-time per critical section.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/ProfileOpenMP.pm

import re
import os
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from XML.Utils                                    import xml_to_dict
from Galacticus.Build.SourceTree                  import walk_tree
from Galacticus.Build.SourceTree.Process          import register_process
from Galacticus.Build.SourceTree.Parse.ModuleUses import add_uses


def _profile_enabled():
    return '-DOMPPROFILE' in os.environ.get('GALACTICUS_FCFLAGS', '').split()


def _load_critical_sections():
    """Load `$BUILDPATH/openMPCriticalSections.xml` and return
    `{name: {'id': <int>}}`.
    """
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        raise RuntimeError(
            "process_profile_openmp: BUILDPATH environment variable is not set")
    path = os.path.join(build_path, 'openMPCriticalSections.xml')
    if not os.path.exists(path):
        raise RuntimeError(
            "process_profile_openmp: critical section enumeration file does not exist: "
            + path)
    root = ET.parse(path).getroot()
    # Perl's XML::Simple groups `<critical name=…/>` children by name; do the
    # same here so `descriptor['critical']['my_section']` works as in Perl.
    descriptor = xml_to_dict(root, keyed_tags={'critical'})
    return descriptor.get('critical', {}) or {}


_OMP_CRITICAL_RE = re.compile(r'^\s*!\$omp\s+critical\s*\(([a-z0-9_]+)\)', re.IGNORECASE)


def _enclosing_subprogram(node):
    cur = node
    while cur is not None and cur.get('type') not in ('function', 'subroutine'):
        if cur.get('parent') is None:
            return None
        cur = cur.get('parent')
    return cur


def _enclosing_module_name(node):
    cur = node
    while cur is not None:
        if cur.get('type') == 'module':
            return cur.get('name', '')
        if cur.get('parent') is None:
            return ''
        cur = cur.get('parent')
    return ''


def process_profile_openmp(tree, options):
    """Mirrors Profile_OpenMP() from ProfileOpenMP.pm."""
    if not _profile_enabled():
        return
    descriptor = _load_critical_sections()

    module_uses = {
        'moduleUse': {
            'OpenMP_Utilities_Data': {'intrinsic': False, 'all': True},
            'OMP_Lib':               {'intrinsic': False, 'all': True},
        },
        'moduleOrder': ['OpenMP_Utilities_Data', 'OMP_Lib'],
    }

    for node in walk_tree(tree):
        if node.get('type') != 'code':
            continue
        content = node.get('content', '')
        new_lines = []
        for line in content.splitlines(keepends=True):
            m = _OMP_CRITICAL_RE.match(line)
            if not m:
                new_lines.append(line)
                continue
            name = m.group(1).lower()
            if name not in descriptor:
                new_lines.append(line)
                continue
            subprogram = _enclosing_subprogram(node)
            if subprogram is None or subprogram.get('type') not in ('function', 'subroutine'):
                new_lines.append(line)
                continue
            # Modules used by the profiling support must not pull those very
            # modules back in — skip the substitution in those cases.
            if _enclosing_module_name(node) == 'iso_varying_string':
                new_lines.append(line)
                continue
            add_uses(subprogram, module_uses)
            section_id = descriptor[name].get('id', 0)
            new_lines.append(
                "ompProfileTimeWaitStart=OMP_Get_WTime()\n"
                + line
                + "ompProfileTimeWaitEnd=OMP_Get_WTime()\n"
                "ompProfileTimeWaitEnd=ompProfileTimeWaitEnd-ompProfileTimeWaitStart\n"
                + f"criticalSectionWaitTime({section_id})=criticalSectionWaitTime({section_id})+ompProfileTimeWaitEnd\n"
            )
        node['content'] = ''.join(new_lines)


register_process('profileOpenMP', process_profile_openmp)
