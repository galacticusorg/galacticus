#!/usr/bin/env python3
"""Generate `Makefile_Module_Dependencies`: per-module and per-submodule
Makefile rules that tell the Galacticus build system how to produce `.mod`
and `.smod` files from their owning object files.

For every `.f` / `.f90` source file under `<installDir>/source/` (one level
deep), this script:

* Extracts every `module <name>` and `submodule (<parent>[:<grand>]) <name>`
  declaration (skipping text inside `!!{…!!}` LaTeX and `!![…!!]` XML blocks).
* For every `functionClass` directive in the file, walks each instance file
  listed in `directiveLocations.xml` to discover the concrete derived types
  and synthesises one functionClass-submodule entry per derived type.
* Follows `include '<leaf>'` statements that reference files under the
  top-level `source/` directory, scanning them for the same constructs.
* Emits Makefile rules for each discovered module / submodule target: the
  `<mod>.mod`, `<mod>@<sub>.smod`, `<src>.o`, `<src>.p.F90`, `<mod>.mod.d`,
  `<mod>.mod.gv`, `<src>.m`, `<src>.smod` dependencies the main build
  graph expects.

Results are cached in `Makefile_Module_Dependencies.blob` (pickle) so
subsequent runs only rescan the files whose mtimes have advanced past the
cache's own mtime.

Mirrors scripts/build/moduleDependencies.pl.
Andrew Benson (ported to Python 2026).
"""

import os
import pickle
import re
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from Galacticus.Build.Directives import extract_directives
from Galacticus.Build.SourceTree import parse_file, walk_tree
from List.ExtraUtils             import as_array
from XML.Utils                   import xml_to_dict


# ---------------------------------------------------------------------------
# Line-level regexes
# ---------------------------------------------------------------------------

_XML_OPEN_RE     = re.compile(r'^\s*!!\[')
_XML_CLOSE_RE    = re.compile(r'^\s*!!\]')
_LATEX_OPEN_RE   = re.compile(r'^\s*!!\{')
_LATEX_CLOSE_RE  = re.compile(r'^\s*!!\}')

_STRIP_COMMENT_RE = re.compile(r'\s*(!.*)?$')

_MODULE_LINE_RE    = re.compile(r'^\s*module\s+([a-zA-Z0-9_]+)$', re.IGNORECASE)
_SUBMODULE_LINE_RE = re.compile(
    r'^\s*submodule\s*\(\s*([a-z0-9_]+)(?::([a-z0-9_]+))?\s*\)'
    r'\s*([a-zA-Z0-9_]+)$',
    re.IGNORECASE,
)
_INCLUDE_LINE_RE   = re.compile(r"""include\s+'(\w+)'""", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Cache / helpers
# ---------------------------------------------------------------------------

def _file_identifier(path):
    """Perl `(my $id = $path) =~ s/\\//_/g; $id =~ s/^\\._??//;`."""
    return re.sub(r'^\._?', '', path.replace('/', '_'))


def _load_cache(blob_path):
    if not os.path.exists(blob_path):
        return {}, None
    try:
        with open(blob_path, 'rb') as fh:
            cache = pickle.load(fh)
    except (pickle.UnpicklingError, EOFError, AttributeError, ValueError,
            ImportError, ModuleNotFoundError):
        return {}, None
    if not isinstance(cache, dict):
        return {}, None
    return cache, os.stat(blob_path).st_mtime


def _load_directive_locations(path):
    if not os.path.exists(path):
        return {}
    return xml_to_dict(ET.parse(path).getroot())


def _source_directory_descriptors(source_root, build_path):
    """Return the list of `{path, leaf}` descriptors -- the top-level `source/`
    directory with leaf "", plus every immediate subdirectory with
    `leaf="<subdir>/"`.  Also `mkdir -p` the matching build subdirectories so
    emitted object paths can land there.
    """
    descriptors = [{'path': source_root, 'leaf': ''}]
    if not os.path.exists(source_root):
        return descriptors
    for entry in sorted(os.listdir(source_root)):
        if entry in ('.', '..'):
            continue
        full = os.path.join(source_root, entry)
        if os.path.isdir(full):
            descriptors.append({'path': full, 'leaf': entry + '/'})
    for d in descriptors:
        os.makedirs(os.path.join(build_path, d['leaf']), exist_ok=True)
    return descriptors


def _list_source_files(directory):
    try:
        names = os.listdir(directory)
    except OSError:
        sys.exit("moduleDependencies.py: can not open the source directory: "
                 + directory)
    return sorted(
        n for n in names
        if re.search(r'\.(f|f90)$', n, re.IGNORECASE) and not n.startswith('.#')
    )


# ---------------------------------------------------------------------------
# functionClass submodule discovery
# ---------------------------------------------------------------------------

def _find_function_class_submodules(directive_name, function_class_file,
                                    locations):
    """Return the list of submodule records produced by one `functionClass`
    directive.  Mirrors moduleDependencies.pl:97-133.
    """
    # Regex is dynamic: the class name is substituted in.
    type_re = re.compile(
        r'^\s*type\s*'
        r'(?:,\s*(?:abstract|public|private|extends\s*\(([a-zA-Z0-9_]+)\))\s*)*'
        r'(?:::)?\s*' + re.escape(directive_name) + r'([a-z0-9_]+)\s*$',
        re.IGNORECASE,
    )

    tree = parse_file(function_class_file)
    classes = {}
    for node in walk_tree(tree):
        ntype = node.get('type')
        if ntype == directive_name:
            name = (node.get('directive') or {}).get('name')
            if name:
                classes.setdefault(name, {})['name'] = name
        elif ntype == 'type':
            m = type_re.match(node.get('opener') or '')
            if m and m.group(1) is not None:
                class_name = directive_name + m.group(2)
                classes.setdefault(class_name, {})['extends'] = m.group(1)

    submodules = []
    leaf = re.sub(r'^.*/(.*)\.F90$', r'\1', function_class_file)
    for class_name, info in classes.items():
        own_name   = info.get('name') or class_name
        extends    = info.get('extends')
        submodule_name = own_name + '_'
        if extends == directive_name + 'Class':
            parent = None
        elif extends is not None:
            parent = extends + '_'
        else:
            parent = None
        submodules.append({
            'name':          submodule_name,
            'fileName':      leaf,
            'extends':       parent,
            'source':        function_class_file,
            'functionClass': True,
        })
    return submodules


# ---------------------------------------------------------------------------
# Per-file Fortran scan
# ---------------------------------------------------------------------------

def _scan_file(stack, entry, source_root):
    """Drive the include stack rooted at `stack`, populating `entry` (which
    may be a plain dict; expected keys populated: `files`, `modulesProvided`,
    `submodulesProvided`).
    """
    while stack:
        file_path = stack.pop()
        in_xml   = False
        in_latex = False
        try:
            fh = open(file_path, 'r', errors='replace')
        except OSError:
            sys.exit("moduleDependencies.py: can not open input file: "
                     + file_path)
        try:
            for line in fh:
                # Leave XML/LaTeX blocks on their closing markers.
                if _XML_CLOSE_RE.match(line):
                    in_xml = False
                if _LATEX_CLOSE_RE.match(line):
                    in_latex = False
                if in_xml or in_latex:
                    continue

                stripped = _STRIP_COMMENT_RE.sub('', line)

                m = _MODULE_LINE_RE.match(stripped)
                if m:
                    entry.setdefault('modulesProvided', []).append(
                        m.group(1).lower() + '.mod',
                    )

                m = _SUBMODULE_LINE_RE.match(stripped)
                if m:
                    parent_mod = m.group(1).lower() + '.mod'
                    submodule_name = m.group(3).lower()
                    leaf = re.sub(r'^.*/([^/]+)\.F90$', r'\1', file_path)
                    entry.setdefault('submodulesProvided', []).append({
                        'moduleName': parent_mod,
                        'submodule':  {
                            'name':          submodule_name,
                            'fileName':      leaf,
                            'source':        file_path,
                            'extends':       None,
                            'functionClass': False,
                        },
                    })

                m = _INCLUDE_LINE_RE.search(stripped)
                if m:
                    include_path = os.path.join(source_root, m.group(1))
                    stack.append(include_path)
                    entry.setdefault('files', []).append(include_path)

                # Enter XML/LaTeX blocks on their opening markers (check
                # AFTER processing the line so `!![` itself is not mistaken
                # for code).
                if _XML_OPEN_RE.match(line):
                    in_xml = True
                if _LATEX_OPEN_RE.match(line):
                    in_latex = True
        finally:
            fh.close()


# ---------------------------------------------------------------------------
# Submodule-by-module map
# ---------------------------------------------------------------------------

def _build_submodule_map(modules_per_file):
    """Return a dict mapping lowercase `<module>.mod` to the ordered list of
    submodule records associated with it.  Combines both functionClass-
    synthesised submodules and Fortran `submodule (…)` statements.
    """
    submodules = {}
    for file_id, entry in modules_per_file.items():
        if not entry.get('submodules'):
            continue
        provided = entry.get('modulesProvided') or []
        if len(provided) != 1:
            raise RuntimeError(
                "moduleDependencies.py: submodules associated with multiple "
                "modules"
            )
        module_file_name = provided[0]
        module_name = re.sub(r'\.mod$', '', module_file_name)
        submodules[module_name.lower() + '.mod'] = list(entry['submodules'])

    for file_id, entry in modules_per_file.items():
        for rec in entry.get('submodulesProvided') or []:
            module_file_name = rec['moduleName']
            submodules.setdefault(module_file_name, []).append(rec['submodule'])

    return submodules


# ---------------------------------------------------------------------------
# Makefile emission
# ---------------------------------------------------------------------------

def _write_makefile(path, modules_per_file, submodules_by_module, work_dir):
    """Write `Makefile_Module_Dependencies`.  Mirrors the output loop at
    moduleDependencies.pl:207-300.
    """
    with open(path, 'w') as mk:
        for file_identifier in sorted(modules_per_file):
            entry  = modules_per_file[file_identifier]
            mods   = entry.get('modulesProvided') or []
            if not mods:
                continue
            source_file_name = entry['sourceFileName']
            leaf             = entry['sourceDirectoryDescriptor']['leaf']
            object_file_name = re.sub(
                r'\.(f|f90)$', '.o', source_file_name,
                flags=re.IGNORECASE,
            )
            obj_dep_file     = re.sub(r'\.o$', '.d', object_file_name)

            for module_file_name in mods:
                mod_key   = module_file_name.lower()
                submodule_list = submodules_by_module.get(mod_key, [])
                mod_path  = work_dir + module_file_name
                leaf_obj  = work_dir + leaf + object_file_name

                # <build>/foo.mod: <build>/<leaf>/foo.o   (+conditional rebuild)
                mk.write(f"{mod_path}: {leaf_obj}\n")
                mk.write(f"\t@if [ ! -f {mod_path} ]; then \\\n")
                mk.write(f"\t  rm {leaf_obj} ; \\\n")
                mk.write(f"\t  $(MAKE) {leaf_obj} ; \\\n")
                mk.write("\tfi\n\n")

                module_name = re.sub(r'\.mod$', '', module_file_name)

                for submodule in submodule_list:
                    sub_name  = submodule['name'].lower()
                    sub_file  = submodule['fileName']
                    smod_path = (
                        f"{work_dir}{module_name}@{sub_name}.smod"
                    )
                    sub_obj = f"{work_dir}{sub_file}.o"

                    mk.write(f"{smod_path}: {sub_obj}\n")
                    mk.write(f"\t@if [ ! -f {smod_path} ]; then \\\n")
                    mk.write(f"\t  rm {sub_obj} ; \\\n")
                    mk.write(f"\t  $(MAKE) {sub_obj} ; \\\n")
                    mk.write("\tfi\n\n")

                    if submodule['extends'] is not None:
                        matches = [
                            s for s in submodule_list
                            if s['name'] == submodule['extends']
                        ]
                        if not matches:
                            print("no matching submodule found:")
                            print(
                                f"\t'{submodule['name']}' extends "
                                f"'{submodule['extends']}'"
                            )
                            print("\tavailable submodules are:")
                            for avail in submodule_list:
                                print(f"\t\t'{avail['name']}'")
                            sys.exit(
                                "ERROR: moduleDependencies.py: "
                                "no matching submodule found"
                            )
                        if len(matches) > 1:
                            print(
                                "ERROR: class '" + submodule['name']
                                + "' is an extension of class '"
                                + submodule['extends']
                                + "' which is multiply defined in files:"
                            )
                            for m in matches:
                                print("\t" + m['source'])
                            sys.exit(
                                "moduleDependencies.py: multiple matching "
                                "submodules found"
                            )
                        depends_on = matches[0]['fileName'] + '.o'
                    else:
                        depends_on = leaf + object_file_name

                    mk.write(f"{sub_obj}: {work_dir}{depends_on}\n\n")

                    if submodule['functionClass']:
                        preprocessed = re.sub(
                            r'\.F90$', '.p.F90', source_file_name,
                        )
                        p_file = (
                            f"{work_dir}{sub_file}.p.F90"
                        )
                        p_src  = f"{work_dir}{leaf}{preprocessed}"
                        mk.write(
                            f"{p_file}.up: {p_src} {submodule['source']}\n"
                        )
                        mk.write(f"\t@if [ ! -f {p_file} ]; then \\\n")
                        mk.write(f"\t  rm {p_src} ; \\\n")
                        mk.write(f"\t  $(MAKE) {p_src} ; \\\n")
                        mk.write("\tfi\n")
                        mk.write(f"{p_file}: {p_file}.up\n\n")

                # <build>/foo.mod.d : aggregates the per-source .d files.
                dep_targets = ' '.join(
                    work_dir + s['fileName'] + '.d'
                    for s in submodule_list
                )
                mod_d    = mod_path + '.d'
                leaf_dep = work_dir + leaf + obj_dep_file
                mk.write(f"{mod_d}: {leaf_dep} {dep_targets}\n")
                mk.write(f"\t@echo {leaf_obj} > {mod_d}~\n")
                mk.write(f"\t@cat {leaf_dep} >> {mod_d}~\n")
                for submodule in submodule_list:
                    sub_file = submodule['fileName']
                    mk.write(f"\t@echo {work_dir}{sub_file}.o >> {mod_d}~\n")
                    mk.write(f"\t@cat {work_dir}{sub_file}.d >> {mod_d}~\n")
                mk.write(f"\t@if cmp -s {mod_d} {mod_d}~ ; then \\\n")
                mk.write(f"\t rm {mod_d}~ ; \\\n")
                mk.write("\telse \\\n")
                mk.write(f"\t mv {mod_d}~ {mod_d} ; \\\n")
                mk.write("\tfi\n\n")

                # <build>/foo.mod.gv : GraphViz source-dependency rule.
                gv_src = work_dir + leaf + source_file_name + '.gv'
                mk.write(
                    f"{mod_path}.gv: {leaf_dep} {gv_src}\n"
                )
                mk.write(
                    f"\t@echo {leaf}{source_file_name} > "
                    f"{work_dir}{mod_key}.gv\n"
                )

            # Module-list file: `<src>.m`.
            modules_list_file_name = re.sub(
                r'\.(f|f90)$', '.m', source_file_name,
                flags=re.IGNORECASE,
            )
            modules_list_path = (
                work_dir + leaf + modules_list_file_name
            )
            mk.write(f"{modules_list_path}:\n")
            director = '> '
            for module_file_name in mods:
                mk.write(
                    f"\t@echo {work_dir}{module_file_name} "
                    f"{director} {modules_list_path}\n"
                )
                director = '>>'
                sub_list = submodules_by_module.get(module_file_name.lower(), [])
                if sub_list:
                    smod_name = re.sub(r'\.mod$', '.smod', module_file_name)
                    mk.write(
                        f"\t@echo {work_dir}{smod_name} "
                        f"{director} {modules_list_path}\n"
                    )
            mk.write("\n")

            # Per-submodule `.m` files.
            for module_file_name in mods:
                module_name = re.sub(r'\.mod$', '', module_file_name)
                for submodule in submodules_by_module.get(
                        module_file_name.lower(), []):
                    sub_name = submodule['name'].lower()
                    sub_file = submodule['fileName']
                    sub_m    = f"{work_dir}{sub_file}.m"
                    mk.write(f"{sub_m}:\n")
                    mk.write(
                        f"\t@echo {work_dir}{module_name}@"
                        f"{sub_name}.smod > {sub_m}\n\n"
                    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv):
    if len(argv) != 2:
        print("Usage: moduleDependencies.py <sourceDirectory>",
              file=sys.stderr)
        sys.exit(1)
    install_directory = argv[1]
    build_path = os.environ['BUILDPATH']
    work_dir   = build_path + '/'
    blob_path  = work_dir + 'Makefile_Module_Dependencies.blob'

    locations = _load_directive_locations(
        work_dir + 'directiveLocations.xml',
    )

    source_root = os.path.join(install_directory, 'source')
    descriptors = _source_directory_descriptors(source_root, build_path)

    modules_per_file, cache_mtime = _load_cache(blob_path)
    have_per_file = cache_mtime is not None

    # Build the unstripped identifier list used for add/remove detection
    # (matches the Perl's first `@fileIdentifiers` loop).
    unstripped_identifiers = []
    for desc in descriptors:
        for name in _list_source_files(desc['path']):
            unstripped_identifiers.append(
                (desc['path'] + '/' + name).replace('/', '_'),
            )
    unstripped_set = set(unstripped_identifiers)

    force_rescan = False
    if have_per_file:
        if any(fid not in modules_per_file for fid in unstripped_identifiers):
            force_rescan = True
        if any(fid not in unstripped_set for fid in modules_per_file):
            force_rescan = True

    # Per-file scan.
    for desc in descriptors:
        for name in _list_source_files(desc['path']):
            file_path      = desc['path'] + '/' + name
            file_identifier = _file_identifier(file_path)

            rescan = True
            if have_per_file and file_identifier in modules_per_file:
                tracked = modules_per_file[file_identifier].get('files') or []
                stale = any(
                    os.path.exists(t)
                    and os.stat(t).st_mtime > cache_mtime
                    for t in tracked
                )
                rescan = bool(stale)
            if not (rescan or force_rescan):
                continue

            modules_per_file.pop(file_identifier, None)

            # functionClass -> synthesise submodule records.
            function_classes = extract_directives(file_path, 'functionClass')
            submodules_for_file = []
            for fc in function_classes:
                for fc_file in as_array(
                    (locations.get(fc['name']) or {}).get('file'),
                ):
                    submodules_for_file.extend(
                        _find_function_class_submodules(
                            fc['name'], fc_file, locations,
                        )
                    )

            entry = modules_per_file.setdefault(file_identifier, {})
            if submodules_for_file:
                entry['submodules'] = submodules_for_file
            entry.setdefault('files', []).append(file_path)
            entry['sourceFileName']           = name
            entry['sourceDirectoryDescriptor'] = desc
            entry.setdefault('modulesProvided', [])

            _scan_file([file_path], entry, source_root)

    # Build the submodule-by-module map.
    submodules_by_module = _build_submodule_map(modules_per_file)

    # Emit the Makefile and persist the cache.
    _write_makefile(
        work_dir + 'Makefile_Module_Dependencies',
        modules_per_file,
        submodules_by_module,
        work_dir,
    )

    with open(blob_path, 'wb') as fh:
        pickle.dump(modules_per_file, fh, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    main(sys.argv)
