#!/usr/bin/env python3
# Locate source files which have dependencies on modules.
# Python port of scripts/build/useDependencies.pl.
# Andrew Benson (ported to Python 2026)

import glob
import os
import pickle
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from Galacticus.Build.Directives import extract_directives, xml_to_dict_simple
from List.ExtraUtils import as_array, hash_list, smart_push


# External modules (ignored for dependency analysis of the source code).
EXTERNAL_MODULES = {
    'omp_lib', 'hdf5', 'h5tb', 'h5lt', 'h5global', 'h5fortran_types',
    'fox_common', 'fox_dom', 'fox_wxml', 'fox_utils', 'mpi', 'mpi_f08',
}

# Modules that require a library to be linked (module name -> library name).
MODULE_LIBRARIES = {
    'nearest_neighbors':   'ANN',
    'points_convex_hull':  'qhullcpp',
    'fftw3':               'fftw3',
    'fox_common':          'FoX_common',
    'fox_dom':             'FoX_dom',
    'fox_wxml':            'FoX_wxml',
    'fox_utils':           'FoX_utils',
    'hdf5':                'hdf5_fortran',
    'h5tb':                'hdf5hl_fortran',
    'vectors':             'blas',
    'models_likelihoods':  'matheval',
    'input_parameters':    'matheval',
    'interface_gsl':       'gsl',
    'output_versioning':   'git2',
}

# C includes that require a library to be linked (include name -> library).
INCLUDE_LIBRARIES = {
    'crypt': 'crypt',
}


def _load_xml_simple(path):
    """Load an XML file and convert it to a nested dict using XML::Simple's
    default auto-folding semantics (the caller expects e.g.
    `locations[name]['file']` to be a string for one file or a list for many,
    and `state_storables['functionClasses']` to be a dict keyed by class
    name)."""
    if not os.path.exists(path):
        return None
    tree = ET.parse(path)
    return xml_to_dict_simple(tree.getroot())


def _parse_preprocessor_directives(build_path):
    """Extract -D<SYMBOL> flags that are in force when compiling, gathered
    from: the top-level Makefile and any `Makefile*` under $BUILDPATH
    (respecting ifdef/ifeq/endif stacks); GALACTICUS_FCFLAGS; and the C
    compiler's built-in predefined macros."""
    directives = []
    makefiles = ['Makefile'] + sorted(glob.glob(os.path.join(build_path, 'Makefile*')))
    fcflags_re = re.compile(r'^\s*FCFLAGS\s*\+??=')
    dmacro_re = re.compile(r'-D([0-9A-Z]+)')
    for makefile_name in makefiles:
        if not os.path.exists(makefile_name):
            continue
        conditions_stack = [1]
        with open(makefile_name, 'r', encoding='utf-8', errors='replace') as handle:
            for line in handle:
                if line.startswith('ifdef '):
                    name = line[len('ifdef '):].strip().split()[0] if line[len('ifdef '):].strip() else ''
                    conditions_stack.append(1 if name in os.environ else 0)
                elif line.startswith('ifeq '):
                    conditions_stack.append(1)
                elif line.startswith('endif'):
                    if len(conditions_stack) > 1:
                        conditions_stack.pop()
                if fcflags_re.match(line) and all(conditions_stack):
                    for token in line.split():
                        m = dmacro_re.match(token)
                        if m:
                            directives.append(m.group(1))
    fcflags = os.environ.get('GALACTICUS_FCFLAGS', '')
    for token in fcflags.split():
        m = dmacro_re.match(token)
        if m:
            directives.append(m.group(1))
    compiler = os.environ.get('CCOMPILER', 'gcc')
    try:
        output = subprocess.run(
            [compiler, '-dM', '-E', '-'],
            input='', capture_output=True, text=True, check=False,
        ).stdout
        for line in output.splitlines():
            parts = line.split()
            if len(parts) >= 2:
                directives.append(parts[1])
    except (FileNotFoundError, OSError):
        pass
    return directives


def _identifier_for(path):
    ident = path.replace('/', '_')
    return re.sub(r'^\._?', '', ident)


def _collect_source_files(root_source_directory, build_path):
    """Find all f/f90/c/cpp/inc files under <root>/source/* and under
    $BUILDPATH/libgalacticus/. Returns a list of dicts with fileName,
    fullPathFileName and subDirectoryName fields."""
    source_subdirs = [os.path.join(root_source_directory, 'source')]
    first = source_subdirs[0]
    if os.path.isdir(first):
        for name in sorted(os.listdir(first)):
            path = os.path.join(first, name)
            if os.path.isdir(path) and not name.startswith('.'):
                source_subdirs.append(path)
    source_subdirs.append(os.path.join(build_path, 'libgalacticus'))

    ext_re = re.compile(r'\.(f|f90|c|cpp|inc)$', re.IGNORECASE)
    root_prefix = root_source_directory + '/source'
    files = []
    for subdir in source_subdirs:
        if not os.path.isdir(subdir):
            continue
        if subdir.startswith(root_prefix):
            sub_directory_name = subdir[len(root_prefix):].lstrip('/')
        elif subdir.startswith(build_path):
            sub_directory_name = subdir[len(build_path):].lstrip('/')
        else:
            sub_directory_name = ''
        for name in sorted(os.listdir(subdir)):
            if name.startswith('.#'):
                continue
            if not ext_re.search(name):
                continue
            files.append({
                'fileName': name,
                'fullPathFileName': os.path.join(subdir, name),
                'subDirectoryName': sub_directory_name,
            })
    return files


def main(root_source_directory):
    build_path = os.environ['BUILDPATH']
    work_directory = build_path.rstrip('/') + '/'

    locations = _load_xml_simple(os.path.join(work_directory, 'directiveLocations.xml'))
    state_storables = _load_xml_simple(os.path.join(work_directory, 'stateStorables.xml'))

    preprocessor_directives = _parse_preprocessor_directives(build_path)

    blob_path = os.path.join(work_directory, 'Makefile_Use_Dependencies.blob')
    uses_per_file = {}
    update_time = None
    if os.path.exists(blob_path):
        try:
            with open(blob_path, 'rb') as handle:
                uses_per_file = pickle.load(handle)
            update_time = os.path.getmtime(blob_path)
        except Exception:
            uses_per_file = {}
            update_time = None

    source_files = _collect_source_files(root_source_directory, build_path)

    # Detect file additions / removals - either forces a full rescan.
    identifiers = {_identifier_for(f['fullPathFileName']) for f in source_files}
    force_rescan = any(ident not in uses_per_file for ident in identifiers)
    if not force_rescan:
        for ident in uses_per_file:
            if ident in ('eventHooksManager',):
                continue
            if ident not in identifiers:
                force_rescan = True
                break

    event_hook_modules = []

    for source_file in source_files:
        full_path = source_file['fullPathFileName']
        file_identifier = _identifier_for(full_path)

        # Cache check: only rescan this file if anything in its recorded
        # dependency set is newer than the cached blob (or a full rescan
        # was forced by file additions / removals).
        rescan = True
        if not force_rescan and update_time is not None and file_identifier in uses_per_file:
            rescan = False
            for f in uses_per_file[file_identifier].get('files', []):
                try:
                    if os.stat(f).st_mtime > update_time:
                        rescan = True
                        break
                except OSError:
                    pass
        if not rescan and not force_rescan:
            continue

        uses_per_file.pop(file_identifier, None)
        entry = uses_per_file.setdefault(file_identifier, {})
        entry['files'] = [full_path]
        entry['modulesUsed'] = []
        entry['dependenciesExplicit'] = []
        entry['modulesProvided'] = {}
        entry['submodules'] = []
        entry['libraryDependencies'] = {}

        files_to_scan = [full_path]

        # Extract the directive kinds we care about from the main file.
        directives = {
            kind: extract_directives(full_path, kind)
            for kind in (
                'functionClass', 'inputParameter', 'enumeration',
                'eventHook', 'eventHookStatic', 'eventHookManager',
                'functionsGlobal',
            )
        }

        _process_directives(
            source_file, directives, entry, uses_per_file, event_hook_modules,
            files_to_scan, locations, state_storables, work_directory,
        )

        _scan_files(
            files_to_scan, entry, directives, root_source_directory,
            work_directory, preprocessor_directives, locations,
        )

    _aggregate_event_hooks(uses_per_file, event_hook_modules)
    submodules_map = _build_submodules_map(uses_per_file, work_directory)
    _write_makefile(
        source_files, uses_per_file, submodules_map, work_directory,
    )

    with open(blob_path, 'wb') as handle:
        pickle.dump(uses_per_file, handle)


def _aggregate_event_hooks(uses_per_file, event_hook_modules):
    """Merge the collected event_hook_modules into the file that declared the
    eventHooksManager directive. Mirrors the Perl block that populates
    usesPerFile{eventHooksFileIdentifier}{modulesUsed}."""
    manager = uses_per_file.get('eventHooksManager')
    if manager is None:
        return
    manager_modules = manager.setdefault('modules', [])
    manager_modules.extend(event_hook_modules)
    manager['modules'] = sorted(set(manager_modules))
    file_identifier = manager.get('fileIdentifier')
    if file_identifier is None:
        return
    record = uses_per_file.setdefault(file_identifier, {})
    modules_used = record.setdefault('modulesUsed', [])
    modules_used.extend(manager['modules'])
    record['modulesUsed'] = sorted(set(modules_used))


def _build_submodules_map(uses_per_file, work_directory):
    """Build a mapping from module-file (absolute) to the list of submodule
    names associated with it. Two sources feed this map: functionClass
    submodules (each functionClass-defining module gets an implementation
    submodule per instance) and explicit `submodule(parent:child) name`
    declarations."""
    submodules = {}
    for file_identifier, record in uses_per_file.items():
        if file_identifier == 'eventHooksManager':
            continue
        file_submodules = record.get('submodules')
        if not file_submodules:
            continue
        modules_provided = list(record.get('modulesProvided', {}).keys())
        if len(modules_provided) != 1:
            raise SystemExit(
                "useDependencies.py: submodules [" + ', '.join(file_submodules)
                + "] associated with multiple modules ["
                + ', '.join(modules_provided) + ']'
            )
        parent_key = work_directory + modules_provided[0].lower()
        submodules[parent_key] = list(file_submodules)
    for file_identifier, record in uses_per_file.items():
        if file_identifier == 'eventHooksManager':
            continue
        for submodule in record.get('submodulesProvided', []):
            parent_key = work_directory + submodule['moduleFileName']
            submodules.setdefault(parent_key, []).append(submodule['submoduleName'])
    return submodules


def _write_makefile(source_files, uses_per_file, submodules_map, work_directory):
    """Emit Makefile_Use_Dependencies with the object-file, library-file,
    dependency-file and GraphViz rules expected by the build."""
    target = work_directory + 'Makefile_Use_Dependencies'
    with open(target, 'w', encoding='utf-8') as out:
        for source_file in source_files:
            full_path = source_file['fullPathFileName']
            file_identifier = _identifier_for(full_path)
            record = uses_per_file.get(file_identifier)
            if record is None:
                continue
            object_file_name = re.sub(
                r'\.(f|f90|c|cpp)$', '.o', source_file['fileName'],
                flags=re.IGNORECASE,
            )
            dependency_file_name = re.sub(
                r'\.(f|f90|c|cpp|inc)$', '.d', source_file['fileName'],
                flags=re.IGNORECASE,
            )
            sub = source_file['subDirectoryName']
            work_subdir = work_directory + (sub + '/' if sub else '')

            # Library file rule (.fl): echoes library names, one per line.
            if not object_file_name.endswith('.Inc'):
                library_file_name = re.sub(r'\.o$', '.fl', object_file_name)
                lib_path = work_subdir + library_file_name
                out.write(lib_path + ':\n')
                redirect = '>'
                for lib in sorted(record.get('libraryDependencies', {})):
                    out.write('\t@echo ' + lib + ' ' + redirect + ' ' + lib_path + '\n')
                    redirect = '>>'

            modules_used = record.get('modulesUsed', [])
            explicit = record.get('dependenciesExplicit', [])
            if not (modules_used or explicit):
                continue

            # Drop self-provided modules from the modulesUsed set and sort.
            provided = record.get('modulesProvided', {})
            filtered_modules = sorted(
                set(m for m in modules_used if m not in provided)
            )
            record['modulesUsed'] = filtered_modules

            # Expand submodules into .smod references for any module used
            # that has submodules.
            submodules_used = []
            for module_used in filtered_modules:
                sub_names = submodules_map.get(module_used)
                if not sub_names:
                    continue
                module_base = re.sub(r'\.mod$', '', module_used)
                for sub_name in sub_names:
                    submodules_used.append(
                        module_base + '@' + sub_name.lower() + '.smod'
                    )

            obj_path = work_subdir + object_file_name
            dep_path = work_subdir + dependency_file_name
            omp_obj = work_directory + 'utility.OpenMP.workaround.o'

            out.write(
                obj_path + ': ' + omp_obj + ' '
                + ' '.join(list(filtered_modules) + list(submodules_used) + list(explicit))
                + ' Makefile\n'
            )

            dependencies_used = [
                m + '.d' if m.endswith('.mod') else m for m in filtered_modules
            ]
            explicit_d = [re.sub(r'\.o$', '.d', e) for e in explicit]
            out.write(
                dep_path + ': ' + ' '.join(dependencies_used + explicit_d) + '\n'
            )
            out.write('\t@echo ' + obj_path + ' > ' + dep_path + '~\n')
            out.write('\t@echo ' + omp_obj + ' >> ' + dep_path + '~\n')
            for explicit_dep in explicit:
                explicit_d_file = re.sub(r'\.o$', '.d', explicit_dep)
                prefix = '' if '/' in explicit_dep else work_directory
                out.write(
                    '\t@cat ' + prefix + explicit_d_file + ' >> ' + dep_path + '~\n'
                )
            for dep_used in dependencies_used:
                out.write('\t@cat ' + dep_used + ' >> ' + dep_path + '~\n')
            out.write('\t@sort -u ' + dep_path + '~ -o ' + dep_path + '~\n')
            out.write(
                '\t@if cmp -s ' + dep_path + ' ' + dep_path + '~ ; then \\\n'
            )
            out.write('\t rm ' + dep_path + '~ ; \\\n')
            out.write('\telse \\\n')
            out.write('\t mv ' + dep_path + '~ ' + dep_path + ' ; \\\n')
            out.write('\tfi\n\n')

            # GraphViz rule (.gv) for dependency visualisation. The Perl
            # maps from modulesUsed (not dependenciesUsed); the `.d` -> `.gv`
            # substitution is effectively a no-op because modulesUsed are
            # `.mod` paths. Preserved verbatim for parity.
            graphvizes_used = [
                m + '.gv' if m.endswith('.d') else m for m in filtered_modules
            ]
            gv_file_name = source_file['fileName'] + '.gv'
            gv_path = work_subdir + gv_file_name
            label = sub + source_file['fileName']
            out.write(
                gv_path + ': ' + dep_path + ' ' + ' '.join(graphvizes_used) + '\n'
            )
            out.write('\t@echo \\"' + label + '\\" > ' + gv_path + '\n')
            for explicit_dep in explicit:
                explicit_d_file = re.sub(r'\.o$', '.d', explicit_dep)
                prefix = '' if '/' in explicit_dep else work_directory
                out.write(
                    "\t@awk '{print \"\\\"" + label + "\\\" -> \\\"\"$$1\"\\\"\"}' "
                    + prefix + explicit_d_file + ' >> ' + gv_path + '\n'
                )
            for gv_used in graphvizes_used:
                out.write(
                    "\t@awk '{print \"\\\"" + label + "\\\" -> \\\"\"$$1\"\\\"\"}' "
                    + gv_used + ' >> ' + gv_path + '\n'
                )
                out.write(
                    "\t@cat `awk '{print \"" + work_directory
                    + "\"$$1\".gv\"}' " + gv_used + '` >> ' + gv_path + '\n'
                )
            out.write('\t@sort -u ' + gv_path + ' -o ' + gv_path + '\n\n')


def _process_directives(source_file, directives, entry, uses_per_file,
                        event_hook_modules, files_to_scan, locations,
                        state_storables, work_directory):
    """Translate per-file directive declarations into module and library
    dependencies. Mirrors the directive-handling block of the original Perl
    useDependencies.pl (roughly lines 155-323)."""

    modules_used = entry['modulesUsed']
    library_deps = entry['libraryDependencies']

    # functionClass: pull the implementation files into the scan stack and
    # add the function_classes module as a dependency.
    if directives['functionClass']:
        for function_class in directives['functionClass']:
            name = function_class.get('name')
            if name and locations:
                entry_locations = locations.get(name, {})
                smart_push(files_to_scan, entry_locations.get('file'))
        modules_used.append(work_directory + 'function_classes.mod')

    # eventHook: accumulate events_hooks.mod and any imported module names.
    if directives['eventHook']:
        modules_used.append(work_directory + 'events_hooks.mod')
        for event_hook in directives['eventHook']:
            imports = event_hook.get('import')
            if not imports:
                continue
            module_block = imports.get('module')
            if not module_block:
                continue
            if isinstance(module_block, dict) and 'name' in module_block:
                modules = as_array(module_block)
            elif isinstance(module_block, dict):
                modules = hash_list(module_block, key_as='name')
            else:
                modules = as_array(module_block)
            for module in modules:
                name = module.get('name') if isinstance(module, dict) else None
                if name:
                    event_hook_modules.append(
                        work_directory + name.lower() + '.mod'
                    )

    # eventHookManager: stash a cross-file reference so the aggregation step
    # can inject all event hook modules into the declaring file.
    if directives['eventHookManager']:
        sub = source_file['subDirectoryName']
        work_subdir = work_directory + (sub + '/' if sub else '')
        object_name = re.sub(
            r'\.(f|f90|c|cpp)$', '.o',
            work_subdir + source_file['fileName'],
            flags=re.IGNORECASE,
        )
        uses_per_file.setdefault('eventHooksManager', {})
        uses_per_file['eventHooksManager']['objectFileName'] = object_name
        uses_per_file['eventHooksManager']['fileIdentifier'] = _identifier_for(
            source_file['fullPathFileName']
        )

    # functionsGlobal: two flavours.
    if directives['functionsGlobal']:
        has_pointers = any(d.get('type') == 'pointers' for d in directives['functionsGlobal'])
        if has_pointers:
            modules_used.append(work_directory + 'error.mod')
            modules_used.append(work_directory + 'input_parameters.mod')
            files = locations.get('functionGlobal', {}).get('file') if locations else None
            for fg_file in as_array(files):
                for fg_directive in extract_directives(fg_file, set_root_element_type='functionGlobal'):
                    module_field = fg_directive.get('module')
                    if module_field is None:
                        continue
                    for module in as_array(module_field):
                        if not isinstance(module, str):
                            continue
                        m = re.match(r'^([a-zA-Z0-9_]+)', module)
                        if not m:
                            continue
                        module_name = m.group(1).lower()
                        if module_name == 'iso_c_binding':
                            continue
                        modules_used.append(work_directory + module_name + '.mod')
        has_establish = any(d.get('type') == 'establish' for d in directives['functionsGlobal'])
        if has_establish:
            files = locations.get('functionGlobal', {}).get('file') if locations else None
            for fg_file in as_array(files):
                module_name = _find_containing_module(fg_file, state_storables)
                if module_name is None:
                    raise SystemExit(
                        "useDependencies.py: unable to locate containing "
                        "module for global function in file '" + fg_file + "'"
                    )
                modules_used.append(work_directory + module_name.lower() + '.mod')

    # eventHookStatic: look up the module that holds the statically-hooked
    # event in each of the files listed for that event name.
    if directives['eventHookStatic']:
        for static in directives['eventHookStatic']:
            name = static.get('name')
            if not name:
                continue
            files = locations.get(name, {}).get('file') if locations else None
            for hook_file in as_array(files):
                module_name = _find_containing_module(hook_file, state_storables)
                if module_name is None:
                    raise SystemExit(
                        "useDependencies.py: unable to locate containing "
                        "module for static event '" + name
                        + "' in file '" + hook_file + "'"
                    )
                modules_used.append(work_directory + module_name.lower() + '.mod')

    # Implicit module dependencies driven by directive kinds.
    if directives['functionClass'] or directives['inputParameter']:
        modules_used.append(work_directory + 'input_parameters.mod')
    if any(
        d.get('encodeFunction') == 'yes' and 'errorValue' not in d
        or d.get('decodeFunction') == 'yes' and 'errorValue' not in d
        for d in directives['enumeration']
    ):
        modules_used.append(work_directory + 'error.mod')
    if directives['enumeration']:
        modules_used.append(work_directory + 'enumerations.mod')

    # Modules referenced by the method definitions of each functionClass.
    for function_class in directives['functionClass']:
        methods = function_class.get('method')
        if not methods:
            continue
        # If method is a single dict with a 'name' field it's a single method;
        # otherwise the dict is keyed by method name.
        if isinstance(methods, dict) and 'name' in methods:
            method_iter = [methods]
        elif isinstance(methods, dict):
            method_iter = list(methods.values())
        else:
            method_iter = methods  # list
        for method in method_iter:
            modules_field = method.get('modules')
            if modules_field is None:
                continue
            if isinstance(modules_field, dict):
                module_names = [k.lower() for k in sorted(modules_field.keys())]
            else:
                module_names = modules_field.lower().split()
            for module_name in module_names:
                if module_name in EXTERNAL_MODULES:
                    continue
                modules_used.append(work_directory + module_name + '.mod')


# Line-level regex patterns used by _scan_files.
_LIB_COMMENT_RE   = re.compile(r'^\s*!;\s*([a-zA-Z0-9_]+)\s*$')
_C_INCLUDE_RE     = re.compile(r'^\s*\#include\s+<([a-zA-Z0-9_]+)\.h>')
_PP_IFDEF_RE      = re.compile(r'^\#ifdef\s+([0-9A-Za-z_]+)\s*$')
_PP_IF_RE         = re.compile(r'^\#if\s')
_PP_IFNDEF_RE     = re.compile(r'^\#ifndef\s+([0-9A-Z_]+)\s*$')
_PP_ENDIF_RE      = re.compile(r'^\#endif\s*$')
_PP_ELSE_RE       = re.compile(r'^\#else\s*$')
_USE_RE           = re.compile(
    r'^\s*(!\$\s)??\s*use\s*(?:::|\s)\s*([a-zA-Z0-9_]+)',
    re.IGNORECASE,
)
_OMP_PARALLEL_RE  = re.compile(r'^\s*\!\$omp\s+parallel', re.IGNORECASE)
_OMP_CRITICAL_RE  = re.compile(r'^\s*\!\$omp\s+critical\s*\([a-z0-9_]+\)', re.IGNORECASE)
_EXPLICIT_DEP_RE  = re.compile(r'^\s*(?:\!|//):\s*(.*)$')
_PROGRAM_RE       = re.compile(r'^\s*program\s', re.IGNORECASE)
_MODULE_DECL_RE   = re.compile(r'^\s*module\s+([a-zA-Z0-9_]+)')
_SUBMODULE_RE     = re.compile(
    r'^\s*submodule\s*\(\s*([a-z0-9_]+)(?::([a-z0-9_]+))??\s*\)\s*([a-zA-Z0-9_]+)',
    re.IGNORECASE,
)
_INCLUDE_LINE_RE  = re.compile(r'^\s*include\s+[\'"]([\w\.\-]+)[\'"]', re.IGNORECASE)
_XML_START_RE     = re.compile(r'^\s*!!\[')
_XML_END_RE       = re.compile(r'^\s*!!\]')
_LATEX_START_RE   = re.compile(r'^\s*!!\{')
_LATEX_END_RE     = re.compile(r'^\s*!!\}')
_TYPE_INC_RE      = re.compile(r'(.*)\.type\.inc$')


def _scan_files(files_to_scan, entry, directives, root_source_directory,
                work_directory, preprocessor_directives, locations):
    """Scan every file on the stack line-by-line, growing the stack as
    include statements are discovered, and populate the entry with modules
    used / provided / submodules / library dependencies / explicit
    dependencies."""

    preprocessor_set = set(preprocessor_directives)
    root_source_subdir = os.path.join(root_source_directory, 'source') + '/'

    modules_used = entry['modulesUsed']
    modules_provided = entry['modulesProvided']
    submodules = entry['submodules']
    library_deps = entry['libraryDependencies']
    dependencies_explicit = entry['dependenciesExplicit']
    submodules_provided = entry.setdefault('submodulesProvided', [])

    while files_to_scan:
        full_path = files_to_scan.pop()
        # Ensure the file exists (a build-generated include file may need to
        # be made first).
        if not os.path.exists(full_path):
            leaf_match = re.search(r'/([\w\.]+?)$', full_path)
            leaf = leaf_match.group(1) if leaf_match else full_path
            subprocess.run(['make', leaf], check=False)
        if full_path.lower().endswith('.cpp'):
            library_deps['stdc++'] = 1

        conditional_stack = []
        conditionally_compile = True
        in_xml = False
        in_latex = False

        try:
            handle = open(full_path, 'r', encoding='utf-8', errors='replace')
        except OSError as err:
            raise SystemExit(
                "useDependencies.py: can't open input file: " + str(err)
            )
        with handle:
            for line in handle:
                if _XML_END_RE.match(line):
                    in_xml = False
                if _LATEX_END_RE.match(line):
                    in_latex = False
                if in_xml or in_latex:
                    continue

                m = _LIB_COMMENT_RE.match(line)
                if m:
                    library_deps[m.group(1)] = 1

                m = _C_INCLUDE_RE.match(line)
                if m:
                    include_file = m.group(1).lower()
                    library = INCLUDE_LIBRARIES.get(include_file)
                    if library is not None and conditionally_compile:
                        library_deps[library] = 1

                if line.startswith('#'):
                    pm = _PP_IFDEF_RE.match(line)
                    if pm:
                        conditional_stack.append({'name': pm.group(1), 'state': 1})
                    elif _PP_IF_RE.match(line):
                        conditional_stack.append({'name': 'conditional', 'state': 1})
                    else:
                        pm = _PP_IFNDEF_RE.match(line)
                        if pm:
                            conditional_stack.append({'name': pm.group(1), 'state': 0})
                        elif _PP_ENDIF_RE.match(line):
                            if conditional_stack:
                                conditional_stack.pop()
                        elif _PP_ELSE_RE.match(line):
                            if conditional_stack:
                                conditional_stack[-1]['state'] = 1 - conditional_stack[-1]['state']
                    conditionally_compile = True
                    for cond in conditional_stack:
                        present = cond['name'] in preprocessor_set
                        active = cond['state'] if present else 1 - cond['state']
                        if active == 0:
                            conditionally_compile = False
                            break

                if conditionally_compile:
                    m = _USE_RE.match(line)
                    if m:
                        used_module = m.group(2).lower()
                        library = MODULE_LIBRARIES.get(used_module)
                        if library is not None:
                            library_deps[library] = 1
                        if used_module not in EXTERNAL_MODULES:
                            modules_used.append(
                                work_directory + used_module + '.mod'
                            )
                    if _OMP_PARALLEL_RE.match(line):
                        modules_used.append(work_directory + 'events_filters.mod')
                    if _OMP_CRITICAL_RE.match(line):
                        modules_used.append(work_directory + 'openmp_utilities_data.mod')
                    m = _EXPLICIT_DEP_RE.match(line)
                    if m:
                        dependencies_explicit.extend(m.group(1).split())
                    if _PROGRAM_RE.match(line):
                        modules_used.append(work_directory + 'iso_varying_string.mod')
                    m = _MODULE_DECL_RE.match(line)
                    if m:
                        module_name = m.group(1)
                        if module_name != 'procedure':
                            modules_provided[module_name + '.mod'] = 1
                            # A functionClass-derived submodule is generated
                            # for every implementation, tracked per-class.
                            if directives['functionClass'] and locations:
                                for fc in directives['functionClass']:
                                    fc_name = fc.get('name')
                                    if not fc_name:
                                        continue
                                    files = locations.get(fc_name, {}).get('file')
                                    for fc_file in as_array(files):
                                        for instance in extract_directives(fc_file, fc_name):
                                            instance_name = instance.get('name')
                                            if instance_name:
                                                submodules.append(instance_name + '_')
                    m = _SUBMODULE_RE.match(line)
                    if m:
                        parent = m.group(1).lower()
                        child = m.group(3).lower()
                        submodules_provided.append({
                            'submoduleName': child,
                            'moduleFileName': parent + '.mod',
                        })
                    m = _INCLUDE_LINE_RE.match(line)
                    if m:
                        preprocessed_include = m.group(1)
                        raw_include = re.sub(r'\.inc$', '.Inc', preprocessed_include)
                        raw_path = root_source_subdir + raw_include
                        pp_path = work_directory + preprocessed_include
                        if os.path.exists(raw_path):
                            files_to_scan.append(raw_path)
                            entry['files'].append(raw_path)
                        elif os.path.exists(pp_path):
                            files_to_scan.append(pp_path)
                            entry['files'].append(pp_path)
                        else:
                            tm = _TYPE_INC_RE.match(preprocessed_include)
                            if tm and locations:
                                base = tm.group(1)
                                loc_files = locations.get(base, {}).get('file')
                                smart_push(files_to_scan, loc_files)
                                smart_push(entry['files'], loc_files)

                if _XML_START_RE.match(line):
                    in_xml = True
                if _LATEX_START_RE.match(line):
                    in_latex = True


def _find_containing_module(source_path, state_storables):
    """Scan a source file for the first `module <name>` declaration, falling
    back to the functionClass's declared module in stateStorables.xml."""
    module_name = None
    function_class_name = None
    module_decl = re.compile(r'^\s*module\s+([a-zA-Z0-9_]+)')
    module_procedure = re.compile(r'^\s*module\s+procedure\s+([a-zA-Z0-9_]+)')
    xml_open = re.compile(r'^\s*<(\S+)')
    try:
        with open(source_path, 'r', encoding='utf-8', errors='replace') as handle:
            for line in handle:
                m = module_decl.match(line)
                if m and not module_procedure.match(line):
                    module_name = m.group(1)
                    break
                m = xml_open.match(line)
                if m:
                    element_name = m.group(1)
                    class_key = element_name + 'Class'
                    function_classes = (
                        state_storables.get('functionClasses') if state_storables else None
                    )
                    if function_classes and class_key in function_classes:
                        function_class_name = class_key
    except OSError:
        return None
    if module_name is None and function_class_name is not None:
        function_classes = state_storables.get('functionClasses') if state_storables else {}
        entry = function_classes.get(function_class_name, {})
        if isinstance(entry, dict):
            module_name = entry.get('module')
    return module_name


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit('Usage: useDependencies.py <sourcedir>')
    main(sys.argv[1])
