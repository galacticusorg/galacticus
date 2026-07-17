#!/usr/bin/env python3
"""Generate `Makefile_Use_Dependencies`: per-source-file dependency rules
covering module `use`s, library linking, OpenMP helpers, `eventHook` /
`functionClass` wiring, and GraphViz source-tree visualisation.

Follows the same preprocessor-conditional and directive-aware scanning
logic as the Perl original, including detection of `-D<NAME>` flags in
`FCFLAGS` (across active `ifdef`/`ifeq` branches of every `Makefile*`)
plus any macros defined via `$(GALACTICUS_FCFLAGS)` or produced by the
configured C compiler's `-dM -E` output.

Mirrors scripts/build/useDependencies.pl.
Andrew Benson (ported to Python 2026).
"""

import os
import pickle
import re
import subprocess
import sys
import xml.etree.ElementTree as ET


from Galacticus.Build.Directives   import extract_directives
from Galacticus.Build.FileChanges  import update as file_changes_update
from Galacticus.Build.ParallelScan import scan as parallel_scan
from List.ExtraUtils              import as_array, hash_list, smart_push
from XML.Utils                    import xml_to_dict
from Galacticus.Build.Libraries import (
    EXTERNAL_MODULES,
    MODULE_LIBRARIES,
    INCLUDE_LIBRARIES,
)
from Galacticus.Build.ScanCache import (
    file_identifier as _file_identifier,
    load_cache      as _load_cache,
    prune           as _prune_cache,
)


# Version stamp for the Makefile_Use_Dependencies.blob cache. Bump when the
# scan rules change so stale-rule entries are discarded.
_BLOB_VERSION = 2

# Directives consulted per source file (module-level so the parallel worker can
# see it). Order matters only in that it is fixed.
_DIRECTIVE_NAMES = (
    'functionClass', 'inputParameter', 'enumeration',
    'eventHook', 'eventHookStatic', 'eventHookManager',
    'functionsGlobal', 'objectDestructor',
)

# Read-only context shared with forked scan workers (inherited via copy-on-write,
# so it is not pickled per task). Populated in main() before the scan.
_WORKER = {}


def _scan_one(task):
    """Worker: do the full per-file scan for one source file and return
    `(file_identifier, entry, event_hook_modules, manager)`.

    Mirrors the body of main()'s per-file loop exactly, but captures the two
    cross-file side effects -- the ordered `event_hook_modules` pool and the
    singleton `eventHooksManager` record -- into locals so they can be merged
    deterministically (in file order) by the caller instead of mutating shared
    state from a worker.
    """
    file_identifier, sf = task
    file_path = sf['fullPathFileName']

    entry = {
        'files':                [file_path],
        'modulesUsed':          [],
        'dependenciesExplicit': [],
        'modulesProvided':      {},
        'submodules':           [],
        'submodulesProvided':   [],
        'libraryDependencies':  {},
    }

    # One read/parse of the file for all directive names (extract_directives
    # performs a full file read per call, so per-name calls would cost
    # len(_DIRECTIVE_NAMES) reads). rootElementType is popped so the per-name
    # directive dicts are identical to those a single-name call returns.
    directives = {name: [] for name in _DIRECTIVE_NAMES}
    for directive in extract_directives(file_path, _DIRECTIVE_NAMES,
                                        set_root_element_type=True):
        directives[directive.pop('rootElementType')].append(directive)

    file_names_to_process = [file_path]
    event_hook_modules = []
    manager = {}   # captures any uses_per_file['eventHooksManager'] write

    _apply_directive_requirements(
        entry, sf, directives, _WORKER['locations'], _WORKER['state_storables'],
        _WORKER['work_dir'], file_identifier, event_hook_modules,
        manager, file_names_to_process,
    )
    _scan_source_file(
        entry, file_names_to_process, sf, directives,
        _WORKER['locations'], _WORKER['root_source_dir'], _WORKER['work_dir'],
        _WORKER['preprocessor_set'],
    )
    return file_identifier, entry, event_hook_modules, manager


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------

def _load_xml(path):
    if not os.path.exists(path):
        return None
    return xml_to_dict(ET.parse(path).getroot())


# ---------------------------------------------------------------------------
# Preprocessor directive collection
# ---------------------------------------------------------------------------

_IFDEF_RE   = re.compile(r'^\s*ifdef\s+([A-Za-z_][A-Za-z0-9_]*)')
_IFNDEF_RE  = re.compile(r'^\s*ifndef\s+([A-Za-z_][A-Za-z0-9_]*)')
_IFEQ_RE    = re.compile(r'^\s*if(?:n?)eq\s')
_ELSE_RE    = re.compile(r'^\s*else\b(.*)$')
_ENDIF_RE   = re.compile(r'^\s*endif\b')
_FCFLAGS_RE = re.compile(r'^\s*FCFLAGS\s*\+??=')
_DASH_D_RE  = re.compile(r'-D([0-9A-Za-z_]+)')


def _collect_preprocessor_directives(build_path):
    """Return the list of preprocessor macro names active for the build.

    Parses every `Makefile*` (under the project root and under
    `$BUILDPATH`) for `-D<NAME>` flags on active `FCFLAGS += …` lines,
    honouring `ifdef`/`ifeq`/`endif` conditionals the same way the Perl
    original does.  Falls back to the C compiler's `-dM -E` output for
    any toolchain-provided macros, and merges `GALACTICUS_FCFLAGS` env
    var contents if set.
    """
    directives = []
    makefiles  = ['Makefile']
    try:
        for name in sorted(os.listdir(build_path)):
            if name.startswith('Makefile'):
                makefiles.append(os.path.join(build_path, name))
    except OSError:
        pass

    for makefile in makefiles:
        if not os.path.exists(makefile):
            continue
        # Stack entries: True (branch known active), False (known
        # inactive), None (unevaluable — `ifeq`/`ifneq` conditions over make
        # variables, and any `else <condition>` chain). A -D flag is
        # collected unless a *definitely false* branch encloses it: for
        # unevaluable conditions every branch is collected, deliberately
        # over-approximating the macro set (over-depending is benign;
        # under-depending silently drops dependencies). `ifdef`/`ifndef` are
        # approximated against the environment. Every opener pushes and
        # every `endif` pops, so an unrecognized condition can no longer
        # desynchronize the stack (previously `ifneq` did not push but its
        # `endif` still popped, deactivating or reactivating outer branches
        # at random).
        conditions_stack = [True]
        with open(makefile, 'r', errors='replace') as fh:
            for line in fh:
                m = _IFDEF_RE.match(line)
                if m:
                    conditions_stack.append(m.group(1) in os.environ)
                    continue
                m = _IFNDEF_RE.match(line)
                if m:
                    conditions_stack.append(m.group(1) not in os.environ)
                    continue
                if _IFEQ_RE.match(line):
                    conditions_stack.append(None)
                    continue
                m = _ELSE_RE.match(line)
                if m:
                    if len(conditions_stack) > 1:
                        taken = conditions_stack.pop()
                        if m.group(1).strip():
                            # `else if…`: a chained condition — unevaluable
                            # in general, so collect from it.
                            conditions_stack.append(None)
                        elif taken is None:
                            conditions_stack.append(None)
                        else:
                            conditions_stack.append(not taken)
                    continue
                if _ENDIF_RE.match(line):
                    if len(conditions_stack) > 1:
                        conditions_stack.pop()
                    continue
                if (_FCFLAGS_RE.match(line)
                        and all(c is not False for c in conditions_stack)):
                    for tok in line.split():
                        m = _DASH_D_RE.match(tok)
                        if m:
                            directives.append(m.group(1))

    fcflags_env = os.environ.get('GALACTICUS_FCFLAGS')
    if fcflags_env:
        for tok in fcflags_env.split():
            m = _DASH_D_RE.match(tok)
            if m:
                directives.append(m.group(1))

    compiler = os.environ.get('CCOMPILER', 'gcc')
    try:
        result = subprocess.run(
            [compiler, '-dM', '-E', '-'],
            input='', text=True, capture_output=True, check=False,
        )
        for line in result.stdout.splitlines():
            cols = line.split()
            if len(cols) >= 2:
                directives.append(cols[1])
    except OSError:
        pass

    return directives


# ---------------------------------------------------------------------------
# Source-file discovery
# ---------------------------------------------------------------------------

_SOURCE_SUFFIX_RE = re.compile(r'\.(f|f90|c|cpp|inc)$', re.IGNORECASE)


# ---------------------------------------------------------------------------
# Directive-based special-case handling
# ---------------------------------------------------------------------------

_MODULE_LINE_FIRST_RE = re.compile(r'^\s*module\s+([a-zA-Z0-9_]+)')
_MODULE_PROC_RE       = re.compile(r'^\s*module\s+procedure\s+([a-zA-Z0-9_]+)')
_DIRECTIVE_OPEN_RE    = re.compile(r'^\s*<(\S+)')


def _find_containing_module(file_path, state_storables):
    """Walk `file_path` line-by-line looking for the first `module <name>`
    declaration (not `module procedure`).  If none is found, fall back to
    the first `<X>` directive whose `XClass` is registered as a functionClass
    in `state_storables` and use its `module` attribute.

    Mirrors the two `while (my $line = <$file>)` blocks at
    useDependencies.pl:232-252 and 259-279.  Returns the bare module name
    or `None` if no module could be located.
    """
    module_name        = None
    function_class_key = None

    with open(file_path, 'r', errors='replace') as fh:
        for line in fh:
            m1 = _MODULE_LINE_FIRST_RE.match(line)
            m2 = _MODULE_PROC_RE      .match(line)
            if m1 and not m2:
                module_name = m1.group(1)
                break
            m3 = _DIRECTIVE_OPEN_RE.match(line)
            if m3:
                element = m3.group(1)
                candidate_key = element + 'Class'
                if (state_storables
                        and candidate_key in (
                            state_storables.get('functionClasses') or {}
                        )):
                    function_class_key = candidate_key

    if module_name is None and function_class_key is not None:
        fc = state_storables['functionClasses'][function_class_key]
        if isinstance(fc, dict):
            module_name = fc.get('module')
    return module_name


def _functionclass_state_storables_map(state_storables):
    """XML::Simple under `KeyAttr` implicit keying gives us
    `state_storables['functionClasses']` as either a list of dicts (default
    xml_to_dict behaviour) or, when the downstream code keys by name, a
    dict.  Normalise to `{name -> {module, ...}}` for easy lookup.
    """
    if not state_storables:
        return {}
    fcs = state_storables.get('functionClasses')
    if isinstance(fcs, dict):
        return fcs
    if isinstance(fcs, list):
        out = {}
        for entry in fcs:
            if isinstance(entry, dict) and 'name' in entry:
                out[entry['name']] = entry
        return out
    return {}


def _modules_from_eventhook(event_hook):
    """Return every module-name dict an `eventHook` directive's `import` block
    references.  Mirrors useDependencies.pl:198-208.
    """
    imp = event_hook.get('import')
    if not isinstance(imp, dict):
        return []
    module_field = imp.get('module')
    if module_field is None:
        return []
    if isinstance(module_field, dict) and 'name' in module_field:
        return as_array(module_field)
    if isinstance(module_field, dict):
        return hash_list(module_field, key_as='name')
    # List of dicts (already normalised by xml_to_dict).
    return as_array(module_field)


def _method_iter(function_class):
    """Iterate over a functionClass directive's methods.  Mirrors
    useDependencies.pl:298:
        exists($method->{'name'}) ? $method : map { $method->{$k} } keys
    """
    method = function_class.get('method')
    if method is None:
        return []
    if isinstance(method, dict):
        if 'name' in method:
            return [method]
        return list(method.values())
    if isinstance(method, list):
        return method
    return []


def _method_module_names(method):
    """Return the list of module names referenced by a single method.

    Mirrors useDependencies.pl:298-321.  `modules` can be:
      * an attribute-form scalar `modules="A B C"` -- space-separated names;
      * a dict from one `<modules>` block with `<name>` (+ `<only>`) children
        -- Python xml_to_dict representation;
      * a list of such dicts when the method has multiple `<modules>` blocks;
      * an XML::Simple-style dict keyed by module name (legacy Perl
        representation preserved for compatibility).
    """
    mods = method.get('modules')
    if mods is None:
        return []
    # Attribute-form scalar.
    if isinstance(mods, str):
        return [tok.lower() for tok in mods.split()]
    # One `<modules>` block -> dict.  If it describes one module directly
    # (has a `name` key) extract that name; otherwise fall back to
    # XML::Simple auto-keyed-by-name form.
    if isinstance(mods, dict):
        if 'name' in mods:
            name = mods['name']
            if isinstance(name, str):
                return [name.lower()]
            return []
        return [name.lower() for name in sorted(mods.keys())]
    # Multiple `<modules>` siblings -> list of dicts (each with a `name` key).
    if isinstance(mods, list):
        names = []
        for item in mods:
            if isinstance(item, dict):
                name = item.get('name')
                if isinstance(name, str):
                    names.append(name.lower())
            elif isinstance(item, str):
                names.append(item.lower())
        return names
    return []


# Per-process memos for _apply_directive_requirements: every file carrying a
# `functionsGlobal` (or `eventHookStatic`) directive needs the same
# information about the functionGlobal location files, which previously was
# re-derived — including full re-reads of those files — once per requesting
# source file. Workers are forked processes, so each builds these at most
# once per run.
_FG_POINTER_MODULES    = None
_CONTAINING_MODULE_MEMO = {}


def _functionglobal_pointer_modules(locations):
    """Module names imported by `functionsGlobal type="pointers"` users: the
    (deduplicated-in-order) modules declared by every `functionGlobal`
    directive, excluding iso_c_binding."""
    global _FG_POINTER_MODULES
    if _FG_POINTER_MODULES is None:
        modules = []
        for fg_file in as_array(
            (locations.get('functionGlobal') or {}).get('file')
            if locations else None,
        ):
            for d in extract_directives(fg_file, 'functionGlobal'):
                module_val = d.get('module')
                if module_val is None:
                    continue
                for module in as_array(module_val):
                    name = str(module).strip()
                    m = re.match(r'^([a-zA-Z0-9_]+)', name)
                    if m and m.group(1).lower() != 'iso_c_binding':
                        modules.append(m.group(1).lower())
        _FG_POINTER_MODULES = modules
    return _FG_POINTER_MODULES


def _find_containing_module_memo(file_name, fc_map):
    """Memoized _find_containing_module (the fc_map is identical across all
    calls within one run)."""
    if file_name not in _CONTAINING_MODULE_MEMO:
        _CONTAINING_MODULE_MEMO[file_name] = _find_containing_module(
            file_name, {'functionClasses': fc_map},
        )
    return _CONTAINING_MODULE_MEMO[file_name]


def _apply_directive_requirements(entry, source_file, directives,
                                  locations, state_storables,
                                  work_dir, file_identifier,
                                  event_hook_modules, uses_per_file,
                                  file_names_to_process):
    """Consume the directives dict and stage all the implicit module /
    library dependencies they imply.

    Mirrors useDependencies.pl:181-323.  Mutates `entry` (the per-file
    record), `event_hook_modules`, `uses_per_file`, and
    `file_names_to_process` in place.
    """
    fc_map = _functionclass_state_storables_map(state_storables)

    # functionClass: pull instance source files onto the scan stack, and
    # always require function_classes.mod.
    if directives['functionClass']:
        for fc in directives['functionClass']:
            smart_push(
                file_names_to_process,
                (locations.get(fc['name']) or {}).get('file')
                if locations else None,
            )
        entry['modulesUsed'].append(work_dir + 'function_classes.mod')

    # eventHook: track a cross-file pool of imported modules that will be
    # injected into the eventHooksManager-owning file during aggregation.
    if directives['eventHook']:
        entry['modulesUsed'].append(work_dir + 'events_hooks.mod')
        for event_hook in directives['eventHook']:
            for module in _modules_from_eventhook(event_hook):
                name = (module.get('name') or '').lower()
                if name:
                    event_hook_modules.append(work_dir + name + '.mod')

    # eventHookManager: remember which file carries the manager definition.
    if directives['eventHookManager']:
        sub_dir   = source_file['subDirectoryName']
        work_sub  = work_dir + (sub_dir + '/' if sub_dir else '')
        obj_file  = re.sub(
            r'\.(f|f90|c|cpp)$', '.o', source_file['fileName'],
            flags=re.IGNORECASE,
        )
        uses_per_file.setdefault('eventHooksManager', {})['objectFileName'] = (
            work_sub + obj_file
        )
        uses_per_file['eventHooksManager']['fileIdentifier'] = file_identifier

    # functionsGlobal: type=pointers and type=establish impose different
    # module requirements (see useDependencies.pl:215-253).
    if directives['functionsGlobal']:
        has_pointers  = any(
            d.get('type') == 'pointers'  for d in directives['functionsGlobal']
        )
        has_establish = any(
            d.get('type') == 'establish' for d in directives['functionsGlobal']
        )

        if has_pointers:
            entry['modulesUsed'].append(work_dir + 'error.mod')
            entry['modulesUsed'].append(work_dir + 'input_parameters.mod')
            for module_name in _functionglobal_pointer_modules(locations):
                entry['modulesUsed'].append(
                    work_dir + module_name + '.mod'
                )

        if has_establish:
            for fg_file in as_array(
                (locations.get('functionGlobal') or {}).get('file')
                if locations else None,
            ):
                module_name = _find_containing_module_memo(fg_file, fc_map)
                if module_name is None:
                    sys.exit(
                        "useDependencies.py: unable to locate containing "
                        f"module for global function in file '{fg_file}'"
                    )
                entry['modulesUsed'].append(
                    work_dir + module_name.lower() + '.mod'
                )

    # eventHookStatic: pull in the module that *declares* each static event.
    if directives['eventHookStatic']:
        for ehs in directives['eventHookStatic']:
            ehs_files = as_array(
                (locations.get(ehs['name']) or {}).get('file')
                if locations else None,
            )
            for ehs_file in ehs_files:
                module_name = _find_containing_module_memo(ehs_file, fc_map)
                if module_name is None:
                    sys.exit(
                        "useDependencies.py: unable to locate containing "
                        f"module for static event '{ehs['name']}' in file "
                        f"'{ehs_file}'"
                    )
                entry['modulesUsed'].append(
                    work_dir + module_name.lower() + '.mod'
                )

    # functionClass / inputParameter → input_parameters.mod.
    if directives['functionClass'] or directives['inputParameter']:
        entry['modulesUsed'].append(work_dir + 'input_parameters.mod')

    # enumeration: errorless encode/decode forms need error.mod; any
    # enumeration directive needs enumerations.mod.
    if any(
        e.get('encodeFunction') == 'yes' and 'errorValue' not in e
        for e in directives['enumeration']
    ) or any(
        e.get('decodeFunction') == 'yes' and 'errorValue' not in e
        for e in directives['enumeration']
    ):
        entry['modulesUsed'].append(work_dir + 'error.mod')
    if directives['enumeration']:
        entry['modulesUsed'].append(work_dir + 'enumerations.mod')

    # objectDestructor: the directive processor (ObjectBuilder.py) injects a
    # `use :: Error` into every file containing this directive (for the
    # negative-reference-count Error_Report call), so such files depend on
    # error.mod. Most files acquire this transitively via their other uses, but
    # a file whose only other dependencies are themselves error.mod-free (e.g.
    # objects.pool.F90, which uses just function_classes.mod) needs it stated
    # explicitly to be ordered correctly in a parallel build.
    if directives['objectDestructor']:
        entry['modulesUsed'].append(work_dir + 'error.mod')

    # functionClass methods: each method may declare module dependencies
    # that this file inherits.
    for fc in directives['functionClass']:
        for method in _method_iter(fc):
            for module_name in _method_module_names(method):
                if module_name in EXTERNAL_MODULES:
                    continue
                entry['modulesUsed'].append(
                    work_dir + module_name + '.mod'
                )


# ---------------------------------------------------------------------------
# Per-file scan (preprocessor state machine + Fortran pattern extraction)
# ---------------------------------------------------------------------------

_XML_OPEN_RE    = re.compile(r'^\s*!!\[')
_XML_CLOSE_RE   = re.compile(r'^\s*!!\]')
_LATEX_OPEN_RE  = re.compile(r'^\s*!!\{')
_LATEX_CLOSE_RE = re.compile(r'^\s*!!\}')

_LIB_HINT_RE        = re.compile(r'^\s*!;\s*([a-zA-Z0-9_]+)\s*$')
_C_INCLUDE_RE       = re.compile(r'^\s*#include\s+<([a-zA-Z0-9_]+)\.h>')
_PP_IFDEF_RE        = re.compile(r'^#ifdef\s+([0-9A-Za-z_]+)\s*$')
_PP_IFNDEF_RE       = re.compile(r'^#ifndef\s+([0-9A-Za-z_]+)\s*$')
_PP_IF_RE           = re.compile(r'^#if\s')
_PP_ELIF_RE         = re.compile(r'^#elif\s')
_PP_ENDIF_RE        = re.compile(r'^#endif\s*$')
_PP_ELSE_RE         = re.compile(r'^#else\s*$')
_USE_LINE_RE        = re.compile(
    r'^\s*(!\$\s)??\s*use\s*(::|\s)\s*([a-zA-Z0-9_]+)',
    re.IGNORECASE,
)
_OMP_PARALLEL_RE    = re.compile(r'^\s*!\$omp\s+parallel', re.IGNORECASE)
_OMP_CRITICAL_RE    = re.compile(
    r'^\s*!\$omp\s+critical\s*\([a-z0-9_]+\)', re.IGNORECASE,
)
_EXPLICIT_DEPS_RE   = re.compile(r'^\s*(?:!|//):\s*(.*)$')
_PROGRAM_LINE_RE    = re.compile(r'^\s*program\s', re.IGNORECASE)
_MODULE_DECL_RE     = re.compile(r'^\s*module\s+([a-zA-Z0-9_]+)')
_SUBMODULE_DECL_RE  = re.compile(
    r'^\s*submodule\s*\(\s*([a-z0-9_]+)(?::([a-z0-9_]+))?\s*\)\s*([a-zA-Z0-9_]+)',
    re.IGNORECASE,
)
_INCLUDE_LINE_RE    = re.compile(
    # The include name may be a hierarchical path (e.g.
    # 'objects/nodes/components/disk/standard/bound_functions.inc'), so allow
    # '/' in addition to word characters, '.' and '-'. Without '/', includes of
    # generated component code in the hierarchical source tree are not followed,
    # so the `use` statements they contain (e.g. of component *_data modules)
    # are missed and the resulting object dependency (.d) lists are incomplete.
    r"""^\s*include\s+(['"])([\w./\-]+)(['"])""",
    re.IGNORECASE,
)


def _update_conditional_compile(stack, preprocessor_directives):
    """Return True when every conditional block on `stack` is active.

    Entries flagged `unknown` (a `#if <expression>` or a chain containing
    `#elif`) are treated as active regardless of `#else` flips — every
    branch of an unevaluable conditional is scanned. Mirrors the foreach
    loop at useDependencies.pl:374-383, with the unknown-entry extension.
    """
    for entry in stack:
        if entry.get('unknown'):
            continue
        name_defined = entry['name'] in preprocessor_directives
        active = entry['state'] if name_defined else 1 - entry['state']
        if active == 0:
            return False
    return True


def _scan_source_file(sources_entry, file_names_to_process, source_file,
                     directives, locations, root_source_dir, work_dir,
                     preprocessor_directives_set):
    """Drain `file_names_to_process`, performing the line-by-line scan that
    populates `sources_entry`.  Mirrors useDependencies.pl:324-461.
    """
    while file_names_to_process:
        full_path = file_names_to_process.pop()

        if re.search(r'\.cpp$', full_path, re.IGNORECASE):
            sources_entry['libraryDependencies']['stdc++'] = True

        preprocessor_stack    = []
        conditionally_compile = True
        in_xml                = False
        in_latex              = False

        try:
            fh = open(full_path, 'r', errors='replace')
        except OSError:
            sys.exit("Can't open input file: " + full_path)
        try:
            for line in fh:
                if _XML_CLOSE_RE.match(line):
                    in_xml = False
                if _LATEX_CLOSE_RE.match(line):
                    in_latex = False
                if in_xml or in_latex:
                    continue

                m = _LIB_HINT_RE.match(line)
                if m:
                    sources_entry['libraryDependencies'][m.group(1)] = True

                m = _C_INCLUDE_RE.match(line)
                if m:
                    include_file = m.group(1).lower()
                    if (include_file in INCLUDE_LIBRARIES
                            and conditionally_compile):
                        sources_entry['libraryDependencies'][
                            INCLUDE_LIBRARIES[include_file]
                        ] = True

                # Preprocessor lines.
                if line.startswith('#'):
                    m = _PP_IFDEF_RE.match(line)
                    if m:
                        preprocessor_stack.append({
                            'name': m.group(1), 'state': 1,
                        })
                    elif _PP_IF_RE.match(line):
                        # A general `#if <expression>` cannot be evaluated
                        # here; treat every branch of it as active so its
                        # `use` statements are still recorded (over-depending
                        # is benign; the previous behavior treated the block
                        # as inactive and silently dropped dependencies).
                        preprocessor_stack.append({
                            'name': None, 'state': 1, 'unknown': True,
                        })
                    else:
                        m = _PP_IFNDEF_RE.match(line)
                        if m:
                            preprocessor_stack.append({
                                'name': m.group(1), 'state': 0,
                            })
                        elif _PP_ENDIF_RE.match(line):
                            if preprocessor_stack:
                                preprocessor_stack.pop()
                        elif _PP_ELIF_RE.match(line):
                            # A chain containing `#elif` can no longer be
                            # tracked by defined-ness flips — degrade the
                            # whole chain to unevaluable (all branches
                            # active).
                            if preprocessor_stack:
                                preprocessor_stack[-1]['unknown'] = True
                        elif _PP_ELSE_RE.match(line):
                            if preprocessor_stack:
                                preprocessor_stack[-1]['state'] = (
                                    1 - preprocessor_stack[-1]['state']
                                )
                    conditionally_compile = _update_conditional_compile(
                        preprocessor_stack, preprocessor_directives_set,
                    )

                if conditionally_compile:
                    m = _USE_LINE_RE.match(line)
                    if m:
                        used = m.group(3).lower()
                        if used in MODULE_LIBRARIES:
                            sources_entry['libraryDependencies'][
                                MODULE_LIBRARIES[used]
                            ] = True
                        if used not in EXTERNAL_MODULES:
                            sources_entry['modulesUsed'].append(
                                work_dir + used + '.mod'
                            )

                    if _OMP_PARALLEL_RE.match(line):
                        sources_entry['modulesUsed'].append(
                            work_dir + 'events_filters.mod'
                        )

                    if _OMP_CRITICAL_RE.match(line):
                        sources_entry['modulesUsed'].append(
                            work_dir + 'openmp_utilities_data.mod'
                        )

                    m = _EXPLICIT_DEPS_RE.match(line)
                    if m:
                        sources_entry['dependenciesExplicit'].extend(
                            m.group(1).split()
                        )

                    if _PROGRAM_LINE_RE.match(line):
                        sources_entry['modulesUsed'].append(
                            work_dir + 'iso_varying_string.mod'
                        )

                    m = _MODULE_DECL_RE.match(line)
                    if m:
                        module_name = m.group(1)
                        if module_name != 'procedure':
                            # Store in the same `<build>/<lower>.mod` form
                            # that `modulesUsed` uses, so the self-reference
                            # filter in _write_makefile can actually match.
                            # The Perl original stored `ModuleName.mod`
                            # here -- case preserved, no path -- which
                            # meant the filter was a no-op; fixed in the
                            # Python port.
                            sources_entry['modulesProvided'][
                                work_dir + module_name.lower() + '.mod'
                            ] = True
                            if directives['functionClass']:
                                for fc in directives['functionClass']:
                                    instance_files = as_array(
                                        (locations.get(fc['name']) or {})
                                            .get('file')
                                        if locations else None,
                                    )
                                    for instance_file in instance_files:
                                        for d in extract_directives(
                                            instance_file, fc['name'],
                                        ):
                                            sources_entry['submodules'].append(
                                                d['name'] + '_'
                                            )

                    m = _SUBMODULE_DECL_RE.match(line)
                    if m:
                        sources_entry['submodulesProvided'].append({
                            'submoduleName':  m.group(3).lower(),
                            'moduleFileName': m.group(1).lower() + '.mod',
                        })

                    m = _INCLUDE_LINE_RE.match(line)
                    if m:
                        pre   = m.group(2)
                        raw   = re.sub(r'\.inc$', '.Inc', pre)
                        raw_path = os.path.join(
                            root_source_dir, 'source', raw,
                        )
                        pre_path = work_dir + pre
                        if os.path.exists(raw_path):
                            file_names_to_process.append(raw_path)
                            sources_entry['files'].append(raw_path)
                        elif os.path.exists(pre_path):
                            file_names_to_process.append(pre_path)
                            sources_entry['files'].append(pre_path)
                        else:
                            m2 = re.match(r'(.*)\.type\.inc', pre)
                            if m2:
                                extra = as_array(
                                    (locations.get(m2.group(1)) or {})
                                        .get('file')
                                    if locations else None,
                                )
                                file_names_to_process.extend(extra)
                                sources_entry['files'].extend(extra)

                if _XML_OPEN_RE.match(line):
                    in_xml = True
                if _LATEX_OPEN_RE.match(line):
                    in_latex = True
        finally:
            fh.close()


# ---------------------------------------------------------------------------
# Source-file discovery
# ---------------------------------------------------------------------------

def _source_files_to_process(root_source_dir, build_path):
    """Return the ordered list of source-file descriptors to scan.

    Walks `<root>/source`, all of its subdirectories (recursively), and
    `$BUILDPATH/libgalacticus`.  Each descriptor is a dict with
    `fileName`, `fullPathFileName`, and `subDirectoryName` (relative to
    `<root>/source` or `$BUILDPATH`).
    """
    root_source = os.path.join(root_source_dir, 'source')
    directories = []
    if os.path.exists(root_source):
        for dirpath, dirnames, _ in os.walk(root_source):
            dirnames[:] = sorted(d for d in dirnames if not d.startswith('.'))
            directories.append(dirpath)
    if not directories:
        directories = [root_source]
    # The libgalacticus subdirectory contains auto-generated Fortran wrappers produced by libraryInterfaces.py; it only exists
    # for library builds. Skip it gracefully when absent so that non-library builds don't fail here.
    libgalacticus_dir = os.path.join(build_path, 'libgalacticus')
    if os.path.isdir(libgalacticus_dir):
        directories.append(libgalacticus_dir)

    # Perl strips `<root>/source` or `$BUILDPATH` prefixes to derive the
    # per-file `subDirectoryName` field.  Matching regex:
    #   s/^(<root>/source|<BUILDPATH>)/?//
    prefix_re = re.compile(
        r'^(' + re.escape(root_source) + r'|' + re.escape(build_path) + r')/?'
    )

    results = []
    for directory in directories:
        try:
            names = os.listdir(directory)
        except OSError:
            sys.exit("useDependencies.py: can not open the source directory: "
                     + directory)
        sub_dir = prefix_re.sub('', directory)
        for name in sorted(names):
            if not _SOURCE_SUFFIX_RE.search(name) or name.startswith('.#'):
                continue
            results.append({
                'fileName':         name,
                'fullPathFileName': os.path.join(directory, name),
                'subDirectoryName': sub_dir,
            })
    return results


# ---------------------------------------------------------------------------
# Aggregation + submodule map
# ---------------------------------------------------------------------------

def _finalise_event_hooks_manager(uses_per_file, event_hook_modules):
    """Inject accumulated event-hook module deps into the file that carries
    the `eventHookManager` directive (if any).  Mirrors useDependencies.pl:
    464-470.
    """
    manager = uses_per_file.get('eventHooksManager')
    if not manager:
        return
    modules = sorted(set(event_hook_modules))
    manager['modules'] = modules
    fid = manager.get('fileIdentifier')
    if not fid or fid not in uses_per_file:
        return
    entry = uses_per_file[fid]
    entry.setdefault('modulesUsed', [])
    entry['modulesUsed'].extend(modules)
    entry['modulesUsed'] = sorted(set(entry['modulesUsed']))


def _build_submodule_map(uses_per_file, work_dir):
    """Return `{work_dir+lc(mod): [submodule names]}` covering both
    functionClass-synthesised submodules and Fortran `submodule (…)`
    statements.  Mirrors useDependencies.pl:472-489.
    """
    submodules = {}
    for file_id, entry in uses_per_file.items():
        if file_id == 'eventHooksManager':
            continue
        subs = entry.get('submodules') or []
        if not subs:
            continue
        provided = list((entry.get('modulesProvided') or {}).keys())
        if len(provided) != 1:
            raise RuntimeError(
                f"useDependencies.py: submodules [{', '.join(subs)}] "
                f"associated with multiple modules [{', '.join(provided)}]"
            )
        # `modulesProvided` keys are already stored in `<build>/<lower>.mod`
        # form, so use directly (this is the format the submodule-map
        # consumers look up against `modulesUsed`).
        submodules[provided[0]] = list(subs)

    for file_id, entry in uses_per_file.items():
        if file_id == 'eventHooksManager':
            continue
        for sub in entry.get('submodulesProvided') or []:
            key = work_dir + sub['moduleFileName']
            submodules.setdefault(key, []).append(sub['submoduleName'])
    return submodules


# ---------------------------------------------------------------------------
# Makefile emission
# ---------------------------------------------------------------------------

def _write_makefile(path, uses_per_file, source_files, submodules,
                    work_dir):
    """Emit `Makefile_Use_Dependencies`.  Mirrors useDependencies.pl:
    491-557.
    """
    with open(path, 'w') as mk:
        for sf in source_files:
            file_identifier = _file_identifier(sf['fullPathFileName'])
            entry = uses_per_file.get(file_identifier)
            if not entry:
                continue

            object_file_name = re.sub(
                r'\.(f|f90|c|cpp)$', '.o', sf['fileName'],
                flags=re.IGNORECASE,
            )
            dependency_file_name = re.sub(
                r'\.(f|f90|c|cpp|inc)$', '.d', sf['fileName'],
                flags=re.IGNORECASE,
            )
            sub_dir  = sf['subDirectoryName']
            work_sub = work_dir + (sub_dir + '/' if sub_dir else '')

            # Library sidecar (.fl): write its content directly, only-if-changed
            # (preserving the mtime -- and so sparing dependent .o targets -- when
            # the library set is unchanged), and also emit the traditional echo
            # rule. The emitted rules have no prerequisites, so once a .fl exists
            # its recipe can never run again: without the direct write here, a
            # change in a source file's library dependencies (or any corruption of
            # an existing .fl) would never be repaired short of a clean build. The
            # echo rule is kept as the recovery path for a manually-deleted .fl.
            if not object_file_name.endswith('.Inc'):
                library_file_name = re.sub(r'\.o$', '.fl', object_file_name)
                fl_path = work_sub + library_file_name
                fl_dir = os.path.dirname(fl_path)
                if fl_dir:
                    os.makedirs(fl_dir, exist_ok=True)
                fl_tmp = fl_path + '.tmp'
                with open(fl_tmp, 'w') as fl:
                    for lib in sorted((entry.get('libraryDependencies') or {}).keys()):
                        fl.write(lib + '\n')
                file_changes_update(fl_path, fl_tmp)
                mk.write(f"{fl_path}:\n")
                redirect = '>'
                for lib in sorted((entry.get('libraryDependencies') or {}).keys()):
                    mk.write(f"\t@echo {lib} {redirect} {fl_path}\n")
                    redirect = '>>'

            modules_used  = entry.get('modulesUsed') or []
            explicit_deps = entry.get('dependenciesExplicit') or []
            if not (modules_used or explicit_deps):
                continue

            # Filter out self-provided modules, sort+uniq.
            provided = entry.get('modulesProvided') or {}
            modules_used = sorted({m for m in modules_used if m not in provided})
            entry['modulesUsed'] = modules_used

            # Any submodule entries associated with modules we use land
            # here as `<mod>@<sub>.smod` tokens.
            submodules_used = []
            for module_used in modules_used:
                subs = submodules.get(module_used)
                if not subs:
                    continue
                module_name = re.sub(r'\.mod$', '', module_used)
                for sub in subs:
                    submodules_used.append(f"{module_name}@{sub.lower()}.smod")

            obj_path = work_sub + object_file_name
            dep_path = work_sub + dependency_file_name

            # Object-file rule.
            deps_list = (
                [work_dir + 'utility/OpenMP/workaround.o']
                + modules_used + submodules_used + explicit_deps
            )
            mk.write(
                f"{obj_path}: " + ' '.join(deps_list) + ' Makefile\n'
            )

            # Dependency-file (`.d`) aggregation rule.
            dependencies_used = [
                m + '.d' if m.endswith('.mod') else m for m in modules_used
            ]
            explicit_d = [
                re.sub(r'\.o$', '.d', d) for d in explicit_deps
            ]
            mk.write(
                f"{dep_path}: " +
                ' '.join(dependencies_used + explicit_d) + '\n'
            )
            mk.write(f"\t@echo {obj_path} > {dep_path}~\n")
            mk.write(
                f"\t@echo {work_dir}utility/OpenMP/workaround.o >> "
                f"{dep_path}~\n"
            )
            for dep_explicit in explicit_deps:
                dep_d = re.sub(r'\.o$', '.d', dep_explicit)
                prefix = '' if '/' in dep_explicit else work_dir
                mk.write(f"\t@cat {prefix}{dep_d} >> {dep_path}~\n")
            for dep in dependencies_used:
                mk.write(f"\t@cat {dep} >> {dep_path}~\n")
            mk.write(f"\t@sort -u {dep_path}~ -o {dep_path}~\n")
            mk.write(f"\t@if cmp -s {dep_path} {dep_path}~ ; then \\\n")
            mk.write(f"\t rm {dep_path}~ ; \\\n")
            mk.write("\telse \\\n")
            mk.write(f"\t mv {dep_path}~ {dep_path} ; \\\n")
            mk.write("\tfi\n\n")

            # GraphViz rule.
            graphvizes_used = [
                m + '.gv' if m.endswith('.d') else m for m in modules_used
            ]
            gv_path = work_sub + sf['fileName'] + '.gv'
            mk.write(f"{gv_path}: {dep_path} " + ' '.join(graphvizes_used)
                     + '\n')
            label = sub_dir + sf['fileName']
            mk.write(f"\t@echo \\\"{label}\\\" > {gv_path}\n")
            for dep_explicit in explicit_deps:
                dep_d = re.sub(r'\.o$', '.d', dep_explicit)
                prefix = '' if '/' in dep_explicit else work_dir
                mk.write(
                    "\t@awk '{print \"\\\"" + label +
                    "\\\" -> \\\"\"$$1\"\\\"\"}' " +
                    f"{prefix}{dep_d} >> {gv_path}\n"
                )
            for gv_used in graphvizes_used:
                mk.write(
                    "\t@awk '{print \"\\\"" + label +
                    "\\\" -> \\\"\"$$1\"\\\"\"}' " +
                    f"{gv_used} >> {gv_path}\n"
                )
                mk.write(
                    "\t@cat `awk '{print \"" + work_dir +
                    "\"$$1\".gv\"}' " + gv_used + f"` >> {gv_path}\n"
                )
            mk.write(f"\t@sort -u {gv_path} -o {gv_path}\n\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv):
    if len(argv) != 2:
        print("Usage: useDependencies.py <sourceDirectory>", file=sys.stderr)
        sys.exit(1)
    root_source_dir = argv[1]
    build_path      = os.environ['BUILDPATH']
    work_dir        = build_path + '/'

    locations       = _load_xml(work_dir + 'directiveLocations.xml')
    state_storables = _load_xml(work_dir + 'stateStorables.xml')

    preprocessor_directives = _collect_preprocessor_directives(build_path)

    uses_per_file, cache_mtime = _load_cache(
        work_dir + 'Makefile_Use_Dependencies.blob',
    )
    # Discard caches written under older scan rules (the conditional-
    # compilation handling changed what a scan records); a mismatch costs
    # one full rescan.
    if uses_per_file.pop('blobVersion', None) != _BLOB_VERSION:
        uses_per_file, cache_mtime = {}, None
    have_per_file = cache_mtime is not None

    source_files = _source_files_to_process(root_source_dir, build_path)

    # Cache-invalidation: add/remove detection (and the stale-entry prune
    # below). Identifiers are canonicalized via `file_identifier()` so they
    # use the same key form as the cache entries.
    current_ids = [
        _file_identifier(sf['fullPathFileName']) for sf in source_files
    ]
    current_id_set = set(current_ids)
    force_rescan = False
    if have_per_file:
        if any(fid not in uses_per_file for fid in current_ids):
            force_rescan = True
        if any(fid not in current_id_set
               for fid in uses_per_file
               if fid != 'eventHooksManager'):
            force_rescan = True

    preprocessor_set = frozenset(preprocessor_directives)
    event_hook_modules = []

    # Decide which files need a (re)scan, preserving `source_files` order so the
    # merge below reproduces the serial accumulation order exactly.
    scan_list = []
    for sf in source_files:
        file_identifier = _file_identifier(sf['fullPathFileName'])
        rescan = True
        if have_per_file and file_identifier in uses_per_file:
            tracked = uses_per_file[file_identifier].get('files') or []
            stale = any(
                os.path.exists(t) and os.stat(t).st_mtime > cache_mtime
                for t in tracked
            )
            rescan = bool(stale)
        if not (rescan or force_rescan):
            continue
        scan_list.append((file_identifier, sf))

    # Scan the files concurrently, then merge results in scan-list order. The
    # per-file `entry`, the ordered `event_hook_modules` pool, and the singleton
    # `eventHooksManager` record are reconstructed exactly as a serial run would.
    _WORKER.update({
        'locations':        locations,
        'state_storables':  state_storables,
        'work_dir':         work_dir,
        'root_source_dir':  root_source_dir,
        'preprocessor_set': preprocessor_set,
    })
    for file_identifier, entry, file_event_hook_modules, manager in parallel_scan(
            scan_list, _scan_one, 'useDependencies.py'):
        uses_per_file.pop(file_identifier, None)
        uses_per_file[file_identifier] = entry
        event_hook_modules.extend(file_event_hook_modules)
        if 'eventHooksManager' in manager:
            uses_per_file['eventHooksManager'] = manager['eventHooksManager']

    # Drop cache entries for source files that no longer exist, or their
    # `use` and submodule rules would be re-emitted into
    # Makefile_Use_Dependencies forever (phantom dependencies on deleted
    # files). The `eventHooksManager` key is bookkeeping, not a per-file
    # entry, so it is preserved.
    _prune_cache(uses_per_file, current_id_set,
                 reserved={'eventHooksManager'})

    _finalise_event_hooks_manager(uses_per_file, event_hook_modules)

    submodule_map = _build_submodule_map(uses_per_file, work_dir)

    _write_makefile(
        work_dir + 'Makefile_Use_Dependencies',
        uses_per_file, source_files, submodule_map, work_dir,
    )

    with open(work_dir + 'Makefile_Use_Dependencies.blob', 'wb') as fh:
        uses_per_file['blobVersion'] = _BLOB_VERSION
        pickle.dump(uses_per_file, fh, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    main(sys.argv)
