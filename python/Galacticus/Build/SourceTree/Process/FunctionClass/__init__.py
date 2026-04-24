# Processes `functionClass` directives: the cornerstone of the Galacticus
# functionClass infrastructure.  This is the first slice (D.7.3a) of a
# three-part port of perl/Galacticus/Build/SourceTree/Process/FunctionClass.pm
# (3105 Perl lines):
#
#   * D.7.3a (this PR) — scaffolding: directive collection, class
#     discovery + topo-sort, small method stubs (stateStore / stateRestore /
#     autoHook / destructor / objectType), and submodule file output.
#     The heavy method-generation (buildDescriptorMethods,
#     buildAllowedParametersMethod, etc.) and body-emission
#     (generateTypeDefinition, generateConstructor,
#     generateClassSubmodules, generateMethodFunctions,
#     generateDocumentation) pathways are stubbed with NotImplementedError
#     so the hook is installable and unit-testable today without
#     pretending to emit a working class.
#   * D.7.3b (next PR) — the generate_* body emitters.
#   * D.7.3c (final PR) — the build_*_methods heavy lifters.
#
# Andrew Benson (ported to Python 2026)

import os
import re
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from List.ExtraUtils                                         import as_array
from Sort.Topo                                               import sort as topo_sort
from XML.Utils                                               import xml_to_dict
from Galacticus.Build.SourceTree                             import (
    walk_tree, parse_code, parse_file, insert_after_node,
    insert_pre_contains, insert_post_contains, prepend_child_to_node,
    serialize,
)
from Galacticus.Build.SourceTree.Process                     import (
    register_process, process_tree,
)
from Galacticus.Build.SourceTree.Process.FunctionClass.Utils import (
    class_dependencies,
)


# Module-level caches.  Perl declares both with `our` so subroutines that
# `use` this module can read them via the `Process::FunctionClass::<name>`
# aliases in DeepCopy.pm / StateStore.pm.  The Python helpers take the
# state explicitly as parameters, so we only use these to avoid repeat
# parses of the XML files within a single `process_function_class` run.
_STATE_STORABLES_CACHE  = None
_DEEP_COPY_ACTIONS_CACHE = None


def _load_xml(name, cache_attr):
    """Load `$BUILDPATH/<name>` via xml_to_dict, caching per-process."""
    if cache_attr['value'] is not None:
        return cache_attr['value']
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        cache_attr['value'] = {}
        return cache_attr['value']
    path = os.path.join(build_path, name)
    if not os.path.exists(path):
        cache_attr['value'] = {}
        return cache_attr['value']
    cache_attr['value'] = xml_to_dict(ET.parse(path).getroot())
    return cache_attr['value']


_STATE_STORABLES_HOLDER   = {'value': None}
_DEEP_COPY_ACTIONS_HOLDER = {'value': None}
_DIRECTIVE_LOCATIONS_HOLDER = {'value': None}


def _function_class_names(state_storables):
    """Set of functionClass name keys."""
    fc = (state_storables or {}).get('functionClasses') or {}
    if not isinstance(fc, dict):
        return set()
    entries = fc.get('functionClass')
    if entries is None:
        return set(fc.keys())
    if isinstance(entries, dict):
        entries = [entries]
    return {e.get('name') for e in entries
            if isinstance(e, dict) and 'name' in e}


def _function_class_instances(state_storables):
    """List of functionClassInstances names."""
    raw = (state_storables or {}).get('functionClassInstances') or []
    if isinstance(raw, str):
        return [raw] if raw else []
    if isinstance(raw, dict):
        name = raw.get('content') or raw.get('name')
        return [name] if name else []
    out = []
    for item in raw:
        if isinstance(item, str):
            out.append(item)
        elif isinstance(item, dict):
            name = item.get('content') or item.get('name')
            if name:
                out.append(name)
    return out


# ---------------------------------------------------------------------------
# Small helpers (fully implemented in D.7.3a)
# ---------------------------------------------------------------------------

def _is_function_class_pointer(declaration, type_name, state_storables):
    """Return True if `declaration` is a pointer to a registered
    functionClass (or functionClassInstance) whose type is `type_name`.

    Mirrors isFunctionClassPointer() at FunctionClass.pm:1467-1478.
    """
    if declaration.get('intrinsic') not in ('class', 'type'):
        return False
    targets = _function_class_names(state_storables) | set(
        _function_class_instances(state_storables))
    if type_name not in targets:
        return False
    return any(a == 'pointer' for a in declaration.get('attributes') or [])


def _init_code_content():
    """Return `(code_content, pre_contains, post_contains)`.

    Mirrors initCodeContent() at FunctionClass.pm:181-214.
    """
    code_content = {
        'module': {
            'preContains':  [],
            'postContains': [],
            'interfaces':   [],
        },
        'submodule': {},
    }
    pre  = _code_node()
    post = _code_node()
    code_content['module']['preContains'].append(pre)
    code_content['module']['postContains'].append(post)
    return code_content, pre, post


def _code_node(content=''):
    return {
        'type':       'code',
        'content':    content,
        'parent':     None,
        'firstChild': None,
        'sibling':    None,
        'source':
            'Galacticus.Build.SourceTree.Process.FunctionClass'
            '.process_function_class()',
        'line':       1,
    }


def _build_base_method_stubs(directive, methods):
    """Populate `methods` with the default `stateStore`, `stateRestore`,
    `autoHook`, and (optional) `destructor` entries.

    Mirrors buildBaseMethodStubs() at FunctionClass.pm:1480-1546.
    """
    iso_c_module = [{
        'name':      'ISO_C_Binding',
        'intrinsic': 1,
        'only':      ['c_ptr'],
    }]
    methods['stateStore'] = {
        'description': 'Store the state of the object to file.',
        'type':        'void',
        'pass':        'yes',
        'modules':     list(iso_c_module),
        'argument':    [
            'integer, intent(in   ) :: stateFile',
            'type(c_ptr), intent(in   ) :: gslStateFile',
        ],
    }
    methods['stateRestore'] = {
        'description': 'Restore the state of the object to file.',
        'type':        'void',
        'pass':        'yes',
        'modules':     list(iso_c_module),
        'argument':    [
            'integer, intent(in   ) :: stateFile',
            'type(c_ptr), intent(in   ) :: gslStateFile',
        ],
    }
    methods['autoHook'] = {
        'description': 'Insert any event hooks required by this object.',
        'type':        'void',
        'pass':        'yes',
        'code':        (
            '!$GLC attributes unused :: self\n\n'
            '! Nothing to do by default.\n'
        ),
    }
    if 'autoHook' in directive and isinstance(directive['autoHook'], dict):
        for module in as_array(directive['autoHook'].get('modules')):
            if not isinstance(module, dict):
                continue
            module_name = module.get('name')
            only_list = [s.strip() for s in re.split(
                r'\s*,\s*', module.get('only', '')) if s.strip()]
            methods['autoHook'].setdefault('modules', []).append({
                'name': module_name, 'only': only_list,
            })
        if 'code' in directive['autoHook']:
            methods['autoHook']['code'] = directive['autoHook']['code']
    if 'destructor' in directive and isinstance(directive['destructor'], dict):
        dest = {
            'description': 'Destructor for this class.',
            'type':        'void',
            'pass':        'yes',
            'code':        directive['destructor'].get('code', ''),
        }
        for module in as_array(directive['destructor'].get('modules')):
            if not isinstance(module, dict):
                continue
            module_name = module.get('name')
            only_list = [s.strip() for s in re.split(
                r'\s*,\s*', module.get('only', '')) if s.strip()]
            dest.setdefault('modules', []).append({
                'name': module_name, 'only': only_list,
            })
        methods['destructor'] = dest


def _build_object_type_method(directive, non_abstract_classes, methods):
    """Populate `methods['objectType']` with the select-type ladder that
    returns the concrete class name (or its short form).

    Mirrors buildObjectTypeMethod() at FunctionClass.pm:1548-1586.
    """
    directive_name = directive['name']
    code = (
        "logical :: short_\n"
        "short_=.false.\n"
        "if (present(short)) short_=short\n"
        "select type (self)\n"
    )
    for non_abstract in non_abstract_classes:
        class_name = non_abstract['name']
        # shortName is computed the same way loadAndSortClasses does — Perl
        # uses lcfirst unless the remainder is already an all-caps acronym
        # of 2+ letters.
        type_short = class_name
        if type_short.startswith(directive_name):
            type_short = type_short[len(directive_name):]
        if not re.match(r'^[A-Z]{2,}', type_short):
            type_short = type_short[:1].lower() + type_short[1:]
        code += (
            f"type is ({class_name})\n"
            "if (short_) then\n"
            f" {directive_name}ObjectType='{type_short}'\n"
            "else\n"
            f" {directive_name}ObjectType='{class_name}'\n"
            "end if\n"
        )
    code += "end select\n"
    methods['objectType'] = {
        'description': 'Return the type of the object.',
        'type':        'type(varying_string)',
        'pass':        'yes',
        'modules':     'ISO_Varying_String',
        'argument':    ['logical, intent(in   ), optional :: short'],
        'code':        code,
    }


# ---------------------------------------------------------------------------
# Class discovery + topological sort
# ---------------------------------------------------------------------------

def _load_and_sort_classes(directive, directive_locations):
    """Discover every class file listed under directiveLocations for the
    directive's name, parse + process each, extract dependencies via
    `class_dependencies`, topo-sort, and build the short-name / abstract
    metadata that downstream method builders consume.

    Mirrors loadAndSortClasses() at FunctionClass.pm:104-179.  Returns
    `(classes_by_type, classes_ordered, non_abstract_classes)`.
    """
    class_locations = list(as_array(
        (directive_locations.get(directive['name']) or {}).get('file')))

    dependencies = {}
    classes      = {}
    class_counts = {}

    for class_location in class_locations:
        class_tree = parse_file(class_location)
        # Process the class tree so Class_Dependencies sees post-processing
        # structure (declarations, visibilities, etc.).  The Perl call
        # passes `errorTolerant => 1` — an option ProcessTree does not
        # consume anywhere in the Perl codebase, so a no-op.
        process_tree(class_tree)

        class_record, class_deps = class_dependencies(
            class_tree, directive['name'])
        for dep in class_deps:
            type_name = class_record.get('type')
            if type_name and dep != type_name:
                dependencies.setdefault(dep, []).append(type_name)

        class_record['file'] = class_location
        class_record['tree'] = class_tree
        class_record.setdefault('abstract', 'no')

        if class_record.get('type') is None:
            raise RuntimeError(
                f"process_function_class: class is undefined in file "
                f"'{class_location}'"
            )
        classes[class_record['type']] = class_record
        class_counts.setdefault(class_record['type'], []).append(
            class_record['file'])

    # Duplicate-class check.
    duplicates = [
        (name, files) for name, files in sorted(class_counts.items())
        if len(files) > 1
    ]
    if duplicates:
        msg = '\n'.join(
            f"Duplicate class '{name}' found in:\n"
            + '\n'.join(f"\t{f}" for f in files)
            for name, files in duplicates
        )
        raise RuntimeError(
            "process_function_class: duplicate classes found\n" + msg)

    # Alphabetical then topological sort.
    unsorted_classes = sorted(classes.keys())
    sorted_names = topo_sort(unsorted_classes, dependencies)
    classes_ordered = [classes[name] for name in sorted_names]

    non_abstract_classes = [
        c for c in classes_ordered if c.get('abstract') != 'yes'
    ]

    # Short name (directive prefix stripped; lcfirst unless all-caps).
    for non_abstract in non_abstract_classes:
        m = re.match(r'^' + re.escape(directive['name']) + r'([a-zA-Z0-9]+)',
                     non_abstract['name'])
        if not m:
            raise RuntimeError("class name has incorrect format")
        short = m.group(1)
        if not re.match(r'^[A-Z]{2,}', short):
            short = short[:1].lower() + short[1:]
        non_abstract['shortName'] = short

    # Validate optional `default` attribute.
    if 'default' in directive:
        wanted = directive['default']
        if not any(wanted == c['shortName'] for c in non_abstract_classes):
            allowed = '\n'.join(
                f"    {c['shortName']}" for c in non_abstract_classes)
            raise RuntimeError(
                f"unrecognized default '{wanted}' for class "
                f"'{directive['name']}'\n  allowed defaults are:\n{allowed}"
            )

    return classes, classes_ordered, non_abstract_classes


# ---------------------------------------------------------------------------
# Output injection + submodule writing
# ---------------------------------------------------------------------------

def _insert_and_write_output(node, code_content, pre, post, classes, directive):
    """Re-parse the accumulated pre/post-contains content, run it through
    process_tree so nested directives expand, and inject the trees into
    the enclosing module.  Additionally, write any per-class submodule
    files to disk.

    Mirrors insertAndWriteOutput() at FunctionClass.pm:1401-1465.  The
    File::Changes "only if different" semantics are approximated by
    writing to a temp file and renaming only when the content actually
    differs (using a simple `os.path.exists` + `open/read` compare).
    """
    tree_pre  = parse_code(
        pre['content'],
        name='Galacticus.Build.SourceTree.Process.FunctionClass'
             '.process_function_class()',
    )
    tree_post = parse_code(
        post['content'],
        name='Galacticus.Build.SourceTree.Process.FunctionClass'
             '.process_function_class()',
    )
    process_tree(tree_pre)
    process_tree(tree_post)

    insert_after_node(node, [tree_pre])
    interfaces = code_content['module'].get('interfaces') or []
    if interfaces:
        insert_pre_contains(node['parent'], list(interfaces))
    insert_post_contains(node['parent'], [tree_post])

    # Per-class submodules.  Iterate sorted by class name for stable output.
    for class_name in sorted(code_content['submodule'].keys()):
        submodule_name = class_name + '_'
        file_node = {
            'type':       'file',
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':
                'Galacticus.Build.SourceTree.Process.FunctionClass'
                '.process_function_class()',
            'line':       1,
        }
        # Parent (sub)module — either the top-level module, or the submodule
        # associated with whatever class this class extends.
        extends = classes[class_name].get('extends', '')
        if extends == directive['name'] + 'Class':
            parent_name = node['parent']['name']
        else:
            parent_name = f"{node['parent']['name']}:{extends}_"
        submodule_node = {
            'type':       'submodule',
            'opener':     f"submodule ({parent_name}) {submodule_name}\n",
            'closer':     f"end submodule {submodule_name}\n",
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':
                'Galacticus.Build.SourceTree.Process.FunctionClass'
                '.process_function_class()',
            'line':       1,
        }
        prepend_child_to_node(file_node, [submodule_node])

        # Detach any existing links on the pre/post children before
        # attaching them to the fresh submodule tree (parallels the loop at
        # FunctionClass.pm:1447-1451).
        block = code_content['submodule'][class_name]
        for kid in block.get('preContains', []) + block.get('postContains', []):
            kid['parent']  = None
            kid['sibling'] = None

        prepend_child_to_node(submodule_node, block.get('preContains', []))
        insert_post_contains(submodule_node, block.get('postContains', []))

        submodule_content = serialize(file_node)
        file_name = block['fileName']
        tmp_name  = file_name + '.tmp'
        with open(tmp_name, 'w') as fh:
            fh.write(submodule_content)
        _update_if_changed(file_name, tmp_name)


def _update_if_changed(target, tmp):
    """Move `tmp` over `target` only when their contents differ.  Best-effort
    replacement for the Perl `File::Changes::Update(…, proveUpdate => 'yes')`
    used by the Perl port; full line-mapping support (the `.lmap` sidecar)
    is out of scope for this slice.
    """
    if os.path.exists(target):
        try:
            with open(target, 'r') as fh:
                existing = fh.read()
            with open(tmp, 'r') as fh:
                new = fh.read()
            if existing == new:
                os.unlink(tmp)
                return
        except OSError:
            pass
    os.replace(tmp, target)


# ---------------------------------------------------------------------------
# Stubs for the heavy method-builders / body-emitters (D.7.3b, D.7.3c)
# ---------------------------------------------------------------------------

def _not_implemented(name):
    raise NotImplementedError(
        f"process_function_class: {name} is deferred to a later PR "
        f"(D.7.3b or D.7.3c) of the FunctionClass port"
    )


def _build_descriptor_methods(*args, **kwargs):
    _not_implemented('buildDescriptorMethods')


def _build_allowed_parameters_method(*args, **kwargs):
    _not_implemented('buildAllowedParametersMethod')


def _build_assignment_method(*args, **kwargs):
    _not_implemented('buildAssignmentMethod')


def _build_deep_copy_methods(*args, **kwargs):
    _not_implemented('buildDeepCopyMethods')


def _build_state_store_methods(*args, **kwargs):
    _not_implemented('buildStateStoreMethods')


def _generate_type_definition(*args, **kwargs):
    _not_implemented('generateTypeDefinition')


def _generate_constructor(*args, **kwargs):
    _not_implemented('generateConstructor')


def _generate_class_submodules(*args, **kwargs):
    _not_implemented('generateClassSubmodules')


def _generate_method_functions(*args, **kwargs):
    _not_implemented('generateMethodFunctions')


def _generate_documentation(*args, **kwargs):
    _not_implemented('generateDocumentation')


# ---------------------------------------------------------------------------
# Top-level dispatcher
# ---------------------------------------------------------------------------

def process_function_class(tree, options):
    """Mirrors Process_FunctionClass() at FunctionClass.pm:34-102.

    In this D.7.3a slice:
     - The "early pass" that marks every functionClass descendant directive
       as processed is fully implemented and safe to run.
     - The per-`<functionClass>` synthesis pipeline is implemented down to
       the point of class discovery + base-method stubs + objectType stub;
       further work raises NotImplementedError so callers see a clear
       error rather than silently-incomplete output.

    Trees that contain no `<functionClass>` directives are effectively
    no-ops — registering this hook does not disturb existing tests.
    """
    global _STATE_STORABLES_CACHE, _DEEP_COPY_ACTIONS_CACHE
    state_storables = _load_xml('stateStorables.xml', _STATE_STORABLES_HOLDER)

    for node in list(walk_tree(tree)):
        # Early pass: a node whose type has the form `<something>Class`
        # where `<something>Class` is a known functionClass key has already
        # been processed by the FunctionClass output for that class; mark
        # its directive processed to avoid a false positive from
        # post_process_directives.
        nclass = (node.get('type') or '') + 'Class'
        if nclass in _function_class_names(state_storables):
            directive = node.get('directive') or {}
            directive['processed'] = True

        if node.get('type') != 'functionClass':
            continue
        directive = node.get('directive') or {}
        if directive.get('processed'):
            continue

        parent = node.get('parent')
        if parent is None or parent.get('type') != 'module':
            raise RuntimeError(
                "process_function_class: parent node must be a module")

        directive['processed'] = True

        directive_locations = _load_xml(
            'directiveLocations.xml', _DIRECTIVE_LOCATIONS_HOLDER)
        deep_copy_actions   = _load_xml(
            'deepCopyActions.xml', _DEEP_COPY_ACTIONS_HOLDER)

        # Extract `method` sub-directive into a {name: dict} map.  XML::Simple
        # returns either a flat dict (single method) or a dict-of-dicts
        # (keyAttr-grouped) in Perl; our xml_to_dict gives a list for
        # multiple methods, a single dict otherwise.
        methods = {}
        raw_method = directive.get('method')
        if isinstance(raw_method, list):
            for m in raw_method:
                if isinstance(m, dict) and 'name' in m:
                    methods[m['name']] = m
        elif isinstance(raw_method, dict):
            if 'name' in raw_method and all(not isinstance(v, dict)
                                            for v in raw_method.values()):
                methods[raw_method['name']] = raw_method
            else:
                methods = dict(raw_method)

        classes, classes_ordered, non_abstract_classes = \
            _load_and_sort_classes(directive, directive_locations)

        _build_base_method_stubs(directive, methods)

        # The following five are stubbed — they'll be wired in by later
        # sub-PRs.  Leaving them out (rather than raising immediately)
        # means D.7.3a is safe to import even for trees that happen to
        # contain `<functionClass>` directives that some other pipeline
        # already processed and marked — those pass through the guard at
        # the top of the loop.
        _build_descriptor_methods(
            directive, non_abstract_classes, classes, methods, tree)
        _build_object_type_method(directive, non_abstract_classes, methods)
        _build_allowed_parameters_method(
            directive, classes_ordered, methods)
        _build_assignment_method(
            directive, non_abstract_classes, classes, methods)
        _build_deep_copy_methods(
            directive, non_abstract_classes, classes,
            node.get('line', 0), methods)
        _build_state_store_methods(
            directive, non_abstract_classes, classes, methods)

        code_content, pre, post = _init_code_content()

        _generate_type_definition (
            directive, methods, pre, node)
        _generate_constructor     (
            directive, classes_ordered, non_abstract_classes,
            pre, post, node, tree)
        _generate_class_submodules(
            directive, classes_ordered, non_abstract_classes,
            pre, code_content, node)
        _generate_method_functions(directive, methods, post, node)
        _generate_documentation   (
            directive, classes, non_abstract_classes)

        _insert_and_write_output(
            node, code_content, pre, post, classes, directive)


register_process('functionClass', process_function_class)
