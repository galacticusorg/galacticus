"""Processes `functionClass` directives: the cornerstone of the Galacticus
functionClass infrastructure.  This is the first slice (D.7.3a) of a
three-part port:

  * D.7.3a (this PR) — scaffolding: directive collection, class
    discovery + topo-sort, small method stubs (stateStore / stateRestore /
    autoHook / destructor / objectType), and submodule file output.
    The heavy method-generation (buildDescriptorMethods,
    buildAllowedParametersMethod, etc.) and body-emission
    (generateTypeDefinition, generateConstructor,
    generateClassSubmodules, generateMethodFunctions,
    generateDocumentation) pathways are stubbed with NotImplementedError
    so the hook is installable and unit-testable today without
    pretending to emit a working class.
  * D.7.3b (next PR) — the generate_* body emitters.
  * D.7.3c (final PR) — the build_*_methods heavy lifters.

Andrew Benson (ported to Python 2026)
"""

import os
import re
import xml.etree.ElementTree as ET


from List.ExtraUtils                                         import as_array
from Sort.Topo                                               import sort as topo_sort
from XML.Utils                                               import xml_to_dict
from Galacticus.Build.StateStorables                         import (
    function_class_names    as _shared_function_class_names,
    function_class_instances as _shared_function_class_instances,
)
from Galacticus.Build.SourceTree                             import (
    walk_tree, parse_code, parse_file, insert_after_node,
    insert_pre_contains, insert_post_contains, prepend_child_to_node,
    serialize, set_visibility,
)
from Galacticus.Build.SourceTree.Process                     import (
    register_process, process_tree,
)
from Galacticus.Build.SourceTree.Parse.Declarations          import (
    parse_declaration, build_declarations,
)
from Galacticus.Build.SourceTree.Parse.ModuleUses            import add_uses, _as_entry_list
from Galacticus.Build.SourceTree.Parse.Visibilities          import update_visibilities
from Galacticus.Build.SourceTree.Process.FunctionClass.Utils import (
    class_dependencies,
)
from Galacticus.Build.SourceTree.Process.SourceIntrospection import location
from Galacticus.Build.FileChanges                            import update as file_changes_update


# ---------------------------------------------------------------------------
# Small formatting helpers
# ---------------------------------------------------------------------------

def _format_variable_definitions(declarations):
    """Render a list of declaration dicts into Fortran declaration lines.

    The output is deliberately simple and unaligned — we only need
    functional equivalence here, and a simple
        `<intrinsic>(<type>), <attr>, <attr> :: <var>, <var>\\n`
    per declaration is enough to keep the downstream Fortran parser /
    compiler happy.  Matches the subset of features
    generateClassSubmodules actually passes through
    (intrinsic, type, attributes, variables).
    """
    out = ''
    for decl in declarations:
        if not isinstance(decl, dict):
            continue
        intrinsic  = decl.get('intrinsic') or ''
        type_text  = decl.get('type')
        attributes = decl.get('attributes') or []
        variables  = decl.get('variables') or []
        line = intrinsic
        if type_text is not None:
            # type in declarations is stored with parens stripped; re-wrap.
            has_parens = type_text.startswith('(') and type_text.endswith(')')
            line += type_text if has_parens else f'({type_text})'
        for attr in attributes:
            line += f', {attr}'
        line += ' :: ' + ', '.join(variables) + '\n'
        out += line
    return out


def _source_digest_binding(name):
    """Return the C-binding declaration for the per-class MD5 symbol."""
    return (
        f'   character(C_Char), dimension(23), bind(C, name="{name}MD5") '
        f':: {name}5\n'
    )


def _ucfirst(s):
    return s[:1].upper() + s[1:] if s else s


def _lcfirst(s):
    return s[:1].lower() + s[1:] if s else s


def _short_name(full_name, directive_name):
    """Strip the directive's name prefix and apply the same `lcfirst-unless-
    all-caps` convention as loadAndSortClasses / buildObjectTypeMethod /
    generateDocumentation."""
    short = full_name
    if short.startswith(directive_name):
        short = short[len(directive_name):]
    if not re.match(r'^[A-Z]{2,}', short):
        short = _lcfirst(short)
    return short


def _xml_escape(s):
    """Escape `&`, `<`, `>`, `"` for XML attribute / text contexts."""
    return (
        str(s)
        .replace('&', '&amp;')
        .replace('<', '&lt;')
        .replace('>', '&gt;')
        .replace('"', '&quot;')
    )


_LATEX_ESCAPES = {
    '\\': r'\textbackslash{}',
    '_':  r'\_',
    '&':  r'\&',
    '%':  r'\%',
    '$':  r'\$',
    '#':  r'\#',
    '{':  r'\{',
    '}':  r'\}',
    '~':  r'\textasciitilde{}',
    '^':  r'\textasciicircum{}',
}


def _latex_encode(s):
    """Minimal LaTeX special-character escape covering the subset of
    characters this pipeline needs.
    """
    return ''.join(_LATEX_ESCAPES.get(c, c) for c in str(s))


# Module-level caches.  The helpers take the state explicitly as
# parameters, so these exist only to avoid repeat parses of the XML files
# within a single `process_function_class` run.
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
    return _shared_function_class_names(state_storables)


def _function_class_instances(state_storables):
    return _shared_function_class_instances(state_storables)


# ---------------------------------------------------------------------------
# Small helpers (fully implemented in D.7.3a)
# ---------------------------------------------------------------------------

def _is_function_class_pointer(declaration, type_name, state_storables):
    """Return True if `declaration` is a pointer to a registered
    functionClass (or functionClassInstance) whose type is `type_name`.
    """
    if declaration.get('intrinsic') not in ('class', 'type'):
        return False
    targets = _function_class_names(state_storables) | set(
        _function_class_instances(state_storables))
    if type_name not in targets:
        return False
    return any(a == 'pointer' for a in declaration.get('attributes') or [])


def _init_code_content():
    """Return `(code_content, pre_contains, post_contains)`."""
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
        # shortName is computed the same way loadAndSortClasses does —
        # lowercase the first letter unless the remainder is already an
        # all-caps acronym of 2+ letters.
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

    Returns `(classes_by_type, classes_ordered, non_abstract_classes)`.
    """
    class_locations = list(as_array(
        (directive_locations.get(directive['name']) or {}).get('file')))

    dependencies = {}
    classes      = {}
    class_counts = {}

    for class_location in class_locations:
        class_tree = parse_file(class_location)
        # Process the class tree so Class_Dependencies sees post-processing
        # structure (declarations, visibilities, etc.).
        process_tree(class_tree)

        class_record, class_deps = class_dependencies(
            class_tree, directive['name'])
        for dep in class_deps:
            type_name = class_record.get('type')
            if type_name and dep != type_name:
                # Sort.Topo's predecessor convention: `dependencies[X] = [Y, …]`
                # means X must come *after* Y — i.e. Y is X's predecessor.
                # `class_dependencies` returns the parents that this class
                # extends (and any cross-referenced sibling types), so each
                # `dep` must be emitted *before* `type_name`.  Earlier this
                # appended `type_name` under `dep`, which produced the
                # opposite order: a child (e.g. cosmologyFunctionsMatterDarkEnergy)
                # would be emitted before its parent
                # (cosmologyFunctionsMatterLambda), which gfortran rejects.
                dependencies.setdefault(type_name, []).append(dep)

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

    Files are updated only if their content actually changed: write to a
    temp file and rename only when the content differs (using a simple
    `os.path.exists` + `open/read` compare).
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
        # attaching them to the fresh submodule tree.
        block = code_content['submodule'][class_name]
        for kid in block.get('preContains', []) + block.get('postContains', []):
            kid['parent']  = None
            kid['sibling'] = None

        prepend_child_to_node(submodule_node, block.get('preContains', []))
        insert_post_contains(submodule_node, block.get('postContains', []))

        # Write the submodule with its `.lmap` line-number-mapping sidecar
        # and `.up` sentinel, exactly as preprocess.py does for directly
        # preprocessed files — postprocess.py needs the map to translate
        # compiler diagnostics back to original source lines, and the
        # Makefile's rebuild logic relies on the sentinel.
        submodule_content, submodule_mappings = serialize(
            file_node, annotate=True, strip_mappings=True)
        file_name = block['fileName']
        tmp_name  = file_name + '.tmp'
        with open(tmp_name, 'w') as fh:
            fh.write(submodule_content)
        file_changes_update(file_name, tmp_name, prove_update=True)
        with open(file_name + '.lmap', 'w') as fh:
            fh.write(submodule_mappings)


# ---------------------------------------------------------------------------
# Stubs for the heavy method-builders / body-emitters (D.7.3b, D.7.3c)
# ---------------------------------------------------------------------------

def _levenshtein(a, b):
    """Iterative Wagner-Fischer edit distance.  Used by
    buildDescriptorMethods when it suggests a "did you mean" alternative
    for mis-matched parameter names.
    """
    if a == b:
        return 0
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        curr = [i] + [0] * len(b)
        for j, cb in enumerate(b, 1):
            curr[j] = min(
                prev[j] + 1,
                curr[j - 1] + 1,
                prev[j - 1] + (0 if ca == cb else 1),
            )
        prev = curr
    return prev[-1]


def _descriptor_discover_class(non_abstract_class, directive, classes,
                               state_storables):
    """Walk one class's tree (and its parent classes, and the directive's
    base-class `<data>` declarations) to collect:
      - `potential_names`: classification of every declared member into
                           `parameters`, `enumerations`, `statefulTypes`,
                           `objects`, `linkedListObjects`, `linkedLists`
                           (delegated to `potential_descriptor_parameters`).
      - `descriptor_parameters`: the actual parameter/enumeration/
                                 statefulType/object/linkedList references
                                 found inside the matching parameters-
                                 constructor function via
                                 `<inputParameter>` / `<objectBuilder>`
                                 directives.
      - `sub_parameters`: dict of `name → {parent, source}` for every
                          `X = Y%subParameters(...)` line in constructor
                          code bodies.
      - `declaration_matches` (bool): the constructor has
                           `type(inputParameters) :: parameters`.
      - `parent_constructor_used` (bool): the constructor assigns
                           `result%extensionOf = …`.
      - `extension_of`: the parent-class type extracted from the type's
                           `extends(...)` opener.
      - `failure_message`: list of human-readable reasons the auto
                           descriptor can't be built.
      - `supported` (int): 1 OK, <1 an explicit failure code.
    """
    from Galacticus.Build.SourceTree.Process.FunctionClass.Descriptor import (
        potential_descriptor_parameters,
    )
    from Galacticus.Build.SourceTree.Process.FunctionClass.Utils import (
        trimlc, striplc, strip_variable_name,
    )
    from Galacticus.Build.SourceTree.Parse.Declarations import parse_declaration
    from Galacticus.Build.FortranUtils import get_fortran_line
    import io

    potential_names        = {}
    descriptor_parameters  = {}
    sub_parameters         = {}
    declaration_matches    = False
    parent_constructor_used = False
    extension_of           = None
    failure_message        = []
    supported              = 1

    non_abstract_class['hasCustomDescriptor'] = False

    # --- potentialNames pass: walk the class and parent-class type bodies.
    cls = non_abstract_class
    while cls is not None:
        node = (cls.get('tree') or {}).get('firstChild')
        while node is not None and (
                node.get('type') != 'type'
                or node.get('name') != cls.get('name')):
            node = node.get('sibling')
        if node is None:
            break
        if cls is non_abstract_class:
            m = re.search(
                r',\s*extends\s*\(\s*([a-zA-Z0-9_]+)\s*\)',
                node.get('opener') or '')
            if m:
                extension_of = m.group(1)

        walker = node.get('firstChild')
        while walker is not None:
            if walker.get('type') == 'declaration':
                potential_descriptor_parameters(
                    walker.get('declarations'),
                    non_abstract_class, cls, state_storables, potential_names)
            # See note at line 877: in our parse `contains` has no children;
            # post-contains members are siblings, so a sibling-walk reaches
            # them naturally.
            walker = walker.get('sibling')

        cls = (None if cls.get('extends') == directive['name']
               else classes.get(cls.get('extends')))

    # Directive-level `<data>` declarations.
    for data in as_array(directive.get('data')):
        declaration_source = None
        if isinstance(data, dict):
            if data.get('scope') == 'self':
                declaration_source = data.get('content')
        else:
            declaration_source = data
        if declaration_source is None:
            continue
        declaration = parse_declaration(declaration_source)
        if declaration is None:
            raise RuntimeError(
                "process_function_class: unable to parse variable declaration")
        potential_descriptor_parameters(
            declaration, non_abstract_class, None, state_storables,
            potential_names)

    # --- Locate the parameters constructor.
    node = (non_abstract_class.get('tree') or {}).get('firstChild')
    while node is not None and (
            node.get('type') != 'interface'
            or node.get('name') != non_abstract_class['name']):
        node = node.get('sibling')
    if node is None:
        return (potential_names, descriptor_parameters, sub_parameters,
                declaration_matches, parent_constructor_used,
                extension_of, failure_message, supported)

    constructors = []
    inner = node.get('firstChild')
    while inner is not None:
        if inner.get('type') == 'moduleProcedure':
            constructors.extend(inner.get('names') or [])
        inner = inner.get('sibling')

    # --- Walk the class tree looking for each constructor function.
    node = (non_abstract_class.get('tree') or {}).get('firstChild')
    while node is not None:
        if (node.get('type') == 'function'
                and node.get('name') in constructors):
            opener = node.get('opener') or ''
            name   = node['name']
            sig_re = (
                r'^\s*(recursive\s+)??function\s+' + re.escape(name)
                + r'\s*\(\s*parameters\s*\)'
            )
            if re.match(sig_re, opener):
                result_m = re.search(
                    r'result\s*\(\s*([a-zA-Z0-9_]+)\s*\)\s*$', opener)
                result = result_m.group(1) if result_m else name

                sub_iter = walk_tree(node)
                next(sub_iter)
                for cnode in sub_iter:
                    ctype = cnode.get('type')
                    if ctype == 'declaration':
                        for declaration in cnode.get('declarations') or []:
                            if (declaration.get('intrinsic') == 'type'
                                    and trimlc(declaration.get('type'))
                                        == 'inputparameters'
                                    and 'parameters'
                                        in (declaration.get('variables')
                                            or [])):
                                declaration_matches = True
                    elif ctype == 'code':
                        code_stream = io.StringIO(cnode.get('content') or '')
                        while True:
                            raw_line, processed_line, _ = get_fortran_line(
                                code_stream)
                            # get_fortran_line returns ('','','') at EOF
                            # rather than None, so detect that explicitly.
                            if not raw_line and not processed_line:
                                break
                            # Subparameter sources.
                            sm = re.match(
                                r'^\s*([a-zA-Z0-9_]+)\s*=\s*([a-zA-Z0-9_]+)'
                                r'\s*%\s*subParameters\s*\(',
                                processed_line)
                            if sm:
                                sub_parameters[sm.group(1)] = {
                                    'parent': sm.group(2),
                                    'source': processed_line,
                                }
                            # Parent constructor usage.
                            if extension_of and re.match(
                                    r'^\s*' + re.escape(result)
                                    + r'\s*%\s*' + re.escape(extension_of)
                                    + r'\s*=',
                                    processed_line):
                                parent_constructor_used = True
                    elif ctype == 'inputParameter':
                        d = cnode.setdefault('directive', {})
                        if 'source' not in d:
                            supported = -4
                            failure_message.append(
                                "unsourced parameters not supported")
                            continue
                        if 'name' not in d:
                            continue
                        # Resolve the internal name.
                        if 'variable' in d:
                            vm = re.match(r'(.*)%(.*)', d['variable'])
                            if vm:
                                obj  = vm.group(1)
                                elem = vm.group(2)
                                if obj.lower() == result.lower():
                                    raw_name = elem
                                else:
                                    raw_name = obj
                            else:
                                raw_name = re.sub(r'\(.+\)$', '',
                                                   d['variable'])
                        else:
                            raw_name = d['name']
                        name_lc    = raw_name.lower()
                        # Try parameters.
                        param_vars = [
                            v for decl in potential_names.get(
                                'parameters', [])
                            for v in decl.get('variableNames') or []]
                        if any(v.lower() == name_lc for v in param_vars):
                            descriptor_parameters.setdefault(
                                'parameters', []).append({
                                    'name':      raw_name,
                                    'inputName': d['name'],
                                    'source':    d['source'],
                                })
                        elif any(v.lower() == name_lc + '_'
                                 for v in param_vars):
                            descriptor_parameters.setdefault(
                                'parameters', []).append({
                                    'name':      raw_name + '_',
                                    'inputName': d['name'],
                                    'source':    d['source'],
                                })
                        elif any(
                                v.lower() == name_lc
                                for decl in potential_names.get(
                                    'enumerations', [])
                                for v in decl.get('variableNames') or []):
                            descriptor_parameters.setdefault(
                                'enumerations', []).append({
                                    'name':      raw_name,
                                    'inputName': d['name'],
                                    'source':    d['source'],
                                })
                        elif any(
                                v == name_lc
                                for decl in potential_names.get(
                                    'statefulTypes', [])
                                for v in decl.get('variables') or []):
                            descriptor_parameters.setdefault(
                                'statefulTypes', []).append({
                                    'name':      raw_name,
                                    'inputName': d['name'],
                                    'source':    d['source'],
                                })
                        else:
                            supported = -1
                            msg = (f"could not find a matching internal "
                                   f"variable for parameter [{raw_name}]")
                            if param_vars:
                                distances = [
                                    _levenshtein(name_lc, p.lower())
                                    for p in param_vars]
                                best = min(distances)
                                idx  = distances.index(best)
                                guess = re.sub(r'_', '', param_vars[idx])
                                msg += f" - did you mean [{guess}]"
                            failure_message.append(msg)
                    elif ctype == 'objectBuilder':
                        d = cnode.setdefault('directive', {})
                        if 'source' not in d:
                            supported = -6
                            failure_message.append(
                                "unsourced objects not supported")
                            continue
                        obj_name = re.sub(
                            r'([a-zA-Z0-9_]+\s*%\s*)?([a-zA-Z0-9_]+).*',
                            r'\2', d.get('name') or '')
                        obj_name = re.sub(r'\s', '', obj_name)
                        if any(o == obj_name.lower()
                               for o in potential_names.get('objects', [])):
                            descriptor_parameters.setdefault(
                                'objects', []).append({
                                    'name':   obj_name,
                                    'source': d['source'],
                                })
                        elif obj_name in potential_names.get(
                                'linkedListObjects', []):
                            descriptor_parameters.setdefault(
                                'linkedLists', []).append(
                                potential_names.get('linkedLists', {})
                                .get(obj_name))
                        else:
                            supported = -5
                            failure_message.append(
                                f"could not find a matching internal object "
                                f"for object [{obj_name}]")
        # In our parse `contains`
        # is a self-closing sibling marker — every post-contains procedure
        # is its sibling, so a plain sibling-walk visits them naturally.
        node = node.get('sibling')

    # --- Validate sub-parameter hierarchy.
    for sp_name in sorted(sub_parameters.keys()):
        parent = sub_parameters[sp_name]['parent']
        if parent not in sub_parameters and parent != 'parameters':
            supported = -7
            failure_message.append("subparameter hierarchy failure")

    return (potential_names, descriptor_parameters, sub_parameters,
            declaration_matches, parent_constructor_used,
            extension_of, failure_message, supported)


def _build_descriptor_methods(directive, non_abstract_classes, classes,
                              methods, tree, state_storables):
    """Populate `methods['descriptor']` (and in step 5c `hashedDescriptor`)
    with the auto-descriptor code that serialises the class's parameter
    constructor into an inputParameters tree.

    This
    sub-step (5a) implements the discovery pass and registers a skeleton
    descriptor method with the outer select-type frame; the per-class
    descriptor body (5b) and the hashedDescriptor companion (5c) follow.
    """
    descriptor_code             = ''
    descriptor_modules          = {'Input_Parameters': True}
    add_sub_parameters          = {}
    add_label                   = False
    rank_maximum                = 0
    descriptor_used             = False
    file_modification_added     = False
    descriptor_linked_list_vars = []

    descriptor_code += "logical :: includeFileModificationTimes_\n"
    descriptor_code += "if (present(includeFileModificationTimes)) then\n"
    descriptor_code += (
        " includeFileModificationTimes_=includeFileModificationTimes\n")
    descriptor_code += "else\n"
    descriptor_code += " includeFileModificationTimes_=.false.\n"
    descriptor_code += "end if\n"
    descriptor_code += "select type (self)\n"

    for non_abstract in non_abstract_classes:
        (potential_names, descriptor_parameters, sub_parameters,
         declaration_matches, parent_constructor_used, extension_of,
         failure_message, supported) = _descriptor_discover_class(
             non_abstract, directive, classes, state_storables)

        from Galacticus.Build.SourceTree.Process.FunctionClass.LinkedList \
            import auto_descriptor_linked_list
        from Galacticus.Build.SourceTree.Process.FunctionClass.Utils \
            import declaration_rank

        short = non_abstract['name']
        if short.startswith(directive['name']):
            short = short[len(directive['name']):]
        if not re.match(r'^[A-Z]{2,}', short):
            short = short[:1].lower() + short[1:]
        label = short

        descriptor_code += f"type is ({non_abstract['name']})\n"

        if non_abstract.get('hasCustomDescriptor'):
            loc_node = non_abstract.get('node') or {}
            loc_expr = location(loc_node, loc_node.get('line', 0))
            descriptor_code += (
                f" call Error_Report('custom descriptor exists - "
                f"this should not happen'//{loc_expr})\n"
            )
            descriptor_modules['Error'] = True
        else:
            if (declaration_matches
                    and (supported == 1
                         or 'descriptorSpecial' in non_abstract)):
                descriptor_used = True
                descriptor_code += " if (present(includeClass)) then\n"
                descriptor_code += "  includeClass_=includeClass\n"
                descriptor_code += " else\n"
                descriptor_code += "  includeClass_=.true.\n"
                descriptor_code += " end if\n"
                descriptor_code += (
                    f" if (includeClass_) call descriptor%addParameter"
                    f"('{directive['name']}','{label}')\n"
                )
                if descriptor_parameters:
                    add_sub_parameters['parameters'] = True
                    descriptor_code += (
                        f"parameters=descriptor%subparameters"
                        f"('{directive['name']}')\n"
                    )
                    for sp_name in sorted(sub_parameters.keys()):
                        add_sub_parameters[sp_name] = True
                        descriptor_code += sub_parameters[sp_name]['source']

                    for parameter in descriptor_parameters.get(
                            'parameters') or []:
                        for declaration in potential_names.get(
                                'parameters', []):
                            if parameter['name'].lower() not in (
                                    declaration.get('variables') or []):
                                continue
                            intrinsic = declaration.get('intrinsic') or ''
                            fmt       = None
                            function  = None
                            is_logical = False
                            if intrinsic == 'type':
                                function  = 'char'
                            elif intrinsic == 'logical':
                                add_label  = True
                                is_logical = True
                            elif intrinsic == 'double precision':
                                add_label  = True
                                fmt = 'e17.10'
                            elif intrinsic == 'integer':
                                add_label  = True
                                fmt = 'i17'
                            elif intrinsic == 'character':
                                function  = 'trim'
                            rank = declaration_rank(declaration)
                            if rank > 0:
                                if rank > rank_maximum:
                                    rank_maximum = rank
                                descriptor_code += "parameterValues=''\n"
                                for i in range(1, rank + 1):
                                    descriptor_code += (
                                        " parameterValues=parameterValues"
                                        "//'['\n"
                                    )
                                    descriptor_code += (
                                        f"do i{i}=lbound(self%"
                                        f"{parameter['name']},dim={i}),"
                                        f"ubound(self%{parameter['name']},"
                                        f"dim={i})\n"
                                    )
                                idx = ",".join(
                                    f"i{i}" for i in range(1, rank + 1))
                                if function:
                                    descriptor_code += (
                                        f" parameterValues="
                                        f"parameterValues//{function}"
                                        f"(self%{parameter['name']}"
                                        f"({idx}))\n"
                                    )
                                else:
                                    if is_logical:
                                        descriptor_code += (
                                            f"if (self%{parameter['name']}"
                                            f"({idx})) then\n"
                                        )
                                        descriptor_code += (
                                            "  parameterLabel='true'\n")
                                        descriptor_code += "else\n"
                                        descriptor_code += (
                                            "  parameterLabel='false'\n")
                                        descriptor_code += "end if\n"
                                    else:
                                        descriptor_code += (
                                            f"write (parameterLabel,"
                                            f"'({fmt})') self%"
                                            f"{parameter['name']}({idx})\n"
                                        )
                                    descriptor_code += (
                                        " parameterValues=parameterValues"
                                        "//trim(adjustl(parameterLabel))\n"
                                    )
                                descriptor_code += (
                                    f" if (i{rank} /= size(self%"
                                    f"{parameter['name']},dim={rank})) "
                                    f"parameterValues=parameterValues//','\n"
                                )
                                for i in range(1, rank + 1):
                                    descriptor_code += "end do\n"
                                    descriptor_code += (
                                        " parameterValues=parameterValues"
                                        "//']'\n"
                                    )
                                    if i != 1:
                                        descriptor_code += (
                                            f" if (i{i - 1} /= size(self%"
                                            f"{parameter['name']},"
                                            f"dim={i - 1})) "
                                            f"parameterValues="
                                            f"parameterValues//','\n"
                                        )
                                descriptor_code += (
                                    f"call {parameter['source']}"
                                    f"%addParameter('"
                                    f"{parameter['inputName']}',"
                                    f"char(parameterValues))\n"
                                )
                            else:
                                if function:
                                    descriptor_code += (
                                        f"call {parameter['source']}"
                                        f"%addParameter('"
                                        f"{parameter['inputName']}',"
                                        f"{function}(self%"
                                        f"{parameter['name']}))\n"
                                    )
                                else:
                                    if is_logical:
                                        descriptor_code += (
                                            f"if (self%{parameter['name']})"
                                            f" then\n"
                                        )
                                        descriptor_code += (
                                            "  parameterLabel='true'\n")
                                        descriptor_code += "else\n"
                                        descriptor_code += (
                                            "  parameterLabel='false'\n")
                                        descriptor_code += "end if\n"
                                    else:
                                        descriptor_code += (
                                            f"write (parameterLabel,"
                                            f"'({fmt})') self%"
                                            f"{parameter['name']}\n"
                                        )
                                    descriptor_code += (
                                        f"call {parameter['source']}"
                                        f"%addParameter('"
                                        f"{parameter['inputName']}',"
                                        f"trim(adjustl(parameterLabel)))\n"
                                    )

                    for parameter in descriptor_parameters.get(
                            'enumerations') or []:
                        for declaration in potential_names.get(
                                'enumerations', []):
                            if parameter['name'].lower() not in (
                                    declaration.get('variables') or []):
                                continue
                            fmt = 'i17'
                            add_label = True
                            rank = declaration_rank(declaration)
                            if rank > 0:
                                if rank > rank_maximum:
                                    rank_maximum = rank
                                descriptor_code += "parameterValues=''\n"
                                for i in range(1, rank + 1):
                                    descriptor_code += (
                                        " parameterValues=parameterValues"
                                        "//'['\n"
                                    )
                                    descriptor_code += (
                                        f"do i{i}=lbound(self%"
                                        f"{parameter['name']},dim={i}),"
                                        f"ubound(self%{parameter['name']},"
                                        f"dim={i})\n"
                                    )
                                idx = ",".join(
                                    f"i{i}" for i in range(1, rank + 1))
                                descriptor_code += (
                                    f"write (parameterLabel,'({fmt})') "
                                    f"self%{parameter['name']}({idx})%ID\n"
                                )
                                descriptor_code += (
                                    " parameterValues=parameterValues"
                                    "//trim(adjustl(parameterLabel))\n"
                                )
                                descriptor_code += (
                                    f" if (i{rank} /= size(self%"
                                    f"{parameter['name']},dim={rank})) "
                                    f"parameterValues=parameterValues//','\n"
                                )
                                for i in range(1, rank + 1):
                                    descriptor_code += "end do\n"
                                    descriptor_code += (
                                        " parameterValues=parameterValues"
                                        "//']'\n"
                                    )
                                    if i != 1:
                                        descriptor_code += (
                                            f" if (i{i - 1} /= size(self%"
                                            f"{parameter['name']},"
                                            f"dim={i - 1})) "
                                            f"parameterValues="
                                            f"parameterValues//','\n"
                                        )
                                descriptor_code += (
                                    f"call {parameter['source']}"
                                    f"%addParameter('"
                                    f"{parameter['inputName']}',"
                                    f"char(parameterValues))\n"
                                )
                            else:
                                descriptor_code += (
                                    f"write (parameterLabel,'({fmt})') "
                                    f"self%{parameter['name']}%ID\n"
                                )
                                descriptor_code += (
                                    f"call {parameter['source']}"
                                    f"%addParameter('"
                                    f"{parameter['inputName']}',"
                                    f"trim(adjustl(parameterLabel)))\n"
                                )

                    for parameter in descriptor_parameters.get(
                            'statefulTypes') or []:
                        for declaration in potential_names.get(
                                'statefulTypes', []):
                            if parameter['name'].lower() not in (
                                    declaration.get('variables') or []):
                                continue
                            t = declaration.get('type') or ''
                            fmt = None
                            is_logical = False
                            if t == 'statefulInteger':
                                fmt = 'i17'
                            elif t == 'statefulDouble':
                                fmt = 'e17.10'
                            elif t == 'statefulLogical':
                                is_logical = True
                            else:
                                raise RuntimeError("unknown stateful-type")
                            add_label = True
                            rank = declaration_rank(declaration)
                            if rank > 0:
                                if rank > rank_maximum:
                                    rank_maximum = rank
                                descriptor_code += "parameterValues=''\n"
                                for i in range(1, rank + 1):
                                    descriptor_code += (
                                        " parameterValues=parameterValues"
                                        "//'['\n"
                                    )
                                    descriptor_code += (
                                        f"do i{i}=lbound(self%"
                                        f"{parameter['name']},dim={i}),"
                                        f"ubound(self%{parameter['name']},"
                                        f"dim={i})\n"
                                    )
                                idx = ",".join(
                                    f"i{i}" for i in range(1, rank + 1))
                                descriptor_code += (
                                    f"if (self%{parameter['name']}"
                                    f"({idx})%isSet) then\n"
                                )
                                if is_logical:
                                    descriptor_code += (
                                        f"if (self%{parameter['name']}"
                                        f"({idx})%value) then\n"
                                    )
                                    descriptor_code += (
                                        "  parameterLabel='true'\n")
                                    descriptor_code += "else\n"
                                    descriptor_code += (
                                        "  parameterLabel='false'\n")
                                    descriptor_code += "end if\n"
                                else:
                                    descriptor_code += (
                                        f"write (parameterLabel,"
                                        f"'({fmt})') self%"
                                        f"{parameter['name']}({idx})"
                                        f"%value\n"
                                    )
                                descriptor_code += " else\n"
                                descriptor_code += (
                                    "  parameterLabel='?'\n")
                                descriptor_code += " end if\n"
                                descriptor_code += (
                                    " parameterValues=parameterValues//"
                                    "trim(adjustl(parameterLabel))\n"
                                )
                                descriptor_code += (
                                    f" if (i{rank} /= size(self%"
                                    f"{parameter['name']},dim={rank})) "
                                    f"parameterValues=parameterValues//','\n"
                                )
                                for i in range(1, rank + 1):
                                    descriptor_code += "end do\n"
                                    descriptor_code += (
                                        " parameterValues=parameterValues"
                                        "//']'\n"
                                    )
                                    if i != 1:
                                        descriptor_code += (
                                            f" if (i{i - 1} /= size(self%"
                                            f"{parameter['name']},"
                                            f"dim={i - 1})) "
                                            f"parameterValues="
                                            f"parameterValues//','\n"
                                        )
                                descriptor_code += (
                                    f"call {parameter['source']}"
                                    f"%addParameter('"
                                    f"{parameter['inputName']}',"
                                    f"char(parameterValues))\n"
                                )
                            else:
                                descriptor_code += (
                                    f"if (self%{parameter['name']}%isSet) "
                                    f"then\n"
                                )
                                if is_logical:
                                    descriptor_code += (
                                        f"if (self%{parameter['name']}"
                                        f"%value) then\n"
                                    )
                                    descriptor_code += (
                                        "  parameterLabel='true'\n")
                                    descriptor_code += "else\n"
                                    descriptor_code += (
                                        "  parameterLabel='false'\n")
                                    descriptor_code += "end if\n"
                                else:
                                    descriptor_code += (
                                        f"write (parameterLabel,"
                                        f"'({fmt})') self%"
                                        f"{parameter['name']}%value\n"
                                    )
                                descriptor_code += " else\n"
                                descriptor_code += (
                                    "  parameterLabel='?'\n")
                                descriptor_code += " end if\n"
                                descriptor_code += (
                                    f"call {parameter['source']}"
                                    f"%addParameter('"
                                    f"{parameter['inputName']}',"
                                    f"trim(adjustl(parameterLabel)))\n"
                                )

                    for obj in descriptor_parameters.get('objects') or []:
                        descriptor_code += (
                            f"if (associated(self%{obj['name']})) "
                            f"call self%{obj['name']}%descriptor"
                            f"(parameters,includeClass=.true.,"
                            f"includeFileModificationTimes="
                            f"includeFileModificationTimes)\n"
                        )
                    for linked in descriptor_parameters.get(
                            'linkedLists') or []:
                        linked_code, linked_module = (
                            auto_descriptor_linked_list(
                                linked, descriptor_linked_list_vars))
                        descriptor_code += linked_code
                        if linked_module:
                            descriptor_modules[linked_module] = True

                if parent_constructor_used:
                    descriptor_code += (
                        f"call self%{extension_of}%descriptor"
                        f"(descriptor,includeClass=.false.,"
                        f"includeFileModificationTimes="
                        f"includeFileModificationTimes)\n"
                    )
            elif (not declaration_matches
                    and 'descriptorSpecial' not in non_abstract):
                raise RuntimeError(
                    f"Automatic descriptor can not be built for class "
                    f"'{non_abstract['name']}': "
                    f"parameter-based constructor not found"
                )
            elif (supported != 1
                    and 'descriptorSpecial' not in non_abstract):
                raise RuntimeError(
                    f"Automatic descriptor can not be built for class "
                    f"'{non_abstract['name']}' because:\n   "
                    + "\n   ".join(failure_message)
                )

        # Run-time file dependencies.
        class_chain_rank_max = 0
        cls = non_abstract
        while cls is not None:
            if 'runTimeFileDependencies' in cls:
                if not file_modification_added:
                    descriptor_code = (
                        "integer :: status\n"
                        "character(len=30) :: timeModification\n"
                        "integer :: countRunTimeFileDependency\n"
                        "type(varying_string) :: "
                        "fileDependencyParameterName\n"
                        + descriptor_code
                    )
                    descriptor_modules['File_Utilities']  = True
                    descriptor_modules['String_Handling'] = True
                    descriptor_modules['Error']           = True
                    file_modification_added = True
                descriptor_code += (
                    "if (includeFileModificationTimes_) then\n"
                    "countRunTimeFileDependency=0\n"
                )
                paths = (
                    cls['runTimeFileDependencies']
                    .get('paths', '')).split()
                for path in paths:
                    rank = 0
                    for declaration in potential_names.get(
                            'parameters', []):
                        if path.lower() in (
                                declaration.get('variables') or []):
                            rank = declaration_rank(declaration)
                    if rank > class_chain_rank_max:
                        class_chain_rank_max = rank
                    loc_node = non_abstract.get('node') or {}
                    introspection = location(loc_node,
                                             loc_node.get('line', 0))
                    if rank > 0:
                        selector = ('('
                                    + ",".join(
                                        f"i{i}"
                                        for i in range(1, rank + 1))
                                    + ')')
                        for i in range(1, rank + 1):
                            descriptor_code += (
                                f"do i{i}=lbound(self%{path},dim={i}),"
                                f"ubound(self%{path},dim={i})\n"
                            )
                    else:
                        selector = ''
                    descriptor_code += (
                        f"timeModification="
                        f"File_Modification_Time(self%{path}{selector},"
                        f"status)\n"
                        f"if (status == errorStatusSuccess) then\n"
                        f" countRunTimeFileDependency="
                        f"countRunTimeFileDependency+1\n"
                        f" fileDependencyParameterName="
                        f"var_str(\"runTimeFileDependency\")//"
                        f"countRunTimeFileDependency\n"
                        f" call descriptor%addParameter("
                        f"char(fileDependencyParameterName),"
                        f"char(self%{path}{selector}//\": \"//"
                        f"trim(timeModification)))\n"
                        f"else if (status /= errorStatusNotExist) then\n"
                        f" call Error_Report('unable to get file "
                        f"modification time'//{introspection})\n"
                        f"end if\n"
                    )
                    if rank > 0:
                        for _ in range(rank):
                            descriptor_code += "end do\n"
                descriptor_code += "end if\n"
            cls = (None if cls.get('extends') == directive['name']
                   else classes.get(cls.get('extends')))
        # Remember whether the parameter-loop section needs `parameterValues`
        # before merging the file-dependency rank into `rank_maximum` —
        # parameterValues is only used for parameter loops.
        parameter_values_used = rank_maximum > 0
        # Loop indices for the runtime-file-dependency section and the
        # parameter-loop section share the names i1, i2, …, so we want a
        # single `integer :: i1, …, iN` declaration covering the larger of
        # the two ranks.  Emitting both separately produces duplicate
        # `integer :: i1` declarations and gfortran rejects with "Symbol
        # 'i1' at (1) already has basic type of INTEGER".  Fold the
        # file-dependency rank into the running max so the parameter-loop
        # block below emits the combined declaration.
        rank_maximum = max(rank_maximum, class_chain_rank_max)

        if 'descriptorSpecial' in non_abstract:
            descriptor_code += (
                f" call self%{non_abstract['descriptorSpecial']}"
                f"(parameters)\n"
            )

    descriptor_code += "end select\n"
    if descriptor_linked_list_vars:
        descriptor_code = (
            _format_variable_definitions(descriptor_linked_list_vars)
            + descriptor_code
        )

    if descriptor_used:
        descriptor_code = "logical :: includeClass_\n" + descriptor_code
    else:
        descriptor_code = (
            " !$GLC attributes unused :: descriptor, includeClass\n"
            + descriptor_code
        )
    if add_sub_parameters:
        descriptor_code = (
            "type(inputParameters) :: "
            + ",".join(sorted(add_sub_parameters.keys())) + "\n"
            + descriptor_code
        )
    if add_label:
        descriptor_code = (
            "character(len=18) :: parameterLabel\n" + descriptor_code)
    if rank_maximum > 0:
        prefix = (
            "integer :: "
            + ",".join(f"i{i}" for i in range(1, rank_maximum + 1))
            + "\n"
        )
        if parameter_values_used:
            prefix += "type(varying_string) :: parameterValues\n"
        descriptor_code = prefix + descriptor_code

    methods['descriptor'] = {
        'description': ('Return an input parameter list descriptor which '
                        'could be used to recreate this object.'),
        'type':        'void',
        'pass':        'yes',
        'modules':     ' '.join(sorted(descriptor_modules.keys())),
        'argument':    [
            'type(inputParameters), intent(inout) :: descriptor',
            'logical, intent(in   ), optional :: includeClass, '
            'includeFileModificationTimes',
        ],
        'code':        descriptor_code,
    }

    # --- hashedDescriptor ---
    directive_name = directive['name']
    hashed = (
        "logical                        :: includeSourceDigest_\n"
        "type   (inputParameters)       :: descriptor\n"
        "type   (varying_string )       :: descriptorString\n"
        "!   Workaround starts here.\n"
        "! type   (varying_string ), save :: descriptorStringPrevious, "
        "hashedDescriptorPrevious\n"
        "! !$omp threadprivate(descriptorStringPrevious,"
        "hashedDescriptorPrevious)\n"
        "! Workaround ends here.\n"
        "descriptor=inputParameters()\n"
        "! Disable live nodeLists in FoX as updating these nodeLists leads "
        "to memory leaks.\n"
        "call setLiveNodeLists(descriptor%document%document,.false.)\n"
        "call self%descriptor(descriptor,includeClass=.true.,"
        "includeFileModificationTimes=includeFileModificationTimes)\n"
        "descriptorString=descriptor%serializeToString()\n"
        "call descriptor%destroy()\n"
        "if (present(includeSourceDigest)) then\n"
        " includeSourceDigest_=includeSourceDigest\n"
        "else\n"
        " includeSourceDigest_=.false.\n"
        "end if\n"
        "if (includeSourceDigest_) then\n"
        "select type (self)\n"
    )
    for non_abstract in non_abstract_classes:
        hashed += (
            f"type is ({non_abstract['name']})\n"
            f"descriptorString=descriptorString//\":sourceDigest\\{{\"//"
            f"String_C_To_Fortran({non_abstract['name']}5)//\"\\}}\"\n"
        )
    hashed += (
        "end select\n"
        "end if\n"
        "!   Workaround starts here.\n"
        "!   if (descriptorString /= descriptorStringPrevious) then\n"
        "!      descriptorStringPrevious=         descriptorString\n"
        "!      hashedDescriptorPrevious=Hash_MD5(descriptorString)\n"
        "!   end if\n"
        f"!   {directive_name}HashedDescriptor=hashedDescriptorPrevious\n"
        f"   {directive_name}HashedDescriptor=Hash_MD5(descriptorString)\n"
        "! Workaround ends here.\n"
    )
    methods['hashedDescriptor'] = {
        'description': ('Return a hash of the descriptor for this object, '
                        'optionally include the source code digest in the '
                        'hash.'),
        'type':        'type(varying_string)',
        'pass':        'yes',
        'modules':     ('ISO_Varying_String String_Handling Input_Parameters '
                        'Hashes_Cryptographic FoX_DOM'),
        'argument':    [
            'logical, intent(in   ), optional :: includeSourceDigest, '
            'includeFileModificationTimes'
        ],
        'code':        hashed,
    }


def _build_allowed_parameters_method(directive, classes_ordered, methods):
    """Populate `methods['allowedParameters']` with the nested
    select-type ladder that announces every parameter name a class
    accepts via its inputParameter / objectBuilder directives.
    """
    from Galacticus.Build.SourceTree.Process.FunctionClass.LinkedList import (
        allowed_parameters_linked_list,
    )
    from Galacticus.Build.SourceTree.Process.FunctionClass.Utils import (
        trimlc, striplc, strip_variable_name,
    )
    from Galacticus.Build.FortranUtils import get_fortran_line
    import io

    directive_name = directive['name']
    allowed_parameters = {}
    parameters_present = False

    # --- Discovery pass: walk every class's constructors. ---
    for class_rec in classes_ordered:
        class_name = class_rec['name']
        allowed_parameters[class_name] = {
            'declarationMatches': False,
            'parameters':         {},
        }

        # Find the interface node named after the class.
        node = (class_rec.get('tree') or {}).get('firstChild')
        while node is not None and (
                node.get('type') != 'interface'
                or node.get('name') != class_name):
            node = node.get('sibling')
        if node is None:
            continue

        # Collect constructor names (from moduleProcedure children).
        constructors = []
        inner = node.get('firstChild')
        while inner is not None:
            if inner.get('type') == 'moduleProcedure':
                constructors.extend(inner.get('names') or [])
            inner = inner.get('sibling')

        # Walk the whole class tree looking for a function named
        # <constructor>(parameters).
        node = (class_rec.get('tree') or {}).get('firstChild')
        while node is not None:
            if (node.get('type') == 'function'
                    and node.get('name') in constructors):
                opener = node.get('opener') or ''
                name   = node['name']
                # Match `function NAME(parameters)` with an optional
                # `recursive` prefix. The mandatory whitespace lives INSIDE
                # the `(recursive\s+)?` group so non-recursive openers (which
                # have no leading word and no whitespace) match too — the
                # previous form `(recursive)??\s+function` required whitespace
                # before `function` unconditionally and therefore missed every
                # non-recursive constructor, leaving the `objects` accumulator
                # empty and the generated `allowedParameters` method without
                # its `if (associated(self%X_)) call self%X_%allowedParameters
                # (allowedParameters,'parameters',.true.)` lines for nested
                # object pointers.
                sig_re = (
                    r'^\s*(recursive\s+)?function\s+' + re.escape(name)
                    + r'\s*\(\s*parameters\s*\)'
                )
                if re.match(sig_re, opener):
                    result_m = re.search(
                        r'result\s*\(\s*([a-zA-Z0-9_]+)\s*\)\s*$', opener)
                    result = result_m.group(1) if result_m else name

                    # Walk the function subtree.
                    sub_iter = walk_tree(node)
                    next(sub_iter)  # skip the function node itself
                    for cnode in sub_iter:
                        ctype = cnode.get('type')
                        if ctype == 'declaration':
                            for declaration in (
                                    cnode.get('declarations') or []):
                                if (declaration.get('intrinsic') == 'type'
                                        and trimlc(declaration.get('type'))
                                            == 'inputparameters'
                                        and 'parameters'
                                            in (declaration.get('variables')
                                                or [])):
                                    allowed_parameters[class_name][
                                        'declarationMatches'] = True
                        elif ctype == 'code':
                            # Rewrite `result%X=Y(parameters)` lines to
                            # bracket them with DsblVldtn increments.
                            new_content = ''
                            modified = False
                            code = io.StringIO(cnode.get('content') or '')
                            while True:
                                raw_line, processed_line, _ = get_fortran_line(code)
                                # get_fortran_line returns ('','','') at EOF,
                                # not None.
                                if not raw_line and not processed_line:
                                    break
                                pm = re.match(
                                    r'^\s*' + re.escape(result)
                                    + r'%([a-zA-Z0-9_]+)\s*=([a-zA-Z0-9_]+)'
                                    + r'\(\s*parameters\s*\)',
                                    processed_line)
                                if pm:
                                    allowed_parameters[class_name][
                                        'classParent'] = pm.group(1)
                                    new_content += (
                                        f"{directive_name}DsblVldtn="
                                        f"{directive_name}DsblVldtn+1\n"
                                    )
                                    new_content += raw_line
                                    new_content += (
                                        f"{directive_name}DsblVldtn="
                                        f"{directive_name}DsblVldtn-1\n"
                                    )
                                    modified = True
                                else:
                                    new_content += raw_line
                            if modified:
                                cnode['content'] = new_content
                        elif ctype == 'inputParameter':
                            src = (cnode.setdefault('directive', {})).get('source')
                            iname = (cnode.setdefault('directive', {})).get('name')
                            if src is not None and iname is not None:
                                slot = (allowed_parameters[class_name]
                                        ['parameters']
                                        .setdefault(src, {})
                                        .setdefault('all', []))
                                slot.append(iname)
                        elif ctype == 'objectBuilder':
                            src = (cnode.setdefault('directive', {})).get('source')
                            if src is None:
                                continue
                            d = cnode.setdefault('directive', {})
                            param_name = d.get('parameterName') or d.get('class')
                            slot_src = (allowed_parameters[class_name]
                                        ['parameters']
                                        .setdefault(src, {}))
                            slot_src.setdefault('all',     []).append(param_name)
                            slot_src.setdefault('classes', []).append(param_name)
                            # Scan the class's type block for a matching
                            # pointer member and collect its variable name.
                            type_node = (class_rec.get('tree') or {}).get(
                                'firstChild')
                            while type_node is not None:
                                if (type_node.get('type') == 'type'
                                        and (type_node.get('name') or '').lower()
                                            == class_name.lower()):
                                    child = type_node.get('firstChild')
                                    while child is not None:
                                        if child.get('type') == 'declaration':
                                            for dec in (
                                                    child.get('declarations')
                                                    or []):
                                                if (dec.get('intrinsic') == 'class'
                                                        and trimlc(dec.get('type'))
                                                            == trimlc(
                                                                d.get('class'))
                                                            + 'class'):
                                                    slot_objects = (
                                                        slot_src.setdefault(
                                                            'objects', []))
                                                    for v in (
                                                            dec.get('variables')
                                                            or []):
                                                        # `v` arrives from
                                                        # `parse_declaration`'s
                                                        # `keep_qualifiers=True`
                                                        # mode and so still
                                                        # carries any
                                                        # `=>null()` or
                                                        # `=initial` tail
                                                        # (e.g.
                                                        # `galacticfilter_=>null()`).
                                                        # `strip_variable_name`
                                                        # removes that tail
                                                        # so the comparison
                                                        # against the
                                                        # directive's `name`
                                                        # attribute (which is
                                                        # the bare member
                                                        # name) actually
                                                        # matches — without
                                                        # this strip, every
                                                        # `<objectBuilder>`
                                                        # falls through and
                                                        # `objects` stays
                                                        # empty, dropping the
                                                        # nested
                                                        # `if (associated(self%X_))
                                                        # call self%X_%allowedParameters(…)`
                                                        # forwarding lines.
                                                        bare = strip_variable_name(v)
                                                        cmp_direct = (
                                                            bare.lower()
                                                            == striplc(
                                                                d.get('name')))
                                                        cmp_qual = (
                                                            (result
                                                             + '%' + bare).lower()
                                                            == striplc(
                                                                d.get('name')))
                                                        if cmp_direct or cmp_qual:
                                                            slot_objects.append(bare)
                                        child = child.get('sibling')
                                    break
                                type_node = type_node.get('sibling')

            # Same `contains` walk as in `_descriptor_discover_class`: in
            # our parse `contains` has no children, so just walk siblings —
            # the post-contains procedures are siblings of the marker.
            node = node.get('sibling')

    # --- Emission pass. ---
    allowed_parameters_linked_list_variables = []
    code = "select type (self)\n"
    modules = {'ISO_Varying_String': True}

    for class_rec in classes_ordered:
        class_name = class_rec['name']
        if not allowed_parameters[class_name]['declarationMatches']:
            continue
        code += f"type is ({class_name})\n"
        climb_name = class_name
        while climb_name is not None:
            entry = allowed_parameters.get(climb_name) or {}
            for src in sorted((entry.get('parameters') or {}).keys()):
                slot = entry['parameters'][src]
                code += "  if (objectsOnly) then\n"
                classes_list = slot.get('classes') or []
                if classes_list:
                    parameters_present = True
                    code += f"   if (sourceName == '{src}') then\n"
                    code += "     countNew=0\n"
                    code += "     if (allocated(allowedParameters)) then\n"
                    for entry_name in classes_list:
                        code += "       isNew=.true.\n"
                        code += "       do j=1,size(allowedParameters)\n"
                        code += (
                            f"          if (allowedParameters(j) "
                            f"== '{entry_name}') then\n")
                        code += "             isNew=.false.\n"
                        code += "             exit\n"
                        code += "          end if\n"
                        code += "       end do\n"
                        code += "       if (isNew) countNew=countNew+1\n"
                    code += "       if (countNew > 0) then\n"
                    code += (
                        "         call move_alloc("
                        "allowedParameters,allowedParametersTmp)\n")
                    code += (
                        "         allocate(allowedParameters("
                        "size(allowedParametersTmp)+countNew))\n")
                    code += (
                        "         allowedParameters("
                        "1:size(allowedParametersTmp))=allowedParametersTmp\n")
                    code += "         deallocate(allowedParametersTmp)\n"
                    for entry_name in classes_list:
                        code += "       isNew=.true.\n"
                        code += (
                            "       do j=1,"
                            "size(allowedParameters)-countNew\n")
                        code += (
                            f"          if (allowedParameters(j) "
                            f"== '{entry_name}') then\n")
                        code += "             isNew=.false.\n"
                        code += "             exit\n"
                        code += "          end if\n"
                        code += "       end do\n"
                        code += "       if (isNew) then\n"
                        code += "           countNew=countNew-1\n"
                        code += (
                            f"           allowedParameters("
                            f"size(allowedParameters)-countNew)"
                            f"='{entry_name}'\n")
                        code += "         end if\n"
                    code += "       end if\n"
                    code += "     else\n"
                    code += (
                        f"       allocate(allowedParameters("
                        f"{len(classes_list)}))\n")
                    for i, entry_name in enumerate(classes_list):
                        code += (
                            f"       allowedParameters({i + 1})"
                            f"='{entry_name}'\n")
                    code += "     end if\n"
                    code += "   end if\n"
                code += "  else\n"
                all_list = slot.get('all') or []
                if all_list:
                    parameters_present = True
                    code += f"   if (sourceName == '{src}') then\n"
                    code += "     countNew=0\n"
                    code += "     if (allocated(allowedParameters)) then\n"
                    for entry_name in all_list:
                        code += "       isNew=.true.\n"
                        code += "       do j=1,size(allowedParameters)\n"
                        code += (
                            f"          if (allowedParameters(j) "
                            f"== '{entry_name}') then\n")
                        code += "             isNew=.false.\n"
                        code += "             exit\n"
                        code += "          end if\n"
                        code += "       end do\n"
                        code += "       if (isNew) countNew=countNew+1\n"
                    code += "       if (countNew > 0) then\n"
                    code += (
                        "         call move_alloc("
                        "allowedParameters,allowedParametersTmp)\n")
                    code += (
                        "         allocate(allowedParameters("
                        "size(allowedParametersTmp)+countNew))\n")
                    code += (
                        "         allowedParameters("
                        "1:size(allowedParametersTmp))=allowedParametersTmp\n")
                    code += "         deallocate(allowedParametersTmp)\n"
                    for entry_name in all_list:
                        code += "       isNew=.true.\n"
                        code += (
                            "       do j=1,"
                            "size(allowedParameters)-countNew\n")
                        code += (
                            f"          if (allowedParameters(j) "
                            f"== '{entry_name}') then\n")
                        code += "             isNew=.false.\n"
                        code += "             exit\n"
                        code += "          end if\n"
                        code += "       end do\n"
                        code += "       if (isNew) then\n"
                        code += "           countNew=countNew-1\n"
                        code += (
                            f"           allowedParameters("
                            f"size(allowedParameters)-countNew)"
                            f"='{entry_name}'\n")
                        code += "         end if\n"
                    code += "       end if\n"
                    code += "     else\n"
                    code += (
                        f"       allocate(allowedParameters("
                        f"{len(all_list)}))\n")
                    for i, entry_name in enumerate(all_list):
                        code += (
                            f"       allowedParameters({i + 1})"
                            f"='{entry_name}'\n")
                    code += "     end if\n"
                    code += "   end if\n"
                code += "  end if\n"

                # Only emit objects/linked-list for the starting class.
                if climb_name == class_name:
                    if slot.get('objects'):
                        parameters_present = True
                    for obj_name in slot.get('objects') or []:
                        code += (
                            f"  if (associated(self%{obj_name})) "
                            f"call self%{obj_name}%allowedParameters"
                            f"(allowedParameters,'{src}',.true.)\n"
                        )
                    linked_code, linked_module = (
                        allowed_parameters_linked_list(
                            class_rec,
                            allowed_parameters_linked_list_variables, src))
                    code += linked_code
                    if linked_module:
                        modules[linked_module] = True
            climb_name = entry.get('classParent')
    code += "end select\n"

    if parameters_present:
        code = (
            "type   (varying_string), allocatable, dimension(:) :: "
            "allowedParametersTmp\n"
            + code)
        code = (
            "integer                                            :: "
            "countNew, j\n"
            + code)
        code = (
            "logical                                            :: isNew\n"
            + code)
    else:
        code = (
            "!$GLC attributes unused :: self, allowedParameters, "
            "sourceName\n"
        )
    code = (_format_variable_definitions(
                allowed_parameters_linked_list_variables)
            + code)

    methods['allowedParameters'] = {
        'description': (
            'Return a list of parameter names allowed for this object.'),
        'type':        'void',
        'recursive':   'yes',
        'pass':        'yes',
        'modules':     ' '.join(modules.keys()),
        'argument':    [
            'type     (varying_string), dimension(:), allocatable, '
            'intent(inout) :: allowedParameters',
            'character(len=*         )                           , '
            'intent(in   ) :: sourceName',
            'logical                                             , '
            'intent(in   ) :: objectsOnly',
        ],
        'code':        code,
    }


def _build_assignment_method(directive, non_abstract_classes, classes,
                             methods, state_storables):
    """Populate `methods['assignment(=)']` with the defined-assignment
    operator overload.
    """
    from Galacticus.Build.SourceTree.Process.FunctionClass.DeepCopy import (
        generate_assignment_allocatable_code,
    )
    from Galacticus.Build.SourceTree.Process.FunctionClass.LinkedList import (
        assigner_linked_list,
    )
    from Galacticus.Build.SourceTree.Parse.Declarations import parse_declaration

    assignment = {'code': ''}
    rank_maximum_ref = [0]
    assigner_modules = {'Error': True}
    assigner_linked_list_variables = []

    assignment['code'] += "select type (self)\n"
    for non_abstract in non_abstract_classes:
        assignment['code'] += f"type is ({non_abstract['name']})\n"
        assignment['code'] += "  select type (from)\n"
        assignment['code'] += f"  type is ({non_abstract['name']})\n"

        cls = non_abstract
        while cls is not None:
            node = (cls.get('tree') or {}).get('firstChild')
            while node is not None and (
                    node.get('type') != 'type'
                    or node.get('name') != cls.get('name')):
                node = node.get('sibling')
            if node is None:
                break

            child = node.get('firstChild')
            while child is not None:
                if child.get('type') == 'declaration':
                    for declaration in as_array(child.get('declarations')):
                        if not isinstance(declaration, dict):
                            continue
                        # Skip type-bound markers — `generic`/`final` and
                        # bare `procedure ::` bindings — none of which are
                        # data members.  A `procedure(intf), pointer ::
                        # foo` is *also* an `intrinsic == 'procedure'`
                        # declaration but IS a data member that must be
                        # copied (with `=>`); that case falls through to
                        # the assignment loop below.
                        intrinsic_decl = declaration.get('intrinsic')
                        decl_attrs     = declaration.get('attributes') or []
                        if intrinsic_decl in ('generic', 'final'):
                            continue
                        if (intrinsic_decl == 'procedure'
                                and not any(a == 'pointer' for a in decl_attrs)):
                            continue
                        attributes = decl_attrs
                        is_pointer     = any(a == 'pointer'     for a in attributes)
                        is_allocatable = any(a == 'allocatable' for a in attributes)
                        assigner = '=>' if is_pointer else '='
                        type_text = declaration.get('type') or ''
                        intrinsic = declaration.get('intrinsic') or ''
                        if intrinsic in ('class', 'type'):
                            type_text = type_text.strip()
                        reference_count = _is_function_class_pointer(
                            declaration, type_text, state_storables)
                        variables = declaration.get('variables') or []
                        for obj in variables:
                            from Galacticus.Build.SourceTree.Process.\
                                FunctionClass.Utils import strip_variable_name
                            name = strip_variable_name(obj)
                            allocatable_this = is_allocatable
                            allocated = 'allocated'
                            if ('assignment' in cls
                                    and isinstance(cls.get('assignment'), dict)
                                    and 'forceArrayAssign'
                                        in cls['assignment']
                                    and name in cls['assignment']
                                        ['forceArrayAssign'].split()):
                                allocatable_this = True
                                allocated = 'associated'
                            if allocatable_this:
                                generate_assignment_allocatable_code(
                                    assignment, declaration, name, allocated,
                                    rank_maximum_ref)
                            elif ('linkedList' in cls
                                    and name in (
                                        cls['linkedList'].get('variable')
                                        or '').split()):
                                pass  # handled later via assigner_linked_list
                            else:
                                assignment['code'] += (
                                    f"    self%{name}{assigner}from%{name}\n"
                                )
                                force_reference_count = (
                                    isinstance(cls.get('assignment'), dict)
                                    and isinstance(
                                        cls['assignment'].get('functionClass'),
                                        dict)
                                    and any(
                                        v.lower() == name.lower()
                                        for v in (cls['assignment']
                                                  ['functionClass']
                                                  .get('variables') or '')
                                        .split()
                                    )
                                )
                                if force_reference_count:
                                    assignment['code'] += (
                                        f"    select type (object_ => "
                                        f"self%{name})\n"
                                        f"    class is (functionClass)\n"
                                        f"    "
                                        f"{'if (associated(object_)) ' if is_pointer else ''}"
                                        f"call object_%referenceCountIncrement()\n"
                                        f"    end select\n"
                                    )
                                elif reference_count:
                                    assignment['code'] += (
                                        f"    "
                                        f"{'if (associated(self%'+name+')) ' if is_pointer else ''}"
                                        f"call self%{name}%referenceCountIncrement()\n"
                                    )
                child = child.get('sibling')

            if 'linkedList' in cls:
                linked_code, linked_module = assigner_linked_list(
                    cls['linkedList'], assigner_linked_list_variables)
                assignment['code'] += linked_code
                if linked_module:
                    assigner_modules[linked_module] = True

            cls = (None if cls.get('extends') == directive['name']
                   else classes.get(cls.get('extends')))

        assignment['code'] += "  class default\n"
        loc_node = non_abstract.get('node') or {}
        loc_expr = location(loc_node, loc_node.get('line', 0))
        assignment['code'] += (
            f"    call Error_Report('self and from types do not match'//"
            f"{loc_expr})\n"
        )
        assignment['code'] += "  end select\n"
    assignment['code'] += "end select\n"

    # Base-class `<data>` declarations.
    for data in as_array(directive.get('data')):
        declaration_source = None
        if isinstance(data, dict):
            if data.get('scope') == 'self':
                declaration_source = data.get('content')
        else:
            declaration_source = data
        if declaration_source is None:
            continue
        declaration = parse_declaration(declaration_source)
        if declaration is None:
            continue
        attributes = declaration.get('attributes') or []
        is_pointer     = any(a == 'pointer'     for a in attributes)
        is_allocatable = any(a == 'allocatable' for a in attributes)
        assigner = '=>' if is_pointer else '='
        type_text = declaration.get('type') or ''
        intrinsic = declaration.get('intrinsic') or ''
        if intrinsic in ('class', 'type'):
            type_text = type_text.strip()
        reference_count = _is_function_class_pointer(
            declaration, type_text, state_storables)
        for obj in declaration.get('variables') or []:
            from Galacticus.Build.SourceTree.Process.FunctionClass.Utils \
                import strip_variable_name
            name = strip_variable_name(obj)
            if is_allocatable:
                generate_assignment_allocatable_code(
                    assignment, declaration, name, 'allocated',
                    rank_maximum_ref)
            else:
                assignment['code'] += (
                    f"    self%{name}{assigner}from%{name}\n"
                )
            if reference_count:
                assignment['code'] += (
                    f"    "
                    f"{'if (associated(self%'+name+')) ' if is_pointer else ''}"
                    f"call self%{name}%referenceCountIncrement()\n"
                )

    # functionClass base-class fields.
    assignment['code'] += "self%isDefaultOfClass=from%isDefaultOfClass\n"
    assignment['code'] += "self%referenceCount=from%referenceCount\n"
    assignment['code'] += "return\n"

    # Prepend required variable declarations.
    if rank_maximum_ref[0] > 0:
        assignment['code'] = (
            'integer :: '
            + ','.join(f"i{i}__" for i in range(1, rank_maximum_ref[0] + 1))
            + '\n' + assignment['code']
        )
    if assigner_linked_list_variables:
        assignment['code'] = (
            _format_variable_definitions(assigner_linked_list_variables)
            + assignment['code']
        )

    methods['assignment(=)'] = {
        'description': 'Assign the object.',
        'type':        'void',
        'recursive':   'yes',
        'pass':        'yes',
        'selfIntent':  'out',
        'modules':     ' '.join(sorted(assigner_modules.keys())),
        'argument':    [
            f"class({directive['name']}Class), intent(in   ) :: from"
        ],
        'code':        assignment['code'],
    }


def _build_deep_copy_methods(directive, non_abstract_classes, classes,
                             line_number, methods,
                             state_storables, deep_copy_actions):
    """Populate `methods` with `deepCopy`, `deepCopy_`, `deepCopyReset`,
    `deepCopyFinalize`.
    """
    from Galacticus.Build.SourceTree.Process.FunctionClass.DeepCopy   import (
        deep_copy_declarations,
    )
    from Galacticus.Build.SourceTree.Process.FunctionClass.LinkedList import (
        deep_copy_linked_list,
    )
    from Galacticus.Build.SourceTree.Parse.Declarations import parse_declaration

    deep_copy = {
        'rankMaximum':        0,
        'needReferenceCount': 0,
        'code':               '',
        'resetCode':          '',
        'finalizeCode':       '',
        'assignments':        '',
        'modules':            {},
        'resetModules':       {},
        'finalizeModules':    {},
    }
    linked_list_variables          = []
    linked_list_reset_variables    = []
    linked_list_finalize_variables = []

    deep_copy['resetCode']    += "self%copiedSelf => null()\n"
    deep_copy['resetCode']    += "select type (self)\n"
    # NOTE (issue #695, D4): deepCopyFinalize deliberately does NOT null copiedSelf. A generated recursive shim
    # (<name>Recursive) resolves its weak recursiveSelf pointer to the copy of the real object via
    # recursiveSelf%copiedSelf, and that resolution runs during deepCopyFinalize. Because finalize is top-down and
    # the real object is an ancestor of the shim, nulling copiedSelf here would clear it before the shim could read
    # it. copiedSelf is only ever read inside a deepCopy operation, which always begins with deepCopyReset nulling
    # it, so leaving it stale between operations is harmless.
    deep_copy['finalizeCode'] += "select type (self)\n"
    deep_copy['code']         += "select type (self)\n"

    for non_abstract in non_abstract_classes:
        cls = non_abstract
        deep_copy['assignments'] = ''
        deep_copy['resetCode']    += f"type is ({non_abstract['name']})\n"
        deep_copy['finalizeCode'] += f"type is ({non_abstract['name']})\n"
        found_deep_copy_names = []

        node = None
        while cls is not None:
            node = (cls.get('tree') or {}).get('firstChild')
            while node is not None and (
                    node.get('type') != 'type'
                    or node.get('name') != cls.get('name')):
                node = node.get('sibling')
            if node is None:
                break

            linked_code, linked_reset, linked_finalize, linked_module = (
                deep_copy_linked_list(
                    cls, non_abstract, linked_list_variables,
                    linked_list_reset_variables,
                    linked_list_finalize_variables))
            deep_copy['assignments']  += linked_code
            deep_copy['resetCode']    += linked_reset
            deep_copy['finalizeCode'] += linked_finalize
            if linked_module:
                deep_copy['modules'][linked_module]         = True
                deep_copy['resetModules'][linked_module]    = True
                deep_copy['finalizeModules'][linked_module] = True

            dc = cls.get('deepCopy') or {}
            ignore_block = dc.get('ignore') if isinstance(dc, dict) else None
            ignore = []
            if isinstance(ignore_block, dict) and 'variables' in ignore_block:
                ignore = [
                    v.strip() for v in re.split(
                        r'\s*,\s*', ignore_block['variables'])
                    if v.strip()
                ]

            child = node.get('firstChild')
            while child is not None:
                if child.get('type') == 'declaration':
                    deep_copy_declarations(
                        cls, non_abstract, child,
                        child.get('declarations'),
                        ignore, line_number, deep_copy,
                        found_deep_copy_names,
                        state_storables, deep_copy_actions)
                child = child.get('sibling')

            cls = (None if cls.get('extends') == directive['name']
                   else classes.get(cls.get('extends')))

        # Base-class <data> declarations.
        for data in as_array(directive.get('data')):
            declaration_source = None
            if isinstance(data, dict):
                if data.get('scope') == 'self':
                    declaration_source = data.get('content')
            else:
                declaration_source = data
            if declaration_source is None:
                continue
            declaration = parse_declaration(declaration_source)
            if declaration is None:
                continue
            deep_copy_declarations(
                cls, non_abstract, node, [declaration], [],
                line_number, deep_copy, found_deep_copy_names,
                state_storables, deep_copy_actions)

        deep_copy['code'] += f"type is ({non_abstract['name']})\n"
        deep_copy['code'] += "select type (destination)\n"
        deep_copy['code'] += f"type is ({non_abstract['name']})\n"
        deep_copy['code'] += "destination=self\n"
        deep_copy['code'] += deep_copy['assignments']
        deep_copy['code'] += "class default\n"
        loc_node = non_abstract.get('node') or {}
        loc_expr = location(loc_node, loc_node.get('line', 0))
        deep_copy['code'] += (
            f"call Error_Report('destination and source types do not match'//"
            f"{loc_expr})\n"
        )
        deep_copy['code'] += "end select\n"

        deep_copy['modules']['Error'] = True

        # Validate explicit deepCopy functionClass variables.
        cls = non_abstract
        while cls is not None:
            dc = cls.get('deepCopy') or {}
            fc_block = dc.get('functionClass') if isinstance(dc, dict) else None
            if isinstance(fc_block, dict) and 'variables' in fc_block:
                for variable in re.split(r'\s*,\s*',
                                          fc_block['variables']):
                    variable = variable.strip()
                    if not variable:
                        continue
                    if variable.lower() not in {
                            v.lower() for v in found_deep_copy_names}:
                        raise RuntimeError(
                            f"Error: unable to find variable '{variable}' "
                            f"marked for deep copy in class '{cls['name']}'"
                        )
            cls = (None if cls.get('extends') == directive['name']
                   else classes.get(cls.get('extends')))

    deep_copy['code']         += "end select\n"
    deep_copy['resetCode']    += "end select\n"
    deep_copy['finalizeCode'] += "end select\n"
    deep_copy['code']         += "call destination%referenceCountReset()\n"
    deep_copy['code']         += (
        "destination%stateOperationID=0_c_size_t\n")

    # Prepend variable declarations.
    deep_copy['code']         = (
        _format_variable_definitions(linked_list_variables)
        + deep_copy['code'])
    deep_copy['resetCode']    = (
        _format_variable_definitions(linked_list_reset_variables)
        + deep_copy['resetCode'])
    deep_copy['finalizeCode'] = (
        _format_variable_definitions(linked_list_finalize_variables)
        + deep_copy['finalizeCode'])
    if deep_copy['rankMaximum'] > 0:
        deep_copy['code'] = (
            'integer :: '
            + ','.join(f"i{i}" for i in range(1, deep_copy['rankMaximum'] + 1))
            + '\n' + deep_copy['code']
        )
    if deep_copy['needReferenceCount']:
        deep_copy['code'] = (
            'integer :: referenceCount__\n' + deep_copy['code']
        )

    directive_name = directive['name']
    methods['deepCopy'] = {
        'description': ('Perform a deep copy of the object. This is a wrapper '
                        'around the actual deep-copy code.'),
        'type':        'void',
        'recursive':   'yes',
        'pass':        'yes',
        'selfTarget':  'yes',
        'argument':    [
            f"class({directive_name}Class), intent(inout) :: destination"
        ],
        'code':        "call self%deepCopy_(destination)",
    }
    methods['deepCopy_'] = {
        'description': 'Perform a deep copy of the object.',
        'type':        'void',
        'recursive':   'yes',
        'pass':        'yes',
        'modules':     ' '.join(sorted(deep_copy['modules'].keys())),
        'argument':    [
            f"class({directive_name}Class), intent(inout) :: destination"
        ],
        'code':        deep_copy['code'],
    }
    methods['deepCopyReset'] = {
        'description': ('Reset deep copy pointers in this object and any '
                        'objects that it uses.'),
        'type':        'void',
        'recursive':   'yes',
        'pass':        'yes',
        'code':        deep_copy['resetCode'],
    }
    methods['deepCopyFinalize'] = {
        'description': ('Finalize a deep copy in this object and any '
                        'objects that it uses.'),
        'type':        'void',
        'recursive':   'yes',
        'pass':        'yes',
        'code':        deep_copy['finalizeCode'],
    }
    if deep_copy['resetModules']:
        methods['deepCopyReset']['modules'] = ' '.join(
            sorted(deep_copy['resetModules'].keys()))
    if deep_copy['finalizeModules']:
        methods['deepCopyFinalize']['modules'] = ' '.join(
            sorted(deep_copy['finalizeModules'].keys()))


def _build_state_store_methods(directive, non_abstract_classes, classes,
                               methods, state_storables):
    """Populate `methods` with `stateStore`, `stateStore_`, `stateRestore`,
    `stateRestore_`.
    """
    from Galacticus.Build.SourceTree.Process.FunctionClass.StateStore  import (
        state_store_variables, state_store_explicit_function,
    )
    from Galacticus.Build.SourceTree.Process.FunctionClass.LinkedList  import (
        state_store_linked_list,
    )
    from Galacticus.Build.SourceTree.Parse.Declarations import parse_declaration

    state_stores = {
        'stateFileUsed':          False,
        'gslStateFileUsed':       False,
        'rankMaximum':            0,
        'allocatablesFound':      False,
        'explicitFunctionsFound': False,
        'dimensionalsFound':      False,
        'labelUsed':              False,
        'stateStoreModules':      {
            'Display':            True,
            'ISO_Varying_String': True,
            'String_Handling':    True,
            'ISO_C_Binding':      True,
        },
        'stateRestoreModules':    {
            'Display':            True,
            'ISO_Varying_String': True,
            'String_Handling':    True,
            'ISO_C_Binding':      True,
        },
    }
    state_store_code   = ''
    state_restore_code = ''
    state_linked_list_variables = []
    output_unused_variables = []
    input_unused_variables  = []

    state_store_code   += "position=FTell(stateFile)\n"
    state_restore_code += "position=FTell(stateFile)\n"
    directive_name = directive['name']
    state_store_code += (
        f"call displayIndent(var_str('storing state for \"{directive_name}\" "
        f"[position: ')//position//']',verbosity=verbosityLevelWorking)\n"
    )
    state_restore_code += (
        f"call displayIndent(var_str('restoring state for \"{directive_name}\" "
        f"[position: ')//position//']',verbosity=verbosityLevelWorking)\n"
    )
    state_store_code   += "select type (self)\n"
    state_restore_code += "select type (self)\n"

    for non_abstract in non_abstract_classes:
        state_store = {
            'staticVariables':      [],
            'outputCode':           '',
            'inputCode':            '',
            'excludes':             [],
            'hasCustomStateStore':  False,
            'hasCustomStateRestore': False,
        }
        state_store_code   += f"type is ({non_abstract['name']})\n"
        state_restore_code += f"type is ({non_abstract['name']})\n"

        # NOTE (issue #695): a recursive="yes" class's shim is now a distinct
        # generated type (<name>Recursive) that overrides stateStore/stateRestore
        # to forward to recursiveSelf, so the concrete class no longer needs an
        # isRecursive short-circuit here.

        state_store_code   += (
            "if (self%stateOperationID == stateOperationID) then\n")
        state_store_code   += (
            " call displayUnindent('skipping - already stored',"
            "verbosity=verbosityLevelWorking)\n")
        state_store_code   += " return\n"
        state_store_code   += "end if\n"
        state_store_code   += "self%stateOperationID=stateOperationID\n"
        state_restore_code += (
            "if (self%stateOperationID == stateOperationID) then\n")
        state_restore_code += (
            " call displayUnindent('skipping - already restored',"
            "verbosity=verbosityLevelWorking)\n")
        state_restore_code += " return\n"
        state_restore_code += "end if\n"
        state_restore_code += "self%stateOperationID=stateOperationID\n"
        state_store_code   += (
            f" call displayMessage('object type \"{non_abstract['name']}\"',"
            "verbosity=verbosityLevelWorking)\n")
        state_restore_code += (
            f" call displayMessage('object type \"{non_abstract['name']}\"',"
            "verbosity=verbosityLevelWorking)\n")

        explicit_names_found = []

        cls = non_abstract
        while cls is not None:
            node = (cls.get('tree') or {}).get('firstChild')
            while node is not None and (
                    node.get('type') != 'type'
                    or node.get('name') != cls.get('name')):
                node = node.get('sibling')
            if node is None:
                break

            # Variables to exclude from state store/restore.
            exclude_block = (
                (cls.get('stateStorable') or {}).get('exclude')
                if isinstance(cls.get('stateStorable'), dict) else None)
            if isinstance(exclude_block, dict) \
                    and 'variables' in exclude_block:
                state_store['excludes'] = [
                    v.strip() for v in re.split(
                        r'\s*,\s*', exclude_block['variables'])
                    if v.strip()
                ]
            else:
                state_store['excludes'] = []

            # Walk declarations (descending into `contains` blocks).
            walker = node.get('firstChild')
            while walker is not None:
                if walker.get('type') == 'declaration':
                    state_store_variables(
                        state_stores, state_store, cls,
                        walker.get('declarations'),
                        explicit_names_found, state_storables)
                walker = (walker.get('firstChild')
                          if walker.get('type') == 'contains'
                          else walker.get('sibling'))

            # Linked-list store/restore.
            ll_input, ll_output, ll_module = state_store_linked_list(
                cls, non_abstract, state_linked_list_variables)
            state_store['inputCode']  += ll_input
            state_store['outputCode'] += ll_output
            if ll_module:
                state_stores['stateStoreModules'][ll_module] = True

            # Explicit stateStore function (at class level).
            ss_block = non_abstract.get('stateStore') or {}
            if (isinstance(ss_block, dict)
                    and isinstance(ss_block.get('stateStore'), dict)
                    and 'restore' in ss_block['stateStore']):
                state_stores['explicitFunctionsFound'] = True
            explicit_in, explicit_out, explicit_modules = (
                state_store_explicit_function(non_abstract))
            state_store['inputCode']  += explicit_in
            state_store['outputCode'] += explicit_out
            for module in sorted(explicit_modules.keys()):
                state_stores['stateStoreModules'][module]   = True
                state_stores['stateRestoreModules'][module] = True

            cls = (None if cls.get('extends') == directive['name']
                   else classes.get(cls.get('extends')))

        # Directive-level excludes (NB: this deliberately resets excludes
        # AFTER the per-class walk — so any base-class <data> processed below
        # uses the directive-level excludes list).
        ss_block = directive.get('stateStorable') or {}
        if isinstance(ss_block, dict) \
                and isinstance(ss_block.get('exclude'), dict) \
                and 'variables' in ss_block['exclude']:
            state_store['excludes'] = [
                v.strip() for v in re.split(
                    r'\s*,\s*', ss_block['exclude']['variables'])
                if v.strip()
            ]
        else:
            state_store['excludes'] = []

        # Base-class <data> declarations.
        for data in as_array(directive.get('data')):
            declaration_source = None
            if isinstance(data, dict):
                if data.get('scope') == 'self':
                    declaration_source = data.get('content')
            else:
                declaration_source = data
            if declaration_source is None:
                continue
            declaration = parse_declaration(declaration_source)
            if declaration is None:
                raise RuntimeError(
                    "process_function_class: unable to parse variable "
                    "declaration")
            state_store_variables(
                state_stores, state_store, None,
                [declaration], explicit_names_found, state_storables)

        # Validate every variable named in <stateStorable><functionClass
        # variables="..."/> was found.
        cls = non_abstract
        while cls is not None:
            ss = cls.get('stateStorable') or {}
            fc_block = ss.get('functionClass') if isinstance(ss, dict) else None
            if isinstance(fc_block, dict) and 'variables' in fc_block:
                explicit_lc = {v.lower() for v in explicit_names_found}
                for variable in re.split(
                        r'\s*,\s*', fc_block['variables']):
                    variable = variable.strip()
                    if not variable:
                        continue
                    if variable.lower() not in explicit_lc:
                        raise RuntimeError(
                            f"Error: unable to find variable "
                            f"'{variable}' marked as state storable in "
                            f"class '{cls['name']}'"
                        )
            cls = (None if cls.get('extends') == directive['name']
                   else classes.get(cls.get('extends')))

        if state_store['staticVariables']:
            state_stores['stateFileUsed'] = True

        if state_store['hasCustomStateStore']:
            loc_node = non_abstract.get('node') or {}
            loc_expr = location(loc_node, loc_node.get('line', 0))
            state_store_code += (
                " call Error_Report('custom state store function exists - "
                f"this should not happen'//{loc_expr})\n"
            )
            state_stores['stateStoreModules']['Error'] = True
        else:
            for v in state_store['staticVariables']:
                state_stores['labelUsed'] = True
                state_store_code += (
                    " if (displayVerbosity() >= verbosityLevelWorking) then\n"
                )
                state_store_code += "   write (label,'(i16)') 0\n"
                state_store_code += (
                    f"  call displayMessage('storing \"{v}\" with size '//"
                    "trim(adjustl(label))//' bytes')\n"
                )
                state_store_code += " end if\n"
            if state_store['staticVariables']:
                state_store_code += (
                    " write (stateFile) "
                    + ", &\n  & ".join(
                        f"self%{v}" for v in state_store['staticVariables'])
                    + "\n"
                )
            state_store_code += state_store['outputCode']

        if state_store['hasCustomStateRestore']:
            loc_node = non_abstract.get('node') or {}
            loc_expr = location(loc_node, loc_node.get('line', 0))
            state_restore_code += (
                " call Error_Report('custom state restore function exists - "
                f"this should not happen'//{loc_expr})\n"
            )
            state_stores['stateRestoreModules']['Error'] = True
        else:
            for v in state_store['staticVariables']:
                state_restore_code += (
                    f" call displayMessage('restoring \"{v}\"',"
                    "verbosity=verbosityLevelWorking)\n"
                )
            if state_store['staticVariables']:
                state_restore_code += (
                    " read (stateFile) "
                    + ", &\n  & ".join(
                        f"self%{v}" for v in state_store['staticVariables'])
                    + "\n"
                )
            state_restore_code += state_store['inputCode']

    state_store_code   += "end select\n"
    state_store_code   += (
        "call displayUnindent('done',verbosity=verbosityLevelWorking)\n")
    state_store_code   += "return\n"
    state_restore_code += "end select\n"
    state_restore_code += (
        "call displayUnindent('done',verbosity=verbosityLevelWorking)\n")
    state_restore_code += "return\n"

    if not state_stores['gslStateFileUsed']:
        output_unused_variables.append('gslStateFile')
        input_unused_variables .append('gslStateFile')
    if not state_stores['stateFileUsed']:
        output_unused_variables.append('stateFile')
        input_unused_variables .append('stateFile')

    prefix_store = ''
    if state_stores['rankMaximum'] > 0:
        prefix_store = (
            ' integer :: '
            + ', '.join(
                f"i{i}" for i in range(1, state_stores['rankMaximum'] + 1))
            + '\n'
        )
    if output_unused_variables:
        prefix_store += (
            ' !$GLC attributes unused :: '
            + ', '.join(output_unused_variables) + '\n'
        )
    state_store_code = prefix_store + state_store_code

    prefix_restore = ''
    if state_stores['rankMaximum'] > 0:
        prefix_restore = (
            ' integer :: '
            + ', '.join(
                f"i{i}" for i in range(1, state_stores['rankMaximum'] + 1))
            + '\n'
        )
    if input_unused_variables:
        prefix_restore += (
            ' !$GLC attributes unused :: '
            + ', '.join(input_unused_variables) + '\n'
        )
    state_restore_code = prefix_restore + state_restore_code

    if state_stores['allocatablesFound']:
        prefix = ''
        if state_stores['dimensionalsFound']:
            prefix += (
                'integer(c_size_t), allocatable, dimension(:) '
                ':: storedShape\n')
        prefix += ('logical                                      '
                   ':: wasAllocated\n')
        state_restore_code = prefix + state_restore_code
    if state_stores['explicitFunctionsFound']:
        state_restore_code = 'logical :: wasAssociated\n' + state_restore_code

    if state_stores['labelUsed']:
        state_store_code = ' character(len=16) :: label\n' + state_store_code

    state_store_code   = ' integer(c_size_t) :: position\n' + state_store_code
    state_restore_code = ' integer(c_size_t) :: position\n' + state_restore_code
    state_store_code   = (
        _format_variable_definitions(state_linked_list_variables)
        + state_store_code)
    state_restore_code = (
        _format_variable_definitions(state_linked_list_variables)
        + state_restore_code)

    args = [
        'integer, intent(in   ) :: stateFile',
        'type(c_ptr), intent(in   ) :: gslStateFile',
        'integer(c_size_t), intent(in   ) :: stateOperationID',
    ]
    methods['stateStore'] = {
        'description': 'Store the state of this object to file.',
        'type':        'void',
        'pass':        'yes',
        'argument':    list(args),
        'code':        "call self%stateStore_(stateFile,gslStateFile,"
                       "stateOperationID)",
    }
    methods['stateStore_'] = {
        'description': 'Store the state of this object to file.',
        'type':        'void',
        'pass':        'yes',
        'modules':     ' '.join(
            sorted(state_stores['stateStoreModules'].keys())),
        'argument':    list(args),
        'code':        state_store_code,
    }
    methods['stateRestore'] = {
        'description': 'Restore the state of this object from file.',
        'type':        'void',
        'pass':        'yes',
        'argument':    list(args),
        'code':        "call self%stateRestore_(stateFile,gslStateFile,"
                       "stateOperationID)",
    }
    methods['stateRestore_'] = {
        'description': 'Restore the state of this object from file.',
        'type':        'void',
        'pass':        'yes',
        'modules':     ' '.join(
            sorted(state_stores['stateRestoreModules'].keys())),
        'argument':    list(args),
        'code':        state_restore_code,
    }


def _generate_type_definition(directive, methods, pre, node):
    """Emit the `type, extends(...) :: <name>Class` block with methods
    table, generic interfaces, destructor binding, and the module-level
    data declarations.
    """
    parent = node['parent']
    directive_name = directive['name']

    set_visibility(parent, directive_name + 'Class', 'public')
    set_visibility(parent, directive_name,           'public')

    extends = directive.get('extends', 'functionClass')

    pre['content'] += f'   type, extends({extends}) :: {directive_name}Class\n'
    pre['content'] += '    private\n'
    pre['content'] += '    integer(c_size_t) :: stateOperationID=0\n'
    pre['content'] += (
        f'    class({directive_name}Class), public, pointer :: '
        'copiedSelf => null()\n'
    )
    add_uses(parent, {
        'moduleUse': {
            'Function_Classes': {'intrinsic': False, 'all': True},
            'ISO_C_Binding':    {'intrinsic': True,  'all': True},
        },
        'moduleOrder': ['Function_Classes', 'ISO_C_Binding'],
    })

    # Per-directive `<data>` entries (self-scope go in the type, module-
    # scope deferred until after `end type`).
    for data in as_array(directive.get('data')):
        if isinstance(data, dict):
            if data.get('scope', 'self') == 'self':
                pre['content'] += (data.get('content') or '') + '\n'
        else:
            pre['content'] += str(data) + '\n'

    # contains + methods metadata.
    pre['content'] += '    contains\n'
    pre['content'] += '    !![\n'
    pre['content'] += '    <methods>\n'

    generics = {}
    for method_name in sorted(methods):
        if method_name == 'destructor':
            continue
        method = methods[method_name]
        argument_list = _method_arguments_to_latex(method)

        pre['content'] += f'     <method method="{method_name}">\n'
        pre['content'] += '      <description>\n'
        description = method.get('description') or ''
        for line in description.split('\n'):
            pre['content'] += '       ' + _xml_escape(line) + '\n'
        pre['content'] += '      </description>\n'
        pre['content'] += '     </method>\n'

        for generic in as_array(directive.get('generic')):
            if not isinstance(generic, dict):
                continue
            generic_methods = list(as_array(generic.get('method')))
            if method_name in generic_methods:
                g_entry = generics.setdefault(generic['name'], {
                    'name':         generic['name'],
                    'description':  [],
                    'argumentList': [],
                })
                g_entry['type'] = method.get('type')
                g_entry['description'].append(description)
                g_entry['argumentList'].append(argument_list)

    for gname in sorted(generics):
        g = generics[gname]
        pre['content'] += f'    <method method="{g["name"]}">\n'
        pre['content'] += '       <description>\n'
        pre['content'] += '        ' + ' | '.join(g['description']) + '\n'
        pre['content'] += '       </description>\n'
        pre['content'] += '    </method>\n'

    pre['content'] += '    </methods>\n'
    pre['content'] += '    !!]\n'

    # Type-bound procedure bindings (skip destructor, assignment(=)).
    for method_name in sorted(methods):
        if method_name in ('destructor', 'assignment(=)'):
            continue
        method = methods[method_name]
        if 'function' in method:
            function_name = method['function']
        else:
            function_name = (directive_name + _ucfirst(method_name)
                             + ('' if 'code' in method else '__'))
        pre['content'] += (
            f'    procedure :: {method_name} => {function_name}\n'
        )
    pre['content'] += f'procedure :: {directive_name}Assignment\n'

    # Generic interfaces.
    for generic in as_array(directive.get('generic')):
        if not isinstance(generic, dict):
            continue
        gmethods = list(as_array(generic.get('method')))
        pre['content'] += (
            f'    generic :: {generic["name"]} => {", ".join(gmethods)}\n'
        )
    pre['content'] += (
        f'    generic :: assignment(=) => {directive_name}Assignment\n'
    )

    if 'destructor' in methods:
        pre['content'] += f'   final :: {directive_name}Destructor\n'
    pre['content'] += f'   end type {directive_name}Class\n\n'

    # Module-scope `<data scope="module">` entries (after `end type`).
    for data in as_array(directive.get('data')):
        if isinstance(data, dict) and data.get('scope') == 'module':
            content = data.get('content') or ''
            pre['content'] += content + '\n'
            if data.get('threadprivate') == 'yes':
                m = re.search(r'::\s*(.*)$', content)
                if m:
                    vars_text = m.group(1)
                    names = [re.sub(r'\s*=.*$', '', v).strip()
                             for v in re.split(r'\s*,\s*', vars_text)
                             if v.strip()]
                    pre['content'] += (
                        f'   !$omp threadprivate({",".join(names)})\n'
                    )

    # State variable for input-parameter validation.
    pre['content'] += f'   integer :: {directive_name}DsblVldtn=0\n'
    pre['content'] += f'   !$omp threadprivate({directive_name}DsblVldtn)\n'
    pre['content'] += (
        f'   !$GLC ignore unused :: {directive_name}DsblVldtn\n'
    )


def _method_arguments_to_latex(method):
    """Render each of a method's argument declarations to the LaTeX
    notation FunctionClass embeds into the <methods>/<description>
    metadata.  We lean on our existing `parse_declaration` rather than
    indexing into INTRINSIC_DECLARATIONS directly.
    """
    args = list(as_array(method.get('argument') or []))
    if not args:
        return ''
    parts = []
    for arg in args:
        decl = parse_declaration(arg)
        if decl is None:
            continue
        intrinsic = decl['intrinsic']
        type_val  = decl.get('type')
        attrs     = decl.get('attributes') or []
        for variable in decl.get('variables') or []:
            piece = (
                r'\textcolor{red}{\textless ' + _latex_encode(intrinsic)
            )
            if type_val is not None:
                piece += _latex_encode(type_val)
            piece += r'\textgreater} ' + _latex_encode(variable)
            for attr in attrs:
                if attr == 'intent(in)':
                    piece += r'\argin'
                elif attr == 'intent(out)':
                    piece += r'\argout'
                elif attr == 'intent(inout)':
                    piece += r'\arginout'
            parts.append(piece)
    return ','.join(parts)


def _generate_constructor(directive, classes_ordered, non_abstract_classes,
                          pre, post, node, tree):
    """Emit the `<name>CnstrctrPrmtrs` parameter-driven constructor +
    associated interface + the recursive-build state variables when any
    class opts into recursion.
    """
    directive_name = directive['name']
    parent = node['parent']
    loc_expr = location(node, node.get('line', 0))

    pre['content'] += f'   interface {directive_name}\n'
    pre['content'] += f'    module procedure {directive_name}CnstrctrPrmtrs\n'
    pre['content'] += f'   end interface {directive_name}\n'

    allow_recursion = any(
        c.get('recursive') == 'yes' for c in classes_ordered
    )
    if allow_recursion:
        # A recursive="yes" class detects a bounded construction cycle -- a re-entrant build that re-discovers the
        # node currently under construction -- by querying the shared object-build stack (see the short-circuit
        # below), rather than via the per-family thread-private RecursiveBuildNode/RecursiveBuildObject module
        # variables used previously. The object under construction is recorded on the stack entry by the factory
        # (Input_Parameters_Build_Stack_Object_Set) after allocation but before dispatch. See issue #695.
        add_uses(parent, {
            'moduleUse': {
                'Input_Parameters': {'intrinsic': False, 'all': True},
            },
            'moduleOrder': ['Input_Parameters'],
        })

    # Append directive name to the per-file `.p` parameter-names file. This must
    # be written to $BUILDPATH/<path relative to source/>.p, mirroring the
    # hierarchical source tree, because parameterDependencies.py reads it from
    # exactly that location. The tree node's `name` is only the basename, so use
    # its `source` (the full path passed to parse_file); using the basename
    # would collapse every `_class.F90` into one colliding `_class.p` and lose
    # all functionClass parameter names from `knownParameterNames`.
    if tree.get('type') == 'file':
        build_path  = os.environ.get('BUILDPATH')
        source_path = tree.get('source') or tree.get('name', '')
        rel         = source_path
        exec_path   = os.environ.get('GALACTICUS_EXEC_PATH')
        if exec_path:
            src_root = os.path.abspath(os.path.join(exec_path, 'source'))
            ap       = os.path.abspath(source_path)
            if ap.startswith(src_root + os.sep):
                rel = os.path.relpath(ap, src_root)
        rel = re.sub(r'^(?:\./)?source/', '', rel)
        m = re.match(r'(.+)\.F90$', rel)
        if build_path and m:
            p_path = os.path.join(build_path, m.group(1) + '.p')
            os.makedirs(os.path.dirname(p_path) or '.', exist_ok=True)
            with open(p_path, 'a') as fh:
                fh.write(directive_name + '\n')

    recursive_prefix = 'recursive ' if allow_recursion else ''
    post['content'] += (
        f'   {recursive_prefix}function '
        f'{directive_name}CnstrctrPrmtrs'
        '(parameters,copyInstance,parameterName) result(self)\n'
    )
    post['content'] += '      !!{\n'
    post['content'] += (
        f'      Return a pointer to a newly created \\mono{{{directive_name}}} '
        'object as specified by the provided parameters.\n'
    )
    post['content'] += '      !!}\n'
    post['content'] += (
        '      use :: Input_Parameters  , only : inputParameter         , '
        'inputParameters        , Input_Parameters_Build_Stack_Push, '
        'Input_Parameters_Build_Stack_Pop, Input_Parameters_Build_Stack_Object_Set, '
        'Input_Parameters_Build_Stack_Recursive_Object\n'
        '      use :: Locks  , only : ompLock\n'
        '      use :: Error  , only : Error_Report\n'
        '      use :: ISO_Varying_String, only : varying_string         , '
        'char           , trim, operator(//), operator(==), assignment(=)\n'
        '      implicit none\n'
        f'      class    ({directive_name}Class), pointer :: self\n'
        '      type     (inputParameters), intent(inout)           :: parameters\n'
        '      integer                   , intent(in   ), optional :: copyInstance\n'
        '      character(len=*          ), intent(in   ), optional :: parameterName\n'
        '      type     (inputParameters)                          :: subParameters\n'
    )
    if 'default' in directive:
        post['content'] += (
            '      type     (inputParameter ), pointer                 :: parameterNode\n'
            '      type     (ompLock        ), save                    :: addLock\n'
            '      logical                   , save                    :: addLockInitialized=.false.\n'
            '      logical                                             :: needLock\n'
        )
    if allow_recursion:
        # Weak pointer to the object already under construction for the re-discovered node, returned by the
        # object-build-stack query used to detect a bounded construction cycle (issue #695). A recursive="yes"
        # class implies a <default>, so `parameterNode` above is always declared alongside this.
        post['content'] += (
            '      class    (*               ), pointer                 :: recursiveObject\n'
        )
    post['content'] += (
        '      type     (varying_string )                          '
        ':: message      , instanceName, parameterName_\n'
        '      integer                                             '
        ':: copyInstance_\n\n'
        '      if (present(parameterName)) then\n'
        '        parameterName_=parameterName\n'
        '      else\n'
        f"        parameterName_='{directive_name}'\n"
        '      end if\n'
        '      if (present(copyInstance)) then\n'
        '        copyInstance_=copyInstance\n'
        '      else\n'
        '        copyInstance_=1\n'
        '      end if\n'
    )

    if 'default' in directive:
        default = directive['default']
        default_uc = _ucfirst(default)
        # Find the matching default class entry.
        target_name = directive_name + default_uc
        match = next(
            (c for c in non_abstract_classes if c['name'] == target_name),
            None,
        )
        post['content'] += (
            '      if (.not.addLockInitialized) then\n'
            f'      !$omp critical (addLockInitialize{default_uc})\n'
            '          if (.not.addLockInitialized) then\n'
            '          addLock=ompLock()\n'
            '          addLockInitialized=.true.\n'
            '      end if\n'
            f'      !$omp end critical (addLockInitialize{default_uc})\n'
            '      end if\n'
            '      needLock=.not.addLock%ownedByThread()\n'
            '      if (needLock) call addLock%set()\n'
        )
        post['content'] += (
            f"      if (parameterName_ == '{directive_name}' "
            f".and. copyInstance_ == 1 "
            f".and. .not.parameters%isPresent(char(parameterName_))) then\n"
            f"        call parameters%addParameter"
            f"('{directive_name}','{default}')\n"
            f"        parameterNode => parameters%node"
            f"('{directive_name}',requireValue=.true.)\n"
        )
        post['content'] += (
            '        subParameters=parameters%subParameters(char(parameterName_))\n'
            f'        allocate({target_name} :: self)\n'
        )
        post['content'] += (
            '        call Input_Parameters_Build_Stack_Push(subParameters%parameters,'
            f"'{directive_name}',{'.true.' if allow_recursion else '.false.'},{loc_expr})\n"
        )
        if match is not None and match.get('recursive') == 'yes':
            # Record the in-progress object on the build-stack entry so a re-entrant build that re-discovers this
            # node can retrieve it and return a shim wired to it (issue #695).
            post['content'] += (
                '        call Input_Parameters_Build_Stack_Object_Set(self)\n'
            )
        post['content'] += (
            '        select type (self)\n'
            f'          type is ({target_name})\n'
            f'            self={target_name}(subParameters)\n'
            '         end select\n'
            '        call Input_Parameters_Build_Stack_Pop()\n'
        )
        post['content'] += '         call parameterNode%objectSet(self)\n'
        post['content'] += '      else\n'

    if allow_recursion:
        # Detect a bounded construction cycle: query the object-build stack for an in-progress object being built
        # for this same node and class. If found, this build has re-discovered the node currently under
        # construction, so return the generated shim type -- carrying only a weak recursiveSelf pointer back to the
        # in-progress object -- instead of building again. See issue #695.
        shim_type = directive_name + 'Recursive'
        post['content'] += (
            '        parameterNode => parameters%node'
            '(char(parameterName_),requireValue=.true.)\n'
            f'        recursiveObject => Input_Parameters_Build_Stack_Recursive_Object'
            f"(parameterNode,'{directive_name}')\n"
            '        if (associated(recursiveObject)) then\n'
            f'           allocate({shim_type} :: self)\n'
            '           select type (self)\n'
            f'           type is ({shim_type})\n'
            '              select type (recursiveObject)\n'
            f'              class is ({directive_name}Class)\n'
            '                 self%recursiveSelf => recursiveObject\n'
            '              end select\n'
            '           end select\n'
            '           if (needLock) call addLock%unset()\n'
            '           return\n'
            '        end if\n'
        )

    # Main build ladder.
    post['content'] += (
        '      call parameters%value(char(parameterName_),'
        'instanceName,copyInstance=copyInstance_)\n'
        '      subParameters=parameters%subParameters'
        '(char(parameterName_),copyInstance=copyInstance_)\n'
        '      call Input_Parameters_Build_Stack_Push(subParameters%parameters,'
        f"'{directive_name}',{'.true.' if allow_recursion else '.false.'},{loc_expr})\n"
        '      select case (char(instanceName))\n'
    )
    for c in non_abstract_classes:
        name = _short_name(c['name'], directive_name)
        post['content'] += f"     case ('{name}')\n"
        post['content'] += f"        allocate({c['name']} :: self)\n"
        if c.get('recursive') == 'yes':
            # Record the in-progress object on the (already-pushed) build-stack entry so a re-entrant build that
            # re-discovers this node can retrieve it and return a shim wired to it (issue #695).
            post['content'] += (
                '        call Input_Parameters_Build_Stack_Object_Set(self)\n'
            )
        post['content'] += (
            '        select type (self)\n'
            f"          type is ({c['name']})\n"
            f"            self={c['name']}(subParameters)\n"
            '         end select\n'
        )
    post['content'] += '      case default\n'
    post['content'] += (
        "         message='Unrecognized type \"'//trim(instanceName)//'\" "
        "Available options are:'\n"
    )
    class_names = sorted([c['name'] for c in non_abstract_classes])
    for cn in class_names:
        sn = _short_name(cn, directive_name)
        post['content'] += (
            f"         message=message//char(10)//'   -> {sn}'\n"
        )
    post['content'] += (
        f'         call Error_Report(message//{loc_expr})\n'
        '      end select\n'
        '      call Input_Parameters_Build_Stack_Pop()\n'
    )
    if 'default' in directive:
        post['content'] += '      end if\n'
        post['content'] += '      if (needLock) call addLock%unset()\n'
    post['content'] += '      return\n'
    post['content'] += (
        f'   end function {directive_name}CnstrctrPrmtrs\n\n'
    )


def _generate_class_submodules(directive, classes_ordered, non_abstract_classes,
                               pre, code_content, node):
    """Populate `code_content['submodule'][<class>]` and append module-level
    interfaces/declarations for every class in the hierarchy.

    Walks
    each class's parsed tree, partitioning its content between:
      - the per-class submodule (interior functions, private declarations,
        submodule-scoped variables),
      - the parent module's interface list (type definitions, type-bound
        procedure interfaces, visibility statements, public declarations,
        module-scoped variables),
      - module-level module-use statements (forwarded via `add_uses`).
    """
    directive_name = directive['name']
    parent = node['parent']
    build_path = os.environ.get('BUILDPATH') or ''

    for class_record in classes_ordered:
        set_visibility(parent, class_record['type'], 'public')

        # Submodule output path: swap the source file's `.F90` for `.p.F90`
        # and place it under BUILDPATH, preserving the file's path *relative
        # to the source tree* so the hierarchical layout is mirrored (e.g.
        # `source/cosmology/parameters/simple.F90` ->
        # `$BUILDPATH/cosmology/parameters/simple.p.F90`).  Using only the
        # basename would flatten every implementation into `$BUILDPATH/` and
        # collide / mis-locate them once the source tree is hierarchical.
        class_file = class_record.get('file') or ''
        exec_path  = os.environ.get('GALACTICUS_EXEC_PATH')
        rel = None
        if exec_path:
            source_root = os.path.abspath(os.path.join(exec_path, 'source'))
            abs_class   = os.path.abspath(class_file)
            if abs_class.startswith(source_root + os.sep):
                rel = os.path.relpath(abs_class, source_root)
        if rel is None:
            rel = os.path.basename(class_file)
        rel = re.sub(r'\.F90$', '.p.F90', rel)
        fn  = os.path.join(build_path, rel) if build_path else rel
        if build_path:
            os.makedirs(os.path.dirname(fn), exist_ok=True)

        block = {
            'fileName':    fn,
            'preContains':  [_code_node_submodule()],
            'postContains': [_code_node_submodule()],
            'interfaces':   [],
        }
        code_content['submodule'][class_record['type']] = block

        submodule_pre  = block['preContains']
        submodule_post = block['postContains']

        module_scoped     = []
        module_symbols    = []
        module_use_nodes  = []

        class_tree = class_record.get('tree')
        if class_tree is None:
            continue

        class_node = class_tree.get('firstChild')
        contained = False
        while class_node is not None:
            if class_node.get('type') == 'contains':
                # `contains` is a self-closing sibling marker in our parse —
                # post-contains procedures are siblings of it, not children.
                contained = True
                class_node = class_node.get('sibling')
                continue
            if contained:
                submodule_post.append(class_node)
                if class_node.get('type') in ('function', 'subroutine'):
                    # Only emit an interface if this function's lowercased
                    # name is in the interfaces list we accumulated from
                    # the pre-contains walk.
                    name_lc = (class_node.get('name') or '').lower()
                    if name_lc in block['interfaces']:
                        _emit_submodule_function_interface(
                            class_node, code_content, module_symbols,
                            module_use_nodes,
                        )
            else:
                _handle_pre_contains_node(
                    class_node, directive, code_content, block, parent,
                    module_scoped, module_symbols, module_use_nodes,
                    submodule_pre,
                )
            class_node = class_node.get('sibling')

        # Filter imported symbols on each captured moduleUse to just those
        # required for module-scope interfaces, then hand them to add_uses.
        module_symbols = sorted({s.lower() for s in module_symbols})
        for mu_node in module_use_nodes:
            uses = mu_node.get('moduleUse') or {}
            for module_name in sorted(uses.keys()):
                kept_entries = []
                for entry in _as_entry_list(uses[module_name]):
                    if entry.get('all'):
                        kept_entries.append(entry)
                        continue
                    only = entry.get('only', {}) or {}
                    kept = {
                        sym: flag for sym, flag in only.items()
                        if sym.lower() in module_symbols
                    }
                    if kept:
                        entry['only'] = kept
                        kept_entries.append(entry)
                if kept_entries:
                    uses[module_name] = kept_entries
                else:
                    del uses[module_name]
            if uses:
                add_uses(parent, {
                    'moduleUse':   uses,
                    'moduleOrder': sorted(uses.keys()),
                })

    # Emit source-digest bindings for every non-abstract class.
    for nac in non_abstract_classes:
        pre['content'] += _source_digest_binding(nac['name'])

    # ISO_C_Binding is required for those bindings.
    add_uses(parent, {
        'moduleUse': {
            'ISO_C_Binding': {
                'intrinsic': True,
                'only':      {'C_Char': True},
            },
        },
        'moduleOrder': ['ISO_C_Binding'],
    })


def _code_node_submodule():
    """Code node factory used specifically by generate_class_submodules's
    per-class pre/postContains slots (source attribution differs slightly
    from the top-level one).
    """
    return {
        'type':       'code',
        'content':    '',
        'parent':     None,
        'firstChild': None,
        'sibling':    None,
        'source':
            'Galacticus.Build.SourceTree.Process.FunctionClass'
            '.process_function_class()',
        'line':       1,
    }


def _emit_submodule_function_interface(class_node, code_content,
                                       module_symbols, module_use_nodes):
    """For a contained function/subroutine whose name needs a module
    interface, build the interface nodes and append them to
    `code_content['module']['interfaces']`.  Also capture the type names
    used in the interface into `module_symbols` so add_uses can later
    import them.
    """
    interface_opener = _code_node_submodule()
    interface_opener['content'] = 'interface\n'
    interface_opener['content'] += 'module ' + (class_node.get('opener') or '')
    interface_nodes = [interface_opener]

    opener = class_node.get('opener') or ''
    m = re.search(r'result\s*\(\s*([a-zA-Z0-9_]+)\s*\)', opener)
    return_name = m.group(1) if m else class_node.get('name')

    function_node = class_node.get('firstChild')
    while function_node is not None:
        if function_node.get('type') == 'declaration':
            local_declarations = []
            for declaration in function_node.get('declarations') or []:
                variables_lc = [
                    v.lower() for v in (declaration.get('variables') or [])
                ]
                attributes = declaration.get('attributes') or []
                has_intent = any(
                    re.match(r'intent\s*\(\s*(in|out|inout)\s*\)', a,
                             re.IGNORECASE)
                    for a in attributes
                )
                external = any(
                    re.match(r'external', a, re.IGNORECASE)
                    for a in attributes
                )
                if ((return_name or '').lower() in variables_lc
                        or has_intent):
                    # Function-result + intent() declarations always go into
                    # the interface.
                    if (function_node.get('firstChild') is None
                            or function_node['firstChild'].get('type') != 'code'):
                        raise RuntimeError(
                            'process_function_class: expected a "code" '
                            'node as first child of function declaration'
                        )
                    iface_decl = _code_node_submodule()
                    iface_decl['source'] = function_node['firstChild'].get('source')
                    iface_decl['line']   = declaration.get('line', 1)
                    iface_decl['content'] = _format_variable_definitions(
                        [declaration])
                    interface_nodes.append(iface_decl)
                    _capture_declaration_type(declaration, module_symbols)
                elif external or declaration.get('intrinsic') == 'procedure':
                    # External/procedure declarations split into module
                    # vs submodule scope based on whether the variable is
                    # a dummy argument of the enclosing unit.
                    _split_external_or_procedure(
                        declaration, function_node, local_declarations,
                        interface_nodes, module_symbols,
                    )
                else:
                    local_declarations.append(declaration)

            function_node['declarations'] = local_declarations
            build_declarations(function_node)
        elif function_node.get('type') == 'moduleUse':
            first_child = function_node.get('firstChild') or {}
            content_node = _code_node_submodule()
            content_node['content'] = first_child.get('content') or ''
            from copy import deepcopy
            module_use_node = {
                'type':       'moduleUse',
                'parent':     None,
                'sibling':    None,
                'firstChild': content_node,
                'source':
                    'Galacticus.Build.SourceTree.Process.FunctionClass'
                    '.process_function_class()',
                'line':       1,
                'moduleUse':  deepcopy(function_node.get('moduleUse') or {}),
            }
            content_node['parent'] = module_use_node
            module_use_nodes.append(module_use_node)
        function_node = function_node.get('sibling')

    interface_closer = _code_node_submodule()
    interface_closer['content'] = (
        (class_node.get('closer') or '') + 'end interface\n'
    )
    class_node['opener'] = (
        f'module procedure {class_node.get("name")}\n'
    )
    class_node['closer'] = (
        f'end procedure {class_node.get("name")}\n'
    )
    interface_nodes.append(interface_closer)
    code_content['module']['interfaces'].extend(interface_nodes)


def _capture_declaration_type(declaration, module_symbols):
    """Extract the type name from a declaration (stripping `len=…=` prefixes)
    and append it to module_symbols unless purely numeric.
    """
    t = declaration.get('type')
    if t is None:
        return
    t = re.sub(r'^.+=', '', t)
    t = re.sub(r'\s', '', t)
    if re.match(r'^\d+$', t):
        return
    module_symbols.append(t)


def _split_external_or_procedure(declaration, function_node,
                                 local_declarations, interface_nodes,
                                 module_symbols):
    """Split an `external` or `procedure` declaration between module scope
    (dummy-argument matches → interface) and submodule scope (every other
    variable → local declarations).
    """
    from Fortran.Utils import UNIT_OPENERS
    from copy import deepcopy
    parent_unit = function_node.get('parent') or {}
    parent_type = parent_unit.get('type')
    spec        = UNIT_OPENERS.get(parent_type) or {}
    regex       = spec.get('regex')
    arg_idx     = spec.get('arguments')
    arguments = []
    if regex is not None and arg_idx is not None:
        m = regex.match(parent_unit.get('opener') or '')
        if m is not None:
            groups = m.groups()
            if arg_idx < len(groups) and groups[arg_idx]:
                arguments = [a.strip().lower()
                             for a in re.split(r'\s*,\s*', groups[arg_idx])
                             if a.strip()]

    module_scope    = deepcopy(declaration)
    submodule_scope = deepcopy(declaration)
    module_scope['variables']     = []
    module_scope['variableNames'] = []
    submodule_scope['variables']     = []
    submodule_scope['variableNames'] = []

    variables      = declaration.get('variables')      or []
    variable_names = declaration.get('variableNames')  or []
    for i, var_name in enumerate(variable_names):
        if var_name.lower() in arguments:
            module_scope['variables'].append(variables[i])
            module_scope['variableNames'].append(var_name)
        else:
            submodule_scope['variables'].append(variables[i])
            submodule_scope['variableNames'].append(var_name)

    if submodule_scope['variables']:
        local_declarations.append(submodule_scope)
    if module_scope['variables']:
        iface_decl = _code_node_submodule()
        iface_decl['content'] = _format_variable_definitions([module_scope])
        interface_nodes.append(iface_decl)
        _capture_declaration_type(declaration, module_symbols)


def _handle_pre_contains_node(class_node, directive, code_content, block,
                              parent, module_scoped, module_symbols,
                              module_use_nodes, submodule_pre):
    """Dispatch a single pre-contains class-tree node to the right output
    partition.
    """
    ntype = class_node.get('type')
    directive_name = directive['name']

    if ntype == 'scoping':
        d = class_node.setdefault('directive', {})
        inner = d.get('module') if isinstance(d, dict) else None
        if isinstance(inner, dict) and 'variables' in inner:
            module_scoped.extend(
                v.strip()
                for v in re.split(r'\s*,\s*', inner['variables'])
                if v.strip()
            )
    elif ntype == 'moduleUse':
        add_uses(parent, class_node)
    elif ntype in ('type', 'interface', directive_name):
        code_content['module']['interfaces'].append(class_node)
        if ntype == 'interface':
            inode = class_node.get('firstChild')
            while inode is not None:
                if inode.get('type') == 'moduleProcedure':
                    block['interfaces'].extend(
                        n.lower() for n in inode.get('names') or []
                    )
                inode = inode.get('sibling')
        elif ntype == 'type':
            _collect_type_bindings_and_module_symbols(
                class_node, block, module_symbols)
    elif ntype == 'visibility':
        vis = class_node.get('visibility') or {}
        public_names = sorted((vis.get('public') or {}).keys())
        block['interfaces'].extend(n.lower() for n in public_names)
        if 'private' in vis:
            del vis['private']
        update_visibilities(class_node)
        code_content['module']['interfaces'].append(class_node)
    elif ntype == 'declaration':
        _partition_declaration(class_node, code_content, module_scoped,
                               module_symbols, submodule_pre)
    else:
        submodule_pre.append(class_node)


def _collect_type_bindings_and_module_symbols(class_node, block,
                                              module_symbols):
    """For a `type` node, accumulate the names of bound procedures (→
    block['interfaces']) and the types referenced pre-contains
    (→ module_symbols).
    """
    post_contains = False
    type_node = class_node.get('firstChild')
    while type_node is not None:
        if type_node.get('type') == 'contains':
            # Self-closing sibling marker — set the flag and skip past.
            post_contains = True
            type_node = type_node.get('sibling')
            continue
        if post_contains:
            if type_node.get('type') == 'declaration':
                for declaration in type_node.get('declarations') or []:
                    for variable in declaration.get('variables') or []:
                        m = re.search(r'=>([a-z0-9_]+)', variable)
                        if m:
                            block['interfaces'].append(m.group(1))
                        else:
                            block['interfaces'].append(variable)
            elif type_node.get('type') == 'code':
                content = type_node.get('content') or ''
                m = re.search(r'final\s*::\s*([a-zA-Z0-9_]+)', content)
                if m:
                    block['interfaces'].append(m.group(1).lower())
        else:
            if type_node.get('type') == 'declaration':
                for declaration in type_node.get('declarations') or []:
                    _capture_declaration_type(declaration, module_symbols)
        type_node = type_node.get('sibling')


def _partition_declaration(class_node, code_content, module_scoped,
                           module_symbols, submodule_pre):
    """Split a `declaration` node's variables into module-scope (public or
    listed in `<scoping>`) vs submodule-scope, rebuilding both sides.
    """
    from copy import deepcopy
    submodule_declarations = []
    for declaration in class_node.get('declarations') or []:
        module_scope_used = False
        attributes = declaration.get('attributes') or []
        is_public = any(a.lower() == 'public' for a in attributes)

        if is_public:
            new_node = _new_declaration_node([declaration])
            code_content['module']['interfaces'].append(new_node)
            build_declarations(new_node)
            module_scope_used = True
        else:
            module_variables = []
            submodule_variables = []
            variables = declaration.get('variables') or []
            variable_names = declaration.get('variableNames') or []
            for i, var in enumerate(variables):
                vname = variable_names[i] if i < len(variable_names) else ''
                if any(s.lower() == vname.lower() for s in module_scoped):
                    module_variables.append((var, vname))
                else:
                    submodule_variables.append((var, vname))

            if module_variables:
                module_decl = deepcopy(declaration)
                module_decl['variables']     = [v[0] for v in module_variables]
                module_decl['variableNames'] = [v[1] for v in module_variables]
                new_node = _new_declaration_node([module_decl])
                code_content['module']['interfaces'].append(new_node)
                build_declarations(new_node)
                module_scope_used = True

            if submodule_variables:
                submodule_decl = deepcopy(declaration)
                submodule_decl['variables']     = [v[0] for v in submodule_variables]
                submodule_decl['variableNames'] = [v[1] for v in submodule_variables]
                # Strip public/private from submodule decls — submodule vars
                # are private by definition and can't carry such attributes.
                submodule_decl['attributes'] = [
                    a for a in submodule_decl.get('attributes') or []
                    if a not in ('public', 'private')
                ]
                submodule_declarations.append(submodule_decl)

        if module_scope_used:
            _capture_declaration_type(declaration, module_symbols)

    class_node['declarations'] = submodule_declarations
    build_declarations(class_node)
    submodule_pre.append(class_node)


def _new_declaration_node(declarations):
    """Build a fresh declaration node with a code firstChild ready for
    build_declarations to populate.
    """
    new_node = {
        'type':       'declaration',
        'sibling':    None,
        'parent':     None,
        'source':
            'Galacticus.Build.SourceTree.Process.FunctionClass'
            '.process_function_class()',
        'line':       1,
        'declarations': list(declarations),
    }
    new_node['firstChild'] = {
        'type':       'code',
        'content':    '',
        'sibling':    None,
        'parent':     new_node,
        'firstChild': None,
        'source':     new_node['source'],
        'line':       new_node['line'],
    }
    return new_node


def _generate_method_functions(directive, methods, post, node):
    """Emit the default-implementation body for every method that has no
    explicit `function=` override.  Each stub either runs the
    directive-supplied `<code>` block, or calls Error_Report with a
    "null method" message that names the object's concrete class at
    runtime.
    """
    directive_name = directive['name']
    loc_expr = location(node, node.get('line', 0))

    for method_name in sorted(methods):
        method = methods[method_name]
        if 'function' in method:
            continue

        arguments = list(as_array(method.get('argument') or []))
        pass_self = method.get('pass', 'yes') == 'yes'

        argument_list = ''
        argument_code = ''
        unused_variables = ['self']
        if pass_self:
            intrinsic = 'type' if method_name == 'destructor' else 'class'
            intent    = method.get('selfIntent', 'inout')
            piece = f'      {intrinsic}({directive_name}Class), intent({intent})'
            if method.get('selfTarget') == 'yes':
                piece += ', target'
            piece += ' :: self\n'
            argument_code += piece

        separator = ''
        for arg in arguments:
            # Extract variable list after `::`.
            parts = arg.split('::', 1)
            variables = parts[1].strip() if len(parts) == 2 else ''
            is_openmp = re.match(r'^\s*!\$', arg) is not None
            if is_openmp:
                argument_list += ' &\n!$ & ' + separator + variables + ' &\n& '
            else:
                argument_list += separator + variables
            argument_code += '      ' + arg + '\n'
            separator = ','
            decl = parse_declaration(arg)
            if decl is not None:
                unused_variables.extend(decl.get('variables') or [])

        if method_name == 'assignment(=)':
            non_operator_name = 'assignment'
        else:
            non_operator_name = method_name
        extension = '' if 'code' in method else '__'
        function_name = directive_name + _ucfirst(non_operator_name) + extension
        recursive_prefix = (
            'recursive ' if method.get('recursive') == 'yes' else ''
        )
        elemental_prefix = (
            'elemental ' if method.get('elemental') == 'yes' else ''
        )

        method_type = method.get('type', 'void')
        self_line = ''
        type_prefix = ''
        if method_type == 'void':
            category = 'subroutine'
        elif method_type.startswith('class'):
            category = 'function'
            self_line = (
                f'      {method_type}, pointer :: {function_name}\n'
            )
        elif method_type.startswith('type') or ',' in method_type:
            category = 'function'
            self_line = f'      {method_type} :: {function_name}\n'
        else:
            category = 'function'
            type_prefix = method_type + ' '

        post['content'] += (
            f'   {recursive_prefix}{elemental_prefix}{type_prefix}{category} '
            f'{function_name}(self'
        )
        if argument_list:
            post['content'] += ',' + argument_list
        post['content'] += ')\n'
        post['content'] += '      !!{\n'
        post['content'] += (
            f'      Default implementation of the \\mono{{{method_name}}} '
            f'method for the \\mono{{{directive_name}}} class.\n'
        )
        post['content'] += '      !!}\n'

        if 'code' in method:
            modules = method.get('modules')
            if modules is not None:
                module_list = _normalise_modules(modules)
                for mod in module_list:
                    prefix = '!$ ' if 'OMP_Lib' in mod['name'] else ''
                    only_clause = ''
                    if mod.get('only'):
                        only_clause = (
                            ', only : ' + ','.join(as_array(mod['only']))
                        )
                    post['content'] += (
                        f'      {prefix}use {mod["name"]}{only_clause}\n'
                    )
        else:
            post['content'] += '      use Error             , only : Error_Report\n'
            post['content'] += '      use ISO_Varying_String, only : char\n'
        post['content'] += '      implicit none\n'
        post['content'] += argument_code
        post['content'] += self_line
        post['content'] += (
            f'!$GLC attributes unused :: {", ".join(unused_variables)}\n'
        )

        if 'code' in method:
            body = method['code']
            body = '      ' + body.replace('\n', '\n      ')
            post['content'] += body + '\n'
        else:
            post['content'] += (
                f"      call Error_Report('this is a null method - initialize the "
                f"{directive_name} object before use and/or check that the \"'//"
                f"char(self%objectType())//'\" class implements this method'//"
                f'{loc_expr})\n'
            )
            if category == 'function':
                # Avoid warnings about unset function values.
                post['content'] += f'      {function_name}='
                set_value = None
                if method_type.startswith('class'):
                    set_value = '> null()'
                elif re.match(r'^type\s*\(\s*(.*)\s*\)', method_type):
                    m = re.match(r'^type\s*\(\s*(.*)\s*\)', method_type)
                    if ',pointer' in method_type.replace(' ', '') \
                            or ', pointer' in method_type.lower():
                        set_value = '>null()'
                    else:
                        set_value = (m.group(1) if m else '') + '()'
                elif method_type.startswith('integer'):
                    set_value = '0'
                elif method_type.startswith('double precision'):
                    set_value = '0.0d0'
                elif method_type.startswith('logical'):
                    set_value = '.false.'
                if set_value is None:
                    raise RuntimeError(
                        f"generate_method_functions: do not know how to set "
                        f"'{method_type}'"
                    )
                post['content'] += set_value + '\n'
            post['content'] += '      return\n'
        post['content'] += f'   end {category} {function_name}\n\n'


# Framework methods that the recursive shim handles explicitly (i.e. NOT by
# forwarding to recursiveSelf via the physics-method path). Everything else in
# the family's `methods` map is a physics method that the shim forwards.
_SHIM_FRAMEWORK_METHODS = {
    'stateStore', 'stateRestore', 'stateStore_', 'stateRestore_',
    'autoHook', 'deepCopy', 'deepCopy_', 'deepCopyReset', 'deepCopyFinalize',
    'descriptor', 'hashedDescriptor', 'allowedParameters', 'objectType',
}

def _shim_forward_argument_names(method):
    """Return the comma-separated actual-argument list for forwarding a call
    of `method` to recursiveSelf, reconstructed from the argument
    declarations exactly as _generate_method_functions builds its
    argument_list (so optional arguments propagate their presence)."""
    argument_list = ''
    separator = ''
    for arg in as_array(method.get('argument') or []):
        parts = arg.split('::', 1)
        variables = parts[1].strip() if len(parts) == 2 else ''
        if re.match(r'^\s*!\$', arg) is not None:
            argument_list += ' &\n!$ & ' + separator + variables + ' &\n& '
        else:
            argument_list += separator + variables
        separator = ','
    return argument_list


def _generate_recursive_shim(directive, methods, non_abstract_classes,
                             pre, post, node):
    """Emit the generated shim type `<name>Recursive` for a functionClass
    family that contains a recursive="yes" class (issue #695).

    The shim is a lightweight stand-in returned by the factory when a bounded
    construction cycle re-discovers the object currently under construction. It
    holds only a weak `recursiveSelf` back-pointer to that object and forwards
    every method to it; deepCopy, stateStore, descriptor, autoHook and
    allowedParameters are handled centrally here rather than by ~100 lines of
    hand-written per-class boilerplate.
    """
    parent = node['parent']
    directive_name = directive['name']
    shim_type = directive_name + 'Recursive'
    loc_expr = location(node, node.get('line', 0))

    set_visibility(parent, shim_type, 'public')

    physics = [
        m for m in sorted(methods)
        if m not in _SHIM_FRAMEWORK_METHODS
        and m not in ('destructor', 'assignment(=)')
    ]

    # ---- type definition ----
    pre['content'] += f'   type, extends({directive_name}Class) :: {shim_type}\n'
    pre['content'] += (
        '     !!{\n'
        f'     A generated shim for the \\refClass{{{directive_name}Class}} '
        'family. Returned by the factory when a bounded construction cycle\n'
        '     re-enters the object currently under construction; forwards every '
        'method to the (weakly-referenced) real object. See issue \\#695.\n'
        '     !!}\n'
    )
    pre['content'] += (
        f'     class({directive_name}Class), pointer :: recursiveSelf => null()\n'
    )
    pre['content'] += '     logical :: parentDeferred = .false.\n'
    pre['content'] += '    contains\n'
    for m in physics:
        fn = shim_type + _ucfirst(m)
        pre['content'] += f'     procedure :: {m} => {fn}\n'
    pre['content'] += (
        f'     procedure :: allowedParameters => {shim_type}AllowedParameters\n'
        f'     procedure :: autoHook          => {shim_type}AutoHook\n'
        f'     procedure :: deepCopy          => {shim_type}DeepCopy\n'
        f'     procedure :: deepCopyReset     => {shim_type}DeepCopyReset\n'
        f'     procedure :: deepCopyFinalize  => {shim_type}DeepCopyFinalize\n'
        f'     procedure :: descriptor        => {shim_type}Descriptor\n'
        f'     procedure :: hashedDescriptor  => {shim_type}HashedDescriptor\n'
        f'     procedure :: objectType        => {shim_type}ObjectType\n'
        f'     procedure :: stateStore        => {shim_type}StateStore\n'
        f'     procedure :: stateRestore      => {shim_type}StateRestore\n'
        f'     procedure :: isRecursiveShim   => {shim_type}IsRecursiveShim\n'
    )
    pre['content'] += f'   end type {shim_type}\n\n'

    # ---- physics-method forwarders ----
    for m in physics:
        method = methods[m]
        fn = shim_type + _ucfirst(m)
        arg_names = _shim_forward_argument_names(method)
        # Rebuild the declaration block (self + arguments) with the shim's own
        # passed-object type, matching the base binding's self intent/target so
        # the override interface is compatible.
        self_intent = method.get('selfIntent', 'inout')
        self_target = ', target' if method.get('selfTarget') == 'yes' else ''
        argument_code = (
            f'      class({shim_type}), intent({self_intent}){self_target} :: self\n'
        )
        for arg in as_array(method.get('argument') or []):
            argument_code += '      ' + arg + '\n'
        method_type = method.get('type', 'void')
        recursive_prefix = 'recursive ' if method.get('recursive') == 'yes' else ''
        elemental_prefix = 'elemental ' if method.get('elemental') == 'yes' else ''
        self_line = ''
        type_prefix = ''
        if method_type == 'void':
            category = 'subroutine'
        elif method_type.startswith('class'):
            category = 'function'
            self_line = f'      {method_type}, pointer :: {fn}\n'
        elif method_type.startswith('type') or ',' in method_type:
            category = 'function'
            self_line = f'      {method_type} :: {fn}\n'
        else:
            category = 'function'
            type_prefix = method_type + ' '
        header = f'   {recursive_prefix}{elemental_prefix}{type_prefix}{category} {fn}(self'
        if arg_names:
            header += ',' + arg_names
        header += ')\n'
        post['content'] += header
        post['content'] += '      implicit none\n'
        post['content'] += argument_code
        post['content'] += self_line
        call_args = f'({arg_names})' if arg_names else '()'
        if category == 'subroutine':
            post['content'] += (
                f'      call self%recursiveSelf%{m}{call_args}\n'
            )
        else:
            post['content'] += (
                f'      {fn}=self%recursiveSelf%{m}{call_args}\n'
            )
        post['content'] += f'   end {category} {fn}\n\n'

    # ---- framework overrides ----
    # allowedParameters: no-op (the shim owns no child objects; the real object
    # already contributed its parameters when it was built).
    post['content'] += (
        f'   recursive subroutine {shim_type}AllowedParameters(self,allowedParameters,sourceName,objectsOnly)\n'
        '      use ISO_Varying_String, only : varying_string\n'
        '      implicit none\n'
        f'      class    ({shim_type}   ), intent(inout)                            :: self\n'
        '      type     (varying_string), dimension(:), allocatable, intent(inout) :: allowedParameters\n'
        '      character(len=*         )                           , intent(in   ) :: sourceName\n'
        '      logical                                             , intent(in   ) :: objectsOnly\n'
        '      !$GLC attributes unused :: self, allowedParameters, sourceName, objectsOnly\n'
        '      return\n'
        f'   end subroutine {shim_type}AllowedParameters\n\n'
    )
    # autoHook: no-op. The real object already hooked itself when it was built;
    # the objectBuilder template calls autoHook on the shim too, so forwarding
    # or a default hook would double-register on calculationResetEvent.
    post['content'] += (
        f'   subroutine {shim_type}AutoHook(self)\n'
        '      implicit none\n'
        f'      class({shim_type}), intent(inout) :: self\n'
        '      !$GLC attributes unused :: self\n'
        '      return\n'
        f'   end subroutine {shim_type}AutoHook\n\n'
    )
    # isRecursiveShim: identifies this object as a shim so the objectBuilder
    # template skips caching it in the parameter node in place of the real
    # object (issue #695).
    post['content'] += (
        f'   logical function {shim_type}IsRecursiveShim(self)\n'
        '      implicit none\n'
        f'      class({shim_type}), intent(in   ) :: self\n'
        '      !$GLC attributes unused :: self\n'
        f'      {shim_type}IsRecursiveShim=.true.\n'
        f'   end function {shim_type}IsRecursiveShim\n\n'
    )
    # objectType / descriptor / hashedDescriptor / stateStore / stateRestore:
    # forward to the real object.
    post['content'] += (
        f'   function {shim_type}ObjectType(self,short)\n'
        '      use ISO_Varying_String, only : varying_string\n'
        '      implicit none\n'
        f'      type (varying_string)                            :: {shim_type}ObjectType\n'
        f'      class({shim_type}   ), intent(inout)             :: self\n'
        '      logical               , intent(in   ), optional :: short\n'
        f'      {shim_type}ObjectType=self%recursiveSelf%objectType(short)\n'
        f'   end function {shim_type}ObjectType\n\n'
    )
    post['content'] += (
        f'   subroutine {shim_type}Descriptor(self,descriptor,includeClass,includeFileModificationTimes)\n'
        '      use Input_Parameters, only : inputParameters\n'
        '      implicit none\n'
        f'      class({shim_type}    ), intent(inout)           :: self\n'
        '      type (inputParameters), intent(inout)           :: descriptor\n'
        '      logical               , intent(in   ), optional :: includeClass, includeFileModificationTimes\n'
        '      call self%recursiveSelf%descriptor(descriptor,includeClass,includeFileModificationTimes)\n'
        f'   end subroutine {shim_type}Descriptor\n\n'
    )
    post['content'] += (
        f'   function {shim_type}HashedDescriptor(self,includeSourceDigest,includeFileModificationTimes)\n'
        '      use ISO_Varying_String, only : varying_string\n'
        '      implicit none\n'
        f'      type (varying_string)                            :: {shim_type}HashedDescriptor\n'
        f'      class({shim_type}   ), intent(inout)             :: self\n'
        '      logical               , intent(in   ), optional :: includeSourceDigest, includeFileModificationTimes\n'
        f'      {shim_type}HashedDescriptor=self%recursiveSelf%hashedDescriptor(includeSourceDigest,includeFileModificationTimes)\n'
        f'   end function {shim_type}HashedDescriptor\n\n'
    )
    for store in ('stateStore', 'stateRestore'):
        post['content'] += (
            f'   subroutine {shim_type}{_ucfirst(store)}(self,stateFile,gslStateFile,stateOperationID)\n'
            '      use, intrinsic :: ISO_C_Binding, only : c_ptr, c_size_t\n'
            '      implicit none\n'
            f'      class  ({shim_type}), intent(inout) :: self\n'
            '      integer             , intent(in   ) :: stateFile\n'
            '      type   (c_ptr      ), intent(in   ) :: gslStateFile\n'
            '      integer(c_size_t   ), intent(in   ) :: stateOperationID\n'
            f'      call self%recursiveSelf%{store}(stateFile,gslStateFile,stateOperationID)\n'
            f'   end subroutine {shim_type}{_ucfirst(store)}\n\n'
        )
    # deepCopy: the shim is freshly mold-allocated by its owner (base fields at
    # default, referenceCount=0 -- matching what the generated assignment's
    # intent(out) reset produces for any copied object), so we only wire the
    # weak recursiveSelf back-pointer to the copy of the real object. If the
    # real object has already been copied in this operation its copiedSelf is
    # set; otherwise we defer resolution to deepCopyFinalize (see D4 note in
    # _build_deep_copy_methods: copiedSelf is deliberately NOT nulled during
    # finalize, so it survives for this resolution).
    post['content'] += (
        f'   recursive subroutine {shim_type}DeepCopy(self,destination)\n'
        '      use, intrinsic :: ISO_C_Binding, only : c_size_t\n'
        '      implicit none\n'
        f'      class({shim_type}      ), intent(inout), target :: self\n'
        f'      class({directive_name}Class), intent(inout)     :: destination\n'
        '      select type (destination)\n'
        f'      type is ({shim_type})\n'
        '         if (associated(self%recursiveSelf)) then\n'
        '            if (associated(self%recursiveSelf%copiedSelf)) then\n'
        '               destination%recursiveSelf  => self%recursiveSelf%copiedSelf\n'
        '               destination%parentDeferred =  .false.\n'
        '            else\n'
        '               destination%recursiveSelf  => self%recursiveSelf\n'
        '               destination%parentDeferred =  .true.\n'
        '            end if\n'
        '         end if\n'
        # Match the generated deepCopy: a freshly mold-allocated copy has
        # referenceCount=0, so reset it to 1 (a live reference held by the
        # owner) and clear the state-store id -- otherwise the owner's
        # objectDestructor decrements it to -1 and aborts.
        '         call destination%referenceCountReset()\n'
        '         destination%stateOperationID=0_c_size_t\n'
        '      end select\n'
        f'   end subroutine {shim_type}DeepCopy\n\n'
    )
    post['content'] += (
        f'   recursive subroutine {shim_type}DeepCopyReset(self)\n'
        '      implicit none\n'
        f'      class({shim_type}), intent(inout) :: self\n'
        '      self%copiedSelf => null()\n'
        f'   end subroutine {shim_type}DeepCopyReset\n\n'
    )
    # deepCopyFinalize: resolve a deferred parent. recursiveSelf currently
    # points at the ORIGINAL real object; its copiedSelf now holds the copy.
    post['content'] += (
        f'   recursive subroutine {shim_type}DeepCopyFinalize(self)\n'
        '      use Error, only : Error_Report\n'
        '      implicit none\n'
        f'      class({shim_type}), intent(inout) :: self\n'
        '      if (self%parentDeferred) then\n'
        '         if (associated(self%recursiveSelf).and.associated(self%recursiveSelf%copiedSelf)) then\n'
        '            self%recursiveSelf  => self%recursiveSelf%copiedSelf\n'
        '            self%parentDeferred =  .false.\n'
        '         else\n'
        f"            call Error_Report('recursive shim parent was not copied'//{loc_expr})\n"
        '         end if\n'
        '      end if\n'
        f'   end subroutine {shim_type}DeepCopyFinalize\n\n'
    )


def _normalise_modules(modules):
    """Normalise the `modules` entry on a method dict into a list of
    `{name, only}` shapes.
    """
    if isinstance(modules, list):
        return [m for m in modules if isinstance(m, dict) and 'name' in m]
    if isinstance(modules, dict):
        if 'name' in modules and all(not isinstance(v, dict)
                                     for v in modules.values()):
            return [modules]
        return [
            dict(attrs, name=name) if isinstance(attrs, dict)
            else {'name': name}
            for name, attrs in modules.items()
        ]
    if isinstance(modules, str):
        return [{'name': tok} for tok in modules.split() if tok]
    return []


def process_function_class(tree, options):
    """Process `functionClass` directives in the tree.

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
            directive = node.setdefault('directive', {})
            directive['processed'] = True

        if node.get('type') != 'functionClass':
            continue
        directive = node.setdefault('directive', {})
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

        # Extract `method` sub-directive into a {name: dict} map.
        # xml_to_dict gives a list for multiple methods, a single dict
        # otherwise.
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
            directive, non_abstract_classes, classes, methods, tree,
            state_storables)
        _build_object_type_method(directive, non_abstract_classes, methods)
        _build_allowed_parameters_method(
            directive, classes_ordered, methods)
        _build_assignment_method(
            directive, non_abstract_classes, classes, methods,
            state_storables)
        _build_deep_copy_methods(
            directive, non_abstract_classes, classes,
            node.get('line', 0), methods,
            state_storables, deep_copy_actions)
        _build_state_store_methods(
            directive, non_abstract_classes, classes, methods,
            state_storables)

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

        # Emit the generated shim type for families that opt in to bounded
        # construction recursion (issue #695).
        if any(c.get('recursive') == 'yes' for c in classes_ordered):
            _generate_recursive_shim(
                directive, methods, non_abstract_classes, pre, post, node)

        _insert_and_write_output(
            node, code_content, pre, post, classes, directive)


register_process('functionClass', process_function_class)
