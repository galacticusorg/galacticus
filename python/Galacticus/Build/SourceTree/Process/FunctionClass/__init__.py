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
    serialize, set_visibility,
)
from Galacticus.Build.SourceTree.Process                     import (
    register_process, process_tree,
)
from Galacticus.Build.SourceTree.Parse.Declarations          import (
    parse_declaration, build_declarations, declaration_exists, get_declaration,
)
from Galacticus.Build.SourceTree.Parse.ModuleUses            import add_uses
from Galacticus.Build.SourceTree.Parse.Visibilities          import update_visibilities
from Galacticus.Build.SourceTree.Process.FunctionClass.Utils import (
    class_dependencies,
)
from Galacticus.Build.SourceTree.Process.SourceIntrospection import location


# ---------------------------------------------------------------------------
# Small formatting helpers (local equivalents of tiny Perl utilities)
# ---------------------------------------------------------------------------

def _format_variable_definitions(declarations):
    """Render a list of declaration dicts into Fortran declaration lines.

    Lightweight stand-in for Perl `Fortran::Utils::Format_Variable_Definitions`,
    whose output is a Text::Table-aligned column layout.  We only need
    functional equivalence here — a simple
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
    """Return the C-binding declaration for the per-class MD5 symbol.

    Mirrors SourceDigest::Binding() at SourceDigest.pm:60-65.
    """
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
    """Minimal LaTeX special-character escape.  Matches the subset of
    Perl `LaTeX::Encode::latex_encode` that FunctionClass.pm relies on.
    """
    return ''.join(_LATEX_ESCAPES.get(c, c) for c in str(s))


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


def _build_assignment_method(directive, non_abstract_classes, classes,
                             methods, state_storables):
    """Populate `methods['assignment(=)']` with the defined-assignment
    operator overload.

    Mirrors buildAssignmentMethod() at FunctionClass.pm:2537-2683.
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


def _build_deep_copy_methods(*args, **kwargs):
    _not_implemented('buildDeepCopyMethods')


def _build_state_store_methods(*args, **kwargs):
    _not_implemented('buildStateStoreMethods')


def _generate_type_definition(directive, methods, pre, node):
    """Emit the `type, extends(...) :: <name>Class` block with methods
    table, generic interfaces, destructor binding, and the module-level
    data declarations.

    Mirrors generateTypeDefinition() at FunctionClass.pm:216-417.
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
    metadata.  Mirrors the nested-regex loop at lines 267-308 of
    generateTypeDefinition; we lean on our existing `parse_declaration`
    rather than the Perl-style INTRINSIC_DECLARATIONS indices.
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

    Mirrors generateConstructor() at FunctionClass.pm:419-593.
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
        pre['content'] += (
            f'   type(inputParameter), pointer :: '
            f'{directive_name}RecursiveBuildNode => null()\n'
            f'   class({directive_name}Class), pointer :: '
            f'{directive_name}RecursiveBuildObject => null()\n'
            f'   !$omp threadprivate({directive_name}RecursiveBuildNode,'
            f'{directive_name}RecursiveBuildObject)\n\n'
        )
        add_uses(parent, {
            'moduleUse': {
                'Input_Parameters': {'intrinsic': False, 'all': True},
            },
            'moduleOrder': ['Input_Parameters'],
        })

    # Append directive name to the per-file `.p` parameter-names file.
    if tree.get('type') == 'file':
        build_path = os.environ.get('BUILDPATH')
        tree_name = tree.get('name', '')
        m = re.match(r'(.+)\.F90$', tree_name)
        if build_path and m:
            with open(os.path.join(build_path, m.group(1) + '.p'), 'a') as fh:
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
        'inputParameters\n'
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
        if match is not None and match.get('recursive') == 'yes':
            post['content'] += (
                f'        {directive_name}RecursiveBuildNode   => parameterNode\n'
                f'        {directive_name}RecursiveBuildObject => self\n'
            )
        post['content'] += (
            '        select type (self)\n'
            f'          type is ({target_name})\n'
            f'            self={target_name}(subParameters)\n'
            '         end select\n'
        )
        if match is not None and match.get('recursive') == 'yes':
            post['content'] += (
                f'        {directive_name}RecursiveBuildNode   => null()\n'
                f'        {directive_name}RecursiveBuildObject => null()\n'
            )
        post['content'] += '         call parameterNode%objectSet(self)\n'
        post['content'] += '      else\n'

    if allow_recursion:
        post['content'] += (
            '        parameterNode => parameters%node'
            '(char(parameterName_),requireValue=.true.)\n'
            f'        if (associated(parameterNode,'
            f'{directive_name}RecursiveBuildNode)) then\n'
        )
        for c in non_abstract_classes:
            if c.get('recursive') != 'yes':
                continue
            class_name = c['name']
            post['content'] += (
                f'           select type ({directive_name}RecursiveBuildObject)\n'
                f'              type is ({class_name})\n'
                f'              allocate({class_name} :: self)\n'
                '              select type (self)\n'
                f'              type is ({class_name})\n'
                '                 self%isRecursive=.true.\n'
                f'                 self%recursiveSelf => '
                f'{directive_name}RecursiveBuildObject\n'
                '              end select\n'
                '           end select\n'
            )
        post['content'] += (
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
        '      select case (char(instanceName))\n'
    )
    for c in non_abstract_classes:
        name = _short_name(c['name'], directive_name)
        post['content'] += f"     case ('{name}')\n"
        post['content'] += f"        allocate({c['name']} :: self)\n"
        if c.get('recursive') == 'yes':
            post['content'] += (
                f"        {directive_name}RecursiveBuildNode   => parameterNode\n"
                f"        {directive_name}RecursiveBuildObject => self\n"
            )
        post['content'] += (
            '        select type (self)\n'
            f"          type is ({c['name']})\n"
            f"            self={c['name']}(subParameters)\n"
            '         end select\n'
        )
        if c.get('recursive') == 'yes':
            post['content'] += (
                f"        {directive_name}RecursiveBuildNode   => null()\n"
                f"        {directive_name}RecursiveBuildObject => null()\n"
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

    Mirrors generateClassSubmodules() at FunctionClass.pm:595-1024.  Walks
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

        # Submodule output path: swap the source file's basename `.F90`
        # for `.p.F90` and place it under BUILDPATH.
        class_file = class_record.get('file') or ''
        base = os.path.basename(class_file)
        base = re.sub(r'\.F90$', '.p.F90', base)
        fn = os.path.join(build_path, base) if build_path else base

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
                class_node = class_node.get('firstChild')
                contained = True
                if class_node is None:
                    break
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
                entry = uses[module_name]
                if entry.get('all'):
                    continue
                only = entry.get('only', {}) or {}
                kept = {
                    sym: flag for sym, flag in only.items()
                    if sym.lower() in module_symbols
                }
                entry['only'] = kept
                if not kept:
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

    Mirrors the branch at FunctionClass.pm:657-813.
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
    variable → local declarations).  Mirrors lines 706-759 of
    generateClassSubmodules.
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
    partition.  Mirrors the large if/elif ladder at lines 815-985 of
    generateClassSubmodules.
    """
    ntype = class_node.get('type')
    directive_name = directive['name']

    if ntype == 'scoping':
        d = class_node.get('directive') or {}
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
    (→ module_symbols).  Mirrors lines 843-880 of generateClassSubmodules.
    """
    post_contains = False
    type_node = class_node.get('firstChild')
    while type_node is not None:
        if type_node.get('type') == 'contains':
            post_contains = True
            type_node = type_node.get('firstChild')
            if type_node is None:
                break
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

    Mirrors lines 889-980 of generateClassSubmodules.
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

    Mirrors generateMethodFunctions() at FunctionClass.pm:1026-1167.
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


def _normalise_modules(modules):
    """Normalise the `modules` entry on a method dict into a list of
    `{name, only}` shapes.  Mirrors the `reftype($method->{'modules'})`
    branch at generateMethodFunctions.pm:1109-1124.
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


def _generate_documentation(directive, classes, non_abstract_classes):
    """Write a LaTeX `doc/physics/<descriptive>.tex` file describing every
    class in the functionClass hierarchy.

    Mirrors generateDocumentation() at FunctionClass.pm:1169-1399.  The
    Perl version uses `system("mkdir -p doc/physics")` + a hard-coded
    write; we do the same, relative to the current working directory,
    which matches the Perl build-system assumption.
    """
    directive_name   = directive['name']
    descriptive_name = directive.get('descriptiveName', directive_name)

    doc = (
        f'\\section{{{descriptive_name}}}\\label{{phys:{directive_name}}}'
        f'\\hyperdef{{physics}}{{{directive_name}}}{{}}\n\n'
    )
    if 'description' in directive:
        doc += directive['description'] + '\n\n'
    if 'default' in directive:
        default_target = directive_name + _ucfirst(directive['default'])
        doc += (
            f'Default implementation: \\refPhysics{{{default_target}}}\n\n'
        )
    else:
        doc += 'No default implementation\n\n'

    for class_name in sorted(classes.keys(), key=str.lower):
        class_record = classes[class_name]
        suffix = class_record['name']
        if suffix.startswith(directive_name):
            suffix = suffix[len(directive_name):]
        if not re.match(r'^[A-Z]{2}', suffix):
            suffix = _lcfirst(suffix)
        doc += (
            f'\\subsection{{\\mono{{{suffix}}}}}'
            f'\\label{{phys:{class_record["name"]}}}'
            f'\\hyperdef{{physics}}{{{class_record["name"]}}}{{}}\n\n'
        )
        doc += (class_record.get('description') or '') + '\n\n'
        if ('default' in directive
                and directive_name + _ucfirst(directive['default'])
                    == class_record['name']):
            doc += '\\noindent \\textbf{(Default)}\n\n'
        doc += (
            f'\\noindent \\emph{{Implemented by}} '
            f'\\refClass{{{class_record["name"]}}}\n'
        )

        # Walk the tree for this class to extract constructor parameter /
        # objectBuilder information.  The level of detail here matches
        # Perl's: we walk interface/moduleProcedure pairs, then hunt for
        # matching function nodes to find inputParameter / objectBuilder
        # directives inside.
        tree = class_record.get('tree')
        if tree is None:
            continue

        constructors = _collect_constructor_names(tree, class_name)
        parameters, objects = _collect_doc_parameters_and_objects(
            tree, constructors, class_name, class_record, classes, directive,
        )
        if parameters:
            doc += (
                '\n\n\\noindent\\emph{Parameters}\n'
                '\\begin{description}\n'
                + '\n'.join(parameters)
                + '\n\\end{description}\n'
            )
        if objects:
            sorted_objects = sorted(objects)
            doc += (
                '\n\\noindent\\emph{Classes used}\n\n'
                '\\begin{tabular}{ll}\n'
            )
            for i in range(0, len(sorted_objects), 2):
                doc += f'\\refPhysics{{{sorted_objects[i]}}}'
                if i + 1 < len(sorted_objects):
                    doc += f' & \\refPhysics{{{sorted_objects[i + 1]}}}'
                doc += '\\\\\n'
            doc += '\\end{tabular}\n\n'

    out_dir = os.path.join('doc', 'physics')
    os.makedirs(out_dir, exist_ok=True)
    file_base = re.sub(r'\s+', '_', descriptive_name.lower())
    out_path = os.path.join(out_dir, f'{file_base}.tex')
    with open(out_path, 'w') as fh:
        fh.write(doc)


def _collect_constructor_names(tree, class_name):
    """Find the `interface <class>` block in the class tree and return the
    list of module-procedure names it declares.  Mirrors the first walk
    inside generateDocumentation() at lines 1196-1207.
    """
    node = tree.get('firstChild')
    while node is not None and not (
            node.get('type') == 'interface'
            and node.get('name') == class_name):
        node = node.get('sibling')
    if node is None:
        return []
    constructors = []
    child = node.get('firstChild')
    while child is not None:
        if child.get('type') == 'moduleProcedure':
            constructors.extend(child.get('names') or [])
        child = child.get('sibling')
    return constructors


def _collect_doc_parameters_and_objects(
        tree, constructors, class_name, class_record, classes, directive):
    """For each constructor function in `tree`, walk its body collecting
    inputParameter directives (→ `parameters` list, as LaTeX-encoded
    `\\item` strings) and objectBuilder directives (→ `objects` list of
    class names).  Mirrors lines 1209-1376 of generateDocumentation().
    """
    parameters = []
    objects    = []
    directive_name = directive['name']

    node = tree.get('firstChild')
    while node is not None:
        if (node.get('type') == 'function'
                and node.get('name') in constructors):
            for cnode in walk_tree(node):
                if cnode is node:
                    continue
                if cnode.get('type') == 'inputParameter':
                    entry = _format_input_parameter_doc(
                        cnode, node, class_record, classes,
                        directive_name,
                    )
                    if entry is not None:
                        parameters.append(entry)
                elif cnode.get('type') == 'objectBuilder':
                    objects.append(
                        (cnode.get('directive') or {}).get('class'))
        # Descend into `contains` blocks — the Perl idiom
        # `$node = $node->{'type'} eq "contains" ? $node->{'firstChild'}
        #                                         : $node->{'sibling'};`
        if node.get('type') == 'contains':
            node = node.get('firstChild')
        else:
            node = node.get('sibling')
    return parameters, [o for o in objects if o]


def _format_input_parameter_doc(cnode, function_node, class_record,
                                classes, directive_name):
    """Render one inputParameter directive into a LaTeX `\\item[...]` line
    for the class's documentation entry.
    """
    cdir = cnode.get('directive') or {}
    variable_name = cdir.get('variable') or cdir.get('name') or ''
    variable_name = re.sub(r'\(.+\)$', '', variable_name)

    declaration = None
    m = re.match(
        r'([a-zA-Z0-9_]+)(\s*\(\s*[a-zA-Z0-9_:,]\s*\)\s*)??'
        r'%([a-zA-Z0-9_]+)',
        variable_name,
    )
    if m:
        object_name = m.group(1)
        member      = m.group(3)
        if object_name in ('self', function_node.get('name')):
            # Walk up the extends chain looking for the member declaration.
            cursor = class_record
            while cursor is not None:
                type_node = _find_type_node_named(cursor.get('tree'),
                                                  cursor.get('name'))
                if type_node is not None and declaration_exists(
                        type_node, member):
                    declaration = get_declaration(type_node, member)
                if declaration is not None:
                    break
                extends = cursor.get('extends')
                if extends is None or extends == directive_name:
                    break
                cursor = classes.get(extends)
        else:
            # Non-self object — handle stateful{Double,Integer,Logical}.
            if declaration_exists(function_node, object_name):
                d_tmp = get_declaration(function_node, object_name)
                if (d_tmp.get('intrinsic') == 'type'
                        and d_tmp.get('type') in (
                            'statefulDouble', 'statefulInteger',
                            'statefulLogical')):
                    declaration = dict(d_tmp)
                    if d_tmp['type'] == 'statefulDouble':
                        declaration['intrinsic'] = 'double precision'
                    elif d_tmp['type'] == 'statefulInteger':
                        declaration['intrinsic'] = 'integer'
                    else:
                        declaration['intrinsic'] = 'logical'
                    declaration['type'] = None
    else:
        if declaration_exists(function_node, variable_name):
            declaration = get_declaration(function_node, variable_name)

    if declaration is None:
        if 'type' in cdir and 'cardinality' in cdir:
            declaration = {
                'parameterType':        cdir['type'],
                'parameterCardinality': cdir['cardinality'],
            }
        else:
            # Bail out with a clear message — matches Perl's `die('abort')`.
            raise RuntimeError(
                f"generate_documentation: unable to find parameter variable "
                f"declaration for \"{cdir.get('name')}\" in class "
                f"\"{class_record.get('name')}\""
            )

    # Determine type + cardinality.
    if 'parameterType' in declaration:
        type_name = declaration['parameterType']
    else:
        intrinsic = declaration.get('intrinsic')
        if intrinsic == 'double precision':
            type_name = 'real'
        elif intrinsic == 'integer':
            type_name = 'integer'
        elif intrinsic == 'logical':
            type_name = 'boolean'
        elif intrinsic == 'character':
            type_name = 'string'
        elif (intrinsic == 'type'
              and declaration.get('type') == 'varying_string'):
            type_name = 'string'
        else:
            raise RuntimeError(
                'generate_documentation: unable to determine parameter type')

    if 'parameterCardinality' in declaration:
        cardinality = declaration['parameterCardinality']
    else:
        attrs = declaration.get('attributes') or []
        dims = [a for a in attrs if a.startswith('dimension')]
        if dims:
            d = re.match(r'^dimension\s*\(\s*(.+?)\s*\)', dims[0])
            if not d:
                raise RuntimeError(
                    'generate_documentation: unable to parse dimension attribute')
            shape_parts = [s.strip() for s in d.group(1).split(',')]
            cardinality_max = 1
            for part in shape_parts:
                if part.isdigit():
                    cardinality_max *= int(part)
                else:
                    cardinality_max = '*'
                    break
        else:
            cardinality_max = 1
        if cardinality_max == '*':
            cardinality = ('0..*' if 'defaultValue' in cdir else '1..*')
        else:
            if 'defaultValue' in cdir:
                cardinality = f'0,{cardinality_max}'
            else:
                cardinality = cardinality_max

    description = (
        f'\\item[\\mono{{[{_latex_encode(cdir.get("name", ""))}]}}] '
        f'({type_name}; {cardinality}) '
    )
    if 'defaultValue' in cdir:
        value = _latex_encode(cdir['defaultValue'])
        if value == '.true.':
            value = 'true'
        elif value == '.false.':
            value = 'false'
        if type_name == 'real':
            value = re.sub(r'(\d)d([\+\-0-9])', r'\1e\2', value)
        if type_name == 'integer':
            value = value.replace(r'\_c\_size\_t', '')
        if type_name == 'string':
            m2 = re.match(r"^var\\_str\(['\"](.*)['\"]\)$", value)
            if m2:
                value = m2.group(1)
        default_source = (
            '; ' + cdir['defaultSource']
            if 'defaultSource' in cdir else ''
        )
        description += (
            f' \\{{\\mono{{{value}}}{default_source}\\}} '
        )
    description += cdir.get('description', '')
    return description


def _find_type_node_named(tree, class_name):
    """Walk `tree` looking for a `type` node with the given name.  Used by
    _format_input_parameter_doc to crawl up the extends chain."""
    if tree is None:
        return None
    node = tree.get('firstChild')
    while node is not None:
        if (node.get('type') == 'type'
                and node.get('name') == class_name):
            return node
        node = node.get('sibling')
    return None


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
            directive, non_abstract_classes, classes, methods,
            state_storables)
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
