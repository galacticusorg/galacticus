# Processes `classDocumentation` directives (only during a documentation
# build — gated on $GALACTICUS_BUILD_DOCS=yes): walks every `type` node in
# the current tree (plus, for `objects.nodes.F90`, the associated
# `.Inc` include files), builds a per-class record of methods + their
# bound functions + argument declarations, and dumps the result to
# `$BUILDPATH/<basename>.classes.xml`.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/ClassDocumentation.pm

import copy
import os
import re
import sys
import xml.etree.ElementTree as ET


from List.ExtraUtils                     import as_array
from XML.Utils                           import xml_to_dict
from Fortran.Utils                       import UNIT_OPENERS
from Galacticus.Build.Directives         import extract_directives
from Galacticus.Build.SourceTree         import walk_tree, parse_file
from Galacticus.Build.SourceTree.Parse.Declarations import parse_declaration
from Galacticus.Build.SourceTree.Process import (
    register_postprocess, is_inner_process_tree_call,
)


# Classes already emitted across prior hook invocations — matches Perl's
# module-level `@outputPrevious`.
_OUTPUT_PREVIOUS = []


def _build_tree_list(tree):
    """For the `objects.nodes.F90` tree, also parse the include files that
    contribute type-bound procedure definitions via `<component>` directives.
    Mirrors the include-file hack at ClassDocumentation.pm:33-49.
    """
    trees = [tree]
    if tree.get('name') != 'objects.nodes.F90':
        return trees
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        return trees
    loc_path = os.path.join(build_path, 'directiveLocations.xml')
    if not os.path.exists(loc_path):
        return trees
    directive_locations = xml_to_dict(ET.parse(loc_path).getroot())
    include_files = ['objects.nodes.components.Inc']
    component_block = directive_locations.get('component') or {}
    for component_file in as_array(component_block.get('file')):
        for d in extract_directives(component_file, 'component'):
            functions = d.get('functions')
            if isinstance(functions, str) and functions.endswith('.inc'):
                include_files.append(functions.replace('.inc', '.Inc'))
    for include_file in include_files:
        include_path = os.path.join(build_path, include_file)
        if os.path.exists(include_path):
            trees.append(parse_file(include_path))
    return trees


def _extract_function_arguments(function_node):
    """Given a function or subroutine opener, return its argument-name list
    (with `self` excluded).  Uses UNIT_OPENERS' own regex — not a parallel
    split — so the result matches what the Perl port does.
    """
    ntype = function_node.get('type')
    spec  = UNIT_OPENERS.get(ntype) or {}
    regex = spec.get('regex')
    arg_idx = spec.get('arguments')
    if regex is None or arg_idx is None:
        return []
    m = regex.match(function_node.get('opener') or '')
    if not m:
        return []
    groups = m.groups()
    if arg_idx >= len(groups) or groups[arg_idx] is None:
        return []
    args_text = groups[arg_idx]
    return [a.strip() for a in re.split(r'\s*,\s*', args_text)
            if a.strip() and a.strip() != 'self']


# ---------------------------------------------------------------------------
# Type-node harvesting
# ---------------------------------------------------------------------------

def _iter_contained_children(type_node):
    """Yield children of `type_node` that live after the `contains` marker.

    Mirrors the inner `if ($child->{'type'} eq "contains")` loop at
    ClassDocumentation.pm:61-78.  Our `contains` is a self-closing sibling
    marker so the iteration is simpler than Perl's (which nests under the
    contains node).
    """
    child = type_node.get('firstChild')
    contained = False
    while child is not None:
        if child.get('type') == 'contains':
            contained = True
            child = child.get('sibling')
            continue
        if contained:
            yield child
        child = child.get('sibling')


def _populate_class_descriptions(type_node, class_record):
    """Fill `descriptions` (from `<methods>` directives) and `functions`
    (declaration records) into `class_record` from the contains-area children
    of `type_node`.  Matches the assignments at ClassDocumentation.pm:69-74.
    """
    for child in _iter_contained_children(type_node):
        ntype = child.get('type')
        if ntype == 'methods':
            directive = child.get('directive') or {}
            methods = list(as_array(directive.get('method')))
            # Normalise the method-name key.  `<method name="…">` (used by
            # most parent classes) and `<method method="…">` (used by most
            # child classes) both need to surface as `method:` for the
            # downstream consumer (`scripts/doc/extractData.py`) to find
            # the method when it walks the inheritance chain to remove
            # already-described "missing" methods.  Perl's XML::Simple
            # would collapse multiple `<method name="…">` children into a
            # name-keyed hash — our `xml_to_dict` keeps them as a list of
            # `{'name': ..., ...}` dicts, so we explicitly promote `name`
            # to `method` here to give both forms the same shape.
            for m in methods:
                if isinstance(m, dict) and 'method' not in m and 'name' in m:
                    m['method'] = m['name']
            class_record.setdefault('descriptions', []).extend(methods)
        elif ntype == 'declaration':
            class_record.setdefault('functions', []).append(
                child.get('declarations') or []
            )


def _normalise_directive_method(method):
    """Reshape a `<method>` directive dict so its `type` and `arguments`
    fields match the structure the doc consumer expects.

    The `<functionClass>` directive's children carry their type/argument
    information as raw Fortran-source strings (e.g.
    `type(keplerOrbit)`, `type(treeNode), intent(inout) :: node, host`),
    while `extractData.py`'s `declaration_builder` requires a parsed
    declaration dict (`{intrinsic, type, attributes, variables}`) — the
    same shape `_process_function` produces for `<type>`-derived
    methods.  Convert the strings here so both code paths feed the
    consumer the same dict shape.

    * `method['type']` (string): wrap in `:: __dummy__`, parse as a
      declaration, return the dict.  When the directive has no `<type>`
      (i.e. the method is subroutine-like), set the literal string
      `'subroutine'` — the consumer recognises that form.
    * `method['argument']` (list of strings, each a declaration line):
      parse each, then re-package as
      `arguments=[{'argument': [parsed, ...]}]` plus a parallel
      `argumentList=['name1, name2, ...']` of just the variable names —
      again matching the shape `_process_function` produces.
    """
    if not isinstance(method, dict):
        return

    raw_type = method.get('type')
    if isinstance(raw_type, str) and raw_type.strip():
        stripped = raw_type.strip()
        if stripped.lower() == 'void':
            # FunctionClass auto-generates void-return method stubs
            # (autoHook, descriptor, deepCopy, …) with `type: 'void'` —
            # the doc consumer's `declaration_builder` only recognises
            # the literal string `'subroutine'` for the void-return
            # case, so translate.
            method['type'] = 'subroutine'
        else:
            # `parse_declaration` requires a `::` separator and at least one
            # variable name.  Append a synthetic name we discard afterwards.
            parsed = parse_declaration(stripped + ' :: __dummy__')
            if parsed is not None:
                parsed['variables']     = []
                parsed['variableNames'] = []
                method['type'] = parsed
            # else: leave the raw string in place — the consumer will raise,
            # but at least we've not actively lost information.
    elif raw_type is None:
        method['type'] = 'subroutine'

    raw_args = method.get('argument')
    if raw_args is not None:
        arg_strings = list(as_array(raw_args)) if not isinstance(raw_args, list) \
            else list(raw_args)
        parsed_args = []
        names_collected = []
        for arg in arg_strings:
            if not isinstance(arg, str):
                continue
            parsed = parse_declaration(arg.strip())
            if parsed is None:
                continue
            parsed_args.append(parsed)
            names_collected.extend(parsed.get('variableNames') or [])
        if parsed_args:
            method.setdefault('arguments',
                              []).append({'argument': parsed_args})
            method.setdefault('argumentList',
                              []).append(', '.join(names_collected))


def _populate_class_from_function_class_directive(node, classes):
    """Synthesise a class record for the abstract `<name>Class` type that
    FunctionClass auto-generates from a `<functionClass>` directive.

    The `<functionClass>` directive at module scope carries the parent
    class's methods directly as `<method name="…">` children.  Promote the
    `name` attribute on each method to a `method` key (so the rest of the
    pipeline — including the doc consumer at
    `scripts/doc/extractData.py`'s inheritance walk — can find them by
    method name), parse each method's `<type>` / `<argument>` strings
    into the declaration-dict shape the consumer expects (see
    `_normalise_directive_method`), assemble them into a `<descriptions>`
    list, and stash the synthesised class under `<base_name>Class` in
    the classes dict.
    """
    directive = node.get('directive') or {}
    base_name = directive.get('name')
    if not base_name:
        return
    class_name = base_name + 'Class'
    class_record = classes.setdefault(class_name, {'name': class_name})

    # FunctionClass auto-generates `type, abstract, extends(functionClass) ::
    # <name>Class`, so the synthesised record needs to advertise that
    # parentage too.  The doc consumer walks the `extends` chain to set
    # `isFunctionClass = True` on every descendant, which gates the
    # filtering of auto-generated methods (autoHook, descriptor, deepCopy,
    # …) from the "missing method descriptions" warning.  Without this,
    # every concrete child class shows the full set of auto-generated
    # method names in the docs build's warning output.
    class_record.setdefault('extends', 'functionClass')

    # Deep-copy each method dict before mutating it.  These dicts live in
    # the original `<functionClass>` directive and are consumed later by
    # FunctionClass — which expects `method['type']` to remain the raw
    # Fortran-source string (`type(keplerOrbit)`, `double precision`, …)
    # so it can do its own startswith/string checks.  The doc-side
    # normalisation (`_normalise_directive_method`) replaces that string
    # with a parsed declaration dict, so without the deep-copy
    # FunctionClass crashes with
    # `AttributeError: 'dict' object has no attribute 'startswith'`.
    methods = [copy.deepcopy(m) if isinstance(m, dict) else m
               for m in as_array(directive.get('method'))]
    for m in methods:
        if isinstance(m, dict) and 'method' not in m and 'name' in m:
            m['method'] = m['name']
        _normalise_directive_method(m)
    class_record.setdefault('descriptions', []).extend(methods)

    # Abstract base classes have no procedure bindings to resolve and no
    # missing methods of their own.  Initialise the fields explicitly so
    # `_dict_to_xml_element` emits an empty `<missingMethods/>` /
    # `<genericUses/>` rather than silently omitting them — extractData.py
    # tolerates either, but matching the shape of `<type>`-derived class
    # records keeps the consumer's logic uniform.
    class_record.setdefault('missingMethods', '')
    class_record.setdefault('genericUses',    '')


def _resolve_method_bindings(class_record):
    """For each method description, resolve the bound Fortran function name(s)
    by scanning the class's collected function declarations.

    Mirrors the nested loop at ClassDocumentation.pm:85-124 — handles
    `deferred` attributes (recorded as the interface name), `=>` procedure
    bindings, and generic-binding indirection.
    """
    descriptions = class_record.get('descriptions') or []
    functions_lists = class_record.get('functions') or []

    for method in descriptions:
        target = (method.get('method') or '').lower()
        for functions in functions_lists:
            for function in as_array(functions):
                variables = function.get('variables') or []
                if not variables:
                    continue
                variable_name = re.sub(r'\s*=>.*', '', variables[0])
                if variable_name != target:
                    continue
                # Matching function — interface for deferred, else parse
                # the `=>` procedure list.
                if 'deferred' in (function.get('attributes') or []):
                    method['interface'] = function.get('type')
                    continue
                # Look for `=>boundFunction[,…]` ANYWHERE in the variable
                # text — `variables[0]` is something like
                # `parametersselect=>jiang2014parametersselect`, not a string
                # that starts with `=>`, so `re.match` (anchored) would
                # always miss.  The Perl original used `m//` which is
                # unanchored; mirror that with `re.search`.  When this
                # check fails the method's `boundFunctions` never gets
                # populated, and later `_process_function` can't attach a
                # type — leaving `extractData.py` to print
                # "missing function type" warnings for every type-bound
                # method in the class.
                m = re.search(r'=>\s*([a-zA-Z0-9_,]+)', variables[0])
                if not m:
                    continue
                for variable in variables:
                    fm = re.match(
                        r'\s*([a-zA-Z0-9_.()+\-*/=<>]+)=>', variable)
                    if not fm:
                        continue
                    function_name = variable[fm.end():]
                    if function.get('intrinsic') == 'generic':
                        function_name = _resolve_generic_binding(
                            function_name, functions_lists)
                    method.setdefault('boundFunctions', []).append(
                        function_name)


def _resolve_generic_binding(function_name, functions_lists):
    """Resolve a generic-binding target by scanning all siblings.  Follows
    `methodName => boundFunction` chains and falls back to a sibling's
    `deferred` interface type when present.  Mirrors ClassDocumentation.pm:
    101-114.
    """
    for functions in functions_lists:
        for sibling in as_array(functions):
            variables = sibling.get('variables') or []
            if not variables:
                continue
            first = variables[0]
            m = re.match(r'(\S+)\s*=>\s*(\S+)', first)
            if m:
                method_name = m.group(1)
                bound = m.group(2)
                if method_name.lower() == function_name.lower():
                    return bound
            elif (first.lower() == function_name.lower()
                  and sibling.get('type') is not None
                  and 'deferred' in (sibling.get('attributes') or [])):
                return sibling['type']
    return function_name


def _compute_missing_and_generic(class_record):
    """Populate the `missingMethods` and `genericUses` lists for a class.

    Mirrors ClassDocumentation.pm:125-152.
    """
    descriptions = class_record.get('descriptions') or []
    described_methods = {
        (d.get('method') or '').lower() for d in descriptions
    }
    functions_lists = class_record.get('functions') or []

    missing = []
    generic = []
    for functions in functions_lists:
        for function in as_array(functions):
            variables = function.get('variables') or []
            if not variables:
                continue
            function_name = re.sub(r'=>.*', '', variables[0])
            if (function.get('intrinsic') == 'final'
                    or '=>' not in variables[0]
                    or function_name.lower() in described_methods):
                continue
            used_in_generic = _is_used_in_generic(function_name, functions_lists)
            if used_in_generic:
                generic.append(function_name)
            else:
                missing.append(function_name)
    class_record['missingMethods'] = ' '.join(missing)
    class_record['genericUses']    = ' '.join(generic)


def _is_used_in_generic(function_name, functions_lists):
    target = function_name.lower()
    for functions in functions_lists:
        for function_generic in as_array(functions):
            if function_generic.get('intrinsic') != 'generic':
                continue
            for variable in function_generic.get('variables') or []:
                bound = re.sub(r'.*=>\s*', '', variable)
                if bound == target:
                    return True
    return False


# ---------------------------------------------------------------------------
# Function processing (argument type extraction)
# ---------------------------------------------------------------------------

def _process_function(function_node, classes):
    """Mirrors processFunction() at ClassDocumentation.pm:183-258.

    For every method description whose `interface` or `boundFunctions`
    names this function, attach the argument declarations and — for
    function-typed methods — the return-type declaration.
    """
    function_name = function_node.get('name')
    arguments = _extract_function_arguments(function_node)

    for class_record in classes.values():
        for method in class_record.get('descriptions') or []:
            in_interface = (method.get('interface')
                            and method['interface'].lower()
                            == (function_name or '').lower())
            in_bindings  = any(
                b.lower() == (function_name or '').lower()
                for b in method.get('boundFunctions') or []
            )
            if not (in_interface or in_bindings):
                continue

            method['type'] = function_node.get('type')
            method.setdefault('argumentList', []).append(
                ', '.join(arguments))

            function_return_name = None
            if method['type'] == 'function':
                m = re.search(
                    r'\sresult\s*\(\s*([a-zA-Z0-9_]+)\s*\)',
                    function_node.get('opener') or '',
                )
                function_return_name = m.group(1) if m else function_name

            argument_declarations = []
            function_content = function_node.get('firstChild')
            while function_content is not None:
                if function_content.get('type') == 'declaration':
                    for declaration in function_content.get('declarations') or []:
                        matched_names = [
                            a for a in arguments
                            if a in (declaration.get('variableNames') or [])
                        ]
                        if matched_names:
                            declaration_copy = copy.deepcopy(declaration)
                            if declaration_copy.get('type') is None:
                                declaration_copy['type'] = ''
                            declaration_copy['variableNames'] = list(matched_names)
                            declaration_copy['variables']     = list(matched_names)
                            argument_declarations.append(declaration_copy)
                        # Return-type detection.
                        if (function_return_name
                                and any(v.lower()
                                        == function_return_name.lower()
                                        for v in (declaration.get('variableNames')
                                                  or []))):
                            return_type = copy.deepcopy(declaration)
                            if return_type.get('type') is None:
                                return_type['type'] = ''
                            return_type['variableNames'] = [function_return_name]
                            return_type['variables']     = [function_return_name]
                            method['type'] = return_type
                function_content = function_content.get('sibling')

            method.setdefault('arguments', []).append(
                {'argument': argument_declarations})

            # If this is still a bare `function` string (no declaration found
            # for the return value), fall back to the opener's intrinsic/kind.
            if method['type'] == 'function':
                spec = UNIT_OPENERS.get('function') or {}
                m = spec.get('regex').match(function_node.get('opener') or '')
                if m:
                    intrinsic = (m.groups()[spec['intrinsic']]
                                 if spec.get('intrinsic') is not None else '')
                    kind_idx  = spec.get('kind')
                    kind      = (m.groups()[kind_idx]
                                 if kind_idx is not None
                                 and kind_idx < len(m.groups()) else None)
                    method['type'] = {
                        'intrinsic': intrinsic or '',
                        'type':      kind or '',
                    }
                else:
                    raise RuntimeError(
                        "process_class_documentation: can not determine return "
                        "type for function")


# ---------------------------------------------------------------------------
# Output (XML dump)
# ---------------------------------------------------------------------------

def _dict_to_xml_element(tag, data):
    """Render `data` as nested XML under <tag>.

    Mimics XML::Simple's `XMLout(..., NoAttr => 1)` output — every key of a
    dict becomes a child element, lists become repeated children.
    """
    root = ET.Element(tag)
    _fill(root, data)
    return root


def _fill(elem, data):
    if isinstance(data, dict):
        for key, value in data.items():
            if isinstance(value, list):
                for item in value:
                    child = ET.SubElement(elem, key)
                    _fill(child, item)
            elif isinstance(value, dict):
                child = ET.SubElement(elem, key)
                _fill(child, value)
            else:
                child = ET.SubElement(elem, key)
                if value is not None:
                    child.text = str(value)
    elif isinstance(data, list):
        for item in data:
            child = ET.SubElement(elem, 'item')
            _fill(child, item)
    else:
        if data is not None:
            elem.text = str(data)


def _write_classes_xml(tree, classes):
    """Write `$BUILDPATH/<basename>.classes.xml` for this tree's .F90/.Inc
    source, skipping the `objects.nodes.components.Inc` special case.
    Mirrors ClassDocumentation.pm:173-180.
    """
    build_path = os.environ.get('BUILDPATH')
    name       = os.path.basename(tree.get('name') or '')
    if not build_path or not classes:
        return
    m = re.match(r'(.+)\.(F90|Inc)$', name)
    if not m:
        return
    if name == 'objects.nodes.components.Inc':
        return
    out_path = os.path.join(build_path, f"{m.group(1)}.classes.xml")
    root = _dict_to_xml_element('classes', classes)
    ET.ElementTree(root).write(out_path, encoding='utf-8', xml_declaration=True)
    _OUTPUT_PREVIOUS.extend(sorted(classes.keys()))


# ---------------------------------------------------------------------------
# Main dispatcher
# ---------------------------------------------------------------------------

def process_class_documentation(tree, options):
    """Mirrors Process_ClassDocumentation() from ClassDocumentation.pm."""
    if os.environ.get('GALACTICUS_BUILD_DOCS') != 'yes':
        return

    # Skip when called from a nested `process_tree` invocation — those are
    # `Generics._expand_subtree` and `FunctionClass._insert_and_write_output`
    # processing isolated subtrees that contain only PART of what we need
    # for documentation extraction (e.g. a generic-cloned `polyRankdouble`
    # type without its associated `doubleRank` / `doubleShape` function
    # nodes — those land in *separate* sibling clones).  Running here would
    # write a `<basename>.classes.xml` containing each clone's incomplete
    # record, only to be overwritten by the next clone's incomplete record;
    # waiting for the OUTER call lets us walk the fully-assembled tree and
    # match every method to its bound function in a single sweep.
    if is_inner_process_tree_call():
        return

    trees = _build_tree_list(tree)

    classes = {}
    function_list = []

    # Outer-tree source: child class type nodes that FunctionClass inserts
    # into THIS parent's tree (via `insert_pre_contains` of the
    # `code_content['module']['interfaces']` list it accumulated from
    # parsing each child class file) carry the CHILD source filename in
    # their `source` field, NOT the parent's.  Their function bodies
    # stay in separate submodule files, so walking those types here
    # would record them with empty `<arguments>` / `<argumentList>` /
    # `<type>` (no function nodes available to feed
    # `_process_function`).  Each child file gets its own preprocess.py
    # invocation that writes a complete record; let that own the
    # documentation for the child types.
    outer_source = (tree.get('source') or tree.get('name') or '')
    def _is_inserted_from_other_source_file(node):
        """True for a `type` node whose `source` is a *different real
        Fortran source file*.  Synthetic sources from FunctionClass /
        Generics inner parses (those starting with `Galacticus.`) and
        nodes with no `source` attribution still pass through, since
        their function bodies live in this same outer tree."""
        src = node.get('source') or ''
        if not src or src == outer_source:
            return False
        if src.startswith('Galacticus.'):
            return False
        return src.endswith('.F90') or src.endswith('.Inc')

    for t in trees:
        for node in walk_tree(t):
            ntype = node.get('type')
            if ntype == 'type':
                if _is_inserted_from_other_source_file(node):
                    continue
                class_name = node.get('name')
                class_record = classes.setdefault(class_name, {
                    'name': class_name,
                })
                _populate_class_descriptions(node, class_record)
                opener = node.get('opener') or ''
                m = re.search(r',\s*extends\s*\(\s*([a-zA-Z0-9_]+)\s*\)', opener)
                if m:
                    class_record['extends'] = m.group(1)
                _resolve_method_bindings(class_record)
                _compute_missing_and_generic(class_record)
            elif ntype == 'functionClass':
                # Synthesise a class record for the abstract `*Class` type
                # that FunctionClass auto-generates from this directive.
                # ClassDocumentation runs before FunctionClass in the topo
                # order, so the generated `type, abstract :: virialOrbitClass`
                # never reaches a `type` node walk in this pass — but the
                # `<functionClass>` directive carries every method we'd need
                # in its `<method name="…">` children, so we can populate the
                # parent class's `descriptions` directly.  Without this
                # entry the doc consumer's inheritance walk
                # (`scripts/doc/extractData.py`) cannot resolve the parent's
                # methods, and EVERY derived class lists them in
                # `<missingMethods>`.
                _populate_class_from_function_class_directive(node, classes)
            elif ntype in ('subroutine', 'function'):
                function_list.append(node)

    for function in function_list:
        _process_function(function, classes)

    # Scrub in-memory-only helpers and already-reported classes.
    for class_record in list(classes.values()):
        class_record.pop('functions', None)
        name = class_record.get('name')
        if name and name in _OUTPUT_PREVIOUS:
            classes.pop(name, None)

    _write_classes_xml(tree, classes)


# Register as a POSTPROCESS hook so we run AFTER every process hook
# (functionClass, generics, eventHooks, stateStorable, deepCopyActions, …)
# has finished expanding the tree.  Generic templates clone each sibling
# subtree independently, and FunctionClass inserts its auto-generated
# `<base>Class` type via `insert_after_node` — both expansions only become
# visible together in the OUTER tree once their respective process hooks
# have completed, so a process-time hook running in topo order would walk
# an incomplete tree.  The `is_inner_process_tree_call()` guard above
# additionally suppresses the postprocess pass during the nested
# `process_tree` calls those hooks make on isolated subtrees.
register_postprocess('classDocumentation', process_class_documentation)
