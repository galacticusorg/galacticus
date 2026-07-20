"""Provides a tree-based parser for Galacticus Fortran source files.

Andrew Benson (ported to Python 2026)
"""

import re
import os


from Galacticus.Build.FortranUtils import get_fortran_line
from Fortran.Utils import UNIT_OPENERS, UNIT_CLOSERS, check_no_parameterized_derived_type
from Galacticus.Build.SourceTree.Parse.Declarations import parse_declaration

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_file(filename):
    """Read a Fortran source file and return the root AST node."""
    with open(filename, 'r', errors='replace') as fh:
        content = fh.read()
    return parse_code(content, name=_name_relative_to_source(filename),
                      source=filename)


def _name_relative_to_source(filename):
    """Return `filename` relative to its innermost `source/` directory.

    The `source/` tree is organized hierarchically, so a bare basename is
    ambiguous (many subdirectories contain e.g. an `_class.F90`).  The file
    node's `name` — reported by `{introspection:location}` expansions and
    build error messages — therefore carries the path below `source/`
    (e.g. `intergalactic_medium/state/internal.F90`).  Files that do not
    live under a `source/` directory (generated includes under $BUILDPATH,
    test fixtures) keep their basename.
    """
    parts = os.path.normpath(filename).split(os.sep)
    for i in range(len(parts) - 2, -1, -1):
        if parts[i] == 'source':
            return '/'.join(parts[i + 1:])
    return os.path.basename(filename)


def parse_code(code, name='<string>', source=None, instrument=True):
    """Build an AST from a Fortran source string.

    Used by Process/Generics when it serializes a macro-expanded subtree and
    needs to re-parse it from its textual form.

    `instrument=True` (the default) tags every `{introspection:location}`
    placeholder with the line number on which it appears; downstream
    `process_source_introspection` then expands those into a full Fortran
    expression naming the surrounding scope.  Re-parses of generic-expanded
    or otherwise synthesised content should pass `instrument=False` so the
    line numbers aren't tagged a second time.
    """
    if instrument:
        from Galacticus.Build.SourceTree.Process.SourceIntrospection import (
            instrument as _instrument,
        )
        code = _instrument(code)
    root = {
        'type':       'file',
        'name':       name,
        'content':    code,
        'parent':     None,
        'firstChild': None,
        'sibling':    None,
        'source':     source if source is not None else name,
        'line':       0,
    }
    _build_tree(root)
    return root


def walk_tree(node):
    """Depth-first generator over all nodes in the tree.

    Yields every node exactly once
    in pre-order (parent before children).
    """
    yield node
    child = node.get('firstChild')
    while child:
        yield from walk_tree(child)
        child = child.get('sibling')


def children(node):
    """Return a list of direct child nodes."""
    result = []
    child = node.get('firstChild')
    while child:
        result.append(child)
        child = child.get('sibling')
    return result


# ---------------------------------------------------------------------------
# Internal helpers — tree construction
# ---------------------------------------------------------------------------

def _link_children(parent, child_list):
    """Attach child_list nodes as children of parent, linking siblings."""
    for i, child in enumerate(child_list):
        child['parent'] = parent
        child['sibling'] = child_list[i + 1] if i + 1 < len(child_list) else None
    parent['firstChild'] = child_list[0] if child_list else None


def _build_tree(root):
    """Parse root['content'] into child AST nodes, then recurse, running
    the directive/moduleUse/declaration parse passes over the result.
    """
    # Step 1: comment out embedded LaTeX (`!!{ ... !!}`) and XML (`!![ ... !!]`)
    # blocks by prefixing each body line with "!< ".
    # The markers themselves remain unmodified (they are
    # already valid Fortran comments because they begin with "!!"), but the
    # body lines need the "!< " prefix to be valid Fortran when serialized
    # straight to the compiler.  We store the commented content in code nodes
    # so that serialize() emits compilable Fortran.
    root['content'] = _comment_embedded(root['content'])

    # Step 2: split content into structural units.
    unit_children = _parse_units(root)
    _link_children(root, unit_children)

    # Step 3: run the parse passes over the whole tree.
    _pass_directives(root)
    _pass_module_uses(root)
    _pass_declarations(root)
    _pass_visibilities(root)
    _pass_openmp(root)

def _comment_embedded(content):
    """Comment out embedded LaTeX and XML blocks so that they are not processed as Fortran code"""

    import io
    fh = io.StringIO(content)

    code_commented = ""
    in_LaTeX = False
    in_XML = False
    for line in fh:
        # Detect the end of a LaTeX section and change state.
        if re.match(r'^\s*!!\}',line):
            in_LaTeX = False
        # Detect the end of an XML section and change state.
        if re.match(r'^\s*!!\]',line):
            in_XML = False
        # Comment out LaTeX and XML, unless it is already commented out.
        if ( in_LaTeX or in_XML ) and not re.match(r'^\s*\!<',line):
            line = "!< "+line
        code_commented = code_commented + line
        # Detect the start of a LaTeX section and change state.
        if re.match(r'^\s*!!\{',line):
            in_LaTeX = True
        # Detect the start of an XML section and change state.
        if re.match(r'^\s*!!\[',line):
            in_XML = True
    return code_commented
        
def _make_code_node(content, source, line, parent=None):
    return {
        'type':       'code',
        'content':    content,
        'parent':     parent,
        'firstChild': None,
        'sibling':    None,
        'source':     source,
        'line':       line,
    }


def _parse_units(parent, first_line=None):
    """Split parent['content'] into child nodes (clean implementation).

    parent['content'] is expected to already have any embedded LaTeX/XML
    blocks commented out by _comment_embedded() — every line stored in a
    code/opener/closer field is the same line we will later serialize, so
    the round-trip output is valid Fortran.

    `first_line` is the absolute (1-based) line number, in the original
    source file, of the first line of parent['content'].  When omitted it
    is derived from parent['line'], which for a root file node (line 0)
    yields line 1.  Every child node's `line` is set to the absolute
    1-based line number of its first line — these feed the `.lmap`
    line-number mappings and the `{introspection:location}` expansions,
    so they must survive arbitrarily deep nesting.
    """
    content   = parent.get('content', '')
    source    = parent.get('source', 'unknown')
    if first_line is None:
        first_line = parent.get('line', 0) + 1

    import io
    fh = io.StringIO(content)

    # Each stack entry: [node_dict, inner_lines_list].  inner_lines_list
    # holds (absolute_line, raw_text) tuples for buffered source lines and
    # ('\x00NODE\x00', node) sentinels for pre-built child nodes.
    stack         = []
    top_children  = []
    raw_code_buf  = []
    raw_code_line = first_line
    current_line  = first_line   # absolute line number of the next unread line

    def flush_code(destination_list):
        nonlocal raw_code_buf
        if raw_code_buf:
            destination_list.append(
                _make_code_node(''.join(raw_code_buf), source, raw_code_line))
            raw_code_buf  = []

    while True:
        raw_line, processed_line, _ = get_fortran_line(fh)
        if not raw_line and not processed_line:
            break

        # Reject parameterized derived types loudly: the parser and downstream
        # generators cannot represent type parameters and would silently drop a
        # PDT from generated code.  See issue #114.
        check_no_parameterized_derived_type(processed_line)

        n_newlines    = raw_line.count('\n')
        start_line    = current_line
        line_after    = current_line + n_newlines

        # ---- check for a closer matching the innermost open unit ----
        if stack:
            top_node = stack[-1][0]
            closer_re = UNIT_CLOSERS.get(top_node['type'])
            if closer_re and closer_re.match(processed_line):
                top_node['closer'] = raw_line
                # Recurse into the accumulated inner content.
                # inner_lines may contain pre-built sentinel nodes interspersed
                # with raw line chunks, so use _children_from_mixed_lines.
                inner_lines = stack[-1][1]
                if inner_lines:
                    inner_children = _children_from_mixed_lines(inner_lines, top_node)
                    _link_children(top_node, inner_children)
                stack.pop()
                # Deliver the closed node to the new top scope.
                if stack:
                    # Flush raw code inside outer unit, then add closed node.
                    # (raw_code_buf accumulates top-level code; it would be
                    # empty here in practice, but handle defensively.)
                    if raw_code_buf:
                        stack[-1][1].append((raw_code_line, ''.join(raw_code_buf)))
                        raw_code_buf  = []
                    # The closed node goes into the outer unit's inner lines as
                    # a pre-built node object (we mark it with a sentinel key).
                    stack[-1][1].append(('\x00NODE\x00', top_node))
                else:
                    flush_code(top_children)
                    top_children.append(top_node)
                current_line = line_after
                continue

        # ---- check for an opener ----
        matched_opener = False
        for unit_type, spec in UNIT_OPENERS.items():
            m = spec['regex'].match(processed_line)
            if not m:
                continue

            # Disambiguate: "module procedure/function/subroutine" is NOT a module.
            if unit_type == 'module' and re.match(
                    r'^\s*module\s+(?:procedure|function|subroutine)\b',
                    processed_line, re.IGNORECASE):
                continue

            # Extract the unit name.  A negative `name_idx` means "no name"
            # (used by self-closing markers like `contains`).
            name_idx  = spec.get('unit_name', 0)
            groups    = m.groups()
            if name_idx < 0 or name_idx >= len(groups) or not groups[name_idx]:
                unit_name = ''
            else:
                unit_name = groups[name_idx].strip()

            # Flush pending raw code to the appropriate destination.
            if stack:
                if raw_code_buf:
                    stack[-1][1].append((raw_code_line, ''.join(raw_code_buf)))
                    raw_code_buf  = []
            else:
                flush_code(top_children)

            # Build the new node.
            node = {
                'type':       unit_type,
                'name':       unit_name,
                'opener':     raw_line,
                'parent':     None,
                'firstChild': None,
                'sibling':    None,
                'source':     source,
                'line':       start_line,
            }

            # moduleProcedure and `contains` are self-closing: the former has
            # no `end procedure` counterpart in the files we parse, and the
            # latter marks the transition from declarations to subprograms
            # inside a module/function/subroutine (its scope closes with the
            # enclosing unit, which we do not need to model structurally —
            # treating `contains` as a sibling marker keeps parsing simple and
            # serializes identically).
            if unit_type in ('moduleProcedure', 'contains'):
                if unit_type == 'moduleProcedure':
                    names = [n.strip() for n in re.split(r'\s*,\s*', unit_name) if n.strip()]
                    node['names'] = names
                if stack:
                    stack[-1][1].append(('\x00NODE\x00', node))
                else:
                    top_children.append(node)
            else:
                stack.append((node, []))

            matched_opener = True
            current_line   = line_after
            break

        if matched_opener:
            continue

        # ---- plain code line ----
        if stack:
            stack[-1][1].append((start_line, raw_line))
        else:
            if not raw_code_buf:
                raw_code_line = start_line
            raw_code_buf.append(raw_line)

        current_line = line_after

    # Close any units left open (e.g. partial files, included fragments).
    while stack:
        top_node   = stack[-1][0]
        inner_lines = stack[-1][1]
        if inner_lines:
            inner_children = _children_from_mixed_lines(inner_lines, top_node)
            _link_children(top_node, inner_children)
        stack.pop()
        flush_code(top_children)
        top_children.append(top_node)

    flush_code(top_children)

    # Second pass: resolve sentinel NODE entries inside inner_lines and rebuild
    # children lists that contain pre-built nodes.
    # (These were injected when a sub-unit was closed while an outer unit was open.)
    _resolve_sentinels(top_children)

    return top_children


def _resolve_sentinels(node_list):
    """Walk node_list and fix up any inner buffers that contain sentinel nodes."""
    # No-op: sentinel resolution is handled inline by _children_from_mixed_lines.
    pass


def _children_from_mixed_lines(inner_lines, parent):
    """Build child node list from inner_lines.

    inner_lines is a list where each element is either:
      - an (absolute_line, text) 2-tuple — a raw Fortran source chunk
        (already commented for embedded LaTeX/XML by _comment_embedded())
        together with the 1-based line number of its first line in the
        original source file, or
      - a ('\x00NODE\x00', node) 2-tuple — a pre-built child node injected
        when a sub-unit closed while an outer unit was still open.

    Runs of raw chunks are joined and re-parsed via _parse_units(), passing
    the absolute line number of the run's first line so that the resulting
    child nodes carry correct original-source line numbers; pre-built nodes
    pass through unchanged.
    """
    children = []
    line_buf = []
    run_line = None
    source   = parent.get('source', 'unknown')

    def _flush():
        nonlocal run_line
        if not line_buf:
            return
        dummy = {
            'content': ''.join(line_buf),
            'source':  source,
            'line':    run_line,
        }
        children.extend(_parse_units(dummy, first_line=run_line))
        line_buf.clear()
        run_line = None

    for item in inner_lines:
        if not (isinstance(item, tuple) and len(item) == 2):
            raise TypeError(
                f"_children_from_mixed_lines: expected an (line, text) or "
                f"sentinel tuple, got {type(item).__name__!r}: {item!r}"
            )
        if item[0] == '\x00NODE\x00':
            _flush()
            children.append(item[1])
        else:
            if not line_buf:
                run_line = item[0]
            line_buf.append(item[1])
    _flush()
    return children


# ---------------------------------------------------------------------------
# Parse pass 1: directives
# ---------------------------------------------------------------------------



def _pass_directives(tree):
    """Walk code nodes extracting `!![…!!]` XML directive blocks into directive nodes.

    Delegates to Galacticus.Build.SourceTree.Parse.Directives.parse_directives().
    Lazy import avoids a circular dependency at module-load time.
    """
    from Galacticus.Build.SourceTree.Parse.Directives import parse_directives
    parse_directives(tree)


def replace_node(old_node, new_nodes):
    """Replace old_node in the tree with new_nodes."""
    if not new_nodes:
        return
    parent = old_node.get('parent')
    if parent is None:
        return

    # Find old_node in parent's child list.
    prev = None
    child = parent.get('firstChild')
    while child:
        if child is old_node:
            break
        prev  = child
        child = child.get('sibling')

    if child is None:
        return  # old_node not found

    # Wire new_nodes into the sibling chain.
    for i, n in enumerate(new_nodes):
        n['parent'] = parent
        n['sibling'] = new_nodes[i + 1] if i + 1 < len(new_nodes) else old_node.get('sibling')

    if prev is None:
        parent['firstChild'] = new_nodes[0]
    else:
        prev['sibling'] = new_nodes[0]


def insert_before_node(node, new_nodes):
    """Insert new_nodes as siblings immediately before node in the parent's child list."""
    parent = node.get('parent')
    if parent is None:
        raise ValueError("insert_before_node: cannot insert before a root node")
    for i, n in enumerate(new_nodes):
        n['parent']  = parent
        n['sibling'] = new_nodes[i + 1] if i + 1 < len(new_nodes) else node
    if parent.get('firstChild') is node:
        parent['firstChild'] = new_nodes[0]
    else:
        prev = parent['firstChild']
        while prev.get('sibling') is not node:
            prev = prev['sibling']
        prev['sibling'] = new_nodes[0]


def insert_after_node(node, new_nodes):
    """Insert new_nodes as siblings immediately after node in the parent's child list."""
    parent = node.get('parent')
    if parent is None:
        raise ValueError("insert_after_node: cannot insert after a root node")
    tail = node.get('sibling')
    for i, n in enumerate(new_nodes):
        n['parent']  = parent
        n['sibling'] = new_nodes[i + 1] if i + 1 < len(new_nodes) else tail
    node['sibling'] = new_nodes[0]


def _find_child_by_type(node, child_type):
    """Return the first direct child of `node` with the given type, or None."""
    child = node.get('firstChild')
    while child is not None:
        if child.get('type') == child_type:
            return child
        child = child.get('sibling')
    return None


def _last_child(node):
    """Return the last direct child of `node`, or None if it has none."""
    child = node.get('firstChild')
    if child is None:
        return None
    while child.get('sibling') is not None:
        child = child['sibling']
    return child


def insert_pre_contains(node, new_nodes):
    """Insert new_nodes before the `contains` child of node.

    If no `contains` exists, the new nodes are appended at the end of
    node's children.
    """
    contains_node = _find_child_by_type(node, 'contains')
    if contains_node is not None:
        insert_before_node(contains_node, new_nodes)
        return
    last = _last_child(node)
    if last is None:
        _link_children(node, new_nodes)
    else:
        insert_after_node(last, new_nodes)


def insert_post_contains(node, new_nodes):
    """Insert new_nodes as the first children placed after the `contains`
    marker of node.

    Auto-creates a `contains` node if the parent does not already have one.
    In our representation `contains` is a self-closing sibling marker, so
    "prepend as children of contains" becomes "insert immediately after
    contains" — the serialized output is identical either way.
    """
    contains_node = _find_child_by_type(node, 'contains')
    if contains_node is None:
        contains_node = {
            'type':       'contains',
            'opener':     'contains\n',
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':     node.get('source', 'unknown'),
            'line':       node.get('line', 0),
        }
        last = _last_child(node)
        if last is None:
            _link_children(node, [contains_node])
        else:
            insert_after_node(last, [contains_node])
    insert_after_node(contains_node, new_nodes)


def prepend_child_to_node(node, new_nodes):
    """Insert new_nodes as the first children of node."""
    first_child = node.get('firstChild')
    if first_child is None:
        _link_children(node, new_nodes)
    else:
        insert_before_node(first_child, new_nodes)


def set_visibility(node, unit_name, visibility):
    """Ensure `unit_name` is listed in node's `public` or `private` visibility.

    Auto-creates
    a `visibility` child node if one does not exist, placing it after any
    existing `moduleUse` child — Fortran requires visibility statements to
    follow `use` statements — and
    regenerates its `firstChild` code content with the sorted `public` /
    `private` lists.
    """
    if visibility not in ('public', 'private'):
        raise ValueError(
            f"set_visibility: visibility must be 'public' or 'private', "
            f"got {visibility!r}")

    visibility_node = _find_child_by_type(node, 'visibility')
    if visibility_node is None:
        visibility_node = {
            'type':       'visibility',
            'visibility': {},
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':     node.get('source', 'unknown'),
            'line':       node.get('line', 0),
        }
        visibility_node['firstChild'] = {
            'type':       'code',
            'content':    '',
            'parent':     visibility_node,
            'firstChild': None,
            'sibling':    None,
            'source':     visibility_node['source'],
            'line':       visibility_node['line'],
        }
        # Find the last moduleUse child; visibilities must come after it.
        child = node.get('firstChild')
        module_use_node = None
        while child is not None:
            if child.get('type') == 'moduleUse':
                module_use_node = child
            child = child.get('sibling')
        if module_use_node is not None:
            insert_after_node(module_use_node, [visibility_node])
        else:
            prepend_child_to_node(node, [visibility_node])

    vis_dict = visibility_node.setdefault('visibility', {})
    vis_dict.setdefault(visibility, {})[unit_name] = True

    # Rebuild the visibility code in fixed private-then-public order so the
    # regenerated code is deterministic.
    content = ''
    for level in ('private', 'public'):
        entries = vis_dict.get(level)
        if entries is None:
            continue
        content += "  " + level
        if entries:
            content += " :: " + ", ".join(sorted(entries.keys()))
        content += "\n"
    visibility_node['firstChild']['content'] = content + "\n"


# ---------------------------------------------------------------------------
# Parse pass 2: module uses
# ---------------------------------------------------------------------------

def _pass_module_uses(tree):
    """Walk code nodes extracting 'use' statements into moduleUse nodes.

    Delegates to Galacticus.Build.SourceTree.Parse.ModuleUses.parse_module_uses(),
    which produces the full-featured node structure (moduleOrder, openMP,
    conditions, firstChild with raw text).  The lazy import avoids a circular
    dependency at module-load time.
    """
    from Galacticus.Build.SourceTree.Parse.ModuleUses import parse_module_uses
    parse_module_uses(tree)


# ---------------------------------------------------------------------------
# Parse pass 3: declarations
# ---------------------------------------------------------------------------

def _pass_declarations(tree):
    """Walk code nodes extracting declaration lines into declaration nodes."""
    nodes_to_replace = []

    for node in walk_tree(tree):
        if node['type'] != 'code':
            continue

        content   = node.get('content', '')
        new_nodes = []
        code_buf  = []
        decl_buf  = []
        decls     = []
        implicit_none = False
        code_line     = node['line']   # absolute line of first buffered code line
        decl_line     = node['line']   # absolute line of first buffered decl line
        current_line  = node['line']   # absolute line of the next unread line

        def flush_code():
            if code_buf:
                new_nodes.append(_make_code_node(
                    ''.join(code_buf), node['source'], code_line))
                code_buf.clear()

        def flush_decls():
            nonlocal implicit_none
            if decl_buf or decls:
                dn = {
                    'type':         'declaration',
                    'implicitNone': implicit_none,
                    'declarations': list(decls),
                    'parent':       None,
                    'firstChild':   None,
                    'sibling':      None,
                    'source':       node['source'],
                    'line':         decl_line,
                }
                # Store raw text in firstChild so serialize() can reconstruct
                # the original source.
                dn['firstChild'] = {
                    'type':       'code',
                    'content':    ''.join(decl_buf),
                    'parent':     dn,
                    'sibling':    None,
                    'firstChild': None,
                    'source':     node['source'],
                    'line':       decl_line,
                }
                new_nodes.append(dn)
                decl_buf.clear()
                decls.clear()
                implicit_none = False

        import io
        fh = io.StringIO(content)
        while True:
            raw_line, processed_line, _ = get_fortran_line(fh)
            if not raw_line and not processed_line:
                break

            start_line = current_line

            is_decl = False
            if re.match(r'^\s*implicit\s+none\s*$', processed_line, re.IGNORECASE):
                is_decl       = True
                implicit_none = True

            decl = parse_declaration(processed_line)
            if decl:
                is_decl = True
                decl['line'] = start_line

            if is_decl:
                flush_code()
                if not decl_buf:
                    decl_line = start_line
                decl_buf.append(raw_line)
                if decl:
                    decls.append(decl)
            else:
                flush_decls()
                if not code_buf:
                    code_line = start_line
                code_buf.append(raw_line)

            current_line = start_line + raw_line.count('\n')

        flush_decls()
        flush_code()

        if not (len(new_nodes) == 1 and new_nodes[0].get('type') == 'code'
                and new_nodes[0].get('content') == content):
            nodes_to_replace.append((node, new_nodes))

    for old_node, new_nodes in nodes_to_replace:
        replace_node(old_node, new_nodes)


# ---------------------------------------------------------------------------
# Parse pass 4: visibilities
# ---------------------------------------------------------------------------

def _pass_visibilities(tree):
    """Walk code nodes extracting public/private statements into visibility nodes.

    Delegates to Galacticus.Build.SourceTree.Parse.Visibilities.parse_visibilities().
    Lazy import avoids a circular dependency at module-load time.
    """
    from Galacticus.Build.SourceTree.Parse.Visibilities import parse_visibilities
    parse_visibilities(tree)


# ---------------------------------------------------------------------------
# Parse pass 5: OpenMP
# ---------------------------------------------------------------------------

def _pass_openmp(tree):
    """Walk code nodes extracting !$omp parallel directives into openMP nodes.

    Delegates to Galacticus.Build.SourceTree.Parse.OpenMP.parse_openmp().
    Lazy import avoids a circular dependency at module-load time.
    """
    from Galacticus.Build.SourceTree.Parse.OpenMP import parse_openmp
    parse_openmp(tree)


# ---------------------------------------------------------------------------
# Serialization
# ---------------------------------------------------------------------------

def serialize(node, annotate=False, strip_mappings=False):
    """Reconstruct Fortran source text from an AST node and its siblings.

    Parameters
    ----------
    node : dict
        Root AST node.
    annotate : bool, default False
        When True, emit `!--> <origLine> <outLine> "<source>"` line-number
        mapping comments ahead of each node's serialized content.  When
        False, the serialized output is pure Fortran.  (The default is False
        because most callers just want plain source text.)
    strip_mappings : bool, default False
        When True, the annotations are removed from the returned source and
        collected into a second string.  In that mode the function returns a
        `(source, mappings)` 2-tuple.  When False (the default), the function
        returns a single `source` string — matching the behaviour used by
        `Generics` / `FunctionClass` callers that do not need mappings.

    The algorithm is a sibling-chain walk that recurses into firstChild:
      - code nodes      → emit content directly
      - all other nodes → emit opener (if any) + serialize(firstChild) + closer (if any)
    """
    result, mappings = _serialize(node, annotate=annotate, strip_mappings=strip_mappings)
    if strip_mappings:
        return result, mappings
    return result


_MAPPING_RE = re.compile(r'^!-->\s+(\d+)\s+(\d+)\s+"(.+)"')


def _serialize(node, annotate, strip_mappings):
    """Internal helper: always returns (source, mappings) and tracks the
    cumulative output line number across the sibling chain.
    """
    serialization = ''
    mappings      = ''
    line_number   = 0
    current       = node

    while current:
        # Emit line-mapping annotation for this node.  In strip mode the
        # recorded output line is 1-based (the node's first line in the
        # stripped output), matching the `line_number + 1` used when inner
        # annotations are collected below — postprocess.py's remap
        # arithmetic needs every entry on the same basis.
        if annotate and 'source' in current and 'line' in current:
            if strip_mappings:
                mappings += f'!--> {current["line"]} {line_number + 1} "{current["source"]}"\n'
            else:
                serialization += f'!--> {current["line"]} {line_number} "{current["source"]}"\n'
                line_number   += 1

        # Serialize this node's own content.
        if current.get('type') == 'code':
            node_text = current.get('content') or ''
        else:
            node_text = current.get('opener') or ''
            child     = current.get('firstChild')
            if child:
                # Children never strip mappings themselves — we do that at the
                # outermost level as we fold their text back in.
                child_text, _ = _serialize(child, annotate=annotate, strip_mappings=False)
                node_text += child_text
            node_text += current.get('closer') or ''

        if strip_mappings:
            stripped = ''
            for line in node_text.splitlines(keepends=True):
                m = _MAPPING_RE.match(line)
                if m:
                    mappings += f'!--> {m.group(1)} {line_number + 1} "{m.group(3)}"\n'
                else:
                    line_number += 1
                    stripped    += line
            node_text = stripped
        else:
            line_number += node_text.count('\n')

        serialization += node_text
        current        = current.get('sibling')

    return serialization, mappings


def analyze_tree(tree, options=None):
    """Run every registered analyze hook on the tree.

    The analyze
    hook registry lives in Galacticus.Build.SourceTree.Process so that Analyze
    submodules can register themselves at import time; if no analyze modules
    have been imported, this is a no-op.
    """
    from Galacticus.Build.SourceTree.Process import ANALYZE_HOOKS
    if options is None:
        options = {}
    for name in sorted(ANALYZE_HOOKS.keys()):
        ANALYZE_HOOKS[name](tree, options)
    return tree
