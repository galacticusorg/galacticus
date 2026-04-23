# Provides a tree-based parser for Galacticus Fortran source files.
# Andrew Benson (ported to Python 2026)
#
# Mirrors the subset of perl/Galacticus/Build/SourceTree.pm used by the
# Galacticus build scripts, together with the Parse::Directives and
# Parse::ModuleUses parsers that libraryInterfaces.py requires.

import re
import os
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from build.fortran_utils import get_fortran_line
from Fortran.Utils import UNIT_OPENERS, UNIT_CLOSERS
from Galacticus.Build.SourceTree.Parse.Declarations import parse_declaration

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_file(filename):
    """Read a Fortran source file and return the root AST node.

    Mirrors Perl Galacticus::Build::SourceTree::ParseFile().
    """
    with open(filename, 'r', errors='replace') as fh:
        content = fh.read()
    return parse_code(content, name=os.path.basename(filename), source=filename)


def parse_code(code, name='<string>', source=None):
    """Build an AST from a Fortran source string.

    Mirrors Perl Galacticus::Build::SourceTree::ParseCode(code, fileName).
    Used by Process/Generics when it serializes a macro-expanded subtree and
    needs to re-parse it from its textual form.
    """
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

    Mirrors the Perl Walk_Tree loop idiom.  Yields every node exactly once
    in pre-order (parent before children).
    """
    yield node
    child = node.get('firstChild')
    while child:
        yield from walk_tree(child)
        child = child.get('sibling')


def children(node):
    """Return a list of direct child nodes.

    Mirrors Perl Galacticus::Build::SourceTree::Children().
    """
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
    """Parse root['content'] into child AST nodes, then recurse.

    Mirrors Perl Build_Children + Parse_Unit + the directive/moduleUse/
    declaration parse hooks.
    """
    # Step 1: comment out LaTeX and XML blocks so that _parse_units() is not confused
    # by documentation text that could look like structural Fortran keywords.  Keep the
    # original content so _parse_units() can store unmodified raw lines in code nodes,
    # which is required for correct round-trip serialization.
    original_content   = root['content']
    root['content']    = _comment_embedded(original_content)
    root['_orig_content'] = original_content   # parallel reader used by _parse_units

    # Step 2: split raw content into structural units.
    unit_children = _parse_units(root)
    del root['_orig_content']
    _link_children(root, unit_children)

    # Step 3: run the parse passes over the whole tree (matches the set of
    # parseHooks registered in perl/Galacticus/Build/SourceTree.pm).
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


def _parse_units(parent):
    """Split parent['content'] into child nodes (clean implementation)."""
    content          = parent.get('content', '')
    # _orig_content is set by _build_tree() before calling _comment_embedded(); for
    # recursively-created dummy parents (from _children_from_mixed_lines), the content
    # is already the original because we accumulate from fh_orig below.
    original_content = parent.get('_orig_content', content)
    source           = parent.get('source', 'unknown')
    base_line        = parent.get('line',   0)

    import io
    fh      = io.StringIO(content)           # processed (possibly commented) — for detection
    fh_orig = io.StringIO(original_content)  # original — for raw line storage

    # Each stack entry: [node_dict, inner_lines_list]
    stack         = []
    top_children  = []
    raw_code_buf  = []
    raw_code_line = base_line
    current_line  = base_line

    def flush_code(destination_list):
        nonlocal raw_code_buf, raw_code_line
        if raw_code_buf:
            destination_list.append(
                _make_code_node(''.join(raw_code_buf), source, raw_code_line))
            raw_code_buf  = []
            raw_code_line = current_line

    while True:
        raw_line,      processed_line, _ = get_fortran_line(fh)
        raw_orig_line, _,              _ = get_fortran_line(fh_orig)
        if not raw_line and not processed_line:
            break

        n_newlines    = raw_line.count('\n')
        line_after    = current_line + n_newlines

        # ---- check for a closer matching the innermost open unit ----
        if stack:
            top_node = stack[-1][0]
            closer_re = UNIT_CLOSERS.get(top_node['type'])
            if closer_re and closer_re.match(processed_line):
                top_node['closer'] = raw_orig_line   # store original
                # Recurse into the accumulated inner content.
                # inner_lines may contain pre-built sentinel nodes interspersed
                # with raw string chunks, so use _children_from_mixed_lines.
                inner_lines = stack[-1][1]
                if inner_lines:
                    inner_children = _children_from_mixed_lines(inner_lines, top_node)
                    _link_children(top_node, inner_children)
                stack.pop()
                # Deliver the closed node to the new top scope.
                if stack:
                    # Flush raw code inside outer unit, then add closed node.
                    # (raw_code_buf is accumulates top-level code; it would be
                    # empty here in practice, but handle defensively.)
                    if raw_code_buf:
                        s = ''.join(raw_code_buf)
                        stack[-1][1].append((s, s))   # same for modified/original
                        raw_code_buf  = []
                        raw_code_line = line_after
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

            # Extract the unit name.
            name_idx  = spec.get('unit_name', 0)
            groups    = m.groups()
            unit_name = (groups[name_idx].strip()
                         if name_idx < len(groups) and groups[name_idx]
                         else '')

            # Flush pending raw code to the appropriate destination.
            if stack:
                if raw_code_buf:
                    stack[-1][1].append(''.join(raw_code_buf))
                    raw_code_buf  = []
                    raw_code_line = line_after
            else:
                flush_code(top_children)

            # Build the new node.
            node = {
                'type':       unit_type,
                'name':       unit_name,
                'opener':     raw_orig_line,   # store original (unmodified) opener
                'parent':     None,
                'firstChild': None,
                'sibling':    None,
                'source':     source,
                'line':       current_line,
            }

            # moduleProcedure is self-closing (no "end procedure" counterpart
            # in the files we parse).
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
        # Store (modified, original) pair so that _children_from_mixed_lines can
        # pass modified content to _parse_units for structural detection while
        # keeping original content for code-node storage.
        if stack:
            stack[-1][1].append((raw_line, raw_orig_line))
        else:
            raw_code_buf.append(raw_orig_line)

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
    """Build child node list from inner_lines (mixed line-pairs and sentinel tuples).

    inner_lines is a list where each element is either:
      - a (modified_str, orig_str) 2-tuple — a raw Fortran source pair, or
      - a ('\x00NODE\x00', node)   2-tuple — a pre-built child node.

    Line-pair runs are joined separately: modified content is passed as 'content'
    (used for structural unit detection) and original content as '_orig_content'
    (used for code-node raw-text storage).  This mirrors the parallel-reader
    approach in _parse_units at the top level.
    """
    children = []
    mod_buf  = []
    orig_buf = []
    source   = parent.get('source', 'unknown')
    line     = parent.get('line',   0)

    def _flush():
        if not mod_buf:
            return
        mod_text  = ''.join(mod_buf)
        orig_text = ''.join(orig_buf)
        dummy = {
            'content':       mod_text,
            '_orig_content': orig_text,
            'source':        source,
            'line':          line,
        }
        children.extend(_parse_units(dummy))
        mod_buf.clear()
        orig_buf.clear()

    for item in inner_lines:
        if isinstance(item, tuple) and item[0] == '\x00NODE\x00':
            _flush()
            children.append(item[1])  # pre-built node
        else:
            # Must be a (modified_str, orig_str) 2-tuple; plain strings are not
            # accepted — they would silently index as character sequences.
            if (
                not isinstance(item, tuple)
                or len(item) != 2
                or not isinstance(item[0], str)
                or not isinstance(item[1], str)
            ):
                raise TypeError(
                    f"_children_from_mixed_lines: expected a (modified_str, orig_str) "
                    f"2-tuple, got {type(item).__name__!r}: {item!r}"
                )
            mod_buf.append(item[0])
            orig_buf.append(item[1])
    _flush()
    return children


# ---------------------------------------------------------------------------
# Parse pass 1: directives  (mirrors Parse::Directives)
# ---------------------------------------------------------------------------



def _pass_directives(tree):
    """Walk code nodes extracting `!![…!!]` XML directive blocks into directive nodes.

    Delegates to Galacticus.Build.SourceTree.Parse.Directives.parse_directives().
    Lazy import avoids a circular dependency at module-load time.
    """
    from Galacticus.Build.SourceTree.Parse.Directives import parse_directives
    parse_directives(tree)


def replace_node(old_node, new_nodes):
    """Replace old_node in the tree with new_nodes.

    Mirrors Perl Galacticus::Build::SourceTree::ReplaceNode().
    """
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
    """Insert new_nodes as siblings immediately before node in the parent's child list.

    Mirrors Perl Galacticus::Build::SourceTree::InsertBeforeNode().
    """
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
    """Insert new_nodes as siblings immediately after node in the parent's child list.

    Mirrors Perl Galacticus::Build::SourceTree::InsertAfterNode().
    """
    parent = node.get('parent')
    if parent is None:
        raise ValueError("insert_after_node: cannot insert after a root node")
    tail = node.get('sibling')
    for i, n in enumerate(new_nodes):
        n['parent']  = parent
        n['sibling'] = new_nodes[i + 1] if i + 1 < len(new_nodes) else tail
    node['sibling'] = new_nodes[0]


# ---------------------------------------------------------------------------
# Parse pass 2: module uses  (mirrors Parse::ModuleUses)
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
# Parse pass 3: declarations  (mirrors Parse::Declarations)
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

        def flush_code():
            if code_buf:
                new_nodes.append(_make_code_node(
                    ''.join(code_buf), node['source'], node['line']))
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
                    'line':         node['line'],
                }
                # Store raw text in firstChild so serialize() can reconstruct
                # the original source — mirrors Perl Parse_Declarations behaviour.
                dn['firstChild'] = {
                    'type':       'code',
                    'content':    ''.join(decl_buf),
                    'parent':     dn,
                    'sibling':    None,
                    'firstChild': None,
                    'source':     node['source'],
                    'line':       node['line'],
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

            is_decl = False
            if re.match(r'^\s*implicit\s+none\s*$', processed_line, re.IGNORECASE):
                is_decl       = True
                implicit_none = True

            decl = parse_declaration(processed_line)
            if decl:
                is_decl = True

            if is_decl:
                flush_code()
                decl_buf.append(raw_line)
                if decl:
                    decls.append(decl)
            else:
                flush_decls()
                code_buf.append(raw_line)

        flush_decls()
        flush_code()

        if not (len(new_nodes) == 1 and new_nodes[0].get('type') == 'code'
                and new_nodes[0].get('content') == content):
            nodes_to_replace.append((node, new_nodes))

    for old_node, new_nodes in nodes_to_replace:
        replace_node(old_node, new_nodes)


# ---------------------------------------------------------------------------
# Parse pass 4: visibilities  (mirrors Parse::Visibilities)
# ---------------------------------------------------------------------------

def _pass_visibilities(tree):
    """Walk code nodes extracting public/private statements into visibility nodes.

    Delegates to Galacticus.Build.SourceTree.Parse.Visibilities.parse_visibilities().
    Lazy import avoids a circular dependency at module-load time.
    """
    from Galacticus.Build.SourceTree.Parse.Visibilities import parse_visibilities
    parse_visibilities(tree)


# ---------------------------------------------------------------------------
# Parse pass 5: OpenMP  (mirrors Parse::OpenMP)
# ---------------------------------------------------------------------------

def _pass_openmp(tree):
    """Walk code nodes extracting !$omp parallel directives into openMP nodes.

    Delegates to Galacticus.Build.SourceTree.Parse.OpenMP.parse_openmp().
    Lazy import avoids a circular dependency at module-load time.
    """
    from Galacticus.Build.SourceTree.Parse.OpenMP import parse_openmp
    parse_openmp(tree)


# ---------------------------------------------------------------------------
# Serialization  (mirrors Perl Galacticus::Build::SourceTree::Serialize)
# ---------------------------------------------------------------------------

def serialize(node):
    """Reconstruct Fortran source text from an AST node and its siblings.

    Mirrors Perl Galacticus::Build::SourceTree::Serialize(node, annotate => 0).

    The algorithm is a sibling-chain walk that recurses into firstChild:
      - code nodes      → emit content directly
      - all other nodes → emit opener (if any) + serialize(firstChild) + closer (if any)

    This works for every node type produced by the parser:
      - structural units (module, subroutine, …): have opener/closer
      - moduleUse nodes: firstChild holds the raw/reformatted use text
      - declaration nodes: firstChild holds the raw declaration text
      - directive nodes: firstChild holds the raw !![…!!] block text
    """
    result  = ""
    current = node
    while current:
        if current.get('type') == 'code':
            result += current.get('content') or ''
        else:
            result += current.get('opener') or ''
            child = current.get('firstChild')
            if child:
                result += serialize(child)
            result += current.get('closer') or ''
        current = current.get('sibling')
    return result
