"""Resolve Galacticus parameter files into a clean, ready-to-run form.

This is the Python re-implementation of the file-level transformations that
Galacticus applies when it reads a parameter file, so that a fully-resolved file
can be handed to ``Galacticus.exe`` (see the implementation plan in
``~/.claude/plans/galacticus-parameter-resolver.md``).

Stage 1 (this module, initial version) implements the two transformations that
Galacticus bakes eagerly into the DOM and that the
``--output-processed-parameters`` oracle can certify:

* **XInclude** -- ``<xi:include href=... xpointer="xpointer(parameters/*)"/>`` is
  replaced by the element children of the referenced document's ``<parameters>``
  root (Galacticus' own limited XInclude; mirrors ``source/utility/IO/XML.F90``).
* **Change files** -- the ``<changes>`` operation set
  (``source/utility/input_parameters.F90``) applied, in command-line order, each
  ``<change>`` in document order, to the post-XInclude tree.

NOT handled here (by design -- see the plan):

* ``id``/``idRef`` are validated but NEVER dereferenced (they are runtime
  pointers -> a shared object; dereferencing would create separate copies).
  Reference validation + conditionals are Stage 2.
* ``=[...]`` math expressions are left INTACT for Galacticus to evaluate at
  runtime (deferred).

Processing order matches Galacticus: XInclude -> changes.

Uses ``lxml`` (already a dependency of the ``galacticus`` launcher).
"""

import copy
import os
import re

from lxml import etree


_XINCLUDE_NS = 'http://www.w3.org/2001/XInclude'
_XI_INCLUDE = f'{{{_XINCLUDE_NS}}}include'
_XI_FALLBACK = f'{{{_XINCLUDE_NS}}}fallback'
_MAX_INCLUDE_DEPTH = 64

# Galacticus accepts only `xpointer(<root>/*)` (select the element children of the
# referenced document's root) or no xpointer (the whole root element).
_XPOINTER_CHILDREN_RE = re.compile(r'^xpointer\(\s*/?(\w+)/\*\s*\)$')


class ResolveError(Exception):
    """A fatal resolution error (mirrors a Galacticus ``Error_Report``)."""


# --- load / save ------------------------------------------------------------

def load(path):
    """Parse a parameter or changes file into an lxml ``ElementTree``.

    Comments and whitespace are preserved so the emitted file stays readable.
    """
    parser = etree.XMLParser(remove_blank_text=False, remove_comments=False,
                             resolve_entities=False)
    return etree.parse(path, parser)


def to_bytes(tree):
    """Serialize a resolved tree to bytes, with an XML declaration."""
    root = tree.getroot() if isinstance(tree, etree._ElementTree) else tree
    return etree.tostring(root, encoding='UTF-8', xml_declaration=True)


def _is_element(node):
    """True for true element nodes (excludes comments and processing instructions)."""
    return isinstance(node.tag, str)


def _element_children(node):
    """Element children only -- mirrors Fortran ``XML_Get_Child_Elements``."""
    return [child for child in node if _is_element(child)]


# --- path resolution --------------------------------------------------------

def resolve_path(node, path, root):
    """Return the list of nodes matching ``path`` (Galacticus legacy semantics).

    Galacticus' internal XPath treats a *bare* leading step (not ``.`` / ``..``)
    as ABSOLUTE from the ``<parameters>`` root, whereas standard XPath is
    context-relative.  We preserve that for back-compat while letting richer,
    standard XPath (leading ``/`` or ``//``, predicates, functions, axes -- all
    provided by lxml) pass through unchanged.  An empty path selects ``node``.
    """
    text = path.strip()
    if text == '':
        return [node]
    first_step = text.split('/', 1)[0]
    if first_step in ('.', '..') or text.startswith('/'):
        return node.xpath(text)           # relative to node, or absolute document
    return root.xpath(text)               # legacy bare path: absolute from root


def _first(node, path, root):
    matches = resolve_path(node, path, root)
    return matches[0] if matches else None


# --- XInclude ---------------------------------------------------------------

def _locate_include(href, base_dir):
    """Locate an xi:include target file, or ``None``.

    Resolves STRICTLY relative to the including file (absolute hrefs verbatim),
    exactly as Galacticus does (``source/utility/IO/XML.F90``). We deliberately
    do NOT add the ``GALACTICUS_EXEC_PATH`` leniency that the *validator* uses for
    template files: the resolver's output must match what Galacticus itself would
    produce, so an include Galacticus cannot find must fail here too.
    """
    if os.path.isabs(href):
        return href if os.path.exists(href) else None
    candidate = os.path.normpath(os.path.join(base_dir, href))
    return candidate if os.path.exists(candidate) else None


def _xinclude_fallback_nodes(include_element):
    fallback = include_element.find(_XI_FALLBACK)
    return [copy.deepcopy(child) for child in fallback] if fallback is not None else []


def _xpointer_select(sub_root, xpointer, href):
    """Nodes a ``<xi:include>`` contributes, per its xpointer."""
    if xpointer is None:
        return [copy.deepcopy(sub_root)]            # whole referenced root
    match = _XPOINTER_CHILDREN_RE.match(xpointer.strip())
    if not match:
        raise ResolveError(f"unsupported xi:include xpointer '{xpointer}' "
                           f"(href '{href}')")
    if sub_root.tag != match.group(1):
        raise ResolveError(
            f"xi:include xpointer '{xpointer}' does not match referenced root "
            f"element '{sub_root.tag}' (href '{href}')")
    return [copy.deepcopy(child) for child in sub_root]


def _resolve_include(include_element, base_dir, seen, depth):
    href = include_element.get('href')
    if not href:
        raise ResolveError("xi:include without an href is not supported")
    target = _locate_include(href, base_dir)
    if target is None:
        nodes = _xinclude_fallback_nodes(include_element)
        if nodes or include_element.find(_XI_FALLBACK) is not None:
            return nodes
        raise ResolveError(f"xi:include href '{href}' does not resolve "
                           f"(relative to {base_dir} or the repository root)")
    if target in seen:
        raise ResolveError(f"xi:include cycle detected at '{href}'")
    sub_root = load(target).getroot()
    _expand_xincludes(sub_root, os.path.dirname(target), seen | {target}, depth + 1)
    return _xpointer_select(sub_root, include_element.get('xpointer'), href)


def _expand_xincludes(element, base_dir, seen=frozenset(), depth=0):
    """Recursively replace ``<xi:include>`` descendants of ``element`` in place."""
    if depth > _MAX_INCLUDE_DEPTH:
        raise ResolveError("xi:include nesting too deep")
    for child in list(element):
        if child.tag == _XI_INCLUDE:
            nodes = _resolve_include(child, base_dir, seen, depth)
            index = element.index(child)
            element.remove(child)
            for offset, node in enumerate(nodes):
                element.insert(index + offset, node)
        elif _is_element(child):
            _expand_xincludes(child, base_dir, seen, depth)


# --- change files -----------------------------------------------------------

def _changes_root(tree):
    """The first ``<changes>`` element in a change-file tree (Fortran uses
    ``XML_Get_First_Element_By_Tag_Name(changesDoc,'changes')``)."""
    root = tree.getroot()
    if root.tag == 'changes':
        return root
    found = next(root.iter('changes'), None)
    if found is None:
        raise ResolveError("change file has no <changes> element")
    return found


def _apply_change(root, change):
    """Apply one ``<change>`` element to the parameter tree ``root``."""
    change_type = change.get('type')
    path = change.get('path')
    if change_type is None:
        raise ResolveError("`change` element must have the `type` attribute")
    if path is None:
        raise ResolveError("`change` element must have the `path` attribute")

    exists = path == '' or bool(resolve_path(root, path, root))
    if exists:
        target = root if path == '' else _first(root, path, root)
        if change_type == 'replaceOrAppend':
            change_type = 'replace'
    elif change_type == 'replaceOrAppend':
        # Path absent: append to the parent (strip the last path step).
        parent_path = path.rsplit('/', 1)[0] if '/' in path else ''
        target = root if parent_path == '' else _first(root, parent_path, root)
        if target is None:
            raise ResolveError(f"path '{parent_path}' does not exist")
        change_type = 'append'
    else:
        raise ResolveError(f"path '{path}' does not exist")

    if change_type == 'remove':
        target.getparent().remove(target)

    elif change_type == 'update':
        if target.get('value') is None:
            raise ResolveError(
                "can not update the `value` in a parameter that has no `value`")
        if change.get('value') is None:
            raise ResolveError("`change` element must have the `value` attribute")
        base = target.get('value') if change.get('append') == 'true' else ''
        target.set('value', base + change.get('value'))

    elif change_type == 'append':
        for new_node in _element_children(change):
            target.append(copy.deepcopy(new_node))

    elif change_type in ('insertBefore', 'insertAfter', 'replace'):
        parent = target.getparent()
        index = parent.index(target)
        for offset, new_node in enumerate(_element_children(change)):
            imported = copy.deepcopy(new_node)
            if change_type == 'insertAfter':
                parent.insert(index + 1 + offset, imported)
            else:
                parent.insert(index + offset, imported)
        if change_type == 'replace':
            parent.remove(target)

    elif change_type == 'replaceWith':
        target_path = change.get('target')
        if target_path is None:
            raise ResolveError(
                '`change` element with `type="replaceWith"` must have the '
                '`target` attribute')
        source = _first(root, target_path, root)
        if source is None:
            raise ResolveError(f"target path '{target_path}' does not exist")
        parent = target.getparent()
        index = parent.index(target)
        parent.insert(index, copy.deepcopy(source))
        parent.remove(target)

    elif change_type == 'encapsulate':
        parent = target.getparent()
        index = parent.index(target)
        wrapper = None
        for offset, new_node in enumerate(_element_children(change)):
            imported = copy.deepcopy(new_node)
            parent.insert(index + offset, imported)
            if wrapper is None:
                wrapper = imported
        if wrapper is None:
            raise ResolveError("`encapsulate` change provides no wrapping element")
        parent.remove(target)
        wrapper.append(target)

    else:
        raise ResolveError(f"unknown change type `{change_type}`")


def apply_changes(root, change_files):
    """Apply each change file (in order), each ``<change>`` in document order."""
    for change_file in change_files:
        change_tree = load(change_file)
        _expand_xincludes(change_tree.getroot(), os.path.dirname(change_file) or '.')
        for change in _element_children(_changes_root(change_tree)):
            if change.tag == 'change':
                _apply_change(root, change)


# --- top-level --------------------------------------------------------------

def resolve_tree(tree, base_dir, change_files=()):
    """Resolve an lxml parameter ``tree`` in place: XInclude then changes."""
    root = tree.getroot()
    _expand_xincludes(root, base_dir or '.')
    if change_files:
        apply_changes(root, change_files)
    return tree


def resolve_file(path, change_files=(), output=None):
    """Resolve ``path`` (+ optional change files); optionally write ``output``.

    Returns the resolved lxml ``ElementTree``.
    """
    tree = load(path)
    resolve_tree(tree, os.path.dirname(os.path.abspath(path)), change_files)
    if output is not None:
        with open(output, 'wb') as handle:
            handle.write(to_bytes(tree))
    return tree
