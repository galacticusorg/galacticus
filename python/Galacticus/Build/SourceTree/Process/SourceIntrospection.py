"""Expands `{introspection:location[:compact]}` placeholders in Fortran source
into an expression that reports the enclosing file/module/subroutine chain
and line number at runtime.

Andrew Benson (ported to Python 2026)
"""

import re


from Galacticus.Build.SourceTree         import walk_tree
from Galacticus.Build.SourceTree.Process import register_process


def instrument(code):
    """Tag every `{introspection:location[:compact]}` placeholder with the
    source line number on which it appears.

    Intended to run
    once at parse time, before the tree is built: the line numbers baked in
    here are consumed later by `process_source_introspection` to emit the
    actual location expression.
    """
    out = []
    for i, line in enumerate(code.splitlines(keepends=True), start=1):
        out.append(
            re.sub(
                r'\{(introspection:location(?::compact)?)\}',
                lambda m: '{' + m.group(1) + ':' + str(i) + '}',
                line,
            )
        )
    return ''.join(out)


def location(node, line_number, compact=False):
    """Build the Fortran expression that names every structural ancestor of
    `node` plus `line_number`.

    The emitted text is a `//`-concatenated expression that is safe to
    splice into any string context in Fortran.
    """
    if compact:
        expr = "''"
    else:
        expr = "char(10)//' Occurred at:'"

    branch = node
    while branch is not None:
        btype = branch.get('type')
        if btype in ('file', 'module', 'function', 'subroutine'):
            bname = branch.get('name', '')
            if compact:
                expr += "//'" + btype + "(" + bname + ")'"
            else:
                pad = ' ' * (10 - len(btype))
                expr += "//char(10)//'   " + pad + btype + ":" + bname + "'"
        if 'directive' in branch and not compact:
            expr += "//char(10)//'    directive:" + (btype or '') + "'"
        branch = branch.get('parent')

    if compact:
        expr += "//':" + str(line_number) + "'"
    else:
        expr += "//'   [line " + str(line_number) + "]'"
    return expr


def process_source_introspection(tree, options):
    """Replace tagged `{introspection:location:NNN}` placeholders with the
    result of `location(node, NNN)`.
    """
    pattern = re.compile(r'\{introspection:location(:compact)?:(\d+)\}')

    for node in walk_tree(tree):
        if node.get('type') != 'code':
            continue
        content = node.get('content')
        if content is None:
            raise RuntimeError(
                "process_source_introspection: code node has no content")
        if 'introspection:location' not in content:
            continue

        new_lines = []
        for line in content.splitlines(keepends=True):
            def _substitute(m, _line_node=node):
                compact = bool(m.group(1))
                line_no = int(m.group(2))
                return location(_line_node, line_no, compact=compact)
            new_lines.append(pattern.sub(_substitute, line))
        node['content'] = ''.join(new_lines)


register_process('sourceIntrospection', process_source_introspection)
