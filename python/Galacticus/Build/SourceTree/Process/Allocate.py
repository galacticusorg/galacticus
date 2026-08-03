"""Processes `allocate` directives: emits an `allocate(var(…))` statement
whose rank is inferred from the variable's declaration and whose shape is
drawn from either a `shape=` or `size=` argument on the directive.

Andrew Benson (ported to Python 2026)
"""



from Galacticus.Build.SourceTree                    import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process            import register_process
from Galacticus.Build.SourceTree.Parse.Declarations import get_declaration


def _declaration_rank(parent, variable):
    """Return the rank (number of dimensions) of `variable` declared under
    `parent`, or 0 if the variable has no `dimension(...)` attribute.
    """
    declaration = get_declaration(parent, variable)
    for attr in declaration.get('attributes') or []:
        if attr.startswith('dimension'):
            # Count commas inside dimension(...): rank = comma count + 1.
            return attr.count(',') + 1
    return 0


def _declaration_has_generic_placeholder(parent, variable):
    """Return True if `variable`'s declaration carries an unresolved generic-
    template placeholder (e.g. `{Type¦rank}`) — in which case the rank can
    not yet be determined and we must defer allocate processing until the
    post-generics rerun on the cloned subtree.

    Galacticus's generic templates use the broken-pipe character `¦`
    (U+00A6) as their separator, which never appears naturally in Fortran
    source — its presence anywhere in the declaration's structured fields
    means a placeholder hasn't been substituted yet.
    """
    declaration = get_declaration(parent, variable)
    fields = list(declaration.get('attributes') or [])
    if declaration.get('type'):
        fields.append(declaration['type'])
    fields.extend(declaration.get('variables') or [])
    return any(isinstance(f, str) and '¦' in f for f in fields)


def _inside_generic_template(node):
    """Return True if `node` is inside a construct whose opener carries an
    unresolved generic-template placeholder (e.g. a subroutine named
    `Foo{Type¦label}`) — in which case declarations inside it may not yet be
    parseable at all (a declaration whose *type* is a placeholder, such as
    `{Type¦intrinsic} :: x`, is not recognized as a declaration until the
    placeholder is substituted).
    """
    p = node.get('parent')
    while p is not None:
        opener = p.get('opener')
        if opener and '¦' in opener:
            return True
        p = p.get('parent')
    return False


def process_allocate(tree, options):
    """Process `allocate` directives in the tree."""
    for node in walk_tree(tree):
        if node.get('type') != 'allocate':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue

        variable = directive['variable']
        if 'rank' in directive:
            rank = int(directive['rank'])
        else:
            # Defer if the declaration is still a generic template — the
            # `{Type¦rank}` placeholder will become a literal `, dimension(:)`
            # (or an empty string) when generics expansion clones the
            # subroutine into one copy per instance.  Leaving the directive
            # un-flagged means generics' clone-then-reparse step preserves
            # its `!![…!!]` markers, and the inner `process_tree` call on
            # each cloned subtree re-invokes us with a fully-resolved
            # declaration.
            try:
                deferred = _declaration_has_generic_placeholder(node['parent'], variable)
            except RuntimeError:
                # A declaration whose type is itself a placeholder (e.g.
                # `{Type¦intrinsic} :: x`) is not recognizable as a declaration
                # at all until generics expansion; defer if we are inside an
                # unexpanded generic template, otherwise the declaration is
                # genuinely missing.
                deferred = _inside_generic_template(node)
                if not deferred:
                    raise
            if deferred:
                continue
            rank = _declaration_rank(node['parent'], variable)

        directive['processed'] = True

        if 'shape' in directive:
            source = directive['shape']
            indices = "" if rank == 0 else "(" + ",".join(
                f"{source}({i})" for i in range(1, rank + 1)) + ")"
        elif 'size' in directive:
            source = directive['size']
            indices = "" if rank == 0 else "(" + ",".join(
                f"size({source},dim={i})" for i in range(1, rank + 1)) + ")"
        else:
            raise RuntimeError("process_allocate: no source given for allocation")

        code  = "! Auto-generated allocation\n"
        code += f"allocate({variable}{indices})\n"
        code += "! End auto-generated allocation\n"

        insert_after_node(node, [{
            'type':       'code',
            'content':    code,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':     node.get('source', 'unknown'),
            'line':       node.get('line', 0),
        }])


register_process('allocate', process_allocate, before=['generics'])
