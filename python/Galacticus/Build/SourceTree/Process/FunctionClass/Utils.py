"""Small utility helpers shared by the functionClass pipeline.

Andrew Benson (ported to Python 2026)
"""

import re


from Galacticus.Build.SourceTree import walk_tree


def _walk_forward(node):
    """Yield `node` and every subsequent node in depth-first order.

    Unlike `walk_tree(subtree)`, this walk does not stop at the end of the
    starting node's subtree — it continues through siblings and
    up-then-siblings all the way through the rest of the whole tree.  We
    materialise this by starting at the top of the tree, collecting every
    node into a list, and iterating from the given node's position onward.
    """
    root = node
    while root.get('parent') is not None:
        root = root['parent']
    sequence = list(walk_tree(root))
    try:
        idx = sequence.index(node)
    except ValueError:
        idx = 0
    for i in range(idx, len(sequence)):
        yield sequence[i]


def class_dependencies(class_node, directive_name):
    """Walk forward from `class_node` until we reach the first `type` node
    whose name is `<directive_name><suffix>`, capturing:
      - a `class` record with the directive attributes, its `extends` parent,
        and its own type name;
      - the ordered list of type-dependency names the class pulls in (its
        parent type plus any `class(…)` / `type(…)` members whose type is
        itself an instance of the same functionClass family).
    """
    class_record = {}
    dependencies = []
    pattern = re.compile(
        r'^\s*type\s*(,\s*abstract\s*|,\s*public\s*|,\s*private\s*'
        r'|,\s*extends\s*\(([a-zA-Z0-9_]+)\)\s*)*'
        r'(?:::)?\s*' + re.escape(directive_name) + r'([a-z0-9_]+)\s*$',
        re.IGNORECASE,
    )

    for node in _walk_forward(class_node):
        # Directive node — copy every attribute.
        if node.get('type') == directive_name:
            class_record['node'] = node
            for key, value in sorted((node.setdefault('directive', {})).items()):
                class_record[key] = value
            continue

        if node.get('type') == 'type':
            m = pattern.match(node.get('opener') or '')
            if not m or m.group(2) is None:
                continue
            class_record['extends'] = m.group(2)
            class_record['type']    = directive_name + m.group(3)
            dependencies.append(m.group(2))
            # Pull in cross-references: `class(Xfoo)` / `type(Xfoo)` members
            # whose type is another member of the same functionClass family.
            child = node.get('firstChild')
            while child is not None:
                if child.get('type') == 'declaration':
                    for declaration in child.get('declarations') or []:
                        dep_type = (declaration.get('type') or '').strip()
                        if (declaration.get('intrinsic') in ('class', 'type')
                                and dep_type.startswith(directive_name)
                                and not re.search(r'Class\s*$', dep_type)):
                            dependencies.append(dep_type)
                child = child.get('sibling')
            break

    # Deduplicate and sort the dependency list.
    return class_record, sorted(set(dependencies))


def latex_breakable(text):
    """Insert a LaTeX soft hyphen between every lower→upper case transition."""
    return re.sub(r'([a-z])([A-Z])', r'\1\\-\2', text or '')


def trimlc(text):
    """Lowercase, then strip leading+trailing whitespace."""
    return (text or '').lower().strip()


def striplc(text):
    """Lowercase, then strip ALL whitespace."""
    return re.sub(r'\s', '', (text or '').lower())


def lctrim(text):
    """Strip trailing whitespace, then lowercase."""
    return re.sub(r'\s*$', '', (text or '')).lower()


def strip_variable_name(text):
    """Return the leading [a-zA-Z0-9_]+ identifier, dropping everything after
    it (array indices, `=` initializers, etc.).
    """
    m = re.match(r'^([a-zA-Z0-9_]+)', text or '')
    return m.group(1) if m else ''


def declaration_rank(declaration):
    """Return the rank of `declaration` (number of commas in its
    `dimension(…)` attribute, +1).  Scalars → 0.
    """
    attributes = declaration.get('attributes') or []
    if not any(re.match(r'^dimension\s*\(', a) for a in attributes):
        return 0
    dim_text = ''.join(
        m.group(1)
        for m in (re.match(r'^dimension\s*\(([a-zA-Z0-9_,:\s]+)\)', a)
                  for a in attributes)
        if m
    )
    return dim_text.count(',') + 1
