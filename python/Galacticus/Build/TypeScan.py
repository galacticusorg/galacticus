"""Shared helpers for scanning derived-type definitions in parsed source trees.

Used by the deepCopyActions and stateStorables catalog scripts, which both
need the map of derived types (name, parent, abstractness) defined in a file
together with the AST nodes carrying a given directive. Historically each
script carried its own copy (a legacy of the Perl-to-Python port).
"""

import re

from Galacticus.Build.SourceTree import walk_tree

__all__ = ['TYPE_OPENER_RE', 'collect_types_and_directives', 'inherits_from']


# Matches a derived-type definition opener, capturing (abstract, extends
# parent, type name).
TYPE_OPENER_RE = re.compile(
    r'^\s*type\s*'
    r'(?:,\s*(?:(abstract)|public|private|extends\s*\(([a-zA-Z0-9_]+)\))\s*)*'
    r'(?:::)?\s*([a-zA-Z0-9_]+)\s*$',
    re.IGNORECASE,
)


def collect_types_and_directives(tree, directive_type):
    """Walk `tree`, returning ``(type-dict, directive-node-list)``.

    type-dict: ``{type_name: {'extends': parent-or-None, 'abstract': bool}}``
    directive-node-list: every AST node whose ``type`` attribute equals
    `directive_type`.

    Raises ``RuntimeError`` if a type-definition opener cannot be parsed
    (which would silently corrupt the inheritance map).
    """
    classes    = {}
    directives = []
    for node in walk_tree(tree):
        if node.get('type') == directive_type:
            directives.append(node)
        if node.get('type') == 'type':
            m = TYPE_OPENER_RE.match(node.get('opener') or '')
            if not m:
                raise RuntimeError(
                    "TypeScan: unable to parse type definition opener: "
                    f"{node.get('opener')!r}"
                )
            abstract, extends, name = m.group(1), m.group(2), m.group(3)
            classes[name] = {
                'extends':  extends,
                'abstract': abstract is not None,
            }
    return classes, directives


def inherits_from(classes, class_name, base_class):
    """True if `class_name` derives (directly or transitively) from
    `base_class` within the `classes` map built by
    :func:`collect_types_and_directives`.
    """
    current = class_name
    while current is not None:
        if current == base_class:
            return True
        current = classes.get(current, {}).get('extends')
    return False
