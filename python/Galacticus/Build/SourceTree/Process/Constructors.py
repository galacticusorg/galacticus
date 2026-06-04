"""Processes `constructorAssign` directives: for each named variable, assigns
it to `self` (or pointer-assigns when prefixed with `*` / `*/`), allocates
allocatable-with-deferred-shape members, and bumps the reference count of
functionClass pointer members so the construction chain is correct.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/SourceTree/Process/Constructors.pm
"""

import os
import re
import xml.etree.ElementTree as ET


from XML.Utils                                      import xml_to_dict
from Galacticus.Build.StateStorables                import (
    function_class_names    as _shared_function_class_names,
    function_class_instances as _shared_function_class_instances,
)
from Galacticus.Build.SourceTree                    import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process            import register_process
from Galacticus.Build.SourceTree.Parse.Declarations import get_declaration


_STATE_STORABLES = None


def _state_storables():
    """Lazily load `$BUILDPATH/stateStorables.xml`, caching at module scope.

    Mirrors the `unless ($stateStorables)` idiom in Constructors.pm:47-49.
    """
    global _STATE_STORABLES
    if _STATE_STORABLES is not None:
        return _STATE_STORABLES
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        raise RuntimeError(
            "process_constructors: BUILDPATH is not set")
    path = os.path.join(build_path, 'stateStorables.xml')
    if not os.path.exists(path):
        raise RuntimeError(
            f"process_constructors: stateStorables.xml not found at {path}")
    _STATE_STORABLES = xml_to_dict(ET.parse(path).getroot())
    return _STATE_STORABLES


def _function_class_type_set(state_storables):
    """Return the set of lowercased type names considered functionClasses.

    Mirrors `keys(%{$stateStorables->{'functionClasses'}})` plus the
    `functionClassInstances` list in Constructors.pm:106.
    """
    return (
        {n.lower() for n in _shared_function_class_names(state_storables)}
        | {n.lower() for n in _shared_function_class_instances(state_storables)}
    )


def _return_value_label(parent):
    """Return the function's result name.

    Perl uses the `result(...)` identifier from the function opener if
    present, otherwise the function name itself (Constructors.pm:53-58).
    """
    opener = parent.get('opener') or ''
    m = re.search(r'result\s*\(\s*([a-zA-Z0-9_]+)\s*\)\s*$', opener)
    if m:
        return m.group(1)
    return parent.get('name', '')


_VARIABLE_TOKEN_RE = re.compile(r'^(\*?)(/?)([a-zA-Z0-9_]+)')


def process_constructors(tree, options):
    """Mirrors Process_Constructors() from Constructors.pm."""
    function_class_types = None  # lazy-loaded on first hit

    for node in walk_tree(tree):
        if node.get('type') != 'constructorAssign':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue

        parent = node.get('parent') or {}
        if parent.get('type') not in ('function', 'moduleProcedure'):
            # Walk to the root for the error message, matching Perl.
            root = node
            while root.get('parent') is not None:
                root = root['parent']
            raise RuntimeError(
                f"process_constructors: parent node must be a function "
                f"(current parent is '{parent.get('type')}') in file "
                f"'{root.get('name')}'"
            )

        if function_class_types is None:
            function_class_types = _function_class_type_set(_state_storables())

        allocate = directive.get('allocate', 'yes')
        return_value = _return_value_label(parent)
        directive['processed'] = True

        lines = "  ! Auto-generated constructor assignment\n"

        variables_text = (directive.get('variables') or '').strip()
        for token in [v for v in re.split(r'\s*,\s*', variables_text) if v]:
            m = _VARIABLE_TOKEN_RE.match(token)
            if not m:
                raise RuntimeError("process_constructors: syntax error")
            is_pointer    = m.group(1) == '*'
            do_count      = m.group(2) != '/'
            argument_name = m.group(3)
            assigner      = ' => ' if is_pointer else '='

            default_match = re.search(r'=\s*(.+)', token)
            has_default = default_match is not None
            default     = default_match.group(1) if has_default else None

            # Perl always calls GetDeclaration here when the parent is a
            # function; for moduleProcedure parents no declaration is fetched
            # (the generated code is intentionally missing the optional /
            # allocatable handling in that case — matches Perl).
            declaration = None
            if parent.get('type') != 'moduleProcedure':
                declaration = get_declaration(parent, argument_name)

            attributes = (declaration or {}).get('attributes') or []
            optional = (f"if (present({argument_name})) "
                        if any(a == 'optional' for a in attributes) else '')

            # Auto-allocate for deferred-shape allocatable members when
            # requested.  Matches Constructors.pm:82-87.
            if allocate == 'yes' and declaration is not None:
                dim_matches = [
                    re.search(r'dimension\s*\(([:,]+)\)', a)
                    for a in attributes
                ]
                dim_inner = ''.join(dm.group(1) for dm in dim_matches if dm)
                rank = dim_inner.count(':')
                if dim_inner and rank > 0:
                    sizes = ','.join(
                        f"size({argument_name},dim={i})" for i in range(1, rank + 1)
                    )
                    lines += (
                        f"   {optional}allocate({return_value}%{argument_name}"
                        f"({sizes}))\n"
                    )

            if optional == '':
                lines += f"   {return_value}%{argument_name}{assigner}{argument_name}\n"
            elif has_default:
                lines += (
                    f"   {optional} then\n"
                    f"{return_value}%{argument_name}{assigner}{argument_name}\n"
                    f"else\n"
                    f"{return_value}%{argument_name}{assigner}{default}\n"
                    f"end if\n"
                )
            else:
                lines += (
                    f"   {optional}{return_value}%{argument_name}"
                    f"{assigner}{argument_name}\n"
                )

            # Reference-count bump for functionClass pointer members.
            if (declaration is not None
                    and declaration.get('intrinsic') in ('type', 'class')
                    and is_pointer
                    and do_count):
                lowered_type = re.sub(r'\s', '', (declaration.get('type') or '').lower())
                if lowered_type in function_class_types:
                    if optional:
                        lines += f"   {optional} then\n"
                    if is_pointer:
                        lines += (
                            f"   if (associated({return_value}%{argument_name})) "
                        )
                    lines += (
                        f" call {return_value}%{argument_name}%referenceCountIncrement()\n"
                    )
                    if optional:
                        lines += "   end if\n"
                elif declaration.get('type') == '*':
                    if optional:
                        lines += f"   {optional} then\n"
                    if is_pointer:
                        lines += (
                            f"   if (associated({return_value}%{argument_name})) then\n"
                        )
                    lines += (
                        f"select type(s__ => {return_value}%{argument_name})\n"
                        "class is (functionClass)\n"
                        " call s__%referenceCountIncrement()\n"
                        "end select\n"
                    )
                    if is_pointer:
                        lines += "   end if\n"
                    if optional:
                        lines += "   end if\n"

        lines += "  ! End auto-generated constructor assignment.\n\n"

        insert_after_node(node, [{
            'type':       'code',
            'content':    lines,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':
                'Galacticus.Build.SourceTree.Process.Constructors.process_constructors()',
            'line':       1,
        }])


register_process('constructors', process_constructors)
