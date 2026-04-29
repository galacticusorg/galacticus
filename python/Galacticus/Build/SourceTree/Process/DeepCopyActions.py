# Processes `deepCopyActions` directives: for each directive naming a base
# class, walks every derived (non-abstract) type in the same tree and
# emits an `XDeepCopyActions(self)` subroutine whose `select type` branches
# call the requested `setTo` / `methodCall` actions declared per level of
# the inheritance chain.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/DeepCopyActions.pm

import os
import re
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from List.ExtraUtils                     import as_array
from Galacticus.Build.SourceTree         import (
    walk_tree, parse_code, children, insert_post_contains,
)
from Galacticus.Build.SourceTree.Process import register_process, process_tree


_TYPE_OPENER_RE = re.compile(
    r'^\s*type\s*'
    r'((?:,\s*(?:abstract|public|private|extends\s*\([a-zA-Z0-9_]+\))\s*)*)'
    r'(?:::)?\s*([a-zA-Z0-9_]+)\s*$',
    re.IGNORECASE,
)
_EXTENDS_ATTR_RE = re.compile(r'extends\(([a-zA-Z0-9_]+)\)', re.IGNORECASE)


def _parse_type_opener(opener):
    """Return (type_name, extends_name_or_None, abstract_bool) parsed from a
    `type [, attr, ... ::] name` opener line, or None if unparseable.

    Mirrors the regex branch at DeepCopyActions.pm:41-55.  Openers containing
    `{…}` generic placeholders are ignored (Perl's explicit carve-out).
    """
    if '{' in opener:
        return None
    m = _TYPE_OPENER_RE.match(opener)
    if not m:
        raise RuntimeError(
            "process_deep_copy_actions: unable to parse type definition opener")
    attrs_text = m.group(1).strip().strip(',').strip()
    type_name  = m.group(2)
    extends    = None
    abstract   = False
    for attr in re.split(r'\s*,\s*', attrs_text) if attrs_text else []:
        if attr == 'abstract':
            abstract = True
        em = _EXTENDS_ATTR_RE.match(attr)
        if em:
            extends = em.group(1)
    return type_name, extends, abstract


def _emit_deep_copy_action(class_name, classes, directive):
    """Render the generated `XDeepCopyActions` subroutine text.

    Mirrors DeepCopyActions.pm:76-135.
    """
    out  = f"subroutine {class_name}DeepCopyActions(self)\n"
    out += " !!{\n"
    out += " Perform actions needed for deep copy of this object.\n"
    out += " !!}\n"
    out += " implicit none\n"
    out += f" class({class_name}), intent(inout) :: self\n"
    out += " select type (self)\n"

    # Emit one `type is (Sub)` branch per non-abstract class descended from
    # the target class.  The Perl code walks classes in sorted-key order.
    for sub_name in sorted(classes.keys()):
        info = classes[sub_name]
        if info['abstract']:
            continue
        # Is `sub_name` derived (transitively) from `class_name`?
        matches = False
        cursor  = sub_name
        while cursor is not None:
            if cursor == class_name:
                matches = True
                break
            cursor = classes.get(cursor, {}).get('extends')
        if not matches:
            continue

        out += f" type is ({sub_name})\n"

        # For each ancestor level (including sub_name itself), pull setTo /
        # methodCall actions declared at that level, provided the target
        # class has any declarations in its body (matches Perl's
        # `$classNode->{'type'} eq "declaration"` guard).
        method_calls = []
        cursor = sub_name
        while cursor is not None:
            class_node = classes[cursor]['node']
            scope_directive = directive.get(cursor) if isinstance(directive.get(cursor), dict) else None

            child = class_node.get('firstChild')
            while child is not None:
                if child.get('type') == 'declaration' and scope_directive is not None:
                    for set_to in as_array(scope_directive.get('setTo')):
                        variables = re.split(
                            r'\s*,\s*', set_to.get('variables', ''))
                        state = set_to.get('state', '')
                        for v in variables:
                            out += f" self%{v}={state}\n"
                    for method_call in as_array(scope_directive.get('methodCall')):
                        args = method_call.get('arguments', '')
                        method_calls.append(
                            f" call self%{method_call.get('method', '')}({args})"
                        )
                child = child.get('sibling')
            cursor = classes[cursor].get('extends')

        if method_calls:
            out += "\n".join(method_calls) + "\n"

    out += "  end select\n"
    out += "  return\n"
    out += f"end subroutine {class_name}DeepCopyActions\n"
    return out


def _emit_type_binding(class_name):
    """Render the `procedure :: deepCopyActions => XDeepCopyActions` binding
    plus its `<methods>` metadata directive.
    """
    return (
        "    !![\n"
        "    <methods>\n"
        "     <method method=\"deepCopyActions\" description=\"Perform actions "
        "needed for deep copy of this object.\"/>\n"
        "    </methods>\n"
        "    !!]\n"
        f"    procedure :: deepCopyActions => {class_name}DeepCopyActions\n"
    )


def _insert_parsed(parent, source_text, run_process_tree=False):
    sub_tree = parse_code(source_text, name='DeepCopyActions')
    if run_process_tree:
        process_tree(sub_tree)
    kids = children(sub_tree)
    for k in kids:
        k['parent'] = None
    insert_post_contains(parent, kids)


def process_deep_copy_actions(tree, options):
    """Mirrors Process_DeepCopyActions() from DeepCopyActions.pm."""
    directive_nodes = []
    classes         = {}

    # Pass 1: collect directives and type nodes in a single walk.
    for node in walk_tree(tree):
        ntype = node.get('type')
        if ntype == 'deepCopyActions':
            directive = node.setdefault('directive', {})
            directive['processed'] = True
            directive_nodes.append(node)
            continue
        if ntype == 'type':
            parsed = _parse_type_opener(node.get('opener', ''))
            if parsed is None:
                continue
            type_name, extends, abstract = parsed
            classes[type_name] = {
                'node':     node,
                'extends':  extends,
                'abstract': abstract,
            }

    for directive_node in directive_nodes:
        directive = directive_node['directive']
        class_name = directive.get('class')
        if class_name not in classes:
            raise RuntimeError(
                f"process_deep_copy_actions: class '{class_name}' not found")

        action_text  = _emit_deep_copy_action(class_name, classes, directive)
        binding_text = _emit_type_binding(class_name)

        # Type-binding into the base class's type node.
        _insert_parsed(classes[class_name]['node'], binding_text)

        # Generated subroutine into the directive's parent module/container,
        # after being run through process_tree so any inner directives get
        # expanded.
        _insert_parsed(
            directive_node['parent'], action_text, run_process_tree=True)


register_process('deepCopyActions', process_deep_copy_actions)
