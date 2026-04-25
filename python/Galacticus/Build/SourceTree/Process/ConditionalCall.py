# Processes `conditionalCall` directives: emits the full 2^N if-then ladder
# over the set of argument conditions, calling the underlying routine once
# per boolean combination with only the arguments whose condition is true.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Process/ConditionalCall.pm

import os
import re
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from List.ExtraUtils                                import as_array
from Galacticus.Build.SourceTree                    import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process            import register_process
from Galacticus.Build.SourceTree.Parse.Declarations import (
    add_declarations, declaration_exists,
)


def _normalize_arguments(raw):
    """Return a dict-of-dicts keyed by argument name.

    Our xml_to_dict (unlike Perl's XML::Simple with its default KeyAttr)
    does not group same-tag children by their `name` attribute.  Accept
    whichever shape it returned:

      - list of dicts (multiple `<argument …/>`)          → key by name
      - dict with no nested dicts (single `<argument …/>`) → wrap as {name: dict}
      - dict already keyed by name (dict-of-dicts)        → pass through
    """
    if isinstance(raw, list):
        return {a['name']: a for a in raw}
    if isinstance(raw, dict):
        nested = any(isinstance(v, dict) for v in raw.values())
        if nested:
            return raw
        if 'name' not in raw:
            raise RuntimeError(
                "conditionalCall: malformed `argument` directive (no `name`)")
        return {raw['name']: raw}
    raise RuntimeError(
        "conditionalCall: unexpected type for `argument`: "
        + type(raw).__name__)


def _argument_condition(argument):
    """Return the Fortran expression representing this argument's condition.

    Mirrors ConditionalCall.pm:42-49.  Prefers `parameterPresent` (which
    becomes a call to the parameter's `isPresent` method) over a literal
    `condition` attribute.
    """
    if 'parameterPresent' in argument:
        param_name = argument.get('parameterName', argument['name'])
        return f"{argument['parameterPresent']}%isPresent('{param_name}')"
    if 'condition' in argument:
        return argument['condition']
    raise RuntimeError("conditionalCall: no condition specified")


def process_conditional_call(tree, options):
    """Mirrors Process_ConditionalCall() from ConditionalCall.pm."""
    for node in walk_tree(tree):
        if node.get('type') != 'conditionalCall':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        directive['processed'] = True

        if 'argument' not in directive:
            raise RuntimeError(
                "conditionalCall: at least one argument must be given")
        arguments = _normalize_arguments(directive['argument'])

        # Attach each argument's resolved condition text and build the unique
        # condition list.  The ordering here (sorted by condition string) is
        # the ordering Perl uses for id assignment below.
        for arg in arguments.values():
            arg['conditionActual'] = _argument_condition(arg)

        unique_conditions = sorted({a['conditionActual'] for a in arguments.values()})
        condition_id = {c: i + 1 for i, c in enumerate(unique_conditions)}

        # Emit logical declarations for any `condition<i>__` that don't yet exist.
        needed_vars = [
            f'condition{i}__' for i in range(1, len(unique_conditions) + 1)
            if not declaration_exists(node['parent'], f'condition{i}__')
        ]
        variables_to_add = []
        if needed_vars:
            variables_to_add.append({
                'intrinsic':     'logical',
                'type':          None,
                'openMP':        False,
                'attributes':    [],
                'variables':     list(needed_vars),
                'variableNames': list(needed_vars),
            })

        code = "! Auto-generated conditional function call.\n"
        for cond in unique_conditions:
            code += f"condition{condition_id[cond]}__={cond}\n"

        n = len(unique_conditions)
        for state in range(2 ** n):
            bits = [(state >> (n - 1 - j)) & 1 for j in range(n)]   # MSB-first
            conjuncts = [
                ('' if b else '.not.') + f'condition{j + 1}__'
                for j, b in enumerate(bits)
            ]
            code += "if (" + " .and. ".join(conjuncts) + ") then\n"
            for call_template in as_array(directive.get('call')):
                code += _render_call(call_template, arguments, bits, condition_id)
            code += "end if\n"

        insert_after_node(node, [{
            'type':       'code',
            'content':    code,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':     node.get('source', 'unknown'),
            'line':       node.get('line', 0),
        }])
        if variables_to_add:
            add_declarations(node['parent'], variables_to_add)


def _render_call(call_template, arguments, bits, condition_id):
    """Render one `<call>` template for one boolean state `bits`.

    Mirrors ConditionalCall.pm:74-100: lines without `{conditions}` are copied
    verbatim; the `{conditions}` line becomes the call with only those
    optional arguments whose condition is true in this state.  Raises if the
    template contains no `{conditions}` marker at all.
    """
    out = ""
    found = False
    for line in call_template.splitlines(keepends=True):
        m = re.match(r'^(.*)\{conditions\}(.*)$', line, re.DOTALL)
        if not m:
            out += line
            continue
        found = True
        pre  = m.group(1)
        post = m.group(2)
        prefix = "" if pre.rstrip('\n').endswith('(') else ","
        optional_args = [
            f"{arg['name']}={arg['value']}"
            for arg in arguments.values()
            if bits[condition_id[arg['conditionActual']] - 1] == 1
        ]
        out += pre
        if optional_args:
            out += prefix + ",".join(optional_args)
        out += post
        if not post.endswith('\n'):
            out += '\n'
    if not found:
        raise RuntimeError(
            "conditionalCall: syntax error in call element (missing `{conditions}`)")
    # The XML parser strips trailing whitespace from the `<call>` text body,
    # so a template that ends with `)` (no trailing newline) leaves `out`
    # without one too — the caller then appends `end if\n` directly to the
    # closing parenthesis: `…)end if`.  Make sure rendered call text always
    # ends with a newline.
    if not out.endswith('\n'):
        out += '\n'
    return out


register_process('conditionalCall', process_conditional_call)
