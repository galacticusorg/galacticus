"""Processes `inputParameter` directives: emits the
`source%value('name', var, defaultValue=…, writeOutput=…)` call that binds
a Fortran variable to a named run-time parameter, and ensures the enclosing
subprogram imports `Input_Parameters`. For each bound declared by a
`<minimum>` or `<maximum>` element, it also emits a call to
`inputParameterBoundCheck` (module `Input_Parameters_Bounds`), which reports
a fatal error if the value read violates the bound.

Andrew Benson (ported to Python 2026)
"""

import re

from Galacticus.Build.SourceTree                  import walk_tree, insert_after_node
from Galacticus.Build.SourceTree.Process          import register_process
from Galacticus.Build.SourceTree.Parse.ModuleUses import add_uses


def _bound(raw):
    """Parse a `<minimum>`/`<maximum>` directive value into `(text, inclusive)`.

    `raw` is the `xml_to_dict` value: a plain string (`<minimum>0.0</minimum>`), or a dict when the
    element carries attributes (`<minimum inclusive="false">0.0</minimum>`).
    """
    if isinstance(raw, dict):
        text      = str(raw.get('content', '')).strip()
        inclusive = raw.get('inclusive', 'true') != 'false'
    else:
        text      = str(raw).strip()
        inclusive = True
    return text, inclusive


def _double_literal(text, parameter):
    """Return `text`, a numeric bound, as a Fortran double precision literal."""
    try:
        float(text.lower().replace('d', 'e'))
    except ValueError:
        raise RuntimeError(
            f"process_input_parameters: bound '{text}' of parameter '{parameter}' is not a number")
    if re.search(r'[dD]', text):
        return text
    if re.search(r'[eE]', text):
        return re.sub(r'[eE]', 'd', text)
    return text + 'd0'


def process_input_parameters(tree, options):
    """Process `inputParameter` directives in the tree."""
    for node in walk_tree(tree):
        if node.get('type') != 'inputParameter':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        directive['processed'] = True

        code  = "  ! Auto-generated input parameter\n"
        if 'name' in directive:
            source        = directive['source']
            parameter     = directive['name']
            # Add delimiters unless the name is actually a function call.
            parameter_arg = parameter if '(' in parameter else f"'{parameter}'"
            variable      = directive.get('variable', directive['name'])
            code += f"  call {source}%value({parameter_arg},{variable}"
            if 'defaultValue' in directive:
                code += f",defaultValue={directive['defaultValue']}"
            if 'writeOutput' in directive:
                write_flag = ('.false.' if directive['writeOutput'] == 'no'
                              else '.true.')
                code += f",writeOutput={write_flag}"
            code += ")\n"
            # Enforce any declared bounds.
            bounds = [
                (element, *_bound(directive[element]))
                for element in ('minimum', 'maximum') if element in directive
            ]
            for element, text, inclusive in bounds:
                code += (
                    f"  call inputParameterBoundCheck({source},{parameter_arg},{variable},"
                    f"{_double_literal(text, parameter)},'{text}',"
                    f"isMinimum={'.true.' if element == 'minimum' else '.false.'},"
                    f"isInclusive={'.true.' if inclusive else '.false.'})\n"
                )
            if bounds:
                add_uses(node['parent'], {
                    'moduleUse': {
                        'Input_Parameters_Bounds': {
                            'intrinsic': False,
                            'only':      {'inputParameterBoundCheck': True},
                        },
                    },
                    'moduleOrder': ['Input_Parameters_Bounds'],
                })
        code += "  ! End auto-generated input parameter\n\n"

        insert_after_node(node, [{
            'type':       'code',
            'content':    code,
            'parent':     None,
            'firstChild': None,
            'sibling':    None,
            'source':
                'Galacticus.Build.SourceTree.Process.InputParameter.process_input_parameters()',
            'line':       1,
        }])

        add_uses(node['parent'], {
            'moduleUse': {
                'Input_Parameters': {
                    'intrinsic': False,
                    'all':       True,
                },
            },
            'moduleOrder': ['Input_Parameters'],
        })


register_process('inputParameters', process_input_parameters)
