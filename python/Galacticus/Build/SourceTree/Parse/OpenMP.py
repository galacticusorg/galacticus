# Contains a Python module which implements parsing of OpenMP directives.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Parse/OpenMP.pm

import re
import io
import os
import sys


from build.fortran_utils import get_fortran_line

_OMP_PARALLEL_RE = re.compile(
    r'^\s*!\$omp\s+(end\s+)?parallel\s+(do|workshare)?\s?(.*)',
    re.IGNORECASE,
)

# Options which attach to the composite (iterator) directive in a combined
# `parallel do / workshare`.  Matches @iteratorOptions in OpenMP.pm:23.
_ITERATOR_OPTIONS = ('schedule',)


def _extract_balanced_parens(s):
    """If s starts with optional whitespace then '(', return (content, remainder)
    where content is the balanced parenthesised substring (including the outer
    parens) and remainder is the rest of s with that content stripped.  Otherwise
    return ('', s).

    Replaces the Perl `$RE{balanced}{-parens=>'()'}` pattern from Regexp::Common.
    """
    i = 0
    while i < len(s) and s[i].isspace():
        i += 1
    if i >= len(s) or s[i] != '(':
        return '', s
    depth = 0
    j = i
    while j < len(s):
        c = s[j]
        if c == '(':
            depth += 1
        elif c == ')':
            depth -= 1
            if depth == 0:
                return s[i:j+1], s[j+1:]
        j += 1
    # Unbalanced — treat as no match.
    return '', s


def _parse_options(options_text):
    """Parse the trailing option text of an !$omp parallel[...] directive into
    a list of {'name': str, 'content': str} dicts.  Mirrors the option loop at
    OpenMP.pm:47-64.
    """
    options_text = options_text.lstrip()
    results = []
    while options_text:
        m = re.match(r'^([a-zA-Z_][a-zA-Z_0-9]*)(.*)$', options_text, re.DOTALL)
        if not m:
            break
        name = m.group(1)
        remainder = m.group(2)
        content, _ = _extract_balanced_parens(remainder)
        # Strip the name and any immediately following whitespace.
        options_text = re.sub(r'^' + re.escape(name) + r'\s*', '', options_text)
        if content:
            options_text = options_text[len(content):] if options_text.startswith(content) else options_text
        options_text = options_text.lstrip()
        results.append({'name': name, 'content': content})
    return results


def parse_openmp(tree):
    """Walk the tree extracting `!$omp parallel[...]` directives into openMP nodes.

    Mirrors Parse_OpenMP() from perl/Galacticus/Build/SourceTree/Parse/OpenMP.pm.
    """
    from Galacticus.Build.SourceTree import walk_tree, replace_node, _make_code_node

    nodes_to_replace = []

    for node in walk_tree(tree):
        if node.get('type') != 'code':
            continue

        content = node.get('content', '')
        source  = node.get('source', 'unknown')
        line_no = node.get('line', 0)

        new_nodes_combined = []
        raw_code_buf       = []
        raw_code_line      = line_no
        current_line       = line_no

        fh = io.StringIO(content)

        while True:
            raw_line, processed_line, _ = get_fortran_line(fh)
            if not raw_line and not processed_line:
                break

            n_newlines = raw_line.count('\n')
            line_after = current_line + n_newlines

            m = _OMP_PARALLEL_RE.match(processed_line)
            if m:
                is_closer = m.group(1) is not None
                composite = m.group(2)      # "do" / "workshare" / None
                options   = m.group(3) or ''
                omp_options = _parse_options(options)

                parallel_opts  = [o for o in omp_options
                                  if o['name'] not in _ITERATOR_OPTIONS]
                composite_opts = [o for o in omp_options
                                  if o['name'] in _ITERATOR_OPTIONS]

                if composite is not None:
                    code_parallel  = "!$omp " + ("end " if is_closer else "") + "parallel"
                    for opt in parallel_opts:
                        code_parallel += " " + opt['name'] + opt['content']
                    code_parallel += "\n"
                    code_composite = "!$omp " + ("end " if is_closer else "") + composite
                    for opt in composite_opts:
                        code_composite += " " + opt['name'] + opt['content']
                    code_composite += "\n"
                else:
                    code_parallel  = raw_line
                    code_composite = None

                new_nodes = []

                node_parallel = {
                    'type':       'openMP',
                    'name':       'parallel',
                    'isCloser':   is_closer,
                    'options':    parallel_opts,
                    'parent':     None,
                    'firstChild': None,
                    'sibling':    None,
                    'source':     source,
                    'line':       current_line,
                }
                node_parallel['firstChild'] = {
                    'type':       'code',
                    'content':    code_parallel,
                    'parent':     node_parallel,
                    'sibling':    None,
                    'firstChild': None,
                    'source':     source,
                    'line':       current_line,
                }
                new_nodes.append(node_parallel)

                if composite is not None:
                    node_composite = {
                        'type':       'openMP',
                        'name':       composite,
                        'isCloser':   is_closer,
                        'options':    composite_opts,
                        'parent':     None,
                        'firstChild': None,
                        'sibling':    None,
                        'source':     source,
                        'line':       current_line,
                    }
                    node_composite['firstChild'] = {
                        'type':       'code',
                        'content':    code_composite,
                        'parent':     node_composite,
                        'sibling':    None,
                        'firstChild': None,
                        'source':     source,
                        'line':       current_line,
                    }
                    if is_closer:
                        new_nodes.insert(0, node_composite)
                    else:
                        new_nodes.append(node_composite)

                if raw_code_buf:
                    code_node = _make_code_node(
                        ''.join(raw_code_buf), source, raw_code_line)
                    new_nodes.insert(0, code_node)
                    raw_code_buf = []
                    raw_code_line = line_after

                new_nodes_combined.extend(new_nodes)
            else:
                raw_code_buf.append(raw_line)

            current_line = line_after

        if raw_code_buf:
            new_nodes_combined.append(_make_code_node(
                ''.join(raw_code_buf), source, raw_code_line))
            raw_code_buf = []

        single_code = (
            len(new_nodes_combined) == 1
            and new_nodes_combined[0].get('type') == 'code'
        )
        if not single_code and new_nodes_combined:
            nodes_to_replace.append((node, new_nodes_combined))

    for old_node, new_node_list in nodes_to_replace:
        replace_node(old_node, new_node_list)


def update(node):
    """Regenerate the `!$omp …` line in node['firstChild']['content'] from
    node's name / isCloser / options.

    Mirrors OpenMP.pm:174-179.
    """
    content = "!$omp " + ("end " if node.get('isCloser') else "") + node['name']
    for opt in node.get('options', []):
        content += " " + opt['name'] + opt.get('content', '')
    content += "\n"
    if node.get('firstChild') is None:
        node['firstChild'] = {
            'type':       'code',
            'content':    content,
            'parent':     node,
            'sibling':    None,
            'firstChild': None,
            'source':     node.get('source', 'unknown'),
            'line':       node.get('line', 0),
        }
    else:
        node['firstChild']['content'] = content


def copyin(node, variable_names):
    """Add `copyin(...)` variables to an openMP node, deduplicating existing entries.

    Mirrors OpenMP.pm:181-208.
    """
    options = node.setdefault('options', [])
    copyin_opt = next((o for o in options if o['name'] == 'copyin'), None)
    if copyin_opt is None:
        copyin_opt = {'name': 'copyin', 'content': '()'}
        options.append(copyin_opt)

    inner = copyin_opt.get('content', '()').strip()
    m = re.match(r'^\((.*)\)$', inner, re.DOTALL)
    inner = m.group(1) if m else ''
    names = [n.strip() for n in re.split(r'\s*,\s*', inner) if n.strip()]

    for name in variable_names:
        if name not in names:
            names.append(name)

    copyin_opt['content'] = '(' + ','.join(names) + ')'
    update(node)
