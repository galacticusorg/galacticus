# Contains a Python module which implements parsing of module uses in the Galacticus preprocessor system.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Parse/ModuleUses.pm

import re
import io
import copy
import os
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

from build.fortran_utils import get_fortran_line

_MODULE_USE_RE = re.compile(
    r'^\s*(!\$)?\s*use\s*(,\s*(intrinsic))?\s*(::)?\s*([a-zA-Z0-9_]+)'
    r'\s*(,\s*only\s*:)?\s*([a-zA-Z0-9_()/=*\-+.,\s]+)?\s*$',
    re.IGNORECASE,
)


def parse_module_uses(tree):
    """Walk the tree parsing Fortran 'use' statements into moduleUse nodes.

    Mirrors Parse_ModuleUses() from perl/Galacticus/Build/SourceTree/Parse/ModuleUses.pm.
    Called by _pass_module_uses() in Galacticus.Build.SourceTree.
    """
    from Galacticus.Build.SourceTree import walk_tree, replace_node, _make_code_node

    nodes_to_replace = []

    for node in walk_tree(tree):
        if node.get('type') != 'code':
            continue

        content = node.get('content', '')
        source  = node.get('source', 'unknown')
        line_no = node.get('line', 0)

        new_nodes       = []
        raw_code_buf    = []
        raw_use_buf     = []
        raw_pp_buf      = []   # preprocessor lines between/before use lines
        module_uses     = {}
        module_order    = []
        pp_stack        = []   # list of {'name': str, 'invert': bool}

        raw_code_line   = line_no
        raw_use_line    = line_no
        current_line    = line_no

        fh = io.StringIO(content)

        def _flush_code():
            nonlocal raw_code_buf, raw_code_line, raw_use_line
            if raw_pp_buf:
                raw_code_buf.extend(raw_pp_buf)
                raw_pp_buf.clear()
            if raw_code_buf:
                new_nodes.append(_make_code_node(
                    ''.join(raw_code_buf), source, raw_code_line))
                raw_code_buf.clear()
                raw_code_line = current_line
                raw_use_line  = current_line

        def _flush_uses():
            nonlocal module_uses, module_order, raw_use_buf, raw_code_line, raw_use_line
            if not module_uses:
                return
            mu_node = {
                'type':        'moduleUse',
                'moduleUse':   module_uses,
                'moduleOrder': module_order,
                'parent':      None,
                'firstChild':  None,
                'sibling':     None,
                'source':      source,
                'line':        raw_use_line,
            }
            child = _make_code_node(''.join(raw_use_buf), source, raw_use_line)
            child['parent'] = mu_node
            mu_node['firstChild'] = child
            new_nodes.append(mu_node)
            module_uses   = {}
            module_order  = []
            raw_use_buf   = []
            raw_code_line = current_line
            raw_use_line  = current_line

        while True:
            raw_line, processed_line, _ = get_fortran_line(fh)
            if not raw_line and not processed_line:
                # EOF — flush whatever is pending
                if raw_pp_buf:
                    raw_code_buf.extend(raw_pp_buf)
                    raw_pp_buf.clear()
                _flush_code()
                _flush_uses()
                break

            n_newlines   = raw_line.count('\n')
            current_line = current_line + n_newlines

            m = _MODULE_USE_RE.match(processed_line)
            is_use = bool(m)

            if is_use:
                # Flush any preceding plain code first.
                _flush_code()
                # Absorb any buffered preprocessor lines into the use block.
                if raw_pp_buf:
                    raw_use_buf.extend(raw_pp_buf)
                    raw_pp_buf.clear()
                raw_use_buf.append(raw_line)

                is_openmp    = bool(m.group(1))
                is_intrinsic = bool(m.group(3))
                mod_name     = m.group(5)
                only_text    = m.group(7)

                if mod_name not in module_uses:
                    if is_intrinsic:
                        module_order.insert(0, mod_name)
                    else:
                        module_order.append(mod_name)

                entry = module_uses.setdefault(mod_name, {
                    'openMP':    is_openmp,
                    'intrinsic': is_intrinsic,
                })
                entry['openMP']    = is_openmp
                entry['intrinsic'] = is_intrinsic

                if pp_stack:
                    entry['conditions'] = copy.deepcopy(pp_stack)

                if only_text and 'all' not in entry:
                    only_text = only_text.strip()
                    only_dict = entry.setdefault('only', {})
                    for sym in re.split(r'\s*,\s*', only_text):
                        sym = re.sub(r'\s', '', sym)
                        if sym:
                            only_dict[sym] = True
                else:
                    entry['all'] = True
                    entry.pop('only', None)

            elif processed_line.startswith('#'):
                # Preprocessor directive.
                if module_uses:
                    # Inside a use block — keep it in the use buffer.
                    is_use = True
                    raw_use_buf.append(raw_line)
                else:
                    raw_pp_buf.append(raw_line)

                # Update the preprocessor stack.
                pm = re.match(r'^#ifdef\s+([A-Z0-9_]+)', processed_line)
                if pm:
                    pp_stack.append({'name': pm.group(1), 'invert': False})
                pm = re.match(r'^#ifndef\s+([A-Z0-9_]+)', processed_line)
                if pm:
                    pp_stack.append({'name': pm.group(1), 'invert': True})
                if re.match(r'^#endif', processed_line) and pp_stack:
                    pp_stack.pop()
                if re.match(r'^#else', processed_line) and pp_stack:
                    top = pp_stack.pop()
                    top['invert'] = not top['invert']
                    pp_stack.append(top)

            else:
                # Plain code line — close any open use block first.
                _flush_uses()
                if raw_pp_buf:
                    raw_code_buf.extend(raw_pp_buf)
                    raw_pp_buf.clear()
                raw_code_buf.append(raw_line)

        # Only replace if the node actually changed.
        single_code = (
            len(new_nodes) == 1
            and new_nodes[0].get('type') == 'code'
            and new_nodes[0].get('content') == content
        )
        if not single_code and new_nodes:
            nodes_to_replace.append((node, new_nodes))

    for old_node, new_node_list in nodes_to_replace:
        replace_node(old_node, new_node_list)


def update_uses(uses_node):
    """Regenerate formatted Fortran 'use' statements from structured node data.

    Mirrors UpdateUses() from perl/Galacticus/Build/SourceTree/Parse/ModuleUses.pm.
    Rewrites uses_node['firstChild']['content'] in place.
    """
    from build.fortran_utils import get_fortran_line

    # Determine indentation from the existing raw content.
    indent = ""
    raw_content = uses_node['firstChild']['content'] if uses_node.get('firstChild') else ""
    fh = io.StringIO(raw_content)
    while indent == "":
        raw_line, processed_line, _ = get_fortran_line(fh)
        if not raw_line and not processed_line:
            break
        m = re.match(r'^(\s*)', raw_line)
        if m:
            indent = m.group(1)

    module_use  = uses_node.get('moduleUse', {})
    module_order = uses_node.get('moduleOrder', list(module_use.keys()))

    openmp_any   = any(v.get('openMP')    for v in module_use.values())
    intrinsic_any = any(v.get('intrinsic') for v in module_use.values())

    name_len_max = max((len(n) for n in module_use), default=0)

    # Compute 4-column max widths for 'only' symbols.
    col_max = [0, 0, 0, 0]
    for mod_name in module_use:
        only = module_use[mod_name].get('only', {})
        if only:
            for i, sym in enumerate(sorted(only.keys())):
                j = i % 4
                if len(sym) > col_max[j]:
                    col_max[j] = len(sym)

    content = ""
    for mod_name in module_order:
        entry = module_use.get(mod_name, {})

        # Emit preprocessor conditions.
        for cond in entry.get('conditions', []):
            directive = "#ifndef" if cond['invert'] else "#ifdef"
            content += directive + " " + cond['name'] + "\n"

        use_line = indent

        if entry.get('openMP'):
            use_line += "!$ "
        elif openmp_any:
            use_line += "   "

        use_line += "use"

        if entry.get('intrinsic'):
            use_line += ", intrinsic"
        elif intrinsic_any:
            use_line += "           "

        use_line += " :: " + mod_name

        only = entry.get('only', {})
        symbol_count = len(only)
        if symbol_count > 0:
            use_line += " " * (name_len_max - len(mod_name)) + ", only : "
            offset_len = len(use_line)
            symbols = sorted(only.keys())
            remaining = symbol_count
            for i, sym in enumerate(symbols):
                j = i % 4
                use_line += sym
                remaining -= 1
                if remaining > 0:
                    use_line += " " * (col_max[j] - len(sym))
                    use_line += ", "
                    if j == 3:
                        use_line += "&\n"
                        if entry.get('openMP'):
                            use_line += "!$"
                        use_line += "          &" + " " * (offset_len - 11)

        use_line += "\n"
        content += use_line

        # Close preprocessor conditions.
        for _ in entry.get('conditions', []):
            content += "#endif\n"

    if uses_node.get('firstChild') is None:
        uses_node['firstChild'] = {
            'type':       'code',
            'content':    content,
            'parent':     uses_node,
            'sibling':    None,
            'firstChild': None,
            'source':     uses_node.get('source', 'unknown'),
            'line':       uses_node.get('line', 0),
        }
    else:
        uses_node['firstChild']['content'] = content


def add_uses(node, module_uses_node):
    """Add module uses to an existing node, creating a moduleUse child if needed.

    Mirrors AddUses() from perl/Galacticus/Build/SourceTree/Parse/ModuleUses.pm.

    Parameters
    ----------
    node : dict
        An AST node whose children will receive the new module uses.
    module_uses_node : dict
        A moduleUse-structured dict with 'moduleUse' and 'moduleOrder' keys.
    """
    from Galacticus.Build.SourceTree import insert_before_node

    # Find an existing moduleUse child.
    uses_node = None
    child = node.get('firstChild')
    while child:
        if child.get('type') == 'moduleUse':
            uses_node = child
        child = child.get('sibling')

    if uses_node is None:
        # Create a new, empty moduleUse node.
        uses_node = {
            'type':        'moduleUse',
            'moduleUse':   {},
            'moduleOrder': [],
            'parent':      None,
            'sibling':     None,
            'source':      'Galacticus.Build.SourceTree.Parse.ModuleUses.add_uses()',
            'line':        1,
        }
        uses_node['firstChild'] = {
            'type':       'code',
            'content':    '',
            'parent':     uses_node,
            'sibling':    None,
            'firstChild': None,
            'source':     uses_node['source'],
            'line':       1,
        }
        first_child = node.get('firstChild')
        if first_child is None:
            uses_node['parent'] = node
            uses_node['sibling'] = None
            node['firstChild'] = uses_node
        else:
            insert_before_node(first_child, [uses_node])

    # Merge module uses.
    new_module_use   = module_uses_node.get('moduleUse', {})
    new_module_order = module_uses_node.get('moduleOrder', sorted(new_module_use.keys()))

    for mod_name in new_module_order:
        new_entry = new_module_use[mod_name]

        if mod_name not in uses_node['moduleOrder']:
            if new_entry.get('intrinsic'):
                uses_node['moduleOrder'].insert(0, mod_name)
            else:
                uses_node['moduleOrder'].append(mod_name)

        existing = uses_node['moduleUse'].setdefault(mod_name, {})
        existing['openMP']    = new_entry.get('openMP', False)
        existing['intrinsic'] = new_entry.get('intrinsic', False)
        if 'conditions' in new_entry:
            existing['conditions'] = copy.deepcopy(new_entry['conditions'])

        if 'all' not in existing:
            if new_entry.get('all'):
                existing['all'] = True
                existing.pop('only', None)
            else:
                only = existing.setdefault('only', {})
                for sym in new_entry.get('only', {}):
                    only[sym] = True

    update_uses(uses_node)
