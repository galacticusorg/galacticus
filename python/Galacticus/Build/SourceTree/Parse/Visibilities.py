"""Contains a Python module which implements parsing of visibilities in the Galacticus preprocessor system.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/SourceTree/Parse/Visibilities.pm
"""

import re
import io


from Galacticus.Build.FortranUtils import get_fortran_line

_VISIBILITY_RE = re.compile(r'^\s*(public|private)\s*(::)?\s*(.*?)\s*$', re.IGNORECASE)


def parse_visibilities(tree):
    """Walk the tree parsing Fortran visibility statements into visibility nodes.

    Mirrors Parse_Visibilities() from perl/Galacticus/Build/SourceTree/Parse/Visibilities.pm.
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
        raw_vis_buf     = []
        visibilities    = {}

        raw_code_line   = line_no
        raw_vis_line    = line_no
        current_line    = line_no

        fh = io.StringIO(content)

        def _flush_code():
            nonlocal raw_code_buf, raw_code_line, raw_vis_line
            if raw_code_buf:
                new_nodes.append(_make_code_node(
                    ''.join(raw_code_buf), source, raw_code_line))
                raw_code_buf.clear()
                raw_code_line = current_line
                raw_vis_line  = current_line

        def _flush_vis():
            nonlocal visibilities, raw_vis_buf, raw_code_line, raw_vis_line
            if not visibilities and not raw_vis_buf:
                return
            vis_node = {
                'type':       'visibility',
                'visibility': visibilities,
                'parent':     None,
                'firstChild': None,
                'sibling':    None,
                'source':     source,
                'line':       raw_vis_line,
            }
            child = _make_code_node(''.join(raw_vis_buf), source, raw_vis_line)
            child['parent'] = vis_node
            vis_node['firstChild'] = child
            new_nodes.append(vis_node)
            visibilities  = {}
            raw_vis_buf   = []
            raw_code_line = current_line
            raw_vis_line  = current_line

        while True:
            raw_line, processed_line, _ = get_fortran_line(fh)
            if not raw_line and not processed_line:
                _flush_code()
                _flush_vis()
                break

            n_newlines = raw_line.count('\n')
            line_after = current_line + n_newlines

            m = _VISIBILITY_RE.match(processed_line)
            if m:
                _flush_code()
                raw_vis_buf.append(raw_line)
                level = m.group(1).lower()
                names = m.group(3)
                # Mirror Perl: each line clears then repopulates the bucket for
                # this visibility level.
                visibilities[level] = {}
                if names:
                    for sym in re.split(r'\s*,\s*', names):
                        sym = sym.strip()
                        if sym:
                            visibilities[level][sym] = True
            else:
                _flush_vis()
                raw_code_buf.append(raw_line)

            current_line = line_after

        single_code = (
            len(new_nodes) == 1
            and new_nodes[0].get('type') == 'code'
            and new_nodes[0].get('content') == content
        )
        if not single_code and new_nodes:
            nodes_to_replace.append((node, new_nodes))

    for old_node, new_node_list in nodes_to_replace:
        replace_node(old_node, new_node_list)


def update_visibilities(visibilities_node):
    """Regenerate Fortran visibility statements from a visibility node.

    Mirrors UpdateVisibilities() from Parse/Visibilities.pm.  Rewrites
    visibilities_node['firstChild']['content'] in place.
    """
    content = "! Galacticus::Build::SourceTree::Parse::Visibilities(): updated\n"
    vis = visibilities_node.get('visibility', {}) or {}
    for level in ("public", "private"):
        names = vis.get(level) or {}
        if not names:
            continue
        content += level + " :: " + ", ".join(sorted(names.keys())) + "\n"
    if visibilities_node.get('firstChild') is None:
        visibilities_node['firstChild'] = {
            'type':       'code',
            'content':    content,
            'parent':     visibilities_node,
            'sibling':    None,
            'firstChild': None,
            'source':     visibilities_node.get('source', 'unknown'),
            'line':       visibilities_node.get('line', 0),
        }
    else:
        visibilities_node['firstChild']['content'] = content
