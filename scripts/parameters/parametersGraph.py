#!/usr/bin/env python3
"""Map the structure of a Galacticus model: the objects its parameter file
selects, the class defaults it leaves unstated, and what each builds.

Usage:
    parametersGraph.py <parameterFile> [--changes <file> ...] [--catalog <file>]
                       [--json <file>] [--html <file>] [--fragment]
                       [--data-url <url>]

``--json <file>`` writes the graph itself: nodes (objects), edges (which object
builds which) and groups (physics subsystems). ``--html <file>`` writes a
self-contained, interactive map of it, to open in a browser. ``--fragment``
writes that page without its ``<html>`` wrapper, for embedding in another page.
With ``--changes``, what the change files added, moved, altered or removed is
marked, against the model without them.
``--data-url <url>`` has the page load the graph from that URL (the file written
by ``--json``, served beside it) instead of carrying the graph within itself.
With neither output named, a summary is printed.

The parameter catalog is rebuilt from the source tree (a few seconds) unless
``--catalog <file>`` names a `parameters.catalog.json` to reuse.

Andrew Benson (2026).
"""

import argparse
import collections
import html
import json
import os
import sys

# Make the in-tree Python packages importable when run directly.
_REPOSITORY = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
sys.path.insert(0, os.path.join(_REPOSITORY, 'python'))

from Galacticus.Parameters import resolve
from Galacticus.Parameters.graph import build_graph, compare_graphs

_TEMPLATE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'parametersGraphTemplate.html')
_DATA_MARKER = '/*MODEL_GRAPH*/null'
_URL_MARKER = '/*MODEL_URL*/null'
_BODY_MARKER = '<!--BODY-->'
_TITLE_MARKER = '<title>Galacticus Model Map</title>'


def _load_catalog(path):
    if path is not None:
        with open(path) as handle:
            return json.load(handle)
    # Set, not setdefault: the directive parser finds its schemas through this.
    os.environ['GALACTICUS_EXEC_PATH'] = _REPOSITORY
    from Galacticus.Parameters.catalog import build_catalog
    return build_catalog(os.path.join(_REPOSITORY, 'source'), log=lambda message: None)


def _render(graph, fragment, data_url=None):
    with open(_TEMPLATE) as handle:
        template = handle.read()
    if any(marker not in template for marker in (_DATA_MARKER, _URL_MARKER, _BODY_MARKER, _TITLE_MARKER)):
        raise RuntimeError(f'{_TEMPLATE} lacks one of its markers')
    # Title the page for its model, so that the maps of several can be told apart.
    model = ' + '.join(os.path.splitext(part)[0]
                       for part in [graph['model']['name'] or ''] + graph['model'].get('changes', []) if part)
    if model:
        template = template.replace(_TITLE_MARKER, f'<title>{html.escape(model)} Model Map</title>')
    if data_url is not None:
        page = template.replace(_URL_MARKER, json.dumps(data_url))
    else:
        # "</" would end the script element which carries the data.
        data = json.dumps(graph, separators=(',', ':')).replace('</', '<\\/')
        page = template.replace(_DATA_MARKER, data)
    if fragment:
        return page.replace(_BODY_MARKER, '')
    head, body = page.split(_BODY_MARKER, 1)
    return ('<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
            f'{head}</head>\n<body>\n{body}</body>\n</html>\n')


def _summary(graph):
    nodes, edges = graph['nodes'], graph['edges']
    titles = {group['id']: group['title'] for group in graph['groups']}
    lines = [
        f"{graph['model']['name']}",
        f"  objects           : {len(nodes)}",
        f"    set in the file : {sum(1 for node in nodes if node['origin'] == 'file')}",
        f"    class defaults  : {sum(1 for node in nodes if node['origin'] == 'default')}",
        f"    not seen built  : {sum(1 for node in nodes if not node['reached'])}",
    ] + ([
        f"    changed         : {sum(1 for node in nodes if node.get('changed'))}",
        f"    removed         : {len(graph['model'].get('removed', []))}",
    ] if graph['model'].get('changes') else []) + [
        f"  links             : {len(edges)}",
    ]
    counts = collections.Counter(node['group'] for node in nodes)
    for group in graph['groups']:
        if counts[group['id']]:
            lines.append(f"  {titles[group['id']]:<36}: {counts[group['id']]}")
    return '\n'.join(lines)


def main(argv):
    parser = argparse.ArgumentParser(
        prog='parametersGraph.py',
        description='Map the structure of the model described by a Galacticus parameter file.')
    parser.add_argument('parameterFile')
    parser.add_argument('--changes', action='append', default=[], metavar='FILE',
                        help='a change file to apply to the parameter file (repeatable)')
    parser.add_argument('--catalog', metavar='FILE',
                        help='a parameters.catalog.json to reuse, rather than rebuilding it')
    parser.add_argument('--json', metavar='FILE', help='write the graph as JSON')
    parser.add_argument('--html', metavar='FILE', help='write an interactive map of the graph')
    parser.add_argument('--fragment', action='store_true',
                        help='write the map without its <html> wrapper, for embedding')
    parser.add_argument('--data-url', metavar='URL',
                        help='have the map load the graph from this URL, rather than carry it')
    arguments = parser.parse_args(argv[1:])

    try:
        tree = resolve.resolve_file(arguments.parameterFile, change_files=arguments.changes)
    except (OSError, resolve.ResolveError) as error:
        print(f'parametersGraph.py: {error}', file=sys.stderr)
        return 1
    catalog = _load_catalog(arguments.catalog)
    source_root = os.path.join(_REPOSITORY, 'source')
    name = os.path.basename(arguments.parameterFile)
    graph = build_graph(tree.getroot(), catalog, source_root=source_root, name=name)
    if arguments.changes:
        # Mark what the change files altered, against the model without them.
        base = build_graph(resolve.resolve_file(arguments.parameterFile).getroot(), catalog,
                           source_root=source_root, name=name)
        compare_graphs(graph, base)
        graph['model']['changes'] = [os.path.basename(path) for path in arguments.changes]

    if arguments.json:
        with open(arguments.json, 'w') as handle:
            json.dump(graph, handle, indent=1)
            handle.write('\n')
    if arguments.html:
        with open(arguments.html, 'w') as handle:
            handle.write(_render(graph, arguments.fragment, arguments.data_url))
    if not (arguments.json or arguments.html):
        print(_summary(graph))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
