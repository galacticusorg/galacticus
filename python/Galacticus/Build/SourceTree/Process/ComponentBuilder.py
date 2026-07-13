"""Processes the `componentBuilder` directive: builds the node-component class
hierarchy at preprocess time and grafts it in place of the directive.

Andrew Benson (2026)

This replaces the retired `component` include-directive path (the old
`scripts/build/buildCode.py` + `Galacticus.Build.Hooks` mechanism).  A single
`<componentBuilder/>` directive lives in `source/objects/nodes/_class.F90`;
when encountered it:

  1. Reads `directiveLocations.xml` for every source file carrying a
     `<component>` directive.
  2. Extracts, XSD-validates, and parses each `<component>` body through the
     existing `Galacticus.Build.Components` generator (unchanged), then runs
     its phased `generate_output` to synthesize the Fortran source.
  3. Parses the synthesized source into a sub-tree and runs the full
     `process_tree` pipeline over it — expanding any directives the generator
     embedded in its own output (e.g. `<allocate>`, `<methods>`,
     `<optionalArgument>`) exactly as `buildCode.py` did before.
  4. Grafts the fragment's children in place of the directive node.

The generated content already carries its own module `contains` marker, so the
hand-written procedures that follow the directive in `_class.F90` remain
post-`contains` — reproducing the layout the compiler-`include` produced.

The heavy `Galacticus.Build.Components` package (~50 modules) is imported
*lazily*, inside the handler, only on first encounter of a `componentBuilder`
node.  `Process/all.py` (and `scripts/build/preprocess.py`) import this module
for every one of the ~2000 preprocess jobs in a clean build, so importing the
generator eagerly here would load those 50 modules needlessly for every job
that has no `componentBuilder` directive.

Mirrors the cross-file gathering pattern already used by
`Process/FunctionClass`, `Process/EventHooks`, etc.
"""

import os
import xml.etree.ElementTree as ET


from List.ExtraUtils                     import as_array
from XML.Utils                           import xml_to_dict
from Galacticus.Build.Directives         import extract_directives
from Galacticus.Build.SourceTree         import (
    walk_tree, parse_code, children, replace_node,
)
from Galacticus.Build.SourceTree.Process import register_process, process_tree


# `XMLin($xml, ForceArray => ["data","property","binding"])` in the Perl-era
# driver (buildCode.py:91): the component generator expects these tags as lists
# even when only a single instance is present.
_FORCE_ARRAY = {'data', 'property', 'binding'}

_DIRECTIVE_LOCATIONS = None


def _load_directive_locations():
    """Return the parsed `directiveLocations.xml`, cached across nodes."""
    global _DIRECTIVE_LOCATIONS
    if _DIRECTIVE_LOCATIONS is not None:
        return _DIRECTIVE_LOCATIONS
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        raise RuntimeError("process_component_builder: BUILDPATH is not set")
    path = os.path.join(build_path, 'directiveLocations.xml')
    _DIRECTIVE_LOCATIONS = xml_to_dict(ET.parse(path).getroot())
    return _DIRECTIVE_LOCATIONS


def _build_component_content():
    """Run the component generator over every `<component>` directive and
    return the synthesized Fortran source string (`build['content']`).

    Mirrors the driver loop of the retired `buildCode.py`: validate + parse
    each directive, then run the whole-build `generate_output`.  The generator
    package is imported here (lazily) rather than at module scope.
    """
    # Lazy, first-use import of the (heavy) component generator package.
    import Galacticus.Build.Components as Components

    directive_locations = _load_directive_locations()
    files = list(as_array((directive_locations.get('component') or {}).get('file')))

    build = {}
    for file_name in files:
        for document in extract_directives(
            file_name, 'component',
            force_array=_FORCE_ARRAY, include_raw_xml=True,
        ):
            raw_xml = document.pop('rawXML')
            Components.validate(raw_xml, file_name)
            build['currentDocument'] = document
            Components.parse_directive(build)

    Components.generate_output(build)
    return build['content']


def process_component_builder(tree, options):
    """Mirrors the `component` build, now driven from the SourceTree pipeline.

    Every preprocess invocation walks the tree here, but the (heavy) generator
    machinery only runs when a `componentBuilder` node is actually present —
    i.e. only when preprocessing `source/objects/nodes/_class.F90`.
    """
    for node in list(walk_tree(tree)):
        if node.get('type') != 'componentBuilder':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        directive['processed'] = True

        content = _build_component_content()

        # Parse the synthesized source and run the full process pipeline over
        # it, exactly as buildCode.py did (buildCode.py:194-199), so directives
        # the generator embedded in its own output are expanded and marked
        # processed.  The parsed fragment's directive nodes come back
        # `processed`, so the outer walk (and later hooks) skip them.
        build_path = os.environ['BUILDPATH']
        fragment_name = os.path.join(build_path, 'objects.nodes.components.p.Inc')
        sub_tree = parse_code(content, name=fragment_name)
        process_tree(sub_tree)

        kids = children(sub_tree)
        for kid in kids:
            kid['parent'] = None
        replace_node(node, kids)


register_process('componentBuilder', process_component_builder)
