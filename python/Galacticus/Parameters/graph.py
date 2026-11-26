"""Build the object graph of a Galacticus model from its parameter file.

Andrew Benson (2026)

A model is a set of ``functionClass`` objects, each of which builds (through
``<objectBuilder>`` directives) the other objects it depends on. Only some of
those objects are named in the parameter file: the rest are class defaults.
This module combines a resolved parameter tree with the parameter catalog
(`Galacticus.Parameters.catalog`) to recover the whole graph:

* a **node** for each object: one per parameter element which selects an
  implementation, and one per class default built at the root;
* an **edge** from each object to each object it builds, resolved as
  ``objectBuilder`` does (the nearest element of that name in the object's own
  parameters or an enclosing element, else the class default);
* a **group** for each node, a physics subsystem assigned from the path of the
  source file which implements it.

Node components (``componentDisk`` and so on) are not ``functionClass`` objects,
so are not in the catalog. Where a source tree is given, they are added from
the ``<component>``, ``<objectBuilder>`` and ``<inputParameter>`` directives of
`source/objects/nodes/components`.

A few objects are built by a hand-written call rather than an
``<objectBuilder>`` (through ``Functions_Global``, where a direct dependence
would make a cycle of modules). Those declare what they build in the
``<allowedParameters>`` of their class directive, which is read here too.

An object is resolved as ``objectBuilder`` resolves it:

* one named for its class is looked for in the parameters of the object which
  builds it, then in each enclosing element, and failing that is the class
  default, built once at the root (so seeing what is at the root, whatever
  built it);
* one given a parameter name of its own has no class default: where no element
  of that name is found it is not built, unless its builder carries an inline
  ``<default>``. Those of node components are followed; those of the few
  ``functionClass`` implementations which have one are not in the catalog, so
  are not.

An implementation which extends another builds what that one builds only where
its parameter constructor calls that of the other. Most do not: they build
their own objects and call the internal constructor. So the objects of an
extended implementation are followed only where that call is found in the
source.

The graph is static: it is what the directives say will be built. Objects
constructed in hand-written code which declares nothing are not seen.
"""

import os
import re

from lxml import etree

from Galacticus.Parameters.query import resolved_parameters

__all__ = ['GROUPS', 'build_graph', 'compare_graphs', 'group_of']

# Physics subsystems, in the order a model is read: from the cosmology, through
# the growth of structure and the evolution of galaxies, to what is output.
GROUPS = [
    ('cosmology',          'Cosmology',
     'The background universe: its contents and expansion history.'),
    ('structure',          'Structure formation',
     'Linear theory and halo statistics: power spectra, growth, mass functions.'),
    ('trees',              'Merger trees',
     'How the merger history of each halo is built or read.'),
    ('engine',             'Evolution engine',
     'The task, and the machinery which evolves each tree through time.'),
    ('components',         'Galaxy components',
     'What each node carries: halo, hot gas, disk, spheroid, black hole.'),
    ('halos',              'Dark matter halos',
     'Halo sizes, density profiles, concentrations and spins.'),
    ('satellites',         'Subhalos and mergers',
     'Orbits, tidal effects and the merging of satellites with their hosts.'),
    ('cgm',                'Circumgalactic gas and cooling',
     'Accretion onto halos, the hot atmosphere, and how it cools and inflows.'),
    ('galaxy',             'Galactic structure',
     'Mass distributions, sizes and dynamical instabilities of galaxies.'),
    ('starFormation',      'Star formation and feedback',
     'Turning gas into stars, and the outflows which stars drive.'),
    ('blackHoles',         'Black holes',
     'Seeding, accretion, winds and mergers of supermassive black holes.'),
    ('stellarPopulations', 'Stellar populations and dust',
     'Stellar evolution, spectra, luminosities and dust.'),
    ('radiation',          'Intergalactic medium and radiation',
     'The state of intergalactic gas and the radiation fields which act on it.'),
    ('output',             'Outputs and analyses',
     'When and what is written, and the analyses computed from it.'),
    ('numerics',           'Numerics and statistics',
     'Random numbers, distributions and other shared numerical tools.'),
]

# The group of an implementation is the first of these patterns to match the
# path of its source file, relative to `source/`. The node operators apply the
# physics of every subsystem to a node, so each is filed with its physics.
_GROUP_RULES = [(re.compile(pattern), group) for pattern, group in (
    (r'^nodes/operators/physics/black_holes',                                         'blackHoles'),
    (r'^nodes/operators/physics/(CGM|circumgalactic_medium|cooling|ram_pressure)',    'cgm'),
    (r'^nodes/operators/physics/(star_formation|stellar_feedback|nuclear_star_clusters|impulsive_outflow)',
                                                                                      'starFormation'),
    (r'^nodes/operators/physics/(satellite|tidal|position|subsubhalo)',               'satellites'),
    (r'^nodes/operators/satellite',                                                   'satellites'),
    (r'^nodes/operators/physics/(dark_matter|halo_|Bertschinger)',                    'halos'),
    (r'^nodes/operators/physics/(bar_instability|inclination|misalignment|disk)',     'galaxy'),
    (r'^nodes/operators/(analyses|output|metadata)',                                  'output'),
    (r'^nodes/operators',                                                             'engine'),
    (r'^nodes/property_extractor',                                                    'output'),
    (r'^merger_trees/(evolve|evolver|initialize|node_evolver|node_merger)',           'engine'),
    (r'^merger_trees/outputter',                                                      'output'),
    (r'^merger_trees',                                                                'trees'),
    (r'^cosmology',                                                                   'cosmology'),
    (r'^(structure_formation|halo_model|dark_matter/)',                               'structure'),
    (r'^(dark_matter_halos|dark_matter_profiles)',                                    'halos'),
    (r'^(satellites|tidal_stripping)',                                                'satellites'),
    (r'^(hot_halo|cooling|accretion/|circumgalactic_medium|chemical|atomic|ram_pressure_stripping)',
                                                                                      'cgm'),
    (r'^(galactic_structure|galactic_dynamics|mass_distributions)',                   'galaxy'),
    (r'^(star_formation|stellar_feedback|nuclear_star_clusters|stellar_astrophysics/feedback)',
                                                                                      'starFormation'),
    (r'^(black_holes|accretion_disks)',                                               'blackHoles'),
    (r'^(stellar_populations|stellar_astrophysics|dust|hiiRegions)',                  'stellarPopulations'),
    (r'^(intergalactic_medium|radiation|radiative_transfer)',                         'radiation'),
    (r'^(tasks|universe|meta|computational_domains)',                                 'engine'),
    (r'^(output|galactic/filters|geometry|statistics/mass_function)',                 'output'),
)]

_COMPONENTS_DIRECTORY = os.path.join('objects', 'nodes', 'components')


def group_of(source_file):
    """Return the group id for an implementation in `source_file` (a path
    relative to `source/`)."""
    path = (source_file or '').replace(os.sep, '/')
    if path.startswith('./'):
        path = path[2:]
    if path.startswith('source/'):
        path = path[len('source/'):]
    for pattern, group in _GROUP_RULES:
        if pattern.search(path):
            return group
    return 'numerics'


def _upper_first(text):
    return text[:1].upper() + text[1:]


def _plain(text):
    """Reduce a reStructuredText description to one line of plain text."""
    if not text:
        return None
    text = re.sub(r':cite:\w*:`([^`]*)`', r'[\1]', text)
    text = re.sub(r':[\w-]+:`([^`<]*?)\s*(<[^`]*>)?`', r'\1', text)
    text = text.replace('``', '')
    return ' '.join(text.split())


def _plain_default(value):
    """Render a Fortran default as a user would write it in a parameter file."""
    if value is None:
        return None
    value = str(value).strip()
    match = re.match(r"^var_str\(\s*'(.*)'\s*\)$", value)
    if match:
        return match.group(1)
    if value.lower() in ('.true.', '.false.'):
        return value.strip('.').lower()
    return re.sub(r'(?<=[\d.])[dD](?=[+-]?\d)', 'e', value)


def _elements(element):
    """The element children of `element` (comments and the like excluded)."""
    return [child for child in element if isinstance(child.tag, str)]


def _scan_components(source_root):
    """Return ``{(class, name): component}`` from the node component sources,
    each with its ``isDefault`` flag, the component it ``extends``, the objects
    it builds, and the parameters it reads."""
    from Galacticus.Build.Directives import extract_directives

    directory = os.path.join(source_root, _COMPONENTS_DIRECTORY)
    per_file = {}
    for path, directories, names in os.walk(directory):
        directories.sort()
        for name in sorted(names):
            if not name.endswith('.F90'):
                continue
            file_name = os.path.join(path, name)
            directives = extract_directives(
                file_name, ('component', 'objectBuilder', 'inputParameter'),
                set_root_element_type=True)
            if directives:
                per_file[file_name] = directives
    components = {}
    owner_of_file = {}
    owners_in_directory = {}
    for file_name, directives in per_file.items():
        for directive in directives:
            if directive.get('rootElementType') != 'component':
                continue
            key = (directive.get('class'), directive.get('name'))
            extends = directive.get('extends') or {}
            components[key] = {
                'class':       key[0],
                'name':        key[1],
                'isDefault':   str(directive.get('isDefault')).lower() == 'true',
                'extends':     (extends.get('class'), extends.get('name')) if extends else None,
                'description': _plain(directive.get('description')),
                'sourceFile':  os.path.relpath(file_name, source_root).replace(os.sep, '/'),
                'objects':     [],
                'parameters':  [],
            }
            owner_of_file[file_name] = key
            owners_in_directory.setdefault(os.path.dirname(file_name), []).append(key)
    for file_name, directives in per_file.items():
        owner = owner_of_file.get(file_name)
        if owner is None:
            # A file which declares no component belongs to the one component
            # declared beside it, if there is exactly one.
            owners = owners_in_directory.get(os.path.dirname(file_name), [])
            owner = owners[0] if len(owners) == 1 else None
        if owner is None:
            continue
        component = components[owner]
        for directive in directives:
            kind = directive.get('rootElementType')
            if kind == 'objectBuilder' and directive.get('class'):
                entry = {
                    'class':         directive['class'],
                    'parameterName': directive.get('parameterName') or directive['class'],
                    'own':           directive.get('source') != 'parameters',
                    'default':       directive.get('default'),
                }
                if entry['parameterName'] not in [obj['parameterName'] for obj in component['objects']]:
                    component['objects'].append(entry)
            elif kind == 'inputParameter' and directive.get('name'):
                names = [parameter['name'] for parameter in component['parameters']]
                if directive['name'] not in names:
                    component['parameters'].append({
                        'name':        directive['name'],
                        'default':     _plain_default(directive.get('defaultValue')),
                        'description': _plain(directive.get('description')),
                    })
    return components


def _scan_declared_objects(source_root, catalog):
    """Return ``{type: [object, ...]}`` for the objects which implementations
    build by a hand-written call, and so declare only in the
    ``<allowedParameters>`` of their class directive."""
    from Galacticus.Build.Directives import extract_directives

    function_classes = catalog['functionClasses']
    declared = {}
    by_file = {}
    for type_name, implementation in catalog['implementations'].items():
        if implementation.get('sourceFile'):
            by_file.setdefault(implementation['sourceFile'], []).append(type_name)
    for source_file, type_names in sorted(by_file.items()):
        path = os.path.join(source_root, source_file)
        try:
            with open(path, errors='replace') as handle:
                if '<allowedParameters>' not in handle.read():
                    continue
        except OSError:
            continue
        class_names = {catalog['implementations'][type_name]['functionClass'] for type_name in type_names}
        for directive in extract_directives(path, tuple(sorted(class_names))):
            if directive.get('name') not in type_names:
                continue
            entries = (directive.get('allowedParameters') or {}).get('object') or []
            for entry in entries if isinstance(entries, list) else [entries]:
                if entry.get('name') in function_classes and entry.get('source') == 'parameters':
                    declared.setdefault(directive['name'], []).append({
                        'class': entry['name'], 'parameterName': entry['name'],
                        'repeatable': False, 'source': 'parameters', 'sourceElement': None})
    return declared


def _scan_delegations(source_root, catalog):
    """Return the types whose parameter constructor hands its parameters to
    that of the implementation it extends, and so builds all that it builds."""
    implementations = catalog['implementations']
    delegating = set()
    for type_name, implementation in implementations.items():
        parent = implementation.get('parent')
        if parent not in implementations or not implementation.get('sourceFile'):
            continue
        try:
            with open(os.path.join(source_root, implementation['sourceFile']), errors='replace') as handle:
                text = handle.read()
        except OSError:
            continue
        if re.search(r'=\s*' + re.escape(parent) + r'\s*\(\s*parameters\s*\)', text, re.IGNORECASE):
            delegating.add(type_name)
    return delegating


def _default_element(parameter_name, default):
    """An element holding the inline ``<default>`` of an object builder. It
    is attached to nothing, as the parameters built from such a default are."""
    def fill(element, content):
        for key, value in (content or {}).items():
            if isinstance(value, dict):
                fill(etree.SubElement(element, key), value)
            elif isinstance(value, list):
                for item in value:
                    fill(etree.SubElement(element, key), item if isinstance(item, dict) else {'value': item})
            elif key == 'value':
                element.set('value', str(value))
    element = etree.Element(parameter_name)
    fill(element, (default or {}).get(parameter_name))
    return element


def compare_graphs(graph, base):
    """Mark, in `graph`, what differs from `base` (the same model before its
    change files were applied): each node which `base` lacks, or holds with
    other parameters, gains ``changed``; the nodes of `base` which `graph`
    lacks are listed as ``graph['model']['removed']``. An object which a
    change moves (by wrapping it in another, say) is both."""
    def place(node):
        return (node['class'], node['label'], node['path'])

    def state(node):
        return place(node) + (tuple((p['name'], p['value']) for p in node['parameters'] if p['set']),)

    before = {state(node) for node in base['nodes']}
    after = {place(node) for node in graph['nodes']}
    for node in graph['nodes']:
        if state(node) not in before:
            node['changed'] = True
    graph['model']['removed'] = [
        {'class': node['class'], 'label': node['label'], 'path': node['path']}
        for node in base['nodes'] if place(node) not in after]
    return graph


def build_graph(root, catalog, source_root=None, name=None):
    """Return the object graph of the model whose resolved ``<parameters>``
    element is `root` (see `Galacticus.Parameters.resolve`).

    Parameters
    ----------
    root : lxml.etree._Element
        The resolved ``<parameters>`` root element.
    catalog : dict
        A parameter catalog, as built by `Galacticus.Parameters.catalog`.
    source_root : str, optional
        The `source/` directory. Where given, node components are included, as
        are objects which are built by a hand-written call, and those which an
        implementation builds through the one it extends.
    name : str, optional
        A name for the model, recorded in the result.

    Returns
    -------
    dict
        ``{'model', 'groups', 'nodes', 'edges'}``. Each node has an ``id``, the
        ``class`` and ``label`` (implementation) of its object, its ``group``,
        its ``origin`` (``'file'`` or ``'default'``), whether it is ``reached``
        from the task or a component, and its ``parameters``. Each edge runs
        ``from`` an object ``to`` one it builds, ``as`` a named parameter. A
        node whose element is not named for its class records that ``element``.
    """
    function_classes = catalog['functionClasses']
    implementations = catalog['implementations']
    type_by_label = {}
    for type_name, implementation in implementations.items():
        for label in (implementation['label'], implementation.get('alias')):
            if label:
                type_by_label[(implementation['functionClass'], label)] = type_name

    tree = root.getroottree()
    parent_of = {child: parent for parent in root.iter() for child in parent}
    identified = {(element.tag, element.get('id')): element
                  for element in root.iter()
                  if isinstance(element.tag, str) and element.get('id') is not None}

    declared = _scan_declared_objects(source_root, catalog) if source_root is not None else {}
    delegating = _scan_delegations(source_root, catalog) if source_root is not None else set()

    nodes = {}
    by_id = {}
    order = []
    edges = []
    edge_keys = set()

    def lineage(type_name):
        """An implementation, then each it extends for so long as the
        parameter constructor of one calls that of the next."""
        chain, current = [], type_name
        while current in implementations and current not in chain:
            chain.append(current)
            if current not in delegating:
                break
            current = implementations[current].get('parent')
        return chain

    def schema(type_name):
        """The objects an implementation builds, its own and inherited."""
        objects = []
        for current in lineage(type_name):
            objects.extend(implementations[current].get('objects', []))
            objects.extend(declared.get(current, []))
        return objects

    def dereference(element):
        if element.get('idRef') is not None:
            return identified.get((element.tag, element.get('idRef')))
        return element

    def locate(scope, parameter_name, repeatable=False):
        """The elements named `parameter_name` in `scope` or, failing that, the
        nearest enclosing element: every one of them where the object is
        `repeatable`, else the first. Empty where there is none, which selects
        the class default."""
        current = scope
        while current is not None:
            found = [child for child in _elements(current) if child.tag == parameter_name]
            if found:
                found = found if repeatable else found[:1]
                return [element for element in map(dereference, found) if element is not None]
            current = parent_of.get(current)
        return []

    def path_of(element):
        """A readable path for an element: its ancestors' selections."""
        parts, current = [], element
        while current is not None and current is not root:
            value = current.get('value')
            parts.append(f'{current.tag}[{value}]' if value and len(value) < 48 else current.tag)
            current = parent_of.get(current)
        return ' / '.join(reversed(parts))

    def scalar_values(element, object_names):
        """The scalar parameters set in `element`: its children which are not
        themselves objects."""
        values = {}
        if element is None:
            return values
        for child in _elements(element):
            if child.tag in object_names or child.tag in values:
                continue
            value = child.get('value')
            if value is None and child.text and child.text.strip() and not len(child):
                value = child.text.strip()
            if value is not None:
                values[child.tag] = ' '.join(value.split())
        return values

    def parameters_of(type_name, element):
        objects = schema(type_name)
        object_names = {obj['parameterName'] for obj in objects} | {
            obj['sourceElement'] for obj in objects if obj.get('sourceElement')}
        values = scalar_values(element, object_names | set(function_classes))
        parameters = []
        read_by = lineage(type_name)
        for parameter in resolved_parameters(catalog, type_name):
            if parameter.get('sourceElement') is not None or not parameter.get('name'):
                continue
            if parameter.get('inheritedFrom') and parameter['inheritedFrom'] not in read_by:
                continue
            value = values.pop(parameter['name'], None)
            parameters.append({
                'name':        parameter['name'],
                'value':       value if value is not None else _plain_default(parameter.get('default')),
                'set':         value is not None,
                'description': _plain(parameter.get('description')),
            })
        # Anything else set in the element is read in hand-written code.
        for parameter_name, value in values.items():
            parameters.append({'name': parameter_name, 'value': value, 'set': True, 'description': None})
        return parameters

    def mark_reached(node_id):
        """Record that a node, and so all that it builds, is reached."""
        node = by_id[node_id]
        if node['reached']:
            return
        node['reached'] = True
        for edge in list(edges):
            if edge['from'] == node_id:
                mark_reached(edge['to'])

    def add_edge(source, target, parameter_name):
        if source is None or target is None or (source, target, parameter_name) in edge_keys:
            return
        edge_keys.add((source, target, parameter_name))
        edges.append({'from': source, 'to': target, 'as': parameter_name})

    def visit(class_name, element, reached, key=None):
        """Add the object of `class_name` selected by `element`, or the class
        default where `element` is None, and all that it builds. A `key` marks
        `element` as an inline default, which is in no parameter file. Returns
        the node id."""
        if class_name not in function_classes:
            return None
        label = element.get('value') if element is not None else None
        if label is None:
            label = function_classes[class_name].get('default')
        if label is None:
            return None
        inline = key is not None
        if key is None:
            key = ('default', class_name) if element is None else ('file', tree.getpath(element))
        node = nodes.get(key)
        if node is not None:
            if reached:
                mark_reached(node['id'])
            return node['id']
        type_name = type_by_label.get((class_name, label))
        implementation = implementations.get(type_name, {})
        source_file = implementation.get('sourceFile') or function_classes[class_name].get('sourceFile')
        node = {
            'id':         f'n{len(order)}',
            'kind':       'object',
            'class':      class_name,
            'className':  function_classes[class_name].get('descriptiveName') or class_name,
            'label':      label,
            'type':       type_name,
            'group':      group_of(source_file),
            'origin':     'default' if element is None or inline else 'file',
            'reached':    reached,
            'path':       path_of(element) if element is not None and not inline else None,
            'sourceFile': (source_file or '').replace('./source/', '') or None,
            'parameters': parameters_of(type_name, element) if type_name else [],
        }
        if element is not None and element.tag != class_name:
            node['element'] = element.tag
        nodes[key] = node
        by_id[node['id']] = node
        order.append(node)
        if type_name is None:
            return node['id']
        scope = root if element is None else element
        for obj in schema(type_name):
            build(node['id'], scope, obj, reached)
        return node['id']

    def build(source_id, scope, obj, reached):
        """Resolve the object `obj`, built from the parameters of `scope`."""
        group = obj.get('sourceElement')
        if group:
            # Built from a named group of sub-parameters, where present.
            groups = [child for child in _elements(scope) if child.tag == group]
            scope = groups[0] if groups else scope
        found = locate(scope, obj['parameterName'], obj.get('repeatable', False))
        if not found and obj['parameterName'] != obj['class']:
            # Only an object named for its class has a class default.
            return
        for element in found or [None]:
            add_edge(source_id, visit(obj['class'], element, reached), obj['parameterName'])

    # The task is the root of the model: all else is built on its behalf.
    tasks = [element for element in map(dereference, _elements(root))
             if element is not None and element.tag == 'task']
    for element in tasks or [None]:
        visit('task', element, True)

    # Node components are built for every model which evolves nodes.
    if source_root is not None:
        components = _scan_components(source_root)
        classes = sorted({key[0] for key in components})
        for class_name in classes:
            element_name = 'component' + _upper_first(class_name)
            selected = [child for child in _elements(root) if child.tag == element_name]
            element = dereference(selected[0]) if selected else None
            if element is not None and element.get('value'):
                key = (class_name, element.get('value'))
            else:
                defaults = [k for k, c in components.items() if k[0] == class_name and c['isDefault']]
                key = defaults[0] if defaults else None
            component = components.get(key)
            if component is None or (key[1] == 'null'):
                continue
            # A component inherits what the component it extends builds.
            chain, current = [], key
            while current in components and current not in chain:
                chain.append(current)
                current = components[current]['extends']
            objects, parameters = [], []
            for link in chain:
                objects.extend(o for o in components[link]['objects'] if o not in objects)
                parameters.extend(p for p in components[link]['parameters']
                                  if p['name'] not in [q['name'] for q in parameters])
            object_names = {obj['parameterName'] for obj in objects}
            values = scalar_values(element, object_names | set(function_classes))
            node_parameters = []
            for parameter in parameters:
                value = values.pop(parameter['name'], None)
                node_parameters.append({
                    'name':        parameter['name'],
                    'value':       value if value is not None else parameter['default'],
                    'set':         value is not None,
                    'description': parameter['description'],
                })
            for parameter_name, value in values.items():
                node_parameters.append({'name': parameter_name, 'value': value, 'set': True, 'description': None})
            node = {
                'id':          f'n{len(order)}',
                'kind':        'component',
                'class':       element_name,
                'className':   f'{_upper_first(class_name)} component',
                'label':       key[1],
                'type':        None,
                'group':       'components',
                'origin':      'file' if selected else 'default',
                'reached':     True,
                'path':        path_of(element) if element is not None else None,
                'sourceFile':  component['sourceFile'],
                'description': component['description'],
                'parameters':  node_parameters,
            }
            by_id[node['id']] = node
            order.append(node)
            for obj in objects:
                scope = element if (obj['own'] and element is not None) else root
                found = locate(scope, obj['parameterName'])
                if found:
                    targets = [visit(obj['class'], found[0], True)]
                elif obj['parameterName'] == obj['class']:
                    targets = [visit(obj['class'], None, True)]
                elif obj['default']:
                    targets = [visit(obj['class'], _default_element(obj['parameterName'], obj['default']), True,
                                     key=('inline', node['id'], obj['parameterName']))]
                else:
                    targets = []
                for target in targets:
                    add_edge(node['id'], target, obj['parameterName'])

    # Whatever remains in the file selects an object which nothing above was
    # seen to build: it is used by hand-written code, or not at all.
    parameter_classes = {}
    for implementation in implementations.values():
        for obj in implementation.get('objects', []):
            parameter_classes.setdefault(obj['parameterName'], set()).add(obj['class'])
    for element in root.iter():
        if not isinstance(element.tag, str) or element is root or element.get('idRef') is not None:
            continue
        label = element.get('value')
        if label is None:
            continue
        candidates = [element.tag] if element.tag in function_classes else sorted(
            parameter_classes.get(element.tag, ()))
        for class_name in candidates:
            if (class_name, label) in type_by_label:
                visit(class_name, element, False)
                break

    return {
        'model':  {'name': name},
        'groups': [{'id': group, 'title': title, 'blurb': blurb} for group, title, blurb in GROUPS],
        'nodes':  order,
        'edges':  edges,
    }
