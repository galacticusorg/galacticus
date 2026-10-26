"""Assemble the Galacticus parameter catalog from the Fortran source tree.

Andrew Benson (2026)

The catalog maps every ``functionClass`` implementation to the input
parameters it accepts and the nested objects it builds.  It is assembled in
three phases:

  1. Discover the ``functionClass`` base classes (``<functionClass>`` directives).
  2. Discover the implementations of each base (``<base name="impl"/>``
     registration directives) and the files they live in.
  3. For each implementation, locate its ``*ConstructorParameters`` constructor,
     harvest the ``<inputParameter>`` / ``<objectBuilder>`` directives in that
     constructor's scope, and infer a type for each parameter.

The harvesting reuses the existing SourceTree parser so that parameter→variable
declaration lookups (for type inference) are available, and so the directive
scoping (which constructor a directive belongs to) is exact.

Phases 1 and 2, and the pre-filter that opens phase 3, all need nothing but the
directives of a file, so they share a single read of it (`_scan_directives`);
phase 3 then re-reads only those files that actually register an implementation,
to build a full parse tree. Both passes run over a process pool
(`Galacticus.Build.ParallelScan`), which also sizes itself against make's
jobserver. Results are collected in source-file order, so the catalog is
identical whatever the parallelism.
"""

import os
import re
from collections import Counter

from Galacticus.Build import FileChanges
from Galacticus.Build.Capabilities import parse_requires, parse_withholds
from Galacticus.Build.Directives import extract_directives
from Galacticus.Build.ParallelScan import scan as parallel_scan
from Galacticus.Build.SourceTree import parse_file, walk_tree
from Galacticus.Build.SourceTree.Parse.Declarations import get_declaration
from Galacticus.Parameters.inference import infer_parameter_type, _scalar

_DIRECTIVE_TYPES = ('inputParameter', 'objectBuilder')

# ---------------------------------------------------------------------------
# Read-only data shared with forked workers via copy-on-write. Populated in
# `build_catalog` BEFORE `parallel_scan` is called, so each worker inherits it
# without anything being pickled per task.
# ---------------------------------------------------------------------------
_WORKER = {}


def _iter_source_files(source_root):
    """Yield F90 source files under `source_root`, deterministically ordered."""
    for dirpath, dirnames, filenames in os.walk(source_root):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith('.'))
        for name in sorted(filenames):
            if name.endswith('.F90') and not name.startswith('.#'):
                yield os.path.join(dirpath, name)


def _base_entry(directive, path):
    """Build a phase-1 base-class record from a `functionClass` directive.

    Describes only what the directive itself says. The `implementations` list is
    deliberately NOT included: it is accumulated later, and owning a mutable
    accumulator here would let one base class's labels leak into another's.
    """
    return {
        'default':         _scalar(directive.get('default')),
        'descriptiveName': _scalar(directive.get('descriptiveName')),
        'sourceFile':      path,
    }


def _base_record(entry):
    """A fresh, private base-class record ready to accumulate implementations.

    Copies `entry`, which the caller must be able to reuse unmutated.
    """
    return {**entry, 'implementations': []}


def _enumeration_labels(directive):
    """Extract the label list from an `enumeration` directive, or None."""
    entries = directive.get('entry')
    if entries is None:
        return None
    if isinstance(entries, dict):
        entries = [entries]
    return [e.get('label') for e in entries
            if isinstance(e, dict) and e.get('label')]


def _scan_directives(path):
    """Read one source file once and return everything the early phases need.

    Returns ``(bases, enumerations, named_roots)`` where `bases` and
    `enumerations` are ``(name, value)`` lists in file order, and `named_roots`
    is the set of root element names of every directive carrying a `name`. The
    caller intersects `named_roots` with the discovered base classes to decide,
    without another read, whether this file registers an implementation.
    """
    bases        = []
    enumerations = []
    named_roots  = set()
    for directive in extract_directives(path, '*', set_root_element_type=True):
        root = directive.get('rootElementType')
        name = _scalar(directive.get('name'))
        if name:
            named_roots.add(root)
        if root == 'functionClass':
            if name:
                bases.append((name, _base_entry(directive, path)))
        elif root == 'enumeration':
            labels = _enumeration_labels(directive)
            if name and labels is not None:
                enumerations.append((name, labels))
    return bases, enumerations, named_roots


def discover_base_classes(source_files):
    """Phase 1: return ``{base_name: {default, descriptiveName, sourceFile}}``."""
    bases = {}
    for path in source_files:
        for directive in extract_directives(path, 'functionClass'):
            name = _scalar(directive.get('name'))
            if not name:
                continue
            bases[name] = _base_record(_base_entry(directive, path))
    return bases


def derive_label(base_name, implementation_type):
    """Derive the parameter-selector label for an implementation type.

    Mirrors the canonical ``lcfirst-unless-all-caps`` convention used by
    ``_short_name`` in
    python/Galacticus/Build/SourceTree/Process/FunctionClass/__init__.py
    (and replicated by loadAndSortClasses / buildObjectTypeMethod /
    generateDocumentation): strip the base-class prefix, then lowercase the
    first character *unless* the remainder begins with an all-caps acronym.

    e.g. ``accretionHaloSimple`` -> ``simple`` but
    ``darkMatterProfileDMONFW`` -> ``NFW`` and
    ``darkMatterProfileDMOSIDMCoreNFW`` -> ``SIDMCoreNFW``.
    """
    suffix = implementation_type
    if suffix.startswith(base_name):
        suffix = suffix[len(base_name):]
    if not suffix:
        return implementation_type
    if not re.match(r'^[A-Z]{2,}', suffix):
        suffix = suffix[0].lower() + suffix[1:]
    return suffix


def _module_name(tree):
    for node in walk_tree(tree):
        if node.get('type') == 'module' and node.get('name'):
            return node['name']
    return None


def _find_constructor(function_nodes, label, implementation_count):
    """Return the constructor function node for an implementation, or None."""
    want = (label + 'constructorparameters').lower()
    for node in function_nodes:
        if (node.get('name') or '').lower() == want:
            return node
    # Fallback: a lone implementation in a file with a single *ConstructorParameters.
    candidates = [
        node for node in function_nodes
        if (node.get('name') or '').lower().endswith('constructorparameters')
    ]
    if implementation_count == 1 and len(candidates) == 1:
        return candidates[0]
    return None


# A local parameter handle is bound to a nested element via
#   <var> = <something>%subParameters('elementName' [, ...])
# e.g. `subParameters = parameters%subParameters('componentHotHalo')` or
# `parametersMassDefinitions = parameters%subParameters('massDefinitions')`.
_SUBPARAMETERS_RE = re.compile(
    r"(\w+)\s*=\s*[\w%]+%subParameters\s*\(\s*['\"]([^'\"]+)['\"]")

# An implementation may extend another implementation (not just the base class),
# inheriting its parameters: `type, extends(parentType) :: implType`.
_TYPE_EXTENDS_RE = re.compile(
    r"type\s*,\s*extends\s*\(\s*(\w+)\s*\)\s*::\s*(\w+)", re.IGNORECASE)

# Parameters/sub-parameter wrappers read directly in hand-written constructor
# code (not via `<inputParameter>`/`<objectBuilder>` directives), e.g.
# `parameters%copiesCount('particleProperty')` or `parameters%isPresent('x')`.
_DIRECT_PARAMETER_RE = re.compile(
    r"%(?:value|subParameters|copiesCount|count|isPresent)\s*\(\s*['\"]([^'\"]+)['\"]")


def _parent_types(tree):
    """Map each implementation type to the type it extends (case as written).

    Used to resolve implementation inheritance; the caller decides whether the
    parent is another implementation (inherit its parameters) or the base
    ``<class>Class`` (chain root).
    """
    parents = {}
    for node in walk_tree(tree):
        if node.get('type') != 'type':
            continue
        match = _TYPE_EXTENDS_RE.search(node.get('opener', '') or '')
        if match:
            parents[match.group(2)] = match.group(1)
    return parents


def _direct_parameter_names(scope_node):
    """Names read directly from the parameter tree in hand-written constructor
    code (complements the directive-harvested names; only widens what's
    accepted)."""
    names = set()
    for node in walk_tree(scope_node):
        if node.get('type') != 'code':
            continue
        for match in _DIRECT_PARAMETER_RE.finditer(node.get('content', '') or ''):
            names.add(match.group(1))
    return sorted(names)


# A string parameter that is really an enumeration is passed through
# `enumeration<Name>Encode(char(var))` in the constructor (the generated encoder,
# Enumeration.py).  This links the local variable to its enumeration's name.
_ENUM_ENCODE_RE = re.compile(
    r"enumeration(\w+)Encode\s*\(\s*(?:char\s*\(\s*)?([A-Za-z_]\w*)")


def discover_enumerations(source_files):
    """Return ``{enumerationName: [labels]}`` from all `<enumeration>` directives."""
    enumerations = {}
    for path in source_files:
        for directive in extract_directives(path, 'enumeration'):
            name = _scalar(directive.get('name'))
            if not name:
                continue
            labels = _enumeration_labels(directive)
            if labels is None:
                continue
            enumerations[name] = labels
    return enumerations


# A string default is written in Fortran as `var_str('value')`; the value a user
# actually enters in a parameter file is just the inner string.
_VAR_STR_RE = re.compile(r"""^var_str\s*\(\s*(['"])(.*)\1\s*\)$""", re.IGNORECASE)


def _normalize_default(value):
    """Render a default value as a user would enter it: unwrap `var_str('x')`
    to `x`.  Leaves numeric/boolean literals and expressions unchanged."""
    if value is None:
        return None
    match = _VAR_STR_RE.match(value.strip())
    return match.group(2) if match else value


def _constraint_bound(raw):
    """Parse a `<minimum>`/`<maximum>` directive value into ``{value, inclusive}``.

    `raw` is the `xml_to_dict` value: a plain string (``<minimum>0.0</minimum>``),
    or a dict carrying the `inclusive` attribute
    (``<minimum inclusive="false">0.0</minimum>`` -> ``{'inclusive':..,'content':..}``).
    """
    if raw is None:
        return None
    if isinstance(raw, dict):
        value = _scalar(raw)
        inclusive = raw.get('inclusive', 'true') != 'false'
    else:
        value = str(raw).strip()
        inclusive = True
    if not value:
        return None
    return {'value': value, 'inclusive': inclusive}


def _capture_constraints(directive):
    """Capture optional validate-time constraints (`minimum`/`maximum`/
    `allowedValues`) from an inputParameter directive, or None if it has none."""
    minimum = _constraint_bound(directive.get('minimum'))
    maximum = _constraint_bound(directive.get('maximum'))
    allowed_text = _scalar(directive.get('allowedValues'))
    allowed_values = allowed_text.split() if allowed_text else None
    if minimum is None and maximum is None and not allowed_values:
        return None
    constraints = {}
    if minimum is not None:
        constraints['minimum'] = minimum
    if maximum is not None:
        constraints['maximum'] = maximum
    if allowed_values:
        constraints['allowedValues'] = allowed_values
    return constraints


def _enumeration_links(scope_node):
    """Map a local variable to the enumeration it is encoded into, by scanning
    constructor code for ``enumeration<Name>Encode(char(var))`` calls.  The
    captured name is `_ucfirst`-cased; lower the first letter to recover the
    enumeration directive's name."""
    links = {}
    for node in walk_tree(scope_node):
        if node.get('type') != 'code':
            continue
        for match in _ENUM_ENCODE_RE.finditer(node.get('content', '') or ''):
            captured = match.group(1)
            enum_name = captured[0].lower() + captured[1:]
            links.setdefault(match.group(2), enum_name)
    return links


def _resolve_source_elements(scope_node):
    """Map each local parameter-source variable to the nested element it reads.

    A directive ``<source>X</source>`` where X is the constructor's primary
    ``parameters`` argument reads from the class's own element; any other X is a
    local handle bound to a *nested* element (e.g. ``massDefinitions``).  This
    returns ``{X: elementName}`` for those bound handles so validation knows the
    extra level of nesting; the primary ``parameters`` handle is absent (its
    parameters live directly in the class's element).
    """
    mapping = {}
    for node in walk_tree(scope_node):
        if node.get('type') != 'code':
            continue
        for match in _SUBPARAMETERS_RE.finditer(node.get('content', '') or ''):
            mapping[match.group(1)] = match.group(2)
    return mapping


def _declaration_lookup(constructor_node):
    """Return a ``f(variable) -> declaration_dict`` closure for type inference."""
    if constructor_node is None:
        return lambda variable: None

    def lookup(variable):
        try:
            return get_declaration(constructor_node, variable)
        except Exception:
            return None

    return lookup


def _harvest_parameters(scope_node, declaration_lookup, source_elements,
                        enumeration_links):
    """Collect input parameters and object-builder edges under `scope_node`.

    ``source_elements`` maps a local parameter-source variable to the nested
    element it reads (see `_resolve_source_elements`).  Each parameter/object
    records ``sourceElement`` -- the nested element it lives in, or None when it
    lives directly in the class's own element.  ``enumeration_links`` maps a
    variable to its enumeration name (see `_enumeration_links`); a matched
    parameter records ``enumeration`` so values can be checked against the
    enumeration's labels.
    """
    parameters = []
    objects    = []
    for node in walk_tree(scope_node):
        directive = node.get('directive')
        if not isinstance(directive, dict):
            continue
        if node.get('type') == 'inputParameter':
            inferred = infer_parameter_type(directive, declaration_lookup)
            source = _scalar(directive.get('source'))
            variable = _scalar(directive.get('variable')) or _scalar(directive.get('name'))
            parameters.append({
                'name':          _scalar(directive.get('name')),
                'variable':      _scalar(directive.get('variable')),
                'source':        source,
                'sourceElement': source_elements.get(source),
                'type':          inferred['type'],
                'kind':          inferred['kind'],
                'provenance':    inferred['provenance'],
                'enumeration':   enumeration_links.get(variable),
                'constraints':   _capture_constraints(directive),
                'default':       _normalize_default(_scalar(directive.get('defaultValue'))),
                'cardinality':   _scalar(directive.get('cardinality')),
                'description':   _scalar(directive.get('description')),
            })
        elif node.get('type') == 'objectBuilder':
            class_name = _scalar(directive.get('class'))
            parameter_name = _scalar(directive.get('parameterName')) or class_name
            source = _scalar(directive.get('source'))
            objects.append({
                'class':         class_name,
                'name':          _scalar(directive.get('name')),
                'parameterName': parameter_name,
                'source':        source,
                'sourceElement': source_elements.get(source),
                'repeatable':    _scalar(directive.get('copy')) is not None,
                # Optional arguments of this object's methods which the
                # implementation never supplies (see `Galacticus.Build.Capabilities`).
                'withholds':     [{'method': method, 'argument': argument}
                                  for method, argument in
                                  parse_withholds(_scalar(directive.get('withholds')))],
            })

    # Drop ambiguous enumeration links: when one variable name is shared by
    # several input parameters (e.g. a generic `name` read in many nested blocks,
    # as in taskMergerTreeFileBuilder), a single `Encode` call cannot be
    # attributed to one of them, so attaching the enumeration to all would
    # mis-type the others.  Only keep unambiguous links.
    variable_counts = Counter(
        (p['variable'] or p['name']) for p in parameters
        if (p['variable'] or p['name'])
    )
    for parameter in parameters:
        key = parameter['variable'] or parameter['name']
        if parameter.get('enumeration') and key and variable_counts[key] > 1:
            parameter['enumeration'] = None

    return parameters, objects


def harvest_file(path, base_names, source_root):
    """Phase 3 for one file: return a list of implementation catalog entries."""
    tree = parse_file(path)
    module = _module_name(tree)
    relative = os.path.relpath(path, source_root)

    registrations = [
        node for node in walk_tree(tree)
        if node.get('type') in base_names
        and isinstance(node.get('directive'), dict)
        and _scalar(node['directive'].get('name'))
    ]
    if not registrations:
        return []

    function_nodes = [
        node for node in walk_tree(tree)
        if node.get('type') in ('function', 'subroutine')
    ]
    parents = _parent_types(tree)

    entries = []
    for registration in registrations:
        base = registration['type']
        implementation_type = _scalar(registration['directive']['name'])
        label = derive_label(base, implementation_type)
        # An implementation whose Fortran name had to be abbreviated to stay
        # within the 63-character identifier limit may carry a readable
        # `alias`; both it and the abbreviated label select the implementation.
        alias = _scalar(registration['directive'].get('alias') or '') or None
        constructor = _find_constructor(function_nodes, label, len(registrations))
        scope = constructor if constructor is not None else tree
        parameters, objects = _harvest_parameters(
            scope, _declaration_lookup(constructor),
            _resolve_source_elements(scope), _enumeration_links(scope))
        requires, forwards = parse_requires(registration['directive'])
        entries.append({
            'type':            implementation_type,
            'functionClass':   base,
            'label':           label,
            'alias':           alias,
            'module':          module,
            'sourceFile':      relative,
            'constructorFound': constructor is not None,
            # The type this implementation extends (another implementation type,
            # whose parameters are inherited, or the base `<class>Class` root).
            'parent':          parents.get(implementation_type),
            'parameters':      parameters,
            'objects':         objects,
            # Names read directly in hand-written constructor code (not via
            # directives); accepted in addition to `parameters`/`objects`.
            'directNames':     _direct_parameter_names(scope),
            # Optional method arguments this implementation can not do without,
            # and the objects whose requirements it inherits by passing its own
            # arguments on to them (see `Galacticus.Build.Capabilities`).
            'requires':        [{'method': method, 'argument': argument}
                                for method, argument in requires],
            'forwards':        forwards,
        })
    return entries


def _harvest_worker(path):
    """Phase 3 for one file, in a pool worker: return ``(entries, error)``.

    A parse failure is returned rather than raised: one unparseable file must
    not abort the whole catalog (and would otherwise take the pool down with
    it), so the caller reports it and carries on -- matching the serial
    behavior this replaced.
    """
    try:
        entries = harvest_file(path, _WORKER['base_names'], _WORKER['source_root'])
    except Exception as exc:  # noqa: BLE001 -- reported by the caller.
        return [], str(exc)
    return entries, None


def build_catalog(source_root, log=None, jobs=None):
    """Build the full catalog dict for the source tree rooted at `source_root`.

    Parameters
    ----------
    source_root : str
        Path to the ``source`` directory.
    log : callable, optional
        ``f(message)`` for progress/diagnostics.
    jobs : int, optional
        Worker count for the parallel passes. Defaults to `ParallelScan`'s own
        resolution (``GALACTICUS_BUILD_JOBS``, then make's ``-j``/jobserver).

    Returns
    -------
    dict
        ``{'functionClasses': {...}, 'implementations': {...}}``.
    """
    def _log(message):
        if log is not None:
            log(message)

    source_files = list(_iter_source_files(source_root))

    # --- Phases 1 and 2, plus the phase-3 pre-filter, from one read per file.
    _log(f"scanning {len(source_files)} source files for functionClass bases")
    scans = dict(zip(source_files, parallel_scan(
        source_files, _scan_directives, 'Galacticus.Parameters.catalog', jobs=jobs)))

    bases        = {}
    enumerations = {}
    for path in source_files:                      # file order: later files win,
        file_bases, file_enumerations, _ = scans[path]   # exactly as the serial
        for name, entry in file_bases:                   # loop this replaced.
            bases[name] = _base_record(entry)
        for name, labels in file_enumerations:
            enumerations[name] = labels

    base_names = set(bases)
    _log(f"found {len(base_names)} functionClass base classes")
    _log(f"found {len(enumerations)} enumerations")

    # --- Phase 3. A file needs a full parse only if it registers an
    # implementation, i.e. carries a named directive whose root is a base class.
    registering  = [path for path in source_files if scans[path][2] & base_names]
    harvests     = {}
    to_harvest   = registering
    _log(f"parsing {len(to_harvest)} files that register implementations")

    _WORKER['base_names']  = base_names
    _WORKER['source_root'] = source_root
    for path, result in zip(to_harvest, parallel_scan(
            to_harvest, _harvest_worker, 'Galacticus.Parameters.catalog', jobs=jobs)):
        harvests[path] = result

    implementations = {}
    parsed_files = 0
    for path in registering:
        entries, error = harvests[path]
        if error is not None:
            _log(f"  [parse-error] {os.path.relpath(path, source_root)}: {error}")
            continue
        parsed_files += 1
        for entry in entries:
            implementations[entry['type']] = entry
            base = entry['functionClass']
            if base in bases:
                for selector in (entry['label'], entry.get('alias')):
                    if selector and selector not in bases[base]['implementations']:
                        bases[base]['implementations'].append(selector)

    for base in bases.values():
        base['implementations'].sort()

    _log(f"parsed {parsed_files} implementation files; "
         f"cataloged {len(implementations)} implementations")
    return {
        'functionClasses': bases,
        'implementations': implementations,
        'enumerations':    enumerations,
    }
