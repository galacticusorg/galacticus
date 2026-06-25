"""Validate a Galacticus parameter file against the parameter catalog.

Andrew Benson (2026)

Walks a parameter XML tree hierarchically and reports three classes of problem:

  * **selector**  -- a ``<functionClassBase value="X"/>`` whose ``X`` does not
    name a real implementation of that base (case-sensitive, matching the
    generated `select case` dispatch).  Catches typos, wrong casing
    (``Tinker2008`` vs ``tinker2008``), and non-existent values (``null`` on an
    objectBuilder-built class).
  * **unknown**   -- a scalar parameter, *inside a functionClass element*, whose
    name is not accepted by the selected implementation.
  * **type**      -- a scalar whose literal value does not parse as the
    parameter's inferred type.

Design notes / why this is low-false-positive:

  * Scalar `parameters%value` reads do NOT walk up the parameter tree (unlike
    objectBuilder), so a scalar always lives in its owning object's element.
    We therefore only name-check scalars *within* a functionClass scope; global
    / meta parameters (`formatVersion`, `verbosityLevel`, ...) live at the root
    scope and are never checked, so they cannot be false-flagged.
  * Objects CAN be hoisted to a higher element (objectBuilder walks up), so an
    element child that is itself a functionClass selector is validated and
    recursed into but never flagged as "unexpected" -- hoisted/shared objects
    are legitimate.
  * `sourceElement` wrappers (e.g. ``<massDefinitions>``) are descended into.
  * Values that are not simple literals (expressions, references, lists) are not
    type-checked -- only an unambiguous literal mismatch is reported.
"""

import os
import re
import xml.etree.ElementTree as ET
import collections
from collections import namedtuple

# A single validation finding.  `level` is 'error' (selector/unknown) or
# 'warning' (type, which is softer because value syntax is permissive).
Finding = namedtuple('Finding', ['level', 'kind', 'path', 'message'])

# Element tags that are structural/meta, never functionClass selectors or
# scalar parameters.
_META_TAGS = frozenset({'formatVersion', 'lastModified'})

# --- Migration-aware diagnostics -------------------------------------------
# Galacticus does NOT migrate at run time (project policy: users run
# scripts/aux/parametersMigrate.py and re-verify).  We therefore validate files
# as-written, but use scripts/aux/migrations.xml to turn an "unknown parameter"
# or unrecognised selector value into an actionable "obsolete; renamed to ..."
# message.  Cached per GALACTICUS_EXEC_PATH so tests can override it.
_migrations_cache = {}
_XPATH_LEAF_RE = re.compile(r'([A-Za-z_]\w*)(?:\[[^\]]*\])?\s*$')


def _load_migrations(exec_path):
    """Return ``(name_renames, value_renames)`` parsed from migrations.xml:
    name_renames maps an obsolete parameter (leaf) name to its new name; value
    renames map an obsolete selector value to its new value."""
    name_renames, value_renames = {}, {}
    if not exec_path:
        return name_renames, value_renames
    path = os.path.join(exec_path, 'scripts', 'aux', 'migrations.xml')
    if not os.path.exists(path):
        return name_renames, value_renames
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError:
        return name_renames, value_renames
    for translation in root.iter('translation'):
        xpath = translation.get('xpath', '') or ''
        name_element = translation.find('name')
        if name_element is not None and name_element.get('new'):
            leaf = _XPATH_LEAF_RE.search(xpath)
            if leaf:
                name_renames.setdefault(leaf.group(1), name_element.get('new'))
        for value_element in translation.findall('value'):
            old, new = value_element.get('old'), value_element.get('new')
            if old and new:
                value_renames.setdefault(old, new)
    return name_renames, value_renames


def _migrations():
    exec_path = os.environ.get('GALACTICUS_EXEC_PATH', '')
    if exec_path not in _migrations_cache:
        _migrations_cache[exec_path] = _load_migrations(exec_path)
    return _migrations_cache[exec_path]


def _obsolete_hint(new_name):
    return (f" (obsolete; renamed to '{new_name}' — run "
            f"scripts/aux/parametersMigrate.py to update)")


class _Index:
    """Pre-computed lookups over a catalog for fast validation."""

    def __init__(self, catalog):
        self.function_classes = catalog['functionClasses']
        self.implementations = catalog['implementations']
        self.enumerations = catalog.get('enumerations', {})
        self.base_names = set(self.function_classes)
        # (base, label) -> implementation type name.
        self.type_by_base_label = {}
        for type_name, impl in self.implementations.items():
            self.type_by_base_label[(impl['functionClass'], impl['label'])] = type_name
        self._schema_cache = {}

    def labels_for(self, base):
        return set(self.function_classes.get(base, {}).get('implementations', []))

    def default_label(self, base):
        return self.function_classes.get(base, {}).get('default')

    def implementation(self, base, label):
        return self.implementations.get(self.type_by_base_label.get((base, label)))

    def _inheritance_chain(self, impl_type):
        """Yield `impl_type` and each ancestor implementation it extends, in
        order.  Stops when the parent is not a known implementation (i.e. the
        base ``<class>Class`` root)."""
        seen = set()
        current = impl_type
        while current in self.implementations and current not in seen:
            seen.add(current)
            yield current
            current = self.implementations[current].get('parent')

    def schema(self, impl_type):
        """Return the element-structure schema for an implementation type,
        merged across its inheritance chain: scalars (own element),
        scalar_groups / object_groups (keyed by the nested `sourceElement`),
        objects (own element, keyed by child element name = parameterName), and
        `extras` (names read directly in hand-written code)."""
        if impl_type in self._schema_cache:
            return self._schema_cache[impl_type]
        scalars, scalar_groups = {}, {}
        objects, object_groups = {}, {}
        extras = set()
        # Walk the chain; a subclass implementation accepts its ancestors'
        # parameters too.  Nearer definitions take precedence (set first).
        for type_name in self._inheritance_chain(impl_type):
            impl = self.implementations.get(type_name, {})
            for parameter in impl.get('parameters', []):
                name = parameter.get('name')
                if not name:
                    continue
                element = parameter.get('sourceElement')
                if element is None:
                    scalars.setdefault(name, parameter)
                else:
                    scalar_groups.setdefault(element, {}).setdefault(name, parameter)
            for obj in impl.get('objects', []):
                parameter_name = obj.get('parameterName')
                if not parameter_name:
                    continue
                element = obj.get('sourceElement')
                if element is None:
                    objects.setdefault(parameter_name, obj)
                else:
                    object_groups.setdefault(element, {}).setdefault(parameter_name, obj)
            extras.update(impl.get('directNames', []))
        schema = {
            'scalars': scalars, 'scalar_groups': scalar_groups,
            'objects': objects, 'object_groups': object_groups,
            'extras': extras,
        }
        self._schema_cache[impl_type] = schema
        return schema


def _is_simple_literal(value):
    """True when `value` is a single literal token we can type-check (not an
    expression, reference, list, or path)."""
    v = value.strip()
    if not v or any(c.isspace() for c in v):      # lists / multi-token values
        return False
    if any(c in v for c in '[]%()'):              # references, xpath, calls
        return False
    if '//' in v:                                 # string concatenation
        return False
    return True


def _enumeration_allowed(enum_name, labels):
    """Accepted values for an enumeration: the bare labels and their prefixed
    forms (`<enum><Ucfirst(label)>`), covering both `includesPrefix` modes of
    the generated encoder."""
    allowed = set(labels)
    for label in labels:
        if label:
            allowed.add(enum_name + label[0].upper() + label[1:])
    return allowed


def _to_number(text):
    """Parse a Fortran/plain numeric literal to a float, or None."""
    try:
        return float(text.strip().replace('d', 'e').replace('D', 'e'))
    except (ValueError, AttributeError):
        return None


def _value_matches_type(value, ctype):
    """Return True/False if `value` clearly matches/violates `ctype`, or None
    when we decline to judge (non-literal, or a type we don't check)."""
    if ctype in (None, 'string', 'object', 'unknown'):
        return None
    if not _is_simple_literal(value):
        return None
    v = value.strip()
    if ctype == 'boolean':
        return v.lower() in ('true', 'false')
    if ctype == 'integer':
        try:
            int(v)
            return True
        except ValueError:
            return None        # might be an expression; do not judge
    if ctype == 'real':
        try:
            float(v.replace('d', 'e').replace('D', 'e'))
            return True
        except ValueError:
            return None
    return None


def validate_parameters(root, catalog):
    """Validate a ``<parameters>`` element tree against `catalog`.

    Parameters
    ----------
    root : xml.etree.ElementTree.Element
        The ``<parameters>`` root element.
    catalog : dict
        A loaded `parameters.catalog.json`.

    Returns
    -------
    list[Finding]
    """
    index = _Index(catalog)
    findings = []
    name_renames, value_renames = _migrations()

    def child_value(element):
        return element.get('value')

    def resolve_label(base, element):
        """The implementation label selected for `element` (its `value`, or the
        base default when absent)."""
        value = child_value(element)
        return value if value is not None else index.default_label(base)

    def validate_selector(base, element, path):
        """Check `element`'s `value` names a real implementation of `base`.
        Returns the resolved implementation type (or None)."""
        value = child_value(element)
        labels = index.labels_for(base)
        if value is None:
            label = index.default_label(base)
        elif value in labels:
            label = value
        else:
            # Distinguish wrong-casing from genuinely unknown, for a better message.
            lowered = {l.lower(): l for l in labels}
            if value.lower() in lowered:
                findings.append(Finding(
                    'error', 'selector', path,
                    f"selector value '{value}' should be '{lowered[value.lower()]}' "
                    f"(wrong case) for class '{base}'"))
            else:
                hint = _obsolete_hint(value_renames[value]) if value in value_renames else ""
                findings.append(Finding(
                    'error', 'selector', path,
                    f"selector value '{value}' is not an implementation of "
                    f"class '{base}'{hint}"))
            return None
        return index.type_by_base_label.get((base, label))

    def walk_scope(element, impl_type, path):
        """Validate the children of `element`, whose selected implementation is
        `impl_type` (None at the root scope, where only objects/globals live)."""
        schema = index.schema(impl_type) if impl_type else None
        for child in element:
            if not isinstance(child.tag, str):       # comments / PIs
                continue
            tag = child.tag
            child_path = f"{path}/{tag}"

            if tag in _META_TAGS:
                continue

            # 1. A scalar declared by this implementation (own element).
            if schema and tag in schema['scalars']:
                _check_type(child, schema['scalars'][tag], child_path)
                continue

            # 2. A sourceElement wrapper (e.g. <massDefinitions>): descend and
            #    validate its contents against the grouped scalars/objects.
            if schema and (tag in schema['scalar_groups']
                           or tag in schema['object_groups']):
                walk_group(child, schema['scalar_groups'].get(tag, {}),
                           schema['object_groups'].get(tag, {}), child_path)
                continue

            # 3. An object declared by this implementation under a differently
            #    named element (parameterName != class).
            if schema and tag in schema['objects']:
                base = schema['objects'][tag]['class']
                resolved = validate_selector(base, child, child_path)
                walk_scope(child, resolved, child_path)
                continue

            # 4. A functionClass selector element (own object with
            #    parameterName == class, or a hoisted/shared object). Validate
            #    and recurse, but never flag as unexpected (hoisting is legal).
            if tag in index.base_names:
                resolved = validate_selector(tag, child, child_path)
                walk_scope(child, resolved, child_path)
                continue

            # 5. A name read directly in hand-written constructor code (e.g. a
            #    repeated `<particleProperty>` sub-parameter block). Accept it;
            #    its internal structure is opaque to the catalog, so do not
            #    recurse or flag its children.
            if schema and tag in schema['extras']:
                continue

            # 6. Anything else inside a functionClass scope is an unknown
            #    parameter name.  (At the root scope we do not flag -- globals
            #    and meta parameters legitimately live there.)  Honour
            #    `ignoreWarnings="true"`, Galacticus's own opt-out for
            #    deliberately-unusual parameters (input_parameters.F90).
            if schema is not None and child.get('ignoreWarnings') != 'true':
                hint = _obsolete_hint(name_renames[tag]) if tag in name_renames else ""
                findings.append(Finding(
                    'error', 'unknown', child_path,
                    f"parameter '{tag}' is not accepted by implementation "
                    f"'{impl_type}'{hint}"))

    def walk_group(element, scalars, objects, path):
        """Validate a sourceElement wrapper's children against its grouped
        scalars/objects."""
        for child in element:
            if not isinstance(child.tag, str):
                continue
            tag = child.tag
            child_path = f"{path}/{tag}"
            if tag in _META_TAGS:
                continue
            if tag in scalars:
                _check_type(child, scalars[tag], child_path)
            elif tag in objects:
                base = objects[tag]['class']
                resolved = validate_selector(base, child, child_path)
                walk_scope(child, resolved, child_path)
            elif tag in index.base_names:
                resolved = validate_selector(tag, child, child_path)
                walk_scope(child, resolved, child_path)
            # Grouped wrappers (labels etc.) may also carry plain values we
            # cannot attribute; leave anything else unflagged here.

    def _check_type(element, parameter, path):
        value = child_value(element)
        if value is None or element.get('ignoreWarnings') == 'true':
            return
        # Enumeration-valued parameter: the string must be one of the labels.
        enumeration = parameter.get('enumeration')
        if enumeration and enumeration in index.enumerations:
            if _is_simple_literal(value):
                labels = index.enumerations[enumeration]
                if value.strip() not in _enumeration_allowed(enumeration, labels):
                    findings.append(Finding(
                        'error', 'enumeration', path,
                        f"value '{value}' is not valid for enumeration "
                        f"'{enumeration}' for parameter '{parameter.get('name')}' "
                        f"(allowed: {', '.join(sorted(labels))})"))
            return
        verdict = _value_matches_type(value, parameter.get('type'))
        if verdict is False:
            findings.append(Finding(
                'warning', 'type', path,
                f"value '{value}' is not a valid {parameter.get('type')} for "
                f"parameter '{parameter.get('name')}'"))
        _check_constraints(value, parameter, path)

    def _check_constraints(value, parameter, path):
        """Enforce optional annotated constraints (allowedValues / numeric
        range) on a simple-literal value."""
        constraints = parameter.get('constraints')
        if not constraints or not _is_simple_literal(value):
            return
        literal = value.strip()
        name = parameter.get('name')

        allowed = constraints.get('allowedValues')
        if allowed is not None and literal not in allowed:
            findings.append(Finding(
                'error', 'constraint', path,
                f"value '{value}' is not allowed for parameter '{name}' "
                f"(allowed: {', '.join(allowed)})"))
            return

        minimum = constraints.get('minimum')
        maximum = constraints.get('maximum')
        if minimum is None and maximum is None:
            return
        number = _to_number(literal)
        if number is None:                       # not a numeric literal -> skip
            return
        if minimum is not None:
            bound = _to_number(minimum['value'])
            if bound is not None and (
                    number < bound if minimum['inclusive'] else number <= bound):
                relation = "less than" if minimum['inclusive'] else "less than or equal to"
                findings.append(Finding(
                    'error', 'constraint', path,
                    f"value '{value}' for parameter '{name}' must not be "
                    f"{relation} {minimum['value']}"))
        if maximum is not None:
            bound = _to_number(maximum['value'])
            if bound is not None and (
                    number > bound if maximum['inclusive'] else number >= bound):
                relation = "greater than" if maximum['inclusive'] else "greater than or equal to"
                findings.append(Finding(
                    'error', 'constraint', path,
                    f"value '{value}' for parameter '{name}' must not be "
                    f"{relation} {maximum['value']}"))

    walk_scope(root, None, 'parameters')
    return findings


# --- XInclude expansion ---------------------------------------------------
# Galacticus processes `<xi:include>` at run time, so the tree it actually sees
# is the included files merged in.  We expand includes before validating, so
# included content is checked too, and report an `include` finding when an href
# does not resolve.  Observed usage is uniformly
# `href="..." xpointer="xpointer(parameters/*)"` -- include the children of the
# referenced document's <parameters> root.

_XINCLUDE_NS      = 'http://www.w3.org/2001/XInclude'
_XI_INCLUDE       = f'{{{_XINCLUDE_NS}}}include'
_XI_FALLBACK      = f'{{{_XINCLUDE_NS}}}fallback'
_MAX_INCLUDE_DEPTH = 64
_XPOINTER_CHILDREN_RE = re.compile(r'^xpointer\(\s*/?(\w+)/\*\s*\)$')


def _xinclude_fallback(include_element):
    fallback = include_element.find(_XI_FALLBACK)
    return list(fallback) if fallback is not None else []


def _xpointer_select(sub_root, xpointer, label, findings):
    """Select the nodes a `<xi:include>` contributes, per its xpointer."""
    if xpointer is None:
        return [sub_root]                       # default: the whole referenced root
    match = _XPOINTER_CHILDREN_RE.match(xpointer.strip())
    if match:
        if sub_root.tag == match.group(1):
            return list(sub_root)               # children of e.g. <parameters>
        findings.append(Finding(
            'error', 'include', label,
            f"xi:include xpointer '{xpointer}' does not match referenced root "
            f"element '{sub_root.tag}'"))
        return []
    findings.append(Finding(
        'error', 'include', label,
        f"unsupported xi:include xpointer '{xpointer}'"))
    return []


def _locate_include(href, base_dir):
    """Locate an xi:include target, or None if it cannot be found.

    Tries the standard file-relative path first.  Some parameter files are
    *templates* executed from a different directory (e.g. test runners write a
    copy under testSuite/outputs/ and run it there), so their `../`-relative
    hrefs are relative to the run location, not the source.  As a fallback we
    resolve the path tail against the repository root (GALACTICUS_EXEC_PATH).
    This validates that the referenced content exists -- the actual goal --
    without false-failing on templates, while still catching missing files.
    """
    candidate = os.path.normpath(os.path.join(base_dir, href))
    if os.path.exists(candidate):
        return candidate
    exec_path = os.environ.get('GALACTICUS_EXEC_PATH')
    if exec_path:
        tail = re.sub(r'^(?:\.\./)+', '', href)
        fallback = os.path.normpath(os.path.join(exec_path, tail))
        if os.path.exists(fallback):
            return fallback
    return None


def _resolve_include(include_element, base_dir, findings, label, seen, depth):
    href = include_element.get('href')
    if not href:
        findings.append(Finding('error', 'include', label,
                                "xi:include without an href is not supported"))
        return _xinclude_fallback(include_element)
    target = _locate_include(href, base_dir)
    if target is None:
        findings.append(Finding(
            'error', 'include', label,
            f"xi:include href '{href}' does not resolve "
            f"(relative to {base_dir} or the repository root)"))
        return _xinclude_fallback(include_element)
    if target in seen:
        findings.append(Finding('error', 'include', label,
                                f"xi:include cycle detected at '{href}'"))
        return []
    try:
        sub_root = ET.parse(target).getroot()
    except ET.ParseError as exc:
        findings.append(Finding('error', 'include', label,
                                f"xi:include '{href}' failed to parse: {exc}"))
        return _xinclude_fallback(include_element)
    # Expand nested includes within the referenced file, relative to its own dir.
    _expand_xincludes(sub_root, os.path.dirname(target), findings,
                      os.path.basename(target), seen | {target}, depth + 1)
    return _xpointer_select(sub_root, include_element.get('xpointer'), label, findings)


def _expand_xincludes(element, base_dir, findings, label, seen, depth=0):
    """Recursively replace `<xi:include>` children of `element` in place."""
    if depth > _MAX_INCLUDE_DEPTH:
        findings.append(Finding('error', 'include', label,
                                "xi:include nesting too deep"))
        return
    new_children = []
    for child in list(element):
        if child.tag == _XI_INCLUDE:
            new_children.extend(
                _resolve_include(child, base_dir, findings, label, seen, depth))
        else:
            _expand_xincludes(child, base_dir, findings, label, seen, depth)
            new_children.append(child)
    element[:] = new_children


# --- Reference integrity (anchors) -----------------------------------------
# Galacticus resolves `idRef="X"` at run time to the element with `id="X"` and
# the SAME tag (resolveReferences); an unmatched idRef is a fatal run-time
# error, so we check the same on the XInclude-expanded tree.
#
# NOTE: expression/conditional `[path]` references are intentionally NOT checked.
# The bracket mini-language supports defaults (`[path|0.0]`, so an unresolved
# path is legal), printf-style formats (`[%4.4d|path]`), and relative paths --
# which makes path-resolution checking low-value (defaults legalise misses) and
# error-prone; left for a future expression-aware resolver if needed.


def _check_references(root, findings):
    """Check that every `idRef` anchor resolves to a matching `id`."""
    ids_by_tag = collections.defaultdict(set)
    for element in root.iter():
        if isinstance(element.tag, str) and element.get('id') is not None:
            ids_by_tag[element.tag].add(element.get('id'))

    for element in root.iter():
        if not isinstance(element.tag, str):
            continue
        id_ref = element.get('idRef')
        if id_ref is not None and id_ref not in ids_by_tag.get(element.tag, ()):
            findings.append(Finding(
                'error', 'reference', f"parameters/{element.tag}",
                f"idRef '{id_ref}' on <{element.tag}> has no matching element "
                f"with id='{id_ref}'"))


# --- Structural checks (catalog-independent, opt-in) -----------------------
# These subsume the legacy scripts/aux/validateParameters.py (format 2): a
# top-level parameter must not be duplicated, every parameter element must carry
# a value (attribute, idRef, or children), and an element must not have multiple
# <value> children.  They are OPT-IN (structural=True) because they assume a pure
# parameter file: run blanket over the whole test corpus they false-positive on
# non-parameter XML (merger-tree data files, changes files) and on parameter
# files with embedded tree data.  Intended for curated parameter files (as the
# legacy script was only ever run on parameters.xml + testSuite/.../validation/).

def _check_structure(root, findings):
    top_level = collections.Counter(
        child.tag for child in root if isinstance(child.tag, str))
    for tag, count in top_level.items():
        if count > 1 and tag not in _META_TAGS:
            findings.append(Finding(
                'error', 'duplicate', f"parameters/{tag}",
                f"parameter '{tag}' appears {count} times at the top level; "
                f"it should appear only once"))
    for element in root.iter():
        if (element is root or not isinstance(element.tag, str)
                or element.tag in _META_TAGS or element.tag == 'value'):
            continue                          # `<value>` holders are not parameters
        value_children = [c for c in element if c.tag == 'value']
        if ('value' not in element.attrib and 'idRef' not in element.attrib
                and len(element) == 0):
            findings.append(Finding(
                'error', 'value', f"parameters/{element.tag}",
                f"parameter '{element.tag}' has no value"))
        elif len(value_children) > 1:
            findings.append(Finding(
                'error', 'value', f"parameters/{element.tag}",
                f"parameter '{element.tag}' has multiple values"))


def validate_file(path, catalog=None, structural=False):
    """Parse, XInclude-expand, check references, and validate a parameter file;
    return ``(findings, error)`` where `error` is a parse-error string or None.

    `catalog=None` runs only the catalog-independent checks (XInclude resolution,
    references, and -- if requested -- structural), which is enough to subsume the
    legacy structural validator.  `structural=True` adds the structural checks (no
    duplicate top-level parameters, every parameter has a value, no multiple
    `<value>` children) -- only appropriate for curated, pure parameter files."""
    try:
        tree = ET.parse(path)
    except ET.ParseError as exc:
        return [], f"XML parse error: {exc}"
    root = tree.getroot()
    findings = []
    abspath = os.path.abspath(path)
    _expand_xincludes(root, os.path.dirname(abspath), findings,
                      os.path.basename(path), frozenset({abspath}))
    if structural:
        _check_structure(root, findings)
    _check_references(root, findings)
    if catalog is not None:
        findings.extend(validate_parameters(root, catalog))
    return findings, None
