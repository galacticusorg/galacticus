# Regression test for `_populate_class_descriptions` in
# `Galacticus.Build.SourceTree.Process.ClassDocumentation`.
#
# Bug: `<methods>` directives in parent classes are written as
#
#     <method name="orbit">
#       <description>...</description>
#     </method>
#
# while child classes use
#
#     <method method="parametersSelect" description="..."/>
#
# Our `xml_to_dict` keeps both forms verbatim — the parent's keyword
# attribute lands in the `name` field, the child's in the `method` field.
# Downstream code (`_resolve_method_bindings`,
# `_compute_missing_and_generic`, plus the doc consumer
# `scripts/doc/extractData.py`) only ever looked at `method:`, so the
# parent's methods never made it into the `described_methods` set.  The
# net effect: every parent's bound procedure (orbit,
# densityContrastDefinition, …) appeared in the parent's
# `<missingMethods>`, the consumer's inheritance walk found nothing in
# the parent's descriptions, and the same methods stayed flagged as
# missing in every derived class's `classes.xml`.
#
# Fix: when populating `descriptions`, promote `name` to `method` for
# every method dict that doesn't already carry it.

from Galacticus.Build.SourceTree.Process.ClassDocumentation import (
    _populate_class_descriptions,
)


def _type_node_with_methods(directive_method_payload):
    """Build a minimal `type` node with a single `<methods>` directive in
    its contains area.  `_populate_class_descriptions` walks the
    contains-area children, so we attach a synthetic `contains` marker
    sibling of the `<methods>` node."""
    contains = {
        'type':       'contains',
        'firstChild': None,
        'sibling':    None,
        'parent':     None,
    }
    methods_node = {
        'type':       'methods',
        'directive':  {'method': directive_method_payload},
        'firstChild': None,
        'sibling':    None,
        'parent':     None,
    }
    contains['sibling'] = methods_node
    type_node = {
        'type':       'type',
        'name':       'someClass',
        'firstChild': contains,
        'sibling':    None,
        'parent':     None,
    }
    return type_node


def test_method_with_name_attribute_normalised_to_method_key():
    """`<method name="orbit">` is normalised so the resulting dict has
    `method: 'orbit'` (in addition to `name: 'orbit'`)."""
    type_node = _type_node_with_methods([
        {'name': 'orbit',                       'description': 'Return orbit.'},
        {'name': 'velocityDistributionFunction', 'description': 'Distribution.'},
    ])
    cls = {}
    _populate_class_descriptions(type_node, cls)
    methods = [(d.get('method'), d.get('description'))
               for d in cls.get('descriptions') or []]
    assert methods == [
        ('orbit',                       'Return orbit.'),
        ('velocityDistributionFunction', 'Distribution.'),
    ], methods


def test_method_with_method_attribute_passes_through():
    """`<method method="parametersSelect" …>` is already in the right shape."""
    type_node = _type_node_with_methods([
        {'method': 'parametersSelect', 'description': 'Pick a parameter set.'},
    ])
    cls = {}
    _populate_class_descriptions(type_node, cls)
    methods = [(d.get('method'), d.get('description'))
               for d in cls.get('descriptions') or []]
    assert methods == [('parametersSelect', 'Pick a parameter set.')]


def test_single_method_dict_form():
    """`xml_to_dict` returns a single dict (not a list) when there's only
    one `<method>` child — `as_array` should normalise that, and the
    promotion still applies."""
    type_node = _type_node_with_methods(
        {'name': 'lonelyMethod', 'description': 'Just one.'})
    cls = {}
    _populate_class_descriptions(type_node, cls)
    methods = [(d.get('method'), d.get('description'))
               for d in cls.get('descriptions') or []]
    assert methods == [('lonelyMethod', 'Just one.')]


def test_method_with_both_attributes_keeps_method():
    """If a directive author writes BOTH `method` and `name` (improbable
    but defined behaviour), `method` wins — we don't clobber it."""
    type_node = _type_node_with_methods([
        {'method': 'foo', 'name': 'bar', 'description': '…'},
    ])
    cls = {}
    _populate_class_descriptions(type_node, cls)
    method = (cls['descriptions'][0]).get('method')
    assert method == 'foo', cls['descriptions']


# ---------------------------------------------------------------------------
# functionClass directive synthesis
# ---------------------------------------------------------------------------
#
# Bug: ClassDocumentation runs BEFORE FunctionClass in the topo order, so
# the abstract `<name>Class` type that FunctionClass auto-generates from
# a `<functionClass>` directive is never present in the tree at the time
# ClassDocumentation walks it.  The parent's methods (declared inline as
# `<method name="…">` children of the `<functionClass>` directive) never
# reach any `<descriptions>` block, so the doc consumer's inheritance
# walk in `scripts/doc/extractData.py` cannot resolve them — and EVERY
# derived class winds up listing them in `<missingMethods>`.
#
# Fix: harvest the directive directly when classDocumentation walks the
# tree.

from Galacticus.Build.SourceTree.Process.ClassDocumentation import (
    _populate_class_from_function_class_directive,
)


def test_function_class_directive_creates_parent_class_record():
    """A `<functionClass><name>virialOrbit</name>...</functionClass>` directive
    becomes a `virialOrbitClass` entry whose `descriptions` lists every
    `<method name="…">` child."""
    directive_node = {
        'type':      'functionClass',
        'directive': {
            'name': 'virialOrbit',
            'method': [
                {'name': 'orbit',                       'description': 'Return orbit.'},
                {'name': 'velocityDistributionFunction', 'description': 'Distribution.'},
                {'name': 'energyMean',                  'description': 'Mean energy.'},
            ],
        },
    }
    classes = {}
    _populate_class_from_function_class_directive(directive_node, classes)

    assert 'virialOrbitClass' in classes
    record = classes['virialOrbitClass']
    assert record['name'] == 'virialOrbitClass'
    method_names = [d.get('method') for d in record.get('descriptions') or []]
    assert method_names == ['orbit', 'velocityDistributionFunction', 'energyMean']


def test_function_class_directive_with_single_method():
    """A `<functionClass>` with a single `<method>` child arrives as a dict
    rather than a list — both the directive synthesis AND the
    `name`-to-`method` promotion must work in that case too."""
    directive_node = {
        'type':      'functionClass',
        'directive': {
            'name': 'oneMethod',
            'method': {'name': 'soloMethod', 'description': 'Just one.'},
        },
    }
    classes = {}
    _populate_class_from_function_class_directive(directive_node, classes)

    record = classes.get('oneMethodClass')
    assert record is not None
    assert [d.get('method') for d in record['descriptions']] == ['soloMethod']


def test_function_class_directive_without_name_is_ignored():
    """A malformed directive lacking `<name>` is silently skipped (the
    rest of the pipeline already errors on this — we don't crash on it)."""
    classes = {}
    _populate_class_from_function_class_directive(
        {'type': 'functionClass', 'directive': {}}, classes)
    assert classes == {}


def test_function_class_directive_method_dict_already_has_method_key():
    """If a `<functionClass>` directive author wrote
    `<method method="x" description="…"/>` (the child-class form, rare for
    parents but valid), that `method` key is honoured as-is."""
    directive_node = {
        'type':      'functionClass',
        'directive': {
            'name': 'mixed',
            'method': [
                {'method': 'preNamedX', 'description': '…'},
                {'name':   'fromName',  'description': '…'},
            ],
        },
    }
    classes = {}
    _populate_class_from_function_class_directive(directive_node, classes)
    method_names = [d.get('method') for d in classes['mixedClass']['descriptions']]
    assert method_names == ['preNamedX', 'fromName']
