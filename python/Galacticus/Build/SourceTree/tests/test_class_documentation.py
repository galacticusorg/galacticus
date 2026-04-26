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
