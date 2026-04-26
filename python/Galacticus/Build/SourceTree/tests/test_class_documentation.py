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


# ---------------------------------------------------------------------------
# Method type / argument normalisation
# ---------------------------------------------------------------------------
#
# Bug: the synthesised parent class records produced
# `<descriptions><type>double precision</type>…</descriptions>` — i.e. raw
# Fortran-source strings carried straight through from the directive XML.
# `scripts/doc/extractData.py`'s `declaration_builder` expects either a
# parsed declaration dict (`{intrinsic, type, attributes, …}`) or the
# literal string `'subroutine'`, and crashed on the bare `'double precision'`
# with `ValueError: declaration_builder: unknown type 'double precision'`.
#
# Fix: parse `<type>` / `<argument>` strings into declaration dicts at
# synthesis time so the consumer sees the same shape it gets for type-
# derived methods.


def test_directive_method_type_string_parsed_to_dict():
    """A `<type>type(keplerOrbit)</type>` child becomes
    `{intrinsic: 'type', type: 'keplerOrbit', …}` after synthesis."""
    directive_node = {
        'type':      'functionClass',
        'directive': {
            'name': 'virialOrbit',
            'method': [
                {'name': 'orbit', 'type': 'type(keplerOrbit)',
                 'description': '…'},
            ],
        },
    }
    classes = {}
    _populate_class_from_function_class_directive(directive_node, classes)
    method_type = classes['virialOrbitClass']['descriptions'][0]['type']
    assert isinstance(method_type, dict), method_type
    assert method_type.get('intrinsic') == 'type'
    assert method_type.get('type')      == 'keplerOrbit'


def test_directive_method_double_precision_type_parsed():
    """`double precision` is two-word intrinsic with no kind/type spec."""
    directive_node = {
        'type':      'functionClass',
        'directive': {
            'name': 'orbit',
            'method': [
                {'name': 'energyMean', 'type': 'double precision',
                 'description': '…'},
            ],
        },
    }
    classes = {}
    _populate_class_from_function_class_directive(directive_node, classes)
    method_type = classes['orbitClass']['descriptions'][0]['type']
    assert isinstance(method_type, dict)
    assert method_type.get('intrinsic') == 'double precision'


def test_directive_method_with_no_type_becomes_subroutine():
    """A method with no `<type>` is subroutine-like — the consumer's
    `declaration_builder` handles the literal `'subroutine'` string."""
    directive_node = {
        'type':      'functionClass',
        'directive': {
            'name': 'foo',
            'method': [{'name': 'doSomething', 'description': '…'}],
        },
    }
    classes = {}
    _populate_class_from_function_class_directive(directive_node, classes)
    method_type = classes['fooClass']['descriptions'][0]['type']
    assert method_type == 'subroutine'


def test_directive_method_arguments_parsed_to_declarations():
    """`<argument>type(treeNode), intent(inout) :: node, host</argument>`
    becomes `{argument: [{intrinsic: 'type', type: 'treeNode', ...,
    variableNames: ['node', 'host']}]}`, with a parallel
    `argumentList=['node, host']`."""
    directive_node = {
        'type':      'functionClass',
        'directive': {
            'name': 'orbit',
            'method': [{
                'name':     'orbit',
                'type':     'type(keplerOrbit)',
                'argument': [
                    'type(treeNode), intent(inout) :: node, host',
                    'logical, intent(in) :: acceptUnboundOrbits',
                ],
            }],
        },
    }
    classes = {}
    _populate_class_from_function_class_directive(directive_node, classes)
    method = classes['orbitClass']['descriptions'][0]

    arg_groups = method.get('arguments') or []
    assert len(arg_groups) == 1
    parsed_args = arg_groups[0]['argument']
    assert len(parsed_args) == 2
    assert parsed_args[0]['intrinsic']     == 'type'
    assert parsed_args[0]['type']          == 'treeNode'
    assert parsed_args[0]['variableNames'] == ['node', 'host']
    assert parsed_args[1]['intrinsic']     == 'logical'
    assert parsed_args[1]['variableNames'] == ['acceptUnboundOrbits']

    # And the parallel name list joins all parsed names.
    arg_lists = method.get('argumentList') or []
    assert arg_lists == ['node, host, acceptUnboundOrbits']


def test_directive_method_single_argument_string_handled():
    """`as_array` handles both the list and the scalar-string forms — the
    single-argument case must produce the same output shape."""
    directive_node = {
        'type':      'functionClass',
        'directive': {
            'name': 'foo',
            'method': [{
                'name':     'doIt',
                'argument': 'integer, intent(in) :: x',
            }],
        },
    }
    classes = {}
    _populate_class_from_function_class_directive(directive_node, classes)
    method = classes['fooClass']['descriptions'][0]
    parsed = method['arguments'][0]['argument'][0]
    assert parsed['intrinsic']     == 'integer'
    assert parsed['variableNames'] == ['x']


def test_synthesis_does_not_mutate_directive_method_dicts():
    """ClassDocumentation must NOT mutate the original method dicts in the
    `<functionClass>` directive — FunctionClass reads them later and
    expects `method['type']` to remain the raw Fortran-source string
    (`type(keplerOrbit)`, etc.) so it can do its own `.startswith(…)`
    checks.  An earlier draft that mutated in place crashed FunctionClass
    with `AttributeError: 'dict' object has no attribute 'startswith'`."""
    raw_method = {'name': 'orbit', 'type': 'type(keplerOrbit)'}
    directive_node = {
        'type':      'functionClass',
        'directive': {
            'name':   'virialOrbit',
            'method': [raw_method],
        },
    }
    classes = {}
    _populate_class_from_function_class_directive(directive_node, classes)

    # Synthesised class record gets the parsed-dict form.
    synthesised_type = classes['virialOrbitClass']['descriptions'][0]['type']
    assert isinstance(synthesised_type, dict)

    # Original directive's method dict still has its raw string type.
    assert raw_method['type'] == 'type(keplerOrbit)'
    assert 'method' not in raw_method   # name→method promotion didn't leak
