"""Regression test for `_populate_class_descriptions` in
`Galacticus.Build.SourceTree.Process.ClassDocumentation`.

Bug: `<methods>` directives in parent classes are written as

    <method name="orbit">
      <description>...</description>
    </method>

while child classes use

    <method method="parametersSelect" description="..."/>

Our `xml_to_dict` keeps both forms verbatim — the parent's keyword
attribute lands in the `name` field, the child's in the `method` field.
Downstream code (`_resolve_method_bindings`,
`_compute_missing_and_generic`, plus the doc consumer
`scripts/doc/extractData.py`) only ever looked at `method:`, so the
parent's methods never made it into the `described_methods` set.  The
net effect: every parent's bound procedure (orbit,
densityContrastDefinition, …) appeared in the parent's
`<missingMethods>`, the consumer's inheritance walk found nothing in
the parent's descriptions, and the same methods stayed flagged as
missing in every derived class's `classes.xml`.

Fix: when populating `descriptions`, promote `name` to `method` for
every method dict that doesn't already carry it.
"""

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


# ---------------------------------------------------------------------------
# _resolve_method_bindings — `=>` lookup must be unanchored
# ---------------------------------------------------------------------------
#
# Bug: `_resolve_method_bindings` used `re.match` (anchored at the start
# of the string) to detect a `=>boundFunction` binding in
# `variables[0]`.  But after `parse_declaration` runs on
# `procedure :: parametersSelect => jiang2014ParametersSelect`,
# `variables[0]` is `parametersselect=>jiang2014parametersselect` —
# it starts with the method name, not with `=>`.  So `re.match` always
# returned None, the function never made it into the method's
# `boundFunctions` list, and `_process_function` couldn't attach a
# type.  The doc consumer then printed "missing function type" for
# EVERY type-bound method in EVERY child class.
#
# Fix: switch to `re.search` (unanchored, matching Perl's `m//`).

from Galacticus.Build.SourceTree.Process.ClassDocumentation import (
    _resolve_method_bindings,
)


def test_resolve_method_bindings_finds_arrow_anywhere():
    """A type-bound `procedure :: methodName => boundImpl` lands in the
    parsed `variables[0]` as `methodname=>boundimpl` — `_resolve_method_bindings`
    must find the `=>` even though it's not at the start of the string."""
    class_record = {
        'descriptions': [{'method': 'parametersSelect'}],
        'functions':    [[{
            'intrinsic':  'procedure',
            'attributes': [],
            'variables':  ['parametersselect=>jiang2014parametersselect'],
        }]],
    }
    _resolve_method_bindings(class_record)
    assert class_record['descriptions'][0].get('boundFunctions') \
        == ['jiang2014parametersselect']


def test_resolve_method_bindings_handles_multiple_bindings():
    """`procedure :: foo => a, bar => b` expands to a multi-element
    `variables` list; each `<name>=>` binding should be resolved
    independently for the matching method."""
    class_record = {
        'descriptions': [{'method': 'foo'}, {'method': 'bar'}],
        'functions':    [[
            {'intrinsic': 'procedure', 'attributes': [],
             'variables': ['foo=>fooimpl']},
            {'intrinsic': 'procedure', 'attributes': [],
             'variables': ['bar=>barimpl']},
        ]],
    }
    _resolve_method_bindings(class_record)
    foo_bound = class_record['descriptions'][0].get('boundFunctions')
    bar_bound = class_record['descriptions'][1].get('boundFunctions')
    assert foo_bound == ['fooimpl']
    assert bar_bound == ['barimpl']


def test_resolve_method_bindings_skips_method_with_no_arrow():
    """A `procedure :: name` (no `=>`, deferred-style without a binding)
    leaves `boundFunctions` unset — used by `_compute_missing_and_generic`
    to flag the method as missing."""
    class_record = {
        'descriptions': [{'method': 'foo'}],
        'functions':    [[
            {'intrinsic': 'procedure', 'attributes': [],
             'variables': ['foo']},
        ]],
    }
    _resolve_method_bindings(class_record)
    assert 'boundFunctions' not in class_record['descriptions'][0]


def test_class_documentation_skips_inner_process_tree_calls(tmp_path,
                                                              monkeypatch):
    """ClassDocumentation must skip when invoked from a nested
    `process_tree` call (Generics' clone reparse, FunctionClass'
    insert-and-write-output).  Without the skip, each generic-cloned
    subtree wrote a partial `<basename>.classes.xml` that was then
    overwritten by the next clone — the OUTER pass on the fully-
    expanded tree never got a chance to populate `<descriptions>` and
    `_process_function`-derived `type` fields side by side."""
    import Galacticus.Build.SourceTree.Process as proc_pkg
    from Galacticus.Build.SourceTree.Process.ClassDocumentation import (
        process_class_documentation,
    )

    monkeypatch.setenv('GALACTICUS_BUILD_DOCS', 'yes')
    monkeypatch.setenv('BUILDPATH', str(tmp_path))

    # Build the same `<functionClass>`-bearing tree as
    # `test_class_documentation_runs_at_outer_depth` — at outer depth
    # this writes `utility.foo.classes.xml`; at inner depth (the
    # generics / functionClass clone-reparse case) the function must
    # bail out before opening that file.
    fc_node = {
        'type':       'functionClass',
        'directive':  {
            'name':   'foo',
            'method': [{'name': 'bar', 'description': '…'}],
        },
        'firstChild': None, 'sibling': None, 'parent': None,
    }
    tree = {
        'name':       'utility.foo.F90',
        'firstChild': fc_node,
        'sibling':    None,
        'parent':     None,
        'type':       'file',
    }
    fc_node['parent'] = tree

    # Force the depth flag to "we are nested" — exactly what generics'
    # `_expand_subtree` would do via its own `process_tree(reparsed)`
    # call from within an outer `process_tree`.
    monkeypatch.setattr(proc_pkg, '_PROCESS_TREE_DEPTH', 2)

    process_class_documentation(tree, {})

    # No `<basename>.classes.xml` should have been emitted.
    assert list(tmp_path.iterdir()) == [], list(tmp_path.iterdir())


def test_class_documentation_runs_at_outer_depth(tmp_path, monkeypatch):
    """At the outer level (`_PROCESS_TREE_DEPTH == 1`) classDocumentation
    runs as normal."""
    import Galacticus.Build.SourceTree.Process as proc_pkg
    from Galacticus.Build.SourceTree.Process.ClassDocumentation import (
        process_class_documentation,
    )

    monkeypatch.setenv('GALACTICUS_BUILD_DOCS', 'yes')
    monkeypatch.setenv('BUILDPATH', str(tmp_path))

    # Build a minimal tree with a single `<functionClass>` directive so
    # there's something to write.
    fc_node = {
        'type':       'functionClass',
        'directive':  {
            'name':   'foo',
            'method': [{'name': 'bar', 'description': '…'}],
        },
        'firstChild': None, 'sibling': None, 'parent': None,
    }
    tree = {
        'name':       'utility.foo.F90',
        'firstChild': fc_node,
        'sibling':    None,
        'parent':     None,
        'type':       'file',
    }
    fc_node['parent'] = tree

    monkeypatch.setattr(proc_pkg, '_PROCESS_TREE_DEPTH', 1)
    process_class_documentation(tree, {})

    # The synthesised parent class record was written out.
    out = tmp_path / 'utility.foo.classes.xml'
    assert out.exists(), list(tmp_path.iterdir())
    content = out.read_text()
    assert '<fooClass>' in content


def test_class_documentation_skips_types_from_other_source_files(tmp_path,
                                                                   monkeypatch):
    """FunctionClass inserts each child class's `<type>` node into the
    parent's tree (via `insert_pre_contains` of the
    `code_content['module']['interfaces']` list).  The inserted nodes
    keep the CHILD source filename in their `source` field.  Walking
    them in the parent's outer postprocess pass would emit a
    `<virialOrbitJiang2014>` record under
    `satellites.merging.virial_orbits.classes.xml` — but the
    associated `jiang2014ParametersSelect` function lives in a
    submodule file, never gets walked here, so `_process_function`
    can't attach a type.  Each child file's own preprocess.py
    invocation already produces a complete record in
    `<child>.classes.xml`; let that own the documentation."""
    import Galacticus.Build.SourceTree.Process as proc_pkg
    from Galacticus.Build.SourceTree.Process.ClassDocumentation import (
        process_class_documentation,
    )

    monkeypatch.setenv('GALACTICUS_BUILD_DOCS', 'yes')
    monkeypatch.setenv('BUILDPATH', str(tmp_path))
    monkeypatch.setattr(proc_pkg, '_PROCESS_TREE_DEPTH', 1)
    # Reset the cross-invocation `_OUTPUT_PREVIOUS` accumulator to avoid
    # state leaking from earlier tests in the same pytest session.
    from Galacticus.Build.SourceTree.Process import ClassDocumentation as cd
    monkeypatch.setattr(cd, '_OUTPUT_PREVIOUS', [])

    # A type node whose `source` is a different real source file —
    # exactly the shape FunctionClass produces when it parses each child
    # class file and threads its `<type>` node into the parent tree.
    inserted_child_type = {
        'type':       'type',
        'name':       'virialOrbitJiang2014',
        'opener':     'type, extends(virialOrbitClass) :: virialOrbitJiang2014',
        'source':     'satellites.merging.virial_orbits.Jiang2014.F90',
        'firstChild': None, 'sibling': None, 'parent': None,
    }
    tree = {
        'name':       'satellites.merging.virial_orbits.F90',
        'source':     'satellites.merging.virial_orbits.F90',
        'firstChild': inserted_child_type,
        'sibling':    None,
        'parent':     None,
        'type':       'file',
    }
    inserted_child_type['parent'] = tree

    process_class_documentation(tree, {})

    # The classes.xml file should NOT contain a <virialOrbitJiang2014>
    # entry — that ownership belongs to the child file's own pass.
    out = tmp_path / 'satellites.merging.virial_orbits.classes.xml'
    if out.exists():
        content = out.read_text()
        assert 'virialOrbitJiang2014' not in content, content


def test_class_documentation_keeps_native_types(tmp_path, monkeypatch):
    """A type whose `source` matches the outer tree IS processed."""
    import Galacticus.Build.SourceTree.Process as proc_pkg
    from Galacticus.Build.SourceTree.Process.ClassDocumentation import (
        process_class_documentation,
    )

    monkeypatch.setenv('GALACTICUS_BUILD_DOCS', 'yes')
    monkeypatch.setenv('BUILDPATH', str(tmp_path))
    monkeypatch.setattr(proc_pkg, '_PROCESS_TREE_DEPTH', 1)
    # Reset the cross-invocation `_OUTPUT_PREVIOUS` accumulator to avoid
    # state leaking from earlier tests in the same pytest session.
    from Galacticus.Build.SourceTree.Process import ClassDocumentation as cd
    monkeypatch.setattr(cd, '_OUTPUT_PREVIOUS', [])

    native_type = {
        'type':       'type',
        'name':       'nativeType',
        'opener':     'type :: nativeType',
        'source':     'utility.foo.F90',
        'firstChild': None, 'sibling': None, 'parent': None,
    }
    tree = {
        'name':       'utility.foo.F90',
        'source':     'utility.foo.F90',
        'firstChild': native_type,
        'sibling':    None,
        'parent':     None,
        'type':       'file',
    }
    native_type['parent'] = tree

    process_class_documentation(tree, {})

    out = tmp_path / 'utility.foo.classes.xml'
    assert out.exists()
    assert '<nativeType>' in out.read_text()


def test_class_documentation_keeps_synthetic_source_types(tmp_path,
                                                           monkeypatch):
    """The auto-generated `<base>Class` type that FunctionClass injects
    has source = 'Galacticus.Build.SourceTree.Process.FunctionClass.…'
    — those must STILL be processed, since their function bodies are
    in the same outer tree (FunctionClass also inserts `tree_post`)."""
    import Galacticus.Build.SourceTree.Process as proc_pkg
    from Galacticus.Build.SourceTree.Process.ClassDocumentation import (
        process_class_documentation,
    )

    monkeypatch.setenv('GALACTICUS_BUILD_DOCS', 'yes')
    monkeypatch.setenv('BUILDPATH', str(tmp_path))
    monkeypatch.setattr(proc_pkg, '_PROCESS_TREE_DEPTH', 1)
    # Reset the cross-invocation `_OUTPUT_PREVIOUS` accumulator to avoid
    # state leaking from earlier tests in the same pytest session.
    from Galacticus.Build.SourceTree.Process import ClassDocumentation as cd
    monkeypatch.setattr(cd, '_OUTPUT_PREVIOUS', [])

    autogen_type = {
        'type':       'type',
        'name':       'fooClass',
        'opener':     'type, abstract :: fooClass',
        'source':     'Galacticus.Build.SourceTree.Process.FunctionClass'
                      '.process_function_class()',
        'firstChild': None, 'sibling': None, 'parent': None,
    }
    tree = {
        'name':       'foo.F90',
        'source':     'foo.F90',
        'firstChild': autogen_type,
        'sibling':    None,
        'parent':     None,
        'type':       'file',
    }
    autogen_type['parent'] = tree

    process_class_documentation(tree, {})

    out = tmp_path / 'foo.classes.xml'
    assert out.exists()
    assert '<fooClass>' in out.read_text()


def test_synthesised_class_extends_functionClass():
    """The abstract base record synthesised from a `<functionClass>`
    directive must advertise `extends: 'functionClass'` so the doc
    consumer's `isFunctionClass` chain walk reaches `functionClass`
    when starting from any concrete derived class.  Without this,
    `isFunctionClass` stays False on every descendant and the
    `FUNCTION_CLASS_EXCLUDES` filter at the consumer side never
    triggers — every concrete child class then warns "missing method
    descriptions" for every auto-generated method (autoHook,
    descriptor, deepCopy, …)."""
    directive_node = {
        'type':      'functionClass',
        'directive': {
            'name':   'foo',
            'method': [{'name': 'bar', 'description': '…'}],
        },
    }
    classes = {}
    _populate_class_from_function_class_directive(directive_node, classes)
    assert classes['fooClass'].get('extends') == 'functionClass'


def test_directive_method_with_void_type_translated_to_subroutine():
    """FunctionClass auto-generated method stubs (autoHook, descriptor, …)
    arrive with `<type>void</type>`.  The doc consumer's
    `declaration_builder` only recognises the literal string `'subroutine'`
    for the void-return case, so the synthesiser must translate.

    Without this translation, `parse_declaration` would fail to parse
    `void` as an intrinsic and we'd leave the raw string `'void'` in
    place — `extractData.py` then crashed with
    `ValueError: declaration_builder: unknown type 'void'`."""
    directive_node = {
        'type':      'functionClass',
        'directive': {
            'name': 'myThing',
            'method': [{'name': 'autoHook', 'type': 'void'}],
        },
    }
    classes = {}
    _populate_class_from_function_class_directive(directive_node, classes)
    method_type = classes['myThingClass']['descriptions'][0]['type']
    assert method_type == 'subroutine'


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
