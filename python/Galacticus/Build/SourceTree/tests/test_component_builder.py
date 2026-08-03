"""Tests for `process_component_builder` in
`Galacticus.Build.SourceTree.Process.ComponentBuilder`.

The hook replaces the retired `component` include-directive path (issue #32):
a single `<componentBuilder/>` directive in `source/objects/nodes/_class.F90`
is expanded, at preprocess time, into the synthesized node-component class
hierarchy, which is grafted in place of the directive.

The heavy `Galacticus.Build.Components` generator itself is exercised by
`Galacticus.Build.Components.tests`; here we test the hook's own logic:

* the directive node is replaced by the generated content (`graft`);
* module-level procedures that follow the directive stay after the generated
  `contains` (the generated content brings its own `contains`);
* a tree with no `componentBuilder` node does no work and imports nothing —
  the generator is only touched when the directive is actually present;
* (integration, skipped without a populated `$BUILDPATH`) the real generator
  produces the `nodeComponent` base type.
"""

import os

import pytest

import Galacticus.Build.SourceTree.Process.ComponentBuilder as CB
from Galacticus.Build.SourceTree import parse_code, serialize, walk_tree


_GENERATED_SNIPPET = (
    "  type, public :: testComponent\n"
    "    private\n"
    "  end type testComponent\n"
    "contains\n\n"
    "  subroutine testComponentBuild(self)\n"
    "    implicit none\n"
    "    class(testComponent), intent(inout) :: self\n"
    "    return\n"
    "  end subroutine testComponentBuild\n"
)


def _module_with_directive():
    return (
        "module Test_Nodes\n"
        "  implicit none\n"
        "  private\n"
        "  ! Build node methods.\n"
        "  !![\n"
        "  <componentBuilder/>\n"
        "  !!]\n"
        "  ! A hand-written procedure that must remain post-`contains`.\n"
        "  subroutine handWritten()\n"
        "    implicit none\n"
        "    return\n"
        "  end subroutine handWritten\n"
        "end module Test_Nodes\n"
    )


def test_graft_replaces_directive(monkeypatch):
    """The `<componentBuilder/>` node is replaced by the generated content."""
    monkeypatch.setattr(CB, '_build_component_content',
                        lambda: _GENERATED_SNIPPET)
    # Isolate the graft mechanics from the (separately-tested) full process
    # pipeline, which needs a populated `$BUILDPATH`; our fragment carries no
    # embedded directives, so processing it is a no-op here.
    monkeypatch.setattr(CB, 'process_tree', lambda tree, options=None: tree)

    tree = parse_code(_module_with_directive(), name='test')
    CB.process_component_builder(tree, {})

    out = serialize(tree)
    assert 'type, public :: testComponent' in out
    assert 'subroutine testComponentBuild' in out

    # The directive node is gone from the tree (replaced in place).
    assert not any(
        node.get('type') == 'componentBuilder' for node in walk_tree(tree)
    )


def test_hand_written_procedure_stays_after_generated_contains(monkeypatch):
    """The generated content carries its own module `contains`, so the
    hand-written procedure following the directive ends up after it — the
    layout the compiler `include` produced before the migration."""
    monkeypatch.setattr(CB, '_build_component_content',
                        lambda: _GENERATED_SNIPPET)
    monkeypatch.setattr(CB, 'process_tree', lambda tree, options=None: tree)

    tree = parse_code(_module_with_directive(), name='test')
    CB.process_component_builder(tree, {})
    out = serialize(tree)

    assert 'contains' in out
    assert out.index('contains') < out.index('subroutine handWritten')


def test_no_directive_does_no_work(monkeypatch):
    """A tree with no `componentBuilder` node must not invoke the generator
    (which is why the hook imports `Galacticus.Build.Components` lazily)."""
    def _boom():
        raise AssertionError("generator must not run without a directive")

    monkeypatch.setattr(CB, '_build_component_content', _boom)

    code = (
        "module Plain\n"
        "  implicit none\n"
        "end module Plain\n"
    )
    tree = parse_code(code, name='test')
    before = serialize(tree)
    CB.process_component_builder(tree, {})       # must not raise
    assert serialize(tree) == before


@pytest.mark.skipif(
    not (os.environ.get('BUILDPATH')
         and os.path.exists(os.path.join(
             os.environ.get('BUILDPATH', ''), 'directiveLocations.xml'))),
    reason="requires a populated $BUILDPATH (directiveLocations.xml)",
)
def test_real_generator_builds_node_component():
    """Integration: the real generator synthesizes the `nodeComponent` base
    type and a module-level `contains`."""
    content = CB._build_component_content()
    assert 'type, public :: nodeComponent' in content
    assert '\ncontains\n' in content
