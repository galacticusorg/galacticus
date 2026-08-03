"""Regression tests for ``scripts/doc/extractDocsRST.py``.

Guards the multi-directive-per-block extraction.  A single ``!![ … !!]`` block
routinely declares several ``<inputParameter>`` directives (and, in one case,
two ``<workaround>`` directives); every one must be emitted, not just the
block's first.  A regression here silently dropped ~68% of the parameter docs
(and one compiler workaround) from the rendered ReadTheDocs pages.
"""
import os
import sys

# extractDocsRST self-inserts its own directory on sys.path for its sibling
# imports (latexToRST, extractContributors); do the same so we can import it.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))
import extractDocsRST  # noqa: E402


# One file with a two-parameter block (plus an unrelated <objectBuilder> to
# mimic real source) and a two-workaround block.
_FIXTURE = """\
module Test_Fixture
  !![
  <inputParameter docformat="rst">
    <name>epsilon</name>
    <defaultValue>3.0d0</defaultValue>
    <description>First parameter.</description>
    <source>parameters</source>
  </inputParameter>
  <inputParameter docformat="rst">
    <name>gamma</name>
    <defaultValue>2.5d0</defaultValue>
    <description>Second parameter.</description>
    <source>parameters</source>
  </inputParameter>
  <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
  !!]

  !![
  <workaround type="compiler" PR="110547" url="https://example/110547" docformat="rst">
    <description>First workaround.</description>
  </workaround>
  <workaround type="compiler" PR="110548" url="https://example/110548" docformat="rst">
    <description>Second workaround.</description>
  </workaround>
  !!]
end module Test_Fixture
"""


def _scan(tmp_path):
    (tmp_path / "testFixture.F90").write_text(_FIXTURE)
    return extractDocsRST.scan_source(str(tmp_path))


def test_all_input_parameters_in_a_block_are_extracted(tmp_path):
    _families, _impl, params_by_file, *_rest = _scan(tmp_path)
    # `scan_source` keys this by the file's full path rather than its basename: every
    # pluggable physics class lives in a file named `_class.F90`, so basenames collide
    # across all of them.
    params = params_by_file.get(str(tmp_path / "testFixture.F90"), [])
    assert [p["name"] for p in params] == ["epsilon", "gamma"], (
        "every <inputParameter> in a block must be emitted, not just the first"
    )
    # The non-first parameter must carry its own default/description.
    gamma = next(p for p in params if p["name"] == "gamma")
    assert gamma["default"] == "2.5d0"
    assert gamma["description"].strip() == "Second parameter."


def test_all_workarounds_in_a_block_are_extracted(tmp_path):
    *_head, workarounds = _scan(tmp_path)
    assert sorted(w["pr"] for w in workarounds) == ["110547", "110548"], (
        "every <workaround> in a block must be emitted, not just the first"
    )


# ---------------------------------------------------------------------------
# Method signatures
#
# A `<methods>` block declares a name and a description; the return type and
# arguments are recovered from the procedure each method is bound to, so that a
# class page documents what a method takes and returns rather than just naming
# it.  Blocks that declare a signature themselves (the ones the code generators
# emit) are rendered as declared.
# ---------------------------------------------------------------------------

_METHODS_FIXTURE = """\
module Test_Methods
  !![
  <functionClass docformat="rst">
   <name>testClass</name>
   <descriptiveName>Test Class</descriptiveName>
   <description>A test class.</description>
   <default>example</default>
   <method name="radiusStripped" >
    <description>Return the stripping radius.</description>
    <type>double precision</type>
    <argument>type(treeNode), intent(inout), target :: node</argument>
   </method>
   <method name="reset" >
    <description>Reset the object.</description>
    <type>void</type>
   </method>
  </functionClass>
  !!]

  !![
  <testClass name="testClassExample" docformat="rst">
   <description>An example implementation.</description>
  </testClass>
  !!]

  type :: testClassExample
   contains
     !![
     <methods docformat="rst">
       <method description="Return the timescale." method="timescale"/>
       <method description="Tabulate the timescale." method="tabulate"/>
       <method description="Document a method that no longer exists." method="departed"/>
       <method method="generated" description="A generated method.">
        <type>``double precision``</type>
        <argument>``integer count`` (optional) [in]</argument>
       </method>
     </methods>
     !!]
     procedure :: timescale => exampleTimescale
     procedure :: tabulate  => exampleTabulate
  end type testClassExample

contains

  double precision function exampleTimescale(self,node,factor)
    implicit none
    class           (testClassExample), intent(inout)           :: self
    type            (treeNode        ), intent(inout), target   :: node
    double precision                  , intent(in   ), optional :: factor

    exampleTimescale=0.0d0
    return
  end function exampleTimescale

  subroutine exampleTabulate(self)
    implicit none
    class(testClassExample), intent(inout) :: self

    return
  end subroutine exampleTabulate
end module Test_Methods
"""


def _scan_methods(tmp_path):
    (tmp_path / "testMethods.F90").write_text(_METHODS_FIXTURE)
    return extractDocsRST.scan_source(str(tmp_path))


def _type_methods(tmp_path):
    _families, _impl, _params, methods_by_file, *_rest = _scan_methods(tmp_path)
    return {m["name"]: m
            for m in methods_by_file[str(tmp_path / "testMethods.F90")]}


def test_type_method_signature_is_recovered_from_its_binding(tmp_path):
    methods = _type_methods(tmp_path)
    assert methods["timescale"]["type"] == "``double precision``"
    assert methods["timescale"]["arguments"] == [
        "``type(treeNode) node`` [inout]",
        "``double precision factor`` (optional) [in]",
    ], "the passed object is not an argument; optional arguments are marked"


def test_type_method_bound_to_a_subroutine_has_no_signature(tmp_path):
    methods = _type_methods(tmp_path)
    assert methods["tabulate"]["type"]      is None
    assert methods["tabulate"]["arguments"] == []


def test_type_method_without_a_binding_keeps_its_description(tmp_path):
    """A documented method that no longer exists must not break the page."""
    methods = _type_methods(tmp_path)
    assert methods["departed"]["description"] == \
        "Document a method that no longer exists."
    assert not methods["departed"]["arguments"]


def test_declared_signature_is_used_as_declared(tmp_path):
    """Where a block declares the signature — as the generators' blocks do — it
    is taken as given rather than recovered."""
    methods = _type_methods(tmp_path)
    assert methods["generated"]["type"]      == "``double precision``"
    assert methods["generated"]["arguments"] == \
        ["``integer count`` (optional) [in]"]


def test_rendered_page_shows_signatures(tmp_path):
    families, implementations, params_by_file, methods_by_file, *_rest = \
        _scan_methods(tmp_path)
    page = extractDocsRST.render_family(
        "testClass", families, implementations, params_by_file,
        methods_by_file, {})
    # The interface method's raw declaration is rendered in the notation, not
    # verbatim.
    assert "``radiusStripped`` → ``double precision``" in page
    assert "* ``type(treeNode) node`` [inout]"         in page
    assert "intent(inout), target :: node"         not in page
    # `void` names no return type.
    assert "``reset``\n"     in page
    assert "``reset`` → " not in page
    # An implementation's own methods read the same way.
    assert "``timescale`` → ``double precision``"      in page
    assert "* ``double precision factor`` (optional) [in]" in page
