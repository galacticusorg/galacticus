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
    params = params_by_file.get("testFixture", [])
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
