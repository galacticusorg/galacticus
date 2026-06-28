"""Tests for `Galacticus.Parameters.catalog`.

Covers the two pieces most likely to drift from the rest of the build system:

  1. `derive_label` -- must reproduce the canonical `lcfirst-unless-all-caps`
     selector convention from `FunctionClass/__init__.py` (so `NFW` keeps its
     casing while `Simple` becomes `simple`); a mismatch silently makes catalog
     labels disagree with the strings users write in parameter files.
  2. `harvest_file` end-to-end on a synthetic implementation -- exercises
     registration discovery, constructor location, type inference across all
     three provenance paths, and the `source`->nested-element resolution.
"""

import textwrap

import pytest

from Galacticus.Parameters.catalog import (
    derive_label, _resolve_source_elements, harvest_file,
    discover_enumerations, _enumeration_links, _capture_constraints,
    _normalize_default,
)
from Galacticus.Build.SourceTree import parse_file


@pytest.mark.parametrize("base, implementation_type, expected_label", [
    ('accretionHalo',       'accretionHaloSimple',          'simple'),
    ('cosmologyParameters', 'cosmologyParametersSimple',    'simple'),
    ('darkMatterProfileDMO', 'darkMatterProfileDMONFW',     'NFW'),          # all-caps acronym kept
    ('darkMatterProfileDMO', 'darkMatterProfileDMOSIDMCoreNFW', 'SIDMCoreNFW'),
    ('haloMassFunction',    'haloMassFunctionTinker2008',   'tinker2008'),   # single cap -> lcfirst
])
def test_derive_label(base, implementation_type, expected_label):
    assert derive_label(base, implementation_type) == expected_label


def test_resolve_source_elements():
    """`<var> = parameters%subParameters('elem')` binds a local source handle to
    a nested element; the resolver recovers `{var: elem}`."""
    tree = parse_file_from_string("""\
        subParameters=parameters%subParameters('componentHotHalo')
        parametersMassDefinitions=parameters%subParameters('massDefinitions',requireValue=.false.)
    """)
    mapping = _resolve_source_elements(tree)
    assert mapping['subParameters']             == 'componentHotHalo'
    assert mapping['parametersMassDefinitions'] == 'massDefinitions'


_SNIPPET = textwrap.dedent("""\
    module My_Module
      private
      !![
      <myClass name="myClassSimple">
       <description>A synthetic implementation for testing.</description>
      </myClass>
      !!]
      type, extends(myClassClass) :: myClassSimple
         double precision :: alpha
         logical          :: betaFlag
      end type myClassSimple
    contains

      function simpleConstructorParameters(parameters) result(self)
        type            (myClassSimple   ) :: self
        type            (inputParameters ) :: parameters
        type            (inputParameters ) :: subParameters
        double precision                   :: gamma
        integer                            :: delta
        !![
        <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
        <inputParameter>
          <name>alpha</name>
          <defaultValue>1.0d0</defaultValue>
          <source>parameters</source>
          <description>alpha.</description>
        </inputParameter>
        <inputParameter>
          <name>betaFlag</name>
          <defaultValue>.true.</defaultValue>
          <source>parameters</source>
          <description>beta.</description>
        </inputParameter>
        <inputParameter>
          <name>gamma</name>
          <source>parameters</source>
          <description>gamma (no default; typed from declaration).</description>
        </inputParameter>
        !!]
        subParameters=parameters%subParameters('nestedElem')
        !![
        <inputParameter>
          <name>delta</name>
          <source>subParameters</source>
          <description>delta (nested source).</description>
        </inputParameter>
        !!]
        return
      end function simpleConstructorParameters
    end module My_Module
    """)


def parse_file_from_string(code):
    """Write `code` to a temp file and parse it via the real SourceTree parser."""
    import tempfile
    import os
    fd, path = tempfile.mkstemp(suffix='.F90')
    try:
        with os.fdopen(fd, 'w') as fh:
            fh.write(code)
        return parse_file(path)
    finally:
        os.unlink(path)


def test_discover_enumerations(tmp_path):
    src = tmp_path / "source" / "x"
    src.mkdir(parents=True)
    (src / "e.F90").write_text(
        'module M\n  !![\n  <enumeration>\n   <name>fixedDensityType</name>\n'
        '   <entry label="critical"/>\n   <entry label="mean"/>\n  </enumeration>\n'
        '  !!]\nend module M\n')
    enums = discover_enumerations([str(src / "e.F90")])
    assert enums == {'fixedDensityType': ['critical', 'mean']}


def test_enumeration_links_from_encode_call():
    tree = parse_file_from_string("""\
        densityType=enumerationFixedDensityTypeEncode(char(densityType),includesPrefix=.false.)
        plain=somethingElse(char(other))
    """)
    links = _enumeration_links(tree)
    assert links.get('densityType') == 'fixedDensityType'
    assert 'other' not in links


_SNIPPET_ENUM = textwrap.dedent("""\
    module My_Module
      !![
      <myClass name="myClassEnum"><description>d.</description></myClass>
      !!]
      type, extends(myClassClass) :: myClassEnum
      end type myClassEnum
    contains
      function enumConstructorParameters(parameters) result(self)
        type(myClassEnum    ) :: self
        type(inputParameters) :: parameters
        type(varying_string ) :: densityType
        !![
        <inputParameter><name>densityType</name><source>parameters</source>
          <description>d.</description></inputParameter>
        !!]
        self=myClassEnum(enumerationFixedDensityTypeEncode(char(densityType),includesPrefix=.false.))
        return
      end function enumConstructorParameters
    end module My_Module
    """)


def test_harvest_links_enumeration(tmp_path):
    source_root = tmp_path / "source"
    (source_root / "my").mkdir(parents=True)
    impl = source_root / "my" / "e.F90"
    impl.write_text(_SNIPPET_ENUM)
    entries = harvest_file(str(impl), {"myClass"}, str(source_root))
    params = {p['name']: p for p in entries[0]['parameters']}
    assert params['densityType']['enumeration'] == 'fixedDensityType'


def test_normalize_default():
    assert _normalize_default("var_str('critical')") == 'critical'
    assert _normalize_default('var_str("none")') == 'none'
    assert _normalize_default("var_str( 'x' )") == 'x'
    # Non-var_str defaults pass through unchanged.
    assert _normalize_default('9.97d0') == '9.97d0'
    assert _normalize_default('.true.') == '.true.'
    assert _normalize_default("inputPath(pathTypeDataStatic)//'f.hdf5'") \
        == "inputPath(pathTypeDataStatic)//'f.hdf5'"
    assert _normalize_default(None) is None


def test_capture_constraints():
    # Bare text -> inclusive by default; an `inclusive` attribute is honoured.
    assert _capture_constraints({
        'minimum': '0.0',
        'maximum': {'inclusive': 'false', 'content': '1.0'},
    }) == {
        'minimum': {'value': '0.0', 'inclusive': True},
        'maximum': {'value': '1.0', 'inclusive': False},
    }
    # allowedValues splits on whitespace.
    assert _capture_constraints({'allowedValues': 'fast slow'}) \
        == {'allowedValues': ['fast', 'slow']}
    # No constraint elements -> None.
    assert _capture_constraints({'name': 'x'}) is None


def test_harvest_captures_constraints(tmp_path):
    source_root = tmp_path / "source"
    (source_root / "my").mkdir(parents=True)
    impl = source_root / "my" / "c.F90"
    impl.write_text(textwrap.dedent("""\
        module M
          !![
          <myClass name="myClassBound"><description>d.</description></myClass>
          !!]
          type, extends(myClassClass) :: myClassBound
          end type myClassBound
        contains
          function boundConstructorParameters(parameters) result(self)
            type(myClassBound   ) :: self
            type(inputParameters) :: parameters
            double precision      :: frac
            !![
            <inputParameter><name>frac</name><source>parameters</source>
              <description>d.</description>
              <minimum>0.0</minimum>
              <maximum inclusive="false">1.0</maximum></inputParameter>
            !!]
            return
          end function boundConstructorParameters
        end module M
        """))
    entries = harvest_file(str(impl), {"myClass"}, str(source_root))
    constraints = entries[0]['parameters'][0]['constraints']
    assert constraints['minimum'] == {'value': '0.0', 'inclusive': True}
    assert constraints['maximum'] == {'value': '1.0', 'inclusive': False}


def test_harvest_file_end_to_end(tmp_path):
    source_root = tmp_path / "source"
    (source_root / "my").mkdir(parents=True)
    impl = source_root / "my" / "impl.F90"
    impl.write_text(_SNIPPET)

    entries = harvest_file(str(impl), {"myClass"}, str(source_root))
    assert len(entries) == 1
    entry = entries[0]

    assert entry['type']            == 'myClassSimple'
    assert entry['functionClass']   == 'myClass'
    assert entry['label']           == 'simple'
    assert entry['constructorFound'] is True
    assert entry['sourceFile']      == 'my/impl.F90'

    params = {p['name']: p for p in entry['parameters']}
    assert set(params) == {'alpha', 'betaFlag', 'gamma', 'delta'}

    # alpha: real/double from the default literal.
    assert (params['alpha']['type'], params['alpha']['kind'],
            params['alpha']['provenance']) == ('real', 'double', 'default')
    # betaFlag: boolean from the default literal.
    assert (params['betaFlag']['type'], params['betaFlag']['provenance']) \
        == ('boolean', 'default')
    # gamma: no default -> typed real/double from its `double precision` declaration.
    assert (params['gamma']['type'], params['gamma']['kind'],
            params['gamma']['provenance']) == ('real', 'double', 'declaration')
    # delta: read from the nested `subParameters` handle -> sourceElement resolved.
    assert params['delta']['sourceElement'] == 'nestedElem'
    assert params['delta']['type'] == 'integer'       # from its `integer` declaration

    # objectBuilder edge captured, in the class's own element (no nesting).
    objects = {o['parameterName']: o for o in entry['objects']}
    assert objects['cosmologyFunctions']['class'] == 'cosmologyFunctions'
    assert objects['cosmologyFunctions']['sourceElement'] is None
    assert objects['cosmologyFunctions']['repeatable'] is False
