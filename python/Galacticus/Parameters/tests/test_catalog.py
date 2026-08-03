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

import os
import pickle
import textwrap
import time

import pytest

from Galacticus.Parameters.catalog import (
    derive_label, _resolve_source_elements, harvest_file,
    discover_enumerations, _enumeration_links, _capture_constraints,
    _normalize_default, _scan_directives, build_catalog,
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


# ---------------------------------------------------------------------------
# `_scan_directives` / `build_catalog`: the single-read early phases
# ---------------------------------------------------------------------------

def _write_tree(tmp_path):
    """A miniature source tree: one base class, one enumeration, one
    implementation, and one file that registers nothing."""
    source_root = tmp_path / "source"
    (source_root / "my").mkdir(parents=True)
    (source_root / "my" / "base.F90").write_text(textwrap.dedent("""\
        module Base
          !![
          <functionClass>
           <name>myClass</name>
           <descriptiveName>My Class</descriptiveName>
           <default>enum</default>
          </functionClass>
          <enumeration>
           <name>fixedDensityType</name>
           <entry label="critical"/>
           <entry label="mean"/>
          </enumeration>
          !!]
        end module Base
        """))
    (source_root / "my" / "impl.F90").write_text(_SNIPPET_ENUM)
    (source_root / "my" / "plain.F90").write_text(
        "module Plain\nend module Plain\n")
    return source_root


def test_scan_directives_collects_all_phases_in_one_read(tmp_path):
    """One read must yield the bases, the enumerations, and enough information
    to decide whether the file needs a full parse."""
    source_root = _write_tree(tmp_path)

    bases, enumerations, named_roots = _scan_directives(
        str(source_root / "my" / "base.F90"))
    assert [name for name, _ in bases] == ['myClass']
    assert enumerations == [('fixedDensityType', ['critical', 'mean'])]
    # `base.F90` declares the class but registers no implementation of it.
    assert 'myClass' not in named_roots

    bases, enumerations, named_roots = _scan_directives(
        str(source_root / "my" / "impl.F90"))
    assert bases == [] and enumerations == []
    # `impl.F90` carries `<myClass name="myClassEnum">` -> needs a parse.
    assert 'myClass' in named_roots


def test_scan_directives_on_file_without_directives(tmp_path):
    source_root = _write_tree(tmp_path)
    assert _scan_directives(str(source_root / "my" / "plain.F90")) == (
        [], [], set())


def test_build_catalog_end_to_end(tmp_path):
    source_root = _write_tree(tmp_path)
    catalog = build_catalog(str(source_root), jobs=1)

    assert catalog['enumerations'] == {'fixedDensityType': ['critical', 'mean']}
    assert catalog['functionClasses']['myClass']['implementations'] == ['enum']
    assert catalog['functionClasses']['myClass']['descriptiveName'] == 'My Class'
    implementation = catalog['implementations']['myClassEnum']
    assert implementation['functionClass'] == 'myClass'
    assert implementation['label'] == 'enum'


def test_build_catalog_parallel_matches_serial(tmp_path):
    """The catalog feeds a committed, diffed artefact: it must not depend on the
    worker count."""
    source_root = _write_tree(tmp_path)
    assert build_catalog(str(source_root), jobs=4) == \
           build_catalog(str(source_root), jobs=1)


def test_harvest_worker_contains_parse_failure(tmp_path, monkeypatch):
    """A file the parser chokes on is reported, not raised.

    `ParallelScan` aborts the run if a worker raises, so a single bad file would
    otherwise take the whole catalog down -- where the serial code it replaced
    logged it and carried on. The failure is injected because the real parser is
    forgiving enough that a genuinely unparseable fixture would be contrived.
    """
    import Galacticus.Parameters.catalog as catalog_module

    def _explode(path, base_names, source_root):
        raise RuntimeError("parser exploded")

    monkeypatch.setattr(catalog_module, 'harvest_file', _explode)
    monkeypatch.setitem(catalog_module._WORKER, 'base_names', {'myClass'})
    monkeypatch.setitem(catalog_module._WORKER, 'source_root', str(tmp_path))

    entries, error = catalog_module._harvest_worker(str(tmp_path / "any.F90"))
    assert entries == []
    assert "parser exploded" in error


def test_build_catalog_logs_parse_error_and_keeps_the_rest(tmp_path, monkeypatch):
    """One bad file must not lose the rest of the catalog."""
    import Galacticus.Parameters.catalog as catalog_module
    source_root = _write_tree(tmp_path)
    (source_root / "my" / "broken.F90").write_text(
        _SNIPPET_ENUM.replace('myClassEnum', 'myClassBroken')
                     .replace('enumConstructorParameters',
                              'brokenConstructorParameters'))

    real_harvest = catalog_module.harvest_file

    def _explode_on_broken(path, base_names, source_root_):
        if path.endswith('broken.F90'):
            raise RuntimeError("parser exploded")
        return real_harvest(path, base_names, source_root_)

    monkeypatch.setattr(catalog_module, 'harvest_file', _explode_on_broken)

    messages = []
    # jobs=1: the monkeypatch lives in this process, so it must not be forked away.
    catalog = build_catalog(str(source_root), log=messages.append, jobs=1)

    assert 'myClassEnum' in catalog['implementations']
    assert 'myClassBroken' not in catalog['implementations']
    assert catalog['functionClasses']['myClass']['implementations'] == ['enum']
    assert any('parse-error' in m and 'broken.F90' in m for m in messages)


# ---------------------------------------------------------------------------
# The scan cache
# ---------------------------------------------------------------------------

def _backdate(source_root, seconds=10):
    """Age the source files so they are unambiguously older than a cache written
    now (mtime comparison is `>=`, so equal timestamps force a rescan)."""
    past = time.time() - seconds
    for dirpath, _, filenames in os.walk(str(source_root)):
        for name in filenames:
            os.utime(os.path.join(dirpath, name), (past, past))


def test_cache_matches_uncached(tmp_path):
    """The cache is an optimisation: it must never change the catalog."""
    source_root = _write_tree(tmp_path)
    blob = str(tmp_path / 'cache.blob')

    uncached = build_catalog(str(source_root), jobs=1)
    build_catalog(str(source_root), jobs=1, cache_path=blob)   # populate
    _backdate(source_root)
    warm = build_catalog(str(source_root), jobs=1, cache_path=blob)

    assert warm == uncached


def test_repeated_cached_runs_are_stable(tmp_path):
    """Regression: the cached base-class record must not accumulate labels.

    `build_catalog` appends implementation labels to each base record. When that
    record was owned by the scan cache, every run appended to the *same* list and
    then pickled it, so labels piled up run over run (and survived deleting the
    file that declared them).
    """
    source_root = _write_tree(tmp_path)
    blob = str(tmp_path / 'cache.blob')

    first = build_catalog(str(source_root), jobs=1, cache_path=blob)
    _backdate(source_root)
    for _ in range(3):
        again = build_catalog(str(source_root), jobs=1, cache_path=blob)
        assert again['functionClasses']['myClass']['implementations'] == ['enum']
        assert again == first


def test_cache_detects_modified_file(tmp_path):
    source_root = _write_tree(tmp_path)
    blob = str(tmp_path / 'cache.blob')
    build_catalog(str(source_root), jobs=1, cache_path=blob)

    (source_root / "my" / "impl.F90").write_text(
        _SNIPPET_ENUM.replace('myClassEnum', 'myClassRenamed'))
    catalog = build_catalog(str(source_root), jobs=1, cache_path=blob)

    assert 'myClassRenamed' in catalog['implementations']
    assert 'myClassEnum' not in catalog['implementations']
    assert catalog['functionClasses']['myClass']['implementations'] == ['renamed']


def test_cache_detects_new_and_deleted_files(tmp_path):
    source_root = _write_tree(tmp_path)
    blob = str(tmp_path / 'cache.blob')
    build_catalog(str(source_root), jobs=1, cache_path=blob)
    _backdate(source_root)

    extra = source_root / "my" / "extra.F90"
    extra.write_text(_SNIPPET_ENUM.replace('myClassEnum', 'myClassExtra')
                                  .replace('enumConstructorParameters',
                                           'extraConstructorParameters'))
    catalog = build_catalog(str(source_root), jobs=1, cache_path=blob)
    assert catalog['functionClasses']['myClass']['implementations'] == \
        ['enum', 'extra']

    extra.unlink()
    catalog = build_catalog(str(source_root), jobs=1, cache_path=blob)
    assert catalog['functionClasses']['myClass']['implementations'] == ['enum']
    assert 'myClassExtra' not in catalog['implementations']


def test_cache_drops_entries_for_deleted_files(tmp_path):
    """A deleted file must not linger in the blob forever."""
    source_root = _write_tree(tmp_path)
    blob = str(tmp_path / 'cache.blob')
    extra = source_root / "my" / "extra.F90"
    extra.write_text("module Extra\nend module Extra\n")
    build_catalog(str(source_root), jobs=1, cache_path=blob)
    assert any('extra.F90' in key for key in pickle.load(open(blob, 'rb')))

    extra.unlink()
    build_catalog(str(source_root), jobs=1, cache_path=blob)
    assert not any('extra.F90' in key for key in pickle.load(open(blob, 'rb')))


def test_cache_reuses_harvest_when_an_unrelated_base_is_added(tmp_path,
                                                              monkeypatch):
    """A new functionClass elsewhere must not force the tree to be re-harvested.

    A cached harvest depends on `base_names` only through this file's own
    registrations, so it stays valid while those are unchanged.
    """
    import Galacticus.Parameters.catalog as catalog_module
    source_root = _write_tree(tmp_path)
    blob = str(tmp_path / 'cache.blob')
    build_catalog(str(source_root), jobs=1, cache_path=blob)
    _backdate(source_root)

    (source_root / "my" / "other.F90").write_text(textwrap.dedent("""\
        module Other
          !![
          <functionClass>
           <name>otherClass</name>
           <descriptiveName>Other</descriptiveName>
          </functionClass>
          !!]
        end module Other
        """))

    harvested = []
    real_harvest = catalog_module.harvest_file
    monkeypatch.setattr(catalog_module, 'harvest_file',
                        lambda p, b, s: (harvested.append(p),
                                         real_harvest(p, b, s))[1])
    catalog = build_catalog(str(source_root), jobs=1, cache_path=blob)

    assert 'otherClass' in catalog['functionClasses']
    # impl.F90's registrations are unchanged, so it is not re-parsed.
    assert not any('impl.F90' in p for p in harvested)


def test_corrupt_cache_falls_back_to_full_scan(tmp_path):
    source_root = _write_tree(tmp_path)
    blob = tmp_path / 'cache.blob'
    blob.write_bytes(b'not a pickle at all')
    catalog = build_catalog(str(source_root), jobs=1, cache_path=str(blob))
    assert catalog == build_catalog(str(source_root), jobs=1)


def test_foreign_format_cache_is_discarded(tmp_path):
    """A blob from an older entry layout must be rebuilt, not misread."""
    source_root = _write_tree(tmp_path)
    blob = tmp_path / 'cache.blob'
    with open(blob, 'wb') as fh:
        pickle.dump({'__format__': -1, 'source_my_impl.F90': 'garbage'}, fh)
    catalog = build_catalog(str(source_root), jobs=1, cache_path=str(blob))
    assert catalog == build_catalog(str(source_root), jobs=1)


def test_unwritable_cache_does_not_fail_the_build(tmp_path):
    source_root = _write_tree(tmp_path)
    messages = []
    catalog = build_catalog(str(source_root), jobs=1, log=messages.append,
                            cache_path='/proc/nonexistent/cache.blob')
    assert catalog == build_catalog(str(source_root), jobs=1)
    assert any('cache' in m for m in messages)


def test_cache_invalidated_when_catalog_code_changes(tmp_path, monkeypatch):
    """Editing the code that builds entries must invalidate the blob.

    Otherwise working on the catalog itself silently reuses results computed by
    the previous version -- and `--check`, whose whole job is to answer whether
    the schema is stale, would answer from them.
    """
    import Galacticus.Parameters.catalog as catalog_module
    source_root = _write_tree(tmp_path)
    blob = str(tmp_path / 'cache.blob')
    build_catalog(str(source_root), jobs=1, cache_path=blob)
    _backdate(source_root)

    # Same sources, same mtimes -- but different catalog code.
    monkeypatch.setattr(catalog_module, '_code_fingerprint',
                        lambda: 'a-different-version-of-the-code')
    rescanned = []
    real_scan = catalog_module._scan_directives
    monkeypatch.setattr(catalog_module, '_scan_directives',
                        lambda p: (rescanned.append(p), real_scan(p))[1])

    build_catalog(str(source_root), jobs=1, cache_path=blob)
    assert rescanned, "cache was reused despite the code having changed"


def test_code_fingerprint_is_stable_and_specific():
    from Galacticus.Parameters.catalog import _code_fingerprint
    assert _code_fingerprint() == _code_fingerprint()
    assert len(_code_fingerprint()) == 64      # sha256 hex
