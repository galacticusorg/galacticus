"""Tests for scripts/build/parameterSchema.py (parameter-file XSD generator)."""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
from parameterSchema import build_schema   # noqa: E402

CATALOG = {
    'functionClasses': {
        'accretionHalo':        {'implementations': ['coldMode', 'simple']},
        'criticalOverdensity':  {'implementations': ['fixed']},   # also a parameter -> collision
    },
    'enumerations': {'densityKind': ['critical', 'mean']},
    'implementations': {
        'accretionHaloSimple': {'parameters': [
            {'name': 'alpha', 'enumeration': None},
            {'name': 'densityType', 'enumeration': 'densityKind'},
        ]},
        # makes `criticalOverdensity` also a scalar parameter name
        'haloMassFunctionFixed': {'parameters': [
            {'name': 'criticalOverdensity', 'enumeration': None},
        ]},
    },
}


def test_selector_enum_emitted():
    xsd = build_schema(CATALOG)
    assert '<xs:element name="accretionHalo">' in xsd
    assert '<xs:enumeration value="simple"/>' in xsd
    assert '<xs:enumeration value="coldMode"/>' in xsd


def test_enum_parameter_emitted():
    xsd = build_schema(CATALOG)
    assert '<xs:element name="densityType">' in xsd
    assert '<xs:enumeration value="critical"/>' in xsd
    assert '<xs:enumeration value="mean"/>' in xsd


def test_collision_name_is_unconstrained():
    """A name that is both a functionClass base and a scalar parameter must not
    have its value enumerated (its value may be a label or an ordinary value)."""
    xsd = build_schema(CATALOG)
    assert '<xs:element name="criticalOverdensity">' in xsd
    # 'fixed' is criticalOverdensity's only implementation; it must NOT appear as
    # an enumeration (the element is generic).
    assert '<xs:enumeration value="fixed"/>' not in xsd


def test_tolerated_roots_declared():
    """Roots of non-parameter files that share parameter directories (changes,
    merger trees, ...) are declared so a directory-scoped editor association does
    not flag them."""
    xsd = build_schema(CATALOG)
    for root in ('changes', 'mergerTrees', 'parameterGrid', 'tree'):
        assert f'<xs:element name="{root}">' in xsd


def test_schema_compiles_and_validates():
    etree = pytest.importorskip('lxml.etree')
    schema = etree.XMLSchema(etree.fromstring(build_schema(CATALOG).encode()))
    assert schema.validate(etree.fromstring(
        b'<parameters><accretionHalo value="simple">'
        b'<densityType value="critical"/></accretionHalo></parameters>'))
    # bad selector / bad enum value are rejected
    assert not schema.validate(etree.fromstring(
        b'<parameters><accretionHalo value="bogus"/></parameters>'))
    assert not schema.validate(etree.fromstring(
        b'<parameters><densityType value="bogus"/></parameters>'))
    # collision name accepts an ordinary (numeric) value
    assert schema.validate(etree.fromstring(
        b'<parameters><criticalOverdensity value="1.686"/></parameters>'))
    # a non-parameter file root (e.g. a `<changes>` file co-located with
    # parameter files) is tolerated rather than rejected
    assert schema.validate(etree.fromstring(
        b'<changes><change type="update"/></changes>'))


# ---------------------------------------------------------------------------
# `main`: --check reports drift without writing (used by CI and the git hook)
# ---------------------------------------------------------------------------

@pytest.fixture
def fake_tree(tmp_path, monkeypatch):
    """A source directory whose catalog is stubbed out, so the CLI's write /
    check / exit-code behaviour is exercised without parsing real Fortran."""
    import parameterSchema

    (tmp_path / 'source').mkdir()
    (tmp_path / 'schema').mkdir()
    monkeypatch.setattr(parameterSchema, 'build_catalog',
                        lambda source_root, **keywords: CATALOG)
    # Keep any cache the CLI resolves inside the fixture's own directory.
    monkeypatch.setenv('BUILDPATH', str(tmp_path / 'build'))
    return tmp_path


def _run(tmp_path, *arguments):
    import parameterSchema
    return parameterSchema.main(
        ['parameterSchema.py', str(tmp_path), str(tmp_path / 'out.xsd'),
         *arguments])


def test_writes_schema(fake_tree):
    assert _run(fake_tree) == 0
    assert (fake_tree / 'out.xsd').read_text() == build_schema(CATALOG)


def test_check_passes_when_up_to_date(fake_tree):
    assert _run(fake_tree) == 0
    assert _run(fake_tree, '--check') == 0


def test_check_fails_when_out_of_date(fake_tree, capsys):
    (fake_tree / 'out.xsd').write_text('<xs:schema/>\n')
    assert _run(fake_tree, '--check') == 1
    assert 'out of date' in capsys.readouterr().err


def test_check_fails_when_absent(fake_tree, capsys):
    assert _run(fake_tree, '--check') == 1
    assert 'does not exist' in capsys.readouterr().err


def test_check_does_not_write(fake_tree):
    """The whole point of --check: CI and the pre-commit hook must not have the
    working tree mutated underneath them."""
    stale = '<xs:schema/>\n'
    (fake_tree / 'out.xsd').write_text(stale)
    assert _run(fake_tree, '--check') == 1
    assert (fake_tree / 'out.xsd').read_text() == stale


def test_check_does_not_create_the_file(fake_tree):
    assert _run(fake_tree, '--check') == 1
    assert not (fake_tree / 'out.xsd').exists()


def test_default_output_path(fake_tree):
    import parameterSchema
    assert parameterSchema.main(['parameterSchema.py', str(fake_tree)]) == 0
    assert (fake_tree / 'schema' / 'parameters.xsd').read_text() == \
        build_schema(CATALOG)


def test_cli_wiring(tmp_path, monkeypatch):
    """--jobs and --no-cache must reach `build_catalog`, and the cache must be
    on by default (the pre-commit hook depends on it)."""
    import parameterSchema

    (tmp_path / 'source').mkdir()
    monkeypatch.setenv('BUILDPATH', str(tmp_path / 'build'))
    calls = []

    def _stub_catalog(source_root, **keywords):
        calls.append(keywords)
        return CATALOG

    monkeypatch.setattr(parameterSchema, 'build_catalog', _stub_catalog)

    def run(*arguments):
        calls.clear()
        parameterSchema.main(['parameterSchema.py', str(tmp_path),
                              str(tmp_path / 'out.xsd'), *arguments])
        return calls[0]

    assert run('--jobs', '3')['jobs'] == 3
    assert run()['jobs'] is None                     # left to ParallelScan
    assert run('--no-cache')['cache_path'] is None
    assert run()['cache_path'] == str(
        tmp_path / 'build' / 'parameters.catalog.blob')
