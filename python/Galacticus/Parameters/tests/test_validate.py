"""Tests for `Galacticus.Parameters.validate`.

Drive the validator with a small synthetic catalog (so the cases are explicit)
covering: valid trees, bad selectors (wrong case and non-existent), unknown
parameter names, type mismatches, implementation inheritance, `sourceElement`
wrappers, direct-read `extras`, hoisted/shared objects, and root-scope globals
(which must never be flagged).
"""

import xml.etree.ElementTree as ET

import pytest

import Galacticus.Parameters.validate as validate_module
from Galacticus.Parameters.validate import (
    validate_parameters, validate_file, _check_references, _load_migrations,
)


def _param(name, type_, source_element=None):
    return {'name': name, 'type': type_, 'sourceElement': source_element,
            'kind': None, 'provenance': 'default', 'default': None,
            'variable': None, 'source': 'parameters', 'cardinality': None,
            'description': None}


def _object(class_, parameter_name, source_element=None, repeatable=False):
    return {'class': class_, 'name': parameter_name + '_',
            'parameterName': parameter_name, 'source': 'parameters',
            'sourceElement': source_element, 'repeatable': repeatable}


CATALOG = {
    'functionClasses': {
        'accretionHalo':       {'default': 'simple',
                                'implementations': ['coldMode', 'simple']},
        'cosmologyParameters': {'default': 'simple',
                                'implementations': ['simple']},
        'virialDensityContrast': {'default': 'fixed',
                                  'implementations': ['fixed']},
        'taskX':               {'default': 'builder',
                                'implementations': ['builder']},
    },
    'implementations': {
        'accretionHaloSimple': {
            'functionClass': 'accretionHalo', 'label': 'simple',
            'parent': 'accretionHaloClass',
            'parameters': [_param('alpha', 'real'), _param('flag', 'boolean')],
            'objects': [_object('cosmologyParameters', 'cosmologyParameters')],
            'directNames': [],
        },
        'accretionHaloColdMode': {            # extends accretionHaloSimple
            'functionClass': 'accretionHalo', 'label': 'coldMode',
            'parent': 'accretionHaloSimple',
            'parameters': [_param('beta', 'integer')],
            'objects': [], 'directNames': [],
        },
        'cosmologyParametersSimple': {
            'functionClass': 'cosmologyParameters', 'label': 'simple',
            'parent': 'cosmologyParametersClass',
            'parameters': [_param('OmegaMatter', 'real')],
            'objects': [], 'directNames': [],
        },
        'virialDensityContrastFixed': {
            'functionClass': 'virialDensityContrast', 'label': 'fixed',
            'parent': 'virialDensityContrastClass',
            'parameters': [dict(_param('densityType', 'string'),
                                enumeration='densityKind')],
            'objects': [], 'directNames': [],
        },
        'taskXBuilder': {
            'functionClass': 'taskX', 'label': 'builder',
            'parent': 'taskXClass',
            'parameters': [_param('labels', 'string', source_element='massDefinitions')],
            'objects': [_object('virialDensityContrast', 'virialDensityContrast',
                                source_element='massDefinitions', repeatable=True)],
            'directNames': ['particleProperty'],
        },
    },
    'enumerations': {
        'densityKind': ['critical', 'mean'],
    },
}


def _validate(xml_text):
    return validate_parameters(ET.fromstring(xml_text), CATALOG)


def test_valid_tree_has_no_findings():
    findings = _validate("""
        <parameters>
          <accretionHalo value="simple">
            <alpha value="1.0"/>
            <flag value="true"/>
            <cosmologyParameters value="simple"><OmegaMatter value="0.3"/></cosmologyParameters>
          </accretionHalo>
        </parameters>
    """)
    assert findings == []


def test_selector_wrong_case():
    findings = _validate('<parameters><accretionHalo value="Simple"/></parameters>')
    assert len(findings) == 1
    assert findings[0].kind == 'selector'
    assert "wrong case" in findings[0].message
    assert "'simple'" in findings[0].message


def test_selector_nonexistent():
    findings = _validate('<parameters><accretionHalo value="bogus"/></parameters>')
    assert [f.kind for f in findings] == ['selector']
    assert "not an implementation" in findings[0].message


def test_unknown_parameter_inside_scope():
    findings = _validate(
        '<parameters><accretionHalo value="simple"><nonsense value="1"/>'
        '</accretionHalo></parameters>')
    assert [f.kind for f in findings] == ['unknown']
    assert "nonsense" in findings[0].message


def test_type_mismatch_is_a_warning():
    findings = _validate(
        '<parameters><accretionHalo value="simple"><flag value="3.5"/>'
        '</accretionHalo></parameters>')
    assert [(f.level, f.kind) for f in findings] == [('warning', 'type')]


def test_inherited_parameters_accepted():
    """coldMode extends simple, so alpha/flag (simple's) are accepted alongside
    beta (its own)."""
    findings = _validate("""
        <parameters>
          <accretionHalo value="coldMode">
            <alpha value="1.0"/><flag value="false"/><beta value="2"/>
          </accretionHalo>
        </parameters>
    """)
    assert findings == []


def test_source_element_wrapper_descended():
    """`labels` and repeated `virialDensityContrast` live inside the
    <massDefinitions> wrapper; both must validate without error."""
    findings = _validate("""
        <parameters>
          <taskX value="builder">
            <massDefinitions>
              <labels value="a b"/>
              <virialDensityContrast value="fixed"/>
              <virialDensityContrast value="fixed"/>
            </massDefinitions>
          </taskX>
        </parameters>
    """)
    assert findings == []


def test_direct_read_extras_accepted_and_not_recursed():
    """A direct-read wrapper (particleProperty) is accepted and its children are
    not validated (opaque to the catalog)."""
    findings = _validate("""
        <parameters>
          <taskX value="builder">
            <particleProperty><whatever value="x"/></particleProperty>
          </taskX>
        </parameters>
    """)
    assert findings == []


def test_root_globals_and_meta_not_flagged():
    """Scalars at the root scope (globals/meta) are never name-checked."""
    findings = _validate("""
        <parameters>
          <formatVersion>2</formatVersion>
          <someGlobalParameter value="x"/>
          <cosmologyParameters value="simple"><OmegaMatter value="0.3"/></cosmologyParameters>
        </parameters>
    """)
    assert findings == []


def test_hoisted_shared_object_not_flagged_as_unknown():
    """A functionClass element appearing where the parent impl does not declare
    it (hoisted/shared) is validated and recursed, not flagged."""
    findings = _validate("""
        <parameters>
          <accretionHalo value="simple">
            <alpha value="1.0"/>
            <virialDensityContrast value="fixed"/>
          </accretionHalo>
        </parameters>
    """)
    assert findings == []        # virialDensityContrast is a valid selector, recursed


def test_ignore_warnings_suppresses_unknown():
    """`ignoreWarnings="true"` opts a parameter out of the unknown-name check,
    matching Galacticus's own convention."""
    findings = _validate(
        '<parameters><accretionHalo value="simple">'
        '<nonsense ignoreWarnings="true" value="1"/>'
        '</accretionHalo></parameters>')
    assert findings == []


_XI = 'xmlns:xi="http://www.w3.org/2001/XInclude"'


def test_xinclude_expands_and_passes(tmp_path):
    (tmp_path / "ref.xml").write_text(
        '<parameters><cosmologyParameters value="simple">'
        '<OmegaMatter value="0.3"/></cosmologyParameters></parameters>')
    main = tmp_path / "main.xml"
    main.write_text(
        f'<parameters><xi:include {_XI} href="ref.xml" '
        'xpointer="xpointer(parameters/*)"/></parameters>')
    findings, error = validate_file(str(main), CATALOG)
    assert error is None
    assert findings == []


def test_xinclude_validates_included_content(tmp_path):
    """A bad selector in the *referenced* file must be caught -- proving the
    included content is actually validated."""
    (tmp_path / "ref.xml").write_text(
        '<parameters><cosmologyParameters value="WrongCase"/></parameters>')
    main = tmp_path / "main.xml"
    main.write_text(
        f'<parameters><xi:include {_XI} href="ref.xml" '
        'xpointer="xpointer(parameters/*)"/></parameters>')
    findings, error = validate_file(str(main), CATALOG)
    assert [f.kind for f in findings] == ['selector']


def test_xinclude_unresolved_href_reported(tmp_path):
    main = tmp_path / "main.xml"
    main.write_text(
        f'<parameters><xi:include {_XI} href="missing.xml" '
        'xpointer="xpointer(parameters/*)"/></parameters>')
    findings, error = validate_file(str(main), CATALOG)
    assert error is None
    assert [f.kind for f in findings] == ['include']
    assert "does not resolve" in findings[0].message


def test_xinclude_nested(tmp_path):
    """An include whose target itself includes a third file expands fully."""
    (tmp_path / "inner.xml").write_text(
        '<parameters><cosmologyParameters value="simple">'
        '<OmegaMatter value="0.3"/></cosmologyParameters></parameters>')
    (tmp_path / "mid.xml").write_text(
        f'<parameters><xi:include {_XI} href="inner.xml" '
        'xpointer="xpointer(parameters/*)"/></parameters>')
    main = tmp_path / "main.xml"
    main.write_text(
        f'<parameters><xi:include {_XI} href="mid.xml" '
        'xpointer="xpointer(parameters/*)"/></parameters>')
    findings, error = validate_file(str(main), CATALOG)
    assert error is None
    assert findings == []


def test_xinclude_template_fallback_to_repo_root(tmp_path, monkeypatch):
    """A template file run from elsewhere has `../`-depth relative to the run
    location, not its source; resolution falls back to the repo root so the
    referenced content is still found and validated."""
    reference_dir = tmp_path / "parameters" / "reference"
    reference_dir.mkdir(parents=True)
    (reference_dir / "ref.xml").write_text(
        '<parameters><cosmologyParameters value="simple">'
        '<OmegaMatter value="0.3"/></cosmologyParameters></parameters>')
    monkeypatch.setenv("GALACTICUS_EXEC_PATH", str(tmp_path))
    source_dir = tmp_path / "testSuite" / "parameters"
    source_dir.mkdir(parents=True)
    main = source_dir / "main.xml"
    main.write_text(
        f'<parameters><xi:include {_XI} '
        'href="../../../../parameters/reference/ref.xml" '
        'xpointer="xpointer(parameters/*)"/></parameters>')
    findings, error = validate_file(str(main), CATALOG)
    assert error is None
    assert findings == []


@pytest.mark.parametrize("value, flagged", [
    ('critical',                 False),   # valid bare label
    ('mean',                     False),   # valid bare label
    ('densityKindCritical',      False),   # valid prefixed form
    ('bogus',                    True),    # not a label
    ('Critical',                 True),    # wrong case (matching is case-sensitive)
])
def test_enumeration_value_checked(value, flagged):
    findings = _validate(
        f'<parameters><virialDensityContrast value="fixed">'
        f'<densityType value="{value}"/></virialDensityContrast></parameters>')
    if flagged:
        assert [f.kind for f in findings] == ['enumeration']
    else:
        assert findings == []


def test_enumeration_expression_value_not_judged():
    """A non-literal value (reference/expression) is not enum-checked."""
    findings = _validate(
        '<parameters><virialDensityContrast value="fixed">'
        '<densityType value="[someOtherParameter]"/></virialDensityContrast></parameters>')
    assert findings == []


def test_enumeration_ignore_warnings_suppresses():
    findings = _validate(
        '<parameters><virialDensityContrast value="fixed">'
        '<densityType ignoreWarnings="true" value="bogus"/>'
        '</virialDensityContrast></parameters>')
    assert findings == []


_CONSTRAINT_CATALOG = {
    'functionClasses': {'widget': {'default': 'basic', 'implementations': ['basic']}},
    'implementations': {
        'widgetBasic': {
            'functionClass': 'widget', 'label': 'basic', 'parent': 'widgetClass',
            'parameters': [
                dict(_param('fraction', 'real'),
                     constraints={'minimum': {'value': '0.0', 'inclusive': True},
                                  'maximum': {'value': '1.0', 'inclusive': True}}),
                dict(_param('positive', 'real'),
                     constraints={'minimum': {'value': '0.0', 'inclusive': False}}),
                dict(_param('mode', 'string'),
                     constraints={'allowedValues': ['fast', 'slow']}),
            ],
            'objects': [], 'directNames': [],
        },
    },
    'enumerations': {},
}


def _validate_constraint(inner):
    xml = f'<parameters><widget value="basic">{inner}</widget></parameters>'
    return validate_parameters(ET.fromstring(xml), _CONSTRAINT_CATALOG)


@pytest.mark.parametrize("inner, flagged", [
    ('<fraction value="0.5"/>',  False),   # inside [0,1]
    ('<fraction value="0.0"/>',  False),   # inclusive minimum
    ('<fraction value="1.0"/>',  False),   # inclusive maximum
    ('<fraction value="1.5"/>',  True),    # above maximum
    ('<fraction value="-0.1"/>', True),    # below minimum
    ('<positive value="0.1"/>',  False),   # above exclusive minimum
    ('<positive value="0.0"/>',  True),    # equals exclusive minimum
    ('<mode value="fast"/>',     False),   # allowed value
    ('<mode value="medium"/>',   True),    # not an allowed value
])
def test_constraint_enforcement(inner, flagged):
    findings = _validate_constraint(inner)
    if flagged:
        assert [f.kind for f in findings] == ['constraint']
    else:
        assert findings == []


def test_constraint_skips_expression_value():
    """A non-literal value (reference) is not range-checked."""
    assert _validate_constraint('<fraction value="[otherParam]"/>') == []


def test_constraint_honours_ignore_warnings():
    assert _validate_constraint(
        '<fraction ignoreWarnings="true" value="9.9"/>') == []


def test_constraint_fortran_exponent_value():
    """`d`-exponent literals parse for range checking."""
    assert [f.kind for f in _validate_constraint('<fraction value="2.0d0"/>')] \
        == ['constraint']


def _references(xml_text):
    findings = []
    _check_references(ET.fromstring(xml_text), findings)
    return findings


def test_idref_resolves():
    assert _references(
        '<parameters>'
        '<cosmologyParameters id="c" value="simple"/>'
        '<cosmologyParameters idRef="c"/>'
        '</parameters>') == []


def test_idref_dangling_flagged():
    findings = _references(
        '<parameters><cosmologyParameters idRef="missing"/></parameters>')
    assert [f.kind for f in findings] == ['reference']
    assert "no matching element" in findings[0].message


def test_idref_requires_same_tag():
    """An `id` on a different element tag does not satisfy an `idRef`."""
    findings = _references(
        '<parameters><foo id="c" value="x"/>'
        '<cosmologyParameters idRef="c"/></parameters>')
    assert [f.kind for f in findings] == ['reference']


def test_load_migrations(tmp_path):
    aux = tmp_path / "scripts" / "aux"
    aux.mkdir(parents=True)
    (aux / "migrations.xml").write_text(
        '<migrations><migration commit="x">'
        '<translation xpath="//nodeOperator[@value=\'y\']/massDestruction[@value]">'
        '<name new="massDestructionAbsolute"/></translation>'
        '<translation xpath="//nodePropertyExtractor[@value]">'
        '<value old="descendents" new="descendants"/></translation>'
        '</migration></migrations>')
    name_renames, value_renames = _load_migrations(str(tmp_path))
    assert name_renames["massDestruction"] == "massDestructionAbsolute"
    assert value_renames["descendents"] == "descendants"


def test_unknown_message_includes_migration_hint(tmp_path, monkeypatch):
    aux = tmp_path / "scripts" / "aux"
    aux.mkdir(parents=True)
    (aux / "migrations.xml").write_text(
        '<migrations><migration commit="x">'
        '<translation xpath="//x/oldName[@value]"><name new="newName"/></translation>'
        '</migration></migrations>')
    monkeypatch.setenv("GALACTICUS_EXEC_PATH", str(tmp_path))
    validate_module._migrations_cache.clear()
    try:
        findings = _validate(
            '<parameters><accretionHalo value="simple">'
            '<oldName value="1"/></accretionHalo></parameters>')
    finally:
        validate_module._migrations_cache.clear()
    assert any('obsolete' in f.message and 'newName' in f.message for f in findings)


def test_bad_selector_on_nested_object():
    findings = _validate("""
        <parameters>
          <accretionHalo value="simple">
            <cosmologyParameters value="WrongCase"/>
          </accretionHalo>
        </parameters>
    """)
    assert [f.kind for f in findings] == ['selector']
    assert "cosmologyParameters" in findings[0].message
