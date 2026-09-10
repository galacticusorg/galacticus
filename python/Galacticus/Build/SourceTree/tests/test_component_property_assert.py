"""Tests for `Galacticus.Build.SourceTree.Process.ComponentPropertyAssert`.

The directive is pure code generation, so the tests assert on the emitted
Fortran text and on the module uses wired into the enclosing unit:

  * the accessor test is `.and.`-joined across properties AND requirements;
  * the `AttributeMatch` list is `.intersection.`-joined across properties,
    with all requirements passed as arguments to each term;
  * `Array_Utilities` (which provides `.intersection.`) is imported only when
    there is more than one term to join — importing it for a single-property
    assertion would produce an unused-import warning;
  * `assignTo` emits the test alone (no `Error_Report`, no `Component_List`);
  * `condition` emits the report alone, guarded by the supplied expression;
  * malformed directives raise rather than emitting broken Fortran.
"""

import pytest

from Galacticus.Build.SourceTree.Process.ComponentPropertyAssert import (
    process_component_property_assert,
)


def _build(directive):
    """Minimal AST: a function unit holding a single directive node."""
    parent = {
        'type':       'function',
        'name':       'testConstructor',
        'firstChild': None,
        'sibling':    None,
        'parent':     None,
    }
    node = {
        'type':       'componentPropertyAssert',
        'directive':  dict(directive),
        'firstChild': None,
        'sibling':    None,
        'parent':     parent,
        'source':     '<test>',
        'line':       17,
    }
    parent['firstChild'] = node
    return parent, node


def _emitted(parent):
    """Concatenate the auto-generated code nodes hanging off `parent`."""
    text  = ""
    child = parent['firstChild']
    while child is not None:
        if child.get('type') == 'code':
            text += child.get('content') or ""
        child = child.get('sibling')
    return text


def _uses(parent):
    """Return `{module: set(only-symbols)}` for the unit's moduleUse node."""
    child = parent['firstChild']
    while child is not None:
        if child.get('type') == 'moduleUse':
            out = {}
            for module, entries in child['moduleUse'].items():
                if isinstance(entries, dict):
                    entries = [entries]
                symbols = set()
                for entry in entries:
                    symbols |= set((entry.get('only') or {}).keys())
                out[module] = symbols
            return out
        child = child.get('sibling')
    return {}


def _run(directive):
    parent, node = _build(directive)
    process_component_property_assert(parent, {})
    return parent, node


def _joined(text):
    """Collapse continuation lines so expressions can be matched as one string."""
    return text.replace(' &\n    & ', '').replace('&\n', '')


# ---------------------------------------------------------------------------
# Default (assert) mode
# ---------------------------------------------------------------------------

def test_single_property_emits_test_and_report():
    parent, node = _run({'class':      'darkMatterProfile',
                         'properties': 'scale',
                         'require':    'gettable'})
    text = _joined(_emitted(parent))

    assert node['directive']['processed'] is True
    assert ('if (.not.(defaultDarkMatterProfileComponent%scaleIsGettable()))'
            in text)
    assert "Component_List('darkMatterProfile'," in text
    assert ('defaultDarkMatterProfileComponent%scaleAttributeMatch'
            '(requireGettable=.true.)' in text)
    # No join operator is needed for a lone term.
    assert '.intersection.' not in text


def test_single_property_does_not_import_intersection():
    """`Array_Utilities` is only needed to join two or more match terms."""
    parent, _ = _run({'class':      'darkMatterProfile',
                      'properties': 'scale',
                      'require':    'gettable'})
    uses = _uses(parent)

    assert 'Array_Utilities' not in uses
    assert uses['Galacticus_Nodes'] == {'defaultDarkMatterProfileComponent'}
    assert uses['Error'] == {'Error_Report', 'Component_List'}


def test_multiple_properties_join_test_and_matches():
    parent, _ = _run({'class':      'spheroid',
                      'properties': 'massStellar massGas halfMassRadius',
                      'require':    'gettable'})
    text = _joined(_emitted(parent))
    uses = _uses(parent)

    assert ('if (.not.('
            'defaultSpheroidComponent%massStellarIsGettable()'
            '.and.defaultSpheroidComponent%massGasIsGettable()'
            '.and.defaultSpheroidComponent%halfMassRadiusIsGettable()'
            '))' in text)
    assert text.count('.intersection.') == 2
    assert uses['Array_Utilities'] == {'operator(.intersection.)'}


def test_multiple_requirements_apply_to_every_property():
    """Each requirement contributes its own `Is…()` term, and every
    requirement is passed to each `AttributeMatch` call."""
    parent, _ = _run({'class':      'spin',
                      'properties': 'angularMomentumVector',
                      'require':    'settable gettable'})
    text = _joined(_emitted(parent))

    assert ('defaultSpinComponent%angularMomentumVectorIsGettable()'
            '.and.defaultSpinComponent%angularMomentumVectorIsSettable()'
            in text)
    # Canonical ordering — independent of the order the directive listed them.
    assert ('angularMomentumVectorAttributeMatch'
            '(requireGettable=.true.,requireSettable=.true.)' in text)


def test_properties_accept_comma_separation():
    parent, _ = _run({'class':      'hotHalo',
                      'properties': 'mass, outerRadius',
                      'require':    'gettable'})
    text = _joined(_emitted(parent))

    assert 'defaultHotHaloComponent%massIsGettable()' in text
    assert 'defaultHotHaloComponent%outerRadiusIsGettable()' in text


def test_message_overrides_generated_text():
    parent, _ = _run({'class':      'hotHalo',
                      'properties': 'mass',
                      'require':    'gettable',
                      'message':    'cooling requires a hot halo mass.'})
    text = _emitted(parent)

    assert "'cooling requires a hot halo mass.'" in text
    assert 'must provide' not in text


def test_message_apostrophes_are_escaped():
    """An apostrophe in the message must be doubled or the literal breaks."""
    parent, _ = _run({'class':      'hotHalo',
                      'properties': 'mass',
                      'require':    'gettable',
                      'message':    "the halo's mass is required."})
    text = _emitted(parent)

    assert "'the halo''s mass is required.'" in text


def test_generated_message_names_class_properties_and_requirement():
    parent, _ = _run({'class':      'hotHalo',
                      'properties': 'mass outerRadius',
                      'require':    'gettable'})
    text = _emitted(parent)

    assert 'the "hotHalo" component must provide gettable' in text
    assert '"mass" and "outerRadius" properties' in text


# ---------------------------------------------------------------------------
# assignTo / condition modes
# ---------------------------------------------------------------------------

def test_assign_to_emits_test_only():
    parent, _ = _run({'class':      'disk',
                      'properties': 'radius velocity',
                      'require':    'gettable',
                      'assignTo':   'self%diskSupported'})
    text = _joined(_emitted(parent))
    uses = _uses(parent)

    assert ('self%diskSupported=defaultDiskComponent%radiusIsGettable()'
            '.and.defaultDiskComponent%velocityIsGettable()' in text)
    assert 'Error_Report' not in text
    assert 'Component_List' not in text
    # Only the component object is needed — no Error, no Array_Utilities.
    assert set(uses.keys()) == {'Galacticus_Nodes'}


def test_condition_emits_report_guarded_by_expression():
    parent, _ = _run({'class':      'disk',
                      'properties': 'radius velocity',
                      'require':    'gettable',
                      'condition':  'self%diskSupported'})
    text = _joined(_emitted(parent))

    assert 'if (.not.(self%diskSupported))' in text
    # The accessor test is NOT regenerated — the caller supplied the guard.
    assert 'defaultDiskComponent%radiusIsGettable()' not in text
    # …but the match list still is, since the report needs it.
    assert 'defaultDiskComponent%radiusAttributeMatch' in text


def test_assign_to_and_condition_are_mutually_exclusive():
    with pytest.raises(RuntimeError, match='mutually exclusive'):
        _run({'class':      'disk',
              'properties': 'radius',
              'require':    'gettable',
              'assignTo':   'self%diskSupported',
              'condition':  'self%diskSupported'})


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def test_unknown_requirement_raises():
    with pytest.raises(RuntimeError, match='unknown requirement'):
        _run({'class': 'disk', 'properties': 'radius', 'require': 'writable'})


def test_empty_properties_raises():
    with pytest.raises(RuntimeError, match='`properties` is empty'):
        _run({'class': 'disk', 'properties': '  ', 'require': 'gettable'})


def test_malformed_class_raises():
    with pytest.raises(RuntimeError, match='must be a component class name'):
        _run({'class': 'disk%radius', 'properties': 'radius',
              'require': 'gettable'})


def test_malformed_property_raises():
    with pytest.raises(RuntimeError, match='is not a valid property name'):
        _run({'class': 'disk', 'properties': 'radius()', 'require': 'gettable'})


def test_directive_is_processed_only_once():
    """A second pass must not emit the assertion twice."""
    parent, node = _run({'class':      'disk',
                         'properties': 'radius',
                         'require':    'gettable'})
    process_component_property_assert(parent, {})

    assert _emitted(parent).count('Auto-generated assertion') == 1
