"""Tests for the capability check of `Galacticus.Parameters.validate`: an object
which requires an optional argument that the implementation building it
withholds must be reported, however that object is selected.

Andrew Benson (2026).
"""

import xml.etree.ElementTree as ET

from Galacticus.Parameters.validate import validate_parameters


def _implementation(base, label, objects=(), requires=(), forwards=()):
    return {
        'functionClass': base, 'label': label, 'parent': base + 'Class',
        'parameters': [], 'objects': list(objects), 'directNames': [],
        'requires': [{'method': m, 'argument': a} for m, a in requires],
        'forwards': list(forwards),
    }


def _object(class_, parameter_name, withholds=()):
    return {'class': class_, 'name': parameter_name + '_',
            'parameterName': parameter_name, 'source': 'parameters',
            'sourceElement': None, 'repeatable': False,
            'withholds': [{'method': m, 'argument': a} for m, a in withholds]}


CATALOG = {
    'functionClasses': {
        # The default is the implementation with a requirement, so that a
        # consumer given no mass function at all must also be caught.
        'haloMassFunction': {'default': 'pressSchechter',
                             'implementations': ['multiplier', 'pressSchechter', 'tinker2008']},
        'task':             {'default': 'nodal',
                             'implementations': ['integrator', 'nodal']},
    },
    'implementations': {
        'haloMassFunctionPressSchechter': _implementation(
            'haloMassFunction', 'pressSchechter', requires=[('differential', 'node')]),
        'haloMassFunctionTinker2008': _implementation('haloMassFunction', 'tinker2008'),
        'haloMassFunctionMultiplier': _implementation(
            'haloMassFunction', 'multiplier',
            objects=[_object('haloMassFunction', 'haloMassFunction')],
            forwards=['haloMassFunction_']),
        # Never supplies a node...
        'taskIntegrator': _implementation(
            'task', 'integrator',
            objects=[_object('haloMassFunction', 'haloMassFunction', withholds=[('differential', 'node')])]),
        # ...while this one does.
        'taskNodal': _implementation(
            'task', 'nodal', objects=[_object('haloMassFunction', 'haloMassFunction')]),
    },
    'enumerations': {},
}


def _capability_findings(xml_text):
    findings = validate_parameters(ET.fromstring(xml_text), CATALOG)
    return [f for f in findings if f.kind == 'capability']


def test_compatible_implementation_not_flagged():
    assert _capability_findings("""
        <parameters>
          <task value="integrator"><haloMassFunction value="tinker2008"/></task>
        </parameters>""") == []


def test_explicit_requirement_flagged():
    findings = _capability_findings("""
        <parameters>
          <task value="integrator"><haloMassFunction value="pressSchechter"/></task>
        </parameters>""")
    assert len(findings) == 1
    assert findings[0].level == 'error'
    assert findings[0].path == 'parameters/task/haloMassFunction'
    assert "'node'" in findings[0].message and "'differential'" in findings[0].message


def test_consumer_which_supplies_the_argument_not_flagged():
    assert _capability_findings("""
        <parameters>
          <task value="nodal"><haloMassFunction value="pressSchechter"/></task>
        </parameters>""") == []


def test_requirement_through_a_forwarding_implementation_flagged():
    assert len(_capability_findings("""
        <parameters>
          <task value="integrator">
            <haloMassFunction value="multiplier"><haloMassFunction value="pressSchechter"/></haloMassFunction>
          </task>
        </parameters>""")) == 1
    assert _capability_findings("""
        <parameters>
          <task value="integrator">
            <haloMassFunction value="multiplier"><haloMassFunction value="tinker2008"/></haloMassFunction>
          </task>
        </parameters>""") == []


def test_requirement_of_a_hoisted_object_flagged():
    # The mass function is not inside the task, so objectBuilder finds it in the enclosing element.
    assert len(_capability_findings("""
        <parameters>
          <haloMassFunction value="pressSchechter"/>
          <task value="integrator"/>
        </parameters>""")) == 1


def test_requirement_of_a_default_object_flagged():
    assert len(_capability_findings("""
        <parameters>
          <task value="integrator"/>
        </parameters>""")) == 1


def test_requirement_of_a_default_consumer_checked():
    # The task is itself defaulted at the root (to `nodal`, which supplies a node), and a default
    # multiplier is not selectable here, so nothing is reported - but the walk must not fail.
    assert _capability_findings('<parameters><haloMassFunction value="tinker2008"/></parameters>') == []


def test_requirement_through_an_id_reference_flagged():
    assert len(_capability_findings("""
        <parameters>
          <shared><haloMassFunction id="massFunction" value="pressSchechter"/></shared>
          <task value="integrator"><haloMassFunction idRef="massFunction"/></task>
        </parameters>""")) == 1


def test_self_referencing_forward_terminates():
    # A multiplier with no inner mass function finds itself when objectBuilder walks up - a recursive
    # build which Galacticus rejects at run time. The validator must simply not loop.
    assert _capability_findings("""
        <parameters>
          <task value="integrator"><haloMassFunction value="multiplier"/></task>
        </parameters>""") == []
