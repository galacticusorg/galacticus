"""Tests for `Implementation_Builder`.

FoX's `extractDataContent` writes a two-line message to standard error and then
`stop`s - which exits with status zero - unless it is passed an `iostat`
argument.  A fully-specified merger tree containing one malformed property
therefore used to terminate the run silently, having built no tree and evolved
nothing (issue #1512).  Every generated call must pass `iostat` and report a
non-zero status itself.
"""

from Galacticus.Build.Components.Implementations.CreateDestroy import (
    Implementation_Builder,
)


def _content(properties, *, class_name='satellite', impl_name='orbiting'):
    build  = {}
    member = {'name': impl_name, 'properties': {'property': properties}}
    Implementation_Builder(build, {'name': class_name}, member)
    bound  = build['types']['nodeComponent' + 'Satellite' + 'Orbiting']['boundFunctions']
    assert len(bound) == 1
    return bound[0]['descriptor']


def _property(name, type_, rank):
    return {'name': name, 'data': {'type': type_, 'rank': rank}, 'attributes': {}}


def test_a_status_variable_is_declared():
    descriptor = _content([_property('position', 'double', 1)])
    declared   = [
        variable
        for entry in descriptor['variables']
        for variable in entry['variables']
    ]
    assert 'status' in declared


def test_rank_one_extraction_requests_and_checks_a_status():
    descriptor = _content([_property('position', 'double', 1)])
    content    = descriptor['content']
    assert "call extractDataContent(property,self%positionData(i),iostat=status)" in content
    assert "if (status /= 0) call Node_Component_Builder_Error("                  in content
    assert "'satellite','orbiting','position',1,propertyListLength,status"        in content


def test_rank_zero_extraction_requests_and_checks_a_status():
    descriptor = _content([_property('boundMass', 'double', 0)])
    content    = descriptor['content']
    assert "call extractDataContent(property,self%boundMassData,iostat=status)"    in content
    assert "'satellite','orbiting','boundMass',0,propertyListLength,status"        in content


def test_a_scalar_given_several_elements_is_reported_with_its_name():
    """The count check used to report only 'scalar property must have precisely
    one value', naming neither the property nor the node."""
    content = _content([_property('boundMass', 'double', 0)])['content']
    assert "if (propertyListLength > 1) call Node_Component_Builder_Error(" \
           "'satellite','orbiting','boundMass',0,propertyListLength,0,self%hostNode)" in content


def test_every_extraction_in_a_builder_passes_a_status():
    """No `extractDataContent` call may be left without an `iostat` argument."""
    content = _content([
        _property('position' , 'double' , 1),
        _property('boundMass', 'double' , 0),
        _property('indexHost', 'integer', 0),
    ])['content']
    for line in content.split("\n"):
        if 'extractDataContent(' in line:
            assert 'iostat=status' in line, line
