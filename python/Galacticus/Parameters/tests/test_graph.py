"""Tests for `Galacticus.Parameters.graph` (the object graph of a model)."""

from lxml import etree

from Galacticus.Parameters.graph import _scan_declared_objects, build_graph, group_of


def _object(class_name, parameter_name=None, repeatable=False):
    return {'class': class_name, 'parameterName': parameter_name or class_name,
            'repeatable': repeatable, 'source': 'parameters', 'sourceElement': None}


def _implementation(class_name, label, objects=(), parameters=(), source_file='tasks/x.F90'):
    return {'functionClass': class_name, 'label': label, 'alias': None,
            'parent': class_name + 'Class', 'objects': list(objects),
            'parameters': list(parameters), 'sourceFile': source_file}


CATALOG = {
    'enumerations': {},
    'functionClasses': {
        'task':      {'default': 'evolve', 'descriptiveName': 'Tasks'},
        'operator':  {'default': 'null',   'descriptiveName': 'Operators'},
        'cosmology': {'default': 'simple', 'descriptiveName': 'Cosmology'},
        'unused':    {'default': 'plain',  'descriptiveName': 'Unused'},
    },
    'implementations': {
        'taskEvolve':      _implementation('task', 'evolve', [_object('operator'), _object('cosmology')]),
        'operatorNull':    _implementation('operator', 'null'),
        'operatorMulti':   _implementation('operator', 'multi', [_object('operator', repeatable=True)]),
        'operatorGrowth':  _implementation(
            'operator', 'growth', [_object('cosmology')],
            [{'name': 'rate', 'default': '1.0d0', 'description': 'The :math:`rate`.'}],
            source_file='nodes/operators/physics/star_formation/disks.F90'),
        'cosmologySimple': _implementation('cosmology', 'simple', source_file='cosmology/functions/simple.F90'),
        'cosmologyOther':  _implementation('cosmology', 'other', source_file='cosmology/functions/other.F90'),
        'unusedPlain':     _implementation('unused', 'plain'),
    },
}


def _graph(xml):
    return build_graph(etree.fromstring(xml.encode()), CATALOG, name='test')


def _node(graph, class_name, label):
    found = [node for node in graph['nodes'] if (node['class'], node['label']) == (class_name, label)]
    assert len(found) == 1, found
    return found[0]


def _uses(graph, node):
    by_id = {other['id']: other for other in graph['nodes']}
    return sorted((by_id[edge['to']]['class'], by_id[edge['to']]['label'])
                  for edge in graph['edges'] if edge['from'] == node['id'])


def test_defaults_are_built_for_what_the_file_leaves_out():
    graph = _graph('<parameters/>')
    assert sorted((node['class'], node['label'], node['origin']) for node in graph['nodes']) == [
        ('cosmology', 'simple', 'default'), ('operator', 'null', 'default'), ('task', 'evolve', 'default')]
    assert _uses(graph, _node(graph, 'task', 'evolve')) == [('cosmology', 'simple'), ('operator', 'null')]


def test_nearest_enclosing_element_is_used_and_shared():
    graph = _graph('''<parameters>
      <cosmology value="other"/>
      <operator value="multi">
        <operator value="growth"><rate value="2.0"/></operator>
        <operator value="growth"><cosmology value="simple"/></operator>
      </operator>
    </parameters>''')
    growths = [node for node in graph['nodes'] if node['label'] == 'growth']
    assert len(growths) == 2
    # The first finds the root-level cosmology, which the task shares; the
    # second has one of its own.
    assert _uses(graph, growths[0]) == [('cosmology', 'other')]
    assert _uses(graph, growths[1]) == [('cosmology', 'simple')]
    other = _node(graph, 'cosmology', 'other')
    assert sum(1 for edge in graph['edges'] if edge['to'] == other['id']) == 2
    assert _uses(graph, _node(graph, 'operator', 'multi')) == [('operator', 'growth'), ('operator', 'growth')]
    assert growths[0]['parameters'] == [
        {'name': 'rate', 'value': '2.0', 'set': True, 'description': 'The rate.'}]
    assert growths[1]['parameters'][0]['value'] == '1.0e0' and not growths[1]['parameters'][0]['set']


def test_reference_resolves_to_the_one_object():
    graph = _graph('''<parameters>
      <operator value="multi">
        <operator value="growth"><cosmology value="other" id="c"/></operator>
        <operator value="growth"><cosmology idRef="c"/></operator>
      </operator>
    </parameters>''')
    other = _node(graph, 'cosmology', 'other')
    assert sum(1 for edge in graph['edges'] if edge['to'] == other['id']) == 2


def test_object_which_nothing_builds_is_kept_but_not_reached():
    graph = _graph('<parameters><unused value="plain"/></parameters>')
    assert _node(graph, 'unused', 'plain')['reached'] is False
    assert _node(graph, 'task', 'evolve')['reached'] is True


def test_group_follows_the_source_file():
    assert group_of('cosmology/functions/matter_lambda.F90') == 'cosmology'
    assert group_of('./source/nodes/operators/physics/star_formation/disks.F90') == 'starFormation'
    assert group_of('nodes/operators/multi.F90') == 'engine'
    assert group_of('merger_trees/evolver/standard.F90') == 'engine'
    assert group_of('merger_trees/construct/build.F90') == 'trees'
    assert group_of('numerical/random_numbers/GSL.F90') == 'numerics'
    graph = _graph('<parameters><operator value="growth"/></parameters>')
    assert _node(graph, 'operator', 'growth')['group'] == 'starFormation'


def test_object_built_by_hand_is_read_from_allowed_parameters(tmp_path):
    (tmp_path / 'tasks').mkdir()
    (tmp_path / 'tasks' / 'x.F90').write_text(
        "  !![\n"
        "  <task name=\"taskEvolve\">\n"
        "   <allowedParameters>\n"
        "    <object variable=\"unused_\" name=\"unused\" source=\"parameters\" function=\"f\" module=\"Functions_Global\"/>\n"
        "   </allowedParameters>\n"
        "  </task>\n"
        "  !!]\n")
    declared = _scan_declared_objects(str(tmp_path), CATALOG)
    assert [obj['class'] for obj in declared['taskEvolve']] == ['unused']
    root = etree.fromstring(b'<parameters><unused value="plain"/></parameters>')
    graph = build_graph(root, CATALOG, source_root=str(tmp_path))
    assert _node(graph, 'unused', 'plain')['reached'] is True
    assert ('unused', 'plain') in _uses(graph, _node(graph, 'task', 'evolve'))
