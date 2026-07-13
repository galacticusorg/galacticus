"""Tests for Galacticus.Parameters.resolve (parameter-file resolver, Stage 1)."""

import os

import pytest

from lxml import etree

from Galacticus.Parameters import resolve
from Galacticus.Parameters.resolve import ResolveError


def _tree(xml):
    return etree.ElementTree(etree.fromstring(xml.encode()))


def _root(xml):
    return etree.fromstring(xml.encode())


def _tags(element):
    return [child.tag for child in element if isinstance(child.tag, str)]


# --- path resolution (Stage 0) ---------------------------------------------

def test_resolve_path_legacy_bare_is_absolute_from_root():
    root = _root('<parameters><a><b value="1"/></a><c><b value="2"/></c></parameters>')
    a = root[0]
    # A bare path is absolute-from-root even when evaluated from a nested node.
    matches = resolve.resolve_path(a, 'c/b', root)
    assert [m.get('value') for m in matches] == ['2']


def test_resolve_path_relative_and_self_and_parent():
    root = _root('<parameters><a><b value="1"/></a></parameters>')
    a = root[0]
    assert resolve.resolve_path(a, '.', root) == [a]
    assert resolve.resolve_path(a, './b', root)[0].get('value') == '1'
    assert resolve.resolve_path(a[0], '..', root) == [a]


def test_resolve_path_predicate_and_empty():
    root = _root('<parameters><n value="x"/><n value="y"/></parameters>')
    assert resolve.resolve_path(root, '', root) == [root]
    assert resolve.resolve_path(root, "n[@value='y']", root)[0].get('value') == 'y'


# --- XInclude (Stage 1) -----------------------------------------------------

def _write(tmp_path, name, xml):
    path = tmp_path / name
    path.write_text(xml)
    return str(path)


def test_xinclude_children_of_parameters_root(tmp_path):
    _write(tmp_path, 'cosmology.xml',
           '<parameters><Hubble value="70"/><Omega value="0.3"/></parameters>')
    main = _write(
        tmp_path, 'main.xml',
        '<parameters xmlns:xi="http://www.w3.org/2001/XInclude">'
        '<first value="1"/>'
        '<xi:include href="cosmology.xml" xpointer="xpointer(parameters/*)"/>'
        '<last value="9"/></parameters>')
    root = resolve.resolve_file(main).getroot()
    # Included children are spliced in place, preserving order.
    assert _tags(root) == ['first', 'Hubble', 'Omega', 'last']


def test_xinclude_nested(tmp_path):
    _write(tmp_path, 'inner.xml', '<parameters><deep value="1"/></parameters>')
    _write(tmp_path, 'outer.xml',
           '<parameters xmlns:xi="http://www.w3.org/2001/XInclude">'
           '<mid value="2"/>'
           '<xi:include href="inner.xml" xpointer="xpointer(parameters/*)"/>'
           '</parameters>')
    main = _write(tmp_path, 'main.xml',
                  '<parameters xmlns:xi="http://www.w3.org/2001/XInclude">'
                  '<xi:include href="outer.xml" xpointer="xpointer(parameters/*)"/>'
                  '</parameters>')
    root = resolve.resolve_file(main).getroot()
    assert _tags(root) == ['mid', 'deep']


def test_xinclude_missing_href_raises(tmp_path):
    main = _write(tmp_path, 'main.xml',
                  '<parameters xmlns:xi="http://www.w3.org/2001/XInclude">'
                  '<xi:include href="nope.xml" xpointer="xpointer(parameters/*)"/>'
                  '</parameters>')
    with pytest.raises(ResolveError, match="does not resolve"):
        resolve.resolve_file(main)


def test_xinclude_fallback_used_when_missing(tmp_path):
    main = _write(tmp_path, 'main.xml',
                  '<parameters xmlns:xi="http://www.w3.org/2001/XInclude">'
                  '<xi:include href="nope.xml">'
                  '<xi:fallback><backup value="1"/></xi:fallback>'
                  '</xi:include></parameters>')
    root = resolve.resolve_file(main).getroot()
    assert _tags(root) == ['backup']


# --- changes (Stage 1): one test per operation ------------------------------

def _apply(base_xml, changes_xml, tmp_path):
    """Apply an inline <changes> document to an inline <parameters> document."""
    change_file = _write(tmp_path, 'changes.xml', changes_xml)
    tree = _tree(base_xml)
    resolve.apply_changes(tree.getroot(), [change_file])
    return tree.getroot()


def test_change_remove(tmp_path):
    root = _apply('<parameters><a value="1"/><b value="2"/></parameters>',
                  '<changes><change type="remove" path="a"/></changes>', tmp_path)
    assert _tags(root) == ['b']


def test_change_update(tmp_path):
    root = _apply('<parameters><a value="1"/></parameters>',
                  '<changes><change type="update" path="a" value="9"/></changes>',
                  tmp_path)
    assert root[0].get('value') == '9'


def test_change_update_append(tmp_path):
    root = _apply('<parameters><a value="1 2"/></parameters>',
                  '<changes><change type="update" path="a" value=" 3" '
                  'append="true"/></changes>', tmp_path)
    assert root[0].get('value') == '1 2 3'


def test_change_update_no_value_on_target_raises(tmp_path):
    with pytest.raises(ResolveError, match="no `value`"):
        _apply('<parameters><a><child value="1"/></a></parameters>',
               '<changes><change type="update" path="a" value="9"/></changes>',
               tmp_path)


def test_change_append(tmp_path):
    root = _apply('<parameters><a><x value="1"/></a></parameters>',
                  '<changes><change type="append" path="a">'
                  '<y value="2"/></change></changes>', tmp_path)
    assert _tags(root[0]) == ['x', 'y']


def test_change_insert_before_and_after(tmp_path):
    root = _apply('<parameters><a value="1"/><b value="2"/></parameters>',
                  '<changes>'
                  '<change type="insertBefore" path="b"><before value="0"/></change>'
                  '<change type="insertAfter" path="b"><after value="3"/></change>'
                  '</changes>', tmp_path)
    assert _tags(root) == ['a', 'before', 'b', 'after']


def test_change_replace(tmp_path):
    root = _apply('<parameters><a value="1"/></parameters>',
                  '<changes><change type="replace" path="a">'
                  '<b value="2"/><c value="3"/></change></changes>', tmp_path)
    assert _tags(root) == ['b', 'c']


def test_change_replace_or_append_existing_replaces(tmp_path):
    root = _apply('<parameters><a value="1"/></parameters>',
                  '<changes><change type="replaceOrAppend" path="a">'
                  '<a value="2"/></change></changes>', tmp_path)
    assert _tags(root) == ['a'] and root[0].get('value') == '2'


def test_change_replace_or_append_missing_appends_to_parent(tmp_path):
    root = _apply('<parameters><a><x value="1"/></a></parameters>',
                  '<changes><change type="replaceOrAppend" path="a/y">'
                  '<y value="2"/></change></changes>', tmp_path)
    assert _tags(root[0]) == ['x', 'y']


def test_change_replace_with(tmp_path):
    root = _apply(
        '<parameters><src value="keep"><deep value="d"/></src>'
        '<dst value="drop"/></parameters>',
        '<changes><change type="replaceWith" path="dst" target="src"/></changes>',
        tmp_path)
    # dst is replaced by a deep clone of src.
    assert _tags(root) == ['src', 'src']
    assert root[1].get('value') == 'keep' and _tags(root[1]) == ['deep']


def test_change_encapsulate(tmp_path):
    root = _apply('<parameters><t value="inner"/></parameters>',
                  '<changes><change type="encapsulate" path="t">'
                  '<wrapper value="outer"/></change></changes>', tmp_path)
    assert _tags(root) == ['wrapper']
    wrapper = root[0]
    assert wrapper.get('value') == 'outer' and _tags(wrapper) == ['t']
    assert wrapper[0].get('value') == 'inner'


def test_change_predicate_path(tmp_path):
    root = _apply(
        '<parameters><op value="a"/><op value="b"><inner value="1"/></op></parameters>',
        '<changes><change type="update" path="op[@value=\'b\']/inner" '
        'value="9"/></changes>', tmp_path)
    assert root[1][0].get('value') == '9'


def test_change_unknown_type_raises(tmp_path):
    with pytest.raises(ResolveError, match="unknown change type"):
        _apply('<parameters><a value="1"/></parameters>',
               '<changes><change type="bogus" path="a"/></changes>', tmp_path)


def test_change_missing_path_raises(tmp_path):
    with pytest.raises(ResolveError, match="does not exist"):
        _apply('<parameters><a value="1"/></parameters>',
               '<changes><change type="remove" path="nope"/></changes>', tmp_path)


def _resolve(xml):
    tree = _tree(xml)
    resolve.resolve_tree(tree, base_dir='.')
    return tree.getroot()


# --- references (Stage 2): validate, never dereference -----------------------

def test_check_references_valid():
    root = _root('<parameters><a value="1" id="k"/><a idRef="k"/></parameters>')
    assert resolve.check_references(root) == []


def test_check_references_dangling():
    root = _root('<parameters><a idRef="missing"/></parameters>')
    assert resolve.check_references(root) == [('a', 'missing')]


def test_check_references_requires_same_tag():
    # id 'k' is on <a>, but the idRef is on <b> -> no same-tag match.
    root = _root('<parameters><a id="k"/><b idRef="k"/></parameters>')
    assert resolve.check_references(root) == [('b', 'k')]


def test_resolve_file_dangling_idref_raises(tmp_path):
    main = _write(tmp_path, 'main.xml',
                  '<parameters><a idRef="nope"/></parameters>')
    with pytest.raises(ResolveError, match="no matching element"):
        resolve.resolve_file(main)


def test_references_are_not_dereferenced():
    # The resolved tree keeps id/idRef verbatim (runtime pointers, not copies).
    root = _resolve('<parameters><a value="1" id="k"/><a idRef="k"/></parameters>')
    assert root[0].get('id') == 'k'
    assert root[1].get('idRef') == 'k' and root[1].get('value') is None


# --- conditionals (Stage 2): evaluate + prune -------------------------------

def test_conditional_equal_kept_not_equal_pruned():
    root = _root(
        '<parameters>'
        '<cosmologicalMassVariance value="filteredPower"><sigma_8 value="0.912"/>'
        '</cosmologicalMassVariance>'
        '<active1 value="0.1" active="[cosmologicalMassVariance/sigma_8] != 0.912"/>'
        '<active1 value="0.2" active="[cosmologicalMassVariance/sigma_8] == 0.912"/>'
        '</parameters>')
    resolve.evaluate_conditionals(root)
    kept = [c for c in root if c.tag == 'active1']
    assert len(kept) == 1 and kept[0].get('value') == '0.2'
    assert kept[0].get('active') is None        # marker stripped from survivor


def test_conditional_relative_and_parent_path():
    root = _root('<parameters><wrap><flag value="on"/>'
                 '<dep value="1" active="[../flag] == on"/></wrap></parameters>')
    resolve.evaluate_conditionals(root)
    assert [c.tag for c in root[0]] == ['flag', 'dep']


def test_conditional_dereferences_idref_in_path():
    # Mirrors testsParameters.xml: the path lands on an idRef'd node and follows
    # it to the real definition to read the value.
    root = _resolve(
        '<parameters>'
        '<cosmologyParameters><cosmologicalMassVariance value="x" id="m">'
        '<sigma_8 value="0.9"/></cosmologicalMassVariance></cosmologyParameters>'
        '<cosmologicalMassVariance idRef="m"/>'
        '<active1 value="a" active="[cosmologicalMassVariance/sigma_8] == 0.9"/>'
        '</parameters>')
    assert any(c.tag == 'active1' for c in root)            # kept (0.9 == 0.9)
    assert any(c.get('idRef') == 'm' for c in root)         # idRef preserved


def test_conditional_fixed_point_chain():
    # `c` depends on `b` which depends on `base`; `c` precedes `b` in document
    # order, so a single pass cannot resolve it -- requires the fixed point.
    root = _root('<parameters><base value="go"/>'
                 '<c value="2" active="[b] == 1"/>'
                 '<b value="1" active="[base] == go"/></parameters>')
    resolve.evaluate_conditionals(root)
    assert {c.tag for c in root} == {'base', 'b', 'c'}      # all active, kept


def test_conditional_no_such_child_raises():
    root = _root('<parameters><x value="1" active="[nope] == 1"/></parameters>')
    with pytest.raises(ResolveError, match="no child parameter"):
        resolve.evaluate_conditionals(root)


def test_conditional_cyclic_raises():
    root = _root('<parameters>'
                 '<a value="1" active="[b] == 1"/>'
                 '<b value="1" active="[a] == 1"/></parameters>')
    with pytest.raises(ResolveError, match="cyclic"):
        resolve.evaluate_conditionals(root)


def test_no_resolver_state_leaks_into_output():
    root = _resolve('<parameters><base value="go"/>'
                    '<x value="1" active="[base] == go"/></parameters>')
    assert b'__resolver_active' not in resolve.to_bytes(root)


# --- corpus-grounded end-to-end (real test fixture) -------------------------

_REPO = os.path.abspath(os.path.join(
    os.path.dirname(__file__), os.pardir, os.pardir, os.pardir, os.pardir))
_TESTS_PARAMS = os.path.join(_REPO, 'testSuite', 'parameters', 'testsParameters.xml')


@pytest.mark.skipif(not os.path.exists(_TESTS_PARAMS),
                    reason="testsParameters.xml not found")
def test_resolve_real_tests_parameters():
    root = resolve.resolve_file(_TESTS_PARAMS).getroot()
    active1 = [c for c in root if c.tag == 'active1']
    # sigma_8 is 0.912: only the `== 0.912` variant (value 0.2) survives.
    assert [c.get('value') for c in active1] == ['0.2']
    # no conditional markers remain anywhere; the idRef anchor is preserved.
    assert not any(c.get('active') is not None for c in root.iter()
                   if isinstance(c.tag, str))
    assert any(c.get('idRef') == 'myCMV' for c in root.iter()
              if isinstance(c.tag, str))


def test_changes_applied_in_order(tmp_path):
    # Second change sees the first change's effect (append then update).
    root = _apply('<parameters><a><x value="1"/></a></parameters>',
                  '<changes>'
                  '<change type="append" path="a"><y value="2"/></change>'
                  '<change type="update" path="a/y" value="3"/>'
                  '</changes>', tmp_path)
    assert root[0][1].tag == 'y' and root[0][1].get('value') == '3'
