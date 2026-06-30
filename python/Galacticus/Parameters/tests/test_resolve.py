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


def test_changes_applied_in_order(tmp_path):
    # Second change sees the first change's effect (append then update).
    root = _apply('<parameters><a><x value="1"/></a></parameters>',
                  '<changes>'
                  '<change type="append" path="a"><y value="2"/></change>'
                  '<change type="update" path="a/y" value="3"/>'
                  '</changes>', tmp_path)
    assert root[0][1].tag == 'y' and root[0][1].get('value') == '3'
