# Regression test for `_insert_code_tree` in
# `Galacticus.Build.SourceTree.Process.Enumeration`.
#
# Bug: `_insert_code_tree` parsed its synthetic source text with a
# fixed `name='Enumeration'`, so every node inserted by the enumeration
# emitter (the `enumeration<N>Type` derived type, the
# `enumeration<N>IsEqual` function, validators, encoders, decoders, …)
# carried the literal string `'Enumeration'` in its `source` field.
#
# That string is not a real `.F90` filename, so the
# source-file-skip filter in `classDocumentation` (introduced to keep
# child class types out of the parent's `<basename>.classes.xml`) did
# NOT skip these nodes.  When FunctionClass later threaded a child class
# file's enumerations into a parent's `<module>` interfaces list, the
# parent's `classDocumentation` postprocess pass walked the inserted
# `enumeration<N>Type` even though its `enumeration<N>IsEqual` function
# body had stayed in the child's submodule file — `_process_function`
# couldn't bind, and the doc consumer warned
# "missing function type for method 'operator(==)' of class
# 'enumeration<N>Type'" (treeStatistic, wagner2016*, wittGordon2000*, …).
#
# Fix: parse the synthetic enumeration source text with the OUTER
# source's name, so its nodes inherit the outer file's source field
# and the existing classDocumentation skip works correctly.

import xml.etree.ElementTree as ET
from Galacticus.Build.SourceTree.Process.Enumeration import _insert_code_tree


def _outer_module_node(source):
    """Build a minimal `module` parent whose `source` attribute matches
    `source` — the shape `_insert_code_tree` is called against in the
    real pipeline.  We use a self-closing `contains` marker so the
    `insert_post_contains` inserter has somewhere to anchor."""
    contains = {
        'type':       'contains',
        'firstChild': None,
        'sibling':    None,
        'parent':     None,
        'source':     source,
    }
    parent = {
        'type':       'module',
        'name':       'foo',
        'firstChild': contains,
        'sibling':    None,
        'parent':     None,
        'source':     source,
    }
    contains['parent'] = parent
    return parent


def test_insert_code_tree_inherits_parents_source():
    """Inserted nodes carry the OUTER source name, not 'Enumeration'."""
    parent = _outer_module_node('utility.foo.F90')
    _insert_code_tree(
        parent,
        '\n  type, public :: barType\n  end type barType\n',
        inserter=lambda p, kids: setattr(_inserted_holder, 'kids', kids),
    )
    inserted = _inserted_holder.kids
    assert inserted, "inserter callback didn't capture children"
    assert all(node.get('source') == 'utility.foo.F90' for node in inserted), \
        [node.get('source') for node in inserted]


def test_insert_code_tree_falls_back_when_parent_has_no_source():
    """If the parent has no `source` attribute (defensive), keep the
    legacy `'Enumeration'` placeholder rather than `None`."""
    parent = {
        'type':       'module',
        'name':       'foo',
        'firstChild': None,
        'sibling':    None,
        'parent':     None,
    }
    _insert_code_tree(
        parent,
        '\n  type :: barType\n  end type barType\n',
        inserter=lambda p, kids: setattr(_inserted_holder, 'kids', kids),
    )
    assert all(node.get('source') == 'Enumeration'
               for node in _inserted_holder.kids)


class _inserted_holder:
    """Module-level container so the test's `inserter` callback can
    stash its `kids` argument for the assert step."""
    kids = ()


# ---------------------------------------------------------------------------
# Documentation fragment emission
# ---------------------------------------------------------------------------
#
# The Perl original wrote a `doc/enumerations/definitions/<name>.tex`
# fragment for every `<enumeration>` directive — assembled later by
# `buildDocumentation.sh`'s `ls enumerations/definitions/*.tex` sweep
# into the docs build.  The Python port had a TODO comment but no
# emission, so the docs build assembled an empty list.

from Galacticus.Build.SourceTree.Process.Enumeration import (
    _emit_enumeration_definition_tex,
)


def _enclosed_directive_node(name, description, entries, module_name='Foo',
                              file_name='utility.foo.F90'):
    """Build a tiny `file > module > <enumeration>` chain so the helper
    can find the enclosing module and file via parent walks."""
    file_node = {
        'type': 'file', 'name': file_name, 'parent': None,
        'firstChild': None, 'sibling': None, 'source': file_name,
    }
    module_node = {
        'type': 'module', 'name': module_name, 'parent': file_node,
        'firstChild': None, 'sibling': None, 'source': file_name,
    }
    file_node['firstChild'] = module_node
    enum_node = {
        'type':       'enumeration',
        'directive':  {'name': name, 'description': description,
                       'entry': entries},
        'parent':     module_node,
        'firstChild': None, 'sibling': None,
        'source':     file_name,
    }
    module_node['firstChild'] = enum_node
    return enum_node


def test_enumeration_definition_tex_written(tmp_path, monkeypatch):
    """Each `<enumeration>` directive produces a
    `doc/enumerations/definitions/<name>.tex` fragment with a
    description line and a per-member table."""
    monkeypatch.chdir(tmp_path)

    enum_node = _enclosed_directive_node(
        name='treeStatistic',
        description='Enumeration of tree-walk statistics.',
        entries=[{'label': 'nodeCount'}, {'label': 'endNodeCount'}],
        module_name='Merger_Trees_Operators_Augment',
        file_name='merger_trees.operators.augment.F90',
    )

    _emit_enumeration_definition_tex(
        enum_node, 'treeStatistic',
        'Enumeration of tree-walk statistics.',
        [{'label': 'nodeCount'}, {'label': 'endNodeCount'}])

    out_path = tmp_path / 'doc' / 'enumerations' / 'definitions' \
        / 'treeStatistic.tex'
    assert out_path.exists(), list((tmp_path / 'doc').rglob('*'))

    body = out_path.read_text()
    assert '\\subsection{\\large \\mono{treeStatistic}}' in body
    assert 'Description: & Enumeration of tree-walk statistics.' in body
    assert '\\mono{treeStatisticNodeCount}' in body
    assert '\\mono{treeStatisticEndNodeCount}' in body
    # The cross-reference line is present when an enclosing module is found.
    assert 'Provided by:' in body
    assert 'merger_trees_operators_augment' in body   # encoded basename
    assert 'Merger\\_Trees\\_Operators\\_Augment' in body   # encoded module name


def test_enumeration_definition_tex_without_module(tmp_path, monkeypatch):
    """Enumeration declared at file scope (no enclosing `module`)
    still gets a fragment, just without the cross-reference line."""
    monkeypatch.chdir(tmp_path)

    file_node = {
        'type': 'file', 'name': 'foo.F90', 'parent': None,
        'firstChild': None, 'sibling': None,
    }
    enum_node = {
        'type':       'enumeration',
        'directive':  {'name': 'top', 'entry': [{'label': 'a'}]},
        'parent':     file_node,
        'firstChild': None, 'sibling': None,
    }
    file_node['firstChild'] = enum_node

    _emit_enumeration_definition_tex(
        enum_node, 'top', 'A top-level enum.', [{'label': 'a'}])

    out = (tmp_path / 'doc' / 'enumerations' / 'definitions' / 'top.tex'
           ).read_text()
    assert 'Description: & A top-level enum.' in out
    assert 'Provided by:' not in out
    assert '\\mono{topA}' in out


def test_enumeration_latex_special_characters_escaped(tmp_path, monkeypatch):
    """The label is run through `latex_encode` before substitution so
    underscores and other LaTeX specials don't corrupt the output."""
    monkeypatch.chdir(tmp_path)

    enum_node = _enclosed_directive_node(
        name='underscores',
        description='An & test_with_underscores.',
        entries=[{'label': 'foo_bar'}],
    )
    _emit_enumeration_definition_tex(
        enum_node, 'underscores',
        'An & test_with_underscores.',
        [{'label': 'foo_bar'}])

    body = (tmp_path / 'doc' / 'enumerations' / 'definitions'
            / 'underscores.tex').read_text()
    # The LABEL goes through latex_encode (escapes `_` → `\_`).
    assert '\\mono{underscoresFoo\\_bar}' in body
