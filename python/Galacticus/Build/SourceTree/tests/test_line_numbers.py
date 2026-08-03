"""Regression tests for original-source line-number fidelity.

Bug class: `{introspection:location}` expansions and the `.lmap`
line-number mappings written by preprocess.py must report positions in
the ORIGINAL (unpreprocessed) source.  An earlier draft of the Python
port assigned every nested node the line number of its enclosing unit
(the Perl implementation tracked absolute lines through the recursion),
so runtime error messages named the wrong line and postprocess.py's
compile-diagnostic remapping drifted — in the worst case (missing or
useless map) printing the preprocessed line number in both slots of
"line U [preprocessed line P]".

Also covered here: the file node's `name` must be the path relative to
the `source/` directory — the source tree is hierarchical, so a bare
basename (e.g. `_class.F90`) is ambiguous.
"""

import bisect
import os
import re

from Galacticus.Build.SourceTree import (
    parse_code, parse_file, serialize, walk_tree,
)


_SOURCE = (
    "module foo\n"                                        # line 1
    "  use :: Some_Module\n"                              # line 2
    "  implicit none\n"                                   # line 3
    "  integer :: moduleVariable\n"                       # line 4
    "\n"                                                  # line 5
    "contains\n"                                          # line 6
    "\n"                                                  # line 7
    "  subroutine bar()\n"                                # line 8
    "    implicit none\n"                                 # line 9
    "    integer :: localVariable\n"                      # line 10
    "    !![\n"                                           # line 11
    "    <someDirective name=\"example\"/>\n"             # line 12
    "    !!]\n"                                           # line 13
    "    localVariable=1\n"                               # line 14
    "    call baz()\n"                                    # line 15
    "  end subroutine bar\n"                              # line 16
    "\n"                                                  # line 17
    "  subroutine baz()\n"                                # line 18
    "    implicit none\n"                                 # line 19
    "    write (0,*) 'hello'\n"                           # line 20
    "  end subroutine baz\n"                              # line 21
    "end module foo\n"                                    # line 22
)


def _nodes_by_type(tree, node_type):
    return [n for n in walk_tree(tree) if n.get('type') == node_type]


def test_node_lines_are_absolute_and_one_based():
    tree = parse_code(_SOURCE, name='test.F90', instrument=False)

    module = _nodes_by_type(tree, 'module')[0]
    assert module['line'] == 1

    subroutines = {n['name']: n for n in _nodes_by_type(tree, 'subroutine')}
    assert subroutines['bar']['line'] == 8
    assert subroutines['baz']['line'] == 18

    directive = _nodes_by_type(tree, 'someDirective')[0]
    assert directive['line'] == 12


def test_declaration_lines_are_absolute():
    tree = parse_code(_SOURCE, name='test.F90', instrument=False)
    decl_lines = sorted(
        d['line']
        for n in _nodes_by_type(tree, 'declaration')
        for d in n.get('declarations', [])
    )
    assert decl_lines == [4, 10]


def test_parse_file_name_is_relative_to_source_directory(tmp_path):
    subdir = tmp_path / 'source' / 'intergalactic_medium' / 'state'
    subdir.mkdir(parents=True)
    path = subdir / 'internal.F90'
    path.write_text("module foo\nend module foo\n")
    tree = parse_file(str(path))
    assert tree['name'] == 'intergalactic_medium/state/internal.F90'
    assert tree['source'] == str(path)


def test_parse_file_name_falls_back_to_basename_outside_source(tmp_path):
    path = tmp_path / 'generated.Inc'
    path.write_text("module foo\nend module foo\n")
    tree = parse_file(str(path))
    assert tree['name'] == 'generated.Inc'


def test_lmap_remaps_preprocessed_lines_to_original():
    """Round-trip the postprocess.py remap arithmetic: after inserting
    generated code into the tree, every surviving original line must remap
    from its (shifted) output position back to its original line number."""
    from Galacticus.Build.SourceTree import insert_after_node, _make_code_node

    tree = parse_code(_SOURCE, name='test.F90', source='source/test.F90',
                      instrument=False)

    # Simulate a process hook: inject three generated lines ahead of the
    # first code node inside subroutine bar (source tag ends in '()', the
    # auto-generated marker postprocess.py recognises).
    bar = [n for n in walk_tree(tree) if n.get('name') == 'bar'][0]
    generated = _make_code_node("call generated_1()\ncall generated_2()\n"
                                "call generated_3()\n", 'Test.hook()', 1)
    insert_after_node(bar['firstChild'], [generated])

    stripped, mappings = serialize(tree, annotate=True, strip_mappings=True)

    # Parse the mappings exactly as postprocess.py does.
    lmap = [{'source': 'test.p.F90', 'line': 1, 'lineOriginal': 1}]
    for line in mappings.splitlines():
        m = re.match(r'^!-->\s+(\d+)\s+(\d+)\s+"(.+)"\s*$', line)
        assert m, f"malformed mapping line: {line!r}"
        lmap.append({'source': m.group(3), 'line': int(m.group(1)),
                     'lineOriginal': int(m.group(2))})
    keys = [e['lineOriginal'] for e in lmap]
    assert keys == sorted(keys), "mapping anchors must be ascending"

    original_lines = _SOURCE.splitlines()
    for out_number, text in enumerate(stripped.splitlines(), start=1):
        idx   = bisect.bisect_right(keys, out_number) - 1
        entry = lmap[idx]
        if entry['source'].endswith('()'):
            assert 'generated' in text
            continue
        remapped = out_number - entry['lineOriginal'] + entry['line']
        # Embedded XML body lines gain a `!< ` comment prefix in the output
        # (_comment_embedded) — strip it before comparing against the
        # original text.
        text = re.sub(r'^!<\s?', '', text)
        assert original_lines[remapped - 1] == text, (
            f"output line {out_number} ({text!r}) remapped to original "
            f"line {remapped} ({original_lines[remapped - 1]!r})")
