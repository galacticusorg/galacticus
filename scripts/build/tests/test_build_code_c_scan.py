"""Regression test for buildCode.py's per-line read of C/C++ sources.

`_process_until_include_or_eof` previously read C/C++ files with
`(fh.readline(), fh.readline(), '')`, consuming TWO lines per iteration —
so `raw_line` and `processed_line` were different, consecutive lines and
every other line was skipped. A directive block spanning three consecutive
lines (`!![` / `<tag .../>` / `!!]`) was therefore missed entirely. Latent
in the current tree (the only include directive, `component`, lives in
Fortran), but a real correctness bug for any C/C++ source carrying a
directive.
"""

import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir))
import buildCode as bc  # noqa: E402


def _scan_c(tmp_path, content, directive='myDirective'):
    src = tmp_path / 'probe.c'
    src.write_text(content)
    build = {'codeType': 'c', 'directive': directive,
             'moduleName': '', 'currentFileName': str(src)}
    frame = {'name': str(src), 'position': -1, 'in_xml': False}
    entry = {'files': [str(src)], 'directives': []}
    with open(src) as fh:
        bc._process_until_include_or_eof(fh, frame, build, '.', entry, {})
    return entry['directives']


def test_c_directive_is_found(tmp_path):
    # The directive occupies the three lines that a two-lines-per-iteration
    # read would straddle; a correct single-line read captures it.
    directives = _scan_c(tmp_path,
                         'int x;\n'
                         '!![\n'
                         '<myDirective foo="bar"/>\n'
                         '!!]\n'
                         'int y;\n')
    assert len(directives) == 1
    assert directives[0]['directive'] == {'foo': 'bar'}


def test_c_directive_after_odd_offset_is_found(tmp_path):
    # An extra leading line shifts the directive's parity; it must still be
    # found (the buggy double-read only ever saw one parity of lines).
    directives = _scan_c(tmp_path,
                         'int a;\n'
                         'int b;\n'
                         '!![\n'
                         '<myDirective n="2"/>\n'
                         '!!]\n')
    assert len(directives) == 1
    assert directives[0]['directive'] == {'n': '2'}


def test_c_non_matching_directive_ignored(tmp_path):
    # A directive whose tag is not the one being built is not recorded
    # (confirms the scan still discriminates by tag, not that it records
    # everything it happens to see).
    directives = _scan_c(tmp_path,
                         '!![\n'
                         '<otherThing/>\n'
                         '!!]\n')
    assert directives == []
