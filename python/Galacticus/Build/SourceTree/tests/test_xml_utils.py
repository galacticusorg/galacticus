"""Regression tests for `XML.Utils.xml_to_dict`.

Bug: the original port called `.strip()` on every element's text without
distinguishing single-line from multi-line content.  For directive bodies
carrying Fortran code — `<forEach>…end if\n</forEach>`,
`<call>…\n</call>` — that stripped the trailing newline of the last line.
Downstream emitters concatenate the body with their own `end do\n` /
`end if\n` lines, producing invalid Fortran like `end ifend do`.

Fix: keep multi-line text as-is (with its leading/trailing whitespace,
including the trailing newline) and only treat it as empty when it is
pure whitespace.  Single-line text continues to be stripped — many
directives have human-formatted text content (descriptions, parameter
values) that callers expect normalised.
"""

import xml.etree.ElementTree as ET
from XML.Utils import xml_to_dict


def _parse(xml_text):
    return xml_to_dict(ET.fromstring(xml_text))


def test_multiline_directive_body_keeps_trailing_newline():
    """A `<forEach>` body whose last line is `end if\\n` must NOT have its
    trailing newline stripped."""
    body = "\n  if (.not.passed{index}) then\n  end if\n"
    out = _parse(f'<forEach variable="passed">{body}</forEach>')
    assert isinstance(out, dict)
    assert out['variable'] == 'passed'
    assert out['content'].endswith('\n'), repr(out['content'])
    assert 'end if\n' in out['content']


def test_single_line_text_still_stripped():
    """Single-line text continues to be stripped — preserving the legacy
    behaviour callers rely on."""
    out = _parse('<description>  Foo bar.  </description>')
    assert out == 'Foo bar.'


def test_whitespace_only_multiline_body_treated_as_empty():
    """A multi-line element whose text is all whitespace returns just the
    attribute dict (no `content` key)."""
    out = _parse('<forEach variable="x">\n   \n   \n</forEach>')
    assert out == {'variable': 'x'}


def test_self_closing_element_has_no_content():
    out = _parse('<reference url="http://example.com"/>')
    assert out == {'url': 'http://example.com'}


def test_multiline_text_only_element_returns_string():
    """An element with no attributes/children but multi-line text returns
    the text as-is (with newlines), not a dict."""
    text = "line one\nline two\n"
    out = _parse(f'<note>{text}</note>')
    assert out == text
