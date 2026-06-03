"""Regression tests for `Galacticus.Build.SourceTree.Process.ConditionalCall`.

Covers two bug classes:

  1. The XML parser strips trailing whitespace from `<call>` text bodies, so
     a template ending with `)` carries no trailing newline.  When
     `_render_call` returned that text verbatim the caller appended
     `end if\n` directly to the closing paren, producing `…)end if` — a
     Fortran syntax error.  The fix forces a trailing newline after each
     rendered call so `end if` lands on its own line.

  2. A `<call>` template with NO `{conditions}` marker should raise rather
     than silently emit an unguarded call (which would lose the entire
     point of the directive).
"""

import pytest

from Galacticus.Build.SourceTree.Process.ConditionalCall import _render_call


def test_render_call_appends_trailing_newline():
    """`<call>` bodies stripped of their trailing newline by the XML parser
    must still emit a newline after the rendered call line."""
    call_template = "call myRoutine(self{conditions})"   # no trailing \n
    arguments = {
        'a': {'name': 'a', 'value': 'aVal',  'conditionActual': 'isA'},
        'b': {'name': 'b', 'value': 'bVal',  'conditionActual': 'isB'},
    }
    condition_id = {'isA': 1, 'isB': 2}
    bits = [1, 0]   # only `isA` is true

    out = _render_call(call_template, arguments, bits, condition_id)
    assert out.endswith('\n'), repr(out)
    # The caller appends `end if\n`; the result must put it on its own line.
    combined = out + "end if\n"
    assert ')\nend if' in combined, repr(combined)


def test_render_call_with_no_optional_args_keeps_close_paren():
    call_template = "call myRoutine(self{conditions})"
    arguments = {
        'a': {'name': 'a', 'value': 'aVal', 'conditionActual': 'isA'},
    }
    condition_id = {'isA': 1}
    bits = [0]   # NO optional argument selected

    out = _render_call(call_template, arguments, bits, condition_id)
    # No bare leading comma — the close paren immediately follows `self`.
    assert 'call myRoutine(self)' in out


def test_render_call_with_optional_args_inserts_keyword():
    call_template = "call myRoutine(self{conditions})"
    arguments = {
        'a': {'name': 'a', 'value': 'aVal', 'conditionActual': 'isA'},
        'b': {'name': 'b', 'value': 'bVal', 'conditionActual': 'isB'},
    }
    condition_id = {'isA': 1, 'isB': 2}
    bits = [1, 1]

    out = _render_call(call_template, arguments, bits, condition_id)
    assert 'call myRoutine(self,a=aVal,b=bVal)' in out


def test_render_call_missing_conditions_marker_raises():
    """Forgetting `{conditions}` in the template is a directive author error
    and must be surfaced loudly."""
    call_template = "call myRoutine(self)\n"
    with pytest.raises(RuntimeError, match='missing `{conditions}`'):
        _render_call(call_template, {}, [], {})
