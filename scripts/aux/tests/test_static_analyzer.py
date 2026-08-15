"""Tests for check 5 of scripts/aux/staticAnalyzer.py.

Check 5 requires that an `Error_Report` call whose message is built from a
string literal appends `{introspection:location}`, so that a fatal error says
where it was raised. The cases below pin the two judgements that make the check
usable on this source tree:

  * a message which is a bare variable is *not* flagged, because the location may
    have been appended where that variable was assigned, and
  * Fortran embedded in the `!![ ... !!]` directive blocks is XML, so its
    continuation lines are escaped as `&amp;`; a call spread over several such
    lines must still be joined before the check is applied, or a call which does
    append the location is reported as though it did not.
"""

import os
import subprocess
import sys

_ANALYZER = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, 'staticAnalyzer.py')


def _analyze(tmp_path, source):
    """Run the analyzer over `source`, returning (exitStatus, output)."""
    fileName = tmp_path / 'test.F90'
    fileName.write_text(source)
    result = subprocess.run(
        [sys.executable, _ANALYZER, str(fileName)],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True,
    )
    return result.returncode, result.stdout


def test_literal_message_without_location_is_reported(tmp_path):
    status, output = _analyze(tmp_path, """
  subroutine test()
    implicit none
    call Error_Report('something went wrong')
    return
  end subroutine test
""")
    assert status == 1
    assert 'introspection:location' in output


def test_literal_message_with_location_is_accepted(tmp_path):
    status, output = _analyze(tmp_path, """
  subroutine test()
    implicit none
    call Error_Report('something went wrong'//{introspection:location})
    return
  end subroutine test
""")
    assert status == 0
    assert output == ''


def test_location_on_a_continuation_line_is_accepted(tmp_path):
    """The marker is frequently on a continuation line, so lines must be joined
    before the check is applied."""
    status, output = _analyze(tmp_path, """
  subroutine test()
    implicit none
    call Error_Report(                         &
         &            'something went wrong'// &
         &            {introspection:location}  &
         &           )
    return
  end subroutine test
""")
    assert status == 0
    assert output == ''


def test_variable_message_is_not_reported(tmp_path):
    """The location may have been appended where the variable was assigned, which
    cannot be determined from the call site."""
    status, output = _analyze(tmp_path, """
  subroutine test()
    implicit none
    call Error_Report(message)
    return
  end subroutine test
""")
    assert status == 0
    assert output == ''


def test_xml_escaped_continuation_is_joined(tmp_path):
    """Fortran inside a directive block is XML, so continuations are `&amp;`."""
    status, output = _analyze(tmp_path, """
  !![
  <method name="test">
   <code>
    if (bad)                                      &amp;
       &amp; call Error_Report(                   &amp;
       &amp;                   'something wrong'//&amp;
       &amp;                   {introspection:location} &amp;
       &amp;                  )
   </code>
  </method>
  !!]
""")
    assert status == 0
    assert output == ''


def test_commented_out_call_is_not_reported(tmp_path):
    status, output = _analyze(tmp_path, """
  subroutine test()
    implicit none
    ! call Error_Report('something went wrong')
    return
  end subroutine test
""")
    assert status == 0
    assert output == ''


def test_apostrophe_in_trailing_comment_does_not_defeat_the_exemption(tmp_path):
    """The literal test must apply to the call's arguments, not to the rest of
    the line: an apostrophe in a trailing comment would otherwise make a message
    which is a bare variable look as though it were built from a literal."""
    status, output = _analyze(tmp_path, """
  subroutine test()
    implicit none
    call Error_Report(message) ! we don't know the location here
    return
  end subroutine test
""")
    assert status == 0
    assert output == ''


def test_parenthesis_inside_a_message_does_not_end_the_arguments(tmp_path):
    """A closing parenthesis inside the message must not be mistaken for the end
    of the argument list, which would hide the location that follows it."""
    status, output = _analyze(tmp_path, """
  subroutine test()
    implicit none
    call Error_Report('bad value (see the manual)'//{introspection:location})
    return
  end subroutine test
""")
    assert status == 0
    assert output == ''


def test_location_only_in_a_trailing_comment_is_still_reported(tmp_path):
    """Mentioning the marker in a comment must not satisfy the check."""
    status, output = _analyze(tmp_path, """
  subroutine test()
    implicit none
    call Error_Report('something went wrong') ! {introspection:location}
    return
  end subroutine test
""")
    assert status == 1
    assert 'introspection:location' in output
