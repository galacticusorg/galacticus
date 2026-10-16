"""Tests for scripts/aux/workaroundChecker.py.

The checker looks up the GCC Bugzilla status of every `<workaround>` directive
in the source and fails the `Workarounds` workflow when a bug has been fixed but
its workaround is still present. A fix landing upstream does not always mean the
workaround can be removed at once (we must usually wait for a released compiler
carrying the fix to be available on every platform we build for), so a directive
may carry an `issue` attribute pointing at the issue which tracks its removal.
The cases below pin the behavior that makes that useful:

  * a resolved PR whose workaround is staged behind an *open* issue is reported
    but does not fail the check, while an unstaged one still does, and
  * closing the issue while the workaround remains re-arms the failure --- the
    staging must not silently hide the workaround forever.

Both Bugzilla and the GitHub API are reached through `curl` in a subprocess, so
they are stubbed out here rather than being contacted for real.
"""

import importlib.util
import os
import sys

import pytest

_CHECKER = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, 'workaroundChecker.py')


def _load():
    """Import `workaroundChecker.py` as a module (it has no `.py`-importable name
    on the path, and executes nothing at import time beyond its definitions)."""
    spec = importlib.util.spec_from_file_location('workaroundChecker', _CHECKER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _directive(pr, issue=None):
    attributes = f'type="gfortran" PR="{pr}"'
    if issue is not None:
        attributes += f' issue="{issue}"'
    return ('  !![\n'
            f'  <workaround {attributes} docformat="rst">\n'
            '   <description>\n'
            '   A bug.\n'
            '   </description>\n'
            '  </workaround>\n'
            '  !!]\n')


def _run(monkeypatch, tmp_path, source, statuses, issueStates):
    """Run the checker over `source`, with Bugzilla and GitHub stubbed.

    `statuses` maps PR number to the Bugzilla status to report, `issueStates`
    maps issue URL to the state ("open"/"closed") to report.
    """
    module = _load()
    (tmp_path / 'test.F90').write_text(source)
    monkeypatch.setattr(module, 'DIRECTORIES', (str(tmp_path),))
    monkeypatch.setattr(module.time, 'sleep', lambda _seconds: None)
    monkeypatch.setattr(module, 'issue_state',
                        lambda url: issueStates.get(url, 'UNKNOWN'))

    def _curl(command, **_kwargs):
        pr = command[-1].rsplit('=', 1)[1]
        status = statuses.get(pr, 'NEW')
        return type('Result', (), {
            'returncode': 0,
            'stdout': f'<span id="static_bug_status">{status}\n',
        })()

    monkeypatch.setattr(module.subprocess, 'run', _curl)
    return module.main(), module


def test_resolved_unstaged_workaround_fails(monkeypatch, tmp_path, capsys):
    status, _module = _run(monkeypatch, tmp_path, _directive(88632),
                           {'88632': 'RESOLVED'}, {})
    assert status == 1
    assert '!!! Resolved PRs with workarounds exist !!!' in capsys.readouterr().out


def test_resolved_staged_workaround_does_not_fail(monkeypatch, tmp_path, capsys):
    url = 'https://github.com/galacticusorg/galacticus/issues/1234'
    status, _module = _run(monkeypatch, tmp_path, _directive(88632, 1234),
                           {'88632': 'RESOLVED'}, {url: 'open'})
    output = capsys.readouterr().out
    assert status == 0
    assert 'staged for removal' in output
    assert url in output
    assert '!!!' not in output


def test_staged_workaround_fails_once_its_issue_is_closed(monkeypatch, tmp_path,
                                                          capsys):
    url = 'https://github.com/galacticusorg/galacticus/issues/1234'
    status, _module = _run(monkeypatch, tmp_path, _directive(88632, 1234),
                           {'88632': 'RESOLVED'}, {url: 'closed'})
    output = capsys.readouterr().out
    assert status == 1
    assert 'issue is CLOSED but workaround remains' in output


def test_unresolved_workaround_never_fails(monkeypatch, tmp_path):
    status, _module = _run(monkeypatch, tmp_path, _directive(88632),
                           {'88632': 'ASSIGNED'}, {})
    assert status == 0


@pytest.mark.parametrize('issue,expected', [
    ('1234',  'https://github.com/galacticusorg/galacticus/issues/1234'),
    ('#1234', 'https://github.com/galacticusorg/galacticus/issues/1234'),
    ('https://github.com/other/repo/issues/7',
     'https://github.com/other/repo/issues/7'),
    # Attribute values in directives are XML-escaped, as they are read back out
    # of the XML directive blocks - the checker must unescape them.
    ('https:&#x2F;&#x2F;github.com&#x2F;other&#x2F;repo&#x2F;issues&#x2F;7',
     'https://github.com/other/repo/issues/7'),
])
def test_issue_attribute_forms(issue, expected):
    module = _load()
    assert module.issue_url(module.attribute(f'issue="{issue}"', 'issue')) == expected


def test_directive_without_pr_is_ignored(monkeypatch, tmp_path):
    source = ('  !![\n'
              '  <workaround type="gfortran" url="https://example.com/" docformat="rst">\n'
              '   <description>\n'
              '   A bug with no PR.\n'
              '   </description>\n'
              '  </workaround>\n'
              '  !!]\n')
    status, module = _run(monkeypatch, tmp_path, source, {}, {})
    assert status == 0
    assert module.WORKAROUNDS == {}


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
