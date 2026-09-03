"""Unit tests for the Windows/WSL delegation layer.

Network-free and Windows-free: ``wsl.exe`` is never run. Every test replaces
:func:`galacticus_launcher.wsl._run` (the single choke point through which
``wsl.exe`` is invoked) with a fake keyed on the arguments, and forces the host
detection as needed.
"""

import pytest

pytest.importorskip("platformdirs")
pytest.importorskip("requests")

from galacticus_launcher import cli, wsl


LIST_OUTPUT = """  NAME              STATE           VERSION
* docker-desktop    Running         2
  Ubuntu            Stopped         2
  Debian            Stopped         1
"""


class FakeWSL:
    """Stand-in for ``wsl._run``: records calls and answers from a script."""

    def __init__(self, responses):
        # responses: list of (predicate(args) -> bool, Result) checked in order.
        self.responses = responses
        self.calls = []

    def __call__(self, args, *, capture=True):
        self.calls.append((list(args), capture))
        for predicate, result in self.responses:
            if predicate(args):
                # A callable response computes its Result from the arguments.
                return result(args) if callable(result) else result
        raise AssertionError(f"unexpected wsl.exe invocation: {args}")


def _fake(monkeypatch, responses, *, windows=True, have_wsl=True):
    fake = FakeWSL(responses)
    monkeypatch.setattr(wsl, "_run", fake)
    monkeypatch.setattr(wsl, "is_windows", lambda: windows)
    monkeypatch.setattr(wsl.shutil, "which", lambda name: "C:\\Windows\\wsl.exe" if have_wsl else None)
    return fake


# --- decoding wsl.exe output -----------------------------------------------

def test_decode_utf16_and_utf8():
    utf16 = "\ufeffDefault Version: 2\r\n".encode("utf-16-le")
    assert wsl.decode(utf16) == "Default Version: 2\n"
    assert wsl.decode("/mnt/c/Users/x\n".encode("utf-8")) == "/mnt/c/Users/x\n"
    assert wsl.decode(b"") == ""


# --- distribution list & status --------------------------------------------

def test_parse_distribution_list():
    found = wsl.parse_distribution_list(LIST_OUTPUT)
    assert [(d.name, d.version, d.default) for d in found] == [
        ("docker-desktop", 2, True), ("Ubuntu", 2, False), ("Debian", 1, False)]


def test_status_skips_container_runtime_distribution(monkeypatch):
    _fake(monkeypatch, [(lambda a: a[:2] == ["--list", "--verbose"],
                         wsl.Result(0, LIST_OUTPUT))])
    current = wsl.status()
    assert (current.state, current.distribution, current.version) == (wsl.STATE_READY, "Ubuntu", 2)


def test_status_prefers_default_distribution(monkeypatch):
    listing = "  NAME    STATE    VERSION\n  Ubuntu  Stopped  2\n* Debian  Running  2\n"
    _fake(monkeypatch, [(lambda a: a[0] == "--list", wsl.Result(0, listing))])
    assert wsl.status().distribution == "Debian"


def test_status_no_distribution(monkeypatch):
    _fake(monkeypatch, [(lambda a: a[0] == "--list", wsl.Result(1, "no distributions")),
                        (lambda a: a[0] == "--status", wsl.Result(0, "Default Version: 2"))])
    assert wsl.status().state == wsl.STATE_NO_DISTRIBUTION


def test_status_not_installed(monkeypatch):
    _fake(monkeypatch, [(lambda a: a[0] == "--list", wsl.Result(1, "")),
                        (lambda a: a[0] == "--status", wsl.Result(1, "not installed"))])
    assert wsl.status().state == wsl.STATE_NOT_INSTALLED


def test_status_unsupported(monkeypatch):
    _fake(monkeypatch, [], windows=False)
    assert wsl.status().state == wsl.STATE_UNSUPPORTED
    fake = _fake(monkeypatch, [], have_wsl=False)
    assert wsl.status().state == wsl.STATE_UNSUPPORTED
    assert fake.calls == []   # nothing is run when wsl.exe is absent


def test_hint():
    assert "virtualization" in wsl.hint("Error code: Wsl/Service/0x80370102")
    assert wsl.hint("all fine") is None


# --- command construction --------------------------------------------------

def test_launcher_command_forwards_arguments_verbatim():
    command = wsl.launcher_command(["run", "a b.xml", "--dry-run"])
    assert command[:2] == ["/bin/sh", "-c"]
    assert wsl.VENV in command[2] and '"$@"' in command[2]
    assert command[3:] == ["galacticus", "run", "a b.xml", "--dry-run"]


def test_run_builds_wsl_invocation(monkeypatch):
    fake = _fake(monkeypatch, [(lambda a: True, wsl.Result(0, ""))])
    wsl.run("Ubuntu", ["/bin/true"], user="root", capture=False)
    assert fake.calls == [(["-d", "Ubuntu", "-u", "root", "--exec", "/bin/true"], False)]


def test_environment_forwards_release_overrides(monkeypatch):
    monkeypatch.setenv("GALACTICUS_RELEASE_TAG", "v1.2.3")
    monkeypatch.delenv("GALACTICUS_DATASETS_REF", raising=False)
    monkeypatch.setenv("WSLENV", "FOO/p")
    assert wsl._environment()["WSLENV"] == "FOO/p:GALACTICUS_RELEASE_TAG"
    monkeypatch.delenv("GALACTICUS_RELEASE_TAG")
    monkeypatch.delenv("WSLENV")
    assert "WSLENV" not in wsl._environment()


# --- path translation ------------------------------------------------------

def _wslpath_fake(args):
    """Answer the batched `wslpath` call: prefix every path with /wsl."""
    paths = args[args.index("sh") + 1:]
    return wsl.Result(0, "".join(f"/wsl{p}\n" for p in paths))


def test_translate_arguments(monkeypatch, tmp_path):
    fake = _fake(monkeypatch, [(lambda a: "--exec" in a, _wslpath_fake)])
    parameters = tmp_path / "model.xml"
    parameters.write_text("<parameters/>")
    output = tmp_path / "out.xml"            # does not exist yet
    argv = ["run", str(parameters), "missing.xml", "--dry-run", "-o", str(output),
            "--", "--tools"]
    translated = wsl.translate_arguments("Ubuntu", argv)
    assert translated == ["run", f"/wsl{parameters}", "missing.xml", "--dry-run",
                          "-o", f"/wsl{tmp_path}/out.xml", "--", "--tools"]
    # One wslpath invocation for both paths; the output file goes via its parent.
    assert len(fake.calls) == 1
    assert fake.calls[0][0][-2:] == [str(parameters), str(tmp_path)]


def test_translate_arguments_without_paths_runs_nothing(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)   # so the bundled example path does not exist here
    fake = _fake(monkeypatch, [])
    assert wsl.translate_arguments("Ubuntu", ["info"]) == ["info"]
    assert wsl.translate_arguments("Ubuntu", ["run", "parameters/quickTest.xml"]) == \
        ["run", "parameters/quickTest.xml"]
    assert fake.calls == []


def test_translate_paths_reports_failure(monkeypatch):
    _fake(monkeypatch, [(lambda a: True, wsl.Result(1, "wslpath: boom"))])
    with pytest.raises(RuntimeError, match="boom"):
        wsl.translate_paths("Ubuntu", ["C:\\x"])


# --- provisioning ----------------------------------------------------------

def test_is_release():
    assert wsl.is_release("1.2.3")
    assert not wsl.is_release("0.0.0")
    assert not wsl.is_release("0.0.0.post1.dev395")
    assert not wsl.is_release("1.2")


def test_pip_target(monkeypatch, tmp_path):
    assert wsl.pip_target("Ubuntu", version="1.2.3") == ["galacticus==1.2.3"]
    monkeypatch.setattr(wsl, "source_tree", lambda: None)
    assert wsl.pip_target("Ubuntu", version="0.0.0") == ["galacticus"]
    monkeypatch.setattr(wsl, "translate_paths", lambda d, p: [f"/wsl{x}" for x in p])
    assert wsl.pip_target("Ubuntu", version="0.0.0", tree=tmp_path) == ["-e", f"/wsl{tmp_path}"]


def test_source_tree_detects_checkout():
    # These tests run from the repository, so the launcher is a source tree.
    tree = wsl.source_tree()
    assert tree is not None and (tree / "pyproject.toml").is_file()


def test_installed_checks_version(monkeypatch):
    monkeypatch.setattr(wsl, "__version__", "1.2.3")
    _fake(monkeypatch, [(lambda a: True, wsl.Result(0, "galacticus 1.2.3\n"))])
    assert wsl.installed("Ubuntu")
    _fake(monkeypatch, [(lambda a: True, wsl.Result(0, "galacticus 1.0.0\n"))])
    assert not wsl.installed("Ubuntu")
    _fake(monkeypatch, [(lambda a: True, wsl.Result(127, ""))])
    assert not wsl.installed("Ubuntu")


def test_provision_falls_back_to_apt(monkeypatch):
    monkeypatch.setattr(wsl, "pip_target", lambda d: ["galacticus==1.2.3"])
    attempts = []

    def fake_run(distribution, command, *, user=None, capture=True):
        attempts.append((user, command[-1]))
        # First venv attempt fails (no python3-venv); apt and the retry succeed.
        if user is None and len(attempts) == 1:
            return wsl.Result(1, "")
        return wsl.Result(0, "")

    monkeypatch.setattr(wsl, "run", fake_run)
    assert wsl.provision("Ubuntu", log=lambda *a, **k: None) == 0
    assert [user for user, _ in attempts] == [None, "root", None]
    assert "python3-venv" in attempts[1][1]


# --- dispatch --------------------------------------------------------------

def test_dispatch_not_ready(monkeypatch, capsys):
    monkeypatch.setattr(wsl, "status",
                        lambda: wsl.Status(wsl.STATE_NOT_INSTALLED, None, None, ""))
    assert wsl.dispatch(["run", "x.xml"]) == 2
    assert "install-wsl" in capsys.readouterr().err


def test_dispatch_rejects_wsl1(monkeypatch, capsys):
    monkeypatch.setattr(wsl, "status",
                        lambda: wsl.Status(wsl.STATE_READY, "Ubuntu", 1, ""))
    assert wsl.dispatch(["run", "x.xml"]) == 2
    assert "WSL 2" in capsys.readouterr().err


def test_dispatch_runs_inside_distribution(monkeypatch):
    monkeypatch.setattr(wsl, "status",
                        lambda: wsl.Status(wsl.STATE_READY, "Ubuntu", 2, ""))
    monkeypatch.setattr(wsl, "installed", lambda d: True)
    monkeypatch.setattr(wsl, "translate_arguments", lambda d, argv: ["run", "/mnt/c/m.xml"])
    calls = []

    def fake_run(distribution, command, *, user=None, capture=True):
        calls.append((distribution, command, capture))
        return wsl.Result(7, "")

    monkeypatch.setattr(wsl, "run", fake_run)
    assert wsl.dispatch(["run", "C:\\m.xml"]) == 7
    (distribution, command, capture), = calls
    assert distribution == "Ubuntu" and capture is False
    assert command == wsl.launcher_command(["run", "/mnt/c/m.xml"])


def test_dispatch_provisions_when_missing(monkeypatch):
    monkeypatch.setattr(wsl, "status",
                        lambda: wsl.Status(wsl.STATE_READY, "Ubuntu", 2, ""))
    monkeypatch.setattr(wsl, "installed", lambda d: False)
    monkeypatch.setattr(wsl, "provision", lambda d, log: 3)
    assert wsl.dispatch(["info"]) == 3


# --- CLI routing -----------------------------------------------------------

def test_cli_forwards_everything_on_windows(monkeypatch):
    seen = []
    monkeypatch.setattr(wsl, "is_windows", lambda: True)
    monkeypatch.setattr(wsl, "dispatch", lambda argv: seen.append(argv) or 7)
    assert cli.main(["model.xml", "changes.xml", "--dry-run"]) == 7
    assert seen == [["run", "model.xml", "changes.xml", "--dry-run"]]
    assert cli.main(["clean", "--dry-run"]) == 7
    assert seen[-1] == ["clean", "--dry-run"]


def test_cli_install_wsl_routes_to_install(monkeypatch):
    seen = {}
    monkeypatch.setattr(wsl, "install", lambda **kw: seen.update(kw) or 0)
    assert cli.main(["install-wsl", "--yes", "--no-tools"]) == 0
    assert seen == {"yes": True, "tools": False}


def test_cli_install_wsl_is_noop_off_windows(capsys):
    assert cli.main(["install-wsl"]) == 0
    assert "only needed on Windows" in capsys.readouterr().out


# --- `install-wsl` flow ----------------------------------------------------

def _status(state, distribution=None, version=None):
    return wsl.Status(state, distribution, version, "")


def test_install_flow_installs_feature_then_asks_for_reboot(monkeypatch, capsys):
    monkeypatch.setattr(wsl, "is_windows", lambda: True)
    monkeypatch.setattr(wsl, "status", lambda: _status(wsl.STATE_NOT_INSTALLED))
    monkeypatch.setattr(wsl, "enable_wsl", lambda: 0)
    assert wsl.install(yes=True) == 0
    assert "Reboot" in capsys.readouterr().out


def test_install_flow_refuses_without_terminal(monkeypatch, capsys):
    monkeypatch.setattr(wsl, "is_windows", lambda: True)
    monkeypatch.setattr(wsl, "status", lambda: _status(wsl.STATE_NOT_INSTALLED))
    monkeypatch.setattr(wsl.sys.stdin, "isatty", lambda: False)
    monkeypatch.setattr(wsl, "enable_wsl", lambda: pytest.fail("must not install"))
    assert wsl.install(yes=False) == 1
    assert "--yes" in capsys.readouterr().err


def test_install_flow_installs_distribution_then_provisions(monkeypatch):
    monkeypatch.setattr(wsl, "is_windows", lambda: True)
    states = iter([_status(wsl.STATE_NO_DISTRIBUTION),
                   _status(wsl.STATE_READY, "Ubuntu", 2)])
    monkeypatch.setattr(wsl, "status", lambda: next(states))
    monkeypatch.setattr(wsl, "install_distribution", lambda name=None: 0)
    steps = []
    monkeypatch.setattr(wsl, "provision", lambda d, log: steps.append("provision") or 0)

    def fake_run(distribution, command, *, user=None, capture=True):
        steps.append(command[-2:])
        return wsl.Result(0, "")

    monkeypatch.setattr(wsl, "run", fake_run)
    assert wsl.install(yes=True, tools=False, log=lambda *a, **k: None) == 0
    assert steps == ["provision", ["install", "--no-tools"]]
