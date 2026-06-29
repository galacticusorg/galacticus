"""Tests for the launcher's `resolve` subcommand and validate-on-resolved wiring.

These exercise how `galacticus_launcher` drives `Galacticus.Parameters.resolve`
and `.validate` (Stages 4-5 of the parameter-resolver work). They import the
in-repo `Galacticus.Parameters` tree via a stub install whose `exec_path` is the
repository root.
"""

import json
from pathlib import Path

import pytest

from galacticus_launcher import cli, download, paths
from galacticus_launcher import validate as launcher_validate

REPO = Path(__file__).resolve().parents[3]

# Minimal catalog: `accretionHalo` accepts only `simple`/`coldMode`.
MIN_CATALOG = {
    "functionClasses": {
        "accretionHalo": {"default": "simple",
                          "implementations": ["coldMode", "simple"]},
    },
    "implementations": {
        "accretionHaloSimple": {
            "functionClass": "accretionHalo", "label": "simple",
            "parent": "accretionHaloClass",
            "parameters": [], "objects": [], "directNames": [],
        },
    },
    "enumerations": {},
}


def _install():
    return paths.Install(
        source=paths.SOURCE_ENVIRONMENT, tag=None, exec_path=REPO,
        data_path=None, tools_path=None, dynamic_path=None,
        binary=None, assets=None)


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


# --- run argument splitting -------------------------------------------------

def test_split_change_files():
    assert cli._split_change_files(["c1.xml", "c2.xml", "--", "-x"]) \
        == (["c1.xml", "c2.xml"], ["-x"])
    assert cli._split_change_files(["c1.xml", "--opt"]) == (["c1.xml"], ["--opt"])
    assert cli._split_change_files(["--opt", "v"]) == ([], ["--opt", "v"])
    assert cli._split_change_files([]) == ([], [])


# --- resolve subcommand -----------------------------------------------------

def test_resolve_subcommand_expands_xinclude(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "resolve", lambda *a, **k: _install())
    _write(tmp_path, "frag.xml",
           '<parameters><accretionHalo value="simple"/></parameters>')
    main = _write(tmp_path, "main.xml",
                  '<parameters xmlns:xi="http://www.w3.org/2001/XInclude">'
                  '<xi:include href="frag.xml" xpointer="xpointer(parameters/*)"/>'
                  '</parameters>')
    out = tmp_path / "resolved.xml"
    assert cli.main(["resolve", str(main), "-o", str(out)]) == 0
    assert "xi:include" not in out.read_text()
    assert "accretionHalo" in out.read_text()


def test_resolve_subcommand_conditionals(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "resolve", lambda *a, **k: _install())
    text = ('<parameters>'
            '<cosmologicalMassVariance value="filteredPower"><sigma_8 value="0.9"/>'
            '</cosmologicalMassVariance>'
            '<active1 value="x" active="[cosmologicalMassVariance/sigma_8] != 0.9"/>'
            '<active1 value="y" active="[cosmologicalMassVariance/sigma_8] == 0.9"/>'
            '</parameters>')
    main = _write(tmp_path, "main.xml", text)
    out = tmp_path / "resolved.xml"
    assert cli.main(["resolve", str(main), "-o", str(out)]) == 0
    pruned = out.read_text()
    assert 'value="y"' in pruned and 'value="x"' not in pruned   # inactive pruned
    assert "active=" not in pruned                               # marker stripped
    # --no-conditionals keeps them.
    out2 = tmp_path / "resolved2.xml"
    assert cli.main(["resolve", str(main), "-o", str(out2), "--no-conditionals"]) == 0
    assert out2.read_text().count("active1") == 2


def test_resolve_subcommand_error_returns_1(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(paths, "resolve", lambda *a, **k: _install())
    main = _write(tmp_path, "main.xml",
                  '<parameters><accretionHalo idRef="missing"/></parameters>')
    out = tmp_path / "resolved.xml"
    assert cli.main(["resolve", str(main), "-o", str(out)]) == 1
    assert "cannot resolve" in capsys.readouterr().err


# --- validate on the resolved tree (Stage 4) --------------------------------

def test_validate_checks_xincluded_content(tmp_path, monkeypatch):
    # A bad selector inside an XIncluded fragment must be caught -- proving
    # validation runs on the RESOLVED tree, not the raw file.
    catalog = _write(tmp_path, "catalog.json", json.dumps(MIN_CATALOG))
    monkeypatch.setenv("GALACTICUS_PARAMETER_CATALOG", str(catalog))
    _write(tmp_path, "frag.xml",
           '<parameters><accretionHalo value="bogus"/></parameters>')
    main = _write(tmp_path, "main.xml",
                  '<parameters xmlns:xi="http://www.w3.org/2001/XInclude">'
                  '<xi:include href="frag.xml" xpointer="xpointer(parameters/*)"/>'
                  '</parameters>')
    result = launcher_validate.validate(str(main), _install())
    assert result.method == "catalog"
    assert not result.ok
    assert any(f.kind == "selector" for f in result.findings)


def test_validate_surfaces_resolve_error(tmp_path, monkeypatch):
    catalog = _write(tmp_path, "catalog.json", json.dumps(MIN_CATALOG))
    monkeypatch.setenv("GALACTICUS_PARAMETER_CATALOG", str(catalog))
    main = _write(tmp_path, "main.xml",
                  '<parameters><accretionHalo value="simple"/>'
                  '<accretionHalo idRef="missing"/></parameters>')
    result = launcher_validate.validate(str(main), _install())
    assert not result.ok
    assert any(f.kind == "resolve" for f in result.findings)


def test_validate_ok_on_resolved_good_file(tmp_path, monkeypatch):
    catalog = _write(tmp_path, "catalog.json", json.dumps(MIN_CATALOG))
    monkeypatch.setenv("GALACTICUS_PARAMETER_CATALOG", str(catalog))
    main = _write(tmp_path, "main.xml",
                  '<parameters><accretionHalo value="simple"/></parameters>')
    result = launcher_validate.validate(str(main), _install())
    assert result.ok and result.method == "catalog"


# --- catalog generation at provision time (Stage 4 delivery: option B) -------

def _install_at(exec_path):
    return paths.Install(
        source=paths.SOURCE_MANAGED, tag="v1.0.0", exec_path=exec_path,
        data_path=None, tools_path=None, dynamic_path=None, binary=None, assets=None)


def _fake_generator(exec_path):
    gen = exec_path / "scripts" / "build" / "parameterCatalog.py"
    gen.parent.mkdir(parents=True, exist_ok=True)
    gen.write_text("# stub generator\n")
    return gen


def test_provision_catalog_generates_when_missing(tmp_path, monkeypatch):
    _fake_generator(tmp_path)
    captured = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        Path(cmd[3]).write_text("{}")        # the generator writes the catalog
    monkeypatch.setattr(download.subprocess, "run", fake_run)

    assert download._provision_catalog(_install_at(tmp_path), force=False,
                                       log=lambda *a: None) is True
    assert (tmp_path / "parameters.catalog.json").is_file()
    assert captured["cmd"][2:] == [str(tmp_path), str(tmp_path / "parameters.catalog.json")]


def test_provision_catalog_skips_when_present(tmp_path, monkeypatch):
    _fake_generator(tmp_path)
    (tmp_path / "parameters.catalog.json").write_text("{}")
    ran = []
    monkeypatch.setattr(download.subprocess, "run", lambda *a, **k: ran.append(a))
    assert download._provision_catalog(_install_at(tmp_path), force=False,
                                       log=lambda *a: None) is False
    assert ran == []                         # not regenerated
    # ... but force=True regenerates.
    monkeypatch.setattr(download.subprocess, "run",
                        lambda cmd, **k: Path(cmd[3]).write_text("{}"))
    assert download._provision_catalog(_install_at(tmp_path), force=True,
                                       log=lambda *a: None) is True


def test_provision_catalog_without_source_is_noop(tmp_path):
    assert download._provision_catalog(_install_at(tmp_path), force=False,
                                       log=lambda *a: None) is False


def test_find_catalog_locations(tmp_path, monkeypatch):
    monkeypatch.delenv("GALACTICUS_PARAMETER_CATALOG", raising=False)
    install = _install_at(tmp_path)
    assert launcher_validate.find_catalog(install) is None
    # Build-dir location (where `make parameters-catalog` writes).
    build = tmp_path / "work" / "build"
    build.mkdir(parents=True)
    (build / "parameters.catalog.json").write_text("{}")
    assert launcher_validate.find_catalog(install) == build / "parameters.catalog.json"
    # Exec-path root (managed provision) takes precedence over the build dir.
    root_catalog = tmp_path / "parameters.catalog.json"
    root_catalog.write_text("{}")
    assert launcher_validate.find_catalog(install) == root_catalog
    # Explicit env override beats both.
    override = tmp_path / "custom.json"
    override.write_text("{}")
    monkeypatch.setenv("GALACTICUS_PARAMETER_CATALOG", str(override))
    assert launcher_validate.find_catalog(install) == override


def test_provision_catalog_failure_is_best_effort(tmp_path, monkeypatch):
    _fake_generator(tmp_path)

    def boom(*a, **k):
        raise download.subprocess.CalledProcessError(1, "parameterCatalog.py")
    monkeypatch.setattr(download.subprocess, "run", boom)
    messages = []
    assert download._provision_catalog(_install_at(tmp_path), force=False,
                                       log=messages.append) is False
    assert any("could not generate" in m for m in messages)
