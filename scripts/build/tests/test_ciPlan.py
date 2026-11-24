"""Tests for the CI/CD workflow planner, and of the workflow's wiring to it."""

import glob
import json
import os
import re
import sys
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import ciPlan  # noqa: E402

REPOSITORY = ciPlan.REPOSITORY

# Jobs which follow the gate rather than being checked by it, with the job each must wait for.
AFTER_GATE = {"Submit-Notarization-MacOS": "CI-Gate", "Deploy": "CI-Gate", "Update-Release-Tag": "Deploy"}
# Jobs which only report, and whose failure must not block a merge.
NOT_GATED = set(AFTER_GATE) | {"Validation-Report"}

# Unit test programs which no job runs. These predate the check below; a new test program must be added to a job instead.
UNLISTED_TEST_PROGRAMS = {
    "tests.output_analysis_capabilities.exe",
    "tests.resource_manager.exe",
    "tests.splashback_radii.exe",
}


@pytest.fixture(name="workflow", scope="module")
def fixture_workflow():
    return ciPlan.Workflow()


def jobs_of(outputs):
    return json.loads(outputs["jobs"])


def skipped_of(outputs):
    return json.loads(outputs["skipEntries"])


def pull_request(workflow, changed, **options):
    return ciPlan.plan(workflow, "pull_request", ref="refs/pull/1/merge", changed=changed, run_id="7", **options)[0]


# --- The workflow's wiring -------------------------------------------------------------------------------------------------------


def test_every_job_is_gated_by_the_plan(workflow):
    for name in workflow.names():
        job = workflow.jobs[name]
        assert "Plan" in workflow.needs(name), f"{name} must need the Plan job"
        assert f"contains(fromJSON(needs.Plan.outputs.jobs), '{name}')" in str(job.get("if", "")), (
            f"{name} must run only if the plan lists it"
        )


def test_the_gate_needs_every_job_which_builds_or_tests(workflow):
    gate = workflow.jobs["CI-Gate"]
    # The gate must run when jobs it needs fail or are skipped, so its condition needs a status function to displace the implicit
    # success(). It must not use always(), which keeps a superseded run alive until the gate is given a runner.
    assert gate["if"] == "${{ !cancelled() || github.event_name == 'merge_group' }}"
    expected = {"Plan"} | (set(workflow.names()) - NOT_GATED)
    assert set(gate["needs"]) == expected


def test_jobs_after_the_gate_wait_for_it_explicitly(workflow):
    for name, predecessor in AFTER_GATE.items():
        condition = workflow.jobs[name]["if"]
        assert predecessor in workflow.needs(name)
        # The implicit success() would skip these whenever an earlier job in the chain was skipped.
        assert "!cancelled()" in condition
        assert f"needs.{predecessor}.result == 'success'" in condition


def test_jobs_refer_only_to_jobs_which_they_need(workflow):
    """A reference to `needs.<job>` for a job which is not needed evaluates to nothing, silently."""
    for name, job in workflow.jobs.items():
        needed = set(workflow.needs(name))
        referenced = set(re.findall(r"needs\.([A-Za-z0-9_-]+)\.(?:outputs|result)", json.dumps(job)))
        assert referenced <= needed, f"{name} refers to {sorted(referenced - needed)} without needing them"
        for needs in needed:
            assert needs in workflow.jobs, f"{name} needs {needs}, which does not exist"


def test_matrices_of_tests_can_skip_entries(workflow):
    for name in workflow.names():
        job = workflow.jobs[name]
        parameters = job.get("with") or {}
        target = parameters.get("file", parameters.get("executable"))
        if isinstance(target, str) and target.startswith("${{ matrix."):
            key = re.fullmatch(r"\$\{\{\s*matrix\.(\w+)\s*\}\}", target).group(1)
            assert parameters.get("enabled") == f"${{{{ !contains(fromJSON(needs.Plan.outputs.skipEntries), matrix.{key}) }}}}", name


def test_the_tables_name_only_what_exists(workflow):
    names = set(workflow.names())
    for group in ciPlan.PR_GROUPS.values():
        assert set(group) <= names
    assert set(ciPlan.DEPLOY_JOBS) <= names
    assert set(ciPlan.PULL_REQUEST_JOBS) <= names
    entries = {entry for name in names for entry in workflow.entries(name)}
    assert set(ciPlan.PR_DEFERRED_ENTRIES) <= entries


def test_every_unit_test_program_is_run_exactly_once(workflow):
    programs = set()
    for file in glob.glob(str(REPOSITORY / "source" / "tests" / "**" / "*.F90"), recursive=True):
        if re.search(r"^\s*program\s+\w+", Path(file).read_text(errors="replace"), re.IGNORECASE | re.MULTILINE):
            programs.add("tests." + os.path.relpath(file, REPOSITORY / "source" / "tests")[:-4].replace(os.sep, ".") + ".exe")
    # The jobs which test the default build must not run a test twice. (The jobs which test other builds may repeat one.)
    default = [executable for name in ("Codes-Quick", "Codes", "Parameters", "IO-HDF5", "Crash", "Fail") for executable in workflow.entries(name)]
    assert len(default) == len(set(default)), "a unit test is run by more than one job"
    others = [executable for name in ("Codes-MPI", "Codes-Debugging", "Codes-MacOS") for executable in workflow.entries(name)]
    assert set(default + others) <= programs, "a job runs a unit test which has no source"
    assert programs - set(default + others) == UNLISTED_TEST_PROGRAMS


def test_unit_tests_run_singly_have_an_artifact_each(workflow):
    """A job which runs one unit test downloads that executable alone, so the build must upload it by that name."""
    uploads = {
        step["with"]["name"]: step["with"]["path"]
        for step in workflow.jobs["Build-Executable-Linux-Non-Static"]["steps"]
        if "upload-artifact" in step.get("uses", "")
    }
    expected = set()
    for name in ("Codes", "Parameters", "IO-HDF5", "Crash", "Fail"):
        artifact = workflow.jobs[name]["with"]["artifact"]
        for executable in workflow.entries(name):
            assert artifact.replace("${{ matrix.executable }}", executable) == "test-exec-" + executable, name
            assert uploads.get("test-exec-" + executable) == executable, f"the build does not upload {executable} on its own"
            expected.add("test-exec-" + executable)
    assert {name for name in uploads if name.startswith("test-exec-")} == expected, "the build uploads an executable no job asks for"
    # The job which runs all of the quick tests takes the whole set.
    assert workflow.jobs["Codes-Quick"]["with"]["artifact"] == "test-execs"
    assert uploads["test-execs"] == "tests.*.exe"


# --- Path rules ------------------------------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "pattern, path, matches",
    [
        ("source/**", "source/a/b.F90", True),
        ("source/**", "sources/a.F90", False),
        ("python/**/tests/**", "python/Galacticus/Parameters/tests/test_query.py", True),
        ("python/**/tests/**", "python/tests/test_x.py", True),
        ("scripts/build/ghPages*.yml", "scripts/build/ghPagesMetricsManifest.yml", True),
        ("scripts/build/ghPages*.yml", "scripts/build/sub/ghPagesX.yml", False),
        ("Makefile", "Makefile", True),
        ("Makefile", "docs/Makefile", False),
    ],
)
def test_glob_to_regex(pattern, path, matches):
    assert bool(ciPlan.glob_to_regex(pattern).match(path)) is matches


@pytest.mark.parametrize(
    "path, selection",
    [
        ("docs/manuals/developer-guide/coding.rst", ()),
        ("README.md", ()),
        (".github/workflows/prChecks.yml", ()),
        ("source/dark_matter_profiles_DMO/NFW.F90", ("tooling", "unit", "model")),
        ("Makefile", ciPlan.FULL),
        ("aux/dependencies.yml", ciPlan.FULL),
        ("scripts/build/preprocess.py", ciPlan.FULL),
        ("scripts/build/ciPlan.py", ("tooling",)),
        ("scripts/build/tests/test_ciPlan.py", ("tooling",)),
        ("python/Galacticus/Build/SourceTree/process.py", ciPlan.FULL),
        ("python/Galacticus/Build/SourceTree/tests/test_enumeration.py", ("tooling",)),
        ("python/Galacticus/Parameters/query.py", ("tooling", "model")),
        (".github/workflows/cicd.yml", ciPlan.FULL),
        (".github/workflows/testModel.yml", ciPlan.FULL),
        (".github/actions/benchmark/action.yml", ciPlan.FULL),
        ("parameters/quickTest.xml", ("tooling",)),
        ("testSuite/test-bar-instability.py", ("tooling", "model")),
    ],
)
def test_classify(path, selection):
    assert ciPlan.classify(path) == selection


# --- Plans -----------------------------------------------------------------------------------------------------------------------


def test_draft_pull_request_runs_nothing(workflow):
    outputs = pull_request(workflow, ["source/a.F90"], draft=True)
    assert outputs["tier"] == "skip" and jobs_of(outputs) == []


def test_documentation_only_pull_request_runs_nothing(workflow):
    outputs = pull_request(workflow, ["docs/index.rst", "README.md"])
    assert outputs["tier"] == "skip" and jobs_of(outputs) == []


def test_source_change_runs_the_subset(workflow):
    outputs = pull_request(workflow, ["source/dark_matter_profiles_DMO/NFW.F90", "docs/index.rst"])
    jobs = set(jobs_of(outputs))
    assert outputs["tier"] == "pr"
    # Both default builds are brought in by the tests which need them, and nothing else is built.
    assert {name for name in jobs if name.startswith("Build-")} == {"Build-Executable-Linux", "Build-Executable-Linux-Non-Static"}
    assert {"Test-Preprocessor", "Codes-Quick", "Model-Integration", "Miscellaneous", "Regressions", "Watch-Updates"} <= jobs
    # Slow and specialized jobs are left to the merge queue, and nothing is published.
    assert not jobs & {"Codes", "Codes-MPI", "Codes-MacOS", "Validate", "Validation-Report", "Merger-Tree-Builder", "Python-Interface"}
    assert not jobs & set(ciPlan.PUBLISH_JOBS)
    assert skipped_of(outputs) == sorted(ciPlan.PR_DEFERRED_ENTRIES)
    assert outputs["deploy"] == "false"


def test_python_test_change_needs_no_build(workflow):
    outputs = pull_request(workflow, ["python/Galacticus/Parameters/tests/test_query.py"])
    assert set(jobs_of(outputs)) == {"Test-Preprocessor", "Python-Utils-Tests", "Watch-Updates"}


@pytest.mark.parametrize("changed", [["Makefile"], ["source/a.F90", "scripts/build/preprocess.py"], [".github/workflows/cicd.yml"], None])
def test_unpredictable_changes_run_everything(workflow, changed):
    outputs = pull_request(workflow, changed)
    jobs = set(jobs_of(outputs))
    assert outputs["tier"] == "full"
    assert jobs == set(workflow.names()) - set(ciPlan.PUBLISH_JOBS)
    assert skipped_of(outputs) == []


def test_label_runs_everything(workflow):
    outputs = pull_request(workflow, ["source/a.F90"], labels=["bug", ciPlan.LABEL_FULL])
    assert outputs["tier"] == "full" and "Validation-Report" in jobs_of(outputs)


def test_merge_group_runs_everything_but_publishes_nothing(workflow):
    outputs = ciPlan.plan(workflow, "merge_group", ref="refs/heads/gh-readonly-queue/master/pr-1-abc", run_id="7")[0]
    jobs = set(jobs_of(outputs))
    assert outputs["tier"] == "full"
    assert jobs == set(workflow.names()) - set(ciPlan.PUBLISH_JOBS) - set(ciPlan.PULL_REQUEST_JOBS)
    assert outputs["deploy"] == "false" and outputs["artifactRunId"] == "7" and skipped_of(outputs) == []


def test_push_of_a_tested_commit_only_deploys(workflow):
    outputs = ciPlan.plan(workflow, "push", ref="refs/heads/master", tested_run_id="42", run_id="7")[0]
    assert outputs["tier"] == "deploy"
    assert set(jobs_of(outputs)) == set(ciPlan.DEPLOY_JOBS)
    assert outputs["deploy"] == "true" and outputs["artifactRunId"] == "42"


def test_push_of_an_untested_commit_runs_everything_then_deploys(workflow):
    outputs = ciPlan.plan(workflow, "push", ref="refs/heads/master", run_id="7")[0]
    assert outputs["tier"] == "full"
    assert set(jobs_of(outputs)) == set(workflow.names()) - set(ciPlan.PULL_REQUEST_JOBS)
    assert outputs["deploy"] == "true" and outputs["artifactRunId"] == "7"


def test_dispatch_deploys_only_from_master(workflow):
    onMaster = ciPlan.plan(workflow, "workflow_dispatch", ref="refs/heads/master", run_id="7")[0]
    onBranch = ciPlan.plan(workflow, "workflow_dispatch", ref="refs/heads/feature", run_id="7")[0]
    subset = ciPlan.plan(workflow, "workflow_dispatch", ref="refs/heads/master", tier="pr", run_id="7")[0]
    assert onMaster["tier"] == onBranch["tier"] == "full"
    assert "Deploy" in jobs_of(onMaster) and "Deploy" not in jobs_of(onBranch)
    assert subset["tier"] == "pr" and subset["deploy"] == "false" and "Deploy" not in jobs_of(subset)


def test_changed_test_runs_alone(workflow):
    outputs = pull_request(workflow, ["testSuite/test-bar-instability.py"])
    jobs = set(jobs_of(outputs))
    assert outputs["tier"] == "pr"
    assert jobs == {"Test-Preprocessor", "Python-Utils-Tests", "Build-Executable-Linux", "Miscellaneous", "Watch-Updates"}
    enabled = set(workflow.entries("Miscellaneous")) - set(skipped_of(outputs))
    assert enabled == {"test-bar-instability.py"}


def test_changed_slow_test_runs_despite_being_deferred(workflow):
    slow = "test-checkpointing.py"
    assert slow in ciPlan.PR_DEFERRED_ENTRIES
    outputs = pull_request(workflow, ["testSuite/" + slow])
    assert slow not in skipped_of(outputs) and "Miscellaneous" in jobs_of(outputs)


def test_changed_test_of_another_build_brings_that_build(workflow):
    outputs = pull_request(workflow, ["testSuite/test-stateRestore_MPI.py"])
    assert {"MPI-State-Restore", "Build-Executables-MPI-Linux"} <= set(jobs_of(outputs))
    assert "Build-Executable-Linux" not in jobs_of(outputs)


def test_unidentifiable_test_suite_change_runs_every_model_test(workflow):
    outputs = pull_request(workflow, ["testSuite/aFileWhichNoTestMentions.dat"])
    assert set(ciPlan.PR_GROUPS["model"]) <= set(jobs_of(outputs))
    assert skipped_of(outputs) == sorted(ciPlan.PR_DEFERRED_ENTRIES)


def test_test_suite_change_with_source_change_is_not_selective(workflow):
    outputs = pull_request(workflow, ["testSuite/test-bar-instability.py", "source/a.F90"])
    assert {"Codes-Quick", "Regressions", "Methods"} <= set(jobs_of(outputs))
    assert skipped_of(outputs) == sorted(ciPlan.PR_DEFERRED_ENTRIES)


def test_closure_brings_needed_jobs(workflow):
    assert workflow.closure(["Python-Interface"]) == ["Build-Library-Linux", "Python-Interface"]


# --- The gate --------------------------------------------------------------------------------------------------------------------


def results(**jobs):
    return {name.replace("_", "-"): {"result": result} for name, result in jobs.items()}


def test_gate_passes_when_planned_jobs_succeed_and_others_are_skipped():
    assert ciPlan.gate(results(Plan="success", A="success", B="skipped"), ["A"]) == []


def test_gate_passes_a_skipped_run():
    assert ciPlan.gate(results(Plan="success", A="skipped"), []) == []


@pytest.mark.parametrize("result", ["failure", "cancelled"])
def test_gate_fails_on_failure_or_cancellation(result):
    assert ciPlan.gate(results(Plan="success", A=result), ["A"])
    # Even for a job which was not planned.
    assert ciPlan.gate(results(Plan="success", A=result), [])


def test_gate_fails_when_a_planned_job_did_not_run():
    problems = ciPlan.gate(results(Plan="success", Build="failure", Test="skipped"), ["Build", "Test"])
    assert len(problems) == 2


def test_gate_fails_when_the_plan_failed():
    assert ciPlan.gate(results(Plan="failure", A="skipped"), [])
    assert ciPlan.gate(results(Plan="skipped", A="skipped"), [])


def test_gate_command(monkeypatch, capsys):
    monkeypatch.delenv("GITHUB_STEP_SUMMARY", raising=False)
    monkeypatch.setenv("CI_NEEDS", json.dumps(results(Plan="success", A="success")))
    monkeypatch.setenv("CI_JOBS", '["A"]')
    assert ciPlan.main(["gate"]) == 0
    assert capsys.readouterr().out.startswith("SUCCESS")
    monkeypatch.setenv("CI_NEEDS", json.dumps(results(Plan="success", A="skipped")))
    assert ciPlan.main(["gate"]) == 1
    assert capsys.readouterr().out.startswith("FAILED")
    # With no plan to check against, a run in which nothing failed still does not pass.
    monkeypatch.setenv("CI_NEEDS", json.dumps(results(Plan="failure", A="skipped")))
    monkeypatch.setenv("CI_JOBS", "")
    assert ciPlan.main(["gate"]) == 1


def test_plan_command_reads_changed_files(monkeypatch, tmp_path, capsys):
    changed = tmp_path / "changed"
    changed.write_text("docs/index.rst\0source/a b.F90\0")
    output = tmp_path / "output"
    monkeypatch.setenv("GITHUB_OUTPUT", str(output))
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(tmp_path / "summary"))
    assert ciPlan.main(["plan", "--event", "pull_request", "--labels", "null", "--changed-files", str(changed), "--run-id", "7"]) == 0
    written = dict(line.split("=", 1) for line in output.read_text().splitlines())
    assert written["tier"] == "pr" and written["artifactRunId"] == "7"
    assert "Codes-Quick" in json.loads(written["jobs"])
    assert "`pr` tier" in (tmp_path / "summary").read_text()
    # A list of changed files which is missing means that the changes are unknown, not that there are none.
    assert ciPlan.main(["plan", "--event", "pull_request", "--changed-files", str(tmp_path / "missing")]) == 0
    assert "tier=full" in capsys.readouterr().out
