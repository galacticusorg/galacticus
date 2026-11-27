"""Tests for the CI/CD workflow planner, and of the workflow's wiring to it."""

import glob
import json
import os
import re
import subprocess
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


def test_a_pass_of_the_full_suite_is_recorded_under_the_name_which_the_plan_looks_for(workflow):
    record = "ci-full-passed-"
    steps = workflow.jobs["CI-Gate"]["steps"]
    check = next(index for index, step in enumerate(steps) if "ciPlan.py gate" in str(step.get("run", "")))
    uploads = [(index, step) for index, step in enumerate(steps) if "upload-artifact" in str(step.get("uses", ""))]
    assert len(uploads) == 1
    index, upload = uploads[0]
    # The record must follow the check (a step runs only if those before it succeeded), and be made only by the full suite.
    assert index > check
    assert upload["if"] == "needs.Plan.outputs.tier == 'full'"
    assert upload["with"]["name"] == record + "${{ needs.Plan.outputs.content }}"
    plan = workflow.jobs["Plan"]
    assert plan["outputs"]["content"].split() == "${{ steps.plan.outputs.content }}".split()
    lookup = next(step for step in plan["steps"] if step.get("id") == "tested")
    # The key is made by the planner itself, so that what is recorded and what is looked up can not come to differ.
    assert "key=$(python scripts/build/ciPlan.py key)" in lookup["run"]
    assert f"artifacts?name={record}${{key}}" in lookup["run"]
    # A record made by a pull request from a fork, which runs its own copy of the workflow, must not be believed.
    assert ".workflow_run.repository_id == .workflow_run.head_repository_id" in lookup["run"]
    assert "--tested-run-id" in str(plan["steps"][-1]["run"])
    assert '--content       "${{ steps.tested.outputs.key }}"' in str(plan["steps"][-1]["run"])


def downloads_of(job):
    """Return `(artifact, parameters)` for each step of a job which downloads an artifact, directly or by the benchmark action."""
    found = []
    for step in job.get("steps", []):
        uses, parameters = str(step.get("uses", "")), step.get("with") or {}
        if "download-artifact" in uses:
            found.append((parameters.get("name") or parameters.get("pattern"), parameters))
        elif uses.endswith("/actions/benchmark"):
            found.append((parameters.get("artifact") or f"{parameters['type']}-{parameters['suffix']}", parameters))
    return found


def test_the_deploy_tier_can_find_every_artifact_which_it_downloads(workflow):
    """In the `deploy` tier only `DEPLOY_JOBS` run, so an artifact must come either from one of those or from the tested run."""
    deployed = [name for name in workflow.jobs if name in ciPlan.DEPLOY_JOBS]
    built = {(step.get("with") or {}).get("name") for name in deployed for step in workflow.jobs[name].get("steps", [])
             if "upload-artifact" in str(step.get("uses", ""))}
    fromTestedRun = "${{ needs.Plan.outputs.tier == 'deploy' && needs.Plan.outputs.resultsRunId || '' }}"
    counted = {True: 0, False: 0}
    for name in deployed:
        for artifact, parameters in downloads_of(workflow.jobs[name]):
            other = "run-id" in parameters
            counted[other] += 1
            if not other:
                assert artifact in built, f"{name} downloads {artifact} from its own run, where no job of the deploy tier makes it"
                assert "github-token" not in parameters
                continue
            assert artifact not in built, f"{name} downloads {artifact} from the tested run, but a job of the deploy tier makes it"
            if "download-artifact" in str(parameters) or "github-token" in parameters:
                # Reading another run's artifacts needs a token, which is passed only when another run is read.
                assert parameters["run-id"] == "${{ needs.Plan.outputs.resultsRunId }}", (name, artifact)
                assert parameters["github-token"] == "${{ needs.Plan.outputs.tier == 'deploy' && github.token || '' }}", (name, artifact)
            else:
                assert parameters["run-id"] == fromTestedRun, (name, artifact)
    assert counted[True] and counted[False]


def test_the_builds_of_the_deploy_tier_need_only_each_other(workflow):
    for name in ciPlan.DEPLOY_BUILD_JOBS:
        assert set(workflow.needs(name)) <= {"Plan"} | set(ciPlan.DEPLOY_BUILD_JOBS), name
    # `Validate` publishes the build profile, and in the `deploy` tier nothing else that it needs waits for that build.
    assert "Build-Executable-Linux" in workflow.needs("Validate")
    assert "needs.Plan.outputs.tier == 'deploy'" in workflow.jobs["Validate"]["if"]
    attribute = next(step for step in workflow.jobs["Validate"]["steps"] if "ciPlan.py attribute" in str(step.get("run", "")))
    assert attribute["if"] == "needs.Plan.outputs.tier == 'deploy'"


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
    """A job which runs one unit test executable downloads that executable alone, so the build must upload it by that name."""
    uploads = {
        step["with"]["name"]: step["with"]["path"]
        for step in workflow.jobs["Build-Executable-Linux-Non-Static"]["steps"]
        if "upload-artifact" in step.get("uses", "")
    }
    expected = set()
    for name in workflow.names():
        job = workflow.jobs[name]
        parameters = job.get("with") or {}
        artifact = parameters.get("artifact", "")
        if artifact == "test-execs":
            # Only the job which runs all of the quick tests takes the whole set, which is several gigabytes.
            assert name == "Codes-Quick", f"{name} downloads every test executable"
        if not artifact.startswith("test-exec-"):
            continue
        executable = parameters["executable"]
        if executable == "${{ matrix.executable }}":
            executables = job["strategy"]["matrix"]["executable"]
        else:
            executables = [executable]
        for executable in executables:
            assert artifact.replace("${{ matrix.executable }}", executable) == "test-exec-" + executable, name
            assert uploads.get("test-exec-" + executable) == executable, f"the build does not upload {executable} on its own"
            expected.add("test-exec-" + executable)
    assert {"Codes", "Parameters", "IO-HDF5", "Crash", "Fail", "Equilibrium-Iteration-Limit"} <= {
        name for name in workflow.names() if (workflow.jobs[name].get("with") or {}).get("artifact", "").startswith("test-exec-")
    }
    assert {name for name in uploads if name.startswith("test-exec-")} == expected, "the build uploads an executable no job asks for"
    assert uploads["test-execs"] == "tests.*.exe"


def test_a_job_which_runs_a_test_file_with_a_named_executable_is_one_entry(workflow):
    assert workflow.entries("Equilibrium-Iteration-Limit") == ["testSuite/test-equilibrium-iteration-limit.py"]


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
        ("Dockerfile", ciPlan.FULL),
        ("pyproject.toml", ("tooling", "model")),
        ("pytest.ini", ("tooling",)),
        ("scripts/doc/extractContributors.py", ("tooling",)),
        ("scripts/aux/parametersMigrate.py", ("tooling", "model")),
        ("scripts/analysis/someAnalysis.py", ciPlan.FULL),
        ("constraints/pipelines/darkMatter/haloMassFunctionBase.xml", ciPlan.FULL),
        ("tutorials/library/executeNotebooks.py", ciPlan.FULL),
        ("AGENTS.md", ()),
        ("docs/agents/ci.md", ()),
        ("aux/words.dict", ()),
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
    # Slow and specialized jobs are left to the full suite, and nothing is published.
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


MERGE_GROUP = "refs/heads/gh-readonly-queue/master/pr-1-abc"


def test_merge_group_of_tested_content_runs_nothing(workflow):
    outputs, reasons = ciPlan.plan(workflow, "merge_group", ref=MERGE_GROUP, tested_run_id="42", content="abc123", run_id="7")
    assert outputs["tier"] == "skip" and jobs_of(outputs) == [] and outputs["deploy"] == "false"
    assert "abc123" in reasons[0] and "42" in reasons[0]


def test_merge_group_of_untested_content_is_refused(workflow):
    """The merge queue never runs the full suite - it could not finish in the time the queue allows - so it must not pass."""
    outputs, reasons = ciPlan.plan(workflow, "merge_group", ref=MERGE_GROUP, content="abc123", run_id="7")
    assert outputs["tier"] == "blocked" and jobs_of(outputs) == [] and outputs["deploy"] == "false"
    assert ciPlan.LABEL_FULL in reasons[0]


def test_plan_command_fails_for_a_refused_merge_group(monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("GITHUB_OUTPUT", str(tmp_path / "output"))
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(tmp_path / "summary"))
    arguments = ["plan", "--event", "merge_group", "--ref", MERGE_GROUP, "--content", "abc123", "--run-id", "7"]
    assert ciPlan.main(arguments) == 1
    assert "::error" in capsys.readouterr().out
    assert "`blocked` tier" in (tmp_path / "summary").read_text()
    assert ciPlan.main(arguments + ["--tested-run-id", "42"]) == 0
    # The failure of the plan is what fails the gate.
    assert ciPlan.gate(results(Plan="failure", A="skipped"), [])


def test_only_the_full_suite_makes_content_tested(workflow):
    """A merge group or a push is believed to be tested only on the word of the lookup, never by default."""
    for event, ref in (("merge_group", MERGE_GROUP), ("push", "refs/heads/master")):
        assert ciPlan.plan(workflow, event, ref=ref, tested_run_id="", run_id="7")[0]["tier"] in ("blocked", "full")


def test_push_of_tested_content_rebuilds_and_deploys(workflow):
    outputs = ciPlan.plan(workflow, "push", ref="refs/heads/master", tested_run_id="42", content="abc123", run_id="7")[0]
    assert outputs["tier"] == "deploy"
    assert set(jobs_of(outputs)) == set(ciPlan.DEPLOY_JOBS)
    # What is published is built by this run; the results published with it are those of the run which tested the content.
    assert set(ciPlan.DEPLOY_BUILD_JOBS) <= set(jobs_of(outputs))
    assert outputs["deploy"] == "true" and outputs["resultsRunId"] == "42" and outputs["content"] == "abc123"


def test_push_of_an_untested_commit_runs_everything_then_deploys(workflow):
    outputs = ciPlan.plan(workflow, "push", ref="refs/heads/master", run_id="7")[0]
    assert outputs["tier"] == "full"
    assert set(jobs_of(outputs)) == set(workflow.names()) - set(ciPlan.PULL_REQUEST_JOBS)
    assert outputs["deploy"] == "true" and outputs["resultsRunId"] == "7"


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
    assert written["tier"] == "pr" and written["resultsRunId"] == "7"
    assert "Codes-Quick" in json.loads(written["jobs"])
    assert "`pr` tier" in (tmp_path / "summary").read_text()
    # A list of changed files which is missing means that the changes are unknown, not that there are none.
    assert ciPlan.main(["plan", "--event", "pull_request", "--changed-files", str(tmp_path / "missing")]) == 0
    assert "tier=full" in capsys.readouterr().out


# --- Attributing published results ------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("variable", ["ANALYSES_DATA", "PONOSV_DATA", "BARYONICSUPPRESSION_DATA"])
def test_attribute_replaces_only_the_commit(variable):
    original = {"repoUrl": "https://example.org", "commit": {"id": "old", "message": "Merge a into b"}, "results": [{"name": "é", "value": 1.5}]}
    text = f"window.{variable} = " + json.dumps(original, indent=4, ensure_ascii=False)
    commit = {"id": "new", "message": "merged"}
    attributed = ciPlan.attribute(text, commit)
    assert attributed.startswith(f"window.{variable} = ")
    assert json.loads(attributed.split(" = ", 1)[1]) == {**original, "commit": commit}
    # Attributing a file to the commit which it already describes changes nothing, so the format is that of the original.
    assert ciPlan.attribute(text, original["commit"]) == text


@pytest.mark.parametrize("text", ["{}", "window.DATA = []", 'window.DATA = {"results": []}', 'window.DATA = {"commit": "abc"}'])
def test_attribute_refuses_what_it_does_not_recognize(text):
    with pytest.raises(ValueError):
        ciPlan.attribute(text, {"id": "new"})


def git(repository, *arguments):
    environment = dict(os.environ, GIT_AUTHOR_NAME="A. Author", GIT_AUTHOR_EMAIL="author@example.org",
                       GIT_COMMITTER_NAME="C. Committer", GIT_COMMITTER_EMAIL="committer@example.org",
                       GIT_AUTHOR_DATE="2026-01-02T03:04:05+00:00", GIT_COMMITTER_DATE="2026-02-03T04:05:06-08:00")
    return subprocess.run(["git", "-C", str(repository), "-c", "commit.gpgsign=false", "-c", "core.hooksPath=/dev/null", *arguments],
                          capture_output=True, text=True, check=True, env=environment).stdout.strip()


def test_head_commit_describes_a_commit_as_the_validation_scripts_do(tmp_path):
    git(tmp_path, "init", "--quiet")
    (tmp_path / "file").write_text("content")
    git(tmp_path, "add", "file")
    git(tmp_path, "commit", "--quiet", "-m", "A subject\n\nAnd a body.")
    commit = ciPlan.head_commit(tmp_path)
    revision = git(tmp_path, "rev-parse", "HEAD")
    assert commit == {
        # The validation scripts record the author, but the time at which the commit was made.
        "author": {"name": "A. Author", "email": "author@example.org"},
        "id": revision,
        "message": "A subject\n\nAnd a body.\n",
        "timestamp": "2026-02-03 04:05:06-08:00",
        "url": "https://github.com/galacticusorg/galacticus/commit/" + revision,
    }


def test_attribute_command_rewrites_files_in_place(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(ciPlan, "head_commit", lambda: {"id": "new"})
    file = tmp_path / "results_a.json"
    file.write_text('window.ANALYSES_DATA = {"commit": {"id": "old"}, "results": []}', encoding="utf-8")
    assert ciPlan.main(["attribute", str(file)]) == 0
    assert json.loads(file.read_text(encoding="utf-8").split(" = ", 1)[1]) == {"commit": {"id": "new"}, "results": []}
    assert "new" in capsys.readouterr().out


# --- The key which identifies tested content --------------------------------------------------------------------------------------


@pytest.fixture(name="checkout")
def fixture_checkout(tmp_path):
    git(tmp_path, "init", "--quiet")
    for name, text in {"source/a.F90": "program a", "docs/index.rst": "Title", "Makefile": "all:", "README.md": "Read me"}.items():
        (tmp_path / name).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / name).write_text(text)
    git(tmp_path, "add", ".")
    git(tmp_path, "commit", "--quiet", "-m", "Initial")
    return tmp_path


def change(checkout, name, text, executable=False):
    file = checkout / name
    file.parent.mkdir(parents=True, exist_ok=True)
    file.write_text(text)
    if executable:
        file.chmod(0o755)
    git(checkout, "add", ".")
    git(checkout, "commit", "--quiet", "-m", "Change " + name)
    return ciPlan.content_key(checkout)


def test_content_key_ignores_files_which_can_not_affect_the_build_or_the_tests(checkout):
    key = ciPlan.content_key(checkout)
    assert re.fullmatch(r"[0-9a-f]{40}", key)
    assert change(checkout, "docs/index.rst", "A new title") == key
    assert change(checkout, "docs/manuals/new.rst", "A new page") == key
    assert change(checkout, "AGENTS.md", "Guidance") == key
    # The key is of what is committed, whatever lies in the working tree.
    (checkout / "source/a.F90").write_text("program uncommitted")
    assert ciPlan.content_key(checkout) == key


@pytest.mark.parametrize("name, text", [("source/a.F90", "program b"), ("source/new.F90", "program a"), ("Makefile", "all: a"),
                                        ("testSuite/test-new.py", "pass"), ("constraints/pipelines/a.xml", "<a/>")])
def test_content_key_changes_with_any_file_which_can(checkout, name, text):
    key = ciPlan.content_key(checkout)
    assert change(checkout, name, text) != key


def test_content_key_changes_with_a_rename_a_removal_or_a_mode(checkout):
    keys = [ciPlan.content_key(checkout)]
    keys.append(change(checkout, "source/a.F90", "program a", executable=True))
    git(checkout, "mv", "source/a.F90", "source/b.F90")
    git(checkout, "commit", "--quiet", "-m", "Rename")
    keys.append(ciPlan.content_key(checkout))
    git(checkout, "rm", "--quiet", "source/b.F90")
    git(checkout, "commit", "--quiet", "-m", "Remove")
    keys.append(ciPlan.content_key(checkout))
    assert len(set(keys)) == len(keys)
    # The same content has the same key, however it was arrived at.
    change(checkout, "source/a.F90", "program a")
    (checkout / "source/a.F90").chmod(0o644)
    assert change(checkout, "docs/index.rst", "Restored") == keys[0]


def test_every_file_which_a_job_of_the_suite_names_is_part_of_the_content(workflow):
    """A test file which matched no path rule would be left out of the key, and so a change to it would never be tested."""
    for name in workflow.names():
        for entry in workflow.entries(name):
            path = ciPlan.entry_path(entry)
            if path:
                assert ciPlan.classify(path), f"{path}, run by {name}, matches no path rule"


def test_key_command(capsys):
    assert ciPlan.main(["key"]) == 0
    assert capsys.readouterr().out.strip() == ciPlan.content_key()
