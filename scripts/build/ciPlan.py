#!/usr/bin/env python3
"""Plan which jobs of the CI/CD workflow (``.github/workflows/cicd.yml``) run for an event, and gate its result.

The full suite costs over a hundred runner-hours, so it is not run on every push to every pull request. Instead each run of the
workflow is given a *tier* by its ``Plan`` job, and every other job runs only if this script lists it:

``full``
    Every job. Used for merge groups (the merge queue tests each group against the exact state of ``master`` it will land on), for
    pull requests which change the build system or carry the ``ci:full`` label, for manual dispatches, and for any push to
    ``master`` which did not come through the merge queue.
``pr``
    A subset chosen from the files a pull request changes: the two default Linux builds, the quick unit tests, and the model tests
    which run in a few minutes each. Tests too slow or too specialized for that (other build configurations, macOS, MPI,
    validation, benchmarks, profiling) are left to the merge queue.
``deploy``
    A push to ``master`` of a commit which the merge queue has already tested. Nothing is rebuilt or retested: the deployment jobs
    run using the artifacts of the merge group's run.
``skip``
    Nothing: a draft pull request, or one which changes no file that can affect the build or the tests.

The ``CI-Gate`` job then checks the outcome against the plan. It is the one status check that needs to be required by branch
protection, and it reports for every run - including a skipped one - so that a required check is never left pending.

What a pull request runs is decided by the tables below, not by the workflow file. Three properties are deliberate:

* Selection fails safe. A test added to the matrix of a job in the `pr` tier runs for pull requests unless it is listed here
  as deferred; a changed file which matches no rule selects nothing only because every path that can affect the build or the
  tests has a rule; and a list of changed files which could not be determined gets the full suite. (A *job* which this script
  does not know about is not in the `pr` tier: it runs in the full suite, and so in every merge group, until it is added to a
  group below.)
* A selected job brings the jobs it ``needs`` with it, so a test can never be planned without its build.
* A job which was planned but did not run fails the gate, as does one which failed or was canceled.

Run with::

    ciPlan.py plan --event <event> [--ref <ref>] [--draft] [--labels <json>] [--changed-files <file>] [--tier <tier>]
                   [--tested-run-id <id>] [--run-id <id>]
    ciPlan.py gate [--needs <json>] [--jobs <json>]
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[2]
WORKFLOW = REPOSITORY / ".github" / "workflows" / "cicd.yml"

TIERS = ("skip", "pr", "full", "deploy")

# The label which runs the full suite on a pull request.
LABEL_FULL = "ci:full"

# Jobs which are not gated by the plan: the plan itself, and the gate which checks it.
UNGATED_JOBS = ("Plan", "CI-Gate")

# Jobs which publish something, and so run only for a push to (or a dispatch on) `master`.
PUBLISH_JOBS = ("Submit-Notarization-MacOS", "Deploy", "Update-Release-Tag")

# The jobs of the `deploy` tier. `Build-Docker-Linux` builds its image from source and pushes it, and `Validate` publishes the
# validation results to `gh-pages`; neither can be replaced by an artifact of the merge group's run.
DEPLOY_JOBS = ("Build-Docker-Linux", "Validate") + PUBLISH_JOBS

# Jobs which comment on a pull request, and so run only for one. `Validation-Report` needs the validation results, which only
# the full suite produces.
PULL_REQUEST_JOBS = {"Watch-Updates": ("pr", "full"), "Validation-Report": ("full",)}

# The jobs of the `pr` tier, in groups which the path rules below select. The builds are not listed: they are brought in by the
# `needs` of whichever of these jobs are selected.
PR_GROUPS = {
    # Pure-Python checks: the preprocessor regression tests, parameter file validation, and the schema.
    "tooling": ("Test-Preprocessor", "Python-Utils-Tests"),
    # Unit tests of the non-static build. `Codes-Quick` runs every unit test which takes less than a few minutes in one job; the
    # few slow ones are the `Codes` matrix, which is left to the merge queue.
    "unit": ("Codes-Quick", "Parameters", "IO-HDF5", "Crash", "Fail"),
    # Model tests of the static build.
    "model": (
        "Cosmology",
        "Mass-Conservation",
        "Model-Integration",
        "Merger-Tree-Write",
        "Miscellaneous",
        "Miscellaneous-Unexpected",
        "Miscellaneous-IRATE",
        "Miscellaneous-Deep",
        "Miscellaneous-Models",
        "Outputs",
        "Emulator-Design",
        "Emulator-Predict",
        "Regressions",
        "Methods",
    ),
}

# Matrix entries of `pr` tier jobs which are left to the merge queue because they are slow: each of these took over five minutes
# to run (not counting job setup) when the tiers were introduced. Listing an entry here only defers it - it still runs in every
# merge group - so the list errs on the side of brevity, and an entry which is not listed runs on pull requests.
PR_DEFERRED_ENTRIES = (
    # Miscellaneous
    "test-constraint-deterministic-spins.py",
    "test-checkpointing.py",
    "test-impulsiveHeating.py",
    "test-output-tree-contiguousity.py",
    "test-concentration-Ludlow2016.py",
    "test-dust-emission-energy-balance.py",
    "test-merger-tree-export.sh",
    "test-concentration-Johnson2021.py",
    "test-prompt-cusps-unformed.py",
    "test-SIDM-parametric.py",
    "test-parallelTreeBuild.py",
    "test-luminosities-from-SFH-lightcone.py",
    "test-transfer-functions-axion.py",
    "test-power-spectrum-normalizer-amplitude-scalar.py",
    "test-luminosities-from-SED-lightcone.py",
    "test-promptCuspNFW.py",
    "test-pruneLightcone.py",
    "test-satellite-distance-minimum.py",
    "test-postprocess-outputs.py",
    "test-massProfileSatelliteBoundMass.py",
    "test-decayingDM-shell-crossing-capture.py",
    "test-treeFilterLabels.py",
    # Regressions
    "testSuite/regressions/preInfallEvolutionNotOnMainBranch.xml",
    "testSuite/regressions/particulate.py",
    "testSuite/regressions/heatingMonotonicWeakShellCrossingRadiusRoot.xml",
    # Miscellaneous-Models
    "testSuite/parameters/tidalStrippingGalacticComponentsTruncated.xml",
    "testSuite/parameters/blackHoleBinarySeparationGrowthRateChandrasekharIntegral.xml",
    "testSuite/parameters/tidalStrippingGalacticComponentsHeated.xml",
    "testSuite/parameters/tidalStrippingGalacticComponents.xml",
    "testSuite/parameters/mergerTreeEvolverThreaded.xml",
)

# What a changed file selects. The first rule whose pattern matches a file decides for that file; a file matching no rule can
# not affect the build or the tests (documentation, for the most part) and selects nothing. `FULL` selects the full suite, and is
# used where the effect of a change can not be predicted from its path: the build system, the code generator, the dependencies,
# and the workflows themselves. In a pattern `*` matches within one path component and `**` matches across them.
FULL = "full"
PATH_RULES = (
    # Scripts which live with the build system but take no part in the build.
    ("scripts/build/tests/**", ("tooling",)),
    ("scripts/build/ciPlan.py", ("tooling",)),
    ("scripts/build/buildGhPagesIndex.py", ("tooling",)),
    ("scripts/build/buildPrValidationReport.py", ("tooling",)),
    ("scripts/build/prValidationReportTemplate.html", ("tooling",)),
    ("scripts/build/ghPages*.yml", ("tooling",)),
    # The build system, the code generator, and the dependencies.
    ("Makefile", FULL),
    ("aux/dependencies.yml", FULL),
    ("scripts/build/**", FULL),
    ("python/Galacticus/Build/**/tests/**", ("tooling",)),
    ("python/Galacticus/Build/**", FULL),
    # The workflows which make up the suite.
    (".github/workflows/cicd.yml", FULL),
    (".github/workflows/testModel.yml", FULL),
    (".github/workflows/testCode.yml", FULL),
    (".github/workflows/profile.yml", FULL),
    (".github/actions/**", FULL),
    # Source code.
    ("source/**", ("tooling", "unit", "model")),
    # Python: tests of the Python modules need no build; the modules themselves are used by the model tests.
    ("python/**/tests/**", ("tooling",)),
    ("python/**", ("tooling", "model")),
    # Parameter files and their schema are validated by the tooling jobs.
    ("schema/**", ("tooling",)),
    ("parameters/**", ("tooling",)),
    ("parameters.xml", ("tooling",)),
    # The test suite. A pull request which changes nothing else runs only the tests it affects - see `select_affected()`.
    ("testSuite/**", ("tooling", "model")),
)


def glob_to_regex(pattern):
    """Translate a path pattern in which `*` matches within a path component and `**` matches across components."""
    regex = ""
    position = 0
    while position < len(pattern):
        if pattern.startswith("**/", position):
            regex += "(?:.*/)?"
            position += 3
        elif pattern.startswith("**", position):
            regex += ".*"
            position += 2
        elif pattern[position] == "*":
            regex += "[^/]*"
            position += 1
        elif pattern[position] == "?":
            regex += "[^/]"
            position += 1
        else:
            regex += re.escape(pattern[position])
            position += 1
    return re.compile(regex + r"\Z")


_PATH_RULES = tuple((glob_to_regex(pattern), selection) for pattern, selection in PATH_RULES)


def classify(path):
    """Return what a changed file selects: `FULL`, a tuple of `PR_GROUPS` names, or an empty tuple for a file which matters not."""
    for regex, selection in _PATH_RULES:
        if regex.match(path):
            return selection
    return ()


class Workflow:
    """The parts of the workflow which planning needs: each job's `needs`, and the test files of the jobs which run them."""

    def __init__(self, path=WORKFLOW):
        import yaml  # Imported here so that `gate`, which needs no workflow, runs without PyYAML.

        self.jobs = yaml.safe_load(Path(path).read_text())["jobs"]

    def names(self):
        return [name for name in self.jobs if name not in UNGATED_JOBS]

    def needs(self, name):
        needs = self.jobs[name].get("needs", [])
        return [needs] if isinstance(needs, str) else list(needs)

    def closure(self, names):
        """Return the named jobs together with every job that they need, in workflow order."""
        selected = set()
        pending = list(names)
        while pending:
            name = pending.pop()
            if name in selected or name in UNGATED_JOBS:
                continue
            selected.add(name)
            pending.extend(self.needs(name))
        return [name for name in self.jobs if name in selected]

    def entries(self, name):
        """Return the test files (or test executables) which a job runs: one per matrix entry, or the single one it names."""
        job = self.jobs[name]
        parameters = job.get("with") or {}
        target = parameters.get("file", parameters.get("executable"))
        if not isinstance(target, str):
            return []
        matrixed = re.fullmatch(r"\$\{\{\s*matrix\.(\w+)\s*\}\}", target)
        if not matrixed:
            # A job which runs unit tests (and so names no test file) may name several executables, separated by white space.
            return [target] if "file" in parameters else target.split()
        values = ((job.get("strategy") or {}).get("matrix") or {}).get(matrixed.group(1), [])
        return [value for value in values if isinstance(value, str)] if isinstance(values, list) else []


def entry_path(entry, repository=REPOSITORY):
    """Return the path of a test file in the repository. Files are named relative to the root, or to `testSuite/`."""
    for candidate in (entry, "testSuite/" + entry):
        if (Path(repository) / candidate).is_file():
            return candidate
    return None


def select_affected(workflow, changed, repository=REPOSITORY):
    """Select the tests affected by changes to files of the test suite.

    A test is affected by a changed file if it *is* that file, or if its script mentions the file by name or (for a file in a
    subdirectory of `testSuite/`, such as a set of parameter files which a script reads as a whole) mentions that subdirectory.
    Every job which runs a test is considered, not only those of the `pr` tier, so that a pull request which changes (say) an MPI
    test runs that test. Returns a `(jobs, skipped entries)` pair, or `None` if some changed file affects no test that could be
    identified - a shared helper, say - in which case the caller must fall back to running every model test.
    """
    tests = {}
    for name in workflow.names():
        for entry in workflow.entries(name):
            path = entry_path(entry, repository)
            text = (Path(repository) / path).read_text(errors="replace") if path else ""
            tests[(name, entry)] = (path, text)
    affected = set()
    for file in changed:
        mentions = [Path(file).name]
        parts = Path(file).parts
        if len(parts) > 3 and parts[0] == "testSuite":
            mentions.append("/".join(parts[1:-1]) + "/")
        found = {key for key, (path, text) in tests.items() if path == file or any(mention in text for mention in mentions)}
        if not found:
            return None
        affected |= found
    jobs = sorted({name for name, _ in affected})
    # An entry may be run by more than one job, and is enabled or disabled by its name alone, so it is skipped only if no
    # selected job runs it.
    enabled = {entry for _, entry in affected}
    skipped = sorted({entry for name in jobs for entry in workflow.entries(name)} - enabled)
    return jobs, skipped


def plan(workflow, event, ref="", draft=False, labels=(), changed=None, tier=None, tested_run_id="", run_id="",
         repository=REPOSITORY):
    """Return the plan for an event as a dictionary of the `Plan` job's outputs, plus a list of reasons for the summary.

    `changed` is the list of files changed by a pull request, or `None` if it could not be determined.
    """
    reasons = []
    deploying = event in ("push", "workflow_dispatch") and ref == "refs/heads/master"
    groups = set()
    suite = []
    selective = False
    if event == "workflow_dispatch":
        tier = tier or "full"
        reasons.append(f"Manually dispatched with tier `{tier}`.")
        if tier == "pr":
            groups = set(PR_GROUPS)
    elif event == "push":
        if tested_run_id:
            tier = "deploy"
            reasons.append(f"This commit was tested by the merge queue in run {tested_run_id}; deploying its artifacts.")
        else:
            tier = "full"
            reasons.append("This commit was not tested by the merge queue (or its artifacts are gone); running the full suite.")
    elif event == "merge_group":
        tier = "full"
        reasons.append("Merge group: running the full suite against the state of the branch this group will land on.")
    elif event == "pull_request":
        if draft:
            tier = "skip"
            reasons.append("Draft pull request: nothing is run until it is marked ready for review.")
        elif LABEL_FULL in labels:
            tier = "full"
            reasons.append(f"The pull request carries the `{LABEL_FULL}` label.")
        elif changed is None:
            tier = "full"
            reasons.append("The files changed by this pull request could not be determined; running the full suite.")
        else:
            decisive = {file: classify(file) for file in changed}
            forcing = sorted(file for file, selection in decisive.items() if selection == FULL)
            if forcing:
                tier = "full"
                reasons.append("Files whose effect can not be predicted from their path were changed: "
                               + ", ".join(f"`{file}`" for file in forcing[:10]) + (", ..." if len(forcing) > 10 else "") + ".")
            else:
                for selection in decisive.values():
                    groups.update(selection)
                if groups:
                    tier = "pr"
                    relevant = [file for file, selection in decisive.items() if selection]
                    # A pull request which touches only the test suite (and files which need no build) runs just the tests
                    # that it affects.
                    suite = [file for file in relevant if file.startswith("testSuite/")]
                    selective = bool(suite) and all(file.startswith("testSuite/") or decisive[file] == ("tooling",)
                                                    for file in relevant)
                    reasons.append("Pull request subset selected by changed paths: " + ", ".join(sorted(groups)) + ".")
                else:
                    tier = "skip"
                    reasons.append("No changed file can affect the build or the tests.")
    else:
        tier = "full"
        reasons.append(f"Unrecognized event `{event}`; running the full suite.")
    if tier not in TIERS:
        raise ValueError(f"unknown tier '{tier}'")

    skipped = []
    if tier == "skip":
        jobs = []
    elif tier == "deploy":
        jobs = [name for name in workflow.names() if name in DEPLOY_JOBS]
    elif tier == "full":
        jobs = workflow.names()
    else:
        selected = [name for group in sorted(groups) for name in PR_GROUPS[group] if name in workflow.jobs]
        skipped = sorted(PR_DEFERRED_ENTRIES)
        if selective:
            affected = select_affected(workflow, suite, repository)
            if affected is None:
                reasons.append("A changed test suite file affects no identifiable test; running every pull request model test.")
            else:
                models, skipped_affected = affected
                selected = [name for name in selected if name not in PR_GROUPS["model"]] + models
                # A test which the pull request itself changes runs even if it is otherwise deferred as slow.
                changed_entries = {entry for name in models for entry in workflow.entries(name)
                                   if entry_path(entry, repository) in suite}
                skipped = sorted((set(skipped) | set(skipped_affected)) - changed_entries)
                reasons.append("Only the test suite changed: running the " + str(len(models)) + " job(s) which it affects.")
        jobs = workflow.closure(selected)
    if tier in ("pr", "full"):
        jobs = [name for name in jobs if name not in PUBLISH_JOBS or (tier == "full" and deploying)]
        jobs = [name for name in jobs
                if name not in PULL_REQUEST_JOBS or (event == "pull_request" and tier in PULL_REQUEST_JOBS[name])]
        jobs += [name for name, tiers in PULL_REQUEST_JOBS.items()
                 if event == "pull_request" and tier in tiers and name in workflow.jobs and name not in jobs]
    outputs = {
        "tier": tier,
        "jobs": json.dumps(jobs),
        "skipEntries": json.dumps(skipped),
        "deploy": "true" if deploying and tier in ("full", "deploy") else "false",
        # The run whose artifacts the deployment jobs download: the merge group's run in the `deploy` tier, else this run.
        "artifactRunId": str(tested_run_id if tier == "deploy" else run_id),
    }
    return outputs, reasons


def gate(needs, planned):
    """Check the results of the jobs which the gate needs against the plan. Returns a list of problems, empty for a pass."""
    problems = []
    if (needs.get("Plan") or {}).get("result") != "success":
        problems.append("the `Plan` job did not succeed, so nothing can be trusted to have run")
    for name, job in sorted(needs.items()):
        result = job.get("result")
        if result in ("failure", "cancelled"):
            problems.append(f"`{name}` {'failed' if result == 'failure' else 'was canceled'}")
        elif result == "skipped" and name in planned:
            problems.append(f"`{name}` was planned but did not run")
    return problems


def write_outputs(outputs):
    target = os.environ.get("GITHUB_OUTPUT")
    lines = "".join(f"{key}={value}\n" for key, value in outputs.items())
    if target:
        with open(target, "a") as file:
            file.write(lines)
    sys.stdout.write(lines)


def write_summary(text):
    target = os.environ.get("GITHUB_STEP_SUMMARY")
    if target:
        with open(target, "a") as file:
            file.write(text)


def read_changed_files(path):
    """Read a NUL-separated list of file names, as written by `git diff --name-only -z`. Returns `None` if there is no list."""
    if not path or not Path(path).is_file():
        return None
    return [name for name in Path(path).read_text(errors="surrogateescape").split("\0") if name]


def main(arguments=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)
    planner = commands.add_parser("plan", help="decide which jobs run")
    planner.add_argument("--event", required=True, help="the name of the event which triggered the workflow")
    planner.add_argument("--ref", default="", help="the ref the workflow runs on")
    planner.add_argument("--draft", default="false", help="'true' if the pull request is a draft")
    planner.add_argument("--labels", default="[]", help="the labels of the pull request, as a JSON list")
    planner.add_argument("--changed-files", help="a file holding the NUL-separated names of the files changed")
    planner.add_argument("--tier", default="", help="the tier requested by a manual dispatch")
    planner.add_argument("--tested-run-id", default="", help="the merge group run which tested this commit, if any")
    planner.add_argument("--run-id", default="", help="the ID of this run")
    planner.add_argument("--workflow", default=str(WORKFLOW), help="the workflow file to plan for")
    gater = commands.add_parser("gate", help="check the results of a run against its plan")
    gater.add_argument("--needs", default=os.environ.get("CI_NEEDS"),
                       help="the `needs` context of the gate job, as JSON (default: the CI_NEEDS environment variable)")
    gater.add_argument("--jobs", default=os.environ.get("CI_JOBS"),
                       help="the planned jobs, as a JSON list (default: the CI_JOBS environment variable)")
    options = parser.parse_args(arguments)

    if options.command == "plan":
        labels = json.loads(options.labels or "[]") or []
        outputs, reasons = plan(
            Workflow(options.workflow),
            options.event,
            ref=options.ref,
            draft=options.draft == "true",
            labels=labels,
            changed=read_changed_files(options.changed_files),
            tier=options.tier or None,
            tested_run_id=options.tested_run_id,
            run_id=options.run_id,
        )
        write_outputs(outputs)
        jobs = json.loads(outputs["jobs"])
        skipped = json.loads(outputs["skipEntries"])
        summary = [f"## CI plan: `{outputs['tier']}` tier", ""] + [f"- {reason}" for reason in reasons] + [""]
        summary.append(f"{len(jobs)} job(s) planned" + (f", with {len(skipped)} matrix entries left to the merge queue" if skipped else "") + ".")
        if jobs and outputs["tier"] != "full":
            summary += ["", "<details><summary>Planned jobs</summary>", ""] + [f"- `{name}`" for name in jobs] + ["", "</details>"]
        write_summary("\n".join(summary) + "\n")
        return 0

    if not options.needs:
        parser.error("the results to check must be given by --needs or CI_NEEDS")
    needs = json.loads(options.needs)
    planned = json.loads(options.jobs or "[]") if (needs.get("Plan") or {}).get("result") == "success" else []
    problems = gate(needs, planned)
    ran = sorted(name for name, job in needs.items() if job.get("result") == "success" and name != "Plan")
    if problems:
        write_summary("## :x: CI gate failed\n\n" + "".join(f"- {problem}\n" for problem in problems))
        for problem in problems:
            print("FAILED: " + problem.replace("`", ""))
        return 1
    write_summary(f"## :white_check_mark: CI gate passed\n\n{len(ran)} job(s) ran and succeeded.\n")
    print(f"SUCCESS: {len(ran)} job(s) ran and succeeded")
    return 0


if __name__ == "__main__":
    sys.exit(main())
