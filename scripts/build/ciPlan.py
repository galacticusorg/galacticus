#!/usr/bin/env python3
"""Plan which jobs of the CI/CD workflow (``.github/workflows/cicd.yml``) run for an event, and gate its result.

The full suite costs over a hundred runner-hours, so it is not run on every push to every pull request. Instead each run of the
workflow is given a *tier* by its ``Plan`` job, and every other job runs only if this script lists it:

``full``
    Every job. Used for pull requests which change the build system or carry the ``ci:full`` label, for manual dispatches, and
    for any push to ``master`` of content which has not already passed the full suite.
``pr``
    A subset chosen from the files a pull request changes: the two default Linux builds, the quick unit tests, and the model tests
    which run in a few minutes each. Tests too slow or too specialized for that (other build configurations, macOS, MPI,
    validation, benchmarks, profiling) are left to the full suite.
``deploy``
    A push to ``master`` of content which has already passed the full suite. Nothing is retested: what is published is rebuilt
    from the commit on ``master`` (an executable records the commit it was built from) and published along with the validation
    results, benchmarks and profiles of the run which tested that content.
``skip``
    Nothing: a draft pull request, one which changes no file that can affect the build or the tests, or a merge group whose
    content has already passed the full suite.
``blocked``
    Nothing, and the run fails: a merge group whose content has not passed the full suite.

"Content" here means the state of every file which can affect the build or the tests - every file other than those which the
path rules below say can not - and is identified by a digest of those files, its *key* (see ``content_key()``). A run of the full suite which
passes records the key of what it tested (the ``CI-Gate`` job uploads an artifact named for it), and the ``Plan`` job of a
later run looks that record up. For a pull request the content is that of the pull request merged into ``master`` as it was
when the run began; the merge queue forms the same content for that pull request for as long as nothing which can affect the
build or the tests has reached ``master`` in the meantime. So the full suite is run on a pull request (or on a rollup of
several), not by the merge queue: the queue allows a required check at most six hours to report, which the full suite can not
promise, and so a merge group is only ever verified against the record. A group whose content has no record is refused rather
than tested, which keeps untested content off ``master`` by construction.

A change to files which the rules ignore - documentation, for the most part - leaves the key as it was. A pull request which
makes only such changes therefore runs nothing, passes the merge queue on the record of the content already on ``master``, and
does not disturb the record of a pull request which is waiting to be queued. A file which matches no rule at all is taken to
affect everything: it is part of the content, and a change to it selects the full suite.

Content which has passed the full suite is not tested again by a pull request either. A push to a pull request which changes
only ignored files, after its full run has passed, runs nothing.

The ``CI-Gate`` job then checks the outcome against the plan. It is the one status check that needs to be required by branch
protection, and it reports for every run - including a skipped one - so that a required check is never left pending.

What a pull request runs is decided by the tables below, not by the workflow file. Three properties are deliberate:

* Selection fails safe. A test added to the matrix of a job in the `pr` tier runs for pull requests unless it is listed here
  as deferred; a changed file which matches no rule selects the full suite; and a list of changed files which could not be
  determined gets the full suite. (A *job* which this script
  does not know about is not in the `pr` tier: it runs in the full suite, which everything passes before it reaches
  ``master``, until it is added to a group below.)
* A pull request which adds a test runs it. A change to the workflow selects the full suite, with one exception: where the
  workflow differs from that of ``master`` only in the tests which its jobs run, the pull request stays in the `pr` tier and
  runs the tests which it adds as well - whichever job they were added to, and whether or not they are listed as deferred.
  Otherwise every pull request which came with a new test would run the full suite at each push.
* A selected job brings the jobs it ``needs`` with it, so a test can never be planned without its build.
* A job which was planned but did not run fails the gate, as does one which failed or was canceled.

Run with::

    ciPlan.py plan --event <event> [--ref <ref>] [--draft] [--labels <json>] [--changed-files <file>] [--tier <tier>]
                   [--base-workflow <file>] [--tested-run-id <id>] [--lookup-failed] [--fork] [--content <key>]
                   [--run-id <id>]
    ciPlan.py key
    ciPlan.py gate [--needs <json>] [--jobs <json>]
    ciPlan.py attribute <results file> ...
"""

import argparse
import copy
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[2]
WORKFLOW_PATH = ".github/workflows/cicd.yml"
WORKFLOW = REPOSITORY / WORKFLOW_PATH

TIERS = ("skip", "pr", "full", "deploy", "blocked")

# The label which runs the full suite on a pull request.
LABEL_FULL = "ci:full"

# Jobs which are not gated by the plan: the plan itself, and the gate which checks it.
UNGATED_JOBS = ("Plan", "CI-Gate")

# Jobs which publish something, and so run only for a push to (or a dispatch on) `master`.
PUBLISH_JOBS = ("Submit-Notarization-MacOS", "Deploy", "Update-Release-Tag")

# The jobs which build what is published. An executable records the commit it was built from, so these are built again from the
# commit on `master` even where its content has already been tested as part of a pull request.
DEPLOY_BUILD_JOBS = (
    # Not a build, but the job which makes the parameter catalog that is published - which takes it a few minutes.
    "Test-Preprocessor",
    "Build-Executable-Linux",
    "Build-Library-Linux",
    "Build-Tools",
    "Package-Datasets",
    "Build-Executable-MacOS-M1",
    "Build-Library-MacOS",
    "Build-Docker-Linux",
)

# The jobs of the `deploy` tier: those builds, `Validate` (which publishes to `gh-pages` the validation results and benchmarks of
# the run which tested this content), and the jobs which publish.
DEPLOY_JOBS = DEPLOY_BUILD_JOBS + ("Validate",) + PUBLISH_JOBS

# Jobs which comment on a pull request, and so run only for one. `Validation-Report` needs the validation results, which only
# the full suite produces.
PULL_REQUEST_JOBS = {"Watch-Updates": ("pr", "full"), "Validation-Report": ("full",)}

# The jobs of the `pr` tier, in groups which the path rules below select. The builds are not listed: they are brought in by the
# `needs` of whichever of these jobs are selected.
PR_GROUPS = {
    # Pure-Python checks: the preprocessor regression tests, parameter file validation, and the schema.
    "tooling": ("Test-Preprocessor", "Python-Utils-Tests"),
    # Unit tests of the non-static build. `Codes-Quick` runs every unit test which takes less than a few minutes in one job; the
    # few slow ones are the `Codes` matrix, which is left to the full suite.
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

# Matrix entries of `pr` tier jobs which are left to the full suite because they are slow: each of these took over five minutes
# to run (not counting job setup) when the tiers were introduced. Listing an entry here only defers it - it still runs in the
# full suite, which everything passes before it reaches `master` - so the list errs on the side of brevity, and an entry which is
# not listed runs on pull requests.
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

# What a changed file selects. The first rule whose pattern matches a file decides for that file. `FULL` selects the full suite,
# and is used where the effect of a change can not be predicted from its path: the build system, the code generator, the
# dependencies, and the workflows themselves. `IGNORE` is for a file which can not affect the build or the tests
# (documentation, for the most part): it selects nothing, and it is left out of the key which identifies tested content (see
# `content_key()`), so a change to it is never tested - by a pull request, by the merge queue, or on `master`. A file which
# matches no rule is treated as `FULL`, so that a new kind of file is tested until someone decides otherwise; a test checks that
# every file in the repository has a rule, so that the decision is made when the file is added. In a pattern `*` matches within
# one path component and `**` matches across them.
FULL = "full"
IGNORE = ()
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
    # The workflows which make up the suite, and the image which one of its jobs builds.
    # (A change to this one which only adds tests to its jobs, or removes them, is treated more leniently: see `plan()`.)
    (WORKFLOW_PATH, FULL),
    (".github/workflows/testModel.yml", FULL),
    (".github/workflows/testCode.yml", FULL),
    (".github/workflows/profile.yml", FULL),
    (".github/actions/**", FULL),
    ("Dockerfile", FULL),
    # The configuration of the Python package and of its tests.
    ("pyproject.toml", ("tooling", "model")),
    ("pytest.ini", ("tooling",)),
    # Other scripts. Those for the documentation and for development are used by no model; the developer tools are used by
    # some model tests (parameter migration, for one); the rest (analysis, emulation, and so on) are used by validations and
    # tests which the `pr` tier does not run.
    ("scripts/doc/**", ("tooling",)),
    ("scripts/aux/**", ("tooling", "model")),
    ("scripts/**", FULL),
    # Constraint pipelines and tutorials are run by validations and tests which the `pr` tier does not run.
    ("constraints/**", FULL),
    ("tutorials/**", FULL),
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
    # Files which can not affect the build or the tests. The workflows and actions which the suite uses have rules above, and
    # a test checks that it uses no other: the rest of `.github/` is other workflows, and templates.
    (".github/**", IGNORE),
    ("docs/**", IGNORE),
    ("plots/**", IGNORE),
    ("aux/words.aff", IGNORE),
    ("aux/words.dict", IGNORE),
    (".agents/**", IGNORE),
    (".claude/**", IGNORE),
    (".devcontainer/**", IGNORE),
    (".vscode/**", IGNORE),
    ("AGENTS.md", IGNORE),
    ("CITATION.cff", IGNORE),
    ("CONTRIBUTING.md", IGNORE),
    ("COPYING", IGNORE),
    ("README.md", IGNORE),
    (".fortls", IGNORE),
    (".gitattributes", IGNORE),
    (".gitconfig", IGNORE),
    (".gitignore", IGNORE),
    (".readthedocs.yaml", IGNORE),
    (".yamllint", IGNORE),
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


def rule_for(path):
    """Return the pattern of the rule which decides for a file, or `None` if no rule matches it."""
    for (pattern, _), (regex, _) in zip(PATH_RULES, _PATH_RULES):
        if regex.match(path):
            return pattern
    return None


def classify(path):
    """Return what a changed file selects: `FULL`, a tuple of `PR_GROUPS` names, or `IGNORE` for a file which matters not.

    A file which matches no rule selects `FULL`: what is not known to be harmless is tested.
    """
    for regex, selection in _PATH_RULES:
        if regex.match(path):
            return selection
    return FULL


class Workflow:
    """The parts of the workflow which planning needs: each job's `needs`, and the test files of the jobs which run them."""

    def __init__(self, path=WORKFLOW, document=None):
        if document is None:
            import yaml  # Imported here so that `gate`, which needs no workflow, runs without PyYAML.

            document = yaml.safe_load(Path(path).read_text())
        self.document = document
        self.jobs = document["jobs"]

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

    @staticmethod
    def _target(job):
        """Return where a job names the tests which it runs, as `(the parameter naming them, the matrix key if any)`."""
        parameters = job.get("with") or {}
        field = "file" if "file" in parameters else "executable"
        target = parameters.get(field)
        if not isinstance(target, str):
            return None, None
        matrixed = re.fullmatch(r"\$\{\{\s*matrix\.(\w+)\s*\}\}", target)
        return field, matrixed.group(1) if matrixed else None

    def entries(self, name):
        """Return the test files (or test executables) which a job runs: one per matrix entry, or the single one it names."""
        job = self.jobs[name]
        field, key = self._target(job)
        if field is None:
            return []
        if key is None:
            # A job which runs unit tests (and so names no test file) may name several executables, separated by white space.
            target = job["with"][field]
            return [target] if field == "file" else target.split()
        values = ((job.get("strategy") or {}).get("matrix") or {}).get(key, [])
        return [value for value in values if isinstance(value, str)] if isinstance(values, list) else []

    def without_entries(self):
        """Return the workflow's document with the tests which its jobs run removed: everything `entries()` does not report."""
        document = copy.deepcopy(self.document)
        for job in document["jobs"].values():
            field, key = self._target(job)
            if field is None:
                continue
            if key is None:
                job["with"][field] = ""
            else:
                matrix = (job.get("strategy") or {}).get("matrix")
                if isinstance(matrix, dict) and isinstance(matrix.get(key), list):
                    matrix[key] = [value for value in matrix[key] if not isinstance(value, str)]
        return document


def added_entries(base, workflow):
    """Compare a workflow with the version which a pull request is based on.

    Returns `{job: [the tests added to it]}` if the two differ only in the tests which their jobs run - so an empty dictionary
    if they differ only in comments or layout, or only by tests removed - and `None` if they differ in any other way.
    """
    if base.without_entries() != workflow.without_entries():
        return None
    added = {}
    for name in workflow.names():
        before = set(base.entries(name))
        new = [entry for entry in workflow.entries(name) if entry not in before]
        if new:
            added[name] = new
    return added


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


def content_key(repository=REPOSITORY, revision="HEAD"):
    """Return the key which identifies the content of a commit: a digest of every file which can affect the build or the tests.

    The digest is of the names, modes and contents (as git's object names) of the files which the path rules do not ignore, so
    two commits have the same key exactly when they differ, if at all, only in ignored files.
    """
    listing = subprocess.run(["git", "-C", str(repository), "ls-tree", "-r", "-z", "--full-tree", revision],
                             capture_output=True, check=True).stdout
    digest = hashlib.sha256()
    for entry in listing.split(b"\0"):
        if entry and classify(entry.split(b"\t", 1)[1].decode("utf-8", "surrogateescape")):
            digest.update(entry + b"\0")
    return digest.hexdigest()[:40]


def plan(workflow, event, ref="", draft=False, labels=(), changed=None, tier=None, tested_run_id="", content="", run_id="",
         repository=REPOSITORY, base=None, lookup_failed=False, fork=False):
    """Return the plan for an event as a dictionary of the `Plan` job's outputs, plus a list of reasons for the summary.

    `changed` is the list of files changed by a pull request, or `None` if it could not be determined, and `base` is the
    workflow as it is on the branch which the pull request is based on, or `None` if that could not be read. `content` is the
    key of the content which the run is for (see `content_key()`), and `tested_run_id` the run of the full suite which passed
    for that content, if there is one. `lookup_failed` says that the record of such a run could not be looked for at all, which
    is not the same as there being none. `fork` says that a pull request comes from a fork of the repository.
    """
    key = content
    content = f"this content (key `{key}`)" if key else "this content"
    reasons = []
    deploying = event in ("push", "workflow_dispatch") and ref == "refs/heads/master"
    groups = set()
    suite = []
    selective = False
    added = {}
    if event == "workflow_dispatch":
        tier = tier or "full"
        reasons.append(f"Manually dispatched with tier `{tier}`.")
        if tier == "pr":
            groups = set(PR_GROUPS)
    elif event == "push":
        if tested_run_id:
            tier = "deploy"
            reasons.append(f"The full suite passed for {content} in run {tested_run_id}: rebuilding what is published from this "
                           "commit, and publishing it with the results of that run.")
        elif lookup_failed:
            tier = "full"
            reasons.append(f"The record of a run of the full suite for {content} could not be looked for (the GitHub API "
                           "returned an error); running the full suite.")
        else:
            tier = "full"
            reasons.append(f"No run of the full suite has passed for {content} (or its record is gone); running the full suite.")
    elif event == "merge_group":
        if tested_run_id:
            tier = "skip"
            reasons.append(f"Merge group: the full suite passed for {content} in run {tested_run_id}, so there is nothing to run.")
        elif lookup_failed:
            tier = "blocked"
            reasons.append(f"Merge group: the record of a run of the full suite for {content} could not be looked for (the GitHub "
                           "API returned an error). This says nothing about the pull request: add it to the merge queue again.")
        else:
            tier = "blocked"
            reasons.append(f"Merge group: no run of the full suite has passed for {content}, and the merge queue does not run "
                           "it. Merge `master` into the pull request (or into a rollup of several), run the full suite on it "
                           f"(the `{LABEL_FULL}` label), and queue it again once that has passed. If the full suite did pass "
                           "for this pull request, then `master` has since gained a change which can affect the build or the "
                           "tests.")
    elif event == "pull_request":
        if draft:
            tier = "skip"
            reasons.append("Draft pull request: nothing is run until it is marked ready for review.")
        elif tested_run_id:
            tier = "skip"
            reasons.append(f"The full suite passed for {content} in run {tested_run_id}, so there is nothing to run.")
        elif LABEL_FULL in labels:
            tier = "full"
            reasons.append(f"The pull request carries the `{LABEL_FULL}` label.")
        elif changed is None:
            tier = "full"
            reasons.append("The files changed by this pull request could not be determined; running the full suite.")
        else:
            decisive = {file: classify(file) for file in changed}
            if decisive.get(WORKFLOW_PATH) == FULL and base is not None:
                # A workflow which differs only in the tests that its jobs run does not need the full suite: only those tests.
                tests = added_entries(base, workflow)
                if tests is not None:
                    added = tests
                    decisive[WORKFLOW_PATH] = ()
                    count = sum(len(entries) for entries in added.values())
                    reasons.append("The workflow changed only in the tests which its jobs run"
                                   + (f": the {count} added will be run (" + ", ".join(f"`{entry}`" for entries in added.values()
                                                                                       for entry in entries) + ")." if count
                                      else ", and none was added."))
            forcing = sorted(file for file, selection in decisive.items() if selection == FULL)
            if forcing:
                tier = "full"
                reasons.append("Files whose effect can not be predicted from their path were changed: "
                               + ", ".join(f"`{file}`" for file in forcing[:10]) + (", ..." if len(forcing) > 10 else "") + ".")
            else:
                for selection in decisive.values():
                    groups.update(selection)
                if groups or added:
                    tier = "pr"
                    relevant = [file for file, selection in decisive.items() if selection]
                    # A pull request which touches only the test suite (and files which need no build) runs just the tests
                    # that it affects.
                    suite = [file for file in relevant if file.startswith("testSuite/")]
                    selective = bool(suite) and all(file.startswith("testSuite/") or decisive[file] == ("tooling",)
                                                    for file in relevant)
                    if groups:
                        reasons.append("Pull request subset selected by changed paths: " + ", ".join(sorted(groups)) + ".")
                else:
                    tier = "skip"
                    reasons.append("No changed file can affect the build or the tests.")
    else:
        tier = "full"
        reasons.append(f"Unrecognized event `{event}`; running the full suite.")
    if tier not in TIERS:
        raise ValueError(f"unknown tier '{tier}'")
    if event == "pull_request" and fork and tier != "skip":
        reasons.append("This pull request is from a fork, so no run of it can make a record which the merge queue will believe"
                       + (", this run included" if tier == "full" else "") + ". To land it, merge its branch into a branch of "
                       f"this repository (a rollup) and run the full suite on a pull request for that branch (the `{LABEL_FULL}` "
                       "label).")

    skipped = []
    if tier in ("skip", "blocked"):
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
                # A test which the pull request affects - which it changes, or a file of which it changes - runs even if it
                # is otherwise deferred as slow.
                affected_entries = {entry for name in models for entry in workflow.entries(name)} - set(skipped_affected)
                skipped = sorted((set(skipped) | set(skipped_affected)) - affected_entries)
                reasons.append("Only the test suite changed: running the " + str(len(models)) + " job(s) which it affects.")
        if added:
            # Tests which the pull request adds to the workflow run, in whichever job they were added to. A job brought in
            # for that alone runs nothing else. A test is enabled or disabled by its name alone, in every job which names it,
            # so one is skipped only if no selected job is to run it.
            wanted = {entry for name in workflow.closure(selected) for entry in workflow.entries(name)} - set(skipped)
            wanted |= {entry for entries in added.values() for entry in entries}
            selected = selected + [name for name in added if name not in selected]
            named = {entry for name in workflow.closure(selected) for entry in workflow.entries(name)}
            skipped = sorted((set(skipped) | named) - wanted)
        # A job which runs tests, but every one of whose tests is skipped, would run nothing: leave it out, with its builds.
        idle = [name for name in selected if workflow.entries(name) and not set(workflow.entries(name)) - set(skipped)]
        if idle:
            selected = [name for name in selected if name not in idle]
            reasons.append("Left out, as every test which they run is left to the full suite: "
                           + ", ".join(f"`{name}`" for name in idle) + ".")
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
        # The run whose validation results, benchmarks and profiles are published: in the `deploy` tier the run which tested
        # this content, else this run.
        "resultsRunId": str(tested_run_id if tier == "deploy" else run_id),
        # The key of the content which this run is for, under which `CI-Gate` records a pass of the full suite.
        "content": key,
    }
    return outputs, reasons


def gate(needs, planned):
    """Check the results of the jobs which the gate needs against the plan. Returns a list of problems, empty for a pass."""
    problems = []
    if (needs.get("Plan") or {}).get("result") != "success":
        problems.append("the `Plan` job did not succeed (its summary says why), so nothing can be trusted to have run")
    for name, job in sorted(needs.items()):
        result = job.get("result")
        if result in ("failure", "cancelled"):
            problems.append(f"`{name}` {'failed' if result == 'failure' else 'was canceled'}")
        elif result == "skipped" and name in planned:
            problems.append(f"`{name}` was planned but did not run")
    return problems


# What a file of validation results begins with: it is JavaScript, which assigns the results to a variable.
RESULTS_PREFIX = re.compile(r"window\.\w+ = ")


def head_commit(repository=REPOSITORY):
    """Describe the commit which is checked out, in the form in which a file of validation results describes a commit."""
    fields = subprocess.run(["git", "-C", str(repository), "log", "-1", "--format=%H%x00%an%x00%ae%x00%cI%x00%B"],
                            capture_output=True, text=True, check=True).stdout.split("\0")
    revision, name, email, committed, message = fields
    from datetime import datetime

    return {
        "author": {"name": name, "email": email},
        "id": revision,
        # A commit message as `git log` prints it ends with one more newline than the message itself.
        "message": message[:-1] if message.endswith("\n") else message,
        "timestamp": str(datetime.fromisoformat(committed)),
        "url": "https://github.com/galacticusorg/galacticus/commit/" + revision,
    }


def attribute(text, commit):
    """Return the text of a file of validation results with the commit which it describes replaced by `commit`.

    The results published by a `deploy` tier run are those of the run which tested the same content as part of a pull request,
    and so describe the commit which that run had checked out - a merge of the pull request made for the test, which is on no
    branch. They are attributed instead to the commit on `master` which has that content.
    """
    prefix = RESULTS_PREFIX.match(text)
    if not prefix:
        raise ValueError("not a file of validation results")
    data = json.loads(text[prefix.end():])
    if not isinstance(data, dict) or not isinstance(data.get("commit"), dict):
        raise ValueError("the results describe no commit")
    data["commit"] = commit
    return prefix.group(0) + json.dumps(data, indent=4, ensure_ascii=False)


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
    planner.add_argument("--base-workflow", help="a file holding the workflow as it is on the branch a pull request is based on")
    planner.add_argument("--tier", default="", help="the tier requested by a manual dispatch")
    planner.add_argument("--tested-run-id", default="", help="the run of the full suite which passed for this content, if any")
    planner.add_argument("--lookup-failed", default="false", help="'true' if the record of such a run could not be looked for")
    planner.add_argument("--fork", default="false", help="'true' if the pull request comes from a fork of the repository")
    planner.add_argument("--content", default="", help="the key of the content which this run is for")
    planner.add_argument("--run-id", default="", help="the ID of this run")
    planner.add_argument("--workflow", default=str(WORKFLOW), help="the workflow file to plan for")
    gater = commands.add_parser("gate", help="check the results of a run against its plan")
    gater.add_argument("--needs", default=os.environ.get("CI_NEEDS"),
                       help="the `needs` context of the gate job, as JSON (default: the CI_NEEDS environment variable)")
    gater.add_argument("--jobs", default=os.environ.get("CI_JOBS"),
                       help="the planned jobs, as a JSON list (default: the CI_JOBS environment variable)")
    attributer = commands.add_parser("attribute", help="attribute files of validation results to the commit which is checked out")
    attributer.add_argument("files", nargs="+", help="the files of results to rewrite in place")
    commands.add_parser("key", help="print the key which identifies the content of the commit which is checked out")
    options = parser.parse_args(arguments)

    if options.command == "key":
        print(content_key())
        return 0

    if options.command == "attribute":
        commit = head_commit()
        for name in options.files:
            file = Path(name)
            file.write_text(attribute(file.read_text(encoding="utf-8"), commit), encoding="utf-8")
            print(f"{name}: attributed to {commit['id']}")
        return 0

    if options.command == "plan":
        labels = json.loads(options.labels or "[]") or []
        base = None
        if options.base_workflow and Path(options.base_workflow).is_file():
            try:
                base = Workflow(options.base_workflow)
            except Exception as error:  # Whatever is wrong with it, the workflows can not be compared: the full suite is run.
                print(f"The workflow of the base branch could not be read ({error}).", file=sys.stderr)
        outputs, reasons = plan(
            Workflow(options.workflow),
            options.event,
            ref=options.ref,
            draft=options.draft == "true",
            labels=labels,
            changed=read_changed_files(options.changed_files),
            tier=options.tier or None,
            tested_run_id=options.tested_run_id,
            content=options.content,
            run_id=options.run_id,
            base=base,
            lookup_failed=options.lookup_failed == "true",
            fork=options.fork == "true",
        )
        write_outputs(outputs)
        jobs = json.loads(outputs["jobs"])
        skipped = json.loads(outputs["skipEntries"])
        summary = [f"## CI plan: `{outputs['tier']}` tier", ""] + [f"- {reason}" for reason in reasons] + [""]
        summary.append(f"{len(jobs)} job(s) planned" + (f", with {len(skipped)} matrix entries left to the full suite" if skipped else "") + ".")
        if jobs and outputs["tier"] != "full":
            summary += ["", "<details><summary>Planned jobs</summary>", ""] + [f"- `{name}`" for name in jobs] + ["", "</details>"]
        write_summary("\n".join(summary) + "\n")
        if options.fork == "true" and options.event == "pull_request" and outputs["tier"] != "skip":
            print("::warning title=Pull request from a fork::" + reasons[-1].replace("`", ""))
        if outputs["tier"] == "blocked":
            # Failing here fails the gate, which is what removes the pull request from the merge queue.
            for reason in reasons:
                print("::error title=Not tested::" + reason.replace("`", ""))
            return 1
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
