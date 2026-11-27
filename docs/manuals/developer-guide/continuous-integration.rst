Testing and Continuous Integration
==================================

This page describes how Galacticus is tested, both locally and through the
automated continuous-integration (CI) pipeline that runs on every pull request.

Test suites
-----------

Galacticus has two complementary sets of tests.

The model regression suite (``testSuite/``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``testSuite/`` directory contains the regression tests that exercise the
compiled model — building small models, checking physics results against known
answers, validating file formats, and so on. They are orchestrated by
``testSuite/test-all.py``, which runs every ``testSuite/test-*.py`` script in
turn and writes a log for each to ``testSuite/outputs/``:

.. code-block:: bash

   export GALACTICUS_EXEC_PATH=`pwd`
   python3 testSuite/test-all.py

An individual test can be run directly, which is usually what you want while
iterating on a specific area:

.. code-block:: bash

   python3 testSuite/test-Python-interface.py

Most of these tests require a built ``Galacticus.exe`` and the run-time datasets
(``GALACTICUS_DATA_PATH``); see :doc:`../user-guide/installation/index`.

Three scripts are deliberately not run by any CI workflow, and are run manually
(or through ``test-all.py``) only: ``test-regressions.py``, which is superseded
by the ``testSuite/regressions/`` matrix, and
``validate-concentrationConvergence.py`` and ``validate-mergerTreeConvergence.py``,
which are long validation runs submitting jobs through a queue manager.

.. _test-outcome-markers:

Test outcome markers
~~~~~~~~~~~~~~~~~~~~~~

A test script reports its outcome by printing a marker at the start of a line,
and then exiting with status zero:

``FAILED: <what went wrong>``
   The test failed.

``SUCCESS: <what passed>``
   The test passed.

``SKIPPED: <why>``
   A prerequisite could not be met — too few MPI processes, a build without
   debugging enabled, a missing executable — so the test ran nothing at all.
   ``test-all.py`` reports these as a third outcome, and ``testModel.yml``
   surfaces them in the job summary, so that a misconfigured runner cannot
   silently stop testing.

Both harnesses judge a test by scanning its log for ``FAIL`` as a substring, and
both also fail a test whose script exits non-zero. An uncaught exception is
therefore detected — but as a traceback rather than as a readable message, which
is why failures are reported through a marker and a zero exit instead. The one
failure that genuinely disappears is a path which exits zero without printing a
marker, so every early-exit guard (a missing input file, a model that did not
run) must print ``FAILED:`` first.

Because the match is on the substring ``FAIL`` anywhere in the log, a *passing*
run must never emit it: keep ``FAILED`` at the start of its line, and do not
print counts or summaries containing it. For the same reason ``SKIPPED`` is
reserved for a whole-test skip; a sub-check which is skipped should say so in
plain words.

The convention is enforced by ``scripts/aux/auditTestMarkers.py``, which audits
every ``test-*.py`` and ``validate-*.py`` script in well under a second:

.. code-block:: bash

   ./scripts/aux/auditTestMarkers.py           # report
   ./scripts/aux/auditTestMarkers.py --check   # exit 1 if anything is reported

The *Validate-Test-Markers* pull-request check runs it with ``--check``. Unlike
the source-formatting job this check is blocking: the suite is fully compliant,
so anything it reports is new drift rather than a pre-existing condition. A
string literal which must contain a marker word without being a marker (a
``grep`` pattern, say) is exempted by ending its line with the pragma
``# markers: exempt``.

.. _test-independent-references:

Writing an independent reference for a test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A test that compares Galacticus against an independent calculation (a Python
reference script, say) can disagree for reasons that look exactly like a
physics error but have nothing to do with the physics. Their signature is a
disagreement that does not shrink when the numerical tolerances are tightened.
The usual causes are:

* **Unit definitions.** Copy them from Galacticus; do not choose your own. In
  Galacticus a Gyr is :math:`10^9` *sidereal* years
  (``3.15581497635456e16`` s), not Julian years (``3.1556952e16`` s), a
  difference of :math:`3.8\times 10^{-5}` which enters squared wherever a rate is
  converted from km/s/Mpc to Gyr\ :sup:`-1`. An Mpc is
  ``3.08567758135e22`` m. Both are defined in
  ``source/numerical/constants/astronomical.F90``. Measured constants (the
  gravitational constant, atomic data) may legitimately come from published standard
  values, which keeps the reference independent; a unit definition may not.

* **Epoch.** A virial radius defined by a ``virialDensityContrast`` relative to
  the *mean* matter density depends on the expansion factor at the node's time.
  At :math:`t=13.8` Gyr under ``matterLambda`` the expansion factor is 0.99997,
  so the mean density is :math:`10^{-4}` higher than today and the virial radius
  :math:`3.3\times 10^{-5}` smaller; under ``staticUniverse`` the expansion
  factor is exactly 1.

* **Initial conditions.** Take them from the merger tree, not from
  recomputing them. A hand-written tree states, for example, a satellite's position
  and velocity, and the model integrates *those*. Recomputing "the same"
  quantity in the reference can give a slightly different value (if the tree
  carried a ``matterLambda`` virial radius while the model ran a static
  universe, say), and the difference grows along the orbit.

* **Dark matter fraction.** When baryonic components are present, the dark
  matter halo holds :math:`1-\Omega_\mathrm{b}/\Omega_\mathrm{M}` of the node
  mass, not all of it: an NFW reference normalized to the full basic mass is too
  heavy by :math:`\Omega_\mathrm{b}/\Omega_\mathrm{M}` (about 16%). Check the
  enclosed mass in the test itself, with the ``massProfile`` extractor (e.g.
  ``virialRadius:all:all:0.2``).

When a disagreement will not shrink as a tolerance is tightened, suspect an
input mismatch before the physics, and print the inputs the two calculations
share to find it.

The Python unit tests (``pytest``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The build system's Python modules (under ``python/``) and the documentation
tooling (under ``scripts/doc/``) have unit tests that are collected by ``pytest``
(the discovered paths are configured in ``pytest.ini``):

.. code-block:: bash

   export GALACTICUS_EXEC_PATH=`pwd`
   python -m pytest -v

These run quickly and do not require a built executable. Note that a bare
``pytest`` does **not** run the ``testSuite/`` model regression tests described
above — those are a separate suite.

The CI pipeline
---------------

The CI workflows live in ``.github/workflows/``. The two most relevant to
contributors are described below; the others are scheduled (cron) or utility
jobs (for example ``linkCheck.yml``, ``bibliographyUpdater.yml``,
``cloudyTableUpdate.yml``, ``slocReport.yml``, and ``notarizeMacOS.yml``).
One, ``largerRunnerBuild.yml``, runs only when dispatched manually or called by
another workflow: it builds on a GitHub-hosted larger runner, which has the
memory for a build at high parallelism but, unlike the standard runners, is
billed by the minute. ``libraryRelease.yml`` calls it for each versioned
release, to build the library with link-time optimization - which needs more
memory than a standard runner has - and publish it to that release.

Pull-request checks (``prChecks.yml``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These fast checks run on every pull request, mostly over the files changed in
the PR:

* **Validate-YAML / Validate-XML / Validate-Python-Scripts** — syntax and
  well-formedness of changed YAML, XML, and Python files.
* **Fortran-Static-Analysis** — ``scripts/aux/staticAnalyzer.py`` checks for
  common Fortran issues (empty constructors/destructors, duplicate variables in
  ``constructorAssign`` directives, destructors and internal constructors which
  are never bound, and similar).
* **Validate-Docstrings-RST** — ``scripts/doc/convertDocstringsToRST.py --check``
  ensures every embedded docstring is reStructuredText (no old-style LaTeX).
* **Validate-Test-Markers** — ``scripts/aux/auditTestMarkers.py --check`` audits
  the ``testSuite`` scripts against the outcome-marker convention described in
  :ref:`test-outcome-markers`.
* **Spell-Check-RST** — builds the documentation with ``sphinxcontrib-spelling``
  and reports possible misspellings as an advisory PR comment; add legitimate
  technical terms to ``aux/words.dict``.

Build, test, and benchmark (``cicd.yml``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main pipeline builds Galacticus on Linux (and macOS), runs the preprocessor
and Python regression tests, and executes model and benchmark runs. It makes use
of the reusable ``testCode.yml`` and ``testModel.yml`` workflows and compares
outputs against references (``hdf5Diff.yml``). This pipeline is heavier and takes
considerably longer than the pull-request checks.

The full suite costs over a hundred hours of runner time, so it is not run for
every push to every pull request. The first job of each run, **Plan**, gives the
run a *tier*, and every other job runs only if the plan lists it:

* **full** — every job. Used for a pull request which changes the build system,
  the code generator, the dependencies, or the workflows themselves (other than
  by adding tests to the jobs of ``cicd.yml``: see below), for a pull request
  carrying the ``ci:full`` label, for a manual dispatch, and for a push to
  ``master`` of content which has not already passed the full suite.
* **pr** — a subset chosen from the files which a pull request changes. A change
  to the source runs the two default Linux builds, the unit tests, and the model
  tests which take no more than a few minutes each. A change to the test suite
  alone runs just the tests which it affects. Other build configurations, macOS,
  MPI, validation, benchmarks, and profiling are left to the full suite.
* **deploy** — a push to ``master`` of content which has already passed the full
  suite. Nothing is tested again. What is published is built again from the
  commit on ``master``, because an executable records the commit from which it
  was built, and is published together with the validation results, benchmarks,
  and profiles of the run which tested that content.
* **skip** — nothing: a draft pull request, one which changes no file that can
  affect the build or the tests, or a merge group whose content has already
  passed the full suite.
* **blocked** — nothing, and the run fails: a merge group whose content has not
  passed the full suite.

The last job, **CI-Gate**, checks the outcome against the plan: it passes only if
no job failed or was canceled and every planned job actually ran. It is the one
status check which branch protection should require, and it reports for every
run — including a skipped one.

Getting a change onto ``master``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Everything passes the full suite before it reaches ``master``, but the merge
queue does not run it. The queue allows a required check at most six hours to
report, and the full suite takes longer than that whenever runners are scarce.
Instead, a run of the full suite which passes records the *content* which it
tested, and the merge queue admits only a merge group whose content has such a
record.

Content here means the state of every file which can affect the build or the
tests: every file which matches one of the path rules in
``scripts/build/ciPlan.py``. It is identified by a digest of those files (print
it for a checkout with ``scripts/build/ciPlan.py key``). Two commits which
differ only in other files — documentation, for the most part — have the same
content. For a pull request, the content tested is that of the pull request
merged into ``master`` as ``master`` was when the run began. The merge queue
forms that same content for the pull request for as long as nothing which can
affect the build or the tests has reached ``master`` in the meantime.

So, to land one pull request, or several together:

#. Let each pull request pass its own checks in the **pr** tier.
#. Choose what will run the full suite. For a single pull request, give it the
   ``ci:full`` label. For several, make a new branch from ``master``, merge
   each of their branches into it, and open a pull request for that branch with
   the ``ci:full`` label. The label is read when a run begins: add it before
   opening the pull request, or add it and then push.
#. When that run has passed, add the pull request to the merge queue. Its merge
   group is verified against the record in about a minute. Merging the combined
   pull request with a merge commit also marks the pull requests which it
   combines as merged.
#. The push to ``master`` which follows runs the **deploy** tier.

A pull request which changes no file that can affect the build or the tests
needs none of this. Its content is the content already on ``master``, which has
a record, so it can be queued as soon as its own checks have passed, and it
does not disturb a pull request which is waiting to be queued.

If a change which can affect the build or the tests reached ``master`` after
the run of the full suite began, the merge group is of content which no run has
tested: its run fails, and the pull request is removed from the queue. Merge
``master`` into the branch, let the full suite run again, and queue it once
more. Two such pull requests queued at once fare the same way: the second is of
content which includes the first.

A record is an artifact of the run which made it, and expires with that run's
other artifacts. Only a record made by a run in this repository is believed: a
pull request from a fork runs its own copy of the workflow. A branch of this
repository could likewise change the workflow to record a pass which it had not
earned, so a pull request which changes ``cicd.yml`` deserves the same review as
any other change to what guards ``master``.

Changing what runs
^^^^^^^^^^^^^^^^^^

The rules which choose a tier and a subset are tables in
``scripts/build/ciPlan.py``, which is the place to change what a pull request
runs. Its path rules must cover every file which can affect the build or the
tests, because a file which matches no rule is also left out of the content
described above, and so a change to it is never tested. When a new kind of file
comes to be used by the build or by a test, give it a rule.

A pull request which adds a test to ``cicd.yml`` — an executable in the list of
``Codes-Quick``, or an entry in the matrix of a job — does not thereby run the
full suite. The plan compares the workflow with that of ``master``, and where
the two differ only in the tests which their jobs run, the pull request stays in
the **pr** tier and runs the tests which it adds as well: whichever job they were
added to, and whether or not they are listed as deferred. A job brought in for
an added test alone runs nothing else. Any other change to the workflow, however
small, runs the full suite.

When adding a job to ``cicd.yml``, make it need the ``Plan`` job, give it the
same ``if:`` condition as its neighbors (with its own name), and add it to the
``needs`` of ``CI-Gate``. The tests in ``scripts/build/tests/test_ciPlan.py``
fail if any of these is missed. A new job runs only in the full tier unless it is
added to a group in ``ciPlan.py``. A new unit test executable is added to the
list of the ``Codes-Quick`` job, which runs all of the quick unit tests one after
another in a single job; a unit test which takes several minutes or more belongs
instead in the matrix of the ``Codes`` job, where each gets a job of its own. A
test added to that matrix also needs a step in the
``Build-Executable-Linux-Non-Static`` job which uploads its executable as an
artifact of its own, so that its job need not download every test executable.

A test which is known to fail by chance can be given more than one attempt with
the ``attempts`` input of ``testModel.yml``. A test which passes only on a later
attempt is reported as a warning on the run, so that it remains visible.

.. _git-hooks:

Git hooks
---------

Galacticus shares a set of `git hooks <https://github.com/galacticusorg/gitHooks>`_
across its repositories. Installing them is optional but recommended: they catch
many of the same problems as the CI checks above, but *before* a commit is made.
The hooks include:

* **pre-commit** — a dispatcher (``pre-commit``) that runs every executable in
  ``pre-commit.d/`` in order. The main check script (``00-galacticus``) performs
  Fortran static analysis, validates ``.bib``, XML, YAML and Python files,
  spell-checks documentation, compiles embedded XML directives, flags leftover
  debugging statements, and verifies that the generated parameter-file schema
  (``schema/parameters.xsd``) is up to date with the sources — regenerating it
  only when a staged change could affect it (a Fortran source, the schema
  generator, the catalog builder, or the schema itself). A second script
  (``01-claude-review``) optionally reviews the staged diff with ``claude``.
* **commit-msg** / **prepare-commit-msg** — enforce and (optionally, with
  ``claude``) pre-fill
  `Conventional Commits <https://www.conventionalcommits.org>`_ messages.
* **pre-push** — asks for confirmation before pushing directly to
  ``master``/``main``.

To install them, clone the ``gitHooks`` repository and symlink it in place of
your clone's ``.git/hooks`` directory (the ``pre-commit`` dispatcher expects to
find ``pre-commit.d/`` under ``.git/hooks``):

.. code-block:: bash

   git clone https://github.com/galacticusorg/gitHooks.git /path/to/gitHooks
   cd /path/to/your/galacticus
   rm -rf .git/hooks
   ln -s /path/to/gitHooks .git/hooks

The two hooks that call ``claude`` spend the committer's own Claude usage, so
they are opt-in: set ``GALACTICUS_CLAUDE_HOOKS=1`` (or
``git config galacticus.claudeHooks true``) to enable them.

The hooks degrade gracefully when optional tools (for example ``hunspell``,
``yamllint``, ``lxml``, or ``claude``) are not installed; see the
`gitHooks README <https://github.com/galacticusorg/gitHooks/blob/main/README.md>`_
for the full list of checks and configuration options.

Before opening a pull request
-----------------------------

At a minimum, build Galacticus, run the quick test, and run the Python unit
tests; if your change touches the model, run the relevant ``testSuite/test-*.py``
script(s) too. See the repository's `CONTRIBUTING.md
<https://github.com/galacticusorg/galacticus/blob/master/CONTRIBUTING.md>`_ for
the full checklist.

.. _manual-sec-mergingPullRequests:

Merging pull requests
---------------------

Pull requests are merged with a **merge commit**. Do not squash-merge or
rebase-merge them.

This is not merely a matter of taste in history: parameter file migrations
(:galacticus-ref:`parameterMigrations`) are keyed by the full hash of the commit
that introduced each change, and ``parametersMigrate.py`` matches those hashes
against the repository's git ancestry by exact string equality. Squash- and
rebase-merging both rewrite commit hashes, which would leave every affected
``migrations.xml`` entry pointing at a commit that no longer exists in the
history. Such a migration matches nothing and is silently skipped — no error is
raised, and parameter files quietly stop migrating correctly.

Preserving commit hashes on merge is therefore a correctness requirement, and
applies to every pull request, since a later branch may add a migration keyed to
a commit merged earlier.
