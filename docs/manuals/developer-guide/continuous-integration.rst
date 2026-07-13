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

Pull-request checks (``prChecks.yml``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These fast checks run on every pull request, mostly over the files changed in
the PR:

* **Validate-YAML / Validate-XML / Validate-Python-Scripts** — syntax and
  well-formedness of changed YAML, XML, and Python files.
* **Fortran-Static-Analysis** — ``scripts/aux/staticAnalyzer.py`` checks for
  common Fortran issues (empty constructors/destructors, duplicate variables in
  ``constructorAssign`` directives, and similar).
* **Validate-Docstrings-RST** — ``scripts/doc/convertDocstringsToRST.py --check``
  ensures every embedded docstring is reStructuredText (no old-style LaTeX).
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
  (``01-claude-review``) optionally reviews the staged diff.
* **commit-msg** / **prepare-commit-msg** — enforce and pre-fill
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
