Versions and Releases
=====================

Semantic Versioning
-------------------

Galacticus follows the concepts of `semantic versioning <https://semver.org/>`_, such that version numbers are of the form ``X.Y.Z`` where ``X``, ``Y``, and ``Z`` are integers which are incremented according to the following rules:

#. The major version, ``X``, is incremented when incompatible `API <https://en.wikipedia.org/wiki/API>`_ changes are made;
#. The minor version, ``Y``, is incremented when adding functionality that is backward compatible; and,
#. The patch version, ``Z``, is incremented when a backwards compatible bug fix is made.

The "API" for Galacticus consists of:

#. Inputs:

   #. The parameter file passed to Galacticus on the command line;
   #. Any data files (e.g. merger tree files) used by Galacticus.

#. Outputs:

   #. The main output HDF5 file;
   #. Any auxiliary files that are output.

Changes to the format or syntax of any of these files are considered to be incompatible API changes. Note that additions to the formats that don't break backward compatibility (e.g. adding a new class with new parameters) *are* backwards compatible and so do not require a major version increment.

Versions are implemented through GitHub's "`release <https://docs.github.com/en/github/administering-a-repository/releasing-projects-on-github/about-releases>`_" mechanism, and so are based on ``git`` "`tags <https://git-scm.com/book/en/v2/Git-Basics-Tagging>`_". Associated with each release are a statically-linked binary executable and documentation. Additionally a Docker image is built for each version, which can be retrieved using ``docker pull galacticusorg/galacticus:X.Y.Z``.

The Galacticus ``datasets`` repo has corresponding versions - it's recommended to using matching versions of the `galacticus <https://github.com/galacticusorg/galacticus>`_ and `datasets <https://github.com/galacticusorg/galacticus>`_ repos.

When to Create a Release
------------------------

Day-to-day development is delivered continuously through the rolling
``bleeding-edge`` release (updated by CI on every push to ``master``), so a
tagged release is a deliberate, stable checkpoint rather than something cut for
every change. Create one when:

#. a backward-compatible **bug fix** is ready that users should be able to pin to
   — increment the patch version ``Z``;
#. new **backward-compatible functionality** has landed (e.g. a new class,
   parameter, or output) — increment the minor version ``Y`` and reset ``Z`` to
   ``0``; or
#. an **incompatible API change** has been made to the parameter-file syntax,
   input data formats, or output HDF5 format (see the API definition above) —
   increment the major version ``X`` and reset ``Y`` and ``Z`` to ``0``.

A practical cadence is to tag a minor release periodically once a useful batch of
features has accumulated and the test suite is green, and to tag a patch release
promptly when a fix matters to users. Always tag from a ``master`` commit whose
CI is passing, so the binaries on ``bleeding-edge`` (which the versioned release
mirrors) correspond to the tagged code.

Determining the Exact Version Used
----------------------------------

You can determine the exact version (i.e. the Git revision hash) of a Galacticus executable by using the `report <https://galacticus.readthedocs.io/en/latest/manuals/developer-guide/index.html>`_ task. To do this, simply run the parameter file `report.xml <https://github.com/galacticusorg/galacticus/blob/master/parameters/report.xml>`_:

.. code-block:: bash

   ./Galacticus.exe parameters/report.xml

You will see output similar to:

.. code-block:: text

                 ##
      ####        #                  #
     #   #        #             #
    #       ###   #  ###   ### ###  ##   ### ## ##   ##
    #       #  #  #  #  # #  #  #    #  #  #  #  #  #
    #   ###  ###  #   ### #     #    #  #     #  #   #
     #   #  #  #  #  #  # #     #    #  #     #  #    #
      ####  #### ### ####  ###   ## ###  ###   #### ##

    © 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
      2017, 2018, 2019, 2020, 2021, 2022
      - Andrew Benson

   MM: -> Begin task: report
   MM:     This is Galacticus: revision fb96dd63391b52e801e8c1059e6ece2ca433efd6   (branch: master; build time: Sun Oct 23 16:13:16 UTC 2022)
   MM:     Built with: :GSL_version[2.6]:FoX_version[4.1.2]:HDF5_version[1.8.9]:FCCOMPILER[gfortran]:PREPROCESSOR[cpp]:CCOMPILER[gcc]:CPPCOMPILER[g++]:FCFLAGS[-ffree-line-length-none -frecursive -DBUILDPATH='./work/build' -J./work/build/moduleBuild/ -I./work/build/ -fintrinsic-modules-path /home/abenson/Galacticus/Tools/finclude -fintrinsic-modules-path /home/abenson/Galacticus/Tools/include -fintrinsic-modules-path /home/abenson/Galacticus/Tools/include/gfortran -fintrinsic-modules-path /home/abenson/Galacticus/Tools/lib/gfortran/modules -L/home/abenson/Galacticus/Tools/lib -L/home/abenson/Galacticus/Tools/lib64 -pthread -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -fdump-core -O3 -ffinite-math-only -fno-math-errno -fopenmp -g -DPROCPS -DOFDUNAVAIL -DFFTW3AVAIL -DANNAVAIL -DQHULLAVAIL -DMATHEVALAVAIL]:FCFLAGS_NOOPT[-ffree-line-length-none -frecursive -DBUILDPATH='./work/build' -J./work/build/moduleBuild/ -I./work/build/ -fintrinsic-modules-path /home/abenson/Galacticus/Tools/finclude -fintrinsic-modules-path /home/abenson/Galacticus/Tools/include -fintrinsic-modules-path /home/abenson/Galacticus/Tools/include/gfortran -fintrinsic-modules-path /home/abenson/Galacticus/Tools/lib/gfortran/modules -L/home/abenson/Galacticus/Tools/lib -L/home/abenson/Galacticus/Tools/lib64 -pthread -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -fdump-core -g]:CFLAGS[-fopenmp -DBUILDPATH='./work/build' -I./source/ -I./work/build/ -I/home/abenson/Galacticus/Tools/include -g -DOFDLOCKS -DPROCPS -DOFDUNAVAIL]:CPPFLAGS[-fopenmp -DBUILDPATH='./work/build' -I./source/ -I./work/build/ -I/home/abenson/Galacticus/Tools/include -I/home/abenson/Galacticus/Tools/include/libqhullcpp -g -DOFDLOCKS -DPROCPS -DOFDUNAVAIL -DANNAVAIL -DQHULLAVAIL -DMATHEVALAVAIL]:FCCOMPILER_VERSION[...]
   MM: <- Done task: report

The ``fb96dd63391b52e801e8c1059e6ece2ca433efd6`` in the above is the Git revision hash from which this copy of Galacticus was built.

Release Steps
-------------

The binary executables (Linux and macOS), the run-time tools, and the Docker
image are all built and published automatically by the CI/CD pipeline
(``.github/workflows/cicd.yml``) on every push to ``master``, which keeps the
rolling ``bleeding-edge`` release up to date. Cutting a versioned release is
therefore lightweight:

#. Ensure all CI checks are green on the ``master`` commit you intend to tag.
   This includes that both the MPI and non-MPI builds compile cleanly without
   errors or warnings, and — importantly — that the CI/CD run for that commit has
   finished and uploaded the binaries and tools to ``bleeding-edge``.
#. Tag that commit with a semantic-version tag and push it::

      git tag vX.Y.Z
      git push origin vX.Y.Z

   The ``Publish`` workflow (``.github/workflows/publish.yml``) then builds and
   publishes the ``galacticus`` Python package to PyPI and creates the versioned
   GitHub release, mirroring the current ``bleeding-edge`` binary and tools
   assets onto it (see `Publishing the Python Package to PyPI`_ below). This is
   what makes ``pip install galacticus`` resolve to the new version.
#. Create the corresponding
   `datasets <https://github.com/galacticusorg/datasets>`_ release so that
   matching versions of the code and data are available.

.. note::

   Because the versioned release **mirrors the assets currently on**
   ``bleeding-edge``, only tag a commit once that commit's CI/CD run has finished
   and refreshed ``bleeding-edge`` — otherwise the versioned release would pick up
   stale binaries. The publish workflow guards against this by checking that the
   ``bleeding-edge`` tag points at the tagged commit, and fails fast if not.

.. note::

   The LaTeX/PDF manuals have been retired; documentation is published
   automatically on `ReadTheDocs <https://galacticus.readthedocs.io/>`_. The
   statically-linked Linux and macOS binaries, the run-time tool archives, the
   Docker image, and the archives of external build dependencies are all produced
   and uploaded automatically by CI (or maintained out of band), so none of them
   need to be assembled or relinked by hand as part of a release.

Publishing the Python Package to PyPI
-------------------------------------

``pip install galacticus`` installs the launcher described in
:doc:`../user-guide/installation/pip`, which downloads the pre-built executable,
datasets, and tools on first use. Publishing a new version to
`PyPI <https://pypi.org/project/galacticus/>`_ is automated by the
``.github/workflows/publish.yml`` workflow and is triggered simply by pushing a
semantic-version tag:

.. code-block:: bash

   git tag v1.2.3
   git push origin v1.2.3

On a version tag the workflow:

#. builds the source distribution and wheel — the version is derived from the
   git tag by `setuptools_scm <https://setuptools-scm.readthedocs.io/>`_, so the
   ``X.Y.Z`` you tag is exactly the version on PyPI;
#. publishes them to PyPI using `Trusted Publishing
   <https://docs.pypi.org/trusted-publishers/>`_ (OIDC), so no API token is
   stored in the repository; and
#. mirrors the current ``bleeding-edge`` binary and tools assets onto the
   versioned GitHub release, so the launcher fetches artefacts pinned to that
   version.

Because the launcher resolves the GitHub release tag from the installed package
version, the binary assets **must** be present on the ``vX.Y.Z`` release; the
workflow handles this automatically, but if you create a release by hand, attach
the binary and tools assets too.

One-time setup
~~~~~~~~~~~~~~

Before the first publish, register a PyPI (and TestPyPI) *Trusted Publisher* for
this repository, with the workflow filename ``publish.yml`` and environment
``pypi``. To dry-run the pipeline without touching the real index, run the
``Publish`` workflow manually (``workflow_dispatch``) with the ``testpypi`` input
set to ``true``; this builds and uploads to TestPyPI only.
