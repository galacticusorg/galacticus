Installing with ``pip``
=======================

The quickest way to get a working Galacticus for running models is to install
the ``galacticus`` package from `PyPI <https://pypi.org/project/galacticus/>`_::

   pip install galacticus

This installs a small launcher (the ``galacticus`` command). The first time you
run a model, the launcher downloads the pre-built executable, the run-time
``datasets``, and the pre-built ``tools`` for your platform, and sets the
required environment variables for you — so there is nothing else to configure.

Pre-built binaries are available for Linux (x86-64), macOS (Intel x86-64), and
macOS (Apple Silicon). On other platforms — including native Windows — there is
no pre-built binary; build :doc:`from source <source-linux>` instead (on Windows,
use WSL 2 and the Linux build).

.. note::

   The pre-built macOS binaries are compiled on a recent macOS and will only run
   on that version or newer. The launcher checks this before running: if your
   macOS is too old it stops with a clear message (rather than a cryptic
   ``dyld`` error) telling you to upgrade macOS or build from source.

Running a model
---------------

.. code-block:: bash

   galacticus run parameters/quickTest.xml

On first use you will see the launcher fetch the executable, datasets, and tools;
subsequent runs reuse the cached copies. ``galacticus run`` validates the
parameter file before dispatching it; pass ``--no-validate`` to skip that, and
any other arguments (e.g. ``--dry-run``) are passed straight through to the
executable. ``galacticus <file>`` is shorthand for ``galacticus run <file>``.

The model writes its output (by default ``galacticus.hdf5``) to the current
directory, exactly as the executable does when run directly. See
:doc:`../running` for what to do with the output.

Commands
--------

``galacticus install``
   Download (or complete) the install without running a model.

``galacticus update``
   Re-download the install for the current package version.

``galacticus validate <file>``
   Validate a parameter file without running it.

``galacticus clean``
   Purge the regenerable data cache so it cannot grow without bound. Use
   ``--older-than N`` to remove only files older than ``N`` days, ``--all`` to
   empty the cache, ``--dry-run`` to report how much would be freed without
   deleting, and ``--prune-installs`` to also remove superseded per-version
   installs. ``clean`` never removes the executable, datasets, or tools.

``galacticus info``
   Show the resolved install, the environment variables it sets, and the current
   cache size.

Where things are stored
-----------------------

The launcher keeps two separate locations (resolved via
`platformdirs <https://pypi.org/project/platformdirs/>`_):

* a **durable data directory** (``user_data_dir``) for the executable, datasets,
  and pre-built tools — managed by ``install``/``update``; and
* a **cache directory** (``user_cache_dir``) for regenerable data such as
  transfer functions and stellar-population spectra — safe to delete, and the
  only thing ``galacticus clean`` ever touches.

Tools are deliberately kept in the durable directory (via
``GALACTICUS_TOOLS_PATH``): a binary-only install has no compilers, so losing the
pre-built tools to a cache purge would leave Galacticus unable to rebuild them.

Using an existing build
-----------------------

The launcher also works as a front-end to a Galacticus you built yourself. If
``GALACTICUS_EXEC_PATH`` and ``GALACTICUS_DATA_PATH`` are already set and the
executable is present, or if ``GALACTICUS_HOME`` points at a build/clone tree
containing ``Galacticus.exe``, the launcher uses that install and skips all
downloads. Run ``galacticus info`` to see which install is in effect.

.. note::

   The launcher fetches assets from the GitHub release matching the installed
   package version (development installs track the rolling ``bleeding-edge``
   release). For a versioned release the run-time datasets are pinned to a
   specific ``datasets`` commit, recorded on the release, so a given package
   version always installs the same data. ``GALACTICUS_RELEASE_TAG`` and
   ``GALACTICUS_DATASETS_REF`` override the release tag and datasets ref
   respectively.
