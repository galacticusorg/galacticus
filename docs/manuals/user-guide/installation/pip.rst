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

On first use you will see the launcher fetch the executable, datasets, and tools,
with a progress bar for each download and for unpacking each archive; the
downloads run concurrently, and each large one is split across several
connections, so the whole set arrives at roughly the speed of your link rather
than of one stream. Subsequent runs reuse the cached copies. ``galacticus run`` validates the parameter file
before dispatching it; pass ``--no-validate`` to skip that, and any other
arguments (e.g. ``--dry-run``) are passed straight through to the executable.
``galacticus <file>`` is shorthand for ``galacticus run <file>``. If you do not
need the run-time tools, ``galacticus install --no-tools`` makes that first-use
download substantially smaller and faster — see :ref:`pip-no-tools`.

The bundled example parameter files (such as ``parameters/quickTest.xml``)
resolve against the install, so the command above works from any directory --
you do not need to ``cd`` into the install tree. A relative path that exists
in your current directory always takes precedence, so your own parameter files
are found first.

The model writes its output (by default ``galacticus.hdf5``) to the current
directory, exactly as the executable does when run directly. See
:doc:`../running` for what to do with the output.

Commands
--------

``galacticus install``
   Download (or complete) the install without running a model. Pass
   ``--no-tools`` to skip the pre-built tools archive — see
   :ref:`pip-no-tools` below.

``galacticus update``
   Re-download the install for the current package version.

``galacticus validate <file> [change files...]``
   Validate a parameter file without running it. Validation is performed on the
   *resolved* tree (XInclude, any change files, and ``active=`` conditionals are
   applied first), so it checks the structure Galacticus will actually build.

``galacticus resolve <file> [change files...] -o <output>``
   Apply the file-level transformations Galacticus performs when reading a
   parameter file — XInclude expansion, change-file application, and ``active=``
   conditional evaluation/pruning — and write a single, clean, self-contained
   parameter file to ``<output>``. Math expressions (``=[...]``) and ``id``/
   ``idRef`` anchors are left intact for Galacticus to handle at run time. Pass
   ``--no-conditionals`` to leave conditionals in place, or ``--validate`` to
   validate the result. This needs no executable or download.

   This is the recommended way to use the launcher with **MPI**: resolve once,
   then launch the (unchanged) executable under ``mpirun`` on the resolved file —

   .. code-block:: bash

      galacticus resolve model.xml changes.xml -o resolved.xml
      mpirun -n 16 Galacticus.exe resolved.xml

   Do **not** run ``mpirun galacticus run …``: that would resolve and launch once
   per rank. (For a single-process run, ``galacticus run --resolve <file>``
   resolves to a temporary file and runs that.)

``galacticus clean``
   Purge the regenerable data cache so it cannot grow without bound. Use
   ``--older-than N`` to remove only files older than ``N`` days, ``--all`` to
   empty the cache, ``--dry-run`` to report how much would be freed without
   deleting, and ``--prune-installs`` to also remove superseded per-version
   installs. ``clean`` never removes the executable, datasets, or tools.

``galacticus info``
   Show the resolved install, the environment variables it sets, and the current
   cache size.

.. _pip-no-tools:

Installing without the pre-built tools
--------------------------------------

The ``tools`` archive (CAMB, CLASS, Cloudy, FSPS, …) is by far the largest thing
the launcher downloads, and many models never touch it — a model reading a
tabulated transfer function or a pre-computed stellar population, for example,
needs none of it. Skip it for a faster, smaller install::

   galacticus install --no-tools

The choice is remembered, so a later ``galacticus run`` will not silently
download the archive you just declined. Add the tools at any time with::

   galacticus install --tools

A model that *does* need a tool will fail without it: a binary-only install has
no compilers, so Galacticus cannot build the missing tool itself. ``galacticus
info`` reports whether tools are installed or were skipped.

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

For such a build, catalog-aware validation needs the parameter catalog. A managed
install downloads a pre-built one from the release (falling back to generating it
if the release publishes none), but a source build has to make its own: generate
it once with ``make parameters-catalog`` (it is written to the build tree where
the launcher looks for it). Without it, ``galacticus validate`` falls back to the
executable's ``--dry-run``.

.. note::

   The parameter catalog is not the same thing as ``schema/parameters.xsd``,
   which *is* committed and so arrives with the source archive. That is a
   deliberately lax schema for editor assistance; validation of a parameter
   file against the parameters each selected implementation actually accepts
   uses the catalog, which is derived from the source and published as a
   release asset.

.. note::

   The launcher fetches assets from the GitHub release matching the installed
   package version (development installs track the rolling ``bleeding-edge``
   release). The run-time datasets are installed from a snapshot published on
   the release and pinned to the ``datasets`` commit it was taken from, so a
   given release always installs the same data — the commit is recorded on the
   release as ``datasets.ref``. ``GALACTICUS_RELEASE_TAG`` and
   ``GALACTICUS_DATASETS_REF`` override the release tag and datasets ref
   respectively; asking for a specific ref fetches that commit from the
   ``datasets`` repository instead of using the release snapshot, which is
   slower but is how you track ``master`` or test an unreleased data change.

.. note::

   ``GALACTICUS_DOWNLOAD_CONNECTIONS`` sets how many byte-range requests a
   single large asset is split across (default 4). Set it to ``1`` if a
   network or proxy handles ranged requests badly.
