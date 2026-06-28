Installation
============

Galacticus is designed to run on Linux, with experimental support for macOS.
There are several ways to install it, listed here from easiest to most involved:

* **pip** — the quickest way to get running for end users. ``pip install
  galacticus`` installs a launcher that downloads the right pre-built binary,
  datasets, and tools on first use and sets the environment up for you; see
  :doc:`pip`.
* **Pre-compiled binary** — download and configure the binary yourself; no
  compilation required. Available for :doc:`Linux <binary>` and
  :doc:`macOS <binary-macos>`.
* **Container** — a ready-to-use Docker image with Galacticus already installed;
  see :doc:`container`.
* **From source** — needed if you want to modify or extend the code. Step-by-step
  instructions are provided for :doc:`Linux <source-linux>` and
  :doc:`macOS <source-macos>`. For an automated source build you can also try the
  `installation script <https://github.com/galacticusorg/installationscripts/wiki>`_.

Whichever route you choose, Galacticus needs a set of run-time **datasets**
(downloaded separately) and a small number of environment variables to be set.
These are covered in the instructions for each route and in :doc:`../running`.

Prerequisites
-------------

The pre-compiled binary and container routes bundle everything that is needed, so
the following are required only when **building from source**:

Compiler
   A recent GCC is required for the Fortran 2003/2008 features that Galacticus
   uses: ``gfortran`` ≥ 16. Earlier versions will not compile Galacticus; later
   versions may work but are less extensively tested.

Required libraries
   `GSL <http://www.gnu.org/software/gsl/>`_, ``zlib``,
   `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_,
   `FoX <https://github.com/galacticusorg/fox>`_, and BLAS. HDF5 and FoX must be
   built with the same GCC version used to build Galacticus, to ensure module
   compatibility.

Optional libraries
   `FFTW3 <http://www.fftw.org/>`_ (survey-geometry calculations), ANN (analysis
   of N-body point data, and using MCMC posteriors as priors), and
   ``libmatheval`` (derived-parameter expressions in posterior simulations).
   These are detected at build time; if they are not found, the code that uses
   them is simply not compiled.

Python
   Python ≥ 3.9, used by the build system to preprocess the source code. Install
   the bundled ``galacticus`` package and its dependencies with a single editable
   install from the repository root::

      pip install -e .

   See :doc:`source-linux` for details and the optional ``[emulation]`` and
   ``[test]`` extras.

.. toctree::
   :maxdepth: 1
   :caption: Installation routes

   pip
   binary
   binary-macos
   container
   source-linux
   source-macos
   advanced-options
