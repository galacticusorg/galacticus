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

Changes to the format or syntax of any of these files are considered to incompatible API changes. Note that additions to the formats that don't break backward compatibility (e.g. adding a new class with new parameters) *are* backwards compatible and so do not require a major version increment.

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

#. Both MPI and non-MPI builds must compile cleanly without errors or warnings;
#. Create releases on GitHub:

   #. `galacticus <https://github.com/galacticusorg/galacticus>`_, with assets:

      #. Statically-linked binary, ``galacticus.exe``.

      .. note::

         The LaTeX/PDF manuals have been retired. Documentation is now published automatically on `ReadTheDocs <https://galacticus.readthedocs.io/>`_, so there are no documentation PDFs to build or attach as release assets.

   #. `datasets <https://github.com/galacticusorg/datasets>`_ release in GitHub;
   #. `galacticusDockerBuildEnv <https://github.com/galacticusorg/galacticusDockerBuildEnv>`_

#. Publish the ``galacticus`` Python package to PyPI (see `Publishing the Python Package to PyPI`_ below). This is what makes ``pip install galacticus`` resolve to the new version.

#. Add archives of external dependencies to the GitHub release. Currently these are:

   #. `CAMB <https://github.com/cmbant/CAMB>`_

      #. plus `forutils <https://github.com/cmbant/forutils>`_

   #. `RecFast <https://www.astro.ubc.ca/people/scott/recfast.html>`_
   #. `FSPS <https://github.com/cconroy20/fsps>`_
   #. `Cloudy <https://gitlab.nublado.org/cloudy/cloudy/-/wikis/home>`_
   #. `FFTLog <http://jila.colorado.edu/~ajsh/FFTLog/fftlog.tgz>`_
   #. `FoX <https://github.com/andreww/fox>`_
   #. `ANN <http://www.cs.umd.edu/~mount/ANN>`_
   #. `libmatheval <https://github.com/galacticusorg/libmatheval>`_
   #. `FFTW3 <https://www.fftw.org/>`_
   #. `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_
   #. `gfortran <https://gfortran.meteodat.ch/>`_
   #. `GSL <https://www.gnu.org/software/gsl/>`_
   #. `Guile <https://www.gnu.org/software/guile/>`_

Docker images can be built locally following the instructions `here <https://github.com/galacticusorg/galacticus/wiki/Building-Docker-Images>`_.

Docker images can be exported as follows:

.. code-block:: bash

   sudo docker pull ghcr.io/galacticusorg/buildenv:latest
   sudo docker run --name buildenvV1.0.0 galacticusorg/buildenv
   sudo docker export --output="/home/abensonca/Scratch/galacticusorg_buildenv_v1.0.0.tar" buildenvV1.0.0
   sudo docker stop buildenvV1.0.0
   sudo bzip2 galacticusorg_buildenv_v1.0.0.tar

MacOS Binaries
~~~~~~~~~~~~~~

Compilation instructions for MacOS are given `here <https://github.com/galacticusorg/galacticus/wiki/Installation-from-source-on-MacOS>`_. To create statically-linked binaries on MacOS proceed as follows:

#. Compile Galacticus as usual.
#. Run ``./Galacticus.exe parameters/buildTools.xml`` to build all run-time tools (RecFast, FSPS, CAMB, Class, Cloudy).
#. Since MacOS does not really support static-linking we need to do some manual re-linking to making statically linked binaries. To make this easier a script ``./scripts/build/staticRelinker.py`` is provided. Static executables can be created using the following commands:

   .. code-block:: bash

      cd ~/galacticus
      ~/galacticus/scripts/build/staticRelinker.py gfortran-16 `cat ./work/build/Galacticus.d` ./work/build/Galacticus.parameters.o ./work/build/Galacticus.md5s.o -o Galacticus.exe -ffree-line-length-none -frecursive -DBUILDPATH=\'./work/build\' -J./work/build/moduleBuild/ -I./work/build/ -fintrinsic-modules-path /usr/local/include -L/usr/local/lib -L/opt/local/lib -fintrinsic-modules-path /usr/local/finclude -fintrinsic-modules-path /usr/local/include/gfortran -fintrinsic-modules-path /usr/local/include -fintrinsic-modules-path /usr/lib/gfortran/modules -fintrinsic-modules-path /usr/include/gfortran -fintrinsic-modules-path /usr/include -fintrinsic-modules-path /usr/finclude -fintrinsic-modules-path /usr/lib64/gfortran/modules -fintrinsic-modules-path /usr/lib64/openmpi/lib -pthread -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -fdump-core -O3 -ffinite-math-only -fno-math-errno -fopenmp -g -DOFDUNAVAIL -DFFTW3AVAIL -DANNAVAIL -DMATHEVALAVAIL `./scripts/build/libraryDependencies.py Galacticus.exe -ffree-line-length-none -frecursive -DBUILDPATH=\'./work/build\' -J./work/build/moduleBuild/ -I./work/build/ -fintrinsic-modules-path /usr/local/include -L/usr/local/lib -L/opt/local/lib -fintrinsic-modules-path /usr/local/finclude -fintrinsic-modules-path /usr/local/include/gfortran -fintrinsic-modules-path /usr/local/include -fintrinsic-modules-path /usr/lib/gfortran/modules -fintrinsic-modules-path /usr/include/gfortran -fintrinsic-modules-path /usr/include -fintrinsic-modules-path /usr/finclude -fintrinsic-modules-path /usr/lib64/gfortran/modules -fintrinsic-modules-path /usr/lib64/openmpi/lib -pthread -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -fdump-core -O3 -ffinite-math-only -fno-math-errno -fopenmp -g -DOFDUNAVAIL -DFFTW3AVAIL -DANNAVAIL -DMATHEVALAVAIL`

      cd ~/datasets/dynamic/RecFast
      ~/galacticus/scripts/build/staticRelinker.py gfortran-16 recfast.for -o recfast.exe -O3 -ffixed-form -ffixed-line-length-none

      cd ~/datasets/dynamic/CAMB-1.3.2/fortran
      ~/galacticus/scripts/build/staticRelinker.py gfortran-16 -cpp -Ofast -fopenmp -fintrinsic-modules-path /usr/local/include -L/usr/local/lib -L/opt/local/lib -JRelease -IRelease/ -I"Users/ajb/datasets/dynamic/CAMB-1.3.2/fortran/../forutils/Release/" Release/inidriver.o Release/libcamb.a  -L"/Users/ajb/datasets/dynamic/CAMB-1.3.2/fortran/../forutils/Release/" -lforutils -o camb

      cd ~/datasets/dynamic/class_public-3.0.2
      ~/galacticus/scripts/build/staticRelinker.py gcc-16 -O3 -fopenmp -g -fPIC -o class build/growTable.o build/dei_rkck.o build/sparse.o build/evolver_rkck.o build/evolver_ndf15.o build/arrays.o build/parser.o build/quadrature.o build/hyperspherical.o build/common.o build/trigonometric_integrals.o build/input.o build/background.o build/thermodynamics.o build/perturbations.o build/primordial.o build/fourier.o build/transfer.o build/harmonic.o build/lensing.o build/distortions.o build/wrap_recfast.o build/injection.o build/noninjection.o build/hyrectools.o build/helium.o build/hydrogen.o build/history.o build/wrap_hyrec.o build/energy_injection.o build/output.o build/class.o -lm

      cd ~/datasets/dynamic/fsps-3.2/src
      ~/galacticus/scripts/build/staticRelinker.py gfortran-16 -O3 -cpp -fPIC -mcmodel=medium  -o autosps.exe autosps.o sps_vars.o sps_utils.o compsp.o csp_gen.o galacticus_IMF.o ssp_gen.o getmags.o locate.o funcint.o sps_setup.o pz_convol.o get_tuniv.o intsfwght.o imf.o imf_weight.o add_dust.o getspec.o sbf.o add_bs.o mod_hb.o add_remnants.o getindx.o smoothspec.o mod_gb.o add_nebular.o write_isochrone.o sfhstat.o linterp.o tsum.o add_agb_dust.o linterparr.o ztinterp.o vacairconv.o igm_absorb.o get_lumdist.o attn_curve.o sfh_weight.o sfhlimit.o sfhinfo.o setup_tabular_sfh.o agn_dust.o

      cd ~/datasets/dynamic/c23.01/source
      ~/galacticus/scripts/build/staticRelinker.py g++-16 -O3 -ftrapping-math -fnop-math-errno -ftree-vectorize -Wall -g -o cloudy.exe maincl.o -L. -lcloudy

#. Package the executables (and data):

   .. code-block:: bash

      zip -r galacticusExecutablesMacOS.zip galacticus/Galacticus.exe datasets/dynamic/fsps-3.2 datasets/dynamic/class_public-3.0.2 datasets/dynamic/CAMB-1.3.2 datasets/dynamic/RecFast datasets/dynamic/c17.02

#. Add ``galacticusExecutablesMacOS.zip`` to the GitHub release.

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
