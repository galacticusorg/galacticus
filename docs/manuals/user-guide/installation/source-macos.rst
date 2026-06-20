Installing from Source (macOS)
==============================

.. note::

   These installation instructions are still in beta-testing. Please report failures on the `discussion forums <https://github.com/galacticusorg/galacticus/discussions>`_.

A beta-release of an installation script (that automates the process of installing Galacticus from source) for MacOS is available in the `installationScripts <https://github.com/galacticusorg/installationscripts>`_ repo. You can find instructions for using this script `here <https://github.com/galacticusorg/installationscripts/wiki>`_.

.. note::

   **GCC 16** is the minimum supported compiler version. These instructions install GCC and the build dependencies through `MacPorts <https://www.macports.org/>`_, so the compilers are named ``gcc-mp-16``, ``g++-mp-16``, and ``gfortran-mp-16`` throughout. If you install GCC through a different package manager (e.g. Homebrew), substitute the corresponding executable names (``gcc-16``, ``g++-16``, ``gfortran-16``).

Install Xcode Command Line Tools
--------------------------------

If you don't already have the Xcode command line tools installed, install them now:

.. code-block:: bash

   if [[ ! $(xcode-select -p) ]]; then
       xcode-select --install
   fi
   export PATH=$PATH:/opt/local/bin:/usr/local/bin

Install MacPorts
----------------

Download the appropriate MacPorts for your system. The following will detect your OS version and download and install the correct version.

.. code-block:: bash

   os_ver=$(sw_vers -productVersion)
   IFS='.' read -r -a ver <<< "$os_ver"
   if [[ "${ver}" -eq 11 ]]; then
       macportsversion=2.7.1
       macportsbase=2.7.1-11-BigSur
   elif [[ "${ver}" -eq 12 ]]; then
       macportsversion=2.9.1
       macportsbase=2.9.1-12-Monterey
   elif [[ "${ver}" -eq 13 ]]; then
       macportsversion=2.9.1
       macportsbase=2.9.1-13-Ventura
   elif [[ "${ver}" -eq 14 ]]; then
       macportsversion=2.9.1
       macportsbase=2.9.1-14-Sonoma
   else
       echo Unknown MacOS version: ${os_ver}
       exit 1
   fi
   curl -L https://github.com/macports/macports-base/releases/download/v${macportsversion}/MacPorts-${macportsbase}.pkg --output MacPorts-${macportsbase}.pkg
   sudo installer -pkg ./MacPorts-${macportsbase}.pkg -target /
   rm ./MacPorts-${macportsbase}.pkg

Install GCC
-----------

GCC 16 is required to build Galacticus. Install it through MacPorts:

.. code-block:: bash

   sudo port install gcc16

This provides the ``gcc-mp-16``, ``g++-mp-16``, and ``gfortran-mp-16`` executables used in the remaining steps.

.. note::

   GCC 16 is very recent and may not yet be packaged for your system. If ``port install gcc16`` cannot find it, install the most recent GCC that is available (adjusting the ``-16`` suffixes in the commands below to match), or build GCC from source as described in the `Linux instructions <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/installation/source-linux.html>`_.

Install Guile
-------------

.. code-block:: bash

   sudo port install guile-3.0

Install GSL
-----------

.. code-block:: bash

   sudo port install gsl

Install ``libmatheval``
-----------------------

.. code-block:: bash

   curl -L https://github.com/galacticusorg/libmatheval/releases/download/latest/libmatheval-1.1.13.tar.gz --output libmatheval-1.1.13.tar.gz
   tar xvfz libmatheval-1.1.13.tar.gz
   cd libmatheval-1.1.13
   sed -E -i~ s/"#undef HAVE_SCM_T_BITS"/"#define HAVE_SCM_T_BITS 1"/ config.h.in
   CC=gcc-mp-16 ./configure --prefix=/usr/local
   make -j
   sudo make install
   cd ..
   rm -rf libmatheval-1.1.13.tar.gz libmatheval-1.1.13

Install QHull
~~~~~~~~~~~~~

.. code-block:: bash

   curl -L http://www.qhull.org/download/qhull-2020-src-8.0.2.tgz --output qhull-2020-src-8.0.2.tgz
   tar xvfz qhull-2020-src-8.0.2.tgz
   cd qhull-2020.2
   make -j CC=gcc-mp-16 CXX=g++-mp-16
   sudo make install
   cd ..
   rm -rf qhull-2020-src-8.0.2.tgz qhull-2020.2

Install HDF5
------------

.. code-block:: bash

   curl -L https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_5/downloads/hdf5-1.14.5.tar.gz --output hdf5-1.14.5.tar.gz
   tar -vxzf hdf5-1.14.5.tar.gz
   cd hdf5-1.14.5
   if [[ "${ver}" -eq 13 ]]; then
      # For MacOS 13 force use of the classic linker as the new linker does not support the '-commons' option - see https://trac.macports.org/ticket/68194#comment:15
      CC=gcc-mp-16 CXX=g++-mp-16 FC=gfortran-mp-16 LDFLAGS=-Wl,-ld_classic ./configure --prefix=/usr/local --enable-fortran --enable-build-mode=production
   else
      CC=gcc-mp-16 CXX=g++-mp-16 FC=gfortran-mp-16                         ./configure --prefix=/usr/local --enable-fortran --enable-build-mode=production
   fi
   make -j3
   sudo make install
   cd ..
   rm -rf hdf5-1.14.5 hdf5-1.14.5.tar.gz

Install FoX
-----------

.. code-block:: bash

   curl -L https://github.com/galacticusorg/fox/archive/refs/tags/v4.1.3.tar.gz --output FoX-4.1.3.tar.gz
   tar xvfz FoX-4.1.3.tar.gz
   cd fox-4.1.3
   FC=gfortran-mp-16 ./configure --prefix=/usr/local
   make -j
   sudo make install
   cd ..
   rm -rf fox-4.1.3 FoX-4.1.3.tar.gz

Install FFTW3
-------------

.. code-block:: bash

   curl -L ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz --output fftw-3.3.4.tar.gz
   tar xvfz fftw-3.3.4.tar.gz
   cd fftw-3.3.4
   F77=gfortran-mp-16 CC=gcc-mp-16 ./configure --prefix=/usr/local
   make -j
   sudo make install
   cd ..
   rm -rf fftw-3.3.4 fftw-3.3.4.tar.gz

Install ANN
-----------

.. code-block:: bash

   curl -L http://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz --output ann_1.1.2.tar.gz
   tar xvfz ann_1.1.2.tar.gz
   cd ann_1.1.2
   sed -E -i~ s/"C\+\+ = g\+\+"/"C\+\+ = g\+\+\-mp\-16"/ Make-config
   make macosx-g++
   sudo cp bin/* /usr/local/bin/.
   sudo cp lib/* /usr/local/lib/.
   sudo cp -R include/* /usr/local/include/.

Install Python dependencies
---------------------------

Galacticus uses Python during compilation to preprocess the source code, and ships various Python scripts under ``scripts/`` for analysis and pipeline tasks. Python ≥ 3.9 is required.

The Python package ``galacticus`` (declared in the top-level ``pyproject.toml``) installs the modules under ``python/`` onto the Python import path together with their third-party dependencies — ``numpy``, ``scipy``, ``h5py``, ``lxml``, ``matplotlib``, ``requests``, ``PyYAML``, ``GitPython``, ``PyPDF2``, ``termcolor`` — via a single editable install run from the cloned repository's root:

.. code-block:: bash

   pip3 install -e .

This is sufficient for compiling Galacticus and running its supporting scripts. Two optional extras are also available:

* ``pip3 install -e '.[emulation]'`` — adds dependencies needed by the emulator pipelines under ``scripts/emulation/``.
* ``pip3 install -e '.[test]'`` — adds ``pytest`` for running the Python unit-test suite.

Install Galacticus
------------------

.. code-block:: bash

   git clone https://github.com/galacticusorg/galacticus.git
   cd galacticus
   export GALACTICUS_EXEC_PATH=`pwd`
   export FCCOMPILER=gfortran-mp-16
   export CCOMPILER=gcc-mp-16
   export CPPCOMPILER=g++-mp-16
   export LIBRARY_PATH=/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk/usr/lib
   export GALACTICUS_FCFLAGS="-fintrinsic-modules-path /usr/local/include -fintrinsic-modules-path /usr/local/finclude -L/usr/local/lib -L/opt/local/lib -I/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk/usr/include"
   if [[ "${ver}" -eq 13 ]]; then
       export GALACTICUS_FCFLAGS="$GALACTICUS_FCFLAGS -Wl,-ld_classic"
   fi
   export GALACTICUS_CFLAGS="-I/usr/local/include -I/opt/local/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk/usr/include"
   export GALACTICUS_CPPFLAGS="-I/usr/local/include -I/opt/local/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk/usr/include"
   make -j Galacticus.exe

.. note::

   The macOS SDK version in these paths (``MacOSX15.4.sdk``) must be changed to match the SDK actually installed on your system. You can list the available SDKs under ``/Library/Developer/CommandLineTools/SDKs/``.

Finally
-------

Follow the instructions to download run-time `datasets <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/installation/binary.html>`_ used by Galacticus. You may want to set the above environment variable exports in your ``~/.zshenv`` so they are set automatically.
