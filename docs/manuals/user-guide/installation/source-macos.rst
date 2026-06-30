Installing from Source (macOS)
==============================

.. note::

   These installation instructions are still in beta-testing. Please report failures on the `discussion forums <https://github.com/galacticusorg/galacticus/discussions>`_.

A beta-release of an installation script (that automates the process of installing Galacticus from source) for MacOS is available in the `installationScripts <https://github.com/galacticusorg/installationscripts>`_ repo. You can find instructions for using this script `here <https://github.com/galacticusorg/installationscripts/wiki>`_.

.. note::

   **GCC 16** is the minimum supported compiler version. These instructions install GCC and the build dependencies through `Homebrew <https://brew.sh/>`_, so the compilers are named ``gcc-16``, ``g++-16``, and ``gfortran-16`` throughout. If you install GCC through a different package manager, substitute the corresponding executable names.

Install Xcode Command Line Tools
--------------------------------

If you don't already have the Xcode command line tools installed, install them now:

.. code-block:: bash

   if [[ ! $(xcode-select -p) ]]; then
       xcode-select --install
   fi
   export PATH=$PATH:$(brew --prefix)/bin:/usr/local/bin

Install Homebrew
----------------

If you don't already have `Homebrew <https://brew.sh/>`_ installed, install it now:

.. code-block:: bash

   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

Also record your macOS major version. It is used for a linker workaround in the HDF5 and Galacticus build steps below:

.. code-block:: bash

   os_ver=$(sw_vers -productVersion)
   IFS='.' read -r -a ver <<< "$os_ver"

Install GCC
-----------

GCC 16 (the minimum supported version) is required to build Galacticus. Install the GCC compilers through Homebrew:

.. code-block:: bash

   brew install gcc

Homebrew names the executables with their major version — ``gcc-16``, ``g++-16``, and ``gfortran-16`` — and these names are used throughout the remaining steps.

.. note::

   ``brew install gcc`` installs the current GCC release. If Homebrew installs a major version newer than 16, adjust the ``-16`` suffixes in the commands below to match (and ensure the version is at least 16).

Install Guile
-------------

.. code-block:: bash

   brew install guile

Install GSL
-----------

.. code-block:: bash

   brew install gsl

Install ``libmatheval``
-----------------------

.. code-block:: bash

   curl -L https://github.com/galacticusorg/libmatheval/releases/download/latest/libmatheval-1.1.13.tar.gz --output libmatheval-1.1.13.tar.gz
   tar xvfz libmatheval-1.1.13.tar.gz
   cd libmatheval-1.1.13
   sed -E -i~ s/"#undef HAVE_SCM_T_BITS"/"#define HAVE_SCM_T_BITS 1"/ config.h.in
   CC=gcc-16 ./configure --prefix=/usr/local
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
   make -j CC=gcc-16 CXX=g++-16
   sudo make install
   cd ..
   rm -rf qhull-2020-src-8.0.2.tgz qhull-2020.2

Install HDF5
------------

HDF5 2.x is built with CMake (the Autotools ``./configure`` build system was removed in HDF5 2.0). Install CMake via Homebrew if you do not already have it:

.. code-block:: bash

   brew install cmake

.. code-block:: bash

   curl -L https://github.com/HDFGroup/hdf5/releases/download/2.1.0/hdf5-2.1.0.tar.gz --output hdf5-2.1.0.tar.gz
   tar -vxzf hdf5-2.1.0.tar.gz
   cd hdf5-2.1.0
   cmakeExtraFlags=()
   if [[ "${ver}" -eq 13 ]]; then
      # For MacOS 13 force use of the classic linker as the new linker does not support the '-commons' option - see https://trac.macports.org/ticket/68194#comment:15
      cmakeExtraFlags=(-DCMAKE_EXE_LINKER_FLAGS=-Wl,-ld_classic -DCMAKE_SHARED_LINKER_FLAGS=-Wl,-ld_classic -DCMAKE_MODULE_LINKER_FLAGS=-Wl,-ld_classic)
   fi
   cmake -S . -B build \
       -DCMAKE_INSTALL_PREFIX=/usr/local \
       -DCMAKE_BUILD_TYPE=Release \
       -DCMAKE_C_COMPILER=gcc-16 \
       -DCMAKE_CXX_COMPILER=g++-16 \
       -DCMAKE_Fortran_COMPILER=gfortran-16 \
       -DHDF5_BUILD_FORTRAN=ON \
       -DHDF5_BUILD_HL_LIB=ON \
       -DHDF5_ENABLE_DEPRECATED_SYMBOLS=OFF \
       "${cmakeExtraFlags[@]}"
   cmake --build build -j3
   sudo cmake --install build
   cd ..
   rm -rf hdf5-2.1.0 hdf5-2.1.0.tar.gz

.. note::

   Galacticus builds against either HDF5 1.14 or HDF5 2.x without modification, so if you already have a working HDF5 1.14 installation you may continue to use it. Support for HDF5 1.x will, however, eventually be deprecated, so HDF5 2.x is recommended for new installations.

Install FoX
-----------

.. code-block:: bash

   curl -L https://github.com/galacticusorg/fox/archive/refs/tags/v4.1.4.tar.gz --output FoX-4.1.4.tar.gz
   tar xvfz FoX-4.1.4.tar.gz
   cd fox-4.1.4
   FC=gfortran-16 ./configure --prefix=/usr/local
   make -j
   sudo make install
   cd ..
   rm -rf fox-4.1.4 FoX-4.1.4.tar.gz

Install FFTW3
-------------

.. code-block:: bash

   curl -L ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz --output fftw-3.3.4.tar.gz
   tar xvfz fftw-3.3.4.tar.gz
   cd fftw-3.3.4
   F77=gfortran-16 CC=gcc-16 ./configure --prefix=/usr/local
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
   sed -E -i~ s/"C\+\+ = g\+\+"/"C\+\+ = g\+\+\-16"/ Make-config
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
   export FCCOMPILER=gfortran-16
   export CCOMPILER=gcc-16
   export CPPCOMPILER=g++-16
   export LIBRARY_PATH=/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk/usr/lib
   export GALACTICUS_FCFLAGS="-fintrinsic-modules-path /usr/local/include -fintrinsic-modules-path /usr/local/finclude -L/usr/local/lib -L$(brew --prefix)/lib -I/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk/usr/include"
   if [[ "${ver}" -eq 13 ]]; then
       export GALACTICUS_FCFLAGS="$GALACTICUS_FCFLAGS -Wl,-ld_classic"
   fi
   export GALACTICUS_CFLAGS="-I/usr/local/include -I$(brew --prefix)/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk/usr/include"
   export GALACTICUS_CPPFLAGS="-I/usr/local/include -I$(brew --prefix)/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk/usr/include"
   make -j Galacticus.exe

.. note::

   The macOS SDK version in these paths (``MacOSX15.4.sdk``) must be changed to match the SDK actually installed on your system. You can list the available SDKs under ``/Library/Developer/CommandLineTools/SDKs/``.

Finally
-------

Follow the instructions to download run-time `datasets <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/installation/binary.html>`_ used by Galacticus. You may want to set the above environment variable exports in your ``~/.zshenv`` so they are set automatically.
