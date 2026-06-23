Installing from Source (Linux)
==============================

Notes
-----

These instructions are intended for installation of Galacticus from source on a Linux machine. For instructions for MacOS see `MacOS <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/installation/source-macos.html>`_.

Before you begin, you might want to consider using the `installation script <https://github.com/galacticusorg/installationscripts/wiki>`_ instead which attempts to automate the process of installing all dependencies and compiling Galacticus for you.

If you have root access (or if your sysadmin is friendly) some of Galacticus dependencies (those marked as "may be available via your distro's package manager") might be available through your distro's package manager. Where possible, installation via the package manager is a good option. Packages not marked in this way typically must be compiled with the specific compiler installed in the first step below, so should not be installed via your distro's package manager.

Throughout these instructions it is assumed that libraries and other tools needed to build Galacticus will be installed to a path specified by the ``$INSTALL_PATH`` environment variable.

Requirements
------------

Note that, if you choose to install any of the following libraries and tools as root, by using the ``sudo`` command for example, be aware that this can change your environment resulting in installation to locations you do not intend. To avoid this you can force certain environment variables to be set, e.g.:

.. code-block:: bash

   sudo PATH="$PATH" HOME="$HOME" LD_LIBRARY_PATH="$LD_LIBRARY_PATH" make install

will run the command ``make install`` as root, but ensuring that the relevant environment variables (``HOME``, ``PATH``, and ``LD_LIBRARY_PATH``) are preserved.

Compilers
~~~~~~~~~

Galacticus is usually built using the GCC compiler. Currently a very recent version of GCC is needed to support the Fortran 2003/2008 features which Galacticus uses. It is recommended to use `GCC 16 <https://gcc.gnu.org/gcc-16/criteria.html>`_ - earlier versions are not compatible with Galacticus and will not compile it. Later versions may also work, but have not been extensively tested. Currently, GCC 16 can be compiled from source, as follows:

.. code-block:: bash

   wget https://ftp.gnu.org/gnu/gcc/gcc-16.1.0/gcc-16.1.0.tar.gz
   tar xvf gcc-16.1.0.tar.gz
   cd gcc-16.1.0
   ./contrib/download_prerequisites
   cd ..
   mkdir gcc-build
   cd gcc-build
   ../gcc/configure --prefix=$INSTALL_PATH --enable-languages=c,c++,fortran --disable-multilib
   make -j8
   make -j8 install

The ``-j8`` in the last two commands will use 8 cores to build GCC - if you have fewer or more cores you can change this number as appropriate - the more cores you use the faster the build will happen.

Libraries
~~~~~~~~~

There are several required (and some optional) libraries needed for Galacticus. The lists below give the preferred version for each library. The HDF5, and FoX libraries should be compiled using the same version of GCC as used to compile Galacticus to ensure module file compatibility - instructions are given for these libraries, and ``$INSTALL_PATH`` should be added to the start of your ``PATH`` and ``LD_LIBRARY_PATH`` environment variables to ensure that the correct GCC is used, for example:

.. code-block:: bash

   export PATH=$INSTALL_PATH/bin:$PATH
   export LD_LIBRARY_PATH=$INSTALL_PATH/lib64:$INSTALL_PATH/lib:$LD_LIBRARY_PATH

Required
^^^^^^^^

* GSL (may be available via your distro's package manager) `v2.6 <http://mirror.rit.edu/gnu/gsl/gsl-2.6.tar.gz>`_

.. code-block:: bash

   wget "http://mirror.rit.edu/gnu/gsl/gsl-2.6.tar.gz"
   tar xvfz gsl-2.6.tar.gz
   cd  gsl-2.6
   ./configure --prefix=$INSTALL_PATH
   make
   make check
   make install

* zlib `v1.3.1 <https://zlib.net/zlib-1.3.2.tar.gz>`_

.. code-block:: bash

   wget https://zlib.net/zlib-1.3.2.tar.gz
   tar -vxzf zlib-1.3.2.tar.gz
   cd zlib-1.3.2
   ./configure --prefix=$INSTALL_PATH
   make
   make install

* HDF5 `v1.14.5 <https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_5/downloads/hdf5-1.14.5.tar.gz>`_

.. code-block:: bash

   wget https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_5/downloads/hdf5-1.14.5.tar.gz
   tar -vxzf hdf5-1.14.5.tar.gz
   cd hdf5-1.14.5
   F9X=gfortran ./configure --prefix=$INSTALL_PATH --enable-fortran --enable-build-mode=production
   make
   make install

* FoX `v4.1.3 <https://github.com/galacticusorg/fox/archive/refs/tags/v4.1.3.tar.gz>`_

.. code-block:: bash

   wget https://github.com/galacticusorg/fox/archive/refs/tags/v4.1.3.tar.gz -O fox-4.1.3.tar.gz
   tar -vxzf fox-4.1.3.tar.gz
   cd fox-4.1.3
   FC=gfortran ./configure --prefix=$INSTALL_PATH
   make
   make install

* BLAS (may be available via your distro's package manager)

.. code-block:: bash

   wget http://www.netlib.org/blas/blas-3.11.0.tgz
   tar -vxzf blas-3.11.0.tgz
   cd BLAS-3.11.0
   make
   cp -f blas_LINUX.a $INSTALL_PATH/lib/libblas.a

Optional
^^^^^^^^

The following libraries are optional. Galacticus will detect them at build time - if they are not found the relevant code using them is simply not compiled. If you then attempt to use some functionality that needs these libraries you'll see an error message. In general they are required only for unusual applications of Galacticus, so you can probably ignore them

* FFTW `v3.3.4 <http://www.fftw.org/fftw-3.3.4.tar.gz>`_

  * Used only for calculations related to survey geometries - unless you know for sure that you want to do such calculations with Galacticus you do not need to install FFTW

* ann (may be available via your distro's package manager) `v1.1.2 <http://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz>`_

  * Used only for analysis of N-body point data and when using MCMC chain posteriors as priors on parameters - if you don't plan to do either of these you do not need to install ann.

* matheval (may be available via your distro's package manager) `v1.1.13 <https://github.com/galacticusorg/libmatheval/releases/download/latest/libmatheval-1.1.13.tar.gz>`_

  * Used only for evaluating of derived parameter expressions when performing simulations of the model posterior - if you don't plan to do this you do not need to install matheval

Python dependencies
~~~~~~~~~~~~~~~~~~~~

Galacticus uses Python during compilation to preprocess the source code, and ships a number of Python scripts under ``scripts/`` for analysis and pipeline tasks.  Python ≥ 3.9 is required.

The Python package ``galacticus`` (declared in the top-level ``pyproject.toml``) installs the modules under ``python/`` onto the Python import path together with their third-party dependencies — ``numpy``, ``scipy``, ``h5py``, ``lxml``, ``matplotlib``, ``requests``, ``PyYAML``, ``GitPython``, ``PyPDF2``, ``termcolor`` — via a single editable install run from the repository root:

.. code-block:: bash

   pip install -e .

This is sufficient for compiling Galacticus and running its supporting scripts.

Two optional extras are also available:

* ``pip install -e '.[emulation]'`` — adds the heavier dependencies needed by the emulator pipelines under ``scripts/emulation/`` (TensorFlow, lenstronomy, colossus, samana, tf-keras, tqdm, mcfit).
* ``pip install -e '.[test]'`` — adds ``pytest`` for running the Python unit-test suite.

If ``pip`` is not available system-wide, exporting ``PYTHONPATH`` before invoking ``make`` is also sufficient to put the Galacticus modules on the import path (the third-party dependencies must then be available in the active Python environment by some other means):

.. code-block:: bash

   export PYTHONPATH=$(pwd)/python:$PYTHONPATH

Compiling Galacticus
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone git@github.com:galacticusorg/galacticus.git
   cd galacticus
   export GALACTICUS_EXEC_PATH=`pwd`
   make Galacticus.exe

To build with MPI parallelism replace the ``make`` command in the above with:

.. code-block:: bash

   make GALACTICUS_BUILD_OPTION=MPI Galacticus.exe

Galacticus can take a long time to compile, so you may want to do a parallel make. For example:

.. code-block:: bash

   make -j8 Galacticus.exe

will use up to 8 simultaneous jobs while building Galacticus - if your machine has sufficient cores you can use an even larger number.

Download required datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download and unpack the datasets needed at runtime:

.. code-block:: bash

   wget https://github.com/galacticusorg/datasets/archive/master.zip -O datasets.zip
   unzip datasets.zip

and set the ``GALACTICUS_DATA_PATH`` environment variable to the path to these datasets, e.g.:

.. code-block:: bash

   export GALACTICUS_DATA_PATH=/path/to/datasets/folder

Running a minimal example
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A simple test example can be run using:

.. code-block:: bash

   ./Galacticus.exe parameters/quickTest.xml

If successful a ``galacticus.hdf5`` file will be produced containing a small number of galaxies.
