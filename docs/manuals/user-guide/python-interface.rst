.. _manual-sec-pythonInterface:

Python Interface
================

Galacticus exposes a subset of its functionality through a shared library, ``libgalacticus.so``, together with a generated Python module, ``galacticus.py``, that provides a ``ctypes``-based interface to that library. This chapter describes how to build and install the library, how to use it from Python, and how to diagnose common problems.

.. _manual-sec-pythonInterface-build:

Building and Installing the Python Library
------------------------------------------

.. _manual-sec-pythonInterface-build-prerequisites:

Prerequisites
~~~~~~~~~~~~~

Building the Python library requires the same toolchain used to build the main Galacticus executable (a Fortran 95 compiler, a C compiler, and the :term:`GNU` Make build system) plus the following additional requirements:

* Python 3 (:math:`\geq` 3.9), together with the Python packages used by the build-time preprocessing and code-generation scripts (such as ``scripts/build/preprocess.py`` and ``scripts/build/codeDirectivesParse.py``).  These are declared in the top-level ``pyproject.toml`` and installed in one step from the repository root with:

  .. code-block:: none

     pip install -e .

  This editable install places the modules under ``python/`` on the Python import path and installs ``numpy``, ``scipy``, ``h5py``, ``lxml``, ``matplotlib``, and the other declared dependencies.  Alternatively, exporting ``PYTHONPATH=$(pwd)/python`` before invoking ``make`` is sufficient for the build itself (the third-party dependencies must then be available in the active Python environment by some other means).

No additional Python packages (e.g., ``numpy``) are required merely to *use* the library---the generated module relies only on the Python standard library's ``ctypes`` module.

.. _manual-sec-pythonInterface-build-build:

Build
~~~~~

The shared library and its Python interface module are built by invoking ``make`` with the ``GALACTICUS_BUILD_OPTION=lib`` flag and specifying the ``libgalacticus.so`` target:

.. code-block:: none

   make GALACTICUS_BUILD_OPTION=lib libgalacticus.so

This build step performs the following actions:

#. The Python script ``scripts/build/libraryInterfaces.py`` is run. It reads the list of classes to be exposed from ``source/libraryClasses.xml``, generates C-binding Fortran wrapper code in ``${BUILDPATH}/libgalacticus/``, and writes the Python interface module to ``galacticus.py`` in the Galacticus root directory.
#. The wrapper code is preprocessed and compiled.
#. Everything is linked into the shared library ``libgalacticus.so`` in the Galacticus root directory.

Both ``galacticus.py`` and ``libgalacticus.so`` are placed directly in the Galacticus root directory immediately after the build completes.

.. _manual-sec-pythonInterface-build-install:

Installation
~~~~~~~~~~~~

After building, the two artifacts must be placed in a directory layout that ``galacticus.py`` expects. The standard layout used by the Galacticus release bundles is:

.. code-block:: none

   <prefix>/
     galacticus/
       lib/
         libgalacticus.so
         <runtime dependency .so files>
       python/
         galacticus.py

To create this layout from the build directory, run:

.. code-block:: none

   mkdir -p galacticus/lib galacticus/python
   mv galacticus.py     galacticus/python/
   mv libgalacticus.so  galacticus/lib/

.. _manual-sec-pythonInterface-build-pythonpath:

Making the Module Importable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Python must be able to find ``galacticus.py``. Set the ``PYTHONPATH`` environment variable to the ``galacticus/python`` directory before running any Python scripts that import the module:

.. code-block:: none

   export PYTHONPATH=/path/to/prefix/galacticus/python:$PYTHONPATH

where ``/path/to/prefix`` is the directory that contains the ``galacticus/`` tree created in Section :galacticus-ref:`pythonInterface:build:install`.

The module also requires that ``galacticus/lib/libgalacticus.so`` can be found *relative to the current working directory* at import time. Therefore, Python scripts using the library should be run from ``/path/to/prefix``:

.. code-block:: none

   cd /path/to/prefix
   python3 myScript.py

.. _manual-sec-pythonInterface-usage:

Using the Python Interface
--------------------------

.. _manual-sec-pythonInterface-usage-import:

Importing the Module
~~~~~~~~~~~~~~~~~~~~

Once the prerequisites in Section :galacticus-ref:`pythonInterface:build:pythonpath` are satisfied, the module is imported with a standard Python import:

.. code-block:: none

   import galacticus

Importing the module automatically loads ``libgalacticus.so`` via ``ctypes`` and calls the internal Galacticus initialisation routine ``libGalacticusInitL``, which sets up event hooks and the HDF5 access lock.

.. _manual-sec-pythonInterface-usage-classes:

Class Mapping
~~~~~~~~~~~~~

For each Galacticus *functionClass* listed in ``source/libraryClasses.xml``, the generated module contains:

* An abstract base class bearing the functionClass name (e.g. ``galacticus.cosmologyParameters``).
* One concrete subclass for each concrete implementation that is not explicitly excluded in ``libraryClasses.xml`` (e.g. ``galacticus.cosmologyParametersSimple``).

Concrete classes are instantiated by calling their constructor with the same arguments that the corresponding Galacticus Fortran constructor expects, converted to their Python equivalents (floats for double-precision reals, booleans for logicals, integers for integers, and other class instances for object arguments). Methods exposed on the Fortran side are available as Python methods of the same name.

The set of currently exposed functionClasses is defined in ``source/libraryClasses.xml`` and includes, among others: ``cosmologyParameters``, ``cosmologyFunctions``, ``powerSpectrumPrimordial``, ``transferFunction``, ``linearGrowth``, ``cosmologicalMassVariance``, ``criticalOverdensity``, ``virialDensityContrast``, and ``haloMassFunction``.

.. _manual-sec-pythonInterface-usage-example:

Basic Example
~~~~~~~~~~~~~

The following example constructs a simple flat :math:`\Lambda`\ CDM cosmology and evaluates some basic quantities:

.. code-block:: none

   import galacticus

   # Cosmological parameters: OmegaMatter, OmegaBaryon, OmegaDarkEnergy,
   #                          temperatureCMB, HubbleConstant
   cospar = galacticus.cosmologyParametersSimple(0.3, 0.045, 0.7, 2.78, 70.0)
   print("Omega_Matter =", cospar.OmegaMatter())
   print("Hubble constant =", cospar.HubbleConstant(), "km/s/Mpc")

   # Cosmological functions
   cosfun = galacticus.cosmologyFunctionsMatterLambda(cospar)
   t0 = cosfun.cosmicTime(1.0)
   print("Age of the universe =", t0, "Gyr")
   print("Expansion factor at t =", t0, "Gyr:", cosfun.expansionFactor(t0))

.. _manual-sec-pythonInterface-usage-workflow:

Typical Workflows
~~~~~~~~~~~~~~~~~

A common use-case is to build up a chain of Galacticus objects that mirrors the parameter-file object hierarchy. For example, to compute the halo mass function:

.. code-block:: none

   import galacticus

   cospar  = galacticus.cosmologyParametersSimple(0.3, 0.045, 0.7, 2.78, 70.0)
   cosfun  = galacticus.cosmologyFunctionsMatterLambda(cospar)
   dmpart  = galacticus.darkMatterParticleCDM()
   transfn = galacticus.transferFunctionCAMB(
                 dmpart, cospar, cosfun, 0.0, 0)
   lingrow = galacticus.linearGrowthCollisionlessMatter(cospar, cosfun)
   psprim  = galacticus.powerSpectrumPrimordialPowerLaw(
                 0.965, 0.0, 0.0, 1.0, False)
   pswin   = galacticus.powerSpectrumWindowFunctionTopHat(cospar)
   pstran  = galacticus.powerSpectrumPrimordialTransferredSimple(
                 psprim, transfn, lingrow)
   cmv     = galacticus.cosmologicalMassVarianceFilteredPower(
                 sigma8=0.8, tolerance=1.0e-4, toleranceTopHat=1.0e-4,
                 nonMonotonicIsFatal=True, monotonicInterpolation=False,
                 truncateAtParticleHorizon=False,
                 cosmologyParameters_=cospar, cosmologyFunctions_=cosfun,
                 linearGrowth_=lingrow,
                 powerSpectrumPrimordialTransferred_=pstran,
                 powerSpectrumWindowFunction_=pswin)
   crover  = galacticus.criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(
                 lingrow, cosfun, cmv, dmpart, True)
   hmf     = galacticus.haloMassFunctionShethTormen(
                 cospar, cmv, crover, 0.707, 0.3, 0.322183)
   m = 1.0e12  # Solar masses
   print("dn/dln(M) at M=1e12 Msun, z=0:", hmf.differential(13.8, m) * m)

.. _manual-sec-pythonInterface-usage-limitations:

Notes and Limitations
~~~~~~~~~~~~~~~~~~~~~

* **Object lifetimes:** Galacticus objects created via the Python interface are owned by the Python garbage collector. Destructors are called automatically when objects go out of scope, but care should be taken not to hold references to a child object longer than its parent---for example, do not use a ``cosmologyFunctions`` object after the ``cosmologyParameters`` object it was constructed from has been garbage-collected.
* **Limited class coverage:** Only the functionClasses listed in ``source/libraryClasses.xml`` are exposed. Classes that rely on unlimited polymorphic constructor arguments, or that have been explicitly excluded, may not be available.
* **No direct HDF5 I/O:** The Python interface does not currently support reading or writing Galacticus output files; use ``h5py`` or similar tools for that purpose.
* **Single Galacticus instance:** Only one ``libgalacticus.so`` instance can be loaded per Python process, because the library maintains global Fortran state.

.. _manual-sec-pythonInterface-troubleshooting:

Troubleshooting
---------------

.. _manual-sec-pythonInterface-troubleshooting-moduleNotFound:

``ModuleNotFoundError: No module named 'galacticus'``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Python cannot locate ``galacticus.py``. Ensure that ``PYTHONPATH`` includes the directory containing ``galacticus.py`` (i.e. ``galacticus/python/`` relative to the install prefix):

.. code-block:: none

   export PYTHONPATH=/path/to/prefix/galacticus/python:$PYTHONPATH

Then verify:

.. code-block:: none

   python3 -c "import galacticus; print('OK')"

.. _manual-sec-pythonInterface-troubleshooting-cannotOpen:

``OSError`` / ``libgalacticus.so: cannot open shared object file``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``ctypes`` call inside ``galacticus.py`` constructs the library path as:

.. code-block:: none

   os.path.join(os.getcwd(), "galacticus/lib/libgalacticus.so")

This means Python must be run from the directory that *contains* the ``galacticus/`` tree, not from inside it:

.. code-block:: none

   # Correct: run from the install prefix
   cd /path/to/prefix
   python3 myScript.py

   # Incorrect: run from inside galacticus/
   cd /path/to/prefix/galacticus
   python3 ../myScript.py   # will fail to find the library

.. _manual-sec-pythonInterface-troubleshooting-deps:

Missing Runtime Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If ``libgalacticus.so`` loads but then immediately raises an ``OSError`` about a missing ``.so`` file (e.g. ``libgfortran.so.5: cannot open shared object file``), the runtime library dependencies were not copied during packaging. Re-run the packaging step described in Section :galacticus-ref:`pythonInterface:build:install`, making sure that all non-system libraries are included. To inspect what is needed:

.. code-block:: none

   ldd /path/to/prefix/galacticus/lib/libgalacticus.so

Any library listed as *not found* must be located and copied into ``galacticus/lib/``.

Alternatively, set ``LD_LIBRARY_PATH`` to include the directory containing the missing library. For example, on a system with GCC 12 installed on an x86_64 architecture:

.. code-block:: none

   export LD_LIBRARY_PATH=/usr/lib/gcc/x86_64-linux-gnu/12:$LD_LIBRARY_PATH

Replace ``12`` with your installed GCC version (e.g. ``11``, ``13``) and adjust the architecture component (e.g. ``aarch64-linux-gnu``) to match your system.

.. _manual-sec-pythonInterface-troubleshooting-abi:

Wrong Python Version or ABI Mismatch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``galacticus.py`` module is a pure Python 3 script that uses only the standard library; it has no compiled extension module and therefore no Python ABI requirement of its own. However, if the *Fortran runtime* linked into ``libgalacticus.so`` was built against a different version of ``libgfortran`` than the one on the system, you may see symbol errors at load time. In that case, ensure that the Fortran compiler version used to build ``libgalacticus.so`` matches what is installed. Inspect the required ``libgfortran`` version with:

.. code-block:: none

   objdump -p /path/to/prefix/galacticus/lib/libgalacticus.so \
       | grep GFORTRAN

.. _manual-sec-pythonInterface-troubleshooting-symbol:

Runtime Errors: Symbol Not Found
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An ``AttributeError`` or ``OSError`` of the form *symbol not found* when calling a method usually means the ``libgalacticus.so`` on disk is older than ``galacticus.py``, or vice versa (e.g. after a partial rebuild). Rebuild the library completely and re-install both ``galacticus.py`` and ``libgalacticus.so`` together:

.. code-block:: none

   make GALACTICUS_BUILD_OPTION=lib libgalacticus.so
   mv galacticus.py    /path/to/prefix/galacticus/python/
   mv libgalacticus.so /path/to/prefix/galacticus/lib/

Always keep ``galacticus.py`` and ``libgalacticus.so`` in sync; they are generated together and must be deployed together.

.. _manual-sec-pythonInterface-troubleshooting-init:

Segmentation Fault or Fortran Initialization Error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If Python crashes with a segmentation fault immediately after ``import galacticus``, the most likely cause is that the Fortran runtime was not properly initialized. This can happen if ``libgalacticus.so`` is loaded more than once in the same process, or if the ``GALACTICUS_EXEC_PATH`` and ``GALACTICUS_DATA_PATH`` environment variables are not set. Set these before running:

.. code-block:: none

   export GALACTICUS_EXEC_PATH=/path/to/galacticus/source
   export GALACTICUS_DATA_PATH=/path/to/galacticus/datasets
