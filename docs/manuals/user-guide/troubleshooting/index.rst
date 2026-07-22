Troubleshooting
===============

This page contains information to help you when Galacticus fails, or doesn't do what you want.

* For errors that occur when you are trying to compile Galacticus from source, see Compile-Time Error Messages.
* For errors that occur when you are running Galacticus, see Run-Time Error Messages.
* For debugging help as you are developing Galacticus, see Debugging.

Beyond actual error messages, sometimes Galacticus may be working (i.e. running without any error message), but not performing as you want. The remainder of this page describes some common problems that you may encounter when running Galacticus and offers advice about how to resolve them.

Warnings about parameters
-------------------------

When Galacticus reads your parameter file it may issue warnings if there is something that it suspects is a problem (but which isn't definitely an error). These warnings will appear in the output log, and are highlighted by the text ``WARNING``. It is recommended that you should check these warnings messages carefully - usually they identify a real problem in your parameter file. If you are sure that the parameter identified is correct as you intended you can add an attribute ``ignoreWarnings="yes"`` to the parameter in question to suppress the error message.

``unrecognized parameter``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you see a message similar to:

.. code-block:: text

   unrecognized parameter [neutrinoMassSum in parameters/transferFunction/] (did you mean [neutrinoMassSummed]?)

this indicates that the named parameter, ``neutrinoMassSum`` in this case, was not recognized as a valid parameter name. Galacticus tells you where in your parameter file this parameter can be found (inside the ``parameters/transferFunction/`` section in this case), and also suggests what it thinks you might have intended - in this case the parameter name you probably wanted was ``neutrinoMassSummed``. If you correct the name of the variable, the warning message will no longer appear.

.. _manual-sec-migratingParameterFiles:

Migrating parameter files to new versions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The names and allowed values of parameters often change between versions of Galacticus. To permit easy and error-free migration between versions a script is provided to translate parameter files from earlier to later versions. To migrate a parameter file simply use:

.. code-block:: bash

   ./scripts/aux/parametersMigrate.py parameters.xml newParameters.xml

By default, this script will translate from the previous to the current version of Galacticus. If your parameter file contains a ``version`` element then this will be used to determine which version of Galacticus the parameter file was constructed for. The migration script will then migrate the parameter file through all intermediate versions to bring it into compliance with the current version. You can also specify input and output versions directly:

.. code-block:: bash

   ./scripts/aux/parametersMigrate.py parameters.xml newParameters.xml --inputVersion 0.9.0 --outputVersion 0.9.3

will convert ``parameters.xml`` from version 0.9.0 syntax to version 0.9.3 syntax.

``multiple copies of parameter``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you see a message similar to:

.. code-block:: text

   multiple copies of parameter [cosmologyFunctions] present - only the first will be utilized

this indicates that the named parameter, ``cosmologyFunctions`` in this case, appears more than once in the same section of the parameter file. In most cases, Galacticus expects only a single copy of a parameter (exceptions to this are inside of specific classes which allow multiple copies of certain parameters, e.g. the ``<nodeOperator value="multi">`` parameter allows multiple ``<nodeOperator....>`` parameter inside of it - these will not cause warning messages). In this case, it is warning you that is is using the first copy, and will ignore the additional copies.

Often, this can happen if you added a parameter to a parameter file without realizing that it was already present. You should check your parameter file carefully, and keep only the copy of the named parameter that you actually want - otherwise Galacticus might not be using the parameter settings that you intended.

Non-deterministic per-tree random number sequences
--------------------------------------------------

Each merger tree in Galacticus is assigned its own random number generator (which is used for building trees via extended Press-Schechter theory, sampling subhalo orbital parameters from distribution functions, etc.). The random seed for this random number generator is set to ``seed``+``treeIndex``, where ``seed`` is typically a value set via the parameter file - for example, for the commonly-used GSL random number generator class a parameter file might contain:

.. code-block:: xml

   <randomNumberGenerator value="GSL">
     <seed value="219"/>
   </randomNumberGenerator>

In the above, ``treeIndex`` is the numerical index of the tree. This ensures that each tree has a unique random number sequence, determined by the ``treeIndex``. As such, the sequence of random numbers for each tree will be identical every time a model is run (allowing for deterministic behavior).

Many random number generators (including the ``GSL`` class) allow for the ``seed`` to be offset by either the OpenMP thread number, the MPI process number, or both. This can be necessary in other applications (e.g. MCMC simulations) where each thread/process is required to produce different sequences of random numbers. However, for tree building this results in non-deterministic behavior - if you re-run the same model, a given tree may be processed by a different thread/process and so will receive a different ``seed``, resulting in completely different tree structure and evolution.

Potentially worse, it's possible in this case for two or more trees to be run with exactly the same seed, resulting in identical (or, at least, strongly correlated) trees.

If such behavior is detected, you will see a warning message:

.. code-block:: text

   WARNING: per-tree random number sequences may not be deterministic

For the ``GSL`` random number generator class, deterministic behavior is the default. If you see the above message you can enforce deterministic behavior by explicitly disabling offsetting of the seed by OpenMP thread number or MPI process number:

.. code-block:: xml

   <randomNumberGenerator value="GSL">
     <seed value="219"/>
     <mpiRankOffset value="false"/>
     <ompThreadOffset value="false"/>
   </randomNumberGenerator>

Low CPU utilization with large numbers of output redshifts
----------------------------------------------------------

If a model is failing to make use of the majority of the available CPU cycles, and you are running a model with a large number of output redshifts the problem may be that I/O to disk is limiting the rate at which merger trees can be processed. I/O occurs through the HDF5 library which provides caching functionality. Therefore, this problem can often be mitigated by expanding the size of the HDF5 library's cache. Galacticus allows you to do this using a set of input parameters:

* ``[hdf5CacheElementsCount]``: HDF5 limits the number of objects that it will store in its cache. Increasing this number will allow more data to be cached and potentially make disk I/O more efficient. We have had good results by setting this number to some factor (e.g. 2) times the product of the number of output redshifts and the number of properties being output in each snapshot.

* ``[hdf5CacheSizeBytes]``: HDF5 also limits the size of the cache in bytes. We have had good results setting this to a factor of a few above the product of, the chunk size (see below) and the size of each output property (8 bytes).

* ``[hdf5SieveBufferSize]``: HDF5 uses a sieve buffer to speed reading/writing of partial datasets. Increasing the buffer size (specified here in bytes) may improve I/O performance. We have had good results using a value of 512Kb.

* ``[hdf5UseLatestFormat]``: Normally, HDF5 selects an internal file format to used based on maximum portability. If you set this option to ``true`` HDF5 will instead use the latest format that it supports—typically allowing it to employ optimizations unavailable to older versions of HDF5. Note that this can make the output file unreadable by older versions of the HDF5 library.

Additionally, you can ensure that compression is switched off in the HDF5 output by setting ``[hdf5CompressionLevel]``\ =-1. Finally, adjusting the HDF5 chunk size via the ``[hdf5ChunkSize]`` parameter may make for more efficient I/O. HDF5 datasets are read/written in chunks of this size. Increasing the size may improve I/O performance.

.. toctree::
   :maxdepth: 1

   compile-time-errors
   run-time-errors
   debugging
   tree-deadlocks
   profiling
