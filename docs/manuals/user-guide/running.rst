Running Galacticus
==================

If you have not yet installed Galacticus you should follow the instructions in :doc:`installation/index` to do so.

Setting the Environment
-----------------------

Before running Galacticus you will need to set two environment variables which specify where the Galacticus source code and datasets can be found. First, the environment variable ``GALACTICUS_EXEC_PATH`` should be set to the full path to the build directory. Second, the environment variable ``GALACTICUS_DATA_PATH`` should be set to the full path to the ``datasets`` directory which was created when you installed Galacticus.

Running Galacticus
------------------

Galacticus is run using

.. code-block:: none

    Galacticus.exe <parameterFile>

where ``parameterFile`` is the name of the file containing parameter settings for Galacticus. Galacticus will display messages indicating its progress as it runs.

A simple example, useful to check that everything is working as expected for you, is to do:

.. code-block:: none

    Galacticus.exe $GALACTICUS_EXEC_PATH/parameters/quickTest.xml

This will run a small model. You should expect to see output which looks something like this:

.. code-block:: none

                 ##
      ####        #                  #
     #   #        #             #
    #       ###   #  ###   ### ###  ##   ### ## ##   ##
    #       #  #  #  #  # #  #  #    #  #  #  #  #  #
    #   ###  ###  #   ### #     #    #  #     #  #   #
     #   #  #  #  #  #  # #     #    #  #     #  #    #
      ####  #### ### ####  ###   ## ###  ###   #### ##

     2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
     2017, 2018, 2019, 2020, 2021, 2022
      - Andrew Benson

   MM: -> Begin task: merger tree evolution
    0:     -> Evolving tree number 1 {1}
    9:     -> Evolving tree number 3 {1}
    7:     -> Evolving tree number 4 {1}
    3:     -> Evolving tree number 2 {1}
    0:         Output tree data at t=  13.79 Gyr
    1:     -> Evolving tree number 5 {1}
   11:     -> Evolving tree number 6 {1}
    0:     <- Finished tree
    9:         Output tree data at t=  13.79 Gyr
    9:     <- Finished tree
    3:         Output tree data at t=  13.79 Gyr
    3:     <- Finished tree
    7:         Output tree data at t=  13.79 Gyr
    7:     <- Finished tree
    1:         Output tree data at t=  13.79 Gyr
    1:     <- Finished tree
   11:         Output tree data at t=  13.79 Gyr
   11:     <- Finished tree
   MM: <- Done task: merger tree evolution

In this simple model, Galacticus is told to build six merger trees, and form galaxies within them, outputting the results at :math:`z=0`. After the initial "Galacticus" banner is displayed there are several messages shown. Each line begins with either "``MM:``", or a number. Galacticus runs in parallel using the available cores on your computer---these prefixes on each line tell you which parallel thread is sending the message. "``MM:``" means that the message is from the main thread (and that Galacticus is currently running in a serial part of the code), while a numerical prefixes gives the number of that parallel thread reporting the message (starting from 0).

The first message, from the main thread, tells you that the "merger tree evolution" task has begun---in this task, Galacticus forms galaxies within a set of merger trees. After that initial message you'll see that Galacticus enters parallel calculation and messages start being reported from each parallel thread. The threads report when they begin evolving a merger tree (in this case merger trees are numbered consecutively, starting from 1), when they are outputting the data from that tree (at a time of :math:`13.79` Gyr, which corresponds to :math:`z=0` in this model), and when they have finished the tree. Once all trees are finished there is a final message from the main thread (parallel calculation is over at this point) indicating that the merger tree evolution task is finished. Galacticus then exits, since it has finished all of the work we asked it to do.

The results will have been written to an output file named ``galacticus.hdf5``, which you can see by doing:

.. code-block:: none

   ls -l galacticus.hdf5

which will show something similar to:

.. code-block:: none

   -rw-r--r-- 1 abenson users 568601 Oct 15 20:58 galacticus.hdf5

Galacticus outputs use the `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ format---libraries for accessing HDF5 files exist in all major languages, including Python\ [#]_. We will explore the structure of the output file in more detail in Section :galacticus-ref:`outputFile`, but for now we'll simply explore some of the key features using the command line.

The content of the output file can be explored using the "``h5ls``" tool, for example:

.. code-block:: none

   $ h5ls galacticus.hdf5
   Build                    Group
   Outputs                  Group
   Parameters               Group
   Version                  Group

The file contains several groups---the data corresponding to the galaxies that were formed can be found in the ``Outputs`` group:

.. code-block:: none

   $ h5ls galacticus.hdf5/Outputs
   Output1                  Group

In this case there is only a single output, corresponding to :math:`z=0`, we can explore that group using:

.. code-block:: none

   $ h5ls galacticus.hdf5/Outputs/Output1
   mergerTreeCount          Dataset {6/Inf}
   mergerTreeIndex          Dataset {6/Inf}
   mergerTreeSeed           Dataset {6/Inf}
   mergerTreeStartIndex     Dataset {6/Inf}
   mergerTreeWeight         Dataset {6/Inf}
   nodeData                 Group

In this output group we find some datasets which give information on the merger trees that were used, plus another group, "``nodeData``", inside of which the galaxy data is stored:

.. code-block:: none

   $ h5ls galacticus.hdf5/Outputs/Output1/nodeData
   basicMass                Dataset {7/Inf}
   basicTimeLastIsolated    Dataset {7/Inf}
   blackHoleCount           Dataset {7/Inf}
   blackHoleMass            Dataset {7/Inf}
   blackHoleSpin            Dataset {7/Inf}
   darkMatterProfileScale   Dataset {7/Inf}
   diskAbundancesGasMetals  Dataset {7/Inf}
   diskAbundancesStellarMetals Dataset {7/Inf}
   diskAngularMomentum      Dataset {7/Inf}
   diskMassGas              Dataset {7/Inf}
   diskMassStellar          Dataset {7/Inf}
   diskRadius               Dataset {7/Inf}
   diskVelocity             Dataset {7/Inf}
   hotHaloAbundancesMetals  Dataset {7/Inf}
   hotHaloAngularMomentum   Dataset {7/Inf}
   hotHaloMass              Dataset {7/Inf}
   hotHaloOuterRadius       Dataset {7/Inf}
   hotHaloOutflowedAbundancesMetals Dataset {7/Inf}
   hotHaloOutflowedAngularMomentum Dataset {7/Inf}
   hotHaloOutflowedMass     Dataset {7/Inf}
   hotHaloUnaccretedMass    Dataset {7/Inf}
   nodeIndex                Dataset {7/Inf}
   nodeIsIsolated           Dataset {7/Inf}
   parentIndex              Dataset {7/Inf}
   satelliteBoundMass       Dataset {7/Inf}
   satelliteIndex           Dataset {7/Inf}
   satelliteMergeTime       Dataset {7/Inf}
   siblingIndex             Dataset {7/Inf}
   spheroidAbundancesGasMetals Dataset {7/Inf}
   spheroidAbundancesStellarMetals Dataset {7/Inf}
   spheroidAngularMomentum  Dataset {7/Inf}
   spheroidMassGas          Dataset {7/Inf}
   spheroidMassStellar      Dataset {7/Inf}
   spheroidRadius           Dataset {7/Inf}
   spheroidVelocity         Dataset {7/Inf}
   spinSpin                 Dataset {7/Inf}

Each dataset here is an array containing the named property of each galaxy formed in the model. In this tiny example model only 7 galaxies were formed. To see the values of each property we can do:

.. code-block:: none

   $ h5ls -d galacticus.hdf5/Outputs/Output1/nodeData/diskMassStellar
   diskMassStellar          Dataset {7/Inf}
       Data:
           (0) 0, 719592.675733525, 0, 36664762.1643418, 75533178.4182319, 0, 764729109.215361

which lists the mass of stars in each galaxy disk (in units of :math:`\mathrm{M}_\odot`) (note that some of them are zero---these halos in the merger tree either formed no galaxy, or formed a galaxy with no disk component).

These data can be extracted and analyzed using any software or language that supports reading HDF5 files. The `Dendros <https://github.com/galacticusorg/dendros>`_ package provides ready-made tools for analyzing and plotting Galacticus output, including on-the-fly analyses, MCMC chain diagnostics, and posterior corner plots.

Dry Runs
~~~~~~~~

You can tell Galacticus to parse your parameter file, report any warnings, and write the parameters to the output HDF5 file, but then do nothing else (i.e. don't actually run the model) by adding the ``--dry-run`` option, for example:

.. code-block:: none

    Galacticus.exe parameters.xml --dry-run

This can be useful to check that your parameter file is valid, and allow you to explore the values of any parameters that were set to defaults (as these will have been output to the HDF5 file).

Run-time Estimation and Progress Reporting
------------------------------------------

When evolving a set of merger trees (the ``evolveForests`` task), Galacticus can report on the progress of the run and estimate how long it will take. This is useful both for monitoring a long-running job and, when a cost model is supplied, for estimating the total run time up-front so that an appropriately-sized cluster job can be requested.

Progress reports
~~~~~~~~~~~~~~~~~

By default, Galacticus reports whole-run progress during merger tree evolution. The reporting cadence is controlled by the ``[timeIntervalReportProgress]`` parameter of the ``evolveForests`` task, which gives the minimum wall-clock interval (in seconds; default ``60``) between reports; setting it to zero or a negative value disables progress reporting. For example, to report every ten seconds:

.. code-block:: xml

    <task value="evolveForests">
      <timeIntervalReportProgress value="10.0"/>
    </task>

Three kinds of message are produced:

* a **start-of-run** line, giving the number of trees to be processed and (if a cost model is available; see below) an estimate of the total run time;
* **throttled progress** reports during the run, giving the number of trees processed, the percentage complete, the elapsed time, and an estimated time remaining;
* an **end-of-run summary**, giving the total wall time, total core-time, and the mean and maximum per-tree processing times---phrased to be directly useful when sizing a subsequent job.

The estimated time remaining is *self-calibrating*: it is continually corrected using the ratio of actual to predicted processing time measured over the trees completed so far in the current run. It therefore converges to an accurate value as the run proceeds, even for a model it has never seen before. When trees are built in order of descending mass (the default), the most expensive trees complete first, so the estimate converges quickly.

Under MPI, the progress report is made by the master process only, and the tree count shown is the globally-consistent count of forests claimed across all processes.

Cost models: estimating run time before the run starts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Without any prior information, no estimate of the *total* run time can be made at the start of a run (progress is instead reported purely in terms of the number of trees completed). To obtain a start-of-run estimate---for example, to decide how much wall time to request for a cluster job---you can supply a *cost model* calibrated from a previous run of the same (or a similar) model. Using such a model is a two-step process.

**Step 1: calibrate.** Enable the ``treeProcessingTimer`` merger tree operator in a (typically shorter, or representative) run:

.. code-block:: xml

    <mergerTreeOperator value="treeProcessingTimer"/>

This records the processing time of each tree and, at the end of the run, fits the (base-10 logarithm of the) processing time as a quadratic function of the (base-10 logarithm of the) tree mass, and separately of the tree node count. The fit coefficients, residual scatter, and range of validity are written into the ``metaData/treeTiming`` group of that run's output file, alongside the raw per-tree timing data.

**Step 2: estimate.** In a subsequent run, supply that output file as the cost model via the ``file`` :galacticus-class:`metaTreeProcessingTime` class:

.. code-block:: xml

    <metaTreeProcessingTime value="file">
      <fileName value="previousRun.hdf5"/>
    </metaTreeProcessingTime>

Galacticus will then report a start-of-run estimate of the total core-time and the projected wall time given the number of workers, and will use the cost model to weight the running time-remaining estimate by the *predicted* work of each tree (which is much more accurate early in the run than a simple count of trees, especially when a few large trees dominate the cost). If a maximum wall time (``[walltimeMaximum]``) is set and the projected time exceeds it, a warning is issued.

Because the time to process a tree depends on the physics of the specific model (as well as on mass, resolution, number of outputs, and so on), the estimate is necessarily *conditional on the calibration model*; this is noted in the reported message.

For merger trees that are **built** on-the-fly, the root masses are known before evolution, so the mass-based fit is used. For merger trees that are **read** from files, the root masses are not known until a tree has been evolved, but the number of nodes in each tree *is* known from the file metadata; in this case the node-count-based fit is used automatically (the ``file`` cost model reads it from the same output file).

The ``file`` cost model also accepts the legacy XML format (a ``<timing><fit>`` document of mass-based coefficients); an ``.hdf5`` or ``.h5`` file name selects the HDF5 format described above. Fit coefficients can be produced, or combined across several output files, using the ``scripts/analysis/treeProcessingTimeFit.py`` helper script (see :doc:`analysis`).

Estimate-only (run-time sizing) mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To obtain the start-of-run estimate *without* actually evolving any trees---for instance, when sizing a job request---set the ``[estimateRunTimeOnly]`` parameter of the ``evolveForests`` task:

.. code-block:: xml

    <task value="evolveForests">
      <estimateRunTimeOnly value="true"/>
    </task>

Galacticus will compute and report the estimate and then exit. This requires a cost model and tree census to produce a useful estimate. (This differs from the ``--dry-run`` command-line option described above, which parses the parameter file and exits without attempting any estimate.)

.. [#] `h5py <https://www.h5py.org/>`_
