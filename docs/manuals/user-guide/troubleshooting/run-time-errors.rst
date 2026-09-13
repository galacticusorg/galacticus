Run-Time Error Messages
=======================

As you are running Galacticus it's quite likely that you'll eventually run into an error. This page lists common types of such run-time errors, along with guidance on what they mean and how to fix them.

For *compile-time* errors (i.e. things that go wrong when you are building the ``Galacticus.exe`` executable) you should refer to :doc:`compile-time-errors`. For other problems (Galacticus is not performing as you expect, or doing what you want) see the general :doc:`index`.

Reading a fatal error message
-----------------------------

A deliberate fatal error raised by Galacticus has a fixed layout:

.. code-block:: text

   Fatal error:
   unable to allocate `particleID` (1000000000000 elements)
   HELP: the run needs more memory than is available to it. Reduce the size of the problem (for
   example, the number of trees or particles being processed), reduce the number of OpenMP threads,
   since each holds its own copy of much of the state, or run where more memory is available
    Occurred at:
      subroutine:randomImport
            file:nBody/import/random.F90   [line 159]
     => Error occurred in master thread

The first line after ``Fatal error:`` is the message itself. Many messages are followed by a ``HELP:`` line (shown in green on a terminal) which states the most likely fix; try that first. The ``Occurred at:`` block names the procedure, module, and source file that raised the error, and the final line identifies the OpenMP thread (and, under MPI, the process and host). A stack trace follows the message; see :doc:`debugging` for how to read it. When reporting a problem, include the whole block.

Error Message Diagnosis and Reporting
-------------------------------------

If a Galacticus run fails, this flowchart can help to diagnose a failed Galacticus run. You can click the highlighted "See here" boxes to find solutions. (Note that you may need to ctrl-click these links for them to open.) If the flowchart doesn't allow you to resolve the problem, you can ask questions from the community in our `discussion forums <https://github.com/galacticusorg/galacticus/discussions>`_, or open a `bug report <https://github.com/galacticusorg/galacticus/issues/new?assignees=abensonca&labels=bug&projects=&template=Bug-Report.yml&title=%5BBug%5D%3A+>`_.

.. mermaid::

   flowchart TB
     Start{{Error message includes the phrase 'Fatal error:'?}}
     Start --yes--> DynamicNotExistQ
     Start --no--> System
     DynamicNotExistQ{{"Error message includes the phrase
                       'dataset 'xyz' does not exist in 'datasets/dynamic/..'?"}}
     DynamicNotExistQ --yes--> DynamicNotExist
     DynamicNotExistQ --no--> ParameterEmptyValueQ
     DynamicNotExist(See here)
     style DynamicNotExist fill:#74c7db
     click DynamicNotExist href "#dataset-xyz-does-not-exist-in-datasets-dynamic" "dummy"
     ParameterEmptyValueQ{{"Error message includes the phrase
                       'empty value in parameter [xyz]'?"}}
     ParameterEmptyValueQ --yes--> ParameterEmptyValue
     ParameterEmptyValueQ --no--> RamPressureComponentQ
     ParameterEmptyValue(See here)
     style ParameterEmptyValue fill:#74c7db
     click ParameterEmptyValue href "#empty-value-in-parameter" "dummy"
     RamPressureComponentQ{{"Error message includes the phrase
                       'only ＂xyz＂ components are supported by the ＂abc＂ ramPressureStripping class'?"}}
     RamPressureComponentQ --yes--> RamPressureComponent
     RamPressureComponentQ --no--> DataFileQ
     RamPressureComponent(See here)
     style RamPressureComponent fill:#74c7db
     click RamPressureComponent href "#inconsistent-assumptions-for-ram-pressure-models" "dummy"
     DataFileQ{{"Error message includes the phrase
                       'Unable to find data file'?"}}
     DataFileQ --yes--> DataFile
     DataFileQ --no--> CloseFileQ
     DataFile(See here)
     style DataFile fill:#74c7db
     click DataFile href "#unable-to-find-data-file" "dummy"
     CloseFileQ{{"Error message includes the phrase
                       'unable to close file object '/dev/shm/glcTmpPar...'?"}}
     CloseFileQ --yes--> CloseFile
     CloseFileQ --no--> AllocateQ
     CloseFile(See here)
     style CloseFile fill:#74c7db
     click CloseFile href "#temporary-parameter-files-fill-dev-shm" "dummy"
     AllocateQ{{"Error message includes the phrase
                       'unable to allocate'?"}}
     AllocateQ --yes--> Allocate
     AllocateQ --no--> RecursiveQ
     Allocate(See here)
     style Allocate fill:#74c7db
     click Allocate href "#unable-to-allocate-memory" "dummy"
     RecursiveQ{{"Error message includes the phrase
                       'composites a member of its own class'?"}}
     RecursiveQ --yes--> Recursive
     RecursiveQ --no--> DeadlockQ
     Recursive(See here)
     style Recursive fill:#74c7db
     click Recursive href "#recursive-object-construction" "dummy"
     DeadlockQ{{"Error message includes the phrase
                       'merger tree appears to be deadlocked'?"}}
     DeadlockQ --yes--> Deadlock
     DeadlockQ --no--> HistoryGridQ
     Deadlock(See here)
     style Deadlock fill:#74c7db
     click Deadlock href "#merger-tree-appears-to-be-deadlocked" "dummy"
     HistoryGridQ{{"Error message includes the phrase
                       'serialized history increment cannot extend the time grid'?"}}
     HistoryGridQ --yes--> HistoryGrid
     HistoryGridQ --no--> DownloadQ
     HistoryGrid(See here)
     style HistoryGrid fill:#74c7db
     click HistoryGrid href "#star-formation-history-cannot-extend-its-time-grid" "dummy"
     DownloadQ{{"Error message includes the phrase
                       'failed to download from'?"}}
     DownloadQ --yes--> Download
     DownloadQ --no--> RadiusZeroQ
     Download(See here)
     style Download fill:#74c7db
     click Download href "#failed-to-download-a-file" "dummy"
     RadiusZeroQ{{"Error message includes the phrase
                       'Radius specifier evaluates to a radius of zero'?"}}
     RadiusZeroQ --yes--> RadiusZero
     RadiusZero(See here)
     style RadiusZero fill:#74c7db
     click RadiusZero href "#radius-specifier-evaluates-to-a-radius-of-zero" "dummy"
     System{{Error message includes the phrase:}}
     System --Floating point exception--> FPE
     System --Segmentation fault--> SegFault
     System --Bus error--> Bus
     System --Illegal instruction--> Illegal
     System --Command terminated by signal 9--> Signal9
     FPE(See here)
     style FPE fill:#74c7db
     click FPE href "#floating-point-errors-and-segfaults" "dummy"
     SegFault(See here)
     style SegFault fill:#74c7db
     click SegFault href "#floating-point-errors-and-segfaults" "dummy"
     Bus(See here)
     style Bus fill:#74c7db
     click Bus href "#bus-errors-illegal-instructions-and-signal-9" "dummy"
     Illegal(See here)
     style Illegal fill:#74c7db
     click Illegal href "#bus-errors-illegal-instructions-and-signal-9" "dummy"
     Signal9(See here)
     style Signal9 fill:#74c7db
     click Signal9 href "#bus-errors-illegal-instructions-and-signal-9" "dummy"

Floating point errors and Segfaults
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A floating point exception or segmentation fault ("segfault") indicate a bug in Galacticus. These should never happen! If they do, please file a bug report using this `form <https://github.com/galacticusorg/galacticus/issues/new?assignees=abensonca&labels=bug&projects=&template=Bug-Report.yml&title=%5BBug%5D%3A+>`_, providing as much information as you can so that we can diagnose and fix the problem.

Before reporting, check whether the stack trace names a known case:

* ``lossConeTabulate`` (``satellites/merging/virial_orbits/loss_cone.F90``): the ``lossCone`` virial orbit class can raise an invalid floating point operation when its orbital tabulation is first built at a very early epoch (redshift above about 8), because the environmental boost factor becomes 0/0 for the most massive tabulated hosts. This is tracked as `issue #1426 <https://github.com/galacticusorg/galacticus/issues/1426>`_; until it is fixed, the workaround is to avoid trees whose progenitors request an orbit at such early times (for example by raising the mass resolution), or to use a different ``virialOrbit`` class.

Bus errors, Illegal instructions, and Signal 9
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bus errors, Illegal instructions, and Signal 9 can sometimes occur if Galacticus runs out of memory. To see if this is a likely cause of the error, first look at the output log from Galacticus - you should see some reports of memory usage, which will look something like this:

.. code-block:: text

   Memory usage: 222.029 MB / 62.813 GB

The first number if the memory currently used, the second is the total memory available on the system. If the first number is close to the second, an out of memory error is likely.

But, keep in mind that out of memory errors can occur even if the reported memory usage is lower than the memory available. Other processes running on the system will be using some of that memory and, if you're running an a compute cluster, your job may have limited memory assigned to it. In that case, try adjusting your job submission parameters to request more memory.

Also keep in mind that Galacticus reports memory usage only periodically, so actual usage may be greater than the last reported value. You can look at peak memory usage using the ``time`` command. If you normally run a model using, e.g.:

.. code-block:: console

   ./Galacticus.exe parameters.xml

try instead:

.. code-block:: console

   /usr/bin/time -v ./Galacticus.exe parameters.xml

When this model finishes (or fails with a bus error) you'll see a report like this:

.. code-block:: text

           Command being timed: "./Galacticus.exe parameters.xml"
           User time (seconds): 31.99
           System time (seconds): 3.24
           Percent of CPU this job got: 515%
           Elapsed (wall clock) time (h:mm:ss or m:ss): 0:06.83
           Average shared text size (kbytes): 0
           Average unshared data size (kbytes): 0
           Average stack size (kbytes): 0
           Average total size (kbytes): 0
           Maximum resident set size (kbytes): 4046928
           Average resident set size (kbytes): 0
           Major (requiring I/O) page faults: 1
           Minor (reclaiming a frame) page faults: 273030
           Voluntary context switches: 409900
           Involuntary context switches: 12038
           Swaps: 0
           File system inputs: 1824
           File system outputs: 432
           Socket messages sent: 0
           Socket messages received: 0
           Signals delivered: 0
           Page size (bytes): 4096
           Exit status: 0

The line:

.. code-block:: text

           Maximum resident set size (kbytes): 4046928

reports the peak memory used. (Note that some older versions of the ``time`` command `incorrectly report <https://access.redhat.com/errata/RHBA-2015:0710.html>`_ the memory use as four times higher than it actually is...)

If an out of memory issue does not seem to be the cause of the error please file a bug report using this `form <https://github.com/galacticusorg/galacticus/issues/new?assignees=abensonca&labels=bug&projects=&template=Bug-Report.yml&title=%5BBug%5D%3A+>`_, providing as much information as you can so that we can diagnose and fix the problem.

``dataset 'xyz' does not exist in 'datasets/dynamic/.....``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When you run Galacticus, lots of files get created in the ``datasets/dynamic/`` path - these are often tabulations of functions that are expensive to evaluate (so they are computed once, saved to file, and then quickly re-read on future runs).

Error messages which state:

.. code-block:: text

   dataset 'xyz' does not exist in 'datasets/dynamic/.....

usually mean that the creation of one of these files was interrupted - leaving the file in a non-usable state. This could happen because some other error occurred, or because you killed a model.

If you see this error, the first thing to try is to simply remove the file mentioned in the error message, and then try running your model again. All files in ``datasets/dynamic/`` will be re-made as needed - so deleting them will just cause Galacticus to remake the file - hopefully successfully this time.

If the same error occurs again, please file a bug report using this `form <https://github.com/galacticusorg/galacticus/issues/new?assignees=abensonca&labels=bug&projects=&template=Bug-Report.yml&title=%5BBug%5D%3A+>`_, providing as much information as you can so that we can diagnose and fix the problem.

Empty value in parameter
------------------------

If you receive an error message of the form:

.. code-block:: text

   empty value in parameter [xyz]

this means that, in your parameter file, the named parameter (``xyz``) has no content in its ``value`` attribute. For example:

.. code-block:: text

   <xyz value=""/>

This is invalid - all parameters with a ``value`` attribute *must* have a non-empty content for that value. To resolve this error, set the appropriate value in the named parameter.

Inconsistent assumptions for ram pressure models
------------------------------------------------

If you receive an error message of:

.. code-block:: text

   only "disk" components are supported by the "simpleCylindrical" ramPressureStripping class

or

.. code-block:: text

   only "spheroid" components are supported by the "simpleSpherical" ramPressureStripping class

you are most likely attempting to use a ram pressure stripping model with a symmetry mismatched to the component to which it is being applied. For example, galactic disks have cylindrical symmetry - if you try to apply a ram pressure stripping model which assumes spherical symmetry to a disk component, you'll get an error like the above.

The recommended approach is to include the relevant ``ramPressureStripping`` parameter inside the corresponding ``nodeOperator``. For example:

.. code-block:: text

   <nodeOperator value="ramPressureMassLossSpheroids">
     <ramPressureStripping value="simpleSpherical">
       <rateFractionalMaximum value="10.0"/>
       <beta value="1"/>
     </ramPressureStripping>
   </nodeOperator>
   <nodeOperator value="ramPressureMassLossDisks">
     <ramPressureStripping value="simpleCylindrical">
       <rateFractionalMaximum value="10.0"/>
       <beta value="1"/>
     </ramPressureStripping>
   </nodeOperator>

This ensures that the ram pressure stripping model of the correct symmetry is found and utilized by the correct ``nodeOperator``.

Unable to find data file
------------------------

An error of the form:

.. code-block:: text

   Fatal error:
   Unable to find data file "./static/foo/bar.xml"

means that Galacticus is looking for a data file (that it needs to read at run time), but is not finding it. Most commonly this is because you do not have the environment variable ``GALACTICUS_DATA_PATH`` set. This environment variable must be set to the path where you installed the `datasets <https://github.com/galacticusorg/datasets>`_ repo. You can do this using, e.g.:

.. code-block:: console

   export GALACTICUS_DATA_PATH=/path/to/where/you/installed/datasets

(putting in the correct path of course).

Temporary parameter files fill ``/dev/shm``
-------------------------------------------

If every run fails immediately after reading the parameter file with:

.. code-block:: text

   Fatal error:
   unable to close file object '/dev/shm/glcTmpPar.145209.1'
    Occurred at:
      subroutine:IO_HDF5_Finalize_Shared
          module:IO_HDF5
            file:utility/IO/HDF5/_module.F90   [line 748]
   HDF5: infinite loop closing library

the shared-memory file system ``/dev/shm`` is full. Each run writes a small temporary HDF5 file, ``/dev/shm/glcTmpPar.<pid>.<n>``, holding the parameters to be recorded in the output. Versions of Galacticus before August 2026 never removed this file (`issue #1391 <https://github.com/galacticusorg/galacticus/issues/1391>`_), so after several hundred runs on one machine, or far fewer inside a Docker container (where ``/dev/shm`` defaults to 64 MB), the file system fills and every later run fails at startup. The failure is caused by *previous* runs, so it can look like a regression introduced by whatever you changed last.

Remove the orphaned files, keeping any that belong to a model still running:

.. code-block:: bash

   for f in /dev/shm/glcTmpPar.*; do
       p=$(basename "$f" | cut -d. -f2)
       [ -d "/proc/$p" ] || rm -f "$f"
   done

Inside Docker, also start the container with a larger ``--shm-size``. Current versions of Galacticus remove the file as soon as it is opened, so the problem does not recur once you update. (On macOS the file is placed in ``/tmp`` instead, where the same leak was less likely to be noticed.)

Unable to allocate memory
-------------------------

If you see:

.. code-block:: text

   unable to allocate `position` (3000000000 elements)
   HELP: the run needs more memory than is available to it. Reduce the size of the problem (for
   example, the number of trees or particles being processed), reduce the number of OpenMP threads,
   since each holds its own copy of much of the state, or run where more memory is available

Galacticus asked for more memory than the system (or your batch job) would give it. The message names the array and its size, which tells you what is driving the requirement: particle arrays when reading N-body data, node arrays when reading large merger tree files, and so on. The remedies are those in the ``HELP`` line, and the discussion of memory limits under `Bus errors, Illegal instructions, and Signal 9`_ applies here too. Not every allocation is guarded in this way, so an out-of-memory condition can also appear as one of those signals rather than as this message.

Recursive object construction
-----------------------------

If you see:

.. code-block:: text

   a [darkMatterProfileDMO] composites a member of its own class but no such [darkMatterProfileDMO]
   was provided explicitly - this would lead to an infinite recursive build; provide a
   [darkMatterProfileDMO] explicitly (see issue 397)

you have selected a *decorator* implementation: a class that wraps another implementation of the same class and modifies its results. Examples are ``darkMatterProfileDMO value="heated"`` (which heats another dark matter profile), ``criticalOverdensity value="renormalize"``, ``cosmologicalMassVariance value="scaled"``, and ``mergerTreeConstructor value="filter"``. A decorator needs to be told which implementation it wraps. If you do not nest one inside it, Galacticus looks for the nearest ``darkMatterProfileDMO`` in the parameter file, finds the decorator itself, and stops rather than build it forever.

Nest the inner implementation inside the decorator:

.. code-block:: xml

   <darkMatterProfileDMO value="heated">
     <darkMatterProfileDMO value="NFW"/>
   </darkMatterProfileDMO>

Some recursive references are legitimate and are handled automatically since `pull request #1264 <https://github.com/galacticusorg/galacticus/pull/1264>`_ (for example a class that needs a pointer back to the object which contains it); this message is raised only where the recursion would be infinite.

Merger tree appears to be deadlocked
------------------------------------

If you see:

.. code-block:: text

   merger tree appears to be deadlocked (see preceding report) - check timestep criteria

no node in a tree can take a step forward, usually because a timestep criterion or a satellite merging time has been set inconsistently. The report printed before the error lists the state of every node. See :doc:`tree-deadlocks` for how to read that report and find the offending node.

Star formation history cannot extend its time grid
--------------------------------------------------

If you see:

.. code-block:: text

   serialized history increment cannot extend the time grid
      subroutine:History_Increment_Serialized
          module:Histories
            file:objects/history.F90

you are recording star formation histories with ``<starFormationHistory value="adaptive"/>`` in a model whose output times are widely separated (for example a set of outputs near redshift 0 together with a set near redshift 2). The adaptive class builds its time grid from the first output and cannot extend it far enough to reach the later ones. This is tracked as `issue #1441 <https://github.com/galacticusorg/galacticus/issues/1441>`_. Until it is resolved, run the widely separated groups of outputs as separate models, each with only the outputs it needs.

Failed to download a file
-------------------------

If you see:

.. code-block:: text

   failed to download from "https://example.org/some/file.tar.gz"

Galacticus needed a file it does not ship (an external tool such as CAMB, FSPS, or RecFast, or a data file) and could not fetch it. The downloader tries ``wget`` and then falls back to ``curl``, retrying each a few times, so first check that at least one of them is installed and that the machine you are running on has outbound network access. Compute nodes on clusters often do not; in that case run the model once on a login node so that the downloaded tools and the files under ``datasets/dynamic/`` are in place, then submit the job. If your site requires a proxy, set the usual ``https_proxy`` environment variable, which both downloaders honor.

Two failure modes are worth knowing about. Some hosts refuse ``wget`` outright, answering every request with a ``404`` because of the way it negotiates TLS, even though the file exists; the fallback to ``curl`` handles this, but only if ``curl`` is installed (`pull request #1376 <https://github.com/galacticusorg/galacticus/pull/1376>`_). And the original RecFast download site is unreliable, so a mirror is tried as well (`pull request #1390 <https://github.com/galacticusorg/galacticus/pull/1390>`_). If a download fails persistently for a reason outside your control, you can fetch the file by other means and place it at the path Galacticus expects, which is given in the log immediately before the error.

Radius specifier evaluates to a radius of zero
----------------------------------------------

If you see:

.. code-block:: text

   Radius specifier evaluates to a radius of zero:

followed by the specifier with the offending part highlighted, a radius given in your parameter file (for example one used to define where a property is measured) has evaluated to zero. Most often the specifier refers to a component, such as a spheroid or a nuclear star cluster, that is absent or empty in the galaxy in question, so its scale radius is zero. Either choose a specifier based on a component that is always present, or use one of the extractors which define a value at zero radius (the ``projectedMass`` extractor, for example, does). The ``HELP:`` line links to the documentation of the radius specifier syntax.
