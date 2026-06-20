Run-Time Error Messages
=======================

As you are running Galacticus it's quite likely that you'll eventually run into an error. This page lists common types of such run-time errors, along with guidance on what they mean and how to fix them.

For *compile-time* errors (i.e. things that go wrong when you are actually running the compiled executable, ``Galacticus.exe``) you should refer to this `page <https://github.com/galacticusorg/galacticus/wiki/Troubleshooting:-Compile-Time-Error-Messages>`_. For other problems (Galacticus is not performing as you expect, or doing what you want) see the general troubleshooting `page <https://github.com/galacticusorg/galacticus/wiki/Troubleshooting>`_.

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
     click DynamicNotExist href "https://github.com/galacticusorg/galacticus/wiki/Troubleshooting:-Run-Time-Error-Messages#dataset-xyz-does-not-exist-in-datasetsdynamic" "dummy"
     ParameterEmptyValueQ{{"Error message includes the phrase
                       'empty value in parameter [xyz]'?"}}
     ParameterEmptyValueQ --yes--> ParameterEmptyValue
     ParameterEmptyValueQ --no--> RamPressureComponentQ
     ParameterEmptyValue(See here)
     style ParameterEmptyValue fill:#74c7db
     click ParameterEmptyValue href "https://github.com/galacticusorg/galacticus/wiki/Troubleshooting:-Run-Time-Error-Messages%3A-Run-Time-Error-Messages#empty-value-in-parameter" "dummy"
     RamPressureComponentQ{{"Error message includes the phrase
                       'only ＂xyz＂ components are supported by the ＂abc＂ ramPressureStripping class'?"}}
     RamPressureComponentQ --yes--> RamPressureComponent
     RamPressureComponent(See here)
     style RamPressureComponent fill:#74c7db
     click RamPressureComponent href "https://github.com/galacticusorg/galacticus/wiki/Troubleshooting%3A-Run-Time-Error-Messages#inconsistent-assumptions-for-ram-pressure-models" "dummy"
     System{{Error message includes the phrase:}}
     System --Floating point exception--> FPE
     System --Segmentation fault--> SegFault
     System --Bus error--> Bus
     System --Illegal instruction--> Illegal
     System --Command terminated by signal 9--> Signal9
     FPE(See here)
     style FPE fill:#74c7db
     click FPE href "https://github.com/galacticusorg/galacticus/wiki/Troubleshooting:-Run-Time-Error-Messages#floating-point-errors-and-segfaults" "dummy"
     SegFault(See here)
     style SegFault fill:#74c7db
     click SegFault href "https://github.com/galacticusorg/galacticus/wiki/Troubleshooting:-Run-Time-Error-Messages#floating-point-errors-and-segfaults" "dummy"
     Bus(See here)
     style Bus fill:#74c7db
     click Bus href "https://github.com/galacticusorg/galacticus/wiki/Troubleshooting%3A-Run-Time-Error-Messages#bus-errors-illegal-instructions-and-signal-9" "dummy"
     Illegal(See here)
     style Illegal fill:#74c7db
     click Illegal href "https://github.com/galacticusorg/galacticus/wiki/Troubleshooting%3A-Run-Time-Error-Messages#bus-errors-illegal-instructions-and-signal-9" "dummy"
     Signal9(See here)
     style Signal9 fill:#74c7db
     click Signal9 href "https://github.com/galacticusorg/galacticus/wiki/Troubleshooting%3A-Run-Time-Error-Messages#bus-errors-illegal-instructions-and-signal-9" "dummy"

Floating point errors and Segfaults
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A floating point exception or segmentation fault ("segfault") indicate a bug in Galacticus. These should never happen! If they do, please file a bug report using this `form <https://github.com/galacticusorg/galacticus/issues/new?assignees=abensonca&labels=bug&projects=&template=Bug-Report.yml&title=%5BBug%5D%3A+>`_, providing as much information as you can so that we can diagnose and fix the problem.

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
