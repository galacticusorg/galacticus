Debugging
=========

This page provides guidance on how to debug errors and problems that you may encounter while developing Galacticus (adding new features, modifying existing classes etc.)

Breakdown of an Error Message
-----------------------------

If something goes wrong at run-time, you'll likely see an error message. In this section we'll break down each part of that message to understand what it means.

The Stack Trace
~~~~~~~~~~~~~~~

Galacticus outputs a `stack trace <https://en.wikipedia.org/wiki/Stack_trace>`_ - essentially a list of which functions Galacticus had progressed through to at the time of the error, along with the names of the files in which those functions are called, and the line number in those files. It might look something like this:

.. code-block:: text

   #0 0x7bcd10 in signal_handler_sigsegv
       at ./work/build/error.p.F90:441
   #1 0x7f5c47be63ff in ???
   #2 0x9d283e in __iso_varying_string_MOD_op_assign_vs_vs
       at ./work/build/ISO_Varying_String/iso_varying_string.p.F90:544
   #3 0xf56467 in __node_property_extractors_MOD_newextractorconstructorinternal
       at ./work/build/nodes.property_extractor.my_new_extractor.p.F90:410
   #4 0xf57598 in __node_property_extractors_MOD_newextractconstructorparameters
       at ./work/build/nodes.property_extractor.my_new_extractor.p.F90:228

The numbers at the start of each line (``#0``, ``#1``, etc.) enumerate the depth through the function calls. The top line, ``#0``, is the current function, ``#1`` is the function that called ``#0`` and so on. Usually the first function listed, ``#0``, is just an error handler function (in this case ``signal_handler_sigsegv()`` which handles segfault errors), and often ``#1`` is some unknown function internal to the OS (``???``) so we ignore it.

After that, we have functions of interest in the Galacticus source code. In this example, we were working on a new ``nodePropertyExtractor`` in a file we named ``nodes.property_extractor.my_new_extractor.F90``. This is first mentioned in ``#3``, so that is most likely where the problem occurs.

To understand the stack trace it's important to know that Galacticus first pre-processes any Fortran file from the ``source/`` directory (manipulating the code in a variety of ways), and emits it to a similarly named file in ``work/build/`` with an extension ``.p.F90`` (so, ``source/foo.F90`` is pre-processed to become ``work/build/foo.p.F90``). You'll see that the file names in the stack trace refer to these pre-processed files. The line numbers also refer to lines in the pre-processed file, *not* in your original source file.

So:

.. code-block:: text

   #3 0xf56467 in __node_property_extractors_MOD_newextractorconstructorinternal
       at ./work/build/nodes.property_extractor.my_new_extractor.p.F90:410

means that you should open up the file ``./work/build/nodes.property_extractor.my_new_extractor.p.F90`` and look at line 410 to find the cause of the error.
