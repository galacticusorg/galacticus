Profiling Galacticus Runs
=========================

If Galacticus is running too slowly (particularly if you modified the code and find that it now runs more slowly than before) the best thing to do is to profile the code to see exactly where it is slowing down. There are different ways to do this, with varying degrees of information returned. The following is a brief guide to different profiling approaches.

Simple profiling using ``gdb``
------------------------------

This is a crude, but sometimes effective approach. Simply run your model in the ``gdb`` debugger:

.. code-block:: bash

   gdb --args ./Galacticus.exe myParameters.xml

then type ``run`` to start the model. Every few seconds use ``ctrl+c`` to interrupt the model, and then use ``info threads`` to see which function each thread is currently in. You can also use ``thread X`` followed by ``where`` to switch to a specific numbered thread ``X`` and see the entire stack trace of where that thread is in the code. Then use ``cont`` to continue the model, and repeat.

Sometimes you'll find that there are one or two functions that most of the threads are always in - these are the slow functions, so you can then try to figure out how to speed them up. This is obviously a very crude approach - but it's easy to try and sometimes gives useful insights.

Profiling using ``gprof``
-------------------------

For this, you need to recompile your code with some special options. Set these by doing:

.. code-block:: bash

   export GALACTICUS_FCFLAGS="$GALACTICUS_FCFLAGS -pg"
   export GALACTICUS_CFLAGS="$GALACTICUS_CFLAGS -pg"
   export GALACTICUS_CPPFLAGS="$GALACTICUS_CPPFLAGS -pg"

and then do a full recompile. The ``-pg`` option makes the compiler add extra instructions to track how long each function is used. Then, run the model as usual. Once it's finished you'll get an extra output file called ``gmon.out``. You can then do:

.. code-block:: bash

   gprof ./Galaticus.exe gmon.out > profile.log

which will create a file ``profile.log`` containing profiling information. You can find information `here <https://ftp.gnu.org/old-gnu/Manuals/gprof-2.9.1/html_chapter/gprof_5.html>`_ on how to interpret that output.

Using ``valgrind``
------------------

Compile as normal (i.e. not with the ``-pg`` option from above). Then run using:

.. code-block:: bash

   valgrind --tool=callgrind ./Galacticus.exe myParameters.xml

This will run the model, but it will be slow (often 20-30 times slower than normal). Once it finishes (or if you get bored of waiting and kill it with ``ctrl+c``) you'll get a file ``callgrind.out.<pid>`` where ``<pid>`` is the process ID of the job. The best way to view the output is using ``kcachegrind`` on Linux or ``qcachegrind`` on Mac OS which you can install using:

.. code-block:: bash

   sudo apt install kcachegrind

on Ubuntu or similar Linux distros, or

.. code-block:: bash

   brew install qcachegrind

on Mac OS.

Open the ``callgrind.out.<pid>`` file in ``kcachegrind``/``qcachegrind`` and you'll get a view of function calls, how many times they were called, what fraction of the total time was spent in each one, etc. More details on how to interpret the information can be found `here <https://kcachegrind.github.io/html/Documentation.html>`_.
