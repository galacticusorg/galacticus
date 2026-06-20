Advanced Compilation Options
============================

Debugging Options
-----------------

Thread Analyzing
~~~~~~~~~~~~~~~~~

Before running Galacticus inside a thread analyzer (such as `DRD <https://valgrind.org/docs/manual/drd-manual.html>`_ or `Helgrind <https://valgrind.org/docs/manual/hg-manual.html>`_) you should:

* Configure GCC with the ``--disable-linux-futex`` option and recompile it - this is necessary since those thread analyzers don't correctly handle Linux futexes.

* Compile Galacticus with the ``-DTHREADSAFEIO`` option, which forces I/O operations to be compiled inside OpenMP ``critical`` sections - this is necessary because thread analyzers (`rightly <https://gcc.gnu.org/bugzilla/show_bug.cgi?id=92836>`_ or `wrongly <https://gcc.gnu.org/bugzilla/show_bug.cgi?id=50175>`_) claim that ``gfortran`` I/O operations are not thread safe.

Selecting Components to Build
-----------------------------

In Galacticus, the content of each node of a merger tree is represented by a collection of "components" (e.g. a galactic disk, a supermassive black hole, etc.). Normally, Galacticus is built with support for all defined components. However, this means that (even if you set some component to "null" at run time) there is some memory usage associated with the component in each node. For memory-heavy situations (such as running very high resolution merger trees) this extra memory can be a problem.

Therefore, it is possible to specify a list of "active" components at build time. Support for other, inactive, components is then not built, reducing the memory overhead at runtime. The list of active components is specified as a space-separated list of component names in the ``GALACTICUS_ACTIVE_COMPONENTS`` environment variable. For example:

.. code-block:: bash

   export GALACTICUS_ACTIVE_COMPONENTS="basic darkMatterProfile satellite"

would cause Galacticus to be built with support for only the ``basic``, ``darkMatterProfile``, and ``satellite`` components.

.. note::

   If any physical process attempts to access an unsupported component at run time, Galacticus will exit with an error message.
