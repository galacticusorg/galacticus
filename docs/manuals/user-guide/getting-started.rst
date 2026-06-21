Getting Started
===============

This page is a quick on-ramp for new users: from installing Galacticus to
running your first model and looking at its output.

1. Install Galacticus
---------------------

Choose an installation route in :doc:`installation/index`. The quickest is the
:doc:`pre-compiled binary <installation/binary>`; if you intend to modify the
code, build :doc:`from source <installation/source-linux>` instead.

2. Set up the environment
-------------------------

Galacticus needs two environment variables — ``GALACTICUS_EXEC_PATH`` (the
Galacticus directory) and ``GALACTICUS_DATA_PATH`` (the run-time datasets) — as
described in :doc:`running`. The run-time datasets are downloaded separately; see
the installation instructions for your platform.

3. Run your first model
-----------------------

Galacticus is run by passing it a *parameter file*:

.. code-block:: bash

   Galacticus.exe parameters/quickTest.xml

``quickTest.xml`` is a small, fast model included with the source. A successful
run writes its results to an HDF5 file named ``galacticus.hdf5`` in the working
directory. :doc:`running` walks through this run step by step and shows how to
explore the output.

What is a parameter file?
-------------------------

Everything about a Galacticus model — the physics included, the cosmology, and
what to output — is controlled by an XML *parameter file*. A minimal example
looks like:

.. code-block:: xml

   <parameters>
     <outputFileName value="myModel.hdf5"/>
     <!-- ... further parameters ... -->
   </parameters>

Each element sets one parameter. Many parameters select a *class* implementation
for a piece of physics — for example, which halo mass function or cooling model
to use. The :doc:`advanced` chapter describes parameter files in depth:
conditional parameters, math expressions, sub-parameters, and how the output is
structured. The full list of available parameters for every class is documented
in the :doc:`physics reference <../../physics/index>`.

4. Analyze the output
---------------------

Galacticus output is standard HDF5 and can be read from any language. The
companion `Dendros <https://github.com/galacticusorg/dendros>`_ package provides
ready-made analysis and plotting tools, including on-the-fly analyses, MCMC chain
diagnostics, and posterior corner plots (``pip install dendros``).

Next steps
----------

* :doc:`running` — a detailed walk-through of running a model and reading its
  output.
* :doc:`advanced` — parameter files, output structure, and advanced usage.
* :doc:`troubleshooting/index` — what to do if something goes wrong.
* A collection of
  `tutorials <https://github.com/galacticusorg/galacticus/wiki/Tutorials>`_ is
  available on the wiki.
