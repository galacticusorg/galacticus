Output
======

Galacticus writes one :term:`HDF5` file per run. At each output time the `outputter <https://galacticus.readthedocs.io/en/latest/physics/mergerTreeOutputter.html>`_ walks every node of a tree and writes the properties selected by the `property extractors <https://galacticus.readthedocs.io/en/latest/physics/nodePropertyExtractor.html>`_ into the ``Outputs/Output<N>/nodeData`` group, one dataset per property, with the trees concatenated and indexed by ``mergerTreeIndex``, ``mergerTreeStartIndex``, and ``mergerTreeCount``. A `galactic filter <https://galacticus.readthedocs.io/en/latest/physics/galacticFilter.html>`_ can restrict which nodes are written, for example to central galaxies above a stellar mass.

.. mermaid::

   flowchart LR
      Node([Node])
      Extractor[<a href='https://galacticus.readthedocs.io/en/latest/physics/nodePropertyExtractor.html' style='text-decoration: none'>Property Extractors</a>]
      Filter[<a href='https://galacticus.readthedocs.io/en/latest/physics/galacticFilter.html' style='text-decoration: none'>Galactic Filter</a>]
      Times[<a href='https://galacticus.readthedocs.io/en/latest/physics/outputTimes.html' style='text-decoration: none'>Output Times</a>]
      Outputter[<a href='https://galacticus.readthedocs.io/en/latest/physics/mergerTreeOutputter.html' style='text-decoration: none'>Outputter</a>]
      File[(HDF5 file)]
      Node --> Filter
      Filter --> Extractor
      Times --> Outputter
      Extractor --> Outputter
      Outputter --> File

There are more than 160 property extractors, from the basic node indices, masses, and positions through luminosities in any filter, star formation histories, density profiles at chosen radii, and merger histories. The ``multi`` extractor combines any number of them, and each extractor declares what it produces per node: a scalar, a tuple of related scalars, an array, or a list of arrays. Alongside ``nodeData``, the file records every parameter the model read (including defaults) in the ``Parameters`` group, the version and build information, and the results of any on-the-fly analyses. The ``standard`` outputter is the usual choice; alternatives write the full evolver state for postprocessing, or nothing at all when only analyses are wanted.
