N-body Simulation Analysis
==========================

Galacticus can act on N-body simulation output directly, without any galaxy formation. The ``NBodyAnalyze`` task reads particle and halo data through an `importer <https://galacticus.readthedocs.io/en/latest/physics/nbodyImporter.html>`_ (GADGET binary and HDF5, Rockstar, IRATE, Millennium CSV, or a random realization for testing), then passes the data through a chain of `operators <https://galacticus.readthedocs.io/en/latest/physics/nbodyOperator.html>`_, each of which adds, transforms, filters, or writes properties. Results are written back into the HDF5 file when the source is HDF5, or exported.

.. mermaid::

   flowchart LR
      Importer[<a href='https://galacticus.readthedocs.io/en/latest/physics/nbodyImporter.html' style='text-decoration: none'>Importer</a>]
      Filter[Filter operators]
      Compute[Property operators]
      Statistics[Statistics operators]
      Export[Export operators]
      Data([nBodyData])
      Importer --> Data
      Data --> Filter
      Filter --> Compute
      Compute --> Statistics
      Statistics --> Export
      click Filter href "https://galacticus.readthedocs.io/en/latest/physics/nbodyOperator.html"
      click Compute href "https://galacticus.readthedocs.io/en/latest/physics/nbodyOperator.html"
      click Statistics href "https://galacticus.readthedocs.io/en/latest/physics/nbodyOperator.html"
      click Export href "https://galacticus.readthedocs.io/en/latest/physics/nbodyOperator.html"

More than fifty operators exist. Filters select particles or halos by box, sphere, ID, property range, or convex hull, and remove contaminated regions of zoom simulations. Property operators compute self-bound masses (with a Barnes-Hut tree), angular momenta and spins, energy tensors and axis ratios, concentrations, formation and last-major-merger times, environmental overdensities, and convex-hull volumes. Statistics operators build mass functions, subhalo mass and radius functions, concentration distributions, and their covariances, and the results feed the same likelihood classes used for galaxy constraints, so an N-body halo mass function can constrain a cosmological or dark matter parameter in exactly the way a galaxy stellar mass function constrains a feedback parameter.
