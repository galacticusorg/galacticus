Analysis Tools
==============

A collection of tools useful for analyzing Galacticus data:

* `Dendros <https://github.com/galacticusorg/dendros>`_ - The companion analysis and visualization package for Galacticus. Dendros provides automated plotting of Galacticus on-the-fly analyses, MCMC chain diagnostics (e.g. Gelman-Rubin convergence, acceptance rate, posterior predictive checks), maximum-likelihood model extraction, and posterior visualization (kernel density estimates, corner plots). Available on `PyPI <https://pypi.org/project/dendros/>`_:

  .. code-block:: bash

     pip install dendros

* `SubScript <https://github.com/cgannonucm/SubScript>`_ - Utility functions for analyzing subhalo distributions, focusing on the Galacticus output format. The goal of this package is to facilitate quick statistical analysis of subhalo distributions across multiple trees. Written by `Charles Gannon <https://github.com/cgannonucm>`_.

* ``scripts/analysis/treeProcessingTimeFit.py`` - Fits a cost model for merger tree processing time from the per-tree timing data recorded by the ``treeProcessingTimer`` merger tree operator (in the ``metaData/treeTiming`` group of a Galacticus output file). Timing data from several output files may be combined into a single fit. The resulting coefficients can be written to an XML file (consumable by the ``file`` :galacticus-class:`metaTreeProcessingTime` class) and/or written back into a Galacticus HDF5 file. This is useful for building a cost model from multiple runs, or for (re-)fitting after the fact; a single run of the ``treeProcessingTimer`` operator performs the equivalent fit automatically. See :doc:`running` for how the resulting cost model is used to estimate run times. For example, to combine the timing data from several runs into a cost model:

  .. code-block:: bash

     scripts/analysis/treeProcessingTimeFit.py run1.hdf5 run2.hdf5 run3.hdf5 --outputHDF5 costModel.hdf5
