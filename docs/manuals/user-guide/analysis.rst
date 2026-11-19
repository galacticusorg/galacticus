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

* ``scripts/analysis/mergerTreeVariance.py`` - Decomposes the tree-to-tree variance of counts of progenitors and subhalos in merger trees into contributions from individual splits, :math:`\mathrm{Var}(N)=\langle\sum\Delta^2\rangle`. It reads trees exported (before evolution) by the ``export`` merger tree operator, measures the expected count of a single branch, :math:`n_1(M,\omega)`, from the trees themselves, and checks the decomposition and its supporting identities for three count statistics. The :galacticus-class:`mergerTreeOperatorBranchingVariance` merger tree operator makes the same measurement on trees in memory, without exporting them, in two runs over a set of trees: a ``table`` run accumulates the :math:`n_1` tables, and a ``jumps`` run evaluates the jumps at every split using them. The script is the reference implementation of the operator, and reproduces it on the same trees when given the operator's settings: the same statistics, thresholds, target redshift, and bin widths, ``--tableMethod linear`` and ``--crossProduct`` (the operator's estimators), ``--foldOffset`` set to the base random number seed of the run (so that trees are assigned to the same folds, by seed), and ``--cumulativeOmegaTable`` set to the output of the ``jumps`` run (from which the script takes the operator's relation between redshift and :math:`\omega`, and its grid for the cumulative variance). ``scripts/analysis/mergerTreeVariancePlot.py`` makes diagnostic plots of the script's output. For example, for trees exported by the ``table`` run of ``testSuite/parameters/mergerTreeBranchingVarianceTable.xml`` (whose base seed is 220), with the output of the matching ``jumps`` run in ``jumps.hdf5``:

  .. code-block:: bash

     python3 scripts/analysis/mergerTreeVariance.py trees.hdf5 --statistics ABC --massThresholds 1.0e10,3.0e10 --redshiftTarget 1.0 \
         --binWidthMass 0.0125 --binWidthOmega 0.025 --tableMethod linear --crossProduct --foldOffset 220 \
         --cumulativeOmegaTable jumps.hdf5 --outputFileName variance.hdf5
     python3 scripts/analysis/mergerTreeVariancePlot.py variance.hdf5 --outputFileName variance.pdf
