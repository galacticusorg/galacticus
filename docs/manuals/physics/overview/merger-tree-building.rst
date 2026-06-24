Merger Tree Building
====================

Below is a flowchart indicating the ingredients of Galacticus merger tree building algorithm. (Galacticus is highly modular - many different ingredients can be included and excluded - this is intended just as a typical example.)

.. mermaid::

   flowchart LR
      Builder[<a href='https://galacticus.readthedocs.io/en/latest/physics/mergerTreeBuilder.html' style='text-decoration: none'>Builder</a>]
      Masses[<a href='https://galacticus.readthedocs.io/en/latest/physics/mergerTreeBuildMassDistribution.html' style='text-decoration: none'>Mass Distribution</a>]
      Branching[<a href='https://galacticus.readthedocs.io/en/latest/physics/mergerTreeBranchingProbability.html' style='text-decoration: none'>Branching Distribution</a>]
      Excursion[<a href='https://galacticus.readthedocs.io/en/latest/physics/excursionSetFirstCrossing.html' style='text-decoration: none'>Excursion Set Solver</a>]
      Controller[<a href='https://galacticus.readthedocs.io/en/latest/physics/mergerTreeBuildController.html' style='text-decoration: none'>Controller</a>]
      Variance["<a href='https://galacticus.readthedocs.io/en/latest/physics/cosmologicalMassVariance.html' style='text-decoration: none'>σ(M)</a>"]
      Critical["<a href='https://galacticus.readthedocs.io/en/latest/physics/criticalOverdensity.html' style='text-decoration: none'>δ<sub>c</sub>(t)</a>"]
      Cosmology[<a href='https://galacticus.readthedocs.io/en/latest/physics/cosmologyParameters.html' style='text-decoration: none'>Cosmology</a>]
      Resolution[<a href='https://galacticus.readthedocs.io/en/latest/physics/mergerTreeMassResolution.html' style='text-decoration: none'>Resolution</a>]
      Cosmology --> Builder
      Masses --> Builder
      Controller --> Builder
      Branching --> Controller
      Excursion -.-> Branching
      Variance --> Builder
      Critical --> Builder
      Resolution --> Controller
