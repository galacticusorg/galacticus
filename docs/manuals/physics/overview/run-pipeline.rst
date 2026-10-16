Anatomy of a Run
================

Everything Galacticus does is a `task <https://galacticus.readthedocs.io/en/latest/physics/task.html>`_. The default task, ``evolveForests``, is the galaxy formation pipeline: build or read merger trees, apply any tree operators, evolve every tree with the engine, write the requested outputs, and finalize the on-the-fly analyses. Other tasks build external tools and tabulations, analyze N-body simulations, run Bayesian parameter estimation, or perform radiative transfer, and several tasks can be chained with the ``multi`` task. Below is the order of operations for the default task; each box links to the class that controls it.

.. mermaid::

   flowchart LR
      Task[<a href='https://galacticus.readthedocs.io/en/latest/physics/task.html' style='text-decoration: none'>Task</a>]
      Constructor[<a href='https://galacticus.readthedocs.io/en/latest/physics/mergerTreeConstructor.html' style='text-decoration: none'>Tree Constructor</a>]
      Filter[<a href='https://galacticus.readthedocs.io/en/latest/physics/mergerTreeFilter.html' style='text-decoration: none'>Tree Filter</a>]
      Operators[<a href='https://galacticus.readthedocs.io/en/latest/physics/mergerTreeOperator.html' style='text-decoration: none'>Tree Operators</a>]
      Evolver[<a href='https://galacticus.readthedocs.io/en/latest/physics/mergerTreeEvolver.html' style='text-decoration: none'>Evolver</a>]
      Times[<a href='https://galacticus.readthedocs.io/en/latest/physics/outputTimes.html' style='text-decoration: none'>Output Times</a>]
      Outputter[<a href='https://galacticus.readthedocs.io/en/latest/physics/mergerTreeOutputter.html' style='text-decoration: none'>Outputter</a>]
      Analyses[<a href='https://galacticus.readthedocs.io/en/latest/physics/outputAnalysis.html' style='text-decoration: none'>Analyses</a>]
      Task --> Constructor
      Constructor --> Filter
      Filter --> Operators
      Operators --> Evolver
      Times --> Evolver
      Evolver --> Outputter
      Evolver --> Analyses
      Operators --> Analyses

Trees that fail the `filter <https://galacticus.readthedocs.io/en/latest/physics/mergerTreeFilter.html>`_ are skipped entirely, which is how a run can be restricted to trees whose base node passes a `galactic filter <https://galacticus.readthedocs.io/en/latest/physics/galacticFilter.html>`_ (a mass range, say), or to a chosen set of tree indices. `Tree operators <https://galacticus.readthedocs.io/en/latest/physics/mergerTreeOperator.html>`_ act on whole trees at four points: before construction, before initialization, before evolution, and after evolution. Examples are pruning branches below a mass, perturbing masses, augmenting trees with unresolved substructure, exporting trees to other formats, and computing tree-level statistics such as conditional mass functions. The `output times <https://galacticus.readthedocs.io/en/latest/physics/outputTimes.html>`_ set both when properties are recorded and where the evolver must stop to record them.
