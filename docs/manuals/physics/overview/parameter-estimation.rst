Parameter Estimation
====================

The ``posteriorSample`` task fits model parameters to data. A set of `model parameters <https://galacticus.readthedocs.io/en/latest/physics/modelParameter.html>`_, each with a prior and an optional mapping (for example logarithmic), defines the space. A `simulation <https://galacticus.readthedocs.io/en/latest/physics/posteriorSampleSimulation.html>`_ moves a population of chains through that space: the default is differential evolution MCMC, with annealed, tempered, and stochastic variants, a particle swarm, and a grid for exhaustive scans. At every step each chain writes a parameter file, runs a Galacticus model, and evaluates a `likelihood <https://galacticus.readthedocs.io/en/latest/physics/posteriorSampleLikelihood.html>`_.

.. mermaid::

   flowchart LR
      Parameters[<a href='https://galacticus.readthedocs.io/en/latest/physics/modelParameter.html' style='text-decoration: none'>Model Parameters</a>]
      Initialize[<a href='https://galacticus.readthedocs.io/en/latest/physics/posteriorSampleStateInitialize.html' style='text-decoration: none'>Initial State</a>]
      State[<a href='https://galacticus.readthedocs.io/en/latest/physics/posteriorSampleState.html' style='text-decoration: none'>Chain State</a>]
      Simulation[<a href='https://galacticus.readthedocs.io/en/latest/physics/posteriorSampleSimulation.html' style='text-decoration: none'>Simulation</a>]
      Model([Galacticus model])
      Analyses[<a href='https://galacticus.readthedocs.io/en/latest/physics/outputAnalysis.html' style='text-decoration: none'>Analyses</a>]
      Likelihood[<a href='https://galacticus.readthedocs.io/en/latest/physics/posteriorSampleLikelihood.html' style='text-decoration: none'>Likelihood</a>]
      Convergence[<a href='https://galacticus.readthedocs.io/en/latest/physics/posteriorSampleConvergence.html' style='text-decoration: none'>Convergence</a>]
      Parameters --> Initialize
      Initialize --> State
      State --> Simulation
      Simulation --> Model
      Model --> Analyses
      Analyses --> Likelihood
      Likelihood --> Simulation
      Simulation --> Convergence

The ``galaxyPopulation`` likelihood runs a full model and sums the log-likelihoods of its on-the-fly analyses, so a constraint is assembled from the same analysis objects used for output. Other likelihoods fit halo mass functions, spin distributions, projected correlation functions, or spectral energy distributions, or combine independent likelihoods. Chains start from the prior, a Latin hypercube, a Gaussian around a known maximum, or a previous run's state file, and `convergence <https://galacticus.readthedocs.io/en/latest/physics/posteriorSampleConvergence.html>`_ is judged by the Gelman-Rubin statistic by default. Chains run in parallel under MPI, one model per process, and every proposed and accepted state is logged for later analysis.
