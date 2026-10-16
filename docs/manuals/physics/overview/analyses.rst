On-the-fly Analyses
===================

Rather than writing every galaxy and computing statistics afterward, Galacticus can accumulate an observable while it runs. An `analysis <https://galacticus.readthedocs.io/en/latest/physics/outputAnalysis.html>`_ receives each node at each output time, extracts a property, transforms it to the observed frame, weights it, bins it, and, once all trees on all processes are done, normalizes the result and writes it to the output file together with its covariance. Because each analysis also returns a log-likelihood against the observational data it was built for, the same objects drive Bayesian parameter estimation.

.. mermaid::

   flowchart LR
      Node([Node])
      Filter[<a href='https://galacticus.readthedocs.io/en/latest/physics/galacticFilter.html' style='text-decoration: none'>Galactic Filter</a>]
      Extractor[<a href='https://galacticus.readthedocs.io/en/latest/physics/nodePropertyExtractor.html' style='text-decoration: none'>Property Extractor</a>]
      PropertyOp[<a href='https://galacticus.readthedocs.io/en/latest/physics/outputAnalysisPropertyOperator.html' style='text-decoration: none'>Property Operator</a>]
      WeightOp[<a href='https://galacticus.readthedocs.io/en/latest/physics/outputAnalysisWeightOperator.html' style='text-decoration: none'>Weight Operator</a>]
      Survey[<a href='https://galacticus.readthedocs.io/en/latest/physics/surveyGeometry.html' style='text-decoration: none'>Survey Geometry</a>]
      DistOp[<a href='https://galacticus.readthedocs.io/en/latest/physics/outputAnalysisDistributionOperator.html' style='text-decoration: none'>Distribution Operator</a>]
      Normalizer[<a href='https://galacticus.readthedocs.io/en/latest/physics/outputAnalysisDistributionNormalizer.html' style='text-decoration: none'>Normalizer</a>]
      Analysis[<a href='https://galacticus.readthedocs.io/en/latest/physics/outputAnalysis.html' style='text-decoration: none'>Analysis</a>]
      Likelihood[<a href='https://galacticus.readthedocs.io/en/latest/physics/posteriorSampleLikelihood.html' style='text-decoration: none'>Likelihood</a>]
      Node --> Filter
      Filter --> Extractor
      Extractor --> PropertyOp
      PropertyOp --> DistOp
      Survey --> WeightOp
      WeightOp --> Analysis
      DistOp --> Analysis
      Analysis --> Normalizer
      Normalizer --> Likelihood

`Property operators <https://galacticus.readthedocs.io/en/latest/physics/outputAnalysisPropertyOperator.html>`_ take logarithms, convert to magnitudes or to observed-frame quantities using the cosmology, and add systematic offsets. `Weight operators <https://galacticus.readthedocs.io/en/latest/physics/outputAnalysisWeightOperator.html>`_ apply :math:`1/V_\mathrm{max}` volumes and completeness from a `survey geometry <https://galacticus.readthedocs.io/en/latest/physics/surveyGeometry.html>`_, of which about twenty are built in (SDSS, GAMA, PRIMUS, ULTRAVISTA, ZFOURGE, ALFALFA, and Local Group census footprints among them). `Distribution operators <https://galacticus.readthedocs.io/en/latest/physics/outputAnalysisDistributionOperator.html>`_ convolve the binned result with measurement errors, incompleteness, or lensing magnification, and `normalizers <https://galacticus.readthedocs.io/en/latest/physics/outputAnalysisDistributionNormalizer.html>`_ turn counts into number densities per dex or per unit volume. Around seventy analyses exist: stellar mass and luminosity functions at many redshifts, size-mass and black hole scaling relations, HI mass functions, correlation functions, subhalo mass functions and radial distributions, and Local Group satellite statistics.
