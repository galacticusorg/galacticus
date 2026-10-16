Stellar Populations and Dust
============================

Stars formed at each step are treated as a simple stellar population: coeval, with a single metallicity, and distributed in mass according to an `initial mass function <https://galacticus.readthedocs.io/en/latest/physics/initialMassFunction.html>`_. Integrating the `stellar astrophysics <https://galacticus.readthedocs.io/en/latest/physics/stellarAstrophysics.html>`_ tables of lifetime, ejected mass, and yield over that mass function gives the `population <https://galacticus.readthedocs.io/en/latest/physics/stellarPopulation.html>`_'s recycled fraction, metal yield, and energy input as functions of age. `Type Ia supernovae <https://galacticus.readthedocs.io/en/latest/physics/supernovaeTypeIa.html>`_ add iron-peak elements on a delay-time distribution. Whether these returns are applied instantaneously or spread over stellar lifetimes is the choice of `stellar population properties <https://galacticus.readthedocs.io/en/latest/physics/stellarPopulationProperties.html>`_ class.

.. mermaid::

   flowchart LR
      IMF[<a href='https://galacticus.readthedocs.io/en/latest/physics/initialMassFunction.html' style='text-decoration: none'>Initial Mass Function</a>]
      Astro[<a href='https://galacticus.readthedocs.io/en/latest/physics/stellarAstrophysics.html' style='text-decoration: none'>Stellar Astrophysics</a>]
      SNIa[<a href='https://galacticus.readthedocs.io/en/latest/physics/supernovaeTypeIa.html' style='text-decoration: none'>Type Ia Supernovae</a>]
      Population[<a href='https://galacticus.readthedocs.io/en/latest/physics/stellarPopulation.html' style='text-decoration: none'>Stellar Population</a>]
      Properties[<a href='https://galacticus.readthedocs.io/en/latest/physics/stellarPopulationProperties.html' style='text-decoration: none'>Recycling and Yields</a>]
      Feedback[<a href='https://galacticus.readthedocs.io/en/latest/physics/stellarFeedback.html' style='text-decoration: none'>Energy Input</a>]
      Spectra[<a href='https://galacticus.readthedocs.io/en/latest/physics/stellarPopulationSpectra.html' style='text-decoration: none'>Spectra</a>]
      Postprocess[<a href='https://galacticus.readthedocs.io/en/latest/physics/stellarPopulationSpectraPostprocessor.html' style='text-decoration: none'>IGM and Birth Clouds</a>]
      Dust[<a href='https://galacticus.readthedocs.io/en/latest/physics/dustAttenuation.html' style='text-decoration: none'>Dust Attenuation</a>]
      Curve[<a href='https://galacticus.readthedocs.io/en/latest/physics/dustExtinctionCurve.html' style='text-decoration: none'>Extinction Curve</a>]
      Luminosity([Luminosities])
      IMF --> Population
      Astro --> Population
      SNIa --> Population
      Population --> Properties
      Population --> Feedback
      IMF --> Spectra
      Spectra --> Postprocess
      Curve --> Dust
      Postprocess --> Dust
      Dust --> Luminosity

Luminosities come from `spectra <https://galacticus.readthedocs.io/en/latest/physics/stellarPopulationSpectra.html>`_ tabulated by FSPS (or read from file) as a function of age and metallicity, convolved with the star formation history and integrated through filter response curves. `Postprocessors <https://galacticus.readthedocs.io/en/latest/physics/stellarPopulationSpectraPostprocessor.html>`_ multiply in effects between the stars and the observer that do not depend on the galaxy's own dust: absorption by the intergalactic medium (:galacticus-class:`stellarPopulationSpectraPostprocessorInoue2014`, :galacticus-class:`stellarPopulationSpectraPostprocessorMadau1995`, :galacticus-class:`stellarPopulationSpectraPostprocessorMeiksin2006`), Lyman continuum suppression, and birth-cloud attenuation of young populations. The galaxy's own dust is applied last by a `dust attenuation <https://galacticus.readthedocs.io/en/latest/physics/dustAttenuation.html>`_ model, from a simple screen scaled with the metal surface density to the inclination-dependent radiative transfer atlas of :cite:t:`ferrara_atlas_1999`, using a chosen `extinction curve <https://galacticus.readthedocs.io/en/latest/physics/dustExtinctionCurve.html>`_ (:galacticus-class:`dustExtinctionCurveCalzetti2000`, :galacticus-class:`dustExtinctionCurveCardelli1989`, and others). The same dust model is available inside analyses through a property operator, so an observed luminosity function can be attenuated consistently with the written luminosities.
