Intergalactic Medium and Reionization
=====================================

The thermal and ionization history of the :term:`IGM` enters galaxy formation in two places: it sets the ultraviolet background that photo-ionizes and heats gas, and it raises the Jeans mass so that low-mass halos accrete less than the cosmic baryon fraction. The `IGM state <https://galacticus.readthedocs.io/en/latest/physics/intergalacticMediumState.html>`_ class provides the neutral and ionized fractions of hydrogen and helium, the temperature, and the electron-scattering optical depth as functions of time, either from RecFast, from a file, from an instantaneous reionization prescription, or evolved internally by the ``intergalacticMediumStateEvolve`` universe operator from the ionizing emissivity of the model's own galaxies.

.. mermaid::

   flowchart LR
      State[<a href='https://galacticus.readthedocs.io/en/latest/physics/intergalacticMediumState.html' style='text-decoration: none'>IGM State</a>]
      Radiation[<a href='https://galacticus.readthedocs.io/en/latest/physics/radiationField.html' style='text-decoration: none'>Radiation Field</a>]
      Filtering[<a href='https://galacticus.readthedocs.io/en/latest/physics/intergalacticMediumFilteringMass.html' style='text-decoration: none'>Filtering Mass</a>]
      Accretion[<a href='https://galacticus.readthedocs.io/en/latest/physics/accretionHalo.html' style='text-decoration: none'>Halo Accretion</a>]
      Cooling[<a href='https://galacticus.readthedocs.io/en/latest/physics/coolingFunction.html' style='text-decoration: none'>Cooling Function</a>]
      Environment[<a href='https://galacticus.readthedocs.io/en/latest/physics/haloEnvironment.html' style='text-decoration: none'>Environment</a>]
      Galaxy([Galaxy])
      State --> Filtering
      State --> Radiation
      Filtering --> Accretion
      Radiation --> Cooling
      Environment --> Accretion
      Accretion --> Galaxy
      Cooling --> Galaxy

The `filtering mass <https://galacticus.readthedocs.io/en/latest/physics/intergalacticMediumFilteringMass.html>`_ of :cite:t:`gnedin_effect_2000` converts the IGM temperature history into the halo mass below which accretion is suppressed, and the `halo accretion <https://galacticus.readthedocs.io/en/latest/physics/accretionHalo.html>`_ classes apply that suppression, optionally with a dependence on the large-scale `environment <https://galacticus.readthedocs.io/en/latest/physics/haloEnvironment.html>`_ that reionized earlier or later. `Radiation fields <https://galacticus.readthedocs.io/en/latest/physics/radiationField.html>`_ (the cosmic microwave background, tabulated or internally computed ultraviolet backgrounds, and stellar fields) are integrated against atomic cross-sections to give photo-ionization and photo-heating rates for the cooling and chemistry of circumgalactic gas.
