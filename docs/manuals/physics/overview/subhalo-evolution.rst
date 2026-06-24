Subhalo Evolution
=================

Below is a flowchart indicating the ingredients of Galacticus subhalo evolution physics model. (Galacticus is highly modular - many different ingredients can be included and excluded - this is intended just as a typical example.)

.. mermaid::

   flowchart LR
      Subhalo
      Host[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#massDistribution' style='text-decoration: none'>Host potential</a>]
      Friction[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#satelliteDynamicalFriction' style='text-decoration: none'>Dynamical friction</a>]
      Heat[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#satelliteTidalHeatingRate' style='text-decoration: none'>Tidal heating</a>]
      Strip[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#satelliteTidalStripping' style='text-decoration: none'>Tidal stripping</a>]
      Host -->|<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#virialOrbit' style='text-decoration: none'>orbit</a>| Subhalo
      Friction -->|<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#satelliteDynamicalFriction' style='text-decoration: none'>orbit decay</a>| Subhalo
      Heat -->|<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#darkMatterProfileHeating' style='text-decoration: none'>expansion</a>| Subhalo
      Strip -->|<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#satelliteTidalStripping' style='text-decoration: none'>mass loss</a>| Subhalo
