Star Formation
==============

Below is a flowchart indicating the ingredients of Galacticus star formation model. (Galacticus is highly modular - many different ingredients can be included and excluded - this is intended just as a typical example.)

.. mermaid::

   flowchart LR
      Disk[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#starFormationRateDisks' style='text-decoration: none'>Disk SFR</a>]
      Spheroid[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#starFormationRateSpheroids' style='text-decoration: none'>Spheroid SFR</a>]
      SFR[["SFR: ψ⭑"]]
      Timescale[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#starFormationTimescale' style='text-decoration: none'>Timescale</a>]
      Surface["<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#starFormationRateSurfaceDensityDisks' style='text-decoration: none'>SFR surface density: Σ⭑</a>"]
      SFR --> Disk
      SFR --> Spheroid
      Timescale -.-> SFR
      Surface -.-> SFR
