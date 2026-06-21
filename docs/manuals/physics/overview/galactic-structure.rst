Galactic Structure
==================

Below is a flowchart indicating the ingredients of Galacticus galactic structure model. (Galacticus is highly modular - many different ingredients can be included and excluded - this is intended just as a typical example.)

.. mermaid::

   flowchart LR
      Galaxy[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>Galaxy structure</a>]
      DMO[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>Dark matter-only profile</a>]
      DM[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>Dark matter profile</a>]
      Solver[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>Structure solver</a>]
      DMO --> DM
      DM --> Solver
      Solver --> Galaxy
      Galaxy --> Solver
