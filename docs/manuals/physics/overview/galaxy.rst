Galaxy Physics
==============

Below is a flowchart indicating the physical components and processes that typically occur in the evolution of a galaxy in Galacticus. (Galacticus is highly modular - many different components and processes can be included and excluded - this is intended just as a typical example.)

.. mermaid::

   flowchart LR
       subgraph Disk
        direction TB
        DiskISM[ISM]
        DiskStars[Stars]
       end
       subgraph Spheroid
        direction LR
        SpheroidISM[ISM]
        SpheroidStars[Stars]
       end
       subgraph Environment
        direction LR
        CGM[<a href='cgm.html' style='text-decoration: none'>CGM</a>]
        IGM
       end
       CGM -->|<a href='cgm-cooling.html' style='text-decoration: none'>cooling</a>| DiskISM
       DiskISM -->|<a href='star-formation.html' style='text-decoration: none'>star formation</a>| DiskStars
       DiskStars -->|<a href='https://galacticus.readthedocs.io/en/latest/physics/stellarPopulationProperties.html' style='text-decoration: none'>recycling</a>| DiskISM
       SpheroidISM -->|<a href='star-formation.html' style='text-decoration: none'>star formation</a>| SpheroidStars
       SpheroidStars -->|<a href='https://galacticus.readthedocs.io/en/latest/physics/stellarPopulationProperties.html' style='text-decoration: none'>recycling</a>| SpheroidISM
       Disk -->|<a href='https://galacticus.readthedocs.io/en/latest/physics/galacticDynamicsBarInstability.html' style='text-decoration: none'>instability</a>| Spheroid
       DiskISM -->|<a href='outflows.html' style='text-decoration: none'>outflow</a>| Environment
       SpheroidISM -->|<a href='outflows.html' style='text-decoration: none'>outflow</a>| Environment
