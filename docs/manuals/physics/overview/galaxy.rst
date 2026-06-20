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
        CGM[<a href='https://github.com/galacticusorg/galacticus/wiki/CGM-Physics' style='text-decoration: none'>CGM</a>]
        IGM
       end
       CGM -->|<a href='https://github.com/galacticusorg/galacticus/wiki/CGM-Cooling-Physics' style='text-decoration: none'>cooling</a>| DiskISM
       DiskISM -->|<a href='https://github.com/galacticusorg/galacticus/wiki/Star-Formation-Physics' style='text-decoration: none'>star formation</a>| DiskStars
       DiskStars -->|<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>recycling</a>| DiskISM
       SpheroidISM -->|<a href='https://github.com/galacticusorg/galacticus/wiki/Star-Formation-Physics' style='text-decoration: none'>star formation</a>| SpheroidStars
       SpheroidStars -->|<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>recycling</a>| SpheroidISM
       Disk -->|<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>instability</a>| Spheroid
       DiskISM -->|<a href='https://github.com/galacticusorg/galacticus/wiki/Outflow-Physics' style='text-decoration: none'>outflow</a>| Environment
       SpheroidISM -->|<a href='https://github.com/galacticusorg/galacticus/wiki/Outflow-Physics' style='text-decoration: none'>outflow</a>| Environment
