Circumgalactic Medium
=====================

Below is a flowchart indicating the physical components and processes that typically occur in the evolution of the circumgalactic medium. (Galacticus is highly modular - many different components and processes can be included and excluded - this is intended just as a typical example.)

.. mermaid::

   flowchart LR
   	IGM
       CGM
   	Outflowed
   	Galaxy
   	IGM -->|<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#accretionHalo' style='text-decoration: none'>accretion</a>| CGM
   	Outflowed -->|<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#hotHaloOutflowReincorporation' style='text-decoration: none'>reincorportation</a>| CGM
   	CGM -->|<a href='https://github.com/galacticusorg/galacticus/wiki/CGM-Cooling-Physics' style='text-decoration: none'>cooling</a>| Galaxy
   	Galaxy -- outflow --> Outflowed
   	Galaxy -- outflow --> IGM
