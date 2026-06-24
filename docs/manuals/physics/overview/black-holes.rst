Supermassive Black Holes
========================

Below is a flowchart indicating the physical components and processes that typically occur in the evolution of supermassive black holes. (Galacticus is highly modular - many different components and processes can be included and excluded - this is intended just as a typical example.)

.. mermaid::

   flowchart LR
       IGM
       CGM[<a href='https://github.com/galacticusorg/galacticus/wiki/CGM-Physics' style='text-decoration: none'>CGM</a>]
       Spheroid
       Disk[Accretion disk]
       SMBH
       CGM -- accretion --> SMBH
       Spheroid -- accretion --> SMBH
       Disk -->|<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#nodeOperatorBlackHolesWinds' style='text-decoration: none'>energy input</a>| Spheroid -- outflow --> CGM
       Disk -->|<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html#blackHoleCGMHeating' style='text-decoration: none'>energy input</a>| CGM -- outflow --> IGM
