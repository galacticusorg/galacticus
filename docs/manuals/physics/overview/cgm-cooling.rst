CGM Cooling
===========

Below is a flowchart indicating the ingredients of Galacticus CGM cooling model. (Galacticus is highly modular - many different ingredients can be included and excluded - this is intended just as a typical example.)

.. mermaid::

   flowchart LR
      Galaxy
      Cold[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>Cold Mode Inflow</a>]
      Lambda[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>Cooling Function</a>]
      Radius[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>Cooling Radius</a>]
      Rate[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>Cooling Rate</a>]
      Time[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>Cooling Time</a>]
      RadiusFreefall[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>Freefall Radius</a>]
      Available[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>Time Available Cooling]
      AvailableFreefall[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>Time Available Freefall]
      Infall[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>Infall Radius</a>]
      Angular[<a href='https://galacticus.readthedocs.io/en/latest/physics/index.html' style='text-decoration: none'>Angular Momentum Content</a>]
      Lambda --> Time
      Time --> Radius
      Available --> Radius
      Radius --> Infall
      AvailableFreefall --> RadiusFreefall
      RadiusFreefall --> Infall
      Infall --> Rate
      Rate --> Galaxy
      Cold --> Galaxy
      Angular --> Galaxy
