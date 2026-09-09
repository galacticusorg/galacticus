Alternative Dark Matter
=======================

Cold dark matter is the default, but Galacticus can follow warm, fuzzy, self-interacting, and decaying dark matter through the whole pipeline. The `dark matter particle <https://galacticus.readthedocs.io/en/latest/physics/darkMatterParticle.html>`_ class declares the candidate and its properties (a thermal relic mass, an axion mass, a self-interaction cross-section, a decay lifetime and kick velocity). The `transfer function <https://galacticus.readthedocs.io/en/latest/physics/transferFunction.html>`_ then carries the small-scale cutoff (Bode 2001 for thermal warm dark matter, the ETHOS fitting form, Murgia 2017 or Passaglia 2022 for fuzzy dark matter, or a direct axionCAMB computation), and the `critical overdensity <https://galacticus.readthedocs.io/en/latest/physics/criticalOverdensity.html>`_ and `merger tree branching <https://galacticus.readthedocs.io/en/latest/physics/mergerTreeBranchingProbability.html>`_ classes have implementations that account for the suppressed collapse of small halos so that tree building reproduces the reduced abundance of low-mass progenitors.

.. mermaid::

   flowchart LR
      Particle[<a href='https://galacticus.readthedocs.io/en/latest/physics/darkMatterParticle.html' style='text-decoration: none'>Dark Matter Particle</a>]
      Transfer[<a href='https://galacticus.readthedocs.io/en/latest/physics/transferFunction.html' style='text-decoration: none'>Transfer Function</a>]
      Critical[<a href='https://galacticus.readthedocs.io/en/latest/physics/criticalOverdensity.html' style='text-decoration: none'>Critical Overdensity</a>]
      Branching[<a href='https://galacticus.readthedocs.io/en/latest/physics/mergerTreeBranchingProbability.html' style='text-decoration: none'>Branching Probability</a>]
      Profile[<a href='https://galacticus.readthedocs.io/en/latest/physics/darkMatterProfileDMO.html' style='text-decoration: none'>Halo Profile</a>]
      Heating[<a href='https://galacticus.readthedocs.io/en/latest/physics/darkMatterProfileHeating.html' style='text-decoration: none'>Profile Heating</a>]
      Trees([Merger trees])
      Subhalos([Subhalo evolution])
      Particle --> Transfer
      Transfer --> Critical
      Critical --> Branching
      Branching --> Trees
      Particle --> Profile
      Heating --> Profile
      Profile --> Subhalos

Halo structure follows: the `dark matter-only profile <https://galacticus.readthedocs.io/en/latest/physics/darkMatterProfileDMO.html>`_ class provides soliton-cored profiles for fuzzy dark matter, isothermal and parametric cores for self-interacting dark matter, profiles depleted by decays, and profiles `heated <https://galacticus.readthedocs.io/en/latest/physics/darkMatterProfileHeating.html>`_ by tidal shocks, impulsive outflows, decay kicks, or two-body relaxation. Because subhalo evolution (see :doc:`subhalo-evolution`) reads the profile through the same interface, tidal stripping and dynamical friction respond to the modified structure without further changes. Dedicated node operators handle the model-specific physics that has no cold dark matter analog, such as gravothermal core collapse in self-interacting models and mass loss through decays.
