About Galacticus
================

At its core, Galacticus is a semi-analytic model of galaxy formation. It solves equations describing how galaxies evolve in a merging hierarchy of dark matter halos in a dark matter-dominated universe. Galacticus has much in common with other semi-analytic models, such as the range of physical processes included and the type of quantities that it can predict. Galacticus also provides a wide variety of other functionality, such as computing halo mass functions, power spectra, analyzing particle simulations, and performing :term:`MCMC` simulations.

Analysis and visualization of Galacticus outputs (including on-the-fly analysis plots, :term:`MCMC` chain diagnostics, and posterior corner plots) is provided by the companion `Dendros <https://github.com/galacticusorg/dendros>`_ package, which is available on `PyPI <https://pypi.org/project/dendros/>`_.

In designing Galacticus our main goal was to make the code flexible, modular and easily extensible. Much greater priority was placed on making the code easy to use and modify than on making it fast. We believe that a modular and extensible nature is crucial as galaxy formation is an evolving science. In particular, key design features are:

Extensible implementations for all functions:
   Essentially all functions within Galacticus are designed to be extensible following an Object Oriented methodology, meaning that you can write your own version and insert it into Galacticus easily. For example, suppose you want to use an improved functional form for the :term:`CDM` halo mass function. You would simply write a new ``haloMassFunction`` class that computes this mass function, decorate it with a short directive (see Section :galacticus-ref:`CodeDirectives`) which explains to the build system how to insert this class into Galacticus. A recompile of the code will then incorporate your new function.

Extensible components for tree nodes:
   The basic structure in Galacticus is a merger tree, which consists of a linked tree of nodes (each corresponding a dark matter halo and its content) which have various properties. Galacticus works by evolving the nodes forward in time subject to a collection of differential equations and other rules. Each node can contain an arbitrary number of *components*. A component may be a dark matter halo, a galactic disk, a black hole etc. Each component may have an arbitrary number of *properties* (some of which may be evolving, others of which can be fixed). Galacticus makes it easy to add additional components. For example, suppose you wanted to add a "stellar halo" components (consisting of stars stripped from satellite galaxies). To do this, you would write a module which specifies the following for this component:

   * Properties (their names, types, and ranks);
   * Functions describing the differential equations which govern the evolution of the properties;
   * Functions describing how the component responds to various events (e.g. the node becoming a satellite, a galaxy-galaxy merger, etc.);
   * "Pipes" which allow for flows of mass/energy/etc. from one component to another.

   Short directives embedded in this module explain to the Galacticus build system how to incorporate the new component. A recompile will then build your new component into Galacticus. Typically, a new component can be created quickly by copying an existing one and modifying it as necessary. Furthermore, multiple implementations of a component are allowed. For example, Galacticus contains a component which tracks the scale length of the dark matter halo. You could add a new component which additionally tracks the axis ratios of the (now triaxial) halo. A simple input parameter then allows you to select which implementation will be used in a given run.

Centralized ODE solver:
   Galacticus evolves nodes in merger trees by calling an ODE solver which integrates forward in time to solve for the evolution of the properties of each component in a node. This means that you do not need to provide explicit solutions for ODEs (in many cases such solutions are not available anyway) and timestepping is automatically handled to achieve a specified level of precision. The ODE solver allows for the evolution to be interrupted. A component may trigger an interrupt at any time and may do so for a number of reasons. A typical use is to actually create a component within a given node---for example when gas first begins to cool and inflow in a node the disk component must be created. Other uses include interrupting evolution when a merging event occurs.

Getting Galacticus
------------------

Galacticus is available in many different forms: a precompiled binary, a `Docker <https://www.docker.com/>`_ image, and the full source code. See :doc:`installation/index` for downloads and installation instructions.

License
-------

Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022 Andrew Benson `<abenson@carnegiescience.edu> <mailto:abenson@carnegiescience.edu>`_

Galacticus is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Galacticus is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Galacticus.  If not, see `<http://www.gnu.org/licenses/> <http://www.gnu.org/licenses/>`_.
