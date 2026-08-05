.. _manual-sec-Components:

Node Components
===============

A :term:`node` in a merger tree is built from a set of *components*, each of which
holds the state associated with one part of the physical system---the galactic
disk, the circumgalactic medium, the central black hole, the halo's orbit within
its host, and so on. Components are the model's state variables: they define
which properties exist, what their types are, how they are output, and (for
evolvable properties) how their absolute tolerance scales are set.

Each component *class* (``disk``, ``hotHalo``, ``spin``, …) may have several
*implementations*, differing in how much detail they track. Exactly one
implementation of each class is active in any given model, selected by an input
parameter named after the class:

.. code-block:: xml

   <componentDisk value="standard"/>

Parameters belonging to a component implementation are given inside that
element, forming a sub-parameter namespace:

.. code-block:: xml

   <componentDisk value="standard">
     <radiusStructureSolver value="1.0"/>
   </componentDisk>

In addition to the implementations listed below, every component class has a
``null`` implementation. Selecting this---which has no properties and does not
respond to any events---effectively switches off the relevant component class.
This is safe only if no other active implementation expects to get or set
properties of that class. Where a class has no implementation marked as the
default, ``null`` is used unless a choice is made explicitly.

Some implementations *extend* another implementation of the same class, adding
properties to it and inheriting the rest. For example, ``hotHalo``
``coldMode`` extends ``hotHalo`` ``standard``, adding separate cold-mode mass,
abundance, and angular momentum reservoirs.

.. _manual-sec-ComponentsStateAndEvolution:

State and Evolution
-------------------

Components store state; they do not, in general, evolve it. Almost all of the
physics of galaxy formation---star formation, stellar and AGN feedback,
gas accretion and cooling, bar instabilities, ram pressure and tidal stripping,
satellite orbital integration, halo angular momentum and profile evolution---is
implemented by :galacticus-class:`nodeOperatorClass` objects, which contribute
terms to the rates of change of component properties during differential
evolution, and which respond to events such as node mergers and promotions.

This division matters when reading a model's parameter file, or when looking for
where a piece of physics lives: with few exceptions, a component implementation
supplies *what is tracked*, and the selected set of node operators supplies *how
it changes*. No component implementation currently defines a differential
evolution task of its own. Documentation for the individual node operators, and
for the pluggable physics classes they use, is generated from the source and can
be found under the physics class reference.

The tasks components *do* retain are: initialization and thread
initialization, setting ODE absolute tolerance scales, pre-evolution and
post-timestep fix-ups (such as zeroing small negative masses), declaring which
properties are inactive, plausibility checks, participation in the galactic
structure solver, and state store/retrieve for checkpointing. A few also attach
to events---most commonly ``satelliteMerger``, to move mass between components
when galaxies merge, as dictated by the selected
:galacticus-class:`mergerMassMovementsClass`.

.. _manual-sec-ComponentClasses:

Component Classes
-----------------

The full list of component classes, their implementations, and the properties
and parameters each carries is generated directly from the source and is given
in the :doc:`/nodeComponents`. The classes currently defined are ``basic``,
``blackHole``, ``darkMatterProfile``, ``disk``, ``hotHalo``, ``NSC`` (nuclear
star cluster), ``position``, ``satellite``, ``spheroid``, and ``spin``.

.. note::

   The ``spin`` component class tracks halo *angular momentum*, not the
   dimensionless spin parameter :math:`\lambda`. Models which assign spins from
   a distribution, propagate them along branches, or build them up from merging
   progenitors do so through node operators---see, for example,
   :galacticus-class:`nodeOperatorHaloAngularMomentumRandom` and
   :galacticus-class:`nodeOperatorHaloAngularMomentumVitvitska2002`.


Properties and Deferred Rates
-----------------------------

Each property declared by a component carries attributes controlling whether it
can be read, written, and evolved as part of the :term:`ODE` system, together
with its output units and description. These are declared by the
``<component>`` directives in ``source/objects/nodes/components/``, from which
both the ``nodeComponent`` class hierarchy (see the developer guide) and the
:doc:`/nodeComponents` are generated at build time.

A few properties are *virtual* and have *deferred* rate functions. These act as
connection points between components: one component (or a node operator) pushes
a rate through the property, and the component which owns it decides what to do
with the material. The most important are:

``outflowingMass``, ``outflowingAngularMomentum``, ``outflowingAbundances``
   Declared by the ``hotHalo`` component. Galactic components which expel gas
   send that mass---plus metals and angular momentum---through these properties,
   where it is received into the circumgalactic medium. The ``verySimpleDelayed``
   implementation redirects them into its outflowed reservoir instead.

``massSink``
   Declared by the ``hotHalo`` component. Removes gas, with proportionate
   angular momentum and elements, from the hot gas halo.

``massGasSink``
   Declared by the ``spheroid`` and ``NSC`` components. Removes gas, with
   proportionate angular momentum and elements, from that component.

``energyGasInput``
   Declared by the ``spheroid`` component. Energy sent through this property is
   added to the gas of the spheroid and drives an outflow at a rate

   .. math::

      \dot{M}_\mathrm{outflow, spheroid} = \beta_\mathrm{spheroid, energy} {\dot{E}_\mathrm{gas, spheroid} \over V_\mathrm{spheroid}^2},

   where :math:`\beta_\mathrm{spheroid, energy}=`\ ``[componentSpheroid/efficiencyEnergeticOutflow]``.
   Input energy should be in units of :math:`M_\odot` km\ :math:`^2`
   s\ :math:`^{-2}` Gyr\ :math:`^{-1}` and must be positive---energy cannot be
   removed from the gas by this route.

.. _manual-sec-ComponentsIntegralEvolution:

Inactive Properties and Integral Evolution
------------------------------------------

Stellar luminosities are expensive to carry through the :term:`ODE` solver, and
their rates of change do not feed back on any other property. They can therefore
be treated as *inactive*: solved for by direct integration across each timestep,
rather than as part of the differential evolution system. This is enabled by
``[componentDisk/inactiveLuminositiesStellar]`` and
``[componentSpheroid/inactiveLuminositiesStellar]``. The ``standard``
``satellite`` implementation offers the same treatment for its bound mass via
``[componentSatellite/inactiveBoundMass]``.

For the spheroid the integral is straightforward. Across a timestep :math:`t_i`
to :math:`t_{i+1}`,

.. math::

   \Delta L_\mathrm{spheroid, stars} = \int_{t_{i}}^{t_{i+1}} \mathrm{d}t \, \Psi \mathcal{L}_\lambda(t_0-t,Z_\mathrm{gas}),

which is just the usual production of starlight by star formation, with
:math:`\mathcal{L}_\lambda(t,Z)` the luminosity-to-mass ratio in the relevant
band for a stellar population of age :math:`t` and metallicity :math:`Z`. Since
:math:`M_\mathrm{stars, spheroid, formed}(t)` has already been solved for during
differential evolution, the substitution :math:`\Psi = \dot{M}_\mathrm{stars,
spheroid, formed}(t)` avoids recomputing the star formation rate.

The disk is more involved, because material---and the starlight associated with
it---is transferred out of the disk by bar instabilities. The disk therefore
tracks a quantity :math:`f_\mathrm{retained}`, the fraction of mass retained in
the disk under transfer processes whose rate is proportional to the current mass,
evolved as

.. math::

   \dot{f}_\mathrm{retained} = - {f_\mathrm{retained} \over \tau_\mathrm{bar}},

where :math:`\tau_\mathrm{bar}` is the bar instability timescale (see
:galacticus-class:`galacticDynamicsBarInstabilityClass`). Its initial value is
arbitrary, as only ratios of the quantity are used. The luminosity changes are
then

.. math::

   \Delta L_\mathrm{disk, stars} = \int_{t_{i}}^{t_{i+1}} \mathrm{d}t \left[ \Psi \mathcal{L}_\lambda(t_0-t,Z_\mathrm{gas}) {f_\mathrm{retained}(t_{i+1}) \over f_\mathrm{retained}(t)} + L_\mathrm{disk, stars}(t_i) {\dot{f}_\mathrm{retained}(t) \over f_\mathrm{retained}(t_{i})} \right],

and

.. math::

   \Delta L_\mathrm{spheroid, stars} = \int_{t_{i}}^{t_{i+1}} \mathrm{d}t \left[ \Psi \mathcal{L}_\lambda(t_0-t,Z_\mathrm{gas}) \left( 1 - {f_\mathrm{retained}(t_{i+1}) \over f_\mathrm{retained}(t)} \right) - L_\mathrm{disk, stars}(t_i) {\dot{f}_\mathrm{retained}(t) \over f_\mathrm{retained}(t_{i})} \right].

The first term in the disk integral is the usual production of starlight,
reduced by a factor :math:`f_\mathrm{retained}(t_{i+1}) / f_\mathrm{retained}(t)`
to account for that fraction which will still be in the disk at the end of the
timestep. The second term accounts for transfer of starlight produced at
:math:`t<t_i`, since it integrates to :math:`L_\mathrm{disk, stars}(t_i)
(f_\mathrm{retained}(t_{i+1}) / f_\mathrm{retained}(t_{i}) - 1)`, which is the
change in the starlight of the disk over the timestep. The spheroid integral
simply reflects that starlight not retained in the disk is transferred to the
spheroid.

.. _manual-sec-ComponentsStructure:

Galactic Structure
------------------

The ``standard`` disk and spheroid implementations participate in the galactic
structure solver, which determines their radii. Both use a *dimensionless* mass
distribution object, scaled by the component's current mass and length scale, so
any implementation of :galacticus-class:`massDistributionClass` with the
appropriate symmetry may be used. The distribution is chosen with
``[componentDisk/massDistributionDisk]`` (which must be cylindrically symmetric;
default ``exponentialDisk``) and
``[componentSpheroid/massDistributionSpheroid]`` (which must be spherically
symmetric; default ``hernquist``).

Disk radii
~~~~~~~~~~

The radial size of the disk is found by solving for equilibrium---that is, the
radius at which the angular momentum of material is sufficient to provide
rotational support---at the radius ``[componentDisk/radiusStructureSolver]``,
given in units of the disk scale length. Converting from the mean specific
angular momentum of the disk to the specific angular momentum at that radius
uses the ratio ``[componentDisk/ratioAngularMomentumSolverRadius]``. Assuming a
flat rotation curve this ratio is :math:`I_1/I_2`, where

.. math::

   I_n = \int_0^\infty \Sigma(R) R^n \mathrm{d}R

and :math:`\Sigma(R)` is the disk surface density profile; this is the default,
computed from the selected mass distribution. Where either moment is infinite a
default of :math:`1/2` is used instead.

Setting ``[componentDisk/structureSolverUseCole2000Method]``\ :math:`=`\
``true`` alters this behavior to match the structure solver of
:cite:t:`cole_hierarchical_2000`, in which adiabatic contraction of the dark
matter halo is solved for assuming that the disk has a *spherical* mass
distribution. In this case the specific angular momentum passed to the structure
solver is reduced by a term accounting for the difference between the rotation
curve of a thin disk and that of a spherical mass distribution of the same mass,

.. math::

   j \rightarrow \left[ j^2 - \left( V_\mathrm{disk}^2(R_\mathrm{s}) R_\mathrm{s} - M_\mathrm{disk}(<R_\mathrm{s}) \right) \mathrm{G} M_\mathrm{disk} r \right]^{1/2},

where :math:`V_\mathrm{disk}` and :math:`M_\mathrm{disk}(<R)` are evaluated for
the dimensionless mass distribution at the structure solver radius
:math:`R_\mathrm{s}`. The result is truncated at zero to trap the imaginary
values which can arise while the solver explores the allowed range of radii.
Note that in this case---as in :cite:author:`cole_hierarchical_2000`
:cite:year:`cole_hierarchical_2000`---the resulting disk will not precisely
satisfy :math:`j(r) = r V_\mathrm{c}(r)`, where :math:`V_\mathrm{c}(r)` is the
net rotation curve.

Spheroid angular momentum
~~~~~~~~~~~~~~~~~~~~~~~~~

The spheroid tracks a *pseudo*-angular momentum---effectively the angular
momentum the spheroid would have were it rotationally supported rather than
pressure supported. The parameter
``[componentSpheroid/ratioAngularMomentumScaleRadius]`` controls the ratio of
the specific pseudo-angular momentum at the scale radius to the mean specific
pseudo-angular momentum. Its default is :math:`I_2/I_3`, where

.. math::

   I_n = \int_0^\infty \rho(r) r^n \mathrm{d}r

and :math:`\rho(r)` is the spheroid density profile, which is appropriate for a
flat rotation curve. For some profiles---the Hernquist profile, for
example---one or both of :math:`I_2` and :math:`I_3` is infinite; in such cases
a default of :math:`1/2` is used. If a finite truncation radius or a different
rotation curve is assumed this ratio may be finite, and the parameter allows
these assumptions to be controlled directly.
