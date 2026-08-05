Definitions and Conventions Used in Galacticus
==============================================

Galacticus adopts various definitions and conventions internally. These are explained below.

Halo Masses and Dark Matter Mass
--------------------------------

Halo masses require some care in specifying exactly what mass they represent due to the way in which merger trees are typically created. For example, when merger trees are extracted from N-body simulations, those simulations frequently represent *all* matter as collisionless. That is, the simulation contains a density :math:`\Omega_\mathrm{M}=\Omega_\mathrm{DM}+\Omega_\mathrm{b}` which is the sum of dark and baryonic matter densities, but all of this mass is represented as collisionless particles. Similarly, masses in merger trees built through Monte Carlo techniques typically represent all mass as collisionless.

The exact way in which masses within Galacticus are defined and used in specified in the following subsections.

Masses in the Basic Component
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``basic`` component (see :ref:`manual-sec-ComponentClasses`) tracks the mass of each halo as defined in the merger tree. As such, it should be considered to be the mass which the halo would have if baryonic matter behaved just as dark matter. Note that these masses are inclusive of subhalos---that is, the mass of a host halo includes the mass of all of its subhalos.

Dark Matter Profiles
~~~~~~~~~~~~~~~~~~~~

The dark matter profile functions (see :galacticus-class:`darkMatterProfileDMO`) return masses and densities etc. which are normalized to match the mass of the ``basic`` component at the virial radius of the halo. As such, their returned values should be considered to represent the case where baryonic matter behaves as dark matter. This is a convention, and is useful for calculations of large scale structure for example.

.. _manual-sec-massDistributions:

Mass Distributions
~~~~~~~~~~~~~~~~~~

All structural properties of a node---densities, enclosed masses, gravitational potentials, rotation curves, surface densities, and so on---are obtained from :galacticus-class:`massDistributionClass` objects. A ``massDistribution`` is a fully general description of the spatial distribution of some mass, providing methods such as ``density()``, ``densitySphericalAverage()``, ``massEnclosedBySphere()``, ``massEnclosedByCylinder()``, ``potential()``, ``rotationCurve()``, ``surfaceDensity()``, ``acceleration()``, and ``tidalTensor()``, together with inverse functions such as ``radiusEnclosingMass()`` and ``radiusEnclosingDensity()``. Kinematic properties (e.g. ``velocityDispersion1D()``, obtained by solving the Jeans equation) are provided by an associated ``kinematicsDistribution`` object, retrieved via the ``kinematicsDistribution()`` method.

Each node component which carries mass supplies its own mass distribution, tagged with a *component type* and a *mass type*. For example, the ``standard`` disk component supplies two distributions---one of ``componentTypeDisk``/``massTypeStellar`` and one of ``componentTypeDisk``/``massTypeGaseous``---while the ``scale`` dark matter profile component supplies the dark matter halo distribution.

The mass distribution of an entire node is obtained through the ``massDistribution`` method of ``treeNode``:

.. code-block:: fortran

   class(massDistributionClass), pointer :: massDistribution_

   massDistribution_ => node%massDistribution(componentType,massType)

This returns a :galacticus-class:`massDistributionComposite`---the superposition of the individual component distributions---restricted to those which match the requested ``componentType`` and ``massType``. Both arguments are optional and default to ``componentTypeAll`` and ``massTypeAll`` respectively, so that ``node%massDistribution()`` returns the total mass distribution of the node. (The object returned is reference counted, and so must be released, via ``<objectDestructor>``, once it is no longer needed.) When a component contributes no mass matching the request the result is a :galacticus-class:`massDistributionZero`, so callers need not special-case empty components.

Component and Mass Types
^^^^^^^^^^^^^^^^^^^^^^^^

The allowed ``componentType`` values are ``all``, ``disk``, ``spheroid``, ``nuclearStarCluster``, ``hotHalo``, ``coldHalo``, ``darkHalo``, ``blackHole``, and ``darkMatterOnly``; the allowed ``massType`` values are ``all``, ``dark``, ``baryonic``, ``galactic``, ``gaseous``, ``stellar``, and ``blackHole``. (These are the same enumerations used by the :ref:`radius specifiers <manual-sec-radiusSpecifiers>`.)

Selection is hierarchical---``massTypeBaryonic`` matches both gaseous and stellar mass, while ``massTypeGalactic`` matches gaseous, stellar, and black hole mass, but only in the disk, spheroid, and black hole components (i.e. it excludes the circumgalactic medium). Thus, for example,

.. code-block:: fortran

   massDistribution_ => node%massDistribution(massType=massTypeStellar)

returns the combined stellar mass distribution of the disk, spheroid, and nuclear star cluster, whereas

.. code-block:: fortran

   massDistribution_ => node%massDistribution(componentTypeDisk,massTypeGaseous)

returns just the gaseous distribution of the disk.

Dark Matter Mass Convention
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The distribution returned for the dark matter halo depends on which component type was requested, and this is where the distinction between "all matter" and "dark matter only" (see above) is made:

* ``componentTypeDarkHalo``---and therefore also ``componentTypeAll``---returns the distribution built by the :galacticus-class:`darkMatterProfile` class. This represents the *dark matter* content of the halo alone. The :galacticus-class:`darkMatterProfileDarkMatterOnly` implementation, for instance, takes the profile supplied by :galacticus-class:`darkMatterProfileDMO` (which, as described above, is normalized to the ``basic`` component mass and so treats baryons as dark matter) and rescales its mass by a factor

  .. math::

     f_\mathrm{DM} = {\Omega_\mathrm{M}-\Omega_\mathrm{b} \over \Omega_\mathrm{M}},

  leaving only the dark matter part of the profile. Baryonic contributions to the mass/density/etc. are then supplied by the mass distributions of the components which represent those baryons (disk, spheroid, hot halo, etc.), so that the composite ``componentTypeAll`` distribution accounts for all of the mass in the node without double counting it. Other implementations of :galacticus-class:`darkMatterProfile`---such as the default, :galacticus-class:`darkMatterProfileAdiabaticGnedin2004`---apply the same rescaling, but additionally modify the *shape* of the profile in response to the baryons (here, via adiabatic contraction).

* ``componentTypeDarkMatterOnly`` returns the unmodified :galacticus-class:`darkMatterProfileDMO` distribution---that is, the profile normalized to the full ``basic`` component mass, with no :math:`f_\mathrm{DM}` rescaling and no baryonic back-reaction. This is the appropriate choice when comparing with dark matter-only N-body simulations, or when a quantity is calibrated against such simulations. Note that this component type is *not* included in ``componentTypeAll``: no other component responds to it, so a ``componentTypeDarkMatterOnly`` request yields a distribution containing only the dark matter-only halo profile.

Satellite Virial Orbits
~~~~~~~~~~~~~~~~~~~~~~~

These functions (see :galacticus-class:`virialOrbit`) typically use the ``basic`` component mass in determining parameters of an orbit, since they are typically calibrated to simulations of collisionless matter only.

Satellite Merging Timescales
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These functions (see :galacticus-class:`satelliteMergingTimescales`) typically use the ``basic`` component mass in determining parameters of an orbit, since they are typically calibrated to simulations of collisionless matter only.

Dynamical Friction
~~~~~~~~~~~~~~~~~~

These functions (see :galacticus-class:`satelliteDynamicalFriction`) evaluate densities through the relevant mass distribution, and so correctly account for the fraction of the ``basic`` component mass which is in the form of dark matter.

Galactic Structure Radius Solvers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These functions (see :galacticus-class:`galacticStructureSolver`) determine the radii of galactic components (such as disk and spheroid), typically by iteratively seeking a solution in which their angular momenta and radii are consistent (assuming rotational support) with the net gravitational potential of the entire system (galaxy plus dark matter halo).

Luminosity Units
----------------

Galaxy luminosities are output in the :term:`AB magnitude` system, such that a luminosity of :math:`1` corresponds to an object of :math:`0^\mathrm{th}` absolute magnitude in the :term:`AB magnitude` system. This implies that the luminosities are in units of :math:`4.4659\times 10^{13}` W/Hz.

.. _manual-sec-GalacticusVelocityDefinitions:

Peculiar Velocities
-------------------

Velocities in Galacticus are always *physical* velocities. When reading merger tree properties (including velocities) from file it is often convenient to store velocities without the Hubble flow contribution, as "peculiar velocities", in the file---see the :ref:`forest halos group <forest-halos-group>` for how to specify whether or not  the velocities included in the file include the Hubble flow or not.

If peculiar velocities are stored it is important to use the same definition of peculiar velocity as is used by Galacticus. Defining :math:`t` to be physical time and :math:`\mathbf{x}` to be comoving position, Galacticus uses the conventional definition of peculiar velocity in a cosmological context, namely that it is the deviation of the physical velocity from the Hubble flow. Physical coordinates are given by :math:`\mathbf{r} = a\mathbf{x}`, so the peculiar velocity is

.. math::

   \mathbf{v}_\mathrm{pec} \equiv {\mathrm{d} \mathbf{r} \over \mathrm{d} t} - H \mathbf{r} = a {\mathrm{d} \mathbf{x}\over\mathrm{d} t} = {\mathrm{d}\mathbf{x}\over\mathrm{d}\eta},

where :math:`\mathrm{d}\eta = \mathrm{d}t/a` is conformal time.

Gravitational Potentials
------------------------

Gravitational potentials are measured in velocity units (i.e. km\ :math:`^2`/s\ :math:`^2`).

No single, global convention is imposed on the arbitrary constant offset in the potential. Instead, each :galacticus-class:`massDistributionClass` fixes its own zero point. Most adopt the usual choice of :math:`\Phi(r) \rightarrow 0` as :math:`r \rightarrow \infty` (e.g. :galacticus-class:`massDistributionNFW`), but this is not always possible---the potential of :galacticus-class:`massDistributionIsothermal`, for example, diverges as :math:`r \rightarrow \infty` and so is instead measured relative to that distribution's own reference length. The ``potential()`` method of a :galacticus-class:`massDistributionComposite` simply sums the potentials of its constituent distributions, each retaining its own zero point.

Absolute potentials are therefore not meaningful to compare between different mass distributions. Physically meaningful quantities should instead be computed from *differences* in the potential, for which the ``potentialDifference()`` method is provided:

.. code-block:: fortran

   potential = massDistribution_%potentialDifference(coordinates1,coordinates2)

Where a particular zero point is required, it is established by the calculation which needs it, rather than by the mass distribution. In particular, ``keplerOrbit`` objects (see :galacticus-class:`virialOrbit`) assume that the potential at the virial radius of the halo is :math:`\Phi(r_\mathrm{virial})=-V_\mathrm{virial}^2`, which is consistent with the potential at the virial radius of the halo considered as a point mass, as is used in Keplerian orbit calculations. That offset is applied by the ``Satellite_Orbits`` module, which evaluates the potential relative to the virial radius using ``potentialDifference()`` and then subtracts :math:`\mathrm{G}M_\mathrm{virial}/r_\mathrm{virial}`.

.. _manual-sec-radiusSpecifiers:

Radius Specifiers
-----------------

Several :galacticus-class:`nodePropertyExtractorClass`\ s extract properties at one or more radii in a node. Galacticus provides a flexible way of specifying such radii as fractions of physical quantities such as the virial radius, disk half-mass radius, etc.

Each radius specifier should take the form:

.. code-block:: none

     radiusType:componentType:massType:radius

The elements of this colon-separated specifier determine the radius at which a property is computed, which components/mass types should be counted, and whether baryonic loading of the halo should be accounted for. The elements have the following meaning:

``radius``
   the numerical value of the radius at which to compute the property (with units specified by the ``radiusType`` element);

``radiusType``
   specifies the units of the ``radius`` element---valid options are ``diskRadius``, ``hotHaloOuterRadius``, ``diskHalfMassRadius``, ``spheroidRadius``, ``spheroidHalfMassRadius``, ``darkMatterScaleRadius``, ``virialRadius``, ``solitonRadiusCore``, ``solitonRadiusSoliton``, just ``radius`` (which implies radii are given in units of Mpc), ``galacticMassFraction{<fraction>}``, ``stellarMassFraction{<fraction>}``, or  ``galacticLightFraction{<fraction>}{<luminosity>}``, where the final three form specify a radius containing a fixed ``<fraction>`` of the galactic or stellar mass, or light respectively (for the case of galactic light, ``<luminosity>`` specifies the band, e.g. ``SDSS_r:rest:z0.0000``);

   The ``solitonRadiusCore`` and ``solitonRadiusSoliton`` options express radii in units of the core and transition radii of the :term:`FDM` soliton respectively. These are meaningful only in models which grow solitons---see :galacticus-class:`darkMatterProfileDMOSolitonNFW`. Where they are not (because the model contains no soliton-forming dark matter profile, or because no soliton formed in the halo in question), the radius is *undefined*: no property is computed there, and a large negative sentinel value is written to the output instead. Filter such entries out before analysis;

``componentType``
   specifies which components of the node should be counted---allowed values are ``all``, ``disk``, ``spheroid``, ``hotHalo``, ``darkHalo``, and ``blackHole``;

``massType``
   specifies which types of mass should be counted---allowed values are ``all``, ``dark``, ``baryonic``, ``galactic``, ``gaseous``, ``stellar``, and ``blackHole``.

Zero radii
^^^^^^^^^^

A radius specifier can evaluate to zero---for example ``diskRadius:all:all:1.0`` in a node whose disk has zero radius. This is common: it happens whenever the component on which the radius is based is absent or empty in the node in question, and is quite distinct from an *undefined* radius (see ``solitonRadiusCore`` above), which is reported with a sentinel value instead.

Such a radius is perfectly well defined, but not every property can be evaluated there:

* enclosed masses (:galacticus-class:`nodePropertyExtractorMassProfile`), projected masses (:galacticus-class:`nodePropertyExtractorProjectedMass`), rotation curves, and velocity dispersions all exist at zero radius, and are reported there---note that the enclosed and projected masses are the mass of any central point mass, such as a black hole, and not necessarily zero;
* densities (:galacticus-class:`nodePropertyExtractorDensityProfile`, :galacticus-class:`nodePropertyExtractorDensityDMOProfile`) exist only if the mass distribution being evaluated has a finite central density---which is true of a disk, but not of an :term:`NFW` halo;
* projected densities (:galacticus-class:`nodePropertyExtractorProjectedDensity`) are never evaluated at zero radius, as the line of sight integral is performed in the logarithm of radius.

Where the property does not exist, the model stops with an error naming the offending specifier and the node in which it was evaluated. To carry on regardless, set ``zeroRadiusIsFatal`` to ``false`` in the property extractor concerned: the undefined-value sentinel is then written in place of the property, and should be filtered out before analysis. The radius itself is still reported (as zero), since---unlike an undefined radius---it is perfectly well defined.

Half-mode and quarter-mode masses
---------------------------------

For non-cold dark matter models with a matter power spectrum suppressed on small scales, a characteristic suppression scale can be defined by comparing the linear transfer function with that in the\ :term:`CDM` model. A popular definition is the half-mode (quarter-mode) mass, which corresponds to the wavenumber at which the transfer function is suppressed by a factor of two (four) relative to a :term:`CDM` transfer function.
