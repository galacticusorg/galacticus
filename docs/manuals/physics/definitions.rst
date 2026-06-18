Definitions and Conventions Used in Galacticus
==============================================

Galacticus adopts various definitions and conventions internally. These are explained below.

Halo Masses and Dark Matter Mass
--------------------------------

Halo masses require some care in specifying exactly what mass they represent due to the way in which merger trees are typically created. For example, when merger trees are extracted from N-body simulations, those simulations frequently represent *all* matter as collisionless. That is, the simulation contains a density :math:`\Omega_\mathrm{M}=\Omega_\mathrm{DM}+\Omega_\mathrm{b}` which is the sum of dark and baryonic matter densities, but all of this mass is represented as collisionless particles. Similarly, masses in merger trees built through Monte Carlo techniques typically represent all mass as collisionless.

The exact way in which masses within Galacticus are defined and used in specified in the following subsections.

Masses in the Basic Component
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``basic`` component (see Section :galacticus-ref:`ComponentBasicProperties`) tracks the mass of each halo as defined in the merger tree. As such, it should be considered to be the mass which the halo would have if baryonic matter behaved just as dark matter. Note that these masses are inclusive of subhalos---that is, the mass of a host halo includes the mass of all of its subhalos.

Dark Matter Profiles
~~~~~~~~~~~~~~~~~~~~

The dark matter profile functions (see :galacticus-class:`darkMatterProfileDMO`) return masses and densities etc. which are normalized to match the mass of the ``basic`` component at the virial radius of the halo. As such, their returned values should be considered to represent the case where baryonic matter behaves as dark matter. This is a convention, and is useful for calculations of large scale structure for example.

Galactic Structure Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The various galactic structure functions assume that the masses/densities/etc. reported by the dark matter profile functions should be scaled by a factor :math:`(\Omega_\mathrm{M}-\Omega_\mathrm{b})/\Omega_\mathrm{b}` to leave only the dark matter part of the profile. Baryonic contributions to the mass/density/etc. will be provided by the components representing those mass distributions.

Satellite Virial Orbits
~~~~~~~~~~~~~~~~~~~~~~~

These functions (see :galacticus-class:`virialOrbit`) typically use the ``basic`` component mass in determining parameters of an orbit, since they are typically calibrated to simulations of collisionless matter only.

Satellite Merging Timescales
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These functions (see :galacticus-class:`satelliteMergingTimescales`) typically use the ``basic`` component mass in determining parameters of an orbit, since they are typically calibrated to simulations of collisionless matter only.

Dynamical Friction
~~~~~~~~~~~~~~~~~~

These functions (see :galacticus-class:`satelliteDynamicalFriction`) evaluate densities through the relevant galactic structure function, and so correctly account for the fraction of the ``basic`` component mass which is in the form of dark matter.

Galactic Structure Radius Solvers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These functions :galacticus-class:`galacticStructureSolver`) determine the radii of galactic components (such as disk and spheroid), typically by iteratively seeking a solution in which their angular momenta and radii are consistent (assuming rotational support) with the net gravitational potential of the entire system (galaxy plus dark matter halo).

Luminosity Units
----------------

Galaxy luminosities are output in the :term:`AB magnitude` system, such that a luminosity of :math:`1` corresponds to an object of :math:`0^\mathrm{th}` absolute magnitude in the :term:`AB magnitude` system. This implies that the luminosities are in units of :math:`4.4659\times 10^{13}` W/Hz.

.. _manual-sec-GalacticusVelocityDefinitions:

Peculiar Velocities
-------------------

Velocities in Galacticus are always *physical* velocities. When reading merger tree properties (including velocities) from file it is often convenient to store velocities without the Hubble flow contribution, as "peculiar velocities", in the file---see `here <https://github.com/galacticusorg/galacticus/wiki/Merger-Tree-File-Format#forest-halos-group>`_ for how to specify whether or not  the velocities included in the file include the Hubble flow or not.

If peculiar velocities are stored it is important to use the same definition of peculiar velocity as is used by Galacticus. Defining :math:`t` to be physical time and :math:`\mathbf{x}` to be comoving position, Galacticus uses the conventional definition of peculiar velocity in a cosmological context, namely that it is the deviation of the physical velocity from the Hubble flow. Physical coordinates are given by :math:`\mathbf{r} = a\mathbf{x}`, so the peculiar velocity is

.. math::

   \mathbf{v}_\mathrm{pec} \equiv {\mathrm{d} \mathbf{r} \over \mathrm{d} t} - H \mathbf{r} = a {\mathrm{d} \mathbf{x}\over\mathrm{d} t} = {\mathrm{d}\mathbf{x}\over\mathrm{d}\eta},

where :math:`\mathrm{d}\eta = \mathrm{d}t/a` is conformal time.

Gravitational Potentials
------------------------

Gravitational potentials are measured in velocity units (i.e. km\ :math:`^2`/s\ :math:`^2`), and the arbitrary constant offset is chosen such that the total gravitational potential in any halo at the virial radius is :math:`\Phi(r_\mathrm{virial})=-V_\mathrm{virial}^2`. This choice is made for two reasons:

#. some mass distributions used have potentials which diverge as :math:`r\rightarrow\infty`, so the usual choice of :math:`\Phi(r) \rightarrow 0` as :math:`r \rightarrow \infty` is not applicable;
#. this choice is consistent with the potential at the virial radius of the halo considered as a point mass as is used in Keplerian orbit calculations.

Note that the choice of constant offset for the potential of any mass distribution or galactic component is irrelevant---the galactic structure function which computes potential will ensure that the potential is always offset to match the definition given above.

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
   specifies the units of the ``radius`` element---valid options are ``diskRadius``, ``hotHaloOuterRadius``, ``diskHalfMassRadius``, ``spheroidRadius``, ``spheroidHalfMassRadius``, ``darkMatterScaleRadius``, ``virialRadius``, just ``radius`` (which implies radii are given in units of Mpc), ``galacticMassFraction{<fraction>}``, ``stellarMassFraction{<fraction>}``, or  ``galacticLightFraction{<fraction>}{<luminosity>}``, where the final three form specify a radius containing a fixed ``<fraction>`` of the galactic or stellar mass, or light respectively (for the case of galactic light, ``<luminosity>`` specifies the band, e.g. ``SDSS_r:rest:z0.0000``);

``componentType``
   specifies which components of the node should be counted---allowed values are ``all``, ``disk``, ``spheroid``, ``hotHalo``, ``darkHalo``, and ``blackHole``;

``massType``
   specifies which types of mass should be counted---allowed values are ``all``, ``dark``, ``baryonic``, ``galactic``, ``gaseous``, ``stellar``, and ``blackHole``.

Half-mode and quarter-mode masses
---------------------------------

For non-cold dark matter models with a matter power spectrum suppressed on small scales, a characteristic suppression scale can be defined by comparing the linear transfer function with that in the\ :term:`CDM` model. A popular definition is the half-mode (quarter-mode) mass, which corresponds to the wavenumber at which the transfer function is suppressed by a factor of two (four) relative to a :term:`CDM` transfer function.
