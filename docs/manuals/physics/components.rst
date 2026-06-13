.. _manual-sec-Components:

Node Components
===============

In addition to the implementations described here, each component class has a "``null``" implementation. Selecting this implementation---which has no properties and does not respond to any events---effectively switches off the relevant component class. Of course, this is safe only if none of the other active implementations expect to get or set properties of the component class (or if they rely on a sensible implementation of that class).

(Supermassive) Black Hole
-------------------------

"Standard" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The standard black hole implementation defines the following properties:

``mass``
   The mass of the black hole: :math:`M_\bullet` ``[blackHoleMass]``.

``spin``
   The spin of the black hole, :math:`j_\bullet` ``[blackHoleSpin]``.

``radialPosition``
   The radial position of the black hole: :math:`r_\bullet` ``[blackHoleRadialPosition]``

Initialization
^^^^^^^^^^^^^^

Black holes are not initialized, they are created (with a seed mass given by ``blackHoleSeedMass`` and zero spin) as needed.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

In the standard black implementation the mass and spin evolve as:

.. math::

   \dot{M}_\bullet &=& (1-\epsilon_\mathrm{radiation}-\epsilon_\mathrm{jet}) \dot{M}_0 \\
   \dot{j}_\bullet &=& \dot{j}(M_\bullet,j_\bullet,\dot{M}_0),

where :math:`\dot{M}_0` is the rest mass accretion rate, :math:`\epsilon_\mathrm{radiation}` is the radiative efficiency of the accretion flow feeding the black hole, :math:`\epsilon_\mathrm{jet}` is the efficiency with which accretion power is converted to jet power and :math:`\dot{j}(M_\bullet,j_\bullet,\dot{M}_0)` is the spin-up function of that accretion flow (see :galacticus-class:`accretionDisks`). The rest mass accretion rate is computed assuming Bondi-Hoyle-Lyttleton accretion from the spheroid gas reservoir (with an assumed temperature of ``[bondiHoyleAccretionTemperatureSpheroid]``) enhanced by a factor of ``[bondiHoyleAccretionEnhancementSpheroid]`` and from the host halo (with whatever temperature that hot halo temperature profile specifies; see :galacticus-class:`hotHaloTemperatureProfile`) enhanced by a factor of ``[bondiHoyleAccretionEnhancementHotHalo]``. For accretion from the hot halo, the Bondi radius is limited to the outer radius of the hot halo. Additionally, the accretion rate is limited to:

.. math::

   \dot{M}_\mathrm{0,~hot~halo,~maximum} = M_\mathrm{hot}/\tau_\mathrm{sound~crossing},

where :math:`\tau_\mathrm{sound~crossing}=r_\mathrm{hot~halo~outer}/c_\mathrm{s}` where :math:`r_\mathrm{hot~halo~outer}` is the outer radius of the hot halo and :math:`c_\mathrm{s}` is the speed of sound in the hot halo.

If ``[bondiHoyleAccretionHotModeOnly]``\ :math:`=`\ ``true`` then the accretion occurs only from that fraction of the hot halo gas which was accreted in the "hot mode", otherwise accretion occurs from the entirety of the hot halo reservoir. In the first case a simple estimate of the hot mode fraction is made:

.. math::

   f_\mathrm{hot} = \left\{ \begin{array}{ll} 1 & \hbox{ if } x < 0.9 \\ y(x)^2[2 y(x) - 3]+1  & \hbox{ if } 0.9 \le x \le 1.0 \\ 0 & \hbox{ if } x > 1.0, \end{array} \right.

where :math:`x = r_\mathrm{cool}/r_\mathrm{virial}` and :math:`y(x)=[x-0.9]/[1.0-0.9]`.

The rest mass accretion rate is removed (as a mass sink) from the spheroid and hot halo components appropriately. The black hole is assumed to cause feedback in two ways:

Radio-mode
   If ``[blackHoleHeatsHotHalo]``\ :math:`=`\ ``true`` then any jet power from the black hole-accretion disk system (see :galacticus-class:`accretionDisks`) is included in the hot halo heating rate providing that the halo is in the slow cooling regime\ [#]_ (i.e. if the cooling radius is smaller than the virial radius; see, for example, :cite:author:`benson_cold_2010` :cite:year:`benson_cold_2010`);

Quasar-mode
   A mechanical wind luminosity of :cite:p:`ostriker_momentum_2010`

   .. math::

      L_\mathrm{wind} = \epsilon_{\bullet, \mathrm{wind}} H(\epsilon_\mathrm{radiation},1,s) \dot{M}_0 \clight^2,

   where :math:`\epsilon_{\bullet \mathrm{wind}}=`\ ``[blackHoleWindEfficiency]`` is the black hole wind efficiency, :math:`s`\ =\ ``blackHoleWindEfficiencyScalesWithRadiativeEfficiency``, and

   .. math::

      H(a,b,c) = \left\{\begin{array}{ll}a & \hbox{ if } c=\hbox\mono{true} \\b & \hbox{ if } c=\hbox\mono{false},\end{array}\right.

   is added to the gas :term:`component` of the spheroid (which, presumably, will respond with an outflow for example---see Section :galacticus-ref:`ComponentSpheroid` for details of how specific implementations of the spheroid component respond to the addition of energy) if and only if the wind pressure (at the spheroid characteristic radius) is less than the typical thermal pressure in the spheroid gas :cite:p:`ciotti_feedbackcentral_2009`, i.e.

   .. math::

      P_\mathrm{wind} &<& P_\mathrm{ISM} \nonumber \\
      \frac{1}{2}\rho_\mathrm{wind} V_\mathrm{wind}^2 &<& {3 \mathrm{k_B} T_\mathrm{ISM} \langle \rho_\mathrm{ISM}\rangle \over 2 m_\mathrm{H}}.

   Since :math:`\Omega r^2 \rho_\mathrm{wind} V_\mathrm{wind}^3 = L_\mathrm{wind}` where :math:`\Omega` is the solid angle of the wind flow, this can be rearranged to give :math:`\langle\rho_\mathrm{ISM}\rangle > \rho_\mathrm{wind, critical}` where

   .. math::

      \rho_\mathrm{wind,critical} = {2 m_\mathrm{H} L_\mathrm{wind} \over 3 \Omega r^2 V_\mathrm{wind} \mathrm{k_B} T_\mathrm{ISM}}.

   This critical wind density is computed at the characteristic radius of the spheroid, :math:`r_\mathrm{spheroid}`, assuming :math:`V_\mathrm{wind}=10^4`\ km/s, :math:`T_\mathrm{ISM}=10^4`\ K and :math:`\Omega=\pi`, and the :term:`ISM` density is approximated by

   .. math::

      \langle\rho_\mathrm{ISM}\rangle = {3 M_\mathrm{gas, spheroid} \over 4 \pi} r_\mathrm{spheroid}^3.

   For numerical ease, the fraction, :math:`f_\mathrm{wind}`, of the wind luminosity added to the spheroid is adjusted smoothly through the :math:`\rho_\mathrm{ISM}\approx\rho_\mathrm{wind,critical}` region according to

   .. math::

      f_\mathrm{wind} = \left\{ \begin{array}{ll} 0 & \hbox{ if } x < 0, \\ 3x^2-2x^3 & \hbox{ if } 0 \le x \le 1, \\ 1 & \hbox{ if } x > 1, \end{array} \right.

   where :math:`x=\rho_\mathrm{ISM}/\rho_\mathrm{wind,critical}-1/2`.

The radial position, :math:`r_\bullet`, evolves according to the selected radial migration model (see :galacticus-class:`blackHoleBinarySeparationGrowthRate`).

Interactions between black hole triplets are accounted for if ``[tripleBlackHoleInteraction]``\ :math:`=`\ ``true`` (and if at least three black holes exist within the :term:`node` of course). In this case the triple is treated as consisting of an inner binary (assumed to be the central black hole and the black hole closest to it) and a third, singleton black hole. When the tertiary black hole reaches a separation of

.. math::

   a_\mathrm{h}= {\mathrm{G} (M_{\bullet, 1} + M_\mathrm{bullet, 2}) \over 4 \sigma^2}

it is assumed to undergo a triple interaction with the binary. Once a triple interaction occurs, no further triple interaction for the specific tertiary black hole can occur unless the host galaxy merges with another galaxy, at which point the black holes from the merging galaxy are eligible for another triple interaction in their new host.

The logic of what happens in a triple black hole interaction is taken from :cite:t:`volonteri_assembly_2003`. Labeling the central black hole as :math:`1`, its binary partner as :math:`2` and the tertiary black hole as :math:`3`, and defining

.. math::

   q_3 = {M_{\bullet, 3} \over M_{\bullet, 1} + M_\mathrm{\bullet, 2} },

then if :math:`q_3 \le 2` then, if :math:`M_{\bullet, 3} \le M_\mathrm{\bullet, 2}` we set

.. math::

   a_3 = {a_2 \over 1 + 0.4 q_3},

and define

.. math::

   E_\mathrm{bind} = {\mathrm{G} M_{\bullet, 3} M_{\bullet, 1} \over a_3},

and

.. math::

   \Delta K =0.4 q_3 E_\mathrm{bind},

:math:`i=3` and :math:`j=2`.

Otherwise if :math:`q_3 \le 2` and :math:`M_{\bullet, 3} > M_\mathrm{\bullet, 2}` we set

.. math::

   a_3 = { a_3 \over 1 +0.4 q_3},

and define

.. math::

   E_\mathrm{bind} = {\mathrm{G} M_{\bullet, 2} M_{\bullet, 1} \over a_2},

and

.. math::

   \Delta K =0.4 q_3 E_\mathrm{bind},

:math:`i=2` and :math:`j=3`.

Finally, if :math:`q_3 > 2`, then we set

.. math::

   a_3 =0.53 a_3,

and define

.. math::

   E_\mathrm{bind} = {\mathrm{G} M_{\bullet, 2} M_{\bullet, 1} \over a_2},

and

.. math::

   \Delta K =0.9 q_3 E_\mathrm{bind},

:math:`i=2`, :math:`j=3`.

Black hole :math:`i` is identified as the "ejected" hole, with black hole :math:`j` becoming the new binary member. Therefore

.. math::

   M_{\bullet, \mathrm ejected} = M_{\bullet, i}.

and

.. math::

   M_\mathrm{binary} = M_\mathrm{\bullet, j} + M_\mathrm{\bullet, 1}.

The imparted velocities of these two systems are

.. math::

   V_\mathrm{ejected} = \left[ {2 \Delta K \over (1+M_{\bullet, \mathrm ejected}/M_\mathrm{binary} ) M_{\bullet, \mathrm ejected} }\right]^{1/2}

and

.. math::

   V_\mathrm{binary} = \left[ {2 \Delta K \over (1+M_\mathrm{binary} /M_{\bullet, \mathrm ejected}) M_\mathrm{binary}}\right]^{1/2}.

If

.. math::

   \frac{1}{2} V_\mathrm{ejected|binary}^2 + \Phi(a_\mathrm{ejected|binary}) \ge 0

for either velocity, then that system is ejected from the node. Ejected black holes are removed from the node. If the binary is ejected the central black hole is replaced with a "null", zero mass placeholder.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None.\\

*Satellite merging:* The black holes in the two merging galaxies can be instantaneously merged, or taken at an initial separation (see :galacticus-class:`blackHoleBinaryInitialSeparation`), it is then evolved until reaching zero separation whereupon it is assumed to undergo merger. Properties are computed using the selected black hole binary merger model (see :galacticus-class:`blackHoleBinaryMerger`). In addition, the recoil velocity of the new black hole due to gravitational wave emission is computed using the selected method (see :galacticus-class:`blackHoleBinaryRecoil`), and if greater than the potential at the center of the galaxy, is assumed to have escaped the galaxy. Black holes which escape the galaxy are simply discarded and no longer tracked. For computational purposes, they are replaced with a "null", zero mass black hole at the center of the galaxy. If any other black hole comes within a distance

.. math::

   a_\mathrm{h} = {\mathrm{G} M_\bullet \over 4 \sigma^2},

where :math:`\sigma` is approximated to be the virial velocity of the dark matter halo, it is promoted to being the new "central" black hole of the node.\\

*Node promotion:* None.\\

Additional Output
^^^^^^^^^^^^^^^^^

If the ``[blackHoleOutputAccretion]`` input parameter is set to true, then rest mass accretion rate (in :math:`M_\odot` Gyr\ :math:`^{-1}`), jet power (in :math:`M_\odot` km\ :math:`^2` s\ :math:`^{-1}` Gyr\ :math:`^{-1}`) and radiative efficiency of the black hole\ [#]_ are output as ``blackHoleAccretionRate``, ``blackHoleJetPower`` and ``blackHoleRadiativeEfficiency`` respectively.

If the ``[blackHoleOutputData]`` input parameter is set to true, then the Masses (in :math:`M_\odot`), Spins (for now just a scalar with no direction), final Radius (in :math:`Mpc`), timescales (in :math:`Gyr`) until merger, accretion rates (in :math:`M_\odot per Gyr`) and radiative Efficiencies of all the black holes in the galaxy are given as outputs in the ``blackHole`` section of the output hdf5. This also saves the tree :term:`node` and merger tree index for further use when using the data.

The outputs of mergers are also automatically saved, as outputs in the ``blackHoleMergers`` section of the output hdf5. Those outputs are the time at which mergers happened and the mass ratio between  the two merging black holes.

"Simple" Implementation
~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The simple black hole implementation defines the following property:

``mass``
   The mass of the black hole: :math:`M_\bullet` ``[blackHoleMass]``.

Initialization
^^^^^^^^^^^^^^

Black holes are not initialized, they are created (with a seed mass given by ``blackHoleSeedMass``) as needed.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

In the simple black hole implementation the mass evolves as:

.. math::

   \dot{M}_\bullet &=& (1-\epsilon_\mathrm{wind}) \epsilon_\mathrm{BH} \dot{M}_{\star,\mathrm{spheroid}} \\

where :math:`\epsilon_\mathrm{BH}` is the ratio of rates at which the black hole and stellar spheroid grow. The black hole is assumed to cause feedback in two ways:

Radio-mode
   If ``[blackHoleHeatsHotHalo]``\ :math:`=`\ ``true`` and ``[blackHoleAccretesFromHotHalo]``\ :math:`=`\ ``false`` then a power :math:`\epsilon_\mathrm{heat} \epsilon_\mathrm{BH} \dot{M}_{\star,\mathrm{spheroid}} \mathrm{c}^2` where :math:`\epsilon_\mathrm{heat}=`\ ``[blackHoleHeatingEfficiency]`` is included in the hot halo heating rate providing that the halo is in the slow cooling regime (i.e. if the cooling radius is smaller than the virial radius; see, for example, :cite:author:`benson_cold_2010` :cite:year:`benson_cold_2010`) and the accretion rate onto the black hole is reduced by :math:`\epsilon_\mathrm{heat} \epsilon_\mathrm{BH} \dot{M}_{\star,\mathrm{spheroid}}`. If ``[blackHoleHeatsHotHalo]``\ :math:`=`\ ``true`` and ``[blackHoleAccretesFromHotHalo]``\ :math:`=`\ ``true`` then a power :math:`\epsilon_\mathrm{heat} \dot{M}_\mathrm{Eddington} \mathrm{c}^2` is included in the hot halo heating rate providing that the halo is in the slow cooling regime and the accretion rate onto the black hole is increased\ [#]_ by :math:`\dot{M}_\mathrm{Eddington} \epsilon_\mathrm{heat} (1-\epsilon_\mathrm{jet})/\epsilon_\mathrm{jet}`, where :math:`\epsilon_\mathrm{jet}=`\ ``[blackHoleJetEfficiency]``;

Quasar-mode
   A mechanical wind luminosity of :cite:p:`ostriker_momentum_2010`

   .. math::

      L_\mathrm{wind} = \epsilon_{\bullet, wind} \dot{M}_0 \clight^2,

   where :math:`\epsilon_{\bullet wind}=`\ ``[blackHoleWindEfficiency]`` is the black hole wind efficiency, is added to the gas :term:`component` of the spheroid (which, presumably, will respond with an outflow for example).

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None.\\

*Satellite merging:* The black holes in the two merging galaxies are instantaneously merged. Properties are computed using the selected black hole binary merger model (see :galacticus-class:`blackHoleBinaryMerger`.\\

*Node promotion:* None.\\

Additional Output
^^^^^^^^^^^^^^^^^

If the ``[blackHoleOutputAccretion]`` input parameter is set to true, then rest mass accretion rate (in :math:`M_\odot` Gyr\ :math:`^{-1}`) is output as ``blackHoleAccretionRate``.

Hot Halo
--------

"Very Simple" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The very simple hot halo implementation defines the following properties:

``mass``
   The mass of gas in the hot halo: :math:`M_\mathrm{hot}` ``[hotHaloMass]``.

and the following pipes:

``hotHaloCoolingMass``
   The net cooling rate of gas mass is sent through this pipe. Any :term:`component` may claim this pipe and connect to it, allowing it to receive the cooling gas.

``outflowingMass``
   Galactic components that wish to expel gas due to an outflow can send that mass  through this pipe, where it will be received into the hot halo component.

Initialization
^^^^^^^^^^^^^^

At initialization, nodes are assigned a mass of gas equal to their own mass, minus the mass of any progenitors, multiplied by :math:`\Omega_\mathrm{b}/\Omega_\mathrm{Matter}`.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

In the very simple hot halo implementation the hot gas mass and heavy element mass(es) evolves as:

.. math::

   \dot{M}_\mathrm{hot} = - \dot{M}_\mathrm{cooling} + \dot{M}_\mathrm{outflow},

where :math:`\dot{M}_\mathrm{cooling}` is the rate of mass loss from the hot halo due to cooling (see :galacticus-class:`coolingRate`. In the above :math:`\dot{M}_\mathrm{outflow}` is the net rate of outflow from any components in the node. For satellite galaxies, the outflow is instead directed to the hot halo of the host :term:`node`.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* Any hot gas from the merging halo is transferred to its host halo.\\

*Satellite merging:* Any hot halo of the satellite :term:`node` is added to that of the host :term:`node` and the hot halo :term:`component` removed from the satellite node.\\

*Node promotion:* Any hot halo of the parent :term:`node` is added to that of the :term:`node` prior to promotion.\\

*Halo formation:* None.\\

"Very Simple Delayed" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The delayed very simple hot halo implementation extends the "very simple" implementation by adding a reservoir for outflowed gas from which the hot halo is gradually replenished. It defines the following properties:

``outflowedMass``
   The mass of gas in the outflowed reservoir of the hot halo: :math:`M_\mathrm{outflowed}` ``[hotHaloOutflowedMass]``.

Initialization
^^^^^^^^^^^^^^

At initialization, nodes are assigned zero outflowed mass.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

This implementation steals the ``outflowingMass`` pipe from the very simple component and redirects it to the outflowed mass reservoir. The resulting rates of change of hot and outflowed masses are then:

.. math::

   \dot{M}_\mathrm{hot}      &=& - \dot{M}_\mathrm{cooling} + \dot{M}_\mathrm{reincorporation}, \nonumber \\
   \dot{M}_\mathrm{outflowed} &=& + \dot{M}_\mathrm{outflow} - \dot{M}_\mathrm{reincorporation}, \nonumber \\

where :math:`\dot{M}_\mathrm{reincorporation}` is the reincorporation rate of outflowed gas (see Section ).

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* Any outflowed gas from the merging halo is transferred to its host halo.\\

*Satellite merging:* Any outflowed halo of the satellite :term:`node` is added to that of the host :term:`node` and the outflowed mass :term:`component` removed from the satellite node.\\

*Node promotion:* Any outflowed gas of the parent :term:`node` is added to that of the :term:`node` prior to promotion.\\

*Halo formation:* None.\\

"Standard" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The standard hot halo implementation defines the following properties:

``unaccretedMass``
   The mass of gas which could have accreted onto the halo if it always accreted baryons and dark matter in the universal proportion, but which failed to do so (e.g. perhaps due to being photoheated to a high temperature and so being resistant to accretion into shallow potential wells): :math:`M_\mathrm{failed}`.

``mass``
   The mass of gas in the hot halo: :math:`M_\mathrm{hot}` ``[hotHaloMass]``.

``angularMomentum``
   The angular momentum of the gas in the hot halo, :math:`J_\mathrm{hot}` ``[hotHaloAngularMomentum]``.

``abundances``
   The mass(es) of heavy elements in gas in the hot halo, :math:`M_{Z, \mathrm{hot}}` ``[hotHalo{abundanceName}]``.

``outflowedMass``
   The mass of gas from outflows in the hot halo: :math:`M_\mathrm{outflowed}` ``[hotHaloOutflowedMass]``.

``outflowedAngularMomentum``
   The angular momentum of the outflowed gas in the hot halo, :math:`J_\mathrm{outflowed}` ``[hotHaloOutflowedAngularMomentum]``.

``outflowedAbundances``
   The mass(es) of heavy elements in outflowed gas, :math:`M_{Z, \mathrm{outflowed}}` ``[hotHaloOutflowed{abundanceName}]``.

``chemicals``
   The mass(es) of molecules in the hot gas, :math:`M_\mathrm{chemical}` ``[hotHaloChemicals{chemicalName}]``.

``outerRadius``
   The outer boundary radius of the hot halo: :math:`r_\mathrm{hot, outer}` ``[hotHaloOuterRadius]``.

``strippedMass``
   The mass of gas which has been stripped from the hot halo (by ram pressure or tidal forces for example): :math:`M_\mathrm{hot, stripped}`. This property is computed only if ``[hotHaloTrackStrippedGas]``\ :math:`=`\ ``true``.

``strippedAbundances``
   The mass(es) of heavy elements in gas that has been stripped from the hot halo (by ram pressure or tidal forces for example), :math:`M_{Z, \mathrm{hot, stripped}}`. These properties are computed only if ``[hotHaloTrackStrippedGas]``\ :math:`=`\ ``true``.

and the following pipes:

``heatSource``
   Energy sent through this pipe is added to the hot halo and used to offset the cooling rate (see below; heat pushed should be in units if :math:`M_\odot` (km/s)\ :math:`^2` Gyr\ :math:`^{-1}`).

``cooling[Mass|Angular_Momentum|Abundances]_To``
   The net cooling rate of gas mass (and metal content and magnitude of angular momentum) is sent through this pipe. Any :term:`component` may claim this pipe and connect to it, allowing it to receive the cooling gas.

``outflowing[Mass|AngularMomentum|Abundances]_To``
   Galactic components that wish to expel gas due to an outflow can send that mass (plus metals and angular momentum) through this pipe, where it will be received into the hot halo component.

``massSink``
   Removes gas (and proportionate amounts of angular momentum and elements) from the hot gas halo.

Initialization
^^^^^^^^^^^^^^

At initialization, any nodes with no children are assigned a hot halo mass, and failed accreted mass as dictated by the baryonic accretion method (see :galacticus-class:`accretionHalo`) and angular momentum based on the accreted mass and the halo spin parameter.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

In the standard hot halo implementation the hot gas mass and heavy element mass(es) evolves as:

.. math::

   \dot{M}_\mathrm{failed} &=& \dot{M}_\mathrm{failed~accretion} \\
   \dot{M}_\mathrm{hot} &=& \dot{M}_\mathrm{accretion} - \dot{M}_\mathrm{cooling} + \dot{M}_\mathrm{outflow,return} - \dot{M}_\mathrm{expelled} - \dot{M}_\mathrm{hot, stripped}, \\
   \dot{M}_{Z, \mathrm{hot}} &=& - \dot{M}_\mathrm{cooling} {M_{Z, \mathrm{hot}}\over M_\mathrm{hot}} + \dot{M}_{Z, \mathrm{outflow,return}} - \dot{M}_{Z,\mathrm{expelled}} - \dot{M}_{Z, \mathrm{hot, stripped}}, \\
   \dot{M}_\mathrm{chemical} &=& - [\dot{M}_\mathrm{cooling} + \dot{M}_\mathrm{expelled} + \dot{M}_\mathrm{hot, stripped}] {M_\mathrm{chemical}\over M_\mathrm{hot}} + f_\mathrm{chemical,outflow} \dot{M}_\mathrm{outflow,return} \nonumber \\
   & & + \dot{M}_\mathrm{chemical,reactions}, \\
   \dot{r}_\mathrm{hot, outer} &=& \left\{ \begin{array}{ll} {r_\mathrm{rp} - r_\mathrm{hot, outer} \over \tau_\mathrm{dynamical}} & \hbox{ if } r_\mathrm{rp} < r_\mathrm{hot, outer} \\
   0 & \hbox{ otherwise.} \end{array} \right. \\
   \dot{M}_\mathrm{hot, stripped} &=& -4 \pi \rho_\mathrm{hot}(r_\mathrm{hot, outer}) r_\mathrm{hot, outer}^2 \dot{r}_\mathrm{hot, outer} + \dot{M}_\mathrm{outflows} f_\mathrm{outflow, stripped} \\
   \dot{M}_{Z, \mathrm{hot, stripped}} &=& -4 \pi \rho_\mathrm{hot}(r_\mathrm{hot, outer}) r_\mathrm{hot, outer}^2 \dot{r}_\mathrm{hot, outer} (M_{Z, \mathrm{hot}} / M_\mathrm{hot}) \nonumber \\
   & & + \dot{M}_{Z, \mathrm{outflows}} f_\mathrm{outflow, stripped}

where :math:`r_\mathrm{rp}` is the ram pressure stripping radius as computed by the ``hotHaloRamPressureStripping`` method (see :galacticus-class:`hotHaloRamPressureStripping`), :math:`\dot{M}_\mathrm{accretion}` is the rate of growth of the hot :term:`component` due to accretion from the :term:`IGM` and :math:`\dot{M}_\mathrm{failed~accretion}` is the rate of failed accretion from the :term:`IGM` (these may include a :term:`component` due to transfer of mass from the failed to accreted reservoirs) and :math:`\dot{M}_\mathrm{cooling}` is the rate of mass loss from the hot halo due to cooling (see :galacticus-class:`coolingRate`---cooling rates are computed using the current :term:`node` if ``[hotHaloCoolingFromNode]``\ :math:`=`\ ``current node`` or from the formation :term:`node` if that parameter is set to ``formation node``) minus any heating rate defined as

.. math::

   \dot{M}_\mathrm{heating} = \dot{E}_\mathrm{input} / V_\mathrm{virial}^2,

where :math:`\dot{E}_\mathrm{input}` is the rate at which energy is being sent through the "energy input" pipe and :math:`V_\mathrm{virial}` is the virial velocity of the halo. The net cooling rate is never allowed to drop below zero. If the mass heating rate exceeds the mass cooling rate and ``[hotHaloExcessHeatDrivesOutflow]``\ :math:`=`\ ``false`` then the excess energy is not used and :math:`\dot{M}_\mathrm{expelled}=0`. Alternatively, if ``[hotHaloExcessHeatDrivesOutflow]``\ :math:`=`\ ``true`` then

.. math::

   \dot{M}_\mathrm{expelled} = \left\{ \begin{array}{ll} \alpha_\mathrm{expel} M_\mathrm{hot}/\tau_\mathrm{dynamical} & \hbox{ if } \dot{M}_\mathrm{heating} - \dot{M}_\mathrm{cool} > \alpha_\mathrm{expel} M_\mathrm{hot}/\tau_\mathrm{dynamical} \\ \dot{M}_\mathrm{heating} - \dot{M}_\mathrm{cool}, & \hbox{ otherwise,} \end{array} \right.

where :math:`\dot{M}_\mathrm{cool}` is the intrinsic cooling rate in the halo (i.e. the cooling rate in the absence of any heating) and :math:`\alpha_\mathrm{expel}=`\ ``[hotHaloExpulsionRateMaximum]`` limits the maximum rate at which mass can be expelled from the halo.

In the above, :math:`f_\mathrm{chemical,return}` if the mass fraction of each chemical species in the outflowed gas and is assumed to be equal to that given by the atomic ionization state functions (see :galacticus-class:`chemicalState`) at the virial temperature and mean density of the halo. Finally, :math:`\dot{M}_\mathrm{chemical,reactions}` represents the rate of change of masses of chemical species due to chemical and atomic processes and is computed using the chemical rates functions (see :galacticus-class:`chemicalReactionRate`). The angular momentum of the hot gas evolves as:

.. math::

   \dot{J}_\mathrm{hot} = \dot{M}_\mathrm{accretion} {\dot{J}_\mathrm{node} \over \dot{M}_\mathrm{node}} - \dot{M}_\mathrm{cooling} r_\mathrm{cool} V_\mathrm{rotate} + \dot{J}_\mathrm{outflow,return} - \dot{M}_\mathrm{expelled} {J_\mathrm{hot} \over M_\mathrm{hot}},

where :math:`\dot{M}_\mathrm{node}` and :math:`\dot{J}_\mathrm{node}` are defined in Section :galacticus-ref:`ComponentBasicProperties`. For the outflowed components:

.. math::

   \dot{M}_\mathrm{outflowed} &=& - \dot{M}_\mathrm{outflow,return} + \dot{M}_\mathrm{outflows} (1-f_\mathrm{outflow, stripped}), \\
   \dot{M}_{Z, \mathrm{outflowed}} &=& - \dot{M}_{Z, \mathrm{outflow,return}} + \dot{M}_{Z, \mathrm{outflows}} (1-f_\mathrm{outflow, stripped}), \\

and:

.. math::

   \dot{J}_\mathrm{outflowed} = - \dot{J}_\mathrm{outflow,return} + \dot{J}_\mathrm{outflows}.

In the above

.. math::

   \dot{M}|\dot{M}_Z|\dot{J}_\mathrm{outflow,return} = \alpha_\mathrm{outflow~return~rate} {M|M_Z|J_\mathrm{outflowed}\over \tau_\mathrm{dynamical, halo}},

where :math:`\alpha_\mathrm{outflow~return~rate}=(`\ ``hotHaloOutflowReturnRate``) is an input parameter controlling the rate at which gas flows from the outflowed to hot reservoirs, and :math:`\dot{M}|\dot{M}_Z|\dot{J}_\mathrm{outflows}` are the net rates of outflow from any components in the node.

In the above, :math:`f_\mathrm{outflow, stripped}` is the fraction of outflowing material assumed to be stripped from the halo. The is computed following the algorithm of :cite:t:`font_colours_2008`, namely

.. math::

   f_\mathrm{outflow, stripped} = \epsilon_\mathrm{strip} {M_\mathrm{hot, outer} \over  M_\mathrm{hot, virial}},

where :math:`\epsilon_\mathrm{strip}=`\ ``[hotHaloOutflowStrippingEfficiency]`` is an input parameter, :math:`M_\mathrm{hot, outer}` is the mass of hot gas contained within the outer radius of the hot halo and :math:`M_\mathrm{hot, virial}` is the mass of hot gas that would be present if the hot halo extended to the virial radius (i.e. if no stripping had occurred).

A fraction :math:`1-`\ ``[hotHaloAngularMomentumLossFraction]`` of the cooling angular momentum rate, :math:`\dot{M}_\mathrm{cooling} r_\mathrm{cool} V_\mathrm{rotate}`, is sent through the ``Hot_Halo_Cooling_Angular_Momentum`` pipe.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* If the ``starveSatellites`` parameter is true, then any hot halo properties of the minor :term:`node` are added to those of the major :term:`node` and the hot halo :term:`component` removed from the minor node. Additionally in this case, any material outflowed or stripped from the the satellite galaxy to its hot halo is transferred to the hot halo of the host dark matter halo after each timestep. (Alternatively, if ``starveSatellitesOutflowed``\ :math:`=`\ ``true`` then only the outflowed and stripped gas is transferred to the host halo---the main hot gas reservoir is left in place.) If stripped mass is being tracked (i.e. if ``[hotHaloTrackStrippedGas]``\ :math:`=`\ ``true``) then any stripped mass is transferred from the satellite galaxy to the hot halo of the host dark matter halo after each timestep. If ``[hotHaloNodeMergerLimitBaryonFraction]``\ :math:`=`\ ``true`` then the hot gas content of the merged node is limited such that the total baryon content of the node (including satellites) does not exceed the universal baryon fraction, if possible. Any gas removed to enforce this limit is placed into the unaccreted gas reservoir, from which is may eventually be reaccreted.\\

*Satellite merging:* If the ``starveSatellites`` parameter is false, then any hot halo properties of the satellite :term:`node` are added to those of the host :term:`node` and the hot halo :term:`component` removed from the satellite node.\\

*Node promotion:* Any hot halo properties of the parent :term:`node` are added to those of the :term:`node` prior to promotion.\\

*Halo formation:* If ``[hotHaloOutflowReturnOnFormation]``\ :math:`=`\ ``true`` then all outflowed gas is returned to the hot gas reservoir on `halo formation events <https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Development.pdf\#sec.HaloFormationEvents>`_.\\

"Outflow Tracking" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The outflow tracking hot halo implementation extends the "standard" implementation by adding the following properties:

``trackedOutflowMass``
   The mass of gas in the hot halo which arrived there directly via outflow: :math:`M_\mathrm{outflow, track}` ``[hotHaloTrackedOutflowMass]``.

``trackedOutflowAbundances``
   The mass of elements in the hot halo which arrived there directly via outflow: :math:`M_{Z, \mathrm{outflow, track}}` ``[hotHaloTrackedOutflowAbundances]``.

Initialization
^^^^^^^^^^^^^^

Outflowed masses and element masses are initialized to zero.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

The tracked outflow masses evolve according to:

.. math::

   \dot{M}|\dot{M}_{Z_\mathrm{outflow, track}} &=& \alpha_\mathrm{outflow~return~rate} {M|M_{Z_\mathrm{outflowed}}\over \tau_\mathrm{dynamical, halo}} -  M|M_{Z_\mathrm{outflow, track}} \dot{M}_\mathrm{expelled} / M \\

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None.\\

*Satellite merging:* None.\\

*Node promotion:* None.\\

*Halo formation:* None.\\

Galactic Disk
-------------

"Very Simple" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This implementation assumes a disk with no structural properties---it consists of just gas and stellar masses.

Properties
^^^^^^^^^^

The very simple galactic disk implementation defines the following properties:

``massGas``
   The mass of gas in the disk: :math:`M_\mathrm{disk, gas}` [``diskMassGas``];

``massStellar``
   The mass of stars in the disk: :math:`M_\mathrm{disk, stars}` [``diskMassStellar``].

Initialization
^^^^^^^^^^^^^^

No initialization is performed---disks are created as needed.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

In the very simple galactic disk implementation the gas mass evolves as:

.. math::

   \dot{M}_\mathrm{disk, gas} = \dot{M}_\mathrm{cooling} - \dot{M}_\mathrm{outflow, disk} - \dot{M}_\mathrm{stars, disk},

where the rate of change of stellar mass is

.. math::

   \dot{M}_\mathrm{stars, disk} = \Psi - \dot{R},

where :math:`\dot{R}` is the rate of mass recycling from stars, and

.. math::

   \Psi = {M_\mathrm{disk, gas} \over \tau_\mathrm{disk, star~formation}}

with :math:`\tau_\mathrm{disk, star~formation}` being the greater of the star formation timescale and :math:`\Gamma_\mathrm{disk, star formation, minimum} \tau_\mathrm{dyn}`, where :math:`\tau_\mathrm{dyn}` is the dynamical time of the halo, and :math:`\Gamma_\mathrm{disk, star formation, minimum}=`\ ``[diskStarFormationTimescaleMinimum]``. The outflow rate, :math:`\dot{M}_\mathrm{outflow, disk}`, is computed for the current star formation rate and gas properties by the prescriptions for non-expulsive supernova feedback (see :galacticus-class:`stellarFeedbackOutflows`), but is limited to a maximum of :math:`M_\mathrm{disk, gas}/ \Gamma_\mathrm{disk, outflow, minimum} \tau_\mathrm{dyn}`, where :math:`\Gamma_\mathrm{disk, outflow, minimum}=`\ ``[diskOutflowTimescaleMinimum]``. This outflow is piped to the hot halo component.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None\\

*Satellite merging:* Disks may be destroyed (or, potentially, created or otherwise modified) as the result of a satellite merging event, as dictated by the selected merger remnant mass movement method (see :galacticus-class:`mergerMassMovements`).\\

*Node promotion:* None\\

"Very Simple Size" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This implementation extends the ``verySimple`` disk class by adding a half-mass radius property. No other assumptions about the mass distribution are made.

Properties
^^^^^^^^^^

The very simple size galactic disk implementation defines the following additional properties:

``radius``
   The half-mass radius of the disk: :math:`r_\mathrm{disk, 1/2}` [``diskRadius``].

Initialization
^^^^^^^^^^^^^^

No initialization is performed---disks are created as needed.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

N/A---disk radii are computed using the selected galactic structure solver.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None\\

*Satellite merging:* None\\

*Node promotion:* None\\

.. _manual-sec-DiskStandard:

"Standard" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~

This implementation assumes a disk with a cylindrical symmetry and a flattened profile in which stars trace gas. Currently, one option for the density profile is allowed\ [#]_:

``exponentialDisk``
   Assumes an exponential profile, :math:`\rho(R,z) = \rho_0 \exp(-r/r_\mathrm{s}) \hbox{sech}^2(z/z_\mathrm{s})`, for the disk :term:`component` of a galaxy. The thickness, :math:`z_\mathrm{s}/r_\mathrm{s}` is set by the parameter ``[heightToRadialScaleDisk]``.

Properties
^^^^^^^^^^

The standard galactic disk implementation defines the following properties:

``massGas``
   The mass of gas in the disk: :math:`M_\mathrm{disk, gas}` [``diskMassGas``].

``abundancesGas``
   The mass of elements in the gaseous disk: :math:`M_{Z, \mathrm{disk, gas}}` [``diskAbundancesGas{abundanceName}``].

``massStellar``
   The mass of stars in the disk: :math:`M_\mathrm{disk, stars}` [``diskMassStellar``].

``abundancesStellar``
   The mass of elements in the stellar disk: :math:`M_{Z, \mathrm{disk, stars}}` [``diskAbundancesStellar{abundanceName}``].

``luminositiesStellar``
   The luminosities (in multiple bands) of the stellar disk: :math:`L_\mathrm{disk, stars}` [``diskLuminositiesStellar{luminosityName}``].

``angularMomentum``
   The angular momentum of the disk, :math:`J_\mathrm{disk}` [``diskAngularMomentum``].

``radius``
   The radial scale length of the disk, :math:`R_\mathrm{disk}` [``diskRadius``].

``velocity``
   The circular velocity of the disk at :math:`R_\mathrm{disk}`, :math:`V_\mathrm{disk}` [``diskVelocity``].

``massStellarFormed``
   The total mass of stars ever formed in the disk: :math:`M_\mathrm{disk, stars, formed}` .

``fractionMassRetained``
   The fraction of mass retained in the disk as a result of transfer processes in which the transfer is proportional to the current mass: :math:`f_\mathrm{retained}`.

Initialization
^^^^^^^^^^^^^^

No initialization is performed---disks are created as needed.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

In the standard galactic disk implementation the gas mass evolves as:

.. math::

   \dot{M}_\mathrm{disk, gas} = \dot{M}_\mathrm{cooling} - \dot{M}_\mathrm{outflow, disk} - \dot{M}_\mathrm{stars, disk} - {M_\mathrm{disk, gas}\over \tau_\mathrm{bar}} - \dot{M}_\mathrm{ram pressure} - {M_\mathrm{disk, gas} \over M_\mathrm{disk, gas} + M_\mathrm{disk, stars}} \dot{M}_\mathrm{tidal},

where the rate of change of stellar mass is

.. math::

   \dot{M}_\mathrm{stars, disk} &=& \Psi - \dot{R} - {M_\mathrm{stars,disk} \over \tau_\mathrm{bar}} - {M_\mathrm{disk, stars} \over M_\mathrm{disk, gas} + M_\mathrm{disk, stars}} \dot{M}_\mathrm{tidal}, \\
   \dot{M}_\mathrm{stars, disk, formed } &=& \Psi,

with

.. math::

   \Psi = {M_\mathrm{disk, gas} \over \tau_\mathrm{disk, star~formation}}

with :math:`\tau_\mathrm{disk, star~formation}` being the star formation timescale and :math:`\dot{R}` is the rate of mass recycling from stars and :math:`\tau_\mathrm{bar}` is a bar instability timescale (see :galacticus-class:`galacticDynamicsBarInstability`). If ``[diskStarFormationInSatellites]``\ :math:`=`\ ``true`` then the star formation rate is forced to zero in satellite galaxies. The mass removed from the disk by the bar instability mechanism is added to the active spheroid component. Element abundances (including total metals) evolve according to:

.. math::

   \dot{M}_{Z, \mathrm{disk, gas}} = \dot{M}_{Z \mathrm{cooling}} - \dot{M}_{Z, \mathrm{outflow, disk}} - \dot{M}_{Z, \mathrm{stars, disk}} + \dot{y} - {M_{Z, \mathrm{disk, gas}} \over M_\mathrm{disk, gas}} \dot{M}_\mathrm{ram pressure} - {M_{Z, \mathrm{disk, gas}} \over M_\mathrm{disk, gas} + M_\mathrm{disk, stars}} \dot{M}_\mathrm{tidal},

and

.. math::

   \dot{M}_{Z, \mathrm{stars, disk}} = \Psi {M_{Z, \mathrm{disk, gas}} \over M_\mathrm{disk, gas}} - \dot{R}_Z - {M_{Z, \mathrm{disk, stars}} \over M_\mathrm{disk, gas} + M_\mathrm{disk, stars}} \dot{M}_\mathrm{tidal}

where :math:`\dot{y}` is the rate of element yield from stars and :math:`\dot{R}_Z` is the rate of element recycling. The angular momentum evolves as:

.. math::

   \dot{J}_\mathrm{disk} = \dot{J}_\mathrm{cooling} - \left[ \dot{M}_\mathrm{outflow, disk} + {M_\mathrm{disk, gas}  + M_\mathrm{disk, stars} \over \tau_\mathrm{bar}} + \dot{M}_\mathrm{ram pressure} + \dot{M}_\mathrm{tidal}\right] {J_\mathrm{disk} \over M_\mathrm{disk, gas}}.

The outflow rate, :math:`\dot{M}_\mathrm{outflow, disk}`, is computed for the current star formation rate and gas properties by the stellar properties subsystem (see :galacticus-class:`stellarPopulationProperties`) and prescriptions for expulsive and non-expulsive supernova feedback (see :galacticus-class:`stellarFeedbackOutflows` and :galacticus-class:`stellarFeedbackOutflows` respectively), but is not allowed to exceed :math:`M_\mathrm{gas, disk}/ \alpha_\mathrm{outflow minimum, disk} \tau_\mathrm{disk, dynamical}`, where :math:`\tau_\mathrm{disk, dynamical}=R_\mathrm{disk}/V_\mathrm{disk}` is the dynamical time of the disk and :math:`\alpha_\mathrm{outflow minimum, disk}=`\ ``[diskOutflowTimescaleMinimum]`` is the shortest timescale (in units of the dynamical timescale) on which gas can be removed from the disk. This limit prevents the disk being depleted on arbitrarily short timescales. The non-expulsive :term:`component` of the outflow is piped to the hot halo component.  The ram pressure and tidal mass loss rates, :math:`\dot{M}_\mathrm{ram pressure}` and :math:`\dot{M}_\mathrm{tidal}`, are computed using the selected methods (see :galacticus-class:`hotHaloRamPressureTimescale` and :galacticus-class:`tidalStripping` respectively). Finally, stellar luminosities evolve according to:

.. math::

   \dot{L}_\mathrm{disk, stars} = \Psi \mathcal{L}_\lambda(t_0-t,Z_\mathrm{gas})  - {L_\mathrm{disk, stars} \over \tau_\mathrm{bar}} - {L_\mathrm{disk, stars} \over M_\mathrm{disk, gas} + M_\mathrm{disk, stars}} \dot{M}_\mathrm{tidal},

where :math:`\mathcal{L}_\lambda(t,Z)` is the luminosity-to-mass ratio in the relevant band for a stellar population of age :math:`t` and metallicity :math:`Z`. The fraction of mass retained in the disk, :math:`f_\mathrm{retained}`, is given by:

.. math::

   \dot{f}_\mathrm{retained} = - {f_\mathrm{retained} \over \tau_\mathrm{bar}},

where the initial value of :math:`f_\mathrm{retained}` is arbitrary as we will be interested only in ratios of the quantity.

Integral Evolution
^^^^^^^^^^^^^^^^^^

The standard disk component supports solving for the stellar luminosities as inactive variables, if the parameter ``[diskLuminositiesStellarInactive]``\ :math:`=`\ ``true``. Across a timestep :math:`t_i` to :math:`t_{i+1}` the change in luminosity can be written as the following integral:

.. math::

   \Delta L_\mathrm{disk, stars} &=& \int_{t_{i}}^{t_{i+1}} \mathrm{d}t \Psi \mathcal{L}_\lambda(t_0-t,Z_\mathrm{gas}) {f_\mathrm{retained}(t_{i+1}) \over f_\mathrm{retained}(t)}  + L_\mathrm{disk, stars}(t_i) {\dot{f}_\mathrm{retained}(t) \over f_\mathrm{retained}(t_{i})}, \\
   \Delta L_\mathrm{spheroid, stars} &=& \int_{t_{i}}^{t_{i+1}} \mathrm{d}t  \Psi \mathcal{L}_\lambda(t_0-t,Z_\mathrm{gas}) \left( 1 - {f_\mathrm{retained}(t_{i+1}) \over f_\mathrm{retained}(t)} \right) - L_\mathrm{disk, stars}(t_i) {\dot{f}_\mathrm{retained}(t) \over f_\mathrm{retained}(t_{i})}. \\

The first term in the integral for disk stellar luminosity is just the usual production of starlight due to star formation, reduced by a factor :math:`f_\mathrm{retained}(t_{i+1}) / f_\mathrm{retained}(t)` to account for that which will be retained in the disk at the end of the timestep. The second term accounts for transfer of starlight produced at :math:`t<t_i`, since it integrates to :math:`L_\mathrm{disk, stars}(t_i) (f_\mathrm{retained}(t_{i+1}) / f_\mathrm{retained}(t_{i}) - 1)` which is the change in the starlight of the disk over the timestep. The functions :math:`f_\mathrm{retained}(t)` and :math:`M_\mathrm{stars, disk, formed}(t)` have already been solved for during the differential evolution phase. Therefore, :math:`f_\mathrm{retained}(t)` and :math:`\dot{f}_\mathrm{retained}(t)` are available for this calculation. Furthermore, we can make the replacement :math:`\Psi = \dot{M}_\mathrm{stars, disk, formed}(t)` in the above, allowing us to avoid recomputing :math:`\Psi` directly in these integral. The integral for the spheroid component simply reflects that starlight not retained in the disk is transferred to the spheroid.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None\\

*Satellite merging:* Disks may be destroyed (or, potentially, created or otherwise modified) as the result of a satellite merging event, as dictated by the selected merger remnant mass movement method (see :galacticus-class:`mergerMassMovements`).\\

*Node promotion:* None\\

Additional Output
^^^^^^^^^^^^^^^^^

If the ``[diskOutputStarFormationRate]`` input parameter is set to true, then the instantaneous star formation rate in the disk (in units of :math:`M_\odot` Gyr\ :math:`^{-1}`) will be included in the output, as ``diskStarFormationRate``.

Structure
^^^^^^^^^

The radial size of the disk is found solving for equilibrium (i.e. the radius is such that the angular momentum of material at that radius is sufficient to provide rotational support) at the specified ``[diskStructureSolverRadius]`` which is given in units of the disk scale length. In converting from the mean specific angular momentum of the disk to the angular momentum at that radius, a flat rotation curve is assumed, i.e.:

.. math::

   j(r)/\langle j \rangle = r V \left/ {\int_0^\infty 2 \pi r^\prime \Sigma(r^\prime) r^\prime V \mathrm{d} r^\prime \over \int_0^\infty 2 \pi r^\prime \Sigma(r^\prime) \mathrm{d} r^\prime} \right. .

The option ``[diskRadiusSolverCole2000Method]``, if set to ``true``, alters this behavior to match that of the structure solver used by :cite:t:`cole_hierarchical_2000`, in which adiabatic contraction of the dark matter halo is solved for assuming that the disk has a spherical mass distribution. The specific angular momentum passed to the structure solver will be modified as follows in this case:

.. math::

   j(r) \rightarrow \left[ j^2(r) - \left( V_\mathrm{disk}^2(r) r^2 - \mathrm{G} M_\mathrm{disk}(<r) r \right) \right]^{1/2},

where :math:`V_\mathrm{disk}` is the rotation curve in the plane of the disk. This adjustment accounts for the difference between a thin disk and spherical mass distribution. Note that in this case (as in :cite:author:`cole_hierarchical_2000` :cite:year:`cole_hierarchical_2000`) the resulting disk will not precisely satisfy :math:`j(r) = r V_\mathrm{c}(r)` where :math:`V_\mathrm{c}(r)` is the net rotation curve.

.. _manual-sec-ComponentSpheroid:

Galactic Spheroid
-----------------

"Very Simple" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This implementation assumes a spheroid consisting of just gas and stellar masses, plus a half-mass radius (with no other assumption about structure being made).

Properties
^^^^^^^^^^

The very simple galactic spheroid implementation defines the following properties:

``massGas``
   The mass of gas in the spheroid: :math:`M_\mathrm{spheroid, gas}` [``spheroidMassGas``];

``massStellar``
   The mass of stars in the spheroid: :math:`M_\mathrm{spheroid, stars}` [``spheroidMassStellar``].;

``radius``
   The half-mass radius of the spheroid: :math:`r_\mathrm{spheroid, 1/2}` [``spheroidRadius``].

Initialization
^^^^^^^^^^^^^^

No initialization is performed---spheroids are created as needed.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

Spheroid radii are computed using the selected galactic structure solver. In the very simple galactic spheroid implementation the gas mass evolves as:

.. math::

   \dot{M}_\mathrm{spheroid, gas} = \dot{M}_\mathrm{cooling} - \dot{M}_\mathrm{outflow, spheroid} - \dot{M}_\mathrm{stars, spheroid},

where the rate of change of stellar mass is

.. math::

   \dot{M}_\mathrm{stars, spheroid} = \Psi - \dot{R},

where :math:`\dot{R}` is the rate of mass recycling from stars, and

.. math::

   \Psi = {M_\mathrm{spheroid, gas} \over \tau_\mathrm{spheroid, star~formation}}

with :math:`\tau_\mathrm{spheroid, star~formation}` being the greater of the star formation timescale and :math:`\Gamma_\mathrm{spheroid, star formation, minimum} \tau_\mathrm{dyn}`, where :math:`\tau_\mathrm{dyn}` is the dynamical time of the halo, and :math:`\Gamma_\mathrm{spheroid, star formation, minimum}=`\ ``[spheroidStarFormationTimescaleMinimum]``. The outflow rate, :math:`\dot{M}_\mathrm{outflow, spheroid}`, is computed for the current star formation rate and gas properties by the prescriptions for non-expulsive supernova feedback (see :galacticus-class:`stellarFeedbackOutflows`), but is limited to a maximum of :math:`M_\mathrm{spheroid, gas}/ \Gamma_\mathrm{spheroid, outflow, minimum} \tau_\mathrm{dyn}`, where :math:`\Gamma_\mathrm{spheroid, outflow, minimum}=`\ ``[spheroidOutflowTimescaleMinimum]``. This outflow is piped to the hot halo component.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None\\

*Satellite merging:* Spheroids may be destroyed (or, potentially, created or otherwise modified) as the result of a satellite merging event, as dictated by the selected merger remnant mass movement method (see :galacticus-class:`mergerMassMovements`).\\

*Node promotion:* None\\

"Standard" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~

The standard spheroid implementation assumes a spheroid density profile described by a single length scale in which stars trace gas. Currently, two options for the density profile are allowed\ [#]_:

``hernquist``
   Assumes a Hernquist profile :cite:p:`hernquist_analytical_1990` for the spheroidal :term:`component` of a galaxy.

``sersic``
   Assumes a Sérsic profile (:cite:author:`sersic_influence_1963` :cite:year:`sersic_influence_1963`; see also :cite:author:`mazure_exact_2002` :cite:year:`mazure_exact_2002`) for the spheroidal :term:`component` of a galaxy in which stars trace gas. The projected density profile of the spheroid is given by:

   .. math::

      \Sigma(R) \propto \exp\left(-b_\mathrm{n} R^{1/n} \right),

   where the Sérsic index, :math:`n=`\ ``[spheroidSersicIndex]`` and the coefficient :math:`b_\mathrm{n}=2.303(0.8689 n-0.1447)` :cite:t:`wadadekar_two-dimensional_1999`. The 3D density distribution for a given :math:`n` is inferred by solving the relevant inverse Abel integral.

Properties
^^^^^^^^^^

The standard galactic spheroid implementation defines the following properties:

``masGas``
   The mass of gas in the spheroid: :math:`M_\mathrm{spheroid, gas}` [``spheroidMassGas``].

``abundancesGas``
   The mass of elements in the gaseous spheroid: :math:`M_{Z, \mathrm{spheroid, gas}}` [``spheroidAbundancesGas{abundanceName}``].

``massStellar``
   The mass of stars in the spheroid: :math:`M_\mathrm{spheroid, stars}` [``spheroidMassStellar``].

``abundancesStellar``
   The mass of elements in the stellar spheroid: :math:`M_{Z, \mathrm{spheroid, stars}}` [``spheroidAbundancesStellar{abundanceName}``].

``luminositiesStellar``
   The luminosities (in multiple bands) of the stellar spheroid: :math:`L_\mathrm{spheroid, stars}` [``spheroidLuminositiesStellar{luminosityName}``].

``angularMomentum``
   The pseudo-angular momentum\ [#]_ of the spheroid, :math:`J_\mathrm{spheroid}` [``spheroidAngularMomentum``]. The parameter ``[spheroidAngularMomentumAtScaleRadius]`` controls the ratio of the specific pseudo-angular momentum at the scale radius of the standard spheroid to the mean specific pseudo-angular momentum. By default, this parameter is set to ``[spheroidAngularMomentumAtScaleRadius]`` :math:`= I_2/I_3`, where

   .. math::

      I_n = \int_0^\infty \rho(r) r^n \mathrm{d}r,

   and :math:`\rho(r)` is the spheroid density profile, which is appropriate for a flat rotation curve. In some cases (e.g. the Hernquist profile) one or both of :math:`I_2` and :math:`I_3` can be infinite. In such cases ``[spheroidAngularMomentumAtScaleRadius]`` :math:`=0.5` is assumed by default. If a finite truncation radius is assumed, or a different rotation curve is assumed, this ratio may be finite. The ``[spheroidAngularMomentumAtScaleRadius]`` parameter allows control over these assumptions.

``radius``
   The radial scale length of the spheroid, :math:`r_\mathrm{spheroid}` [``spheroidRadius``].

``velocity``
   The circular velocity of the spheroid at :math:`r_\mathrm{spheroid}`, :math:`V_\mathrm{spheroid}` [``spheroidVelocity``].

``massStellarFormed``
   The total mass of stars ever formed in the spheroid: :math:`M_\mathrm{spheroid, stars, formed}`.

and the following pipes:

``energyInput``
   Energy sent through this pipe is added to the gas of the spheroid and will result in an outflow (see below). Input energy should be in units of :math:`M_\odot` km\ :math:`^2` s\ :math:`^{-2}` Gyr\ :math:`^{-1}` and must be positive (energy cannot be removed from the gas via this pipe).

``massGasSink``
   Removes gas (and proportionate amounts of angular momentum and elements) from the spheroid gas. Removed mass should be in units of :math:`M_\odot` and must be positive (a negative mass sink would add mass to the spheroid which is not allowed via this pipe).

Initialization
^^^^^^^^^^^^^^

No initialization is performed---spheroids are created as needed.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

In the standard galactic spheroid implementation the gas mass evolves as\ [#]_:

.. math::

   \dot{M}_\mathrm{spheroid, gas} = - \dot{M}_\mathrm{outflow, spheroid} - \dot{M}_\mathrm{stars, spheroid} - \dot{M}_\mathrm{ram pressure} - {M_\mathrm{spheroid, gas} \over M_\mathrm{spheroid, gas} + M_\mathrm{spheroid, stars}} \dot{M}_\mathrm{tidal},

where the rate of change of stellar mass is

.. math::

   \dot{M}_\mathrm{stars, spheroid} &=& \Psi - \dot{R} - {M_\mathrm{spheroid, stars} \over M_\mathrm{spheroid, gas} + M_\mathrm{spheroid, stars}} \dot{M}_\mathrm{tidal}, \\
   \dot{M}_\mathrm{stars, disk, formed } &=& \Psi,

with

.. math::

   \Psi = {M_\mathrm{spheroid, gas} \over \tau_\mathrm{spheroid, star~formation}}

with :math:`\tau_\mathrm{spheroid, star~formation}` being the star formation timescale and :math:`\dot{R}` is the rate of mass recycling from stars. If ``[spheroidStarFormationInSatellites]``\ :math:`=`\ ``true`` then the star formation rate is forced to zero in satellite galaxies. Element abundances (including total metals) evolve according to:

.. math::

   \dot{M}_{Z, \mathrm{spheroid, gas}} &=& - \dot{M}_{Z, \mathrm{outflow, spheroid}} - \dot{M}_{Z, \mathrm{stars, spheroid}} + \dot{y} - {M_{Z, \mathrm{spheroid, gas}} \over M_\mathrm{spheroid, gas}} \dot{M}_\mathrm{ram pressure} \nonumber \\
   & & - {M_{Z, \mathrm{spheroid, gas}} \over M_\mathrm{spheroid, gas} + M_\mathrm{spheroid, stars}} \dot{M}_\mathrm{tidal},

and

.. math::

   \dot{M}_{Z, \mathrm{stars, spheroid}} = \Psi {M_{Z, \mathrm{spheroid, gas}} \over M_\mathrm{spheroid, gas}} - \dot{R}_Z - {M_{Z, \mathrm{spheroid, stars}} \over M_\mathrm{spheroid, gas} + M_\mathrm{spheroid, stars}} \dot{M}_\mathrm{tidal}

where :math:`\dot{y}` is the rate of element yield from stars and :math:`\dot{R}_Z` is the rate of element recycling. The angular momentum evolves as:

.. math::
   :label: eq-SpheroidStandardAngularMomentumEvolution

   \dot{J}_\mathrm{spheroid} = - (\dot{M}_\mathrm{outflow, spheroid} + \dot{M}_\mathrm{ram pressure}+\dot{M}_\mathrm{tidal}) {J_\mathrm{spheroid} \over M_\mathrm{spheroid, gas} + M_\mathrm{spheroid, stars}} + |\mathcal{T}| ( M_\mathrm{spheroid, gas} + M_\mathrm{spheroid, stars} ) R_\mathrm{spheroid}^2.

The outflow rate, :math:`\dot{M}_\mathrm{outflow, spheroid}`, is computed for the current star formation rate and gas properties by the stellar properties subsystem (see :galacticus-class:`stellarPopulationProperties`) and prescriptions for expulsive and non-expulsive supernova feedback (see :galacticus-class:`stellarFeedbackOutflows` and :galacticus-class:`stellarFeedbackOutflows` respectively), with an additional contribution given by

.. math::

   \dot{M}_\mathrm{outflow, spheroid} = \beta_\mathrm{spheroid, energy} {\dot{E}_\mathrm{gas, spheroid} \over V_\mathrm{spheroid}^2}

where :math:`\beta_\mathrm{spheroid, energy}=`\ ``[spheroidEnergeticOutflowMassRate]`` is an input parameter, and :math:`\dot{E}_\mathrm{gas,spheroid}` is any input energy sent through the ``Tree_Node_Spheroid_Gas_Energy_Input`` pipe, but is not allowed to exceed :math:`M_\mathrm{gas, spheroid}/ \alpha_\mathrm{outflow minimum, spheroid} \tau_\mathrm{spheroid, dynamical}`, where :math:`\tau_\mathrm{spheroid, dynamical}=R_\mathrm{spheroid}/V_\mathrm{spheroid}` is the dynamical time of the spheroid and :math:`\alpha_\mathrm{outflow minimum, spheroid}=`\ ``[spheroidOutflowTimescaleMinimum]`` is the shortest timescale (in units of the dynamical timescale) on which gas can be removed from the spheroid. This limit prevents the spheroid being depleted on arbitrarily short timescales. The non-expulsive :term:`component` of the outflow is piped to the hot halo component. The ram pressure and tidal mass loss rates, :math:`\dot{M}_\mathrm{ram pressure}` and :math:`\dot{M}_\mathrm{tidal}`, are computed using the selected methods (see :galacticus-class:`hotHaloRamPressureTimescale` and :galacticus-class:`tidalStripping` respectively). The final term in :eq:`eq-SpheroidStandardAngularMomentumEvolution` accounts for tidal heating of the spheroid due to the tidal field, :math:`\mathcal{T}`.

Finally, stellar luminosities evolve according to:

.. math::

   \dot{L}_\mathrm{spheroid, stars} = \Psi \mathcal{L}_\lambda(t_0-t,Z_\mathrm{gas}) - {L_\mathrm{spheroid, stars} \over M_\mathrm{spheroid, gas} + M_\mathrm{spheroid, stars}} \dot{M}_\mathrm{tidal},

where :math:`\mathcal{L}_\lambda(t,Z)` is the luminosity-to-mass ratio in the relevant band for a stellar population of age :math:`t` and metallicity :math:`Z`.

Integral Evolution
^^^^^^^^^^^^^^^^^^

The standard spheroid component supports solving for the stellar luminosities as inactive variables, if the parameter ``[spheroidLuminositiesStellarInactive]``\ :math:`=`\ ``true``. Across a timestep :math:`t_i` to :math:`t_{i+1}` the change in luminosity can be written as the following integral:

.. math::

   \Delta L_\mathrm{spheroid, stars} = \int_{t_{i}}^{t_{i+1}} \mathrm{d}t \Psi \mathcal{L}_\lambda(t_0-t,Z_\mathrm{gas}),

which is just the usual production of starlight due to star formation. As the function :math:`M_\mathrm{stars, spheroid, formed}(t)` has already been solved for during the differential evolution phase we can make the replacement :math:`\Psi = \dot{M}_\mathrm{stars, disk, formed}(t)` in the above, allowing us to avoid recomputing :math:`\Psi` directly in these integral.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None\\

*Satellite merging:* Spheroids may be created as the result of a satellite merging event, as dictated by the selected merger remnant mass movement model (see :galacticus-class:`mergerMassMovements`).\\

*Node promotion:* None.\\

Additional Output
^^^^^^^^^^^^^^^^^

If the ``[spheroidOutputStarFormationRate]`` input parameter is set to true, then the instantaneous star formation rate in the spheroid (in units of :math:`M_\odot` Gyr\ :math:`^{-1}`) will be included in the output, as ``spheroidStarFormationRate``.

.. _manual-sec-ComponentBasicProperties:

Basic Properties
----------------

Basic properties are the total mass of a :term:`node` and the cosmic time at which it currently exists.

"Non-evolving" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The non-evolving basic properties implementation defines the following properties:

``mass``
   The total mass of the node: :math:`M_\mathrm{node}` [``basicMass``].

``time``
   The time at which the :term:`node` is defined: :math:`t_\mathrm{node}`.

``timeLastIsolated``
   The time at which the :term:`node` was last an isolated halo (i.e. not a subhalo): [``basicTimeLastIsolated``].

Initialization
^^^^^^^^^^^^^^

All basic properties are required to be initialized by the merger tree construction routine.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

Properties are evolved according to:

.. math::

   \dot{M}_\mathrm{node} &=& 0 \\
   \dot{t}_\mathrm{node} &=& 1.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None.\\

*Satellite merging:* None.\\

*Node promotion:* :math:`M_\mathrm{node}` is updated to the :term:`node` mass of the parent prior to promotion.\\

"Standard Implementation
~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The standard basic properties implementation defines the following properties:

``mass``
   The total mass of the node: :math:`M_\mathrm{node}` [``basicMass``].

``time``
   The time at which the :term:`node` is defined: :math:`t_\mathrm{node}`.

``timeLastIsolated``
   The time at which the :term:`node` was last an isolated halo (i.e. not a subhalo): [``basicTimeLastIsolated``].

Initialization
^^^^^^^^^^^^^^

All basic properties are required to be initialized by the merger tree construction routine.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

Properties are evolved according to:

.. math::

   \dot{M}_\mathrm{node} &=& \left\{\begin{array}{ll}{M_\mathrm{node, parent} - M_\mathrm{node} \over t_\mathrm{node, parent} - t_\mathrm{node}} & \hbox{ if primary progenitor} \\
   0 & \hbox{ otherwise}, \end{array} \right. \\
   \dot{t}_\mathrm{node} &=& 1,

where the "parent" subscript indicates a property of the parent :term:`node` in the merger tree.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None.\\

*Satellite merging:* None.\\

*Node promotion:* :math:`M_\mathrm{node}` is updated to the :term:`node` mass of the parent prior to promotion.\\

"Standard-Extended Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The standard-extended basic properties implementation extends the standard implementation and defines the following additional properties:

``massBertschinger``
   The "Bertschinger" mass of the node, defined as the mass enclosing a density contrast equal to that for the spherical collapse model: :math:`M_\mathrm{node, Bertschinger}` [``basicMassBertschinger``].

``accretionRateBertschinger``
   The accretion rate of "Bertschinger" mass onto the node: :math:`\dot{M}_\mathrm{node, Bertschinger}` [``basicAccretionRateBertschinger``].

``radiusTurnaround``
   The turnaround radius corresponding to the current "Bertschinger"mass: :math:`R_\mathrm{ta}` [``basicRadiusTurnaround``]

Initialization
^^^^^^^^^^^^^^

The Bertschinger mass and accretion rate for each node are computed from :math:`M_\mathrm{basic}` and the dark matter density profile. The turnaround radius is computed from the virial radius (under the appropriate definition for the Bertschinger mass) and the ratio of turnaround to virial radius from spherical collapse models. The specific spherical collapse model to use is determined by the ``[nodeComponentBasicExtendedSphericalCollapseType]`` parameter---a value of "``matterLambda``" causes the solution for matter+cosmological constant universes to be used, while a value of "``matterDarkEnergy``" causes the solution for matter+dark energy universes to be used.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

None.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* The accretion rate of Bertschinger mass is set to zero.\\

*Satellite merging:* None.\\

*Node promotion:* :math:`\dot{M}_\mathrm{node, Bertschinger}` is updated to the :term:`node` Bertschinger accretion rate of the parent prior to promotion.\\

.. _manual-sec-ComponentPosition:

Position
--------

The position :term:`component` implements the position and velocity of each galaxy. See Section :galacticus-ref:`GalacticusVelocityDefinitions` for important notes on velocity definitions in Galacticus.

"cartesian" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The Cartesian position implementation defines the following properties:

``position``
   The 3-D position of the node: :math:`\mathbf{x}` [``positionPosition[X|Y|Z]``].

``velocity``
   The 3-D velocity of the node: :math:`\mathbf{v}` [``positionVelocity[X|Y|Z]``].

``positionHistory``
   The history of the node's position in 6-D phase space, usually used for satellite nodes.

Initialization
^^^^^^^^^^^^^^

None---all properties are assumed to have been preset, usually by the merger tree construction routine.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

None. Positions and velocities do not evolve for a given node. When output, if a 6-D position history is available than the position and velocity from the history entry closest to the output time will be used\ [#]_.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* If ``positionsPresetSatelliteToHost``\ :math:`=`\ ``true`` then the position and velocity of the satellite node is set equal to that of the host node (this is useful if position data is not available for orphaned halos for example, which would otherwise remain fixed in physical coordinates at their last known position---the position/velocity will also be updated to that of the new host each time a satellite's host changes), otherwise, none.\\

*Satellite merging:* None.\\

*Node promotion:* The position and velocity are updated to those of the parent node.\\

Satellite Orbit
---------------

This :term:`component` tracks the orbital properties of subhalos.

"Preset" Implementation
~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The preset satellite orbit implementation defines the following properties:

``mergeTime``
   The time until the satellite will merge with its host: :math:`t_\mathrm{satellite, merge}` [``satelliteMergeTime``].

``timeOfMerging``
   The cosmological time at which the satellite will merge with its host: :math:`T_\mathrm{satellite, merge}`.

``boundMass``
   The remaining, total bound mass of the satellite (this property is read only---it is determined from the ``boundMassHistory`` property).

``boundMassHistory``
   A history time-series of the total bound mass of the satellite.

``virialOrbit``
   The orbit (a `keplerOrbit <https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Development.pdf\#sec.KeplerOrbits>`_ object) of the satellite at virial orbit crossing.

Note that the ``mergeTime`` and ``timeOfMerging`` effectively provide the same information. For that reason, setting one of them will automatically set the other accordingly.

Initialization
^^^^^^^^^^^^^^

None. This method assumes that merging times and bound mass histories will be set externally (usually when the merger tree is constructed).

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

None.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None.\\

*Satellite merging:* None.\\

*Node promotion:* None.\\

"Merge Time" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The "merge time" satellite orbit implementation defines the following properties:

``timeOfMerging``
   The time at which the satellite will merge with its host: :math:`t_\mathrm{satellite, merge}`.

"Standard" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The standard satellite orbit implementation extends the "merge time" implementation and defines the following additional properties:

``boundMass``
   The remaining, total bound mass of the satellite: :math:`M_\mathrm{node,bound}` [``satelliteBoundMass``].

``virialOrbit``
   The orbit (returned as a `keplerOrbit <https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Development.pdf\#sec.KeplerOrbits>`_ object) of the satellite at the point of virial radius crossing.

.. _manual-sec-ComponentSatelliteOrbiting:

"Orbiting" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The orbiting satellite orbit implementation defines the following properties:

``position``
   The 3-dimensional position of the satellite relative to its host: :math:`\mathbf{r}=(x,y,z)`.

``velocity``
   The 3-dimensional velocity of the satellite relative to its host: :math:`\mathbf{v}=(v_x,v_y,v_z)`.

``mergeTime``
   The time until the satellite will merge with its host: :math:`t_\mathrm{satellite, merge}` [``satelliteMergeTime``].

``boundMass``
   The remaining, total bound mass of the satellite: :math:`M_\mathrm{node,bound}` [``satelliteBoundMass``].

``virialOrbit``
   The orbit (returned as a `keplerOrbit <https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Development.pdf\#sec.KeplerOrbits>`_ object) of the satellite at the point of virial radius crossing.

``tidalTensorPathIntegrated``
   The time integral of the tidal tensor along the orbit of the satellite from initialization: :math:`G_{ij}`.

``tidalHeatingNormalized``
   The tidal heating energy per radius squared deposited into the satellite: :math:`Q_\mathrm{tidal}`.

Initialization
^^^^^^^^^^^^^^

The satellite position, if not already assigned, is selected such that its magnitude is set according to its virial orbit parameters and its direction is chosen at random assuming an isotropic distribution.  The satellite velocity is similarly selected.  The integrated tidal tensor, :math:`G_{ij}`, is initialized to the null tensor, while :math:`Q_\mathrm{tidal}` is initialized to zero.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

Properties are evolved according to:

.. math::

   \dot{\mathbf{r}}&=&\mathbf{v} \\
   \dot{\mathbf{v}}&=&-\frac{G M_\mathrm{host}(<r)\mathbf{r}}{r^3}\left(1+\frac{M_\mathrm{node,bound}}{M_\mathrm{host}(<r)}\right)+\mathbf{a}_\mathrm{DF} \\
   \dot{G}_{ij}&=&g_{ij}-G_{ij}/T_\mathrm{orb},

with :math:`r=|\mathbf{r}|`, :math:`\mathbf{a}_{DF}` set to the rate given by the ``satelliteDynamicalFriction`` (see :galacticus-class:`satelliteDynamicalFriction`), :math:`g_{ij}` being the tidal tensor, :math:`T_\mathrm{orb}` being the satellite's orbital period, :math:`\dot{M}_\mathrm{node,bound}` set to the rate given by the ``satelliteTidalStripping`` (see :galacticus-class:`satelliteTidalStripping`), and :math:`\dot{Q}_\mathrm{tidal}` set to the rate given by the ``satelliteTidalHeating`` (see :galacticus-class:`satelliteTidalHeatingRate`). Note that in the evolution of :math:`G_{ij}` a decay term is artificially introduced such that the integration of :math:`g_{ij}` is effectively over just the previous orbit (i.e. over the previous tidal shock).

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* The :term:`component` is created and the time to merging is given a value of -1 to indicate that the satellite is not about to merge as unmerged.  Once the satellite satisfies one of the two merging conditions:

.. math::

   r&<&R_\mathrm{host}+R_\mathrm{satellite} \\
   M_\mathrm{node,bound}&<&f_\mathrm{d} M_\mathrm{node,basic},

with :math:`R` being the node's half-mass radius, :math:`M_\mathrm{node,basic}` being the node's initial mass, and :math:`f_\mathrm{d}=`\ ``[satelliteOrbitingDestructionMassFraction]``, the node is considered merged and the time to merging is set to zero. The bound mass is set to the current total mass of the node. A virial orbit is selected using the ``virialOrbits`` (see :galacticus-class:`virialOrbit`). \\

*Satellite merging:* None.\\

*Node promotion:* Not applicable (component only exists for satellite nodes).\\

.. _manual-sec-DarkMatterHaloSpinComponent:

Dark Matter Halo Spin
---------------------

"Random" Implementation
~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The random dark matter halo spin implementation defines the following properties:

``spin``
   The spin parameter of the halo: :math:`\lambda` [``spinSpin``].

Initialization
^^^^^^^^^^^^^^

The spin parameter of each node, if not already assigned, is selected at random from a distribution of spin parameters. This value is assigned to the earliest progenitor of the halo traced along its primary branch. The value is then propagated forward along the primary branch until the :term:`node` mass exceeds that of the :term:`node` for which the spin was selected by a factor of ``[randomSpinResetMassFactor]``, at which point a new spin is selected at random, and the process repeated until the end of the branch is reached.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

The spin parameter does not evolve.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None.\\

*Satellite merging:* None.\\

*Node promotion:* The spin is updated to equal that of the parent node. (The two will differ only if this is a case where the new halo :term:`node` was sufficiently more massive than the :term:`node` for which a spin was last selected that a new spin value was chosen.)\\

"Preset" Implementation
~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The preset dark matter halo spin implementation defines the following properties:

``spin``
   The spin parameter of the halo: :math:`\lambda` [``spinSpin``].

``spinGrowthRate``
   The growth rate spin parameter of the halo (in units of Gyr\ :math:`^{-1}`).

Initialization
^^^^^^^^^^^^^^

The spin parameter of each :term:`node` is assumed to have been preset prior to merger tree initialization. The growth rate is computed assuming linear growth with time along each branch.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

The spin parameter evolves linearly with time between :term:`node` and parent node.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None.\\

*Satellite merging:* None.\\

*Node promotion:* The spin and growth rate are updated to equal those of the parent node.\\

"Preset3D" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~

The type extends the ``preset`` spin component implementation, and so provides all properties of that implementation plus those described below.

Properties
^^^^^^^^^^

The preset dark matter halo spin implementation defines the following properties:

``spinVector``
   The spin vector of the halo: :math:`\lambda` [``spinSpinVector``].

``spinVectorGrowthRate``
   The growth rate of the spin vector of the halo (in units of Gyr\ :math:`^{-1}`).

Initialization
^^^^^^^^^^^^^^

The spin vector of each :term:`node` is assumed to have been preset prior to merger tree initialization. The growth rate is computed assuming linear growth with time along each branch.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

The spin parameter evolves linearly with time between :term:`node` and parent node.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None.\\

*Satellite merging:* None.\\

*Node promotion:* The spin vector and growth rate are updated to equal those of the parent node.\\

"Vitvitska Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The Vitvitska dark matter halo spin implementation follows the :cite:t:`vitvitska_origin_2002` model for halo spins, and defines the following properties:

``spin``
   The spin parameter of the halo: :math:`\lambda` [``spinSpin``].

Initialization
^^^^^^^^^^^^^^

The spin parameters of nodes having no children are drawn at random from the selected spin distribution function. For other nodes, their angular momentum is set equal to the sum of the spin angular momentum of their primary progenitor and the orbital angular momenta of all non-primary progenitors (which are determined from the orbital parameters) modulated by a factor dependent on the mass ratio of the merging halos:

.. math::

   \mathbf{J}^\mathrm{(spin)}_\mathrm{parent} = \mathbf{J}^\mathrm{(spin)}_{1} + \sum_{i=2}^N \left(1+{m_i \over m_{1}}\right)^{-\alpha} \mathbf{J}^\mathrm{(orbit)}_{i}.

where :math:`\alpha=`\ ``[spinVitvitskaMergerRatioExponent]``.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

The spin parameter does not evolve.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None.\\

*Satellite merging:* None.\\

*Node promotion:* None.

.. _manual-sec-DarkMatterProfileComponent:

Dark Matter Profile
-------------------

This :term:`component` stores dynamic properties associated with dark matter halo density profiles.

.. _manual-sec-DarkMatterProfileScale:

"Scale" Implementation
~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The scale dark matter profile implementation defines the following properties:

``scale``
   The scale length of the density profile [``darkMatterProfileScale``];

``scaleGrowthRate``
   The growth rate of the scale length of the density profile.

Initialization
^^^^^^^^^^^^^^

The scale length of each node, if not already assigned, is assigned using the concentration parameter function (see :galacticus-class:`darkMatterProfileConcentration`), but is not allowed to drop below ``[darkMatterProfileMinimumConcentration]``, such that the scale length is equal to the virial radius divided by that concentration.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

The scale radius grows linearly with time to interpolate between the scale radii of child and parent halos.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None.\\

*Satellite merging:* None.\\

*Node promotion:* None.\\

.. _manual-sec-DarkMatterProfileScalePreset:

"Scale Preset" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The "scale preset" dark matter profile implementation defines the following properties:

``scale``
   The scale length of the density profile [``darkMatterProfileScale``];

``scaleGrowthRate``
   The growth rate of the scale length of the density profile.

Initialization
^^^^^^^^^^^^^^

The scale length of each node is assumed to be assigned during tree construction.

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

The scale radius grows linearly with time to interpolate between the scale radii of child and parent halos.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None.\\

*Satellite merging:* None.\\

*Node promotion:* None.\\

.. _manual-sec-DarkMatterProfileScaleShape:

"Scale+Shape" Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Properties
^^^^^^^^^^

The scale\ :math:`+`\ shape dark matter profile implementation defines the following properties:

``scale``
   The scale length of the density profile [``darkMatterProfileScale``];

``scaleGrowthRate``
   The growth rate of the scale length of the density profile.

``shape``
   A shape parameter describing the density profile [``darkMatterProfileShape``];

``shapeGrowthRate``
   The growth rate of the shape parameter of the density profile.

Initialization
^^^^^^^^^^^^^^

The scale length of each node, if not already assigned, is assigned using the concentration parameter function (see :galacticus-class:`darkMatterProfileConcentration`), but is not allowed to drop below ``[darkMatterProfileMinimumConcentration]``, such that the scale length is equal to the virial radius divided by that concentration. The shape parameter of each :term:`node` is assigned using the dark matter profile shape model (see :galacticus-class:`darkMatterProfileShape`).

Differential Evolution
^^^^^^^^^^^^^^^^^^^^^^

The scale radius and shape parameter of each node grow linearly with time to interpolate between the scale radii and shape parameters of child and parent halos.

Event Evolution
^^^^^^^^^^^^^^^

*Node mergers:* None.\\

*Satellite merging:* None.\\

*Node promotion:* None.\\

.. [#] Specifically, the jet power multiplied by :math:`f_\mathrm{hot} [(M_\mathrm{hot}/M_\mathrm{total}) (\Omega_\mathrm{M}/\Omega_\mathrm{b})]^2` is added to the hot halo heating rate. The dependence on the gas fraction in the hot halo ensures that the heating rate goes smoothly to zero as the hot halo becomes depleted of gas.
.. [#] Technically of the black hole plus accretion disk system.
.. [#] Note that mass is not removed from the hot halo to compensate, since the accretion rate is independent of the hot halo mass this could lead to negative mass in the halo.
.. [#] The disk density distribution is handled internally using a :galacticus-class:`massDistributionClass` object. As such, any mass distribution implemented as an extension of the ``massDistribution`` class (and which is described by a single length scale) could be trivially added to the standard disk component.
.. [#] The spheroid density distribution is handled internally using a :galacticus-class:`massDistributionClass` object. As such, any mass distribution implemented as an extension of the ``massDistribution`` class (and which is described by a single length scale) could be trivially added to the standard spheroid component.
.. [#] Effectively the angular momentum that the spheroid would have, were it rotationally supported rather than pressure supported.
.. [#] There may be an additional contribution to the mass and angular momentum rates of change in the spheroid due to material transferred from the disk :term:`component` via the bar instability mechanism (see Section :galacticus-ref:`DiskStandard`). This is not included here as it is not intrinsic to this specific spheroid implementation---it is handled explicitly by the disk :term:`component` and so applies equally to any spheroid :term:`component` implementation.
.. [#] While interpolation could be used this is usually a bad idea. For nodes that are satellites in a halo for example, no simple interpolation algorithm can correctly account for the complex orbital dynamics by which the position and velocity is actually evolving.
