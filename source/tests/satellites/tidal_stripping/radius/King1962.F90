!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Contains a program to test the :cite:t:`king_structure_1962` satellite tidal radius against independently computed reference
  values.
  !!}

program Test_Satellite_Tidal_Stripping_Radius_King1962
  !!{RST
  Test the :galacticus-class:`satelliteTidalStrippingRadiusKing1962` class against reference values computed independently of
  Galacticus by the ``referenceKing1962.py`` script in the `galacticusDevTools
  <https://github.com/galacticusorg/galacticusDevTools>`_ repository.

  An NFW satellite orbits within an NFW host. For a tidal field which is stretching, the tidal radius, :math:`r_\mathrm{t}`, is
  the solution of

  .. math::

     \frac{\mathrm{G} M_\mathrm{sat}(<r_\mathrm{t})}{r_\mathrm{t}^3} = \gamma_\mathrm{c} \omega^2 - \frac{\mathrm{d}^2\Phi}{\mathrm{d}R^2},

  with :math:`\mathrm{d}^2\Phi/\mathrm{d}R^2 = 4 \pi \mathrm{G} \rho_\mathrm{host}(R) - 2 \mathrm{G} M_\mathrm{host}(<R)/R^3` for
  a spherical host, which reduces to the Jacobi radius for a point-mass host and satellite on a circular orbit. The cases tested
  cover tangential, radial, and general orbits, and several values of :math:`\gamma_\mathrm{c}`. A radial orbit and a
  tangential orbit with :math:`\gamma_\mathrm{c}=0` must give the same radius, since in both the centrifugal term vanishes.

  Where no tidal radius exists (here forced by using a null tidal field with :math:`\gamma_\mathrm{c}=0`) the radius enclosing
  the bound dark matter mass of the satellite, capped at its virial radius, must be returned. Bound masses below, between, and
  above the dark matter and total virial masses of the satellite are tested. Finally, a satellite with zero bound mass must have
  zero tidal radius.
  !!}
  use :: Cosmology_Functions            , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters           , only : cosmologyParametersSimple
  use :: Dark_Matter_Halo_Scales        , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Display                        , only : displayVerbositySet                               , verbosityLevelStandard
  use :: Events_Hooks                   , only : eventsHooksInitialize
  use :: Functions_Global_Utilities     , only : Functions_Global_Set
  use :: Galacticus_Nodes               , only : mergerTree                                        , nodeClassHierarchyFinalize         , nodeClassHierarchyInitialize, nodeComponentBasic          , &
       &                                         nodeComponentDarkMatterProfile                    , nodeComponentSatellite             , treeNode
  use :: Input_Parameters               , only : inputParameters
  use :: Node_Components                , only : Node_Components_Initialize                        , Node_Components_Thread_Initialize  , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Satellite_Tidal_Stripping_Radii, only : satelliteTidalStrippingRadiusKing1962
  use :: Satellites_Tidal_Fields        , only : satelliteTidalFieldNull                           , satelliteTidalFieldSphericalSymmetry
  use :: Unit_Tests                     , only : Assert                                            , Unit_Tests_Begin_Group             , Unit_Tests_End_Group        , Unit_Tests_Finish
  use :: Virial_Density_Contrast        , only : virialDensityContrastBryanNorman1998
  implicit none
  type            (inputParameters                                   )                 :: parameters
  type            (cosmologyParametersSimple                         ), pointer        :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                    ), pointer        :: cosmologyFunctions_
  type            (virialDensityContrastBryanNorman1998              ), pointer        :: virialDensityContrast_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition), pointer        :: darkMatterHaloScale_
  type            (satelliteTidalFieldSphericalSymmetry              ), pointer        :: satelliteTidalFieldSpherical_
  type            (satelliteTidalFieldNull                           ), pointer        :: satelliteTidalFieldNull_
  type            (satelliteTidalStrippingRadiusKing1962             ), pointer        :: radiusGammaUnity_                          , radiusGammaZero_                          , &
       &                                                                                  radiusGammaHalf_                           , radiusNoTide_
  type            (mergerTree                                        )                 :: tree
  type            (treeNode                                          ), pointer        :: node                                       , nodeHost
  class           (nodeComponentBasic                                ), pointer        :: basic                                      , basicHost
  class           (nodeComponentDarkMatterProfile                    ), pointer        :: darkMatterProfile                          , darkMatterProfileHost
  class           (nodeComponentSatellite                            ), pointer        :: satellite
  ! Halo properties. These, and the cosmology in the parameter file, must match those assumed by the reference script.
  double precision                                                    , parameter      :: massHost                   =1.0d12             , radiusScaleHost          =2.5d-2          , & ! [M☉], [Mpc]
       &                                                                                  massSatellite              =1.0d10             , radiusScaleSatellite     =4.0d-3            ! [M☉], [Mpc]
  ! Orbital positions [Mpc] and velocities [km/s] relative to the host, for each case.
  integer                                                             , parameter      :: countCases                 =5
  double precision                                                    , dimension(3,5) :: position                   =reshape([0.10d0,0.00d0, 0.00d0, &
       &                                                                                                                       0.10d0,0.00d0, 0.00d0, &
       &                                                                                                                       0.10d0,0.00d0, 0.00d0, &
       &                                                                                                                       0.03d0,0.00d0, 0.00d0, &
       &                                                                                                                       0.10d0,0.12d0,-0.05d0],[3,5]), &
       &                                                                                  velocity                   =reshape([  0.0d0,150.0d0,0.0d0, &
       &                                                                                                                     150.0d0,  0.0d0,0.0d0, &
       &                                                                                                                       0.0d0,150.0d0,0.0d0, &
       &                                                                                                                       0.0d0,200.0d0,0.0d0, &
       &                                                                                                                     100.0d0,120.0d0,0.0d0],[3,5])
  character       (len=48                                            ), dimension(5)   :: labelCase
  ! Reference values computed by referenceKing1962.py. Virial radii [Mpc].
  double precision                                                    , parameter      :: radiusVirialHostReference  =2.6384248500d-01   , radiusVirialSatelliteReference=5.6843140239d-02
  ! Tidal radii [Mpc] for each case.
  double precision                                                    , dimension(5)   :: radiusTidalReference       =[1.5084654228d-02,2.0431773802d-02,2.0431773802d-02,4.0112720182d-03,3.0954794353d-02]
  ! Radii [Mpc] where no tidal radius exists, for bound masses of 0.5, 0.9, and 1.1 times the satellite mass. Since the dark
  ! matter fraction is 0.844, these span the cases where the bound dark matter mass is below the dark matter virial mass, where
  ! it is below that mass but the total bound mass is above it, and where it exceeds that mass (and so is capped at the virial
  ! radius).
  double precision                                                    , dimension(3)   :: fractionMassBound          =[0.5d0,0.9d0,1.1d0]
  double precision                                                    , dimension(3)   :: radiusFallbackReference    =[1.8197174841d-02,4.6176596479d-02,5.6843140239d-02]
  character       (len=48                                            ), dimension(3)   :: labelFallback
  ! Tolerances. The tidal radius is found by the root finder in "radiusEnclosingDensityNumerical", which has a relative tolerance
  ! of 10⁻³, so we allow twice that. The radius enclosing a given mass in an NFW profile is found analytically (via the Lambert W
  ! function), as are the virial radii, so the only differences there arise from finding the expansion factor corresponding to
  ! the time of the nodes.
  double precision                                                    , parameter      :: toleranceTidal             =2.0d-3             , toleranceAnalytic        =1.0d-5
  double precision                                                    , dimension(5)   :: radiusTidal
  double precision                                                    , dimension(3)   :: radiusFallback
  integer                                                                              :: i

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Satellite tidal stripping radius: King (1962)")
  ! Labels for each case.
  labelCase    (1)='tangential orbit, R=100 kpc, γ=1'
  labelCase    (2)='radial orbit, R=100 kpc, γ=1'
  labelCase    (3)='tangential orbit, R=100 kpc, γ=0'
  labelCase    (4)='tangential orbit, R=30 kpc, γ=1'
  labelCase    (5)='general orbit, R=164 kpc, γ=½'
  labelFallback(1)='bound dark mass below dark virial mass'
  labelFallback(2)='bound mass between dark and total virial mass'
  labelFallback(3)='bound mass above total virial mass'
  parameters=inputParameters('testSuite/parameters/satellites/tidalStrippingRadiusKing1962.xml')
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Build the objects required by the tidal radius models. These must match those built from the parameter file by the dark matter
  ! profile component.
  allocate(cosmologyParameters_         )
  allocate(cosmologyFunctions_          )
  allocate(virialDensityContrast_       )
  allocate(darkMatterHaloScale_         )
  allocate(satelliteTidalFieldSpherical_)
  allocate(satelliteTidalFieldNull_     )
  allocate(radiusGammaUnity_            )
  allocate(radiusGammaZero_             )
  allocate(radiusGammaHalf_             )
  allocate(radiusNoTide_                )
  cosmologyParameters_         =cosmologyParametersSimple                         (OmegaMatter=0.3153d0,OmegaBaryon=0.0493d0,OmegaDarkEnergy=0.6847d0,temperatureCMB=2.72548d0,HubbleConstant=67.36d0)
  cosmologyFunctions_          =cosmologyFunctionsMatterLambda                    (cosmologyParameters_=cosmologyParameters_)
  virialDensityContrast_       =virialDensityContrastBryanNorman1998              (allowUnsupportedCosmology=.false.,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_)
  darkMatterHaloScale_         =darkMatterHaloScaleVirialDensityContrastDefinition(cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,virialDensityContrast_=virialDensityContrast_)
  satelliteTidalFieldSpherical_=satelliteTidalFieldSphericalSymmetry              (factorBoost=1.0d0,darkMatterHaloScale_=darkMatterHaloScale_)
  satelliteTidalFieldNull_     =satelliteTidalFieldNull                           ()
  radiusGammaUnity_            =satelliteTidalStrippingRadiusKing1962             (efficiencyCentrifugal=1.0d0,applyPreInfall=.false.,cosmologyParameters_=cosmologyParameters_,darkMatterHaloScale_=darkMatterHaloScale_,satelliteTidalField_=satelliteTidalFieldSpherical_)
  radiusGammaZero_             =satelliteTidalStrippingRadiusKing1962             (efficiencyCentrifugal=0.0d0,applyPreInfall=.false.,cosmologyParameters_=cosmologyParameters_,darkMatterHaloScale_=darkMatterHaloScale_,satelliteTidalField_=satelliteTidalFieldSpherical_)
  radiusGammaHalf_             =satelliteTidalStrippingRadiusKing1962             (efficiencyCentrifugal=0.5d0,applyPreInfall=.false.,cosmologyParameters_=cosmologyParameters_,darkMatterHaloScale_=darkMatterHaloScale_,satelliteTidalField_=satelliteTidalFieldSpherical_)
  radiusNoTide_                =satelliteTidalStrippingRadiusKing1962             (efficiencyCentrifugal=0.0d0,applyPreInfall=.false.,cosmologyParameters_=cosmologyParameters_,darkMatterHaloScale_=darkMatterHaloScale_,satelliteTidalField_=satelliteTidalFieldNull_     )
  ! Build a host and a satellite. The satellite is attached to the host as a satellite (not as a progenitor), so that it is
  ! neither on the main branch nor the primary progenitor of its host, and so that the host is the node with which it merges.
  nodeHost              => treeNode(hostTree=tree)
  node                  => treeNode(hostTree=tree)
  tree    %nodeBase       => nodeHost
  nodeHost%firstSatellite => node
  node    %parent         => nodeHost
  basicHost             => nodeHost%basic            (autoCreate=.true.)
  basic                 => node    %basic            (autoCreate=.true.)
  darkMatterProfileHost => nodeHost%darkMatterProfile(autoCreate=.true.)
  darkMatterProfile     => node    %darkMatterProfile(autoCreate=.true.)
  satellite             => node    %satellite        (autoCreate=.true.)
  call basicHost            %massSet     (massHost                             )
  call basicHost            %timeSet     (cosmologyFunctions_%cosmicTime(1.0d0))
  call basic                %massSet     (massSatellite                        )
  call basic                %timeSet     (cosmologyFunctions_%cosmicTime(1.0d0))
  call darkMatterProfileHost%scaleSet    (radiusScaleHost                      )
  call darkMatterProfile    %scaleSet    (radiusScaleSatellite                 )
  call satellite            %boundMassSet(massSatellite                        )
  ! Check that the virial radii match those of the reference calculation. This is not a test of the tidal radius model, but
  ! confirms that the halos here are those assumed by the reference script.
  call Assert('host virial radius matches reference'     ,darkMatterHaloScale_%radiusVirial(nodeHost),radiusVirialHostReference     ,relTol=toleranceAnalytic)
  call Assert('satellite virial radius matches reference',darkMatterHaloScale_%radiusVirial(node    ),radiusVirialSatelliteReference,relTol=toleranceAnalytic)
  ! Evaluate tidal radii for each orbit.
  do i=1,countCases
     call satellite%positionSet(position(:,i))
     call satellite%velocitySet(velocity(:,i))
     select case (i)
     case (3)
        radiusTidal(i)=radiusGammaZero_ %radius(node)
     case (5)
        radiusTidal(i)=radiusGammaHalf_ %radius(node)
     case default
        radiusTidal(i)=radiusGammaUnity_%radius(node)
     end select
     call Assert('tidal radius: '//trim(labelCase(i)),radiusTidal(i),radiusTidalReference(i),relTol=toleranceTidal)
  end do
  ! With no centrifugal term the tidal radius must not depend on the direction of the orbital velocity.
  call Assert('radial orbit and γ=0 agree',radiusTidal(2),radiusTidal(3),relTol=1.0d-6)
  ! Evaluate radii where no tidal radius exists.
  do i=1,size(fractionMassBound)
     call satellite%boundMassSet(fractionMassBound(i)*massSatellite)
     radiusFallback(i)=radiusNoTide_%radius(node)
     call Assert('no tidal radius: '//trim(labelFallback(i)),radiusFallback(i),radiusFallbackReference(i),relTol=toleranceAnalytic)
  end do
  ! A satellite with no bound mass has zero tidal radius.
  call satellite%boundMassSet(0.0d0)
  call Assert('zero bound mass gives zero tidal radius',radiusGammaUnity_%radius(node),0.0d0)
  ! Clean up.
  call node    %destroy()
  call nodeHost%destroy()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
end program Test_Satellite_Tidal_Stripping_Radius_King1962
