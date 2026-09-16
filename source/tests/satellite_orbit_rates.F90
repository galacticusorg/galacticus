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

!!{RST
Contains a program which tests the rates driving satellite orbital evolution against independent reference values.
!!}

!+    Contributions to this file made by: Andrew Benson, Claude.

program Test_Satellite_Orbit_Rates
  !!{RST
  Tests the two rates which drive satellite orbital evolution against reference values computed independently by the
  ``satelliteOrbitRates.py`` script in the `galacticusDevTools
  &lt;https://github.com/galacticusorg/galacticusDevTools&gt;`_ repository: the :cite:t:`chandrasekhar_dynamical_1943` dynamical
  friction acceleration, and the :cite:t:`zentner_physics_2005` tidal mass loss rate together with the
  :cite:t:`king_structure_1962` tidal radius it is built on.

  The rates are compared pointwise, at a grid of phase-space configurations, rather than through an integrated orbit. Decay
  time and bound mass are the *outcome* of these rates, and a comparison of those alone could not distinguish a wrong rate
  from a rate correctly computed but wrongly assembled into the orbital differential equations. Verifying the rates first
  separates the two, as was done for the stellar population yields.

  The grid varies the orbital radius, which sets the host density and velocity dispersion the friction term sees and the
  tidal field the stripping term sees; the direction of the velocity, which separates the angular and radial orbital
  frequencies that :cite:t:`zentner_physics_2005` takes the larger of, and moves the Chandrasekhar velocity ratio across the
  range where its error function factor varies; and the satellite's mass and concentration. One configuration places the
  satellite deep inside its own virial radius, where the extended-mass suppression factor in the Chandrasekhar integral is
  genuinely below unity - about 0.37 - rather than saturated, which nothing else in the grid exercises. Note that at fixed
  concentration both rates are exactly linear in the satellite mass, so the configurations varying it alone check that
  proportionality and nothing more.

  The satellite's bound mass is always the total mass of its own profile. Pairing a reduced bound mass with an unstripped
  profile is a state the code never produces, and it would drive the mass outside the tidal radius to zero, making the
  comparison one of zero against zero.
  !!}
  use :: Display                       , only : displayMessage                   , displayVerbositySet               , verbosityLevelStandard
  use :: Error                         , only : Error_Handler_Register
  use :: Events_Hooks                  , only : eventsHooksInitialize
  use :: Functions_Global_Utilities    , only : Functions_Global_Set
  use :: Dark_Matter_Halo_Scales       , only : darkMatterHaloScaleClass
  use :: Galacticus_Nodes              , only : nodeClassHierarchyFinalize       , nodeClassHierarchyInitialize      , nodeComponentBasic                 , nodeComponentDarkMatterProfile, &
       &                                        nodeComponentSatellite           , treeNode
  use :: Input_Parameters              , only : inputParameters
  use :: ISO_Varying_String            , only : varying_string                   , assignment(=)                     , var_str
  use :: Node_Components               , only : Node_Components_Initialize       , Node_Components_Thread_Initialize , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Satellite_Dynamical_Friction  , only : satelliteDynamicalFrictionClass
  use :: Satellite_Tidal_Stripping     , only : satelliteTidalStrippingClass
  use :: Satellite_Tidal_Stripping_Radii, only : satelliteTidalStrippingRadiusClass
  use :: Unit_Tests                    , only : Assert                           , Unit_Tests_Begin_Group            , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  ! The host halo, shared by every configuration.
  double precision, parameter :: massHost               =1.0d12, concentrationHost=10.0d0
  double precision, parameter :: timeNode               =13.8d0
  ! The reference configurations and rates, emitted by `satelliteOrbitRates.py --fortran`.
  integer         , parameter :: countConfigurations=14
  double precision, dimension(countConfigurations), parameter :: radius                 =[ 6.218943119192954d-02, 6.218943119192954d-02, 6.218943119192954d-02, 1.554735779798238d-01, 1.554735779798238d-01, 1.554735779798238d-01, 3.109471559596477d-01, 3.109471559596477d-01, 3.109471559596477d-01, 1.554735779798238d-01, 1.554735779798238d-01, 1.554735779798238d-01, 1.554735779798238d-01, 3.109471559596477d-02]
  double precision, dimension(countConfigurations), parameter :: velocityRadial         =[ 0.000000000000000d+00,-5.880586697325706d+01,-1.058505605518627d+02, 0.000000000000000d+00,-5.880586697325706d+01,-1.058505605518627d+02, 0.000000000000000d+00,-5.880586697325706d+01,-1.058505605518627d+02,-5.880586697325706d+01,-5.880586697325706d+01,-5.880586697325706d+01,-5.880586697325706d+01,-5.880586697325706d+01]
  double precision, dimension(countConfigurations), parameter :: velocityTangential     =[ 1.176117339465141d+02, 5.880586697325706d+01, 2.352234678930283d+01, 1.176117339465141d+02, 5.880586697325706d+01, 2.352234678930283d+01, 1.176117339465141d+02, 5.880586697325706d+01, 2.352234678930283d+01, 5.880586697325706d+01, 5.880586697325706d+01, 5.880586697325706d+01, 5.880586697325706d+01, 5.880586697325706d+01]
  double precision, dimension(countConfigurations), parameter :: massSatellite          =[ 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+09, 1.000000000000000d+11, 1.000000000000000d+11]
  double precision, dimension(countConfigurations), parameter :: concentrationSatellite =[ 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 5.000000000000000d+00, 3.000000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01]
  double precision, dimension(countConfigurations), parameter :: accelerationXReference =[-0.000000000000000d+00, 8.227487739763317d+00, 1.210102576731724d+01,-0.000000000000000d+00, 1.214489492257159d+00, 1.677253745658934d+00,-0.000000000000000d+00, 2.562005033353899d-01, 3.246185578498261d-01, 1.214489492257159d+00, 1.214489492257159d+00, 1.214489492257159d-01, 1.214489492257159d+01, 1.249521497639164d+02]
  double precision, dimension(countConfigurations), parameter :: accelerationYReference =[-1.236727818950592d+01,-8.227487739763317d+00,-2.689116837181610d+00,-1.672708177997145d+00,-1.214489492257159d+00,-3.727230545908742d-01,-3.136121117902790d-01,-2.562005033353899d-01,-7.213745729996136d-02,-1.214489492257159d+00,-1.214489492257159d+00,-1.214489492257159d-01,-1.214489492257159d+01,-1.249521497639164d+02]
  double precision, dimension(countConfigurations), parameter :: radiusTidalReference   =[ 1.104123197924320d-02, 1.350429546012691d-02, 1.454231114437747d-02, 2.458983649053660d-02, 2.892694227370823d-02, 3.066739554481213d-02, 4.540098877351536d-02, 5.367071859024051d-02, 5.710513523656589d-02, 2.551288381211287d-02, 3.014988927732068d-02, 1.342669722482055d-02, 6.232120791102678d-02, 1.737763807376192d-02]
  double precision, dimension(countConfigurations), parameter :: rateMassLossReference  =[-3.518323768396382d+10,-2.615997176117954d+10,-2.331310955015278d+10,-9.123778069701977d+09,-6.441304614718612d+09,-5.628203204272208d+09,-1.810087941830191d+09,-8.498977712976755d+08,-5.691448590360527d+08,-9.800858222280577d+09,-5.025796852263382d+09,-6.441304614718614d+08,-6.441304614718614d+10,-5.343634066009216d+11]
  ! The assertion tolerances, justified rather than tuned, and different for the two rates because their floors differ.
  !
  ! The tidal radius is found by the root finder in `radiusEnclosingDensityNumerical`, whose relative tolerance is 10^-3, so
  ! nothing downstream of it can agree better than that; 2 x 10^-3 is allowed, following the same reasoning as the King (1962)
  ! test in PR #1499. Measured here: 2.1 x 10^-4.
  !
  ! The mass loss rate is built on that radius and inherits its error *amplified*, because the mass outside the tidal radius,
  ! `boundMass - M(<r_tidal)`, is a difference of two comparable masses. The amplification is largest where the satellite is
  ! least stripped - at an orbital radius of one host virial radius the tidal radius approaches the satellite's own virial
  ! radius and the difference becomes small - and reaches a factor of about 8.5 across this grid, turning a radius difference of
  ! 1.9 x 10^-4 into 1.6 x 10^-3. With the solver's own 10^-3 that bounds the rate at roughly 8.5 x 10^-3; 5 x 10^-3 is allowed,
  ! which keeps a factor of three over the largest value measured.
  !
  ! Note that this sensitivity is a consequence of correcting the mass normalization: while `massOuter` carried a spurious
  ! floor of the baryon fraction of the bound mass it was never small, and so was never delicate.
  !
  ! The dynamical friction acceleration does not go through that root finder. Its floor is the quadrature of the isotropic
  ! Jeans integral giving the host velocity dispersion, and it is held to a tighter tolerance so that the assertion keeps its
  ! teeth. Measured here: 6.7 x 10^-5.
  double precision, parameter :: toleranceAcceleration  =5.0d-4, toleranceTidal=2.0d-3, toleranceRateMassLoss=5.0d-3
  class           (darkMatterHaloScaleClass          ), pointer :: darkMatterHaloScale_
  class           (satelliteDynamicalFrictionClass   ), pointer :: satelliteDynamicalFriction_
  class           (satelliteTidalStrippingClass      ), pointer :: satelliteTidalStripping_
  class           (satelliteTidalStrippingRadiusClass), pointer :: satelliteTidalStrippingRadius_
  type            (treeNode                          ), pointer :: nodeHost                          , nodeSatellite
  class           (nodeComponentBasic                ), pointer :: basicHost                         , basicSatellite
  class           (nodeComponentDarkMatterProfile    ), pointer :: profileHost                       , profileSatellite
  class           (nodeComponentSatellite            ), pointer :: satellite
  type            (inputParameters                   )          :: parameters
  character       (len=128                           )          :: message
  integer                                                       :: iConfiguration
  double precision                                              :: radiusVirialHost                  , radiusVirialSatellite      , &
       &                                                           differenceMaximum                 , differenceAccelerationX    , &
       &                                                           differenceAccelerationY           , differenceRadiusTidal      , &
       &                                                           differenceRateMassLoss            , radiusTidal                , &
       &                                                           rateMassLoss
  double precision                                , dimension(3) :: acceleration

  call displayVerbositySet              (verbosityLevelStandard)
  call Error_Handler_Register           (                      )
  parameters=inputParameters(var_str('testSuite/parameters/satelliteOrbitRates.xml'))
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <objectBuilder class="darkMatterHaloScale"           name="darkMatterHaloScale_"           source="parameters"/>
  <objectBuilder class="satelliteDynamicalFriction"    name="satelliteDynamicalFriction_"    source="parameters"/>
  <objectBuilder class="satelliteTidalStripping"       name="satelliteTidalStripping_"       source="parameters"/>
  <objectBuilder class="satelliteTidalStrippingRadius" name="satelliteTidalStrippingRadius_" source="parameters"/>
  !!]
  call Unit_Tests_Begin_Group("Satellite orbital evolution rates")
  differenceMaximum=0.0d0
  do iConfiguration=1,countConfigurations
     ! Build a fresh host and satellite for every configuration. A node memoizes the mass distributions built from its
     ! components, so reusing one and resetting its properties would leave the previous configuration's distribution in place.
     nodeHost         => treeNode                  (                 )
     nodeSatellite    => treeNode                  (                 )
     basicHost        => nodeHost     %basic            (autoCreate=.true.)
     profileHost      => nodeHost     %darkMatterProfile(autoCreate=.true.)
     basicSatellite   => nodeSatellite%basic            (autoCreate=.true.)
     profileSatellite => nodeSatellite%darkMatterProfile(autoCreate=.true.)
     satellite        => nodeSatellite%satellite        (autoCreate=.true.)
     ! Attach the satellite to its host, so that `mergesWith` and `isSatellite` resolve.
     nodeHost     %firstSatellite => nodeSatellite
     nodeSatellite%parent         => nodeHost
     call basicHost     %massSet(massHost                        )
     call basicHost     %timeSet(timeNode                        )
     call basicSatellite%massSet(massSatellite (iConfiguration)  )
     call basicSatellite%timeSet(timeNode                        )
     radiusVirialHost     =darkMatterHaloScale_%radiusVirial(nodeHost     )
     radiusVirialSatellite=darkMatterHaloScale_%radiusVirial(nodeSatellite)
     call profileHost     %scaleSet(radiusVirialHost     /concentrationHost                    )
     call profileSatellite%scaleSet(radiusVirialSatellite/concentrationSatellite(iConfiguration))
     ! Place the satellite. The orbit lies in the x-y plane, with the radial direction along x.
     call satellite%boundMassSet(massSatellite(iConfiguration))
     call satellite%positionSet ([radius        (iConfiguration),0.0d0                            ,0.0d0])
     call satellite%velocitySet ([velocityRadial(iConfiguration),velocityTangential(iConfiguration),0.0d0])
     ! Evaluate the rates.
     acceleration=satelliteDynamicalFriction_   %acceleration (nodeSatellite)
     radiusTidal =satelliteTidalStrippingRadius_%radius       (nodeSatellite)
     rateMassLoss=satelliteTidalStripping_      %massLossRate (nodeSatellite)
     ! Report the fractional difference of every comparison, so that a reader can see how much of the tolerance is in use. A
     ! test passing only because its tolerance is loose otherwise looks identical to one passing sharply.
     differenceAccelerationX=abs(acceleration(1)/accelerationXReference(iConfiguration)-1.0d0)
     differenceAccelerationY=abs(acceleration(2)/accelerationYReference(iConfiguration)-1.0d0)
     differenceRadiusTidal  =abs(radiusTidal    /radiusTidalReference  (iConfiguration)-1.0d0)
     differenceRateMassLoss =abs(rateMassLoss   /rateMassLossReference (iConfiguration)-1.0d0)
     ! The x-acceleration is exactly zero for circular orbits, where the fractional difference is not defined; skip it there.
     if (accelerationXReference(iConfiguration) /= 0.0d0) &
          & differenceMaximum=max(differenceMaximum,differenceAccelerationX)
     differenceMaximum=max(differenceMaximum,differenceAccelerationY,differenceRadiusTidal,differenceRateMassLoss)
     write (message,'(a,i0,a,e12.5,a,e12.5,a,e12.5)') 'configuration ',iConfiguration,': fractional difference, a ', &
          & differenceAccelerationY,', r_tidal ',differenceRadiusTidal,', dM/dt ',differenceRateMassLoss
     call displayMessage(trim(message))
     write (message,'(a,i0)') 'dynamical friction acceleration, radial, configuration ',iConfiguration
     call Assert(trim(message),acceleration(1),accelerationXReference(iConfiguration),absTol=1.0d-12,relTol=toleranceAcceleration)
     write (message,'(a,i0)') 'dynamical friction acceleration, tangential, configuration ',iConfiguration
     call Assert(trim(message),acceleration(2),accelerationYReference(iConfiguration),relTol=toleranceAcceleration)
     write (message,'(a,i0)') 'tidal radius, configuration ',iConfiguration
     call Assert(trim(message),radiusTidal    ,radiusTidalReference  (iConfiguration),relTol=toleranceTidal)
     write (message,'(a,i0)') 'tidal mass loss rate, configuration ',iConfiguration
     call Assert(trim(message),rateMassLoss   ,rateMassLossReference (iConfiguration),relTol=toleranceRateMassLoss)
     ! Detach before destroying, so that the host does not attempt to destroy the satellite a second time.
     nodeHost     %firstSatellite => null()
     nodeSatellite%parent         => null()
     call nodeSatellite%destroy()
     call nodeHost     %destroy()
     deallocate(nodeSatellite)
     deallocate(nodeHost     )
  end do
  write (message,'(a,e12.5)') 'largest fractional difference over all configurations: ',differenceMaximum
  call displayMessage(trim(message))
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  !![
  <objectDestructor name="darkMatterHaloScale_"          />
  <objectDestructor name="satelliteDynamicalFriction_"   />
  <objectDestructor name="satelliteTidalStripping_"      />
  <objectDestructor name="satelliteTidalStrippingRadius_"/>
  !!]
end program Test_Satellite_Orbit_Rates
