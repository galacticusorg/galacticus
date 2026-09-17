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
Contains a program which tests the satellite tidal heating rate against independent reference values.
!!}

!+    Contributions to this file made by: Andrew Benson, Claude.

program Test_Satellite_Tidal_Heating_Rate
  !!{RST
  Tests the :cite:t:`gnedin_tidal_1999` satellite tidal heating rate against reference values computed independently by the
  ``satelliteTidalHeatingRate.py`` script in the `galacticusDevTools
  &lt;https://github.com/galacticusorg/galacticusDevTools&gt;`_ repository.

  The rate is compared pointwise, with the path-integrated tidal tensor set directly, rather than through an integrated orbit.
  The accumulation of that path integral along an orbit is a separate step, and a comparison of the accumulated heating alone
  could not distinguish a wrong rate from a correct rate wrongly accumulated.

  The satellite is placed off the coordinate axes, and every path integral has non-zero off-diagonal elements, so that both the
  orientation of the host's tidal tensor and the double contraction of the two tensors are exercised; on an axis the tidal
  tensor is diagonal and would test neither. One configuration on the x-axis is kept as the simplest case. The radius, speed,
  satellite concentration and bound mass move the adiabatic correction :math:`[1+(\omega\tau)^2]^{-\gamma}` from near unity to
  :math:`\sim 10^{-5}`, so that its exponent is tested and not merely its presence. Bound masses below the basic mass exercise
  the minimum taken in the satellite's half mass, a tiny bound mass exercises the fallback to the virial orbital frequency, and a
  path integral anti-aligned with the tidal field exercises the clamping of a negative rate to zero.

  The satellite's basic mass is not varied: with a fixed virial density contrast the internal orbital frequency, and so the rate,
  is independent of it, and such configurations would repeat another exactly.
  !!}
  use :: Display                   , only : displayMessage                , displayVerbositySet              , verbosityLevelStandard
  use :: Error                     , only : Error_Handler_Register
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass
  use :: Galacticus_Nodes          , only : nodeClassHierarchyFinalize    , nodeClassHierarchyInitialize     , nodeComponentBasic                 , nodeComponentDarkMatterProfile, &
       &                                    nodeComponentSatellite        , treeNode
  use :: Input_Parameters          , only : inputParameters
  use :: ISO_Varying_String        , only : var_str
  use :: Node_Components           , only : Node_Components_Initialize    , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Satellite_Tidal_Heating   , only : satelliteTidalHeatingRateClass
  use :: Tensors                   , only : tensorRank2Dimension3Symmetric
  use :: Unit_Tests                , only : Assert                        , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  ! The host halo, shared by every configuration.
  double precision, parameter :: massHost           =1.0d12, concentrationHost=10.0d0
  double precision, parameter :: timeNode           =13.8d0
  ! The reference configurations and rates, emitted by `satelliteTidalHeatingRate.py --fortran`. Positions are in Mpc, speeds in
  ! km/s, path-integrated tidal tensors in (km/s/Mpc)² Gyr, and heating rates in (km/s/Mpc)² Gyr⁻¹.
  integer         , parameter :: countConfigurations=13
  double precision, dimension(countConfigurations), parameter :: positionX              =[ 1.795194348578051d-02, 5.385583045734153d-02, 1.795194348578051d-01, 2.826842234878456d-02, 2.826842234878456d-02, 9.328103463593093d-02, 2.826842234878456d-02, 2.826842234878456d-02, 2.826842234878456d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02]
  double precision, dimension(countConfigurations), parameter :: positionY              =[ 1.795194348578051d-02, 5.385583045734153d-02, 1.795194348578051d-01,-4.711403724797426d-02,-4.711403724797426d-02, 0.000000000000000d+00,-4.711403724797426d-02,-4.711403724797426d-02,-4.711403724797426d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02]
  double precision, dimension(countConfigurations), parameter :: positionZ              =[ 1.795194348578051d-02, 5.385583045734153d-02, 1.795194348578051d-01, 7.538245959675884d-02, 7.538245959675884d-02, 0.000000000000000d+00, 7.538245959675884d-02, 7.538245959675884d-02, 7.538245959675884d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02, 5.385583045734153d-02]
  double precision, dimension(countConfigurations), parameter :: speed                  =[ 2.352273917617680d+02, 1.764205438213260d+02, 5.880684794044201d+01, 7.056821752853041d+01, 4.704547835235361d+02, 1.764205438213260d+02, 1.764205438213260d+02, 1.764205438213260d+02, 1.764205438213260d+02, 1.764205438213260d+02, 1.764205438213260d+02, 1.764205438213260d+02, 1.764205438213260d+02]
  double precision, dimension(countConfigurations), parameter :: tensorPathIntegratedXX =[ 1.430775831020361d+05, 1.430775831020361d+05, 1.430775831020361d+05,-4.292327493061082d+04,-4.292327493061082d+04, 1.430775831020361d+05,-1.430775831020361d+05,-1.430775831020361d+05,-4.292327493061082d+04, 4.292327493061082d+04, 4.292327493061082d+04, 4.292327493061082d+04,-1.430775831020361d+05]
  double precision, dimension(countConfigurations), parameter :: tensorPathIntegratedXY =[ 5.723103324081444d+04, 5.723103324081444d+04, 5.723103324081444d+04, 1.144620664816289d+05, 1.144620664816289d+05, 5.723103324081444d+04,-5.723103324081444d+04,-5.723103324081444d+04, 1.144620664816289d+05,-1.144620664816289d+05,-1.144620664816289d+05,-1.144620664816289d+05,-5.723103324081444d+04]
  double precision, dimension(countConfigurations), parameter :: tensorPathIntegratedXZ =[-4.292327493061082d+04,-4.292327493061082d+04,-4.292327493061082d+04,-7.153879155101805d+04,-7.153879155101805d+04,-4.292327493061082d+04, 4.292327493061082d+04, 4.292327493061082d+04,-7.153879155101805d+04, 7.153879155101805d+04, 7.153879155101805d+04, 7.153879155101805d+04, 4.292327493061082d+04]
  double precision, dimension(countConfigurations), parameter :: tensorPathIntegratedYY =[-7.153879155101805d+04,-7.153879155101805d+04,-7.153879155101805d+04,-8.584654986122165d+04,-8.584654986122165d+04,-7.153879155101805d+04, 7.153879155101805d+04, 7.153879155101805d+04,-8.584654986122165d+04, 8.584654986122165d+04, 8.584654986122165d+04, 8.584654986122165d+04, 7.153879155101805d+04]
  double precision, dimension(countConfigurations), parameter :: tensorPathIntegratedYZ =[ 2.861551662040722d+04, 2.861551662040722d+04, 2.861551662040722d+04,-1.001543081714253d+05,-1.001543081714253d+05, 2.861551662040722d+04,-2.861551662040722d+04,-2.861551662040722d+04,-1.001543081714253d+05, 1.001543081714253d+05, 1.001543081714253d+05, 1.001543081714253d+05,-2.861551662040722d+04]
  double precision, dimension(countConfigurations), parameter :: tensorPathIntegratedZZ =[-7.153879155101805d+04,-7.153879155101805d+04,-7.153879155101805d+04, 1.287698247918325d+05, 1.287698247918325d+05,-7.153879155101805d+04, 7.153879155101805d+04, 7.153879155101805d+04, 1.287698247918325d+05,-1.287698247918325d+05,-1.287698247918325d+05,-1.287698247918325d+05, 7.153879155101805d+04]
  double precision, dimension(countConfigurations), parameter :: massSatellite          =[ 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10]
  double precision, dimension(countConfigurations), parameter :: massBound              =[ 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 1.000000000000000d+10, 6.000000000000000d+09, 3.000000000000000d+09, 5.000000000000000d+08, 1.000000000000000d+03, 1.000000000000000d+10]
  double precision, dimension(countConfigurations), parameter :: concentrationSatellite =[ 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 5.000000000000000d+00, 3.000000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01, 1.500000000000000d+01]
  double precision, dimension(countConfigurations), parameter :: rateHeatingReference   =[ 7.352871526342838d+05, 4.081543134638479d+04, 3.857328174759763d-01, 7.268422955010923d+03, 2.459168536178509d+05, 2.040771567319241d+05, 3.089602274814853d+05, 6.794783392315524d+04, 1.838192993311568d+04, 8.259565157905137d+02, 7.302336898509945d+00, 1.462768047066803d+05, 0.000000000000000d+00]
  ! The assertion tolerance, justified rather than tuned.
  !
  ! Every ingredient of the rate that depends only on the geometry - the tidal tensor's orientation and contraction, and the
  ! half-mass radius as a fraction of the virial radius - was measured to agree with the reference to rounding. What remains
  ! is a difference of 7.5 x 10⁻⁶ in the virial radii themselves (a mean density differing by 2.3 x 10⁻⁵), not traced further.
  ! The adiabatic correction amplifies a fractional error in the internal orbital frequency by up to 2γ x²/(1+x²), approaching
  ! 5 at the largest x = ωτ in the grid, so that the rate inherits up to about 6 x 10⁻⁵. Measured: 4.7 x 10⁻⁵. 2 x 10⁻⁴ is
  ! allowed. Any error in the form of the rate - a missing factor of 2π in x, off-diagonal elements counted once, a wrong
  ! exponent - changes it by far more.
  double precision                                , parameter :: toleranceRate           =2.0d-4
  class           (darkMatterHaloScaleClass      ), pointer   :: darkMatterHaloScale_
  class           (satelliteTidalHeatingRateClass), pointer   :: satelliteTidalHeatingRate_
  type            (treeNode                      ), pointer   :: nodeHost                  , nodeSatellite
  class           (nodeComponentBasic            ), pointer   :: basicHost                 , basicSatellite
  class           (nodeComponentDarkMatterProfile), pointer   :: profileHost               , profileSatellite
  class           (nodeComponentSatellite        ), pointer   :: satellite
  type            (inputParameters               )            :: parameters
  character       (len=128                       )            :: message
  integer                                                     :: iConfiguration
  double precision                                            :: radiusVirialHost          , radiusVirialSatellite, &
       &                                                         differenceMaximum         , differenceRate       , &
       &                                                         rateHeating

  call displayVerbositySet              (verbosityLevelStandard)
  call Error_Handler_Register           (                      )
  parameters=inputParameters(var_str('testSuite/parameters/satelliteTidalHeatingRate.xml'))
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <objectBuilder class="darkMatterHaloScale"       name="darkMatterHaloScale_"       source="parameters"/>
  <objectBuilder class="satelliteTidalHeatingRate" name="satelliteTidalHeatingRate_" source="parameters"/>
  !!]
  call Unit_Tests_Begin_Group("Satellite tidal heating rate")
  differenceMaximum=0.0d0
  do iConfiguration=1,countConfigurations
     ! Build a fresh host and satellite for every configuration. A node memoizes the mass distributions built from its
     ! components, so reusing one and resetting its properties would leave the previous configuration's distribution in place.
     nodeHost         => treeNode                       (                 )
     nodeSatellite    => treeNode                       (                 )
     basicHost        => nodeHost     %basic            (autoCreate=.true.)
     profileHost      => nodeHost     %darkMatterProfile(autoCreate=.true.)
     basicSatellite   => nodeSatellite%basic            (autoCreate=.true.)
     profileSatellite => nodeSatellite%darkMatterProfile(autoCreate=.true.)
     satellite        => nodeSatellite%satellite        (autoCreate=.true.)
     ! Attach the satellite to its host, so that `mergesWith` and `isSatellite` resolve.
     nodeHost     %firstSatellite => nodeSatellite
     nodeSatellite%parent         => nodeHost
     call basicHost     %massSet(massHost                     )
     call basicHost     %timeSet(timeNode                     )
     call basicSatellite%massSet(massSatellite(iConfiguration))
     call basicSatellite%timeSet(timeNode                     )
     radiusVirialHost     =darkMatterHaloScale_%radiusVirial(nodeHost     )
     radiusVirialSatellite=darkMatterHaloScale_%radiusVirial(nodeSatellite)
     call profileHost     %scaleSet(radiusVirialHost     /concentrationHost                    )
     call profileSatellite%scaleSet(radiusVirialSatellite/concentrationSatellite(iConfiguration))
     ! Place the satellite. Only the speed enters the rate, so the velocity is set along the y-axis throughout.
     call satellite%boundMassSet                (massBound(iConfiguration))
     call satellite%positionSet                 ([positionX(iConfiguration),positionY(iConfiguration),positionZ(iConfiguration)])
     call satellite%velocitySet                 ([0.0d0                    ,speed    (iConfiguration),0.0d0                    ])
     call satellite%tidalTensorPathIntegratedSet(                                           &
          &                                      tensorRank2Dimension3Symmetric(            &
          &                                                                     tensorPathIntegratedXX(iConfiguration), &
          &                                                                     tensorPathIntegratedXY(iConfiguration), &
          &                                                                     tensorPathIntegratedXZ(iConfiguration), &
          &                                                                     tensorPathIntegratedYY(iConfiguration), &
          &                                                                     tensorPathIntegratedYZ(iConfiguration), &
          &                                                                     tensorPathIntegratedZZ(iConfiguration)  &
          &                                                                    )            &
          &                                     )
     ! Evaluate the rate.
     rateHeating=satelliteTidalHeatingRate_%heatingRate(nodeSatellite)
     ! Report the fractional difference of every comparison, so that a reader can see how much of the tolerance is in use. A
     ! test passing only because its tolerance is loose otherwise looks identical to one passing sharply. The rate clamped to zero
     ! has no fractional difference, and is compared absolutely below.
     if (rateHeatingReference(iConfiguration) > 0.0d0) then
        differenceRate   =abs(rateHeating/rateHeatingReference(iConfiguration)-1.0d0)
        differenceMaximum=max(differenceMaximum,differenceRate)
        write (message,'(a,i0,a,e12.5,a,e12.5)') 'configuration ',iConfiguration,': dQ/dt ',rateHeating,', fractional difference ',differenceRate
     else
        write (message,'(a,i0,a,e12.5)'        ) 'configuration ',iConfiguration,': dQ/dt ',rateHeating
     end if
     call displayMessage(trim(message))
     write (message,'(a,i0)') 'tidal heating rate, configuration ',iConfiguration
     call Assert(trim(message),rateHeating,rateHeatingReference(iConfiguration),absTol=1.0d-30,relTol=toleranceRate)
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
  <objectDestructor name="darkMatterHaloScale_"      />
  <objectDestructor name="satelliteTidalHeatingRate_"/>
  !!]
end program Test_Satellite_Tidal_Heating_Rate
