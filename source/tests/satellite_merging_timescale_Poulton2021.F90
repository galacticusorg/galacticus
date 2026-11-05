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
Contains a program which tests the :cite:t:`poulton_extracting_2021` satellite merging timescale.
!!}

!+    Contributions to this file made by: Claude.

program Test_Satellite_Merging_Timescale_Poulton2021
  !!{RST
  Tests the :cite:t:`poulton_extracting_2021` satellite merging timescale against a direct evaluation of their equations 2--8
  for an NFW host. The reference evaluation takes the host mass enclosed within the orbital radius from the analytic NFW
  enclosed mass, capped at the virial mass beyond the virial radius (their equation 2), and computes from it, and the reduced
  mass of the satellite, the orbital energy, angular momentum, eccentricity, pericentric distance, and dynamical time. Orbits
  both inside and outside of the virial radius are tested, so that the cap is exercised, as are unbound orbits, which receive
  the same fitted, finite, merging time.
  !!}
  use :: Display                         , only : displayVerbositySet                   , verbosityLevelStandard
  use :: Error                           , only : Error_Handler_Register
  use :: Events_Hooks                    , only : eventsHooksInitialize
  use :: Functions_Global_Utilities      , only : Functions_Global_Set
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMOClass
  use :: Galacticus_Nodes                , only : nodeClassHierarchyFinalize            , nodeClassHierarchyInitialize     , nodeComponentBasic                 , nodeComponentDarkMatterProfile, &
       &                                          treeNode
  use :: Input_Parameters                , only : inputParameters
  use :: ISO_Varying_String              , only : var_str
  use :: Kepler_Orbits                   , only : keplerOrbit
  use :: Node_Components                 , only : Node_Components_Initialize            , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal        , MpcPerKmPerSToGyr
  use :: Satellite_Merging_Timescales    , only : satelliteMergingTimescalesPoulton2021
  use :: Unit_Tests                      , only : Assert                                , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish, &
       &                                          compareGreaterThan
  implicit none
  ! The host halo, and the parameters of the Poulton et al. (2021) fit (their favored values).
  double precision                                       , parameter    :: massHost                    = 1.0d12                  , concentrationHost  =10.0d0, &
       &                                                                   timeNode                    =13.8d00
  double precision                                       , parameter    :: A                           = 5.5d00                  , c                  = 0.2d0, &
       &                                                                   bInterior                   =-0.5d00                  , bExterior          =-1.0d0
  ! Orbital configurations: radius (in units of the virial radius), radial and tangential velocities (in units of the virial
  ! velocity), and satellite mass. The third set of velocities gives an unbound orbit at every radius.
  integer                                                , parameter    :: countRadii                  =4                        , countVelocities    =3     , &
       &                                                                   countMasses                 =2
  double precision, dimension(countRadii     )           , parameter    :: radiusFractional            =[0.3d0,0.7d0,1.5d0,2.0d0]
  double precision, dimension(countVelocities)           , parameter    :: velocityRadialFractional    =[-0.5d0,-0.9d0,-1.6d0]   , &
       &                                                                   velocityTangentialFractional=[+0.6d0,+0.3d0,+1.2d0]
  double precision, dimension(countMasses    )           , parameter    :: massSatellite               =[1.0d10,1.0d11]
  class           (darkMatterHaloScaleClass             ), pointer      :: darkMatterHaloScale_
  class           (darkMatterProfileDMOClass            ), pointer      :: darkMatterProfileDMO_
  type            (satelliteMergingTimescalesPoulton2021)               :: satelliteMergingTimescales_
  type            (treeNode                             ), pointer      :: nodeHost                                              , nodeSatellite
  class           (nodeComponentBasic                   ), pointer      :: basicHost                                             , basicSatellite
  class           (nodeComponentDarkMatterProfile       ), pointer      :: profileHost
  type            (inputParameters                      )               :: parameters
  type            (keplerOrbit                          )               :: orbit
  character       (len=128                              )               :: message
  integer                                                               :: iRadius                                               , iVelocity                 , &
       &                                                                   iMass
  double precision                                                      :: radiusVirial                                          , velocityVirial            , &
       &                                                                   radius                                                , velocityRadial            , &
       &                                                                   velocityTangential                                    , massEnclosed              , &
       &                                                                   massReduced                                           , energyOrbital             , &
       &                                                                   angularMomentum                                       , eccentricity              , &
       &                                                                   radiusPericenter                                      , timeDynamical             , &
       &                                                                   timeMergingReference                                  , timeMerging               , &
       &                                                                   radiusScaleFree

  call displayVerbositySet              (verbosityLevelStandard)
  call Error_Handler_Register           (                      )
  parameters=inputParameters(var_str('testSuite/parameters/satelliteOrbitRates.xml'))
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
  <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
  !!]
  satelliteMergingTimescales_=satelliteMergingTimescalesPoulton2021(A,c,bInterior,bExterior,darkMatterHaloScale_,darkMatterProfileDMO_)
  call Unit_Tests_Begin_Group("Poulton et al. (2021) satellite merging timescale")
  do iMass=1,countMasses
     ! Build a fresh host and satellite, with the satellite attached to its host.
     nodeHost       => treeNode                       (                 )
     nodeSatellite  => treeNode                       (                 )
     basicHost      => nodeHost     %basic            (autoCreate=.true.)
     profileHost    => nodeHost     %darkMatterProfile(autoCreate=.true.)
     basicSatellite => nodeSatellite%basic            (autoCreate=.true.)
     nodeHost     %firstSatellite => nodeSatellite
     nodeSatellite%parent         => nodeHost
     call basicHost     %massSet(massHost            )
     call basicHost     %timeSet(timeNode            )
     call basicSatellite%massSet(massSatellite(iMass))
     call basicSatellite%timeSet(timeNode            )
     radiusVirial  =darkMatterHaloScale_%radiusVirial  (nodeHost)
     velocityVirial=darkMatterHaloScale_%velocityVirial(nodeHost)
     call profileHost%scaleSet(radiusVirial/concentrationHost)
     do iRadius=1,countRadii
        do iVelocity=1,countVelocities
           radius            =radiusFractional            (iRadius  )*radiusVirial
           velocityRadial    =velocityRadialFractional    (iVelocity)*velocityVirial
           velocityTangential=velocityTangentialFractional(iVelocity)*velocityVirial
           orbit=keplerOrbit()
           call orbit%massesSet            (massSatellite(iMass),massHost)
           call orbit%radiusSet            (radius                       )
           call orbit%velocityRadialSet    (velocityRadial               )
           call orbit%velocityTangentialSet(velocityTangential           )
           timeMerging=satelliteMergingTimescales_%timeUntilMerging(nodeSatellite,orbit)
           ! Reference evaluation of Poulton et al. (2021) equations 2-8.
           radiusScaleFree =concentrationHost*min(radius,radiusVirial)/radiusVirial
           massEnclosed    =+massHost                                                                   &
                &           *(log(1.0d0+radiusScaleFree  )-radiusScaleFree  /(1.0d0+radiusScaleFree  )) &
                &           /(log(1.0d0+concentrationHost)-concentrationHost/(1.0d0+concentrationHost))
           massReduced     =+massSatellite(iMass)*massEnclosed/(massSatellite(iMass)+massEnclosed)
           energyOrbital   =+0.5d0*massReduced*(velocityRadial**2+velocityTangential**2) &
                &           -gravitationalConstant_internal*massEnclosed*massSatellite(iMass)/radius
           angularMomentum =+massReduced*radius*velocityTangential
           eccentricity    =+sqrt(1.0d0+2.0d0*energyOrbital*angularMomentum**2/(gravitationalConstant_internal*massEnclosed*massSatellite(iMass))**2/massReduced)
           radiusPericenter=+angularMomentum**2/(1.0d0+eccentricity)/gravitationalConstant_internal/massEnclosed/massSatellite(iMass)/massReduced
           timeDynamical   =+sqrt(radius**3/gravitationalConstant_internal/massEnclosed)*MpcPerKmPerSToGyr
           if (radius < radiusVirial) then
              timeMergingReference=A*timeDynamical*(radius/radiusVirial)**bInterior*(radiusPericenter/radiusVirial)**c
           else
              timeMergingReference=A*timeDynamical*(radius/radiusVirial)**bExterior*(radiusPericenter/radiusVirial)**c
           end if
           write (message,'(a,f3.1,a,i1,a,e7.1)') 'r/R_vir=',radiusFractional(iRadius),', velocity ',iVelocity,', M_sat=',massSatellite(iMass)
           call Assert(trim(message),timeMerging,timeMergingReference,relTol=1.0d-6)
           ! Confirm that the third set of velocities does give an unbound orbit, so that the unbound case is tested.
           if (iVelocity == 3) call Assert(trim(message)//' (unbound)',eccentricity,1.0d0,compare=compareGreaterThan)
        end do
     end do
     ! Detach before destroying, so that the host does not attempt to destroy the satellite a second time.
     nodeHost     %firstSatellite => null()
     nodeSatellite%parent         => null()
     call nodeSatellite%destroy()
     call nodeHost     %destroy()
     deallocate(nodeSatellite)
     deallocate(nodeHost     )
  end do
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  !![
  <objectDestructor name="darkMatterHaloScale_" />
  <objectDestructor name="darkMatterProfileDMO_"/>
  !!]
end program Test_Satellite_Merging_Timescale_Poulton2021
