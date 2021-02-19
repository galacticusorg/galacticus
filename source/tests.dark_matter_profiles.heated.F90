!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a program which tests the \glspl{node} implementation.

program Test_Dark_Matter_Profiles_Heated
  !% Tests the heated dark matter profile implementation. An isothermal dark matter halo is used, since analytic solutions are
  !% available for this case. Specifically, the initial radius in the unheated profile is given by $r_\mathrm{i}(r) =
  !% (r_\mathrm{h}^4/4 r^2 + r_\mathrm{h}^2)^{1/2}-r_\mathrm{h}^2/2 r$ where $r$ is the radius in the heated profile, and
  !% $r_\mathrm{h}=(\mathrm{G} M_\mathrm{v} / 2 Q r_\mathrm{v})^{1/2}$ is a characteristic heating radius. Here, $M_\mathrm{v}$,
  !% and $r_\mathrm{v}$ are the virial mass and radius of the halo respectively, and $Q r_\mathrm{i}^2$ is the specific heat input
  !% to the density profile, with $Q$ assumed to be a constant (as expected for tidal heating). Assuming no shell crossing, the
  !% enclosed mass in the final profile is simply $M(r) = M_\mathrm{v} r_\mathrm{i}(r)/r_\mathrm{v}$, from which the density of
  !% the final profile is found as $\rho(r) = (4 \pi r^2)^{-1} \mathrm{d} M(r) / \mathrm{d} r$.
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScale            , darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMOHeated     , darkMatterProfileDMOIsothermal, darkMatterProfileHeatingTidal
  use :: Dark_Matter_Profiles_Generic    , only : nonAnalyticSolversFallThrough
  use :: Display                         , only : displayVerbositySet            , verbosityLevelStandard
  use :: Events_Hooks                    , only : eventsHooksInitialize
  use :: Galacticus_Nodes                , only : nodeClassHierarchyFinalize     , nodeClassHierarchyInitialize  , nodeComponentBasic           , nodeComponentSatellite, &
          &                                       treeNode
  use :: ISO_Varying_String              , only : assignment(=)                  , varying_string
  use :: Input_Parameters                , only : inputParameters
  use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
  use :: Numerical_Constants_Math        , only : Pi
  use :: Unit_Tests                      , only : Assert                         , Unit_Tests_Begin_Group        , Unit_Tests_End_Group         , Unit_Tests_Finish
  implicit none
  double precision                                , parameter    :: time                           =13.8d00
  double precision                                , parameter    :: massVirial                     = 1.0d10
  double precision                                , parameter    :: heatingSpecific                = 1.0d06
  double precision                                , parameter    :: coefficientSecondOrder         = 0.0d00
  double precision                                , parameter    :: correlationVelocityRadius      =-1.0d00
  logical                                         , parameter    :: velocityDispersionApproximate  =.true.
  class           (nodeComponentBasic            ), pointer      :: basic
  class           (nodeComponentSatellite        ), pointer      :: satellite
  class           (darkMatterHaloScaleClass      ), pointer      :: darkMatterHaloScale_
  double precision                                , dimension(3) :: radiusVirialFractional         =[0.1d0,0.5d0,1.0d0]
  type            (treeNode                      )               :: node
  type            (varying_string                )               :: parameterFile
  type            (inputParameters               )               :: parameters
  type            (darkMatterProfileDMOIsothermal)               :: darkMatterProfileDMOIsothermal_
  type            (darkMatterProfileDMOHeated    )               :: darkMatterProfileDMOHeated_
  type            (darkMatterProfileHeatingTidal )               :: darkMatterProfileHeatingTidal_
  double precision                                               :: radiusVirial                                       , radiusHeated         , &
       &                                                            density                                            , densityAnalytic      , &
       &                                                            radiusInitial                                      , radiusInitialAnalytic, &
       &                                                            massEnclosed                                       , massEnclosedAnalytic , &
       &                                                            radius
  integer                                                        :: i
  character       (len=5                         )               :: radiusLabel

  ! Initialize.
  call displayVerbositySet(verbosityLevelStandard)
  parameterFile='testSuite/parameters/darkMatterProfileHeated.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  ! Initialize event hooks.
  call eventsHooksInitialize()
  call Unit_Tests_Begin_Group("Heated dark matter profiles")
  call nodeClassHierarchyInitialize(parameters)
  ! Create the dark matter profiles.
  darkMatterHaloScale_            => darkMatterHaloScale           (                                                                                                                                               )
  darkMatterProfileHeatingTidal_  =  darkMatterProfileHeatingTidal (coefficientSecondOrder       ,correlationVelocityRadius                                                                                        )
  darkMatterProfileDMOIsothermal_ =  darkMatterProfileDMOIsothermal(                                                                                            darkMatterHaloScale_                               )
  darkMatterProfileDMOHeated_     =  darkMatterProfileDMOHeated    (nonAnalyticSolversFallThrough,velocityDispersionApproximate,darkMatterProfileDMOIsothermal_,darkMatterHaloScale_,darkMatterProfileHeatingTidal_)
  ! Set up the node.
  basic     => node%basic    (autoCreate=.true.)
  satellite => node%satellite(autoCreate=.true.)
  call basic    %massSet                  (massVirial     )
  call basic    %timeSet                  (time           )
  call basic    %timeLastIsolatedSet      (time           )
  call satellite%tidalHeatingNormalizedSet(heatingSpecific)
  radiusVirial=darkMatterHaloScale_%virialRadius(node)
  ! Compute the characteristic radius for heating.
  radiusHeated=sqrt(gravitationalConstantGalacticus*massVirial/2.0d0/heatingSpecific/radiusVirial)
  ! Compute initial radius, enclosed mass, and density for a variety of radii and compare to the analytic solutions.
  do i=1,size(radiusVirialFractional)
     write (radiusLabel,'(f4.2)') radiusVirialFractional(i)
     radius               =+radiusVirial                                           &
          &                *radiusVirialFractional                   (     i     )
     radiusInitial        =+darkMatterProfileDMOHeated_%radiusInitial(node,radius)
     massEnclosed         =+darkMatterProfileDMOHeated_%enclosedMass (node,radius)
     density              =+darkMatterProfileDMOHeated_%density      (node,radius)
     radiusInitialAnalytic=+sqrt(                            &
          &                      +  radiusHeated         **4 &
          &                      /4.0d0                      &
          &                      /  radius               **2 &
          &                      +  radiusHeated         **2 &
          &                     )                            &
          &                -        radiusHeated         **2 &
          &                /2.0d0                            &
          &                /        radius
     massEnclosedAnalytic =+massVirial                       &
          &                *        radiusInitialAnalytic    &
          &                /        radiusVirial
     densityAnalytic      =+massVirial                       &
          &                /4.0d0                            &
          &                /Pi                               &
          &                /        radiusVirial             &
          &                /        radius               **2 &
          &                *(                                &
          &                  -      radiusHeated         **4 &
          &                  /4.0d0                          &
          &                  /      radius               **3 &
          &                  /sqrt(                          &
          &                        +radiusHeated         **4 &
          &                        /4.0d0                    &
          &                        /radius               **2 &
          &                        +radiusHeated         **2 &
          &                       )                          &
          &                  +      radiusHeated         **2 &
          &                  /2.0d0                          &
          &                  /      radius               **2 &
          &                 )
     call Assert('r='//trim(radiusLabel)//'rᵥ; initial radius',radiusInitial,radiusInitialAnalytic,relTol=1.0d-6)
     call Assert('r='//trim(radiusLabel)//'rᵥ; enclosed mass' ,massEnclosed ,massEnclosedAnalytic ,relTol=1.0d-6)
     call Assert('r='//trim(radiusLabel)//'rᵥ; density'       ,density      ,densityAnalytic      ,relTol=1.0d-6)
  end do
  call nodeClassHierarchyFinalize()
  call Unit_Tests_End_Group      ()
  call Unit_Tests_Finish         ()
end program Test_Dark_Matter_Profiles_Heated
