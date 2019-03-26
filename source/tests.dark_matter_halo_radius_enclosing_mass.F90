!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!+    Contributions to this file made by: Xiaolong Du, Andrew Benson.

!% Contains a program which tests the calculation of dark matter halo radius enclosing a given mass.

program Test_Dark_Matter_Halo_Radius_Enclosing_Mass
  !% Tests the calculation of dark matter halo radius enclosing a given mass.
  use ISO_Varying_String
  use Input_Parameters
  use Galacticus_Nodes       , only : treeNode                           , nodeComponentBasic               , nodeComponentSatellite        , &
       &                              nodeComponentDarkMatterProfile     , nodeClassHierarchyInitialize     , nodeClassHierarchyFinalize
  use Node_Components        , only : Node_Components_Initialize         , Node_Components_Thread_Initialize, Node_Components_Uninitialize  , &
       &                              Node_Components_Thread_Uninitialize
  use Dark_Matter_Halo_Scales
  use Dark_Matter_Profiles
  use Unit_Tests
  use Galacticus_Display
  implicit none
  type            (treeNode                             )               :: node
  class           (nodeComponentBasic                   ), pointer      :: basic
  class           (nodeComponentSatellite               ), pointer      :: satellite
  class           (nodeComponentDarkMatterProfile       ), pointer      :: dmProfile
  class           (darkMatterHaloScaleClass             ), pointer      :: darkMatterHaloScale_
  class           (darkMatterProfileClass               ), pointer      :: darkMatterProfile_
  type            (darkMatterProfileNFW                 ), target       :: darkMatterProfileNFW_
  type            (darkMatterProfileBurkert             ), target       :: darkMatterProfileBurkert_
  type            (darkMatterProfileTruncated           ), target       :: darkMatterProfileTruncated_
  type            (darkMatterProfileTruncatedExponential), target       :: darkMatterProfileTruncatedExponential_
  type            (darkMatterProfileHeated              ), target       :: darkMatterProfileHeated_
  type            (darkMatterProfileHeatingTidal        )               :: darkMatterProfileHeatingTidal_
  double precision                                                      :: time=13.8d0
  double precision                                       , parameter    :: massVirial=1.0d10, concentration=8.0d0
  double precision                                                      :: radiusVirial     , radiusScale
  double precision                                       , dimension(7) :: radiusOverVirialRadius=[0.125d0, 0.250d0, 0.500d0, 1.000d0, 2.000d0, 4.000d0, 8.000d0]
  double precision                                       , dimension(7) :: radius           , radiusRoot
  double precision                                       , dimension(7) :: mass
  logical                                                               :: unimplementedIsFatal=.true.
  double precision                                                      :: radiusFractionalTruncateMinimum=2.0d0, radiusFractionalTruncateMaximum=8.0d0
  double precision                                                      :: radiusFractionalDecay=0.06d0, alpha=1.0d0, beta=3.0d0, gamma=1.0d0
  double precision                                       , parameter    :: heatingSpecific=1.0d06
  type            (varying_string                       )               :: parameterFile
  type            (inputParameters                      )               :: parameters
  integer                                                               :: i, j

  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group('Dark matter halo radius enclosing a give mass')
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/darkMatterHaloRadiusEnclosingMass.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  ! Initialize node components.
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Create the dark matter profiles.
  darkMatterHaloScale_                   => darkMatterHaloScale                  (                                                               )
  darkMatterProfileNFW_                  =  darkMatterProfileNFW                 (                                           darkMatterHaloScale_)
  darkMatterProfileBurkert_              =  darkMatterProfileBurkert             (                                           darkMatterHaloScale_)
  darkMatterProfileTruncated_            =  darkMatterProfileTruncated           (radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum, &
       &                                                                          unimplementedIsFatal,darkMatterProfileNFW_,darkMatterHaloScale_)
  darkMatterProfileTruncatedExponential_ =  darkMatterProfileTruncatedExponential(radiusFractionalDecay, alpha        , beta,gamma               , &
       &                                                                          unimplementedIsFatal,darkMatterProfileNFW_,darkMatterHaloScale_)
  darkMatterProfileHeatingTidal_         =  darkMatterProfileHeatingTidal        (                                                               )
  darkMatterProfileHeated_               =  darkMatterProfileHeated              (unimplementedIsFatal,darkMatterProfileNFW_,darkMatterHaloScale_, &
       &                                                                          darkMatterProfileHeatingTidal_                                 )
  ! Set up the node.
  basic     => node%basic                 (autoCreate=.true.)
  satellite => node%satellite             (autoCreate=.true.)
  dmProfile => node%darkMatterProfile     (autoCreate=.true.)
  call basic    %massSet                  (massVirial       )
  call basic    %timeSet                  (time             )
  call basic    %timeLastIsolatedSet      (time             )
  call satellite%tidalHeatingNormalizedSet(heatingSpecific  )
  ! Get the virial raidus.
  radiusVirial=darkMatterHaloScale_%virialRadius(node)
  ! Compute scale radius.
  radiusScale =radiusVirial/concentration
  call dmProfile%scaleSet(radiusScale)
  ! Test different dark matter profiles.
  radius      =radiusOverVirialRadius*radiusVirial
  do i=1,5
     select case (i)
     case (1)
        call Unit_Tests_Begin_Group('NFW profile'                       )
        darkMatterProfile_ => darkMatterProfileNFW_
     case (2)
        call Unit_Tests_Begin_Group('Burkert profile'                   )
        darkMatterProfile_ => darkMatterProfileBurkert_
     case (3)
        call Unit_Tests_Begin_Group('Truncated profile'                 )
        darkMatterProfile_ => darkMatterProfileTruncated_
     case (4)
        call Unit_Tests_Begin_Group('Exponentially truncated profile'   )
        darkMatterProfile_ => darkMatterProfileTruncatedExponential_
     case default
        call Unit_Tests_Begin_Group('Heated profile'                    )
        darkMatterProfile_ => darkMatterProfileHeated_
     end select
     do j=1,7
        mass      (j)=darkMatterProfile_%enclosedMass       (node, radius(j))
        radiusRoot(j)=darkMatterProfile_%radiusEnclosingMass(node, mass  (j)) 
     end do
     call Assert('radius enclosing a given mass',radius,radiusRoot,relTol=1.0d-6)
     call Unit_Tests_End_Group()
  end do
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  ! Uninitialize node components.
  call nodeClassHierarchyFinalize         ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Test_Dark_Matter_Halo_Radius_Enclosing_Mass
