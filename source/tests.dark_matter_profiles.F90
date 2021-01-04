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

!% Contains a program which tests dark matter profiles.

program Test_Dark_Matter_Profiles
  !% Tests dark matter profiles.
  use :: Cosmology_Functions     , only : cosmologyFunctions            , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScale           , darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOBurkert
  use :: Events_Hooks            , only : eventsHooksInitialize
  use :: Galacticus_Display      , only : Galacticus_Verbosity_Level_Set, verbosityStandard
  use :: Galacticus_Nodes        , only : nodeClassHierarchyFinalize    , nodeClassHierarchyInitialize     , nodeComponentBasic                 , nodeComponentDarkMatterProfile, &
          &                               treeNode
  use :: ISO_Varying_String      , only : varying_string                , assignment(=)
  use :: Input_Parameters        , only : inputParameters
  use :: Node_Components         , only : Node_Components_Initialize    , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests              , only : Assert                        , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (treeNode                      ), pointer      :: node
  class           (nodeComponentBasic            ), pointer      :: basic
  class           (nodeComponentDarkMatterProfile), pointer      :: dmProfile
  class           (cosmologyFunctionsClass       ), pointer      :: cosmologyFunctions_
  class           (darkMatterHaloScaleClass      ), pointer      :: darkMatterHaloScale_
  double precision                                , parameter    :: concentration            =8.0d0                                                          , &
       &                                                            massVirial               =1.0d0
  double precision                                , dimension(7) :: radius                   =[0.125d0, 0.250d0, 0.500d0, 1.000d0, 2.000d0, 4.000d0, 8.000d0]
  double precision                                , dimension(7) :: mass                                                                                     , &
       &                                                            density                                                                                  , &
       &                                                            fourier
  type            (darkMatterProfileDMOBurkert   )               :: darkMatterProfileDMOBurkert_
  type            (varying_string                )               :: parameterFile
  type            (inputParameters               )               :: parameters
  integer                                                        :: i
  double precision                                               :: radiusScale

  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group    ('Dark matter profiles'           )
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/darkMatterProfiles.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Initialize node components.
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Create a node.
  node                         => treeNode                                     (                                )
  ! Create components.
  basic                        => node                       %basic            (autoCreate=.true.               )
  dmProfile                    => node                       %darkMatterProfile(autoCreate=.true.               )
  ! Get required objects.
  cosmologyFunctions_          => cosmologyFunctions                           (                                )
  darkMatterHaloScale_         => darkMatterHaloScale                          (                                )
  darkMatterProfileDMOBurkert_ =  darkMatterProfileDMOBurkert                  (            darkMatterHaloScale_)
  ! Set properties.
  call basic%timeSet     (cosmologyFunctions_%cosmicTime(1.0d0))
  call basic%massSet     (massVirial                           )
  ! Compute scale radius.
  radiusScale                  = +darkMatterHaloScale_       %virialRadius     (            node                ) &
       &                         /concentration
  call dmProfile%scaleSet(radiusScale                          )
  ! Test Burkert profile.
  call Unit_Tests_Begin_Group('Burkert profile')
  do i=1,7
     mass   (i)=darkMatterProfileDMOBurkert_%enclosedMass(node,      radiusScale*radius(i))
     density(i)=darkMatterProfileDMOBurkert_%density     (node,      radiusScale*radius(i))*radiusScale**3
     fourier(i)=darkMatterProfileDMOBurkert_%kSpace      (node,1.0d0/radiusScale/radius(i))
  end do
  call Assert(                        &
       &      'enclosed mass'       , &
       &      mass                  , &
       &      [                       &
       &       4.1583650166653620d-4, &
       &       2.9870571374971090d-3, &
       &       1.8812441757378840d-2, &
       &       8.9614051908432800d-2, &
       &       0.2805458115927276d+0, &
       &       0.5990982283029260d+0, &
       &       1.0000000000000000d+0  &
       &      ]                     , &
       &      relTol=1.0d-6           &
       &     )
  call Assert(                        &
       &      'density'             , &
       &      density               , &
       &      [                       &
       &       4.9082352873191440d-2, &
       &       4.2225259457083800d-2, &
       &       2.9909558782101030d-2, &
       &       1.4020105679109860d-2, &
       &       3.7386948477626290d-3, &
       &       6.5976967901693440d-4, &
       &       9.5863970455452000d-5  &
       &      ]                     , &
       &      relTol=1.0d-6           &
       &     )
  call Assert(                        &
       &      'fourier'             , &
       &      fourier               , &
       &      [                       &
       &       3.2941046717529910d-4, &
       &       6.7507368425877680d-3, &
       &       6.0387631952687390d-2, &
       &       0.2118984282100852d+0, &
       &       0.5171391325515731d+0, &
       &       0.8364966656849830d+0, &
       &       0.9557322027757890d+0  &
       &      ]                     , &
       &      relTol=1.0d-6           &
       &     )
  call Unit_Tests_End_Group       ()
  ! End unit tests.
  call Unit_Tests_End_Group       ()
  call Unit_Tests_Finish          ()
  ! Uninitialize node components.
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Test_Dark_Matter_Profiles
