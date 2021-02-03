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

!% Contains a program which tests the \cite{correa_accretion_2015} concentration-mass relation.

program Test_Correa2015_Concentration
  !% Tests the \cite{correa_accretion_2015} concentration-mass relation.
  use :: Cosmology_Functions                 , only : cosmologyFunctions                 , cosmologyFunctionsClass
  use :: Dark_Matter_Profiles_Concentration  , only : darkMatterProfileConcentration     , darkMatterProfileConcentrationClass
  use :: Display                             , only : displayVerbositySet                , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Galacticus_Function_Classes_Destroys, only : Galacticus_Function_Classes_Destroy
  use :: Galacticus_Nodes                    , only : nodeClassHierarchyInitialize       , nodeComponentBasic                 , treeNode
  use :: ISO_Varying_String                  , only : assignment(=)                      , varying_string
  use :: Input_Parameters                    , only : inputParameters
  use :: Node_Components                     , only : Node_Components_Initialize         , Node_Components_Thread_Initialize  , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests                          , only : Assert                             , Unit_Tests_Begin_Group             , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (treeNode                           ), pointer      :: node
  class           (nodeComponentBasic                 ), pointer      :: basic
  class           (cosmologyFunctionsClass            ), pointer      :: cosmologyFunctions_
  class           (darkMatterProfileConcentrationClass), pointer      :: darkMatterProfileConcentration_
  type            (varying_string                     )               :: parameterFile
  type            (inputParameters                    )               :: parameters
  double precision                                     , dimension(3) :: concentrationTarget            , mass         , &
       &                                                                 redshift                       , concentration
  integer                                                             :: i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Correa et al. 2015 concentration-mass relation")
  ! Test Correa et al. 2015 algorithm.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/Correa2015.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Create a node.
  node  => treeNode      (                 )
  ! Get the basic component.
  basic => node    %basic(autoCreate=.true.)
  ! Get required objects.
  cosmologyFunctions_             => cosmologyFunctions            ()
  darkMatterProfileConcentration_ => darkMatterProfileConcentration()
  ! Specify halo redshifts and masses.
  mass               =[                  &
       &               1.00000000000d12, &
       &               5.64259219617d11, &
       &               2.96360394164d11  &
       &              ]
  redshift           =[                  &
       &               0.00000000000d00, &
       &               1.00000000000d00, &
       &               2.00000000000d00  &
       &              ]
  ! Specify concentration target.
  concentrationTarget=[                  &
       &               8.883912548300d0, &
       &               6.346877625700d0, &
       &               4.979205154090d0  &
       &              ]
  ! Iterate over redshifts.
  do i=1,size(redshift)
     ! Set node properties.
     call basic%massSet(                                                              &
          &                                                               mass    (i) &
          &            )
     call basic%timeSet(                                                              &
          &             cosmologyFunctions_ %cosmicTime                 (             &
          &              cosmologyFunctions_%expansionFactorFromRedshift (            &
          &                                                               redshift(i) &
          &                                                              )            &
          &                                                             )             &
          &            )
     ! Evaluate concentration.
     concentration(i)=darkMatterProfileConcentration_%concentration(node)
  end do
  call Assert('concentration',concentration,concentrationTarget,relTol=2.0d-3)
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call Galacticus_Function_Classes_Destroy()

end program Test_Correa2015_Concentration
