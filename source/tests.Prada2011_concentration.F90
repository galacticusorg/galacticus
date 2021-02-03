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

!% Contains a program which tests the \cite{prada_halo_2011} halo concentration algorithm.

program Test_Prada2011_Concentration
  !% Tests the \cite{prada_halo_2011} halo concentration algorithm. Values of concentration were read from their Figure~12.
  use :: Cosmology_Functions                 , only : cosmologyFunctions                 , cosmologyFunctionsClass
  use :: Cosmology_Parameters                , only : cosmologyParameters                , cosmologyParametersClass           , hubbleUnitsLittleH
  use :: Dark_Matter_Profiles_Concentration  , only : darkMatterProfileConcentration     , darkMatterProfileConcentrationClass
  use :: Display                             , only : displayVerbositySet                , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Galacticus_Function_Classes_Destroys, only : Galacticus_Function_Classes_Destroy
  use :: Galacticus_Nodes                    , only : nodeClassHierarchyInitialize       , nodeComponentBasic                 , treeNode
  use :: ISO_Varying_String                  , only : assignment(=)                      , char                               , varying_string
  use :: Input_Parameters                    , only : inputParameters
  use :: Node_Components                     , only : Node_Components_Initialize         , Node_Components_Thread_Initialize  , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests                          , only : Assert                             , Unit_Tests_Begin_Group             , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (treeNode                           )                                 , pointer :: node
  class           (nodeComponentBasic                 )                                 , pointer :: basic
  class           (cosmologyParametersClass           )                                 , pointer :: cosmologyParameters_
  class           (darkMatterProfileConcentrationClass)                                 , pointer :: darkMatterProfileConcentration_
  type            (varying_string                     )                                           :: message                                                        , parameterFile
  integer                                                                    , parameter          :: massCount                =4
  double precision                                     , dimension(massCount), parameter          :: logMass                  =[11.000d0,12.000d0,13.000d0,14.000d0]
  double precision                                     , dimension(massCount), parameter          :: pradaLogConcentration    =[0.966d0,0.887d0,0.804d0,0.728d0]
  double precision                                     , dimension(massCount)                     :: ourLogConcentration
  class           (cosmologyFunctionsClass )                                            , pointer :: cosmologyFunctions_
  integer                                                                                         :: iMass
  type            (inputParameters         )                                                      :: parameters

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Prada2011 halo concentration algorithm")

  ! Test Prada2011 halo concentration algorithm.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/Prada2011HaloConcentration/testParameters.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)

  ! Create a node.
  node                            => treeNode                            (                 )
  ! Get the basic component.
  basic                           => node                          %basic(autoCreate=.true.)
  ! Get the default cosmology functions object.
  cosmologyFunctions_             => cosmologyFunctions                  (                 )
  ! Get the default concentrations object.
  darkMatterProfileConcentration_ => darkMatterProfileConcentration      (                 )
  ! Get the default cosmology.
  cosmologyParameters_            => cosmologyParameters                 (                 )

  ! Set the time for the node.
  call basic%timeSet(cosmologyFunctions_%cosmicTime(1.00d0))

  ! Loop over halo masses
  do iMass=1,massCount

     ! Set the mass of the original node.
     call basic%massSet(10.0d0**logMass(iMass)/cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH))

     ! Compute and compare concentration at z=0.
     ourLogConcentration(iMass)=log10(darkMatterProfileConcentration_%concentration(node))

  end do

  ! Check that results are as expected.
  message="Halo concentration at z=0"
  call Assert(char(message),ourLogConcentration,pradaLogConcentration,absTol=0.01d0)

  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call Galacticus_Function_Classes_Destroy()

end program Test_Prada2011_Concentration
