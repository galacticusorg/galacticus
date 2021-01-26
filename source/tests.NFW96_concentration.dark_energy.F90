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

!% Contains a program which tests the \cite{navarro_structure_1996} halo concentration algorithm in a dark energy Universe.
!% Comparisons are made to the {\normalfont \ttfamily charden} code written by Julio Navarro.

program Test_NFW96_Concentration_Dark_Energy
  !% Tests the \cite{navarro_structure_1996} halo concentration algorithm in a dark energy Universe. Comparisons are made to the
  !% {\normalfont \ttfamily charden} code written by Julio Navarro.
  use :: Cosmology_Functions                 , only : cosmologyFunctions                 , cosmologyFunctionsClass
  use :: Cosmology_Parameters                , only : cosmologyParameters                , cosmologyParametersClass           , hubbleUnitsLittleH
  use :: Dark_Matter_Profiles_Concentration  , only : darkMatterProfileConcentration     , darkMatterProfileConcentrationClass
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Galacticus_Display                  , only : Galacticus_Verbosity_Level_Set     , verbosityStandard
  use :: Galacticus_Function_Classes_Destroys, only : Galacticus_Function_Classes_Destroy
  use :: Galacticus_Nodes                    , only : nodeClassHierarchyInitialize       , nodeComponentBasic                 , treeNode
  use :: ISO_Varying_String                  , only : varying_string                     , assignment(=)                      , operator(//)                       , char
  use :: Input_Parameters                    , only : inputParameters
  use :: Node_Components                     , only : Node_Components_Initialize         , Node_Components_Thread_Initialize  , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: String_Handling                     , only : operator(//)
  use :: Unit_Tests                          , only : Assert                             , Unit_Tests_Begin_Group             , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (treeNode                           )                         , pointer :: node
  class           (nodeComponentBasic                 )                         , pointer :: basic
  class           (cosmologyParametersClass           )                         , pointer :: cosmologyParameters_
  class           (cosmologyFunctionsClass            )                         , pointer :: cosmologyFunctions_
  class           (darkMatterProfileConcentrationClass)                         , pointer :: darkMatterProfileConcentration_
  integer                                              , dimension(6), parameter          :: chardenLogHaloMass       =[10,11,12,13,14,15]
  double precision                                     , dimension(6), parameter          :: chardenConcentrationZ0   =[10.2700200d00,9.0204391d00,7.8041310d00,6.6154380d00,5.4956946d00,4.4538398d00], chardenConcentrationZ3=[5.8715897d00,5.4417138d00,5.0239682d00,4.6186433d00,4.2366042d00,3.8884208d00]
  type            (varying_string                     )                                   :: message                                                                                                   , parameterFile
  type            (inputParameters                    )                                   :: parameters
  integer                                                                                 :: iMass
  double precision                                                                        :: ourConcentration

  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("NFW96 halo concentration algorithm: dark energy cosmology")

  ! Test NFW96 halo concentration algorithm in a dark energy universe.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/NFW96HaloConcentration/darkEnergy.xml'
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

  ! Get the default cosmology.
  cosmologyParameters_            => cosmologyParameters           ()
  ! Get the default cosmology functions object.
  cosmologyFunctions_             => cosmologyFunctions            ()
  ! Get the default concentrations object.
  darkMatterProfileConcentration_ => darkMatterProfileConcentration()

  ! Loop over halo masses.
  do iMass=1,size(chardenLogHaloMass)

     ! Set the mass of the original node.
     call basic%massSet((10.0d0**chardenLogHaloMass(iMass))/cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH))

     ! Compute and compare concentration at z=0.
     call basic%timeSet(cosmologyFunctions_%cosmicTime(1.00d0))
     ourConcentration=darkMatterProfileConcentration_%concentration(node)
     message="10^"
     message=message//chardenLogHaloMass(iMass)//" M⊙/h halo concentration at z=0"
     call Assert(char(message),ourConcentration,chardenConcentrationZ0(iMass),relTol=0.02d0)

     ! Compute and compare concentration at z=3.
     call basic%timeSet(cosmologyFunctions_%cosmicTime(0.25d0))
     ourConcentration=darkMatterProfileConcentration_%concentration(node)
     message="10^"
     message=message//chardenLogHaloMass(iMass)//" M⊙/h halo concentration at z=3"
     call Assert(char(message),ourConcentration,chardenConcentrationZ3(iMass),relTol=0.01d0)

  end do
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call Galacticus_Function_Classes_Destroy()
end program Test_NFW96_Concentration_Dark_Energy
