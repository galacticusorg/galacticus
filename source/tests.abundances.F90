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

!!{
Contains a program to test abundances objects functions.
!!}

program Test_Abundances
  !!{
  Test abundances objects.
  !!}
  use :: Abundances_Structure      , only : abundances                  , max
  use :: Display                   , only : displayVerbositySet         , verbosityLevelStandard
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: ISO_Varying_String        , only : assignment(=)               , varying_string
  use :: Input_Parameters          , only : inputParameters
  use :: Galacticus_Nodes          , only : nodeClassHierarchyInitialize
  use :: Node_Components           , only : Node_Components_Initialize
  use :: Unit_Tests                , only : Assert                      , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (varying_string )               :: parameterFile
  type            (abundances     )               :: abundances1    , abundances2, abundances3
  double precision                 , dimension(5) :: abundancesArray
  type            (inputParameters)               :: parameters

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Abundances objects functions")

  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/abundances/testAbundances.xml'
  parameters=inputParameters(parameterFile)
  call Functions_Global_Set        (          )
  call nodeClassHierarchyInitialize(parameters)
  call Node_Components_Initialize  (parameters)

  ! Initialize abundances.
  call abundances1%deserialize([1.0d0,2.0d0,3.0d0,4.0d0,5.0d0])
  call abundances2%deserialize([5.0d0,4.0d0,3.0d0,2.0d0,1.0d0])
  abundances3=max(abundances1,abundances2)
  call abundances3%serialize(abundancesArray)
  call Assert('max() function',abundancesArray,[5.0d0,4.0d0,3.0d0,4.0d0,5.0d0])

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Abundances
