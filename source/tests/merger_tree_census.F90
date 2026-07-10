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

program Tests_Merger_Tree_Census
  !!{RST
  Tests of the merger tree census interface (the ``countTrees`` and ``treeMasses`` methods of the ``mergerTreeConstructor`` class),
  which report the number and root masses of the trees that will be constructed. The ``build`` constructor is tested here as it
  knows this information up-front.
  !!}
  use            :: Display                , only : displayVerbositySet    , verbosityLevelStandard
  use            :: Events_Hooks           , only : eventsHooksInitialize
  use            :: Functions_Global_Utilities, only : Functions_Global_Set
  use            :: Galacticus_Nodes       , only : nodeClassHierarchyFinalize, nodeClassHierarchyInitialize
  use            :: Input_Parameters       , only : inputParameters
  use            :: ISO_Varying_String     , only : var_str
  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Merger_Tree_Construction, only : mergerTreeConstructorClass
  use            :: Node_Components         , only : Node_Components_Initialize, Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use            :: Unit_Tests             , only : Assert                 , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (inputParameters          )               :: parameters
  class           (mergerTreeConstructorClass), pointer     :: mergerTreeConstructor_
  integer         (c_size_t                  )              :: countTrees
  double precision                            , allocatable, dimension(:) :: masses
  double precision                            , parameter   :: massTreeMinimum=1.0d10, massTreeMaximum=1.0d13

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group('Merger tree census')
  ! Read in controlling parameters and initialize.
  parameters=inputParameters(var_str('testSuite/parameters/mergerTreeCensus.xml'))
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Build the merger tree constructor.
  !![
  <objectBuilder class="mergerTreeConstructor" name="mergerTreeConstructor_" source="parameters"/>
  !!]
  ! Query the census interface.
  countTrees=mergerTreeConstructor_%countTrees()
  call mergerTreeConstructor_%treeMasses(masses)
  ! The build constructor knows the tree count and masses up-front, so the count should be positive, the masses array should be
  ! allocated with a size matching the count, and all masses should lie within the requested mass range.
  call Assert('tree count is known (positive)'            ,countTrees > 0_c_size_t                     ,.true.)
  call Assert('tree masses array is allocated'            ,allocated(masses)                           ,.true.)
  if (allocated(masses)) then
     call Assert('tree masses array size matches count'   ,size(masses,kind=c_size_t) == countTrees    ,.true.)
     call Assert('all tree masses are within range'       ,all(masses >= massTreeMinimum .and. masses <= massTreeMaximum),.true.)
  end if
  ! Clean up.
  !![
  <objectDestructor name="mergerTreeConstructor_"/>
  !!]
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  call parameters%destroy                 ()
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
end program Tests_Merger_Tree_Census
