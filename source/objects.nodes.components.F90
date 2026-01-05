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
Contains a module which implements top-level functions for node components.
!!}

module Node_Components
  !!{
  Implements top-level functions for node components.
  !!}
  private
  public :: Node_Components_Initialize  , Node_Components_Thread_Initialize  , &
       &    Node_Components_Uninitialize, Node_Components_Thread_Uninitialize

  integer :: initializationCount=0, initializationThreadCount=0
  !$omp threadprivate(initializationThreadCount)

contains

  subroutine Node_Components_Initialize(parameters)
    !!{
    Perform initialization tasks for node components.
    !!}
    use :: Input_Parameters, only : inputParameters
    !![
    <include directive="nodeComponentInitializationTask" type="moduleUse">
    !!]
    include 'node_components.initialize.moduleUse.inc'
    !![
    </include>
    !!]
    implicit none
    type(inputParameters), intent(inout) :: parameters

    if (initializationCount == 0) then
       !![
       <include directive="nodeComponentInitializationTask" type="functionCall" functionType="void">
        <functionArgs>parameters</functionArgs>
       !!]
       include 'node_components.initialize.inc'
       !![
       </include>
       !!]
    end if
    initializationCount=initializationCount+1
    return
  end subroutine Node_Components_Initialize

  subroutine Node_Components_Thread_Initialize(parameters)
    !!{
    Perform per-thread initialization tasks for node components.
    !!}
    use :: Input_Parameters, only : inputParameters
    use :: Events_Hooks    , only : eventsHooksAtLevelToAllLevels   , calculationResetEvent   , openMPThreadBindingAllLevels
    use :: Galacticus_Nodes, only : massDistributionCalculationReset, massDistributionsDestroy, massDistributionsLast
    !![
    <include directive="nodeComponentThreadInitializationTask" type="moduleUse">
    !!]
    include 'node_components.threadInitialize.moduleUse.inc'
    !![
    </include>
    !!]
    implicit none
    type(inputParameters), intent(inout) :: parameters

    if (initializationThreadCount == 0) then
       ! Force all events that attach during this initialization to attach at all levels. This is necessary as the master thread in
       ! an OpenMP parallel section inherits objects from before the parallel region.
       call eventsHooksAtLevelToAllLevels(.true. )
       !![
       <include directive="nodeComponentThreadInitializationTask" type="functionCall" functionType="void">
        <functionArgs>parameters</functionArgs>
       !!]
       include 'node_components.threadInitialize.inc'
       !![
       </include>
       !!]
       ! Attach to an event that will be used to reset massDistributions of treeNodes during evolution.
       call calculationResetEvent%attach(massDistributionsLast,massDistributionCalculationReset,openMPThreadBindingAllLevels,label='massDistribution')
       call massDistributionsDestroy()
       ! Restore event hooking to standard behavior.
       call eventsHooksAtLevelToAllLevels(.false.)
    end if
    initializationThreadCount=initializationThreadCount+1
    return
  end subroutine Node_Components_Thread_Initialize

  subroutine Node_Components_Uninitialize()
    !!{
    Perform uninitialization tasks for node components.
    !!}
    implicit none

    initializationCount=initializationCount-1
    return
  end subroutine Node_Components_Uninitialize

  subroutine Node_Components_Thread_Uninitialize()
    !!{
    Perform per-thread uninitialization tasks for node components.
    !!}
    !![
    <include directive="nodeComponentThreadUninitializationTask" type="moduleUse">
    !!]
    include 'node_components.threadUninitialize.moduleUse.inc'
    !![
    </include>
    !!]
    use :: Events_Hooks    , only : calculationResetEvent
    use :: Galacticus_Nodes, only : massDistributionCalculationReset, massDistributionsLast
    implicit none

    initializationThreadCount=initializationThreadCount-1
    if (initializationThreadCount == 0) then
       !![
       <include directive="nodeComponentThreadUninitializationTask" type="functionCall" functionType="void">
       !!]
       include 'node_components.threadUninitialize.inc'
       !![
       </include>
       !!]
       if (calculationResetEvent%isAttached(massDistributionsLast,massDistributionCalculationReset)) call calculationResetEvent%detach(massDistributionsLast,massDistributionCalculationReset)
    end if
    return
  end subroutine Node_Components_Thread_Uninitialize

end module Node_Components
