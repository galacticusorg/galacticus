!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  use :: Posterior_Sampling_Simulation, only : posteriorSampleSimulation, posteriorSampleSimulationClass

  !![
  <task name="taskPosteriorSample">
   <description>A task which performs sampling from a posterior distribution.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskPosteriorSample
     !!{
     Implementation of a task which performs sampling from a posterior distribution.
     !!}
     private
     class  (posteriorSampleSimulationClass), pointer :: posteriorSampleSimulation_    => null()
     logical                                          :: nodeClassHierarchyInitialized =  .false., initializeNodeClassHierarchy
   contains
     final     ::            posteriorSampleDestructor
     procedure :: perform => posteriorSamplePerform
  end type taskPosteriorSample

  interface taskPosteriorSample
     !!{
     Constructors for the \refClass{taskPosteriorSample} task.
     !!}
     module procedure posteriorSampleConstructorParameters
     module procedure posteriorSampleConstructorInternal
  end interface taskPosteriorSample

contains

  function posteriorSampleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskPosteriorSample} task class which takes a parameter set as input.
    !!}
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize
    use :: Input_Parameters, only : inputParameter              , inputParameters
    use :: Node_Components , only : Node_Components_Initialize  , Node_Components_Thread_Initialize
    implicit none
    type   (taskPosteriorSample           )                :: self
    type   (inputParameters               ), intent(inout) :: parameters
    class  (posteriorSampleSimulationClass), pointer       :: posteriorSampleSimulation_
    type   (inputParameters               ), pointer       :: parametersRoot
    logical                                                :: initializeNodeClassHierarchy

    !![
    <inputParameter>
      <name>initializeNodeClassHierarchy</name>
      <description>If true then initialize the node class hierarchy in the posterior sampling class. This should be set to false if the likelihood function will instead perform this action.</description>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (initializeNodeClassHierarchy) then
       if (associated(parameters%parent)) then
          parametersRoot => parameters%parent
          do while (associated(parametersRoot%parent))
             parametersRoot => parametersRoot%parent
          end do
       else
          parametersRoot => parameters
       end if
       call nodeClassHierarchyInitialize     (parametersRoot)
       call Node_Components_Initialize       (parametersRoot)
       call Node_Components_Thread_Initialize(parametersRoot)
    end if
    !![
    <objectBuilder class="posteriorSampleSimulation" name="posteriorSampleSimulation_" source="parameters"/>
    !!]
    self=taskPosteriorSample(posteriorSampleSimulation_)
    self%initializeNodeClassHierarchy=initializeNodeClassHierarchy
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="posteriorSampleSimulation_"/>
    !!]
    if (initializeNodeClassHierarchy) self%nodeClassHierarchyInitialized=.true.
    return
  end function posteriorSampleConstructorParameters

  function posteriorSampleConstructorInternal(posteriorSampleSimulation_) result(self)
    !!{
    Internal constructor for the \refClass{taskPosteriorSample} task class.
    !!}
    implicit none
    type (taskPosteriorSample           )                        :: self
    class(posteriorSampleSimulationClass), intent(in   ), target :: posteriorSampleSimulation_
    !![
    <constructorAssign variables="*posteriorSampleSimulation_"/>
    !!]

    self%initializeNodeClassHierarchy=.false.
    return
  end function posteriorSampleConstructorInternal

  subroutine posteriorSampleDestructor(self)
    !!{
    Destructor for the \refClass{taskPosteriorSample} task class.
    !!}
    use :: Node_Components, only : Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
    implicit none
    type(taskPosteriorSample), intent(inout) :: self

    !![
    <objectDestructor name="self%posteriorSampleSimulation_"/>
    !!]
    if (self%nodeClassHierarchyInitialized) then
       call Node_Components_Thread_Uninitialize()
       call Node_Components_Uninitialize       ()
    end if
    return
  end subroutine posteriorSampleDestructor

  subroutine posteriorSamplePerform(self,status)
    !!{
    Perform the posterior sampling.
    !!}
    use :: Display         , only : displayIndent     , displayUnindent
    use :: Error, only : errorStatusSuccess
    implicit none
    class  (taskPosteriorSample), intent(inout), target   :: self
    integer                     , intent(  out), optional :: status

    call displayIndent('Begin task: posterior sampling')
    call self%posteriorSampleSimulation_%simulate()
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: posterior sampling' )
    return
  end subroutine posteriorSamplePerform
