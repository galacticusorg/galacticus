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

  !!{
  Implements a class for applying multiple different timestepping criteria.
  !!}

  type, public :: multiMergerTreeEvolveTimestepList
     class(mergerTreeEvolveTimestepClass    ), pointer :: mergerTreeEvolveTimestep_ => null()
     type (multiMergerTreeEvolveTimestepList), pointer :: next                      => null()
  end type multiMergerTreeEvolveTimestepList

  !![
  <mergerTreeEvolveTimestep name="mergerTreeEvolveTimestepMulti">
   <description>A merger tree evolution timestepping class which takes the minimum over multiple other timesteppers.</description>
   <linkedList type="multiMergerTreeEvolveTimestepList" variable="mergerTreeEvolveTimesteps" next="next" object="mergerTreeEvolveTimestep_" objectType="mergerTreeEvolveTimestepClass"/>
   <deepCopy>
     <ignore  variables="mergerTreeEvolveTimestep_"/>
   </deepCopy>
   <stateStorable>
     <exclude variables="mergerTreeEvolveTimestep_"/>
   </stateStorable> 
  </mergerTreeEvolveTimestep>
  !!]
  type, extends(mergerTreeEvolveTimestepClass) :: mergerTreeEvolveTimestepMulti
     !!{
     Implementation of a merger tree evolution timestepping class which takes the minimum over multiple other timesteppers.
     !!}
     private
     type (multiMergerTreeEvolveTimestepList), pointer :: mergerTreeEvolveTimesteps => null()
     class(mergerTreeEvolveTimestepClass    ), pointer :: mergerTreeEvolveTimestep_ => null()
   contains
     final     ::                   multiDestructor
     procedure :: timeEvolveTo   => multiTimeEvolveTo
     procedure :: refuseToEvolve => multiRefuseToEvolve
  end type mergerTreeEvolveTimestepMulti

  interface mergerTreeEvolveTimestepMulti
     !!{
     Constructors for the \refClass{mergerTreeEvolveTimestepMulti} mergerTreeEvolveTimestep.
     !!}
     module procedure multiConstructorParameters
     module procedure multiConstructorInternal
  end interface mergerTreeEvolveTimestepMulti

contains

  function multiConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveTimestepMulti} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (mergerTreeEvolveTimestepMulti    )                :: self
    type   (inputParameters                  ), intent(inout) :: parameters
    type   (multiMergerTreeEvolveTimestepList), pointer       :: mergerTreeEvolveTimestep_
    integer                                                   :: i

    self %mergerTreeEvolveTimesteps => null()
    mergerTreeEvolveTimestep_       => null()
    do i=1,parameters%copiesCount('mergerTreeEvolveTimestep',zeroIfNotPresent=.true.)
       if (associated(mergerTreeEvolveTimestep_)) then
          allocate(mergerTreeEvolveTimestep_%next)
          mergerTreeEvolveTimestep_ => mergerTreeEvolveTimestep_%next
       else
          allocate(self%mergerTreeEvolveTimesteps)
          mergerTreeEvolveTimestep_ => self%mergerTreeEvolveTimesteps
       end if
       !![
       <objectBuilder class="mergerTreeEvolveTimestep" name="mergerTreeEvolveTimestep_%mergerTreeEvolveTimestep_" source="parameters" copy="i"/>
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="mergerTreeEvolveTimestep"/>
    !!]
    return
  end function multiConstructorParameters

  function multiConstructorInternal(mergerTreeEvolveTimesteps) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeEvolveTimestepMulti} merger tree evolution timestep class.
    !!}
    implicit none
    type(mergerTreeEvolveTimestepMulti    )                        :: self
    type(multiMergerTreeEvolveTimestepList), target, intent(in   ) :: mergerTreeEvolveTimesteps
    type(multiMergerTreeEvolveTimestepList), pointer               :: mergerTreeEvolveTimestep_

    self                     %mergerTreeEvolveTimesteps => mergerTreeEvolveTimesteps
    mergerTreeEvolveTimestep_                           => mergerTreeEvolveTimesteps
    do while (associated(mergerTreeEvolveTimestep_))
       !![
       <referenceCountIncrement owner="mergerTreeEvolveTimestep_" object="mergerTreeEvolveTimestep_"/>
       !!]
       mergerTreeEvolveTimestep_ => mergerTreeEvolveTimestep_%next
    end do
    return
  end function multiConstructorInternal

  subroutine multiDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeEvolveTimestepMulti} merger tree evolution timestep class.
    !!}
    implicit none
    type(mergerTreeEvolveTimestepMulti    ), intent(inout) :: self
    type(multiMergerTreeEvolveTimestepList), pointer       :: mergerTreeEvolveTimestep_, mergerTreeEvolveTimestepNext

    if (associated(self%mergerTreeEvolveTimesteps)) then
       mergerTreeEvolveTimestep_ => self%mergerTreeEvolveTimesteps
       do while (associated(mergerTreeEvolveTimestep_))
          mergerTreeEvolveTimestepNext => mergerTreeEvolveTimestep_%next
          !![
          <objectDestructor name="mergerTreeEvolveTimestep_%mergerTreeEvolveTimestep_"/>
          !!]
          deallocate(mergerTreeEvolveTimestep_)
          mergerTreeEvolveTimestep_ => mergerTreeEvolveTimestepNext
       end do
    end if
    return
  end subroutine multiDestructor

  double precision function multiTimeEvolveTo(self,timeEnd,node,task,taskSelf,report,lockNode,lockType)
    !!{
    Perform all mergerTreeEvolveTimesteps.
    !!}
    implicit none
    class           (mergerTreeEvolveTimestepMulti    ), intent(inout), target            :: self
    double precision                                   , intent(in   )                    :: timeEnd
    type            (treeNode                         ), intent(inout), target            :: node
    procedure       (timestepTask                     ), intent(  out), pointer           :: task
    class           (*                                ), intent(  out), pointer           :: taskSelf
    logical                                            , intent(in   )                    :: report
    type            (treeNode                         ), intent(  out), pointer, optional :: lockNode
    type            (varying_string                   ), intent(  out)         , optional :: lockType
    type            (multiMergerTreeEvolveTimestepList)               , pointer           :: mergerTreeEvolveTimestep_
    procedure       (timestepTask                     )               , pointer           :: task_
    class           (*                                )               , pointer           :: taskSelf_
    type            (treeNode                         )               , pointer           :: lockNode_
    type            (varying_string                   ), save                             :: lockType_
    !$omp threadprivate(lockType_)
    double precision                                                                      :: timeEvolveTo

    multiTimeEvolveTo               =  huge(0.0d0)
    task                            => null(     )
    taskSelf                        => null(     )
    if (present(lockNode)) lockNode => null(     )
    if (present(lockType)) lockType = ""
    mergerTreeEvolveTimestep_ => self%mergerTreeEvolveTimesteps
    do while (associated(mergerTreeEvolveTimestep_))
       timeEvolveTo=huge(0.0d0)
       !![
       <conditionalCall>
        <call>timeEvolveTo=mergerTreeEvolveTimestep_%mergerTreeEvolveTimestep_%timeEvolveTo(timeEnd,node,task_,taskSelf_,report{conditions})</call>
        <argument name="lockNode" value="lockNode_" condition="present(lockNode)"/>
        <argument name="lockType" value="lockType_" condition="present(lockType)"/>
       </conditionalCall>
       !!]
       if (timeEvolveTo < multiTimeEvolveTo) then
          self%mergerTreeEvolveTimestep_  => mergerTreeEvolveTimestep_%mergerTreeEvolveTimestep_
          multiTimeEvolveTo               =  timeEvolveTo
          task                            => task_
          taskSelf                        => taskSelf_
          if (present(lockNode)) lockNode => lockNode_
          if (present(lockType)) lockType =  lockType_
       end if
       mergerTreeEvolveTimestep_ => mergerTreeEvolveTimestep_%next
    end do
    return
  end function multiTimeEvolveTo

  logical function multiRefuseToEvolve(self,node)
    !!{
    Refuse to evolve if the timestep is too small.
    !!}
    implicit none
    class(mergerTreeEvolveTimestepMulti), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node

    multiRefuseToEvolve=self%mergerTreeEvolveTimestep_%refuseToEvolve(node)
    return
  end function multiRefuseToEvolve
