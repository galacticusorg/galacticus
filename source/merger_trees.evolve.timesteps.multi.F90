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

  !% Implements a class for applying multiple different timestepping criteria.
  
  type, public :: multiMergerTreeEvolveTimestepList
     class(mergerTreeEvolveTimestepClass    ), pointer :: mergerTreeEvolveTimestep_
     type (multiMergerTreeEvolveTimestepList), pointer :: next                      => null()
  end type multiMergerTreeEvolveTimestepList
  
  !# <mergerTreeEvolveTimestep name="mergerTreeEvolveTimestepMulti" defaultThreadPrivate="yes">
  !#  <description>A merger tree evolution timestepping class which takes the minimum over multiple other timesteppers.</description>
  !# </mergerTreeEvolveTimestep>
  type, extends(mergerTreeEvolveTimestepClass) :: mergerTreeEvolveTimestepMulti
     !% Implementation of a merger tree evolution timestepping class which takes the minimum over multiple other timesteppers.
     private
     type(multiMergerTreeEvolveTimestepList), pointer :: mergerTreeEvolveTimesteps
   contains
     final     ::                       multiDestructor
     procedure :: timeEvolveTo       => multiTimeEvolveTo
     procedure :: deepCopy           => multiDeepCopy
  end type mergerTreeEvolveTimestepMulti

  interface mergerTreeEvolveTimestepMulti
     !% Constructors for the {\normalfont \ttfamily multi} mergerTreeEvolveTimestep.
     module procedure multiConstructorParameters
     module procedure multiConstructorInternal
  end interface mergerTreeEvolveTimestepMulti

contains

  function multiConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily multi} merger tree evolution timestep class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type   (mergerTreeEvolveTimestepMulti      )                :: self
    type   (inputParameters                    ), intent(inout) :: parameters
    type   (multiMergerTreeEvolveTimestepList  ), pointer       :: mergerTreeEvolveTimestep_
    integer                                                     :: i

    self %mergerTreeEvolveTimesteps => null()
    mergerTreeEvolveTimestep_       => null()
    do i=1,parameters%copiesCount('mergerTreeEvolveTimestepMethod',zeroIfNotPresent=.true.)
       if (associated(mergerTreeEvolveTimestep_)) then
          allocate(mergerTreeEvolveTimestep_%next)
          mergerTreeEvolveTimestep_ => mergerTreeEvolveTimestep_%next
       else
          allocate(self%mergerTreeEvolveTimesteps)
          mergerTreeEvolveTimestep_ => self%mergerTreeEvolveTimesteps
       end if
       !# <objectBuilder class="mergerTreeEvolveTimestep" name="mergerTreeEvolveTimestep_%mergerTreeEvolveTimestep_" source="parameters" copy="i"/>
    end do
    return
  end function multiConstructorParameters
  
  function multiConstructorInternal(mergerTreeEvolveTimesteps) result(self)
    !% Internal constructor for the {\normalfont \ttfamily multi} merger tree evolution timestep class.
    implicit none
    type(mergerTreeEvolveTimestepMulti    )                        :: self
    type(multiMergerTreeEvolveTimestepList), target, intent(in   ) :: mergerTreeEvolveTimesteps
    !# <constructorAssign variables="*mergerTreeEvolveTimesteps"/>

    return
  end function multiConstructorInternal
  
  subroutine multiDestructor(self)
    !% Destructor for the {\normalfont \ttfamily multi} merger tree evolution timestep class.
    implicit none
    type(mergerTreeEvolveTimestepMulti    ), intent(inout) :: self
    type(multiMergerTreeEvolveTimestepList), pointer       :: mergerTreeEvolveTimestep_, mergerTreeEvolveTimestepNext

    if (associated(self%mergerTreeEvolveTimesteps)) then
       mergerTreeEvolveTimestep_ => self%mergerTreeEvolveTimesteps
       do while (associated(mergerTreeEvolveTimestep_))
          mergerTreeEvolveTimestepNext => mergerTreeEvolveTimestep_%next
          !# <objectDestructor name="mergerTreeEvolveTimestep_%mergerTreeEvolveTimestep_"/>
          deallocate(mergerTreeEvolveTimestep_      )
          mergerTreeEvolveTimestep_ => mergerTreeEvolveTimestepNext
       end do
    end if
    return 
  end subroutine multiDestructor

  double precision function multiTimeEvolveTo(self,node,task,taskSelf,report,lockNode,lockType)
    !% Perform all mergerTreeEvolveTimesteps.
    implicit none
    class           (mergerTreeEvolveTimestepMulti    ), intent(inout), target  :: self
    type            (treeNode                         ), intent(inout), target  :: node
    procedure       (timestepTask                     ), intent(  out), pointer :: task
    class           (*                                ), intent(  out), pointer :: taskSelf
    logical                                            , intent(in   )          :: report
    type            (treeNode                         ), intent(  out), pointer :: lockNode
    type            (varying_string                   ), intent(  out)          :: lockType
    type            (multiMergerTreeEvolveTimestepList)               , pointer :: mergerTreeEvolveTimestep_
    procedure       (timestepTask                     )               , pointer :: task_
    class           (*                                )               , pointer :: taskSelf_
    type            (treeNode                         )               , pointer :: lockNode_
    type            (varying_string                   )                         :: lockType_
    double precision                                                            :: timeEvolveTo

    multiTimeEvolveTo =  huge(0.0d0)
    task              => null(     )
    taskSelf          => null(     )
    lockNode          => null(     )
    lockType          = ""
    mergerTreeEvolveTimestep_ => self%mergerTreeEvolveTimesteps
    do while (associated(mergerTreeEvolveTimestep_))
       timeEvolveTo=mergerTreeEvolveTimestep_%mergerTreeEvolveTimestep_%timeEvolveTo(node,task_,taskSelf_,report,lockNode_,lockType_)
       if (timeEvolveTo < multiTimeEvolveTo) then
          multiTimeEvolveTo =  timeEvolveTo
          task              => task_
          taskSelf          => taskSelf_
          lockNode          => lockNode_
          lockType          =  lockType_
       end if
       mergerTreeEvolveTimestep_ => mergerTreeEvolveTimestep_%next
    end do
    return
  end function multiTimeEvolveTo

  subroutine multiDeepCopy(self,destination)
    !% Perform a deep copy for the {\normalfont \ttfamily multi} merger tree evolution timestep class.
    use Galacticus_Error
    implicit none
    class(mergerTreeEvolveTimestepMulti    ), intent(inout) :: self
    class(mergerTreeEvolveTimestepClass    ), intent(  out) :: destination
    type (multiMergerTreeEvolveTimestepList), pointer       :: mergerTreeEvolveTimestep_   , mergerTreeEvolveTimestepDestination_, &
         &                                                     mergerTreeEvolveTimestepNew_

    call self%mergerTreeEvolveTimestepClass%deepCopy(destination)
    select type (destination)
    type is (mergerTreeEvolveTimestepMulti)
       destination%mergerTreeEvolveTimesteps => null          ()
       mergerTreeEvolveTimestepDestination_  => null          ()
       mergerTreeEvolveTimestep_             => self%mergerTreeEvolveTimesteps
       do while (associated(mergerTreeEvolveTimestep_))
          allocate(mergerTreeEvolveTimestepNew_)
          if (associated(mergerTreeEvolveTimestepDestination_)) then
             mergerTreeEvolveTimestepDestination_%next       => mergerTreeEvolveTimestepNew_
             mergerTreeEvolveTimestepDestination_            => mergerTreeEvolveTimestepNew_             
          else
             destination          %mergerTreeEvolveTimesteps => mergerTreeEvolveTimestepNew_
             mergerTreeEvolveTimestepDestination_            => mergerTreeEvolveTimestepNew_
          end if
          allocate(mergerTreeEvolveTimestepNew_%mergerTreeEvolveTimestep_,mold=mergerTreeEvolveTimestep_%mergerTreeEvolveTimestep_)
          call mergerTreeEvolveTimestep_%mergerTreeEvolveTimestep_%deepCopy(mergerTreeEvolveTimestepNew_%mergerTreeEvolveTimestep_)
          mergerTreeEvolveTimestep_ => mergerTreeEvolveTimestep_%next
       end do       
    class default
       call Galacticus_Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine multiDeepCopy
