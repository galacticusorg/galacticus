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

  type, public :: multiTaskList
     class(taskClass    ), pointer :: task_ => null()
     type (multiTaskList), pointer :: next  => null()
  end type multiTaskList

  !![
  <task name="taskMulti">
   <description>A task which performs multiple other tasks.</description>
   <linkedList type="multiTaskList" variable="tasks" next="next" object="task_" objectType="taskClass"/>
  </task>
  !!]
  type, extends(taskClass) :: taskMulti
     !!{
     Implementation of a task which performs multiple other tasks.
     !!}
     private
     type(multiTaskList), pointer :: tasks => null()
   contains
     final     ::                       multiDestructor
     procedure :: perform            => multiPerform
     procedure :: requiresOutputFile => multiRequiresOutputFile
  end type taskMulti

  interface taskMulti
     !!{
     Constructors for the \refClass{taskMulti} task.
     !!}
     module procedure multiConstructorParameters
     module procedure multiConstructorInternal
  end interface taskMulti

contains

  function multiConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskMulti} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (taskMulti      )                :: self
    type   (inputParameters), intent(inout) :: parameters
    type   (multiTaskList  ), pointer       :: task_
    integer                                 :: i

    self %tasks => null()
    task_       => null()
    do i=1,parameters%copiesCount('task',zeroIfNotPresent=.true.)
       if (associated(task_)) then
          allocate(task_%next)
          task_ => task_%next
       else
          allocate(self%tasks)
          task_ => self%tasks
       end if
       !![
       <objectBuilder class="task" name="task_%task_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="task"/>
    !!]
    return
  end function multiConstructorParameters

  function multiConstructorInternal(tasks) result(self)
    !!{
    Internal constructor for the \refClass{taskMulti} task class.
    !!}
    implicit none
    type(taskMulti    )                        :: self
    type(multiTaskList), target, intent(in   ) :: tasks
    type(multiTaskList), pointer               :: task_

    self %tasks => tasks
    task_       => tasks
    do while (associated(task_))
       !![
       <referenceCountIncrement owner="task_" object="task_"/>
       !!]
       task_ => task_%next
    end do
    return
  end function multiConstructorInternal

  subroutine multiDestructor(self)
    !!{
    Destructor for the \refClass{taskMulti} task class.
    !!}
    implicit none
    type(taskMulti    ), intent(inout) :: self
    type(multiTaskList), pointer       :: task_, taskNext

    if (associated(self%tasks)) then
       task_ => self%tasks
       do while (associated(task_))
          taskNext => task_%next
          !![
          <objectDestructor name="task_%task_"/>
          !!]
          deallocate(task_)
          task_ => taskNext
       end do
    end if
    return
  end subroutine multiDestructor

  subroutine multiPerform(self,status)
    !!{
    Perform all tasks.
    !!}
    use :: Display, only : displayIndent     , displayUnindent
    use :: Error  , only : errorStatusSuccess
    implicit none
    class  (taskMulti    ), intent(inout), target   :: self
    integer               , intent(  out), optional :: status
    type   (multiTaskList), pointer                 :: task_

    call displayIndent('Begin multiple tasks')
    if (present(status)) status=errorStatusSuccess
    task_ => self%tasks
    do while (associated(task_))
       call task_%task_%perform(status)
       if (present(status) .and. status /= errorStatusSuccess) return
       task_ => task_%next
    end do
    call displayUnindent('Done multiple tasks')
    return
  end subroutine multiPerform

  logical function multiRequiresOutputFile(self)
    !!{
    Returns true if any sub-task requires that the output file be open.
    !!}
    implicit none
    class(taskMulti    ), intent(inout) :: self
    type (multiTaskList), pointer       :: task_

    multiRequiresOutputFile =  .false.
    task_                   => self%tasks
    do while (associated(task_).and..not.multiRequiresOutputFile)
       multiRequiresOutputFile=task_%task_%requiresOutputFile()
       task_ => task_%next
    end do
    return
  end function multiRequiresOutputFile
