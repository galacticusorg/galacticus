!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which defines and keeps track of the current task in {\sc Galacticus}.

module Galacticus_Tasks
  !% Defines and keeps track of the current task in {\sc Galacticus}.
  implicit none
  private
  public :: Galacticus_Task_Do

  integer, parameter :: taskBegin   =0         !  Initial value to indicate the start.     
  integer, parameter :: taskFinished=-1        !  Task that indicates all tasks are done.  
                                                                                        
  integer            :: currentTask =taskBegin !  Variable that tracks the current task.   
                                                                                        
contains

  subroutine Galacticus_Task_Do()
    !% Performs \glc\ tasks.
    !# <include directive="galacticusTask" type="moduleUse">
    include 'galacticus.tasks.task_rules.modules.inc'
    !# </include>
    implicit none
    logical                     :: tasksRemaining=.true.   
    procedure(logical), pointer :: taskFunction  =>null()  
                                                        
    do while (tasksRemaining)
       tasksRemaining=.false.
       !# <include directive="galacticusTask" type="functionCall" functionType="pointer">
       !#  <pointerName>taskFunction</pointerName>
       !#  <onReturn>tasksRemaining=tasksRemaining.or.taskFunction()</onReturn>
       include 'galacticus.tasks.task_rules.inc'
       !# </include>
    end do
    return
  end subroutine Galacticus_Task_Do

end module Galacticus_Tasks
