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

!% Contains a module which implements semaphores.

! Specify an explicit dependence on the semaphores.o object file.
!: ./work/build/semaphores.o

module Semaphores
  !% Implements semaphores.
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: Semaphore_Open, Semaphore_Post_On_Error

  type, public :: semaphore
     type     (c_ptr  ) :: s
     character(len=251) :: name
     integer            :: waitCount
   contains
     procedure :: close  => Semaphore_Close
     procedure :: wait   => Semaphore_Wait
     procedure :: post   => Semaphore_Post
     procedure :: unlink => Semaphore_Unlink
     !@ <objectMethods>
     !@   <object>semaphore</object>
     !@   <objectMethod>
     !@     <method>close</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Close the semaphore.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>wait</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Wait for a semaphore to become available.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>post</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Post to (i.e. release) a semaphore.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>unlink</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Unlink a semaphore.</description>
     !@   </objectMethod>
     !@ </objectMethods>
  end type semaphore

  type :: semaphoreList
     type(semaphore    ), pointer :: self => null()
     type(semaphoreList), pointer :: next => null()
  end type semaphoreList

  type(semaphoreList), pointer :: semaphoreListHead => null()

  interface
     function Semaphore_Open_C(name,initialValue) bind(c,name='Semaphore_Open_C')
       !% Template for a C function that calls {\tt sem\_open()}.
       import
       type     (c_ptr )        :: Semaphore_Open_C
       character(c_char)        :: name
       integer  (c_int ), value :: initialValue       
     end function Semaphore_Open_C
  end interface

  interface
     subroutine Semaphore_Close_C(s) bind(c,name='Semaphore_Close_C')
       !% Template for a C function that calls {\tt sem\_close()}.
       import
       type(c_ptr), value :: s
     end subroutine Semaphore_Close_C
  end interface

  interface
     subroutine Semaphore_Wait_C(s) bind(c,name='Semaphore_Wait_C')
       !% Template for a C function that calls {\tt sem\_wait()}.
       import
       type(c_ptr), value :: s
     end subroutine Semaphore_Wait_C
  end interface

  interface
     subroutine Semaphore_Post_C(s) bind(c,name='Semaphore_Post_C')
       !% Template for a C function that calls {\tt sem\_post()}.
       import
       type(c_ptr), value :: s
     end subroutine Semaphore_Post_C
  end interface

  interface
     subroutine Semaphore_Unlink_C(name) bind(c,name='Semaphore_Unlink_C')
       !% Template for a C function that calls {\tt sem\_unlink()}.
       import
       character(c_char) :: name
     end subroutine Semaphore_Unlink_C
  end interface

contains

  function Semaphore_Open(name,initialValue) result (newSemaphore)
    implicit none
    type     (semaphore    ), pointer       :: newSemaphore
    character(len=*        ), intent(in   ) :: name
    integer                 , intent(in   ) :: initialValue
    type     (semaphoreList), pointer       :: thisSemaphore
    
    allocate(newSemaphore)
    newSemaphore%s        =Semaphore_Open_C(trim(name)//char(0),initialValue)
    newSemaphore%name     =                      name
    newSemaphore%waitCount=0
    ! Maintain a linked list of semaphores.
    if (associated(semaphoreListHead%self)) then
       thisSemaphore => semaphoreListHead
       do while (associated(thisSemaphore%next))
          thisSemaphore => thisSemaphore%next
       end do
       allocate(thisSemaphore%next)
       thisSemaphore      => thisSemaphore%next
       thisSemaphore%self => newSemaphore
       thisSemaphore%next => null()
    else
       semaphoreListHead%self => newSemaphore
       semaphoreListHead%next => null()
    end if
    return
  end function Semaphore_Open

  subroutine Semaphore_Close(self)
    implicit none
    class(semaphore), intent(in   ) :: self

    call Semaphore_Close_C(self%s)
    return
  end subroutine Semaphore_Close

  subroutine Semaphore_Wait(self)
    implicit none
    class(semaphore), intent(inout) :: self

    call Semaphore_Wait_C(self%s)
    !$omp atomic
    self%waitCount=self%waitCount+1
    return
  end subroutine Semaphore_Wait

  subroutine Semaphore_Post(self)
    implicit none
    class(semaphore), intent(inout) :: self

    call Semaphore_Post_C(self%s)
    !$omp atomic
    self%waitCount=self%waitCount-1
    return
  end subroutine Semaphore_Post

  subroutine Semaphore_Unlink(self)
    implicit none
    class(semaphore), intent(in   ) :: self

    call Semaphore_Unlink_C(trim(self%name)//char(0))
    return
  end subroutine Semaphore_Unlink

  subroutine Semaphore_Post_On_Error()
    !% Attempts to post to all open semaphores before exiting the code in an error condition.
    implicit none
    type(semaphoreList), pointer :: thisSemaphore
    
    thisSemaphore => semaphoreListHead
    do while (associated(thisSemaphore))
       do while (thisSemaphore%self%waitCount > 0)
          call thisSemaphore%self%post()       
       end do
       thisSemaphore => thisSemaphore%next
    end do
    return
  end subroutine Semaphore_Post_On_Error

end module Semaphores
