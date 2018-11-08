!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains module which handles hooking of objects into events.

module Events_Hooks
  !% Handles hooking of object function class into events.
  private
  public :: hook

  type :: hook
     !% Class for individual hooked function calls. Stores the function to be called, and the object to be passed as its first
     !% argument.
     class    (*   ), pointer         :: object_   => null()
     procedure(    ), pointer, nopass :: function_ => null()
     class    (hook), pointer         :: next
  end type hook

  type :: eventHook
     !% Class used to define a set of hooked function calls for a given event.
     private
     integer                :: count_ =  0
     type   (hook), pointer :: first_  => null()
   contains
     !@ <objectMethods>
     !@   <object>eventHook</object>
     !@   <objectMethod>
     !@     <method>attach</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless class(*)\textgreater} *object\_\argin, \textcolor{red}{\textless external\textgreater} *function\_\argin</arguments>
     !@     <description>Attach a hook to the event.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>count</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Return a count of the number of hooks into this event.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>first</method>
     !@     <type></type>
     !@     <arguments>\textcolor{red}{\textless *type(hook)\textgreater}</arguments>
     !@     <description>Return a pointer to the first hook into this event.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: attach => eventHookAttach
     procedure :: count  => eventHookCount
     procedure :: first  => eventHookFirst
  end type eventHook
  
  !# <eventHookManager/>
  
contains

  subroutine eventHookAttach(self,object_,function_)
    !% Attach an object to an event hook.
    implicit none
    class   (eventHook), intent(inout)          :: self
    class   (*        ), intent(in   ), target  :: object_
    external                                    :: function_
    type    (hook     )               , pointer :: hook_

    ! Allocate the next entry in our list of hooks.
    if (associated(self%first_)) then
       hook_ => self%first_
       do while (associated(hook_%next))
          hook_ => hook_%next
       end do
       allocate(hook_%next)
       hook_ => hook_%next
    else
       allocate(self%first_)
       hook_ => self%first_
    end if
    ! Create the new hook.
    hook_%object_   => object_
    hook_%function_ => function_
    ! Increment the count of hooks into this event.
    self%count_=self%count_+1
    return
  end subroutine eventHookAttach

  integer function eventHookCount(self)
    !% Return a count of the number of hooks into this event.
    implicit none
    class(eventHook), intent(in   ):: self

    eventHookCount=self%count_
    return
  end function eventHookCount

  function eventHookFirst(self)
    !% Return a pointer to the first hook into this event.
    implicit none
    type (hook     ), pointer      :: eventHookFirst
    class(eventHook), intent(in   ):: self

    eventHookFirst => self%first_
    return
  end function eventHookFirst
  
end module Events_Hooks
