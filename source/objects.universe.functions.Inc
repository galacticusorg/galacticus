!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Defines functions bound to the universe class.

function universeCreateEvent(self) result (newEvent)
  !% Create a new event in a universe.
  implicit none
  class(universe     ), intent(inout) :: self
  type (universeEvent), pointer       :: newEvent, thisEvent
  
  allocate(newEvent)
  nullify(newEvent%next)
  !$omp atomic
  eventID=eventID+1
  newEvent%ID=eventID
  if (associated(self%event)) then
     thisEvent => self%event
     do while (associated(thisEvent%next))
        thisEvent => thisEvent%next
     end do
     thisEvent%next => newEvent
  else
     self%event => newEvent
  end if
  return
end function universeCreateEvent

subroutine universeRemoveEvent(self,event)
  !% Remove an event from {\tt self}.
  implicit none
  class  (universe     ), intent(inout) :: self
  type   (universeEvent), intent(in   ) :: event
  type   (universeEvent), pointer       :: lastEvent, nextEvent, thisEvent
  
  ! Destroy the event.
  thisEvent => self%event
  do while (associated(thisEvent))
     ! Match the event.
     if (thisEvent%ID == event%ID) then
        if (associated(thisEvent,self%event)) then
           self     %event => thisEvent%next
           lastEvent       => self     %event
        else
           lastEvent%next  => thisEvent%next
        end if
        nextEvent => thisEvent%next
        deallocate(thisEvent)
        thisEvent => nextEvent
     else
        lastEvent => thisEvent
        thisEvent => thisEvent%next
     end if
  end do
  return
end subroutine universeRemoveEvent

function universePopTree(self)
  !% Pop a tree from the universe.
  implicit none
  type (mergerTree    ), pointer       :: universePopTree
  class(universe      ), intent(inout) :: self
  type (mergerTreeList), pointer       :: nextTree

  if (associated(self%trees)) then
     universePopTree => self%trees%tree
     nextTree        => self%trees%next
     deallocate(self%trees)
     self%trees      => nextTree
  else
     universePopTree => null()
  end if
  return
end function universePopTree

subroutine universePushTree(self,thisTree)
  !% Pop a tree from the universe.
  implicit none
  class(universe      ), intent(inout)          :: self
  type (mergerTree    ), intent(in   ), pointer :: thisTree
  type (mergerTreeList)               , pointer :: nextTree, newTree

  allocate(newTree)
  newTree%tree => thisTree
  newTree%next => null()
  if (associated(self%trees)) then
     nextTree => self%trees
     do while (associated(nextTree%next))
        nextTree => nextTree%next
     end do
     nextTree%next => newTree
  else
     self%trees => newTree
  end if
  return
end subroutine universePushTree
