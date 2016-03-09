!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Defines functions bound to the merger tree object class.

subroutine Merger_Tree_Destroy(thisTree)
  !% Destroys the entire merger tree.
  implicit none
  class(mergerTree), intent(inout) :: thisTree
  type (treeEvent ), pointer       :: thisEvent, nextEvent 

  select type (thisTree)
  type is (mergerTree)
     ! Destroy all nodes.
     if (associated(thisTree%baseNode)) then
        call thisTree%baseNode%destroyBranch()
        deallocate(thisTree%baseNode)
     end if
     ! Destroy the HDF5 group associated with this tree.
     call thisTree%hdf5Group%destroy()
     ! Destroy any events attached to this tree.
     thisEvent => thisTree%event
     do while (associated(thisEvent))
        nextEvent => thisEvent%next
        deallocate(thisEvent)
        thisEvent => nextEvent
     end do
  end select
  return
end subroutine Merger_Tree_Destroy

function Merger_Tree_Node_Get(thisTree,nodeIndex)
  !% Return a pointer to a node in {\normalfont \ttfamily thisTree} given the index of the node.
  implicit none
  class  (mergerTree    ), intent(in   ), target :: thisTree
  integer(kind=kind_int8), intent(in   )         :: nodeIndex
  type   (mergerTree    ), pointer               :: currentTree
  type   (treeNode      ), pointer               :: Merger_Tree_Node_Get, thisNode

  Merger_Tree_Node_Get => null()
  currentTree => thisTree
  do while (associated(currentTree))
     thisNode => currentTree%baseNode
     do while (associated(thisNode))
        if (thisNode%index() == nodeIndex) then
           Merger_Tree_Node_Get => thisNode
           return
        end if
        thisNode => thisNode%walkTreeWithSatellites()
     end do
     currentTree => currentTree%nextTree
  end do
  return
end function Merger_Tree_Node_Get

function Merger_Tree_Create_Event(self) result (newEvent)
  !% Create a new event in a merger tree.
  implicit none
  class(mergerTree), intent(inout) :: self
  type (treeEvent ), pointer       :: newEvent, thisEvent
  
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
end function Merger_Tree_Create_Event

subroutine Merger_Tree_Remove_Event(self,event)
  !% Removed an event from {\normalfont \ttfamily self}.
  implicit none
  class  (mergerTree), intent(inout) :: self
  type   (treeEvent ), intent(in   ) :: event
  type   (treeEvent ), pointer       :: lastEvent, nextEvent, thisEvent
  
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
end subroutine Merger_Tree_Remove_Event

double precision function Merger_Tree_Earliest_Time(self)
  !% Return the earliest time in a merger tree.
  implicit none
  class           (mergerTree        ), intent(inout), target :: self
  double precision                    , parameter             :: timeInfinity=1.0d30
  type            (mergerTree        ), pointer               :: currentTree
  type            (treeNode          ), pointer               :: thisNode
  class           (nodeComponentBasic), pointer               :: thisBasic
  
  Merger_Tree_Earliest_Time =  timeInfinity
  currentTree               => self
  do while (associated(currentTree))
     thisNode => currentTree%baseNode
     do while (associated(thisNode))
        if (.not.associated(thisNode%firstChild)) then
           thisBasic                 =>                               thisNode %basic()
           Merger_Tree_Earliest_Time =  min(Merger_Tree_Earliest_Time,thisBasic%time ())
        end if
        thisNode => thisNode%walkTreeWithSatellites()
     end do
     currentTree => currentTree%nextTree
  end do
  return
end function Merger_Tree_Earliest_Time

double precision function Merger_Tree_Earliest_Time_Evolving(self)
  !% Return the earliest time in a merger tree.
  implicit none
  class           (mergerTree        ), intent(inout), target :: self
  double precision                    , parameter             :: timeInfinity=1.0d30
  type            (mergerTree        ), pointer               :: currentTree
  type            (treeNode          ), pointer               :: thisNode
  class           (nodeComponentBasic), pointer               :: thisBasic
  
  Merger_Tree_Earliest_Time_Evolving =  timeInfinity
  currentTree                        => self
  do while (associated(currentTree))
     thisNode => currentTree%baseNode
     do while (associated(thisNode))
        if (.not.associated(thisNode%firstChild).and.(associated(thisNode%parent).or..not.associated(thisNode%firstSatellite))) then
           thisBasic                          =>                                        thisNode %basic()
           Merger_Tree_Earliest_Time_Evolving =  min(Merger_Tree_Earliest_Time_Evolving,thisBasic%time ())
        end if
        thisNode => thisNode%walkTreeWithSatellites()
     end do
     currentTree => currentTree%nextTree
  end do
  return
end function Merger_Tree_Earliest_Time_Evolving

double precision function Merger_Tree_Latest_Time(self)
  !% Return the latest time in a merger tree.
  implicit none
  class           (mergerTree        ), intent(inout), target :: self
  type            (mergerTree        ), pointer               :: currentTree
  class           (nodeComponentBasic), pointer               :: baseBasic
  
  Merger_Tree_Latest_Time =  -1.0d0
  currentTree             => self
  do while (associated(currentTree))
     if (associated(currentTree%baseNode)) then
        baseBasic               =>                             currentTree%baseNode%basic()
        Merger_Tree_Latest_Time =  max(Merger_Tree_Latest_Time,baseBasic  %time          ())
     end if
     currentTree => currentTree%nextTree
  end do
  return
end function Merger_Tree_Latest_Time
