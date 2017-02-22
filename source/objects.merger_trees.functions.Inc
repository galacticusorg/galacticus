!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Defines functions bound to the merger tree object class.

subroutine Merger_Tree_Destroy(tree)
  !% Destroys the entire merger tree.
  implicit none
  class(mergerTree), intent(inout) :: tree
  type (treeEvent ), pointer       :: event, eventNext

  select type (tree)
  type is (mergerTree)
     ! Destroy all nodes.
     if (associated(tree%baseNode)) then
        call tree%baseNode%destroyBranch()
        deallocate(tree%baseNode)
     end if
     ! Destroy the HDF5 group associated with this tree.
     call tree%hdf5Group%destroy()
     ! Destroy any events attached to this tree.
     event => tree%event
     do while (associated(event))
        eventNext => event%next
        deallocate(event)
        event => eventNext
     end do
  end select
  return
end subroutine Merger_Tree_Destroy

function Merger_Tree_Node_Get(tree,nodeIndex)
  !% Return a pointer to a node in {\normalfont \ttfamily tree} given the index of the node.
  implicit none
  class  (mergerTree    ), intent(in   ), target :: tree
  integer(kind=kind_int8), intent(in   )         :: nodeIndex
  type   (mergerTree    ), pointer               :: treeCurrent
  type   (treeNode      ), pointer               :: Merger_Tree_Node_Get, node

  Merger_Tree_Node_Get => null()
  treeCurrent => tree
  do while (associated(treeCurrent))
     node => treeCurrent%baseNode
     do while (associated(node))
        if (node%index() == nodeIndex) then
           Merger_Tree_Node_Get => node
           return
        end if
        node => node%walkTreeWithSatellites()
     end do
     treeCurrent => treeCurrent%nextTree
  end do
  return
end function Merger_Tree_Node_Get

function Merger_Tree_Create_Event(self) result (eventNew)
  !% Create a new event in a merger tree.
  implicit none
  class(mergerTree), intent(inout) :: self
  type (treeEvent ), pointer       :: eventNew, event
  
  allocate(eventNew)
  nullify(eventNew%next)
  !$omp atomic
  eventID=eventID+1
  eventNew%ID=eventID
  if (associated(self%event)) then
     event => self%event
     do while (associated(event%next))
        event => event%next
     end do
     event%next => eventNew
  else
     self%event => eventNew
  end if
  return
end function Merger_Tree_Create_Event

subroutine Merger_Tree_Remove_Event(self,event)
  !% Removed an event from {\normalfont \ttfamily self}.
  implicit none
  class  (mergerTree), intent(inout) :: self
  type   (treeEvent ), intent(in   ) :: event
  type   (treeEvent ), pointer       :: eventLast, eventNext, eventCurrent
  
  ! Destroy the event.
  eventLast    => null()
  eventCurrent => self%event
  do while (associated(eventCurrent))
     ! Match the event.
     if (eventCurrent%ID == event%ID) then
        if (associated(eventCurrent,self%event)) then
           self     %event => eventCurrent%next
           eventLast       => self %event
        else
           eventLast%next  => eventCurrent%next
        end if
        eventNext => eventCurrent%next
        deallocate(eventCurrent)
        eventCurrent => eventNext
     else
        eventLast    => eventCurrent
        eventCurrent => eventCurrent%next
     end if
  end do
  return
end subroutine Merger_Tree_Remove_Event

double precision function Merger_Tree_Earliest_Time(self)
  !% Return the earliest time in a merger tree.
  implicit none
  class           (mergerTree        ), intent(inout), target :: self
  double precision                    , parameter             :: timeInfinity=huge(1.0d0)
  type            (mergerTree        ), pointer               :: treeCurrent
  type            (treeNode          ), pointer               :: node
  class           (nodeComponentBasic), pointer               :: basic
  
  Merger_Tree_Earliest_Time =  timeInfinity
  treeCurrent               => self
  do while (associated(treeCurrent))
     node => treeCurrent%baseNode
     do while (associated(node))
        if (.not.associated(node%firstChild)) then
           basic                     =>                               node %basic()
           Merger_Tree_Earliest_Time =  min(Merger_Tree_Earliest_Time,basic%time ())
        end if
        node => node%walkTreeWithSatellites()
     end do
     treeCurrent => treeCurrent%nextTree
  end do
  return
end function Merger_Tree_Earliest_Time

double precision function Merger_Tree_Earliest_Time_Evolving(self)
  !% Return the earliest time in a merger tree.
  implicit none
  class           (mergerTree        ), intent(inout), target :: self
  double precision                    , parameter             :: timeInfinity=huge(1.0d0)
  type            (mergerTree        ), pointer               :: treeCurrent
  type            (treeNode          ), pointer               :: node
  class           (nodeComponentBasic), pointer               :: basic
  
  Merger_Tree_Earliest_Time_Evolving =  timeInfinity
  treeCurrent                        => self
  do while (associated(treeCurrent))
     node => treeCurrent%baseNode
     do while (associated(node))
        if (.not.associated(node%firstChild).and.(associated(node%parent).or..not.associated(node%firstSatellite))) then
           basic                              =>                                        node %basic()
           Merger_Tree_Earliest_Time_Evolving =  min(Merger_Tree_Earliest_Time_Evolving,basic%time ())
        end if
        node => node%walkTreeWithSatellites()
     end do
     treeCurrent => treeCurrent%nextTree
  end do
  return
end function Merger_Tree_Earliest_Time_Evolving

double precision function Merger_Tree_Latest_Time(self)
  !% Return the latest time in a merger tree.
  implicit none
  class           (mergerTree        ), intent(inout), target :: self
  type            (mergerTree        ), pointer               :: treeCurrent
  class           (nodeComponentBasic), pointer               :: baseBasic
  
  Merger_Tree_Latest_Time =  -1.0d0
  treeCurrent             => self
  do while (associated(treeCurrent))
     if (associated(treeCurrent%baseNode)) then
        baseBasic               =>                             treeCurrent%baseNode%basic()
        Merger_Tree_Latest_Time =  max(Merger_Tree_Latest_Time,baseBasic  %time          ())
     end if
     treeCurrent => treeCurrent%nextTree
  end do
  return
end function Merger_Tree_Latest_Time
