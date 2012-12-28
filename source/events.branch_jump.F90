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

!% Contains a module which handles node branch jump events.

module Node_Branch_Jumps
  !% Handles satellite node branch jump events.
  implicit none
  private
  public :: Node_Branch_Jump
  
contains

  subroutine Node_Branch_Jump(thisEvent,thisNode,deadlockStatus)
    !% Moves a satellite node to a different branch of the merger tree.
    use Input_Parameters
    use Galacticus_Nodes
    use Galacticus_Error
    use Merger_Trees
    use Merger_Trees_Evolve_Deadlock_Status
    implicit none
    class  (nodeEvent ),          intent(in   ) :: thisEvent
    type   (treeNode  ), pointer, intent(inout) :: thisNode
    integer            ,          intent(inout) :: deadlockStatus
    type   (nodeEvent ), pointer                :: pairEvent,lastEvent,nextEvent
    type   (treeNode  ), pointer                :: newHost,lastSatellite
    type   (mergerTree)                         :: thisTree
    logical                                     :: pairMatched

    ! If the node is not yet a satellite, it must become one.
    if (.not.thisNode%isSatellite()) call thisTree%mergeNode(thisNode)
    ! Remove the satellite from its current host.
    call thisNode%removeFromHost()
    ! Find the new host and insert the node as a satellite in that new host.
    newHost          => thisEvent%node
    thisNode%sibling => null()
    thisNode%parent  => newHost
    if (associated(newHost%firstSatellite)) then
       lastSatellite                => newHost %lastSatellite()
       lastSatellite%sibling        => thisNode
    else
       newHost      %firstSatellite => thisNode
    end if
    ! Locate the paired event in the host and remove it.
    pairEvent => newHost%event
    lastEvent => newHost%event
    ! Iterate over all events.
    pairMatched=.false.
    do while (associated(pairEvent).and..not.pairMatched)
       ! Match the paired event ID with the current event ID.
       if (pairEvent%ID == thisEvent%ID) then
          pairMatched=.true.
          if (associated(pairEvent,newHost%event)) then
             newHost  %event => pairEvent%next
             lastEvent       => newHost %event
          else
             lastEvent%next  => pairEvent%next
          end if
          nextEvent => pairEvent%next
          deallocate(pairEvent)
          pairEvent => nextEvent
       else
          lastEvent => pairEvent
          pairEvent => pairEvent%next
       end if
    end do
    if (.not.pairMatched) call Galacticus_Error_Report('Node_Branch_Jump','unable to find paired event')
    ! Since we changed the tree, record that the tree is not deadlocked.
    deadlockStatus=isNotDeadlocked
    return
  end subroutine Node_Branch_Jump
  
end module Node_Branch_Jumps
