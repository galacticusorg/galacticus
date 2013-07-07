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

  logical function Node_Branch_Jump(thisEvent,thisNode,deadlockStatus)
    !% Moves a satellite node to a different branch of the merger tree.
    use Input_Parameters
    use Galacticus_Nodes
    use Galacticus_Error
    use Galacticus_Display
    use Merger_Trees_Evolve_Deadlock_Status
    use ISO_Varying_String
    use String_Handling
    implicit none
    class  (nodeEvent     ), intent(in   )          :: thisEvent
    type   (treeNode      ), intent(inout), pointer :: thisNode
    integer                , intent(inout)          :: deadlockStatus
    type   (treeNode      )               , pointer :: lastSatellite , newHost
    type   (varying_string)                         :: message

    ! If the node is not yet a satellite, wait until it is before peforming this task.
    if (.not.thisNode%isSatellite()) then
       Node_Branch_Jump=.false.
       return
    else
       Node_Branch_Jump=.true.
    end if
    ! Report.
    message='Node ['
    message=message//thisNode%index()//'] jumping branch to ['//thisEvent%node%index()//']'
    call Galacticus_Display_Message(message,verbosityInfo)
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
    ! Update the host tree pointer in the node to point to its new tree.
    thisNode%hostTree => newHost%hostTree
    ! Locate the paired event in the host and remove it.
    call newHost%removePairedEvent(thisEvent)
    ! Since we changed the tree, record that the tree is not deadlocked.
    deadlockStatus=isNotDeadlocked
    return
  end function Node_Branch_Jump

end module Node_Branch_Jumps
