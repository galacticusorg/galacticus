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

!% Contains a module which handles node subhalo promotion events.

module Node_Subhalo_Promotions
  !% Handles subhalo promotion events.
  implicit none
  private
  public :: Node_Subhalo_Promotion

contains

  logical function Node_Subhalo_Promotion(thisEvent,thisNode,deadlockStatus)
    !% Promotes a subhalo to be an isolated node.
    use Input_Parameters
    use Galacticus_Nodes
    use Galacticus_Error
    use Merger_Trees_Evolve_Deadlock_Status
    use ISO_Varying_String
    use String_Handling
    use Galacticus_Display
    use Merger_Trees_Evolve_Node
    implicit none
    class  (nodeEvent     ), intent(in   )          :: thisEvent
    type   (treeNode      ), intent(inout), pointer :: thisNode
    integer                , intent(inout)          :: deadlockStatus
    type   (treeNode      )               , pointer :: promotionNode
    type   (mergerTree    )                         :: thisTree
    type   (varying_string)                         :: message

    ! Find the node to promote to.
    promotionNode => thisEvent%node
    ! If the target node has a child, we must wait for that child to be processed before promoting. Note that this should only
    ! happen in cases where the target node was cloned to be its own primary progenitor.
    if (associated(promotionNode%firstChild)) then
       Node_Subhalo_Promotion=.false.
       return
    end if
    ! Report.
    message='Satellite node ['
    message=message//thisNode%index()//'] promoting to isolated node ['//thisEvent%node%index()//']'
    call Galacticus_Display_Message(message,verbosityInfo)
    ! Remove the subhalo from its host.
    call thisNode%removeFromHost()
    ! Make thisNode the primary progenitor of the target node.
    thisNode%parent          => promotionNode
    thisNode%sibling         => null()
    promotionNode%firstChild => thisNode
    ! Promote the halo.
    call Tree_Node_Promote(thisTree,thisNode)
    ! Since we changed the tree, record that the tree is not deadlocked.
    deadlockStatus=isNotDeadlocked
    ! Record that the task was performed.
    Node_Subhalo_Promotion=.true.
    return
  end function Node_Subhalo_Promotion

end module Node_Subhalo_Promotions
