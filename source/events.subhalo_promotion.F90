!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which handles node subhalo promotion events.

module Node_Subhalo_Promotions
  !% Handles subhalo promotion events.
  implicit none
  private
  public :: Node_Subhalo_Promotion

contains

  logical function Node_Subhalo_Promotion(thisEvent,thisNode,deadlockStatus)
    !% Promotes a subhalo to be an isolated node.
    use Galacticus_Nodes
    use Merger_Trees_Evolve_Deadlock_Status
    use ISO_Varying_String
    use String_Handling
    use Galacticus_Display
    use Merger_Trees_Evolve_Node
    implicit none
    class    (nodeEvent                     ), intent(in   )          :: thisEvent
    type     (treeNode                      ), intent(inout), pointer :: thisNode
    integer                                  , intent(inout)          :: deadlockStatus
    type     (treeNode                      )               , pointer :: promotionNode
    class    (nodeComponentBasic            )               , pointer :: parentBasic
    class    (nodeComponentMergingStatistics)               , pointer :: mergingStatistics
    type     (varying_string                )                         :: message
    character(len=12                        )                         :: label
    
    ! Find the node to promote to.
    promotionNode => thisEvent%node
    ! If the target node has a child, we must wait for that child to be processed before promoting. Note that this should only
    ! happen in cases where the target node was cloned to be its own primary progenitor.
    if (associated(promotionNode%firstChild)) then
       Node_Subhalo_Promotion=.false.
       return
    end if
    ! Report.
    if (Galacticus_Verbosity_Level() >= verbosityInfo) then
       write (label,'(f12.6)') thisEvent%time
       message='Satellite node ['
       message=message//thisNode%index()//'] promoting to isolated node ['//thisEvent%node%index()//'] at time '//trim(label)//' Gyr'
       call Galacticus_Display_Message(message)
    end if
    ! Remove the subhalo from its host.
    call thisNode%removeFromHost()
    ! Make thisNode the primary progenitor of the target node.
    thisNode%parent          => promotionNode
    thisNode%sibling         => null()
    promotionNode%firstChild => thisNode
    ! Reset the mass-when-first-isolated property of the merging statistics component if possible.
    mergingStatistics => thisNode%mergingStatistics()
    if (mergingStatistics%massWhenFirstIsolatedIsSettable()) then
       parentBasic => promotionNode%basic()
       call mergingStatistics%massWhenFirstIsolatedSet(parentBasic%mass())
    end if
    ! Promote the halo.
    call Tree_Node_Promote(thisNode)
    ! Since we changed the tree, record that the tree is not deadlocked.
    deadlockStatus=deadlockStatusIsNotDeadlocked
    ! Record that the task was performed.
    Node_Subhalo_Promotion=.true.
    return
  end function Node_Subhalo_Promotion

end module Node_Subhalo_Promotions
