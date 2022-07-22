!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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

!!{
Contains a module which handles node subhalo promotion events.
!!}

module Node_Subhalo_Promotions
  !!{
  Handles subhalo promotion events.
  !!}
  implicit none
  private
  public :: nodeSubhaloPromotionPerform

contains

  logical function nodeSubhaloPromotionPerform(event,node,deadlockStatus)
    !!{
    Promotes a subhalo to be an isolated node.
    !!}
    use :: Display                            , only : displayMessage               , displayVerbosity             , verbosityLevelInfo
    use :: Galacticus_Nodes                   , only : nodeComponentBasic           , nodeEvent                    , treeNode
    use :: ISO_Varying_String                 , only : assignment(=)                , operator(//)                 , varying_string
    use :: Merger_Trees_Evolve_Deadlock_Status, only : deadlockStatusIsNotDeadlocked, enumerationDeadlockStatusType
    use :: String_Handling                    , only : operator(//)
    implicit none
    class    (nodeEvent                    ), intent(in   )          :: event
    type     (treeNode                     ), intent(inout), pointer :: node
    type     (enumerationDeadlockStatusType), intent(inout)          :: deadlockStatus
    type     (treeNode                     )               , pointer :: nodePromotion
    type     (varying_string               )                         :: message
    character(len=12                       )                         :: label

    ! Find the node to promote to.
    nodePromotion => event%node
    ! If the target node has a child, we must wait for that child to be processed before promoting. Note that this should only
    ! happen in cases where the target node was cloned to be its own primary progenitor.
    if (associated(nodePromotion%firstChild)) then
       nodeSubhaloPromotionPerform=.false.
       return
    end if
    ! Report.
    if (displayVerbosity() >= verbosityLevelInfo) then
       write (label,'(f12.6)') event%time
       message='Satellite node ['
       message=message//node%index()//'] promoting to isolated node ['//event%node%index()//'] at time '//trim(label)//' Gyr'
       call displayMessage(message)
    end if
    ! Remove the subhalo from its host.
    call node%removeFromHost  ()
    call node%removeFromMergee()
    ! Make node the primary progenitor of the target node.
    node         %parent     => nodePromotion
    node         %sibling    => null()
    nodePromotion%firstChild => node
    ! Trigger the event.
    !![
    <eventHook name="subhaloPromotion">
     <import>
      <module name="Galacticus_Nodes" symbols="treeNode"/>
     </import>
     <interface>
      type(treeNode), intent(inout), pointer :: node, nodePromotion
     </interface>
     <callWith>node,nodePromotion</callWith>
    </eventHook>
    !!]
    ! Since we changed the tree, record that the tree is not deadlocked.
    deadlockStatus=deadlockStatusIsNotDeadlocked
    ! Record that the task was performed.
    nodeSubhaloPromotionPerform=.true.
    return
  end function nodeSubhaloPromotionPerform

end module Node_Subhalo_Promotions
