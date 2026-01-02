!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  Implements a merger tree operator which prunes all branches which do not contain an ``essential''
  node.
  !!}

  use :: Kind_Numbers, only : kind_int8

  !![
  <mergerTreeOperator name="mergerTreeOperatorPruneNonEssential">
   <description>
    A merger tree operator class which prunes branches that do not directly influence an ``essential'' node. Any branch which
    does not connect to the branch into which the node identified by ID {\normalfont \ttfamily [essentialNodeID]} descends by
    time {\normalfont \ttfamily essentialNodeTime]} will be pruned. Specifying the time is important---if the node is a
    satellite at this time, then the pruning will not remove any progenitors of the parent node in which the essential node
    lives at the specified time.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneNonEssential
     !!{
     A merger tree operator class which prunes branches which do not contain an ``essential'' node.
     !!}
     private
     integer         (kind=kind_int8) :: essentialNodeID
     double precision                 :: essentialNodeTime
   contains
     procedure :: operatePreEvolution => pruneNonEssentialOperatePreEvolution
  end type mergerTreeOperatorPruneNonEssential

  interface mergerTreeOperatorPruneNonEssential
     !!{
     Constructors for the prune-non-essential merger tree operator class.
     !!}
     module procedure pruneNonEssentialConstructorParameters
     module procedure pruneNonEssentialConstructorInternal
  end interface mergerTreeOperatorPruneNonEssential

contains

  function pruneNonEssentialConstructorParameters(parameters) result(self)
    !!{
    Constructor for the prune-non-essential merger tree operator class which takes a parameter set as input.
    !!}
    implicit none
    type            (mergerTreeOperatorPruneNonEssential)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    integer                                                              :: essentialNodeID
    double precision                                                     :: essentialNodeTime

    !![
    <inputParameter>
      <name>essentialNodeID</name>
      <source>parameters</source>
      <description>ID of the essential node to avoid pruning.</description>
    </inputParameter>
    <inputParameter>
      <name>essentialNodeTime</name>
      <source>parameters</source>
      <description>Time of the essential node to avoid pruning.</description>
    </inputParameter>
    !!]
    self=mergerTreeOperatorPruneNonEssential(essentialNodeID,essentialNodeTime)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function pruneNonEssentialConstructorParameters

  function pruneNonEssentialConstructorInternal(essentialNodeID,essentialNodeTime) result(self)
    !!{
    Internal constructor for the prune-non-essential merger tree operator class.
    !!}
    implicit none
    type            (mergerTreeOperatorPruneNonEssential)                :: self
    integer                                              , intent(in   ) :: essentialNodeID
    double precision                                     , intent(in   ) :: essentialNodeTime
    !![
    <constructorAssign variables="essentialNodeID, essentialNodeTime"/>
    !!]
    
   return
  end function pruneNonEssentialConstructorInternal

  subroutine pruneNonEssentialOperatePreEvolution(self,tree)
    !!{
    Perform a prune-non-essential operation on a merger tree.
    !!}
    use :: Galacticus_Nodes              , only : mergerTree                    , nodeComponentBasic             , treeNode
    use :: Merger_Tree_Walkers           , only : mergerTreeWalkerIsolatedNodes
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Clean_Branch, Merger_Tree_Prune_Uniqueify_IDs, Merger_Tree_Prune_Unlink_Parent
    implicit none
    class  (mergerTreeOperatorPruneNonEssential), intent(inout), target  :: self
    type   (mergerTree                         ), intent(inout), target  :: tree
    type   (mergerTree                         )               , pointer :: treeCurrent
    type   (treeNode                           )               , pointer :: node          , nodeEssential, &
         &                                                                  nodePrevious
    class  (nodeComponentBasic                 )               , pointer :: basic
    type   (mergerTreeWalkerIsolatedNodes      )                         :: treeWalker

    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Find the essential node.
       nodeEssential => treeCurrent%getNode(self%essentialNodeID)
       if (associated(nodeEssential)) then
          ! Trace the essential node to the required time.
          basic => nodeEssential%basic()
          do while (basic%time() < self%essentialNodeTime .and. associated(nodeEssential%parent))
             nodeEssential => nodeEssential%parent
             basic         => nodeEssential%basic ()
          end do
          ! Walk the tree, pruning branches.
          treeWalker=mergerTreeWalkerIsolatedNodes(treeCurrent)
          do while (treeWalker%next(node))
             ! Record the parent node to which we will return.
             if     (                                                &
                  &   associated(node%parent)                        &
                  &  .and.                                           &
                  &   .not.                                          &
                  &    (                                             &
                  &      node         %isProgenitorOf(nodeEssential) &
                  &     .or.                                         &
                  &      nodeEssential%isProgenitorOf(node         ) &
                  &    )                                             &
                  & ) then
                ! Return to the previous node.
                call treeWalker%previous(nodePrevious)
                ! Decouple from other nodes.
                call Merger_Tree_Prune_Unlink_Parent(                                                                  &
                     &                               node                                                            , &
                     &                               node%parent                                                     , &
                     &                               .not.                                                             &
                     &                                    (                                                            &
                     &                                      node         %parent%isProgenitorOf(nodeEssential       )  &
                     &                                     .or.                                                        &
                     &                                      nodeEssential       %isProgenitorOf(node         %parent)  &
                     &                                    )                                                          , &
                     &                               .true.                                                            &
                     &                              )
                ! Clean the branch.
                call Merger_Tree_Prune_Clean_Branch(node)
                ! Destroy the branch.
                call node%destroyBranch()
                deallocate(node)
             end if
          end do
       else
          ! Entire tree can be pruned. Destroy all but this base node. (Leaving just
          ! the base node makes the tree inert - i.e. it can not do anything.)
          node => treeCurrent%nodeBase%firstChild
          do while (associated(node))
             nodePrevious => node%sibling
             call Merger_Tree_Prune_Clean_Branch(node)
             call node%destroyBranch()
             deallocate(node)
             node => nodePrevious
          end do
       end if
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    ! Uniqueify nodes.
    call Merger_Tree_Prune_Uniqueify_IDs(tree)
    return
  end subroutine pruneNonEssentialOperatePreEvolution
