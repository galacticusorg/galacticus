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
Contains a module which provides utility functions for pruning branches from merger trees.
!!}

module Merger_Trees_Pruning_Utilities
  implicit none
  private
  public :: Merger_Tree_Prune_Clean_Branch , Merger_Tree_Prune_Unlink_Parent, &
       &    Merger_Tree_Prune_Uniqueify_IDs

contains

  subroutine Merger_Tree_Prune_Clean_Branch(node)
    !!{
    Cleans pointers in a branch about to be pruned to avoid dangling pointer problems during tree evolution.
    !!}
    use :: Galacticus_Nodes   , only : treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodesBranch
    implicit none
    type(treeNode                      ), intent(inout) :: node
    type(treeNode                      ), pointer       :: workNode  , mergeeNode, &
         &                                                 nextMergee
    type(mergerTreeWalkerAllNodesBranch)                :: treeWalker

    ! Walk the branch to be pruned.
    treeWalker=mergerTreeWalkerAllNodesBranch(node)
    do while (treeWalker%next(workNode))
       ! Remove the current node from any merge target it may have.
       call workNode%removeFromMergee()
       ! If the current node has any mergees, unlink them from it.
       if (associated(workNode%firstMergee)) then
          mergeeNode => workNode%firstMergee
          do while (associated(mergeeNode))
             call mergeeNode%removeFromMergee()
             nullify(mergeeNode%mergeTarget)
             nextMergee => mergeeNode%siblingMergee
             nullify(mergeeNode%siblingMergee)
             mergeeNode => nextMergee
          end do
       end if
    end do
    return
  end subroutine Merger_Tree_Prune_Clean_Branch

  subroutine Merger_Tree_Prune_Unlink_Parent(node,parentNode,parentWillBePruned,preservePrimaryProgenitor)
    !!{
    Unlink a parent node from a tree branch which is about to be pruned.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    type   (treeNode          ), intent(inout), pointer :: node              , parentNode
    logical                    , intent(in   )          :: parentWillBePruned, preservePrimaryProgenitor
    type   (treeNode          ), pointer                :: newNode           , workNode
    class  (nodeComponentBasic), pointer                :: newBasic

    ! Check primary progenitor status.
    if (node%isPrimaryProgenitorOf(parentNode)) then
       ! Node is primary progenitor - we must check if the parent will be pruned also.
       if (parentWillBePruned.or..not.preservePrimaryProgenitor) then
          ! Parent will eventually be pruned - simply replace the current first child (about to pruned) with its
          ! sibling. Alternatively, we've been asked to not preserve primary progenitor status by inserting cloned nodes.
          parentNode%firstChild => node%sibling
       else
          ! Parent will not be pruned. Insert a clone of the parent as its own first progenitor to prevent any siblings of the
          ! to-be-pruned branch being misidentified as the primary progenitor.
          allocate(newNode)
          call parentNode%copyNodeTo(newNode)
          newNode%sibling        => node%sibling
          newNode%parent         => parentNode
          parentNode%firstChild  => newNode
          newNode%firstChild     => null()
          newNode%event          => null()
          newNode%firstSatellite => null()
          newNode%firstMergee    => null()
          newNode%mergeTarget    => null()
          newNode%siblingMergee  => null()
          newBasic               => newNode%basic()
          call newBasic%timeSet(newBasic%time()*(1.0d0-1.0d-6))
       end if
    else
       ! Node to be pruned is not the primary progenitor - simply unlink it from its siblings.
       workNode => parentNode%firstChild
       do while (.not.associated(workNode%sibling,node))
          workNode => workNode%sibling
       end do
       workNode%sibling => node%sibling
    end if
    return
  end subroutine Merger_Tree_Prune_Unlink_Parent

  subroutine Merger_Tree_Prune_Uniqueify_IDs(tree)
    !!{
    Ensure that nodes cloned during tree pruning have unique IDs.
    !!}
    use :: Galacticus_Nodes   , only : mergerTree                   , treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    implicit none
    type(mergerTree                   ), target , intent(in   ) :: tree
    type(treeNode                     ), pointer                :: node
    type(mergerTreeWalkerIsolatedNodes)                         :: treeWalker

    ! Iterate over trees.
    treeWalker=mergerTreeWalkerIsolatedNodes(tree,spanForest=.true.)
    do while (treeWalker%next(node))
       if (node%uniqueID() == node%parent%uniqueID()) call node%uniqueIDSet()
    end do
    return
  end subroutine Merger_Tree_Prune_Uniqueify_IDs

end module Merger_Trees_Pruning_Utilities
