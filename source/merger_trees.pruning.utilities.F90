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

  subroutine Merger_Tree_Prune_Unlink_Parent(node,nodeParent,parentWillBePruned,preservePrimaryProgenitor)
    !!{
    Unlink a parent node from a tree branch which is about to be pruned.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    type   (treeNode          ), intent(inout), pointer :: node              , nodeParent
    logical                    , intent(in   )          :: parentWillBePruned, preservePrimaryProgenitor
    type   (treeNode          ), pointer                :: nodeNew           , nodeWork
    class  (nodeComponentBasic), pointer                :: basicNew
    integer                                             :: i
    
    ! Check primary progenitor status.
    if (node%isPrimaryProgenitorOf(nodeParent)) then
       ! Node is primary progenitor - we must check if the parent will be pruned also.
       if (parentWillBePruned.or..not.preservePrimaryProgenitor) then
          ! Parent will eventually be pruned - simply replace the current first child (about to pruned) with its
          ! sibling. Alternatively, we've been asked to not preserve primary progenitor status by inserting cloned nodes.
          nodeParent%firstChild => node%sibling
       else
          ! Parent will not be pruned. Insert a clone of the parent as its own first progenitor to prevent any siblings of the
          ! to-be-pruned branch being misidentified as the primary progenitor.
          allocate(nodeNew)
          call nodeParent%copyNodeTo(nodeNew)
          nodeNew%sibling        => node%sibling
          nodeNew%parent         => nodeParent
          nodeParent%firstChild  => nodeNew
          nodeNew%firstChild     => null()
          nodeNew%event          => null()
          nodeNew%firstSatellite => null()
          nodeNew%firstMergee    => null()
          nodeNew%mergeTarget    => null()
          nodeNew%siblingMergee  => null()
          basicNew               => nodeNew%basic()
          call basicNew%timeSet(basicNew%time()*(1.0d0-1.0d-6))
          if (nodeNew%satelliteCount() > 0) then
             ! Remove any satellite component from the copied node - each branch should have only a single satellite.
             do i=nodeNew%satelliteCount(),1,-1
                call nodeNew%satelliteRemove(i)
             end do
          end if          
       end if
    else
       ! Node to be pruned is not the primary progenitor - simply unlink it from its siblings.
       nodeWork => nodeParent%firstChild
       do while (.not.associated(nodeWork%sibling,node))
          nodeWork => nodeWork%sibling
       end do
       nodeWork%sibling => node%sibling
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
       if (.not.associated(node%parent)) cycle
       if (node%uniqueID() == node%parent%uniqueID()) call node%uniqueIDSet()
    end do
    return
  end subroutine Merger_Tree_Prune_Uniqueify_IDs

end module Merger_Trees_Pruning_Utilities
