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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Implements a clone-pruning operator on merger trees.
!!}

  !![
  <mergerTreeOperator name="mergerTreeOperatorPruneClones" docformat="rst">
   <description>
   Provides a merger tree operator which removes cloned nodes. A clone is identified as a node with no progenitors whose time
   agrees with that of its parent to within a fractional tolerance of :math:`10^{-5}`. Such nodes are inserted by
   :galacticus-class:`mergerTreeConstructorRead` (as a primary progenitor for a node which otherwise has none), and by some merger
   tree pruning operators (to preserve primary progenitor status). Any satellites of the clone are returned to its parent, and
   any nodes which were to merge with the clone will instead merge with its parent.

   A clone is removed only if it is the sole progenitor of its parent. A clone with siblings exists to be the primary progenitor
   of its parent, so that none of those siblings (which are intended to become subhalos) is promoted into the parent; removing it
   would make one of them the primary progenitor, and the tree could then not be evolved correctly. Clones which have events
   attached (e.g. those which act as the target of a subhalo promotion) are also left in place, as they can not be removed
   consistently. Clones become sole progenitors when their siblings are removed, for example by another pruning operator applied
   first.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneClones
     !!{RST
     A clone pruning merger tree operator class.
     !!}
     private
   contains
     procedure :: operatePreEvolution => pruneClonesOperatePreEvolution
  end type mergerTreeOperatorPruneClones

  interface mergerTreeOperatorPruneClones
     !!{RST
     Constructors for the clone pruning merger tree operator class.
     !!}
     module procedure pruneClonesConstructorParameters
  end interface mergerTreeOperatorPruneClones

contains

  function pruneClonesConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the clone pruning merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeOperatorPruneClones)                :: self
    type(inputParameters              ), intent(inout) :: parameters

    self=mergerTreeOperatorPruneClones()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function pruneClonesConstructorParameters

  subroutine pruneClonesOperatePreEvolution(self,tree)
    !!{RST
    Perform a clone pruning operation on a merger tree.
    !!}
    use :: Galacticus_Nodes              , only : mergerTree                    , nodeComponentBasic             , treeNode
    use :: Merger_Tree_Walkers           , only : mergerTreeWalkerIsolatedNodes
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Clean_Branch, Merger_Tree_Prune_Unlink_Parent
    use :: Numerical_Comparison          , only : Values_Agree
    implicit none
    class(mergerTreeOperatorPruneClones), intent(inout), target :: self
    type (mergerTree                   ), intent(inout), target :: tree
    type (treeNode                     ), pointer               :: node      , nodeParent , &
         &                                                         nodeWork  , nodeLast
    class(nodeComponentBasic           ), pointer               :: basic     , basicParent
    type (mergerTreeWalkerIsolatedNodes)                        :: treeWalker
    !$GLC attributes unused :: self

    ! Iterate over nodes in all trees of the forest.
    treeWalker=mergerTreeWalkerIsolatedNodes(tree,spanForest=.true.)
    do while (treeWalker%next(node))
       ! Skip this node if it has no parent, or has progenitors (clones never do).
       if (.not.associated(node%parent).or.associated(node%firstChild)) cycle
       ! Skip nodes which are not the sole progenitor of their parent. Removing such a clone would make one of its siblings the
       ! primary progenitor, although the clone was inserted precisely to prevent this.
       if (.not.associated(node%parent%firstChild,node).or.associated(node%sibling)) cycle
       ! Skip nodes with attached events. These can not be removed consistently as other nodes may hold paired events which
       ! point to them (e.g. a clone inserted as the target for a subhalo promotion).
       if (associated(node%event)) cycle
       ! Test for a clone - a node at the same time as its parent.
       nodeParent  => node      %parent
       basic       => node      %basic ()
       basicParent => nodeParent%basic ()
       if (.not.Values_Agree(basic%time(),basicParent%time(),relTol=1.0d-5)) cycle
       ! Return any satellites of the clone to its parent.
       if (associated(node%firstSatellite)) then
          nodeWork => node%firstSatellite
          do while (associated(nodeWork))
             nodeWork%parent => nodeParent
             nodeWork        => nodeWork%sibling
          end do
          if (associated(nodeParent%firstSatellite)) then
             nodeLast                  => nodeParent%lastSatellite()
             nodeLast  %sibling        => node      %firstSatellite
          else
             nodeParent%firstSatellite => node      %firstSatellite
          end if
          node%firstSatellite => null()
       end if
       ! Any nodes which were to merge with the clone will instead merge with its parent.
       if (associated(node%firstMergee)) then
          nodeWork => node%firstMergee
          do while (associated(nodeWork))
             nodeWork%mergeTarget => nodeParent
             nodeLast             => nodeWork
             nodeWork             => nodeWork  %siblingMergee
          end do
          nodeLast  %siblingMergee => nodeParent%firstMergee
          nodeParent%firstMergee   => node      %firstMergee
          node      %firstMergee   => null()
       end if
       ! Step the walker back to the previously visited node. The walk is depth-first, visiting progenitors before their
       ! descendants, and the clone has no progenitors, so once the clone is removed, the next step from the previous node
       ! reaches the node which would have followed the clone.
       call treeWalker%previous(nodeWork)
       ! Decouple from other nodes.
       call Merger_Tree_Prune_Unlink_Parent(node,nodeParent,.false.,.false.)
       ! Clean the branch.
       call Merger_Tree_Prune_Clean_Branch (node                           )
       ! Destroy the branch.
       call node%destroyBranch()
       deallocate(node)
    end do
    return
  end subroutine pruneClonesOperatePreEvolution
