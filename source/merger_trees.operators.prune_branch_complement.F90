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
  Implements a merger tree operator which prunes all but a specified branch.
  !!}

  use :: Kind_Numbers, only : kind_int8

  !![
  <mergerTreeOperator name="mergerTreeOperatorPruneBranchComplement">
   <description>
    A merger tree operator class which prunes all but the branch starting from {\normalfont \ttfamily [branchNodeID]}.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneBranchComplement
     !!{
     A merger tree operator class which prunes all but a specified branch.
     !!}
     private
     integer(kind=kind_int8) :: branchNodeID
   contains
     procedure :: operatePreEvolution => pruneBranchComplementOperatePreEvolution
  end type mergerTreeOperatorPruneBranchComplement

  interface mergerTreeOperatorPruneBranchComplement
     !!{
     Constructors for the prune-non-essential merger tree operator class.
     !!}
     module procedure pruneBranchComplementConstructorParameters
     module procedure pruneBranchComplementConstructorInternal
  end interface mergerTreeOperatorPruneBranchComplement

contains

  function pruneBranchComplementConstructorParameters(parameters) result(self)
    !!{
    Constructor for the prune-non-essential merger tree operator class which takes a parameter set as input.
    !!}
    implicit none
    type   (mergerTreeOperatorPruneBranchComplement)                :: self
    type   (inputParameters                        ), intent(inout) :: parameters
    integer(kind_int8                              )                :: branchNodeID
    
    !![
    <inputParameter>
      <name>branchNodeID</name>
      <source>parameters</source>
      <variable>branchNodeID</variable>
      <description>ID of the node at the bash of the branch to avoid pruning.</description>
    </inputParameter>
    !!]
    self=mergerTreeOperatorPruneBranchComplement(branchNodeID)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function pruneBranchComplementConstructorParameters

  function pruneBranchComplementConstructorInternal(branchNodeID) result(self)
    !!{
    Internal constructor for the prune-non-essential merger tree operator class.
    !!}
    implicit none
    type   (mergerTreeOperatorPruneBranchComplement)                :: self
    integer(kind_int8                              ), intent(in   ) :: branchNodeID
    !![
    <constructorAssign variables="branchNodeID"/>
    !!]
    
   return
  end function pruneBranchComplementConstructorInternal

  subroutine pruneBranchComplementOperatePreEvolution(self,tree)
    !!{
    Perform a prune-branch-complement operation on a merger tree.
    !!}
    use :: Galacticus_Nodes              , only : mergerTree                    , treeNode, nodeCOmponentBasic
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Clean_Branch
    implicit none
    class(mergerTreeOperatorPruneBranchComplement), intent(inout), target  :: self
    type (mergerTree                             ), intent(inout), target  :: tree
    type (mergerTree                             )               , pointer :: treeCurrent , treeHost  , &
         &                                                                    treeNext
    type (treeNode                               )               , pointer :: node        , nodeBranch, &
         &                                                                    nodePrevious

    ! Iterate over trees.
    treeCurrent => tree
    treeHost    => null()
    do while (associated(treeCurrent))
       ! Find the branch node.
       nodeBranch => treeCurrent%getNode(self%branchNodeID)
       if (associated(nodeBranch)) then
          ! Uncouple the node from any parent.
          if (associated(nodeBranch%parent)) then
             node => nodeBranch%parent%firstChild
             if (associated(node,nodeBranch)) then
                nodeBranch%parent%firstChild => nodeBranch%sibling
             else
                nodePrevious => null()
                do while (associated(node))
                   if (associated(node,nodeBranch)) then
                      nodePrevious%sibling => node%sibling
                      exit
                   end if
                   nodePrevious => node
                   node         => node%sibling
                end do
             end if
             nodeBranch%parent  => null()
             nodeBranch%sibling => null()
          end if
          ! Destroy the original tree.
          call treeCurrent%nodeBase%destroy()
          deallocate(treeCurrent%nodeBase)
          treeCurrent%nodeBase => nodeBranch
          treeHost             => treeCurrent
       end if
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    ! Set the first tree in the forest to be the tree hosting the target branch.
    treeCurrent => tree
    do while (associated(treeCurrent))
       treeNext => treeCurrent%nextTree
       ! Destroy all nodes in the tree if it is not the one hosting the target node.
       if (.not.associated(treeCurrent,treeHost)) call treeCurrent%nodeBase%destroy()
       ! Destroy the tree itself unless it is the initial tree, or the one hosting the target node.
       if (.not.associated(treeCurrent,treeHost).and..not.associated(treeCurrent,tree)) then
          call treeCurrent%destroy()
          deallocate(treeCurrent)
       end if
       treeCurrent => treeNext
    end do
    ! Assign the tree containing the target node to the initial tree in the forest.
    tree%nodeBase     => treeHost%nodeBase
    tree%nextTree     => null()
    tree%firstTree    => tree
    tree%index        =  treeHost%index
    tree%volumeWeight =  treeHost%volumeWeight
    ! Destroy the original host tree (unless it was already the first tree in the forest).
    if (.not.associated(treeHost,tree)) then
       call treeHost%destroy()
       deallocate(treeHost)
    end if
    return
  end subroutine pruneBranchComplementOperatePreEvolution
