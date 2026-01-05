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
  Implements a merger tree operator which prunes tips of branches (i.e. sections from the leaf node to
  the first node with a sibling).
  !!}

  !![
  <mergerTreeOperator name="mergerTreeOperatorPruneBranchTips">
   <description>Complements a merger tree operator which prunes tips of branches (i.e. sections from the leaf node to the first node with a sibling).</description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneBranchTips
     !!{
     A merger tree operator class which prunes tips of branches (i.e. sections from the leaf node to the first node with a
     sibling).
     !!}
     private
   contains
     procedure :: operatePreEvolution => pruneBranchTipsOperatePreEvolution
  end type mergerTreeOperatorPruneBranchTips

  interface mergerTreeOperatorPruneBranchTips
     !!{
     Constructors for the prune-branchTips merger tree operator class.
     !!}
     module procedure pruneBranchTipsConstructorParameters
  end interface mergerTreeOperatorPruneBranchTips

contains

  function pruneBranchTipsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the prune-branchTips merger tree operator class which takes a parameter set as input.
    !!}
    implicit none
    type   (mergerTreeOperatorPruneBranchTips)                :: self
    type   (inputParameters                  ), intent(inout) :: parameters

    self=mergerTreeOperatorPruneBranchTips()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function pruneBranchTipsConstructorParameters

  subroutine pruneBranchTipsOperatePreEvolution(self,tree)
    !!{
    Perform a prune-branchTips operation on a merger tree.
    !!}
    use :: Merger_Tree_Walkers           , only : mergerTreeWalkerIsolatedNodes
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Clean_Branch, Merger_Tree_Prune_Uniqueify_IDs, Merger_Tree_Prune_Unlink_Parent
    implicit none
    class(mergerTreeOperatorPruneBranchTips), intent(inout), target  :: self
    type (mergerTree                       ), intent(inout), target  :: tree
    type (treeNode                         )               , pointer :: node      , nodeWork, nodeDestroy
    type (mergerTreeWalkerIsolatedNodes    )                         :: treeWalker
    !$GLC attributes unused :: self

    treeWalker=mergerTreeWalkerIsolatedNodes(tree,spanForest=.true.)
    do while (treeWalker%next(node))
       if     (                                             &
            &   .not.associated(node%firstChild           ) &
            &  .and.                                        &
            &   .not.associated(node%firstSatellite       ) &
            &  .and.                                        &
            &   .not.associated(node%sibling              ) &
            &  .and.                                        &
            &                   node%isPrimaryProgenitor()  &
            & ) then
          nodeWork => node
          do while (                                                 &
               &     .not.associated(nodeWork%firstSatellite       ) &
               &    .and.                                            &
               &     .not.associated(nodeWork%sibling              ) &
               &    .and.                                            &
               &                     nodeWork%isPrimaryProgenitor()  &
               &    .and.                                            &
               &          associated(nodeWork%parent               ) &
               &   )
             nodeWork => nodeWork%parent
          end do
          if (associated(nodeWork%firstChild)) then
             nodeDestroy => nodeWork%firstChild
             call Merger_Tree_Prune_Unlink_Parent(nodeDestroy,nodeWork,parentWillBePruned=.false.,preservePrimaryProgenitor=.false.)
             call Merger_Tree_Prune_Clean_Branch (nodeDestroy                                                                      )
             call nodeDestroy%destroyBranch()
             deallocate(nodeDestroy)
             call treeWalker%setNode(nodeWork)
          end if
       end if
    end do
    call Merger_Tree_Prune_Uniqueify_IDs(tree)
    return
  end subroutine pruneBranchTipsOperatePreEvolution

