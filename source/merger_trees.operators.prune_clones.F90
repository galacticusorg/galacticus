!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Implements a pruning-by-mass operator on merger trees.
!!}

  !![
  <mergerTreeOperator name="mergerTreeOperatorPruneClones">
   <description>Provides a clone pruning operator on merger trees.</description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneClones
     !!{
     A clone pruning merger tree operator class.
     !!}
     private
   contains
     procedure :: operatePreEvolution => pruneClonesOperatePreEvolution
  end type mergerTreeOperatorPruneClones

  interface mergerTreeOperatorPruneClones
     !!{
     Constructors for the clone pruning merger tree operator class.
     !!}
     module procedure pruneClonesConstructorParameters
  end interface mergerTreeOperatorPruneClones

contains

  function pruneClonesConstructorParameters(parameters) result(self)
    !!{
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
    !!{
    Perform a clone pruning operation on a merger tree.
    !!}
    use :: Galacticus_Nodes              , only : mergerTree                    , nodeComponentBasic             , treeNode
    use :: Merger_Tree_Walkers           , only : mergerTreeWalkerIsolatedNodes
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Clean_Branch, Merger_Tree_Prune_Unlink_Parent
    use :: Numerical_Comparison          , only : Values_Agree
    implicit none
    class(mergerTreeOperatorPruneClones), intent(inout), target :: self
    type (mergerTree                   ), intent(inout), target :: tree
    type (treeNode                     ), pointer               :: node      , nodePrevious
    class(nodeComponentBasic           ), pointer               :: basic     , basicPrevious
    type (mergerTreeWalkerIsolatedNodes)                        :: treeWalker
    !$GLC attributes unused :: self

    ! Iterate over nodes.
    treeWalker=mergerTreeWalkerIsolatedNodes(tree)
    do while (treeWalker%next(node))
       ! Skip this node if it has no parent.
       if (associated(node%parent)) then
          ! Get basic components.
          basic         => node        %basic()
          basicPrevious => nodePrevious%basic()
          if (Values_Agree(basic%time(),basicPrevious%time(),relTol=1.0d-5)) then
             ! Return to the previous node.
             call treeWalker%previous(nodePrevious)
             ! Decouple from other nodes.
             call Merger_Tree_Prune_Unlink_Parent(node,nodePrevious,.false.,.false.)
             ! Clean the branch.
             call Merger_Tree_Prune_Clean_Branch (node                             )
             ! Destroy the branch.
             call node%destroyBranch()
          end if
       end if
    end do
    return
  end subroutine pruneClonesOperatePreEvolution
