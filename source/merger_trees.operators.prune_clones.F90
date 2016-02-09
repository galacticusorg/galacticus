!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a pruning-by-mass operator on merger trees.

  !# <mergerTreeOperator name="mergerTreeOperatorPruneClones">
  !#  <description>Provides a clone pruning operator on merger trees.</description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneClones
     !% A clone pruning merger tree operator class.
     private
   contains
     final     ::            pruneClonesDestructor
     procedure :: operate => pruneClonesOperate
  end type mergerTreeOperatorPruneClones

  interface mergerTreeOperatorPruneClones
     !% Constructors for the clone pruning merger tree operator class.
     module procedure pruneClonesConstructorParameters
     module procedure pruneClonesConstructorInternal
  end interface mergerTreeOperatorPruneClones

contains

  function pruneClonesConstructorParameters(parameters)
    !% Constructor for the clone pruning merger tree operator class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(mergerTreeOperatorPruneClones)                :: pruneClonesConstructorParameters
    type(inputParameters              ), intent(in   ) :: parameters

    return
  end function pruneClonesConstructorParameters

  function pruneClonesConstructorInternal()
    !% Internal constructor for the clone pruning merger tree operator class.
    implicit none
    type(mergerTreeOperatorPruneClones) :: pruneClonesConstructorInternal

    return
  end function pruneClonesConstructorInternal

  elemental subroutine pruneClonesDestructor(self)
    !% Destructor for the merger tree operator function class.
    implicit none
    type(mergerTreeOperatorPruneClones), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine pruneClonesDestructor

  subroutine pruneClonesOperate(self,tree)
    !% Perform a clone pruning operation on a merger tree.
    use Merger_Trees_Pruning_Utilities
    use Numerical_Comparison
    implicit none
    class(mergerTreeOperatorPruneClones), intent(inout)         :: self
    type (mergerTree                   ), intent(inout), target :: tree
    type (treeNode                     ), pointer               :: node       , nodePrevious
    class(nodeComponentBasic           ), pointer               :: basic      , basicPrevious
    type (mergerTree                   ), pointer               :: currentTree

    ! Iterate over trees.
    currentTree => tree
    do while (associated(currentTree))
       ! Get root node of the tree.
       node  => currentTree%baseNode
       ! Walk the tree, pruning clones.
       do while (associated(node))
          ! Skip this node if it has no parent.
          if (associated(node%parent)) then
             ! Record the parent node to which we will return.
             nodePrevious => node%parent
             ! Get basic components.
             basic         => node        %basic()
             basicPrevious => nodePrevious%basic()
             if (Values_Agree(basic%time(),basicPrevious%time(),relTol=1.0d-5)) then
                ! Decouple from other nodes.
                call Merger_Tree_Prune_Unlink_Parent(node,nodePrevious,.false.,.false.)
                ! Clean the branch.
                call Merger_Tree_Prune_Clean_Branch (node                             )
                ! Destroy the branch.
                call currentTree%destroyBranch(node)
                ! Return to parent node.
                node => nodePrevious
             end if
          end if
          node => node%walkTree()
       end do
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
    return
  end subroutine pruneClonesOperate
