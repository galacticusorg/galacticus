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
Implements a pruning operator on merger trees that removes all branches that can not contain any baryons.
!!}

  use :: Accretion_Halos        , only : accretionHalo        , accretionHaloClass
  use :: Virial_Density_Contrast, only : virialDensityContrast, virialDensityContrastClass

  !![
  <mergerTreeOperator name="mergerTreeOperatorPruneBaryons">
   <description>Provides a pruning operator on merger trees that removes all branches that can not contain any baryons.</description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneBaryons
     !!{
     A pruning operator on merger trees that removes all branches that can not contain any baryons.
     !!}
     private
     class(accretionHaloClass        ), pointer :: accretionHalo_         => null()
     class(virialDensityContrastClass), pointer :: virialDensityContrast_ => null()
   contains
     final     ::                        pruneBaryonsDestructor
     procedure :: operatePreEvolution => pruneBaryonsOperatePreEvolution
  end type mergerTreeOperatorPruneBaryons

  interface mergerTreeOperatorPruneBaryons
     !!{
     Constructors for the pruning-by-mass merger tree operator class.
     !!}
     module procedure pruneBaryonsConstructorParameters
     module procedure pruneBaryonsConstructorInternal
  end interface mergerTreeOperatorPruneBaryons

contains

  function pruneBaryonsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the prune-baryons merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (mergerTreeOperatorPruneBaryons)                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(accretionHaloClass            ), pointer       :: accretionHalo_
    class(virialDensityContrastClass    ), pointer       :: virialDensityContrast_

    !![
    <objectBuilder class="accretionHalo"         name="accretionHalo_"         source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self=mergerTreeOperatorPruneBaryons(accretionHalo_,virialDensityContrast_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="accretionHalo_"        />
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function pruneBaryonsConstructorParameters

  function pruneBaryonsConstructorInternal(accretionHalo_,virialDensityContrast_) result(self)
    !!{
    Internal constructor for the prune-baryons merger tree operator class.
    !!}
    implicit none
    type (mergerTreeOperatorPruneBaryons)                        :: self
    class(accretionHaloClass            ), intent(in   ), target :: accretionHalo_
    class(virialDensityContrastClass    ), intent(in   ), target :: virialDensityContrast_
    !![
    <constructorAssign variables="*accretionHalo_, *virialDensityContrast_"/>
    !!]

    return
  end function pruneBaryonsConstructorInternal

  subroutine pruneBaryonsDestructor(self)
    !!{
    Destructor for the merger tree operator function class.
    !!}
    implicit none
    type(mergerTreeOperatorPruneBaryons), intent(inout) :: self

    !![
    <objectDestructor name="self%accretionHalo_"        />
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]
    return
  end subroutine pruneBaryonsDestructor

  subroutine pruneBaryonsOperatePreEvolution(self,tree)
    !!{
    Prune branches from {\normalfont \ttfamily tree}.
    !!}
    use :: Galacticus_Nodes              , only : mergerTree                    , nodeComponentBasic             , treeNode
    use :: Merger_Tree_Walkers           , only : mergerTreeWalkerIsolatedNodes
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Clean_Branch, Merger_Tree_Prune_Uniqueify_IDs, Merger_Tree_Prune_Unlink_Parent
    implicit none
    class           (mergerTreeOperatorPruneBaryons), intent(inout), target :: self
    type            (mergerTree                    ), intent(inout), target :: tree
    type            (treeNode                      ), pointer               :: nodeNext       , node
    class           (nodeComponentBasic            ), pointer               :: basic
    type            (mergerTree                    ), pointer               :: treeCurrent
    type            (mergerTreeWalkerIsolatedNodes )                        :: treeWalker
    logical                                                                 :: didPruning
    double precision                                                        :: densityContrast

    !![
    <workaround type="unknown">
     <description>
       When using the percolation virial density contrast type we run into problems with initialization if we do not
       explicitly cause the percolation tables to be initialized here. The problem is a segmentation fault - I have not
       understood why this happens.
     </description>
    !!]
    ! Compute a density contrast.
    node            => tree                       %nodeBase
    basic           => node                       %basic          (                         )
    densityContrast =  self%virialDensityContrast_%densityContrast(basic%mass(),basic%time())
    !![
    </workaround>
    !!]
    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       didPruning=.true.
       do while (didPruning)
          didPruning=.false.
          ! Get root node of the tree.
          node  => treeCurrent%nodeBase
          basic => node       %basic   ()
          ! Check if the tree has any baryons.
          if (.not.self%accretionHalo_%branchHasBaryons(node)) then
             ! Entire tree can be pruned. Destroy all but this base node. (Leaving just
             ! the base node makes the tree inert - i.e. it can not do anything.)
             node => node%firstChild
             do while (associated(node))
                nodeNext => node%sibling
                call node%destroyBranch()
                deallocate(node)
                node => nodeNext
             end do
          else
             ! Walk the tree, pruning branches.
             treeWalker=mergerTreeWalkerIsolatedNodes(treeCurrent)
             do while (treeWalker%next(node))
                basic => node%basic()
                if (.not.self%accretionHalo_%branchHasBaryons(node).and.node%uniqueID() /= node%parent%uniqueID()) then
                   didPruning=.true.
                   ! Set the tree walker back to the base node so that it exits - we've changed the tree structure so need to
                   ! begin the walk again.
                   call treeWalker%setNode(treeCurrent%nodeBase)
                   ! Decouple from other nodes.
                   call Merger_Tree_Prune_Unlink_Parent(node,node%parent,.not.self%accretionHalo_%branchHasBaryons(node%parent),.true.)
                   ! Clean the branch.
                   call Merger_Tree_Prune_Clean_Branch(node)
                   ! Destroy the branch.
                   call node%destroyBranch()
                   deallocate(node)
                end if
             end do
          end if
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    ! Uniqueify nodes.
    call Merger_Tree_Prune_Uniqueify_IDs(tree)
    return
  end subroutine pruneBaryonsOperatePreEvolution
