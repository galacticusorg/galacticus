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

!!{RST
Implements a pruning-by-filter operator on merger trees.
!!}

!+ Contributions to this file made by: Andrew Benson
!+ This file was generated with assistance from Claude, and reviewed and tested by Andrew Benson.

  use :: Galactic_Filters, only : galacticFilterClass

  !![
  <mergerTreeOperator name="mergerTreeOperatorPruneByFilter" docformat="rst">
   <description>
   A merger tree operator class which prunes branches of merger trees prior to any evolution. Note that the sense of the filter is
   that nodes which *pass* ``[galacticFilter]`` are those which are *removed* from the tree, along with all of their
   progenitors. (To prune the complement---i.e. to retain only those nodes which pass some filter---wrap that filter in a
   :galacticus-class:`galacticFilterNot`.) This is intended for cases in which some earlier stage of tree initialization has
   identified nodes which should not be present in the tree, for example by labeling them (see
   :galacticus-class:`galacticFilterLabelled`).
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneByFilter
     !!{RST
     A pruning-by-filter merger tree operator class.
     !!}
     private
     class  (galacticFilterClass), pointer :: galacticFilter_           => null()
     logical                               :: preservePrimaryProgenitor
   contains
     final     ::                        pruneByFilterDestructor
     procedure :: operatePreEvolution => pruneByFilterOperatePreEvolution
  end type mergerTreeOperatorPruneByFilter

  interface mergerTreeOperatorPruneByFilter
     !!{RST
     Constructors for the prune-by-filter merger tree operator class.
     !!}
     module procedure pruneByFilterConstructorParameters
     module procedure pruneByFilterConstructorInternal
  end interface mergerTreeOperatorPruneByFilter

contains

  function pruneByFilterConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the prune-by-filter merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (mergerTreeOperatorPruneByFilter)                :: self
    type   (inputParameters                ), intent(inout) :: parameters
    class  (galacticFilterClass            ), pointer       :: galacticFilter_
    logical                                                 :: preservePrimaryProgenitor

    !![
    <inputParameter docformat="rst">
      <name>preservePrimaryProgenitor</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, primary progenitor status is preserved even if the primary progenitor is pruned from the tree.
      </description>
    </inputParameter>
    <objectBuilder class="galacticFilter" name="galacticFilter_" source="parameters"/>
    !!]
    self=mergerTreeOperatorPruneByFilter(preservePrimaryProgenitor,galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"/>
    !!]
    return
  end function pruneByFilterConstructorParameters

  function pruneByFilterConstructorInternal(preservePrimaryProgenitor,galacticFilter_) result(self)
    !!{RST
    Internal constructor for the prune-by-filter merger tree operator class.
    !!}
    implicit none
    type   (mergerTreeOperatorPruneByFilter)                        :: self
    class  (galacticFilterClass            ), intent(in   ), target :: galacticFilter_
    logical                                 , intent(in   )         :: preservePrimaryProgenitor
    !![
    <constructorAssign variables="preservePrimaryProgenitor, *galacticFilter_"/>
    !!]

    return
  end function pruneByFilterConstructorInternal

  subroutine pruneByFilterDestructor(self)
    !!{RST
    Destructor for the prune-by-filter merger tree operator class.
    !!}
    implicit none
    type(mergerTreeOperatorPruneByFilter), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine pruneByFilterDestructor

  subroutine pruneByFilterOperatePreEvolution(self,tree)
    !!{RST
    Perform a prune-by-filter operation on a merger tree.
    !!}
    use :: Galacticus_Nodes              , only : mergerTree                    , treeNode
    use :: Merger_Tree_Walkers           , only : mergerTreeWalkerIsolatedNodes
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Clean_Branch, Merger_Tree_Prune_Uniqueify_IDs, Merger_Tree_Prune_Unlink_Parent
    implicit none
    class  (mergerTreeOperatorPruneByFilter), intent(inout), target :: self
    type   (mergerTree                     ), intent(inout), target :: tree
    type   (treeNode                       ), pointer               :: nodeNext   , nodeWork, &
         &                                                             node
    type   (mergerTree                     ), pointer               :: currentTree
    type   (mergerTreeWalkerIsolatedNodes  )                        :: treeWalker
    logical                                                         :: didPruning , parentWillBePruned

    ! Iterate over trees.
    currentTree => tree
    do while (associated(currentTree))
       didPruning=.true.
       do while (didPruning)
          didPruning=.false.
          ! Get the root node of the tree.
          node => currentTree%nodeBase
          if (self%galacticFilter_%passes(node)) then
             ! The base node itself is to be pruned. As we can not remove the base node, destroy all but this base node instead -
             ! leaving just the base node makes the tree inert (i.e. it can not do anything).
             call Merger_Tree_Prune_Clean_Branch(node)
             nodeWork => node%firstChild
             do while (associated(nodeWork))
                nodeNext => nodeWork%sibling
                call Merger_Tree_Prune_Clean_Branch(nodeWork)
                call nodeWork%destroyBranch()
                deallocate(nodeWork)
                nodeWork => nodeNext
             end do
             nullify(node%firstChild)
             nodeWork => node%firstSatellite
             do while (associated(nodeWork))
                nodeNext => nodeWork%sibling
                call Merger_Tree_Prune_Clean_Branch(nodeWork)
                call nodeWork%destroyBranch()
                deallocate(nodeWork)
                nodeWork => nodeNext
             end do
             nullify(node%firstSatellite)
          else
             ! Walk the tree, pruning branches.
             treeWalker=mergerTreeWalkerIsolatedNodes(currentTree)
             do while (treeWalker%next(node))
                if (self%galacticFilter_%passes(node)) then
                   didPruning=.true.
                   ! Set the tree walker back to the base node so that it exits - we've changed the tree structure so need to
                   ! begin the walk again.
                   call treeWalker%setNode(currentTree%nodeBase)
                   ! Determine whether the parent will also be pruned - if it will, there is no need to preserve primary
                   ! progenitor status within it.
                   parentWillBePruned=self%galacticFilter_%passes(node%parent)
                   ! Decouple from other nodes.
                   call Merger_Tree_Prune_Unlink_Parent(node,node%parent,parentWillBePruned,self%preservePrimaryProgenitor)
                   ! Clean the branch.
                   call Merger_Tree_Prune_Clean_Branch (node)
                   ! Destroy the branch.
                   call node%destroyBranch()
                   deallocate(node)
                end if
             end do
          end if
       end do
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
    ! Uniqueify nodes.
    call Merger_Tree_Prune_Uniqueify_IDs(tree)
    return
  end subroutine pruneByFilterOperatePreEvolution
