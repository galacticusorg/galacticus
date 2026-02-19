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
  <mergerTreeOperator name="mergerTreeOperatorPruneByMass">
   <description>
    A merger tree operator class which allows for branches of merger trees to be pruned---i.e. nodes below a specified mass
    limit are removed from the tree prior to any evolution. This can be useful for convergence studies for example. Set
    {\normalfont \ttfamily [massThreshold]} to the desired mass threshold below which nodes will be pruned.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneByMass
     !!{
     A pruning-by-mass merger tree operator class.
     !!}
     private
     double precision :: massThreshold
     logical          :: preservePrimaryProgenitor
   contains
     procedure :: operatePreEvolution => pruneByMassOperatePreEvolution
  end type mergerTreeOperatorPruneByMass

  interface mergerTreeOperatorPruneByMass
     !!{
     Constructors for the pruning-by-mass merger tree operator class.
     !!}
     module procedure pruneByMassConstructorParameters
     module procedure pruneByMassConstructorInternal
  end interface mergerTreeOperatorPruneByMass

contains

  function pruneByMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the prune-by-mass merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeOperatorPruneByMass)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    double precision                                               :: massThreshold
    logical                                                        :: preservePrimaryProgenitor

    !![
    <inputParameter>
      <name>massThreshold</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>Threshold mass below which merger tree branches should be pruned.</description>
    </inputParameter>
    <inputParameter>
      <name>preservePrimaryProgenitor</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, primary progenitor status is preserved even if the primary progenitor is pruned from the tree.</description>
    </inputParameter>
    !!]
    self=mergerTreeOperatorPruneByMass(massThreshold,preservePrimaryProgenitor)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function pruneByMassConstructorParameters

  function pruneByMassConstructorInternal(massThreshold,preservePrimaryProgenitor) result(self)
    !!{
    Internal constructor for the prune-by-mass merger tree operator class.
    !!}
    implicit none
    type            (mergerTreeOperatorPruneByMass)                :: self
    double precision                               , intent(in   ) :: massThreshold
    logical                                        , intent(in   ) :: preservePrimaryProgenitor
    !![
    <constructorAssign variables="massThreshold, preservePrimaryProgenitor"/>
    !!]

    return
  end function pruneByMassConstructorInternal

  subroutine pruneByMassOperatePreEvolution(self,tree)
    !!{
    Perform a prune-by-mass operation on a merger tree.
    !!}
    use :: Galacticus_Nodes              , only : mergerTree                    , nodeComponentBasic             , treeNode
    use :: Merger_Tree_Walkers           , only : mergerTreeWalkerIsolatedNodes
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Clean_Branch, Merger_Tree_Prune_Uniqueify_IDs, Merger_Tree_Prune_Unlink_Parent
    implicit none
    class  (mergerTreeOperatorPruneByMass), intent(inout), target :: self
    type   (mergerTree                   ), intent(inout), target :: tree
    type   (treeNode                     ), pointer               :: nodeNext     , nodeWork   , &
         &                                                           node
    class  (nodeComponentBasic           ), pointer               :: basic        , basicParent
    type   (mergerTree                   ), pointer               :: currentTree
    type   (mergerTreeWalkerIsolatedNodes)                        :: treeWalker
    logical                                                       :: didPruning

    ! Iterate over trees.
    currentTree => tree
    do while (associated(currentTree))
       didPruning=.true.
       do while (didPruning)
          didPruning=.false.
          ! Get root node of the tree.
          node  => currentTree%nodeBase
          basic => node       %basic   ()
          if (basic%mass() < self%massThreshold) then
             ! Entire tree is below threshold. Destroy all but this base node. (Leaving just
             ! the base node makes the tree inert - i.e. it can not do anything.)
             call Merger_Tree_Prune_Clean_Branch (node)
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
                basic => node%basic()
                if (basic%mass() < self%massThreshold) then
                   didPruning=.true.
                   ! Set the tree walker back to the base node so that it exits - we've changed the tree structure so need to
                   ! begin the walk again.
                   call treeWalker%setNode(currentTree%nodeBase)
                   ! Decouple from other nodes.
                   basicParent => node%parent%basic()
                   call Merger_Tree_Prune_Unlink_Parent(node,node%parent,basicParent%mass() < self%massThreshold,self%preservePrimaryProgenitor)
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
  end subroutine pruneByMassOperatePreEvolution
