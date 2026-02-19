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
Implements a merger tree operator that prunes branches below some maximum mass, and all of the tree \emph{after} the final output time.
!!}

  use :: Output_Times, only : outputTimesClass
  
  !![
  <mergerTreeOperator name="mergerTreeOperatorPruneByMassAndTime">
   <description>
    A merger tree operator class that prunes branches below some maximum mass, and all of the tree \emph{after} the final output time.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneByMassAndTime
     !!{
     A merger tree operator that prunes branches below some maximum mass, and all of the tree \emph{after} the final output time.
     !!}
     private
     class           (outputTimesClass), pointer :: outputTimes_  => null()
     double precision                            :: massThreshold
   contains
     final     ::                        pruneByMassAndTimeDestructor
     procedure :: operatePreEvolution => pruneByMassAndTimeOperatePreEvolution
  end type mergerTreeOperatorPruneByMassAndTime

  interface mergerTreeOperatorPruneByMassAndTime
     !!{
     Constructors for the pruning-by-mass merger tree operator class.
     !!}
     module procedure pruneByMassAndTimeConstructorParameters
     module procedure pruneByMassAndTimeConstructorInternal
  end interface mergerTreeOperatorPruneByMassAndTime

contains

  function pruneByMassAndTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the prune-by-mass merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeOperatorPruneByMassAndTime)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (outputTimesClass                    ), pointer       :: outputTimes_
    double precision                                                      :: massThreshold

    !![
    <inputParameter>
      <name>massThreshold</name>
      <source>parameters</source>
      <description>Threshold mass below which merger tree branches should be pruned.</description>
    </inputParameter> 
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    self=mergerTreeOperatorPruneByMassAndTime(massThreshold,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function pruneByMassAndTimeConstructorParameters

  function pruneByMassAndTimeConstructorInternal(massThreshold,outputTimes_) result(self)
    !!{
    Internal constructor for the prune-by-mass merger tree operator class.
    !!}
    implicit none
    type            (mergerTreeOperatorPruneByMassAndTime)                        :: self
    double precision                                      , intent(in   )         :: massThreshold
    class           (outputTimesClass                    ), intent(in   ), target :: outputTimes_
    !![
    <constructorAssign variables="massThreshold, *outputTimes_"/>
    !!]

    return
  end function pruneByMassAndTimeConstructorInternal

  subroutine pruneByMassAndTimeDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeOperatorPruneByMassAndTime} output times class.
    !!}
    implicit none
    type(mergerTreeOperatorPruneByMassAndTime), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"/>
    !!]
    return
  end subroutine pruneByMassAndTimeDestructor

  subroutine pruneByMassAndTimeOperatePreEvolution(self,tree)
    !!{
    Perform a prune-by-mass operation on a merger tree.
    !!}
    use :: Galacticus_Nodes              , only : mergerTree                    , nodeComponentBasic                 , treeNode                       , treeNodeLinkedList
    use :: Merger_Tree_Walkers           , only : mergerTreeWalkerIsolatedNodes , mergerTreeWalkerIsolatedNodesBranch, mergerTreeWalkerAllNodes
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Clean_Branch, Merger_Tree_Prune_Uniqueify_IDs    , Merger_Tree_Prune_Unlink_Parent
    implicit none
    class           (mergerTreeOperatorPruneByMassAndTime), intent(inout), target :: self
    type            (mergerTree                          ), intent(inout), target :: tree
    type            (treeNode                            ), pointer               :: nodeNext    , nodeBranch     , &
         &                                                                           node        , nodeRoot
    class           (nodeComponentBasic                  ), pointer               :: basic       , basicChild     , &
         &                                                                           basicBranch
    type            (mergerTree                          ), pointer               :: currentTree , treeNext
    type            (treeNodeLinkedList                  ), pointer               :: nodeRootHead, nodeRootCurrent, &
         &                                                                           nodeRootNext
    type            (mergerTreeWalkerIsolatedNodes       )                        :: treeWalker
    type            (mergerTreeWalkerAllNodes            )                        :: allWalker
    type            (mergerTreeWalkerIsolatedNodesBranch )                        :: branchWalker
    double precision                                                              :: timeFinal
    logical                                                                       :: pruneBranch , nodesRemain    , &
         &                                                                           keepTree

    ! Make a copy of the root node.
    allocate(nodeRoot)
    call tree%nodeBase%copyNodeTo(nodeRoot)
    ! Get the final output time.
    timeFinal=self%outputTimes_%time(self%outputTimes_%count())    
    ! Walk the trees.
    treeWalker      =  mergerTreeWalkerIsolatedNodes(tree,spanForest=.true.)
    nodesRemain     =  treeWalker%next(node)
    nodeRootHead    => null()
    nodeRootCurrent => null()
    do while (nodesRemain)
       nodesRemain=treeWalker%next(nodeNext)
       ! Skip nodes with no child
       if (associated(node%firstChild)) then
          ! Check for nodes after the final time, with a child before the final time.
          basic      => node           %basic()
          basicChild => node%firstChild%basic()
          if     (                                &
               &   basic     %time() >= timeFinal &
               &  .and.                           &
               &   basicChild%time() <  timeFinal &
               & ) then
             ! Check if any node in this branch was above threshold.
             pruneBranch =.true.
             branchWalker=mergerTreeWalkerIsolatedNodesBranch(node)
             do while (branchWalker%next(nodeBranch))
                basicBranch => nodeBranch%basic()
                if (basicBranch%mass() >= self%massThreshold) then
                   pruneBranch=.false.
                   exit
                end if
             end do
             ! Always decouple the node from the tree.
             if (associated(node%parent)) call Merger_Tree_Prune_Unlink_Parent(node,node%parent,parentWillBePruned=.false.,preservePrimaryProgenitor=.false.)
             node%parent  => null()
             node%sibling => null()
             ! Prune or keep the branch, depending on the masses found.
             if (pruneBranch) then
                ! Branch is to be pruned - clean the branch.
                call Merger_Tree_Prune_Clean_Branch(node)
                ! Destroy the branch.
                call node%destroyBranch()
                deallocate(node)
             else
                ! Branch is to be retained.
                ! Add the node to the list of new root nodes.
                if (associated(nodeRootHead)) then
                   allocate(nodeRootCurrent%next)
                   nodeRootCurrent => nodeRootCurrent%next
                else
                   allocate(nodeRootHead)
                   nodeRootCurrent => nodeRootHead
                end if
                nodeRootCurrent%node => node
             end if
          end if
       end if
       node => nodeNext
    end do
    ! Destroy trees with root nodes that are not on the new list of root nodes.
    currentTree => tree
    do while (associated(currentTree))
       keepTree=.false.
       nodeRootCurrent => nodeRootHead
       do while (associated(nodeRootCurrent))
          if (associated(currentTree%nodeBase,nodeRootCurrent%node)) then
             keepTree=.true.
             exit
          end if
          nodeRootCurrent => nodeRootCurrent%next
       end do
       if (.not.keepTree) then
          ! Branch is to be pruned - clean the branch.
          call Merger_Tree_Prune_Clean_Branch(currentTree%nodeBase)
          ! Destroy the branch.
          call currentTree%nodeBase%destroyBranch()
          deallocate(currentTree%nodeBase)
       end if
       currentTree => currentTree%nextTree
    end do
    ! Destroy the current set of trees.
    currentTree => tree
    do while (associated(currentTree))
       treeNext => currentTree%nextTree
       if (.not.associated(currentTree,tree)) deallocate(currentTree)
       currentTree => treeNext
    end do
    ! Create new trees.
    if (associated(nodeRootHead)) then    
       nodeRootCurrent          => nodeRootHead
       tree           %nextTree => null()
       currentTree              => null()
       do while (associated(nodeRootCurrent))
          if (.not.associated(currentTree)) then
             currentTree => tree
          else
             allocate(currentTree%nextTree)
             currentTree                  => currentTree%nextTree
             currentTree%firstTree        => tree
             currentTree%hostUniverse     => tree%hostUniverse
             currentTree%event            => null()
             currentTree%nextTree         => null()
             currentTree%index            =  tree%index
             currentTree%volumeWeight     =  tree%volumeWeight
             currentTree%initializedUntil =  tree%initializedUntil
             call currentTree%properties%initialize()
             allocate(currentTree%randomNumberGenerator_,mold=tree%randomNumberGenerator_)
             !$omp critical(mergerTreeOperatorMassAndTimeeepCopyReset)
             !![
             <deepCopyReset variables="tree%randomNumberGenerator_"/>
             <deepCopy source="tree%randomNumberGenerator_" destination="currentTree%randomNumberGenerator_"/>
             <deepCopyFinalize variables="currentTree%randomNumberGenerator_"/>
             !!]
             !$omp end critical(mergerTreeOperatorMassAndTimeeepCopyReset)
          end if
          currentTree    %nodeBase => nodeRootCurrent%node
          nodeRootNext             => nodeRootCurrent%next
          deallocate(nodeRootCurrent)
          nodeRootCurrent          => nodeRootNext
       end do
       deallocate(nodeRoot)
    else
       ! No tree remains - insert a single node to make an inert tree.
       tree    %nodeBase   => nodeRoot
       tree    %nextTree   => null()
       tree    %event      => null()
       nodeRoot%parent     => null()
       nodeRoot%firstChild => null()
       nodeRoot%sibling    => null()
       nodeRoot%event      => null()
    end if
    ! Reassign host tree pointers.
    currentTree => tree
    do while (associated(currentTree))
       allWalker=mergerTreeWalkerAllNodes(currentTree,spanForest=.false.)
       do while (allWalker%next(node))
          node%hostTree => currentTree
       end do
       currentTree => currentTree%nextTree
    end do    
    return
  end subroutine pruneByMassAndTimeOperatePreEvolution
