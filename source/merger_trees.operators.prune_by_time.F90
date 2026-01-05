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
  Implements a merger tree operator which prunes branches to end at a fixed time.
  !!}

  !![
  <mergerTreeOperator name="mergerTreeOperatorPruneByTime">
   <description>Provides a merger tree operator which prunes branches to end at a fixed time.</description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneByTime
     !!{
     A merger tree operator class which prunes branches to end at a fixed time.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     double precision                                   :: massMinimum                  , massMaximum     , &
          &                                                timeEarliest                 , redshiftEarliest
   contains
     final     ::                        pruneByTimeDestructor
     procedure :: operatePreEvolution => pruneByTimeOperatePreEvolution
  end type mergerTreeOperatorPruneByTime

  interface mergerTreeOperatorPruneByTime
     !!{
     Constructors for the prune-by-time merger tree operator class.
     !!}
     module procedure pruneByTimeConstructorParameters
     module procedure pruneByTimeConstructorInternal
  end interface mergerTreeOperatorPruneByTime

contains

  function pruneByTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the prune-by-time merger tree operator class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass
    implicit none
    type            (mergerTreeOperatorPruneByTime)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    double precision                                               :: massMinimum        , massMaximum     , &
         &                                                            timeEarliest       , redshiftEarliest

    !![
    <inputParameter>
      <name>redshiftEarliest</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>Redshift at which to truncate merger tree branches.</description>
    </inputParameter>
    <inputParameter>
      <name>massMinimum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>Minimum mass for which to consider merger tree branches for truncation.</description>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <source>parameters</source>
      <defaultValue>huge(0.0d0)</defaultValue>
      <description>Maximum mass for which to consider merger tree branches for truncation.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    timeEarliest=cosmologyFunctions_ %cosmicTime                 (              &
         &        cosmologyFunctions_%expansionFactorFromRedshift (             &
         &                                                         timeEarliest &
         &                                                        )             &
         &                                                       )
    self=mergerTreeOperatorPruneByTime(timeEarliest,massMinimum,massMaximum,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function pruneByTimeConstructorParameters

  function pruneByTimeConstructorInternal(timeEarliest,massMinimum,massMaximum,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the prune-by-time merger tree operator class.
    !!}
    implicit none
    type            (mergerTreeOperatorPruneByTime)                        :: self
    class           (cosmologyFunctionsClass      ), intent(in   ), target :: cosmologyFunctions_
    double precision                               , intent(in   )         :: massMinimum        , massMaximum, &
         &                                                                    timeEarliest
    !![
    <constructorAssign variables="timeEarliest, massMinimum, massMaximum, *cosmologyFunctions_"/>
    !!]
    
    self%redshiftEarliest=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeEarliest))
    return
  end function pruneByTimeConstructorInternal

  subroutine pruneByTimeDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeOperatorPruneByTime} merger tree operator class.
    !!}
    implicit none
    type(mergerTreeOperatorPruneByTime), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine pruneByTimeDestructor

  subroutine pruneByTimeOperatePreEvolution(self,tree)
    !!{
    Perform a prune-by-time operation on a merger tree.
    !!}
    use :: Galacticus_Nodes              , only : mergerTree                    , nodeComponentBasic, treeNode
    use :: Merger_Tree_Walkers           , only : mergerTreeWalkerIsolatedNodes
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Clean_Branch
    implicit none
    class           (mergerTreeOperatorPruneByTime), intent(inout), target  :: self
    type            (mergerTree                   ), intent(inout), target  :: tree
    type            (treeNode                     )               , pointer :: node              , nodeChild   , &
         &                                                                     nodeNew           , nodeNext
    class           (nodeComponentBasic           )               , pointer :: basic             , basicParent , &
         &                                                                     basicChild
    type            (mergerTreeWalkerIsolatedNodes)                         :: treeWalker
    double precision                                                        :: massNow           , massParent  , &
         &                                                                     timeNow           , timeParent  , &
         &                                                                     massAtTimeEarliest
    logical                                                                    nodesRemain

    ! Iterate over nodes.
    treeWalker =mergerTreeWalkerIsolatedNodes(tree,spanForest=.true.)
    nodesRemain=treeWalker%next(nodeNext)
    do while (nodesRemain)
       node        =>                 nodeNext
       nodesRemain =  treeWalker%next(nodeNext)
       ! Skip this node if it is the root node.
       if (associated(node%parent)) then
          ! Get basic components.
          basic       => node       %basic()
          basicParent => node%parent%basic()
          ! Get the time of this node and its parent.
          timeNow   =basic      %time()
          timeParent=basicParent%time()
          ! If the branch from node to parent spans the earliest time, insert a new node at that time.
          if (timeParent >= self%timeEarliest .and. timeNow < self%timeEarliest) then
             ! Get masses of these halos.
             massNow   =basic  %mass()
             massParent=basicParent%mass()
             if (node%isPrimaryProgenitor()) then
                ! Remove the mass in any non-primary progenitors - we don't want to include their mass in the estimated mass
                ! growth rate of this node.
                nodeChild => node%parent%firstChild%sibling
                do while (associated(nodeChild))
                   basicChild => nodeChild%basic()
                   massParent =  massParent-basicChild%mass()
                   nodeChild  => nodeChild%sibling
                end do
             else
                ! Halo is not the primary progenitor of its parent. Assume that its mass does not grow further.
                massParent=massNow
             end if
             ! Determine mass at truncation time.
             massAtTimeEarliest=massNow+(massParent-massNow)*(self%timeEarliest-timeNow)/(timeParent-timeNow)
             ! Check if mass is within range.
             if     (                                        &
                  &   massAtTimeEarliest >= self%massMinimum &
                  &  .and.                                   &
                  &   massAtTimeEarliest <= self%massMaximum &
                  & ) then
                ! Decide whether or not to create a new, interpolated node.                
                if (timeParent == self%timeEarliest) then
                   ! Parent exists precisely at the pruning time. No need to create a new node. Simply uncouple the node to be
                   ! pruned from the parent.
                   if (node%isPrimaryProgenitor()) then
                      node%parent%firstChild  => node%sibling
                   else
                      nodeChild => node%parent%firstChild
                      do while (.not.associated(nodeChild%sibling,node))
                         nodeChild => nodeChild%sibling
                      end do
                      nodeChild%sibling => node%sibling
                  end if
                else
                   ! Create new interpolated node at the pruning time.
                   nodeNew => treeNode(hostTree=node%hostTree)
                   call nodeNew%indexSet(node%index())
                   ! Assign a time and a mass
                   basic => nodeNew%basic(autoCreate=.true.)
                   call basic%timeSet(self%timeEarliest )
                   call basic%massSet(massAtTimeEarliest)
                   ! No child node.
                   nodeNew%firstChild => null()
                   ! Link to parent node.
                   nodeNew%parent     => node%parent
                   ! Link  sibling to current node sibling.
                   nodeNew%sibling    => node%sibling
                   ! Link the parent if necessary.
                   if (node%isPrimaryProgenitor()) then
                      ! Node is the main progenitor of its parent, so simply replace it with the final node in our list.
                      node%parent%firstChild  => nodeNew
                   else
                      ! Node is not the main progenitor of its parent, so find the child node that has it as a sibling.
                      nodeChild => node%parent%firstChild
                      do while (.not.associated(nodeChild%sibling,node))
                         nodeChild => nodeChild%sibling
                      end do
                      nodeChild%sibling => nodeNew
                   end if
                end if
                ! Clean the branch.
                call Merger_Tree_Prune_Clean_Branch(node)
                ! Destroy the branch.
                call node%destroyBranch()
                deallocate(node)
             end if
          end if
       end if
    end do
    return
  end subroutine pruneByTimeOperatePreEvolution
