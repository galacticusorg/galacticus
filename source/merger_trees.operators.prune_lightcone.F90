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
  Implements a prune-by-lightcone operator on merger trees.
  !!}
  
  use :: Geometry_Lightcones           , only : geometryLightconeClass
  use :: Satellite_Oprhan_Distributions, only : satelliteOrphanDistributionClass
  use :: Output_Times                  , only : outputTimesClass

  !![
  <mergerTreeOperator name="mergerTreeOperatorPruneLightcone">
   <description>
    Provides a pruning-by-lightcone operator on merger trees. Trees which have no nodes which lie within the lightcone are
    completely pruned away. If the parameter {\normalfont \ttfamily [splitTrees]} is set to {\normalfont \ttfamily true} then
    any parts of a merger tree which does intersect the light that exist after the latest time at which a constituent node of
    the tree intersects the lightcone will be pruned away also (possibly causing the tree to be split into multiple trees in a
    forest). If the parameter {\normalfont \ttfamily [bufferIsolatedHalos]} is set to {\normalfont \ttfamily true} then, when
    testing whether an isolated halo intersects the lightcone a buffer radius equal in size to the extent of any possible
    orphan galaxies associated with the halo is added around the lightcone---this ensures that if orphan galaxies of the halo
    might possibly intersect the lightcone the halo will not be pruned away.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneLightcone
     !!{
     A pruning-by-lightcone merger tree operator class.
     !!}
     private
     class  (geometryLightconeClass          ), pointer :: geometryLightcone_           => null()
     class  (satelliteOrphanDistributionClass), pointer :: satelliteOrphanDistribution_ => null()
     class  (outputTimesClass                ), pointer :: outputTimes_                 => null()
     logical                                            :: bufferIsolatedHalos                   , positionHistoryAvailable, &
          &                                                splitTrees
   contains
     final     ::                        pruneLightconeDestructor
     procedure :: operatePreEvolution => pruneLightconeOperatePreEvolution
  end type mergerTreeOperatorPruneLightcone

  interface mergerTreeOperatorPruneLightcone
     !!{
     Constructors for the pruning-by-lightcone merger tree operator class.
     !!}
     module procedure pruneLightconeConstructorParameters
     module procedure pruneLightconeConstructorInternal
  end interface mergerTreeOperatorPruneLightcone

contains

  function pruneLightconeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the prune-by-lightcone merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (mergerTreeOperatorPruneLightcone)                :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    class  (geometryLightconeClass          ), pointer       :: geometryLightcone_
    class  (satelliteOrphanDistributionClass), pointer       :: satelliteOrphanDistribution_
    class (outputTimesClass                 ), pointer       :: outputTimes_
    logical                                                  :: bufferIsolatedHalos         , splitTrees

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>bufferIsolatedHalos</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, intersection of a tree with the lightcone will be determined using the positions of non-isolated (a.k.a. ``satellite'') halos, and of isolated halos (a.k.a ``centrals'') with a buffer region (with radius equal to the extent of the orphan satellite distribution---see \refPhysics{satelliteOrphanDistribution}) placed around each such halo, and any intersection of that region with the lightcone is sufficient to prevent pruning of the tree. If this parameter is {\normalfont \ttfamily false} then (unbuffered) positions of all halos are used for determining intersection with the lightcone---this requires complete (i.e. throughout the extent of their existence) knowledge of non-isolated halos prior to application of this operator.</description>
    </inputParameter>
    <inputParameter>
      <name>splitTrees</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, prune away any nodes of the tree that are not needed to determine evolution up to the latest time at which a node is present inside the lightcone. This typically leads to a tree splitting into a forest of trees.</description>
    </inputParameter>
    <objectBuilder class="geometryLightcone"           name="geometryLightcone_"           source="parameters"/>
    <objectBuilder class="satelliteOrphanDistribution" name="satelliteOrphanDistribution_" source="parameters"/>
    <objectBuilder class="outputTimes"                 name="outputTimes_"                 source="parameters"/>
    !!]
    self=mergerTreeOperatorPruneLightcone(geometryLightcone_,satelliteOrphanDistribution_,outputTimes_,bufferIsolatedHalos,splitTrees)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="geometryLightcone_"          />
    <objectDestructor name="satelliteOrphanDistribution_"/>
    <objectDestructor name="outputTimes_"                />
    !!]
    return
  end function pruneLightconeConstructorParameters

  function pruneLightconeConstructorInternal(geometryLightcone_,satelliteOrphanDistribution_,outputTimes_,bufferIsolatedHalos,splitTrees) result(self)
    !!{
    Internal constructor for the prune-by-lightcone merger tree operator class.
    !!}
    implicit none
    type   (mergerTreeOperatorPruneLightcone)                        :: self
    class  (geometryLightconeClass          ), intent(in   ), target :: geometryLightcone_
    class  (satelliteOrphanDistributionClass), intent(in   ), target :: satelliteOrphanDistribution_
    class  (outputTimesClass                ), intent(in   ), target :: outputTimes_
    logical                                  , intent(in   )         :: bufferIsolatedHalos         , splitTrees
    !![
    <constructorAssign variables="*geometryLightcone_, *satelliteOrphanDistribution_, *outputTimes_, bufferIsolatedHalos, splitTrees"/>
    !!]

    call pruneLightconeValidate(self)
    return
  end function pruneLightconeConstructorInternal

  subroutine pruneLightconeDestructor(self)
    !!{
    Destructor for the lightcone merger tree operator function class.
    !!}
    implicit none
    type(mergerTreeOperatorPruneLightcone), intent(inout) :: self

    !![
    <objectDestructor name="self%geometryLightcone_"          />
    <objectDestructor name="self%satelliteOrphanDistribution_"/>
    <objectDestructor name="self%outputTimes_"                />
    !!]
    return
  end subroutine pruneLightconeDestructor

  subroutine pruneLightconeValidate(self)
    !!{
    Validate the lightcone pruning operator.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Component_List          , Error_Report
    use :: Galacticus_Nodes, only : defaultPositionComponent, defaultSatelliteComponent
    implicit none
    class(mergerTreeOperatorPruneLightcone), intent(inout) :: self

    ! If buffering is applied to isolated halos, then satellite halos are to be checked for intersection only for as long as their
    ! position is known. In this case, check if satellite position history can be obtained, and also require that satellite time
    ! of merging can be both read and written.
    if (self%bufferIsolatedHalos) then
       self%positionHistoryAvailable=defaultPositionComponent%positionHistoryIsGettable()
       if     (                                                                                                                               &
            &  .not.(                                                                                                                         &
            &         defaultSatelliteComponent%timeOfMergingIsGettable()                                                                     &
            &        .and.                                                                                                                    &
            &         defaultSatelliteComponent%timeOfMergingIsSettable()                                                                     &
            &       )                                                                                                                         &
            & )                                                                                                                               &
            & call Error_Report                                                                                                               &
            &      (                                                                                                                          &
            &       'buffering isolated halos requires that the timeOfMerging property of the satellite component be gettable and settable'// &
            &       Component_List(                                                                                                           &
            &                      'satellite'                                                                                             ,  &
            &                        defaultSatelliteComponent%timeOfMergingAttributeMatch(requireGettable=.true.)                            &
            &                       .intersection.                                                                                            &
            &                        defaultSatelliteComponent%timeOfMergingAttributeMatch(requireSettable=.true.)                            &
            &                      )                                                                                                       // &
            &       {introspection:location}                                                                                                  &
            &      )
    end if
    return
  end subroutine pruneLightconeValidate

  subroutine pruneLightconeOperatePreEvolution(self,tree)
    !!{
    Perform a prune-by-lightcone operation on a merger tree.
    !!}
    use :: Display                       , only : displayMessage                , displayIndent           , displayUnindent      , verbosityLevelInfo
    use :: Galacticus_Nodes              , only : mergerTree                    , nodeComponentBasic      , nodeComponentPosition, nodeComponentSatellite, &
         &                                        treeNode                      , treeNodeLinkedList
    use :: Histories                     , only : history
    use :: Merger_Tree_Walkers           , only : mergerTreeWalkerIsolatedNodes , mergerTreeWalkerAllNodes
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Clean_Branch
    use :: Satellite_Merging_Timescales  , only : satelliteMergeTimeInfinite
    implicit none
    class           (mergerTreeOperatorPruneLightcone), intent(inout), target :: self
    type            (mergerTree                      ), intent(inout), target :: tree
    type            (treeNode                        ), pointer               :: node                   , nodeNext               , &
         &                                                                       nodeParent             , nodeChild
    type            (mergerTree                      ), pointer               :: treeCurrent            , treeNext
    type            (treeNodeLinkedList              ), pointer               :: forestRootsOriginalHead, forestRootsNewHead     , &
         &                                                                       forestRootsOriginalLast, forestRootsNewLast
    class           (nodeComponentBasic              ), pointer               :: basic                  , basicChild
    class           (nodeComponentPosition           ), pointer               :: position
    class           (nodeComponentSatellite          ), pointer               :: satellite
    type            (history                         )                        :: positionHistory
    type            (mergerTreeWalkerIsolatedNodes   )                        :: treeWalker
    type            (mergerTreeWalkerAllNodes        )                        :: treeWalkerAll
    double precision                                                          :: radiusBuffer           , timeOfMergingCurrent   , &
         &                                                                       timeInLightconeLatest
    logical                                                                   :: nodeIntersectsLightcone, treeIntersectsLightcone, &
         &                                                                       nodeIsNewRoot
    character       (len=16                          )                        :: labelIndex             , labelTime

    write (labelIndex,'(i16)') tree%index
    call displayIndent('Pruning tree ['//trim(adjustl(labelIndex))//'] for lightcone geometry',verbosityLevelInfo)
    ! Set buffer size to zero by default.
    radiusBuffer            =        0.0d0
    ! Initialize the latest time in lightcone to an unphysical value.
    timeInLightconeLatest   =  -huge(0.0d0)
    ! Initialize the state of the tree-in-lightcone condition.
    treeIntersectsLightcone =  .false.
    ! Initialize linked lists of forest root nodes.
    forestRootsNewHead      => null()
    forestRootsOriginalHead => null()
    forestRootsNewLast      => null()
    forestRootsOriginalLast => null()
    ! Iterate over trees.
    nodeParent    => null()
    treeWalkerAll =  mergerTreeWalkerAllNodes(tree,spanForest=.true.)
    do while (treeWalkerAll%next(node))
       ! If buffering isolated halos, set a suitable buffer radius. For isolated halos this is the maximum extent of their
       ! orphan satellite population. For non-isolated halos, no buffer is required. In addition, for non-isolated halos, limit
       ! the time over which they will be checked for intersection with the lightcone to the maximum time for which they have a
       ! defined position.
       if (self%bufferIsolatedHalos) then
          ! Store current parent of the node.
          nodeParent => node%parent
          ! Handle satellites and centrals.
          if (node%isSatellite()) then
             ! No buffer is required for a satellite node.
             radiusBuffer=0.0d0
             ! Find the extent of the position history known for this node. Limit the merging time to the final time for which
             ! position is known (or the current time if no position history is available). Record the current time of merging
             ! so that it can be reset after testing for intersection.
             position             => node    %position       ()
             basic                => node    %basic          ()
             satellite            => node    %satellite      ()
             positionHistory      =  position%positionHistory()
             timeOfMergingCurrent =  satellite%timeOfMerging ()
             if (positionHistory%exists()) then
                call satellite%timeOfMergingSet(positionHistory%time(size(positionHistory%time)))
             else
                call satellite%timeOfMergingSet(basic          %time(                          ))
             end if
          else
             ! For a central - set a suitable buffer.
             radiusBuffer=self%satelliteOrphanDistribution_%extent(node)
          end if
       end if
       ! Test for intersection with the lightcone.
       nodeIntersectsLightcone=self%geometryLightcone_%isInLightcone(node,atPresentEpoch=.false.,radiusBuffer=radiusBuffer)
       ! Reset the time of merging if it was adjusted above.
       if (self%bufferIsolatedHalos) then
          ! Reset time of merging.
          if (node%isSatellite()) call satellite%timeOfMergingSet(timeOfMergingCurrent)
       end if
       ! If intersection with lightcone was detected then this tree can not be pruned - return immediately.
       if (nodeIntersectsLightcone) then
          if (self%splitTrees) then
             ! Trees are to be split - record the latest time at which a node is in the lightcone.
             if (associated(node%parent)) then
                basic => node%parent%basic()
             else
                basic => node       %basic()
             end if
             timeInLightconeLatest  =max(timeInLightconeLatest,basic%time())
             treeIntersectsLightcone=.true.
          else
             ! Trees are not being split - this tree is in the lightcone, so can not be pruned - we can simply return.
             call displayUnindent('done',verbosityLevelInfo)
             return
          end if
       end if
    end do
    if (.not.treeIntersectsLightcone) then
       ! Entire forest is outside lightcone. Destroy all but the base node in the first tree. (Leaving just the base node makes the
       ! tree inert - i.e. it can not do anything.) Destroy any additional trees in the forest.
       node => tree%nodeBase%firstChild
       do while (associated(node))
          nodeNext => node%sibling
          call Merger_Tree_Prune_Clean_Branch(node)
          call node%destroyBranch()
          deallocate(node)
          node => nodeNext
       end do
       treeCurrent          => tree%nextTree
       tree       %nextTree => null()
       do while (associated(treeCurrent))
          treeNext => treeCurrent%nextTree
          ! Destroy the tree.
          call treeCurrent%destroy()
          deallocate(treeCurrent)
          ! Move to the next tree.
          treeCurrent => treeNext
       end do
    else if (self%splitTrees) then
       ! Tree is in lightcone, but we are asked to split it into forests and trim off late-time nodes where possible.
       !! Increase the latest time in the lightcone to the subsequent output time to ensure that the tree will be evolved beyond
       !! this time - necessary in case the actual lightcone crossing occurs slightly after the time at which the node presently
       !! exists.
       if (self%outputTimes_%timeNext(timeInLightconeLatest) > 0.0d0) timeInLightconeLatest=self%outputTimes_%timeNext(timeInLightconeLatest)
       !! Walk the trees, building a list of new root nodes and any original root nodes that we can remove.
       write (labelTime,'(f9.3)') timeInLightconeLatest
       call displayMessage('Splitting tree into forests at time '//trim(adjustl(labelTime))//' Gyr',verbosityLevelInfo)
       treeWalker=mergerTreeWalkerIsolatedNodes(tree,spanForest=.true.)
       do while (treeWalker%next(node))
          nodeIsNewRoot=.false.
          basic => node%basic()
          if (basic%time() == timeInLightconeLatest) then
             ! Node exists precisely at the latest time that we require in the tree. This node becomes a new root node for our new
             ! forest.
             nodeIsNewRoot=.true.
          else if (basic%time() > timeInLightconeLatest) then
             nodeChild => node%firstChild
             do while (associated(nodeChild))                
                basicChild => nodeChild%basic()
                if (basicChild%time() < timeInLightconeLatest) then
                   ! Node exists after the latest time that we need, but has a child which exists before that time. This node
                   ! becomes a new root node for our new forest.
                   nodeIsNewRoot=.true.
                   exit
                end if
                nodeChild => nodeChild%sibling
             end do
          end if
          ! Add this node to a list of root nodes if necessary.
          if (nodeIsNewRoot .or. (.not.associated(node%parent) .and. basic%time() <= timeInLightconeLatest)) then
             ! This is a new root node (created by splitting the tree at the latest time needed for the lightcone), or an original
             ! root node which exists prior to the latest time needed for the lightcone (in which case it must be kept as a root
             ! node as its evolution may affect other trees in the forest).
             if (.not.associated(forestRootsNewHead)) then
                allocate(forestRootsNewHead)
                forestRootsNewHead%node      => node
                forestRootsNewHead%next      => null()
                forestRootsNewLast           => forestRootsNewHead
             else
                allocate(forestRootsNewLast%next)
                forestRootsNewLast           => forestRootsNewLast%next
                forestRootsNewLast%node      => node
                forestRootsNewLast%next      => null()
             end if
          else if (               .not.associated(node%parent) .and. basic%time() >  timeInLightconeLatest ) then
             ! This is a root node in the original forest, and exists after the latest time needed for the lightcone - it can be
             ! removed.             
             if (.not.associated(forestRootsOriginalHead)) then
                allocate(forestRootsOriginalHead)
                forestRootsOriginalHead%node => node
                forestRootsOriginalHead%next => null()
                forestRootsOriginalLast      => forestRootsOriginalHead
             else
                allocate(forestRootsOriginalLast%next)
                forestRootsOriginalLast      => forestRootsOriginalLast%next
                forestRootsOriginalLast%node => node
                forestRootsOriginalLast%next => null()
             end if
          end if
       end do
       ! Remove any merge targets which no longer exist.
       treeWalkerAll=mergerTreeWalkerAllNodes(tree,spanForest=.true.)
       do while (treeWalkerAll%next(node))
          if (.not.associated(node%mergeTarget)) cycle
          satellite => node%satellite()
          if (satellite%timeOfMerging() > timeInLightconeLatest) then
             call node%removeFromMergee()
             node%mergeTarget => null()
             call satellite%timeOfMergingSet(satelliteMergeTimeInfinite)
          end if
       end do
       ! Decouple all new root nodes from their parents and siblings. Also decouple their parent from them (so that we can later
       ! destroy the tree containing the parent without destroy this node too).
       forestRootsNewLast => forestRootsNewHead
       do while (associated(forestRootsNewLast))
          node => forestRootsNewLast%node
          if (associated(node%parent)) then
             if (associated(node%parent%firstChild,node)) then
                node%parent%firstChild => node%sibling
             else
                nodeChild => node%parent%firstChild
                do while (.not.associated(nodeChild%sibling,node))
                   nodeChild => nodeChild%sibling
                end do
                nodeChild%sibling => node%sibling   
             end if
             node%parent  => null()
             node%sibling => null()
          end if
          forestRootsNewLast => forestRootsNewLast%next
       end do
       ! For any of the original root nodes that are not also a new root node, delete their associated tree.
       forestRootsOriginalLast => forestRootsOriginalHead
       do while (associated(forestRootsOriginalLast))
          node                    => forestRootsOriginalLast%node
          basic                   => node                   %basic()
          write (labelIndex,'(i16)' ) node %index()
          write (labelTime ,'(f9.3)') basic%time ()
          call displayMessage('Removing branch from ['//trim(adjustl(labelIndex))//'] at '//trim(adjustl(labelTime))//' Gyr',verbosityLevelInfo)
          call node%destroyBranch()
          deallocate(node)
          forestRootsOriginalLast => forestRootsOriginalLast%next
       end do
       ! Replace the first tree in the forest.
       tree %nodeBase => forestRootsNewHead%node
       basic          => forestRootsNewHead%node%basic()
       write (labelIndex,'(i16)' ) forestRootsNewHead%node%index()
       write (labelTime ,'(f9.3)') basic                  %time ()
       call displayMessage('Create new tree from ['//trim(adjustl(labelIndex))//'] at '//trim(adjustl(labelTime))//' Gyr',verbosityLevelInfo)
       ! Destroy all other trees in the current forest.
       treeCurrent => tree%nextTree
       do while (associated(treeCurrent))
          treeNext             => treeCurrent%nextTree
          treeCurrent%nodeBase => null()
          call treeCurrent%destroy()
          deallocate(treeCurrent)
          treeCurrent          => treeNext
       end do
       tree%nextTree => null()
       ! Create new trees in the forest as needed.
       treeCurrent        => tree
       forestRootsNewLast => forestRootsNewHead%next
       do while (associated(forestRootsNewLast))
          allocate(treeCurrent%nextTree)
          treeCurrent%nextTree%nextTree         => null()
          treeCurrent%nextTree%event            => null()
          treeCurrent%nextTree%nodeBase         => forestRootsNewLast                  %node
          treeCurrent%nextTree%index            =  treeCurrent       %nextTree%nodeBase%index       ()
          treeCurrent%nextTree%firstTree        => tree
          treeCurrent%nextTree%volumeWeight     =  tree                                %volumeWeight
          treeCurrent%nextTree%initializedUntil =  0.0d0
          allocate(treeCurrent%nextTree%randomNumberGenerator_,mold=tree%randomNumberGenerator_)
          !$omp critical(mergerTreeOperatorLightconeDeepCopyReset)
          !![
          <deepCopyReset variables="tree%randomNumberGenerator_"/>
          <deepCopy source="tree%randomNumberGenerator_" destination="treeCurrent%nextTree%randomNumberGenerator_"/>
          <deepCopyFinalize variables="treeCurrent%nextTree%randomNumberGenerator_"/>
          !!]
          !$omp end critical(mergerTreeOperatorLightconeDeepCopyReset)
          call treeCurrent%nextTree%randomNumberGenerator_%seedSet(seed=treeCurrent%nextTree%index,offset=.true.)
          basic => forestRootsNewLast%node%basic()
          write (labelIndex,'(i16)' ) treeCurrent%nextTree%nodeBase%index()
          write (labelTime ,'(f9.3)') basic                        %time ()
          call displayMessage('Create new tree from ['//trim(adjustl(labelIndex))//'] at '//trim(adjustl(labelTime))//' Gyr',verbosityLevelInfo)
          treeCurrent        => treeCurrent       %nextTree          
          forestRootsNewLast => forestRootsNewLast%next
       end do
       ! Reassign host tree pointers in all nodes.
       treeCurrent => tree
       do while (associated(treeCurrent))
          treeWalkerAll=mergerTreeWalkerAllNodes(treeCurrent,spanForest=.false.)
          do while (treeWalkerAll%next(node))
             node%hostTree => treeCurrent
          end do
          treeCurrent => treeCurrent%nextTree
       end do
       ! Clean up forest root node lists.
       forestRootsNewLast => forestRootsNewHead
       do while (associated(forestRootsNewLast))
          forestRootsNewHead => forestRootsNewLast
          forestRootsNewLast => forestRootsNewLast%next
          deallocate(forestRootsNewHead)
       end do
       forestRootsOriginalLast => forestRootsOriginalHead
       do while (associated(forestRootsOriginalLast))
          forestRootsOriginalHead => forestRootsOriginalLast
          forestRootsOriginalLast => forestRootsOriginalLast%next
          deallocate(forestRootsOriginalHead)
       end do
    else
       ! This should not happen.
       call Error_Report('unknown state'//{introspection:location})
    end if
    call displayUnindent('done',verbosityLevelInfo)
    return
  end subroutine pruneLightconeOperatePreEvolution
