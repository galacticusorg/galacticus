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
  Implements an augmenting operator on merger trees.
  !!}
  use            :: Cosmology_Functions  , only : cosmologyFunctionsClass
  use, intrinsic :: ISO_C_Binding        , only : c_size_t
  use            :: Merger_Trees_Builders, only : mergerTreeBuilderClass

  !![
  <mergerTreeOperator name="mergerTreeOperatorAugment">
   <description>Provides a merger tree operator which augments tree resolution by inserting high-resolution branches.</description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorAugment
     !!{
     An augmenting merger tree operator class.
     !!}
     private
     double precision                         , allocatable, dimension(   :) :: timeSnapshots
     integer         (c_size_t               )             , dimension(0:10) :: retryHistogram                        , trialCount
     double precision                                                        :: massCutOff                            , timeEarliest             , &
          &                                                                     toleranceScale                        , massCutOffScaleFactor    , &
          &                                                                     massOvershootScaleFactor
     integer                                                                 :: retryMaximum                          , rescaleMaximum           , &
          &                                                                     attemptsMaximum                       , massCutOffAttemptsMaximum, &
          &                                                                     massOvershootAttemptsMaximum
     logical                                                                 :: performChecks                         , useOneNodeTrees
     class           (mergerTreeBuilderClass ), pointer                      :: mergerTreeBuilder_           => null()
     class           (cosmologyFunctionsClass), pointer                      :: cosmologyFunctions_          => null()
   contains
     !![
     <methods>
       <method description="Build a merger tree starting from the given node." method="buildTreeFromNode" />
       <method description="Determine if a newly built tree is an acceptable match." method="acceptTree" />
       <method description="Graft new branches onto all end-nodes of a newly built tree." method="extendNonOverlapNodes" />
       <method description="Sort child nodes into descending mass order." method="sortChildren" />
       <method description="Reinsert a linked list of non-overlap nodes into their parent tree." method="nonOverlapReinsert" />
     </methods>
     !!]
     final     ::                          augmentDestructor
     procedure :: operatePreEvolution   => augmentOperatePreEvolution
     procedure :: finalize              => augmentFinalize
     procedure :: buildTreeFromNode     => augmentBuildTreeFromNode
     procedure :: acceptTree            => augmentAcceptTree
     procedure :: extendNonOverlapNodes => augmentExtendNonOverlapNodes
     procedure :: sortChildren          => augmentSortChildren
     procedure :: nonOverlapReinsert    => augmentNonOverlapReinsert
  end type mergerTreeOperatorAugment

  interface mergerTreeOperatorAugment
     !!{
     Constructors for the \refClass{mergerTreeOperatorAugment} merger tree operator class.
     !!}
     module procedure augmentConstructorParameters
     module procedure augmentConstructorInternal
  end interface mergerTreeOperatorAugment

  ! A tree index counter used since the tree index is used in initializing the random number
  ! sequence for each tree.
  integer :: indexNewTree=-1

  !![
  <enumeration>
   <name>treeStatistic</name>
   <description>Enumeration of tasks to be performed during a tree walk.</description>
   <entry label="nodeCount"   />
   <entry label="endNodeCount"/>
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>treeBuild</name>
   <description>Enumeration of tree building status.</description>
   <entry label="success"         />
   <entry label="failureTolerance"/>
   <entry label="failureStructure"/>
   <entry label="failureGeneric"  />
  </enumeration>
  !!]

contains

  function augmentConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeOperatorAugment} merger tree operator class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeOperatorAugment)                              :: self
    type            (inputParameters          ), intent(inout)               :: parameters
    double precision                           , allocatable  , dimension(:) :: timeSnapshots
    class           (cosmologyFunctionsClass  ), pointer                     :: cosmologyFunctions_
    class           (mergerTreeBuilderClass   ), pointer                     :: mergerTreeBuilder_
    integer                                                                  :: i
    double precision                                                         :: massCutOff                  , toleranceScale           , &
         &                                                                      massCutOffScaleFactor       , massOvershootScaleFactor
    integer                                                                  :: retryMaximum                , rescaleMaximum           , &
         &                                                                      attemptsMaximum             , massCutOffAttemptsMaximum, &
         &                                                                      massOvershootAttemptsMaximum
    logical                                                                  :: performChecks               , useOneNodeTrees

    !![
    <objectBuilder class="mergerTreeBuilder"  name="mergerTreeBuilder_"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <inputParameter>
      <name>massCutOff</name>
      <source>parameters</source>
      <defaultValue>1.0d10</defaultValue>
      <description>For the {\normalfont \ttfamily augment} operator a description of resolution limit for new trees.</description>
    </inputParameter>
    <inputParameter>
      <name>performChecks</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, perform checks of the augmentation process.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceScale</name>
      <source>parameters</source>
      <defaultValue>0.15d0</defaultValue>
      <description>The tolerance scale used in deciding if a trial tree is an acceptable match.</description>
    </inputParameter>
    <inputParameter>
      <name>retryMaximum</name>
      <source>parameters</source>
      <defaultValue>50</defaultValue>
      <description>The number of tree build attempts to complete before rescaling the tolerance.</description>
    </inputParameter>
    <inputParameter>
      <name>rescaleMaximum</name>
      <source>parameters</source>
      <defaultValue>20</defaultValue>
      <description>The maximum allowed number of tolerance rescalings.</description>
    </inputParameter>
    <inputParameter>
      <name>attemptsMaximum</name>
      <source>parameters</source>
      <defaultValue>10000</defaultValue>
      <description>The maximum allowed number of tree build attempts.</description>
    </inputParameter>
    <inputParameter>
      <name>massCutOffAttemptsMaximum</name>
      <source>parameters</source>
      <defaultValue>50</defaultValue>
      <description>The number of trees with nodes above the mass resolution to allow before adjusting the mass cut-off tolerance.</description>
    </inputParameter>
    <inputParameter>
      <name>massCutOffScaleFactor</name>
      <source>parameters</source>
      <defaultValue>0.05d0</defaultValue>
      <description>The amount by which to increase the mass cut-off scale tolerance after exhausting tree build attempts.</description>
    </inputParameter>
    <inputParameter>
      <name>massOvershootAttemptsMaximum</name>
      <source>parameters</source>
      <defaultValue>50</defaultValue>
      <description>The number of failed trees to allow before increasing the mass of the parent node.</description>
    </inputParameter>
    <inputParameter>
      <name>massOvershootScaleFactor</name>
      <source>parameters</source>
      <defaultValue>0.05d0</defaultValue>
      <description>The amount by which to increase the mass overshoot factor after exhausting tree build attempts.</description>
    </inputParameter>
    <inputParameter>
      <name>useOneNodeTrees</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, trees only consisting of their base node will be augmented.</description>
    </inputParameter>
    !!]
    if (parameters%isPresent('snapshotRedshifts')) then
       allocate(timeSnapshots(parameters%count('snapshotRedshifts')))
       !![
       <inputParameter>
         <name>snapshotRedshifts</name>
         <variable>timeSnapshots</variable>
         <source>parameters</source>
         <description>For {\normalfont \ttfamily augment} description of redshift snapshots.</description>
       </inputParameter>
       !!]
       do i=1,size(timeSnapshots)
          timeSnapshots(i)=cosmologyFunctions_ %cosmicTime                 (                  &
               &            cosmologyFunctions_%expansionFactorFromRedshift (                 &
               &                                                             timeSnapshots(i) &
               &                                                            )                 &
               &                                                           )
       end do
    end if
    !![
    <conditionalCall>
     <call>self=mergerTreeOperatorAugment(massCutOff,performChecks,toleranceScale,retryMaximum,rescaleMaximum,attemptsMaximum,massCutOffAttemptsMaximum,massCutOffScaleFactor,massOvershootAttemptsMaximum,massOvershootScaleFactor,useOneNodeTrees,mergerTreeBuilder_,cosmologyFunctions_{conditions})</call>
     <argument name="timeSnapshots" value="timeSnapshots" condition="parameters%isPresent('snapshotRedshifts')"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeBuilder_" />
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function augmentConstructorParameters

  function augmentConstructorInternal(massCutOff,performChecks,toleranceScale,retryMaximum,rescaleMaximum,attemptsMaximum,massCutOffAttemptsMaximum,massCutOffScaleFactor,massOvershootAttemptsMaximum,massOvershootScaleFactor,useOneNodeTrees,mergerTreeBuilder_,cosmologyFunctions_,timeSnapshots) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeOperatorAugment} merger tree operator class.
    !!}
    use :: Cosmology_Functions, only : cosmologyFunctionsClass
    use :: Sorting            , only : sort
    implicit none
    type            (mergerTreeOperatorAugment)                                        :: self
    double precision                           , intent(in   )                         :: massCutOff                       , toleranceScale           , &
         &                                                                                massCutOffScaleFactor            , massOvershootScaleFactor
    integer                                    , intent(in   )                         :: retryMaximum                     , rescaleMaximum           , &
         &                                                                                attemptsMaximum                  , massCutOffAttemptsMaximum, &
         &                                                                                massOvershootAttemptsMaximum
    double precision                           , intent(in   ), dimension(:), optional :: timeSnapshots
    class           (mergerTreeBuilderClass   ), intent(in   ), pointer                :: mergerTreeBuilder_
    class           (cosmologyFunctionsClass  ), intent(in   ), pointer                :: cosmologyFunctions_
    logical                                    , intent(in   )                         :: performChecks                    , useOneNodeTrees
    double precision                           , parameter                             :: expansionFactorDefault    =0.01d0
    !![
    <constructorAssign variables="massCutOff,performChecks,toleranceScale,retryMaximum,rescaleMaximum,attemptsMaximum,massCutOffAttemptsMaximum,massCutOffScaleFactor,massOvershootAttemptsMaximum,massOvershootScaleFactor,useOneNodeTrees,*mergerTreeBuilder_,*cosmologyFunctions_"/>
    !!]

    if (present(timeSnapshots)) then
       allocate(self%timeSnapshots,mold=timeSnapshots)
       self%timeSnapshots=timeSnapshots
       call sort(self%timeSnapshots)
       self%timeEarliest=min(                                                           &
            &                cosmologyFunctions_%cosmicTime   (expansionFactorDefault), &
            &                self               %timeSnapshots(                     1)  &
            &               )
    else
       self%timeEarliest=    cosmologyFunctions_%cosmicTime   (expansionFactorDefault)
    end if
    self%retryHistogram=0
    self%trialCount    =0
    call self%mergerTreeBuilder_%timeEarliestSet(self%timeEarliest)
    return
  end function augmentConstructorInternal

  subroutine augmentDestructor(self)
    !!{
    Destructor for the augment merger tree operator function class.
    !!}
    implicit none
    type(mergerTreeOperatorAugment), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBuilder_" />
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine augmentDestructor

  subroutine augmentOperatePreEvolution(self,tree)
    !!{
    Augment the resolution of a merger tree by inserting high resolution branches.
    !!}
    use            :: Display            , only : displayIndent                , displayMessage    , displayUnindent, displayVerbosity, &
          &                                       verbosityLevelWorking
    use            :: Galacticus_Nodes   , only : mergerTree                   , nodeComponentBasic, treeNode       , treeNodeList
    use, intrinsic :: ISO_C_Binding      , only : c_size_t
    use            :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    use            :: Sorting            , only : sortIndex
    use            :: String_Handling    , only : operator(//)
    implicit none
    class           (mergerTreeOperatorAugment    ), intent(inout), target       :: self
    type            (mergerTree                   ), intent(inout), target       :: tree
    type            (treeNode                     ), pointer                     :: node
    class           (nodeComponentBasic           ), pointer                     :: basic
    type            (mergerTree                   ), pointer                     :: treeCurrent
    type            (treeNodeList                 ), allocatable  , dimension(:) :: anchorNodes
    double precision                               , allocatable  , dimension(:) :: anchorTimes
    integer         (c_size_t                     ), allocatable  , dimension(:) :: anchorIndex
    type            (varying_string               )                              :: message
    type            (mergerTree                   )                              :: treeBest
    type            (mergerTreeWalkerIsolatedNodes)                              :: treeWalker
    type            (enumerationTreeBuildType     )                              :: treeBuilt
    integer                                                                      :: nodeCount                     , i                             , &
         &                                                                          attemptsRemaining             , rescaleCount                  , &
         &                                                                          massCutoffAttemptsRemaining   , massOvershootAttemptsRemaining, &
         &                                                                          retryCount
    double precision                                                             :: tolerance                     , treeBestWorstFit              , &
         &                                                                          massCutoffScale               , massOvershootScale
    logical                                                                      :: treeBestOverride              , treeNewHasNodeAboveResolution , &
         &                                                                          treeBestHasNodeAboveResolution, newRescale                    , &
         &                                                                          nodeBranches

    ! Iterate over all linked trees in this forest.
    call displayIndent('Augmenting merger tree',verbosityLevelWorking)
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Allocate array of original anchor nodes from which new high-resolution branches will be built.
       nodeCount=augmentTreeStatistics(treeCurrent,treeStatisticNodeCount)
       if (nodeCount > 1 .or. self%useOneNodeTrees) then
          if (displayVerbosity() >= verbosityLevelWorking) then
             message="Number of nodes in tree: "
             message=message//nodeCount
             call displayMessage(message)
          end if
          allocate(anchorNodes(nodeCount))
          allocate(anchorTimes(nodeCount))
          allocate(anchorIndex(nodeCount))
          ! Build pointers to all anchor nodes.
          i         =0
          treeWalker=mergerTreeWalkerIsolatedNodes(treeCurrent)
          do while(treeWalker%next(node))
             i                   =  i               +1
             basic               => node %basic   ()
             anchorNodes(i)%node => node
             anchorTimes(i)      =  basic%time    ()
          end do
          anchorIndex=sortIndex(anchorTimes)
          deallocate(anchorTimes)
          ! Walk the tree.
          i=1
          do while (i <= nodeCount)
             ! Get the node to work with.
             node                           => anchorNodes(anchorIndex(nodeCount-i+1))%node
             basic                          => node                                   %basic()
             ! Initialize the current best-known tree to null.
             treeBestWorstFit               =      3.000d0
             treeBest%nodeBase              => null()
             treeBestOverride               = .false.
             treeBestHasNodeAboveResolution = .false.
             newRescale                     = .false.
             ! Reset all factors used in tree acceptance.
             tolerance                      =  self%toleranceScale
             rescaleCount                   =      0
             retryCount                     =      1
             treeBuilt                      =  treeBuildFailureGeneric
             attemptsRemaining              =  self%attemptsMaximum
             massCutoffScale                =      1.0d0
             massOvershootScale             =      1.0d0
             massCutoffAttemptsRemaining    =  self%massCutOffAttemptsMaximum
             massOvershootAttemptsRemaining =  self%massOvershootAttemptsMaximum
             ! Determine if the node branches.
             nodeBranches=associated(node%firstChild).and.associated(node%firstChild%sibling)
             ! Begin building trees from this node, searching for an acceptable tree.
             if (displayVerbosity() >= verbosityLevelWorking) then
                message="Building tree from node: "
                message=message//node%index()
                call displayIndent(message)
             end if
             do while (                                            &
                  &     treeBuilt         /= treeBuildSuccess      & ! Exit if tree successfully built.
                  &    .and.                                       &
                  &     rescaleCount      <= self%rescaleMaximum   & ! Exit once number of rescalings exceeds maximum allowed.
                  &    .and.                                       &
                  &     attemptsRemaining >  0                     & ! Exit if no more attempts remain.
                  &   )
                treeNewHasNodeAboveResolution=.false.
                treeBuilt                    =self%buildTreeFromNode(                                                 &
                     &                                               node                                           , &
                     &                                               .false.                                        , &
                     &                                               tolerance                                      , &
                     &                                               self%timeEarliest                              , &
                     &                                               treeBest                                       , &
                     &                                               treeBestWorstFit                               , &
                     &                                               treeBestOverride                               , &
                     &                                               massCutoffScale                                , &
                     &                                               massOvershootScale                             , &
                     &                                               treeNewHasNodeAboveResolution                  , &
                     &                                               treeBestHasNodeAboveResolution                 , &
                     &                                               newRescale                                       &
                     &                                              )
                ! Check for exhaustion of retry attempts.
                if (retryCount == self%retryMaximum) then
                   ! Rescale the tolerance to allow less accurate tree matches to be accepted in future.
                   retryCount  =0
                   tolerance   =tolerance   *(1.0d0+self%toleranceScale)
                   ! Check for exhaustion of rescaling attempts.
                   rescaleCount=rescaleCount+1
                   if (rescaleCount > self%rescaleMaximum) call displayMessage('Node build attempts exhausted',verbosityLevelWorking)
                   newRescale = .true.
                else
                   newRescale = .false.
                end if
                ! Increment the retry count in cases where the tree was not accepted due to matching tolerance.
                select case (treeBuilt%ID)
                case (treeBuildSuccess         %ID)
                   retryCount=retryCount-1
                case (treeBuildFailureTolerance%ID)
                   retryCount=retryCount+1
                end select
                ! Decrement the number of attempts remaining and, if the best tree is to be forcibly used, set no attempts remaining.
                attemptsRemaining                                                          =attemptsRemaining-1
                if (treeBestOverride .and. treeBuilt /= treeBuildSuccess) attemptsRemaining=attemptsRemaining+1
                ! If all attempts have been used but no match has been found, insert the best tree on next pass through loop.
                if (attemptsRemaining == 1 .and. associated(treeBest%nodeBase)) treeBestOverride=.true.
                ! If the new tree contains nodes above the mass cut-off, decrement the number of remaining retries.
                if (treeNewHasNodeAboveResolution) then
                   massCutoffAttemptsRemaining=massCutoffAttemptsRemaining-1
                   ! If number of retries is exhausted, adjust the tolerance for declaring nodes to be above the mass cut-off.
                   if (massCutoffAttemptsRemaining == 0) then
                      massCutoffAttemptsRemaining=                self%massCutoffAttemptsMaximum
                      massCutoffScale            =massCutoffScale+self%massCutOffScaleFactor
                   end if
                end if
                ! Increment the mass overshoot budget if attempts are exhausted.
                massOvershootAttemptsRemaining=massOvershootAttemptsRemaining-1
                ! If number of retries is exhausted, adjust the tolerance for declaring nodes to be above the mass cut-off.
                if (massOvershootAttemptsRemaining == 0) then
                   massOvershootAttemptsRemaining=                   self%massOvershootAttemptsMaximum
                   massOvershootScale            =massOvershootScale+self%massOvershootScaleFactor
                end if
             end do
             ! Accumulate the histogram of rescalings.
             if (nodeBranches) then
                self        %retryHistogram(min(rescaleCount,ubound(self%retryHistogram)))= &
                     & +self%retryHistogram(min(rescaleCount,ubound(self%retryHistogram)))  &
                     & +1
                self        %trialCount    (min(rescaleCount,ubound(self%retryHistogram)))= &
                     & +self%trialCount    (min(rescaleCount,ubound(self%retryHistogram)))  &
                     & +max(1,1+retryCount)
             end if
             ! Clean up the best tree if one exists.
             if (associated(treeBest%nodeBase)) then
                call treeBest%nodeBase%destroyBranch()
                deallocate(treeBest%nodeBase)
                treeBest%nodeBase => null()
             end if
             ! Move on to the next node.
             i=i+1
             call displayUnindent('Finished building tree',verbosityLevelWorking)
          end do
          deallocate(anchorNodes)
          deallocate(anchorIndex)
       end if
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    call displayUnindent('done',verbosityLevelWorking)
    return
  end subroutine augmentOperatePreEvolution

  recursive function augmentBuildTreeFromNode(self,node,extendingEndNode,tolerance,timeEarliestIn,treeBest,treeBestWorstFit,treeBestOverride,massCutoffScale,massOvershootScale,treeNewHasNodeAboveResolution,treeBestHasNodeAboveResolution,newRescale)
    use            :: Arrays_Search       , only : searchArrayClosest
    use            :: Error               , only : Error_Report      , errorStatusSuccess
    use            :: Galacticus_Nodes    , only : mergerTree        , nodeComponentBasic, treeNode
    use, intrinsic :: ISO_C_Binding       , only : c_size_t
    use            :: Numerical_Comparison, only : Values_Agree
    use            :: String_Handling     , only : operator(//)
    implicit none
    type            (enumerationTreeBuildType     )                         :: augmentBuildTreeFromNode
    class           (mergerTreeOperatorAugment    ), intent(inout)          :: self
    type            (treeNode                     ), intent(inout), pointer :: node
    double precision                               , intent(in   )          :: timeEarliestIn               , tolerance
    double precision                               , intent(inout)          :: treeBestWorstFit             , massCutoffScale               , &
         &                                                                     massOvershootScale
    logical                                        , intent(in   )          :: extendingEndNode             , treeBestOverride
    logical                                        , intent(inout)          :: treeNewHasNodeAboveResolution, treeBestHasNodeAboveResolution, &
         &                                                                     newRescale
    type            (mergerTree                   ), intent(inout)          :: treeBest
    type            (treeNode                     ), pointer                :: nodeBase                     , childNode                     , &
         &                                                                     primaryProgenitorNode
    class           (nodeComponentBasic           ), pointer                :: basic                        , basicBase                     , &
         &                                                                     basicChild
    type            (mergerTree                   )                         :: newTree
    double precision                                                        :: timeEarliest                 , timeFirstChild                , &
         &                                                                     massInChildren               , newTreeBaseMass
    integer         (c_size_t                     )                         :: timeIndex
    integer                                                                 :: endNodeCount                 , nodeChildCount                , &
         &                                                                     status
    type            (enumerationTreeBuildType     )                         :: treeAccepted
    type            (mergerTreeOperatorPruneByTime)                         :: pruneByTime
    logical                                                                 :: newTreeBest
    type            (varying_string               )                         :: message
    character       (len=16                       )                         :: label
    logical                                                                 :: primaryProgenitorIsClone

    ! Record the primary progenitor of the node, and determine if it is a clone.
    primaryProgenitorNode    => node%firstChild
    primaryProgenitorIsClone =  .false.
    if (associated(primaryProgenitorNode)) then
       ! Detect if the primary progenitor is a clone.
       basic                    => node                 %basic()
       basicChild               => primaryProgenitorNode%basic()
       primaryProgenitorIsClone =  Values_Agree(basic%time(),basicChild%time(),relTol=1.0d-5)
       ! Remove cloned primary progenitor from the tree.
       if (primaryProgenitorIsClone) node%firstChild => primaryProgenitorNode%sibling
    end if
    ! Find the earliest time to which the tree should be built.
    basic => node%basic()
    if (extendingEndNode) then
       ! An end-node is being extended - we can build the tree all the way back to the earliest time.
       timeEarliest =  timeEarliestIn
    else if (associated(node%firstChild)) then
       ! The node has children - set the limit time to the time at which the children exist.
       basicChild   => node      %firstChild%basic()
       timeEarliest =  basicChild%           time ()
    else if (allocated(self%timeSnapshots)) then
       if (size(self%timeSnapshots) == 1) then
          ! Only one snapshot time is given, use the earliest time.
          timeEarliest =  timeEarliestIn
       else
          timeIndex=searchArrayClosest(self%timeSnapshots,basic%time(),tolerance=1.0d-4,status=status)
          if (status /= errorStatusSuccess) then
             write (label,'(f8.5)') basic%time()
             message='Failed to find matching snapshot time for time '//trim(label)//' Gyr.'//char(10)
             write (label,'(f8.5)') self%timeSnapshots(timeIndex)
             message=message//'Closest time was '//trim(label)//' Gyr'
             call Error_Report(message//{introspection:location})
          end if
          if (timeIndex == 1) then
             timeEarliest=timeEarliestIn
          else
             timeEarliest=self%timeSnapshots(timeIndex-1)
          end if
       end if
    else
       timeEarliest=timeEarliestIn
    end if
    ! Check that all children exists at the same time, and accumulate mass in children.
    massInChildren=0.0d0
    if (associated(node%firstChild)) then
       childNode     => node%firstChild
       basic         => childNode%basic()
       timeFirstChild=basic%time()
       do while (associated(childNode))
          basic     => childNode%basic()
          massInChildren=massInChildren+basic%mass()
          if (basic%time() /= timeFirstChild) then
             message="Child node ["
             write (label,'(f8.5)') basic%time()
             message=message                                              // &
                  &  childNode      %index()                              // &
                  &  "] of parent node ["                                 // &
                  &  node           %index()                              // &
                  &  "] does not exist at same time as first child node ["// &
                  &  node%firstChild%index()                              // &
                  &  "]: "                                                // &
                  &  trim(adjustl(label))
             write (label,'(f8.5)') timeFirstChild
             message=message                                                                                               // &
                  &  " Gyr ≠ "                                                                                             // &
                  &  trim(adjustl(label))                                                                                  // &
                  &  " Gyr."                                                                                     //char(10)// &
                  &  "Consider using the 'regridTimes' operator first to force halos onto a fixed array of times"
             call Error_Report(message//{introspection:location})
          end if
          childNode => childNode%sibling
       end do
    end if
    ! Create a new base node, matched to the current node, build a tree from it, and truncate that tree to the desired earliest time.
    pruneByTime                    =  mergerTreeOperatorPruneByTime(                                &
         &                                                          timeEarliest,                   &
         &                                                               0.0d0  ,                   &
         &                                                          huge(0.0d0) ,                   &
         &                                                          self%cosmologyFunctions_        &
         &                                                         )
    nodeBase                       => treeNode                     (node%index(),newTree          )
    basicBase                      => nodeBase%basic               (             autoCreate=.true.)
    basic                          => node    %basic               (                              )
    newTree%nodeBase               => nodeBase
    newTree%randomNumberGenerator_ => node%hostTree%randomNumberGenerator_
    !$omp atomic
    indexNewTree =indexNewTree+1
    newTree%index=indexNewTree
    ! Determine the mass to use for the base mass of the new tree. This will be the mass of the node in the original tree if this
    ! exceeds the scale mass of the child nodes, otherwise, set to the scaled mass of the child nodes.
    newTreeBaseMass=max(basic%mass(),massInChildren)
    call basicBase  %                    timeSet        (basic          %time())
    call basicBase  %                    massSet        (newTreeBaseMass       )
    call self       %mergerTreeBuilder_ %timeEarliestSet(timeEarliest          )
    call self       %mergerTreeBuilder_ %build          (newTree               )
    call pruneByTime%operatePreEvolution                (newTree               )
    ! Assert that the new tree has some branches.
    if (.not.(extendingEndNode.or.associated(newTree%nodeBase%firstChild))) then
       message="proposed tree has no branches - check cut off mass settings"//char(10)
       write (label,'(e16.8)') newTreeBaseMass
       message=message//"   tree root mass = "//trim(label)//char(10)
       write (label,'(e16.8)') self%massCutOff
       message=message//"     cut off mass = "//trim(label)
       call Error_Report(message//{introspection:location})
    end if
    ! Sort children of our node by mass, and gather statistics on number of children and number of end-nodes in the new tree.
    call self%sortChildren(node)
    nodeChildCount=augmentChildCount    (node                             )
    endNodeCount  =augmentTreeStatistics(newTree,treeStatisticEndNodeCount)
    ! Determine if the newly built tree is an acceptable match.
    treeAccepted  =self%acceptTree(                                &
         &                         node                          , &
         &                         newTree                       , &
         &                         nodeChildCount                , &
         &                         extendingEndNode              , &
         &                         tolerance                     , &
         &                         treeBest                      , &
         &                         treeBestWorstFit              , &
         &                         treeBestOverride              , &
         &                         massCutoffScale               , &
         &                         massOvershootScale            , &
         &                         treeNewHasNodeAboveResolution , &
         &                         treeBestHasNodeAboveResolution, &
         &                         newTreeBest                   , &
         &                         primaryProgenitorNode         , &
         &                         primaryProgenitorIsClone        &
         &                         )
    ! Determine whether to use stored best tree, newly created tree, or reject the tree.
    if     (                                          &
         &   associated(treeBest%nodeBase)            &
         &  .and.                                     &
         &   treeAccepted /= treeBuildSuccess         &
         &  .and.                                     &
         &   (                                        &
         &     (                                      &
         &       treeBestWorstFit <= tolerance        &
         &      .and.                                 &
         &       .not.treeBestHasNodeAboveResolution  &
         &     )                                      &
         &     .and.                                  &
         &      .not.newTreeBest                      &
         &     .and.                                  &
         &      newRescale                            &
         &    .or.                                    &
         &     treeBestOverride                       &
         &   )                                        &
         & ) then
       ! The previously stored best matching tree can be used iff:
       !  1) A previously stored best matching tree exists;
       !  2) The newly built tree was not accepted;
       !  3) The previously stored best tree is now within tolerance and has no nodes above the cut-off, or its use is being
       !     forced.
       !
       ! Clean up the newly created tree, and replace it with the best tree.
       if (associated(newTree%nodeBase)) then
          call newTree%nodeBase%destroyBranch()
          deallocate(newTree%nodeBase)
       end if
       newTree%nodeBase => treeBest%nodeBase
       ! Reset the best tree.
       treeBestWorstFit  =  3.0d0
       treeBest%nodeBase => null()
       ! Test for acceptance of the best tree.
       if     (                                                 &
            &   self%acceptTree(                                &
            &                   node                          , &
            &                   newTree                       , &
            &                   nodeChildCount                , &
            &                   extendingEndnode              , &
            &                   tolerance                     , &
            &                   treeBest                      , &
            &                   treeBestWorstFit              , &
            &                   treeBestOverride              , &
            &                   massCutoffScale               , &
            &                   massOvershootScale            , &
            &                   treeNewHasNodeAboveResolution , &
            &                   treeBestHasNodeAboveResolution, &
            &                   newTreeBest                   , &
            &                   primaryProgenitorNode         , &
            &                   primaryProgenitorIsClone        &
            &                  )                                &
            &  ==                                               &
            &   treeBuildSuccess                                &
            & ) then
          augmentBuildTreeFromNode=treeBuildSuccess
       else
          augmentBuildTreeFromNode=treeBuildFailureStructure
       end if
    else if (treeAccepted == treeBuildSuccess) then
       ! If the newly created tree was acceptable, use it.
       augmentBuildTreeFromNode=treeBuildSuccess
       ! Clean up any previously stored best tree.
       if (associated(treeBest%nodeBase)) then
          call treeBest%nodeBase%destroyBranch()
          deallocate(treeBest%nodeBase)
       end if
       treeBestWorstFit  =  3.0d0
       treeBest%nodeBase => null()
    else
       ! The newly created tree was unacceptable, clean it up and return the failure code.
       if (associated(newTree%nodeBase)) then
          call newTree%nodeBase%destroyBranch()
          deallocate(newTree%nodeBase)
       end if
       augmentBuildTreeFromNode=treeAccepted
    end if
    ! Put any cloned progenitor back in.
    if (augmentBuildTreeFromNode /= treeBuildSuccess .and. associated(primaryProgenitorNode)) then
       if (primaryProgenitorIsClone) then
          ! Simply reinsert the clone, but ensure we update its sibling pointer in case the
          ! order of children was altered due to them being sorted by mass.
          primaryProgenitorNode%sibling    => node                 %firstChild
          node                 %firstChild => primaryProgenitorNode
       else
          ! Move the original primary progenitor back into the primary progenitor slot.
          if (.not.associated(node%firstChild,primaryProgenitorNode)) then
             childNode => node%firstChild
             do while (.not.associated(childNode%sibling,primaryProgenitorNode))
                childNode => childNode%sibling
             end do
             ! Modify links to make our non-overlap node the primary progenitor.
             childNode            %sibling    => primaryProgenitorNode%sibling
             primaryProgenitorNode%sibling    => node                 %firstChild
             node                 %firstChild => primaryProgenitorNode
          end if
       end if
    end if
    return
  end function augmentBuildTreeFromNode

  recursive function augmentAcceptTree(self,node,tree,nodeChildCount,extendingEndNode,tolerance,treeBest,treeBestWorstFit,treeBestOverride,massCutoffScale,massOvershootScale,treeNewHasNodeAboveResolution,treeBestHasNodeAboveResolution,newTreeBest,primaryProgenitorNode,primaryProgenitorIsClone)
    !!{
    Determine whether a trial tree is an acceptable match to the original tree structure.
    !!}
    use :: Display            , only : displayMessage               , displayVerbosity  , verbosityLevelWorking
    use :: Galacticus_Nodes   , only : mergerTree                   , nodeComponentBasic, treeNode             , treeNodeList
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    implicit none
    type            (enumerationTreeBuildType     )                                     :: augmentAcceptTree
    class           (mergerTreeOperatorAugment    ), intent(inout)                      :: self
    type            (treeNode                     ), intent(inout)            , pointer :: node                         , primaryProgenitorNode
    type            (treeNode                     )                           , pointer :: nodeCurrent                  , nodePrevious                  , &
         &                                                                                 nodeNonOverlap               , nodeNonOverlapFirst           , &
         &                                                                                 nodeSatellite                , nodeNext
    class           (nodeComponentBasic           )                           , pointer :: basicCurrent                 , basicSort
    type            (mergerTree                   ), intent(inout)            , target  :: tree                         , treeBest
    logical                                        , intent(inout)                      :: treeNewHasNodeAboveResolution, treeBestHasNodeAboveResolution, &
         &                                                                                 newTreeBest
    logical                                        , intent(in   )                      :: treeBestOverride             , extendingEndNode              , &
         &                                                                                 primaryProgenitorIsClone
    double precision                               , intent(inout)                      :: treeBestWorstFit             , massCutoffScale               , &
         &                                                                                 massOvershootScale
    double precision                               , intent(in   )                      :: tolerance
    integer                                        , intent(inout)                      :: nodeChildCount
    type            (treeNodeList                 ), dimension(nodeChildCount)          :: endNodes
    type            (mergerTreeWalkerIsolatedNodes)                                     :: treeWalker
    integer                                                                             :: i                            , j                             , &
         &                                                                                 endNodesSorted
    logical                                                                             :: treeAccepted                 , nodeMassesAgree               , &
         &                                                                                 nodeCurrentBelowAll          , nodesRemain
    double precision                                                                    :: treeCurrentWorstFit
    type            (varying_string               )                                     :: message
    character       (len=12                       )                                     :: label

    ! Initialize.
    treeNewHasNodeAboveResolution =  .false.
    nodeNonOverlapFirst           => null()
    i                             =  1
    endNodesSorted                =  1
    ! Walk through the tree identifying end-nodes.
    treeWalker=mergerTreeWalkerIsolatedNodes(tree)
    nodesRemain=treeWalker%next(nodeNext)
    do while (nodesRemain)
       nodeCurrent =>                 nodeNext
       nodesRemain =  treeWalker%next(nodeNext)
       ! Initialize the current node.
       basicCurrent         => nodeCurrent%basic   ()
       nodeCurrent%hostTree => node       %hostTree
       nodeCurrent%event    => null()
       nodeNonOverlap       => null()
       ! Test for children.
       if (.not.associated(nodeCurrent%firstChild)) then
          ! Store this node in the mass-ordered array of end-nodes such that we build up a list of the most massive overlapping
          ! nodes.
          if (nodeChildCount > 0) then
             ! If insufficient potential overlap nodes have yet to be found, add this node to the list of potential overlap nodes.
             if (endNodesSorted == 1) then
                ! The first node is simply placed into the first element of our end-node array.
                nodeCurrentBelowAll           =  .false.
                endNodes(endNodesSorted)%node => nodeCurrent
             else
                ! Iterate over current list of end-nodes to find insertion point.
                i                  =1
                nodeCurrentBelowAll=.true.
                do while (i < min(endNodesSorted,nodeChildCount))
                   ! Test mass of current node in sorted list.
                   basicSort => endNodes(i)%node%basic()
                   if (basicSort%mass() < basicCurrent%mass()) then
                      ! Current node exceeds the mass of one or more nodes in the current list.
                      nodeCurrentBelowAll=.false.
                      ! Check for a node being pushed off the end of the list.
                      if (endNodesSorted > nodeChildCount) then
                         nodeNonOverlap => endNodes(nodeChildCount)%node
                         ! Test if the node being lost from the end of the list is above the mass resolution.
                         basicSort      => endNodes(nodeChildCount)%node%basic()
                         if (basicSort%mass() > self%massCutOff*massCutoffScale) then
                            if (displayVerbosity() >= verbosityLevelWorking) then
                               write (label,'(e12.6)') basicSort%mass()
                               message="Nonoverlap failure at mass: "//trim(label)
                               call displayMessage(message)
                            end if
                            ! Node being lost from end of list is above the mass resolution - record that we therefore have a new,
                            ! non-overlapping node which is above the mass resolution.
                            treeNewHasNodeAboveResolution=.true.
                         end if
                      end if
                      ! Shift nodes down in the sorted array to make space for the insertion.
                      j=min(endNodesSorted,nodeChildCount)
                      do while (j > i)
                         endNodes(j)%node => endNodes(j-1)%node
                         j                =           j-1
                      end do
                      ! Insert the new node, and trigger exit of the loop.
                      endNodes(i)%node => nodeCurrent
                      exit
                   end if
                   i=i+1
                end do
                ! Handle cases where the current node has lower mass than all other end-nodes currently in the list.
                if (nodeCurrentBelowAll) then
                   if (endNodesSorted <= nodeChildCount) then
                      ! The list of end-nodes is not yet full - insert the node into the next available space.
                      endNodes(endNodesSorted)%node => nodeCurrent
                   else
                      ! The list of end nodes is full - this node can not be inserted into the list.
                      nodeNonOverlap => nodeCurrent
                      if (basicCurrent%mass() > self%massCutOff*massCutoffScale) then
                         if (displayVerbosity() >= verbosityLevelWorking) then
                            write (label,'(e12.6)') basicCurrent%mass()
                            message="Nonoverlap failure at mass: "//trim(label)
                            call displayMessage(message)
                         end if
                         treeNewHasNodeAboveResolution=.true.
                      end if
                   end if
                end if
             end if
          else
             ! No overlap nodes are being sought - so this node is automatically a non-overlap node.
             nodeNonOverlap => nodeCurrent
             if (basicCurrent%mass() > self%massCutOff*massCutoffScale) then
                if (displayVerbosity() >= verbosityLevelWorking) then
                   write (label,'(e12.6)') basicCurrent%mass()
                   message="Nonoverlap failure at mass: "//trim(label)
                   call displayMessage(message)
                end if
                treeNewHasNodeAboveResolution=.true.
             end if
          end if
          endNodesSorted=endNodesSorted+1
       end if
       ! Add the non-overlap node (which is either the current node, or the one it pushed off the end of the most-massive nodes
       ! list) to the list of non-overlap nodes.
       if (.not.extendingEndNode) call augmentNonOverlapListAdd(nodeNonOverlap,nodeNonOverlapFirst)
    end do
    endNodesSorted=endNodesSorted-1
    ! Test for mass-matches between overlap nodes.
    nodeMassesAgree    =.true.
    treeCurrentWorstFit=0.0d0
    if (nodeChildCount > 0 .and. nodeChildCount <= endNodesSorted) then
       i           =  1
       nodeCurrent => node%firstChild
       do while (associated(nodeCurrent).and.nodeMassesAgree)
          if (.not.augmentNodeComparison(nodeCurrent,endNodes(i)%node,tolerance,treeCurrentWorstFit)) nodeMassesAgree=.false.
          i           =  i+1
          nodeCurrent => nodeCurrent%sibling
       end do
    end if
    ! If we have a matching tree, extend the non-overlapping nodes.
    if     (                                      &
         &   nodeChildCount <= endNodesSorted     &
         &  .and.                                 &
         &   (                                    &
         &     .not.treeNewHasNodeAboveResolution &
         &    .or.                                &
         &          treeBestOverride              &
         &   )                                    &
         &  .and.                                 &
         &   nodeMassesAgree                      &
         & ) then
       call self%extendNonOverlapNodes(                     &
            &                          nodeNonOverlapFirst, &
            &                          tolerance          , &
            &                          treeBest           , &
            &                          massCutoffScale    , &
            &                          massOvershootScale   &
            &                         )
    end if
    ! Determine if the tree is accepted.
    if (treeBestOverride .and. associated(tree%nodeBase)) then
       treeAccepted=      nodeMassesAgree                                 &
            &       .and.                                                 &
            &             nodeChildCount                <= endNodesSorted
    else
       treeAccepted=      nodeMassesAgree                                 &
            &       .and.                                                 &
            &             nodeChildCount                <= endNodesSorted &
            &       .and.                                                 &
            &        .not.treeNewHasNodeAboveResolution
    end if
    ! Process the tree depending on acceptance state.
    newTreeBest=.false.
    if (treeAccepted) then
       ! Tree was accepted, insert it into the original tree.
       call augmentSimpleInsert(self,node,tree,endNodes,nodeChildCount,nodeNonOverlapFirst)
       ! Fix up primary progenitor status here.
       if (associated(primaryProgenitorNode)) then
          ! Our node originally had a primary progenitor - so we must reinstate it.
          if (primaryProgenitorIsClone) then
             ! The primary progenitor of the original node was a clone. If we have an available non-overlap node, use that as the
             ! new primary progenitor, otherwise reinstate the clone.
             if (associated(nodeNonOverlapFirst)) then
                ! A non-overlap node is available to use as the primary progenitor. Step up the tree from it forcing it to be the
                ! primary progenitor at each step.
                nodeCurrent => nodeNonOverlapFirst
                do while (.not.associated(nodeCurrent,node))
                   ! If our non-overlap node branch is not the primary progenitor at this point, switch the ordering of
                   ! progenitors to make it so.
                   if (.not.nodeCurrent%isPrimaryProgenitor()) then
                      ! Find the sibling immediately prior to our non-overlap node.
                      nodePrevious => nodeCurrent%parent%firstChild
                      do while (.not.associated(nodePrevious%sibling,nodeCurrent))
                         nodePrevious => nodePrevious%sibling
                      end do
                      ! Modify links to make our non-overlap node the primary progenitor.
                      nodePrevious%sibling           => nodeCurrent%sibling
                      nodeCurrent %sibling           => nodeCurrent%parent %firstChild
                      nodeCurrent %parent%firstChild => nodeCurrent
                   end if
                   nodeCurrent => nodeCurrent%parent
                end do
                ! If the clone primary progenitor contained any satellites, shift them to the parent node.
                if (associated(node%firstSatellite)) then
                   nodeSatellite => node%firstSatellite
                   do while (associated(nodeSatellite%sibling))
                      nodeSatellite => nodeSatellite%sibling
                   end do
                   nodeSatellite%sibling        => primaryProgenitorNode%firstSatellite
                else
                   node         %firstSatellite => primaryProgenitorNode%firstSatellite
                end if
                nodeSatellite => primaryProgenitorNode%firstSatellite
                do while (associated(nodeSatellite))
                   nodeSatellite%parent => node
                   nodeSatellite        => nodeSatellite%sibling
                end do
             else
                ! Reinstate the cloned progenitor.
                primaryProgenitorNode%sibling    => node                 %firstChild
                node                 %firstChild => primaryProgenitorNode
             end if
          else
             ! Ensure that the original primary progenitor remains the primary progenitor.
             nodeCurrent => primaryProgenitorNode
             do while (.not.associated(nodeCurrent,node))
                ! If our node branch is not the primary progenitor at this point, switch the ordering of progenitors to make it
                ! so.
                if (.not.nodeCurrent%isPrimaryProgenitor()) then
                   ! Find the sibling immediately prior to our non-overlap node.
                   nodePrevious => nodeCurrent%parent%firstChild
                   do while (.not.associated(nodePrevious%sibling,nodeCurrent))
                      nodePrevious => nodePrevious%sibling
                   end do
                   ! Modify links to make our non-overlap node the primary progenitor.
                   nodePrevious%sibling           => nodeCurrent%sibling
                   nodeCurrent %sibling           => nodeCurrent%parent %firstChild
                   nodeCurrent %parent%firstChild => nodeCurrent
                end if
                nodeCurrent => nodeCurrent%parent
             end do
          end if
       end if
    else if (                                                      &
         &         nodeChildCount                <= endNodesSorted &
         &   .and.                                                 &
         &    .not.treeNewHasNodeAboveResolution                   &
         &   .and.                                                 &
         &    .not.treeBestOverride                                &
         &  ) then
       ! Tree is structurally acceptable - decide if we want to keep it as the current best tree.
       call self%nonOverlapReinsert(nodeNonOverlapFirst)
       if (treeCurrentWorstFit < treeBestWorstFit) then
          newTreeBest = .true.
          ! Current tree is better than the current best tree. Replace the best tree with the current tree.
          if (associated(treeBest%nodeBase)) then
             call treeBest%nodeBase%destroyBranch()
             deallocate(treeBest%nodeBase)
          end if
          treeBest                      %nodeBase          => tree                         %nodeBase
          tree                          %nodeBase          => null()
          treeBest                      %nodeBase%hostTree => treeBest
          treeBestWorstFit                                 =  treeCurrentWorstFit
          treeBestHasNodeAboveResolution                   =  treeNewHasNodeAboveResolution
          treeWalker                                       =  mergerTreeWalkerIsolatedNodes         (treeBest)
          do while (treeWalker%next(nodeCurrent))
             nodeCurrent%event    => null()
             nodeCurrent%hostTree => treeBest%nodeBase%hostTree
          end do
       end if
    else
       ! Tree is not acceptable or better than the current best tree - destroy it.
       call self%nonOverlapReinsert(nodeNonOverlapFirst)
       if (associated(tree%nodeBase)) then
          call tree%nodeBase%destroyBranch()
          deallocate(tree%nodeBase)
       end if
    end if
    ! Return a suitable status code based on tree acceptance criteria.
    if (treeAccepted) then
      augmentAcceptTree=treeBuildSuccess
    else if (.not.nodeMassesAgree) then
      augmentAcceptTree=treeBuildFailureTolerance
    else
      augmentAcceptTree=treeBuildFailureStructure
    end if
    return
  end function augmentAcceptTree

  subroutine augmentSimpleInsert(self,node,tree,endNodes,nodeChildCount,nodeNonOverlapFirst)
    !!{
    Insert a newly constructed tree into the original tree.
    !!}
    use :: Galacticus_Nodes, only : mergerTree, treeNode, treeNodeList
    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)                            :: self
    type   (treeNode                 ), intent(inout), pointer                   :: node
    type   (mergerTree               ), intent(inout), target                    :: tree
    integer                           , intent(in   )                            :: nodeChildCount
    type   (treeNode                 ), intent(inout), pointer                   :: nodeNonOverlapFirst
    type   (treeNodeList             ), intent(inout), dimension(nodeChildCount) :: endNodes
    type   (treeNode                 ), pointer                                  :: nodeCurrent        , nodeCurrentSibling
    integer                                                                      :: i

    ! Reinsert the non-overlap nodes which have been removed from the tree during acceptance testing.
    call self%nonOverlapReinsert(nodeNonOverlapFirst)
    ! Test if we have overlap nodes.
    if (nodeChildCount > 0) then
       ! Overlap nodes exist - connect them in to the original tree.
       i           =  1
       nodeCurrent => node%firstChild
       do while (associated(nodeCurrent))
          ! For each child of the original node, we first shift the firstChild pointer of the parent to the sibling, then overlap
          ! the corresponding node from the new tree with the child. After the last such overlap, the firstChild pointer of the
          ! parent is null - ready for overlap with the base node of the new tree (see below).
          nodeCurrentSibling            => nodeCurrent%sibling
          node              %firstChild => nodeCurrent%sibling
          call augmentExtendByOverLap(                            &
               &                      endNodes(i)%node          , &
               &                      nodeCurrent               , &
               &                      keepTop           =.true. , &
               &                      exchangeProperties=.false.  &
               &                     )
          nodeCurrent => nodeCurrentSibling
          i           =  i                 +1
       end do
    end if
    ! Overlap the base node of the tree with the node in the original tree.
    call augmentExtendByOverLap(                            &
         &                      node                      , &
         &                      tree%nodeBase             , &
         &                      keepTop           =.false., &
         &                      exchangeProperties=.false.  &
         &                     )
    return
  end subroutine augmentSimpleInsert

  subroutine augmentExtendByOverlap(nodeBottom,nodeTop,keepTop,exchangeProperties)
    !!{
    Conjoin two trees by overlapping the {\normalfont \ttfamily nodeTop} of one tree with the chosen {\normalfont \ttfamily
    nodeBottom} of the other. If {\normalfont \ttfamily keepTop} is {\normalfont \ttfamily true}, {\normalfont \ttfamily
    nodeTop} replaces {\normalfont \ttfamily nodeBottom}, otherwise, {\normalfont \ttfamily nodeBottom} replaces {\normalfont
    \ttfamily nodeTop}. If {\normalfont \ttfamily exchangeProperties} is {\normalfont \ttfamily true}, the mass and time
    information of the deleted node overwrites the mass and time of the retained node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    type   (treeNode          ), intent(inout), pointer :: nodeBottom  , nodeTop
    logical                    , intent(in   )          :: keepTop     , exchangeProperties
    type   (treeNode          )               , pointer :: currentChild, currentSibling
    class  (nodeComponentBasic)               , pointer :: basicBottom , basicTop

    basicBottom => nodeBottom%basic()
    basicTop    => nodeTop   %basic()
    if (keepTop) then
       nodeTop%parent  => nodeBottom%parent
       nodeTop%sibling => nodeBottom%sibling
       if (associated(nodeBottom%parent)) then
          currentSibling => nodeBottom%parent%firstChild
       else
          currentSibling => null()
       end if
       do while (associated(currentSibling))
          if (associated(currentSibling%sibling,nodeBottom)) currentSibling%sibling => nodeTop
          currentSibling => currentSibling%sibling
       end do
       if (associated(nodeBottom%parent)) then
          if (associated(nodeBottom%parent%firstChild, nodeBottom)) nodeBottom%parent%firstChild => nodeTop
       end if
       if (exchangeProperties) then
          call basicTop%timeSet(basicBottom%time())
          call basicTop%massSet(basicBottom%mass())
       end if
       call nodeBottom%destroy()
       deallocate(nodeBottom)
    else
       currentChild          => nodeTop%firstChild
       nodeBottom%firstChild => nodeTop%firstChild
       do while (associated(currentChild))
          currentChild%parent => nodeBottom
          currentChild        => currentChild%sibling
       end do
       if (exchangeProperties) then
          call basicBottom%timeSet(basicTop%time())
          call basicBottom%massSet(basicTop%mass())
       end if
       call nodeTop%destroy()
       deallocate(nodeTop)
    end if
    return
  end subroutine augmentExtendByOverlap

  integer function augmentChildCount(node)
    !!{
    Return a count of the number of child nodes.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    type(treeNode), intent(in   ) :: node
    type(treeNode), pointer       :: childNode

    childNode         => node%firstChild
    augmentChildCount =  0
    do while (associated(childNode))
       if (associated(childNode)) augmentChildCount=augmentChildCount+1
       childNode => childNode%sibling
    end do
    return
  end function augmentChildCount

  subroutine augmentSortChildren(self,node)
    !!{
    Sort the children of the given {\normalfont \ttfamily node} such that they are in descending mass order.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (mergerTreeOperatorAugment), intent(in   )          :: self
    type            (treeNode                 ), intent(inout), target  :: node
    type            (treeNode                 )               , pointer :: nodeCurrent , nodeSort , &
         &                                                                 nodeNext
    class           (nodeComponentBasic       )               , pointer :: basicCurrent, basicSort
    double precision                                                    :: massLargest

    if (.not.associated(node%firstChild)) return
    nodeCurrent          => node%firstChild
    basicCurrent         => nodeCurrent%basic  ()
    nodeNext             => nodeCurrent%sibling
    nodeCurrent %sibling => null()
    massLargest          =  basicCurrent%mass  ()
    do while (associated(nodeNext))
       nodeSort     => node       %firstChild
       nodeCurrent  => nodeNext
       nodeNext     => nodeNext   %sibling
       basicCurrent => nodeCurrent%basic     ()
       if (basicCurrent%mass() > massLargest) then
          node       %firstChild => nodeCurrent
          nodeCurrent%sibling    => nodeSort
          massLargest            =  basicCurrent%mass()
       else
          do while (associated(nodeSort%sibling))
             basicSort => nodeSort%sibling%basic()
             if (basicCurrent%mass() > basicSort%mass()) exit
             nodeSort  => nodeSort%sibling
          end do
          nodeCurrent%sibling => nodeSort   %sibling
          nodeSort   %sibling => nodeCurrent
       end if
    end do
    ! Check results.
    if (self%performChecks) then
       nodeCurrent  => node        %firstChild
       basicCurrent => nodeCurrent %basic     ()
       massLargest  =  basicCurrent%mass      ()
       nodeCurrent  => nodeCurrent %sibling
       do while (associated(nodeCurrent))
          basicCurrent => nodeCurrent%basic()
          if (basicCurrent%mass() > massLargest) call Error_Report('failed to sort child nodes'//{introspection:location})
          massLargest =  basicCurrent%mass   ()
          nodeCurrent => nodeCurrent %sibling
       end do
    end if
    return
  end subroutine augmentSortChildren

  subroutine augmentNonOverlapListAdd(node,listFirstElement)
    !!{
    Add the given node to a linked list of non-overlap nodes in the current trial tree.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type(treeNode), intent(in   ), pointer :: node
    type(treeNode), intent(inout), pointer :: listFirstElement
    type(treeNode)               , pointer :: currentSibling

    ! Return immediately if no node is given.
    if (.not.associated(node)) return
    ! Check if node has a parent node.
    if (associated(node%parent)) then
       ! Unlink the node from its parent.
       currentSibling => node%parent%firstChild
       if (associated(currentSibling, node).and.(associated(node%sibling))) then
          node%parent%firstChild => node%sibling
       else
          do while (associated(currentSibling))
             if (associated(currentSibling%sibling,node)) currentSibling%sibling => node%sibling
             currentSibling => currentSibling%sibling
          end do
       end if
    end if
    ! Insert the node into the list.
    if (.not.associated(listFirstElement)) then
       listFirstElement         => node
       listFirstElement%sibling => null()
    else
       node%sibling     => listFirstElement
       listFirstElement => node
    end if
    return
  end subroutine augmentNonOverlapListAdd

  subroutine augmentNonOverlapReinsert(self,listFirstElement)
    !!{
    Reinsert non-overlap nodes into their tree.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class(mergerTreeOperatorAugment), intent(in   )          :: self
    type (treeNode                 ), intent(in   ), pointer :: listFirstElement
    type (treeNode                 )               , pointer :: nodeCurrent     , listHead

    listHead => listFirstElement
    do while (associated(listHead))
       nodeCurrent => listHead
       listHead    => listHead%sibling
       if (associated(nodeCurrent%parent)) then
          if (associated(nodeCurrent,nodeCurrent%parent%firstChild)) then
             nodeCurrent%sibling => null()
          else
             nodeCurrent%sibling => nodeCurrent%parent%firstChild
          end if
          nodeCurrent%parent%firstChild => nodeCurrent
          call self%sortChildren(nodeCurrent%parent)
       end if
    end do
    return
  end subroutine augmentNonOverlapReinsert

  subroutine augmentExtendNonOverlapNodes(self,nodeNonOverlapFirst,tolerance,treeBest,massCutoffScale,massOvershootScale)
    !!{
    Extend any non-overlap nodes in an accepted tree by growing a new tree from each such node.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : mergerTree  , treeNode
    implicit none
    class           (mergerTreeOperatorAugment), intent(inout)          :: self
    type            (mergerTree               ), intent(inout), target  :: treeBest
    type            (treeNode                 ), intent(in   ), pointer :: nodeNonOverlapFirst
    double precision                           , intent(in   )          :: tolerance
    double precision                           , intent(inout)          :: massCutoffScale        , massOvershootScale
    type            (treeNode                 )               , pointer :: nodeCurrent
    double precision                                                    :: falseWorstFit
    logical                                                             :: falseNewNodeAboveCutoff, falseBestTreeNodeAboveCutoff, &
         &                                                                 falseNewRescale
    type            (enumerationTreeBuildType )                         :: treeStatus

    falseBestTreeNodeAboveCutoff =  .false.
    falseNewNodeAboveCutoff      =  .false.
    falseNewRescale              =  .false.
    falseWorstFit                =  3.0d0
    nodeCurrent                  => nodeNonOverlapFirst
    do while (associated(nodeCurrent))
       treeStatus=self%buildTreeFromNode(                              &
            &                            nodeCurrent                 , &
            &                            .true.                      , &
            &                            tolerance                   , &
            &                            self%timeEarliest           , &
            &                            treeBest                    , &
            &                            falseWorstFit               , &
            &                            .false.                     , &
            &                            massCutoffScale             , &
            &                            massOvershootScale          , &
            &                            falseNewNodeAboveCutoff     , &
            &                            falseBestTreeNodeAboveCutoff, &
            &                            falseNewRescale               &
            &                           )
       if (treeStatus /= treeBuildSuccess) call Error_Report('extension of non-overlap node failed'//{introspection:location})
       nodeCurrent => nodeCurrent%sibling
    end do
    return
  end subroutine augmentExtendNonOverlapNodes

  logical function augmentNodeComparison(nodeNew,nodeOriginal,tolerance,treeCurrentWorstFit)
    !!{
    Compare the masses of an overlap node and the corresponding node in the original tree, testing if they agree to withing
    tolerance.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    type            (treeNode          ), intent(inout), pointer :: nodeNew            , nodeOriginal
    double precision                    , intent(in   )          :: tolerance
    double precision                    , intent(inout)          :: treeCurrentWorstFit
    class           (nodeComponentBasic)               , pointer :: basicNew           , basicOriginal
    double precision                                             :: fitMeasure

    basicNew              => nodeNew     %basic()
    basicOriginal         => nodeOriginal%basic()
    fitMeasure            =  +2.0d0                                     &
         &                   *abs(basicNew%mass()-basicOriginal%mass()) &
         &                   /abs(basicNew%mass()+basicOriginal%mass())
    augmentNodeComparison =  fitMeasure < tolerance
    treeCurrentWorstFit   =  max(treeCurrentWorstFit,fitMeasure)
    return
  end function augmentNodeComparison

  integer function augmentTreeStatistics(tree,desiredOutput)
    !!{
    Walks through tree and quietly collects information specified by {\normalfont \ttfamily desiredOutput} input enumeration and
    returns that information.
    !!}
    use :: Error              , only : Error_Report
    use :: Galacticus_Nodes   , only : mergerTree                   , treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    implicit none
    type   (mergerTree                   ), intent(in   ), target  :: tree
    type   (enumerationTreeStatisticType ), intent(in   )          :: desiredOutput
    type   (treeNode                     )               , pointer :: node
    type   (mergerTreeWalkerIsolatedNodes)                         :: treeWalker
    integer                                                        :: nodeCount    , endNodeCount

    nodeCount   =0
    endNodeCount=0
    treeWalker  =mergerTreeWalkerIsolatedNodes(tree)
    do while (treeWalker%next(node))
       nodeCount                                         =   nodeCount+1
       if (.not.associated(node%firstChild)) endNodeCount=endNodeCount+1
    end do
    ! Return the requested quantity.
    select case (desiredOutput%ID)
    case (treeStatisticNodeCount   %ID)
       augmentTreeStatistics=nodeCount
    case (treeStatisticEndNodeCount%ID)
       augmentTreeStatistics=endNodeCount
    case default
       augmentTreeStatistics=0
       call Error_Report('unknown task requested'//{introspection:location})
    end select
    return
  end function augmentTreeStatistics

  subroutine augmentFinalize(self)
    !!{
    Output augmentation histogram.
    !!}
    use :: Output_HDF5      , only : outputFile
    use :: HDF5_Access      , only : hdf5Access
    use :: IO_HDF5          , only : hdf5Object
    implicit none
    class           (mergerTreeOperatorAugment), intent(inout)               :: self
    integer         (c_size_t                 ), allocatable  , dimension(:) :: retryHistogram        , trialCount
    type            (hdf5Object               )                              :: augmentStatisticsGroup

    ! Output the data.
    !$ call hdf5Access%set()
    ! Check if our output group already exists.
    if (outputFile%hasGroup('augmentStatistics')) then
       ! Our group does exist. Read existing histogram, add them to our own, then write back to file.
       augmentStatisticsGroup=outputFile%openGroup('augmentStatistics','Statistics of merger tree augmentation.',objectsOverwritable=.true.,overwriteOverride=.true.)
       allocate(retryHistogram(size(self%retryHistogram)))
       allocate(trialCount    (size(self%trialCount    )))
       call augmentStatisticsGroup%readDataset('retryHistogram',retryHistogram)
       call augmentStatisticsGroup%readDataset('trialCount'    ,trialCount    )
       self%retryHistogram=self%retryHistogram+retryHistogram
       self%trialCount    =self%trialCount    +trialCount
       deallocate(retryHistogram)
       deallocate(trialCount    )
    else
       ! Our group does not already exist. Simply write the data.
       augmentStatisticsGroup=outputFile%openGroup('augmentStatistics','Statistics of merger tree augmentation.',objectsOverwritable=.true.,overwriteOverride=.true.)
    end if
    call augmentStatisticsGroup%writeDataset(self%retryHistogram,"retryHistogram","Retry histogram []")
    call augmentStatisticsGroup%writeDataset(self%trialCount    ,"trialCount"    ,"Trial counts []"   )
    !$ call hdf5Access%unset()
    return
  end subroutine augmentFinalize
