!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Contains a module which implements an augmenting operator on merger trees.
  use Merger_Trees_Builders

  !# <mergerTreeOperator name="mergerTreeOperatorAugment">
  !#  <description>Provides a merger tree operator which augments tree resolution by inserting high-resolution branches.</description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorAugment
     !% An augmenting merger tree operator class.
     private
     double precision                        , allocatable, dimension (:) :: timeSnapshots
     double precision                                                     :: massResolution    , timeEarliest             , &
          &                                                                  toleranceScale    , massCutOffScaleFactor
     integer                                                              :: retryMaximum      , rescaleMaximum           , &
          &                                                                  attemptsMaximum   , massCutOffAttemptsMaximum
     logical                                                              :: performChecks     , useOneNodeTrees
     class           (mergerTreeBuilderClass), pointer                    :: mergerTreeBuilder_
   contains
     ! AJB : need argument data adding here
     !@ <objectMethods>
     !@   <object>mergerTreeOperatorAugment</object>
     !@   <objectMethod>
     !@     <method>buildTreeFromNode</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Build a merger tree starting from the given node.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>acceptTree</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Determine if a newly built tree is an acceptable match.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>extendNonOverlapNodes</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Graft new branches onto all end-nodes of a newly built tree.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>multiScale</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Do some scaling.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>sortChildren</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Sort child nodes into descending mass order.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>nonOverlapReinsert</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Reinsert a linked list of non-overlap nodes into their parent tree.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                          augmentDestructor
     procedure :: operate               => augmentOperate
     procedure :: buildTreeFromNode     => augmentBuildTreeFromNode
     procedure :: acceptTree            => augmentAcceptTree
     procedure :: extendNonOverlapNodes => augmentExtendNonOverlapNodes
     procedure :: multiScale            => augmentMultiScale
     procedure :: sortChildren          => augmentSortChildren
     procedure :: nonOverlapReinsert    => augmentNonOverlapReinsert
  end type mergerTreeOperatorAugment

  interface mergerTreeOperatorAugment
     !% Constructors for the {\normalfont \ttfamily augment} merger tree operator class.
     module procedure augmentConstructorParameters
     module procedure augmentConstructorInternal
  end interface mergerTreeOperatorAugment

  !# <enumeration>
  !#  <name>treeStatistic</name>
  !#  <description>Enumeration of tasks to be performed during a tree walk.</description>
  !#  <entry label="nodeCount"   />
  !#  <entry label="endNodeCount"/>
  !# </enumeration>

  !# <enumeration>
  !#  <name>treeBuild</name>
  !#  <description>Enumeration of tree building status.</description>
  !#  <entry label="success"         />
  !#  <entry label="failureTolerance"/>
  !#  <entry label="failureStructure"/>
  !#  <entry label="failureGeneric"  />
  !# </enumeration>

contains

  function augmentConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily augment} merger tree operator class which takes a parameter set as input.
    use Input_Parameters2
    use Cosmology_Functions
    use Memory_Management
    use Galacticus_Error
    implicit none
    type            (mergerTreeOperatorAugment)                              :: augmentConstructorParameters
    type            (inputParameters          ), intent(in   )               :: parameters
    double precision                           , allocatable  , dimension(:) :: timeSnapshots
    class           (cosmologyFunctionsClass  ), pointer                     :: cosmologyFunctions_
    class           (mergerTreeBuilderClass   ), pointer                     :: mergerTreeBuilder_ 
    integer                                                                  :: i
    double precision                                                         :: massResolution              , toleranceScale           , &
         &                                                                      massCutOffScaleFactor
    integer                                                                  :: retryMaximum                , rescaleMaximum           , &
         &                                                                      attemptsMaximum             , massCutOffAttemptsMaximum
    logical                                                                  :: performChecks               , useOneNodeTrees
    !# <inputParameterList label="allowedParameterNames" />

    !# <objectBuilder class="mergerTreeBuilder" name="mergerTreeBuilder_" source="parameters"/>
    !# <inputParameter>
    !#   <name>massResolution</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d10</defaultValue>
    !#   <description>For the {\normalfont \ttfamily augment} operator a description of resolution limit for new trees.</description>
    !#   <type>double precision</type>
    !#   <cardinality>0..</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>performChecks</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true, perform checks of the augmentation process.</description>
    !#   <type>logical</type>
    !#   <cardinality>0..</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>toleranceScale</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.925d0</defaultValue>
    !#   <description>The factor by which to scale the tolerance parameter on each rescaling.</description>
    !#   <type>real</type>
    !#   <cardinality>0..</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>retryMaximum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>50</defaultValue>
    !#   <description>The number of tree build attempts to complete before rescaling the tolerance.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>rescaleMaximum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>20</defaultValue>
    !#   <description>The maximum allowed number of tolerance rescalings.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>attemptsMaximum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>10000</defaultValue>
    !#   <description>The maximum allowed number of tree build attempts.</description>
    !#   <type>logical</type>
    !#   <cardinality>0..</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massCutOffAttemptsMaximum</name>
    !#   <source>parameters</source>
    !#   <defaultValue>50</defaultValue>
    !#   <description>The number of trees with nodes above the mass resolution to allow before adjusting the mass cut-off tolerance.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massCutOffScaleFactor</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.05d0</defaultValue>
    !#   <description>The amount by which to increase the mass cut-off scale tolerance after exhausting tree build attempts.</description>
    !#   <type>real</type>
    !#   <cardinality>0..</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>useOneNodeTrees</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true, trees only consisting of their base node will be augmented.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    if (parameters%isPresent('snapshotRedshifts')) then
       allocate(timeSnapshots(parameters%count('snapshotRedshifts')))
       !# <inputParameter>
       !#   <name>snapshotRedshifts</name>
       !#   <variable>timeSnapshots</variable>
       !#   <source>parameters</source>
       !#   <description>For {\normalfont \ttfamily augment} description of redshift snapshots.</description>
       !#   <type>double precision</type>
       !#   <cardinality> 0..</cardinality>
       !# </inputParameter>
       cosmologyFunctions_ => cosmologyFunctions()
       do i=1,size(timeSnapshots)
          timeSnapshots(i)=cosmologyFunctions_ %cosmicTime                 (                  &
               &            cosmologyFunctions_%expansionFactorFromRedshift (                 &
               &                                                             timeSnapshots(i) &
               &                                                            )                 &
               &                                                           )
       end do
       augmentConstructorParameters=augmentConstructorInternal(massResolution,performChecks,toleranceScale,retryMaximum,rescaleMaximum,attemptsMaximum,massCutOffAttemptsMaximum,massCutOffScaleFactor,useOneNodeTrees,mergerTreeBuilder_,timeSnapshots)
    else
       augmentConstructorParameters=augmentConstructorInternal(massResolution,performChecks,toleranceScale,retryMaximum,rescaleMaximum,attemptsMaximum,massCutOffAttemptsMaximum,massCutOffScaleFactor,useOneNodeTrees,mergerTreeBuilder_              )
    end if
    return
  end function augmentConstructorParameters

  function augmentConstructorInternal(massResolution,performChecks,toleranceScale,retryMaximum,rescaleMaximum,attemptsMaximum,massCutOffAttemptsMaximum,massCutOffScaleFactor,useOneNodeTrees,mergerTreeBuilder_,timeSnapshots)
    !% Internal constructor for the {\normalfont \ttfamily augment} merger tree operator class.
    use Memory_Management
    use Cosmology_Functions
    use Sort
    implicit none
    type            (mergerTreeOperatorAugment)                                        :: augmentConstructorInternal
    double precision                           , intent(in   )                         :: massResolution                   , toleranceScale           , &
         &                                                                                massCutOffScaleFactor
    integer                                    , intent(in   )                         :: retryMaximum                     , rescaleMaximum           , &
         &                                                                                attemptsMaximum                  , massCutOffAttemptsMaximum
    double precision                           , intent(in   ), dimension(:), optional :: timeSnapshots
    class           (mergerTreeBuilderClass   ), intent(in   ), pointer                :: mergerTreeBuilder_
    logical                                    , intent(in   )                         :: performChecks                    , useOneNodeTrees
    class           (cosmologyFunctionsClass  )               , pointer                :: cosmologyFunctions_
    double precision                           , parameter                             :: expansionFactorDefault    =0.01d0

    cosmologyFunctions_ => cosmologyFunctions()
    if (present(timeSnapshots)) then
       call Alloc_Array(augmentConstructorInternal%timeSnapshots,shape(timeSnapshots))
       augmentConstructorInternal%timeSnapshots=timeSnapshots
       call Sort_Do(augmentConstructorInternal%timeSnapshots)
       augmentConstructorInternal%timeEarliest=min(                                                                  &
            &                                      cosmologyFunctions_       %cosmicTime   (expansionFactorDefault), &
            &                                      augmentConstructorInternal%timeSnapshots(                     1)  &
            &                                     )
    else
       augmentConstructorInternal%timeEarliest=    cosmologyFunctions_       %cosmicTime   (expansionFactorDefault)
    end if
    augmentConstructorInternal%massResolution            =  massResolution
    augmentConstructorInternal%performChecks             =  performChecks
    augmentConstructorInternal%toleranceScale            =  toleranceScale
    augmentConstructorInternal%retryMaximum              =  retryMaximum
    augmentConstructorInternal%rescaleMaximum            =  rescaleMaximum
    augmentConstructorInternal%attemptsMaximum           =  attemptsMaximum
    augmentConstructorInternal%massCutOffAttemptsMaximum =  massCutOffAttemptsMaximum
    augmentConstructorInternal%massCutOffScaleFactor     =  massCutOffScaleFactor
    augmentConstructorInternal%useOneNodeTrees           =  useOneNodeTrees
    augmentConstructorInternal%mergerTreeBuilder_        => mergerTreeBuilder_
     call augmentConstructorInternal%mergerTreeBuilder_%timeEarliestSet(augmentConstructorInternal%timeEarliest)    
   return
  end function augmentConstructorInternal

  subroutine augmentDestructor(self)
    !% Destructor for the augment merger tree operator function class.
    implicit none
    type(mergerTreeOperatorAugment), intent(inout) :: self

    !# <objectDestructor name="self%mergerTreeBuilder_" />
    return
  end subroutine augmentDestructor

  subroutine augmentOperate(self,tree)
    !% Augment the resolution of a merger tree by inserting high resolution branches.
    use Galacticus_Nodes
    use Galacticus_Display
    use String_Handling
    use Merger_Trees_Pruning_Utilities
    implicit none
    class           (mergerTreeOperatorAugment), intent(inout)                 :: self
    type            (mergerTree               ), intent(inout ), target        :: tree
    type            (treeNode                 ), pointer                       :: node
    class           (nodeComponentBasic       ), pointer                       :: basic
    type            (mergerTree               ), pointer                       :: treeCurrent
    type            (treeNodeList             ), allocatable   , dimension (:) :: anchorNodes
    type            (varying_string           )                                :: message
    type            (mergerTree               )                                :: treeBest
    integer                                                                    :: nodeCount                     , i                            , &
         &                                                                        retryCount                    , treeBuilt                    , &
         &                                                                        attemptsRemaining             , rescaleCount                 , &
         &                                                                        massCutoffAttemptsRemaining
    double precision                                                           :: tolerance                     , treeBestWorstFit             , &
         &                                                                        multiplier                    , constant                     , &
         &                                                                        scalingFactor                 , massCutoffScale
    logical                                                                    :: treeBestOverride              , treeNewHasNodeAboveResolution, &
         &                                                                        treeBestHasNodeAboveResolution, newRescale                   , &
         &                                                                        useOneNodeTrees
    
    ! Iterate over all linked trees in this forest.
    call Galacticus_Display_Indent('Augmenting merger tree',verbosityWorking)
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Allocate array of original anchor nodes from which new high-resolution branches will be built.
       nodeCount=augmentTreeStatistics(treeCurrent,treeStatisticNodeCount)
       if (nodeCount > 1 .or. self%useOneNodeTrees) then
          if (Galacticus_Verbosity_Level() >= verbosityWorking) then
             message="Number of nodes in tree: "
             message=message//nodeCount
             call Galacticus_Display_Message(message)
          end if
          allocate(anchorNodes(nodeCount))
          ! Build pointers to all anchor nodes.
          node => treeCurrent%baseNode
          do i=1,nodeCount
             anchorNodes(i)%node => node
             node                => node%walkTree()
          end do
          ! Walk the tree.
          i=1
          do while (i <= nodeCount)
             ! Get the node to work with.
             node                           => anchorNodes(i)%node
             basic                          => node          %basic()
             ! Initialize the current best-known tree to null.
             treeBestWorstFit               =      3.000d0
             treeBest%baseNode              => null()
             treeBestOverride               = .false.
             treeBestHasNodeAboveResolution = .false.
             newRescale                     = .false.
             ! Reset all factors used in tree acceptance and scaling.
             multiplier                     =      0.000d0
             constant                       =      0.000d0
             scalingFactor                  =      1.000d0
             ! Currently not used: call augmentScaleChildren(self, node, multiplier, constant, scalingFactor)
             tolerance                      =  self%toleranceScale
             rescaleCount                   =      0
             retryCount                     =      1
             treeBuilt                      =  treeBuildFailureGeneric
             attemptsRemaining              =  self%attemptsMaximum
             massCutoffScale                =      1.0d0
             massCutoffAttemptsRemaining    =  self%massCutOffAttemptsMaximum
             ! Begin building trees from this node, searching for an acceptable tree.
             if (Galacticus_Verbosity_Level() >= verbosityWorking) then
                message="Building tree from node: "
                message=message//node%index()
                call Galacticus_Display_Indent(message)
             end if
             do while (                                            &
                  &     treeBuilt         /= treeBuildSuccess      & ! Exit if tree successfully built.
                  &    .and.                                       &
                  &     rescaleCount      <= self%rescaleMaximum   & ! Exit once number of rescalings exceeds maximum allowed.
                  &    .and.                                       &
                  &     attemptsRemaining >  0                     & ! Exit if no more attempts remain.
                  &   )
                treeNewHasNodeAboveResolution=.false.
                treeBuilt                    =self%buildTreeFromNode(                                &
                     &                                               node                          , &
                     &                                               .false.                       , &
                     &                                               tolerance                     , &
                     &                                               self%timeEarliest             , &
                     &                                               treeBest                      , &
                     &                                               treeBestWorstFit              , &
                     &                                               treeBestOverride              , &
                     &                                               multiplier                    , &
                     &                                               constant                      , &
                     &                                               scalingFactor                 , &
                     &                                               massCutoffScale               , &
                     &                                               treeNewHasNodeAboveResolution , &
                     &                                               treeBestHasNodeAboveResolution, &
                     &                                               newRescale                      &
                     &                                              )
                !Check for exhaustion of retry attempts.
                if (retryCount == self%retryMaximum) then
                   ! Rescale the tolerance to allow less accurate tree matches to be accepted in future.
                   retryCount  =0
                   tolerance   =tolerance   *self%toleranceScale
                   ! Check for exhaustion of rescaling attempts.
                   rescaleCount=rescaleCount+1
                   if (rescaleCount > self%rescaleMaximum) call Galacticus_Display_Message('Node build attempts exhausted',verbosityWorking)
                   newRescale = .true.
                else 
                   newRescale = .false.
                end if
                ! Increment the retry count in cases where the tree was not acepted due to matching tolerance.
                select case (treeBuilt)
                case (treeBuildSuccess         )
                   retryCount=retryCount-1
                case (treeBuildFailureTolerance)
                   retryCount=retryCount+1
                end select
                ! Decrement the number of attempts remaining and, if the best tree is to be forcibly used, set no attempts remaining.
                attemptsRemaining                                                          =attemptsRemaining-1
                if (treeBestOverride .and. treeBuilt /= treeBuildSuccess) attemptsRemaining=attemptsRemaining+1
                ! If all attempts have been used but no match has been found, insert the best tree on next pass through loop.
                if (attemptsRemaining == 1 .and. associated(treeBest%baseNode)) treeBestOverride=.true.
                ! If the new tree contains nodes above the mass cut-off, decrement the number of remaining retries.
                if (treeNewHasNodeAboveResolution) then
                   massCutoffAttemptsRemaining=massCutoffAttemptsRemaining-1
                   ! If number of retries is exhausted, adjust the tolerance for declaring nodes to be above the mass cut-off.
                   if (massCutoffAttemptsRemaining == 0) then
                      massCutoffAttemptsRemaining=                self%massCutoffAttemptsMaximum
                      massCutoffScale            =massCutoffScale+self%massCutOffScaleFactor
                   end if
                end if
             end do
             ! Clean up the best tree if one exists.
             if (associated(treeBest%baseNode)) then
                call treeBest%destroyBranch(treeBest%baseNode)
                treeBest%baseNode => null() 
             end if
             ! Move on to the nest node.
             i=i+1
             call Galacticus_Display_Unindent('Finished building tree',verbosityWorking)
          end do
          deallocate(anchorNodes)
       end if
       !Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    call Galacticus_Display_Unindent('done',verbosityWorking)
    return
  end subroutine augmentOperate
  
  recursive integer function augmentBuildTreeFromNode(self,node,extendingEndNode,tolerance,timeEarliestIn,treeBest,treeBestWorstFit,treeBestOverride,multiplier,constant,scalingFactor,massCutoffScale,treeNewHasNodeAboveResolution,treeBestHasNodeAboveResolution, newRescale)
    use, intrinsic :: ISO_C_Binding
    use               Arrays_Search
    use               Galacticus_Nodes
    implicit none
    class           (mergerTreeOperatorAugment    ), intent(inout)         :: self
    type            (treeNode                     ), intent(inout), target :: node
    double precision                               , intent(in   )         :: timeEarliestIn               , tolerance
    double precision                               , intent(inout)         :: treeBestWorstFit             , multiplier                    , &
         &                                                                    constant                     , scalingFactor                 , &
         &                                                                    massCutoffScale
    logical                                        , intent(in   )         :: extendingEndNode
    logical                                        , intent(inout)         :: treeNewHasNodeAboveResolution, treeBestHasNodeAboveResolution, &
         &                                                                    newRescale
    type            (mergerTree                   ), intent(inout)         :: treeBest
    type            (treeNode                     ), pointer               :: baseNode
    class           (nodeComponentBasic           ), pointer               :: basic                        , baseBasic                     , &
         &                                                                    childBasic
    type            (mergerTree                   )                        :: newTree
    double precision                                                       :: timeEarliest
    integer         (c_size_t                     )                        :: timeIndex
    integer                                                                :: endNodeCount                 , nodeChildCount                , &
         &                                                                    i                            , treeAccepted
    type            (mergerTreeOperatorPruneByTime)                        :: pruneByTime
    logical                                                                :: treeBestOverride, newTreeBest

    ! Find the earliest time to which the tree should be built.
    basic => node%basic()
    if (extendingEndNode) then
       ! An end-node is being extended - we can build the tree all the way back to the earliest time.
       timeEarliest =  timeEarliestIn
    else if (associated(node%firstChild)) then
       ! The node has children - set the limit time to the time at which the children exist.
       childBasic   => node      %firstChild%basic()
       timeEarliest =  childBasic%           time ()
    else if (allocated(self%timeSnapshots)) then
       if (size(self%timeSnapshots) == 1) then
          ! Only one snapshot time is given, use the earliest time.
          timeEarliest =  timeEarliestIn
       else
          timeIndex=Search_Array_For_Closest(self%timeSnapshots,basic%time(),tolerance=1.0d-5)
          if (timeIndex == 1) then
             timeEarliest=timeEarliestIn
          else
             timeEarliest=self%timeSnapshots(timeIndex-1)
          end if
       end if
    else
       timeEarliest=timeEarliestIn
    end if
    ! Build trees from the nodes above the mass resolution. 
     ! Create a new base node, matched to the current node, build a tree from it, and truncate that tree to the desired earliest time.
     pruneByTime      =  mergerTreeOperatorPruneByTime(             timeEarliest       )
     baseNode         => treeNode                     (node%index(),newTree            )
     baseBasic        => baseNode%basic               (             autoCreate  =.true.)
     newTree%baseNode => baseNode
     call baseBasic  %                   timeSet        (basic       %time())
     call baseBasic  %                   massSet        (basic       %mass())
     call self       %mergerTreeBuilder_%timeEarliestSet(timeEarliest       ) 
     call self       %mergerTreeBuilder_%build          (newTree            )
     call pruneByTime%operate                           (newTree            )
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
          &                         timeEarliest                  , &
          &                         treeBest                      , &
          &                         treeBestWorstFit              , &
          &                         treeBestOverride              , &
          &                         multiplier                    , &
          &                         constant                      , &
          &                         scalingFactor                 , &
          &                         massCutoffScale               , &
          &                         treeNewHasNodeAboveResolution , &
          &                         treeBestHasNodeAboveResolution, &
          &                         newTreeBest                     &
          &                         )
     ! Determine whether to use stored best tree, newly created tree, or reject the tree.
     if     (                                          &
          &   associated(treeBest%baseNode)            &
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
        if (associated(newTree%baseNode)) call newTree%destroyBranch(newTree%baseNode)
        newTree%baseNode => treeBest%baseNode
        ! Reset the best tree.
        treeBestWorstFit  =  3.0d0
        treeBest%baseNode => null()
        ! Test for acceptance of the best tree.
        if     (                                                 &
             &   self%acceptTree(                                &
             &                   node                          , &
             &                   newTree                       , &
             &                   nodeChildCount                , &
             &                   extendingEndnode              , &
             &                   tolerance                     , &
             &                   timeEarliest                  , &
             &                   treeBest                      , &
             &                   treeBestWorstFit              , &
             &                   treeBestOverride              , &
             &                   multiplier                    , &
             &                   constant                      , &
             &                   scalingFactor                 , &
             &                   massCutoffScale               , &
             &                   treeNewHasNodeAboveResolution , &
             &                   treeBestHasNodeAboveResolution, &
             &                   newTreeBest                     &
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
        if (associated(treeBest%baseNode)) call treeBest%destroyBranch(treeBest%baseNode)
        treeBestWorstFit  =  3.0d0
        treeBest%baseNode => null()
     else
        ! The newly created tree was unacceptable, clean it up and return the failure code.
        if (associated(newTree%baseNode)) call newTree%destroyBranch(newTree%baseNode)
        augmentBuildTreeFromNode=treeAccepted
     end if
    return
  end function augmentBuildTreeFromNode

  recursive integer function augmentAcceptTree(self,node,tree,nodeChildCount,extendingEndNode,tolerance,timeEarliest,treeBest,treeBestWorstFit,treeBestOverride,multiplier,constant,scalingFactor,massCutoffScale,treeNewHasNodeAboveResolution,treeBestHasNodeAboveResolution, newTreeBest)
    use Galacticus_Nodes
    use Galacticus_Display
    use Merger_Trees_Builders
    implicit none
    class           (mergerTreeOperatorAugment), intent(inout)                      :: self
    type            (treeNode                 ), intent(inout)            , target  :: node
    type            (treeNode                 )                           , pointer :: nodeCurrent                  , nodePrevious                  , &
         &                                                                             nodeNonOverlap               , nodeNonOverlapFirst
    class           (nodeComponentBasic       )                           , pointer :: basicCurrent                 , basicSort                     , &
         &                                                                             basicNonOverlap
    type            (mergerTree               ), intent(inout)            , target  :: tree
    type            (mergerTree               ), intent(inout)            , target  :: treeBest
    logical                                    , intent(inout)                      :: treeNewHasNodeAboveResolution, treeBestHasNodeAboveResolution, &
         &                                                                             newTreeBest                  
    logical                                    , intent(in   )                      :: treeBestOverride             , extendingEndNode
    double precision                           , intent(inout)                      :: treeBestWorstFit             , multiplier                    , &
         &                                                                             constant                     , scalingFactor                 , &
         &                                                                             massCutoffScale
    double precision                           , intent(in   )                      :: tolerance                    , timeEarliest
    integer                                    , intent(inout)                      :: nodeChildCount
    type            (treeNodeList             ), dimension(nodeChildCount)          :: endNodes
    integer                                                                         :: i                            , j                             , &
         &                                                                             endNodesSorted
    logical                                                                         :: treeAccepted                 , nodeMassesAgree               , &
         &                                                                             nodeCurrentBelowAll          , treeScalable
    double precision                                                                :: unresolvedMass               , treeMass                      , &
         &                                                                             endNodeMass                  , treeCurrentWorstFit
    type            (varying_string           )                                     :: message
    character       (len=12                   )                                     :: label

    ! Initialize.
    treeNewHasNodeAboveResolution =  .false.
    nodeNonOverlapFirst           => null()
    i                             =  1
    endNodesSorted                =  1
    nodeCurrent                   => tree        %baseNode
    basicCurrent                  => nodeCurrent %basic   ()
    unresolvedMass                =  basicCurrent%mass    ()
    treeMass                      =  basicCurrent%mass    ()
    ! Walk through the tree identifying end-nodes.
    do while (associated(nodeCurrent))
       ! Initialize the current node.
       basicCurrent         => nodeCurrent%basic   ()
       nodeCurrent%hostTree => node       %hostTree
       nodeCurrent%event    => null()
       nodeNonOverlap       => null()
       ! Test for children.
       if (.not.associated(nodeCurrent%firstChild)) then
          ! This is an end-node. Remove its mass from the current unresolved mass.
          unresolvedMass=unresolvedMass-basicCurrent%mass()
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
                         if (basicSort%mass() > self%massResolution*massCutoffScale) then
                            if (Galacticus_Verbosity_Level() >= verbosityWorking) then
                               write (label,'(e12.6)') basicSort%mass()                      
                               message="Nonoverlap failure at mass: "//trim(label)
                               call Galacticus_Display_Message(message)
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
                      if (basicCurrent%mass() > self%massResolution*massCutoffScale) then
                         if (Galacticus_Verbosity_Level() >= verbosityWorking) then
                            write (label,'(e12.6)') basicCurrent%mass()                      
                            message="Nonoverlap failure at mass: "//trim(label)
                            call Galacticus_Display_Message(message)
                         end if
                         treeNewHasNodeAboveResolution=.true.
                      end if
                   end if
                end if
             end if
          else
             ! No overlap nodes are being sought - so this node is automatically a non-overlap node.
             nodeNonOverlap => nodeCurrent
          end if
          endNodesSorted=endNodesSorted+1 
       end if
       ! Walk to the next node.
       nodeCurrent => nodeCurrent%walkTree()
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
          if (.not.augmentNodeComparison(nodeCurrent,endNodes(i)%node,2.0d0-2.0d0*tolerance,treeCurrentWorstFit)) nodeMassesAgree=.false.
          i           =  i+1
          nodeCurrent => nodeCurrent%sibling
       end do
    end if
    ! AJB : not currently used?
    !if (treeNewHasNodeAboveResolution) then
    !  treeNewHasNodeAboveResolution = augmentScaleNodesAboveCutoff(self, node, tree, endNodes, nodeChildCount, nodeNonOverlapFirst)
    !end if 
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
            &                          tree               , &
            &                          nodeNonOverlapFirst, &
            &                          tolerance          , &
            &                          timeEarliest       , &
            &                          treeBest           , &
            &                          massCutoffScale      &
            &                         )
       endNodeMass=0
       i          =1
       do while (i <= nodeChildCount)
          basicCurrent =>             endNodes    (i)%node%basic()
          endNodeMass  =  endNodeMass+basicCurrent   %     mass ()
          i            =  i          +                     1
       end do
       treeScalable=self%multiScale(                     &
            &                       node               , &
            &                       tree               , &
            &                       endNodes           , &
            &                       nodeChildCount     , &
            &                       nodeNonOverlapFirst, &
            &                       tolerance          , &
            &                       timeEarliest       , &
            &                       treeBest           , &
            &                       treeBestWorstFit   , &
            &                       unresolvedMass     , &
            &                       treeMass           , &
            &                       massCutoffScale      &
            &                      )
    else
       ! Tree is not scalable either because of structural difference (mismatch in number of overlap nodes), or because overlap
       ! nodes differ significantly in mass.
       treeScalable=.false.
    end if
    ! Determine if the tree is accepted.
    if (treeBestOverride .and. associated(tree%baseNode)) then
       treeAccepted=      nodeMassesAgree                                 &
            &       .and.                                                 &
            &             nodeChildCount                <= endNodesSorted
    else
       treeAccepted=      nodeMassesAgree                                 &
            &       .and.                                                 &
            &             nodeChildCount                <= endNodesSorted &
            &       .and.                                                 &
            &        .not.treeNewHasNodeAboveResolution                   &
            &       .and.                                                 &
            &             treeScalable
    end if
    ! Process the tree depending on acceptance state.
    newTreeBest = .false.    
    if (treeAccepted) then
       ! Tree was accepted, unscale and insert it into the original tree.
       call augmentResetUniqueIDs (tree)
       call augmentUnscaleChildren(self,node,nodeChildCount,endNodes,multiplier,constant,scalingFactor)
       call augmentResetUniqueIDs (tree)  
       call augmentSimpleInsert   (self,node,tree,endNodes,nodeChildCount,nodeNonOverlapFirst)
       !call augmentSimpleScale(self, node, tree, endNodes, nodeChildCount, nodeNonOverlapFirst)
    else if (nodeChildCount <= endNodesSorted .and. .not.treeBestOverride) then
       ! Tree is structurally acceptable - decide if we want to keep it as the current best tree.
       call self%nonOverlapReinsert(nodeNonOverlapFirst)
       call augmentResetUniqueIDs(tree)
       if (treeCurrentWorstFit < treeBestWorstFit) then
          newTreeBest = .true.
          ! Current tree is better than the current best tree. Replace the best tree with the current tree.
          if (associated(treeBest%baseNode)) call treeBest%destroyBranch(treeBest%baseNode)
          treeBest                      %baseNode          => tree                         %baseNode
          tree                          %baseNode          => null()
          treeBest                      %baseNode%hostTree => treeBest
          treeBestWorstFit                                 =  treeCurrentWorstFit
          treeBestHasNodeAboveResolution                   =  treeNewHasNodeAboveResolution
          nodeCurrent                                      => treeBest                     %baseNode
          do while (associated(nodeCurrent))
             nodeCurrent%event    => null()
             nodeCurrent%hostTree => treeBest  %baseNode %hostTree
             nodeCurrent          => nodeCurrent%walkTree         ()
          end do
       end if
    else
       ! Tree is not acceptable or better than the current best tree - destroy it.
       call self%nonOverlapReinsert(nodeNonOverlapFirst)
       if (associated(tree%baseNode)) call tree%destroyBranch(tree%baseNode)
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
    !% Insert a newly constructed tree into the original tree.
    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)                            :: self
    type   (treeNode                 ), intent(inout)                            :: node
    type   (mergerTree               ), intent(in   ), target                    :: tree
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
         &                      tree%baseNode             , &
         &                      keepTop           =.false., &
         &                      exchangeProperties=.false.  &
         &                     )
    return
  end subroutine augmentSimpleInsert
  
  subroutine augmentExtendByOverlap(nodeBottom,nodeTop,keepTop,exchangeProperties)
    !% Conjoin two trees by overlapping the {\normalfont \ttfamily nodeTop} of one tree with the chosen {\normalfont \ttfamily
    !% nodeBottom} of the other. If {\normalfont \ttfamily keepTop} is {\normalfont \ttfamily true}, {\normalfont \ttfamily
    !% nodeTop} replaces {\normalfont \ttfamily nodeBottom}, otherwise, {\normalfont \ttfamily nodeBottom} replaces {\normalfont
    !% \ttfamily nodeTop}. If {\normalfont \ttfamily exchangeProperties} is {\normalfont \ttfamily true}, the mass and time
    !% information of the deleted node overwrites the mass and time of the retained node.
    use Galacticus_Nodes
    implicit none
    type   (treeNode          ), intent(inout), target  :: nodeBottom  , nodeTop 
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
    end if
    return
  end subroutine augmentExtendByOverlap

  integer function augmentChildCount(node)
    !% Return a count of the number of child nodes.
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
    !% Sort the children of the given {\normalfont \ttfamily node} such that they are in descending mass order.
    use Galacticus_Error
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
          if (basicCurrent%mass() > massLargest) call Galacticus_Error_Report('augmentSortChildren','failed to sort child nodes')
          massLargest =  basicCurrent%mass   ()
          nodeCurrent => nodeCurrent %sibling
       end do
    end if
    return
  end subroutine augmentSortChildren

  subroutine augmentNonOverlapListAdd(node, listFirstElement)
    !% Add the given node to a linked list of non-overlap nodes in the current trial tree.
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
    !% Reinsert non-overlap nodes into their tree.
    implicit none
    class(mergerTreeOperatorAugment), intent(in   )          :: self
    type (treeNode                 ), intent(inout), pointer :: listFirstElement
    type (treeNode                 )               , pointer :: nodeCurrent

    do while (associated(listFirstElement))
       nodeCurrent      => listFirstElement
       listFirstElement => listFirstElement%sibling
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

  subroutine augmentExtendNonOverlapNodes(self,tree,nodeNonOverlapFirst,tolerance,timeEarliest,treeBest,massCutoffScale)
    !% Extend any non-overlap nodes in an accepted tree by growing a new tree from each such node.
    use Galacticus_Error
    use Galacticus_Nodes
    implicit none
    class           (mergerTreeOperatorAugment), intent(inout)          :: self
    type            (mergerTree               ), intent(inout), target  :: tree                   , treeBest
    type            (treeNode                 )               , pointer :: nodeNonOverlapFirst
    double precision                           , intent(in   )          :: tolerance              , timeEarliest
    double precision                           , intent(inout)          :: massCutoffScale 
    type            (treeNode                 )               , pointer :: nodeCurrent            , nodeNonOverlap
    double precision                                                    :: falseWorstFit          , falseMultiplier             , &
         &                                                                 falseConstant          , falseScalingFactor
    logical                                                             :: falseNewNodeAboveCutoff, falseBestTreeNodeAboveCutoff, &
         &                                                                 falseNewRescale
    integer                                                             :: treeStatus
    
    falseBestTreeNodeAboveCutoff =  .false.
    falseNewNodeAboveCutoff      =  .false.
    falseNewRescale              =  .false.
    falseMultiplier              =  0.0d0
    falseConstant                =  0.0d0
    falseScalingFactor           =  1.0d0
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
            &                            falseMultiplier             , &
            &                            falseConstant               , &
            &                            falseScalingFactor          , &
            &                            massCutoffScale             , &
            &                            falseNewNodeAboveCutoff     , &
            &                            falseBestTreeNodeAboveCutoff, &
            &                            falseNewRescale               &
            &                           )
       if (treeStatus /= treeBuildSuccess) call Galacticus_Error_Report('augmentExtendNonOverlapNodes','extension of non-overlap node failed')
       nodeCurrent => nodeCurrent%sibling
    end do
    return
  end subroutine augmentExtendNonOverlapNodes
  
  logical function augmentNodeComparison(nodeNew, nodeOriginal, tolerance, treeCurrentWorstFit)
    !% Compare the masses of an overlap node and the corresponding node in the original tree, testing if they agree to withing
    !% tolerance.
    implicit none   
    type            (treeNode          ), intent(in   ), pointer :: nodeNew            , nodeOriginal
    double precision                    , intent(in   )          :: tolerance
    double precision                    , intent(inout)          :: treeCurrentWorstFit
    class           (nodeComponentBasic)               , pointer :: basicNew           , basicOriginal
    double precision                                             :: thisFit

    basicNew              => nodeNew     %basic()
    basicOriginal         => nodeOriginal%basic()
    thisFit               =  +2.0d0                                     &
         &                   *abs(basicNew%mass()-basicOriginal%mass()) &
         &                   /abs(basicNew%mass()+basicOriginal%mass())
    augmentNodeComparison =  thisFit < tolerance
    treeCurrentWorstFit   =  max(treeCurrentWorstFit,thisFit)
    return
  end function augmentNodeComparison
  
  logical function augmentMultiScale(self,node,tree,endNodes,nodeChildCount,nodeNonOverlapFirst,tolerance,timeEarliest,treeBest,treeBestWorstFit,unresolvedMass,treeMass,massCutoffScale)
    !% 
    use Galacticus_Nodes
    use Merger_Trees_Builders
    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self
    type(treeNode), intent(in   ) :: node
    type (treeNode), pointer :: nodeCurrent, nodePrevious, nodeNonOverlapFirst, nodeCurrentSibling, currentChildNode, currentGrandchild, scaleNode
    type (mergerTree), target :: tree
    type (mergerTree), intent (inout), target :: treeBest
    class (nodeComponentBasic), pointer :: childComponentBasic, basicCurrent, basicSort
    integer, intent(inout) :: nodeChildCount
    integer :: i, retryCount
    type (treeNodeList), dimension(nodeChildCount), intent(inout) :: endNodes
    logical :: treeScaled, lastNodeFound
    double precision, intent(in   ) :: tolerance, timeEarliest
    double precision, intent (inout) :: treeBestWorstFit, unresolvedMass, treeMass, massCutoffScale
    double precision :: massExcess, childNodeMass, currentMass, endNodeMass, falseWorstFit, massDifferenceScaleFactor
  
    falseWorstFit = 3.0
    call self%nonOverlapReinsert(nodeNonOverlapFirst)
    if (nodeChildCount > 0) then
      nodeCurrent => tree%baseNode
      basicCurrent => nodeCurrent%basic()
      i = 1
      currentChildNode => node%firstChild
      massExcess = 0
      endNodeMass = 0
      childNodeMass = 0
      treeMass = basicCurrent%mass()
      do while (associated(currentChildNode))
        nodeCurrent => endNodes(i)%node
        basicCurrent => nodeCurrent%basic()
        childComponentBasic => currentChildNode%basic()
        massExcess = massExcess + basicCurrent%mass() - childComponentBasic%mass()
        endNodeMass = endNodeMass + basicCurrent%mass()
        childNodeMass = childNodeMass + childComponentBasic%mass()
        call nodeCurrent%uniqueIDSet(-nodeCurrent%uniqueID())
        i = i + 1
        currentChildNode => currentChildNode%sibling
      end do
      treeScaled = .true.
      if ( treeScaled) then
        i = 1
        currentChildNode => node%firstChild
        do while (associated(currentChildNode))
          nodeCurrent => endNodes(i)%node
          basicCurrent => nodeCurrent%basic()
          call augmentInsertChildMass(self, tree%baseNode, currentChildNode, nodeCurrent)
          currentChildNode => currentChildNode%sibling
        end do
        if (tree%baseNode%uniqueID() > 0) then
          call tree%baseNode%uniqueIDSet(-tree%baseNode%uniqueID())
        end if
   
      end if
    end if
    treescaled = .true.
    augmentMultiScale = treeScaled
    call augmentResetUniqueIDs(tree)
  end function augmentMultiScale

  subroutine augmentInsertChildMass(self, node, originalChildNode, newChildNode)
    use Galacticus_Nodes
    implicit none
    class (mergerTreeOperatorAugment), intent(inout) :: self
    type (treeNode), pointer :: node, originalChildNode, newChildNode, scaleNode
    class (nodeComponentBasic), pointer :: basicCurrent, childComponentBasic
    double precision :: multiplier, constant, scalingFactor
    double precision :: parentMass, childMass, parentTime, childTime, massDifference

    basicCurrent => node%basic()
    parentMass = basicCurrent%mass()
    parentTime = basicCurrent%time()
    childComponentBasic => originalChildNode%basic()
    childTime = childComponentBasic%time()
    massDifference = childComponentBasic%mass()
    childComponentBasic => newChildNode%basic()
    massDifference = massDifference - childComponentBasic%mass()
    multiplier = (massDifference)/(LOG10(childTime) - LOG10(parentTime))
    constant = -multiplier *LOG10(parentTime)
    scaleNode => newChildNode
    do while (associated(scaleNode))
        basicCurrent => scaleNode%basic()
     !! AJB HACK   call basicCurrent%massSet(basicCurrent%mass() + multiplier*LOG10(basicCurrent%time()) + constant)
        scaleNode => scaleNode%parent
    end do

  end subroutine augmentInsertChildMass

  integer function augmentTreeStatistics(tree,desiredOutput)
    !% Walks through tree and quietly collects information specified by {\normalfont \ttfamily desiredOutput} input enumeration and
    !% returns that information.
    use Galacticus_Nodes
    use Galacticus_Error
    implicit none
    type   (mergerTree), intent(in   ), target  :: tree
    integer            , intent(in   )          :: desiredOutput
    type   (treeNode  )               , pointer :: node
    integer                                     :: nodeCount    , endNodeCount
    
    nodeCount    =  0
    endNodeCount =  0
    node         => tree%baseNode
    do while (associated(node))
       nodeCount                                         =   nodeCount+1
       if (.not.associated(node%firstChild)) endNodeCount=endNodeCount+1
       node => node%walkTree()
    end do
    ! Return the requested quantity.
    select case (desiredOutput)
    case (treeStatisticNodeCount   )
       augmentTreeStatistics=nodeCount
    case (treeStatisticEndNodeCount)
       augmentTreeStatistics=endNodeCount
    case default
       call Galacticus_Error_Report('augmentTreeStatistics','unknown task requested')
    end select
    return
  end function augmentTreeStatistics

  subroutine augmentResetUniqueIDs(tree)
    !% Walk through a given {\normalfont \ttfamily tree} and reset all negative unique IDs back to positive.
    use Galacticus_Nodes
    implicit none
    type(mergerTree), intent(in   ), target  :: tree
    type(treeNode  )               , pointer :: node
 
    node => tree%baseNode
    do while (associated(node)) 
       if (node%uniqueID() < 0) call node%uniqueIDSet(-node%uniqueID())
       node => node%walkTree()
    end do
    return
  end subroutine augmentResetUniqueIDs
  
  subroutine augmentUnscaleChildren (self, node, nodeChildCount, endNodes, multiplier, constant, scalingFactor)
    use Galacticus_Nodes
    implicit none
    class (mergerTreeOperatorAugment), intent(inout) :: self
    type(treeNode), intent(inout) :: node
    type (treeNode), pointer :: currentChild, nodeCurrent
    integer, intent (in   ) :: nodeChildCount
    class (nodeComponentBasic), pointer :: basicCurrent, childComponentBasic
    type (treeNodeList), dimension(nodeChildCount), intent(inout) :: endNodes
    double precision, intent(inout) :: multiplier, constant, scalingFactor
    double precision :: parentMass, childMass, parentTime, childTime
    integer :: i
    
    !Unscales children using parameters saved from augmentScaleChildren function
    if (scalingFactor /= 1.0 .and. multiplier /= 0.0) then
      currentChild => node%firstChild
      childComponentBasic => currentChild%basic()
      childMass = 0
      do while(associated(currentChild)) 
        childComponentBasic => currentChild%basic()
        childMass = childMass + childComponentBasic%mass()
        call childComponentBasic%massSet(childComponentBasic%mass()/scalingFactor)
        currentChild => currentChild%sibling
      end do
      i = 1
      basicCurrent => node%basic()
      parentMass = basicCurrent%mass()
      currentChild => node%firstChild
      do while (associated(currentChild))
        nodeCurrent => endNodes(i)%node
        do while (associated(nodeCurrent))
          basicCurrent = nodeCurrent%basic()
          if (nodeCurrent%uniqueID() > 0) then
            call basicCurrent%massSet((basicCurrent%mass() / parentMass)*(multiplier*LOG10(basicCurrent%time()) + constant))
            call nodeCurrent%uniqueIDSet(-nodeCurrent%uniqueID())
          end if 
          nodeCurrent => nodeCurrent%parent
        end do
        i = i + 1
        currentChild => currentChild%sibling
      end do
    end if  

  end subroutine augmentUnscaleChildren


  !!! CURRENTLY UNUSED !!!

  subroutine augmentSimpleScale(self, node, tree, endNodes, nodeChildCount, nodeNonOverlapFirst)
    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self
    type (treeNode), pointer :: node, nodeCurrent, nodePrevious, nodeNonOverlapFirst, nodeCurrentSibling, nodeSort
    type (mergerTree), target :: tree
    integer, intent(inout) :: nodeChildCount
    integer :: i
    type (treeNodeList), dimension(nodeChildCount), intent(inout) :: endNodes
    logical :: scalableTest
    double precision :: massDifference, nonOverlapMass, scaleFactor
    class (nodeComponentBasic), pointer :: basicCurrent, childComponentBasic

    write (*,*) "Starting Simple Scale"
    massDifference = 0
    i = 1
    if (nodeChildCount > 0) then
      nodeCurrent => node%firstChild
      do while (associated(nodeCurrent))
        basicCurrent => nodeCurrent%basic()
        childComponentBasic => endNodes(i)%node%basic()
        massDifference = massDifference + basicCurrent%mass() - childComponentBasic%mass()
        nodeSort => endNodes(i)%node
        do while (associated(nodeSort))
          call nodeSort%uniqueIDSet(-nodeSort%uniqueID())
          nodeSort => nodeSort%parent
        end do 
        i = i + 1
        nodeCurrent => nodeCurrent%sibling
      end do
      nodeCurrent => nodeNonOverlapFirst
      do while(associated(nodeCurrent))
        basicCurrent => nodeCurrent%basic()
        nonOverlapMass = nonOverlapMass + basicCurrent%mass()
        nodeCurrent => nodeCurrent%sibling
      end do
      if (.not. (nonOverlapMass ==0)) then
        scaleFactor = (nonOverlapMass - massDifference)/nonOverlapMass
      else 
        scaleFactor = 1
      end if
      nodeCurrent => nodeNonOverlapFirst
      do while(associated(nodeCurrent))
        nodeSort => nodeCurrent
        do while(associated(nodeSort))
          basicCurrent => nodeSort%basic()
          if(nodeSort%uniqueID() > 0) then
            call basicCurrent%massSet(scaleFactor * basicCurrent%mass())
          end if
          nodeSort => nodeSort%parent
        end do
        nodeCurrent => nodeCurrent%sibling
      end do
    end if
    call self%nonOverlapReinsert(nodeNonOverlapFirst)

    if (nodeChildCount > 0) then
      i =1
      nodeCurrent => node%firstChild

      do while (associated(nodeCurrent))
        nodeCurrentSibling => nodeCurrent%sibling
        node%firstChild => nodeCurrent%sibling
        call augmentExtendByOverLap(endNodes(i)%node, nodeCurrent,.true., .false.)
        nodeCurrent => nodeCurrentSibling
        i = i + 1
      end do
    else
      nodeCurrent => tree%baseNode
      do while (associated(nodeCurrent))
           nodeCurrent%event => null()
           nodeCurrent%hosttree => node%hostTree
           nodePrevious => nodeCurrent%parent
           nodeCurrent => nodeCurrent%walkTree()
      end do

    end if

    call augmentExtendByOverLap(node, tree%baseNode, .false., .false.)

  end subroutine augmentSimpleScale

  subroutine augmentScaleBranch (self, node, scalingFactor)

    use Galacticus_Nodes
    implicit none
    class (mergerTreeOperatorAugment), intent (inout) :: self
    type (treeNode), pointer :: node, currentChild
    class (nodeComponentBasic), pointer :: thisComponentBasic
    double precision :: scalingFactor
    thisComponentBasic => node%basic()
    call thisComponentBasic%massSet(thisComponentBasic%mass() * scalingFactor)
    currentChild => node%firstChild
    do while (associated(currentChild)) 
      call augmentScaleBranch (self, currentChild, scalingFactor)
      currentChild => currentChild%sibling
    end do

  end subroutine augmentScaleBranch


  subroutine augmentRemoveUnresolvedMass (tree)
    use Galacticus_Nodes
    implicit none 
    type (mergerTree), target :: tree
    class (nodeComponentBasic), pointer :: basicCurrent, childComponentBasic
    type (treeNode), pointer :: nodeCurrent, currentChild, nodePrevious
    double precision :: childMass, unresolvedMass

    nodeCurrent => tree%baseNode
    basicCurrent => nodeCurrent%basic()
    do while (associated (nodeCurrent))
      basicCurrent => nodeCurrent%basic()
      childMass = 0
      unresolvedMass = 0
      currentChild => nodeCurrent%firstChild
      do while(associated(currentChild))
        childComponentBasic => currentChild%basic()
        childMass = childMass + childComponentBasic%mass()
        currentChild => currentChild%sibling
      end do
      if (childMass > 0) then
        unresolvedMass = basicCurrent%mass() - childMass
        call basicCurrent%massSet(basicCurrent%mass() - unresolvedMass)
      end if 
      nodePrevious => nodeCurrent%parent
      nodeCurrent => nodeCurrent%walkTree()
    end do

  end subroutine augmentRemoveUnresolvedMass

  logical function augmentScaleNodesAboveCutoff(self, node, tree, endNodes, nodeChildCount, nodeNonOverlapFirst)
    use Galacticus_Nodes
    implicit none
    class (mergerTreeOperatorAugment), intent (inout) :: self
    type (treeNode), pointer :: node, nodeNonOverlapFirst, currentNonOverlap, currentChild, currentEndNode, scaleNode
    type (mergerTree), target :: tree
    integer, intent(inout) :: nodeChildCount
    type (treeNodeList), dimension(nodeChildCount), intent(inout) :: endNodes
    class (nodeComponentBasic), pointer :: basicCurrent, childComponentBasic
    double precision :: massAboveCutoff, overlapMassDifference, massLargestAboveCutoff, scalingFactor, totalEndNodeMass, scalingMass
    logical :: nodesScalable
    integer :: i

    nodesScalable = .true.
    overlapMassDifference = 0.0
    totalEndNodeMass = 0.0
    i = 1
    currentChild => node%firstChild
    do while (associated(currentChild))
      currentEndNode => endNodes(i)%node
      basicCurrent => currentEndNode%basic()
      childComponentBasic => currentChild%basic()
      overlapMassDifference = overlapMassDifference + basicCurrent%mass() - childComponentBasic%mass()
      totalEndNodeMass = totalEndNodeMass + basicCurrent%mass()
      !do while (associated(currentEndNode))
      !  call currentEndNode%uniqueIDSet(-currentEndNode%uniqueID())
      !  currentEndNode => currentEndNode%parent
      !end do
      i = i + 1
      currentChild => currentChild%sibling
    end do

    massAboveCutoff = 0.0
    massLargestAboveCutoff = 0.0
    currentNonOverlap => nodeNonOverlapFirst
    do while (associated(currentNonOverlap))
      basicCurrent => currentNonOverlap%basic()
      if (basicCurrent%mass() > self%massResolution) then
        massAboveCutoff = massAboveCutoff + basicCurrent%mass()
      end if 
      if (basicCurrent%mass() > massLargestAboveCutoff) then
        massLargestAboveCutoff = basicCurrent%mass()
      end if
      currentNonOverlap => currentNonOverlap%sibling
    end do

    if (overlapMassDifference < 0) then
      nodesScalable = .false.
    else if (overlapMassDifference < massAboveCutoff) then
      if ( massLargestAboveCutoff * (( massAboveCutoff - overlapMassDifference) /  massAboveCutoff) > self%massResolution) then 
        nodesScalable = .false.
      else 
        currentNonOverlap => nodeNonOverlapFirst
        do while (associated(currentNonOverlap))
        basicCurrent => currentNonOverlap%basic()
        if ( basicCurrent%mass() > self%massResolution) then
          scalingMass = basicCurrent%mass() * ( 1 - (massAboveCutoff - overlapMassDifference) / massAboveCutoff)
          scaleNode => currentNonOverlap
          do while (associated(scaleNode))
            basicCurrent => scaleNode%basic()
            call basicCurrent%massSet(basicCurrent%mass() - scalingMass)
            scaleNode => scaleNode%parent
          end do
        end if
        currentNonOverlap => currentNonOverlap%sibling
        end do

        i = 1
        currentChild => node%firstChild
        do while (associated(currentChild))
          currentEndNode => endNodes(i)%node
          basicCurrent => currentEndNode%basic()
          childComponentBasic => currentChild%basic()
          scalingMass = childComponentBasic%mass() - basicCurrent%mass()
          scaleNode => currentEndNode
          do while (associated(scaleNode))
            basicCurrent => scaleNode%basic()
            call basicCurrent%massSet(basicCurrent%mass() + scalingMass)
            scaleNode => scaleNode%parent
          end do
          currentChild => currentChild%sibling
          i = i + 1
        end do
      end if
   
    else 
      nodesScalable = .false.
    end if
     
    !call augmentResetUniqueIDs(tree)
    augmentScaleNodesAboveCutoff = nodesScalable
  end function augmentScaleNodesAboveCutoff



  subroutine augmentScaleChildren(self, node, multiplier, constant, scalingFactor)
    use Galacticus_Nodes
    implicit none
    class (mergerTreeOperatorAugment), intent(inout) :: self
    type (treeNode), pointer :: node, currentChild
    class (nodeComponentBasic), pointer :: basicCurrent, childComponentBasic
    double precision, intent(inout) :: multiplier, constant, scalingFactor
    double precision :: parentMass, childMass, parentTime, childTime

    !If the sum of the child masses of node are above the mass of node
    !ScaleChildren will scale down their masses to equal the mass of node
    !keeping the scaling factor in scalingFactor.
    !The additional mass will be added in logaritmically, so that the total mass
    !of the descendants of each node will follow the relation:
    !Mass = multiplier * log(time) + constant 
    !At time = time_node, Mass = node's mass.
    !At time = time_node%firstChild, Mass = sum of node's childrens' mass.
    basicCurrent => node%basic()
    parentMass = basicCurrent%mass()
    parentTime = basicCurrent%time()
    childMass = 0
    currentChild => node%firstChild
    if (associated(currentChild)) then
      do while (associated(currentChild))
        childComponentBasic => currentChild%basic()
        childMass = childMass + childComponentBasic%mass()
        currentChild => currentChild%sibling
      end do 
      currentChild => node%firstChild
      childComponentBasic => currentChild%basic()
      childTime = childComponentBasic%time()
      if (childMass > parentMass) then
        scalingFactor = parentMass/childMass
      else
        scalingFactor = 1.0
      end if
      constant = parentMass
      multiplier = 0.0
      if (scalingFactor /= 1.0) then
        currentChild => node%firstChild 
        childComponentBasic => currentChild%basic()
        do while (associated(currentChild)) 
          basicCurrent => currentChild%basic()
          call basicCurrent%massSet(scalingFactor*basicCurrent%mass())
          currentChild => currentChild%sibling
        end do
        multiplier = (childMass - parentMass)/ (LOG10(childTime) - LOG10(parentTime))
        constant = parentMass - multiplier* LOG10(parentTime)
      end if
    else
      scalingFactor = 1.0
      constant = parentMass
      multiplier = 0.0
    end if 

  end subroutine augmentScaleChildren

