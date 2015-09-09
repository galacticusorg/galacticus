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
     double precision                                                     :: massResolution    , timeEarliest
     class           (mergerTreeBuilderClass), pointer                    :: mergerTreeBuilder_
   contains
     !@ <objectMethods>
     !@   <object>mergerTreeOperatorAugment</object>
     !@   <objectMethod>
     !@     <method>buildTreeFromNode</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Build a merger tree starting from the given node.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                      augmentDestructor
     procedure :: operate           => augmentOperate
     procedure :: buildTreeFromNode => augmentBuildTreeFromNode
  end type mergerTreeOperatorAugment

  interface mergerTreeOperatorAugment
     !% Constructors for the {\normalfont \ttfamily augment} merger tree operator class.
     module procedure augmentConstructorParameters
     module procedure augmentConstructorInternal
  end interface mergerTreeOperatorAugment

  !# <enumeration>
  !#  <name>walkTask</name>
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
    double precision                                                         :: massResolution
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
    if (parameters%isPresent('snapshotRedshifts')) then
       allocate(timeSnapshots(parameters%count('snapshotRedshifts')))
    else
       call Galacticus_Error_Report('augmentConstructorParameters','parameter [snapshotRedshifts] is required')
    end if
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
    augmentConstructorParameters=augmentConstructorInternal(timeSnapshots,massResolution,mergerTreeBuilder_)
    return
  end function augmentConstructorParameters

  function augmentConstructorInternal(timeSnapshots,massResolution,mergerTreeBuilder_)
    !% Internal constructor for the {\normalfont \ttfamily augment} merger tree operator class.
    use Memory_Management
    use Cosmology_Functions
    use Sort
    implicit none
    type            (mergerTreeOperatorAugment)                              :: augmentConstructorInternal
    double precision                           , intent(in   )               :: massResolution
    double precision                           , intent(in   ), dimension(:) :: timeSnapshots
    class           (mergerTreeBuilderClass   ), intent(in   ), pointer      :: mergerTreeBuilder_
    class           (cosmologyFunctionsClass  )               , pointer      :: cosmologyFunctions_
    double precision                           , parameter                   :: expansionFactorDefault    =0.01d0

    call Alloc_Array(augmentConstructorInternal%timeSnapshots,shape(timeSnapshots))
    augmentConstructorInternal%timeSnapshots      =  timeSnapshots
    augmentConstructorInternal%massResolution     =  massResolution
    augmentConstructorInternal%mergerTreeBuilder_ => mergerTreeBuilder_
    cosmologyFunctions_                           => cosmologyFunctions()
    call Sort_Do(augmentConstructorInternal%timeSnapshots)
    augmentConstructorInternal%timeEarliest=min(                                                                  &
         &                                      cosmologyFunctions_       %cosmicTime   (expansionFactorDefault), &
         &                                      augmentConstructorInternal%timeSnapshots(                     1)  &
         &                                     )
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
    integer                                                                    :: nodeCount                , i                          , &
         &                                                                        retryCount               , retryMaximum               , &
         &                                                                        attemptsRemaining        , rescaleCount               , &
         &                                                                        rescaleMaximum           , treeBuilt                  , &
         &                                                                        massCutoffAttemptsMaximum, massCutoffAttemptsRemaining
    double precision                                                           :: tolerance                , treeBestWorstFit           , &
         &                                                                        multiplier               , constant                   , &
         &                                                                        scalingFactor            , massCutoffScale
    logical                                                                    :: treeBestOverride         , newNodeAboveCutoff         , &
         &                                                                        treeBestNodeAboveCutoff
    
    ! Iterate over all linked trees in this forest.
    call Galacticus_Display_Indent('Augmenting merger tree',verbosityWorking)
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Allocate array of original anchor nodes from which new high-resolution branches will be built.
       nodeCount=augmentSilentWalk(treeCurrent,walkTaskNodeCount)
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
          node                       => anchorNodes(i)%node
          basic                      => node          %basic()
          ! Initialize the current best-known tree to null.
          treeBestWorstFit           =     3.000d0
          treeBest%baseNode          => null()
          treeBestOverride           = .false.
          treeBestNodeAboveCutoff    = .false.
          ! Reset all factors used in tree acceptance and scaling.
          multiplier                 =     0.000d0
          constant                   =     0.000d0
          scalingFactor              =     1.000d0
          ! Currently not used: call augmentScaleChildren(self, node, multiplier, constant, scalingFactor)
          tolerance                  =     0.925d0
          retryMaximum               =    50
          rescaleMaximum             =    20
          rescaleCount               =     0
          retryCount                 =     1
          treeBuilt                  =     0
          attemptsRemaining          = 10000       ! Set remaining number of attempts to maximum allowed.
          massCutoffScale            =     1.0d0
          massCutoffAttemptsMaximum  =    50
          massCutoffAttemptsRemaining= massCutoffAttemptsMaximum
          ! Begin building trees from this node, searching for an acceptable tree.
          if (Galacticus_Verbosity_Level() >= verbosityWorking) then
             message="Building tree from node: "
             message=message//node%index()
             call Galacticus_Display_Indent(message)
          end if
          do while (                                     &
               &     treeBuilt         /= 1              & ! Exit if tree successfully built.
               &    .and.                                &
               &     rescaleCount      <= rescaleMaximum & ! Exit once number of rescalings exceeds maximum allowed.
               &    .and.                                &
               &     attemptsRemaining >  0              & ! Exit if no more attempts remain.
               &   )
             newNodeAboveCutoff=.false.
             treeBuilt         =self%buildTreeFromNode(                         &
                  &                                    node                   , &
                  &                                    .false.                , &
                  &                                    tolerance              , &
                  &                                    self%timeEarliest      , &
                  &                                    treeBest               , &
                  &                                    treeBestWorstFit       , &
                  &                                    treeBestOverride       , &
                  &                                    multiplier             , &
                  &                                    constant               , &
                  &                                    scalingFactor          , &
                  &                                    massCutoffScale        , &
                  &                                    newNodeAboveCutoff     , &
                  &                                    treeBestNodeAboveCutoff  &
                  &                                   )
             ! Check for exhaustion of retry attempts.
             if (retryCount == retryMaximum) then
                ! Rescale the tolerance to allow less accurate tree matches to be accepted in future.
                retryCount  =             0
                tolerance   =tolerance   *0.925d0
                ! Check for exhaustion of rescaling attempts.
                rescaleCount=rescaleCount+1
                if (rescaleCount > rescaleMaximum) call Galacticus_Display_Message('Node build attempts exhausted',verbosityWorking)
             end if
             ! Increment the retry count in cases where the tree was not acepted due to matching tolerance.
             select case (treeBuilt)
             case (treeBuildSuccess         )
                retryCount=retryCount-1
             case (treeBuildFailureTolerance)
                retryCount=retryCount+1
             end select
             ! Decrement the number of attempts remaining and, if the best tree is to be forcibly used, set no attempts remaining.
             attemptsRemaining                      =attemptsRemaining-1
             if (treeBestOverride) attemptsRemaining=                 -1
             ! If all attempts have been used but no match has been foound, insert the best tree on next pass through loop.
             if (attemptsRemaining == 1 .and. associated(treeBest%baseNode))  treeBestOverride=.false. !.true.
             ! If the new tree contains nodes above the mass cut-off, decrement the number of remaining retries.
             if (newNodeAboveCutoff) then
                massCutoffAttemptsRemaining=massCutoffAttemptsRemaining-1
                ! If number of retries is exhausted, adjust the tolerance for declaring nodes to be above the mass cut-off.
                if (massCutoffAttemptsRemaining == 0) then
                   massCutoffAttemptsRemaining=massCutoffAttemptsMaximum
                   massCutoffScale            =massCutoffScale          +0.05d0
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
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
       deallocate(anchorNodes)
    end do
    call Galacticus_Display_Unindent('done',verbosityWorking)
    return
  end subroutine augmentOperate
  
  integer function augmentNoisyWalk(tree,desiredOutput)
    !Walks through tree and prints various types of output.  Will return information about tree based on desiredOutput input string.
    use Galacticus_Nodes
    use Merger_Trees_Pruning_Utilities
    use Input_Parameters
    implicit none
    type   (mergerTree        ), intent(in   ), target :: tree
    type   (treeNode          ), pointer               :: nextNode              , nodePrevious     , &
         &                                                node              , mergeeNode       , &
         &                                                newNode
    class  (nodeComponentBasic), pointer               :: basic    , newBasicComponent, &
         &                                                previousBasicComponent, childBasic, siblingBasic, parentBasic, currentBasic
    type   (mergerTree        ), pointer               :: treeCurrent
    integer :: nodeCount, endNodeCount
    character (len =*) :: desiredOutput
    double precision :: endMass
    nodeCount = 0
    endNodeCount = 0
    endMass = 0
    treeCurrent => tree
    node           => treeCurrent%baseNode
    basic => node%basic()
    do while (associated(node))
      basic => node%basic()
      nodeCount = nodeCount + 1
      if(.not.associated(node%firstChild))  then
        endNodeCount = endNodeCount + 1
        endMass = endMass + basic%mass()
      end if
      write (*,*) 'Node ID --', node%uniqueID()
      !write (*,*) 'Visiting Node with Mass =', basic%mass(), 'At time =', basic%time()
      if(associated(node%firstChild)) then
        childBasic => node%firstChild%basic()
        !write (*,*) 'First Child of Mass = ', childBasic%mass(), ' At time = ', childBasic%time()
      end if
      !if(associated(node%sibling)) then
      !siblingBasic => node%sibling%basic()
      !write (*,*) 'Sibling of Mass = ', siblingBasic%mass(), ' At time = ', siblingBasic%mass()
      !endif
      node => node%walkTree()
    end do
    basic => tree%baseNode%basic()
    if(endMass > basic%mass()) then
      write(*,*) 'End Mass larger than parent mass'
      !write (*,*) 'End Mass -- ', endMass, ' OriginalMass -- ',basic%mass() 
    end if
    if(desiredOutput == 'nodeCount') then
      augmentNoisyWalk = nodeCount
    else if(desiredOutput =='endNodeCount') then
      augmentNoisyWalk = endNodeCount
    else
      augmentNoisyWalk = 0
    end if
  end function augmentNoisyWalk

  integer function augmentSilentWalk(tree,desiredOutput)
    !% Walks through tree and quietly collects information specified by {\normalfont \ttfamily desiredOutput} input enumeration and
    !% returns that information.
    use Galacticus_Nodes
    use Galacticus_Error
    use Merger_Trees_Pruning_Utilities
    implicit none
    type   (mergerTree        ), intent(in   ), target :: tree
    integer                    , intent(in   )         :: desiredOutput
    type   (treeNode          ), pointer               :: nextNode              , nodePrevious     , &
         &                                                node              , mergeeNode       , &
         &                                                newNode
    class  (nodeComponentBasic), pointer               :: basic    , newBasicComponent, &
         &                                                previousBasicComponent
    type   (mergerTree        ), pointer               :: treeCurrent
    integer :: nodeCount, endNodeCount
    double precision :: endNodeMass

    nodeCount = 0
    endNodeCount = 0
    endNodeMass = 0
    treeCurrent => tree
    node           => treeCurrent%baseNode
    basic => node%basic()
    do while (associated(node))
      basic => node%basic()
      nodeCount = nodeCount + 1
      if(.not.associated(node%firstChild))  then
        endNodeCount = endNodeCount + 1
        endNodeMass = endNodeMass + basic%mass()
      end if
      nodePrevious => node%parent
      node => node%walkTree()
    end do
    ! Return the requested quantity.
    select case (desiredOutput)
    case (walkTaskNodeCount)
       augmentSilentWalk=nodeCount
    case (walkTaskEndNodeCount)
       augmentSilentWalk=endNodeCount
    case default
       call Galacticus_Error_Report('augmentSilentWalk','unknown task requested')
    end select
    return
  end function augmentSilentWalk

  subroutine augmentResetUniqueIDs(tree)
    !Walks through tree and returns all negative IDs back to being positive.
    use Galacticus_Nodes
    use Merger_Trees_Pruning_Utilities
    use Input_Parameters
    implicit none
    type   (mergerTree        ), intent(in   ), target :: tree
    type   (treeNode          ), pointer               :: nextNode              , nodePrevious     , &
         &                                                node              , mergeeNode       , &
         &                                                newNode
    class  (nodeComponentBasic), pointer               :: basic    , newBasicComponent, &
         &                                                previousBasicComponent, childBasic, siblingBasic, parentBasic, currentBasic
    type   (mergerTree        ), pointer               :: treeCurrent
    treeCurrent => tree
    node           => treeCurrent%baseNode
    basic => node%basic()
    do while (associated(node)) 
      if (node%uniqueID() < 0) then
        call node%uniqueIDSet(-node%uniqueID())
      end if
      node => node%walkTree()
    end do

  end subroutine augmentResetUniqueIDs


  integer function augmentBuildTreeFromNode(self, node, extendingEndNode, tolerance,timeEarliestIn, treeBest, treeBestWorstFit, treeBestOverride, multiplier, constant, scalingFactor, massCutoffScale, newNodeAboveCutoff, treeBestNodeAboveCutoff)

    use Galacticus_Nodes
    use Merger_Trees_Builders
    use Cosmology_Functions

    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self
    type (treeNode), target:: node
    type (treeNode), pointer :: baseNode, nodePointer
    class (nodeComponentBasic), pointer :: basic
    type (mergerTree) :: newTree
    double precision, intent(in   ) :: timeEarliestIn
    type (mergerTree), intent (inout) :: treeBest
    class (nodeComponentBasic), pointer :: baseBasic, childBasic
    double precision ::  massCutoff, tolerance , timeEarliest
    double precision, intent (inout) :: treeBestWorstFit, multiplier, constant, scalingFactor
    logical :: extendingEndNode
    integer :: noisyCount, endNodeCount, nodeChildCount, i, timeIndex, treeAccepted
    type(mergerTreeOperatorPruneByTime) :: pruneByTime
    class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_
    double precision, intent (inout) :: massCutoffScale
    logical, intent (inout) :: newNodeAboveCutoff
    logical, intent (inout) :: treeBestNodeAboveCutoff
    logical :: treeBestOverride

    basic =>node%basic()
    !Find the next time for timeEarliest.
    if(extendingEndNode) then
      timeEarliest=timeEarliestIn
      !end nodes will be extended to the default min time.
    else 
      if(associated(node%firstChild)) then
        childBasic => node%firstChild%basic()
        timeEarliest = childBasic%time()
        !If there are children, the next time is set to the time of the children.
      else
        if (size(self%timeSnapshots) > 1) then
          i = 1
          timeIndex = -1
          !find the next smallest time in the time snapshots and use that as the timeEarliest.
          do while (i <= size(self%timeSnapshots))
            if((basic%time() <= self%timeSnapshots(i)).or.(basic%time() <= self%timeSnapshots(i)*(1.0 + 1.0d-5)).or.(basic%time() <= self%timeSnapshots(i)*(1.0-1.0d-5 ))) then
              timeIndex = i - 1
              exit
            else 
              i = i + 1
            end if
          end do
          if (timeIndex == -1) then 
            timeIndex = size(self%timeSnapshots)
          end if
          if (timeIndex == 0) then
            timeEarliest = timeEarliestIn
          else 
            timeEarliest = self%timeSnapshots(timeIndex)
          end if
        else 
          !If the size of the snapshot list is reported to be one, revert to child time / default min time determination of earliest time.
          timeEarliest = timeEarliestIn
        end if
      end if
    endif
    massCutoff = self%massResolution
    !Build trees from the nodes above the cutoff. 
    if(basic%mass() > massCutoff) then
      baseNode => treeNode(node%index(), newTree)
      newTree%baseNode => baseNode
      baseBasic => baseNode%basic(autoCreate=.true.)
      call baseBasic%timeSet(basic%time())
      call baseBasic%massSet(basic%mass())
      call self%mergerTreeBuilder_%timeEarliestSet(timeEarliest) 
      call self%mergerTreeBuilder_%build(newTree)
      call pruneByTime%operate(newTree)
      !Set and build tree from the base node of the new tree.
      nodePointer =>node
      endNodeCount = augmentSilentWalk(newTree,walkTaskEndNodeCount)
      call augmentSortChildren(nodePointer)
      nodeChildCount = augmentChildCount(nodePointer)
      !check whether the newly built tree is accepted.
      treeAccepted = augmentAcceptTree(self, nodePointer, newTree, nodeChildCount, extendingEndNode, tolerance, timeEarliest, treeBest, treeBestWorstFit, treeBestOverride, multiplier, constant, scalingFactor, massCutoffScale, newNodeAboveCutoff, treeBestNodeAboveCutoff)
      !check whether currently saved best tree can be accepted now.
      if (((treeBestWorstFit <= tolerance .and. .not. treeBestNodeAboveCutoff) .or. treeBestOverride ).and.(associated(treeBest%baseNode)).and. (.not.(treeAccepted == 1)) ) then
        treeBestWorstFit = 3.0
        if(associated(newTree%baseNode)) then
          call newTree%destroyBranch(newTree%baseNode)
        end if
        newTree%baseNode => treeBest%baseNode
        treeBest%baseNode => null()
        if ( augmentAcceptTree(self, nodePointer, newTree, nodeChildCount, extendingEndnode, tolerance, timeEarliest, treeBest, treeBestWorstFit, treeBestOverride, multiplier, constant, scalingFactor, massCutoffScale, newNodeAboveCutoff, treeBestNodeAboveCutoff) == 1) then
          augmentBuildTreeFromNode = treeBuildSuccess
        else 
          augmentBuildTreeFromNode = treeBuildFailureStructure
        end if


      else if (treeAccepted  == 1 ) then
        augmentBuildTreeFromNode = treeBuildSuccess
        if (associated(treeBest%baseNode)) then
          call treeBest%destroyBranch(treeBest%baseNode)
        end if
        treeBestWorstFit = 3.0
      else
        if (associated(newTree%baseNode)) then
         call newTree%destroyBranch(newTree%baseNode)
        end if
        augmentBuildTreeFromNode = treeAccepted
      end if
    else
      augmentBuildTreeFromNode = treeBuildSuccess

    end if
    !return 1 on successful tree building attempt. 0 on unsuccessful attempt.
  end function augmentBuildTreeFromNode

  subroutine augmentScaleChildren(self, node, multiplier, constant, scalingFactor)
    use Galacticus_Nodes
    implicit none
    class (mergerTreeOperatorAugment), intent(inout) :: self
    type (treeNode), pointer :: node, currentChild
    class (nodeComponentBasic), pointer :: currentComponentBasic, childComponentBasic
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
    currentComponentBasic => node%basic()
    parentMass = currentComponentBasic%mass()
    parentTime = currentComponentBasic%time()
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
          currentComponentBasic => currentChild%basic()
          call currentComponentBasic%massSet(scalingFactor*currentComponentBasic%mass())
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

  subroutine augmentUnscaleChildren (self, node, nodeChildCount, endNodes, multiplier, constant, scalingFactor)
    use Galacticus_Nodes
    implicit none
    class (mergerTreeOperatorAugment), intent(inout) :: self
    type (treeNode), pointer :: node, currentChild, currentNode
    integer, intent (inout) :: nodeChildCount
    class (nodeComponentBasic), pointer :: currentComponentBasic, childComponentBasic
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
      currentComponentBasic => node%basic()
      parentMass = currentComponentBasic%mass()
      currentChild => node%firstChild
      do while (associated(currentChild))
        currentNode => endNodes(i)%node
        do while (associated(currentNode))
          currentComponentBasic = currentNode%basic()
          if (currentNode%uniqueID() > 0) then
            call currentComponentBasic%massSet((currentComponentBasic%mass() / parentMass)*(multiplier*LOG10(currentComponentBasic%time()) + constant))
            call currentNode%uniqueIDSet(-currentNode%uniqueID())
          end if 
          currentNode => currentNode%parent
        end do
        i = i + 1
        currentChild => currentChild%sibling
      end do
    end if  

  end subroutine augmentUnscaleChildren

  integer function augmentAcceptTree(self, node, tree, nodeChildCount, extendingEndNode, tolerance, timeEarliest, treeBest, treeBestWorstFit, treeBestOverride, multiplier, constant, scalingFactor, massCutoffScale, newNodeAboveCutoff, treeBestNodeAboveCutoff)
    use Galacticus_Nodes
    use Merger_Trees_Builders
    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self

    type (treeNode), pointer :: node, currentNode, nodePrevious, currentNodeSibling, nonOverlapNode, firstNonOverlap
    class (nodeComponentBasic), pointer :: currentBasicComponent, sortNodeComponentBasic, nonOverlapComponentBasic
    type (mergerTree), target :: tree
    type (mergerTree), intent (inout), target :: treeBest
    integer :: nodeChildCount, i, j, endNodesSorted, nonChildNodes, retryCount
    type (treeNodeList), dimension(nodeChildCount) :: endNodes
    logical :: treeAccepted, extendingEndNode, endNodeBuildAccept, nodeMassesAgree, currentNodeBelowAll, treeScalable
    double precision :: resolutionLimit, tolerance, massCutoff, timeEarliest, treeCurrentWorstFit, falseWorstFit, unresolvedMass, treeMass, endNodeMass
    double precision, intent(inout) :: treeBestWorstFit, multiplier, constant, scalingFactor, massCutoffScale
    logical, intent (inout) :: newNodeAboveCutoff, treeBestNodeAboveCutoff
    logical :: treeBestOverride
    falseWorstFit = 3.0
    newNodeAboveCutoff = .false.
    massCutoff = self%massResolution
    resolutionLimit = 2.3589041e10
    firstNonOverlap => null()
    nonOverlapNode => null()
    i = 1
    endNodesSorted = 1
    nonChildNodes = 0
    currentNode => tree%baseNode
    currentBasicComponent => currentNode%basic()
    unresolvedMass = currentBasicComponent%mass()
    treeMass = currentBasicComponent%mass()
    !walk through tree and find the end nodes of the tree.
    do while (associated(currentNode))
         currentBasicComponent => currentNode%basic()
         currentNode%event => null()
         currentNode%hostTree => node%hostTree
         if(.not.associated(currentNode%firstChild))  then
           unresolvedMass = unresolvedMass - currentBasicComponent%mass()
           !place end node into endNodes array, corresponding to the largest overlapping nodes.
           if(nodeChildCount > 0) then

             if(endNodesSorted > nodeChildCount) then
               nonChildNodes = nonChildNodes + 1
             end if

             i = 1
             currentNodeBelowAll = .true.
             do while(i <= nodeChildCount)
               if (endNodesSorted > nodeChildCount) then
                 sortNodeComponentBasic => endNodes(i)%node%basic()
                 if (sortNodeComponentBasic%mass() < currentBasicComponent%mass()) then
                   currentNodeBelowAll = .false.
                   j = nodeChildCount
                   nonOverlapNode => endNodes(j)%node
                   sortNodeComponentBasic => endNodes(j)%node%basic()
                   if (sortNodeComponentBasic%mass() > resolutionLimit*massCutoffScale) then 
                     write (*, *) 'NonOverlap Failure at mass --', sortNodeComponentBasic%mass()
                     newNodeAboveCutoff = .true.
                   end if

                   do while (j > i)
                     endNodes(j)%node => endNodes(j-1)%node
                     j = j - 1
                   end do
                   endNodes(i)%node => currentNode
                   i = nodeChildCount
                 end if
               else
                 currentNodeBelowAll = .false.
                 i = endNodesSorted
                 endNodes(i)%node => currentNode
                 i = nodeChildCount
               end if
               i = i + 1
             end do
             if (currentNodeBelowAll) then
               nonOverlapNode => currentNode
               if (currentBasicComponent%mass() > resolutionLimit*massCutoffScale) then
                 write (*,*) 'NonOverlap Failure at mass ==', currentBasicComponent%mass()
                 newNodeAboveCutoff = .true.
               end if
             end if

           else 
             nonOverlapNode => currentNode
             nonChildNodes = nonChildNodes + 1
           end if
           endNodesSorted = endNodesSorted + 1 
         end if
         nodePrevious => currentNode%parent
         currentNode => currentNode%walkTree()
         if (.not.extendingEndNode) then
           call augmentNonOverlapListAdd(nonOverlapNode, firstNonOverlap)
         end if
         nonOverlapNode => null()
    end do
    endNodesSorted = endNodesSorted - 1
    nodeMassesAgree = .true.
    treeCurrentWorstFit = 0.0
    if((nodeChildCount > 0).and.(nodeChildCount <= endNodesSorted)) then
      i = 1
      currentNode => node%firstChild
      do while(associated(currentNode))
        if(.not.augmentNodeComparison(currentNode, endNodes(i)%node, 2.0d0 - 2.0 * tolerance, treeCurrentWorstFit)) then
          nodeMassesAgree = .false.
        end if
        i = i + 1
        currentNode => currentNode%sibling
      end do
    end if

    !if (newNodeAboveCutoff) then
    !  newNodeAboveCutoff = augmentScaleNodesAboveCutoff(self, node, tree, endNodes, nodeChildCount, firstNonOverlap, resolutionLimit)
    !end if 

    if ((nodeChildCount <= endNodesSorted).and.(newNodeAboveCutoff .eqv. .false. .or. treeBestOverride).and.nodeMassesAgree) then
      call augmentExtendNonOverlapNodes(self, tree, firstNonOverlap,massCutoff, tolerance, timeEarliest, treeBest, treeBestWorstFit, massCutoffScale)
      endNodeMass = 0
      i = 1
      do while (i <= nodeChildCount)
        currentBasicComponent => endNodes(i)%node%basic()
        endNodeMass = endNodeMass + currentBasicComponent%mass()
        i = i + 1
      end do
      treeScalable = augmentMultiScale(self, node, tree, endNodes, nodeChildCount, firstNonOverlap, tolerance, timeEarliest, treeBest, treeBestWorstFit, unresolvedMass, treeMass, massCutoffScale)
    else 
      treeScalable = .false.
    end if
    treeAccepted = (nodeChildCount <= endNodesSorted).and.(newNodeAboveCutoff .eqv. .false.).and.treeScalable.and.nodeMassesAgree
    if (treeBestOverride .and. associated(tree%baseNode)) then
      treeAccepted = nodeMassesAgree.and.(nodeChildCount<= endNodesSorted)
    end if

    !call augmentResetUniqueIDs(tree)
    !call augmentUnscaleChildren(self, node, nodeChildCount, endNodes, multiplier, constant, scalingFactor)
    if(treeAccepted) then
      currentBasicComponent => tree%baseNode%basic()
      !write (*,*) 'Building Node ', node%uniqueID(),' at time ', currentBasicComponent%time(), 'to time ', timeEarliest
      call augmentResetUniqueIDs(tree)
      call augmentUnscaleChildren(self, node, nodeChildCount, endNodes, multiplier, constant, scalingFactor)  
      call augmentSimpleInsert(self, node, tree, endNodes, nodeChildCount, firstNonOverlap)
      !call augmentSimpleScale(self, node, tree, endNodes, nodeChildCount, firstNonOverlap)


    else if ((nodeChildCount <= endNodesSorted) .and. .not.treeBestOverride) then
      !write (*,*) 'Updating Best Tree'
      call augmentNonOverlapReinsert(firstNonOverlap)
      if (treeCurrentWorstFit < treeBestWorstFit) then
        if(associated(treeBest%baseNode)) then
          call treeBest%destroyBranch(treeBest%baseNode)
        end if
        treeBest%baseNode => tree%baseNode
        tree%baseNode => null()
        treeBest%baseNode%hostTree => treeBest
        currentNode => treeBest%baseNode
        do while (associated(currentNode))
          currentNode%event => null()
          currentNode%hostTree => treeBest%baseNode%hostTree
          nodePrevious => currentNode%parent
          currentNode => currentNode%walkTree()
        end do
        treeBestWorstFit = treeCurrentWorstFit
        treeBestNodeAboveCutoff = newNodeAboveCutoff
        !write (*,*) 'This Tree is Best Fit With Worst-- ', treeCurrentWorstFit
      end if
    else 
      call augmentNonOverlapReinsert(firstNonOverlap)
      if(associated(tree%baseNode)) then
        call tree%destroyBranch(tree%baseNode)
      end if
    end if
    if (treeAccepted) then    
      augmentAcceptTree = treeBuildSuccess
      if (treeCurrentWorstFit > 0.0 .and. nodeChildCount > 1) then
        open (unit = 55, file = 'toleranceHistogram.py', position = 'append', action = 'write')
        write (55,*) treeCurrentWorstFit, ','
        close (55)
      end if
    else if (.not.nodeMassesAgree) then
      augmentAcceptTree = treeBuildFailureTolerance
    else
      augmentAcceptTree = treeBuildFailureStructure
    end if 

  end function augmentAcceptTree

  subroutine augmentExtendByOverlap(bottomNode, topNode, keepTop, exchangeProperties)
    !This will conjoin two trees by overlapping the topNode of one tree with
    !a chosen bottom node of the other.  If keepTop is true, topNode replaces
    !bottomNode, otherwise, bottomNode replaces topNode.  If exchangeProperties
    !is true, the mass and time information of the deleted node overwrites the 
    !mass and time of the retained node.

    use Galacticus_Nodes

    implicit none
    type (treeNode), pointer :: bottomNode, topNode, currentChild, currentSibling
    logical ::keepTop, exchangeProperties
    class (nodeComponentBasic), pointer :: bottomComponentBasic, topComponentBasic

    bottomComponentBasic => bottomNode%basic()
    topComponentBasic => topNode%basic()

    if (keepTop) then 

      topNode%parent => bottomNode%parent
      topNode%sibling => bottomNode%sibling
      if(associated(bottomNode%parent)) then
        currentSibling => bottomNode%parent%firstChild
      else 
        currentSibling => null()
      end if
      do while(associated(currentSibling))
        if(associated(currentSibling%sibling, bottomNode)) then
          currentSibling%sibling => topNode
        end if 
        currentSibling => currentSibling%sibling
      end do
      if(associated(bottomNode%parent)) then
        if(associated(bottomNode%parent%firstChild, bottomNode)) then
          bottomNode%parent%firstChild => topNode
        end if 
      end if

      if(exchangeProperties) then
        call topComponentBasic%timeSet(bottomComponentBasic%time())
        call topComponentBasic%massSet(bottomComponentBasic%mass())
      end if 

      call bottomNode%destroy()

    else 
      currentChild => topNode%firstChild
      bottomNode%firstChild => topNode%firstChild
      do while(associated(currentChild))
        currentChild%parent => bottomNode
        currentChild => currentChild%sibling
      end do

      if(exchangeProperties) then
        call bottomComponentBasic%timeSet(topComponentBasic%time())
        call bottomComponentBasic%massSet(topComponentBasic%mass())
      end if
      
      call topNode%destroy()
    endif

    return
  end subroutine augmentExtendByOverlap


 integer function augmentChildCount(node)
    type (treeNode), pointer :: node, currentChild
    integer childCount

    currentChild => node%firstChild
    childCount = 0

    do while (associated(currentChild))
      if(associated(currentChild)) then
        childCount = childCount + 1
      end if
      currentChild => currentChild%sibling
    end do
    augmentChildCount = childCount

  end function augmentChildCount

 subroutine augmentSortChildren(node)
    type (treeNode), pointer :: node, currentNode, sortNode, nextNode
    class (nodeComponentBasic), pointer :: currentComponentBasic, sortComponentBasic
    double precision :: largestMass

    if(associated(node%firstChild)) then    
      currentNode => node%firstChild
      currentComponentBasic => currentNode%basic()
      largestMass = currentComponentBasic%mass()
      nextNode => currentNode%sibling
      currentNode%sibling => null()
      do while (associated(nextNode))
        sortNode => node%firstChild
        currentNode => nextNode
        nextNode => nextNode%sibling
        currentComponentBasic => currentNode%basic()
        if (currentComponentBasic%mass() > largestmass) then
          node%firstChild => currentNode
          largestmass = currentComponentBasic%mass()
          currentNode%sibling => sortNode
        else 
          do while (associated(sortNode%sibling))
            sortComponentBasic => sortNode%sibling%basic()
            if (currentComponentBasic%mass() > sortComponentBasic%mass()) then
              exit
            end if
            sortNode => sortNode%sibling
          end do
  
          currentNode%sibling => sortNode%sibling
          sortNode%sibling => currentNode
         end if
      end do
    end if
  end subroutine augmentSortChildren

  subroutine augmentNonOverlapListAdd(node, listFirstElement)

    implicit none
    type (treeNode), pointer :: node, listFirstElement, currentSibling
    if(associated(node)) then
      if(associated(node%parent)) then
        currentSibling => node%parent%firstChild
        if(.not.associated(listFirstElement)) then 
          listFirstElement => node
          if(associated(currentSibling, node).and.(associated(node%sibling))) then
            node%parent%firstChild => node%sibling
          else
            do while (associated(currentSibling))
              if(associated(currentSibling%sibling, node)) then
                currentSibling%sibling => node%sibling
              end if
              currentSibling => currentSibling%sibling
            end do
          end if
          listFirstElement%sibling => null()
        else
          if(associated(currentSibling, node).and.(associated(node%sibling))) then
            node%parent%firstChild => node%sibling
          else
            do while (associated(currentSibling))
              if(associated(currentSibling%sibling, node)) then
                currentSibling%sibling => node%sibling
              end if
              currentSibling => currentSibling%sibling
            end do 
          end if
          node%sibling => listFirstElement
          listFirstElement => node
        end if
      else 
        if(.not.associated(listFirstElement)) then
          listFirstElement => node
          listFirstElement%sibling => null()
        else 
          node%sibling => listFirstElement
          listFirstElement => node
        end if
      end if
    end if
    
  end subroutine augmentNonOverlapListAdd

  subroutine augmentNonOverlapReinsert(listFirstElement)
    implicit none 
    type (treeNode), pointer :: currentNode, listFirstElement
    do while (associated(listFirstElement))
      currentNode => listFirstElement
      listFirstElement => listFirstElement%sibling
      if (associated(currentNode%parent)) then
        if(.not.associated(currentNode, currentNode%parent%firstChild)) then
          currentNode%sibling => currentNode%parent%firstChild
        else
          currentNode%sibling => null()
        end if
        currentNode%parent%firstChild => currentNode
        call augmentSortChildren(currentNode%parent)
      end if
    end do 

  end subroutine augmentNonOverlapReinsert

  subroutine augmentExtendNonOverlapNodes(self, tree, firstNonOverlap,massCutoff, tolerance, timeEarliest, treeBest, treeBestWorstFit, massCutoffScale)
    use Merger_Trees_Builders
    use Galacticus_Nodes
    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self
    type (treeNode), pointer :: currentNode, nonOverlapNode, firstNonOverlap
    class (nodeComponentBasic), pointer :: nonOverlapComponentBasic
    type (mergerTree), target :: tree
    type (mergerTree), intent (inout), target :: treeBest
    integer :: retryCount
    double precision, intent(inout) :: treeBestWorstFit, tolerance, timeEarliest, massCutoff, massCutoffScale 
    double precision ::falseWorstFit, falseMultiplier, falseConstant, falseScalingFactor
    logical :: falseNewNodeAboveCutoff, falseBestTreeNodeAboveCutoff
    falseBestTreeNodeAboveCutoff = .false.
    falseNewNodeAboveCutoff = .false.
    falseMultiplier = 0.0
    falseConstant = 0.0
    falseScalingFactor = 1.0
    falseWorstFit = 3.0
    currentNode => firstNonOverlap
    do while (associated(currentNode))
      nonOverlapNode => currentNode
      nonOverlapComponentBasic => nonOverlapNode%basic()
      currentNode => currentNode%sibling
      if(nonOverlapComponentBasic%mass() > massCutoff) then
        retryCount = 100
        do while ((.not.(self%buildTreeFromNode(nonOverlapNode, .true., tolerance, timeEarliest, treeBest, falseWorstFit, .false., falseMultiplier, falseConstant, falseScalingFactor, massCutoffScale, falseNewNodeAboveCutoff, falseBestTreeNodeAboveCutoff)==1)).and.(retryCount >0))
          retryCount = retryCount - 1
        end do
      end if
    end do
  end subroutine augmentExtendNonOverlapNodes

  logical function augmentNodeComparison(newNode, oldNode, tolerance, treeCurrentWorstFit)

    use Numerical_Comparison

    implicit none   

    type (treeNode), pointer :: newNode, oldNode
    class (nodeComponentBasic), pointer :: newComponentBasic, oldComponentBasic
    double precision :: tolerance, thisFit
    double precision, intent(inout) :: treeCurrentWorstFit

    newComponentBasic => newNode%basic()
    oldComponentBasic => oldNode%basic()
    augmentNodeComparison = Values_Agree(newComponentBasic%mass(), oldComponentBasic%mass(), relTol =tolerance)

    thisFit = 2.0*ABS(newComponentBasic%mass() - oldComponentBasic%mass())/ABS(newComponentBasic%mass() + oldComponentBasic%mass())
    if ( thisFit > treeCurrentWorstFit) then
      treeCurrentWorstFit = thisFit
    end if
  end function augmentNodeComparison

  subroutine augmentSimpleInsert(self, node, tree, endNodes, nodeChildCount, firstNonOverlap)
    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self
    type (treeNode), pointer :: node, currentNode, nodePrevious, firstNonOverlap, currentNodeSibling
    type (mergerTree), target :: tree
    integer, intent(inout) :: nodeChildCount
    integer :: i
    type (treeNodeList), dimension(nodeChildCount), intent(inout) :: endNodes
    logical :: scalableTest

    i = 0
    call augmentNonOverlapReinsert(firstNonOverlap)

    if (nodeChildCount > 0) then
      i =1
      currentNode => node%firstChild

      do while (associated(currentNode))
        currentNodeSibling => currentNode%sibling
        node%firstChild => currentNode%sibling
        call augmentExtendByOverLap(endNodes(i)%node, currentNode,.true., .false.)
        currentNode => currentNodeSibling
        i = i + 1
      end do
    else
      currentNode => tree%baseNode
      do while (associated(currentNode))
           currentNode%event => null()
           currentNode%hosttree => node%hostTree
           nodePrevious => currentNode%parent
           currentNode => currentNode%walkTree()
      end do

    end if
   
    call augmentExtendByOverLap(node, tree%baseNode, .false., .false.)

  end subroutine augmentSimpleInsert

  integer function augmentNegativeChildCount(node)
    use Galacticus_Nodes
    implicit none
    type (treeNode), pointer :: node, currentChild
    integer :: negativeCount
  
    negativeCount = 0
    currentChild => node%firstChild
    do while (associated(currentChild))
      if (currentChild%uniqueID() < 0) then
        negativeCount = negativeCount + 1
      end if
      currentChild => currentChild%sibling
    end do
    augmentNegativeChildCount = negativeCount 

   end function augmentNegativeChildCount

  subroutine augmentSimpleScale(self, node, tree, endNodes, nodeChildCount, firstNonOverlap)
    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self
    type (treeNode), pointer :: node, currentNode, nodePrevious, firstNonOverlap, currentNodeSibling, sortNode
    type (mergerTree), target :: tree
    integer, intent(inout) :: nodeChildCount
    integer :: i
    type (treeNodeList), dimension(nodeChildCount), intent(inout) :: endNodes
    logical :: scalableTest
    double precision :: massDifference, nonOverlapMass, scaleFactor
    class (nodeComponentBasic), pointer :: currentComponentBasic, childComponentBasic

    write (*,*) "Starting Simple Scale"
    massDifference = 0
    i = 1
    if (nodeChildCount > 0) then
      currentnode => node%firstChild
      do while (associated(currentNode))
        currentComponentBasic => currentNode%basic()
        childComponentBasic => endNodes(i)%node%basic()
        massDifference = massDifference + currentComponentBasic%mass() - childComponentBasic%mass()
        sortNode => endNodes(i)%node
        do while (associated(sortNode))
          call sortNode%uniqueIDSet(-sortNode%uniqueID())
          sortNode => sortNode%parent
        end do 
        i = i + 1
        currentNode => currentNode%sibling
      end do
      currentNode => firstNonOverlap
      do while(associated(currentNode))
        currentComponentBasic => currentNode%basic()
        nonOverlapMass = nonOverlapMass + currentComponentBasic%mass()
        currentNode => currentNode%sibling
      end do
      if (.not. (nonOverlapMass ==0)) then
        scaleFactor = (nonOverlapMass - massDifference)/nonOverlapMass
      else 
        scaleFactor = 1
      end if
      currentnode => firstNonOverlap
      do while(associated(currentNode))
        sortNode => currentNode
        do while(associated(sortNode))
          currentComponentBasic => sortNode%basic()
          if(sortNode%uniqueID() > 0) then
            call currentComponentBasic%massSet(scaleFactor * currentComponentBasic%mass())
          end if
          sortNode => sortNode%parent
        end do
        currentNode => currentNode%sibling
      end do
    end if
    call augmentNonOverlapReinsert(firstNonOverlap)

    if (nodeChildCount > 0) then
      i =1
      currentNode => node%firstChild

      do while (associated(currentNode))
        currentNodeSibling => currentNode%sibling
        node%firstChild => currentNode%sibling
        call augmentExtendByOverLap(endNodes(i)%node, currentNode,.true., .false.)
        currentNode => currentNodeSibling
        i = i + 1
      end do
    else
      currentNode => tree%baseNode
      do while (associated(currentNode))
           currentNode%event => null()
           currentNode%hosttree => node%hostTree
           nodePrevious => currentNode%parent
           currentNode => currentNode%walkTree()
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
    class (nodeComponentBasic), pointer :: currentComponentBasic, childComponentBasic
    type (treeNode), pointer :: currentNode, currentChild, nodePrevious
    double precision :: childMass, unresolvedMass

    currentNode => tree%baseNode
    currentComponentBasic => currentNode%basic()
    do while (associated (currentNode))
      currentComponentBasic => currentNode%basic()
      childMass = 0
      unresolvedMass = 0
      currentChild => currentNode%firstChild
      do while(associated(currentChild))
        childComponentBasic => currentChild%basic()
        childMass = childMass + childComponentBasic%mass()
        currentChild => currentChild%sibling
      end do
      if (childMass > 0) then
        unresolvedMass = currentComponentBasic%mass() - childMass
        call currentComponentBasic%massSet(currentComponentBasic%mass() - unresolvedMass)
      end if 
      nodePrevious => currentNode%parent
      currentNode => currentNode%walkTree()
    end do

  end subroutine augmentRemoveUnresolvedMass

  subroutine augmentInsertChildMass(self, node, originalChildNode, newChildNode)
    use Galacticus_Nodes
    implicit none
    class (mergerTreeOperatorAugment), intent(inout) :: self
    type (treeNode), pointer :: node, originalChildNode, newChildNode, scaleNode
    class (nodeComponentBasic), pointer :: currentComponentBasic, childComponentBasic
    double precision :: multiplier, constant, scalingFactor
    double precision :: parentMass, childMass, parentTime, childTime, massDifference

    currentComponentBasic => node%basic()
    parentMass = currentComponentBasic%mass()
    parentTime = currentComponentBasic%time()
    childComponentBasic => originalChildNode%basic()
    childTime = childComponentBasic%time()
    massDifference = childComponentBasic%mass()
    childComponentBasic => newChildNode%basic()
    massDifference = massDifference - childComponentBasic%mass()
    multiplier = (massDifference)/(LOG10(childTime) - LOG10(parentTime))
    constant = -multiplier *LOG10(parentTime)
    scaleNode => newChildNode
    do while (associated(scaleNode))
        currentComponentBasic => scaleNode%basic()
        call currentComponentBasic%massSet(currentComponentBasic%mass() + multiplier*LOG10(currentComponentBasic%time()) + constant)
        scaleNode => scaleNode%parent
    end do

  end subroutine augmentInsertChildMass

  logical function augmentScaleNodesAboveCutoff(self, node, tree, endNodes, nodeChildCount, firstNonOverlap, resolutionLimit)
    use Galacticus_Nodes
    implicit none
    class (mergerTreeOperatorAugment), intent (inout) :: self
    type (treeNode), pointer :: node, firstNonOverlap, currentNonOverlap, currentChild, currentEndNode, scaleNode
    type (mergerTree), target :: tree
    integer, intent(inout) :: nodeChildCount
    type (treeNodeList), dimension(nodeChildCount), intent(inout) :: endNodes
    class (nodeComponentBasic), pointer :: currentComponentBasic, childComponentBasic
    double precision, intent (inout) :: resolutionLimit
    double precision :: massAboveCutoff, overlapMassDifference, largestMassAboveCutoff, scalingFactor, totalEndNodeMass, scalingMass
    logical :: nodesScalable
    integer :: i

    nodesScalable = .true.
    overlapMassDifference = 0.0
    totalEndNodeMass = 0.0
    i = 1
    currentChild => node%firstChild
    do while (associated(currentChild))
      currentEndNode => endNodes(i)%node
      currentComponentBasic => currentEndNode%basic()
      childComponentBasic => currentChild%basic()
      overlapMassDifference = overlapMassDifference + currentComponentBasic%mass() - childComponentBasic%mass()
      totalEndNodeMass = totalEndNodeMass + currentComponentBasic%mass()
      !do while (associated(currentEndNode))
      !  call currentEndNode%uniqueIDSet(-currentEndNode%uniqueID())
      !  currentEndNode => currentEndNode%parent
      !end do
      i = i + 1
      currentChild => currentChild%sibling
    end do

    massAboveCutoff = 0.0
    largestMassAboveCutoff = 0.0
    currentNonOverlap => firstNonOverlap
    do while (associated(currentNonOverlap))
      currentComponentBasic => currentNonOverlap%basic()
      if (currentComponentBasic%mass() > resolutionLimit) then
        massAboveCutoff = massAboveCutoff + currentComponentBasic%mass()
      end if 
      if (currentComponentBasic%mass() > largestMassAboveCutoff) then
        largestMassAboveCutoff = currentComponentBasic%mass()
      end if
      currentNonOverlap => currentNonOverlap%sibling
    end do

    if (overlapMassDifference < 0) then
      nodesScalable = .false.
    else if (overlapMassDifference < massAboveCutoff) then
      if ( largestMassAboveCutoff * (( massAboveCutoff - overlapMassDifference) /  massAboveCutoff) > resolutionLimit) then 
        nodesScalable = .false.
      else 
        currentNonOverlap => firstNonOverlap
        do while (associated(currentNonOverlap))
        currentComponentBasic => currentNonOverlap%basic()
        if ( currentComponentBasic%mass() > resolutionLimit) then
          scalingMass = currentComponentBasic%mass() * ( 1 - (massAboveCutoff - overlapMassDifference) / massAboveCutoff)
          scaleNode => currentNonOverlap
          do while (associated(scaleNode))
            currentComponentBasic => scaleNode%basic()
            call currentComponentBasic%massSet(currentComponentBasic%mass() - scalingMass)
            scaleNode => scaleNode%parent
          end do
        end if
        currentNonOverlap => currentNonOverlap%sibling
        end do

        i = 1
        currentChild => node%firstChild
        do while (associated(currentChild))
          currentEndNode => endNodes(i)%node
          currentComponentBasic => currentEndNode%basic()
          childComponentBasic => currentChild%basic()
          scalingMass = childComponentBasic%mass() - currentComponentBasic%mass()
          scaleNode => currentEndNode
          do while (associated(scaleNode))
            currentComponentBasic => scaleNode%basic()
            call currentComponentBasic%massSet(currentComponentBasic%mass() + scalingMass)
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


  logical function augmentMultiScale(self, node, tree, endNodes, nodeChildCount, firstNonOverlap, tolerance, timeEarliest, treeBest, treeBestWorstFit, unresolvedMass,  treeMass, massCutoffScale)

    use Galacticus_Nodes
    use Merger_Trees_Builders
    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self
    type (treeNode), pointer :: node, currentNode, nodePrevious, firstNonOverlap, currentNodeSibling, currentChildNode,sortNode, currentGrandchild, scaleNode
    type (mergerTree), target :: tree
    type (mergerTree), intent (inout), target :: treeBest
    class (nodeComponentBasic), pointer :: childComponentBasic, currentComponentBasic, sortComponentBasic
    integer, intent(inout) :: nodeChildCount
    integer :: i, retryCount
    type (treeNodeList), dimension(nodeChildCount), intent(inout) :: endNodes
    logical :: treeScaled, lastNodeFound
    double precision, intent (inout) :: tolerance, timeEarliest, treeBestWorstFit, unresolvedMass, treeMass, massCutoffScale
    double precision :: massExcess, childNodeMass, currentMass, endNodeMass, falseWorstFit, massDifferenceScaleFactor
  
    falseWorstFit = 3.0
    call augmentNonOverlapReinsert(firstNonOverlap)
    if (nodeChildCount > 0) then
      currentNode => tree%baseNode
      currentComponentBasic => currentNode%basic()
      i = 1
      currentChildNode => node%firstChild
      massExcess = 0
      endNodeMass = 0
      childNodeMass = 0
      treeMass = currentComponentBasic%mass()
      do while (associated(currentChildNode))
        currentNode => endNodes(i)%node
        currentComponentBasic => currentNode%basic()
        childComponentBasic => currentChildNode%basic()
        massExcess = massExcess + currentComponentBasic%mass() - childComponentBasic%mass()
        endNodeMass = endNodeMass + currentComponentBasic%mass()
        childNodeMass = childNodeMass + childComponentBasic%mass()
        call currentNode%uniqueIDSet(-currentNode%uniqueID())
        i = i + 1
        currentChildNode => currentChildNode%sibling
      end do
      treeScaled = .true.
      if ( treeScaled) then
        i = 1
        currentChildNode => node%firstChild
        do while (associated(currentChildNode))
          currentNode => endNodes(i)%node
          currentComponentBasic => currentNode%basic()
          call augmentInsertChildMass(self, tree%baseNode, currentChildNode, currentNode)
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
