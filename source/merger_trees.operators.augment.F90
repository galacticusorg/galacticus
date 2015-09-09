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
     integer :: augmentTimeCount
     double precision, allocatable, dimension (:) :: augmentTimeSnapshots
     double precision :: augmentResolutionLimit, augmentTimeEarliest
     class(mergerTreeBuilderClass), pointer :: mergerTreeBuilder_
   contains
     final     ::            augmentDestructor
     procedure :: operate => augmentOperate
  end type mergerTreeOperatorAugment

  interface mergerTreeOperatorAugment
     !% Constructors for the {\normalfont \ttfamily augment} merger tree operator class.
     module procedure augmentConstructorParameters
     module procedure augmentConstructorInternal
  end interface mergerTreeOperatorAugment

contains

  function augmentConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily augment} merger tree operator class which takes a parameter set as input.
    use Input_Parameters2
    use Cosmology_Functions
    use Memory_Management
    use Galacticus_Error
    implicit none
    type(mergerTreeOperatorAugment)                :: augmentConstructorParameters
    type(inputParameters              ), intent(in   ) :: parameters
    double precision, allocatable, dimension(:) :: augmentTimeSnapshots
    double precision :: augmentResolutionLimit, augmentTimeEarliest
    integer :: augmentTimeCount
    class (cosmologyFunctionsClass), pointer :: cosmologyFunctions_
    class (mergerTreeBuilderClass), pointer :: mergerTreeBuilder_
    integer :: i
     double precision, parameter :: expansionFactorDefault=0.01d0
   !# <inputParameterList label="allowedParameterNames" />

    !# <objectBuilder class="mergerTreeBuilder" name="mergerTreeBuilder_" source="parameters"/>
    !# <inputParameter>
    !#   <name>augmentResolutionLimit</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d10</defaultValue>
    !#   <description>For the {\normalfont \ttfamily augment} operator a description of resolution limit for new trees.</description>
    !#   <type>double precision</type>
    !#   <cardinality>0..</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>augmentTimeCount</name>
    !#   <source>parameters</source>
    !#   <defaultValue>100</defaultValue>
    !#   <description>Number of points in time to use for building merger trees.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>

    if (parameters%isPresent('augmentSnapshotRedshifts')) then
      allocate(augmentTimeSnapshots(parameters%count('augmentSnapshotRedshifts')))
      if (augmentTimeCount > 1) then
        if(size(augmentTimeSnapshots)/= augmentTimeCount) call Galacticus_Error_Report('augmentConstructorParameters', 'mismatch between [augmentTimeCount] and size of [augmentSnapshotRedshifts]')
      end if
    else
       call Galacticus_Error_Report('augmentConstructorParameters','parameter [augmentTimeSnapshots] is required')
    end if
    !# <inputParameter>
    !#   <name>augmentSnapshotRedshifts</name>
    !#   <variable>augmentTimeSnapshots</variable>
    !#   <source>parameters</source>
    !#   <description>For {\normalfont \ttfamily augment} description of redshift snapshots.</description>
    !#   <type>double precision</type>
    !#   <cardinality> 0..</cardinality>
    !# </inputParameter>
    cosmologyFunctions_ => cosmologyFunctions()
    do i =1,size(augmentTimeSnapshots)
      augmentTimeSnapshots(i) = cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(augmentTimeSnapshots(i)))
    end do

        augmentConstructorParameters%augmentTimeEarliest=min(cosmologyFunctions_%cosmicTime(expansionFactorDefault),augmentTimeSnapshots(1))
    
    augmentConstructorParameters = augmentConstructorInternal(augmentTimeCount, augmentTimeSnapshots, augmentResolutionLimit, augmentTimeEarliest,mergerTreeBuilder_)

    return
  end function augmentConstructorParameters


  function augmentConstructorInternal(augmentTimeCount, augmentTimeSnapshots,augmentResolutionLimit, augmentTimeEarliest,mergerTreeBuilder_)
    !% Internal constructor for the {\normalfont \ttfamily augment} merger tree operator class.
    use Memory_Management
    use Cosmology_Functions
    use Sort
    implicit none
    type            (mergerTreeOperatorAugment)                :: augmentConstructorInternal
    integer, intent (in) :: augmentTimeCount
    double precision                               , intent(in   ) :: augmentResolutionLimit, augmentTimeEarliest
    double precision, intent (in), dimension(:) :: augmentTimeSnapshots
    class(mergerTreeBuilderClass), pointer, intent(in   ) :: mergerTreeBuilder_
    class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_
    double precision, parameter :: expansionFactorDefault=0.01d0

    call Alloc_Array(augmentConstructorInternal%augmentTimeSnapshots, [augmentTimeCount])
    augmentConstructorInternal%augmentTimeSnapshots = augmentTimeSnapshots
    call Sort_Do(augmentConstructorInternal%augmentTimeSnapshots)
    augmentConstructorInternal%augmentResolutionLimit = augmentResolutionLimit

    cosmologyFunctions_ => cosmologyFunctions()
    augmentConstructorInternal%augmentTimeEarliest=min(cosmologyFunctions_%cosmicTime(expansionFactorDefault),augmentConstructorInternal%augmentTimeSnapshots(1))
    augmentConstructorInternal%mergerTreeBuilder_ => mergerTreeBuilder_
    call augmentConstructorInternal%mergerTreeBuilder_%timeEarliestSet(augmentConstructorInternal%augmentTimeEarliest)
    
    return
  end function augmentConstructorInternal

  subroutine augmentDestructor(self)
    !% Destructor for the augment merger tree operator function class.
    implicit none
    type(mergerTreeOperatorAugment), intent(inout) :: self

    !# <objectDestructor name="self%mergerTreeBuilder_" />
    return
  end subroutine augmentDestructor



  subroutine augmentOperate(self, tree)
    !% Walk through nodes of {\normalfont \ttfamily thisTree}.
    use Galacticus_Nodes
    use Merger_Trees_Pruning_Utilities
    use Input_Parameters
    use Memory_Management

    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self
    type   (mergerTree        ), intent(inout ), target :: tree
    type (mergerTree) :: bestTree
    type   (treeNode          ), pointer               :: nextNode              , previousNode     , &
         &                                                thisNode              , mergeeNode       , &
         &                                                newNode
    class  (nodeComponentBasic), pointer               :: thisBasicComponent    , newBasicComponent, &
         &                                                previousBasicComponent
    type   (mergerTree        ), pointer               :: currentTree
    type (treeNodeList), allocatable, dimension (:) :: anchorNodes
    integer :: nodeCount, noisyCount, i, currentNodeCount, retryCount, iterationMax,attemptMax,  rescaleCount, rescaleMax, childCount, treeBuilt
    double precision :: tolerance, bestTreeWorstFit, multiplier, constant, scalingFactor
    double precision :: massCutoffScale
    integer :: massCutoffRescale, massCutoffRetry
    logical :: bestTreeOverride, newNodeAboveCutoff, bestTreeNodeAboveCutoff
    currentTree => tree
    do while (associated(currentTree))

      !Allocate array of original anchor nodes.
      nodeCount = augmentSilentWalk(currentTree, 'nodeCount')
      write (*,*) "Number of Nodes in Tree --", nodeCount
      ALLOCATE (anchorNodes(nodeCount))
      thisNode => currentTree%baseNode
      do i=1,nodeCount
        previousNode => thisNode%parent
        anchorNodes(i)%node => thisNode
        thisBasicComponent => anchorNodes(i)%node%basic()
        thisNode => thisNode%walkTree()
      end do
      i = 1
      thisNode           => anchorNodes(i)%node
      thisBasicComponent => thisNode%basic()
      currentNodeCount = 0
      ! Walk the tree
      do while (i <= nodeCount)
          thisNode => anchorNodes(i)%node
          thisBasicComponent => thisNode%basic()
          currentNodeCount = currentNodeCount + 1
          multiplier = 0.0
          constant = 0.0
          scalingFactor = 1.0
          !call augmentScaleChildren(self, thisNode, multiplier, constant, scalingFactor)
          tolerance = 0.925d0
          iterationMax = 50
          rescaleMax =  20
          rescaleCount = 0
          retryCount = 1
          bestTreeWorstFit = 3.0
          bestTree%baseNode => null()
          treeBuilt = 0
          attemptMax = 10000
          massCutoffScale = 1.0
          massCutoffRescale = 50
          massCutoffRetry = massCutoffRescale
          bestTreeOverride = .false.
          bestTreeNodeAboveCutoff = .false.
          write (*,*) 'Building Tree from Node --', i
          do while (.not.(treeBuilt == 1).and.(rescaleCount <= rescaleMax).and.(attemptMax > 0))
            newNodeAboveCutoff = .false.
            treeBuilt = augmentBuildTreeFrom(self, thisNode, .false., tolerance, self%augmentTimeEarliest, bestTree, bestTreeWorstFit, bestTreeOverride, multiplier, constant, scalingFactor, massCutoffScale, newNodeAboveCutoff, bestTreeNodeAboveCutoff)
            if (retryCount == iterationMax) then
              !Rescale
              retryCount = 0
              tolerance = tolerance*0.925d0 
              rescaleCount = rescaleCount + 1
              if (rescaleCount > rescaleMax) then
                write (*,*) 'Node Build Attempts Exhausted'
              end if              
            end if
            retryCount = retryCount - treeBuilt
            !retryCount only decreased by tolerance failures.  Structural failures do not affect tolerance scaling.
            attemptMax = attemptMax - 1
            if (bestTreeOverride) then
              attemptMax = -1
            end if

            if (attemptMax == 1 .and. associated(bestTree%baseNode)) then
              bestTreeOverride = .false. !.true.
              !If out of attempts and no match within tolerance has been found, insert the best tree on next pass through loop.
            end if

            if (newNodeAboveCutoff) then
              massCutoffRetry = massCutoffRetry - 1
              if (massCutoffRetry == 0) then
                massCutoffRetry = massCutoffRescale
                massCutoffScale = massCutoffScale + 0.05 
              end if 
            end if 
          end do
          i = i + 1
          if (associated(bestTree%baseNode)) then
            call bestTree%destroyBranch(bestTree%baseNode)
            bestTree%baseNode => null() 
          end if
          write (*,*) "Tree Built from Node --", i
      end do
      ! Move to the next tree.
      currentTree => currentTree%nextTree
      DEALLOCATE (anchorNodes)
    end do
    return
  end subroutine augmentOperate

  integer function augmentNoisyWalk(tree,desiredOutput)
    !Walks through tree and prints various types of output.  Will return information about tree based on desiredOutput input string.
    use Galacticus_Nodes
    use Merger_Trees_Pruning_Utilities
    use Input_Parameters
    implicit none
    type   (mergerTree        ), intent(in   ), target :: tree
    type   (treeNode          ), pointer               :: nextNode              , previousNode     , &
         &                                                thisNode              , mergeeNode       , &
         &                                                newNode
    class  (nodeComponentBasic), pointer               :: thisBasicComponent    , newBasicComponent, &
         &                                                previousBasicComponent, childBasic, siblingBasic, parentBasic, currentBasic
    type   (mergerTree        ), pointer               :: currentTree
    integer :: nodeCount, endNodeCount
    character (len =*) :: desiredOutput
    double precision :: endMass
    nodeCount = 0
    endNodeCount = 0
    endMass = 0
    currentTree => tree
    thisNode           => currentTree%baseNode
    thisBasicComponent => thisNode%basic()
    do while (associated(thisNode))
      thisBasicComponent => thisNode%basic()
      nodeCount = nodeCount + 1
      if(.not.associated(thisNode%firstChild))  then
        endNodeCount = endNodeCount + 1
        endMass = endMass + thisBasicComponent%mass()
      end if
      write (*,*) 'Node ID --', thisNode%uniqueID()
      !write (*,*) 'Visiting Node with Mass =', thisBasicComponent%mass(), 'At time =', thisBasicComponent%time()
      if(associated(thisNode%firstChild)) then
        childBasic => thisNode%firstChild%basic()
        !write (*,*) 'First Child of Mass = ', childBasic%mass(), ' At time = ', childBasic%time()
      end if
      !if(associated(thisNode%sibling)) then
      !siblingBasic => thisNode%sibling%basic()
      !write (*,*) 'Sibling of Mass = ', siblingBasic%mass(), ' At time = ', siblingBasic%mass()
      !endif
      thisNode => thisNode%walkTree()
    end do
    thisBasicComponent => tree%baseNode%basic()
    if(endMass > thisBasicComponent%mass()) then
      write(*,*) 'End Mass larger than parent mass'
      !write (*,*) 'End Mass -- ', endMass, ' OriginalMass -- ',thisBasicComponent%mass() 
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
    !Walks through tree and quietly collects information specified by desiredOutput input string and returns that information.
    use Galacticus_Nodes
    use Merger_Trees_Pruning_Utilities
    use Input_Parameters
    implicit none
    type   (mergerTree        ), intent(in   ), target :: tree
    type   (treeNode          ), pointer               :: nextNode              , previousNode     , &
         &                                                thisNode              , mergeeNode       , &
         &                                                newNode
    class  (nodeComponentBasic), pointer               :: thisBasicComponent    , newBasicComponent, &
         &                                                previousBasicComponent
    type   (mergerTree        ), pointer               :: currentTree
    integer :: nodeCount, endNodeCount
    double precision :: endNodeMass
    character (len =*) :: desiredOutput
    nodeCount = 0
    endNodeCount = 0
    endNodeMass = 0
    currentTree => tree
    thisNode           => currentTree%baseNode
    thisBasicComponent => thisNode%basic()
    do while (associated(thisNode))
      thisBasicComponent => thisNode%basic()
      nodeCount = nodeCount + 1
      if(.not.associated(thisNode%firstChild))  then
        endNodeCount = endNodeCount + 1
        endNodeMass = endNodeMass + thisBasicComponent%mass()
      end if
      previousNode => thisNode%parent
      thisNode => thisNode%walkTree()
    end do

    if(desiredOutput == 'nodeCount') then
      augmentSilentWalk = nodeCount
    else if(desiredOutput =='endNodeCount') then
      augmentSilentWalk = endNodeCount
    else
      augmentSilentWalk = 0
    end if
  end function augmentSilentWalk

  subroutine augmentResetUniqueIDs(tree)
    !Walks through tree and returns all negative IDs back to being positive.
    use Galacticus_Nodes
    use Merger_Trees_Pruning_Utilities
    use Input_Parameters
    implicit none
    type   (mergerTree        ), intent(in   ), target :: tree
    type   (treeNode          ), pointer               :: nextNode              , previousNode     , &
         &                                                thisNode              , mergeeNode       , &
         &                                                newNode
    class  (nodeComponentBasic), pointer               :: thisBasicComponent    , newBasicComponent, &
         &                                                previousBasicComponent, childBasic, siblingBasic, parentBasic, currentBasic
    type   (mergerTree        ), pointer               :: currentTree
    currentTree => tree
    thisNode           => currentTree%baseNode
    thisBasicComponent => thisNode%basic()
    do while (associated(thisNode)) 
      if (thisNode%uniqueID() < 0) then
        call thisNode%uniqueIDSet(-thisNode%uniqueID())
      end if
      thisNode => thisNode%walkTree()
    end do

  end subroutine augmentResetUniqueIDs


  integer function augmentBuildTreeFrom(self, thisNode, extendingEndNode, tolerance, timeEarliest, bestTree, bestTreeWorstFit, bestTreeOverride, multiplier, constant, scalingFactor, massCutoffScale, newNodeAboveCutoff, bestTreeNodeAboveCutoff)

    use Galacticus_Nodes
    use Merger_Trees_Builders
    use Cosmology_Functions

    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self
    type (treeNode), target:: thisNode
    type (treeNode), pointer :: baseNode, thisNodePointer
    class (nodeComponentBasic), pointer :: thisBasicComponent
    type (mergerTree) :: newTree
    type (mergerTree), intent (inout) :: bestTree
    class (nodeComponentBasic), pointer :: baseBasic, childBasic
    double precision ::  massCutoff, tolerance, timeEarliest 
    double precision, intent (inout) :: bestTreeWorstFit, multiplier, constant, scalingFactor
    logical :: extendingEndNode
    integer :: noisyCount, endNodeCount, nodeChildCount, i, timeIndex, treeAccepted
    type(mergerTreeOperatorPruneByTime) :: pruneByTime
    class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_
    double precision, parameter :: expansionFactorDefault = 0.01d0
    double precision, intent (inout) :: massCutoffScale
    logical, intent (inout) :: newNodeAboveCutoff
    logical, intent (inout) :: bestTreeNodeAboveCutoff
    logical :: bestTreeOverride

    cosmologyFunctions_ => cosmologyFunctions()
    thisBasicComponent =>thisNode%basic()
    !Find the next time for timeEarliest.
    if(extendingEndNode) then
      timeEarliest=cosmologyFunctions_%cosmicTime(expansionFactorDefault)
      !end nodes will be extended to the default min time.
    else 
      if(associated(thisNode%firstChild)) then
        childBasic => thisNode%firstChild%basic()
        timeEarliest = childBasic%time()
        !If there are children, the next time is set to the time of the children.
      else
        if (self%augmentTimeCount>1) then
          i = 1
          timeIndex = -1
          !find the next smallest time in the time snapshots and use that as the timeEarliest.
          do while (i <= size(self%augmentTimeSnapshots))
            if((thisBasicComponent%time() <= self%augmentTimeSnapshots(i)).or.(thisBasicComponent%time() <= self%augmentTimeSnapshots(i)*(1.0 + 1.0d-5)).or.(thisBasicComponent%time() <= self%augmentTimeSnapshots(i)*(1.0-1.0d-5 ))) then
              timeIndex = i - 1
              exit
            else 
              i = i + 1
            end if
          end do
          if (timeIndex == -1) then 
            timeIndex = size(self%augmentTimeSnapshots)
          end if
          if (timeIndex == 0) then
            timeEarliest = cosmologyFunctions_%cosmicTime(expansionFactorDefault)
          else 
            timeEarliest = self%augmentTimeSnapshots(timeIndex)
          end if
        else 
          !If the size of the snapshot list is reported to be one, revert to child time / default min time determination of earliest time.
          timeEarliest = cosmologyFunctions_%cosmicTime(expansionFactorDefault)
        end if
      end if
    endif
    massCutoff = self%augmentResolutionLimit
    !Build trees from the nodes above the cutoff. 
    if(thisBasicComponent%mass() > massCutoff) then
      baseNode => treeNode(thisNode%index(), newTree)
      newTree%baseNode => baseNode
      baseBasic => baseNode%basic(autoCreate=.true.)
      call baseBasic%timeSet(thisBasicComponent%time())
      call baseBasic%massSet(thisBasicComponent%mass())
      call self%mergerTreeBuilder_%timeEarliestSet(timeEarliest) 
      call self%mergerTreeBuilder_%build(newTree)
      call pruneByTime%operate(newTree)
      !Set and build tree from the base node of the new tree.
      thisNodePointer =>thisNode
      endNodeCount = augmentSilentWalk(newTree, 'endNodeCount')
      call augmentSortChildren(thisNodePointer)
      nodeChildCount = augmentChildCount(thisNodePointer)
      !check whether the newly built tree is accepted.
      treeAccepted = augmentAcceptTree(self, thisNodePointer, newTree, nodeChildCount, extendingEndNode, tolerance, timeEarliest, bestTree, bestTreeWorstFit, bestTreeOverride, multiplier, constant, scalingFactor, massCutoffScale, newNodeAboveCutoff, bestTreeNodeAboveCutoff)
      !check whether currently saved best tree can be accepted now.
      if (((bestTreeWorstFit <= tolerance .and. .not. bestTreeNodeAboveCutoff) .or. bestTreeOverride ).and.(associated(bestTree%baseNode)).and. (.not.(treeAccepted == 1)) ) then
        bestTreeWorstFit = 3.0
        if(associated(newTree%baseNode)) then
          call newTree%destroyBranch(newTree%baseNode)
        end if
        newTree%baseNode => bestTree%baseNode
        bestTree%baseNode => null()
        if ( augmentAcceptTree(self, thisNodePointer, newTree, nodeChildCount, extendingEndnode, tolerance, timeEarliest, bestTree, bestTreeWorstFit, bestTreeOverride, multiplier, constant, scalingFactor, massCutoffScale, newNodeAboveCutoff, bestTreeNodeAboveCutoff) == 1) then
          augmentBuildTreeFrom = 1
        else 
          augmentBuildTreeFrom = 0
        end if


      else if (treeAccepted  == 1 ) then
        augmentBuildTreeFrom = 1
        if (associated(bestTree%baseNode)) then
          call bestTree%destroyBranch(bestTree%baseNode)
        end if
        bestTreeWorstFit = 3.0
      else
        if (associated(newTree%baseNode)) then
         call newTree%destroyBranch(newTree%baseNode)
        end if
        augmentBuildTreeFrom = treeAccepted
      end if
    else
      augmentBuildTreeFrom = 1

    end if
    !return 1 on successful tree building attempt. 0 on unsuccessful attempt.
  end function augmentBuildTreeFrom

  subroutine augmentScaleChildren(self, thisNode, multiplier, constant, scalingFactor)
    use Galacticus_Nodes
    implicit none
    class (mergerTreeOperatorAugment), intent(inout) :: self
    type (treeNode), pointer :: thisNode, currentChild
    class (nodeComponentBasic), pointer :: currentComponentBasic, childComponentBasic
    double precision, intent(inout) :: multiplier, constant, scalingFactor
    double precision :: parentMass, childMass, parentTime, childTime

    !If the sum of the child masses of thisNode are above the mass of thisNode
    !ScaleChildren will scale down their masses to equal the mass of thisNode
    !keeping the scaling factor in scalingFactor.
    !The additional mass will be added in logaritmically, so that the total mass
    !of the descendants of each node will follow the relation:
    !Mass = multiplier * log(time) + constant 
    !At time = time_thisNode, Mass = thisNode's mass.
    !At time = time_thisNode%firstChild, Mass = sum of thisNode's childrens' mass.
    currentComponentBasic => thisNode%basic()
    parentMass = currentComponentBasic%mass()
    parentTime = currentComponentBasic%time()
    childMass = 0
    currentChild => thisNode%firstChild
    if (associated(currentChild)) then
      do while (associated(currentChild))
        childComponentBasic => currentChild%basic()
        childMass = childMass + childComponentBasic%mass()
        currentChild => currentChild%sibling
      end do 
      currentChild => thisNode%firstChild
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
        currentChild => thisNode%firstChild 
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

  subroutine augmentUnscaleChildren (self, thisNode, nodeChildCount, endNodes, multiplier, constant, scalingFactor)
    use Galacticus_Nodes
    implicit none
    class (mergerTreeOperatorAugment), intent(inout) :: self
    type (treeNode), pointer :: thisNode, currentChild, currentNode
    integer, intent (inout) :: nodeChildCount
    class (nodeComponentBasic), pointer :: currentComponentBasic, childComponentBasic
    type (treeNodeList), dimension(nodeChildCount), intent(inout) :: endNodes
    double precision, intent(inout) :: multiplier, constant, scalingFactor
    double precision :: parentMass, childMass, parentTime, childTime
    integer :: i
    !Unscales children using parameters saved from augmentScaleChildren function
    if (scalingFactor /= 1.0 .and. multiplier /= 0.0) then
      currentChild => thisNode%firstChild
      childComponentBasic => currentChild%basic()
      childMass = 0
      do while(associated(currentChild)) 
        childComponentBasic => currentChild%basic()
        childMass = childMass + childComponentBasic%mass()
        call childComponentBasic%massSet(childComponentBasic%mass()/scalingFactor)
        currentChild => currentChild%sibling
      end do
      i = 1
      currentComponentBasic => thisNode%basic()
      parentMass = currentComponentBasic%mass()
      currentChild => thisNode%firstChild
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

  integer function augmentAcceptTree(self, thisNode, tree, nodeChildCount, extendingEndNode, tolerance, timeEarliest, bestTree, bestTreeWorstFit, bestTreeOverride, multiplier, constant, scalingFactor, massCutoffScale, newNodeAboveCutoff, bestTreeNodeAboveCutoff)
    use Galacticus_Nodes
    use Merger_Trees_Builders
    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self

    type (treeNode), pointer :: thisNode, currentNode, previousNode, currentNodeSibling, nonOverlapNode, firstNonOverlap
    class (nodeComponentBasic), pointer :: currentBasicComponent, sortNodeComponentBasic, nonOverlapComponentBasic
    type (mergerTree), target :: tree
    type (mergerTree), intent (inout), target :: bestTree
    integer :: nodeChildCount, i, j, endNodesSorted, nonChildNodes, retryCount
    type (treeNodeList), dimension(nodeChildCount) :: endNodes
    logical :: treeAccepted, extendingEndNode, endNodeBuildAccept, nodeMassesAgree, currentNodeBelowAll, treeScalable
    double precision :: resolutionLimit, tolerance, massCutoff, timeEarliest, currentTreeWorstFit, falseWorstFit, unresolvedMass, treeMass, endNodeMass
    double precision, intent(inout) :: bestTreeWorstFit, multiplier, constant, scalingFactor, massCutoffScale
    logical, intent (inout) :: newNodeAboveCutoff, bestTreeNodeAboveCutoff
    logical :: bestTreeOverride
    falseWorstFit = 3.0
    newNodeAboveCutoff = .false.
    massCutoff = self%augmentResolutionLimit
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
         currentNode%hostTree => thisNode%hostTree
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
         previousNode => currentNode%parent
         currentNode => currentNode%walkTree()
         if (.not.extendingEndNode) then
           call augmentNonOverlapListAdd(nonOverlapNode, firstNonOverlap)
         end if
         nonOverlapNode => null()
    end do
    endNodesSorted = endNodesSorted - 1
    nodeMassesAgree = .true.
    currentTreeWorstFit = 0.0
    if((nodeChildCount > 0).and.(nodeChildCount <= endNodesSorted)) then
      i = 1
      currentNode => thisNode%firstChild
      do while(associated(currentNode))
        if(.not.augmentNodeComparison(currentNode, endNodes(i)%node, 2.0d0 - 2.0 * tolerance, currentTreeWorstFit)) then
          nodeMassesAgree = .false.
        end if
        i = i + 1
        currentNode => currentNode%sibling
      end do
    end if

    !if (newNodeAboveCutoff) then
    !  newNodeAboveCutoff = augmentScaleNodesAboveCutoff(self, thisNode, tree, endNodes, nodeChildCount, firstNonOverlap, resolutionLimit)
    !end if 

    if ((nodeChildCount <= endNodesSorted).and.(newNodeAboveCutoff .eqv. .false. .or. bestTreeOverride).and.nodeMassesAgree) then
      call augmentExtendNonOverlapNodes(self, tree, firstNonOverlap,massCutoff, tolerance, timeEarliest, bestTree, bestTreeWorstFit, massCutoffScale)
      endNodeMass = 0
      i = 1
      do while (i <= nodeChildCount)
        currentBasicComponent => endNodes(i)%node%basic()
        endNodeMass = endNodeMass + currentBasicComponent%mass()
        i = i + 1
      end do
      treeScalable = augmentMultiScale(self, thisNode, tree, endNodes, nodeChildCount, firstNonOverlap, tolerance, timeEarliest, bestTree, bestTreeWorstFit, unresolvedMass, treeMass, massCutoffScale)
    else 
      treeScalable = .false.
    end if
    treeAccepted = (nodeChildCount <= endNodesSorted).and.(newNodeAboveCutoff .eqv. .false.).and.treeScalable.and.nodeMassesAgree
    if (bestTreeOverride .and. associated(tree%baseNode)) then
      treeAccepted = nodeMassesAgree.and.(nodeChildCount<= endNodesSorted)
    end if

    !call augmentResetUniqueIDs(tree)
    !call augmentUnscaleChildren(self, thisNode, nodeChildCount, endNodes, multiplier, constant, scalingFactor)
    if(treeAccepted) then
      currentBasicComponent => tree%baseNode%basic()
      !write (*,*) 'Building Node ', thisNode%uniqueID(),' at time ', currentBasicComponent%time(), 'to time ', timeEarliest
      call augmentResetUniqueIDs(tree)
      call augmentUnscaleChildren(self, thisNode, nodeChildCount, endNodes, multiplier, constant, scalingFactor)  
      call augmentSimpleInsert(self, thisNode, tree, endNodes, nodeChildCount, firstNonOverlap)
      !call augmentSimpleScale(self, thisNode, tree, endNodes, nodeChildCount, firstNonOverlap)


    else if ((nodeChildCount <= endNodesSorted) .and. .not.bestTreeOverride) then
      !write (*,*) 'Updating Best Tree'
      call augmentNonOverlapReinsert(firstNonOverlap)
      if (currentTreeWorstFit < bestTreeWorstFit) then
        if(associated(bestTree%baseNode)) then
          call bestTree%destroyBranch(bestTree%baseNode)
        end if
        bestTree%baseNode => tree%baseNode
        tree%baseNode => null()
        bestTree%baseNode%hostTree => bestTree
        currentNode => bestTree%baseNode
        do while (associated(currentNode))
          currentNode%event => null()
          currentNode%hostTree => bestTree%baseNode%hostTree
          previousNode => currentNode%parent
          currentNode => currentNode%walkTree()
        end do
        bestTreeWorstFit = currentTreeWorstFit
        bestTreeNodeAboveCutoff = newNodeAboveCutoff
        !write (*,*) 'This Tree is Best Fit With Worst-- ', currentTreeWorstFit
      end if
    else 
      call augmentNonOverlapReinsert(firstNonOverlap)
      if(associated(tree%baseNode)) then
        call tree%destroyBranch(tree%baseNode)
      end if
    end if
    if (treeAccepted) then    
      augmentAcceptTree = 1
      if (currentTreeWorstFit > 0.0 .and. nodeChildCount > 1) then
        open (unit = 55, file = 'toleranceHistogram.py', position = 'append', action = 'write')
        write (55,*) currentTreeWorstFit, ','
        close (55)
      end if
    else if (.not.nodeMassesAgree) then
      augmentAcceptTree = -1
    else
      augmentAcceptTree = 0
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


 integer function augmentChildCount(thisNode)
    type (treeNode), pointer :: thisNode, currentChild
    integer childCount

    currentChild => thisNode%firstChild
    childCount = 0

    do while (associated(currentChild))
      if(associated(currentChild)) then
        childCount = childCount + 1
      end if
      currentChild => currentChild%sibling
    end do
    augmentChildCount = childCount

  end function augmentChildCount

 subroutine augmentSortChildren(thisNode)
    type (treeNode), pointer :: thisNode, currentNode, sortNode, nextNode
    class (nodeComponentBasic), pointer :: currentComponentBasic, sortComponentBasic
    double precision :: largestMass

    if(associated(thisNode%firstChild)) then    
      currentNode => thisNode%firstChild
      currentComponentBasic => currentNode%basic()
      largestMass = currentComponentBasic%mass()
      nextNode => currentNode%sibling
      currentNode%sibling => null()
      do while (associated(nextNode))
        sortNode => thisNode%firstChild
        currentNode => nextNode
        nextNode => nextNode%sibling
        currentComponentBasic => currentNode%basic()
        if (currentComponentBasic%mass() > largestmass) then
          thisNode%firstChild => currentNode
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

  subroutine augmentExtendNonOverlapNodes(self, tree, firstNonOverlap,massCutoff, tolerance, timeEarliest, bestTree, bestTreeWorstFit, massCutoffScale)
    use Merger_Trees_Builders
    use Galacticus_Nodes
    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self
    type (treeNode), pointer :: currentNode, nonOverlapNode, firstNonOverlap
    class (nodeComponentBasic), pointer :: nonOverlapComponentBasic
    type (mergerTree), target :: tree
    type (mergerTree), intent (inout), target :: bestTree
    integer :: retryCount
    double precision, intent(inout) :: bestTreeWorstFit, tolerance, timeEarliest, massCutoff, massCutoffScale 
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
        do while ((.not.(augmentBuildTreeFrom(self, nonOverlapNode, .true., tolerance, timeEarliest, bestTree, falseWorstFit, .false., falseMultiplier, falseConstant, falseScalingFactor, massCutoffScale, falseNewNodeAboveCutoff, falseBestTreeNodeAboveCutoff)==1)).and.(retryCount >0))
          retryCount = retryCount - 1
        end do
      end if
    end do
  end subroutine augmentExtendNonOverlapNodes

  logical function augmentNodeComparison(newNode, oldNode, tolerance, currentTreeWorstFit)

    use Numerical_Comparison

    implicit none   

    type (treeNode), pointer :: newNode, oldNode
    class (nodeComponentBasic), pointer :: newComponentBasic, oldComponentBasic
    double precision :: tolerance, thisFit
    double precision, intent(inout) :: currentTreeWorstFit

    newComponentBasic => newNode%basic()
    oldComponentBasic => oldNode%basic()
    augmentNodeComparison = Values_Agree(newComponentBasic%mass(), oldComponentBasic%mass(), relTol =tolerance)

    thisFit = 2.0*ABS(newComponentBasic%mass() - oldComponentBasic%mass())/ABS(newComponentBasic%mass() + oldComponentBasic%mass())
    if ( thisFit > currentTreeWorstFit) then
      currentTreeWorstFit = thisFit
    end if
  end function augmentNodeComparison

  subroutine augmentSimpleInsert(self, thisNode, tree, endNodes, nodeChildCount, firstNonOverlap)
    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self
    type (treeNode), pointer :: thisNode, currentNode, previousNode, firstNonOverlap, currentNodeSibling
    type (mergerTree), target :: tree
    integer, intent(inout) :: nodeChildCount
    integer :: i
    type (treeNodeList), dimension(nodeChildCount), intent(inout) :: endNodes
    logical :: scalableTest

    i = 0
    call augmentNonOverlapReinsert(firstNonOverlap)

    if (nodeChildCount > 0) then
      i =1
      currentNode => thisNode%firstChild

      do while (associated(currentNode))
        currentNodeSibling => currentNode%sibling
        thisNode%firstChild => currentNode%sibling
        call augmentExtendByOverLap(endNodes(i)%node, currentNode,.true., .false.)
        currentNode => currentNodeSibling
        i = i + 1
      end do
    else
      currentNode => tree%baseNode
      do while (associated(currentNode))
           currentNode%event => null()
           currentNode%hosttree => thisNode%hostTree
           previousNode => currentNode%parent
           currentNode => currentNode%walkTree()
      end do

    end if
   
    call augmentExtendByOverLap(thisNode, tree%baseNode, .false., .false.)

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

  subroutine augmentSimpleScale(self, thisNode, tree, endNodes, nodeChildCount, firstNonOverlap)
    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self
    type (treeNode), pointer :: thisNode, currentNode, previousNode, firstNonOverlap, currentNodeSibling, sortNode
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
      currentnode => thisNode%firstChild
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
      currentNode => thisNode%firstChild

      do while (associated(currentNode))
        currentNodeSibling => currentNode%sibling
        thisNode%firstChild => currentNode%sibling
        call augmentExtendByOverLap(endNodes(i)%node, currentNode,.true., .false.)
        currentNode => currentNodeSibling
        i = i + 1
      end do
    else
      currentNode => tree%baseNode
      do while (associated(currentNode))
           currentNode%event => null()
           currentNode%hosttree => thisNode%hostTree
           previousNode => currentNode%parent
           currentNode => currentNode%walkTree()
      end do

    end if

    call augmentExtendByOverLap(thisNode, tree%baseNode, .false., .false.)

  end subroutine augmentSimpleScale

  subroutine augmentScaleBranch (self, thisNode, scalingFactor)

    use Galacticus_Nodes
    implicit none
    class (mergerTreeOperatorAugment), intent (inout) :: self
    type (treeNode), pointer :: thisNode, currentChild
    class (nodeComponentBasic), pointer :: thisComponentBasic
    double precision :: scalingFactor
    thisComponentBasic => thisNode%basic()
    call thisComponentBasic%massSet(thisComponentBasic%mass() * scalingFactor)
    currentChild => thisNode%firstChild
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
    type (treeNode), pointer :: currentNode, currentChild, previousNode
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
      previousNode => currentNode%parent
      currentNode => currentNode%walkTree()
    end do

  end subroutine augmentRemoveUnresolvedMass

  subroutine augmentInsertChildMass(self, thisNode, originalChildNode, newChildNode)
    use Galacticus_Nodes
    implicit none
    class (mergerTreeOperatorAugment), intent(inout) :: self
    type (treeNode), pointer :: thisNode, originalChildNode, newChildNode, scaleNode
    class (nodeComponentBasic), pointer :: currentComponentBasic, childComponentBasic
    double precision :: multiplier, constant, scalingFactor
    double precision :: parentMass, childMass, parentTime, childTime, massDifference

    currentComponentBasic => thisNode%basic()
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

  logical function augmentScaleNodesAboveCutoff(self, thisNode, tree, endNodes, nodeChildCount, firstNonOverlap, resolutionLimit)
    use Galacticus_Nodes
    implicit none
    class (mergerTreeOperatorAugment), intent (inout) :: self
    type (treeNode), pointer :: thisNode, firstNonOverlap, currentNonOverlap, currentChild, currentEndNode, scaleNode
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
    currentChild => thisNode%firstChild
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
        currentChild => thisNode%firstChild
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


  logical function augmentMultiScale(self, thisNode, tree, endNodes, nodeChildCount, firstNonOverlap, tolerance, timeEarliest, bestTree, bestTreeWorstFit, unresolvedMass,  treeMass, massCutoffScale)

    use Galacticus_Nodes
    use Merger_Trees_Builders
    implicit none
    class  (mergerTreeOperatorAugment), intent(inout)         :: self
    type (treeNode), pointer :: thisNode, currentNode, previousNode, firstNonOverlap, currentNodeSibling, currentChildNode,sortNode, currentGrandchild, scaleNode
    type (mergerTree), target :: tree
    type (mergerTree), intent (inout), target :: bestTree
    class (nodeComponentBasic), pointer :: childComponentBasic, currentComponentBasic, sortComponentBasic
    integer, intent(inout) :: nodeChildCount
    integer :: i, retryCount
    type (treeNodeList), dimension(nodeChildCount), intent(inout) :: endNodes
    logical :: treeScaled, lastNodeFound
    double precision, intent (inout) :: tolerance, timeEarliest, bestTreeWorstFit, unresolvedMass, treeMass, massCutoffScale
    double precision :: massExcess, childNodeMass, currentMass, endNodeMass, falseWorstFit, massDifferenceScaleFactor
  
    falseWorstFit = 3.0
    call augmentNonOverlapReinsert(firstNonOverlap)
    if (nodeChildCount > 0) then
      currentNode => tree%baseNode
      currentComponentBasic => currentNode%basic()
      i = 1
      currentChildNode => thisNode%firstChild
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
        currentChildNode => thisNode%firstChild
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
