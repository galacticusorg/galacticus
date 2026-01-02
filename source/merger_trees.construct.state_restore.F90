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
  Implements a merger tree constructor class which constructs a merger tree by restoring state from file.
  !!}

  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  use :: Merger_Tree_Seeds       , only : mergerTreeSeedsClass

  public :: mergerTreeStateStore, mergerTreeStateFromFile

  !![
  <mergerTreeConstructor name="mergerTreeConstructorStateRestored">
   <description>
    A merger tree constructor class which will restore a merger tree whose complete internal state was written to file. It is
    intended primarily for debugging purposes to allow a tree to begin processing just prior to the point of failure. To use
    this method, the following procedure should be followed:
    \begin{enumerate}
     \item Identify a point in the evolution of the tree suitably close to, but before, the point of failure;
     \item Insert appropriate code into \glc\ to have it call the function to store the state of the file and then stop, e.g.:
     \begin{verbatim}
      use :: Merger_Tree_Construction, only : mergerTreeStateStore
      .
      .
      .
      if (&lt;conditions are met&gt;) then
         call mergerTreeStateStore(tree,'storedTree.dat')
         stop 'tree internal state was stored'
      end if
     \end{verbatim}
     \item Run the model ensuring that {\normalfont \ttfamily [stateFileRoot]} is set to a suitable file root name to allow the
     internal state of \glc\ to be stored;
     \item Remove the code inserted above and recompile;
     \item Run \glc\ with an input parameter file identical to the one used previously except with {\normalfont \ttfamily
     [mergerTreeConstruct]}$=${\normalfont \ttfamily stateRestore}, {\normalfont \ttfamily [stateFileRoot]} removed,
     {\normalfont \ttfamily [stateRetrieveFileRoot]} set to the value previously used for {\normalfont \ttfamily
     [stateFileRoot]} and {\normalfont \ttfamily [fileName]}$=${\normalfont \ttfamily storedTree.dat}.
    \end{enumerate}
    This should restore the tree and the internal state of \glc\ precisely from the point where they were saved and produce the
    same subsequent evolution.
    
    Note that currently this method does not support storing and restoring of trees which contain components that have more
    than one instance.
   </description>
   <runTimeFileDependencies paths="fileName"/>
  </mergerTreeConstructor>
  !!]
  type, extends(mergerTreeConstructorClass) :: mergerTreeConstructorStateRestored
     !!{
     A class implementing merger tree construction via restoring state from file.
     !!}
     private
     class(randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
     class(mergerTreeSeedsClass      ), pointer :: mergerTreeSeeds_       => null()
     type (varying_string            )          :: fileName
   contains
     final     ::              stateRestoredDestructor
     procedure :: construct => stateRestoredConstruct
  end type mergerTreeConstructorStateRestored

  interface mergerTreeConstructorStateRestored
     !!{
     Constructors for the \refClass{mergerTreeConstructorStateRestored} merger tree constructor class.
     !!}
     module procedure stateRestoredConstructorParameters
     module procedure stateRestoredConstructorInternal
  end interface mergerTreeConstructorStateRestored

  interface mergerTreeStateFromFile
    module procedure mergerTreeStateFromFileName
    module procedure mergerTreeStateFromFileUnit
  end interface mergerTreeStateFromFile
  
contains

  function stateRestoredConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeConstructorStateRestored} merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (mergerTreeConstructorStateRestored)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(randomNumberGeneratorClass        ), pointer       :: randomNumberGenerator_
    class(mergerTreeSeedsClass              ), pointer       :: mergerTreeSeeds_
    type (varying_string                    )                :: fileName

    !![
    <inputParameter>
      <name>fileName</name>
      <description>The name of the file containing the stored merger tree.</description>
      <defaultValue>var_str('storedTree.dat')</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    <objectBuilder class="mergerTreeSeeds"       name="mergerTreeSeeds_"       source="parameters"/>
    !!]
    self=mergerTreeConstructorStateRestored(fileName,randomNumberGenerator_,mergerTreeSeeds_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    <objectDestructor name="mergerTreeSeeds_"      />
    !!]
    return
  end function stateRestoredConstructorParameters

  function stateRestoredConstructorInternal(fileName,randomNumberGenerator_,mergerTreeSeeds_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeConstructorStateRestored} merger tree operator class.
    !!}
    implicit none
    type (mergerTreeConstructorStateRestored)                        :: self
    type (varying_string                    ), intent(in   )         :: fileName
    class(randomNumberGeneratorClass        ), intent(in   ), target :: randomNumberGenerator_
    class(mergerTreeSeedsClass              ), intent(in   ), target :: mergerTreeSeeds_
    !![
    <constructorAssign variables="fileName, *randomNumberGenerator_, *mergerTreeSeeds_"/>
    !!]

    return
  end function stateRestoredConstructorInternal

  subroutine stateRestoredDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeConstructorStateRestored} merger tree constructor class.
    !!}
    implicit none
    type(mergerTreeConstructorStateRestored), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    <objectDestructor name="self%mergerTreeSeeds_"      />
    !!]
    return
  end subroutine stateRestoredDestructor

  function stateRestoredConstruct(self,treeNumber,finished) result(tree)
    !!{
    Restores the state of a merger tree from file.
    !!}
    use            :: Functions_Global, only : State_Retrieve_
    use            :: Galacticus_Nodes, only : mergerTree
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    implicit none
    type   (mergerTree                        ), pointer       :: tree
    class  (mergerTreeConstructorStateRestored), intent(inout) :: self
    integer(c_size_t                          ), intent(in   ) :: treeNumber
    logical                                    , intent(  out) :: finished

    ! Only one tree to construct.
    if (treeNumber == 1_c_size_t) then
       ! Retrieve stored internal state if possible.
       call State_Retrieve_()
       ! Read the tree(s).
       allocate(tree)
       call mergerTreeStateFromFile(tree,char(self%fileName),self%randomNumberGenerator_,self%mergerTreeSeeds_)
       call self%randomSequenceNonDeterministicWarn(tree)
    else
       nullify(tree)
    end if
    finished=.not.associated(tree)
    return
  end function stateRestoredConstruct

  !![
  <functionGlobal>
    <unitName>mergerTreeStateStore</unitName>
    <type>void</type>
    <module>Galacticus_Nodes, only : mergerTree</module>
    <module>ISO_C_Binding   , only : c_size_t  </module>
    <arguments>type     (mergerTree), intent(in   ), target   :: tree              </arguments>
    <arguments>character(len=*     ), intent(in   )           :: storeFile         </arguments>
    <arguments>integer  (c_size_t  ), intent(in   ), optional :: indexOutput       </arguments>
    <arguments>logical              , intent(in   ), optional :: snapshot  , append</arguments>
  </functionGlobal>
  !!]
  subroutine mergerTreeStateStore(tree,storeFile,indexOutput,snapshot,append)
    !!{
    Store the complete internal state of a merger tree to file.
    !!}
    use            :: Functions_Global   , only : State_Store_
    use            :: Error              , only : Error_Report
    use            :: Galacticus_Nodes   , only : mergerTree                          , nodeEvent, nodeEventBuildFromRaw, treeNode
    use, intrinsic :: ISO_C_Binding      , only : c_size_t
    use            :: Kind_Numbers       , only : kind_int8
    use            :: Merger_Tree_Walkers, only : mergerTreeWalkerAllAndFormationNodes
    implicit none
    type     (mergerTree                          ), intent(in   ), target       :: tree
    character(len=*                               ), intent(in   )               :: storeFile
    integer  (c_size_t                            ), intent(in   ), optional     :: indexOutput
    logical                                        , intent(in   ), optional     :: snapshot         , append
    type     (mergerTree                          ), pointer                     :: treeCurrent
    type     (treeNode                            ), pointer                     :: node
    class    (nodeEvent                           ), pointer                     :: event
    integer  (kind_int8                           ), allocatable  , dimension(:) :: nodeIndices
    integer                                        , allocatable  , dimension(:) :: nodeCountTree
    type     (varying_string                      ), save                        :: storeFilePrevious
    type     (mergerTreeWalkerAllAndFormationNodes)                              :: treeWalker
    integer  (c_size_t                            )                              :: nodeCount
    integer                                                                      :: fileUnit         , eventCount, &
         &                                                                          treeCount        , iTree
    !![
    <optionalArgument name="snapshot" defaultsTo=".true." />
    <optionalArgument name="append"   defaultsTo=".true." />
    !!]

    ! Store internal state.
    if (snapshot_) call State_Store_()
    ! Open an output file. (Append to the old file if the file name has not changed.)
    !$omp critical (mergerTreeStateStore)
    if (append_ .and. trim(storeFile) == storeFilePrevious) then
       open(newunit=fileUnit,file=trim(storeFile),status='old'    ,form='unformatted',access='append')
    else
       storeFilePrevious=trim(storeFile)
       open(newunit=fileUnit,file=trim(storeFile),status='unknown',form='unformatted'                )
    end if
    !$omp end critical (mergerTreeStateStore)
    ! Store the output index if provided.
    if (present(indexOutput)) write (fileUnit) indexOutput
    ! Count trees and check for events attached to trees.
    treeCount   =  0
    treeCurrent => tree
    do while (associated(treeCurrent))
       if (associated(treeCurrent%event)) call Error_Report('tree events not currently supported'//{introspection:location})
       treeCount   =  treeCount           +1
       treeCurrent => treeCurrent%nextTree
    end do
    write (fileUnit) treeCount
    ! Iterate over trees, counting nodes.
    allocate(nodeCountTree(treeCount))
    nodeCountTree =  0
    iTree         =  0
    treeCurrent   => tree
    do while (associated(treeCurrent))
       iTree     =iTree+1
       treeWalker=mergerTreeWalkerAllAndFormationNodes(treeCurrent,spanForest=.false.)
       do while (treeWalker%next(node))
          nodeCountTree(iTree)=nodeCountTree(iTree)+1
       end do
       treeCurrent => treeCurrent%nextTree
    end do
    nodeCount=sum(nodeCountTree)
    write (fileUnit) nodeCountTree
    ! Allocate and populate an array of node indices in the order in which the tree will be traversed.
    allocate(nodeIndices(nodeCount))
    nodeCount =0
    treeWalker=mergerTreeWalkerAllAndFormationNodes(tree,spanForest=.true.)
    do while (treeWalker%next(node))
       nodeCount=nodeCount+1
       nodeIndices(nodeCount)=node%uniqueID()
    end do
    ! Write basic tree information.
    treeCurrent => tree
    do while (associated(treeCurrent))
       write (fileUnit) treeCurrent%index,treeCurrent%volumeWeight,treeCurrent%initializedUntil,treeCurrent%isTreeInitialized,nodeArrayPosition(treeCurrent%nodeBase,nodeIndices)
       treeCurrent => treeCurrent%nextTree
    end do
    ! Output nodes.
    treeWalker=mergerTreeWalkerAllAndFormationNodes(tree,spanForest=.true.)
    ! Iterate over all nodes.
    do while (treeWalker%next(node))
       ! Write all node information.
       ! Indices.
       write (fileUnit) node%index(),node%uniqueID()
       ! Pointers to other nodes.
       write (fileUnit)                                           &
            & nodeArrayPosition(node%parent        ,nodeIndices), &
            & nodeArrayPosition(node%firstChild    ,nodeIndices), &
            & nodeArrayPosition(node%sibling       ,nodeIndices), &
            & nodeArrayPosition(node%firstSatellite,nodeIndices), &
            & nodeArrayPosition(node%mergeTarget   ,nodeIndices), &
            & nodeArrayPosition(node%firstMergee   ,nodeIndices), &
            & nodeArrayPosition(node%siblingMergee ,nodeIndices), &
            & nodeArrayPosition(node%formationNode ,nodeIndices)
       ! Store the node.
       call node%serializeRaw(fileUnit)
       ! Store any events attached to the node.
       eventCount =  0
       event      => node%event
       do while (associated(event))
          eventCount =  eventCount+1
          event      => event%next
       end do
       write (fileUnit) eventCount
       event => node%event
       do while (associated(event))
          call event%serializeRaw(fileUnit)
          if (associated(event%node)) then
             write (fileUnit) nodeArrayPosition(event%node,nodeIndices)
          else
             write (fileUnit) -1
          end if
          event => event%next
       end do
    end do
    close(fileUnit)
    ! Destroy the temporary array of indices.
    deallocate(nodeIndices)
    return

  contains

    integer function nodeArrayPosition(node,nodeIndices)
      !!{
      Returns the position of a node in the output list given its index.
      !!}
      use :: Error             , only : Error_Report
      use :: ISO_Varying_String, only : varying_string
      use :: Kind_Numbers      , only : kind_int8
      use :: String_Handling   , only : operator(//)
      implicit none
      type   (treeNode      ), pointer     , intent(in   ) :: node
      integer(kind_int8     ), dimension(:), intent(in   ) :: nodeIndices
      type   (varying_string)                              :: message
      integer(kind_int8     )                              :: nodeIndex

      if (.not.associated(node)) then
         nodeArrayPosition=-1
      else
         nodeArrayPosition=1
         nodeIndex        =node%uniqueID()
         do while (nodeArrayPosition <= size(nodeIndices) .and. nodeIndices(nodeArrayPosition) /= nodeIndex)
            nodeArrayPosition=nodeArrayPosition+1
         end do
         if (nodeIndices(nodeArrayPosition) /= nodeIndex) then
            message="node ["
            message=message//nodeIndex//"] could not be found in merger tree"
            call Error_Report(message//{introspection:location})
         end if
      end if
      return
    end function nodeArrayPosition

  end subroutine mergerTreeStateStore

  subroutine mergerTreeStateFromFileUnit(tree,fileUnit,randomNumberGenerator_,mergerTreeSeeds_,status)
    !!{
    Read the state of a merger tree from file given the file unit.
    !!}
    use, intrinsic :: ISO_Fortran_Env , only : IOStat_End
    use            :: Error           , only : Error_Report                  , errorStatusSuccess
    use            :: Galacticus_Nodes, only : Galacticus_Nodes_Unique_ID_Set, mergerTree        , nodeEvent, nodeEventBuildFromRaw, &
          &                                    treeNode                      , treeNodeList
    use            :: Kind_Numbers    , only : kind_int8
    use            :: String_Handling , only : operator(//)
    implicit none
    type     (mergerTree                ), intent(inout), target       :: tree
    class    (randomNumberGeneratorClass), intent(inout)               :: randomNumberGenerator_
    class    (mergerTreeSeedsClass      ), intent(inout)               :: mergerTreeSeeds_
    integer                              , intent(in   )               :: fileUnit
    integer                              , intent(  out), optional     :: status
    type     (mergerTree                )               , pointer      :: treeCurrent
    type     (treeNodeList              ), allocatable  , dimension(:) :: nodes
    integer                              , allocatable  , dimension(:) :: nodeCountTree
    class    (nodeEvent                 ), pointer                     :: event                 , eventPrevious
    integer                                                            :: fileStatus            , firstChildIndex   , firstMergeeIndex   , &
         &                                                                firstSatelliteIndex   , formationNodeIndex, iNode              , &
         &                                                                mergeTargetIndex      , nodeArrayIndex    , nodeCount          , &
         &                                                                parentIndex           , siblingIndex      , siblingMergeeIndex , &
         &                                                                treeCount             , iTree             , iNodeTree          , &
         &                                                                eventCount            , iEvent            , eventNodeIndex     , &
         &                                                                status_
    integer  (kind_int8                 )                              :: nodeIndex             , nodeUniqueID      , nodeUniqueIDMaximum
    type     (varying_string            )                              :: message

    if (present(status)) status=errorStatusSuccess
    ! Read number of trees.
    read (fileUnit,iostat=status_) treeCount
    if (status_ /= 0) then
       if (status_ == IOStat_End .and. present(status)) then
          status=status_
          return
       end if
       call Error_Report('failed to read forest from file'//{introspection:location})
    end if
    ! Create trees
    if (treeCount > 1) then
       treeCurrent => tree
       do iTree=2,treeCount
          allocate(treeCurrent%nextTree)
          treeCurrent => treeCurrent%nextTree
       end do
    end if
    ! Read number of nodes.
    allocate(nodeCountTree(treeCount))
    read (fileUnit) nodeCountTree
    nodeCount=sum(nodeCountTree)
    ! Allocate a list of nodes.
    allocate(nodes(nodeCount))
    ! Iterate over trees.
    treeCurrent => tree
    iNode       =  0
    do iTree=1,treeCount
       ! Read basic tree information.
       read (fileUnit,iostat=fileStatus) treeCurrent%index,treeCurrent%volumeWeight,treeCurrent%initializedUntil,treeCurrent%isTreeInitialized,nodeArrayIndex
       ! Create nodes.
       do iNodeTree=1,nodeCountTree(iTree)
          iNode             =  iNode                         +1
          nodes(iNode)%node => treeNode(hostTree=treeCurrent)
       end do
       ! Assign the tree base node.
       treeCurrent%nodeBase => nodes(nodeArrayIndex)%node
       ! Ensure tree events are nullified.
       treeCurrent%event => null()
       ! Set pointer to the first tree in the forest.
       treeCurrent%firstTree => tree
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    ! Loop over all nodes.
    nodeUniqueIDMaximum=-1
    do iNode=1,nodeCount
       ! Read all node information.
       ! Indices.
       read (fileUnit) nodeIndex,nodeUniqueID
       call nodes(iNode)%node%indexSet   (nodeIndex   )
       call nodes(iNode)%node%uniqueIDSet(nodeUniqueID)
       ! Pointers to other nodes.
       read (fileUnit) parentIndex,firstChildIndex,siblingIndex,firstSatelliteIndex,mergeTargetIndex,firstMergeeIndex&
            &,siblingMergeeIndex,formationNodeIndex
       nodes(iNode)%node%parent         => nodePointer(parentIndex        ,nodes)
       nodes(iNode)%node%firstChild     => nodePointer(firstChildIndex    ,nodes)
       nodes(iNode)%node%sibling        => nodePointer(siblingIndex       ,nodes)
       nodes(iNode)%node%firstSatellite => nodePointer(firstSatelliteIndex,nodes)
       nodes(iNode)%node%mergeTarget    => nodePointer(mergeTargetIndex   ,nodes)
       nodes(iNode)%node%firstMergee    => nodePointer(firstMergeeIndex   ,nodes)
       nodes(iNode)%node%siblingMergee  => nodePointer(siblingMergeeIndex ,nodes)
       nodes(iNode)%node%formationNode  => nodePointer(formationNodeIndex ,nodes)
       ! Read the node.
       call nodes(iNode)%node%deserializeRaw(fileUnit)
       ! Find the highest uniqueID.
       nodeUniqueIDMaximum=max(nodeUniqueIDMaximum,nodeUniqueID)
       ! Read any events attached to the node.
       eventPrevious => null()
       read (fileUnit) eventCount
       do iEvent=1,eventCount
          ! Build an event of the correct type from the file.
          event => nodeEventBuildFromRaw(fileUnit)
          ! Read and assign any node pointer that the event may have.
          read (fileUnit) eventNodeIndex
          event%node => nodePointer(eventNodeIndex,nodes)
          ! Link the event to the node.
          if (iEvent == 1) then
             nodes(iNode)%node%event => event
          else
             eventPrevious%next => event
          end if
          eventPrevious => event
          event         => null()
       end do
    end do
    ! Set the global maximum unique ID to the maximum found.
    call Galacticus_Nodes_Unique_ID_Set(nodeUniqueIDMaximum)
    ! Perform sanity checks.
    do iNode=1,nodeCount
       if (associated(nodes(iNode)%node%firstChild)) then
          if (.not.associated(nodes(iNode)%node,nodes(iNode)%node%firstChild%parent)) then
             message="child's parent is not self"
             message=message//char(10)//" -> self                : "//nodes(iNode)%node                  %uniqueID()
             message=message//char(10)//" -> self->child         : "//nodes(iNode)%node%firstChild       %uniqueID()
             message=message//char(10)//" -> self->child->parent : "//nodes(iNode)%node%firstChild%parent%uniqueID()
             call Error_Report(message//{introspection:location})
          end if
       end if
    end do
    ! Destroy the list of nodes.
    deallocate(nodes)
    ! Restart the random number sequence.
    allocate(tree%randomNumberGenerator_,mold=randomNumberGenerator_)
    !$omp critical(mergerTreeStateRestoreDeepCopyReset)
    !![
    <deepCopyReset variables="randomNumberGenerator_"/>
    <deepCopy source="randomNumberGenerator_" destination="tree%randomNumberGenerator_"/>
    <deepCopyFinalize variables="tree%randomNumberGenerator_"/>
    !!]
    !$omp end critical(mergerTreeStateRestoreDeepCopyReset)
    call mergerTreeSeeds_%set(tree)
    return

  contains

    function nodePointer(nodeArrayIndex,nodes)
      !!{
      Return a pointer to a node, given its position in the array of nodes. Return a null pointer if the array index is $-1$.
      !!}
      use :: Galacticus_Nodes, only : treeNode, treeNodeList
      type   (treeNode    ), pointer                     :: nodePointer
      integer                            , intent(in   ) :: nodeArrayIndex
      type   (treeNodeList), dimension(:), intent(in   ) :: nodes

      if (nodeArrayIndex == -1) then
         nodePointer => null()
      else
         nodePointer => nodes(nodeArrayIndex)%node
      end if
      return
    end function nodePointer

  end subroutine mergerTreeStateFromFileUnit

  subroutine mergerTreeStateFromFileName(tree,fileName,randomNumberGenerator_,mergerTreeSeeds_,status,deleteAfterRead)
    !!{
    Read the state of a merger tree from file given a file name.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : mergerTree
    implicit none
    type     (mergerTree                ), intent(inout), target   :: tree
    character(len=*                     ), intent(in   )           :: fileName
    class    (randomNumberGeneratorClass), intent(inout)           :: randomNumberGenerator_
    class    (mergerTreeSeedsClass      ), intent(inout)           :: mergerTreeSeeds_
    integer                              , intent(  out), optional :: status
    logical                              , intent(in   ), optional :: deleteAfterRead
    integer                                                        :: fileUnit              , ioStatus
    !![
    <optionalArgument name="deleteAfterRead" defaultsTo=".false." />
    !!]

    ! Open the file.
    open(newUnit=fileUnit,file=fileName,status='old',form='unformatted',iostat=ioStatus)
    if (ioStatus /= 0) call Error_Report('unable to open file "'//trim(fileName)//'"'//{introspection:location})
    ! Perform the read.
    call mergerTreeStateFromFile(tree,fileUnit,randomNumberGenerator_,mergerTreeSeeds_,status)
    ! Close the file.
    if (deleteAfterRead_) then
       close(fileUnit,status='delete')
    else
       close(fileUnit                )
    end if
    return
  end subroutine mergerTreeStateFromFileName
