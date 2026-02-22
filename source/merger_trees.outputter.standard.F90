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
  Implements the standard class for outputting merger trees.
  !!}

  use :: Cosmology_Functions               , only : cosmologyFunctions   , cosmologyFunctionsClass
  use :: Galactic_Filters                  , only : galacticFilter       , galacticFilterClass
  use :: IO_HDF5                           , only : hdf5Object
  use :: Kind_Numbers                      , only : kind_int8
  use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger, outputPropertyDouble
  use :: Node_Property_Extractors          , only : nodePropertyExtractor, nodePropertyExtractorClass
  use :: Locks                             , only : ompLock

  type outputGroup
     !!{
     Type used for output group information.
     !!}
     logical             :: doubleAttributesWritten, integerAttributesWritten, &
          &                 opened
     type   (hdf5Object) :: hdf5Group              , nodeDataGroup
   contains
     !![
     <methods>
       <method method="assignment(=)" description="Assign outputGroup objects."/>
     </methods>
     !!]
     procedure :: outputGroupAssign
     generic   :: assignment(=)     =>  outputGroupAssign
  end type outputGroup

  ! Parameters controlling the size of output data.
  !![
  <scoping>
   <module variables="standardBufferSizeIncrement"/>
  </scoping>
  !!]
  integer         , parameter :: standardOutputGroupsIncrement=10     , standardBufferSizeIncrement=1024

  ! Lock object used to ensure merger trees are output contiguously.
  type   (ompLock)            :: treeLock
  logical                     :: treeLockInitialized          =.false.
  
  !![
  <mergerTreeOutputter name="mergerTreeOutputterStandard">
   <description>The standard merger tree outputter.</description>
   <stateStorable>
    <restoreTo variables="outputsGroupOpened"                                                                                            state=".false."                    />
    <restoreTo variables="outputGroupsCount  , doublePropertiesWritten, integerPropertiesWritten, doubleBufferCount, integerBufferCount" state="0"                          />
    <restoreTo variables="doublePropertyCount, integerPropertyCount"                                                                     state="-1"                         />
    <restoreTo variables="integerBufferSize  , doubleBufferSize"                                                                         state="standardBufferSizeIncrement"/>
    <exclude   variables="doubleProperty     , integerProperty"                                                                                                             />
   </stateStorable>
  </mergerTreeOutputter>
  !!]
  type, extends(mergerTreeOutputterClass) :: mergerTreeOutputterStandard
     !!{
     Implementation of the standard merger tree outputter.
     !!}
     private
     logical                                                                   :: outputReferences
     type            (varying_string              )                            :: outputsGroupName
     type            (hdf5Object                  )                            :: outputsGroup
     logical                                                                   :: outputsGroupOpened
     integer         (c_size_t                    )                            :: outputGroupsCount       =  0_c_size_t
     integer                                                                   :: doublePropertyCount                  , integerPropertyCount
     integer                                                                   :: doublePropertiesWritten              , integerPropertiesWritten
     integer                                                                   :: doubleBufferCount                    , integerBufferCount
     integer                                                                   :: doubleScalarCount                    , integerScalarCount
     integer                                                                   :: integerBufferSize                    , doubleBufferSize
     type            (outputPropertyInteger       ), allocatable, dimension(:) :: integerProperty
     type            (outputPropertyDouble        ), allocatable, dimension(:) :: doubleProperty
     type            (outputGroup                 ), allocatable, dimension(:) :: outputGroups
     class           (galacticFilterClass         ), pointer                   :: galacticFilter_         => null()
     class           (cosmologyFunctionsClass     ), pointer                   :: cosmologyFunctions_     => null()
     class           (nodePropertyExtractorClass  ), pointer                   :: nodePropertyExtractor_  => null()
   contains
     !![
     <methods>
       <method description="Make an group in the \glc\ file in which to store {\normalfont \ttfamily tree}." method="makeGroup"             />
       <method description="Dump the contents of the integer properties buffer to the \glc\ output file."    method="dumpIntegerBuffer"     />
       <method description="Dump the contents of the double properties buffer to the \glc\ output file."     method="dumpDoubleBuffer"      />
       <method description="Extend the size of the integer buffer."                                          method="extendIntegerBuffer"   />
       <method description="Extend the size of the double buffer."                                           method="extendDoubleBuffer"    />
       <method description="Count up the number of properties that will be output."                          method="propertiesCount"       />
       <method description="Allocate buffers for storage of properties."                                     method="buffersAllocate"       />
       <method description="Set names for the properties."                                                   method="propertyNamesEstablish"/>
       <method description="Create a group in which to store this output."                                   method="outputGroupCreate"     />
       <method description="Perform output of a single node."                                                method="output"                />
     </methods>
     !!]
     final     ::                           standardDestructor
     procedure :: outputTree             => standardOutputTree
     procedure :: outputNode             => standardOutputNode
     procedure :: makeGroup              => standardMakeGroup
     procedure :: dumpIntegerBuffer      => standardDumpIntegerBuffer
     procedure :: extendIntegerBuffer    => standardExtendIntegerBuffer
     procedure :: extendDoubleBuffer     => standardExtendDoubleBuffer
     procedure :: dumpDoubleBuffer       => standardDumpDoubleBuffer
     procedure :: propertiesCount        => standardPropertiesCount
     procedure :: buffersAllocate        => standardBuffersAllocate
     procedure :: propertyNamesEstablish => standardPropertyNamesEstablish
     procedure :: outputGroupCreate      => standardOutputGroupCreate
     procedure :: output                 => standardOutput
  end type mergerTreeOutputterStandard

  interface mergerTreeOutputterStandard
     !!{
     Constructors for the \refClass{mergerTreeOutputterStandard} merger tree outputter.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface mergerTreeOutputterStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeOutputterStandard} merger tree outputter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (mergerTreeOutputterStandard)                :: self
    type   (inputParameters            ), intent(inout) :: parameters
    class  (galacticFilterClass        ), pointer       :: galacticFilter_
    class  (cosmologyFunctionsClass    ), pointer       :: cosmologyFunctions_
    class  (nodePropertyExtractorClass ), pointer       :: nodePropertyExtractor_
    logical                                             :: outputReferences
    type   (varying_string             )                :: outputsGroupName

    !![
    <inputParameter>
      <name>outputsGroupName</name>
      <source>parameters</source>
      <defaultValue>var_str('Outputs')</defaultValue>
      <description>The name of the HDF5 group to which outputs will be written.</description>
    </inputParameter>
    <inputParameter>
      <name>outputReferences</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not references to individual merger tree datasets should be output.</description>
    </inputParameter>
    <objectBuilder class="galacticFilter"        name="galacticFilter_"        source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    !!]
    self=mergerTreeOutputterStandard(outputsGroupName,outputReferences,galacticFilter_,cosmologyFunctions_,nodePropertyExtractor_)
    !![
    <inputParametersValidate source="parameters"   />
    <objectDestructor name="galacticFilter_"       />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="nodePropertyExtractor_"/>
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(outputsGroupName,outputReferences,galacticFilter_,cosmologyFunctions_,nodePropertyExtractor_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeOutputterStandard} merger tree outputter class.
    !!}
    implicit none
    type   (mergerTreeOutputterStandard)                        :: self
    type   (varying_string             ), intent(in   )         :: outputsGroupName
    class  (galacticFilterClass        ), intent(in   ), target :: galacticFilter_
    class  (cosmologyFunctionsClass    ), intent(in   ), target :: cosmologyFunctions_
    class  (nodePropertyExtractorClass ), intent(in   ), target :: nodePropertyExtractor_
    logical                             , intent(in   )         :: outputReferences
    !![
    <constructorAssign variables="outputsGroupName, outputReferences, *galacticFilter_, *cosmologyFunctions_, *nodePropertyExtractor_"/>
    !!]

    self%outputsGroupOpened      =.false.
    self%outputGroupsCount       = 0
    self%doublePropertyCount     =-1
    self%integerPropertyCount    =-1
    self%doubleScalarCount       =-1
    self%integerScalarCount      =-1
    self%doublePropertiesWritten = 0
    self%integerPropertiesWritten= 0
    self%doubleBufferCount       = 0
    self%integerBufferCount      = 0
    self%integerBufferSize       =standardBufferSizeIncrement
    self%doubleBufferSize        =standardBufferSizeIncrement
    allocate(self%integerProperty(0))
    allocate(self% doubleProperty(0))
    !$omp critical(mergerTreeOutputterStandardInitialize)
    if (.not.treeLockInitialized) then
       treeLock           =ompLock()
       treeLockInitialized=.true.
    end if
    !$omp end critical(mergerTreeOutputterStandardInitialize)
    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeOutputterStandard} merger tree outputter class.
    !!}
    implicit none
    type(mergerTreeOutputterStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"       />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%nodePropertyExtractor_"/>
    !!]
    return
  end subroutine standardDestructor

  subroutine standardOutputTree(self,tree,indexOutput,time)
    !!{
    Write properties of nodes in {\normalfont \ttfamily tree} to the \glc\ output file.
    !!}
    use            :: Error              , only : Error_Report
    use            :: Galacticus_Nodes   , only : mergerTree              , nodeComponentBasic, treeNode
    use            :: HDF5_Access        , only : hdf5Access
    use            :: IO_HDF5            , only : hdf5Object
    use, intrinsic :: ISO_C_Binding      , only : c_size_t
    use            :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    implicit none
    class           (mergerTreeOutputterStandard), intent(inout)          :: self
    type            (mergerTree                 ), intent(inout), target  :: tree
    integer         (c_size_t                   ), intent(in   )          :: indexOutput
    double precision                             , intent(in   )          :: time
    type            (treeNode                   )               , pointer :: node
    integer         (kind=hsize_t               ), dimension(1)           :: referenceLength , referenceStart
    class           (nodeComponentBasic         )               , pointer :: basic
    type            (mergerTree                 )               , pointer :: currentTree
    type            (mergerTreeWalkerAllNodes   )                         :: treeWalker
    integer                                                               :: iProperty
    type            (hdf5Object                 )                         :: toDataset

    ! Create an output group.
    call self%outputGroupCreate(indexOutput,time)
    ! Iterate over trees.
    currentTree => tree
    do while (associated(currentTree))
       ! Get the base node of the tree.
       node => currentTree%nodeBase
       ! Skip empty trees.
       if (associated(node)) then
          ! Initialize output buffers.
          ! Count up the number of properties to be output.
          call self%propertiesCount        (time,node)
          ! Ensure buffers are allocated for temporary property storage.
          call self%buffersAllocate        (indexOutput)
          ! Get names for all properties to be output.
          call self%propertyNamesEstablish(time,node)
          ! Loop over all nodes in the tree.
          self%integerPropertiesWritten=0
          self%doublePropertiesWritten =0
          treeWalker=mergerTreeWalkerAllNodes(currentTree)
          do while (treeWalker%next(node))
             ! Get the basic component.
             basic => node%basic()
             ! Perform our output for nodes at the output time.
             if (basic%time() == time) call self%output(node,indexOutput,time)
          end do
          ! Finished output.
          if     (                                                                     &
               &   (                                                                   &
               &     (self%integerPropertyCount > 0 .and. self%integerBufferCount > 0) &
               &    .or.                                                               &
               &     (self% doublePropertyCount > 0 .and. self% doubleBufferCount > 0) &
               &   )                                                                   &
               &  .and.                                                                &
               &   .not.treeLock%ownedByThread()                                       &
               & ) call treeLock%set()
          if         (self%integerPropertyCount > 0 .and. self%integerBufferCount > 0) call self%dumpIntegerBuffer(indexOutput)
          if         (self% doublePropertyCount > 0 .and. self% doubleBufferCount > 0) call self%dumpDoubleBuffer (indexOutput)
          ! Compute the start and length of regions to reference.
          if (.not.treeLock%ownedByThread()) call treeLock%set()
          !$ call hdf5Access%set()
          referenceLength(1)=max(self%integerPropertiesWritten,self%doublePropertiesWritten)
          if      (allocated(self%integerProperty).and.size(self%integerProperty) > 0.and.self%outputGroups(indexOutput)%nodeDataGroup%hasDataset(self%integerProperty(1)%name)) then
             toDataset=self%outputGroups(indexOutput)%nodeDataGroup%openDataset(self%integerProperty(1)%name)
             referenceStart(1)=toDataset%size(1)-referenceLength(1)
          else if (allocated(self% doubleProperty).and.size(self% doubleProperty) > 0.and.self%outputGroups(indexOutput)%nodeDataGroup%hasDataset(self% doubleProperty(1)%name)) then
             toDataset=self%outputGroups(indexOutput)%nodeDataGroup%openDataset(self% doubleProperty(1)%name)
             referenceStart(1)=toDataset%size(1)-referenceLength(1)
          else
             ! Datasets do not yet exist therefore the start reference must be zero.
             referenceStart(1)=0
          end if
          ! Create references to the datasets if requested.
          if (self%outputReferences) then
             ! Ensure that a group has been made for this merger tree.
             call self%makeGroup(currentTree,indexOutput)
             ! Create references for this tree.
             if (self%integerPropertyCount > 0 .and. self%integerPropertiesWritten > 0) then
                do iProperty=1,self%integerPropertyCount
                   toDataset=self%outputGroups(indexOutput)%nodeDataGroup%openDataset(self%integerProperty(iProperty)%name)
                   call currentTree%hdf5Group%createReference1D(toDataset,self%integerProperty(iProperty)%name,referenceStart+1,referenceLength)
                end do
             end if
             if (self%doublePropertyCount > 0  .and. self%doublePropertiesWritten  > 0) then
                do iProperty=1,self%doublePropertyCount
                   toDataset=self%outputGroups(indexOutput)%nodeDataGroup%openDataset(self%doubleProperty(iProperty)%name)
                   call currentTree%hdf5Group%createReference1D(toDataset,self%doubleProperty(iProperty)%name,referenceStart+1,referenceLength)
                end do
             end if
          end if
          ! Store the start position and length of the node data for this tree, along with its volume weight.
          call self%outputGroups(indexOutput)%hdf5Group%writeDataset([currentTree%index]                        ,"mergerTreeIndex"     ,"Index of each merger tree."                                  ,appendTo=.true.)
          call self%outputGroups(indexOutput)%hdf5Group%writeDataset(referenceStart                             ,"mergerTreeStartIndex","Index in nodeData datasets at which each merger tree begins.",appendTo=.true.)
          call self%outputGroups(indexOutput)%hdf5Group%writeDataset(referenceLength                            ,"mergerTreeCount"     ,"Number of nodes in nodeData datasets for each merger tree."  ,appendTo=.true.)
          call self%outputGroups(indexOutput)%hdf5Group%writeDataset([currentTree%volumeWeight                 ],"mergerTreeWeight"    ,"Number density of each tree [Mpc⁻³]."                        ,appendTo=.true.)
          call self%outputGroups(indexOutput)%hdf5Group%writeDataset([currentTree%randomNumberGenerator_%seed()],"mergerTreeSeed"      ,"Random number seed of each merger tree."                     ,appendTo=.true.)
          !$ call hdf5Access%unset()
          if (treeLock%ownedByThread()) call treeLock%unset()
       end if
       ! Skip to the next tree.
       currentTree => currentTree%nextTree
    end do
    return
  end subroutine standardOutputTree

  subroutine standardOutputNode(self,node,indexOutput)
    !!{
    Output a single node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class  (mergerTreeOutputterStandard), intent(inout) :: self
    type   (treeNode                   ), intent(inout) :: node
    integer(c_size_t                   ), intent(in   ) :: indexOutput
    class  (nodeComponentBasic         ), pointer       :: basic
    
    ! Create an output group.
    basic => node%basic()
    call self%outputGroupCreate(indexOutput,basic%time())
    ! Initialize output buffers.
    ! Count up the number of properties to be output.
    call self%propertiesCount       (basic%time(),node)
    ! Ensure buffers are allocated for temporary property storage.
    call self%buffersAllocate       (indexOutput      )
    ! Get names for all properties to be output.
    call self%propertyNamesEstablish(basic%time(),node)
    self%integerPropertiesWritten=0
    self%doublePropertiesWritten =0
    ! Perform our output.
    call self%output(node,indexOutput,basic%time())
    ! Finished output.
    if     (                                                                     &
         &   (                                                                   &
         &     (self%integerPropertyCount > 0 .and. self%integerBufferCount > 0) &
         &    .or.                                                               &
         &     (self% doublePropertyCount > 0 .and. self% doubleBufferCount > 0) &
         &   )                                                                   &
         &  .and.                                                                &
         &   .not.treeLock%ownedByThread()                                       &
         & ) call treeLock%set()
    if (self%integerPropertyCount > 0 .and. self%integerBufferCount > 0) call self%dumpIntegerBuffer(indexOutput)
    if (self% doublePropertyCount > 0 .and. self% doubleBufferCount > 0) call self%dumpDoubleBuffer (indexOutput)
    if (treeLock%ownedByThread()) call treeLock%unset()
    return
  end subroutine standardOutputNode

  subroutine standardOutput(self,node,indexOutput,time)
    !!{
    Output the provided node.
    !!}
    use :: Calculations_Resets     , only : Calculations_Reset
    use :: Multi_Counters          , only : multiCounter
    use :: Node_Property_Extractors, only : elementTypeDouble         , elementTypeInteger       , nodePropertyExtractorIntegerScalar, nodePropertyExtractorIntegerTuple, &
         &                                  nodePropertyExtractorMulti, nodePropertyExtractorNull, nodePropertyExtractorScalar       , nodePropertyExtractorTuple       , &
         &                                  nodePropertyExtractorArray, nodePropertyExtractorList, nodePropertyExtractorList2D       , nodePropertyExtractorIntegerList
    use :: Poly_Ranks              , only : polyRankInteger           , polyRankDouble           , assignment(=)
    implicit none
    class           (mergerTreeOutputterStandard), intent(inout)                   :: self
    type            (treeNode                   ), intent(inout)                   :: node
    integer         (c_size_t                   ), intent(in   )                   :: indexOutput
    double precision                             , intent(in   )                   :: time
    integer         (kind_int8                  ), allocatable  , dimension(:    ) :: integerTuple
    double precision                             , allocatable  , dimension(:    ) :: doubleTuple
    integer         (kind_int8                  ), allocatable  , dimension(:,:  ) :: integerArray
    double precision                             , allocatable  , dimension(:,:  ) :: doubleArray
    double precision                             , allocatable  , dimension(:,:,:) :: doubleArray2D
    type            (polyRankInteger            ), allocatable  , dimension(:    ) :: integerProperties
    type            (polyRankDouble             ), allocatable  , dimension(:    ) :: doubleProperties
    integer                                      , allocatable  , dimension(:    ) :: doubleRanks      , integerRanks
    integer         (c_size_t                   ), allocatable  , dimension(:    ) :: shape_
    integer                                                                        :: doubleProperty   , integerProperty, &
         &                                                                            i
    logical                                                                        :: nodePassesFilter
    type            (multiCounter               )                                  :: instance

    ! Reset calculations (necessary in case the last node to be evolved is the first one we output, in which case
    ! calculations would not be automatically reset because the node unique ID will not have changed).
    call Calculations_Reset (node)
    ! Test whether this node passes all output filters.
    nodePassesFilter=self%galacticFilter_%passes(node)
    if (nodePassesFilter) then
       ! Ensure output buffers are allocated.
       do i=1,self%integerScalarCount
          if (.not.allocated(self%integerProperty(i)%scalar)) allocate(self%integerProperty(i)%scalar(self%integerBufferSize))
       end do
       do i=1,self%doubleScalarCount
          if (.not.allocated(self%doubleProperty (i)%scalar)) allocate(self%doubleProperty (i)%scalar(self%doubleBufferSize))
       end do
       ! Initialize the instance counter.
       instance=multiCounter([1_c_size_t])
       call self%nodePropertyExtractor_%addInstances(node,instance)
       do while (instance%increment())
          if (self%integerPropertyCount > 0) then
             integerProperty=0
             self%integerBufferCount=self%integerBufferCount+1
          end if
          if (self%doublePropertyCount > 0) then
             doubleProperty=0
             self%doubleBufferCount=self%doubleBufferCount+1
          end if
          ! Populate the output buffers with properties. We first populate with any "extra" properties that may be
          ! being computed, and then call the standard treeNode output method to populate with all "standard"
          ! properties.
          call node%output(integerProperty,self%integerBufferCount,self%integerProperty,doubleProperty,self%doubleBufferCount,self%doubleProperty,time,instance)
          ! Handle any extracted properties.
          select type (extractor_ => self%nodePropertyExtractor_)
          type  is (nodePropertyExtractorNull         )
             ! Null extractor - simply ignore.
             class is (nodePropertyExtractorScalar       )
                ! Scalar property extractor - extract and store the value.
             if    (.not.allocated(self%doubleProperty (doubleProperty +1)%scalar)) allocate(self%doubleProperty(doubleProperty +1)%scalar(                      self%doubleBufferSize))
             self   %doubleProperty (doubleProperty +1)%scalar(self%doubleBufferCount )=extractor_      %extract     (                  node     ,instance)
             doubleProperty                                                            =+doubleProperty                                                     &
                  &                                                                     +1
             class is (nodePropertyExtractorTuple        )
                ! Tuple property extractor - extract and store the values.
             doubleTuple =extractor_%extract       (node,time,instance)
             do i=1,+extractor_%elementCount(                  time)
                if (.not.allocated(self%doubleProperty (doubleProperty +i)%scalar)) allocate(self%doubleProperty (doubleProperty +i)%scalar(                      self%doubleBufferSize))
                self%doubleProperty (doubleProperty +i)%scalar(self%doubleBufferCount )=doubleTuple (  i)
             end do
             deallocate(doubleTuple )
             doubleProperty                                                            =+doubleProperty                                                     &
                  &                                                                     +extractor_     %elementCount(                       time         )
             class is (nodePropertyExtractorIntegerScalar)
                ! Integer scalar property extractor - extract and store the value.
             if    (.not.allocated(self%integerProperty(integerProperty+1)%scalar)) allocate(self%integerProperty(integerProperty+1)%scalar(                      self%integerBufferSize))
             self   %integerProperty(integerProperty+1)%scalar(self%integerBufferCount)=extractor_      %extract     (                  node,time,instance)
             integerProperty                                                           =+integerProperty                                                    &
                  &                                                                     +1
             class is (nodePropertyExtractorIntegerTuple )
                ! Integer tuple property extractor - extract and store the values.
             integerTuple=extractor_%extract       (node,time,instance)
             do i=1,extractor_%elementCount(                   time)
                if (.not.allocated(self%integerProperty(integerProperty+i)%scalar)) allocate(self%integerProperty(integerProperty+i)%scalar(                      self%integerBufferSize))
                self%integerProperty(integerProperty+i)%scalar(self%integerBufferCount)=integerTuple(  i)
             end do
             deallocate(integerTuple)
             integerProperty                                                           =+integerProperty                                                    &
                  &                                                                     +extractor_     %elementCount(                       time         )
             class is (nodePropertyExtractorArray        )
                ! Array property extractor - extract and store the values.
             doubleArray =extractor_%extract       (node,time,instance)
             do i=1,+extractor_%elementCount(                  time)
                if (     allocated(self%doubleProperty (doubleProperty +i)%scalar))                          deallocate(self%doubleProperty (doubleProperty +i)%scalar)
                if (     allocated(self%doubleProperty (doubleProperty +i)%rank1 )) then 
                   if (size(self%doubleProperty (doubleProperty +i)%rank1,dim=1) /= size(doubleArray,dim=1)) deallocate(self%doubleProperty (doubleProperty +i)%rank1 )
                end if
                if (.not.allocated(self%doubleProperty (doubleProperty +i)%rank1)) allocate(self%doubleProperty (doubleProperty +i)%rank1(size(doubleArray,dim=1),self%doubleBufferSize))
                self%doubleProperty (doubleProperty +i)%rank1(:,self%doubleBufferCount)=doubleArray (:,i)
             end do
             deallocate(doubleArray )
             doubleProperty                                                            =+doubleProperty                                                     &
                  &                                                                     +extractor_     %elementCount(                       time         )
             class is (nodePropertyExtractorList         )
                ! List property extractor - extract and store the values.
             doubleArray =extractor_%extract       (node     ,instance)
             do i=1,+extractor_%elementCount(                      )
                if (     allocated(self%doubleProperty (doubleProperty +i)%scalar     ))                          deallocate(self%doubleProperty (doubleProperty +i)%scalar)
                if (     allocated(self%doubleProperty (doubleProperty +i)%rank1      ))                          deallocate(self%doubleProperty (doubleProperty +i)%rank1 )
                if (.not.allocated(self%doubleProperty (doubleProperty +i)%rank1VarLen)) allocate(self%doubleProperty (doubleProperty +i)%rank1VarLen(self%doubleBufferSize))
                if (associated(self%doubleProperty (doubleProperty +i)%rank1VarLen (self%doubleBufferCount )%row)) then
                   if (size(self%doubleProperty (doubleProperty +i)%rank1VarLen (self%doubleBufferCount )%row) /= size(doubleArray,dim=1)) deallocate(self%doubleProperty (doubleProperty +i)%rank1VarLen (self%doubleBufferCount )%row)
                end if
                if (.not.associated(self%doubleProperty (doubleProperty +i)%rank1VarLen (self%doubleBufferCount )%row)) then
                   allocate(self%doubleProperty (doubleProperty +i)%rank1VarLen (self%doubleBufferCount )%row(size(doubleArray,dim=1)))
                end if
                self%doubleProperty (doubleProperty +i)%rank1VarLen(self%doubleBufferCount)%row=doubleArray(:,i)
             end do
             deallocate(doubleArray )
             doubleProperty                                                            =+doubleProperty                                                     &
                  &                                                                     +extractor_     %elementCount(                                    )
             class is (nodePropertyExtractorIntegerList  )
                ! Integer list property extractor - extract and store the values.
             integerArray=extractor_%extract       (node     ,instance)
             do i=1,+extractor_%elementCount(                      )
                if (     allocated(self%integerProperty(integerProperty +i)%scalar     ))                          deallocate(self%integerProperty(integerProperty+i)%scalar)
                if (     allocated(self%integerProperty(integerProperty +i)%rank1      ))                          deallocate(self%integerProperty(integerProperty+i)%rank1 )
                if (.not.allocated(self%integerProperty(integerProperty +i)%rank1VarLen)) allocate(self%integerProperty (integerProperty +i)%rank1VarLen(self%integerBufferSize))
                if (associated(self%integerProperty(integerProperty+i)%rank1VarLen(self%integerBufferCount)%row)) then
                   if (size(self%integerProperty(integerProperty+i)%rank1VarLen(self%integerBufferCount)%row) /= size(integerArray,dim=1)) deallocate(self%integerProperty(integerProperty+i)%rank1VarLen(self%integerBufferCount)%row)
                end if
                if (.not.associated(self%integerProperty(integerProperty+i)%rank1VarLen(self%integerBufferCount)%row)) then
                   allocate(self%integerProperty(integerProperty+i)%rank1VarLen(self%integerBufferCount)%row(size(integerArray,dim=1)))
                end if
                self%integerProperty(integerProperty+i)%rank1VarLen(self%integerBufferCount)%row=integerArray(:,i)
             end do
             deallocate(integerArray)
             integerProperty                                                           =+integerProperty                                                    &
                  &                                                                     +extractor_     %elementCount(                                    )
             class is (nodePropertyExtractorList2D       )
                ! 2D list property extractor - extract and store the values.
             doubleArray2D =extractor_%extract       (node     ,instance)
             do i=1,+extractor_%elementCount(                      )
                if (     allocated(self%doubleProperty (doubleProperty +i)%scalar     ))                          deallocate(self%doubleProperty (doubleProperty +i)%scalar     )
                if (     allocated(self%doubleProperty (doubleProperty +i)%rank1      ))                          deallocate(self%doubleProperty (doubleProperty +i)%rank1      )
                if (     allocated(self%doubleProperty (doubleProperty +i)%rank1VarLen))                          deallocate(self%doubleProperty (doubleProperty +i)%rank1VarLen)
                if (.not.allocated(self%doubleProperty (doubleProperty +i)%rank2VarLen)) allocate(self%doubleProperty (doubleProperty +i)%rank2VarLen(self%doubleBufferSize))
                if (associated(self%doubleProperty (doubleProperty +i)%rank2VarLen (self%doubleBufferCount )%row)) then
                   if     (                                                                                                                            &
                        &   size(self%doubleProperty (doubleProperty +i)%rank2VarLen (self%doubleBufferCount )%row,dim=1) /= size(doubleArray2D,dim=1) &
                        &  .or.                                                                                                                        &
                        &   size(self%doubleProperty (doubleProperty +i)%rank2VarLen (self%doubleBufferCount )%row,dim=2) /= size(doubleArray2D,dim=2) &
                        & ) deallocate(self%doubleProperty (doubleProperty +i)%rank2VarLen (self%doubleBufferCount )%row)
                end if
                if (.not.associated(self%doubleProperty (doubleProperty +i)%rank2VarLen (self%doubleBufferCount )%row)) then
                   allocate(self%doubleProperty (doubleProperty +i)%rank2VarLen (self%doubleBufferCount )%row(size(doubleArray2D,dim=1),size(doubleArray2D,dim=2)))
                end if
                self%doubleProperty (doubleProperty +i)%rank2VarLen(self%doubleBufferCount)%row=doubleArray2D(:,:,i)
             end do
             deallocate(doubleArray2D )
             doubleProperty                                                            =+doubleProperty                                                     &
                  &                                                                     +extractor_     %elementCount(                                    )
             class is (nodePropertyExtractorMulti        )
                ! Multi property extractor - extract and store the values.
             doubleProperties =extractor_%extractDouble (node,time,instance,doubleRanks)
             do i=1,extractor_%elementCount(elementTypeDouble ,time)
                select case (doubleRanks(i))
                case (0)
                   ! Scalar property.
                   if (     allocated(self%doubleProperty (doubleProperty +i)%rank1 ))            deallocate(self%doubleProperty (doubleProperty +i)%rank1 )
                   if (.not.allocated(self%doubleProperty (doubleProperty +i)%scalar)) then
                      allocate(self%doubleProperty(doubleProperty+i)%scalar(          self%doubleBufferSize))
                   end if
                   self%doubleProperty (doubleProperty +i)%scalar(  self%doubleBufferCount )=doubleProperties(i)
                case (1)
                   ! Rank-1 array property.
                   if (     allocated(self%doubleProperty(doubleProperty +i)%scalar))             deallocate(self%doubleProperty (doubleProperty +i)%scalar)
                   if (     allocated(self%doubleProperty(doubleProperty +i)%rank1 )) then
                      shape_=doubleProperties(i)%shape()
                      if (size(self%doubleProperty (doubleProperty +i)%rank1,dim=1) /= shape_(1)) deallocate(self%doubleProperty (doubleProperty +i)%rank1 )
                      deallocate(shape_)
                   end if
                   if (.not.allocated(self%doubleProperty(doubleProperty +i)%rank1 )) then
                      shape_=doubleProperties (i)%shape()
                      allocate(self%doubleProperty (doubleProperty +i)%rank1(shape_(1),self%doubleBufferSize))
                      deallocate(shape_)
                   end if
                   self%doubleProperty (doubleProperty +i)%rank1 (:,self%doubleBufferCount )=doubleProperties(i)
                case (-1)
                   ! Rank-1 list property
                   if (     allocated(self%doubleProperty(doubleProperty +i)%scalar     ))             deallocate(self%doubleProperty (doubleProperty +i)%scalar)
                   if (     allocated(self%doubleProperty(doubleProperty +i)%rank1      ))             deallocate(self%doubleProperty (doubleProperty +i)%rank1 )
                   if (.not.allocated(self%doubleProperty(doubleProperty +i)%rank1VarLen)) allocate(self%doubleProperty (doubleProperty +i)%rank1VarLen(self%doubleBufferSize))
                   if (associated(self%doubleProperty (doubleProperty +i)%rank1VarLen (self%doubleBufferCount )%row)) then
                      shape_=doubleProperties(i)%shape()
                      if (size(self%doubleProperty (doubleProperty +i)%rank1VarLen (self%doubleBufferCount )%row) /= shape_(1)) deallocate(self%doubleProperty (doubleProperty +i)%rank1VarLen (self%doubleBufferCount )%row)
                      deallocate(shape_)
                   end if
                   if (.not.associated(self%doubleProperty (doubleProperty +i)%rank1VarLen (self%doubleBufferCount )%row)) then
                      shape_=doubleProperties(i)%shape()
                      allocate(self%doubleProperty (doubleProperty +i)%rank1VarLen (self%doubleBufferCount )%row(shape_(1)))
                      deallocate(shape_)
                   end if
                   self%doubleProperty (doubleProperty +i)%rank1VarLen (self%doubleBufferCount )%row=doubleProperties(i)
                case (-2)
                   ! Rank-2 list property
                   if (     allocated(self%doubleProperty(doubleProperty +i)%scalar     ))             deallocate(self%doubleProperty (doubleProperty +i)%scalar     )
                   if (     allocated(self%doubleProperty(doubleProperty +i)%rank1      ))             deallocate(self%doubleProperty (doubleProperty +i)%rank1      )
                   if (     allocated(self%doubleProperty(doubleProperty +i)%rank1VarLen))             deallocate(self%doubleProperty (doubleProperty +i)%rank1VarLen)
                   if (.not.allocated(self%doubleProperty(doubleProperty +i)%rank2VarLen)) allocate(self%doubleProperty (doubleProperty +i)%rank2VarLen(self%doubleBufferSize))
                   if (associated(self%doubleProperty (doubleProperty +i)%rank2VarLen (self%doubleBufferCount )%row)) then
                      shape_=doubleProperties(i)%shape()
                      if     (                                                                                                            &
                           &   size(self%doubleProperty (doubleProperty +i)%rank2VarLen (self%doubleBufferCount )%row,dim=1) /= shape_(1) &
                           &  .or.                                                                                                        &
                           &   size(self%doubleProperty (doubleProperty +i)%rank2VarLen (self%doubleBufferCount )%row,dim=2) /= shape_(2) &
                           & ) deallocate(self%doubleProperty (doubleProperty +i)%rank2VarLen (self%doubleBufferCount )%row)
                      deallocate(shape_)
                   end if
                   if (.not.associated(self%doubleProperty (doubleProperty +i)%rank2VarLen (self%doubleBufferCount )%row)) then
                      shape_=doubleProperties(i)%shape()
                      allocate(self%doubleProperty (doubleProperty +i)%rank2VarLen (self%doubleBufferCount )%row(shape_(1),shape_(2)))
                      deallocate(shape_)
                   end if
                   self%doubleProperty (doubleProperty +i)%rank2VarLen (self%doubleBufferCount )%row=doubleProperties(i)
                case default
                   call Error_Report('unsupported rank for output property'//{introspection:location})
                end select
             end do
             deallocate(doubleProperties)
             doubleProperty                                                            =+doubleProperty                                                     &
                  &                                                                     +extractor_     %elementCount(elementTypeDouble,     time         )
             integerProperties =extractor_%extractInteger (node,time,instance,integerRanks)
             do i=1,extractor_%elementCount(elementTypeInteger ,time)
                select case (integerRanks(i))
                case (0)
                   ! Scalar property.
                   if (     allocated(self%integerProperty(integerProperty+i)%rank1 ))            deallocate(self%integerProperty(integerProperty+i)%rank1 )
                   if (.not.allocated(self%integerProperty(integerProperty+i)%scalar)) then
                      allocate(self%integerProperty(integerProperty+i)%scalar(          self%integerBufferSize))
                   end if
                   self%integerProperty (integerProperty+i)%scalar(self%integerBufferCount  )=integerProperties(i)
                case (1)
                   ! Rank-1 array property.
                   if (     allocated(self%integerProperty(integerProperty+i)%scalar))             deallocate(self%integerProperty(integerProperty+i)%scalar)
                   if (     allocated(self%integerProperty(integerProperty+i)%rank1 )) then
                      shape_=integerProperties(i)%shape()
                      if (size(self%integerProperty(integerProperty+i)%rank1,dim=1) /= shape_(1)) deallocate(self%integerProperty(integerProperty+i)%rank1)
                      deallocate(shape_)
                   end if
                   if (.not.allocated(self%integerProperty(integerProperty+i)%rank1)) then
                      shape_=integerProperties(i)%shape()
                      allocate(self%integerProperty(integerProperty+i)%rank1(shape_(1),self%integerBufferSize))
                      deallocate(shape_)
                   end if
                   self%integerProperty (integerProperty+i)%rank1 (:,self%integerBufferCount)=integerProperties(i)
                case (-1)
                   ! Rank-1 list property
                   if (     allocated(self%integerProperty(integerProperty+i)%scalar     ))             deallocate(self%integerProperty (integerProperty+i)%scalar)
                   if (     allocated(self%integerProperty(integerProperty+i)%rank1      ))             deallocate(self%integerProperty (integerProperty+i)%rank1 )
                   if (.not.allocated(self%integerProperty(integerProperty+i)%rank1VarLen)) allocate(self%integerProperty (integerProperty+i)%rank1VarLen(self%integerBufferSize))
                   if (associated(self%integerProperty (integerProperty+i)%rank1VarLen (self%integerBufferCount )%row)) then
                      shape_=integerProperties(i)%shape()
                      if (size(self%integerProperty (integerProperty+i)%rank1VarLen(self%integerBufferCount )%row) /= shape_(1)) deallocate(self%integerProperty(integerProperty+i)%rank1VarLen(self%integerBufferCount)%row)
                      deallocate(shape_)
                   end if
                   if (.not.associated(self%integerProperty (integerProperty+i)%rank1VarLen(self%integerBufferCount)%row)) then
                      shape_=integerProperties(i)%shape()
                      allocate(self%integerProperty (integerProperty+i)%rank1VarLen(self%integerBufferCount)%row(shape_(1)))
                      deallocate(shape_)
                   end if
                   self%integerProperty (integerProperty+i)%rank1VarLen (self%integerBufferCount )%row=integerProperties(i)
                case default
                   call Error_Report('unsupported rank for output property'//{introspection:location})
                end select
             end do
             deallocate(integerProperties)
             integerProperty                                                           =+integerProperty                                                    &
                  &                                                                    +extractor_     %elementCount(elementTypeInteger,     time         )
             class default
             call Error_Report('unsupported property extractor class'//{introspection:location})
          end select
          ! If buffer is full, extend it.
          if (self%integerBufferCount == self%integerBufferSize) call self%extendIntegerBuffer()
          if (self% doubleBufferCount == self% doubleBufferSize) call self%extendDoubleBuffer ()
       end do
    end if
    ! Perform any extra output tasks.
    !![
    <eventHook name="mergerTreeExtraOutput">
      <callWith>node,indexOutput,node%hostTree,nodePassesFilter,treeLock</callWith>
    </eventHook>  
    !!]
    return
  end subroutine standardOutput

  subroutine standardMakeGroup(self,tree,indexOutput)
    !!{
    Make an group in the \glc\ file in which to store {\normalfont \ttfamily tree}.
    !!}
    use            :: Galacticus_Nodes                , only : mergerTree
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Numerical_Constants_Astronomical, only : megaParsec
    use            :: String_Handling                 , only : operator(//)
    implicit none
    class  (mergerTreeOutputterStandard), intent(inout) :: self
    type   (mergerTree                 ), intent(inout) :: tree
    integer(c_size_t                   ), intent(in   ) :: indexOutput
    type   (varying_string             )                :: description, groupName

    ! Create a name for the group.
    groupName='mergerTree'
    groupName=groupName//tree%index
    ! Create a comment for the group.
    description='Data for nodes within merger tree ID='
    description=description//tree%index
    ! Create a group for the tree.
    tree%hdf5Group=self%outputGroups(indexOutput)%hdf5Group%openGroup(char(groupName),char(description))
    ! Add the merger tree weight to the group.
    call tree%hdf5Group%writeAttribute(tree%volumeWeight  ,"volumeWeight"         )
    call tree%hdf5Group%writeAttribute(1.0d0/megaParsec**3,"volumeWeightUnitsInSI")
    return
  end subroutine standardMakeGroup

  subroutine standardDumpIntegerBuffer(self,indexOutput)
    !!{
    Dump the contents of the integer properties buffer to the \glc\ output file.
    !!}
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class  (mergerTreeOutputterStandard), intent(inout) :: self
    integer(c_size_t                   ), intent(in   ) :: indexOutput
    integer                                             :: iProperty  , iMetaDatum
    type   (hdf5Object                 )                :: dataset    

    ! Write integer data from the buffer.
    if (self%integerPropertyCount > 0) then
       !$ call hdf5Access%set()
       do iProperty=1,self%integerPropertyCount
          if      (allocated(self%integerProperty(iProperty)%scalar     )) then
          call self%outputGroups(indexOutput)%nodeDataGroup%writeDataset(                                                                             &
               &                                                                  self%integerProperty(iProperty)%scalar (1:self%integerBufferCount), &
               &                                                                  self%integerProperty(iProperty)%name                              , &
               &                                                                  self%integerProperty(iProperty)%comment                           , &
               &                                                         appendTo=.true.                                                              &
               &                                                        )
          else if (allocated(self%integerProperty(iProperty)%rank1      )) then
             call self%outputGroups(indexOutput)%nodeDataGroup%writeDataset(                                                                                          &
                  &                                                                         self%integerProperty(iProperty)%rank1      (:,1:self%integerBufferCount), &
                  &                                                                         self%integerProperty(iProperty)%name                                    , &
                  &                                                                         self%integerProperty(iProperty)%comment                                 , &
                  &                                                         appendTo       =.true.                                                                  , &
                  &                                                         appendDimension=2                                                                         &
                  &                                                        )
          else if (allocated(self%integerProperty(iProperty)%rank1VarLen)) then
             call self%outputGroups(indexOutput)%nodeDataGroup%writeDataset(                                                                                          &
                  &                                                                         self%integerProperty(iProperty)%rank1VarLen(  1:self%integerBufferCount), &
                  &                                                                         self%integerProperty(iProperty)%name                                    , &
                  &                                                                         self%integerProperty(iProperty)%comment                                 , &
                  &                                                         appendTo       =.true.                                                                    &
                  &                                                        )
          end if
          if (.not.self%outputGroups(indexOutput)% integerAttributesWritten.and. self%integerProperty(iProperty)%unitsInSI /= 0.0d0) then
             dataset=self%outputGroups(indexOutput)%nodeDataGroup%openDataset(self%integerProperty(iProperty)%name)
             call dataset%writeAttribute(self%integerProperty(iProperty)%unitsInSI,"unitsInSI")
             do iMetaDatum=1,self%integerProperty(iProperty)%metaDataRank0%size()
                call     dataset%writeAttribute(self%integerProperty(iProperty)%metaDataRank0%value(iMetaDatum),char(self%integerProperty(iProperty)%metaDataRank0%key(iMetaDatum)))
             end do
             do iMetaDatum=1,self%integerProperty(iProperty)%metaDataRank1%size()
                call     dataset%writeAttribute(self%integerProperty(iProperty)%metaDataRank1%value(iMetaDatum),char(self%integerProperty(iProperty)%metaDataRank1%key(iMetaDatum)))
             end do
          end if
       end do
       self%integerPropertiesWritten=self%integerPropertiesWritten+self%integerBufferCount
       self%integerBufferCount=0
       self%outputGroups(indexOutput)%integerAttributesWritten=.true.
       !$ call hdf5Access%unset()
    end if
    return
  end subroutine standardDumpIntegerBuffer

  subroutine standardDumpDoubleBuffer(self,indexOutput)
    !!{
    Dump the contents of the double precision properties buffer to the \glc\ output file.
    !!}
    use            :: HDF5_Access  , only : hdf5Access
    use            :: IO_HDF5      , only : hdf5Object
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class  (mergerTreeOutputterStandard), intent(inout) :: self
    integer(c_size_t                   ), intent(in   ) :: indexOutput
    integer                                             :: iProperty  , iMetaDatum
    type   (hdf5Object                 )                :: dataset
    
    ! Write double data from the buffer.
    if (self%doublePropertyCount > 0) then
       !$ call hdf5Access%set()
       do iProperty=1,self%doublePropertyCount
          if      (allocated(self%doubleProperty(iProperty)%scalar     )) then
             call self%outputGroups(indexOutput)%nodeDataGroup%writeDataset(                                                                                        &
                  &                                                                         self%doubleProperty(iProperty)%scalar     (  1:self%doubleBufferCount), &
                  &                                                                         self%doubleProperty(iProperty)%name                                   , &
                  &                                                                         self%doubleProperty(iProperty)%comment                                , &
                  &                                                         appendTo       =.true.                                                                  &
                  &                                                        )
          else if (allocated(self%doubleProperty(iProperty)%rank1      )) then
             call self%outputGroups(indexOutput)%nodeDataGroup%writeDataset(                                                                                        &
                  &                                                                         self%doubleProperty(iProperty)%rank1      (:,1:self%doubleBufferCount), &
                  &                                                                         self%doubleProperty(iProperty)%name                                   , &
                  &                                                                         self%doubleProperty(iProperty)%comment                                , &
                  &                                                         appendTo       =.true.                                                                , &
                  &                                                         appendDimension=2                                                                       &
                  &                                                        )
          else if (allocated(self%doubleProperty(iProperty)%rank1VarLen)) then
             call self%outputGroups(indexOutput)%nodeDataGroup%writeDataset(                                                                                        &
                  &                                                                         self%doubleProperty(iProperty)%rank1VarLen(  1:self%doubleBufferCount), &
                  &                                                                         self%doubleProperty(iProperty)%name                                   , &
                  &                                                                         self%doubleProperty(iProperty)%comment                                , &
                  &                                                         appendTo       =.true.                                                                  &
                  &                                                        )
          else if (allocated(self%doubleProperty(iProperty)%rank2VarLen)) then
             call self%outputGroups(indexOutput)%nodeDataGroup%writeDataset(                                                                                        &
                  &                                                                         self%doubleProperty(iProperty)%rank2VarLen(  1:self%doubleBufferCount), &
                  &                                                                         self%doubleProperty(iProperty)%name                                   , &
                  &                                                                         self%doubleProperty(iProperty)%comment                                , &
                  &                                                         appendTo       =.true.                                                                  &
                  &                                                        )
          end if
          if (.not.self%outputGroups(indexOutput)% doubleAttributesWritten.and. self%doubleProperty(iProperty)%unitsInSI /= 0.0d0) then
             dataset=self%outputGroups(indexOutput)%nodeDataGroup%openDataset(self%doubleProperty(iProperty)%name)
             call        dataset%writeAttribute(self%doubleProperty(iProperty)%unitsInSI       ,"unitsInSI"        )
             do iMetaDatum=1,self%doubleProperty(iProperty)%metaDataRank0%size()
                call     dataset%writeAttribute(self%doubleProperty(iProperty)%metaDataRank0%value(iMetaDatum),char(self%doubleProperty(iProperty)%metaDataRank0%key(iMetaDatum)))
             end do
             do iMetaDatum=1,self%doubleProperty(iProperty)%metaDataRank1%size()
                call     dataset%writeAttribute(self%doubleProperty(iProperty)%metaDataRank1%value(iMetaDatum),char(self%doubleProperty(iProperty)%metaDataRank1%key(iMetaDatum)))
             end do
             if (allocated(self%doubleProperty(iProperty)%rank1Descriptors     ) .and. size(self%doubleProperty(iProperty)%rank1Descriptors     ) > 0) then
                call self%outputGroups(indexOutput)%nodeDataGroup%writeDataset(                                                                                       &
                     &                                                              self%doubleProperty(iProperty)%rank1Descriptors                                 , &
                     &                                                         trim(self%doubleProperty(iProperty)%name                  )//"Columns"               , &
                     &                                                         trim(self%doubleProperty(iProperty)%comment               )//" (column descriptions)"  &
                     &                                                        )
             end if
             if (allocated(self%doubleProperty(iProperty)%rank1DescriptorValues) .and. size(self%doubleProperty(iProperty)%rank1DescriptorValues) > 0) then
                call self%outputGroups(indexOutput)%nodeDataGroup%writeDataset(                                                                                       &
                     &                                                              self%doubleProperty(iProperty)%rank1DescriptorValues                            , &
                     &                                                         trim(self%doubleProperty(iProperty)%name                  )//"ColumnValues"          , &
                     &                                                         char(self%doubleProperty(iProperty)%rank1DescriptorComment)                            &
                     &                                                        )
                dataset=self%outputGroups(indexOutput)%nodeDataGroup%openDataset(trim(self%doubleProperty(iProperty)%name)//"ColumnValues")
                call dataset%writeAttribute(self%doubleProperty(iProperty)%rank1DescriptorUnitsInSI,"unitsInSI")
             end if
          end if
       end do
       self%doublePropertiesWritten=self%doublePropertiesWritten+self%doubleBufferCount
       self%doubleBufferCount=0
       self%outputGroups(indexOutput)%doubleAttributesWritten=.true.
       !$ call hdf5Access%unset()
    end if
    return
  end subroutine standardDumpDoubleBuffer

  subroutine standardExtendIntegerBuffer(self)
    !!{
    Extend the size of the integer buffer.
    !!}
    use :: IO_HDF5, only : hdf5VarInteger8
    implicit none
    class  (mergerTreeOutputterStandard), intent(inout)                 :: self
    integer(kind_int8                  ), allocatable  , dimension(:  ) :: scalarTemporary
    integer(kind_int8                  ), allocatable  , dimension(:,:) :: rank1Temporary
    type   (hdf5VarInteger8            ), allocatable  , dimension(:  ) :: rank1VarLenTemporary
    integer                                                             :: i                   , j

    do i=1,size(self%integerProperty)
       if (allocated(self%integerProperty(i)%scalar)) then
          call move_alloc(self%integerProperty(i)%scalar,scalarTemporary)
          allocate(self%integerProperty(i)%scalar(                           self%integerBufferSize+standardBufferSizeIncrement))
          self%integerProperty(i)%scalar(  1:self%integerBufferSize)=scalarTemporary
          deallocate(scalarTemporary)
       end if
       if (allocated(self%integerProperty(i)%rank1)) then
          call move_alloc(self%integerProperty(i)%rank1 ,rank1Temporary )
          allocate(self%integerProperty(i)%rank1 (size(rank1Temporary,dim=1),self%integerBufferSize+standardBufferSizeIncrement))
          self%integerProperty(i)%rank1(:,1:self%integerBufferSize)=rank1Temporary
          deallocate(rank1Temporary )
       end if
       if (allocated(self%integerProperty(i)%rank1VarLen)) then
          call move_alloc(self%integerProperty(i)%rank1VarLen,rank1VarLenTemporary)
          allocate(self%integerProperty(i)%rank1VarLen(                           self%integerBufferSize+standardBufferSizeIncrement))
          self%integerProperty(i)%rank1VarLen(  1:self%integerBufferSize)=rank1VarLenTemporary
          ! Nullify row pointers in the temporary array to avoid the associated data being destroyed when our temporary array is
          ! deallocated.
          do j=1,size(rank1VarLenTemporary)
             rank1VarLenTemporary(j)%row => null()
          end do
          deallocate(rank1VarLenTemporary)
       end if
    end do
    self%integerBufferSize=self%integerBufferSize+standardBufferSizeIncrement
    return
  end subroutine standardExtendIntegerBuffer

  subroutine standardExtendDoubleBuffer(self)
    !!{
    Extend the size of the double buffer.
    !!}
    use :: IO_HDF5, only : hdf5VarDouble , hdf5VarDouble2D
    implicit none
    class           (mergerTreeOutputterStandard), intent(inout)                 :: self
    double precision                             , allocatable  , dimension(:  ) :: scalarTemporary
    double precision                             , allocatable  , dimension(:,:) :: rank1Temporary
    type            (hdf5VarDouble              ), allocatable  , dimension(:  ) :: rank1VarLenTemporary
    type            (hdf5VarDouble2D            ), allocatable  , dimension(:  ) :: rank2VarLenTemporary
    integer                                                                      :: i                   , j

    do i=1,size(self%doubleProperty)
       if (allocated(self%doubleProperty(i)%scalar    )) then
          call move_alloc(self%doubleProperty(i)%scalar    ,scalarTemporary)
          allocate(self%doubleProperty(i)%scalar    (                           self%doubleBufferSize+standardBufferSizeIncrement))
          self%doubleProperty(i)%scalar    (  1:self%doubleBufferSize)=scalarTemporary
          deallocate(scalarTemporary)
       end if
       if (allocated(self%doubleProperty(i)%rank1     )) then
          call move_alloc(self%doubleProperty(i)%rank1     ,rank1Temporary )
          allocate(self%doubleProperty(i)%rank1     (size(rank1Temporary,dim=1),self%doubleBufferSize+standardBufferSizeIncrement))
          self%doubleProperty(i)%rank1     (:,1:self%doubleBufferSize)=rank1Temporary
          deallocate(rank1Temporary )
       end if
       if (allocated(self%doubleProperty(i)%rank1VarLen)) then
          call move_alloc(self%doubleProperty(i)%rank1VarLen,rank1VarLenTemporary)
          allocate(self%doubleProperty(i)%rank1VarLen(                           self%doubleBufferSize+standardBufferSizeIncrement))
          self%doubleProperty(i)%rank1VarLen(  1:self%doubleBufferSize)=rank1VarLenTemporary
          ! Nullify row pointers in the temporary array to avoid the associated data being destroyed when our temporary array is
          ! deallocated.
          do j=1,size(rank1VarLenTemporary)
             rank1VarLenTemporary(j)%row => null()
          end do
          deallocate(rank1VarLenTemporary)
       end if
       if (allocated(self%doubleProperty(i)%rank2VarLen)) then
          call move_alloc(self%doubleProperty(i)%rank2VarLen,rank2VarLenTemporary)
          allocate(self%doubleProperty(i)%rank2VarLen(                           self%doubleBufferSize+standardBufferSizeIncrement))
          self%doubleProperty(i)%rank2VarLen(  1:self%doubleBufferSize)=rank2VarLenTemporary
          ! Nullify row pointers in the temporary array to avoid the associated data being destroyed when our temporary array is
          ! deallocated.
          do j=1,size(rank2VarLenTemporary)
             rank2VarLenTemporary(j)%row => null()
          end do
          deallocate(rank2VarLenTemporary)
       end if
    end do
    self%doubleBufferSize=self%doubleBufferSize+standardBufferSizeIncrement
    return
  end subroutine standardExtendDoubleBuffer

  subroutine standardPropertiesCount(self,time,node)
    !!{
    Count up the number of properties that will be output.
    !!}
    use :: Error                   , only : Error_Report
    use :: Galacticus_Nodes        , only : treeNode
    use :: Node_Property_Extractors, only : elementTypeDouble         , elementTypeInteger       , nodePropertyExtractorIntegerScalar, nodePropertyExtractorIntegerTuple, &
         &                                  nodePropertyExtractorMulti, nodePropertyExtractorNull, nodePropertyExtractorScalar       , nodePropertyExtractorTuple       , &
         &                                  nodePropertyExtractorArray, nodePropertyExtractorList, nodePropertyExtractorList2D       , nodePropertyExtractorIntegerList
    implicit none
    class           (mergerTreeOutputterStandard), intent(inout) :: self
    double precision                             , intent(in   ) :: time
    type            (treeNode                   ), intent(inout) :: node

    self%integerPropertyCount=0
    self%doublePropertyCount =0
    call node%outputCount(self%integerPropertyCount,self%doublePropertyCount,time)
    self%integerScalarCount=self%integerPropertyCount
    self% doubleScalarCount=self% doublePropertyCount
    ! Handle extracted properties.
    select type (extractor_ => self%nodePropertyExtractor_)
    type is (nodePropertyExtractorNull)
       ! Null extractor - simply ignore.
    class is (nodePropertyExtractorScalar       )
       ! Scalar property extractor - simply increment the double property output count by one.
       self%doublePropertyCount =self%doublePropertyCount +1
    class is (nodePropertyExtractorTuple        )
       ! Tuple property extractor - increment the double property output count by the number of elements.
       self%doublePropertyCount =self%doublePropertyCount +extractor_%elementCount(time)
    class is (nodePropertyExtractorArray        )
       ! Array property extractor - increment the double property output count by the number of elements.
       self%doublePropertyCount =self%doublePropertyCount +extractor_%elementCount(time)
    class is (nodePropertyExtractorList         )
       ! List property extractor - simply increment the double property output count by one.
       self%doublePropertyCount =self%doublePropertyCount +1
    class is (nodePropertyExtractorList2D       )
       ! 2D list property extractor - simply increment the double property output count by one.
       self%doublePropertyCount =self%doublePropertyCount +1
    class is (nodePropertyExtractorIntegerScalar)
       ! Integer scalar property extractor - simply increment the integer property output count by one.
       self%integerPropertyCount=self%integerPropertyCount+1
    class is (nodePropertyExtractorIntegerTuple )
       ! Integer tuple property extractor - increment the integer property output count by the number of elements.
       self%integerPropertyCount=self%integerPropertyCount+extractor_%elementCount(time)
    class is (nodePropertyExtractorIntegerList  )
       ! Integer list property extractor - simply increment the integer property output count by one.
       self%integerPropertyCount=self%integerPropertyCount+1
    class is (nodePropertyExtractorMulti        )
       ! Multi property extractor - increment double and integer property output counts.
       self%integerPropertyCount=self%integerPropertyCount+extractor_%elementCount(elementTypeInteger,time)
       self% doublePropertyCount=self% doublePropertyCount+extractor_%elementCount(elementTypeDouble ,time)
    class default
       call Error_Report('unsupported property extractor class'//{introspection:location})
    end select
    return
  end subroutine standardPropertiesCount

  subroutine standardBuffersAllocate(self,indexOutput)
    !!{
    Allocate buffers for storage of properties.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class  (mergerTreeOutputterStandard), intent(inout) :: self
    integer(c_size_t                   ), intent(in   ) :: indexOutput

    if (self%integerPropertyCount > 0 .and. (.not.allocated(self%integerProperty) .or. self%integerPropertyCount > size(self%integerProperty))) then
       if (allocated(self%integerProperty)) deallocate(self%integerProperty)
       allocate(self%integerProperty(self%integerPropertyCount))
    end if
    if (self%doublePropertyCount  > 0 .and. (.not.allocated(self%doubleProperty ) .or. self%doublePropertyCount  > size(self%doubleProperty ))) then
       if (allocated(self%doubleProperty )) deallocate(self%doubleProperty )
       allocate(self%doubleProperty (self%doublePropertyCount ))
    end if
    return
  end subroutine standardBuffersAllocate

  subroutine standardPropertyNamesEstablish(self,time,node)
    !!{
    Set names for the properties.
    !!}
    use :: Display                 , only : displayGreen              , displayReset
    use :: Error                   , only : Error_Report
    use :: Galacticus_Nodes        , only : treeNode
    use :: Hashes                  , only : doubleHash                , rank1DoubleHash
    use :: Node_Property_Extractors, only : elementTypeDouble         , elementTypeInteger       , nodePropertyExtractorIntegerScalar, nodePropertyExtractorIntegerTuple, &
         &                                  nodePropertyExtractorMulti, nodePropertyExtractorNull, nodePropertyExtractorScalar       , nodePropertyExtractorTuple       , &
         &                                  nodePropertyExtractorArray, nodePropertyExtractorList, nodePropertyExtractorList2D       , nodePropertyExtractorIntegerList
    implicit none
    class           (mergerTreeOutputterStandard), intent(inout)               :: self
    double precision                             , intent(in   )               :: time
    type            (treeNode                   ), intent(inout)               :: node
    integer                                                                    :: doubleProperty, integerProperty, &
         &                                                                        i
    type            (varying_string             ), allocatable  , dimension(:) :: namesTmp      , descriptionsTmp

    if (allocated(self%integerProperty)) then
       do i=1,size(self%integerProperty)
          if (allocated(self%integerProperty(i)%metaDataRank0)) deallocate(self%integerProperty(i)%metaDataRank0)
          allocate(self%integerProperty(i)%metaDataRank0)
          self%integerProperty(i)%metaDataRank0=doubleHash()
       end do
    end if
    if (allocated(self%doubleProperty )) then
       do i=1,size(self%doubleProperty )
          if (allocated(self%doubleProperty (i)%metaDataRank0)) deallocate(self%doubleProperty (i)%metaDataRank0)
          if (allocated(self%doubleProperty (i)%metaDataRank1)) deallocate(self%doubleProperty (i)%metaDataRank1)
          allocate(self%doubleProperty (i)%metaDataRank0)
          allocate(self%doubleProperty (i)%metaDataRank1)
          self%doubleProperty (i)%metaDataRank0=     doubleHash()
          self%doubleProperty (i)%metaDataRank1=rank1DoubleHash()
       end do
    end if
    integerProperty=0
    doubleProperty =0
    call node%outputNames(integerProperty,self%integerProperty,doubleProperty,self%doubleProperty,time)
    ! Handle extracted properties.
    select type (extractor_ => self%nodePropertyExtractor_)
    type is (nodePropertyExtractorNull)
       ! Null extractor - simply ignore.
    class is (nodePropertyExtractorScalar       )
       ! Scalar property extractor - get the name, description, and units.
       self%doubleProperty (doubleProperty +1                                                                 )%name      =extractor_%name        (                       )
       self%doubleProperty (doubleProperty +1                                                                 )%comment   =extractor_%description (                       )
       self%doubleProperty (doubleProperty +1                                                                 )%unitsInSI =extractor_%unitsInSI   (                       )
       call    extractor_%metaData(node      ,self%doubleProperty (doubleProperty +1)%metaDataRank0,self%doubleProperty (doubleProperty +1)%metaDataRank1)
       doubleProperty =doubleProperty +1
    class is (nodePropertyExtractorTuple        )
       ! Tuple property extractor - get the names, descriptions, and units.
       call extractor_%names       (time,namesTmp       )
       call extractor_%descriptions(time,descriptionsTmp)
       self%doubleProperty (doubleProperty +1:doubleProperty +extractor_%elementCount(                   time))%name      =namesTmp
       self%doubleProperty (doubleProperty +1:doubleProperty +extractor_%elementCount(                   time))%comment   =descriptionsTmp
       self%doubleProperty (doubleProperty +1:doubleProperty +extractor_%elementCount(                   time))%unitsInSI =extractor_%unitsInSI   (                   time)
       do i=1,extractor_%elementCount(                   time)
          call extractor_%metaData(node,i     ,self%doubleProperty (doubleProperty +i)%metaDataRank0,self%doubleProperty (doubleProperty +i)%metaDataRank1)
       end do
       doubleProperty =doubleProperty +extractor_%elementCount(time)
       deallocate(namesTmp       )
       deallocate(descriptionsTmp)
    class is (nodePropertyExtractorArray        )
       ! Array property extractor - get the names, descriptions, and units.
       call extractor_%names       (namesTmp       ,time)
       call extractor_%descriptions(descriptionsTmp,time)
       self%doubleProperty (doubleProperty +1:doubleProperty +extractor_%elementCount(                   time))%name      =namesTmp
       self%doubleProperty (doubleProperty +1:doubleProperty +extractor_%elementCount(                   time))%comment   =descriptionsTmp
       self%doubleProperty (doubleProperty +1:doubleProperty +extractor_%elementCount(                   time))%unitsInSI =extractor_%unitsInSI   (                   time)
       do i=1,extractor_%elementCount(time)
          if (allocated(self%doubleProperty(doubleProperty+i)%rank1Descriptors     )) deallocate(self%doubleProperty(doubleProperty+i)%rank1Descriptors     )
          if (allocated(self%doubleProperty(doubleProperty+i)%rank1DescriptorValues)) deallocate(self%doubleProperty(doubleProperty+i)%rank1DescriptorValues)
          call extractor_%columnDescriptions(self%doubleProperty(doubleProperty+i)%rank1Descriptors,self%doubleProperty(doubleProperty+i)%rank1DescriptorValues,self%doubleProperty(doubleProperty+i)%rank1DescriptorComment,self%doubleProperty(doubleProperty+i)%rank1DescriptorUnitsInSI,time)
       end do
       do i=1,extractor_%elementCount(                   time)
          call extractor_%metaData(node,i     ,self%doubleProperty (doubleProperty +i)%metaDataRank0,self%doubleProperty (doubleProperty +1)%metaDataRank1)
       end do
       doubleProperty =doubleProperty +extractor_%elementCount(time)
       deallocate(namesTmp       )
       deallocate(descriptionsTmp)
    class is (nodePropertyExtractorList         )
       ! List property extractor - get the name, description, and units.
       call extractor_%names       (namesTmp            )
       call extractor_%descriptions(descriptionsTmp     )
       self%doubleProperty (doubleProperty +1:doubleProperty +extractor_%elementCount(                       ))%name      =namesTmp
       self%doubleProperty (doubleProperty +1:doubleProperty +extractor_%elementCount(                       ))%comment   =descriptionsTmp
       self%doubleProperty (doubleProperty +1:doubleProperty +extractor_%elementCount(                       ))%unitsInSI =extractor_%unitsInSI  (                       )
       call    extractor_%metaData(node            ,self%doubleProperty (doubleProperty +1)%metaDataRank0,self%doubleProperty (doubleProperty +1)%metaDataRank1)
       doubleProperty =doubleProperty +1
    class is (nodePropertyExtractorList2D       )
       ! 2D list property extractor - get the name, description, and units.
       call extractor_%names       (namesTmp            )
       call extractor_%descriptions(descriptionsTmp     )
       self%doubleProperty (doubleProperty +1:doubleProperty +extractor_%elementCount(                       ))%name      =namesTmp
       self%doubleProperty (doubleProperty +1:doubleProperty +extractor_%elementCount(                       ))%comment   =descriptionsTmp
       self%doubleProperty (doubleProperty +1:doubleProperty +extractor_%elementCount(                       ))%unitsInSI =extractor_%unitsInSI  (                       )
       call    extractor_%metaData(node  ,time,self%doubleProperty (doubleProperty +1)%metaDataRank0,self%doubleProperty (doubleProperty +1)%metaDataRank1)
       doubleProperty =doubleProperty +1
    class is (nodePropertyExtractorIntegerScalar)
       ! Integer scalar property extractor - get the name, description, and units.
       self%integerProperty(integerProperty+1                                                                 )%name     =extractor_%name        (                       )
       self%integerProperty(integerProperty+1                                                                 )%comment  =extractor_%description (                       )
       self%integerProperty(integerProperty+1                                                                 )%unitsInSI=extractor_%unitsInSI   (                       )
       call    extractor_%metaData(node       ,self%integerProperty(integerProperty+1)%metaDataRank0,self%integerProperty(integerProperty+1)%metaDataRank1)
       integerProperty=integerProperty+1
    class is (nodePropertyExtractorIntegerTuple )
       ! Integer tuple property extractor - get the names, descriptions, and units.
       call extractor_%names       (time,namesTmp       )
       call extractor_%descriptions(time,descriptionsTmp)
       self%integerProperty(integerProperty+1:integerProperty+extractor_%elementCount(                   time))%name     =namesTmp
       self%integerProperty(integerProperty+1:integerProperty+extractor_%elementCount(                   time))%comment  =descriptionsTmp
       self%integerProperty(integerProperty+1:integerProperty+extractor_%elementCount(                   time))%unitsInSI=extractor_%unitsInSI   (                   time)
       do i=1,extractor_%elementCount(                   time)
          call extractor_%metaData(node,i     ,self%integerProperty(integerProperty+i)%metaDataRank0,self%integerProperty(integerProperty+1)%metaDataRank1)
       end do
       integerProperty=integerProperty+extractor_%elementCount(time)
       deallocate(namesTmp       )
       deallocate(descriptionsTmp)
    class is (nodePropertyExtractorIntegerList  )
       ! Integer list property extractor - get the name, description, and units.
       call extractor_%names       (namesTmp            )
       call extractor_%descriptions(descriptionsTmp     )
       self%integerProperty(integerProperty+1:integerProperty+extractor_%elementCount(                       ))%name      =namesTmp
       self%integerProperty(integerProperty+1:integerProperty+extractor_%elementCount(                       ))%comment   =descriptionsTmp
       self%integerProperty(integerProperty+1:integerProperty+extractor_%elementCount(                       ))%unitsInSI =extractor_%unitsInSI  (                       )
       call    extractor_%metaData(node            ,self%integerProperty(integerProperty+1)%metaDataRank0,self%integerProperty(integerProperty+1)%metaDataRank1)
       integerProperty =integerProperty +1
    class is (nodePropertyExtractorMulti        )
       ! Multi property extractor - get the names, descriptions, and units.
       if (extractor_%elementCount(elementTypeDouble ,time) > 0) then
          call extractor_%names       (elementTypeDouble ,time,namesTmp       )
          call extractor_%descriptions(elementTypeDouble ,time,descriptionsTmp)
          self%doubleProperty(doubleProperty +1:doubleProperty +extractor_%elementCount(elementTypeDouble ,time))%name      =namesTmp
          self%doubleProperty(doubleProperty +1:doubleProperty +extractor_%elementCount(elementTypeDouble ,time))%comment   =descriptionsTmp
          self%doubleProperty(doubleProperty +1:doubleProperty +extractor_%elementCount(elementTypeDouble ,time))%unitsInSI =extractor_%unitsInSI   (elementTypeDouble ,time)
          do i=1,extractor_%elementCount(elementTypeDouble,time)
             call extractor_%metaData(node,elementTypeDouble,time,i,self%doubleProperty (doubleProperty +i)%metaDataRank0,self%doubleProperty (doubleProperty +i)%metaDataRank1)
          end do
          do i=1,extractor_%elementCount(elementTypeDouble,time)
             if (allocated(self%doubleProperty(doubleProperty+i)%rank1Descriptors     )) deallocate(self%doubleProperty(doubleProperty+i)%rank1Descriptors     )
             if (allocated(self%doubleProperty(doubleProperty+i)%rank1DescriptorValues)) deallocate(self%doubleProperty(doubleProperty+i)%rank1DescriptorValues)
             call extractor_%columnDescriptions(elementTypeDouble,i,time,self%doubleProperty(doubleProperty+i)%rank1Descriptors,self%doubleProperty(doubleProperty+i)%rank1DescriptorValues,self%doubleProperty(doubleProperty+i)%rank1DescriptorComment,self%doubleProperty(doubleProperty+i)%rank1DescriptorUnitsInSI)
          end do
          doubleProperty =doubleProperty +extractor_%elementCount(elementTypeDouble ,time)
          deallocate(namesTmp       )
          deallocate(descriptionsTmp)
       end if
       if (extractor_%elementCount(elementTypeInteger,time) > 0) then
          call extractor_%names       (elementTypeInteger,time,namesTmp       )
          call extractor_%descriptions(elementTypeInteger,time,descriptionsTmp)
          self%integerProperty(integerProperty+1:integerProperty+extractor_%elementCount(elementTypeInteger,time))%name     =namesTmp
          self%integerProperty(integerProperty+1:integerProperty+extractor_%elementCount(elementTypeInteger,time))%comment  =descriptionsTmp
          self%integerProperty(integerProperty+1:integerProperty+extractor_%elementCount(elementTypeInteger,time))%unitsInSI=extractor_%unitsInSI   (elementTypeInteger,time)
          do i=1,extractor_%elementCount(elementTypeInteger,time)
             call extractor_%metaData(node,elementTypeInteger,time,i,self%integerProperty(integerProperty+i)%metaDataRank0,self%integerProperty(integerProperty+i)%metaDataRank1)
          end do
          integerProperty=integerProperty+extractor_%elementCount(elementTypeInteger,time)
          deallocate(namesTmp       )
          deallocate(descriptionsTmp)
       end if
    class default
       call Error_Report('unsupported property extractor class'//{introspection:location})
    end select
    ! Verify that no names are duplicated.
    if (self%integerPropertyCount > 0) then
       do i=1,self%integerPropertyCount
          if     (                                                                                                                                       &
               &   (self%integerPropertyCount > 1 .and. any(self%integerProperty(i+1:self%integerPropertyCount)%name == self%integerProperty(i)%name))   &
               &  .or.                                                                                                                                   &
               &   (self% doublePropertyCount > 0 .and. any(self% doubleProperty(  1:self% doublePropertyCount)%name == self%integerProperty(i)%name))   &
               & ) call Error_Report(                                                                                                                    &
               &                     "duplicate property name '"//trim(self%integerProperty(i)%name)//"' detected"//char(10)                          // &
               &                     displayGreen()//'  HELP: '//displayReset()                                                                       // &
               &                     "This can happen if you have multiple copies of identical `nodePropertyExtractor` objects in your parameter file"// &
               &                     {introspection:location}                                                                                            &
               &                    )
       end do
    end if
    if (self% doublePropertyCount > 0) then
       do i=1,self% doublePropertyCount
          if     (                                                                                                                                       &
               &   (self% doublePropertyCount > 1 .and. any(self% doubleProperty(i+1:self% doublePropertyCount)%name == self% doubleProperty(i)%name))   &
               & ) call Error_Report(                                                                                                                    &
               &                     "duplicate property name '"//trim(self% doubleProperty(i)%name)//"' detected"//char(10)                          // &
               &                     displayGreen()//'  HELP: '//displayReset()                                                                       // &
               &                     "This can happen if you have multiple copies of identical `nodePropertyExtractor` objects in your parameter file"// &
               &                     {introspection:location}                                                                                            &
               &                    )
       end do
    end if
    return
  end subroutine standardPropertyNamesEstablish

  subroutine standardOutputGroupCreate(self,indexOutput,time)
    !!{
    Create a group in which to store this output.
    !!}
    use            :: Output_HDF5                     , only : outputTo     => outputGroup
    use            :: HDF5_Access                     , only : hdf5Access
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Numerical_Constants_Astronomical, only : gigaYear    , megaParsec
    use            :: String_Handling                 , only : operator(//)
    implicit none
    class           (mergerTreeOutputterStandard), intent(inout)               :: self
    integer         (c_size_t                   ), intent(in   )               :: indexOutput
    double precision                             , intent(in   )               :: time
    type            (outputGroup                ), allocatable  , dimension(:) :: outputGroupsTemporary
    type            (varying_string             )                              :: description          , groupName
    double precision                                                           :: expansionFactor      , distanceComoving
    integer                                                                    :: i
    
    !$ call hdf5Access%set()
    ! Ensure group ID space is large enough.
    if (indexOutput > self%outputGroupsCount) then
       if (allocated(self%outputGroups)) then
          call Move_Alloc(self%outputGroups,outputGroupsTemporary)
          self%outputGroupsCount=max(self%outputGroupsCount+standardOutputGroupsIncrement,(indexOutput/standardOutputGroupsIncrement+1)*standardOutputGroupsIncrement)
          allocate(self%outputGroups(self%outputGroupsCount))
          !![
	  <workaround type="gfortran" PR="46897" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=46897">
	    <seeAlso type="gfortran" PR="57696" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=57696"/>
	    <description>
	      Type-bound defined assignment not done because multiple part array references would occur in intermediate expressions.
	    </description>
	  !!]
          do i=1,size(outputGroupsTemporary)
             self%outputGroups(i)=outputGroupsTemporary(i)
          end do
          !![
	  </workaround>
          !!]
          self%outputGroups(size(outputGroupsTemporary)+1:size(self%outputGroups))%opened                  =.false.
          self%outputGroups(size(outputGroupsTemporary)+1:size(self%outputGroups))%integerAttributesWritten=.false.
          self%outputGroups(size(outputGroupsTemporary)+1:size(self%outputGroups))%doubleAttributesWritten =.false.
          deallocate(outputGroupsTemporary)
       else
          self%outputGroupsCount=max(standardOutputGroupsIncrement,(indexOutput/standardOutputGroupsIncrement+1)*standardOutputGroupsIncrement)
          allocate(self%outputGroups(self%outputGroupsCount))
          self%outputGroups%opened                  =.false.
          self%outputGroups%integerAttributesWritten=.false.
          self%outputGroups%doubleAttributesWritten =.false.
       end if
    end if
    ! Make the enclosing group if it has not been created.
    if (.not.self%outputsGroupOpened) then
       self%outputsGroup      =outputTo%openGroup(char(self%outputsGroupName),'Contains all outputs from Galacticus.')
       self%outputsGroupOpened=.true.
    end if
    !$ call hdf5Access%unset()
    ! Create the group if it has not been created.
    if (.not.self%outputGroups(indexOutput)%opened) then
       ! Create a name for the group.
       groupName='Output'
       groupName=groupName//indexOutput
       ! Create a comment for the group.
       description='Data for output number '
       description=description//indexOutput
       ! Create a group for the tree.
       !$ call hdf5Access%set()
       self%outputGroups(indexOutput)%hdf5Group    =self%outputsGroup                       %openGroup(char(groupName),char(description)                                   )
       self%outputGroups(indexOutput)%nodeDataGroup=self%outputGroups(indexOutput)%hdf5Group%openGroup("nodeData"     ,"Group containing data on all nodes at this output.")
       self%outputGroups(indexOutput)%opened                  =.true.
       self%outputGroups(indexOutput)%integerAttributesWritten=.false.
       self%outputGroups(indexOutput)%doubleAttributesWritten =.false.
       ! Add the time to this group.
       expansionFactor    =+self%cosmologyFunctions_%expansionFactor (time)
       if (expansionFactor > 1.0d0) then
          distanceComoving=-1.0d0
       else if (expansionFactor == 1.0d0) then
          distanceComoving=+0.0d0
       else
          distanceComoving=+self%cosmologyFunctions_%distanceComoving(time)
       end if
       call self%outputGroups(indexOutput)%hdf5Group%writeAttribute(time            ,'outputTime'            )
       call self%outputGroups(indexOutput)%hdf5Group%writeAttribute(gigaYear        ,'timeUnitInSI'          )
       call self%outputGroups(indexOutput)%hdf5Group%writeAttribute(megaParsec      ,'distanceUnitInSI'      )
       call self%outputGroups(indexOutput)%hdf5Group%writeAttribute(expansionFactor ,'outputExpansionFactor' )
       call self%outputGroups(indexOutput)%hdf5Group%writeAttribute(distanceComoving,'outputComovingDistance')
       !$ call hdf5Access%unset()
    end if
    return
  end subroutine standardOutputGroupCreate

  subroutine outputGroupAssign(to,from)
    implicit none
    class(outputGroup), intent(  out) :: to
    class(outputGroup), intent(in   ) :: from

    to%doubleAttributesWritten =from%doubleAttributesWritten
    to%integerAttributesWritten=from%integerAttributesWritten
    to%opened                  =from%opened
    to%hdf5Group               =from%hdf5Group
    to%nodeDataGroup           =from%nodeDataGroup
    return
  end subroutine outputGroupAssign
