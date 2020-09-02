!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implements the standard class for outputting merger trees.

  use :: Cosmology_Functions     , only : cosmologyFunctions   , cosmologyFunctionsClass
  use :: Galactic_Filters        , only : galacticFilter       , galacticFilterClass
  use :: IO_HDF5                 , only : hdf5Object
  use :: Kind_Numbers            , only : kind_int8
  use :: Node_Property_Extractors, only : nodePropertyExtractor, nodePropertyExtractorClass

  type outputGroup
     !% Type used for output group information.
     logical                                        :: doubleAttributesWritten, integerAttributesWritten, &
          &                                            opened
     type   (hdf5Object)                            :: hdf5Group              , nodeDataGroup
     type   (hdf5Object), allocatable, dimension(:) :: integerDataset         , doubleDataset
  end type outputGroup

  ! Parameters controlling the size of output data.
  !# <scoping>
  !#  <module variables="standardNameLengthMax, standardCommentLengthMax, standardBufferSizeIncrement"/>
  !# </scoping>
  integer, parameter :: standardOutputGroupsIncrement=  10, standardNameLengthMax   =256, &
       &                standardBufferSizeIncrement  =1024, standardCommentLengthMax=256

  !# <mergerTreeOutputter name="mergerTreeOutputterStandard">
  !#  <description>The standard merger tree outputter.</description>
  !#  <stateStorable>
  !#   <restoreTo variables="outputsGroupOpened" state=".false."/>
  !#   <restoreTo variables="outputGroupsCount  , doublePropertiesWritten, integerPropertiesWritten, doubleBufferCount   , integerBufferCount"                                                                             state="0"                          />
  !#   <restoreTo variables="doublePropertyCount, integerPropertyCount"                                                                                                                                                    state="-1"                         />
  !#   <restoreTo variables="integerBufferSize  , doubleBufferSize"                                                                                                                                                        state="standardBufferSizeIncrement"/>
  !#   <exclude   variables="integerBuffer      , doubleBuffer           , doublePropertyNames     , integerPropertyNames, doublePropertyComments, integerPropertyComments, doublePropertyUnitsSI, integerPropertyUnitsSI"                                    />
  !#  </stateStorable>
  !# </mergerTreeOutputter>
  type, extends(mergerTreeOutputterClass) :: mergerTreeOutputterStandard
     !% Implementation of the standard merger tree outputter.
     private
     logical                                                                     :: outputReferences
     type            (varying_string              )                              :: outputsGroupName
     type            (hdf5Object                  )                              :: outputsGroup
     logical                                                                     :: outputsGroupOpened
     integer         (c_size_t                    )                              :: outputGroupsCount
     integer                                                                     :: doublePropertyCount              , integerPropertyCount
     integer                                                                     :: doublePropertiesWritten          , integerPropertiesWritten
     integer                                                                     :: doubleBufferCount                , integerBufferCount
     integer                                                                     :: integerBufferSize                , doubleBufferSize
     integer         (kind=kind_int8              ), allocatable, dimension(:,:) :: integerBuffer
     double precision                              , allocatable, dimension(:,:) :: doubleBuffer
     character       (len=standardNameLengthMax   ), allocatable, dimension(:  ) :: doublePropertyNames              , integerPropertyNames
     character       (len=standardCommentLengthMax), allocatable, dimension(:  ) :: doublePropertyComments           , integerPropertyComments
     double precision                              , allocatable, dimension(:  ) :: doublePropertyUnitsSI            , integerPropertyUnitsSI
     type            (outputGroup                 ), allocatable, dimension(:  ) :: outputGroups
     class           (galacticFilterClass         ), pointer                     :: galacticFilter_         => null()
     class           (cosmologyFunctionsClass     ), pointer                     :: cosmologyFunctions_     => null()
     class           (nodePropertyExtractorClass  ), pointer                     :: nodePropertyExtractor_  => null()
   contains
     !@ <objectMethods>
     !@  <object>mergerTreeOutputterStandard</object>
     !@  <objectMethod>
     !@   <method>makeGroup</method>
     !@   <type>void</type>
     !@   <arguments>\textcolor{red}{\textless type(mergerTree)\textgreater} tree\arginout, \textcolor{red}{\textless integer(c\_size\_t)\textgreater} indexOutput\argin</arguments>
     !@   <description>Make an group in the \glc\ file in which to store {\normalfont \ttfamily tree}.</description>
     !@  </objectMethod>
     !@  <objectMethod>
     !@   <method>dumpIntegerBuffer</method>
     !@   <type>void</type>
     !@   <arguments>\textcolor{red}{\textless integer(c\_size\_t)\textgreater} indexOutput\argin</arguments>
     !@   <description>Dump the contents of the integer properties buffer to the \glc\ output file.</description>
     !@  </objectMethod>
     !@  <objectMethod>
     !@   <method>dumpDoubleBuffer</method>
     !@   <type>void</type>
     !@   <arguments>\textcolor{red}{\textless integer(c\_size\_t)\textgreater} indexOutput\argin</arguments>
     !@   <description>Dump the contents of the double properties buffer to the \glc\ output file.</description>
     !@  </objectMethod>
     !@  <objectMethod>
     !@   <method>extendIntegerBuffer</method>
     !@   <type>void</type>
     !@   <arguments></arguments>
     !@   <description>Extend the size of the integer buffer.</description>
     !@  </objectMethod>
     !@  <objectMethod>
     !@   <method>extendDoubleBuffer</method>
     !@   <type>void</type>
     !@   <arguments></arguments>
     !@   <description>Extend the size of the double buffer.</description>
     !@  </objectMethod>
     !@  <objectMethod>
     !@   <method>propertiesCount</method>
     !@   <type>void</type>
     !@   <arguments>\doublezero\ time\argin, \textcolor{red}{\textless type(treeNode)\textgreater} *node\arginout</arguments>
     !@   <description>Count up the number of properties that will be output.</description>
     !@  </objectMethod>
     !@  <objectMethod>
     !@   <method>buffersAllocate</method>
     !@   <type>void</type>
     !@   <arguments>\textcolor{red}{\textless integer(c\_size\_t)\textgreater} indexOutput\argin</arguments>
     !@   <description>Allocate buffers for storage of properties.</description>
     !@  </objectMethod>
     !@  <objectMethod>
     !@   <method>propertyNamesEstablish</method>
     !@   <type>void</type>
     !@   <arguments>\doublezero\ time\argin, \textcolor{red}{\textless type(treeNode)\textgreater} *node\arginout</arguments>
     !@   <description>Set names for the properties.</description>
     !@  </objectMethod>
     !@  <objectMethod>
     !@   <method>outputGroupCreate</method>
     !@   <type>void</type>
     !@   <arguments>\textcolor{red}{\textless integer(c\_size\_t)\textgreater} indexOutput\argin, \doublezero\ time\argin</arguments>
     !@   <description>Create a group in which to store this output.</description>
     !@  </objectMethod>
     !@ </objectMethods>
     final     ::                           standardDestructor
     procedure :: output                 => standardOutput
     procedure :: finalize               => standardFinalize
     procedure :: makeGroup              => standardMakeGroup
     procedure :: dumpIntegerBuffer      => standardDumpIntegerBuffer
     procedure :: extendIntegerBuffer    => standardExtendIntegerBuffer
     procedure :: extendDoubleBuffer     => standardExtendDoubleBuffer
     procedure :: dumpDoubleBuffer       => standardDumpDoubleBuffer
     procedure :: propertiesCount        => standardPropertiesCount
     procedure :: buffersAllocate        => standardBuffersAllocate
     procedure :: propertyNamesEstablish => standardPropertyNamesEstablish
     procedure :: outputGroupCreate      => standardOutputGroupCreate
  end type mergerTreeOutputterStandard

  interface mergerTreeOutputterStandard
     !% Constructors for the {\normalfont \ttfamily standard} merger tree outputter.
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface mergerTreeOutputterStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily standard} merger tree outputter class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (mergerTreeOutputterStandard)                :: self
    type   (inputParameters            ), intent(inout) :: parameters
    class  (galacticFilterClass        ), pointer       :: galacticFilter_
    class  (cosmologyFunctionsClass    ), pointer       :: cosmologyFunctions_
    class  (nodePropertyExtractorClass ), pointer       :: nodePropertyExtractor_
    logical                                             :: outputReferences
    type   (varying_string             )                :: outputsGroupName

    !# <inputParameter>
    !#   <name>outputsGroupName</name>
    !#   <source>parameters</source>
    !#   <defaultValue>var_str('Outputs')</defaultValue>
    !#   <description>The name of the HDF5 group to which outputs will be written.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>outputReferences</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Specifies whether or not references to individual merger tree datasets should be output.</description>
    !# </inputParameter>
    !# <objectBuilder class="galacticFilter"        name="galacticFilter_"        source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    !# <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    self=mergerTreeOutputterStandard(outputsGroupName,outputReferences,galacticFilter_,cosmologyFunctions_,nodePropertyExtractor_)
    !# <inputParametersValidate source="parameters"   />
    !# <objectDestructor name="galacticFilter_"       />
    !# <objectDestructor name="cosmologyFunctions_"   />
    !# <objectDestructor name="nodePropertyExtractor_"/>
    return
  end function standardConstructorParameters

  function standardConstructorInternal(outputsGroupName,outputReferences,galacticFilter_,cosmologyFunctions_,nodePropertyExtractor_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily standard} merger tree outputter class.
    implicit none
    type   (mergerTreeOutputterStandard)                        :: self
    type   (varying_string             ), intent(in   )         :: outputsGroupName
    class  (galacticFilterClass        ), intent(in   ), target :: galacticFilter_
    class  (cosmologyFunctionsClass    ), intent(in   ), target :: cosmologyFunctions_
    class  (nodePropertyExtractorClass ), intent(in   ), target :: nodePropertyExtractor_
    logical                             , intent(in   )         :: outputReferences
    !# <constructorAssign variables="outputsGroupName, outputReferences, *galacticFilter_, *cosmologyFunctions_, *nodePropertyExtractor_"/>

    self%outputsGroupOpened      =.false.
    self%outputGroupsCount       = 0
    self%doublePropertyCount     =-1
    self%integerPropertyCount    =-1
    self%doublePropertiesWritten = 0
    self%integerPropertiesWritten= 0
    self%doubleBufferCount       = 0
    self%integerBufferCount      = 0
    self%integerBufferSize       =standardBufferSizeIncrement
    self%doubleBufferSize        =standardBufferSizeIncrement
    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !% Destructor for the {\normalfont \ttfamily standard} merger tree outputter class.
    implicit none
    type(mergerTreeOutputterStandard), intent(inout) :: self

    call self%finalize()
    !# <objectDestructor name="self%galacticFilter_"       />
    !# <objectDestructor name="self%cosmologyFunctions_"   />
    !# <objectDestructor name="self%nodePropertyExtractor_"/>
    return
  end subroutine standardDestructor

  subroutine standardOutput(self,tree,indexOutput,time,isLastOutput)
    !% Write properties of nodes in {\normalfont \ttfamily tree} to the \glc\ output file.
    use            :: Galacticus_Calculations_Resets, only : Galacticus_Calculations_Reset
    use            :: Galacticus_Error              , only : Galacticus_Error_Report
    use            :: Galacticus_Nodes              , only : mergerTree                   , nodeComponentBasic       , treeNode
    use            :: IO_HDF5                       , only : hdf5Access                   , hdf5Object
    use, intrinsic :: ISO_C_Binding                 , only : c_size_t
    use            :: Merger_Tree_Walkers           , only : mergerTreeWalkerAllNodes
    use            :: Multi_Counters                , only : multiCounter
    use            :: Node_Property_Extractors      , only : elementTypeDouble            , elementTypeInteger       , nodePropertyExtractorIntegerScalar, nodePropertyExtractorIntegerTuple, &
          &                                                  nodePropertyExtractorMulti   , nodePropertyExtractorNull, nodePropertyExtractorScalar       , nodePropertyExtractorTuple
    !# <include directive="mergerTreeOutputTask" type="moduleUse">
    include 'galacticus.output.merger_tree.tasks.modules.inc'
    !# </include>
    implicit none
    class           (mergerTreeOutputterStandard), intent(inout)           :: self
    type            (mergerTree                 ), intent(inout), target   :: tree
    integer         (c_size_t                   ), intent(in   )           :: indexOutput
    double precision                             , intent(in   )           :: time
    logical                                      , intent(in   ), optional :: isLastOutput
    type            (treeNode                   )               , pointer  :: node
    integer         (kind=hsize_t               ), dimension(1)            :: referenceLength , referenceStart
    class           (nodeComponentBasic         )               , pointer  :: basic
    type            (mergerTree                 )               , pointer  :: currentTree
    type            (mergerTreeWalkerAllNodes   )                          :: treeWalker
    integer                                                                :: doubleProperty  , integerProperty, &
         &                                                                    iProperty
    integer         (c_size_t                   )                          :: iGroup
    logical                                                                :: nodePassesFilter
    type            (hdf5Object                 )                          :: toDataset
    type            (multiCounter               )                          :: instance

    ! Main output block.
    !$omp critical(mergerTreeOutputterStandard)
    ! Create an output group.
    call self%outputGroupCreate(indexOutput,time)
    ! Iterate over trees.
    currentTree => tree
    do while (associated(currentTree))
       ! Get the base node of the tree.
       node => currentTree%baseNode
       ! Skip empty trees.
       if (associated(node)) then
          ! Initialize output buffers.
          ! Count up the number of properties to be output.
          call self%propertiesCount        (time,node)
          ! Ensure buffers are allocated for temporary property storage.
          call self%buffersAllocate        (indexOutput      )
          ! Get names for all properties to be output.
          call self%propertyNamesEstablish(time,node)
          ! Loop over all nodes in the tree.
          self%integerPropertiesWritten=0
          self%doublePropertiesWritten =0
          treeWalker=mergerTreeWalkerAllNodes(currentTree)
          do while (treeWalker%next(node))
             ! Get the basic component.
             basic => node%basic()
             if (basic%time() == time) then
                ! Reset calculations (necessary in case the last node to be evolved is the first one we output, in which case
                ! calculations would not be automatically reset because the node unique ID will not have changed).
                call Galacticus_Calculations_Reset (node)
                ! Test whether this node passes all output filters.
                nodePassesFilter=self%galacticFilter_%passes(node)
                if (nodePassesFilter) then
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
                      ! Populate the output buffers with properties. We first populate with any "extra" properites that may be
                      ! being computed, and then call the standard treeNode output method to populate with all "standard"
                      ! properties.
                      !# <include directive="mergerTreeOutputTask" type="functionCall" functionType="void">
                      !#  <functionArgs>node,integerProperty,self%integerBufferCount,self%integerBuffer,doubleProperty,self%doubleBufferCount,self%doubleBuffer,time,instance</functionArgs>
                      include 'galacticus.output.merger_tree.tasks.inc'
                      !# </include>
                      call node%output(integerProperty,self%integerBufferCount,self%integerBuffer,doubleProperty,self%doubleBufferCount,self%doubleBuffer,time,instance)
                      ! Handle any extracted properties.
                      select type (extractor_ => self%nodePropertyExtractor_)
                      type is (nodePropertyExtractorNull)
                         ! Null extractor - simply ignore.
                      class is (nodePropertyExtractorScalar       )
                         ! Scalar property extractor - extract and store the value.
                         self%doubleBuffer (self%doubleBufferCount ,doubleProperty+1                                                                  )=extractor_      %extract       (                  node     ,instance)
                         doubleProperty                                                                                                                =+doubleProperty                                                       &
                              &                                                                                                                         +1
                      class is (nodePropertyExtractorTuple        )
                         ! Tuple property extractor - extract and store the values.
                         self%doubleBuffer (self%doubleBufferCount ,doubleProperty+1 :doubleProperty+extractor_%elementCount(                    time))= extractor_     %extract       (                  node,time,instance)
                         doubleProperty                                                                                                                =+doubleProperty                                                       &
                              &                                                                                                                         +extractor_     %elementCount  (                       time         )
                      class is (nodePropertyExtractorIntegerScalar)
                         ! Integer scalar property extractor - extract and store the value.
                         self%integerBuffer(self%integerBufferCount,integerProperty+1                                                                 )=extractor_      %extract       (                  node,time,instance)
                         integerProperty                                                                                                               =+integerProperty                                                      &
                              &                                                                                                                         +1
                      class is (nodePropertyExtractorIntegerTuple )
                         ! Integer tuple property extractor - extract and store the values.
                         self%integerBuffer(self%integerBufferCount,integerProperty+1:integerProperty+extractor_%elementCount(                   time))= extractor_     %extract       (                  node,time,instance)
                         integerProperty                                                                                                               =+integerProperty                                                      &
                              &                                                                                                                         +extractor_     %elementCount  (                       time         )
                       class is (nodePropertyExtractorMulti        )
                         ! Multi property extractor - extract and store the values.
                         self%doubleBuffer (self%doubleBufferCount ,doubleProperty +1:doubleProperty +extractor_%elementCount(elementTypeDouble ,time))=extractor_      %extractDouble (                  node,time,instance)
                         doubleProperty                                                                                                                =+doubleProperty                                                       &
                              &                                                                                                                         +extractor_     %elementCount  (elementTypeDouble,     time         )
                         self%integerBuffer(self%integerBufferCount,integerProperty+1:integerProperty+extractor_%elementCount(elementTypeInteger,time))=extractor_      %extractInteger(                  node,time,instance)
                         integerProperty                                                                                                               =+integerProperty                                                      &
                              &                                                                                                                         +extractor_     %elementCount  (elementTypeInteger,    time         )
                      class default
                         call Galacticus_Error_Report('unsupported property extractor class'//{introspection:location})
                      end select
                      ! If buffer is full, dump it to file.
                      if (self%integerBufferCount == self%integerBufferSize) call self%extendIntegerBuffer()
                      if (self% doubleBufferCount == self% doubleBufferSize) call self%extendDoubleBuffer ()
                   end do
                end if
             end if
          end do
          ! Finished output.
          if (self%integerPropertyCount > 0 .and. self%integerBufferCount > 0) call self%dumpIntegerBuffer(indexOutput)
          if (self% doublePropertyCount > 0 .and. self% doubleBufferCount > 0) call self%dumpDoubleBuffer (indexOutput)
          ! Compute the start and length of regions to reference.
          !$ call hdf5Access%set()
          referenceLength(1)=max(self%integerPropertiesWritten,self%doublePropertiesWritten)
          if      (allocated(self%integerPropertyNames).and.self%outputGroups(indexOutput)%nodeDataGroup%hasDataset(self%integerPropertyNames(1))) then
             toDataset=self%outputGroups(indexOutput)%nodeDataGroup%openDataset(self%integerPropertyNames(1))
             referenceStart(1)=toDataset%size(1)-referenceLength(1)
             call toDataset%close()
          else if (allocated(self% doublePropertyNames).and.self%outputGroups(indexOutput)%nodeDataGroup%hasDataset(self% doublePropertyNames(1))) then
             toDataset=self%outputGroups(indexOutput)%nodeDataGroup%openDataset(self% doublePropertyNames(1))
             referenceStart(1)=toDataset%size(1)-referenceLength(1)
             call toDataset%close()
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
                   toDataset=self%outputGroups(indexOutput)%nodeDataGroup%openDataset(self%integerPropertyNames(iProperty))
                   call currentTree%hdf5Group%createReference1D(toDataset,self%integerPropertyNames(iProperty),referenceStart+1,referenceLength)
                   call toDataset%close()
                end do
             end if
             if (self%doublePropertyCount > 0  .and. self%doublePropertiesWritten  > 0) then
                do iProperty=1,self%doublePropertyCount
                   toDataset=self%outputGroups(indexOutput)%nodeDataGroup%openDataset(self%doublePropertyNames(iProperty))
                   call currentTree%hdf5Group%createReference1D(toDataset,self%doublePropertyNames(iProperty),referenceStart+1,referenceLength)
                   call toDataset%close()
                end do
             end if
             ! Close the tree group.
             call currentTree%hdf5Group%close()
          end if
          ! Store the start position and length of the node data for this tree, along with its volume weight.
          call self%outputGroups(indexOutput)%hdf5Group%writeDataset([currentTree%index]       ,"mergerTreeIndex"     ,"Index of each merger tree."                                  ,appendTo=.true.)
          call self%outputGroups(indexOutput)%hdf5Group%writeDataset(referenceStart            ,"mergerTreeStartIndex","Index in nodeData datasets at which each merger tree begins.",appendTo=.true.)
          call self%outputGroups(indexOutput)%hdf5Group%writeDataset(referenceLength           ,"mergerTreeCount"     ,"Number of nodes in nodeData datasets for each merger tree."  ,appendTo=.true.)
          call self%outputGroups(indexOutput)%hdf5Group%writeDataset([currentTree%volumeWeight],"mergerTreeWeight"    ,"Number density of each tree [Mpc⁻³]."                        ,appendTo=.true.)
          !$ call hdf5Access%unset()
       end if
       ! Skip to the next tree.
       currentTree => currentTree%nextTree
    end do
    ! Close down if this is the final output.
    if (present(isLastOutput)) then
       if (isLastOutput) then
          ! Close any open output groups.
          !$ call hdf5Access%set()
          do iGroup=1,self%outputGroupsCount
             if (self%outputGroups(iGroup)%opened) then
                if (self%outputGroups(iGroup)%nodeDataGroup%isOpen()) call self%outputGroups(iGroup)%nodeDataGroup%close()
                if (self%outputGroups(iGroup)%hdf5Group    %isOpen()) call self%outputGroups(iGroup)%hdf5Group    %close()
             end if
          end do
          if (self%outputsGroup%isOpen()) call self%outputsGroup%close()
          !$ call hdf5Access%unset()
       end if
    end if
    !$omp end critical(mergerTreeOutputterStandard)
    return
  end subroutine standardOutput

  subroutine standardFinalize(self)
    !% Finalize merger tree output by closing any open groups.
    use :: IO_HDF5, only : hdf5Access
    implicit none
    class  (mergerTreeOutputterStandard), intent(inout) :: self
    integer(c_size_t                   )                :: iGroup, iDataset

    ! Close any open output groups.
    !$ call hdf5Access%set()
    do iGroup=1,self%outputGroupsCount
       if (self%outputGroups(iGroup)%opened) then
          if (allocated(self%outputGroups(iGroup)%integerDataset)) then
             do iDataset=1,size(self%outputGroups(iGroup)%integerDataset,kind=c_size_t)
                if (self%outputGroups(iGroup)%integerDataset(iDataset)%isOpen()) call self%outputGroups(iGroup)%integerDataset(iDataset)%close()
             end do
          end if
          if (allocated(self%outputGroups(iGroup)% doubleDataset)) then
             do iDataset=1,size(self%outputGroups(iGroup)% doubleDataset)
                if (self%outputGroups(iGroup)% doubleDataset(iDataset)%isOpen()) call self%outputGroups(iGroup)% doubleDataset(iDataset)%close()
             end do
          end if
          if (self%outputGroups(iGroup)%nodeDataGroup%isOpen()) call self%outputGroups(iGroup)%nodeDataGroup%close()
          if (self%outputGroups(iGroup)%hdf5Group    %isOpen()) call self%outputGroups(iGroup)%hdf5Group    %close()
       end if
    end do
    if (self%outputsGroup%isOpen()) call self%outputsGroup%close()
    !$ call hdf5Access%unset()
    return
  end subroutine standardFinalize

  subroutine standardMakeGroup(self,tree,indexOutput)
    !% Make an group in the \glc\ file in which to store {\normalfont \ttfamily tree}.
    use            :: Galacticus_Nodes                , only : mergerTree
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Numerical_Constants_Astronomical, only : megaParsec
    use            :: String_Handling                 , only : operator(//)
    implicit none
    class  (mergerTreeOutputterStandard), intent(inout) :: self
    type   (mergerTree                 ), intent(inout) :: tree
    integer(c_size_t                   ), intent(in   ) :: indexOutput
    type   (varying_string             )                :: commentText, groupName

    ! Create a name for the group.
    groupName='mergerTree'
    groupName=groupName//tree%index
    ! Create a comment for the group.
    commentText='Data for nodes within merger tree ID='
    commentText=commentText//tree%index
    ! Create a group for the tree.
    tree%hdf5Group=self%outputGroups(indexOutput)%hdf5Group%openGroup(char(groupName),char(commentText))
    ! Add the merger tree weight to the group.
    call tree%hdf5Group%writeAttribute(tree%volumeWeight  ,"volumeWeight"         )
    call tree%hdf5Group%writeAttribute(1.0d0/megaParsec**3,"volumeWeightUnitsInSI")
    return
  end subroutine standardMakeGroup

  subroutine standardDumpIntegerBuffer(self,indexOutput)
    !% Dump the contents of the integer properties buffer to the \glc\ output file.
    use :: IO_HDF5, only : hdf5Access, hdf5DataTypeInteger8
    implicit none
    class  (mergerTreeOutputterStandard), intent(inout) :: self
    integer(c_size_t                   ), intent(in   ) :: indexOutput
    integer                                             :: iProperty

    ! Write integer data from the buffer.
    if (self%integerPropertyCount > 0) then
       !$ call hdf5Access%set()
       do iProperty=1,self%integerPropertyCount
          if (.not.self%outputGroups(indexOutput)%                                                 integerDataset         (iProperty)%isOpen())  &
               &   self%outputGroups(indexOutput)%                                                 integerDataset         (iProperty)=           &
               &   self%outputGroups(indexOutput)%nodeDataGroup%openDataset(                                                                     &
               &                                                                              self%integerPropertyNames   (iProperty)          , &
               &                                                                              self%integerPropertyComments(iProperty)          , &
               &                                                            datasetDataType  =hdf5DataTypeInteger8                             , &
               &                                                            datasetDimensions=[0_c_size_t]                                     , &
               &                                                            appendTo         =.true.                                             &
               &                                                           )
          call self%outputGroups(indexOutput)%integerDataset(iProperty)%writeDataset(self%integerBuffer(1:self%integerBufferCount,iProperty),self%integerPropertyNames(iProperty),self%integerPropertyComments(iProperty),appendTo=.true.)
          if (.not.self%outputGroups(indexOutput)%integerAttributesWritten.and.self%integerPropertyUnitsSI(iProperty) /= 0.0d0) &
               & call self%outputGroups(indexOutput)%integerDataset(iProperty)%writeAttribute(self%integerPropertyUnitsSI(iProperty),"unitsInSI")
       end do
       self%integerPropertiesWritten=self%integerPropertiesWritten+self%integerBufferCount
       self%integerBufferCount=0
       self%outputGroups(indexOutput)%integerAttributesWritten=.true.
       !$ call hdf5Access%unset()
    end if
    return
  end subroutine standardDumpIntegerBuffer

  subroutine standardDumpDoubleBuffer(self,indexOutput)
    !% Dump the contents of the double precision properties buffer to the \glc\ output file.
    use            :: IO_HDF5      , only : hdf5Access, hdf5DataTypeDouble
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class  (mergerTreeOutputterStandard), intent(inout) :: self
    integer(c_size_t                   ), intent(in   ) :: indexOutput
    integer                                             :: iProperty

    ! Write double data from the buffer.
    if (self%doublePropertyCount > 0) then
       !$ call hdf5Access%set()
       do iProperty=1,self%doublePropertyCount
          if (.not.self%outputGroups(indexOutput)%                                                 doubleDataset         (iProperty)%isOpen())  &
               &   self%outputGroups(indexOutput)%                                                 doubleDataset         (iProperty)=           &
               &   self%outputGroups(indexOutput)%nodeDataGroup%openDataset(                                                                    &
               &                                                                              self%doublePropertyNames   (iProperty)          , &
               &                                                                              self%doublePropertyComments(iProperty)          , &
               &                                                            datasetDataType  =hdf5DataTypeDouble                              , &
               &                                                            datasetDimensions=[0_c_size_t]                                    , &
               &                                                            appendTo         =.true.                                            &
               &                                                           )
          call self%outputGroups(indexOutput)% doubleDataset(iProperty)%writeDataset( self%doubleBuffer(1: self%doubleBufferCount,iProperty), self%doublePropertyNames(iProperty),self%doublePropertyComments(iProperty),appendTo=.true.)
          if (.not.self%outputGroups(indexOutput)% doubleAttributesWritten.and. self%doublePropertyUnitsSI(iProperty) /= 0.0d0) &
               & call self%outputGroups(indexOutput)% doubleDataset(iProperty)%writeAttribute( self%doublePropertyUnitsSI(iProperty),"unitsInSI")
       end do
       self%doublePropertiesWritten=self%doublePropertiesWritten+self%doubleBufferCount
       self%doubleBufferCount=0
       self%outputGroups(indexOutput)%doubleAttributesWritten=.true.
       !$ call hdf5Access%unset()
    end if
    return
  end subroutine standardDumpDoubleBuffer

  subroutine standardExtendIntegerBuffer(self)
    !% Extend the size of the integer buffer.
    use :: Memory_Management, only : allocateArray, deallocateArray
    implicit none
    class  (mergerTreeOutputterStandard), intent(inout)                 :: self
    integer(kind_int8                  ), allocatable  , dimension(:,:) :: integerBufferTemporary
    integer                                                             :: integerPropertyCountMaximum

    integerPropertyCountMaximum=size(self%integerBuffer,dim=2)
    call Move_Alloc   (self%integerBuffer, integerBufferTemporary                                                         )
    call allocateArray(self%integerBuffer,[self%integerBufferSize+standardBufferSizeIncrement,integerPropertyCountMaximum])
    self%integerBuffer(1:self%integerBufferSize,:)=integerBufferTemporary
    self%integerBufferSize=self%integerBufferSize+standardBufferSizeIncrement
    call deallocateArray(                  integerBufferTemporary                                                         )
    return
  end subroutine standardExtendIntegerBuffer

  subroutine standardExtendDoubleBuffer(self)
    !% Extend the size of the double buffer.
    use :: Memory_Management, only : allocateArray, deallocateArray
    implicit none
    class           (mergerTreeOutputterStandard), intent(inout)                 :: self
    double precision                             , allocatable  , dimension(:,:) :: doubleBufferTemporary
    integer                                                                      :: doublePropertyCountMaximum

    doublePropertyCountMaximum=size(self%doubleBuffer,dim=2)
    call Move_Alloc   (self%doubleBuffer, doubleBufferTemporary                                                        )
    call allocateArray(self%doubleBuffer,[self%doubleBufferSize+standardBufferSizeIncrement,doublePropertyCountMaximum])
    self%doubleBuffer(1:self%doubleBufferSize,:)=doubleBufferTemporary
    self%doubleBufferSize=self%doubleBufferSize+standardBufferSizeIncrement
    call deallocateArray(                 doubleBufferTemporary                                                        )
    return
  end subroutine standardExtendDoubleBuffer

  subroutine standardPropertiesCount(self,time,node)
    !% Count up the number of properties that will be output.
    use :: Galacticus_Error        , only : Galacticus_Error_Report
    use :: Galacticus_Nodes        , only : treeNode
    use :: Node_Property_Extractors, only : elementTypeDouble         , elementTypeInteger       , nodePropertyExtractorIntegerScalar, nodePropertyExtractorIntegerTuple, &
          &                                 nodePropertyExtractorMulti, nodePropertyExtractorNull, nodePropertyExtractorScalar       , nodePropertyExtractorTuple
    !# <include directive="mergerTreeOutputPropertyCount" type="moduleUse">
    include 'galacticus.output.merger_tree.property_count.modules.inc'
    !# </include>
    implicit none
    class           (mergerTreeOutputterStandard), intent(inout)          :: self
    double precision                             , intent(in   )          :: time
    type            (treeNode                   ), intent(inout), pointer :: node

    self%integerPropertyCount=0
    self%doublePropertyCount =0
    !# <include directive="mergerTreeOutputPropertyCount" type="functionCall" functionType="void">
    !#  <functionArgs>node,self%integerPropertyCount,self%doublePropertyCount,time</functionArgs>
    include 'galacticus.output.merger_tree.property_count.inc'
    !# </include>
    call node%outputCount(self%integerPropertyCount,self%doublePropertyCount,time)
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
    class is (nodePropertyExtractorIntegerScalar)
       ! Integer scalar property extractor - simply increment the integer property output count by one.
       self%integerPropertyCount=self%integerPropertyCount+1
    class is (nodePropertyExtractorIntegerTuple )
       ! Integer tuple property extractor - increment the integer property output count by the number of elements.
       self%integerPropertyCount=self%integerPropertyCount+extractor_%elementCount(time)
    class is (nodePropertyExtractorMulti        )
       ! Multi proprty extractor - increment double and integer property output counts.
       self%integerPropertyCount=self%integerPropertyCount+extractor_%elementCount(elementTypeInteger,time)
       self% doublePropertyCount=self% doublePropertyCount+extractor_%elementCount(elementTypeDouble ,time)
    class default
       call Galacticus_Error_Report('unsupported property extractor class'//{introspection:location})
    end select
    return
  end subroutine standardPropertiesCount

  subroutine standardBuffersAllocate(self,indexOutput)
    !% Allocate buffers for storage of properties.
    use, intrinsic :: ISO_C_Binding    , only : c_size_t
    use            :: Memory_Management, only : allocateArray, deallocateArray
    implicit none
    class  (mergerTreeOutputterStandard), intent(inout) :: self
    integer(c_size_t                   ), intent(in   ) :: indexOutput

    if (self%integerPropertyCount > 0 .and. (.not.allocated(self%integerBuffer) .or. self%integerPropertyCount > size(self%integerPropertyNames)) ) then
       if (allocated(self%integerBuffer)) then
          call deallocateArray(self%integerBuffer          )
          call deallocateArray(self%integerPropertyNames   )
          call deallocateArray(self%integerPropertyComments)
          call deallocateArray(self%integerPropertyUnitsSI )
       end if
       call allocateArray(self%integerBuffer          ,[self%integerBufferSize,self%integerPropertyCount])
       call allocateArray(self%integerPropertyNames   ,[                       self%integerPropertyCount])
       call allocateArray(self%integerPropertyComments,[                       self%integerPropertyCount])
       call allocateArray(self%integerPropertyUnitsSI ,[                       self%integerPropertyCount])
    end if
    if (self%doublePropertyCount  > 0 .and. (.not.allocated(self%doubleBuffer ) .or. self%doublePropertyCount  > size(self%doublePropertyNames ))) then
       if (allocated(self%doubleBuffer )) then
          call deallocateArray(self%doubleBuffer           )
          call deallocateArray(self%doublePropertyNames    )
          call deallocateArray(self%doublePropertyComments )
          call deallocateArray(self%doublePropertyUnitsSI  )
       end if
       call allocateArray(self%doubleBuffer           ,[self%doubleBufferSize ,self%doublePropertyCount ])
       call allocateArray(self%doublePropertyNames    ,[                       self%doublePropertyCount ])
       call allocateArray(self%doublePropertyComments ,[                       self%doublePropertyCount ])
       call allocateArray(self%doublePropertyUnitsSI  ,[                       self%doublePropertyCount ])
    end if
    ! Allocate datasets.
    if (.not.allocated(self%outputGroups(indexOutput)%integerDataset)) allocate(self%outputGroups(indexOutput)%integerDataset(self%integerPropertyCount))
    if (.not.allocated(self%outputGroups(indexOutput)% doubleDataset)) allocate(self%outputGroups(indexOutput)% doubleDataset( self%doublePropertyCount))
    return
  end subroutine standardBuffersAllocate

  subroutine standardPropertyNamesEstablish(self,time,node)
    !% Set names for the properties.
    use :: Galacticus_Error        , only : Galacticus_Error_Report
    use :: Galacticus_Nodes        , only : treeNode
    use :: Node_Property_Extractors, only : elementTypeDouble         , elementTypeInteger       , nodePropertyExtractorIntegerScalar, nodePropertyExtractorIntegerTuple, &
          &                                 nodePropertyExtractorMulti, nodePropertyExtractorNull, nodePropertyExtractorScalar       , nodePropertyExtractorTuple
    !# <include directive="mergerTreeOutputNames" type="moduleUse">
    include 'galacticus.output.merger_tree.names.modules.inc'
    !# </include>
    implicit none
    class           (mergerTreeOutputterStandard), intent(inout)          :: self
    double precision                             , intent(in   )          :: time
    type            (treeNode                   ), intent(inout), pointer :: node
    integer                                                               :: doubleProperty, integerProperty

    integerProperty=0
    doubleProperty =0
    !# <include directive="mergerTreeOutputNames" type="functionCall" functionType="void">
    !#  <functionArgs>node,integerProperty,self%integerPropertyNames,self%integerPropertyComments,self%integerPropertyUnitsSI,doubleProperty,self%doublePropertyNames,self%doublePropertyComments,self%doublePropertyUnitsSI,time</functionArgs>
    include 'galacticus.output.merger_tree.names.inc'
    !# </include>
    call node%outputNames(integerProperty,self%integerPropertyNames,self%integerPropertyComments,self%integerPropertyUnitsSI,doubleProperty,self%doublePropertyNames,self%doublePropertyComments,self%doublePropertyUnitsSI,time)
    ! Handle extracted properties.
    select type (extractor_ => self%nodePropertyExtractor_)
    type is (nodePropertyExtractorNull)
       ! Null extractor - simply ignore.
    class is (nodePropertyExtractorScalar       )
          ! Scalar property extractor - get the name, description, and units.
       self%doublePropertyNames    (doubleProperty +1                                                                 )=extractor_%name        (                       )
       self%doublePropertyComments (doubleProperty +1                                                                 )=extractor_%description (                       )
       self%doublePropertyUnitsSI  (doubleProperty +1                                                                 )=extractor_%unitsInSI   (                       )
       doubleProperty =doubleProperty +1
    class is (nodePropertyExtractorTuple        )
       ! Tuple property extractor - get the names, descriptions, and units.
       self%doublePropertyNames    (doubleProperty +1:doubleProperty +extractor_%elementCount(                   time))=extractor_%names       (                   time)
       self%doublePropertyComments (doubleProperty +1:doubleProperty +extractor_%elementCount(                   time))=extractor_%descriptions(                   time)
       self%doublePropertyUnitsSI  (doubleProperty +1:doubleProperty +extractor_%elementCount(                   time))=extractor_%unitsInSI   (                   time)
       doubleProperty =doubleProperty +extractor_%elementCount(time)
    class is (nodePropertyExtractorIntegerScalar)
       ! Integer scalar property extractor - get the name, description, and units.
       self%integerPropertyNames   (integerProperty+1                                                                 )=extractor_%name        (                       )
       self%integerPropertyComments(integerProperty+1                                                                 )=extractor_%description (                       )
       self%integerPropertyUnitsSI (integerProperty+1                                                                 )=extractor_%unitsInSI   (                       )
       integerProperty=integerProperty+1
    class is (nodePropertyExtractorIntegerTuple )
       ! Integer tuple property extractor - get the names, descriptions, and units.
       self%integerPropertyNames   (integerProperty+1:integerProperty+extractor_%elementCount(                   time))=extractor_%names       (                   time)
       self%integerPropertyComments(integerProperty+1:integerProperty+extractor_%elementCount(                   time))=extractor_%descriptions(                   time)
       self%integerPropertyUnitsSI (integerProperty+1:integerProperty+extractor_%elementCount(                   time))=extractor_%unitsInSI   (                   time)
       integerProperty=integerProperty+extractor_%elementCount(time)
    class is (nodePropertyExtractorMulti        )
       ! Multi proprty extractor - get the names, descriptions, and units.
       self%doublePropertyNames    (doubleProperty +1:doubleProperty +extractor_%elementCount(elementTypeDouble ,time))=extractor_%names       (elementTypeDouble ,time)
       self%doublePropertyComments (doubleProperty +1:doubleProperty +extractor_%elementCount(elementTypeDouble ,time))=extractor_%descriptions(elementTypeDouble ,time)
       self%doublePropertyUnitsSI  (doubleProperty +1:doubleProperty +extractor_%elementCount(elementTypeDouble ,time))=extractor_%unitsInSI   (elementTypeDouble ,time)
       doubleProperty =doubleProperty +extractor_%elementCount(elementTypeDouble ,time)
       self%integerPropertyNames   (integerProperty+1:integerProperty+extractor_%elementCount(elementTypeInteger,time))=extractor_%names       (elementTypeInteger,time)
       self%integerPropertyComments(integerProperty+1:integerProperty+extractor_%elementCount(elementTypeInteger,time))=extractor_%descriptions(elementTypeInteger,time)
       self%integerPropertyUnitsSI (integerProperty+1:integerProperty+extractor_%elementCount(elementTypeInteger,time))=extractor_%unitsInSI   (elementTypeInteger,time)
       integerProperty=integerProperty+extractor_%elementCount(elementTypeInteger,time)
    class default
       call Galacticus_Error_Report('unsupported property extractor class'//{introspection:location})
    end select
    return
  end subroutine standardPropertyNamesEstablish

  subroutine standardOutputGroupCreate(self,indexOutput,time)
    !% Create a group in which to store this output.
    use            :: Galacticus_HDF5                 , only : galacticusOutputFile
    use            :: IO_HDF5                         , only : hdf5Access
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Memory_Management               , only : Memory_Usage_Record
    use            :: Numerical_Constants_Astronomical, only : gigaYear
    use            :: String_Handling                 , only : operator(//)
    implicit none
    class           (mergerTreeOutputterStandard), intent(inout)               :: self
    integer         (c_size_t                   ), intent(in   )               :: indexOutput
    double precision                             , intent(in   )               :: time
    type            (outputGroup                ), allocatable  , dimension(:) :: outputGroupsTemporary
    type            (varying_string             )                              :: commentText          , groupName

    !$ call hdf5Access%set()
    ! Ensure group ID space is large enough.
    if (indexOutput > self%outputGroupsCount) then
       if (allocated(self%outputGroups)) then
          call Move_Alloc(self%outputGroups,outputGroupsTemporary)
          self%outputGroupsCount=max(self%outputGroupsCount+standardOutputGroupsIncrement,(indexOutput/standardOutputGroupsIncrement+1)*standardOutputGroupsIncrement)
          allocate(self%outputGroups(self%outputGroupsCount))
          self%outputGroups(1:size(outputGroupsTemporary))=outputGroupsTemporary
          self%outputGroups(size(outputGroupsTemporary)+1:size(self%outputGroups))%opened                  =.false.
          self%outputGroups(size(outputGroupsTemporary)+1:size(self%outputGroups))%integerAttributesWritten=.false.
          self%outputGroups(size(outputGroupsTemporary)+1:size(self%outputGroups))%doubleAttributesWritten =.false.
          deallocate(outputGroupsTemporary)
          call Memory_Usage_Record(sizeof(self%outputGroups(1)),blockCount=0)
       else
          self%outputGroupsCount=max(standardOutputGroupsIncrement,(indexOutput/standardOutputGroupsIncrement+1)*standardOutputGroupsIncrement)
          allocate(self%outputGroups(self%outputGroupsCount))
          self%outputGroups%opened                  =.false.
          self%outputGroups%integerAttributesWritten=.false.
          self%outputGroups%doubleAttributesWritten =.false.
          call Memory_Usage_Record(sizeof(self%outputGroups))
       end if
    end if
    ! Make the enclosing group if it has not been created.
    if (.not.self%outputsGroupOpened) then
       self%outputsGroup=galacticusOutputFile%openGroup(char(self%outputsGroupName),'Contains all outputs from Galacticus.')
       self%outputsGroupOpened=.true.
    end if
    !$ call hdf5Access%unset()
    ! Create the group if it has not been created.
    if (.not.self%outputGroups(indexOutput)%opened) then
       ! Create a name for the group.
       groupName='Output'
       groupName=groupName//indexOutput
       ! Create a comment for the group.
       commentText='Data for output number '
       commentText=commentText//indexOutput
       ! Create a group for the tree.
       !$ call hdf5Access%set()
       self%outputGroups(indexOutput)%hdf5Group    =self%outputsGroup                       %openGroup(char(groupName),char(commentText)                                   )
       self%outputGroups(indexOutput)%nodeDataGroup=self%outputGroups(indexOutput)%hdf5Group%openGroup("nodeData"     ,"Group containing data on all nodes at this output.")
       self%outputGroups(indexOutput)%opened                  =.true.
       self%outputGroups(indexOutput)%integerAttributesWritten=.false.
       self%outputGroups(indexOutput)%doubleAttributesWritten =.false.
       ! Add the time to this group.
       call self%outputGroups(indexOutput)%hdf5Group%writeAttribute(time                                          ,'outputTime'           )
       call self%outputGroups(indexOutput)%hdf5Group%writeAttribute(gigaYear                                      ,'timeUnitInSI'         )
       call self%outputGroups(indexOutput)%hdf5Group%writeAttribute(self%cosmologyFunctions_%expansionFactor(time),'outputExpansionFactor')
       !$ call hdf5Access%unset()
    end if
    return
  end subroutine standardOutputGroupCreate
