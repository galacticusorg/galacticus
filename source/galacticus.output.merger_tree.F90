!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements writing a merger tree to the \glc\ output file.

module Galacticus_Output_Merger_Tree
  !% Implements writing a merger tree to the \glc\ output file.
  use, intrinsic :: ISO_C_Binding
  use               ISO_Varying_String
  use               Galacticus_HDF5
  use               Kind_Numbers
  use               Galactic_Filters
  use               Output_Analyses
  implicit none
  private
  public :: Galacticus_Merger_Tree_Output, Galacticus_Merger_Tree_Output_Finalize

  ! Output groups.
  type outputGroup
     !% Type used for output group information.
     logical                                        :: doubleAttributesWritten, integerAttributesWritten, &
          &                                            opened
     type   (hdf5Object)                            :: hdf5Group              , nodeDataGroup
     integer(c_size_t  )                            :: length
     type   (hdf5Object), allocatable, dimension(:) :: integerDataset         , doubleDataset
  end type outputGroup
  
  type            (hdf5Object          )                                         :: outputsGroup
  logical                                                                        :: outputsGroupOpened         =.false.
  type            (outputGroup         )           , allocatable, dimension(:)   :: outputGroups
  integer         (c_size_t            )                                         :: outputGroupsCount          =0
  integer                               , parameter                              :: outputGroupsIncrement      =10

  ! Number of properties of each type.
  integer                                                                        :: doublePropertyCount        =-1     , integerPropertyCount    =-1
  !$omp threadprivate(integerPropertyCount,doublePropertyCount)
  ! Buffers for properties.
  integer                                                                        :: doublePropertiesWritten    =0      , integerPropertiesWritten=0
  integer                                                                        :: doubleBufferCount          =0      , integerBufferCount      =0
  integer                               , parameter                              :: bufferSizeIncrement        =1024   , commentLengthMax        =256, &
       &                                                                            nameLengthMax              =256
  integer                                                                        :: integerBufferSize          =bufferSizeIncrement, doubleBufferSize=bufferSizeIncrement
  integer         (kind=kind_int8      )           , allocatable, dimension(:,:) :: integerBuffer
  double precision                                 , allocatable, dimension(:,:) :: doubleBuffer
  character       (len=nameLengthMax   )           , allocatable, dimension(:)   :: doublePropertyNames                , integerPropertyNames
  character       (len=commentLengthMax)           , allocatable, dimension(:)   :: doublePropertyComments             , integerPropertyComments
  double precision                                 , allocatable, dimension(:)   :: doublePropertyUnitsSI              , integerPropertyUnitsSI
  !$omp threadprivate(integerPropertiesWritten,doublePropertiesWritten,integerBufferCount,doubleBufferCount,integerBuffer)
  !$omp threadprivate(doubleBuffer,integerPropertyNames,doublePropertyNames,integerPropertyUnitsSI,doublePropertyUnitsSI)
  !$omp threadprivate(integerPropertyComments,doublePropertyComments,integerBufferSize,doubleBufferSize)
  ! Flag indicating if module is initialized.
  logical                                                                        :: mergerTreeOutputInitialized=.false.

  ! Flag indicating if merger tree references are to be output.
  logical                                                                        :: outputReferences

  ! Option indicating whether galaxies are to be output.
  logical                                                                        :: outputMergerTrees

  ! Analyses to carry out.
  type            (varying_string      )          , dimension(:),   allocatable  :: analyses
  type            (outputAnalyses      )          , dimension(:),   allocatable  :: outputAnalysisList

  ! Filter to apply to output.
  class           (galacticFilterClass ), pointer                                :: outputFilter => null()
  !$omp threadprivate(outputFilter)

contains

  subroutine Galacticus_Merger_Tree_Output(tree,iOutput,time,isLastOutput)
    !% Write properties of nodes in {\normalfont \ttfamily tree} to the \glc\ output file.
    use, intrinsic :: ISO_C_Binding
    use               Galacticus_Calculations_Resets
    use               Galacticus_Nodes
    use               Galacticus_Output_Open
    use               Input_Parameters2
    use               Galactic_Structure_Radii
    use               Galacticus_Output_Merger_Tree_Data
    use               Multi_Counters
    !# <include directive="mergerTreeOutputTask" type="moduleUse">
    include 'galacticus.output.merger_tree.tasks.modules.inc'
    !# </include>
    !# <include directive="mergerTreeOutputInstances" type="moduleUse">
    include 'galacticus.output.merger_tree.instances.modules.inc'
    !# </include>
    !# <include directive="mergerTreeExtraOutputTask" type="moduleUse">
    include 'galacticus.output.merger_tree.tasks.extra.modules.inc'
    !# </include>
    !# <include directive="mergerTreeAnalysisTask" type="moduleUse">
    include 'galacticus.output.merger_tree.analysis.modules.inc'
    !# </include>
    ! Define two inputParameter ojbects here, even though they are both used to access the same "mergerTreeOutput" parameter
    ! block. This is necessary to trigger automatic finalization of the objects when this function is exited. That finalization
    ! cleans up (e.g. closes the associated HDF5 group). Otherwise we could do this manually, but it's easier to simply have it
    ! done automatically. They are defined on separate lines as, currently, the "objectBuilder" preprocessor directive is
    ! insufficiently intelligent to modify the attributes of these variables if they are on a single line.
    implicit none
    type            (mergerTree        )              , intent(inout), target   :: tree
    integer         (c_size_t          )              , intent(in   )           :: iOutput
    double precision                                  , intent(in   )           :: time
    logical                                           , intent(in   ), optional :: isLastOutput
    type            (treeNode          ), pointer                               :: node                  , nodeNext
    integer         (kind=HSIZE_T      ), dimension(1)                          :: referenceLength       , referenceStart
    class           (nodeComponentBasic), pointer                               :: basic
    type            (mergerTree        ), pointer                               :: currentTree
    integer                                                                     :: doubleProperty        , nodeStatus     , &
         &                                                                         iProperty             , integerProperty, &
         &                                                                         i                     , analysisCount
    integer         (c_size_t          )                                        :: iGroup
    logical                                                                     :: nodePassesFilter      , finished
    type            (hdf5Object        )                                        :: toDataset
    type            (inputParameters   )                                        :: outputParameters
    type            (inputParameters   )                                        :: outputParametersFilter
    type            (multiCounter      )                                        :: instance
    
    ! Initialize if necessary.
    if (.not.mergerTreeOutputInitialized) then
       !$omp critical(Merger_Tree_Output_Initialization)
       if (.not.mergerTreeOutputInitialized) then
          ! Ensure file is open.
          call Galacticus_Output_Open_File
          outputParameters=globalParameters%subParameters('mergerTreeOutput',requirePresent=.false.,requireValue=.false.)
          allocate(analyses(outputParameters%count('analyses',zeroIfNotPresent=.true.)))
          !# <inputParameter>
          !#   <name>outputMergerTrees</name>
          !#   <source>outputParameters</source>
          !#   <defaultValue>.true.</defaultValue>
          !#   <description>Specifies whether or not to output the content of merger trees.</description>
          !#   <type>boolean</type>
          !#   <cardinality>1</cardinality>
          !#   <group>output</group>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>outputReferences</name>
          !#   <source>outputParameters</source>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>Specifies whether or not references to individual merger tree datasets should be output.</description>
          !#   <type>boolean</type>
          !#   <cardinality>1</cardinality>
          !# </inputParameter>
          analysisCount=outputParameters%copiesCount('outputAnalysisMethod',zeroIfNotPresent=.true.)
          if (analysisCount > 0) then
             allocate(outputAnalysisList(analysisCount))
             do i=1,analysisCount
                !# <objectBuilder class="outputAnalysis" name="outputAnalysisList(i)%outputAnalysis" source="outputParameters" copy="i"/>
             end do
          end if
          if (outputParameters%isPresent('analyses')) then
             !# <inputParameter>
             !#   <name>analyses</name>
             !#   <source>outputParameters</source>
             !#   <description>List of analyses to carry out on merger trees.</description>
             !#   <type>string</type>
             !#   <cardinality>1</cardinality>
             !# </inputParameter>
          end if
          ! Flag that the module is now initialized.
          mergerTreeOutputInitialized=.true.
       end if
       !$omp end critical(Merger_Tree_Output_Initialization)
    end if
    ! Get a galactic filter for determining which nodes are to be output. Note that this is done separately from the parameter
    ! reading above, as we need to get a per-thread object here, such that it matches up precisely with any per-thread objects
    ! used later in the output process. Failure to do this can lead to problems. For example, in lightcone output, the output
    ! filter here checks for presence of a node in the lightcone via the lightcone geometry object. Later, the same lightcone
    ! geometry object is used to find the position of the node within the lightcone. If two different objects are used then there
    ! can be numerical errors (e.g. differences in the distance-redshift relation) which cause a node to pass this filter but to
    ! not be found in the lightcone when attempting to find its position.
    if (.not.associated(outputFilter)) then
       outputParametersFilter=globalParameters%subParameters('mergerTreeOutput',requirePresent=.false.,requireValue=.false.)
       !# <objectBuilder class="galacticFilter" name="outputFilter" source="outputParametersFilter"/>          
    end if    
    ! Analysis block.
    ! Iterate over trees.
    currentTree => tree
    do while (associated(currentTree))     
       ! Iterate over nodes.
       node   => currentTree%baseNode
       nodeStatus =  nodeStatusFirst       
       do while (associated(node))
          ! Find the next node.
          nodeNext => node%walkTreeWithSatellites()
          ! Check for final node.
          if (.not.associated(nodeNext)) nodeStatus=nodeStatusLast
          ! Get the basic component.
          basic => node%basic()
          if (basic%time() == time) then
             ! Perform analysis tasks.
             !# <include directive="mergerTreeAnalysisTask" type="functionCall" functionType="void">
             !#  <functionArgs>currentTree,node,nodeStatus,iOutput,analyses</functionArgs>
             include 'galacticus.output.merger_tree.analysis.inc'
             !# </include>
             if (allocated(outputAnalysisList)) then
                do i=1,size(outputAnalysisList)
                   call outputAnalysisList(i)%outputAnalysis%analyze(node,iOutput)
                end do
             end if             
          end if
          ! Move to the next node.
          nodeStatus=nodeStatusNull
          node => nodeNext
       end do
       ! Record end of tree.
       node   => currentTree%baseNode
       nodeStatus =  nodeStatusFinal
       include 'galacticus.output.merger_tree.analysis.inc'
       ! Skip to the next tree.
       currentTree => currentTree%nextTree
    end do
    
    ! Main output block.
    if (outputMergerTrees) then
       !$omp critical(Merger_Tree_Output)
       ! Create an output group.
       call Make_Output_Group(iOutput,time)
       ! Iterate over trees.
       currentTree => tree
       do while (associated(currentTree))          
          ! Get the base node of the tree.
          node => currentTree%baseNode
          ! Skip empty trees.
          if (associated(node)) then          
             ! Initialize output buffers.
             ! Count up the number of properties to be output.
             call Count_Properties        (time,node)          
             ! Ensure buffers are allocated for temporary property storage.
             call Allocate_Buffers        (iOutput      )
             ! Get names for all properties to be output.
             call Establish_Property_Names(time,node)
             ! Loop over all nodes in the tree.
             integerPropertiesWritten=0
             doublePropertiesWritten =0
             finished                =.false.
             do while (.not.finished)
                node => node%walkTreeWithSatellites()
                if (.not.associated(node)) then
                   ! When a null pointer is returned, the full tree has been walked. We then have to
                   ! handle the base node as a special case.
                   node => currentTree%baseNode
                   finished =  .true.
                end if
                ! Get the basic component.
                basic => node%basic()
                if (basic%time() == time) then
                   ! Reset calculations (necessary in case the last node to be evolved is the first one we output, in which case
                   ! calculations would not be automatically reset because the node unique ID will not have changed).
                   call Galacticus_Calculations_Reset (node)
                   ! Ensure that galactic structure is up to date.
                   call Galactic_Structure_Radii_Solve(node)
                   ! Test whether this node passes all output filters.
                   nodePassesFilter=outputFilter%passes(node)
                   if (nodePassesFilter) then
                      ! Initialize the instance counter.
                      instance=multiCounter([1_c_size_t])
                      !# <include directive="mergerTreeOutputInstances" type="functionCall" functionType="void">
                      !#  <functionArgs>node,iOutput,instance</functionArgs>
                      include 'galacticus.output.merger_tree.instances.inc'
                      !# </include>
                      do while (instance%increment())
                         if (integerPropertyCount > 0) then
                            integerProperty=0
                            integerBufferCount=integerBufferCount+1
                         end if
                         if (doublePropertyCount > 0) then
                            doubleProperty=0
                            doubleBufferCount=doubleBufferCount+1
                         end if
                         ! Populate the output buffers with properties. We first populate with any "extra" properites that may be
                         ! being computed, and then call the standard treeNode output method to populate with all "standard"
                         ! properties. This order of operation is chosen because the treeNode output method will perform some clean
                         ! up (e.g. removing stellar luminosities no longer needed after this output). These quantities may be
                         ! required by the "extra" property calculations, thus requiring the extra property calculations to be
                         ! carried out first.                      
                         !# <include directive="mergerTreeOutputTask" type="functionCall" functionType="void">
                         !#  <functionArgs>node,integerProperty,integerBufferCount,integerBuffer,doubleProperty,doubleBufferCount,doubleBuffer,time,instance</functionArgs>
                         include 'galacticus.output.merger_tree.tasks.inc'
                         !# </include>
                         call node%output(integerProperty,integerBufferCount,integerBuffer,doubleProperty,doubleBufferCount,doubleBuffer,time,instance)
                         ! If buffer is full, dump it to file.
                         if (integerBufferCount == integerBufferSize) call Integer_Buffer_Extend()
                         if (doubleBufferCount  ==  doubleBufferSize) call  Double_Buffer_Extend()
                      end do                      
                   end if
                   ! Do any extra output tasks.
                   !# <include directive="mergerTreeExtraOutputTask" type="functionCall" functionType="void">
                   !#  <functionArgs>node,iOutput,currentTree%index,nodePassesFilter</functionArgs>
                   include 'galacticus.output.merger_tree.tasks.extra.inc'
                   !# </include>
                end if
             end do
             ! Finished output.
             if (integerPropertyCount > 0 .and. integerBufferCount > 0) call Integer_Buffer_Dump(iOutput)
             if (doublePropertyCount  > 0 .and. doubleBufferCount  > 0) call  Double_Buffer_Dump(iOutput)
             ! Compute the start and length of regions to reference.
             !$omp critical(HDF5_Access)
             referenceLength(1)=max(integerPropertiesWritten,doublePropertiesWritten)
             referenceStart (1)=outputGroups(iOutput)%length
             ! Create references to the datasets if requested.
             if (outputReferences) then
                ! Ensure that a group has been made for this merger tree.
                call Galacticus_Merger_Tree_Output_Make_Group(currentTree,iOutput)             
                ! Create references for this tree.
                if (integerPropertyCount > 0 .and. integerPropertiesWritten > 0) then
                   do iProperty=1,integerPropertyCount
                      toDataset=outputGroups(iOutput)%nodeDataGroup%openDataset(integerPropertyNames(iProperty))
                      call currentTree%hdf5Group%createReference1D(toDataset,integerPropertyNames(iProperty),referenceStart+1,referenceLength)
                      call toDataset%close()
                   end do
                end if
                if (doublePropertyCount > 0  .and. doublePropertiesWritten  > 0) then
                   do iProperty=1,doublePropertyCount
                      toDataset=outputGroups(iOutput)%nodeDataGroup%openDataset(doublePropertyNames(iProperty))
                      call currentTree%hdf5Group%createReference1D(toDataset,doublePropertyNames(iProperty),referenceStart+1,referenceLength)
                      call toDataset%close()
                   end do
                end if
                ! Close the tree group.
                call currentTree%hdf5Group%close()
             end if
             ! Store the start position and length of the node data for this tree, along with its volume weight.
             call outputGroups(iOutput)%hdf5Group%writeDataset([currentTree%index]       ,"mergerTreeIndex"     ,"Index of each merger tree."                                  ,appendTo=.true.)
             call outputGroups(iOutput)%hdf5Group%writeDataset(referenceStart            ,"mergerTreeStartIndex","Index in nodeData datasets at which each merger tree begins.",appendTo=.true.)
             call outputGroups(iOutput)%hdf5Group%writeDataset(referenceLength           ,"mergerTreeCount"     ,"Number of nodes in nodeData datasets for each merger tree."  ,appendTo=.true.)
             call outputGroups(iOutput)%hdf5Group%writeDataset([currentTree%volumeWeight],"mergerTreeWeight"    ,"Number density of each tree [Mpc⁻³]."                        ,appendTo=.true.)
             ! Increment the number of nodes written to this output group.
             outputGroups(iOutput)%length=outputGroups(iOutput)%length+referenceLength(1)
             !$omp end critical(HDF5_Access)
          end if
          ! Skip to the next tree.
          currentTree => currentTree%nextTree
       end do
       ! Close down if this is the final output.
       if (present(isLastOutput)) then
          if (isLastOutput) then
             ! Close any open output groups.
             !$omp critical(HDF5_Access)
             do iGroup=1,outputGroupsCount
                if (outputGroups(iGroup)%opened) then
                   if (outputGroups(iGroup)%nodeDataGroup%isOpen()) call outputGroups(iGroup)%nodeDataGroup%close()
                   if (outputGroups(iGroup)%hdf5Group    %isOpen()) call outputGroups(iGroup)%hdf5Group    %close()
                end if
             end do
             if (outputsGroup%isOpen()) call outputsGroup%close()
             !$omp end critical(HDF5_Access)
          end if
       end if
       !$omp end critical(Merger_Tree_Output)
    end if
    return
  end subroutine Galacticus_Merger_Tree_Output

  subroutine Galacticus_Merger_Tree_Output_Finalize()
    !% Finalize merger tree output by closing any open groups and finalizing analyses.
    implicit none
    integer(c_size_t) :: iGroup, iDataset
    integer           :: i

    ! Close any open output groups.
    !$omp critical(HDF5_Access)
    do iGroup=1,outputGroupsCount
       if (outputGroups(iGroup)%opened) then
          if (allocated(outputGroups(iGroup)%integerDataset)) then
             do iDataset=1,size(outputGroups(iGroup)%integerDataset,kind=c_size_t)
                if (outputGroups(iGroup)%integerDataset(iDataset)%isOpen()) call outputGroups(iGroup)%integerDataset(iDataset)%close()
             end do
          end if
          if (allocated(outputGroups(iGroup)% doubleDataset)) then
             do iDataset=1,size(outputGroups(iGroup)% doubleDataset)
                if (outputGroups(iGroup)% doubleDataset(iDataset)%isOpen()) call outputGroups(iGroup)% doubleDataset(iDataset)%close()
             end do
          end if
          if (outputGroups(iGroup)%nodeDataGroup%isOpen()) call outputGroups(iGroup)%nodeDataGroup%close()
          if (outputGroups(iGroup)%hdf5Group    %isOpen()) call outputGroups(iGroup)%hdf5Group    %close()
       end if
    end do
    if (outputsGroup%isOpen()) call outputsGroup%close()
    !$omp end critical(HDF5_Access)
    ! Finalize analyses.
    if (allocated(outputAnalysisList)) then
       do i=1,size(outputAnalysisList)
          call outputAnalysisList(i)%outputAnalysis%finalize()
       end do
       deallocate(outputAnalysisList)
    end if
    return
  end subroutine Galacticus_Merger_Tree_Output_Finalize

  subroutine Galacticus_Merger_Tree_Output_Make_Group(tree,iOutput)
    !% Make an group in the \glc\ file in which to store {\normalfont \ttfamily tree}.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Astronomical
    use String_Handling
    use Galacticus_Nodes
    implicit none
    type   (mergerTree    ), intent(inout) :: tree
    integer(c_size_t      ), intent(in   ) :: iOutput
    type   (varying_string)                :: commentText, groupName

    ! Create a name for the group.
    groupName='mergerTree'
    groupName=groupName//tree%index

    ! Create a comment for the group.
    commentText='Data for nodes within merger tree ID='
    commentText=commentText//tree%index

    ! Create a group for the tree.
    tree%hdf5Group=outputGroups(iOutput)%hdf5Group%openGroup(char(groupName),char(commentText))

    ! Add the merger tree weight to the group.
    call tree%hdf5Group%writeAttribute(tree%volumeWeight,"volumeWeight"         )
    call tree%hdf5Group%writeAttribute(1.0d0/megaParsec**3  ,"volumeWeightUnitsInSI")
    return
  end subroutine Galacticus_Merger_Tree_Output_Make_Group

  subroutine Integer_Buffer_Dump(iOutput)
    !% Dump the contents of the integer properties buffer to the \glc\ output file.
    use, intrinsic :: ISO_C_Binding
    implicit none
    integer(c_size_t), intent(in   ) :: iOutput
    integer                          :: iProperty

    ! Write integer data from the buffer.
    if (integerPropertyCount > 0) then
       !$omp critical(HDF5_Access)
       do iProperty=1,integerPropertyCount
          if (.not.outputGroups(iOutput)%                          integerDataset         (iProperty)%isOpen())             &
               &   outputGroups(iOutput)%                          integerDataset         (iProperty)=                      &
               &   outputGroups(iOutput)%nodeDataGroup%openDataset(                                                         &
               &                                                   integerPropertyNames   (iProperty)                     , &
               &                                                   integerPropertyComments(iProperty)                     , &
               &                                                   datasetDataType                   =hdf5DataTypeInteger8, &
               &                                                   datasetDimensions                 =[0_c_size_t]        , &
               &                                                   appendTo                          =.true.                &
               &                                                  )
          call outputGroups(iOutput)%integerDataset(iProperty)%writeDataset(integerBuffer(1:integerBufferCount,iProperty),integerPropertyNames(iProperty) &
               &,integerPropertyComments(iProperty),appendTo=.true.)
          if (.not.outputGroups(iOutput)%integerAttributesWritten.and.integerPropertyUnitsSI(iProperty) /= 0.0d0) &
               & call outputGroups(iOutput)%integerDataset(iProperty)%writeAttribute(integerPropertyUnitsSI(iProperty),"unitsInSI")
       end do
       integerPropertiesWritten=integerPropertiesWritten+integerBufferCount
       integerBufferCount=0
       outputGroups(iOutput)%integerAttributesWritten=.true.
       !$omp end critical(HDF5_Access)
    end if
    return
  end subroutine Integer_Buffer_Dump

  subroutine Integer_Buffer_Extend()
    !% Extend the size of the integer buffer.
    use Memory_Management
    implicit none
    integer(kind=kind_int8), allocatable, dimension(:,:) :: integerBufferTemporary
    integer                                              :: integerPropertyCountMaximum

    integerPropertyCountMaximum=size(integerBuffer,dim=2)
    call Move_Alloc   (integerBuffer, integerBufferTemporary                                            )
    call allocateArray  (integerBuffer,[integerBufferSize+bufferSizeIncrement,integerPropertyCountMaximum])
    integerBuffer(1:integerBufferSize,:)=integerBufferTemporary
    integerBufferSize=integerBufferSize+bufferSizeIncrement
    call deallocateArray(               integerBufferTemporary                                            )
    return
  end subroutine Integer_Buffer_Extend

  subroutine Double_Buffer_Extend()
    !% Extend the size of the double buffer.
    use Memory_Management
    implicit none
    double precision, allocatable, dimension(:,:) :: doubleBufferTemporary
    integer                                       :: doublePropertyCountMaximum

    doublePropertyCountMaximum=size(doubleBuffer,dim=2)
    call Move_Alloc   (doubleBuffer, doubleBufferTemporary                                           )
    call allocateArray  (doubleBuffer,[doubleBufferSize+bufferSizeIncrement,doublePropertyCountMaximum])
    doubleBuffer(1:doubleBufferSize,:)=doubleBufferTemporary
    doubleBufferSize=doubleBufferSize+bufferSizeIncrement
    call deallocateArray(              doubleBufferTemporary                                           )
    return
  end subroutine Double_Buffer_Extend

  subroutine Double_Buffer_Dump(iOutput)
    !% Dump the contents of the double precision properties buffer to the \glc\ output file.
    use, intrinsic :: ISO_C_Binding
    implicit none
    integer(c_size_t  ), intent(in   ) :: iOutput
    integer                            :: iProperty

    ! Write double data from the buffer.
    if (doublePropertyCount > 0) then
       !$omp critical(HDF5_Access)
       do iProperty=1,doublePropertyCount
          if (.not.outputGroups(iOutput)%                           doubleDataset         (iProperty)%isOpen())             &
               &   outputGroups(iOutput)%                           doubleDataset         (iProperty)=                      &
               &   outputGroups(iOutput)%nodeDataGroup%openDataset(                                                         &
               &                                                    doublePropertyNames   (iProperty)                     , &
               &                                                    doublePropertyComments(iProperty)                     , &
               &                                                   datasetDataType                   =hdf5DataTypeDouble  , &
               &                                                   datasetDimensions                 =[0_c_size_t]        , &
               &                                                   appendTo                          =.true.                &
               &                                                  )
          call outputGroups(iOutput)% doubleDataset(iProperty)%writeDataset( doubleBuffer(1: doubleBufferCount,iProperty), doublePropertyNames(iProperty) &
               &, doublePropertyComments(iProperty),appendTo=.true.)
          if (.not.outputGroups(iOutput)% doubleAttributesWritten.and. doublePropertyUnitsSI(iProperty) /= 0.0d0) &
               & call outputGroups(iOutput)% doubleDataset(iProperty)%writeAttribute( doublePropertyUnitsSI(iProperty),"unitsInSI")
       end do
       doublePropertiesWritten=doublePropertiesWritten+doubleBufferCount
       doubleBufferCount=0
       outputGroups(iOutput)%doubleAttributesWritten=.true.
       !$omp end critical(HDF5_Access)
    end if
    return
  end subroutine Double_Buffer_Dump

  subroutine Count_Properties(time,node)
    !% Count up the number of properties that will be output.
    use Galacticus_Nodes
    !# <include directive="mergerTreeOutputPropertyCount" type="moduleUse">
    include 'galacticus.output.merger_tree.property_count.modules.inc'
    !# </include>
    implicit none
    double precision          , intent(in   )          :: time
    type            (treeNode), intent(inout), pointer :: node

    integerPropertyCount=0
    doublePropertyCount =0
    !# <include directive="mergerTreeOutputPropertyCount" type="functionCall" functionType="void">
    !#  <functionArgs>node,integerPropertyCount,doublePropertyCount,time</functionArgs>
    include 'galacticus.output.merger_tree.property_count.inc'
    !# </include>
    call node%outputCount(integerPropertyCount,doublePropertyCount,time)
    return
  end subroutine Count_Properties

  subroutine Allocate_Buffers(iOutput)
    !% Allocate buffers for storage of properties.
    use, intrinsic :: ISO_C_Binding
    use Memory_Management
    implicit none
    integer(c_size_t), intent(in   ) :: iOutput

    if (integerPropertyCount > 0 .and. (.not.allocated(integerBuffer) .or. integerPropertyCount > size(integerPropertyNames)) ) then
       if (allocated(integerBuffer)) then
          call deallocateArray(integerBuffer          )
          call deallocateArray(integerPropertyNames   )
          call deallocateArray(integerPropertyComments)
          call deallocateArray(integerPropertyUnitsSI )
       end if
       call allocateArray(integerBuffer          ,[integerBufferSize,integerPropertyCount])
       call allocateArray(integerPropertyNames                      ,[integerPropertyCount])
       call allocateArray(integerPropertyComments                   ,[integerPropertyCount])
       call allocateArray(integerPropertyUnitsSI                    ,[integerPropertyCount])
    end if
    if (doublePropertyCount  > 0 .and. (.not.allocated(doubleBuffer ) .or. doublePropertyCount  > size(doublePropertyNames ))) then
       if (allocated(doubleBuffer )) then
          call deallocateArray(doubleBuffer           )
          call deallocateArray(doublePropertyNames    )
          call deallocateArray(doublePropertyComments )
          call deallocateArray(doublePropertyUnitsSI  )
       end if
       call allocateArray(doubleBuffer           ,[doubleBufferSize,doublePropertyCount])
       call allocateArray(doublePropertyNames                      ,[doublePropertyCount])
       call allocateArray(doublePropertyComments                   ,[doublePropertyCount])
       call allocateArray(doublePropertyUnitsSI                    ,[doublePropertyCount])
    end if
    ! Allocate datasets.
    if (.not.allocated(outputGroups(iOutput)%integerDataset)) allocate(outputGroups(iOutput)%integerDataset(integerPropertyCount))
    if (.not.allocated(outputGroups(iOutput)% doubleDataset)) allocate(outputGroups(iOutput)% doubleDataset( doublePropertyCount))
    return
  end subroutine Allocate_Buffers

  subroutine Establish_Property_Names(time,node)
    !% Set names for the properties.
    use Galacticus_Nodes
    !# <include directive="mergerTreeOutputNames" type="moduleUse">
    include 'galacticus.output.merger_tree.names.modules.inc'
    !# </include>
    implicit none
    double precision          , intent(in   )          :: time
    type            (treeNode), intent(inout), pointer :: node
    integer                                            :: doubleProperty, integerProperty

    integerProperty=0
    doubleProperty =0
    !# <include directive="mergerTreeOutputNames" type="functionCall" functionType="void">
    !#  <functionArgs>node,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time</functionArgs>
    include 'galacticus.output.merger_tree.names.inc'
    !# </include>
    call node%outputNames(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    return
  end subroutine Establish_Property_Names

  subroutine Make_Output_Group(iOutput,time)
    !% Create a group in which to store this output.
    use, intrinsic :: ISO_C_Binding
    use String_Handling
    use Cosmology_Functions
    use Memory_Management
    use Numerical_Constants_Astronomical
    !# <include directive="outputGroupOutputTask" type="moduleUse">
    include 'galacticus.output.merger_tree.outputGroup.tasks.modules.inc'
    !# </include>
    implicit none
    integer         (c_size_t               ), intent(in   )               :: iOutput
    double precision                         , intent(in   )               :: time
    type            (outputGroup            ), allocatable  , dimension(:) :: outputGroupsTemporary
    class           (cosmologyFunctionsClass), pointer                     :: cosmologyFunctions_
    type            (varying_string         )                              :: commentText              , groupName

    !$omp critical (HDF5_Access)
    ! Ensure group ID space is large enough.
    if (iOutput > outputGroupsCount) then
       if (allocated(outputGroups)) then
          call Move_Alloc(outputGroups,outputGroupsTemporary)
          outputGroupsCount=max(outputGroupsCount+outputGroupsIncrement,(iOutput/outputGroupsIncrement+1)*outputGroupsIncrement)
          allocate(outputGroups(outputGroupsCount))
          outputGroups(1:size(outputGroupsTemporary))=outputGroupsTemporary
          outputGroups(size(outputGroupsTemporary)+1:size(outputGroups))%opened                  =.false.
          outputGroups(size(outputGroupsTemporary)+1:size(outputGroups))%integerAttributesWritten=.false.
          outputGroups(size(outputGroupsTemporary)+1:size(outputGroups))%doubleAttributesWritten =.false.
          outputGroups(size(outputGroupsTemporary)+1:size(outputGroups))%length                  =0
          deallocate(outputGroupsTemporary)
          call Memory_Usage_Record(sizeof(outputGroups(1)),blockCount=0)
       else
          outputGroupsCount=max(outputGroupsIncrement,(iOutput/outputGroupsIncrement+1)*outputGroupsIncrement)
          allocate(outputGroups(outputGroupsCount))
          outputGroups%opened                  =.false.
          outputGroups%integerAttributesWritten=.false.
          outputGroups%doubleAttributesWritten =.false.
          outputGroups%length=0
          call Memory_Usage_Record(sizeof(outputGroups))
       end if
    end if

    ! Make the enclosing group if it has not been created.
    if (.not.outputsGroupOpened) then
       outputsGroup=galacticusOutputFile%openGroup('Outputs','Contains all outputs from Galacticus.')
       outputsGroupOpened=.true.
    end if
    !$omp end critical(HDF5_Access)

    ! Create the group if it has not been created.
    if (.not.outputGroups(iOutput)%opened) then
       ! Create a name for the group.
       groupName='Output'
       groupName=groupName//iOutput
       ! Create a comment for the group.
       commentText='Data for output number '
       commentText=commentText//iOutput

       ! Create a group for the tree.
       !$omp critical(HDF5_Access)
       !@ <outputType>
       !@   <name>nodeData</name>
       !@   <description>A representation of the state of all nodes in the simulation at a given time. It consists of numerous datasets which gives the properties of nodes in all merger trees at that time.</description>
       !@ </outputType>
       outputGroups(iOutput)%hdf5Group    =outputsGroup                   %openGroup(char(groupName),char(commentText))
       outputGroups(iOutput)%nodeDataGroup=outputGroups(iOutput)%hdf5Group%openGroup("nodeData","Group containing data on all nodes at this output.")
       outputGroups(iOutput)%opened                  =.true.
       outputGroups(iOutput)%integerAttributesWritten=.false.
       outputGroups(iOutput)%doubleAttributesWritten =.false.

       ! Get the default cosmology functions object.
       cosmologyFunctions_ => cosmologyFunctions()
       ! Add the time to this group.
       call outputGroups(iOutput)%hdf5Group%writeAttribute(time                                     ,'outputTime'           )
       call outputGroups(iOutput)%hdf5Group%writeAttribute(gigaYear                                 ,'timeUnitInSI'         )
       call outputGroups(iOutput)%hdf5Group%writeAttribute(cosmologyFunctions_%expansionFactor(time),'outputExpansionFactor')
       !$omp end critical(HDF5_Access)

       ! Establish all other properties.
       !# <include directive="outputGroupOutputTask" type="functionCall" functionType="void">
       !#  <functionArgs>outputGroups(iOutput)%hdf5Group,time</functionArgs>
       include 'galacticus.output.merger_tree.outputGroup.tasks.inc'
       !# </include>

    end if
    return
  end subroutine Make_Output_Group

end module Galacticus_Output_Merger_Tree
