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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Implements a task which postprocesses galaxies read back from a Galacticus output file.
!!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Galacticus_Nodes        , only : nodeHierarchyWrapper
  use            :: Node_Components         , only : nodeComponentsWrapper
  use            :: Input_Parameters        , only : inputParameters
  use            :: Merger_Tree_Outputters  , only : mergerTreeOutputterClass
  use            :: Star_Formation_Histories, only : starFormationHistoryClass

  !![
  <task name="taskPostprocessOutputs" docformat="rst">
    <description>
    A task which postprocesses galaxies read back from the output file of a prior Galacticus model. For each galaxy in the file
    named by ``fileName`` a node is rebuilt from the properties which were output for it, and is then passed to the
    :galacticus-class:`mergerTreeOutputterClass` object for output. New properties---spectral energy distributions, emission
    line luminosities, dust attenuated luminosities, and so on---can therefore be computed for the galaxies of an existing
    model, without evolving that model again, by specifying the relevant :galacticus-class:`nodePropertyExtractorClass`
    objects.

    Only what was output can be rebuilt. The properties restored, where datasets for them are present in the file, are:

    * the time (from the ``time`` dataset, as written by :galacticus-class:`nodePropertyExtractorTime`, or, where there is no
      such dataset, the time of the output) and mass of the basic component---since each galaxy in a lightcone has its own
      time, the ``time`` dataset is required for lightcone outputs;
    * the stellar and gas masses, masses of metals in stars and gas, radius, angular momentum, and velocity of the disk and
      spheroid components (only the total mass of metals is restored, so properties which depend on the abundances of
      individual elements can not be computed);
    * the star formation histories of the disk and spheroid components---these are restored by the
      :galacticus-class:`starFormationHistoryClass` object, which must therefore be configured as it was in the model which
      generated the file;
    * the inclination (from the ``inclination`` dataset, in degrees, as written by
      :galacticus-class:`nodePropertyExtractorInclination`), which is then available to a
      :galacticus-class:`galacticInclinationRandom` object---where there is no such dataset, galaxies have zero inclination
      (i.e. are face-on).

    Nodes are rebuilt in isolation: no host halo, satellite status, or position in a merger tree is restored. (So, for
    example, the time at which a node was last isolated is reported as its current time.) A property extractor which requires
    anything else of a node can not be used. The node components must be configured as
    they were in the model which generated the file.

    Since node components are only partially reconstructed, their own output properties would duplicate (or, for properties
    not listed above, misrepresent) those in the original file. Setting ``outputComponentProperties`` to ``false`` in the
    :galacticus-class:`mergerTreeOutputterStandard` outputter restricts output to the properties of the
    ``nodePropertyExtractor``. When postprocessing large numbers of galaxies, setting ``countNodesBuffered`` in that outputter
    to a value of order 1000 also greatly reduces the time spent writing output.

    Galaxies are read from each ``Output`` group within the group named by ``outputsGroupName``. By default all galaxies are
    postprocessed. If ``selectionFileName`` is given then it must name an HDF5 file containing, for each output to be
    postprocessed, a dataset ``[outputsGroupName]/Output[i]/selection`` listing the indices (counting from zero) in the
    ``nodeData`` datasets of the galaxies to postprocess---outputs with no such dataset are skipped.
    Galaxies may be further selected by the :galacticus-class:`galacticFilterClass` object of the outputter.

    Each output written is of the same type, and has the same time, expansion factor, and comoving distance attributes, as
    the output from which its galaxies were read, provided that the outputter writes to a group named ``outputsGroupName``.

    Galaxies are not guaranteed to be output in the order in which they were read. The
    :galacticus-class:`nodePropertyExtractorIndexPostprocessSource` property extractor outputs the index of each galaxy in
    the file from which it was read, allowing each to be matched to its original entry.
    </description>
  </task>
  !!]
  type, extends(taskClass) :: taskPostprocessOutputs
     !!{RST
     Implementation of a task which postprocesses galaxies read back from a Galacticus output file.
     !!}
     private
     type   (inputParameters          ), pointer :: parameters            => null()
     class  (mergerTreeOutputterClass ), pointer :: mergerTreeOutputter_  => null()
     class  (starFormationHistoryClass), pointer :: starFormationHistory_ => null()
     type   (varying_string           )          :: fileName                       , outputsGroupName    , &
          &                                         selectionFileName
     integer(c_size_t                 )          :: countBatch
     integer                                     :: inclinationID                  , indexSourceID
     ! Manager for node class hierarchy and component initialization.
     type   (nodeComponentsWrapper    ), pointer :: nodeComponents_       => null()
     type   (nodeHierarchyWrapper     ), pointer :: nodeHierarchy_        => null()
     type   (resourceManager          )          :: nodeComponentsManager          , nodeHierarchyManager
   contains
     final     ::            postprocessOutputsDestructor
     procedure :: perform => postprocessOutputsPerform
  end type taskPostprocessOutputs

  interface taskPostprocessOutputs
     !!{RST
     Constructors for the :galacticus-class:`taskPostprocessOutputs` task.
     !!}
     module procedure postprocessOutputsConstructorParameters
     module procedure postprocessOutputsConstructorInternal
  end interface taskPostprocessOutputs

  ! Copies of objects used by each thread.
  class(mergerTreeOutputterClass), pointer :: postprocessOutputsOutputter_ => null()
  !$omp threadprivate(postprocessOutputsOutputter_)

  type :: postprocessOutputsColumn
     !!{RST
     Type used to store a batch of a scalar property read from an output file.
     !!}
     logical                                       :: isPresent=.false.
     double precision, allocatable, dimension(:)   :: values
  end type postprocessOutputsColumn

  type :: postprocessOutputsBatch
     !!{RST
     Type used to describe a batch of galaxies to be postprocessed.
     !!}
     integer(c_size_t) :: indexOutput, indexBegin, count
  end type postprocessOutputsBatch

  type :: postprocessOutputsSelection
     !!{RST
     Type used to store the selection of galaxies to be postprocessed in an output.
     !!}
     logical                                      :: isSelected=.false.
     integer(c_size_t), allocatable, dimension(:) :: indices
  end type postprocessOutputsSelection

  ! Indices of the scalar properties which are restored.
  integer, parameter :: postprocessColumnTime             = 1, postprocessColumnMassBasic      = 2, &
       &                postprocessColumnInclination      = 3, postprocessColumnCountNode      = 3
  integer, parameter :: postprocessColumnMassStellar      = 1, postprocessColumnMassGas        = 2, &
       &                postprocessColumnMassMetalsStellar= 3, postprocessColumnMassMetalsGas  = 4, &
       &                postprocessColumnRadius           = 5, postprocessColumnAngularMomentum= 6, &
       &                postprocessColumnVelocity         = 7, postprocessColumnCountComponent = 7
  integer, parameter :: postprocessComponentDisk          = 1, postprocessComponentSpheroid    = 2, &
       &                postprocessComponentCount         = 2
  integer, parameter :: postprocessColumnCount            =+postprocessColumnCountNode      &
       &                                                   +postprocessColumnCountComponent &
       &                                                   *postprocessComponentCount

contains

  function postprocessOutputsConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`taskPostprocessOutputs` task class which takes a parameter set as input.
    !!}
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize
    use :: Node_Components , only : Node_Components_Initialize
    implicit none
    type   (taskPostprocessOutputs   )                        :: self
    type   (inputParameters          ), intent(inout), target :: parameters
    class  (mergerTreeOutputterClass ), pointer               :: mergerTreeOutputter_
    class  (starFormationHistoryClass), pointer               :: starFormationHistory_
    class  (*                        ), pointer               :: dummyPointer_
    type   (inputParameters          ), pointer               :: parametersRoot
    type   (varying_string           )                        :: fileName             , outputsGroupName, &
         &                                                       selectionFileName
    integer(c_size_t                 )                        :: countBatch

    ! Ensure the nodes objects are initialized.
    if (associated(parameters%parent)) then
       parametersRoot => parameters%parent
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       call nodeClassHierarchyInitialize(parametersRoot)
       call Node_Components_Initialize  (parametersRoot)
    else
       parametersRoot => null()
       call nodeClassHierarchyInitialize(parameters    )
       call Node_Components_Initialize  (parameters    )
    end if
    !![
    <inputParameter docformat="rst">
      <name>fileName</name>
      <description>
      The name of the Galacticus output file from which galaxies should be read.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>outputsGroupName</name>
      <defaultValue>var_str('Outputs')</defaultValue>
      <description>
      The name of the group, in the file from which galaxies are read, which contains the outputs.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>selectionFileName</name>
      <defaultValue>var_str('none')</defaultValue>
      <description>
      The name of an HDF5 file containing, for each output to be postprocessed, a dataset
      ``[outputsGroupName]/Output[i]/selection`` listing the indices (counting from zero) of the galaxies to
      postprocess. If set to ``none`` all galaxies in all outputs are postprocessed.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>countBatch</name>
      <defaultValue>1024_c_size_t</defaultValue>
      <description>
      The number of galaxies to read from file, and then postprocess, in each batch.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="mergerTreeOutputter"  name="mergerTreeOutputter_"  source="parameters"/>
    <objectBuilder class="starFormationHistory" name="starFormationHistory_" source="parameters"/>
    !!]
    if (associated(parametersRoot)) then
       self=taskPostprocessOutputs(fileName,outputsGroupName,selectionFileName,countBatch,mergerTreeOutputter_,starFormationHistory_,parametersRoot)
    else
       self=taskPostprocessOutputs(fileName,outputsGroupName,selectionFileName,countBatch,mergerTreeOutputter_,starFormationHistory_,parameters    )
    end if
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    allocate(self%nodeComponents_)
    allocate(self%nodeHierarchy_ )
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi?id=105807" docformat="rst">
      <description>
      ICE when passing a derived type component to a class(*) function argument.
      </description>
    !!]
    dummyPointer_              => self%nodeComponents_
    self%nodeComponentsManager =  resourceManager(dummyPointer_)
    dummyPointer_              => self%nodeHierarchy_
    self%nodeHierarchyManager  =  resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    !![
    <objectDestructor name="mergerTreeOutputter_" />
    <objectDestructor name="starFormationHistory_"/>
    !!]
    return
  end function postprocessOutputsConstructorParameters

  function postprocessOutputsConstructorInternal(fileName,outputsGroupName,selectionFileName,countBatch,mergerTreeOutputter_,starFormationHistory_,parameters) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`taskPostprocessOutputs` task class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (taskPostprocessOutputs   )                        :: self
    type   (varying_string           ), intent(in   )         :: fileName             , outputsGroupName, &
         &                                                       selectionFileName
    integer(c_size_t                 ), intent(in   )         :: countBatch
    type   (inputParameters          ), intent(in   ), target :: parameters
    class  (mergerTreeOutputterClass ), intent(in   ), target :: mergerTreeOutputter_
    class  (starFormationHistoryClass), intent(in   ), target :: starFormationHistory_
    !![
    <constructorAssign variables="fileName, outputsGroupName, selectionFileName, countBatch, *mergerTreeOutputter_, *starFormationHistory_, *parameters"/>
    !!]

    if (countBatch < 1_c_size_t) call Error_Report('`countBatch` must be positive'//{introspection:location})
    !![
    <addMetaProperty component="basic" name="inclination"                id="self%inclinationID" isEvolvable="no"   isCreator="yes"/>
    <addMetaProperty component="basic" name="nodeIndexPostprocessSource" id="self%indexSourceID" type="longInteger" isCreator="yes"/>
    !!]
    return
  end function postprocessOutputsConstructorInternal

  subroutine postprocessOutputsDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`taskPostprocessOutputs` task class.
    !!}
    implicit none
    type(taskPostprocessOutputs), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeOutputter_" />
    <objectDestructor name="self%starFormationHistory_"/>
    !!]
    return
  end subroutine postprocessOutputsDestructor

  subroutine postprocessOutputsPerform(self,status)
    !!{RST
    Postprocess galaxies read back from a Galacticus output file.
    !!}
    use    :: Abundances_Structure            , only : abundances
    use    :: Display                         , only : displayIndent                                   , displayMessage                     , displayUnindent
    use    :: Error                           , only : Error_Report                                    , errorStatusSuccess
    use    :: Galacticus_Nodes                , only : mergerTree                                      , nodeComponentBasic                 , nodeComponentDisk, nodeComponentSpheroid     , &
         &                                             treeNode
    use    :: HDF5                            , only : HSIZE_T
    use    :: HDF5_Access                     , only : hdf5Access
    use    :: Histories                       , only : history
    use    :: IO_HDF5                         , only : hdf5File                                        , hdf5Group                          , hdf5Dataset      , hdf5VarDouble             , &
         &                                             hdf5VarDouble2D
    use    :: ISO_Varying_String              , only : var_str                                         , operator(//)                       , char             , operator(==)
    use    :: Kind_Numbers                    , only : kind_int8
    use    :: Merger_Tree_Outputters          , only : enumerationOutputGroupTypeType                  , enumerationOutputGroupTypeDecode   , outputGroupTypeNode, outputGroupTypeLightcone, &
         &                                             outputGroupTypeSnapshot                         , outputGroupTypeTree
    use    :: Output_HDF5                     , only : outputTo                          => outputGroup
    use    :: Node_Components                 , only : Node_Components_Thread_Initialize               , Node_Components_Thread_Uninitialize
    use    :: Numerical_Constants_Astronomical, only : degreesToRadians
    use    :: String_Handling                 , only : operator(//)
    !$ use :: OMP_Lib                         , only : OMP_Get_Thread_Num
    implicit none
    class           (taskPostprocessOutputs        ), intent(inout), target                                 :: self
    integer                                         , intent(  out), optional                               :: status
    type            (inputParameters               ), allocatable  , save                                   :: parameters
    !$omp threadprivate(parameters)
    logical                                                        , save                                   :: recordParameters
    !$omp threadprivate(recordParameters)
    ! HDF5 objects are allocatable so that they can be released while HDF5 access is locked.
    type            (hdf5File                      ), allocatable                                           :: file              , fileSelection
    type            (hdf5Group                     ), allocatable                                           :: outputsGroup      , outputGroup   , &
         &                                                                                                     selectionsGroup   , selectionGroup, &
         &                                                                                                     outputsGroupTo    , outputGroupTo
    type            (hdf5Group                     ), allocatable  , dimension(:                          ) :: nodeDataGroups
    type            (hdf5Dataset                   ), allocatable                                           :: dataset
    type            (enumerationOutputGroupTypeType), allocatable  , dimension(:                          ) :: typesOutput
    type            (enumerationOutputGroupTypeType)               , dimension(4                          ) :: typesKnown
    type            (varying_string                )                                                        :: typeOutput
    double precision                                                                                        :: attributeValue
    type            (postprocessOutputsSelection   ), allocatable  , dimension(:                          ) :: selections
    type            (postprocessOutputsBatch       ), allocatable  , dimension(:                          ) :: batches           , batchesTmp
    double precision                                , allocatable  , dimension(:                          ) :: timesOutput
    character       (len=32                        )               , dimension(  postprocessColumnCount   ) :: namesColumn
    character       (len=48                        )               , dimension(2,postprocessComponentCount) :: namesHistory
    integer         (c_size_t                      )                                                        :: countOutputs      , indexOutput   , &
         &                                                                                                     countGalaxies     , countBatches  , &
         &                                                                                                     indexBatch        , indexBegin    , &
         &                                                                                                     countGalaxiesTotal
    integer                                                                                                 :: iColumn           , iComponent
    type            (varying_string                )                                                        :: nameOutput        , message
    character       (len=32                        )               , dimension(2                          ) :: prefixComponent
    character       (len=32                        )               , dimension(3                          ) :: namesAttributeTime
    integer                                                                                                 :: iType             , iAttribute

    call displayIndent('Begin task: output postprocess')
    ! Set status to success by default.
    if (present(status)) status=errorStatusSuccess
    ! Construct the names of the datasets from which properties are read.
    prefixComponent                          =['disk    ','spheroid']
    namesColumn(postprocessColumnTime       )= 'time'
    namesColumn(postprocessColumnMassBasic  )= 'basicMass'
    namesColumn(postprocessColumnInclination)= 'inclination'
    do iComponent=1,postprocessComponentCount
       iColumn=postprocessColumnCountNode+(iComponent-1)*postprocessColumnCountComponent
       namesColumn (iColumn+postprocessColumnMassStellar      )=trim(prefixComponent(iComponent))//'MassStellar'
       namesColumn (iColumn+postprocessColumnMassGas          )=trim(prefixComponent(iComponent))//'MassGas'
       namesColumn (iColumn+postprocessColumnMassMetalsStellar)=trim(prefixComponent(iComponent))//'AbundancesStellarMetals'
       namesColumn (iColumn+postprocessColumnMassMetalsGas    )=trim(prefixComponent(iComponent))//'AbundancesGasMetals'
       namesColumn (iColumn+postprocessColumnRadius           )=trim(prefixComponent(iComponent))//'Radius'
       namesColumn (iColumn+postprocessColumnAngularMomentum  )=trim(prefixComponent(iComponent))//'AngularMomentum'
       namesColumn (iColumn+postprocessColumnVelocity         )=trim(prefixComponent(iComponent))//'Velocity'
       namesHistory(1,iComponent                              )=trim(prefixComponent(iComponent))//'StarFormationHistoryMass'
       namesHistory(2,iComponent                              )=trim(prefixComponent(iComponent))//'StarFormationHistoryTimes'
    end do
    ! Open the file, find the outputs which it contains, and divide the galaxies in each into batches.
    countBatches      =0_c_size_t
    countGalaxiesTotal=0_c_size_t
    allocate(batches(64))
    typesKnown        =[outputGroupTypeNode,outputGroupTypeLightcone,outputGroupTypeSnapshot,outputGroupTypeTree]
    !$ call hdf5Access%set()
    allocate(file,fileSelection,outputsGroup,outputGroup,selectionsGroup,selectionGroup,dataset)
    file=hdf5File(char(self%fileName),readOnly=.true.)
    if (.not.file%hasGroup(char(self%outputsGroupName))) call Error_Report('file `'//self%fileName//'` has no `'//self%outputsGroupName//'` group'//{introspection:location})
    outputsGroup=file%openGroup(char(self%outputsGroupName))
    countOutputs=0_c_size_t
    nameOutput  =var_str('Output')//(countOutputs+1_c_size_t)
    do while (outputsGroup%hasGroup(char(nameOutput)))
       countOutputs=countOutputs+1_c_size_t
       nameOutput  =var_str('Output')//(countOutputs+1_c_size_t)
    end do
    if (countOutputs == 0_c_size_t) call Error_Report('no outputs were found in file `'//self%fileName//'`'//{introspection:location})
    allocate(nodeDataGroups(countOutputs))
    allocate(selections    (countOutputs))
    allocate(timesOutput   (countOutputs))
    allocate(typesOutput   (countOutputs))
    if (self%selectionFileName /= 'none') then
       fileSelection=hdf5File(char(self%selectionFileName),readOnly=.true.)
       if (.not.fileSelection%hasGroup(char(self%outputsGroupName))) call Error_Report('selection file `'//self%selectionFileName//'` has no `'//self%outputsGroupName//'` group'//{introspection:location})
       selectionsGroup=fileSelection%openGroup(char(self%outputsGroupName))
    end if
    do indexOutput=1_c_size_t,countOutputs
       nameOutput                 =var_str('Output')//indexOutput
       outputGroup                =outputsGroup%openGroup(char(nameOutput))
       nodeDataGroups(indexOutput)=outputGroup %openGroup('nodeData'      )
       timesOutput   (indexOutput)=-1.0d0
       if (outputGroup%hasAttribute('outputTime')) call outputGroup%readAttribute('outputTime',timesOutput(indexOutput))
       ! Find the type of this output, so that the output written is of the same type.
       typesOutput(indexOutput)=outputGroupTypeNode
       if (outputGroup%hasAttribute('outputType')) then
          call outputGroup%readAttribute('outputType',typeOutput)
          do iType=1,size(typesKnown)
             if (typeOutput == enumerationOutputGroupTypeDecode(typesKnown(iType),includePrefix=.false.)) typesOutput(indexOutput)=typesKnown(iType)
          end do
       end if
       ! Determine how many galaxies are present in this output.
       countGalaxies=0_c_size_t
       do iColumn=1,postprocessColumnCount
          if (nodeDataGroups(indexOutput)%hasDataset(trim(namesColumn(iColumn)))) then
             dataset      =nodeDataGroups(indexOutput)%openDataset(trim(namesColumn(iColumn)))
             countGalaxies=int(dataset%size(1),kind=c_size_t)
             exit
          end if
       end do
       ! Determine which galaxies are to be postprocessed.
       if (self%selectionFileName /= 'none') then
          if (selectionsGroup%hasGroup(char(nameOutput))) then
             selectionGroup=selectionsGroup%openGroup(char(nameOutput))
             if (selectionGroup%hasDataset('selection')) then
                selections(indexOutput)%isSelected=.true.
                call selectionGroup%readDataset('selection',selections(indexOutput)%indices)
                if (size(selections(indexOutput)%indices) > 0) then
                   if     (                                                       &
                        &   any(selections(indexOutput)%indices <  0_c_size_t   ) &
                        &  .or.                                                   &
                        &   any(selections(indexOutput)%indices >= countGalaxies) &
                        & ) call Error_Report('selection for `'//nameOutput//'` must lie within the range of galaxies present'//{introspection:location})
                end if
                countGalaxies=size(selections(indexOutput)%indices,kind=c_size_t)
             else
                countGalaxies=0_c_size_t
             end if
          else
             countGalaxies=0_c_size_t
          end if
       end if
       countGalaxiesTotal=countGalaxiesTotal+countGalaxies
       ! Divide the galaxies into batches.
       indexBegin=0_c_size_t
       do while (indexBegin < countGalaxies)
          countBatches=countBatches+1_c_size_t
          if (countBatches > size(batches)) then
             call move_alloc(batches,batchesTmp)
             allocate(batches(2*size(batchesTmp)))
             batches(1:size(batchesTmp))=batchesTmp
             deallocate(batchesTmp)
          end if
          batches(countBatches)%indexOutput=indexOutput
          batches(countBatches)%indexBegin =indexBegin
          batches(countBatches)%count      =min(self%countBatch,countGalaxies-indexBegin)
          indexBegin                       =indexBegin+batches(countBatches)%count
       end do
    end do
    !$ call hdf5Access%unset()
    message=var_str('postprocessing ')//countGalaxiesTotal//' galaxies from '//countOutputs//' outputs'
    call displayMessage(message)
    ! Begin parallel postprocessing of batches.
    !$omp parallel private(indexBatch)
    allocate(postprocessOutputsOutputter_,mold=self%mergerTreeOutputter_)
    !$omp critical(postprocessOutputsDeepCopy)
    !![
    <deepCopyReset variables="self%mergerTreeOutputter_"/>
    <deepCopy source="self%mergerTreeOutputter_" destination="postprocessOutputsOutputter_"/>
    <deepCopyFinalize variables="postprocessOutputsOutputter_"/>
    !!]
    !$omp end critical(postprocessOutputsDeepCopy)
    !$omp barrier
    ! Call routines to perform initialization which must occur for all threads if run in parallel. The master thread initializes
    ! first, and alone, against a copy of the parameters for which output is enabled, so that parameters read during thread
    ! initialization are recorded in the output file, and so that parameters which take their default values are marked as such
    ! irrespective of the number of threads (see the `taskPostprocessForests` class).
    recordParameters=.true.
    !$ recordParameters=OMP_Get_Thread_Num() == 0
    if (recordParameters) then
       allocate(parameters)
       parameters=inputParameters(self%parameters,noOutput=.false.)
       call Node_Components_Thread_Initialize(parameters)
    end if
    !$omp barrier
    if (.not.recordParameters) then
       allocate(parameters)
       parameters=inputParameters(self%parameters)
       call Node_Components_Thread_Initialize(parameters)
    end if
    !$omp barrier
    !$omp do schedule(dynamic)
    do indexBatch=1_c_size_t,countBatches
       call postprocessOutputsBatchProcess(indexBatch)
    end do
    !$omp end do
    ! Reduce outputs back into the original outputter object.
    call postprocessOutputsOutputter_%reduce(self%mergerTreeOutputter_)
    !![
    <objectDestructor name="postprocessOutputsOutputter_"/>
    !!]
    call Node_Components_Thread_Uninitialize()
    !$omp barrier
    !$omp critical(evolveForestReset)
    call parameters%reset()
    !$omp end critical(evolveForestReset)
    !$omp barrier
    deallocate(parameters)
    !$omp end parallel
    ! Finalize outputs.
    call self%mergerTreeOutputter_%finalize()
    ! The outputter records, as attributes of each output group, the time of the first node which it output to that group. Nodes
    ! are output in an order which depends on the scheduling of threads, and, for a lightcone, each has its own time. Replace these
    ! attributes with those of the output from which the nodes were read. This assumes that the outputter writes to a group of the
    ! same name as that from which the nodes were read.
    namesAttributeTime=[character(len=32) :: 'outputTime','outputExpansionFactor','outputComovingDistance']
    !$ call hdf5Access%set()
    allocate(outputsGroupTo,outputGroupTo)
    if (outputTo%hasGroup(char(self%outputsGroupName))) then
       outputsGroupTo=outputTo%openGroup(char(self%outputsGroupName))
       do indexOutput=1_c_size_t,countOutputs
          nameOutput=var_str('Output')//indexOutput
          if (.not.outputsGroupTo%hasGroup(char(nameOutput))) cycle
          outputGroup  =outputsGroup  %openGroup(char(nameOutput))
          outputGroupTo=outputsGroupTo%openGroup(char(nameOutput))
          do iAttribute=1,size(namesAttributeTime)
             if (.not.outputGroup%hasAttribute(trim(namesAttributeTime(iAttribute)))) cycle
             call outputGroup  %readAttribute (trim(namesAttributeTime(iAttribute)),attributeValue)
             call outputGroupTo%writeAttribute(attributeValue,trim(namesAttributeTime(iAttribute)))
          end do
       end do
    end if
    ! Close the files.
    deallocate(outputGroupTo,outputsGroupTo)
    deallocate(nodeDataGroups)
    deallocate(dataset,selectionGroup,selectionsGroup,outputGroup,outputsGroup,fileSelection,file)
    !$ call hdf5Access%unset()
    call displayUnindent('Done task: output postprocess')
    return

  contains

    subroutine postprocessOutputsBatchProcess(indexBatch)
      !!{RST
      Read, rebuild nodes for, and output, the galaxies of a batch. All working storage is local to this routine, and so is
      private to the thread which calls it.
      !!}
      implicit none
      integer         (c_size_t                ), intent(in   )                                       :: indexBatch
      type            (postprocessOutputsColumn)               , dimension(postprocessColumnCount   ) :: columns
      logical                                                  , dimension(postprocessComponentCount) :: historyIsPresent     , componentIsPresent
      type            (hdf5VarDouble2D         ), allocatable  , dimension(:                        ) :: massesHistoryDisk    , massesHistorySpheroid
      type            (hdf5VarDouble           ), allocatable  , dimension(:                        ) :: timesHistoryDisk     , timesHistorySpheroid
      integer         (c_size_t                ), allocatable  , dimension(:                        ) :: indicesSource
      integer         (kind_int8               ), allocatable  , dimension(:                        ) :: indicesNode
      double precision                                         , dimension(0                        ) :: timesNull
      double precision                                         , dimension(0,0                      ) :: massesNull
      type            (mergerTree              ), pointer                                             :: tree
      type            (treeNode                ), pointer                                             :: node
      class           (nodeComponentBasic      ), pointer                                             :: basic
      class           (nodeComponentDisk       ), pointer                                             :: disk
      class           (nodeComponentSpheroid   ), pointer                                             :: spheroid
      type            (abundances              )                                                      :: abundances_
      type            (history                 )                                                      :: history_
      logical                                                                                         :: indicesNodeArePresent
      integer         (c_size_t                )                                                      :: indexOutput          , i
      integer                                                                                         :: iColumn              , iComponent           , &
           &                                                                                             j

      ! Create a tree to act as host to the nodes which this thread builds.
      allocate(tree)
      tree%index             =  1_kind_int8
      tree%volumeWeight      =  1.0d0
      tree%initializedUntil  =  0.0d0
      tree%isTreeInitialized =  .false.
      tree%firstTree         => tree
      tree%nodeBase          => null()
      tree%event             => null()
      call tree%properties%initialize()
      indexOutput=batches(indexBatch)%indexOutput
      ! Find the indices of the galaxies in this batch.
      if (allocated(indicesSource)) deallocate(indicesSource)
      allocate(indicesSource(batches(indexBatch)%count))
      if (selections(indexOutput)%isSelected) then
         indicesSource=selections(indexOutput)%indices(batches(indexBatch)%indexBegin+1_c_size_t:batches(indexBatch)%indexBegin+batches(indexBatch)%count)
      else
         do i=1_c_size_t,batches(indexBatch)%count
            indicesSource(i)=batches(indexBatch)%indexBegin+i-1_c_size_t
         end do
      end if
      ! Read the properties of the galaxies in this batch.
      !$ call hdf5Access%set()
      do iColumn=1,postprocessColumnCount
         columns(iColumn)%isPresent=nodeDataGroups(indexOutput)%hasDataset(trim(namesColumn(iColumn)))
         if (columns(iColumn)%isPresent) then
            if (selections(indexOutput)%isSelected) then
               call nodeDataGroups(indexOutput)%readDataset(trim(namesColumn(iColumn)),columns(iColumn)%values,readSelection=int (indicesSource   +1_c_size_t ,kind=HSIZE_T)                                                        )
            else
               call nodeDataGroups(indexOutput)%readDataset(trim(namesColumn(iColumn)),columns(iColumn)%values,readBegin    =int([indicesSource(1)+1_c_size_t],kind=HSIZE_T),readCount=int([batches(indexBatch)%count],kind=HSIZE_T))
            end if
         end if
      end do
      indicesNodeArePresent=nodeDataGroups(indexOutput)%hasDataset('nodeIndex')
      if (indicesNodeArePresent) then
         if (selections(indexOutput)%isSelected) then
            call nodeDataGroups(indexOutput)%readDataset('nodeIndex',indicesNode,readSelection=int( indicesSource   +1_c_size_t ,kind=HSIZE_T)                                                        )
         else
            call nodeDataGroups(indexOutput)%readDataset('nodeIndex',indicesNode,readBegin    =int([indicesSource(1)+1_c_size_t],kind=HSIZE_T),readCount=int([batches(indexBatch)%count],kind=HSIZE_T))
         end if
      end if
      do iComponent=1,postprocessComponentCount
         historyIsPresent(iComponent)= nodeDataGroups(indexOutput)%hasDataset(trim(namesHistory(1,iComponent))) &
              &                       .and.                                                                     &
              &                        nodeDataGroups(indexOutput)%hasDataset(trim(namesHistory(2,iComponent)))
      end do
      if (historyIsPresent(postprocessComponentDisk    )) then
         if (selections(indexOutput)%isSelected) then
            call nodeDataGroups(indexOutput)%readDataset(trim(namesHistory(1,postprocessComponentDisk    )),massesHistoryDisk    ,readSelection=int( indicesSource   +1_c_size_t ,kind=HSIZE_T)                                                        )
            call nodeDataGroups(indexOutput)%readDataset(trim(namesHistory(2,postprocessComponentDisk    )),timesHistoryDisk     ,readSelection=int( indicesSource   +1_c_size_t ,kind=HSIZE_T)                                                        )
         else
            call nodeDataGroups(indexOutput)%readDataset(trim(namesHistory(1,postprocessComponentDisk    )),massesHistoryDisk    ,readBegin    =int([indicesSource(1)+1_c_size_t],kind=HSIZE_T),readCount=int([batches(indexBatch)%count],kind=HSIZE_T))
            call nodeDataGroups(indexOutput)%readDataset(trim(namesHistory(2,postprocessComponentDisk    )),timesHistoryDisk     ,readBegin    =int([indicesSource(1)+1_c_size_t],kind=HSIZE_T),readCount=int([batches(indexBatch)%count],kind=HSIZE_T))
         end if
      end if
      if (historyIsPresent(postprocessComponentSpheroid)) then
         if (selections(indexOutput)%isSelected) then
            call nodeDataGroups(indexOutput)%readDataset(trim(namesHistory(1,postprocessComponentSpheroid)),massesHistorySpheroid,readSelection=int( indicesSource   +1_c_size_t ,kind=HSIZE_T)                                                        )
            call nodeDataGroups(indexOutput)%readDataset(trim(namesHistory(2,postprocessComponentSpheroid)),timesHistorySpheroid ,readSelection=int( indicesSource   +1_c_size_t ,kind=HSIZE_T)                                                        )
         else
            call nodeDataGroups(indexOutput)%readDataset(trim(namesHistory(1,postprocessComponentSpheroid)),massesHistorySpheroid,readBegin    =int([indicesSource(1)+1_c_size_t],kind=HSIZE_T),readCount=int([batches(indexBatch)%count],kind=HSIZE_T))
            call nodeDataGroups(indexOutput)%readDataset(trim(namesHistory(2,postprocessComponentSpheroid)),timesHistorySpheroid ,readBegin    =int([indicesSource(1)+1_c_size_t],kind=HSIZE_T),readCount=int([batches(indexBatch)%count],kind=HSIZE_T))
         end if
      end if
      !$ call hdf5Access%unset()
      if (.not.columns(postprocessColumnTime)%isPresent) then
         if (typesOutput(indexOutput) == outputGroupTypeLightcone) &
              & call Error_Report('lightcone output has no `time` dataset---each galaxy in a lightcone has its own time, so the model must output it using a `nodePropertyExtractorTime` object'//{introspection:location})
         if (timesOutput(indexOutput) < 0.0d0) &
              & call Error_Report('output has neither a `time` dataset nor an `outputTime` attribute, so the times of galaxies are unknown'//{introspection:location})
      end if
      ! Determine which components are present.
      do iComponent=1,postprocessComponentCount
         iColumn                       =postprocessColumnCountNode+(iComponent-1)*postprocessColumnCountComponent
         componentIsPresent(iComponent)= any(columns(iColumn+1:iColumn+postprocessColumnCountComponent)%isPresent) &
              &                         .or.                                                                       &
              &                          historyIsPresent(iComponent)
      end do
      ! Build, and output, a node for each galaxy.
      do i=1_c_size_t,batches(indexBatch)%count
         node => treeNode(hostTree=tree)
         if (indicesNodeArePresent) call node%indexSet(indicesNode(i))
         ! Basic component.
         basic => node%basic(autoCreate=.true.)
         if (columns(postprocessColumnTime)%isPresent) then
            call basic%timeSet(columns(postprocessColumnTime)%values(i))
         else
            call basic%timeSet(timesOutput(indexOutput)                )
         end if
         if (columns(postprocessColumnMassBasic  )%isPresent) call basic%massSet                        (columns(postprocessColumnMassBasic       )%values(i))
         if (columns(postprocessColumnInclination)%isPresent) call basic%floatRank0MetaPropertySet      (self%inclinationID,columns(postprocessColumnInclination)%values(i)*degreesToRadians)
         call                                                      basic%longIntegerRank0MetaPropertySet(self%indexSourceID,int(indicesSource(i),kind=kind_int8))
         ! Disk component.
         if (componentIsPresent(postprocessComponentDisk)) then
            j    =  postprocessColumnCountNode+(postprocessComponentDisk-1)*postprocessColumnCountComponent
            disk => node%disk(autoCreate=.true.)
            if (columns(j+postprocessColumnMassStellar      )%isPresent) call disk%massStellarSet     (columns(j+postprocessColumnMassStellar    )%values(i))
            if (columns(j+postprocessColumnMassGas          )%isPresent) call disk%massGasSet         (columns(j+postprocessColumnMassGas        )%values(i))
            if (columns(j+postprocessColumnRadius           )%isPresent) call disk%radiusSet          (columns(j+postprocessColumnRadius         )%values(i))
            if (columns(j+postprocessColumnAngularMomentum  )%isPresent) call disk%angularMomentumSet (columns(j+postprocessColumnAngularMomentum)%values(i))
            if (columns(j+postprocessColumnVelocity         )%isPresent) call disk%velocitySet        (columns(j+postprocessColumnVelocity       )%values(i))
            if (columns(j+postprocessColumnMassMetalsStellar)%isPresent) then
               call abundances_%metallicitySet(columns(j+postprocessColumnMassMetalsStellar)%values(i))
               call disk%abundancesStellarSet(abundances_)
            end if
            if (columns(j+postprocessColumnMassMetalsGas    )%isPresent) then
               call abundances_%metallicitySet(columns(j+postprocessColumnMassMetalsGas    )%values(i))
               call disk%abundancesGasSet    (abundances_)
            end if
            if (historyIsPresent(postprocessComponentDisk)) then
               if (associated(timesHistoryDisk(i)%row) .and. associated(massesHistoryDisk(i)%row)) then
                  call self%starFormationHistory_%reconstruct(node,history_,timesHistoryDisk(i)%row,massesHistoryDisk(i)%row)
               else
                  call self%starFormationHistory_%reconstruct(node,history_,timesNull              ,massesNull              )
               end if
               call disk    %starFormationHistorySet(history_)
               call history_%destroy                (        )
            end if
         end if
         ! Spheroid component.
         if (componentIsPresent(postprocessComponentSpheroid)) then
            j        =  postprocessColumnCountNode+(postprocessComponentSpheroid-1)*postprocessColumnCountComponent
            spheroid => node%spheroid(autoCreate=.true.)
            if (columns(j+postprocessColumnMassStellar      )%isPresent) call spheroid%massStellarSet    (columns(j+postprocessColumnMassStellar    )%values(i))
            if (columns(j+postprocessColumnMassGas          )%isPresent) call spheroid%massGasSet        (columns(j+postprocessColumnMassGas        )%values(i))
            if (columns(j+postprocessColumnRadius           )%isPresent) call spheroid%radiusSet         (columns(j+postprocessColumnRadius         )%values(i))
            if (columns(j+postprocessColumnAngularMomentum  )%isPresent) call spheroid%angularMomentumSet(columns(j+postprocessColumnAngularMomentum)%values(i))
            if (columns(j+postprocessColumnVelocity         )%isPresent) call spheroid%velocitySet       (columns(j+postprocessColumnVelocity       )%values(i))
            if (columns(j+postprocessColumnMassMetalsStellar)%isPresent) then
               call abundances_%metallicitySet(columns(j+postprocessColumnMassMetalsStellar)%values(i))
               call spheroid%abundancesStellarSet(abundances_)
            end if
            if (columns(j+postprocessColumnMassMetalsGas    )%isPresent) then
               call abundances_%metallicitySet(columns(j+postprocessColumnMassMetalsGas    )%values(i))
               call spheroid%abundancesGasSet    (abundances_)
            end if
            if (historyIsPresent(postprocessComponentSpheroid)) then
               if (associated(timesHistorySpheroid(i)%row) .and. associated(massesHistorySpheroid(i)%row)) then
                  call self%starFormationHistory_%reconstruct(node,history_,timesHistorySpheroid(i)%row,massesHistorySpheroid(i)%row)
               else
                  call self%starFormationHistory_%reconstruct(node,history_,timesNull                  ,massesNull                  )
               end if
               call spheroid%starFormationHistorySet(history_)
               call history_%destroy                (        )
            end if
         end if
         call postprocessOutputsOutputter_%outputNode(node,indexOutput,typesOutput(indexOutput))
         call node%destroy()
         deallocate(node)
      end do
      call tree%destroy()
      deallocate(tree)
      return
    end subroutine postprocessOutputsBatchProcess

  end subroutine postprocessOutputsPerform
