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
  Implements a merger tree operator which exports merger trees to
  file.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Merger_Tree_Data_Structure, only : enumerationMergerTreeFormatType

  !![
  <mergerTreeOperator name="mergerTreeOperatorExport">
   <description>
    This operator will export merger trees to a file specified by the {\normalfont \ttfamily [outputFileName]} using the format
    specified by {\normalfont \ttfamily [exportFormat]}. Currently, node indices (plus host indices, which are assumed identical to
    the node indices), descendant indices, masses and redshifts are exported. Positions and velocities are exported if available. If
    {\normalfont \ttfamily IRATE}-format output is requested then ``snapshot'' numbers will be assigned to nodes based on the time
    at which they exist. This usually only makes sense if the nodes are defined on a time grid (i.e. if merger trees were extracted
    from an N-body simulation, or if trees were re-gridded onto such a time grid; see \refPhysics{mergerTreeOperatorRegridTimes}).
    Export happens during the merger tree pre-evolution phase.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorExport
     !!{
     A merger tree operator class which exports merger trees to file.
     !!}
     private
     class  (cosmologyParametersClass       ), pointer :: cosmologyParameters_      => null()
     class  (cosmologyFunctionsClass        ), pointer :: cosmologyFunctions_       => null()
     class  (cosmologicalMassVarianceClass  ), pointer :: cosmologicalMassVariance_ => null()
     type   (varying_string                 )          :: outputFileName
     type   (enumerationMergerTreeFormatType)          :: exportFormat
     logical                                           :: snapshotsRequired
   contains
     final     ::                         exportDestructor
     procedure :: operatePreEvolution  => exportOperatePreEvolution
  end type mergerTreeOperatorExport

  interface mergerTreeOperatorExport
     !!{
     Constructors for the export merger tree operator class.
     !!}
     module procedure exportConstructorParameters
     module procedure exportConstructorInternal
  end interface mergerTreeOperatorExport

contains

  function exportConstructorParameters(parameters) result(self)
    !!{
    Constructor for the export merger tree operator class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameter                   , inputParameters
    use :: Merger_Tree_Data_Structure, only : enumerationMergerTreeFormatEncode
    implicit none
    type (mergerTreeOperatorExport     )                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class(cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class(cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    type (varying_string               )                :: outputFileName            , exportFormat

    !![
    <inputParameter>
      <name>outputFileName</name>
      <source>parameters</source>
      <defaultValue>var_str('galacticusExportedTrees.hdf5')</defaultValue>
      <description>The name of the file to which merger trees should be exported.</description>
    </inputParameter>
    <inputParameter>
      <name>exportFormat</name>
      <source>parameters</source>
      <defaultValue>var_str('galacticus')</defaultValue>
      <description>The output format to use when exporting merger trees.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=exportConstructorInternal(char(outputFileName),enumerationMergerTreeFormatEncode(char(exportFormat),includesPrefix=.false.),cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function exportConstructorParameters

  function exportConstructorInternal(outputFileName,exportFormat,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the export merger tree operator class.
    !!}
    use :: Error                     , only : Error_Report
    use :: Merger_Tree_Data_Structure, only : enumerationMergerTreeFormatIsValid, mergerTreeFormatIrate
    implicit none
    type     (mergerTreeOperatorExport       )                        :: self
    character(len=*                          ), intent(in   )         :: outputFileName
    type     (enumerationMergerTreeFormatType), intent(in   )         :: exportFormat
    class    (cosmologyParametersClass       ), intent(in   ), target :: cosmologyParameters_
    class    (cosmologyFunctionsClass        ), intent(in   ), target :: cosmologyFunctions_
    class    (cosmologicalMassVarianceClass  ), intent(in   ), target :: cosmologicalMassVariance_
    !![
    <constructorAssign variables="outputFileName, exportFormat, *cosmologyParameters_, *cosmologyFunctions_, *cosmologicalMassVariance_"/>
    !!]

    ! Validate the export format.
    if (.not.enumerationMergerTreeFormatIsValid(exportFormat)) call Error_Report('exportFormat is invalid'//{introspection:location})
    ! Construct the object.
    self%snapshotsRequired=(exportFormat == mergerTreeFormatIrate)
    return
  end function exportConstructorInternal

  subroutine exportDestructor(self)
    !!{
    Destructor for the export merger tree operator function class.
    !!}
    implicit none
    type(mergerTreeOperatorExport), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine exportDestructor

  subroutine exportOperatePreEvolution(self,tree)
    !!{
    Output the structure of {\normalfont \ttfamily tree}.
    !!}
    use :: Cosmology_Parameters            , only : hubbleUnitsLittleH
    use :: Dates_and_Times                 , only : Formatted_Date_and_Time
    use :: Galacticus_Nodes                , only : defaultPositionComponent     , mergerTree            , nodeComponentBasic    , nodeComponentPosition      , &
          &                                         treeNode
    use :: Merger_Tree_Data_Structure      , only : mergerTreeData               , metaDataTypeCosmology , metaDataTypeProvenance, propertyTypeDescendantIndex, &
          &                                         propertyTypeHostIndex        , propertyTypeNodeIndex , propertyTypeNodeMass  , propertyTypePositionX      , &
          &                                         propertyTypePositionY        , propertyTypePositionZ , propertyTypeRedshift  , propertyTypeSnapshot       , &
          &                                         propertyTypeTreeIndex        , propertyTypeTreeWeight, propertyTypeVelocityX , propertyTypeVelocityY      , &
          &                                         propertyTypeVelocityZ        , unitsLength           , unitsMass             , unitsVelocity
    use :: Merger_Tree_Walkers             , only : mergerTreeWalkerIsolatedNodes
    use :: Numerical_Constants_Astronomical, only : massSolar                    , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Interpolation         , only : interpolator
    use :: Sorting                         , only : sort
    implicit none
    class           (mergerTreeOperatorExport     ), intent(inout), target         :: self
    type            (mergerTree                   ), intent(inout), target         :: tree
    integer         (kind=size_t                  ), parameter                     :: hdfChunkSize                   =1024
    integer                                        , parameter                     :: hdfCompressionLevel            =   9
    double precision                               , allocatable  , dimension(:  ) :: nodeMass                            , nodeRedshift         , &
         &                                                                            snapshotTime                        , snapshotTimeTemp     , &
         &                                                                            treeWeight
    double precision                               , allocatable  , dimension(:,:) :: nodePosition                        , nodeVelocity
    integer         (kind=kind_int8               ), allocatable  , dimension(:  ) :: descendantIndex                     , nodeIndex            , &
         &                                                                            nodeSnapshot                        , treeIndex
    type            (treeNode                     ), pointer                       :: node
    class           (nodeComponentBasic           ), pointer                       :: basic
    class           (nodeComponentPosition        ), pointer                       :: position
    type            (mergerTree                   ), pointer                       :: treeCurrent
    integer                                        , parameter                     :: snapshotCountIncrement         = 100
    type            (interpolator                 ), allocatable                   :: snapshotInterpolator
    type            (mergerTreeWalkerIsolatedNodes)                                :: treeWalker
    integer                                                                        :: nodeCount                           , snapshotCount
    type            (mergerTreeData               )                                :: mergerTrees

    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Reset the tree data structure.
       call mergerTrees%reset                   (       )
       ! Specify a single tree in the structure.
       call mergerTrees%forestCountSet          (      1)
       ! Specify that we do not want to create individual merger tree reference datasets.
       call mergerTrees%makeReferences          (.false.)
       ! Specify that trees are self-contained (i.e. nodes never move from one tree to another).
       call mergerTrees%setSelfContained        (.true. )
       ! Specify that halo masses do include subhalo contributions.
       call mergerTrees%setIncludesSubhaloMasses(.true. )
       ! Specify units system used.
       call mergerTrees%setUnits(unitsMass    ,unitsInSI=massSolar ,hubbleExponent=0,scaleFactorExponent=0,name="Msolar")
       call mergerTrees%setUnits(unitsLength  ,unitsInSI=megaParsec,hubbleExponent=0,scaleFactorExponent=0,name="Mpc"   )
       call mergerTrees%setUnits(unitsVelocity,unitsInSI=kilo      ,hubbleExponent=0,scaleFactorExponent=0,name="km/s"  )
       ! Set cosmology metadata.
       call mergerTrees%addMetadata(metaDataTypeCosmology ,'OmegaMatter'       ,self%cosmologyParameters_     %OmegaMatter    (                  ))
       call mergerTrees%addMetadata(metaDataTypeCosmology ,'OmegaBaryon'       ,self%cosmologyParameters_     %OmegaBaryon    (                  ))
       call mergerTrees%addMetadata(metaDataTypeCosmology ,'OmegaLambda'       ,self%cosmologyParameters_     %OmegaDarkEnergy(                  ))
       call mergerTrees%addMetadata(metaDataTypeCosmology ,'HubbleParam'       ,self%cosmologyParameters_     %HubbleConstant (hubbleUnitsLittleH))
       call mergerTrees%addMetadata(metaDataTypeCosmology ,'sigma_8'           ,self%cosmologicalMassVariance_%sigma8         (                  ))
       call mergerTrees%addMetadata(metaDataTypeCosmology ,'powerSpectrumIndex',"not specified"                                                   )
       ! Set provenance metadata.
       call mergerTrees%addMetadata(metaDataTypeProvenance,'fileBuiltBy'       ,'Galacticus'                   )
       call mergerTrees%addMetadata(metaDataTypeProvenance,'fileTimestamp'     ,char(Formatted_Date_and_Time()))
       ! Count nodes in the tree.
       nodeCount=0
       treeWalker=mergerTreeWalkerIsolatedNodes(treeCurrent)
       do while (treeWalker%next(node))
          nodeCount=nodeCount+1
       end do
       call mergerTrees%nodeCountSet(nodeCount)
       ! Allocate arrays for serialization.
       allocate(treeIndex      (nodeCount))
       allocate(treeWeight     (nodeCount))
       allocate(nodeIndex      (nodeCount))
       allocate(descendantIndex(nodeCount))
       allocate(nodeMass       (nodeCount))
       allocate(nodeRedshift   (nodeCount))
       if (self%snapshotsRequired                       ) allocate(nodeSnapshot(nodeCount  ))
       if (defaultPositionComponent%positionIsGettable()) then
          allocate(nodePosition(nodeCount,3))
       else
          allocate(nodePosition(        0,0))
       end if
       if (defaultPositionComponent%velocityIsGettable()) then
          allocate(nodeVelocity(nodeCount,3))
       else
          allocate(nodeVelocity(        0,0))
       end if
       ! Find "snapshot" numbers for nodes - relevant only for IRATE output format.
       if (self%snapshotsRequired) then
          allocate(snapshotTime(snapshotCountIncrement))
          node                        => treeCurrent%nodeBase
          basic                       => node       %basic   ()
          snapshotCount               =  1
          snapshotTime(snapshotCount) =  basic      %time    ()
          treeWalker                  =  mergerTreeWalkerIsolatedNodes(treeCurrent)
          do while (treeWalker%next(node))
             basic => node%basic()
             if (all(snapshotTime(1:snapshotCount) /= basic%time())) then
                snapshotCount=snapshotCount+1
                if (snapshotCount > size(snapshotTime)) then
                   call Move_Alloc(snapshotTime,snapshotTimeTemp)
                   allocate(snapshotTime(size(snapshotTimeTemp)+snapshotCountIncrement))
                   snapshotTime(1:size(snapshotTimeTemp))=snapshotTimeTemp
                   deallocate(snapshotTimeTemp)
                end if
                snapshotTime(snapshotCount)=basic%time()
             end if
          end do
          call sort(snapshotTime(1:snapshotCount))
          if (allocated(snapshotInterpolator)) deallocate(snapshotInterpolator)
          allocate(snapshotInterpolator)
          snapshotInterpolator=interpolator(snapshotTime(1:snapshotCount))
          deallocate(snapshotTime)
       else
          snapshotCount=0
       end if
       ! Serialize node data to arrays and write to merger tree data structure.
       treeIndex                =                              treeCurrent%index
       treeWeight               =                              treeCurrent%volumeWeight
       treeWalker               =mergerTreeWalkerIsolatedNodes(treeCurrent             )
       nodeCount                =0
       do while (treeWalker%next(node))
          nodeCount=nodeCount+1
          nodeIndex      (nodeCount)=  node       %index   ()
          descendantIndex(nodeCount)=  node%parent%index   ()
          basic                     => node       %basic   ()
          position                  => node       %position()
          nodeMass       (nodeCount)=                                                                                              basic%mass()
          nodeRedshift   (nodeCount)=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(basic%time()))
          if (defaultPositionComponent%positionIsGettable()) nodePosition(nodeCount,:)=                            position%position()
          if (defaultPositionComponent%velocityIsGettable()) nodeVelocity(nodeCount,:)=                            position%velocity()
          if (self                    %snapshotsRequired   ) nodeSnapshot(nodeCount  )=snapshotInterpolator%locate(basic   %time    (),closest=.true.)
       end do
       call mergerTrees%setProperty(propertyTypeTreeWeight     ,treeWeight     )
       call mergerTrees%setProperty(propertyTypeTreeIndex      ,treeIndex      )
       call mergerTrees%setProperty(propertyTypeNodeIndex      ,nodeIndex      )
       call mergerTrees%setProperty(propertyTypeHostIndex      ,nodeIndex      )
       call mergerTrees%setProperty(propertyTypeDescendantIndex,descendantIndex)
       call mergerTrees%setProperty(propertyTypeNodeMass       ,nodeMass       )
       call mergerTrees%setProperty(propertyTypeRedshift       ,nodeRedshift   )
       if (defaultPositionComponent%positionIsGettable()) then
          call mergerTrees%setProperty(propertyTypePositionX,nodePosition(:,1))
          call mergerTrees%setProperty(propertyTypePositionY,nodePosition(:,2))
          call mergerTrees%setProperty(propertyTypePositionZ,nodePosition(:,3))
       end if
       if (defaultPositionComponent%velocityIsGettable()) then
          call mergerTrees%setProperty(propertyTypeVelocityX,nodeVelocity(:,1))
          call mergerTrees%setProperty(propertyTypeVelocityY,nodeVelocity(:,2))
          call mergerTrees%setProperty(propertyTypeVelocityZ,nodeVelocity(:,3))
       end if
       if (self%snapshotsRequired) call mergerTrees%setProperty(propertyTypeSnapshot,nodeSnapshot)
       ! Write the tree to file.
       call mergerTrees%export(char(self%outputFileName),self%exportFormat,hdfChunkSize,hdfCompressionLevel,append=.true.)
       ! Deallocate arrays.
       deallocate(treeIndex      )
       deallocate(treeWeight     )
       deallocate(nodeIndex      )
       deallocate(descendantIndex)
       deallocate(nodeMass       )
       deallocate(nodeRedshift   )
       if (self%snapshotsRequired                       ) deallocate(nodeSnapshot)
       if (defaultPositionComponent%positionIsGettable()) deallocate(nodePosition)
       if (defaultPositionComponent%velocityIsGettable()) deallocate(nodeVelocity)
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine exportOperatePreEvolution
