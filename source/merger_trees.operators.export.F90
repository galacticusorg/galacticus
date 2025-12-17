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

    Node positions and velocities will be exported if they are available.
    
    If {\normalfont \ttfamily [skipSingleNodeTrees]} is true then trees consisting of a single node (which are typically left over
    from pruning operators, and which are effectively inert) are not output. If {\normalfont \ttfamily [includeScaleRadii]} is
    true then scale radii will be included in the output (if they are available and have been set).  If {\normalfont \ttfamily
    [includeAngularMomenta]} is true then halo angular momenta will be included in the output (if they are available and have been
    set).

    If {\normalfont \ttfamily [includeSubhalos]} is true then subhalos are included in the exported data. Note that if particle
    positions were used to track subhalos after their destruction in an N-body simulation, attempting to output subhalos may lead
    to errors.
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
     logical                                           :: snapshotsRequired                  , includeScaleRadii, &
          &                                               includeAngularMomenta              , includeSubhalos  , &
          &                                               skipSingleNodeTrees
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
    type   (mergerTreeOperatorExport     )                :: self
    type   (inputParameters              ), intent(inout) :: parameters
    class  (cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class  (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class  (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    type   (varying_string               )                :: outputFileName           , exportFormat
    logical                                               :: includeAngularMomenta    , includeScaleRadii  , &
         &                                                   includeSubhalos          , skipSingleNodeTrees

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
    <inputParameter>
      <name>skipSingleNodeTrees</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, trees consisting of a single node (which are inert in Galacticus) are not output.</description>
    </inputParameter>
    <inputParameter>
      <name>includeScaleRadii</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, include scale radii (if available) in the output.</description>
    </inputParameter>
    <inputParameter>
      <name>includeAngularMomenta</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, include angular momenta (if available) in the output.</description>
    </inputParameter>
    <inputParameter>
      <name>includeSubhalos</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, subhalos are included in the exported data. (Note that subhalo export is not supported in cases where particle data was used to track halos after they are destroyed in an N-body simulation.)</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=exportConstructorInternal(char(outputFileName),enumerationMergerTreeFormatEncode(char(exportFormat),includesPrefix=.false.),skipSingleNodeTrees,includeScaleRadii,includeAngularMomenta,includeSubhalos,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function exportConstructorParameters

  function exportConstructorInternal(outputFileName,exportFormat,skipSingleNodeTrees,includeScaleRadii,includeAngularMomenta,includeSubhalos,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the export merger tree operator class.
    !!}
    use :: Error                     , only : Error_Report
    use :: Merger_Tree_Data_Structure, only : enumerationMergerTreeFormatIsValid, mergerTreeFormatIrate
    implicit none
    type     (mergerTreeOperatorExport       )                        :: self
    character(len=*                          ), intent(in   )         :: outputFileName
    type     (enumerationMergerTreeFormatType), intent(in   )         :: exportFormat
    logical                                   , intent(in   )         :: includeScaleRadii        , includeAngularMomenta, &
         &                                                               includeSubhalos          , skipSingleNodeTrees
    class    (cosmologyParametersClass       ), intent(in   ), target :: cosmologyParameters_
    class    (cosmologyFunctionsClass        ), intent(in   ), target :: cosmologyFunctions_
    class    (cosmologicalMassVarianceClass  ), intent(in   ), target :: cosmologicalMassVariance_
    !![
    <constructorAssign variables="outputFileName, exportFormat, skipSingleNodeTrees, includeScaleRadii, includeAngularMomenta, includeSubhalos, *cosmologyParameters_, *cosmologyFunctions_, *cosmologicalMassVariance_"/>
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
    use :: Galacticus_Nodes                , only : defaultPositionComponent         , mergerTree                  , nodeComponentBasic         , nodeComponentPosition       , &
         &                                          defaultDarkMatterProfileComponent, defaultSpinComponent        , treeNode                   , nodeComponentSpin           , &
         &                                          nodeComponentDarkMatterProfile   , nodeComponentSatellite      , nodeEvent                  , nodeEventBranchJump         , &
         &                                          nodeEventSubhaloPromotion        , treeNodeList
    use :: Histories                       , only : history                          , longIntegerHistory
    use :: Merger_Tree_Data_Structure      , only : mergerTreeData                   , metaDataTypeCosmology       , metaDataTypeProvenance     , propertyTypeDescendantIndex , &
          &                                         propertyTypeHostIndex            , propertyTypeNodeIndex       , propertyTypeNodeMass       , propertyTypePositionX       , &
          &                                         propertyTypePositionY            , propertyTypePositionZ       , propertyTypeRedshift       , propertyTypeSnapshot        , &
          &                                         propertyTypeTreeIndex            , propertyTypeTreeWeight      , propertyTypeVelocityX      , propertyTypeVelocityY       , &
          &                                         propertyTypeVelocityZ            , propertyTypeScaleRadius     , propertyTypeAngularMomentum, propertyTypeAngularMomentumX, &
          &                                         propertyTypeAngularMomentumY     , propertyTypeAngularMomentumZ, unitsLength                , unitsMass                   , &
          &                                         unitsVelocity
    use :: Merger_Tree_Walkers             , only : mergerTreeWalkerAllNodes
    use :: Numerical_Constants_Astronomical, only : massSolar                        , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Interpolation         , only : interpolator
    use :: Sorting                         , only : sort                             , sortIndex
    implicit none
    class           (mergerTreeOperatorExport      ), intent(inout), target         :: self
    type            (mergerTree                    ), intent(inout), target         :: tree
    integer         (size_t                        ), parameter                     :: hdfChunkSize             =1024
    integer                                         , parameter                     :: hdfCompressionLevel      =   9
    double precision                                , allocatable  , dimension(:  ) :: nodeMass                      , nodeRedshift         , &
         &                                                                             snapshotTime                  , snapshotTimeTemp     , &
         &                                                                             treeWeight                    , nodeAngularMomentum  , &
         &                                                                             nodeRadiusScale
    double precision                                , allocatable  , dimension(:,:) :: nodePosition                  , nodeVelocity         , &
         &                                                                             nodeAngularMomentum3D
    integer         (kind_int8                     ), allocatable  , dimension(:  ) :: descendantIndex               , nodeIndex            , &
         &                                                                             nodeSnapshot                  , treeIndex            , &
         &                                                                             hostIndex
    type            (treeNode                      ), pointer                       :: node                          , nodeHost 
    class           (nodeComponentBasic            ), pointer                       :: basic                         , basicTarget
    class           (nodeComponentSatellite        ), pointer                       :: satellite                     , satelliteTarget
    class           (nodeComponentPosition         ), pointer                       :: position
    class           (nodeComponentDarkMatterProfile), pointer                       :: darkMatterProfile
    class           (nodeComponentSpin             ), pointer                       :: spin
    type            (mergerTree                    ), pointer                       :: treeCurrent
    class           (nodeEvent                     ), pointer                       :: event
    integer                                         , parameter                     :: snapshotCountIncrement   = 100
    type            (interpolator                  ), allocatable                   :: snapshotInterpolator
    type            (treeNodeList                  ), allocatable  , dimension(:  ) :: nodeEvent_
    double precision                                , allocatable  , dimension(:  ) :: timeEvent
    integer         (c_size_t                      ), allocatable  , dimension(:  ) :: orderEvent
    integer                                                                         :: snapshotCount                 , countEvent
    integer         (kind_int8                     )                                :: nodeIndexMaximum              , nodeIndexOffset      , &
         &                                                                             nodeIndexClone                , descendantIndexEvent
    integer         (c_size_t                      )                                :: indexHistory                  , indexHistoryStart    , &
         &                                                                             nodeCount                     , countHistory
    type            (history                       )                                :: historyMassBound              , historyPosition
    type            (longIntegerHistory             )                               :: historyNodeIndex              , historyNodeIndexTarget
    type            (mergerTreeWalkerAllNodes      )                                :: treeWalker
    type            (mergerTreeData                )                                :: mergerTrees
    logical                                                                         :: positionIsGettable            , velocityIsGettable   , &
         &                                                                             angularMomentumIsGettable     , radiusScaleIsGettable, &
         &                                                                             haveAngularMomentum           , haveRadiusScale      , &
         &                                                                             havePosition                  , haveVelocity         , &
         &                                                                             angularMomentumIsVector       , isClone              , &
         &                                                                             haveMassHistory               , havePositionHistory  , &
         &                                                                             haveNodeIndexHistory          , haveSubhaloPromotion , &
         &                                                                             haveBranchJump
    double precision                                                                :: timeHistory
    
    ! Determine which properties to output.
    positionIsGettable       =defaultPositionComponent         %             positionIsGettable()
    velocityIsGettable       =defaultPositionComponent         %             velocityIsGettable()
    angularMomentumIsGettable=defaultSpinComponent             %      angularMomentumIsGettable() .and. self%includeAngularMomenta
    radiusScaleIsGettable    =defaultDarkMatterProfileComponent%                scaleIsGettable() .and. self%includeScaleRadii
    angularMomentumIsVector  =defaultSpinComponent             %angularMomentumVectorIsGettable()
    havePosition             =positionIsGettable
    haveVelocity             =velocityIsGettable
    haveAngularMomentum      =.false.
    haveRadiusScale          =.false.
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
       ! Count nodes in the tree, and determine if angular momenta and scale radii are available (if we need them).
       nodeCount       =0_c_size_t
       nodeIndexMaximum=0_kind_int8
       treeWalker      =mergerTreeWalkerAllNodes(treeCurrent)
       do while (treeWalker%next(node))
          nodeCount       =nodeCount+1_c_size_t
          nodeIndexMaximum=max(nodeIndexMaximum,node%index())
          ! Determine if we have mass, position, or index history information. If so, use them to determine the number of
          ! descendant subhalos that must be created.
          basic     => node%basic    ()
          satellite => node%satellite()
          position  => node%position ()
          if (self%includeSubhalos) then
             if (satellite%boundMassHistoryIsGettable()) then
                historyMassBound=satellite%boundMassHistory()
                if (allocated(historyMassBound%time)) then
                   if (historyMassBound%time(1) == basic%time()) then
                      nodeCount=nodeCount+size(historyMassBound%time)-1_c_size_t
                   else
                      nodeCount=nodeCount+size(historyMassBound%time)
                   end if
                end if
             else if (satellite%nodeIndexHistoryIsGettable()) then
                historyNodeIndex=satellite%nodeIndexHistory()
                if (allocated(historyNodeIndex%time)) then
                   if (historyNodeIndex%time(1) == basic%time()) then
                      nodeCount=nodeCount+size(historyNodeIndex%time)-1_c_size_t
                   else
                      nodeCount=nodeCount+size(historyNodeIndex%time)
                   end if
                end if
             else if (position%positionHistoryIsGettable()) then
                historyPosition=position%positionHistory()
                if (allocated(historyPosition%time)) then
                   if (historyPosition%time(1) == basic%time()) then
                      nodeCount=nodeCount+size(historyPosition%time)-1_c_size_t
                   else
                      nodeCount=nodeCount+size(historyPosition%time)
                   end if
                end if
             end if
          end if
          ! Determine if we have angular momentum and/or scale information.
          if (angularMomentumIsGettable .and. .not.haveAngularMomentum) then
             spin              => node%spin             ()
             if (spin             %angularMomentum() /= 0.0d0) haveAngularMomentum=.true.
          end if
          if (radiusScaleIsGettable     .and. .not.haveRadiusScale    ) then
             darkMatterProfile => node%darkMatterProfile()
             if (darkMatterProfile%scale          () /= 0.0d0) haveRadiusScale    =.true.
          end if
          ! Check for clones.
          if (node%parent%index() == node%index()) then
             ! Check that this clone has no children.
             if (associated(node%firstChild)) call Error_Report('cloned nodes with children are not supported'//{introspection:location})
          end if
       end do
       call mergerTrees%nodeCountSet(int(nodeCount))
       ! Determine an offset for clone node indices.
       nodeIndexOffset=10_kind_int8**ceiling(log10(dble(nodeIndexMaximum)))       
       ! Allocate arrays for serialization.
       allocate       (treeIndex            (nodeCount  ))
       allocate       (treeWeight           (nodeCount  ))
       allocate       (nodeIndex            (nodeCount  ))
       allocate       (hostIndex            (nodeCount  ))
       allocate       (descendantIndex      (nodeCount  ))
       allocate       (nodeMass             (nodeCount  ))
       allocate       (nodeRedshift         (nodeCount  ))
       if (self%snapshotsRequired) &
            & allocate(nodeSnapshot(nodeCount  ))
       if (     havePosition       ) then
          allocate    (nodePosition         (nodeCount,3))
       else
          allocate    (nodePosition         (        0,0))
       end if
       if (     haveVelocity       ) then
          allocate    (nodeVelocity         (nodeCount,3))
       else
          allocate    (nodeVelocity         (        0,0))
       end if
       if (     haveAngularMomentum) then
          if (angularMomentumIsVector) then
             allocate (nodeAngularMomentum3D(nodeCount,3))
             allocate (nodeAngularMomentum  (        0  ))
          else
             allocate (nodeAngularMomentum3D(        0,0))
             allocate (nodeAngularMomentum  (nodeCount  ))
          end if
       else
          allocate    (nodeAngularMomentum3D(        0,0))
          allocate    (nodeAngularMomentum  (        0  ))
       end if
       if (     haveRadiusScale    ) then
          allocate    (nodeRadiusScale      (nodeCount  ))
       else
          allocate    (nodeRadiusScale      (        0  ))
       end if
       ! Find "snapshot" numbers for nodes - relevant only for IRATE output format.
       if (self%snapshotsRequired) then
          allocate(snapshotTime(snapshotCountIncrement))
          node                        => treeCurrent%nodeBase
          snapshotCount               =  1
          snapshotTime(snapshotCount) =  basic      %time    ()
          treeWalker                  =  mergerTreeWalkerAllNodes(treeCurrent)
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
       treeIndex =                         treeCurrent%index
       treeWeight=                         treeCurrent%volumeWeight
       treeWalker=mergerTreeWalkerAllNodes(treeCurrent             )
       nodeCount =0_c_size_t
       do while (treeWalker%next(node))
          nodeCount=nodeCount+1_c_size_t
          ! Determine if this node is a clone, or has a clone child.
          isClone       = node%index          () == node%parent%index()
          nodeIndexClone=+node%index          ()                        &
               &         +     nodeIndexOffset
          ! Store indices and properties.
          if (isClone) then
             nodeIndex      (nodeCount) =              nodeIndexClone
          else
             nodeIndex      (nodeCount) =  node       %index            ()
          end if
          if (node%isSatellite()) then
             hostIndex      (nodeCount) =  node%parent%index            ()
             descendantIndex(nodeCount) =  -1_kind_int8
          else
             hostIndex      (nodeCount) =  node       %index            ()
             descendantIndex(nodeCount) =  node%parent%index            ()
          end if
          basic                         => node       %basic            ()
          position                      => node       %position         ()
          darkMatterProfile             => node       %darkMatterProfile()
          spin                          => node       %spin             ()
          satellite                     => node       %satellite        ()
          nodeMass          (nodeCount) =                                                                                               basic%mass()
          nodeRedshift      (nodeCount) = self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(basic%time()))
          if (     havePosition       ) nodePosition   (nodeCount,:)=                            position         %position()
          if (     haveVelocity       ) nodeVelocity   (nodeCount,:)=                            position         %velocity()
          if (     haveRadiusScale    ) nodeRadiusScale(nodeCount  )=                            darkMatterProfile%scale   ()
          if (     haveAngularMomentum) then
             if (angularMomentumIsVector) then
                nodeAngularMomentum3D                  (nodeCount,:)=                            spin%angularMomentumVector   ()
             else
                nodeAngularMomentum                    (nodeCount  )=                            spin%angularMomentum   ()
             end if
          end if
          if (self%snapshotsRequired  ) nodeSnapshot   (nodeCount  )=snapshotInterpolator%locate(basic            %time    (),closest=.true.)
          ! Handle subhalo histories (mass, position and velocity). Note that we assume that scale radii and angular momenta do
          ! not change, as Galacticus currently ignores these for subhalo evolution.
          if (self%includeSubhalos) then
             ! First, determine which histories we have available.
             haveMassHistory=.false.
             if (satellite%boundMassHistoryIsGettable()) then
                historyMassBound=satellite%boundMassHistory()
                if (allocated(historyMassBound%time)) haveMassHistory     =.true.
             end if
             haveNodeIndexHistory=.false.
             if (satellite%nodeIndexHistoryIsGettable()) then
                historyNodeIndex=satellite%nodeIndexHistory()
                if (allocated(historyNodeIndex%time)) haveNodeIndexHistory=.true.
             end if
             havePositionHistory=.false.
             if (position % positionHistoryIsGettable()) then
                historyPosition=position%positionHistory()
                if (allocated(historyPosition %time)) havePositionHistory =.true.
             end if
             ! Determine what events (if any) are associated with this subhalo.
             descendantIndexEvent=-1_kind_int8
             if (allocated(timeEvent )) deallocate(timeEvent )
             if (allocated(nodeEvent_)) deallocate(nodeEvent_)
             if (allocated(orderEvent)) deallocate(orderEvent)
             countEvent          =0
             haveSubhaloPromotion=.false.
             haveBranchJump      =.false.
             ! First count the number of branch jumps.
             event => node%event
             do while (associated(event))
                if (associated(event%task)) then
                   select type (event)
                   type is (nodeEventBranchJump)
                      countEvent=countEvent+1
                   end select
                end if
                event => event%next
             end do
             ! Build a list of branch jump events, and record if any subhalo promotion occurs.
             allocate(timeEvent (countEvent))
             allocate(nodeEvent_(countEvent))
             allocate(orderEvent(countEvent))
             countEvent =  0
             event      => node%event
             do while (associated(event))
                if (associated(event%task)) then
                   select type (event)
                   type is (nodeEventBranchJump)
                      ! Branch jump - record the node jumped to, and the time of the jump.
                      countEvent    =countEvent+1
                      haveBranchJump=.true.
                      nodeEvent_(countEvent)%node => event%node
                      timeEvent (countEvent)      =  event%time
                      orderEvent(countEvent)      =  countEvent
                   type is (nodeEventSubhaloPromotion)
                      ! Subhalo promotion - simply set the descendant to the target node.
                      haveSubhaloPromotion=.true.
                      descendantIndexEvent=event%node%index()
                      class default
                      call Error_Report('unknown event type'//{introspection:location})
                   end select
                end if
                event => event%next
             end do
             ! Determine the time ordering of branch jumps.
             orderEvent=sortIndex(timeEvent)
             ! If histories are available, initialize the host node.
             if (haveMassHistory .or. haveNodeIndexHistory .or. havePositionHistory) then
                countEvent=1
                if (node%isSatellite()) then
                   ! The node is an initial satellite.
                   if (node%parent%parent%index() == node%parent%index()) then
                      ! Cloned parent - first move to the cloned parent, then the host of our first step in the history is that
                      ! parent's parent.
                      nodeHost => node%parent%parent%parent
                   else
                      ! Non-cloned parent - the host of our first step in the history is that parent's parent.
                      nodeHost => node%parent%parent
                   end if
                   ! Check for an immediate branch jump.
                   if (haveBranchJump .and. timeEvent(orderEvent(countEvent)) == basic%time()) then
                      ! Branch jump occurs at the initial time. Move the host pointer to the jumped-to branch.
                      nodeHost   => nodeEvent_(orderEvent(countEvent))%node
                      countEvent =  countEvent+1
                      if (countEvent > size(timeEvent)) haveBranchJump=.false.
                   end if
                else
                   ! Node is not an initial satellite. The host of our first step in the history is simply the node's parent.
                   nodeHost => node%parent
                end if
                ! Determine the initial time in the subhalo history, and the number of entries in that history.
                if (haveMassHistory) then
                   timeHistory =     historyMassBound%time(1)
                   countHistory=size(historyMassBound%time   )
                else if (haveNodeIndexHistory) then
                   timeHistory =     historyNodeIndex%time(1)
                   countHistory=size(historyNodeIndex%time   )
                else
                   timeHistory =     historyPosition %time(1)
                   countHistory=size(historyPosition %time   )
                end if
                ! Determine the starting point in the history.
                if (timeHistory == basic%time()) then
                   ! Initial satellite - skip the first entry in the history.
                   indexHistoryStart=2_c_size_t
                else
                   ! Not an initial satellite - begin from the first entry in the history.
                   indexHistoryStart=1_c_size_t
                end if
                ! Iterate over all history entries.
                do indexHistory=1_c_size_t,countHistory+1_c_size_t-indexHistoryStart
                   ! Determine the time at this step in this history.
                   if (haveMassHistory) then
                      timeHistory=historyMassBound%time(indexHistory-1_c_size_t+indexHistoryStart)
                   else if (haveNodeIndexHistory) then
                      timeHistory=historyNodeIndex%time(indexHistory-1_c_size_t+indexHistoryStart)
                   else
                      timeHistory=historyPosition %time(indexHistory-1_c_size_t+indexHistoryStart)
                   end if
                   ! Determine node and descendant indices.
                   if (haveNodeIndexHistory) then
                      ! Node index history is available - simply use the indices from the history.
                      nodeIndex      (nodeCount+indexHistory)=historyNodeIndex%data(indexHistory-1_c_size_t+indexHistoryStart,1)
                      if (indexHistory == countHistory-1_c_size_t+indexHistoryStart) then
                         ! At the end of history - set a null descendant index for now.
                         descendantIndex(nodeCount+indexHistory)=-1_c_size_t
                      else
                         descendantIndex(nodeCount+indexHistory)=historyNodeIndex%data(indexHistory-1_c_size_t+indexHistoryStart+1_c_size_t,1)
                      end if
                   else
                      ! Node index history is not available - generate unique indices for history entries using our offset.
                      nodeIndex      (nodeCount+indexHistory)=nodeIndex(nodeCount)+(indexHistory-1_c_size_t+indexHistoryStart           )*nodeIndexOffset
                      descendantIndex(nodeCount+indexHistory)=nodeIndex(nodeCount)+(indexHistory-1_c_size_t+indexHistoryStart+1_c_size_t)*nodeIndexOffset
                   end if
                   ! Set ths host index to that of the current host.
                   hostIndex      (nodeCount+indexHistory)=nodeHost%index()
                   ! Set node mass.
                   if (haveMassHistory) then
                      ! A mass history is available - use it to set the mass.
                      nodeMass(nodeCount+indexHistory)=historyMassBound%data(indexHistory-1_c_size_t+indexHistoryStart,1)
                   else
                      ! No node mass history is available - assuming no change in mass from the original node.
                      nodeMass(nodeCount+indexHistory)=nodeMass             (nodeCount                                  )
                   end if
                   ! Set node redshift using the current history time.
                   nodeRedshift      (nodeCount+indexHistory) = self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeHistory))
                   ! Set node position and velocity.
                   if (havePositionHistory) then
                      ! A position history is available - use it to set the position and velocity.
                      nodePosition   (nodeCount+indexHistory,:)=historyPosition%data    (indexHistory-1_c_size_t+indexHistoryStart,1:3)
                      nodeVelocity   (nodeCount+indexHistory,:)=historyPosition%data    (indexHistory-1_c_size_t+indexHistoryStart,4:6)
                   else
                      ! No position history is available - assume no change in position from the original node.
                      nodePosition   (nodeCount+indexHistory,:)=position       %position(                                             )
                      nodeVelocity   (nodeCount+indexHistory,:)=position       %velocity(                                             )
                   end if
                   ! Set the scale radius if available. This does not change from the original node.
                   if (haveRadiusScale) nodeRadiusScale(nodeCount+indexHistory)=nodeRadiusScale(nodeCount)
                   ! Set the angular momentum if available. This does not change from the original node.
                   if (haveAngularMomentum) then
                      if (angularMomentumIsVector) then
                         nodeAngularMomentum3D(nodeCount+indexHistory,:)=nodeAngularMomentum3D(nodeCount,:)
                      else
                         nodeAngularMomentum  (nodeCount+indexHistory  )=nodeAngularMomentum  (nodeCount  )
                      end if
                   end if
                   ! Set snapshot index if required.
                   if (self%snapshotsRequired) nodeSnapshot(nodeCount+indexHistory)=snapshotInterpolator%locate(timeHistory,closest=.true.)
                   ! Check for a branch jump at this time step.
                   if (haveBranchJump .and. timeHistory == timeEvent(orderEvent(countEvent))) then
                      ! Move the host to the jumped-to node, and increment the event counter to the next branch jump.
                      nodeHost   => nodeEvent_(orderEvent(countEvent))%node
                      countEvent =  countEvent+1
                      ! If no more branch jumps exist, indicate that.
                      if (countEvent > size(timeEvent)) haveBranchJump=.false.
                   else
                      ! No branch jump occurs at this timestep - simply move the host pointer to its parent.
                      nodeHost => nodeHost%parent
                   end if
                end do
                ! Set the descendant index for the original node.
                if (haveNodeIndexHistory) then
                   ! A node index history is available, so use the first entry from that history.
                   descendantIndex(nodeCount)=historyNodeIndex%data(indexHistoryStart,1)
                else
                   ! No node index history is available, so use the offset index.
                   descendantIndex(nodeCount)=nodeIndex(nodeCount)+nodeIndexOffset
                end if
                ! Increment the node count by the number that were added.
                nodeCount=nodeCount+countHistory-indexHistoryStart+1_c_size_t
                ! Act on any subhalo promotion.             
                if (haveSubhaloPromotion) then
                   ! A subhalo promotion occurs at the end of the subhalo history - set the final descendant to the target node of
                   ! the promotion.
                   descendantIndex(nodeCount)=descendantIndexEvent
                else
                   ! No subhalo promotion occurs - check for any merge target.
                   if (associated(node%mergeTarget)) then
                      ! A merge target is defined - set the descendant index to that of the merge target.
                      basicTarget                => node%mergeTarget%basic()
                      descendantIndex(nodeCount) =  node%mergeTarget%index()
                      if (satellite%timeOfMerging() > basicTarget%time()) then
                         if (associated(node%mergeTarget%parent) .and. node%mergeTarget%parent%index() == node%mergeTarget%index()) then
                            ! Merge target is a clone - satellite information is attached to the target node's parent clone.
                            satelliteTarget => node%mergeTarget%parent%satellite()
                         else
                            ! Merge target is not a clone - satellite information is attached directly to the target node.
                            satelliteTarget => node%mergeTarget       %satellite()
                         end if
                         if (satelliteTarget%nodeIndexHistoryIsGettable()) then
                            historyNodeIndexTarget=satelliteTarget%nodeIndexHistory()
                            if (allocated(historyNodeIndexTarget%time)) then
                               do indexHistory=1,size(historyNodeIndexTarget%time)
                                  if (historyNodeIndexTarget%time(indexHistory) == satellite%timeOfMerging()) then
                                     descendantIndex(nodeCount)=historyNodeIndexTarget%data(indexHistory,1)
                                     exit
                                  end if
                               end do
                            end if
                         end if
                      end if
                   else
                      ! No merge target exists. Either merge targets were not set (in which case the assumption is that the subhalo
                      ! merges with its current host), or the final time was reached, in which case the subhalo has no descendant
                      ! (and neither does its host). In either case, we can therefore simply set the descendant index to the current
                      ! host node index (note that this host node will have been progressed one step beyond the end of the subhalo
                      ! history).
                      descendantIndex(nodeCount)=nodeHost%index()
                   end if
                end if
             else if (associated(node%mergeTarget)) then
                ! Handle immediate subhalo-subhalo mergers.
                basicTarget     => node%mergeTarget%basic   ()
                satelliteTarget => node%mergeTarget%satellite()
                if (satellite%timeOfMerging() > basicTarget%time()) then
                   if (satelliteTarget%nodeIndexHistoryIsGettable()) then
                      historyNodeIndexTarget=satelliteTarget%nodeIndexHistory()
                      do indexHistory=1,size(historyNodeIndexTarget%time)
                         if (historyNodeIndexTarget%time(indexHistory) == satellite%timeOfMerging()) then
                            descendantIndex(nodeCount)=historyNodeIndexTarget%data(indexHistory,1)
                            exit
                         end if
                      end do
                   end if
                end if
             end if
          end if
       end do
       ! Store the merger tree data.
       call       mergerTrees%setProperty(propertyTypeTreeWeight      ,treeWeight                )
       call       mergerTrees%setProperty(propertyTypeTreeIndex       ,treeIndex                 )
       call       mergerTrees%setProperty(propertyTypeNodeIndex       ,nodeIndex                 )
       call       mergerTrees%setProperty(propertyTypeHostIndex       ,hostIndex                 )
       call       mergerTrees%setProperty(propertyTypeDescendantIndex ,descendantIndex           )
       call       mergerTrees%setProperty(propertyTypeNodeMass        ,nodeMass                  )
       call       mergerTrees%setProperty(propertyTypeRedshift        ,nodeRedshift              )
       if (havePosition       ) then
          call    mergerTrees%setProperty(propertyTypePositionX       ,nodePosition         (:,1))
          call    mergerTrees%setProperty(propertyTypePositionY       ,nodePosition         (:,2))
          call    mergerTrees%setProperty(propertyTypePositionZ       ,nodePosition         (:,3))
       end if
       if (haveVelocity       ) then
          call    mergerTrees%setProperty(propertyTypeVelocityX       ,nodeVelocity         (:,1))
          call    mergerTrees%setProperty(propertyTypeVelocityY       ,nodeVelocity         (:,2))
          call    mergerTrees%setProperty(propertyTypeVelocityZ       ,nodeVelocity         (:,3))
       end if
       if (haveRadiusScale    ) then
          call    mergerTrees%setProperty(propertyTypeScaleRadius     ,nodeRadiusScale           )
       end if
       if (haveAngularMomentum) then
          if (angularMomentumIsVector) then
             call mergerTrees%setProperty(propertyTypeAngularMomentumX,nodeAngularMomentum3D(:,1))
             call mergerTrees%setProperty(propertyTypeAngularMomentumY,nodeAngularMomentum3D(:,2))
             call mergerTrees%setProperty(propertyTypeAngularMomentumZ,nodeAngularMomentum3D(:,3))
          else
             call mergerTrees%setProperty(propertyTypeAngularMomentum ,nodeAngularMomentum       )
          end if
       end if
       if (self%snapshotsRequired) then
          call    mergerTrees%setProperty(propertyTypeSnapshot        ,nodeSnapshot              )
       end if
       ! Write the tree to file.
       if (nodeCount > 1_c_size_t .or. .not.self%skipSingleNodeTrees) &
            & call mergerTrees%export(char(self%outputFileName),self%exportFormat,hdfChunkSize,hdfCompressionLevel,append=.true.)
       ! Deallocate arrays.
       deallocate       (treeIndex            )
       deallocate       (treeWeight           )
       deallocate       (nodeIndex            )
       deallocate       (hostIndex            )
       deallocate       (descendantIndex      )
       deallocate       (nodeMass             )
       deallocate       (nodeRedshift         )
       deallocate       (nodePosition         )
       deallocate       (nodeVelocity         )
       deallocate       (nodeRadiusScale      )
       deallocate       (nodeAngularMomentum3D)
       deallocate       (nodeAngularMomentum  )
       if (self%snapshotsRequired)              &
            & deallocate(nodeSnapshot         )
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine exportOperatePreEvolution
