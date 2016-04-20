!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% Contains a module which implements a merger tree operator which exports merger trees to
  !% file.
  
  !# <mergerTreeOperator name="mergerTreeOperatorExport">
  !#  <description>
  !#   A merger tree operator which exports merger trees to file.
  !# </description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorExport
     !% A merger tree operator class which exports merger trees to file.
     private
     type   (varying_string) :: outputFileName
     integer                 :: exportFormat
     logical                 :: snapshotsRequired
   contains
     final     ::             exportDestructor
     procedure :: operate  => exportOperate
  end type mergerTreeOperatorExport
  
  interface mergerTreeOperatorExport
     !% Constructors for the export merger tree operator class.
     module procedure exportConstructorParameters
     module procedure exportConstructorInternal
  end interface mergerTreeOperatorExport

contains

  function exportConstructorParameters(parameters)
    !% Constructor for the export merger tree operator class which takes a
    !% parameter set as input.
    use Merger_Tree_Data_Structure
    use Input_Parameters2
    implicit none
    type   (mergerTreeOperatorExport)                :: exportConstructorParameters
    type   (inputParameters         ), intent(inout) :: parameters
    type   (varying_string          )                   outputFileName             , exportFormatText
    integer                                          :: exportFormat
    !# <inputParameterList label="allowedParameterNames" />

    call parameters%checkParameters(allowedParameterNames)
    !# <inputParameter>
    !#   <name>outputFileName</name>
    !#   <source>parameters</source>
    !#   <defaultValue>var_str('galacticusExportedTrees.hdf5')</defaultValue>
    !#   <description>The name of the file to which merger trees should be exported.</description>
    !#   <type>string</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exportFormat</name>
    !#   <source>parameters</source>
    !#   <variable>exportFormatText</variable>
    !#   <defaultValue>var_str('galacticus')</defaultValue>
    !#   <description>The output format to use when exporting merger trees.</description>
    !#   <type>string</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    ! Validate the output format.
    exportFormat=enumerationMergerTreeFormatEncode(char(exportFormatText),includesPrefix=.false.)
    ! Construct the object.
    exportConstructorParameters=exportConstructorInternal(char(outputFileName),exportFormat)
    return
  end function exportConstructorParameters

  function exportConstructorInternal(outputFileName,exportFormat)
    !% Internal constructor for the export merger tree operator class.
    use Merger_Tree_Data_Structure
    use Galacticus_Error
    implicit none
    type     (mergerTreeOperatorExport)                :: exportConstructorInternal
    character(len=*                   ), intent(in   ) :: outputFileName
    integer                            , intent(in   ) :: exportFormat

    ! Validate the export format.
    if (.not.enumerationMergerTreeFormatIsValid(exportFormat)) call Galacticus_Error_Report('exportConstructorInternal','exportFormat is invalid')
    ! Construct the object.
    exportConstructorInternal%outputFileName   =outputFileName
    exportConstructorInternal%exportFormat     =exportFormat
    exportConstructorInternal%snapshotsRequired=(exportFormat == mergerTreeFormatIrate)
    return
  end function exportConstructorInternal

  elemental subroutine exportDestructor(self)
    !% Destructor for the export merger tree operator function class.
    implicit none
    type(mergerTreeOperatorExport), intent(inout) :: self
    !GCC$ attributes unused :: self
    
    ! Nothing to do.
    return
  end subroutine exportDestructor

  subroutine exportOperate(self,tree)
    !% Output the structure of {\normalfont \ttfamily thisTree}.
    use HDF5
    use Cosmology_Parameters
    use Cosmology_Functions
    use Dates_and_Times
    use Numerical_Constants_Astronomical
    use Numerical_Interpolation
    use Merger_Tree_Data_Structure
    use Galacticus_Nodes
    use Input_Parameters
    use Memory_Management
    use Sort
    use File_Utilities
    use System_Command
    use Cosmological_Mass_Variance
    implicit none
    class           (mergerTreeOperatorExport), intent(inout)                 :: self
    type            (mergerTree              ), intent(inout), target         :: tree
    integer         (kind=size_t             ), parameter                     :: hdfChunkSize                   =1024
    integer                                   , parameter                     :: hdfCompressionLevel            =   9
    double precision                          , allocatable  , dimension(:  ) :: nodeMass                            , nodeRedshift         , &
         &                                                                       snapshotTime                        , snapshotTimeTemp     , &
         &                                                                       treeWeight
    double precision                          , allocatable  , dimension(:,:) :: nodePosition                        , nodeVelocity
    integer         (kind=kind_int8          ), allocatable  , dimension(:  ) :: descendentIndex                     , nodeIndex            , &
         &                                                                       nodeSnapshot                        , treeIndex
    type            (treeNode                ), pointer                       :: node
    class           (nodeComponentBasic      ), pointer                       :: basic
    class           (nodeComponentPosition   ), pointer                       :: position
    type            (mergerTree              ), pointer                       :: treeCurrent
    class           (cosmologyParametersClass), pointer                       :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer                       :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass), pointer                                               :: cosmologicalMassVariance_
    integer                                   , parameter                     :: snapshotCountIncrement         = 100
    integer                                                                   :: nodeCount                           , snapshotCount
    type            (mergerTreeData          )                                :: mergerTrees
    logical                                                                   :: snapshotInterpolatorReset
    type            (fgsl_interp_accel       )                                :: snapshotInterpolatorAccelerator

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
       ! Get the default cosmology.
       cosmologyParameters_      => cosmologyParameters     ()
       cosmologicalMassVariance_ => cosmologicalMassVariance()
       ! Set cosmology metadata.
       call mergerTrees%addMetadata(metaDataTypeCosmology ,'OmegaMatter'       ,cosmologyParameters_     %OmegaMatter    (                  ))
       call mergerTrees%addMetadata(metaDataTypeCosmology ,'OmegaBaryon'       ,cosmologyParameters_     %OmegaBaryon    (                  ))
       call mergerTrees%addMetadata(metaDataTypeCosmology ,'OmegaLambda'       ,cosmologyParameters_     %OmegaDarkEnergy(                  ))
       call mergerTrees%addMetadata(metaDataTypeCosmology ,'HubbleParam'       ,cosmologyParameters_     %HubbleConstant (hubbleUnitsLittleH))
       call mergerTrees%addMetadata(metaDataTypeCosmology ,'sigma_8'           ,cosmologicalMassVariance_%sigma8         (                  ))
       call mergerTrees%addMetadata(metaDataTypeCosmology ,'powerSpectrumIndex',"not specified"                                              )
       ! Set provenance metadata.
       call mergerTrees%addMetadata(metaDataTypeProvenance,'fileBuiltBy'       ,'Galacticus'                   )
       call mergerTrees%addMetadata(metaDataTypeProvenance,'fileTimestamp'     ,char(Formatted_Date_and_Time()))       
       ! Count nodes in the tree.
       nodeCount=0
       node => treeCurrent%baseNode
       do while (associated(node))
          nodeCount=nodeCount+1
          node => node%walkTree()
       end do
       call mergerTrees%nodeCountSet(nodeCount)
       ! Allocate arrays for serialization.
       call Alloc_Array(treeIndex      ,[nodeCount])
       call Alloc_Array(treeWeight     ,[nodeCount])
       call Alloc_Array(nodeIndex      ,[nodeCount])
       call Alloc_Array(descendentIndex,[nodeCount])
       call Alloc_Array(nodeMass       ,[nodeCount])
       call Alloc_Array(nodeRedshift   ,[nodeCount])
       if (self%snapshotsRequired                               ) call Alloc_Array(nodeSnapshot,[nodeCount  ])
       if (defaultPositionComponent%positionIsGettable()) call Alloc_Array(nodePosition,[nodeCount,3])
       if (defaultPositionComponent%velocityIsGettable()) call Alloc_Array(nodeVelocity,[nodeCount,3])
       ! Find "snapshot" numbers for nodes - relevant only for IRATE output format.
       if (self%snapshotsRequired) then
          call Alloc_Array(snapshotTime,[snapshotCountIncrement])
          node           => treeCurrent%baseNode
          basic => node%basic()
          snapshotCount=1
          snapshotTime(snapshotCount)=basic%time()
          do while (associated(node))
             basic => node%basic()
             if (all(snapshotTime(1:snapshotCount) /= basic%time())) then
                snapshotCount=snapshotCount+1
                if (snapshotCount > size(snapshotTime)) then
                   call Move_Alloc(snapshotTime,snapshotTimeTemp)
                   call Alloc_Array(snapshotTime,[size(snapshotTimeTemp)+snapshotCountIncrement])
                   snapshotTime(1:size(snapshotTimeTemp))=snapshotTimeTemp
                   call Dealloc_Array(snapshotTimeTemp)
                end if
                snapshotTime(snapshotCount)=basic%time()
             end if
             node => node%walkTree()
          end do
          call Sort_Do(snapshotTime(1:snapshotCount))
       end if
       ! Get the default cosmology functions object.
       cosmologyFunctions_ => cosmologyFunctions()
       ! Serialize node data to arrays and write to merger tree data structure.
       treeIndex =treeCurrent%index
       treeWeight=treeCurrent%volumeWeight
       nodeCount =0
       node  => treeCurrent%baseNode
       snapshotInterpolatorReset=.true.
       do while (associated(node))
          nodeCount=nodeCount+1
          nodeIndex      (nodeCount)=  node       %index   ()
          descendentIndex(nodeCount)=  node%parent%index   ()
          basic                     => node       %basic   ()
          position                  => node       %position()
          nodeMass       (nodeCount)=                                                                                    basic%mass()
          nodeRedshift   (nodeCount)=cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%expansionFactor(basic%time()))
          if (defaultPositionComponent%positionIsGettable()) nodePosition(nodeCount,:)=position%position()
          if (defaultPositionComponent%velocityIsGettable()) nodeVelocity(nodeCount,:)=position%velocity()
          if (self%snapshotsRequired)                                                                                      &
               & nodeSnapshot(nodeCount)=Interpolate_Locate(                                                               & 
               &                                                    snapshotTime                        (1:snapshotCount), &
               &                                                    snapshotInterpolatorAccelerator                      , &
               &                                                    basic                          %time(               ), &
               &                                            reset  =snapshotInterpolatorReset                            , &
               &                                            closest=.true.                                                 &
               &                                           )
          node => node%walkTree()
       end do
       call Interpolate_Done(interpolationAccelerator=snapshotInterpolatorAccelerator,reset=snapshotInterpolatorReset)
       call mergerTrees%setProperty(propertyTypeTreeWeight     ,treeWeight     )
       call mergerTrees%setProperty(propertyTypeTreeIndex      ,treeIndex      )
       call mergerTrees%setProperty(propertyTypeNodeIndex      ,nodeIndex      )
       call mergerTrees%setProperty(propertyTypeHostIndex      ,nodeIndex      )
       call mergerTrees%setProperty(propertyTypeDescendentIndex,descendentIndex)
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
       !$omp critical (Merger_Tree_Write)
       call mergerTrees%export(char(self%outputFileName),self%exportFormat,hdfChunkSize,hdfCompressionLevel,append=.true.)
       !$omp end critical (Merger_Tree_Write)       
       ! Deallocate arrays.
       call Dealloc_Array(treeIndex      )
       call Dealloc_Array(treeWeight     )
       call Dealloc_Array(nodeIndex      )
       call Dealloc_Array(descendentIndex)
       call Dealloc_Array(nodeMass       )
       call Dealloc_Array(nodeRedshift   )
       if (defaultPositionComponent%positionIsGettable()) call Dealloc_Array(nodePosition)
       if (defaultPositionComponent%velocityIsGettable()) call Dealloc_Array(nodeVelocity)
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do    
    return
  end subroutine exportOperate
