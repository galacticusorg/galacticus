!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which writes merger trees to file.

module Merger_Trees_Write
  !% Writes merger trees to file.
  use ISO_Varying_String
  implicit none
  private
  public :: Merger_Tree_Write
  
  ! Flag indicating if module is initialized.
  logical                 :: moduleInitialized=.false.
  
  ! Flag indicating if first tree has been written.
  logical                 :: firstTreeWritten=.false.

  ! Flag indicating if output is required.
  logical                 :: mergerTreesWrite

  ! Options controlling output.
  type   (varying_string) :: mergerTreeExportFileName,mergerTreeExportOutputFormat

  ! Record of whether snapshots are required.  
  logical                 :: needsSnapshots

contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Write</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Write(thisTree)
    !% Output the structure of {\tt thisTree}.
    use Cosmological_Parameters
    use Cosmology_Functions
    use Dates_and_Times
    use CDM_Power_Spectrum
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    use Numerical_Interpolation
    use Merger_Trees
    use Merger_Tree_Data_Structure
    use Tree_Nodes
    use Input_Parameters
    use Memory_Management
    use Kind_Numbers
    use Sort
    use FGSL
    implicit none
    type   (mergerTree       ), intent(in)                  :: thisTree
    integer                   , parameter                   :: hdfChunkSize=1024, hdfCompressionLevel=9
    double precision          , allocatable, dimension(:  ) :: nodeMass,nodeRedshift
    double precision          , allocatable, dimension(:,:) :: nodePosition,nodeVelocity
    integer(kind=kind_int8   ), allocatable, dimension(:  ) :: descendentIndex,treeIndex,nodeIndex,nodeSnapshot
    type   (treeNode         ), pointer                     :: thisNode
    integer                   , parameter                   :: snapshotCountIncrement=100
    double precision          , allocatable, dimension(:  ) :: snapshotTime,snapshotTimeTemp
    integer                                                 :: nodeCount,snapshotCount
    type   (mergerTreeData   )                              :: mergerTrees
    logical                                                 :: snapshotInterpolatorReset
    type   (fgsl_interp_accel)                              :: snapshotInterpolatorAccelerator

    ! Check if module is initialized.
    if (.not.moduleInitialized) then
       !$omp critical(Merger_Tree_Write_Initialize)
       if (.not.moduleInitialized) then
          ! Get parameter specifying if output is required.
          !@ <inputParameter>
          !@   <name>mergerTreesWrite</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not trees should be written to file.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreesWrite',mergerTreesWrite,defaultValue=.false.)
          ! Get the name of the output file name.
          !@ <inputParameter>
          !@   <name>mergerTreeExportFileName</name>
          !@   <defaultValue>galacticusExportedTrees.hdf5</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the file to which merger trees will be exported.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeExportFileName',mergerTreeExportFileName,defaultValue="galacticusExportedTrees.hdf5")
          ! Get the name of the output file format.
          !@ <inputParameter>
          !@   <name>mergerTreeExportOutputFormat</name>
          !@   <defaultValue>galacticus</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The output format to use when exporting merger trees.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeExportOutputFormat',mergerTreeExportOutputFormat,defaultValue="galacticus")
          needsSnapshots=(mergerTreeExportOutputFormat == "irate")
          ! Flag that module is initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical(Merger_Tree_Write_Initialize)
    end if

    ! Write tree to file.
    if (mergerTreesWrite) then

       ! Reset the tree data structure.
       call mergerTrees%reset                   (       )

       ! Specify a single tree in the structure.
       call mergerTrees%treeCountSet            (      1)

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
       call mergerTrees%addMetadata(metaDataCosmology ,'OmegaMatter'       ,Omega_Matter()                 )
       call mergerTrees%addMetadata(metaDataCosmology ,'OmegaBaryon'       ,Omega_B     ()                 )
       call mergerTrees%addMetadata(metaDataCosmology ,'OmegaLambda'       ,Omega_DE    ()                 )
       call mergerTrees%addMetadata(metaDataCosmology ,'HubbleParam'       ,Little_H_0  ()                 )
       call mergerTrees%addMetadata(metaDataCosmology ,'sigma_8'           ,sigma_8     ()                 )
       call mergerTrees%addMetadata(metaDataCosmology ,'powerSpectrumIndex',"not specified"                )

       ! Set provenance metadata.
       call mergerTrees%addMetadata(metaDataProvenance,'fileBuiltBy'       ,'Galacticus'                   )
       call mergerTrees%addMetadata(metaDataProvenance,'fileTimestamp'     ,char(Formatted_Date_and_Time()))

       ! Count nodes in the tree.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          call thisNode%walkTree(thisNode)
       end do
       call mergerTrees%nodeCountSet(nodeCount)

       ! Allocate arrays for serialization.
       call Alloc_Array(treeIndex      ,[nodeCount])
       call Alloc_Array(nodeIndex      ,[nodeCount])
       call Alloc_Array(descendentIndex,[nodeCount])
       call Alloc_Array(nodeMass       ,[nodeCount])
       call Alloc_Array(nodeRedshift   ,[nodeCount])
       if (needsSnapshots                ) call Alloc_Array(nodeSnapshot,[nodeCount  ])
       if (associated(Tree_Node_Position)) call Alloc_Array(nodePosition,[nodeCount,3])
       if (associated(Tree_Node_Velocity)) call Alloc_Array(nodeVelocity,[nodeCount,3])

       ! Find "snapshot" numbers for nodes - relevant only for IRATE output format.
       if (needsSnapshots) then
          call Alloc_Array(snapshotTime,[snapshotCountIncrement])
          thisNode => thisTree%baseNode
          snapshotCount=1
          snapshotTime(snapshotCount)=Tree_Node_Time(thisNode)
          do while (associated(thisNode))
             if (all(snapshotTime(1:snapshotCount) /= Tree_Node_Time(thisNode))) then
                snapshotCount=snapshotCount+1
                if (snapshotCount > size(snapshotTime)) then
                   call Move_Alloc(snapshotTime,snapshotTimeTemp)
                   call Alloc_Array(snapshotTime,[size(snapshotTimeTemp)+snapshotCountIncrement])
                   snapshotTime(1:size(snapshotTimeTemp))=snapshotTimeTemp
                   call Dealloc_Array(snapshotTimeTemp)
                end if
                snapshotTime(snapshotCount)=Tree_Node_Time(thisNode)
             end if
             call thisNode%walkTree(thisNode)
          end do
          call Sort_Do(snapshotTime(1:snapshotCount))
       end if

       ! Serialize node data to arrays and write to merger tree data structure.
       treeIndex=thisTree%index
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeIndex      (nodeCount)=thisNode           %index()
          descendentIndex(nodeCount)=thisNode%parentNode%index()
          nodeMass       (nodeCount)=Tree_Node_Mass                                                (thisNode)
          nodeRedshift   (nodeCount)=Redshift_From_Expansion_Factor(Expansion_Factor(Tree_Node_Time(thisNode)))
          if (associated(Tree_Node_Position)) call Tree_Node_Position(thisNode,nodePosition(nodeCount,:))
          if (associated(Tree_Node_Velocity)) call Tree_Node_Velocity(thisNode,nodeVelocity(nodeCount,:))
          if (needsSnapshots) then
             nodeSnapshot(nodeCount)=Interpolate_Locate(snapshotCount,snapshotTime,snapshotInterpolatorAccelerator&
                  &,Tree_Node_Time(thisNode),reset=snapshotInterpolatorReset,closest=.true.)
          end if
          call thisNode%walkTree(thisNode)
       end do
       call Interpolate_Done(interpolationAccelerator=snapshotInterpolatorAccelerator,reset=snapshotInterpolatorReset)
       call mergerTrees%setProperty(propertyTypeTreeIndex      ,treeIndex      )
       call mergerTrees%setProperty(propertyTypeNodeIndex      ,nodeIndex      )
       call mergerTrees%setProperty(propertyTypeHostIndex      ,nodeIndex      )
       call mergerTrees%setProperty(propertyTypeDescendentIndex,descendentIndex)
       call mergerTrees%setProperty(propertyTypeNodeMass       ,nodeMass       )
       call mergerTrees%setProperty(propertyTypeRedshift       ,nodeRedshift   )
       if (associated(Tree_Node_Position)) then
          call mergerTrees%setProperty(propertyTypePositionX,nodePosition(:,1))
          call mergerTrees%setProperty(propertyTypePositionY,nodePosition(:,2))
          call mergerTrees%setProperty(propertyTypePositionZ,nodePosition(:,3))
       end if
       if (associated(Tree_Node_Velocity)) then
          call mergerTrees%setProperty(propertyTypeVelocityX,nodeVelocity(:,1))
          call mergerTrees%setProperty(propertyTypeVelocityY,nodeVelocity(:,2))
          call mergerTrees%setProperty(propertyTypeVelocityZ,nodeVelocity(:,3))
       end if
       if (needsSnapshots) call mergerTrees%setProperty(propertyTypeSnapshot,nodeSnapshot)

       ! Write the tree to file.
       !$omp critical (Merger_Tree_Write)
       call mergerTrees%export(char(mergerTreeExportFileName),char(mergerTreeExportOutputFormat),hdfChunkSize,hdfCompressionLevel,append=firstTreeWritten)
       firstTreeWritten=.true.
       !$omp end critical (Merger_Tree_Write)

       ! Deallocate arrays.
       call Dealloc_Array(treeIndex      )
       call Dealloc_Array(nodeIndex      )
       call Dealloc_Array(descendentIndex)
       call Dealloc_Array(nodeMass       )
       call Dealloc_Array(nodeRedshift   )
       if (associated(Tree_Node_Position)) call Dealloc_Array(nodePosition)
       if (associated(Tree_Node_Velocity)) call Dealloc_Array(nodeVelocity)

    end if

    return
  end subroutine Merger_Tree_Write
  
end module Merger_Trees_Write
