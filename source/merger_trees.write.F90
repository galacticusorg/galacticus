!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which writes merger trees to file.

module Merger_Trees_Write
  !% Writes merger trees to file.
  implicit none
  private
  public :: Merger_Tree_Write
  
  ! Flag indicating if module is initialized.
  logical :: moduleInitialized=.false.
  
  ! Flag indicating if first tree has been written.
  logical :: firstTreeWritten=.false.

  ! Flag indicating if output is required.
  logical :: mergerTreesWrite

contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Write</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Write(thisTree)
    !% Output the structure of {\tt thisTree}.
    use Cosmological_Parameters
    use Cosmology_Functions
    use ISO_Varying_String
    use Dates_and_Times
    use CDM_Power_Spectrum
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    use Merger_Trees
    use Merger_Tree_Data_Structure
    use Tree_Nodes
    use Input_Parameters
    use Memory_Management
    use Kind_Numbers
    implicit none
    type   (mergerTree    ), intent(in)                :: thisTree
    integer                , parameter                 :: hdfChunkSize=1024, hdfCompressionLevel=9
    double precision       , allocatable, dimension(:) :: nodeMass,nodeRedshift
    integer(kind=kind_int8), allocatable, dimension(:) :: descendentIndex,treeIndex,nodeIndex
    type   (treeNode      ), pointer                   :: thisNode
    integer                                            :: nodeCount
    type   (mergerTreeData)                            :: mergerTrees
    type   (varyingString )                            :: mergerTreeExportFileName

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
       call mergerTrees%setUnits(unitsMass    ,unitsInSI=massSolar ,hubbleExponent=0,scaleFactorExponent=0)
       call mergerTrees%setUnits(unitsLength  ,unitsInSI=megaParsec,hubbleExponent=0,scaleFactorExponent=0)
       call mergerTrees%setUnits(unitsVelocity,unitsInSI=kilo      ,hubbleExponent=0,scaleFactorExponent=0)

       ! Set cosmology metadata.
       call mergerTrees%addMetadata(metaDataCosmology ,'OmegaMatter'       ,Omega_Matter()                 )
       call mergerTrees%addMetadata(metaDataCosmology ,'OmegaBaryon'       ,Omega_B     ()                 )
       call mergerTrees%addMetadata(metaDataCosmology ,'OmegaLambda'       ,Omega_DE    ()                 )
       call mergerTrees%addMetadata(metaDataCosmology ,'HubbleParam'       ,H_0         ()                 )
       call mergerTrees%addMetadata(metaDataCosmology ,'sigma_8'           ,sigma_8     ()                 )

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
          call thisNode%walkTree(thisNode)
       end do
       call mergerTrees%setProperty(propertyTypeTreeIndex      ,treeIndex      )
       call mergerTrees%setProperty(propertyTypeNodeIndex      ,nodeIndex      )
       call mergerTrees%setProperty(propertyTypeDescendentIndex,descendentIndex)
       call mergerTrees%setProperty(propertyTypeNodeMass       ,nodeMass       )
       call mergerTrees%setProperty(propertyTypeRedshift       ,nodeRedshift   )

       ! Write the tree to file.
       !$omp critical (Merger_Tree_Write)
       call mergerTrees%export(char(mergerTreeExportFileName),hdfChunkSize,hdfCompressionLevel,append=firstTreeWritten)
       firstTreeWritten=.true.
       !$omp end critical (Merger_Tree_Write)

       ! Deallocate arrays.
       call Dealloc_Array(treeIndex      )
       call Dealloc_Array(nodeIndex      )
       call Dealloc_Array(descendentIndex)
       call Dealloc_Array(nodeMass       )
       call Dealloc_Array(nodeRedshift   )

    end if

    return
  end subroutine Merger_Tree_Write
  
end module Merger_Trees_Write
