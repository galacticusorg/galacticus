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


!% Contains a module which outputs the structure of entire merger trees.

module Merger_Tree_Output_Structure
  !% Outputs the structure of entire merger trees.
  implicit none
  private
  public :: Merger_Tree_Structure_Output
  
  ! Flag indicating if module is initialized.
  logical :: structureOutputModuleInitialized=.false.
  
  ! Flag indicating if output is required.
  logical :: mergerTreeStructureOutput

  ! Flag indicating if virial quantities should be included in output.
  logical :: mergerTreeStructureOutputVirialQuantities

  ! HDF5 group index.
  integer :: structureGroupID
  
contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Structure_Output</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Structure_Output(thisTree)
    !% Output the structure of {\tt thisTree}.
    use Merger_Trees
    use Tree_Nodes
    use Input_Parameters
    use Memory_Management
    use Galacticus_HDF5
    use IO_HDF5
    use ISO_Varying_String
    use String_Handling
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    use Kind_Numbers
    !# <include directive="mergerTreeStructureOutputTask" type="moduleUse">
    include 'merger_trees.output_structure.tasks.modules.inc'
    !# </include>
    implicit none
    type(mergerTree),        intent(in)                :: thisTree
    type(treeNode),          pointer                   :: thisNode
    integer(kind=kind_int8), allocatable, dimension(:) :: nodeIndex
    double precision,        allocatable, dimension(:) :: nodeProperty
    type(hdf5Object),        save                      :: structureGroup
    integer                                            :: nodeCount
    type(varying_string)                               :: groupName,groupComment
    type(hdf5Object)                                   :: treeGroup,nodeDataset

    ! Check if module is initialized.
    !$omp critical(structureOutputModuleInitialize)
    if (.not.structureOutputModuleInitialized) then
       ! Get parameter specifying if output is required.
       !@ <inputParameter>
       !@   <name>mergerTreeStructureOutput</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not to output the structure of merger trees prior to evolution.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeStructureOutput',mergerTreeStructureOutput,defaultValue=.false.)
       !@ <inputParameter>
       !@   <name>mergerTreeStructureOutputVirialQuantities</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not to output virial quantities (radius and velocity) when outputting the structure of merger trees prior to evolution.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeStructureOutputVirialQuantities',mergerTreeStructureOutputVirialQuantities,defaultValue=.false.)
       ! Create an output group if necessary.
       !$omp critical(HDF5_Access)
       if (mergerTreeStructureOutput) structureGroup=IO_HDF5_Open_Group(galacticusOutputFile,'mergerTreeStructures','Pre-evolution structures of merger trees.')
       !$omp end critical(HDF5_Access)
       ! Flag that module is initialized.
       structureOutputModuleInitialized=.true.
    end if
    !$omp end critical(structureOutputModuleInitialize)

    ! Output the tree structure history.
    if (mergerTreeStructureOutput) then
       ! Count up number of nodes in the tree.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          call thisNode%walkTree(thisNode)
       end do
       ! Allocate storage space.
       call Alloc_Array(nodeIndex   ,[nodeCount])
       call Alloc_Array(nodeProperty,[nodeCount])
       ! Create a group for this tree structure.
       groupName   ='mergerTree'
       groupName   =groupName//thisTree%index
       groupComment='Pre-evolution structure of merger tree.'
       !$omp critical(HDF5_Access)
       treeGroup   =IO_HDF5_Open_Group(structureGroup,char(groupName),'Pre-evolution structure of merger tree.')
       !$omp end critical(HDF5_Access)

       ! Extract node indices and output to file.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeIndex(nodeCount)=thisNode%index()        
          call thisNode%walkTree(thisNode)
       end do
       !$omp critical(HDF5_Access)
       call treeGroup%writeDataset(nodeIndex,'nodeIndex','Index of the node.')
       !$omp end critical(HDF5_Access)
       
       ! Extract child node indices and output to file.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeIndex(nodeCount)=thisNode%childNode%index()        
          call thisNode%walkTree(thisNode)
       end do
       !$omp critical(HDF5_Access)
       call treeGroup%writeDataset(nodeIndex,'childNodeIndex','Index of the child node.')
       !$omp end critical(HDF5_Access)
       
       ! Extract parent node indices and output to file.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeIndex(nodeCount)=thisNode%parentNode%index()        
          call thisNode%walkTree(thisNode)
       end do
       !$omp critical(HDF5_Access)
       call treeGroup%writeDataset(nodeIndex,'parentNodeIndex','Index of the parent node.')
       !$omp end critical(HDF5_Access)
       
       ! Extract sibling node indices and output to file.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeIndex(nodeCount)=thisNode%siblingNode%index()        
          call thisNode%walkTree(thisNode)
       end do
       !$omp critical(HDF5_Access)
       call treeGroup%writeDataset(nodeIndex,'siblingNodeIndex','Index of the sibling node.')
       !$omp end critical(HDF5_Access)
       
       ! Extract node masses and output to file.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeProperty(nodeCount)=Tree_Node_Mass(thisNode)
          call thisNode%walkTree(thisNode)
       end do
       !$omp critical(HDF5_Access)
       call treeGroup%writeDataset(nodeProperty,'nodeMass','Mass of node.',datasetReturned=nodeDataset)
       call nodeDataset%writeAttribute(massSolar,"unitsInSI")
       call nodeDataset%close()
       !$omp end critical(HDF5_Access)
       
       ! Extract node times and output to file.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeProperty(nodeCount)=Tree_Node_Time(thisNode)
          call thisNode%walkTree(thisNode)
       end do
       !$omp critical(HDF5_Access)
       call treeGroup%writeDataset(nodeProperty,'nodeTime','Time at node.',datasetReturned=nodeDataset)
       call nodeDataset%writeAttribute(gigaYear,"unitsInSI")
       call nodeDataset%close()
       !$omp end critical(HDF5_Access)
       
       ! Check whether output of virial quantities is required.
       if (mergerTreeStructureOutputVirialQuantities) then
          
          ! Extract node virial radii and output to file.
          nodeCount=0
          thisNode => thisTree%baseNode
          do while (associated(thisNode))
             nodeCount=nodeCount+1
             nodeProperty(nodeCount)=Dark_Matter_Halo_Virial_Radius(thisNode)
             call thisNode%walkTree(thisNode)
          end do
          !$omp critical(HDF5_Access)
          call treeGroup%writeDataset(nodeProperty,'nodeVirialRadius','Virial radius of the node [Mpc].',datasetReturned=nodeDataset)
          call nodeDataset%writeAttribute(megaParsec,"unitsInSI")
          call nodeDataset%close()
          !$omp end critical(HDF5_Access)
          
          ! Extract node virial velocity and output to file.
          nodeCount=0
          thisNode => thisTree%baseNode
          do while (associated(thisNode))
             nodeCount=nodeCount+1
             nodeProperty(nodeCount)=Dark_Matter_Halo_Virial_Velocity(thisNode)
             call thisNode%walkTree(thisNode)
          end do
          !$omp critical(HDF5_Access)
          call treeGroup%writeDataset(nodeProperty,'nodeVirialVelocity','Virial velocity of the node [km/s].',datasetReturned=nodeDataset)
          call nodeDataset%writeAttribute(kilo,"unitsInSI")
          call nodeDataset%close()
          !$omp end critical(HDF5_Access)
          
       end if
       
       !$omp critical(HDF5_Access)
       ! Call any subroutines that want to attach data to the merger tree output.
       !# <include directive="mergerTreeStructureOutputTask" type="code" action="subroutine">
       !#  <subroutineArgs>thisTree%baseNode,nodeProperty,treeGroup</subroutineArgs>
       include 'merger_trees.output_structure.tasks.inc'
       !# </include>
       
       ! Close output group.
       call treeGroup%close()
       !$omp end critical(HDF5_Access)
       
       ! Deallocate storage space.
       call Dealloc_Array(nodeIndex   )
       call Dealloc_Array(nodeProperty)
    end if
    
    return
  end subroutine Merger_Tree_Structure_Output
  
end module Merger_Tree_Output_Structure
