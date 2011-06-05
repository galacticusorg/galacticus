!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
    use Tree_Node_Methods
    use Input_Parameters
    use Memory_Management
    use Galacticus_HDF5_Groups
    use HDF5
    use ISO_Varying_String
    use String_Handling
    use Dark_Matter_Halo_Scales
    !# <include directive="mergerTreeStructureOutputTask" type="moduleUse">
    include 'merger_trees.output_structure.tasks.modules.inc'
    !# </include>
    implicit none
    type(mergerTree), intent(in)                :: thisTree
    type(treeNode),   pointer                   :: thisNode
    integer,          allocatable, dimension(:) :: nodeIndex
    double precision, allocatable, dimension(:) :: nodeProperty
    integer                                     :: nodeCount,structureDataID,treeGroupID
    type(varying_string)                        :: groupName,groupComment

    ! Check if module is initialized.
    if (.not.structureOutputModuleInitialized) then
       ! Get parameter specifying if output is required.
       !@ <inputParameter>
       !@   <name>mergerTreeStructureOutput</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not to output the structure of merger trees prior to evolution.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeStructureOutput',mergerTreeStructureOutput,defaultValue=.false.)
       !@ <inputParameter>
       !@   <name>mergerTreeStructureOutputVirialQuantities</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not to output virial quantities (radius and velocity) when outputting the structure of merger trees prior to evolution.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeStructureOutputVirialQuantities',mergerTreeStructureOutputVirialQuantities,defaultValue=.false.)
       ! Create an output group if necessary.
       if (mergerTreeStructureOutput) then
          groupName   ='mergerTreeStructures'
          groupComment='Pre-evolution structures of merger trees.'
          structureGroupID=Galacticus_Output_Make_Group(groupName,groupComment)
       end if
       ! Flag that module is initialized.
       structureOutputModuleInitialized=.true.
    end if

    ! Output the tree structure history.
    if (mergerTreeStructureOutput) then
       ! Count up number of nodes in the tree.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          call thisNode%walkTree()
       end do
       ! Allocate storage space.
       call Alloc_Array(nodeIndex   ,nodeCount,'nodeIndex'   )
       call Alloc_Array(nodeProperty,nodeCount,'nodeProperty')
       ! Create a group for this tree structure.
       groupName   ='mergerTree'
       groupName   =groupName//thisTree%index
       groupComment='Pre-evolution structure of merger tree.'
       treeGroupID =Galacticus_Output_Make_Group(groupName,groupComment,structureGroupID)

       ! Extract node indices and output to file.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeIndex(nodeCount)=thisNode%index()        
          call thisNode%walkTree()
       end do
       structureDataID=0
       call Galacticus_Output_Dataset(treeGroupID,structureDataID,'nodeIndex','Index of the node.',nodeIndex)
    
       ! Extract child node indices and output to file.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeIndex(nodeCount)=thisNode%childNode%index()        
          call thisNode%walkTree()
       end do
       structureDataID=0
       call Galacticus_Output_Dataset(treeGroupID,structureDataID,'childNodeIndex','Index of the child node.',nodeIndex)
    
       ! Extract parent node indices and output to file.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeIndex(nodeCount)=thisNode%parentNode%index()        
          call thisNode%walkTree()
       end do
       structureDataID=0
       call Galacticus_Output_Dataset(treeGroupID,structureDataID,'parentNodeIndex','Index of the parent node.',nodeIndex)
    
       ! Extract sibling node indices and output to file.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeIndex(nodeCount)=thisNode%siblingNode%index()        
          call thisNode%walkTree()
       end do
       structureDataID=0
       call Galacticus_Output_Dataset(treeGroupID,structureDataID,'siblingNodeIndex','Index of the sibling node.',nodeIndex)
    
       ! Extract node masses and output to file.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeProperty(nodeCount)=Tree_Node_Mass(thisNode)
          call thisNode%walkTree()
       end do
       structureDataID=0
       call Galacticus_Output_Dataset(treeGroupID,structureDataID,'nodeMass','Mass of node.',nodeProperty)
    
       ! Extract node times and output to file.
       nodeCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeProperty(nodeCount)=Tree_Node_Time(thisNode)
          call thisNode%walkTree()
       end do
       structureDataID=0
       call Galacticus_Output_Dataset(treeGroupID,structureDataID,'nodeTime','Time at node.',nodeProperty)
    
       ! Check whether output of virial quantities is required.
       if (mergerTreeStructureOutputVirialQuantities) then
          
          ! Extract node virial radii and output to file.
          nodeCount=0
          thisNode => thisTree%baseNode
          do while (associated(thisNode))
             nodeCount=nodeCount+1
             nodeProperty(nodeCount)=Dark_Matter_Halo_Virial_Radius(thisNode)
             call thisNode%walkTree()
          end do
          structureDataID=0
          call Galacticus_Output_Dataset(treeGroupID,structureDataID,'nodeVirialRadius','Virial radius of the node [Mpc].',nodeProperty)
          
          ! Extract node virial velocity and output to file.
          nodeCount=0
          thisNode => thisTree%baseNode
          do while (associated(thisNode))
             nodeCount=nodeCount+1
             nodeProperty(nodeCount)=Dark_Matter_Halo_Virial_Velocity(thisNode)
             call thisNode%walkTree()
          end do
          structureDataID=0
          call Galacticus_Output_Dataset(treeGroupID,structureDataID,'nodeVirialVelocity','Virial velocity of the node [km/s].',nodeProperty)
    
       end if
       
       ! Call any subroutines that want to attach data to the merger tree output.
       !# <include directive="mergerTreeStructureOutputTask" type="code" action="subroutine">
       !#  <subroutineArgs>thisTree%baseNode,nodeProperty,treeGroupID</subroutineArgs>
       include 'merger_trees.output_structure.tasks.inc'
       !# </include>

       ! Deallocate storage space.
       call Dealloc_Array(nodeIndex   )
       call Dealloc_Array(nodeProperty)
    end if

    return
  end subroutine Merger_Tree_Structure_Output
  
end module Merger_Tree_Output_Structure
