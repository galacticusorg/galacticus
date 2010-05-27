!% Contains a module which outputs the structure of entire merger trees.

module Merger_Tree_Output_Structure
  !% Outputs the structure of entire merger trees.
  private
  public :: Merger_Tree_Structure_Output
  
  ! Flag indicating if module is initialized.
  logical :: structureOutputModuleInitialized=.false.
  
  ! Flag indicating if output is required.
  logical :: mergerTreeStructureOutput

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
    
       ! Deallocate storage space.
       call Dealloc_Array(nodeIndex   )
       call Dealloc_Array(nodeProperty)
    end if

    return
  end subroutine Merger_Tree_Structure_Output
  
end module Merger_Tree_Output_Structure
