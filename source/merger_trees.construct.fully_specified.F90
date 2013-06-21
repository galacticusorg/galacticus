!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements building of merger trees using a fully-specified description read from file.

module Merger_Trees_Construct_Fully_Specified
  !% Implements building of merger trees using a fully-specified description read from file.
  use Merger_Trees
  use ISO_Varying_String
  implicit none
  private
  public :: Merger_Tree_Construct_Fully_Specified_Initialize

  ! Record of whether the tree has yet been processed.
  logical                 :: treeProcessed                            =.false. 
  
  ! Name of the file to read the merger tree definition from.
  type   (varying_string) :: mergerTreeConstructFullySpecifiedFileName         
  
contains

  !# <mergerTreeConstructMethod>
  !#  <unitName>Merger_Tree_Construct_Fully_Specified_Initialize</unitName>
  !# </mergerTreeConstructMethod>
  subroutine Merger_Tree_Construct_Fully_Specified_Initialize(mergerTreeConstructMethod,Merger_Tree_Construct)
    !% Initializes the merger tree construction ``fully-specified'' module.
    use Input_Parameters
    implicit none
    type     (varying_string), intent(in   )          :: mergerTreeConstructMethod 
    procedure(              ), intent(inout), pointer :: Merger_Tree_Construct     
    
    ! Check if our method is to be used.
    if (mergerTreeConstructMethod == 'fullySpecified') then
       ! Assign pointer to our merger tree construction subroutine.
       Merger_Tree_Construct => Merger_Tree_Construct_Fully_Specified
       ! Read parameter giving the name of the file to read.
       !@ <inputParameter>
       !@   <name>mergerTreeConstructFullySpecifiedFileName</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the file giving the fully-specified description of the merger tree to process.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeConstructFullySpecifiedFileName',mergerTreeConstructFullySpecifiedFileName)
    end if
    return
  end subroutine Merger_Tree_Construct_Fully_Specified_Initialize
  
  subroutine Merger_Tree_Construct_Fully_Specified(thisTree,skipTree)
    !% Construct a fully-specified merger tree.
    use Galacticus_Nodes
    use FoX_Dom
    use Kind_Numbers
    use Galacticus_Error
    use Memory_Management
    use Galacticus_Display
    implicit none
    type   (mergerTree    )             , intent(inout) :: thisTree                               
    logical                             , intent(in   ) :: skipTree                               
    type   (treeNodeList  ), allocatable, dimension(:)  :: nodeArray                              
    type   (node          ), pointer                    :: doc        , nodeDefinition            
    type   (nodeList      ), pointer                    :: nodes                                  
    integer                                             :: i          , ioErr         , nodeCount 
    integer(kind=kind_int8)                             :: indexValue                             
    logical                                             :: processTree                            
    
    !$omp critical(Merger_Tree_Construct_Fully_Specified_Process)
    ! If the tree is already processed, return.
    processTree=.true.
    if (treeProcessed) processTree=.false.
    ! Mark that the tree has now been processed.
    treeProcessed=.true.
    !$omp end critical(Merger_Tree_Construct_Fully_Specified_Process)
    ! Return if the tree is to be skipped.
    if (skipTree.or..not.processTree) return

    ! Parse the definition file.
    !$omp critical (FoX_DOM_Access)
    doc => parseFile(char(mergerTreeConstructFullySpecifiedFileName),iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('Merger_Tree_Construct_Fully_Specified','unable to read or parse merger tree root mass file')
    ! Get the list of nodes.
    nodes => getElementsByTagname(doc,"node")
    nodeCount=getLength(nodes)
    if (nodeCount <= 0) call Galacticus_Error_Report('Merger_Tree_Construct_Fully_Specified','no nodes were specified')
    ! Create an array of nodes.
    allocate(nodeArray(nodeCount))
    call Memory_Usage_Record(sizeof(nodeArray))
    ! Iterate over nodes.
    do i=1,nodeCount
       ! Create the node.
       call thisTree%createNode(nodeArray(i)%node)
       ! Get the node definition.
       nodeDefinition => item(nodes,i-1)
       ! Assign an index to the node.
       call nodeArray(i)%node%indexSet(Node_Definition_Index(nodeDefinition,'index'))
    end do
    ! Initialize the tree root to null.
    thisTree%index    =  1
    thisTree%baseNode => null()
    ! Begin writing report.
    call Galacticus_Display_Indent('Initial conditions of fully-specified tree',verbosityInfo)
    ! Iterate over nodes.
    do i=1,nodeCount
       ! Get the node definition.
       nodeDefinition => item(nodes,i-1)
       ! Build parent pointers.
       indexValue=Node_Definition_Index(nodeDefinition,'parent'     )
       nodeArray(i)%node%parent     => Node_Lookup(nodeArray,indexValue)
       ! Build child pointers.
       indexValue=Node_Definition_Index(nodeDefinition,'firstChild' )
       nodeArray(i)%node%firstChild => Node_Lookup(nodeArray,indexValue)
       ! Build sibling pointers.
       indexValue=Node_Definition_Index(nodeDefinition,'sibling'    )
       nodeArray(i)%node%sibling    => Node_Lookup(nodeArray,indexValue)
       ! Assign the tree root node if this node has no parent.
       if (.not.associated(nodeArray(i)%node%parent)) then
          if (associated(thisTree%baseNode)) call Galacticus_Error_Report('Merger_Tree_Construct_Fully_Specified','multiple trees are not supported')
          thisTree%baseNode => nodeArray(i)%node
       end if
       ! Build components.
       call nodeArray(i)%node%componentBuilder(nodeDefinition)
       ! Dump the node.
       if (Galacticus_Verbosity_Level() > verbosityInfo) call nodeArray(i)%node%dump()
    end do   
    ! Finished - destroy the XML document.
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    ! Finish writing report.
    call Galacticus_Display_Unindent('done',verbosityInfo)
    ! Destroy the node array.
    deallocate(nodeArray)
    ! Check that we found a root node.
    if (.not.associated(thisTree%baseNode)) call Galacticus_Error_Report('Merger_Tree_Construct_Fully_Specified','no root node was found')
    return
  end subroutine Merger_Tree_Construct_Fully_Specified

  function Node_Definition_Index(nodeDefinition,indexType)
    !% Extract and return an index from a node definition as used when constructing fully-specified merger trees.
    use FoX_Dom
    use Galacticus_Error
    use Galacticus_Nodes
    use Kind_Numbers
    implicit none
    integer  (kind=kind_int8)                         :: Node_Definition_Index 
    type     (node          ), intent(in   ), pointer :: nodeDefinition        
    character(len=*         ), intent(in   )          :: indexType             
    type     (nodeList      )               , pointer :: indexElements         
    type     (node          )               , pointer :: indexElement          
    integer                                           :: indexValue            
    
    ! Find all matching tags.
    indexElements => getElementsByTagname(nodeDefinition,indexType)
    if (getLength(indexElements) < 1 .or. getLength(indexElements) > 1) call Galacticus_Error_Report('Node_Definition_Index','multiple indices specified')
    ! Get the index element.
    indexElement => item(indexElements,0)
    ! Extract the value.
    call extractDataContent(indexElement,indexValue)
    ! Transfer to function result.
    Node_Definition_Index=indexValue
    return
  end function Node_Definition_Index

  function Node_Lookup(nodeArray,indexValue) result (node)
    !% Find the position of a node in the {\tt nodeArray} array given its {\tt indexValue}.
    use Galacticus_Nodes
    use Kind_Numbers
    use Galacticus_Error
    implicit none
    type   (treeNode      ), pointer                     :: node       
    type   (treeNodeList  ), dimension(:), intent(in   ) :: nodeArray  
    integer(kind=kind_int8)              , intent(in   ) :: indexValue 
    integer                                              :: i          
    
    node => null()
    if (indexValue < 0_kind_int8) return
    do i=1,size(nodeArray)
       if (nodeArray(i)%node%index() == indexValue) then
          node => nodeArray(i)%node
          exit
       end if
    end do
    if (.not.associated(node)) call Galacticus_Error_Report('Node_Lookup','unable to find requested node')
    return
  end function Node_Lookup

end module Merger_Trees_Construct_Fully_Specified
