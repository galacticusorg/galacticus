!% Contains a module which implements dumping of a node's data to screen.

module Tree_Nodes_Dump
  !% Implements dumping of a node's data to screen.
  private
  public :: Node_Dump
  
contains

  subroutine Node_Dump(thisNode)
    !% Dump the properties of {\tt thisNode} to screen.
    use Tree_Nodes
    !# <include directive="nodeDumpTask" type="moduleUse">
    include 'objects.tree_node.dump.moduleUse.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    write (0,'(a,i6,a)'   ) 'Dumping data for node [',thisNode%index(),']'
    write (0,'(1x,a20,i6)') 'parent node: ',thisNode%parentNode%index()
    write (0,'(1x,a20,i6)') 'sibling node: ',thisNode%siblingNode%index()
    write (0,'(1x,a20,i6)') 'child node: ',thisNode%childNode%index()
    !# <include directive="nodeDumpTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode</subroutineArgs>
    include 'objects.tree_node.dump.inc'
    !# </include>
    
    return
  end subroutine Node_Dump

end module Tree_Nodes_Dump
