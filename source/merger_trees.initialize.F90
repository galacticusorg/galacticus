!% Contains a module which implements initialization of merger tree structures.

module Merger_Trees_Initialize
  !% Implements initialization of merger tree structures.
  private
  public :: Merger_Tree_Initialize

contains

  subroutine Merger_Tree_Initialize(thisTree)
    !% Walk through all nodes of a tree and call any routines that requested to perform initialization tasks.
    use Merger_Trees
    use Tree_Nodes
    !# <include directive="mergerTreeInitializeTask" type="moduleUse">
    include 'merger_trees.initialize.tasks.modules.inc'
    !# </include>
    implicit none
    type(mergerTree), intent(inout) :: thisTree
    type(treeNode),   pointer       :: thisNode

    if (.not.thisTree%initialized) then
       thisNode => thisTree%baseNode
       do while (associated(thisNode))

          ! Call subroutines to perform any necessary initialization of this node.
          !# <include directive="mergerTreeInitializeTask" type="code" action="subroutine">
          !#  <subroutineArgs>thisNode</subroutineArgs>
          include 'merger_trees.initialize.tasks.inc'
          !# </include>

          call thisNode%walkTree()
       end do
       thisTree%initialized=.true.
    end if

    return
  end subroutine Merger_Tree_Initialize

end module Merger_Trees_Initialize
