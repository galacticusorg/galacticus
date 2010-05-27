!% Contains a module which defines the tree node object and associated methods.

module Tree_Nodes
  !% Defines the tree node object and associated methods.
  use Components
  private
  public :: treeNode, treeNodeList

  type treeNode
     !% The tree node object type.
     integer,         private                   :: nodeIndex
     type(treeNode),  pointer                   :: parentNode,childNode,siblingNode,satelliteNode
     integer,         allocatable, dimension(:) :: componentIndex
     type(component), allocatable, dimension(:) :: components ! memoryManagementIgnore (force memory management system to ignore)
   contains
     ! Node indexing methods.
     !@ <objectMethods>
     !@   <object>treeNode</object>
     !@   <objectMethod>
     !@     <method>index</method>
     !@     <description>Return the index of the node (or $-1$ is the node does not exist).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>indexSet</method>
     !@     <description>Set the index of the node.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                  :: index                  => Tree_Node_Index
     procedure                                  :: indexSet               => Tree_Node_Index_Set
     ! Create/destroy methods.
     !@ <objectMethod>
     !@   <object>treeNode</object>
     !@   <method>destroy</method>
     !@   <description>Destroy the node and all components.</description>
     !@ </objectMethod>
     procedure                                  :: destroy                => Tree_Node_Destroy
     ! Component create/destroy methods.
     !@ <objectMethods>
     !@   <object>treeNode</object>
     !@   <objectMethod>
     !@     <method>componentExists(index)</method>
     !@     <description>Return true if the node has a component of the specified index.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>createComponent(index)</method>
     !@     <description>Create a component of the specified index if one does not already exist.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>destroyComponent(index)</method>
     !@     <description>Destroy the component of the specified index.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>destroyAllComponents</method>
     !@     <description>Destroy all components of the node.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                  :: componentExists        => Tree_Node_Component_Exists
     procedure                                  :: createComponent        => Tree_Node_Allocate_Component
     procedure                                  :: destroyComponent       => Tree_Node_Deallocate_Component
     procedure                                  :: destroyAllComponents   => Tree_Node_Deallocate_All_Components
     ! ! Component property get/set methods.
     ! !@ <objectMethods>
     ! !@   <object>treeNode</object>
     ! !@   <objectMethod>
     ! !@     <method>property(componentIndex,propertyIndex,emptyValue)</method>
     ! !@     <description>Returns the value of the {\tt propertyIndex} property of component {\tt componentIndex}, returning zero (or the optional {\tt emptyValue}) if the component does not exist.</description>
     ! !@   </objectMethod>
     ! !@   <objectMethod>
     ! !@     <method>data(componentIndex,dataIndex,emptyValue)</method>
     ! !@     <description>Returns the value of the {\tt dataIndex} data of component {\tt componentIndex}, returning zero (or the optional {\tt emptyValue}) if the component does not exist.</description>
     ! !@   </objectMethod>
     ! !@ </objectMethods>
     ! procedure :: property => Tree_Node_Component_Property_Get
     ! procedure :: data     => Tree_Node_Component_Data_Get
     ! Tree relational methods.
     !@ <objectMethods>
     !@   <object>treeNode</object>
     !@   <objectMethod>
     !@     <method>isPrimaryProgenitor</method>
     !@     <description>Return true if the node is the primary progenitor of its parent.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isProgenitorOf(descendentNode)</method>
     !@     <description>Return true if the node is a progenitor of {\tt descendentNode} (which may be specified as a node pointer or an index).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isPrimaryProgenitorOf(descendentNode)</method>
     !@     <description>Return true if the node is a primary progenitor of {\tt descendentNode} (which may be specified as a node pointer or an index).</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isOnMainBranch</method>
     !@     <description>Return true if the node is on the main branch of its tree.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                  ::                           Tree_Node_Is_Progenitor_Of_Index
     procedure                                  ::                           Tree_Node_Is_Progenitor_Of_Node 
     procedure                                  ::                           Tree_Node_Is_Primary_Progenitor_Of_Index
     procedure                                  ::                           Tree_Node_Is_Primary_Progenitor_Of_Node
     procedure                                  :: isPrimaryProgenitor    => Tree_Node_Is_Primary_Progenitor
     generic                                    :: isProgenitorOf         => Tree_Node_Is_Progenitor_Of_Index        ,&
          &                                                                  Tree_Node_Is_Progenitor_Of_Node
     generic                                    :: isPrimaryProgenitorOf  => Tree_Node_Is_Primary_Progenitor_Of_Index,&
          &                                                                  Tree_Node_Is_Primary_Progenitor_Of_Node
     procedure                                  :: isOnMainBranch         => Tree_Node_Is_On_Main_Branch
     ! Tree walk methods.
     !@ <objectMethods>
     !@   <object>treeNode</object>
     !@   <objectMethod>
     !@     <method>walkTree(nextNode)</method>
     !@     <description>Walks to the next node in a tree, returning the pointer in {\tt nextNode} if given, otherwise moving the target node to the new node.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>walkTreeWithSatellites(nextNode)</method>
     !@     <description>Walks to the next node in a tree, recursing through satellite nodes, returning the pointer in {\tt nextNode} if given, otherwise moving the target node to the new node.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>walkTreeConstruction(nextNode)</method>
     !@     <description>Walks to the next node in a tree in a manner suitable for tree construction, returning the pointer in {\tt nextNode} if given, otherwise moving the target node to the new node.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>walkBranch(startNode,nextNode)</method>
     !@     <description>Walks to the next node in the branch rooted at {\tt startNode}, returning the pointer in {\tt nextNode} if given, otherwise moving the target node to the new node.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                  ::                           Merger_Tree_Walk_Tree
     procedure                                  ::                           Merger_Tree_Walk_Tree_Same_Node
     generic                                    :: walkTree               => Merger_Tree_Walk_Tree,&
          &                                                                  Merger_Tree_Walk_Tree_Same_Node
     procedure                                  ::                           Merger_Tree_Walk_Tree_With_Satellites
     procedure                                  ::                           Merger_Tree_Walk_Tree_With_Satellites_Same_Node
     generic                                    :: walkTreeWithSatellites => Merger_Tree_Walk_Tree_With_Satellites,&
          &                                                                  Merger_Tree_Walk_Tree_With_Satellites_Same_Node
     procedure                                  ::                           Merger_Tree_Construction_Walk
     procedure                                  ::                           Merger_Tree_Construction_Walk_Same_Node
     generic                                    :: walkTreeConstruction   => Merger_Tree_Construction_Walk,&
          &                                                                  Merger_Tree_Construction_Walk_Same_Node
     procedure                                  :: walkBranch             => Merger_Tree_Walk_Branch
     ! Satellite methods.
     !@ <objectMethods>
     !@   <object>treeNode</object>
     !@   <objectMethod>
     !@     <method>isSatellite</method>
     !@     <description>Returns true if {\tt thisNode} is a satellite.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>removeFromHost</method>
     !@     <description>Removes {\tt thisNode} from its hosts list of satellites.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>lastSatellite</method>
     !@     <description>Returns a pointer to the last attached satellite of {\tt thisNode}.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                  :: isSatellite            => Tree_Node_Is_Satellite
     procedure                                  :: removeFromHost         => Satellite_Remove_from_Host
     procedure                                  :: lastSatellite          => Get_Last_Satellite
  end type treeNode

  type treeNodeList
     !% Type to give a list of treeNodes.
     type(treeNode), pointer :: node
  end type treeNodeList
  
contains

  !! Node indexing.

  integer function Tree_Node_Index(thisNode)
    !% Returns the index of {\tt thisNode}.
    implicit none
    type(treeNode), intent(in), pointer :: thisNode

    if (associated(thisNode)) then
       Tree_Node_Index=thisNode%nodeIndex
    else
       Tree_Node_Index=-1
    end if
    return
  end function Tree_Node_Index

  subroutine Tree_Node_Index_Set(thisNode,index)
    !% Set the index of {\tt thisNode}.
    implicit none
    type(treeNode), intent(inout) :: thisNode
    integer,        intent(in)    :: index

    thisNode%nodeIndex=index
    return
  end subroutine Tree_Node_Index_Set

  !! Creation and destruction.

  subroutine Tree_Node_Destroy(thisNode)
    !% Destroy a node in the tree, along with all components.
    use Memory_Management
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: iComponent

    ! Deallocate list of component indices.
    call Memory_Usage_Record(size(thisNode%componentIndex),addRemove=-4)
    if (allocated(thisNode%componentIndex)) deallocate(thisNode%componentIndex)

    ! Deallocate components.
    if (allocated(thisNode%components)) then
       do iComponent=1,size(thisNode%components)
          if (allocated(thisNode%components(iComponent)%properties)) then
             call Memory_Usage_Record(size(thisNode%components(iComponent)%properties),addRemove=-16)
             deallocate(thisNode%components(iComponent)%properties)
          end if
          if (allocated(thisNode%components(iComponent)%data)) then
             call Memory_Usage_Record(size(thisNode%components(iComponent)%data),addRemove=-8)
             deallocate(thisNode%components(iComponent)%data)
          end if
       end do
       deallocate(thisNode%components)
    end if

    ! Deallocate the tree node object.
    deallocate(thisNode)
    thisNode => null()

    return
  end subroutine Tree_Node_Destroy

  !! Component creation and destruction methods.

  logical function Tree_Node_Component_Exists(thisNode,componentIndex)
    !% Return true if {\tt thisNode} already has a component with index {\tt componentIndex}.
    implicit none
    type(treeNode), pointer, intent(in) :: thisNode
    integer,                 intent(in) :: componentIndex

    Tree_Node_Component_Exists=(thisNode%componentIndex(componentIndex) > 0)
    return
  end function Tree_Node_Component_Exists

  subroutine Tree_Node_Allocate_Component(thisNode,componentIndex,propertyCount,dataCount,historyCount)
    !% Ensure that the component array is allocated with sufficient size.
    use Memory_Management
    implicit none
    integer,         intent(in)                :: componentIndex,propertyCount,dataCount,historyCount
    type(treeNode),  intent(inout)             :: thisNode
    type(component), allocatable, dimension(:) :: tempComponents ! memoryManagementIgnore (force memory management system to ignore)
    integer                                    :: previousSize,thisIndex

    if (thisNode%componentIndex(componentIndex) == -1) then
       if (allocated(thisNode%components)) then
          previousSize=size(thisNode%components)
          call Move_Alloc(thisNode%components,tempComponents)
          allocate(thisNode%components(previousSize+1))
          thisNode%components(1:previousSize)=tempComponents
          deallocate(tempComponents)
       else
          allocate(thisNode%components(1))
       end if
       thisNode%componentIndex(componentIndex)=size(thisNode%components)
       thisIndex=size(thisNode%components)
       thisNode%components(thisIndex)%nextComponentOfType => null()
       if (.not.allocated(thisNode%components(thisIndex)%properties)) then
          if (propertyCount > 0) then
             allocate(thisNode%components(thisIndex)%properties(propertyCount,2))
             call Memory_Usage_Record(propertyCount,addRemove=+16)
             thisNode%components(thisIndex)%properties=0.0d0
          end if
          if (dataCount > 0) then
             allocate(thisNode%components(thisIndex)%data(dataCount))
             call Memory_Usage_Record(dataCount,addRemove=+8)
             thisNode%components(thisIndex)%data=0.0d0
          end if
          if (historyCount > 0) allocate(thisNode%components(thisIndex)%histories(historyCount))
       end if
    end if

    return
  end subroutine Tree_Node_Allocate_Component

  subroutine Tree_Node_Deallocate_Component(thisNode,componentIndex,propertyCount,dataCount)
    !% Ensure that the component array is deallocated.
    use Memory_Management
    implicit none
    integer,         intent(in)                :: componentIndex,propertyCount,dataCount
    type(treeNode),  intent(inout)             :: thisNode
    type(component), allocatable, dimension(:) :: tempComponents ! memoryManagementIgnore (force memory management system to ignore)
    integer                                    :: previousSize,listIndex,timesCount,historyCount,iHistory

    listIndex=thisNode%componentIndex(componentIndex)
    if (listIndex /= -1) then
       ! Count memory that will be freed.
       call Memory_Usage_Record(propertyCount,addRemove=-16)
       call Memory_Usage_Record(dataCount    ,addRemove=-8 )
       do iHistory=1,size(thisNode%components(listIndex)%histories)
          if (allocated(thisNode%components(listIndex)%histories(iHistory)%time)) then
             timesCount  =size(thisNode%components(listIndex)%histories(iHistory)%time      )
             historyCount=size(thisNode%components(listIndex)%histories(iHistory)%data,dim=2)
             call Memory_Usage_Record(timesCount*(1+2*historyCount),addRemove=-8)
          end if
       end do
       previousSize=size(thisNode%components)
       if (previousSize > 1) then
          call Move_Alloc(thisNode%components,tempComponents)
          allocate(thisNode%components(previousSize-1))
          if (listIndex > 1           ) thisNode%components(1             :listIndex-1)=tempComponents(1         &
               &  :listIndex-1)
          if (listIndex < previousSize) thisNode%components(listIndex:previousSize-1  )=tempComponents(listIndex &
               &+1:previousSize    )
          deallocate(tempComponents)
          ! Shift component indices of remaining components by 1.
          where (thisNode%componentIndex > listIndex)
             thisNode%componentIndex=thisNode%componentIndex-1
          end where
       else
          deallocate(thisNode%components)
       end if
       thisNode%componentIndex(componentIndex)=-1
    end if
    return
  end subroutine Tree_Node_Deallocate_Component

  subroutine Tree_Node_Deallocate_All_Components(thisNode)
    !% Ensure that all components of {\tt thisNode} are deallocated.
    use Memory_Management
    implicit none
    type(treeNode),  intent(inout)             :: thisNode
    integer                                    :: listIndex,thisIndex

    ! Deallocate each component and record the memory usage change.
    do thisIndex=1,size(thisNode%componentIndex)
       listIndex=thisNode%componentIndex(thisIndex)
       if (listIndex /= -1) then
          if (allocated(thisNode%components(listIndex)%properties)) then
             call Memory_Usage_Record(size(thisNode%components(listIndex)%properties,dim=1),addRemove=-16)
             deallocate(thisNode%components(listIndex)%properties)
          end if
          if (allocated(thisNode%components(listIndex)%data)) then
             call Memory_Usage_Record(size(thisNode%components(listIndex)%data),addRemove=-8)
             deallocate(thisNode%components(listIndex)%data)
          end if
       end if
    end do
    ! Deallocate the components array.
    deallocate(thisNode%components)
    ! Reset the component index array.
    thisNode%componentIndex=-1
    return
  end subroutine Tree_Node_Deallocate_All_Components

  !! (Tree) relational methods.

  logical function Tree_Node_Is_Primary_Progenitor(thisNode)
    !% Returns true if {\tt thisNode} is the primary progenitor of its parent node.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    if (associated(thisNode%parentNode)) then
       Tree_Node_Is_Primary_Progenitor=associated(thisNode%parentNode%childNode,thisNode)
    else
       Tree_Node_Is_Primary_Progenitor=.false.
    end if
    return
  end function Tree_Node_Is_Primary_Progenitor

  logical function Tree_Node_Is_Primary_Progenitor_Of_Index(thisNode,targetNodeIndex)
    !% Return true if {\tt thisNode} is a progenitor of the node with index {\tt targetNodeIndex}.
    implicit none
    type(treeNode), intent(in), pointer :: thisNode
    integer,        intent(in)          :: targetNodeIndex
    type(treeNode),             pointer :: workNode

    Tree_Node_Is_Primary_Progenitor_Of_Index=.false.
    workNode => thisNode
    do while (associated(workNode))
       if (workNode%index() == targetNodeIndex) then
          Tree_Node_Is_Primary_Progenitor_Of_Index=.true.
          return
       end if
       if (.not.workNode%isPrimaryProgenitor()) return
       workNode => workNode%parentNode
    end do
    return
  end function Tree_Node_Is_Primary_Progenitor_Of_Index

  logical function Tree_Node_Is_Primary_Progenitor_Of_Node(thisNode,targetNode)
    !% Return true if {\tt thisNode} is a progenitor of {\tt targetNode}.
    implicit none
    type(treeNode), intent(in), pointer :: thisNode,targetNode
    type(treeNode),             pointer :: workNode

    Tree_Node_Is_Primary_Progenitor_Of_Node=.false.
    workNode => thisNode
    do while (associated(workNode))
       if (associated(workNode,targetNode)) then
          Tree_Node_Is_Primary_Progenitor_Of_Node=.true.
          return
       end if
       if (.not.workNode%isPrimaryProgenitor()) return
       workNode => workNode%parentNode
    end do
    return
  end function Tree_Node_Is_Primary_Progenitor_Of_Node

  logical function Tree_Node_Is_Progenitor_Of_Index(thisNode,targetNodeIndex)
    !% Return true if {\tt thisNode} is a progenitor of the node with index {\tt targetNodeIndex}.
    implicit none
    type(treeNode), intent(in), pointer :: thisNode
    integer,        intent(in)          :: targetNodeIndex
    type(treeNode),             pointer :: workNode

    Tree_Node_Is_Progenitor_Of_Index=.false.
    workNode => thisNode
    do while (associated(workNode))
       if (workNode%index() == targetNodeIndex) then
          Tree_Node_Is_Progenitor_Of_Index=.true.
          return
       end if
       workNode => workNode%parentNode
    end do
    return
  end function Tree_Node_Is_Progenitor_Of_Index

  logical function Tree_Node_Is_Progenitor_Of_Node(thisNode,targetNode)
    !% Return true if {\tt thisNode} is a progenitor of {\tt targetNode}.
    implicit none
    type(treeNode), intent(in), pointer :: thisNode,targetNode
    type(treeNode),             pointer :: workNode

    Tree_Node_Is_Progenitor_Of_Node=.false.
    workNode => thisNode
    do while (associated(workNode))
       if (associated(workNode,targetNode)) then
          Tree_Node_Is_Progenitor_Of_Node=.true.
          return
       end if
       workNode => workNode%parentNode
    end do
    return
  end function Tree_Node_Is_Progenitor_Of_Node

  logical function Tree_Node_Is_On_Main_Branch(thisNode)
    !% Returns true if {\tt thisNode} is on the main branch.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    type(treeNode),                pointer :: workNode

    Tree_Node_Is_On_Main_Branch=.not.associated(thisNode%parentNode)
    workNode => thisNode
    do while (associated(workNode%parentNode))
       if (.not.workNode%isPrimaryProgenitor()) return
       workNode => workNode%parentNode
    end do
    Tree_Node_Is_On_Main_Branch=.true.  
    return
  end function Tree_Node_Is_On_Main_Branch

  !! Satellite methods.

  logical function Tree_Node_Is_Satellite(thisNode)
    !% Returns true if {\tt thisNode} is a satellite.
    implicit none
    type(treeNode), pointer, intent(in) :: thisNode
    type(treeNode), pointer             :: parentNode,childNode

    parentNode => thisNode%parentNode
    select case (associated(parentNode))
    case (.false.)
       Tree_Node_Is_Satellite=.false.
       return
    case (.true.)
       childNode => parentNode%childNode
       Tree_Node_Is_Satellite=.true.
       do while (associated(childNode))
          if (associated(childNode,thisNode)) then
             Tree_Node_Is_Satellite=.false.
             exit
          end if
          childNode => childNode%siblingNode
       end do
    end select
    return
  end function Tree_Node_Is_Satellite

  subroutine Satellite_Remove_from_Host(satelliteNode)
    !% Remove {\tt satelliteNode} from the linked list of its host node's satellites.
    use Galacticus_Display
    use ISO_Varying_String
    use String_Handling
    implicit none
    type(treeNode),      pointer, intent(in) :: satelliteNode
    type(treeNode),      pointer             :: hostNode,thisNode,previousNode
    type(varying_string)                     :: message

    hostNode => satelliteNode%parentNode
    message='Satellite node ['
    message=message//satelliteNode%index()//'] being removed from host node ['//hostNode%index()//']'
    call Galacticus_Display_Message(message,verbosityInfo)
    if (associated(hostNode%satelliteNode,satelliteNode)) then
       ! This is the first satellite, unlink it, and link to any sibling.
       hostNode%satelliteNode => satelliteNode%siblingNode
    else
       thisNode     => hostNode%satelliteNode
       previousNode => null()
       do while (associated(thisNode))
          if (associated(thisNode,satelliteNode)) then
             ! Found our node, link its older sibling to its younger sibling.
             previousNode%siblingNode => thisNode%siblingNode
             exit
          end if
          previousNode => thisNode
          thisNode => thisNode%siblingNode
       end do
    end if
    return
  end subroutine Satellite_Remove_from_Host

  subroutine Get_Last_Satellite(thisNode,satelliteNode)
    !% Returns a pointer to the final satellite node associated with {\tt thisNode}.
    implicit none
    type(treeNode), intent(in),    pointer :: thisNode
    type(treeNode), intent(inout), pointer :: satelliteNode

    satelliteNode => thisNode%satelliteNode
    do while (associated(satelliteNode%siblingNode))
       satelliteNode => satelliteNode%siblingNode
    end do
    return
  end subroutine Get_Last_Satellite

  !! Tree-walking methods

  subroutine Merger_Tree_Walk_Tree_Same_Node(thisNode)
    !% Simple interface to the \hyperlink{objects.tree_node.F90:tree_nodes:merger_tree_walk_tree}{{\tt
    !% Merger\_Tree\_Walk\_Tree()}} subroutine that does the walk in place, i.e. by updating the input tree node pointer to point
    !% to the next node.
    implicit none
    type (treeNode), pointer, intent(inout) :: thisNode

    call Merger_Tree_Walk_Tree(thisNode,thisNode)
    return
  end subroutine Merger_Tree_Walk_Tree_Same_Node

  subroutine Merger_Tree_Walk_Tree(thisNode,nextNode)
    !% This function provides a mechanism for walking through an entire merger tree. Given a pointer {\tt thisNode}
    !% to a node of the tree, it will return the next node that should be visited in the tree. Thus, if {\tt thisNode} is
    !% initially set to the base of the merger tree and {\tt Merger\_Tree\_Walk()} is called repeatedly it will walk through every node
    !% of the tree. Once the entire tree has been walked, a {\tt null()} pointer will be returned, indicating that there
    !% are no more nodes to walk. Each node will be visited once and once only if the tree is walked in this way.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode,nextNode
    type(treeNode), pointer                :: workNode

    workNode => thisNode

    if (.not.associated(workNode%parentNode)) then
       ! This is the base of the merger tree.
       do while (associated(workNode%childNode))
          workNode => workNode%childNode
       end do
       if (associated(workNode,thisNode)) nullify(workNode)
    else
       if (associated(workNode%siblingNode)) then
          workNode => workNode%siblingNode
          do while (associated(workNode%childNode))
             workNode => workNode%childNode
          end do
       else
          workNode => workNode%parentNode
          if (.not.associated(workNode%parentNode)) workNode => null() ! Terminate when back at tree base.
       end if
    end if
    nextNode => workNode
    return
  end subroutine Merger_Tree_Walk_Tree

  subroutine Merger_Tree_Walk_Tree_With_Satellites_Same_Node(thisNode)
    !% Simple interface to the \hyperlink{objects.tree_node.F90:tree_nodes:merger_tree_walk_tree_with_satellites}{{\tt
    !% Merger\_Tree\_Walk\_Tree\_With \_Satelites()}} subroutine that does the walk in place, i.e. by updating the input tree node
    !% pointer to point to the next node.
    implicit none
    type (treeNode), pointer :: thisNode

    call Merger_Tree_Walk_Tree_With_Satellites(thisNode,thisNode)
    return
  end subroutine Merger_Tree_Walk_Tree_With_Satellites_Same_Node

  subroutine Merger_Tree_Walk_Tree_With_Satellites(thisNode,nextNode)
    !% Merger tree walk function which also descends through satellite nodes.
    implicit none
    type (treeNode), pointer, intent(inout) :: thisNode,nextNode
    type (treeNode), pointer                :: workNode

    workNode => thisNode
    if (.not.associated(workNode%parentNode)) then
       ! This is the base of the merger tree.
       do while (associated(workNode%childNode))
          workNode => workNode%childNode
       end do
       ! Descend through any satellite nodes.
       do while (associated(workNode%satelliteNode))
          workNode => workNode%satelliteNode
       end do
       if (associated(workNode,thisNode)) nullify(workNode)
    else
       if (associated(workNode%siblingNode)) then
          workNode => workNode%siblingNode
          do while (associated(workNode%childNode))
             workNode => workNode%childNode
          end do
          ! Descend through any satellite nodes.
          do while (associated(workNode%satelliteNode))
             workNode => workNode%satelliteNode
          end do
       else
          workNode => workNode%parentNode
          if (.not.associated(workNode%parentNode)) workNode => null() ! Terminate when back at tree base.
       end if
    end if
    nextNode => workNode
    return
  end subroutine Merger_Tree_Walk_Tree_With_Satellites

  subroutine Merger_Tree_Walk_Branch(thisNode,startNode,nextNode)
    !% This function provides a mechanism for walking through the branches of the merger tree. Given a pointer {\tt thisNode}
    !% to a branch of the tree, it will return the next node that should be visited in the tree. Thus, if {\tt thisNode} is
    !% initially set to the base of the merger tree and {\tt Merger\_Tree\_Walk\_Branch()} is called repeatedly it will walk through every node
    !% of the branch. Once the entire branch has been walked, a {\tt null()} pointer will be returned, indicating that there
    !% are no more nodes to walk. Each node will be visited once and once only if the branch is walked in this way.
    implicit none
    type (treeNode), pointer, intent(inout) :: thisNode,startNode,nextNode
    type (treeNode), pointer                :: workNode

    workNode => thisNode

    if (associated(thisNode,startNode)) then
       do while (associated(workNode%childNode))
          workNode => workNode%childNode
       end do
       if (associated(workNode,thisNode)) nullify(workNode)
    else
       if (associated(workNode%siblingNode)) then
          workNode => workNode%siblingNode
          do while (associated(workNode%childNode))
             workNode => workNode%childNode
          end do
       else
          workNode => workNode%parentNode
          if (associated(workNode,startNode)) workNode => null() ! Terminate when back at starting node.
       end if
    end if
    nextNode => workNode
    return
  end subroutine Merger_Tree_Walk_Branch

  subroutine Merger_Tree_Construction_Walk_Same_Node(thisNode)
    !% Simple interface to the \hyperlink{objects.tree_node.F90:tree_nodes:merger_tree_construction_walk}{{\tt Merger\_Tree\_Construction\_Walk()}} subroutine
    !% that does the walk in place, i.e. by updating the input tree node pointer to point to the next node.
    implicit none
    type (treeNode), pointer, intent(inout) :: thisNode

    call Merger_Tree_Construction_Walk(thisNode,thisNode)
    return
  end subroutine Merger_Tree_Construction_Walk_Same_Node

  subroutine Merger_Tree_Construction_Walk(thisNode,nextNode)
    !% This function provides a mechanism for walking through a merger tree that is being built.
    implicit none
    type (treeNode), pointer, intent(inout) :: thisNode,nextNode
    type (treeNode), pointer                :: workNode

    workNode => thisNode

    if (associated(workNode%childNode)) then
       ! Move to the primary child if one exists.
       do while (associated(workNode%childNode))
          workNode => workNode%childNode
       end do
    else
       if (associated(workNode%siblingNode)) then
          workNode => workNode%siblingNode
          do while (associated(workNode%childNode))
             workNode => workNode%childNode
          end do
       else
          do while (associated(workNode))
             if (associated(workNode%parentNode)) then
                workNode => workNode%parentNode
                if (associated(workNode%siblingNode)) then
                   workNode => workNode%siblingNode
                   do while (associated(workNode%childNode))
                      workNode => workNode%childNode
                   end do
                   exit
                end if
             else
                workNode => null()
             end if
          end do
       end if
    end if
    nextNode => workNode
    return
  end subroutine Merger_Tree_Construction_Walk

  !! Component property get methods.

!!! Currently not used as gfortran 4.4 has a bug when using type-bound procedures with optional argumentds. Should be fixed in <gfortran 4.5>.

  ! double precision function Tree_Node_Component_Property_Get(thisNode,componentIndex,propertyIndex,emptyValue)
  !   !% Returns the value of a property identified by index {\tt propertyIndex} in the component identified by {\tt
  !   !% componentIndex}. If the component does not exist then zero (or the optional {\tt emptyValue}) is returned instead.
  !   implicit none
  !   type(treeNode),   intent(inout), pointer  :: thisNode
  !   integer,          intent(in)              :: componentIndex,propertyIndex
  !   double precision, intent(in),    optional :: emptyValue
  !   integer                                   :: thisIndex

  !   if (thisNode%componentExists(componentIndex)) then
  !      thisIndex=thisNode%componentIndex(componentIndex)
  !      Tree_Node_Component_Property_Get=thisNode%components(thisIndex)%properties(propertyIndex,propertyValue)
  !   else
  !      if (present(emptyValue)) then
  !         Tree_Node_Component_Property_Get=emptyValue
  !      else
  !         Tree_Node_Component_Property_Get=0.0d0
  !      end if
  !   end if
  !   return
  ! end function Tree_Node_Component_Property_Get

  ! double precision function Tree_Node_Component_Data_Get(thisNode,componentIndex,dataIndex,emptyValue)
  !   !% Returns the value of a data identified by index {\tt dataIndex} in the component identified by {\tt
  !   !% componentIndex}. If the component does not exist then zero (or the optional {\tt emptyValue}) is returned instead.
  !   implicit none
  !   type(treeNode),   intent(inout), pointer  :: thisNode
  !   integer,          intent(in)              :: componentIndex,dataIndex
  !   double precision, intent(in),    optional :: emptyValue
  !   integer                                   :: thisIndex

  !   if (thisNode%componentExists(componentIndex)) then
  !      thisIndex=thisNode%componentIndex(componentIndex)
  !      Tree_Node_Component_Data_Get=thisNode%components(thisIndex)%data(dataIndex)
  !   else
  !      if (present(emptyValue)) then
  !         Tree_Node_Component_Data_Get=emptyValue
  !      else
  !         Tree_Node_Component_Data_Get=0.0d0
  !      end if
  !   end if
  !   return
  ! end function Tree_Node_Component_Data_Get

end module Tree_Nodes
