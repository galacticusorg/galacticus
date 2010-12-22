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


!% Contains a module which defines the tree node object and associated methods.

module Tree_Nodes
  !% Defines the tree node object and associated methods.
  use Components
  use Kind_Numbers
  use Histories
  private
  public :: treeNode, treeNodeList, Tree_Node_Rate_Rate_Compute_Dummy, Tree_Node_Rate_Adjust_Dummy,&
       & Tree_Node_Rate_Adjust_Array_Dummy, Tree_Node_Rate_Adjust_History_Dummy, Get_Template, Set_Template

  type treeNode
     !% The tree node object type.
     integer(kind=kind_int8), private                   :: nodeIndex,nodeUniqueID
     type(treeNode),          pointer                   :: parentNode,childNode,siblingNode,satelliteNode,mergeNode,mergeeNode,nextMergee
     integer,                 allocatable, dimension(:) :: componentIndex
     type(component),         allocatable, dimension(:) :: components
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
     !@   <objectMethod>
     !@     <method>uniqueID</method>
     !@     <description>Return the unique ID of the node.\label{method:uniqueID}</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>uniqueIDSet</method>
     !@     <description>Set the unique ID of the node.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                  :: index                  => Tree_Node_Index
     procedure                                  :: indexSet               => Tree_Node_Index_Set
     procedure                                  :: uniqueID               => Tree_Node_Unique_ID
     procedure                                  :: uniqueIDSet            => Tree_Node_Unique_ID_Set
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
     !@   <objectMethod>
     !@     <method>mergesWith</method>
     !@     <description>Return a pointer to the node with which this node merges.</description>
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
     procedure                                  :: mergesWith             => Tree_Node_Merge_Node

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
     !@   <objectMethod>
     !@     <method>walkBranchWithSatellites(startNode,nextNode)</method>
     !@     <description>Walks to the next node in the branch rooted at {\tt startNode}, recursing through satellite nodes, returning the pointer in {\tt nextNode} if given, otherwise moving the target node to the new node.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                  ::                             Merger_Tree_Walk_Tree
! <gfortran 4.6> not sure we can support the following since class dummy argument can not be pointer.
!     procedure                                  ::                             Merger_Tree_Walk_Tree_Same_Node
     generic                                    :: walkTree                 => Merger_Tree_Walk_Tree!,&
!          &                                                                    Merger_Tree_Walk_Tree_Same_Node
     procedure                                  ::                             Merger_Tree_Walk_Tree_With_Satellites
! <gfortran 4.6> not sure we can support the following since class dummy argument can not be pointer.
!     procedure                                  ::                             Merger_Tree_Walk_Tree_With_Satellites_Same_Node
     generic                                    :: walkTreeWithSatellites   => Merger_Tree_Walk_Tree_With_Satellites!,&
!          &                                                                    Merger_Tree_Walk_Tree_With_Satellites_Same_Node
     procedure                                  ::                             Merger_Tree_Construction_Walk
! <gfortran 4.6> not sure we can support the following since class dummy argument can not be pointer.
!     procedure                                  ::                             Merger_Tree_Construction_Walk_Same_Node
     generic                                    :: walkTreeConstruction     => Merger_Tree_Construction_Walk!,&
!          &                                                                    Merger_Tree_Construction_Walk_Same_Node
     procedure                                  :: walkBranch               => Merger_Tree_Walk_Branch
     procedure                                  :: walkBranchWithSatellites => Merger_Tree_Walk_Branch_With_Satellites
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
     !@     <method>removeFromMergee</method>
     !@     <description>Removes {\tt thisNode} from its merge target's list of mergees.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>lastSatellite</method>
     !@     <description>Returns a pointer to the last attached satellite of {\tt thisNode}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>earliestProgenitor</method>
     !@     <description>Returns a pointer to the earliest progenitor (along the main branch) of {\tt thisNode}.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                  :: isSatellite            => Tree_Node_Is_Satellite
     procedure                                  :: removeFromHost         => Satellite_Remove_from_Host
     procedure                                  :: removeFromMergee       => Satellite_Remove_from_Mergee
     procedure                                  :: lastSatellite          => Get_Last_Satellite
     procedure                                  :: earliestProgenitor     => Get_Earliest_Progenitor
  end type treeNode

  type treeNodeList
     !% Type to give a list of treeNodes.
     type(treeNode), pointer :: node
  end type treeNodeList

  ! Procedure pointers go here.
  !# <include directive="treeNodeMethodsPointer" type="methods">
  include 'objects.tree_node.methods.inc'
  !# </include>

  ! Pipe procedure pointers go here.
  !# <include directive="treeNodePipePointer" type="pipe">
  include 'objects.tree_node.pipes.inc'
  !# </include>

  abstract interface
     double precision function Get_Template(thisNode)
       import treeNode
       type(treeNode), pointer, intent(inout) :: thisNode
     end function Get_Template
  end interface
  abstract interface
     subroutine Get_Template_Array(thisNode,values)
       import treeNode
       type(treeNode),   pointer,      intent(inout) :: thisNode
       double precision, dimension(:), intent(out)   :: values
     end subroutine Get_Template_Array
  end interface
  abstract interface
     function Get_Template_History(thisNode)
       import treeNode, history
       type(treeNode), pointer, intent(inout) :: thisNode
       type(history)                          :: Get_Template_History
     end function Get_Template_History
  end interface
  abstract interface
     subroutine Set_Template(thisNode,value)
       import treeNode
       type(treeNode),   pointer, intent(inout) :: thisNode
       double precision,          intent(in)    :: value
     end subroutine Set_Template
  end interface
  abstract interface
     subroutine Set_Template_Array(thisNode,values)
       import treeNode
       type(treeNode),   pointer, intent(inout) :: thisNode
       double precision,          intent(in)    :: values(:)
     end subroutine Set_Template_Array
  end interface
  abstract interface
     subroutine Set_Template_History(thisNode,values)
       import treeNode, history
       type(treeNode), pointer, intent(inout) :: thisNode
       type(history),           intent(in)    :: values
     end subroutine Set_Template_History
  end interface
  abstract interface
     subroutine Rate_Adjust_Template(thisNode,interrupt,interruptProcedure,rateAdjustment)
       import treeNode
       type(treeNode),  pointer, intent(inout) :: thisNode
       logical,                  intent(inout) :: interrupt
       procedure(),     pointer, intent(inout) :: interruptProcedure
       double precision,         intent(in)    :: rateAdjustment
     end subroutine Rate_Adjust_Template
  end interface
  abstract interface
     subroutine Rate_Adjust_Template_Array(thisNode,interrupt,interruptProcedure,rateAdjustments)
       import treeNode
       type(treeNode),  pointer, intent(inout) :: thisNode
       logical,                  intent(inout) :: interrupt
       procedure(),     pointer, intent(inout) :: interruptProcedure
       double precision,         intent(in)    :: rateAdjustments(:)
     end subroutine Rate_Adjust_Template_Array
  end interface
  abstract interface
     subroutine Rate_Adjust_Template_History(thisNode,interrupt,interruptProcedure,rateAdjustments)
       import treeNode, history
       type(treeNode),  pointer, intent(inout) :: thisNode
       logical,                  intent(inout) :: interrupt
       procedure(),     pointer, intent(inout) :: interruptProcedure
       type(history),            intent(in)    :: rateAdjustments
     end subroutine Rate_Adjust_Template_History
  end interface
  abstract interface
     subroutine Rate_Compute_Template(thisNode,interrupt,interruptProcedure)
       import treeNode
       type(treeNode), pointer, intent(inout) :: thisNode
       logical,                 intent(inout) :: interrupt
       procedure(),    pointer, intent(inout) :: interruptProcedure
    end subroutine Rate_Compute_Template
  end interface

contains

  integer(kind=kind_int8) function Tree_Node_Unique_ID(thisNode)
    !% Returns the unique ID of {\tt thisNode}.
    implicit none
#ifdef GCC45
    class(treeNode), intent(in), target  :: thisNode
#else
    type(treeNode),  intent(in), pointer :: thisNode
#endif
    type(treeNode),              pointer :: workNode

#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       workNode => thisNode
#ifdef GCC45
    end select
#endif
    if (associated(workNode)) then
       Tree_Node_Unique_ID=workNode%nodeUniqueID
    else
       Tree_Node_Unique_ID=-1
    end if
    return
  end function Tree_Node_Unique_ID

  subroutine Tree_Node_Unique_ID_Set(thisNode,uniqueID)
    !% Set the index of {\tt thisNode}.
    implicit none
#ifdef GCC45
    class(treeNode),         intent(inout) :: thisNode
#else
    type(treeNode),          intent(inout) :: thisNode
#endif
    integer(kind=kind_int8), intent(in)    :: uniqueID

    thisNode%nodeUniqueID=uniqueID
    return
  end subroutine Tree_Node_Unique_ID_Set

  integer(kind=kind_int8) function Tree_Node_Index(thisNode)
    !% Returns the index of {\tt thisNode}.
    implicit none
#ifdef GCC45
    class(treeNode), intent(in), target  :: thisNode
#else
    type(treeNode),  intent(in), pointer :: thisNode
#endif
    type(treeNode),              pointer :: workNode

#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       workNode => thisNode
#ifdef GCC45
    end select
#endif
    if (associated(workNode)) then
       Tree_Node_Index=workNode%nodeIndex
    else
       Tree_Node_Index=-1
    end if
    return
  end function Tree_Node_Index

  subroutine Tree_Node_Index_Set(thisNode,index)
    !% Set the index of {\tt thisNode}.
    implicit none
#ifdef GCC45
    class(treeNode),         intent(inout) :: thisNode
#else
    type(treeNode),          intent(inout) :: thisNode
#endif
    integer(kind=kind_int8), intent(in)    :: index

    thisNode%nodeIndex=index
    return
  end subroutine Tree_Node_Index_Set

  subroutine Tree_Node_Destroy(thisNode)
    !% Destroy a node in the tree, along with all components.
    use Memory_Management
    implicit none
#ifdef GCC45
    class(treeNode), target,  intent(inout) :: thisNode
#else
    type(treeNode),  pointer, intent(inout) :: thisNode
#endif
    integer                                 :: iComponent
    type(treeNode),  pointer                :: workNode

#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       workNode => thisNode
#ifdef GCC45
    end select
#endif

    ! Deallocate list of component indices.
    if (allocated(workNode%componentIndex)) then
       call Memory_Usage_Record(sizeof(workNode%componentIndex),addRemove=-1,memoryType=memoryTypeNodes)
       deallocate(workNode%componentIndex)
    end if

    ! Deallocate components.
    if (allocated(workNode%components)) then
       do iComponent=1,size(workNode%components)
          if (allocated(workNode%components(iComponent)%properties)) then
             call Memory_Usage_Record(sizeof(workNode%components(iComponent)%properties),addRemove=-1,memoryType=memoryTypeNodes)
             deallocate(workNode%components(iComponent)%properties)
          end if
          if (allocated(workNode%components(iComponent)%data)) then
             call Memory_Usage_Record(sizeof(workNode%components(iComponent)%data),addRemove=-1,memoryType=memoryTypeNodes)
             deallocate(workNode%components(iComponent)%data)
          end if
       end do
       deallocate(workNode%components)
       call Memory_Usage_Record(sizeof(workNode%components),addRemove=-1,memoryType=memoryTypeNodes)
    end if

    ! Deallocate the tree node object.
    call Memory_Usage_Record(sizeof(workNode),addRemove=-1,memoryType=memoryTypeNodes)
    deallocate(workNode)
#ifndef GCC45
    ! <gfortran 4.6> This next line does not do anything under gFortran 4.6 since the input object is not a pointer (merely a target).
    workNode => null()
#endif

    return
  end subroutine Tree_Node_Destroy

  logical function Tree_Node_Component_Exists(thisNode,componentIndex)
    !% Return true if {\tt thisNode} already has a component with index {\tt componentIndex}.
    implicit none
#ifdef GCC45
    class(treeNode), intent(in)          :: thisNode
#else
    type(treeNode),  intent(in), pointer :: thisNode
#endif
    integer,         intent(in)          :: componentIndex

    Tree_Node_Component_Exists=(thisNode%componentIndex(componentIndex) > 0)
    return
  end function Tree_Node_Component_Exists

  subroutine Tree_Node_Allocate_Component(thisNode,componentIndex,propertyCount,dataCount,historyCount)
    !% Ensure that the component array is allocated with sufficient size.
    use Memory_Management
    implicit none
    integer,          intent(in)                :: componentIndex,propertyCount,dataCount,historyCount
#ifdef GCC45
    class(treeNode),  intent(inout)             :: thisNode
#else
    type(treeNode),   intent(inout)             :: thisNode
#endif
    type(component),  allocatable, dimension(:) :: tempComponents
    integer                                     :: previousSize,thisIndex

    if (thisNode%componentIndex(componentIndex) == -1) then
       if (allocated(thisNode%components)) then
          previousSize=size(thisNode%components)
          call Move_Alloc(thisNode%components,tempComponents)
          allocate(thisNode%components(previousSize+1))
          thisNode%components(1:previousSize)=tempComponents
          deallocate(tempComponents)
          call Memory_Usage_Record(sizeof(thisNode%components(1)),memoryType=memoryTypeNodes,blockCount=0)
       else
          allocate(thisNode%components(1))
          call Memory_Usage_Record(sizeof(thisNode%components),memoryType=memoryTypeNodes)
       end if
       thisNode%componentIndex(componentIndex)=size(thisNode%components)
       thisIndex=size(thisNode%components)
       thisNode%components(thisIndex)%nextComponentOfType => null()
       if (.not.allocated(thisNode%components(thisIndex)%properties)) then
          if (propertyCount > 0) then
             allocate(thisNode%components(thisIndex)%properties(propertyCount,2))
             call Memory_Usage_Record(sizeof(thisNode%components(thisIndex)%properties),memoryType=memoryTypeNodes)
             thisNode%components(thisIndex)%properties=0.0d0
          end if
          if (dataCount > 0) then
             allocate(thisNode%components(thisIndex)%data(dataCount))
             call Memory_Usage_Record(sizeof(thisNode%components(thisIndex)%data),memoryType=memoryTypeNodes)
             thisNode%components(thisIndex)%data=0.0d0
          end if
          if (historyCount > 0) then
             allocate(thisNode%components(thisIndex)%histories(historyCount))
             call Memory_Usage_Record(sizeof(thisNode%components(thisIndex)%histories),memoryType=memoryTypeNodes)
          end if
       end if
    end if

    return
  end subroutine Tree_Node_Allocate_Component

  subroutine Tree_Node_Deallocate_Component(thisNode,componentIndex,propertyCount,dataCount)
    !% Ensure that the component array is deallocated.
    use Memory_Management
    implicit none
    integer,         intent(in)                :: componentIndex,propertyCount,dataCount
#ifdef GCC45
    class(treeNode), intent(inout)             :: thisNode
#else
    type(treeNode),  intent(inout)             :: thisNode
#endif
    type(component), allocatable, dimension(:) :: tempComponents
    integer                                    :: previousSize,listIndex,timesCount,iHistory

    listIndex=thisNode%componentIndex(componentIndex)
    if (listIndex /= -1) then
       ! Count memory that will be freed.
       call Memory_Usage_Record(sizeof(thisNode%components(listIndex)           ),addRemove=-1,memoryType=memoryTypeNodes)
       call Memory_Usage_Record(sizeof(thisNode%components(listIndex)%properties),addRemove=-1,memoryType=memoryTypeNodes)
       call Memory_Usage_Record(sizeof(thisNode%components(listIndex)%data      ),addRemove=-1,memoryType=memoryTypeNodes)
       if (allocated(thisNode%components(listIndex)%histories)) then
          do iHistory=1,size(thisNode%components(listIndex)%histories)
             call thisNode%components(listIndex)%histories(iHistory)%destroy()
          end do
       end if
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
          call Memory_Usage_Record(sizeof(thisNode%components(1)),addRemove=-1,memoryType=memoryTypeNodes,blockCount=0)
       else
          deallocate(thisNode%components)
          call Memory_Usage_Record(sizeof(thisNode%components),addRemove=-1,memoryType=memoryTypeNodes)
       end if
       thisNode%componentIndex(componentIndex)=-1
    end if
    return
  end subroutine Tree_Node_Deallocate_Component

  subroutine Tree_Node_Deallocate_All_Components(thisNode)
    !% Ensure that all components of {\tt thisNode} are deallocated.
    use Memory_Management
    implicit none
#ifdef GCC45
    class(treeNode), intent(inout)             :: thisNode
#else
    type(treeNode),  intent(inout)             :: thisNode
#endif
    integer                                    :: listIndex,thisIndex,iHistory

    ! Deallocate each component and record the memory usage change.
    do thisIndex=1,size(thisNode%componentIndex)
       listIndex=thisNode%componentIndex(thisIndex)
       if (listIndex /= -1) then
          if (allocated(thisNode%components(listIndex)%properties)) then
             call Memory_Usage_Record(sizeof(thisNode%components(listIndex)%properties),addRemove=-1,memoryType=memoryTypeNodes)
             deallocate(thisNode%components(listIndex)%properties)
          end if
          if (allocated(thisNode%components(listIndex)%data)) then
             call Memory_Usage_Record(sizeof(thisNode%components(listIndex)%data),addRemove=-1,memoryType=memoryTypeNodes)
             deallocate(thisNode%components(listIndex)%data)
          end if
          if (allocated(thisNode%components(listIndex)%histories)) then
             do iHistory=1,size(thisNode%components(listIndex)%histories)
                call thisNode%components(listIndex)%histories(iHistory)%destroy()
             end do
          end if
       end if
    end do
    ! Deallocate the components array.
    call Memory_Usage_Record(sizeof(thisNode%components),addRemove=-1,memoryType=memoryTypeNodes)
    deallocate(thisNode%components)
    ! Reset the component index array.
    thisNode%componentIndex=-1
    return
  end subroutine Tree_Node_Deallocate_All_Components


  logical function Tree_Node_Is_Primary_Progenitor(thisNode)
    !% Returns true if {\tt thisNode} is the primary progenitor of its parent node.
    implicit none
#ifdef GCC45
    class(treeNode), intent(inout), target  :: thisNode
#else
    type(treeNode),  intent(inout), pointer :: thisNode
#endif

#ifdef GCC45
    select type(thisNode)
    type is (treeNode)
#endif
       if (associated(thisNode%parentNode)) then
          Tree_Node_Is_Primary_Progenitor=associated(thisNode%parentNode%childNode,thisNode)
       else
          Tree_Node_Is_Primary_Progenitor=.false.
       end if
#ifdef GCC45
    end select
#endif
    return
  end function Tree_Node_Is_Primary_Progenitor

  logical function Tree_Node_Is_Primary_Progenitor_Of_Index(thisNode,targetNodeIndex)
    !% Return true if {\tt thisNode} is a progenitor of the node with index {\tt targetNodeIndex}.
    implicit none
#ifdef GCC45
    class(treeNode),         intent(in), target  :: thisNode
#else
    type(treeNode),          intent(in), pointer :: thisNode
#endif
    integer(kind=kind_int8), intent(in)          :: targetNodeIndex
    type(treeNode),                      pointer :: workNode

    Tree_Node_Is_Primary_Progenitor_Of_Index=.false.
#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       workNode => thisNode
#ifdef GCC45
    end select
#endif
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
#ifdef GCC45
    class(treeNode), intent(in), target  :: thisNode
#else
    type(treeNode),  intent(in), pointer :: thisNode
#endif
    type(treeNode),  intent(in), pointer :: targetNode
    type(treeNode),              pointer :: workNode

    Tree_Node_Is_Primary_Progenitor_Of_Node=.false.
#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       workNode => thisNode
#ifdef GCC45
    end select
#endif
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
#ifdef GCC45
    class(treeNode),         intent(in), target  :: thisNode
#else
    type(treeNode),          intent(in), pointer :: thisNode
#endif
    integer(kind=kind_int8), intent(in)          :: targetNodeIndex
    type(treeNode),                      pointer :: workNode

    Tree_Node_Is_Progenitor_Of_Index=.false.
#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       workNode => thisNode
#ifdef GCC45
    end select
#endif
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
#ifdef GCC45
    class(treeNode), intent(in), target  :: thisNode
#else
    type(treeNode),  intent(in), pointer :: thisNode
#endif
    type(treeNode),  intent(in), pointer :: targetNode
    type(treeNode),              pointer :: workNode

    Tree_Node_Is_Progenitor_Of_Node=.false.
#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       workNode => thisNode
#ifdef GCC45
    end select
#endif
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
#ifdef GCC45
    class(treeNode), intent(inout), target  :: thisNode
#else
    type(treeNode),  intent(inout), pointer :: thisNode
#endif
    type(treeNode),                 pointer :: workNode

    Tree_Node_Is_On_Main_Branch=.not.associated(thisNode%parentNode)
#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       workNode => thisNode
#ifdef GCC45
    end select
#endif
    do while (associated(workNode%parentNode))
       if (.not.workNode%isPrimaryProgenitor()) return
       workNode => workNode%parentNode
    end do
    Tree_Node_Is_On_Main_Branch=.true.  
    return
  end function Tree_Node_Is_On_Main_Branch


  logical function Tree_Node_Is_Satellite(thisNode)
    !% Returns true if {\tt thisNode} is a satellite.
    implicit none
#ifdef GCC45
    class(treeNode), target, intent(in)  :: thisNode
#else
    type(treeNode),  pointer, intent(in) :: thisNode
#endif
    type(treeNode),  pointer             :: parentNode,childNode,thisNodeActual

#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       thisnodeActual => thisNode
#ifdef GCC45
    end select
#endif
    parentNode => thisNodeActual%parentNode
    select case (associated(parentNode))
    case (.false.)
       Tree_Node_Is_Satellite=.false.
       return
    case (.true.)
       childNode => parentNode%childNode
       Tree_Node_Is_Satellite=.true.
       do while (associated(childNode))
          if (associated(childNode,thisNodeActual)) then
             Tree_Node_Is_Satellite=.false.
             exit
          end if
          childNode => childNode%siblingNode
       end do
    end select
    return
  end function Tree_Node_Is_Satellite

  subroutine Tree_Node_Merge_Node(thisNode,mergesWith)
    !% Returns a pointer to the node with which {\tt thisNode} will merge.
    implicit none
#ifdef GCC45
    class(treeNode), target,  intent(in)    :: thisNode
#else
    type(treeNode),  pointer, intent(in)    :: thisNode
#endif
    type(treeNode),  pointer, intent(inout) :: mergesWith
    type(treeNode),  pointer                :: thisNodeActual

#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       thisNodeActual => thisNode
#ifdef GCC45
    end select
#endif

    ! Check if a specific merge node has been set.
    if (associated(thisNode%mergeNode)) then
       ! One has, so simply return it.
       mergesWith => thisNode%mergeNode
    else
       ! No specific merge node has been set, assume merging with the parent node.
       mergesWith => thisNode%parentNode
    end if
    return
  end subroutine Tree_Node_Merge_Node

  subroutine Satellite_Remove_from_Host(satelliteNode)
    !% Remove {\tt satelliteNode} from the linked list of its host node's satellites.
    use Galacticus_Display
    use ISO_Varying_String
    use String_Handling
    implicit none
#ifdef GCC45
    class(treeNode),     target,  intent(in) :: satelliteNode
#else
    type(treeNode),      pointer, intent(in) :: satelliteNode
#endif
    type(treeNode),      pointer             :: hostNode,thisNode,previousNode,satelliteNodeActual
    type(varying_string)                     :: message

#ifdef GCC45
    select type (satelliteNode)
    type is (treeNode)
#endif
       satelliteNodeActual => satelliteNode
#ifdef GCC45
    end select
#endif

    ! Remove from the parent node satellite list.
    hostNode => satelliteNodeActual%parentNode
    message='Satellite node ['
    message=message//satelliteNodeActual%index()//'] being removed from host node ['//hostNode%index()//']'
    call Galacticus_Display_Message(message,verbosityInfo)
    if (associated(hostNode%satelliteNode,satelliteNodeActual)) then
       ! This is the first satellite, unlink it, and link to any sibling.
       hostNode%satelliteNode => satelliteNodeActual%siblingNode
    else
       thisNode     => hostNode%satelliteNode
       previousNode => null()
       do while (associated(thisNode))
          if (associated(thisNode,satelliteNodeActual)) then
             ! Found our node, link its older sibling to its younger sibling.
             previousNode%siblingNode => thisNode%siblingNode
             exit
          end if
          previousNode => thisNode
          thisNode     => thisNode%siblingNode
       end do
    end if
    return
  end subroutine Satellite_Remove_from_Host

  subroutine Satellite_Remove_from_Mergee(mergeeNode)
    !% Remove {\tt mergeeNode} from the linked list of its host node's satellites.
    use Galacticus_Display
    use ISO_Varying_String
    use String_Handling
    implicit none
#ifdef GCC45
    class(treeNode),     target,  intent(in) :: mergeeNode
#else
    type(treeNode),      pointer, intent(in) :: mergeeNode
#endif
    type(treeNode),      pointer             :: hostNode,thisNode,previousNode,mergeeNodeActual
    type(varying_string)                     :: message

#ifdef GCC45
    select type (mergeeNode)
    type is (treeNode)
#endif
       mergeeNodeActual => mergeeNode
#ifdef GCC45
    end select
#endif

    ! Remove from the mergee list of any merge target.
    if (associated(mergeeNodeActual%mergeNode)) then
       hostNode => mergeeNodeActual%mergeNode
       message='Mergee node ['
       message=message//mergeeNodeActual%index()//'] being removed from merge target ['//hostNode%index()//']'
       call Galacticus_Display_Message(message,verbosityInfo)
       if (associated(hostNode%mergeeNode,mergeeNodeActual)) then
          ! This is the first mergee, unlink it, and link to any sibling.
          hostNode%mergeeNode => mergeeNodeActual%nextMergee
       else
          thisNode     => hostNode%mergeeNode
          previousNode => null()
          do while (associated(thisNode))
             if (associated(thisNode,mergeeNodeActual)) then
                ! Found our node, link its older sibling to its younger sibling.
                previousNode%nextMergee => thisNode%nextMergee
                exit
             end if
             previousNode => thisNode
             thisNode     => thisNode%nextMergee
          end do
       end if
    end if
    return
  end subroutine Satellite_Remove_from_Mergee

  subroutine Get_Last_Satellite(thisNode,satelliteNode)
    !% Returns a pointer to the final satellite node associated with {\tt thisNode}.
    implicit none
#ifdef GCC45
    class(treeNode), intent(in)              :: thisNode
#else
    type(treeNode),  intent(in),     pointer :: thisNode
#endif
    type(treeNode),  intent(inout),  pointer :: satelliteNode

    satelliteNode => thisNode%satelliteNode
    do while (associated(satelliteNode%siblingNode))
       satelliteNode => satelliteNode%siblingNode
    end do
    return
  end subroutine Get_Last_Satellite

! <gfortran 4.6> not sure we can support the following since class dummy argument can not be pointer.
!   subroutine Merger_Tree_Walk_Tree_Same_Node(thisNode)
!     !% Simple interface to the \hyperlink{objects.tree_node.F90:tree_nodes:merger_tree_walk_tree}{{\tt
!     !% Merger\_Tree\_Walk\_Tree()}} subroutine that does the walk in place, i.e. by updating the input tree node pointer to point
!     !% to the next node.
!     implicit none
! #ifdef GCC45
!     class (treeNode), target , intent(inout) :: thisNode
! #else
!     type (treeNode),  pointer, intent(inout) :: thisNode
! #endif
!     type (treeNode),  pointer                :: thisNodeActual

! #ifdef GCC45
!     select type (thisNode)
!     type is (treeNode)
! #endif
!        thisNodeActual => thisNode
! #ifdef GCC45
!     end select
! #endif
!     call Merger_Tree_Walk_Tree(thisNode,thisNodeActual)
!     return
!   end subroutine Merger_Tree_Walk_Tree_Same_Node

  subroutine Get_Earliest_Progenitor(thisNode,progenitorNode)
    !% Returns a pointer to the earliest progenitor of with {\tt thisNode}.
    implicit none
#ifdef GCC45
    class(treeNode), target,  intent(inout) :: thisNode
#else
    type(treeNode),  pointer, intent(inout) :: thisNode
#endif
    type(treeNode),  pointer, intent(inout) :: progenitorNode
    type(treeNode),  pointer                :: thisNodeActual


#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       thisNodeActual => thisNode
#ifdef GCC45
    end select
#endif

    progenitorNode => thisNodeActual
    do while (associated(progenitorNode%childNode))
       progenitorNode => progenitorNode%childNode
    end do
    return
  end subroutine Get_Earliest_Progenitor

  subroutine Merger_Tree_Walk_Tree(thisNode,nextNode)
    !% This function provides a mechanism for walking through an entire merger tree. Given a pointer {\tt thisNode}
    !% to a node of the tree, it will return the next node that should be visited in the tree. Thus, if {\tt thisNode} is
    !% initially set to the base of the merger tree and {\tt Merger\_Tree\_Walk()} is called repeatedly it will walk through every node
    !% of the tree. Once the entire tree has been walked, a {\tt null()} pointer will be returned, indicating that there
    !% are no more nodes to walk. Each node will be visited once and once only if the tree is walked in this way.
    implicit none
#ifdef GCC45
    class(treeNode), target,  intent(inout) :: thisNode
#else
    type(treeNode),  pointer, intent(inout) :: thisNode
#endif
    type(treeNode),  pointer, intent(inout) :: nextNode
    type(treeNode),  pointer                :: workNode,thisNodeActual

#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       thisNodeActual => thisNode
#ifdef GCC45
    end select
#endif
    workNode => thisNodeActual

    if (.not.associated(workNode%parentNode)) then
       ! This is the base of the merger tree.
       do while (associated(workNode%childNode))
          workNode => workNode%childNode
       end do
       if (associated(workNode,thisNodeActual)) nullify(workNode)
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

! <gfortran 4.6> not sure we can support the following since class dummy argument can not be pointer.
!   subroutine Merger_Tree_Walk_Tree_With_Satellites_Same_Node(thisNode)
!     !% Simple interface to the \hyperlink{objects.tree_node.F90:tree_nodes:merger_tree_walk_tree_with_satellites}{{\tt
!     !% Merger\_Tree\_Walk\_Tree\_With \_Satelites()}} subroutine that does the walk in place, i.e. by updating the input tree node
!     !% pointer to point to the next node.
!     implicit none
! #ifdef GCC45
!     class (treeNode), intent(inout)          :: thisNode
! #else
!     type (treeNode),  intent(inout), pointer :: thisNode
! #endif
!     type (treeNode),                 pointer :: thisNodeActual

! #ifdef GCC45
!     select type (thisNode)
!     type is (treeNode)
! #endif
!        thisNodeActual => thisNode
! #ifdef GCC45
!     end select
! #endif
!     call Merger_Tree_Walk_Tree_With_Satellites(thisNode,thisNodeActual)
!     return
!   end subroutine Merger_Tree_Walk_Tree_With_Satellites_Same_Node

  subroutine Merger_Tree_Walk_Tree_With_Satellites(thisNode,nextNode)
    !% Merger tree walk function which also descends through satellite nodes. Note that it is important that the walk descends to
    !% satellites before descending to children: the routines that destroy merger tree branches rely on this since child nodes are
    !% used in testing whether a node is a satellite---if they are destroyed prior to the test being made then problems with
    !% dangling pointers will occur.
    implicit none
#ifdef GCC45
    class (treeNode), target , intent(inout) :: thisNode
#else
    type (treeNode),  pointer, intent(inout) :: thisNode
#endif
    type (treeNode),  pointer, intent(inout) :: nextNode
    type (treeNode),  pointer                :: workNode,thisNodeActual

#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       thisNodeActual => thisNode
#ifdef GCC45
    end select
#endif
    workNode => thisNodeActual
    if (.not.associated(workNode%parentNode)) then
       ! This is the base of the merger tree.
       ! Descend through satellites and children.
       workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       if (associated(workNode,thisNodeActual)) nullify(workNode)
    else
       if (associated(workNode%siblingNode)) then
          workNode => workNode%siblingNode
          ! Descend through satellites and children.
          workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       else
          ! About to move back up the tree. Check if the node we're moving up from is a satellite.
          if (workNode%isSatellite()) then
             ! It is a satellite. Therefore, the parent may have children that have yet to be visited. Check if the parent has
             ! children.
             if (associated(workNode%parentNode%childNode)) then
                ! Parent does have children, so move to the first one.
                workNode => workNode%parentNode%childNode
                ! Descend through satellites and children.
                workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
             else
                ! Parent has no children, so move to the parent.
                workNode => workNode%parentNode
             end if
          else
             ! It is not a satellite, so all satellites and children have been processed.
             workNode => workNode%parentNode
          end if
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
#ifdef GCC45
    class (treeNode), target,  intent(inout) :: thisNode
#else
    type (treeNode),  pointer, intent(inout) :: thisNode
#endif
    type (treeNode),  pointer, intent(inout) :: startNode,nextNode
    type (treeNode),  pointer                :: workNode,thisNodeActual

#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       thisNodeActual => thisNode
#ifdef GCC45
    end select
#endif
    workNode => thisNodeActual

    if (associated(thisNodeActual,startNode)) then
       do while (associated(workNode%childNode))
          workNode => workNode%childNode
       end do
       if (associated(workNode,thisNodeActual)) nullify(workNode)
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

  subroutine Merger_Tree_Walk_Branch_With_Satellites(thisNode,startNode,nextNode)
    !% This function provides a mechanism for walking through the branches of the merger tree. Given a pointer {\tt thisNode} to a
    !% branch of the tree, it will return the next node that should be visited in the tree. Thus, if {\tt thisNode} is initially
    !% set to the base of the merger tree and {\tt Merger\_Tree\_Walk\_Branch()} is called repeatedly it will walk through every
    !% node of the branch. Once the entire branch has been walked, a {\tt null()} pointer will be returned, indicating that there
    !% are no more nodes to walk. Each node will be visited once and once only if the branch is walked in this way. Note that it
    !% is important that the walk descends to satellites before descending to children: the routines that destroy merger tree
    !% branches rely on this since child nodes are used in testing whether a node is a satellite---if they are destroyed prior to
    !% the test being made then problems with dangling pointers will occur.
    implicit none
#ifdef GCC45
    class (treeNode), target,  intent(inout) :: thisNode
#else
    type (treeNode),  pointer, intent(inout) :: thisNode
#endif
    type (treeNode),  pointer, intent(inout) :: startNode,nextNode
    type (treeNode),  pointer                :: workNode,thisNodeActual

#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       thisNodeActual => thisNode
#ifdef GCC45
    end select
#endif
    workNode => thisNodeActual

    if (associated(thisNodeActual,startNode)) then
       ! Descend through satellites and children.
       workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       if (associated(workNode,thisNodeActual)) nullify(workNode)
    else
       if (associated(workNode%siblingNode)) then
          workNode => workNode%siblingNode
          ! Descend through satellites and children.
          workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       else
          ! About to move back up the tree. Check if the node we're moving up from is a satellite.
          if (workNode%isSatellite()) then
             ! It is a satellite. Therefore, the parent may have children that have yet to be visited. Check if the parent
             ! has children.
             if (associated(workNode%parentNode%childNode)) then
                ! Parent does have children, so move to the first one.
                workNode => workNode%parentNode%childNode
                ! Descend through satellites and children.
                workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
             else
                ! Parent has no satellites, so move to the parent.
                workNode => workNode%parentNode
             end if
          else
             ! It is not a satellite, so all satellites and children of the parent must have been processed. Therefore, move to
             ! the parent.
             workNode => workNode%parentNode
          end if
          if (associated(workNode,startNode)) workNode => null() ! Terminate when back at starting node.
       end if
    end if
    nextNode => workNode
    return
  end subroutine Merger_Tree_Walk_Branch_With_Satellites

! <gfortran 4.6> not sure we can support the following since class dummy argument can not be pointer.
!   subroutine Merger_Tree_Construction_Walk_Same_Node(thisNode)
!     !% Simple interface to the \hyperlink{objects.tree_node.F90:tree_nodes:merger_tree_construction_walk}{{\tt Merger\_Tree\_Construction\_Walk()}} subroutine
!     !% that does the walk in place, i.e. by updating the input tree node pointer to point to the next node.
!     implicit none
! #ifdef GCC45
!     class (treeNode), target,  intent(inout) :: thisNode
! #else
!     type (treeNode),  pointer, intent(inout) :: thisNode
! #endif
!     type (treeNode),  pointer                :: thisNodeActual

! #ifdef GCC45
!     select type (thisNode)
!     type is (treeNode)
! #endif
!        thisNodeActual => thisNode
! #ifdef GCC45
!     end select
! #endif

!     call Merger_Tree_Construction_Walk(thisNode,thisNodeActual)
!     return
!   end subroutine Merger_Tree_Construction_Walk_Same_Node

  function Merger_Tree_Walk_Descend_to_Progenitors(thisNode) result (progenitorNode)
    !% Descend to the deepest progenitor (satellites and children) of {\tt thisNode}.
    implicit none
    type(treeNode), intent(in), pointer :: thisNode
    type(treeNode),             pointer :: progenitorNode

    ! Begin at the input node.
    progenitorNode => thisNode
    
    ! Descend through satellites and children.
    do while (associated(progenitorNode%satelliteNode).or.associated(progenitorNode%childNode))
       ! Descend through any satellite nodes.
       do while (associated(progenitorNode%satelliteNode))
          progenitorNode => progenitorNode%satelliteNode
       end do
       ! Descend through any child nodes.
       do while (associated(progenitorNode%childNode))
          progenitorNode => progenitorNode%childNode
       end do
    end do
    return
  end function Merger_Tree_Walk_Descend_to_Progenitors

  subroutine Merger_Tree_Construction_Walk(thisNode,nextNode)
    !% This function provides a mechanism for walking through a merger tree that is being built.
    implicit none
#ifdef GCC45
    class (treeNode), target,  intent(inout) :: thisNode
#else
    type (treeNode),  pointer, intent(inout) :: thisNode
#endif
    type (treeNode),  pointer, intent(inout) :: nextNode
    type (treeNode),  pointer                :: workNode

#ifdef GCC45
    select type (thisNode)
    type is (treeNode)
#endif
       workNode => thisNode
#ifdef GCC45
    end select
#endif

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

  subroutine Tree_Node_Rate_Rate_Compute_Dummy(thisNode,interrupt,interruptProcedure)
    !% A dummy rate compute subroutine that can be pointed to by any sterile methods.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure

    return
  end subroutine Tree_Node_Rate_Rate_Compute_Dummy

  subroutine Tree_Node_Rate_Adjust_Dummy(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% A dummy rate-adjust subroutine that can be pointed to by unconnected pipes.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment

    return
  end subroutine Tree_Node_Rate_Adjust_Dummy

  subroutine Tree_Node_Rate_Adjust_Array_Dummy(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% A dummy rate-adjust subroutine that can be pointed to by unconnected pipes.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment(:)

    return
  end subroutine Tree_Node_Rate_Adjust_Array_Dummy

  subroutine Tree_Node_Rate_Adjust_History_Dummy(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% A dummy rate-adjust subroutine that can be pointed to by unconnected pipes.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure
    type(history),           intent(in)    :: rateAdjustment

    return
  end subroutine Tree_Node_Rate_Adjust_History_Dummy

end module Tree_Nodes
