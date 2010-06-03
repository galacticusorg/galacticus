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






!% Contains a module which implements methods for the tree node object.


module Tree_Node_Methods
  !% Implements methods for the tree node object.
  use Components
  use Tree_Nodes
  use Histories
  private
  public :: Tree_Node_Rate_Rate_Compute_Dummy, Tree_Node_Rate_Adjust_Dummy,&
       & Tree_Node_Rate_Adjust_Array_Dummy, Tree_Node_Rate_Adjust_History_Dummy, Get_Template, Set_Template

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

end module Tree_Node_Methods
