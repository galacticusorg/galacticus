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
