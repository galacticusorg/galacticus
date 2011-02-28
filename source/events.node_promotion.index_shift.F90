!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which optionally shifts the index of a node about to be promoted to its parent node, allowing indices to be
!% tracked along merger trees.

module Node_Promotion_Index_Shifts
  !% Implements optional shifting of the index of a node about to be promoted to its parent node, allowing indices to be
  !% tracked along merger trees.
  private
  public :: Node_Promotion_Index_Shift
  
  ! Flag indicating if this module has been initialized.
  logical :: indexShiftInitialized=.false.
  
  ! Flag indicating whether or not to shift indices.
  logical :: nodePromotionIndexShift

contains
  
  !# <nodePromotionTask>
  !#  <unitName>Node_Promotion_Index_Shift</unitName>
  !# </nodePromotionTask>
  subroutine Node_Promotion_Index_Shift(thisNode)
    !% Shifts the index of {\tt thisNode} to its parent node just prior to promotion, thereby allowing indices to track galaxies
    !% through the tree.
    use Input_Parameters
    use Tree_Nodes
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(treeNode), pointer                :: parentNode

    ! Ensure that the module is initialized.
    !$omp critical (Node_Promotion_Index_Shift_Initialize)
    if (.not.indexShiftInitialized) then
       !@ <inputParameter>
       !@   <name>nodePromotionIndexShift</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether or not the index of a node should be shifted to its parent node prior to promotion.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('nodePromotionIndexShift',nodePromotionIndexShift,defaultValue=.false.)
       ! Flag that the module is now initialized.
       indexShiftInitialized=.true.
    end if
    !$omp end critical (Node_Promotion_Index_Shift_Initialize)
    
    ! Check if we are to perform an index shift.
    if (nodePromotionIndexShift) then
       ! Get the parent node.
       parentNode => thisNode%parentNode
       ! Shift the index from thisNode to the parent node.
       call parentNode%indexSet(thisNode%index())
    end if

    return
  end subroutine Node_Promotion_Index_Shift
  
end module Node_Promotion_Index_Shifts
