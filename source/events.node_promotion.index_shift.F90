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

!% Contains a module which optionally shifts the index of a node about to be promoted to its parent node, allowing indices to be
!% tracked along merger trees.

module Node_Promotion_Index_Shifts
  !% Implements optional shifting of the index of a node about to be promoted to its parent node, allowing indices to be
  !% tracked along merger trees.
  implicit none
  private
  public :: Node_Promotion_Index_Shift

  ! Flag indicating if this module has been initialized.
  logical :: indexShiftInitialized  =.false.

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
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    type(treeNode)               , pointer :: parentNode

    ! Ensure that the module is initialized.
    if (.not.indexShiftInitialized) then
       !$omp critical (Node_Promotion_Index_Shift_Initialize)
       if (.not.indexShiftInitialized) then
          !@ <inputParameter>
          !@   <name>nodePromotionIndexShift</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Specifies whether or not the index of a node should be shifted to its parent node prior to promotion.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('nodePromotionIndexShift',nodePromotionIndexShift,defaultValue=.false.)
          ! Flag that the module is now initialized.
          indexShiftInitialized=.true.
       end if
       !$omp end critical (Node_Promotion_Index_Shift_Initialize)
    end if

    ! Check if we are to perform an index shift.
    if (nodePromotionIndexShift) then
       ! Get the parent node.
       parentNode => thisNode%parent
       ! Shift the index from thisNode to the parent node.
       call parentNode%indexSet(thisNode%index())
    end if

    return
  end subroutine Node_Promotion_Index_Shift

end module Node_Promotion_Index_Shifts
