!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module of node position methods in which properties are null.

module Tree_Node_Methods_Positions_Null
  !% Implements node positions using null data.
  use Tree_Nodes
  use Components
  implicit none
  private
  public :: Tree_Node_Methods_Position_Initialize_Null
  
  ! Define procedure pointers.
  !# <treeNodeMethodsPointer type="array">
  !#  <methodName>Tree_Node_Position</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="array">
  !#  <methodName>Tree_Node_Velocity</methodName>
  !# </treeNodeMethodsPointer>

contains

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Position_Initialize_Null</unitName>
  !#  <optionName default="null">treeNodeMethodPosition</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Position_Initialize_Null(componentOption,componentTypeCount)
    !% Initializes the tree node satellite orbit methods module.
    use ISO_Varying_String
    use Galacticus_Display
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount

    ! Check if this implementation is selected.
    if (componentOption == 'null') then

       ! Display message.
       call Galacticus_Display_Message('Null position method selected',verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Position              => null()
       Tree_Node_Position_Set          => null()
       Tree_Node_Position_Rate_Adjust  => null()
       Tree_Node_Position_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy
       Tree_Node_Velocity              => null()
       Tree_Node_Velocity_Set          => null()
       Tree_Node_Velocity_Rate_Adjust  => null()
       Tree_Node_Velocity_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy
    end if

    return
  end subroutine Tree_Node_Methods_Position_Initialize_Null

end module Tree_Node_Methods_Positions_Null
