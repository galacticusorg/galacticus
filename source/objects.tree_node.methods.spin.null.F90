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






!% A null spin method.

module Tree_Node_Methods_Spin_Null
  !% A null spin method.
  private
  public :: Tree_Node_Methods_Spin_Null_Initialize
  
contains


  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Spin_Null_Initialize</unitName>
  !#  <optionName>treeNodeMethodSpin</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Spin_Null_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node spin null methods module.
    use Tree_Node_Methods
    use ISO_Varying_String
    use Galacticus_Display
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption == 'null') then

       ! Display message.
       message='Null spin method selected'
       call Galacticus_Display_Message(message,verbosityInfo)

    end if
    return
  end subroutine Tree_Node_Methods_Spin_Null_Initialize
  
end module Tree_Node_Methods_Spin_Null
