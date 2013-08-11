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

!% Contains a module which implements a null method for star formation histories.

module Star_Formation_Histories_Null
  !% Implements a null method for star formation histories.
  implicit none
  private
  public :: Star_Formation_Histories_Null_Initialize

contains

  !# <starFormationHistoriesMethod>
  !#  <unitName>Star_Formation_Histories_Null_Initialize</unitName>
  !# </starFormationHistoriesMethod>
  subroutine Star_Formation_Histories_Null_Initialize(starFormationHistoriesMethod,Star_Formation_History_Create_Do&
       &,Star_Formation_History_Scales_Do,Star_Formation_History_Record_Do,Star_Formation_History_Output_Do)
    !% Initializes the metallicity split star formation history module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                    ), intent(in   )          :: starFormationHistoriesMethod
    procedure(Star_Formation_History_Create_Null), intent(inout), pointer :: Star_Formation_History_Create_Do
    procedure(Star_Formation_History_Scales_Null), intent(inout), pointer :: Star_Formation_History_Scales_Do
    procedure(Star_Formation_History_Record_Null), intent(inout), pointer :: Star_Formation_History_Record_Do
    procedure(Star_Formation_History_Output_Null), intent(inout), pointer :: Star_Formation_History_Output_Do

    if (starFormationHistoriesMethod == 'null') then
       ! Associate procedure pointers.
       Star_Formation_History_Create_Do => Star_Formation_History_Create_Null
       Star_Formation_History_Scales_Do => Star_Formation_History_Scales_Null
       Star_Formation_History_Record_Do => Star_Formation_History_Record_Null
       Star_Formation_History_Output_Do => Star_Formation_History_Output_Null
    end if
    return
  end subroutine Star_Formation_Histories_Null_Initialize

  subroutine Star_Formation_History_Create_Null(thisNode,thisHistory)
    !% Create the history required for storing star formation history.
    use Histories
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    type(history ), intent(inout)          :: thisHistory

    ! Do nothing.
    return
  end subroutine Star_Formation_History_Create_Null

  subroutine Star_Formation_History_Record_Null(thisNode,thisHistory,fuelAbundances,starFormationRate)
    !% Record the star formation history for {\tt thisNode}.
    use Histories
    use Galacticus_Nodes
    use Abundances_Structure
    implicit none
    type            (treeNode  ), intent(inout), pointer :: thisNode
    type            (history   ), intent(inout)          :: thisHistory
    type            (abundances), intent(in   )          :: fuelAbundances
    double precision            , intent(in   )          :: starFormationRate

    ! Ensure the history does not exist.
    call thisHistory%destroy()
    return
  end subroutine Star_Formation_History_Record_Null

  subroutine Star_Formation_History_Output_Null(thisNode,nodePassesFilter,thisHistory,iOutput,treeIndex,componentLabel)
    !% Output the star formation history for {\tt thisNode}.
    use Histories
    use Galacticus_Nodes
    use Kind_Numbers
    implicit none
    type     (treeNode      ), intent(inout), pointer :: thisNode
    logical                  , intent(in   )          :: nodePassesFilter
    type     (history       ), intent(inout)          :: thisHistory
    integer                  , intent(in   )          :: iOutput
    integer  (kind=kind_int8), intent(in   )          :: treeIndex
    character(len=*         ), intent(in   )          :: componentLabel

    ! Do nothing.
    return
  end subroutine Star_Formation_History_Output_Null

  subroutine Star_Formation_History_Scales_Null(thisHistory,stellarMass,stellarAbundances)
    !% Set the scalings for error control on the absolute values of star formation histories.
    use Histories
    use Abundances_Structure
    implicit none
    double precision            , intent(in   ) :: stellarMass
    type            (abundances), intent(in   ) :: stellarAbundances
    type            (history   ), intent(inout) :: thisHistory

    ! Do nothing.
    return
  end subroutine Star_Formation_History_Scales_Null

end module Star_Formation_Histories_Null
