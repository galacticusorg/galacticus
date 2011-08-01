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
    use Input_Parameters
    use Numerical_Ranges
    use Memory_Management
    implicit none
    type(varying_string),          intent(in)    :: starFormationHistoriesMethod
    procedure(),          pointer, intent(inout) :: Star_Formation_History_Create_Do,Star_Formation_History_Scales_Do&
         &,Star_Formation_History_Record_Do,Star_Formation_History_Output_Do
    
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
    use Numerical_Ranges
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    type(history),    intent(inout)          :: thisHistory

    ! Do nothing.
    return
  end subroutine Star_Formation_History_Create_Null

  subroutine Star_Formation_History_Record_Null(thisNode,thisHistory,fuelAbundances,starFormationRate)
    !% Record the star formation history for {\tt thisNode}.
    use Histories
    use Numerical_Ranges
    use Tree_Nodes
    use Abundances_Structure
    use Arrays_Search
    implicit none
    type(treeNode),            intent(inout), pointer :: thisNode
    type(history),             intent(inout)          :: thisHistory
    type(abundancesStructure), intent(in)             :: fuelAbundances
    double precision,          intent(in)             :: starFormationRate

    ! Do nothing.
    return
  end subroutine Star_Formation_History_Record_Null

  subroutine Star_Formation_History_Output_Null(thisNode,nodePassesFilter,thisHistory,iOutput,treeIndex,componentLabel)
    !% Output the star formation history for {\tt thisNode}.
    use Histories
    use ISO_Varying_String
    use Galacticus_HDF5
    use IO_HDF5
    use Tree_Nodes
    use String_Handling
    use Kind_Numbers
    implicit none
    type(treeNode),          intent(inout), pointer :: thisNode
    logical,                 intent(in)             :: nodePassesFilter
    type(history),           intent(inout)          :: thisHistory
    integer,                 intent(in)             :: iOutput
    integer(kind=kind_int8), intent(in)             :: treeIndex
    character(len=*),        intent(in)             :: componentLabel
 
    ! Do nothing.
    return
  end subroutine Star_Formation_History_Output_Null

  subroutine Star_Formation_History_Scales_Null(thisHistory,stellarMass,stellarAbundances)
    !% Set the scalings for error control on the absolute values of star formation histories.
    use Histories
    use Stellar_Feedback
    use Abundances_Structure
    use Memory_Management
    implicit none
    double precision,          intent(in)    :: stellarMass
    type(abundancesStructure), intent(in)    :: stellarAbundances
    type(history),             intent(inout) :: thisHistory

    ! Do nothing.
    return
  end subroutine Star_Formation_History_Scales_Null

end module Star_Formation_Histories_Null
