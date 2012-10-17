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

!% Contains a module which implements a null expulsive outflow rate in galactic disks.

module Star_Formation_Expulsive_Feedback_Disks_Null
  !% Implementss a null expulsive outflow rate in galactic disks.
  use Tree_Nodes
  implicit none
  private
  public :: Star_Formation_Expulsive_Feedback_Disks_Null_Initialize
  
contains

  !# <starFormationExpulsiveFeedbackDisksMethod>
  !#  <unitName>Star_Formation_Expulsive_Feedback_Disks_Null_Initialize</unitName>
  !# </starFormationExpulsiveFeedbackDisksMethod>
  subroutine Star_Formation_Expulsive_Feedback_Disks_Null_Initialize(starFormationExpulsiveFeedbackDisksMethod&
       &,Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_Get)
    !% Initializes the ``null'' disk star formation expulsive feedback module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: starFormationExpulsiveFeedbackDisksMethod
    procedure(double precision), pointer, intent(inout) :: Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_Get
    
    if (starFormationExpulsiveFeedbackDisksMethod == 'null') Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_Get => Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_Null
     
    return
  end subroutine Star_Formation_Expulsive_Feedback_Disks_Null_Initialize

  double precision function Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_Null(thisNode,starFormationRate,energyInputRate)
    !% Implements a null expulsive outflow rate for disks.
    use Tree_Nodes
    use Numerical_Constants_Units
    use Stellar_Feedback
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: starFormationRate,energyInputRate

    ! Return a zero outflow rate.
    Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_Null=0.0d0
    return
  end function Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_Null
  
end module Star_Formation_Expulsive_Feedback_Disks_Null
