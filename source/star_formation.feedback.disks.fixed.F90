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

!% Contains a module which implements a fixed outflow rate due to star formation feedback in galactic disks.

module Star_Formation_Feedback_Disks_Fixed
  !% Implements a fixed outflow rate due to star formation feedback in galactic disks.
  implicit none
  private
  public :: Star_Formation_Feedback_Disks_Fixed_Initialize

  ! Parameters of the feedback model.
  double precision :: diskOutflowFraction
  
contains

  !# <starFormationFeedbackDisksMethod>
  !#  <unitName>Star_Formation_Feedback_Disks_Fixed_Initialize</unitName>
  !# </starFormationFeedbackDisksMethod>
  subroutine Star_Formation_Feedback_Disks_Fixed_Initialize(starFormationFeedbackDisksMethod,Star_Formation_Feedback_Disk_Outflow_Rate_Get)
    !% Initializes the ``fixed'' disk star formation feedback module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in   ) :: starFormationFeedbackDisksMethod
    procedure(double precision), pointer, intent(inout) :: Star_Formation_Feedback_Disk_Outflow_Rate_Get
    
    if (starFormationFeedbackDisksMethod == 'fixed') then
       Star_Formation_Feedback_Disk_Outflow_Rate_Get => Star_Formation_Feedback_Disk_Outflow_Rate_Fixed
       ! Get parameters of for the feedback calculation.
       !@ <inputParameter>
       !@   <name>diskOutflowFraction</name>
       !@   <defaultValue>0.01</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The ratio of outflow rate to star formation rate in disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('diskOutflowFraction',diskOutflowFraction,defaultValue=0.01d0)
    end if
    return
  end subroutine Star_Formation_Feedback_Disks_Fixed_Initialize

  double precision function Star_Formation_Feedback_Disk_Outflow_Rate_Fixed(thisNode,starFormationRate,energyInputRate)
    !% Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic disk of {\tt thisNode}. Assumes a
    !% fixed ratio of outflow rate to star formation rate.
    use Tree_Nodes
    use Stellar_Feedback
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: starFormationRate,energyInputRate

    Star_Formation_Feedback_Disk_Outflow_Rate_Fixed=diskOutflowFraction*energyInputRate/feedbackEnergyInputAtInfinityCanonical
    return
  end function Star_Formation_Feedback_Disk_Outflow_Rate_Fixed
  
end module Star_Formation_Feedback_Disks_Fixed
