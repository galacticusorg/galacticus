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

!% Contains a module which implements a ``superwind'' expulsive outflow rate (as in \citep{baugh_can_2005}) due to star
!% formation feedback in galactic disks.

module Star_Formation_Expulsive_Feedback_Disks_Superwind
  !% Implements a``superwind'' outflow rate (as in \citep{baugh_can_2005}) due to star formation feedback in
  !% galactic disks.
  use Tree_Nodes
  implicit none
  private
  public :: Star_Formation_Expulsive_Feedback_Disks_SW_Initialize

  ! Parameters of the feedback model.
  double precision :: diskSuperwindMassLoading,diskSuperwindVelocity
  
contains

  !# <starFormationExpulsiveFeedbackDisksMethod>
  !#  <unitName>Star_Formation_Expulsive_Feedback_Disks_SW_Initialize</unitName>
  !# </starFormationExpulsiveFeedbackDisksMethod>
  subroutine Star_Formation_Expulsive_Feedback_Disks_SW_Initialize(starFormationExpulsiveFeedbackDisksMethod&
       &,Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_Get)
    !% Initializes the ``superwind'' disk star formation expulsive feedback module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: starFormationExpulsiveFeedbackDisksMethod
    procedure(double precision), pointer, intent(inout) :: Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_Get
    
    if (starFormationExpulsiveFeedbackDisksMethod == 'superwind') then
       Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_Get => Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_SW
       ! Get parameters of for the feedback calculation.
       !@ <inputParameter>
       !@   <name>diskSuperwindMassLoading</name>
       !@   <defaultValue>2</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass loading of the disk superwind.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('diskSuperwindMassLoading',diskSuperwindMassLoading,defaultValue=2.0d0)
       !@ <inputParameter>
       !@   <name>diskSuperwindVelocity</name>
       !@   <defaultValue>200 km/s</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The velocity scale of the disk superwind.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('diskSuperwindVelocity',diskSuperwindVelocity,defaultValue=200.0d0)
    end if
    return
  end subroutine Star_Formation_Expulsive_Feedback_Disks_SW_Initialize

  double precision function Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_SW(thisNode,starFormationRate,energyInputRate)
    !% Returns the expulsive outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic disk of {\tt thisNode}. The outflow
    !% rate is given by
    !% \begin{equation}
    !% \dot{M}_{\rm outflow} = f_{\rm SW,0} \left\{ \begin{array}{ll} 1 & \hbox{ if } V_{\rm disk} < V_{\rm disk,SW} \\ (V_{\rm disk,SW}/V_{\rm disk})^2 &  \hbox{ if } V_{\rm disk} \ge V_{\rm disk,SW} \end{array} \right. ,
    !% \end{equation}
    !%  where $V_{\rm disk,SW}=${\tt
    !% [diskSuperwindVelocity]} and $f_{\rm SW,0}=${\tt [diskSuperwindMassLoading]}. Note that the velocity $V_{\rm
    !% disk}$ is whatever characteristic value returned by the disk method. This scaling is functionally similar to
    !% that adopted by \cite{cole_hierarchical_2000} and \cite{baugh_can_2005}, except that they specifically used the
    !% circular velocity at half-mass radius.
    use Tree_Nodes
    use Numerical_Constants_Units
    use Stellar_Feedback
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: starFormationRate,energyInputRate
    double precision                         :: diskVelocity,outflowRateToStarFormationRate

    ! Get disk circular velocity.
    diskVelocity=Tree_Node_Disk_Velocity(thisNode)

    ! Check for zero velocity disk.
    if (diskVelocity <= 0.0d0) then
       Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_SW=0.0d0 ! No well defined answer in this case.
    else
       outflowRateToStarFormationRate=diskSuperwindMassLoading*(diskSuperwindVelocity/max(diskVelocity,diskSuperwindVelocity))**2
       Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_SW=outflowRateToStarFormationRate*energyInputRate &
            &/feedbackEnergyInputAtInfinityCanonical
    end if
    return
  end function Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate_SW
  
end module Star_Formation_Expulsive_Feedback_Disks_Superwind
