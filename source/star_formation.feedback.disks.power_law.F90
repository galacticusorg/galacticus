!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a power-law outflow rate due to star formation feedback in galactic disks.

module Star_Formation_Feedback_Disks_Power_Law
  !% Implements a power-law outflow rate due to star formation feedback in galactic disks.
  implicit none
  private
  public :: Star_Formation_Feedback_Disks_Power_Law_Initialize

  ! Parameters of the feedback model.
  double precision :: diskOutflowExponent, diskOutflowVelocity

contains

  !# <starFormationFeedbackDisksMethod>
  !#  <unitName>Star_Formation_Feedback_Disks_Power_Law_Initialize</unitName>
  !# </starFormationFeedbackDisksMethod>
  subroutine Star_Formation_Feedback_Disks_Power_Law_Initialize(starFormationFeedbackDisksMethod,Star_Formation_Feedback_Disk_Outflow_Rate_Get)
    !% Initializes the ``power law'' disk star formation feedback module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                                     ), intent(in   )          :: starFormationFeedbackDisksMethod
    procedure(Star_Formation_Feedback_Disk_Outflow_Rate_Power_Law), intent(inout), pointer :: Star_Formation_Feedback_Disk_Outflow_Rate_Get

    if (starFormationFeedbackDisksMethod == 'powerLaw') then
       Star_Formation_Feedback_Disk_Outflow_Rate_Get => Star_Formation_Feedback_Disk_Outflow_Rate_Power_Law
       ! Get parameters of for the feedback calculation.
       !@ <inputParameter>
       !@   <name>diskOutflowVelocity</name>
       !@   <defaultValue>250</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The velocity scale at which the \gls{sne}-driven outflow rate equals the star formation rate in disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('diskOutflowVelocity',diskOutflowVelocity,defaultValue=250.0d0)
       !@ <inputParameter>
       !@   <name>diskOutflowExponent</name>
       !@   <defaultValue>3.5</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The velocity scaling of the \gls{sne}-driven outflow rate in disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('diskOutflowExponent',diskOutflowExponent,defaultValue=  3.5d0)
    end if
    return
  end subroutine Star_Formation_Feedback_Disks_Power_Law_Initialize

  double precision function Star_Formation_Feedback_Disk_Outflow_Rate_Power_Law(thisNode,starFormationRate,energyInputRate)
    !% Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic disk of {\tt thisNode}. The outflow
    !% rate is given by
    !% \begin{equation}
    !% \dot{M}_{\rm outflow} = \left({V_{\rm disk,outflow} \over V_{\rm disk}}\right)^{\alpha_{\rm disk,outflow}},
    !% \end{equation}
    !% where $V_{\rm disk,outflow}$(={\tt diskOutflowVelocity}) is the velocity scale at which outflow rate equals star formation
    !% rate and $\alpha_{\rm disk,outflow}$(={\tt diskOutflowExponent}) controls the scaling with velocity. Note that the velocity
    !% $V_{\rm disk}$ is whatever characteristic value returned by the disk method. This scaling is functionally similar to that
    !% adopted by \cite{cole_hierarchical_2000}, but that they specifically used the circular velocity at half-mass radius.
    use Galacticus_Nodes
    use Stellar_Feedback
    implicit none
    type            (treeNode         ), intent(inout), pointer :: thisNode
    class           (nodeComponentDisk)               , pointer :: thisDiskComponent
    double precision                   , intent(in   )          :: energyInputRate  , starFormationRate
    double precision                                            :: diskVelocity     , outflowRateToStarFormationRate

    ! Get the disk.
    thisDiskComponent => thisNode%disk()

    ! Get disk circular velocity.
    diskVelocity=thisDiskComponent%velocity()

    ! Check for zero velocity disk.
    if (diskVelocity <= 0.0d0) then
       Star_Formation_Feedback_Disk_Outflow_Rate_Power_Law=0.0d0 ! No well defined answer in this case.
    else
       outflowRateToStarFormationRate=(diskOutflowVelocity/diskVelocity)**diskOutflowExponent
       Star_Formation_Feedback_Disk_Outflow_Rate_Power_Law=outflowRateToStarFormationRate*energyInputRate&
            &/feedbackEnergyInputAtInfinityCanonical
    end if
    return
  end function Star_Formation_Feedback_Disk_Outflow_Rate_Power_Law

end module Star_Formation_Feedback_Disks_Power_Law
